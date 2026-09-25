!> @file volume_constraint.f90
!! @copyright
!! Copyright (c) 2025, The Neko-TOP Authors
!! All rights reserved.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions
!! are met:
!!
!!   * Redistributions of source code must retain the above copyright
!!     notice, this list of conditions and the following disclaimer.
!!
!!   * Redistributions in binary form must reproduce the above
!!     copyright notice, this list of conditions and the following
!!     disclaimer in the documentation and/or other materials provided
!!     with the distribution.
!!
!!   * Neither the name of the authors nor the names of its
!!     contributors may be used to endorse or promote products derived
!!     from this software without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
!! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
!! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
!! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
!! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
!! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
!! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!! POSSIBILITY OF SUCH DAMAGE.
!
!> Implements the `volume_constraint_t` type.
!! \f$V = \frac{1}{|\Omega_O|}\int_{\Omega_O} \tilde{\rho} d\Omega\f$
!! Either
!! \f$V < V_\text{max} \f$
!! \f$V > V_\text{min} \f$
module volume_constraint
  use constraint, only: constraint_t

  use design, only: design_t
  use brinkman_design, only: brinkman_design_t
  use simulation_m, only: simulation_t
  use mapping_handler, only: mapping_handler_t
  use num_types, only: rp
  use coefs, only: coef_t
  use json_module, only: json_file
  use json_utils, only: json_get, json_get_or_default
  use field, only: field_t
  use scratch_registry, only: neko_scratch_registry
  use neko_config, only: NEKO_BCKND_DEVICE
  use mask_ops, only: mask_exterior_const
  use math, only: glsc2, copy, cmult
  use device_math, only: device_glsc2, device_copy, device_cmult
  use vector_math, only: vector_cmult, vector_cfill
  use math_ext, only: glsc2_mask
  use field_math, only: field_rone, field_copy, field_cmult, field_cfill
  use utils, only: neko_error
  use vector, only: vector_t
  use neko_ext, only: field_to_vector
  implicit none
  private

  !> A constraint on the volume of the design.
  type, public, extends(constraint_t) :: volume_constraint_t
     private

     !> whether it is minimum or maximum volume
     !! is_max = .false.,   ie V > V_min      =>     -V + V_max < 0
     !! is_max = .true. ,   ie V < V_max      =>      V - V_max < 0
     logical :: is_max
     !> Maximum (or minimum) volume
     real(kind=rp) :: limit
     !> Volume of the optimization domain
     real(kind=rp) :: volume_domain
     !> Pointer to the SEM field.
     class(coef_t), pointer :: c_Xh => null()
     !> Mapping cascade
     type(mapping_handler_t) :: mapping
     !> if mapping is needed
     logical :: if_mapping = .false.

   contains

     !> The common constructor using a JSON object.
     procedure, public, pass(this) :: init_json_sim => &
          volume_constraint_init_json_sim
     !> The direct initializer from attributes.
     procedure, public, pass(this) :: init_from_components => &
          volume_constraint_init_from_components
     !> Destructor.
     procedure, public, pass(this) :: free => volume_constraint_free
     !> Computes the source term and adds the result to `fields`.
     procedure, public, pass(this) :: update_value => &
          volume_constraint_update_value
     !> Computes the source term and adds the result to `fields`.
     procedure, public, pass(this) :: update_sensitivity => &
          volume_constraint_update_sensitivity
     !> Get number of log entries
     procedure, public, pass(this) :: get_log_size => &
          volume_constraint_get_log_size
     !> Get header labels for log entries
     procedure, public, pass(this) :: get_log_headers => &
          volume_constraint_get_log_headers
     !> Get values for log entries
     procedure, public, pass(this) :: get_log_values => &
          volume_constraint_get_log_values

     !> Computes the volume of the brinkman_design.
     procedure, private, pass(this) :: compute_volume

  end type volume_constraint_t

contains

  !> The common constructor using a JSON object.
  !! @param this the constraint.
  !! @param json the JSON object.
  !! @param design the design.
  !! @param simulation the simulation.
  subroutine volume_constraint_init_json_sim(this, json, design, simulation)
    class(volume_constraint_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(design_t), intent(in) :: design
    type(simulation_t), target, intent(inout) :: simulation

    character(len=:), allocatable :: mask_name
    character(len=:), allocatable :: name
    logical :: is_max
    real(kind=rp) :: limit

    call json_get_or_default(json, "name", name, "Volume constraint")
    call json_get_or_default(json, "mask_name", mask_name, "")
    call json_get_or_default(json, "is_max", is_max, .false.)
    call json_get(json, "limit", limit)

    call this%init_from_components(design, simulation, name, mask_name, &
         is_max, limit)

    ! check if we have a mapping
    if ('mapping' .in. json) then
       ! Initialize the mapper
       call this%mapping%init_base(this%c_Xh)
       call this%mapping%add(json, 'mapping')
       this%if_mapping = .true.
       ! recompute value after mapping
       call this%update_value(design)
    end if

  end subroutine volume_constraint_init_json_sim

  !> The direct initializer from attributes.
  !! @param this the constraint.
  !! @param design the design.
  !! @param simulation the simulation.
  !! @param name the name of the constraint.
  !! @param mask_name the name of the mask.
  !! @param is_max whether it is a maximum volume constraint.
  !! @param limit The volume limit value.
  subroutine volume_constraint_init_from_components(this, design, simulation, &
       name, mask_name, is_max, limit)
    class(volume_constraint_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(simulation_t), target, intent(inout) :: simulation
    character(len=*), intent(in) :: mask_name
    character(len=*), intent(in) :: name
    logical, intent(in) :: is_max
    real(kind=rp), intent(in) :: limit

    type(field_t), pointer :: work
    integer :: temp_indices(1)

    ! Initialize the base class
    call this%init_base(name, design%size(), mask_name)

    ! Store the attributes
    this%is_max = is_max
    this%limit = limit
    this%c_Xh => simulation%neko_case%fluid%c_Xh

    ! Now we can extract the mask/has_mask from the design
    if (this%has_mask) then

       ! calculate the volume of the optimization domain
       call neko_scratch_registry%request(work, temp_indices(1), .false.)
       call field_rone(work)
       call mask_exterior_const(work, this%mask, 0.0_rp)

       if (NEKO_BCKND_DEVICE .eq. 1) then
          this%volume_domain = device_glsc2(work%x_d, this%c_xh%B_d, &
               work%size())
       else
          this%volume_domain = glsc2_mask(work%x, this%c_Xh%B, &
               design%size(), this%mask%mask%get(), this%mask%size)
       end if

       call neko_scratch_registry%relinquish(temp_indices)
    else
       this%volume_domain = this%c_Xh%volume
    end if

    ! ------------------------------------------------------------------------ !
    ! Initialize the value of constraint

    call this%update_value(design)

    ! ------------------------------------------------------------------------ !
    ! Initialize the sensitivity value
    call vector_cfill(this%sensitivity, -1.0_rp / this%volume_domain)
    if (this%is_max) then
       call vector_cmult(this%sensitivity, -1.0_rp)
    end if

    if (this%has_mask) then
       call mask_exterior_const(this%sensitivity, this%mask, 0.0_rp)
    end if

  end subroutine volume_constraint_init_from_components

  !> Destructor.
  subroutine volume_constraint_free(this)
    class(volume_constraint_t), intent(inout) :: this

    call this%free_base()
  end subroutine volume_constraint_free

  !> The computation of the constraint.
  !! @param this the constraint
  !! @param design the design
  subroutine volume_constraint_update_value(this, design)
    class(volume_constraint_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    real(kind=rp) :: volume

    volume = this%compute_volume(design)

    ! Compute the distance to the target volume
    this%value = this%limit - volume / this%volume_domain

    ! Invert the sign if it is a maximum constraint
    if (this%is_max) this%value = -this%value

  end subroutine volume_constraint_update_value

  !> The computation of the sensitivity.
  !! @param this the constraint
  !! @param design the design
  subroutine volume_constraint_update_sensitivity(this, design)
    class(volume_constraint_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    !> Sensitivity field
    type(field_t), pointer :: unmapped, mapped
    integer :: idx(2)

    if (this%if_mapping) then
       ! Recompute and map backward
       call neko_scratch_registry%request(unmapped, idx(1), .false.)
       call neko_scratch_registry%request(mapped, idx(2), .false.)
       call field_cfill(unmapped, -1.0_rp / this%volume_domain)
       if (this%is_max) then
          call field_cmult(unmapped, -1.0_rp)
       end if
       ! mask
       if (this%has_mask) then
          call mask_exterior_const(unmapped, this%mask, 0.0_rp)
       end if
       ! map backwards
       call this%mapping%apply_backward(mapped, unmapped)
       if (this%has_mask) then
          call mask_exterior_const(mapped, this%mask, 0.0_rp)
       end if
       call field_to_vector(this%sensitivity, mapped)

       call neko_scratch_registry%relinquish(idx)

    else
       ! Sensitivity is just a constant so it should not be updated
    end if

  end subroutine volume_constraint_update_sensitivity

  ! ========================================================================== !
  ! The actual volume computations for different types of designs

  !> Computes the volume of the design.
  !!
  !! Automatically select which design type, or throw an error.
  !! @param this the constraint.
  !! @param design the design.
  function compute_volume(this, design) result(volume)
    class(volume_constraint_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    real(kind=rp) :: volume

    volume = 0.0_rp
    select type (design)
    type is (brinkman_design_t)
       volume = volume_brinkman_design(this, design)

    class default
       call neko_error('Volume constraint only works with brinkman_design')
    end select

  end function compute_volume

  !> Computes the volume of the brinkman_design.
  !! @param this the constraint.
  !! @param design the design.
  function volume_brinkman_design(this, design) result(volume)
    class(volume_constraint_t), intent(inout) :: this
    type(brinkman_design_t), intent(in) :: design
    real(kind=rp) :: volume
    type(field_t), pointer :: work
    type(vector_t), pointer :: values, unmapped_values
    integer :: temp_indices, ind_value, ind_um_value

    call neko_scratch_registry%request(values, ind_value, design%size(), &
         .false.)

    volume = 0.0_rp
    if (this%if_mapping) then
       call neko_scratch_registry%request(unmapped_values, ind_um_value, &
            design%size(), .false.)
       call design%get_values(unmapped_values)
       call this%mapping%apply_forward(values, unmapped_values)
       call neko_scratch_registry%relinquish(ind_um_value)
    else
       call design%get_values(values)
    end if

    if (this%has_mask) then

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call neko_scratch_registry%request(work, temp_indices, .false.)
          call device_copy(work%x_d, values%x_d, design%size())
          call mask_exterior_const(work, this%mask, 0.0_rp)

          volume = device_glsc2(work%x_d, this%c_xh%B_d, design%size())

          call neko_scratch_registry%relinquish(temp_indices)
       else
          volume = glsc2_mask(values%x, this%c_Xh%B, design%size(), &
               this%mask%mask%get(), this%mask%size)
       end if

    else

       if (NEKO_BCKND_DEVICE .eq. 1) then
          volume = device_glsc2(values%x_d, this%c_xh%B_d, design%size())
       else
          volume = glsc2(values%x, this%c_Xh%B, design%size())
       end if

    end if

    call neko_scratch_registry%relinquish(ind_value)

  end function volume_brinkman_design

  !> Return number of log entries for volume constraint.
  !! @param[in] this The constraint object.
  !! @return n Number of log entries.
  function volume_constraint_get_log_size(this) result(n)
    class(volume_constraint_t), intent(in) :: this
    integer :: n

    n = 2
  end function volume_constraint_get_log_size

  !> Populate log header labels for volume constraint.
  !! @param[in] this The constraint object.
  !! @param[out] headers Header labels for each log entry.
  subroutine volume_constraint_get_log_headers(this, headers)
    class(volume_constraint_t), intent(in) :: this
    character(len=*), intent(out) :: headers(:)
    character(len=64) :: prefix

    if (size(headers) .lt. 1) return
    prefix = trim(this%name)
    headers(1) = prefix
    if (size(headers) .lt. 2) return
    headers(2) = trim(prefix) // '.volume'

  end subroutine volume_constraint_get_log_headers

  !> Populate log values for volume constraint.
  !! @param[in] this The constraint object.
  !! @param[out] values Values corresponding to the log headers.
  subroutine volume_constraint_get_log_values(this, values)
    class(volume_constraint_t), intent(in) :: this
    real(kind=rp), intent(out) :: values(:)
    real(kind=rp) :: volume_ratio

    if (size(values) .lt. 1) return
    if (this%is_max) then
       volume_ratio = this%limit + this%value
    else
       volume_ratio = this%limit - this%value
    end if

    values(1) = this%value
    if (size(values) .lt. 2) return
    values(2) = volume_ratio

  end subroutine volume_constraint_get_log_values

end module volume_constraint
