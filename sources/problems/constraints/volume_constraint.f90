! Copyright (c) 2023, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> Implements the `volume_constraint_t` type.
! $V = \frac{1}{|\Omega_O|}\int_{\Omega_O} \tilde{\rho} d\Omega$
! Either
! $V < V_\text{max}$
! $V > V_\text{min}$
module volume_constraint
  use num_types, only : rp
  use field_list, only : field_list_t
  use json_utils, only: json_get, json_get_or_default
  use source_term, only : source_term_t
  use coefs, only : coef_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use utils, only : neko_error
  use field, only: field_t
  use field_math, only: field_col3, field_addcol3
  use user_intf, only: user_t, simulation_component_user_settings
  use json_module, only: json_file
  use steady_simcomp, only: steady_simcomp_t
  use simcomp_executor, only: neko_simcomps
  use fluid_user_source_term, only: fluid_user_source_term_t
  use num_types, only : rp
  use field, only : field_t
  use field_registry, only : neko_field_registry
  use math, only : rzero, copy, chsign
  use device_math, only: device_copy, device_cmult
  use neko_config, only: NEKO_BCKND_DEVICE
  use operators, only: curl, grad
  use scratch_registry, only : neko_scratch_registry
  use constraint, only : constraint_t
  use simulation, only : simulation_t
  use adjoint_scheme, only : adjoint_scheme_t
  use fluid_source_term, only: fluid_source_term_t
  use math, only : glsc2
  use field_math, only: field_rone, field_cmult
  use design, only: design_t
  use topopt_design, only: topopt_design_t
  use mask_ops, only: mask_exterior_const
  use math_ext, only: glsc2_mask
  use json_module, only: json_file
  use json_utils, only: json_get_or_default
  implicit none
  private

  !> A constraint on the volume of the design.
  type, public, extends(constraint_t) :: volume_constraint_t
     private

     !> whether it is minimum or maximum volume
     !! is_max = .false., 	ie V > V_min  		=>		 -V + V_max < 0
     !! is_max = .true. , 	ie V < V_max  		=>		  V - V_max < 0
     logical :: is_max
     !> Maximum (or minimum) volume
     real(kind=rp) :: limit
     !> Volume of the optimization domain
     real(kind=rp) :: volume_domain

     !> Pointer to the SEM field.
     class(coef_t), pointer :: c_Xh => null()

   contains
     !> The common constructor using a JSON object.
     procedure, public, pass(this) :: init_json => volume_constraint_init_json
     !> The direct initializer from attributes.
     procedure, public, pass(this) :: init_from_attributes => &
          volume_constraint_init_attributes
     !> Destructor.
     procedure, public, pass(this) :: free => volume_constraint_free
     !> Computes the source term and adds the result to `fields`.
     procedure, public, pass(this) :: update_value => &
          volume_constraint_update_value
     !> Computes the source term and adds the result to `fields`.
     procedure, public, pass(this) :: update_sensitivity => &
          volume_constraint_update_sensitivity

     !> Computes the volume of the topopt_design.
     procedure, private, pass(this) :: compute_volume

  end type volume_constraint_t

contains

  !> The common constructor using a JSON object.
  !! @param json the JSON object.
  !! @param design the design.
  !! @param simulation the simulation.
  subroutine volume_constraint_init_json(this, json, design, simulation)
    class(volume_constraint_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(design_t), intent(in) :: design
    type(simulation_t), target, intent(inout) :: simulation

    character(len=:), allocatable :: mask_name
    character(len=:), allocatable :: log_name
    logical :: is_max
    real(kind=rp) :: limit

    call json_get_or_default(json, "mask_name", mask_name, "")
    call json_get_or_default(json, "log_name", log_name, "Volume constraint")
    call json_get_or_default(json, "is_max", is_max, .false.)
    call json_get(json, "limit", limit)

    call this%init_from_attributes(design, simulation, log_name, mask_name, &
       is_max, limit)
  end subroutine volume_constraint_init_json

  !> The direct initializer from attributes.
  !! @param design the design.
  !! @param simulation the simulation.
  !! @param mask_name the name of the mask.
  !! @param is_max whether it is a maximum volume constraint.
  !! @param limit The volume limit value.
  subroutine volume_constraint_init_attributes(this, design, simulation, &
       log_name, mask_name, is_max, limit)
    class(volume_constraint_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(simulation_t), target, intent(inout) :: simulation
    character(len=*), intent(in) :: mask_name
    character(len=*), intent(in) :: log_name
    logical, intent(in) :: is_max
    real(kind=rp), intent(in) :: limit

    real(kind=rp) :: volume
    type(field_t), pointer :: work
    integer :: temp_indices(1)

    ! Initialize the base class
    call this%init_base(design%size(), mask_name)

    ! Assign the name to constraint
    this%log_name = log_name 

    ! Store the attributes
    this%is_max = is_max
    this%limit = limit
    this%c_Xh => simulation%neko_case%fluid%c_Xh

    ! Now we can extract the mask/has_mask from the design
    if (this%has_mask) then

       ! calculate the volume of the optimization domain
       call neko_scratch_registry%request_field(work, temp_indices(1))
       call field_rone(work)

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call neko_error('GPU not supported volume constraint')
       else
          this%volume_domain = glsc2_mask(work%x, this%c_Xh%B, &
               design%size(), this%mask%mask, this%mask%size)
       end if

       call neko_scratch_registry%relinquish_field(temp_indices)
    else
       this%volume_domain = this%c_Xh%volume
    end if

    ! ------------------------------------------------------------------------ !
    ! Initialize the value of constraint

    ! Compute the volume of the design
    volume = this%compute_volume(design)

    ! Compute the distance to the target volume
    this%value = this%limit - volume / this%volume_domain

    ! Invert the sign if it is a maximum constraint
    if (this%is_max) this%value = - ( this%value )

    ! ------------------------------------------------------------------------ !
    ! Initialize the sensitivity value

    this%sensitivity = 1.0_rp / this%volume_domain

    ! Invert the sign if it is a maximum constraint
    if (.not. this%is_max) this%sensitivity = (-1.0_rp) * this%sensitivity

    if (this%has_mask) then
       call mask_exterior_const(this%sensitivity, this%mask, 0.0_rp)
    end if

  end subroutine volume_constraint_init_attributes

  !> Destructor.
  subroutine volume_constraint_free(this)
    class(volume_constraint_t), intent(inout) :: this

    call this%free_base()
  end subroutine volume_constraint_free

  !> The computation of the constraint.
  !! @param design the design
  !! @param fluid the fluid scheme
  subroutine volume_constraint_update_value(this, design)
    class(volume_constraint_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    real(kind=rp) :: volume

    volume = this%compute_volume(design)

    ! Compute the distance to the target volume
    this%value = this%limit - volume / this%volume_domain

    ! Invert the sign if it is a maximum constraint
    if (this%is_max) this%value = - ( this%value )

  end subroutine volume_constraint_update_value

  !> The computation of the sensitivity.
  !! @param design the design
  !! @param fluid the fluid scheme
  !! @param adjoint the adjoint scheme
  subroutine volume_constraint_update_sensitivity(this, design)
    class(volume_constraint_t), intent(inout) :: this
    class(design_t), intent(in) :: design

    ! Sensitivity is just a constant so it should not be updated

  end subroutine volume_constraint_update_sensitivity


  ! ========================================================================== !
  ! The actual volume computations for different types of designs


  !> Computes the volume of the design.
  !!
  !! Automatically select which design type, or throw an error.
  !! @param design the design.
  function compute_volume(this, design) result(volume)
    class(volume_constraint_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    real(kind=rp) :: volume

    volume = 0.0_rp
    select type (design)
      type is (topopt_design_t)
       volume = volume_topopt_design(this, design)

      class default
       call neko_error('Volume constraint only works with topopt_design')
    end select

  end function compute_volume

  !> Computes the volume of the topopt_design.
  !! @param design the design.
  function volume_topopt_design(this, design) result(volume)
    class(volume_constraint_t), intent(inout) :: this
    type(topopt_design_t), intent(in) :: design
    real(kind=rp) :: volume

    ! in the future we should be using the mapped design variable
    ! corresponding to this constraint!!!
    if (this%has_mask) then

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call neko_error('GPU not supported volume constraint')
       else
          volume = glsc2_mask(design%design_indicator%x, &
               this%c_Xh%B, design%size(), this%mask%mask, this%mask%size)
       end if

    else

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call neko_error('GPU not supported volume constraint')
       else
          volume = glsc2(design%design_indicator%x, &
               this%c_Xh%B, design%size())
       end if

    end if

  end function volume_topopt_design

end module volume_constraint
