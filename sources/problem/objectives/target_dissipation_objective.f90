! Copyright (c) 2024-25, The Neko Authors
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
!> Implements the `target_dissipation_objective_t` type.
!
module target_dissipation_objective
  use num_types, only: rp
  use field, only: field_t
  use field_math, only: field_col3, field_addcol3, field_copy
  use operators, only: grad
  use adjoint_fluid_pnpn, only: adjoint_fluid_pnpn_t
  use registry, only: neko_registry
  use adjoint_target_dissipation_source_term, only: &
       adjoint_target_dissipation_source_term_t
  use objective, only: objective_t
  use simulation_m, only: simulation_t
  use adjoint_fluid_scheme, only: adjoint_fluid_scheme_t
  use coefs, only: coef_t
  use scratch_registry, only: neko_scratch_registry
  use neko_config, only: NEKO_BCKND_DEVICE
  use math, only: glsc2, copy
  use device_math, only: device_copy, device_glsc2
  use design, only: design_t
  use point_zone, only: point_zone_t
  use mask_ops, only: mask_exterior_const, compute_masked_volume
  use math_ext, only: glsc2_mask
  use utils, only: neko_error
  use json_module, only: json_file
  use json_utils, only: json_get_or_default, json_get
  implicit none
  private

  !> An objective based on total dissipation density
  !! \f$ |\nabla \mathbf{u}|^2 + \chi |\mathbf{u}|^2 \f$.
  type, public, extends(objective_t) :: target_dissipation_objective_t
     private

     !> Pointer to the u field.
     type(field_t), pointer :: u => null()
     !> Pointer to the v field.
     type(field_t), pointer :: v => null()
     !> Pointer to the w field.
     type(field_t), pointer :: w => null()
     !> Pointer to the coefficient field.
     type(coef_t), pointer :: c_Xh => null()
     !> Pointer to Brinkman amplitude (chi).
     type(field_t), pointer :: brinkman_amplitude => null()
     !> Pointer to adjoint u field.
     type(field_t), pointer :: adjoint_u => null()
     !> Pointer to adjoint v field.
     type(field_t), pointer :: adjoint_v => null()
     !> Pointer to adjoint w field.
     type(field_t), pointer :: adjoint_w => null()
     !> Volume of the objective domain.
     real(kind=rp) :: volume
     !> Initial dissipation.
     real(kind=rp) :: initial_dissipation = 0.0_rp
     !> Current dissipation.
     real(kind=rp) :: current_dissipation = 0.0_rp
     !> Current viscous contribution.
     real(kind=rp) :: current_viscous_dissipation = 0.0_rp
     !> Current Brinkman contribution.
     real(kind=rp) :: current_brinkman_dissipation = 0.0_rp
     !> logical for first time
     logical :: is_first_time = .true.
     !> Target fraction of initial dissipation.
     real(kind=rp) :: target_fraction = 0.0_rp

   contains
     !> The common constructor using a JSON object.
     procedure, public, pass(this) :: init_json_sim => &
          target_dissipation_init_json_sim
     !> The direct initializer from attributes.
     procedure, public, pass(this) :: init_from_attributes => &
          target_dissipation_init_attributes
     !> Destructor.
     procedure, public, pass(this) :: free => target_dissipation_free
     !> Computes the value of the objective function.
     procedure, public, pass(this) :: update_value => &
          target_dissipation_update_value
     !> Computes the sensitivity with respect to the coefficient \f$\chi\f$.
     procedure, public, pass(this) :: update_sensitivity => &
          target_dissipation_update_sensitivity
     !> Log sizing and values
     procedure, public, pass(this) :: get_log_size => &
          target_dissipation_get_log_size
     procedure, public, pass(this) :: get_log_headers => &
          target_dissipation_get_log_headers
     procedure, public, pass(this) :: get_log_values => &
          target_dissipation_get_log_values
     !> This is particularly tricky... We're going to overwrite the get value
     !! to do the scaling. So the value we store in the object is the dissipation
     !! to allow it to accumulate, but then the normalization with respect
     !! to the initial only happens on the get value
     procedure, pass(this) :: finalize_value => &
          target_dissipation_finalize_value

  end type target_dissipation_objective_t

contains

  !> The common constructor using a JSON object.
  !! @param this the objective.
  !! @param json the JSON object.
  !! @param design the design.
  !! @param simulation the simulation.
  subroutine target_dissipation_init_json_sim(this, json, design, simulation)
    class(target_dissipation_objective_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(design_t), intent(in) :: design
    type(simulation_t), target, intent(inout) :: simulation

    character(len=:), allocatable :: name
    character(len=:), allocatable :: mask_name
    real(kind=rp) :: weight, target_fraction

    call json_get_or_default(json, "weight", weight, 1.0_rp)
    call json_get_or_default(json, "mask_name", mask_name, "")
    call json_get_or_default(json, "name", name, "Target Dissipation")
    call json_get(json, "target", target_fraction)

    call this%init_from_attributes(design, simulation, weight, name, &
         mask_name, target_fraction)
  end subroutine target_dissipation_init_json_sim

  !> The actual constructor.
  !! @param this the objective.
  !! @param design the design.
  !! @param simulation the simulation.
  !! @param weight the weight of the objective function.
  !! @param name the name of the objective.
  !! @param mask_name the name of the mask.
  !! @param target_fraction target fraction of the initial dissipation.
  subroutine target_dissipation_init_attributes(this, design, simulation, &
       weight, name, mask_name, target_fraction)
    class(target_dissipation_objective_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(simulation_t), target, intent(inout) :: simulation
    real(kind=rp), intent(in) :: weight
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: mask_name
    real(kind=rp), intent(in) :: target_fraction
    type(adjoint_target_dissipation_source_term_t) :: adjoint_forcing

    call this%init_base(name, design%size(), weight, mask_name)

    ! Save the simulation and design
    this%u => neko_registry%get_field('u')
    this%v => neko_registry%get_field('v')
    this%w => neko_registry%get_field('w')
    this%brinkman_amplitude => neko_registry%get_field('brinkman_amplitude')
    this%c_Xh => simulation%fluid%c_Xh
    this%adjoint_u => neko_registry%get_field('u_adj')
    this%adjoint_v => neko_registry%get_field('v_adj')
    this%adjoint_w => neko_registry%get_field('w_adj')
    this%target_fraction = target_fraction

    ! compute the volume of the objective domain
    if (this%has_mask) then
       this%volume = compute_masked_volume(this%mask, this%c_Xh)
    else
       this%volume = this%c_Xh%volume
    end if
    ! append a source term based on the target dissipation
    ! init the adjoint forcing term for the adjoint
    call adjoint_forcing%init_from_components( &
         simulation%adjoint_fluid%f_adj_x, &
         simulation%adjoint_fluid%f_adj_y, &
         simulation%adjoint_fluid%f_adj_z, &
         this%u, this%v, this%w, this%brinkman_amplitude, this%weight, &
         this%mask, this%has_mask, &
         this%c_Xh, this%volume, this%target_fraction, &
         this%current_dissipation, this%initial_dissipation)

    ! append adjoint forcing term based on objective function
    select type (f => simulation%adjoint_fluid)
    type is (adjoint_fluid_pnpn_t)
       call f%source_term%add_source_term(adjoint_forcing)
    end select

  end subroutine target_dissipation_init_attributes

  !> Destructor.
  subroutine target_dissipation_free(this)
    class(target_dissipation_objective_t), intent(inout) :: this
    call this%free_base()

    if (associated(this%u)) nullify(this%u)
    if (associated(this%v)) nullify(this%v)
    if (associated(this%w)) nullify(this%w)
    if (associated(this%brinkman_amplitude)) nullify(this%brinkman_amplitude)
    if (associated(this%c_Xh)) nullify(this%c_Xh)

    if (associated(this%adjoint_u)) nullify(this%adjoint_u)
    if (associated(this%adjoint_v)) nullify(this%adjoint_v)
    if (associated(this%adjoint_w)) nullify(this%adjoint_w)

  end subroutine target_dissipation_free

  !> Compute the objective function.
  !! @param this the objective.
  !! @param design the design.
  subroutine target_dissipation_update_value(this, design)
    class(target_dissipation_objective_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(field_t), pointer :: wo1, wo2, wo3, work
    type(field_t), pointer :: objective_field
    integer :: temp_indices(5)
    integer n
    real(kind=rp) :: viscous_integral, brinkman_integral

    call neko_scratch_registry%request_field(wo1, temp_indices(1), .false.)
    call neko_scratch_registry%request_field(wo2, temp_indices(2), .false.)
    call neko_scratch_registry%request_field(wo3, temp_indices(3), .false.)
    call neko_scratch_registry%request_field(objective_field, temp_indices(4), &
         .false.)
    call neko_scratch_registry%request_field(work, temp_indices(5), .false.)

    ! Compute viscous dissipation contribution |\nabla u|^2.
    call grad(wo1%x, wo2%x, wo3%x, this%u%x, this%c_Xh)
    call field_col3(objective_field, wo1, wo1)
    call field_addcol3(objective_field, wo2, wo2)
    call field_addcol3(objective_field, wo3, wo3)

    call grad(wo1%x, wo2%x, wo3%x, this%v%x, this%c_Xh)
    call field_addcol3(objective_field, wo1, wo1)
    call field_addcol3(objective_field, wo2, wo2)
    call field_addcol3(objective_field, wo3, wo3)

    call grad(wo1%x, wo2%x, wo3%x, this%w%x, this%c_Xh)
    call field_addcol3(objective_field, wo1, wo1)
    call field_addcol3(objective_field, wo2, wo2)
    call field_addcol3(objective_field, wo3, wo3)

    ! Integrate viscous contribution.
    n = wo1%size()
    if (this%has_mask) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          ! note, this could be done more elagantly by writing
          ! device_glsc2_mask
          call field_copy(work, objective_field)
          call mask_exterior_const(work, this%mask, 0.0_rp)
          viscous_integral = device_glsc2(work%x_d, this%c_xh%B_d, n)
       else
          viscous_integral = glsc2_mask(objective_field%x, this%c_Xh%b, &
               n, this%mask%mask%get(), this%mask%size)
       end if
    else
       if (neko_bcknd_device .eq. 1) then
          viscous_integral = device_glsc2(objective_field%x_d, &
               this%c_Xh%b_d, n)
       else
          viscous_integral = glsc2(objective_field%x, this%c_Xh%b, n)
       end if
    end if

    ! Integrate Brinkman contribution \chi * |u|^2.
    call field_col3(work, this%u, this%u)
    call field_addcol3(work, this%v, this%v)
    call field_addcol3(work, this%w, this%w)
    call field_col3(work, work, this%brinkman_amplitude)

    if (this%has_mask) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call mask_exterior_const(work, this%mask, 0.0_rp)
          brinkman_integral = device_glsc2(work%x_d, this%c_xh%B_d, n)
       else
          brinkman_integral = glsc2_mask(work%x, this%c_Xh%b, &
               n, this%mask%mask%get(), this%mask%size)
       end if
    else
       if (neko_bcknd_device .eq. 1) then
          brinkman_integral = device_glsc2(work%x_d, this%c_Xh%b_d, n)
       else
          brinkman_integral = glsc2(work%x, this%c_Xh%b, n)
       end if
    end if

    this%current_viscous_dissipation = 0.5_rp * viscous_integral / this%volume
    this%current_brinkman_dissipation = 0.5_rp * brinkman_integral / this%volume
    this%value = this%current_viscous_dissipation + &
         this%current_brinkman_dissipation

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine target_dissipation_update_value

  subroutine target_dissipation_finalize_value(this)
    class(target_dissipation_objective_t), intent(inout) :: this

    ! Check if it's the first time, and if so, set the initial dissipation
    if (this%is_first_time) then
       this%initial_dissipation = this%value
       this%is_first_time = .false.
    end if
    ! For clarity write the current dissipation clearly
    this%current_dissipation = this%value

    ! Now compute the objective
    this%value = 0.5_rp * (this%current_dissipation / &
         (this%initial_dissipation * this%target_fraction) - 1.0_rp) ** 2

  end subroutine target_dissipation_finalize_value

  !> update_value the sensitivity of the objective function with respect to
  !! \f$\chi\f$
  !! @param this the objective.
  !! @param design the design.
  subroutine target_dissipation_update_sensitivity(this, design)
    class(target_dissipation_objective_t), intent(inout) :: this
    class(design_t), intent(in) :: design

  end subroutine target_dissipation_update_sensitivity

  !> Number of log entries
  function target_dissipation_get_log_size(this) result(n)
    class(target_dissipation_objective_t), intent(in) :: this
    integer :: n

    n = 7
  end function target_dissipation_get_log_size

  !> Header labels for log entries
  subroutine target_dissipation_get_log_headers(this, headers)
    class(target_dissipation_objective_t), intent(in) :: this
    character(len=*), intent(inout) :: headers(:)
    character(len=64) :: prefix

    if (size(headers) .lt. 1) return
    prefix = trim(this%name)
    headers(1) = prefix
    if (size(headers) .lt. 2) return
    headers(2) = trim(prefix) // '.weight'
    if (size(headers) .lt. 3) return
    headers(3) = trim(prefix) // '.current'
    if (size(headers) .lt. 4) return
    headers(4) = trim(prefix) // '.initial'
    if (size(headers) .lt. 5) return
    headers(5) = trim(prefix) // '.ratio'
    if (size(headers) .lt. 6) return
    headers(6) = trim(prefix) // '.viscous'
    if (size(headers) .lt. 7) return
    headers(7) = trim(prefix) // '.brinkman'

  end subroutine target_dissipation_get_log_headers

  !> Values for log entries
  subroutine target_dissipation_get_log_values(this, values)
    class(target_dissipation_objective_t), intent(in) :: this
    real(kind=rp), intent(out) :: values(:)

    if (size(values) .lt. 1) return
    values(1) = this%value
    if (size(values) .lt. 2) return
    values(2) = this%weight
    if (size(values) .lt. 3) return
    values(3) = this%current_dissipation
    if (size(values) .lt. 4) return
    values(4) = this%initial_dissipation
    if (size(values) .lt. 5) return
    values(5) = this%current_dissipation / this%initial_dissipation
    if (size(values) .lt. 6) return
    values(6) = this%current_viscous_dissipation
    if (size(values) .lt. 7) return
    values(7) = this%current_brinkman_dissipation

  end subroutine target_dissipation_get_log_values

end module target_dissipation_objective
