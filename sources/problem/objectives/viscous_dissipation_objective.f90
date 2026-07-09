!> @file viscous_dissipation_objective.f90
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
!> Implements the `viscous_dissipation_objective_t` type.
!
! I promise I'll write this document properly in the future...
!
! But the Borrvall & Petersson (I think) paper had an objective function
! that had 2 terms, dissipation and this term they claimed represented
! out of plane stresses.
! I never really understood that extra term, I also don't think it
! applies to 3D cases, but everyone includes it anyway.
!
! It appears to me to be basically a heuristic penality that targets
! non-binary designs
!
! so let's call
!
! F = \int |\nabla u|^2  + K \int \chi \u^2
!
!      | dissipation |     |"Brinkman dissipation"|
!
! I say "Brinkman dissipation" because they said it came from lubrication
! theory...
! Anyway, we can change all this later (especially the names!)

! If the objective function \frac{\mu}{2} \int |\nabla u|^2,
! the corresponding adjoint forcing is \mu \int \nabla v \cdot \nabla u
!
! for the Brinkman dissipation, the adjoint forcing is \chi u
!
! This has always annoyed me...
! because now I see one objective and one constraint
!
module viscous_dissipation_objective
  use num_types, only: rp
  use field, only: field_t
  use field_math, only: field_col3, field_addcol3, field_cmult, field_add2s2, &
       field_copy
  use operators, only: grad
  use adjoint_fluid_pnpn, only: adjoint_fluid_pnpn_t
  use scratch_registry, only: neko_scratch_registry
  use adjoint_viscous_dissipation_source_term, only: &
       adjoint_viscous_dissipation_source_term_t
  use objective, only: objective_t
  use simulation_m, only: simulation_t
  use adjoint_fluid_scheme, only: adjoint_fluid_scheme_t
  use coefs, only: coef_t
  use registry, only: neko_registry
  use neko_config, only: NEKO_BCKND_DEVICE
  use math, only: glsc2, copy
  use device_math, only: device_copy, device_glsc2
  use design, only: design_t
  use brinkman_design, only: brinkman_design_t
  use point_zone, only: point_zone_t
  use mask_ops, only: mask_exterior_const, compute_masked_volume
  use math_ext, only: glsc2_mask
  use utils, only: neko_error
  use json_module, only: json_file
  use json_utils, only: json_get_or_default
  use continuation_scheduler, only: nekotop_continuation
  implicit none
  private

  !> An objective function corresponding to viscous dissipation
  !! \f$ F =  \int_\Omega \frac{\mu}{2} |\nabla u|^2 d \Omega \f$
  type, public, extends(objective_t) :: viscous_dissipation_objective_t
     private

     !> Pointer to the u field.
     type(field_t), pointer :: u => null()
     !> Pointer to the v field.
     type(field_t), pointer :: v => null()
     !> Pointer to the w field.
     type(field_t), pointer :: w => null()
     !> Pointer to the coefficient field.
     type(coef_t), pointer :: c_Xh => null()
     !> Pointer to adjoint u field.
     type(field_t), pointer :: adjoint_u => null()
     !> Pointer to adjoint v field.
     type(field_t), pointer :: adjoint_v => null()
     !> Pointer to adjoint w field.
     type(field_t), pointer :: adjoint_w => null()
     !> Volume of the objective domain.
     real(kind=rp) :: volume
     !> Dynamic viscosity.
     real(kind=rp) :: viscosity

   contains
     !> The common constructor using a JSON object.
     procedure, public, pass(this) :: init_json_sim => &
          viscous_dissipation_init_json_sim
     !> The direct initializer from attributes.
     procedure, public, pass(this) :: init_from_attributes => &
          viscous_dissipation_init_attributes
     !> Destructor.
     procedure, public, pass(this) :: free => viscous_dissipation_free
     !> Computes the value of the objective function.
     procedure, public, pass(this) :: update_value => &
          viscous_dissipation_update_value
     !> Computes the sensitivity with respect to the coefficient \f$\chi\f$.
     procedure, public, pass(this) :: update_sensitivity => &
          viscous_dissipation_update_sensitivity

  end type viscous_dissipation_objective_t

contains

  !> The common constructor using a JSON object.
  !! @param this the objective.
  !! @param json the JSON object.
  !! @param design the design.
  !! @param simulation the simulation.
  subroutine viscous_dissipation_init_json_sim(this, json, design, simulation)
    class(viscous_dissipation_objective_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(design_t), intent(in) :: design
    type(simulation_t), target, intent(inout) :: simulation

    character(len=:), allocatable :: name
    character(len=:), allocatable :: mask_name
    real(kind=rp) :: weight
    real(kind=rp) :: start_time, end_time

    call nekotop_continuation%json_get_or_register(json, 'weight', &
         this%weight, weight, 1.0_rp)
    call json_get_or_default(json, "mask_name", mask_name, "")
    call json_get_or_default(json, "name", name, "Viscous dissipation")
    call json_get_or_default(json, "start_time", start_time, 0.0_rp)
    call json_get_or_default(json, "end_time", end_time, huge(0.0_rp))

    call this%init_from_attributes(design, simulation, weight, name, &
         mask_name, start_time, end_time)
  end subroutine viscous_dissipation_init_json_sim

  !> The actual constructor.
  !! @param this the objective.
  !! @param design the design.
  !! @param simulation the simulation.
  !! @param weight the weight of the objective function.
  !! @param name the name of the objective.
  !! @param mask_name the name of the mask.
  !! @param start_time start of the integration window.
  !! @param end_time end of the integration window.
  subroutine viscous_dissipation_init_attributes(this, design, simulation, &
       weight, name, mask_name, start_time, end_time)
    class(viscous_dissipation_objective_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(simulation_t), target, intent(inout) :: simulation
    real(kind=rp), intent(in) :: weight
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: mask_name
    real(kind=rp), intent(in) :: start_time
    real(kind=rp), intent(in) :: end_time

    type(adjoint_viscous_dissipation_source_term_t) :: adjoint_forcing

    call this%init_base(name, design%size(), weight, mask_name, &
         start_time, end_time)

    ! Save the simulation and design
    this%u => neko_registry%get_field('u')
    this%v => neko_registry%get_field('v')
    this%w => neko_registry%get_field('w')
    this%c_Xh => simulation%fluid%c_Xh
    this%adjoint_u => neko_registry%get_field('u_adj')
    this%adjoint_v => neko_registry%get_field('v_adj')
    this%adjoint_w => neko_registry%get_field('w_adj')
    this%viscosity = simulation%fluid%mu%x(1,1,1,1)

    ! compute the volume of the objective domain
    if (this%has_mask) then
       this%volume = compute_masked_volume(this%mask, this%c_Xh)
    else
       this%volume = this%c_Xh%volume
    end if

    ! you will need to init this!
    ! append a source term based on the viscous dissipation
    ! init the adjoint forcing term for the adjoint
    call adjoint_forcing%init_from_components( &
         simulation%adjoint_fluid%f_adj_x, &
         simulation%adjoint_fluid%f_adj_y, &
         simulation%adjoint_fluid%f_adj_z, &
         this%u, this%v, this%w, this%weight * this%viscosity, &
         this%mask, this%has_mask, &
         this%c_Xh, this%volume, this%start_time, this%end_time)

    ! append adjoint forcing term based on objective function
    select type (f => simulation%adjoint_fluid)
    type is (adjoint_fluid_pnpn_t)
       call f%source_term%add_source_term(adjoint_forcing)
    end select

  end subroutine viscous_dissipation_init_attributes

  !> Destructor.
  subroutine viscous_dissipation_free(this)
    class(viscous_dissipation_objective_t), intent(inout) :: this
    call this%free_base()

    if (associated(this%u)) nullify(this%u)
    if (associated(this%v)) nullify(this%v)
    if (associated(this%w)) nullify(this%w)
    if (associated(this%c_Xh)) nullify(this%c_Xh)

    if (associated(this%adjoint_u)) nullify(this%adjoint_u)
    if (associated(this%adjoint_v)) nullify(this%adjoint_v)
    if (associated(this%adjoint_w)) nullify(this%adjoint_w)

  end subroutine viscous_dissipation_free

  !> Compute the objective function.
  !! @param this the objective.
  !! @param design the design.
  subroutine viscous_dissipation_update_value(this, design)
    class(viscous_dissipation_objective_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(field_t), pointer :: wo1, wo2, wo3, work
    type(field_t), pointer :: objective_field
    integer :: temp_indices(5)
    integer n

    call neko_scratch_registry%request_field(wo1, temp_indices(1), .false.)
    call neko_scratch_registry%request_field(wo2, temp_indices(2), .false.)
    call neko_scratch_registry%request_field(wo3, temp_indices(3), .false.)
    call neko_scratch_registry%request_field(objective_field, temp_indices(4), &
         .false.)
    call neko_scratch_registry%request_field(work, temp_indices(5), .false.)

    ! update_value the objective function.
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

    ! integrate the field
    n = wo1%size()
    if (this%has_mask) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          ! note, this could be done more elagantly by writing
          ! device_glsc2_mask
          call field_copy(work, objective_field)
          call mask_exterior_const(work, this%mask, 0.0_rp)
          this%value = device_glsc2(work%x_d, this%c_xh%B_d, n)
       else
          this%value = glsc2_mask(objective_field%x, this%c_Xh%b, &
               n, this%mask%mask%get(), this%mask%size)
       end if
    else
       if (neko_bcknd_device .eq. 1) then
          this%value = device_glsc2(objective_field%x_d, &
               this%c_Xh%b_d, n)
       else
          this%value = glsc2(objective_field%x, this%c_Xh%b, n)
       end if
    end if

    this%value = this%value * 0.5_rp * this%viscosity / this%volume

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine viscous_dissipation_update_value

  !> Compute the objective sensitivity with respect to \f$\chi\f$.
  !! @param this the objective.
  !! @param design the design.
  subroutine viscous_dissipation_update_sensitivity(this, design)
    class(viscous_dissipation_objective_t), intent(inout) :: this
    class(design_t), intent(in) :: design

  end subroutine viscous_dissipation_update_sensitivity

end module viscous_dissipation_objective
