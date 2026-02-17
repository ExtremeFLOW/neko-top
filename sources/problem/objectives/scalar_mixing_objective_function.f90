!> @file scalar_mixing_objective_function.f90
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
!> An objective function corresponding to the mixing of a passive scalar
!! \f$ F = \frac{1}{|\Omega_{obj}|}\int_{\Omega_{obj}}
!! \frac{1}{2}\left(\phi - \phi_{ref}\right)^2 d\Omega, \f$
module scalar_mixing_objective
  use num_types, only: rp
  use objective, only: objective_t
  use simulation_m, only: simulation_t
  use design, only: design_t
  use json_module, only: json_file
  use json_utils, only: json_get_or_default
  use field, only: field_t
  use field_math, only: field_copy, field_cadd, field_col2, field_cmult
  use neko_config, only: NEKO_BCKND_DEVICE
  use math, only: glsc2, copy
  use device_math, only: device_glsc2, device_copy
  use mask_ops, only: mask_exterior_const, compute_masked_volume
  use coefs, only: coef_t
  use scratch_registry, only: neko_scratch_registry
  use utils, only: neko_error
  use adjoint_mixing_scalar_source_term, only: &
       adjoint_mixing_scalar_source_term_t
  use neko_ext, only: get_scalar_indicies
  ! delete after the simulation computes u u_adj
  use field_math, only: field_addcol3, field_col3
  implicit none
  private

  !> An objective function corresponding to the mixing of a passive scalar
  !! \f$ F = \frac{1}{|\Omega_{obj}|}\int_{\Omega_{obj}}
  !! \frac{1}{2}\left(\phi - \phi_{ref}\right)^2 d\Omega, \f$
  type, public, extends(objective_t) :: scalar_mixing_objective_t
     private

     !> pointer to the primal passive scalar fields \f$\phi\f$
     type(field_t), pointer :: phi
     !> Target concentration in the optimized region \f$\phi_{ref}\f$
     real(kind=rp) :: phi_ref
     !> Coefficients defined on a given mesh
     class(coef_t), pointer :: coef
     !> Volume of the domain \f$|\Omega_{obj}|\f$
     real(kind=rp) :: domain_volume
     !> name of the scalar field being acted on
     character(len=:), allocatable :: scalar_name

     ! -----------------------------------------------------------------
     ! THESE SHOULD BE DELETED WHEN THE DESIGN UPDATE COMES IN.
     ! These sensitivity should return zero, and the simulation
     ! should compute the u u_adj component
     type(field_t), pointer :: u, v, w, u_adj, v_adj, w_adj

   contains

     !> The common constructor using a JSON object.
     procedure, public, pass(this) :: init_json_sim => &
          scalar_mixing_init_json_sim
     !> The actual constructor.
     procedure, public, pass(this) :: init_from_attributes => &
          scalar_mixing_init_attributes
     !> Destructor.
     procedure, public, pass(this) :: free => scalar_mixing_free
     !> Computes the value of the objective function.
     procedure, public, pass(this) :: update_value => &
          scalar_mixing_update_value
     !> Computes the sensitivity with respect to the coefficient \f$\chi\f$.
     procedure, public, pass(this) :: update_sensitivity => &
          scalar_mixing_update_sensitivity

  end type scalar_mixing_objective_t

contains

  !> The common constructor using a JSON object.
  !! @param this The object
  !! @param json The JSON subdictionary corresponding to your objective
  !! @param design The design
  !! @param simulation The simulation
  subroutine scalar_mixing_init_json_sim(this, json, design, simulation)
    class(scalar_mixing_objective_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(design_t), intent(in) :: design
    type(simulation_t), target, intent(inout) :: simulation
    real(kind=rp) :: weight
    real(kind=rp) :: phi_ref
    character(len=:), allocatable :: name
    character(len=:), allocatable :: mask_name
    character(len=:), allocatable :: scalar_name

    call json_get_or_default(json, "weight", weight, 1.0_rp)
    call json_get_or_default(json, "mask_name", mask_name, "")
    call json_get_or_default(json, "target_concentration", phi_ref, 0.5_rp)
    call json_get_or_default(json, "name", name, "Scalar Mixing")
    call json_get_or_default(json, "scalar_name", scalar_name, "s")

    ! initialize
    call this%init_from_attributes(design, simulation, weight, name, &
         mask_name, phi_ref, scalar_name)
  end subroutine scalar_mixing_init_json_sim

  !> The actual constructor.
  !! @param this The object
  !! @param design the design.
  !! @param simulation the simulation.
  !! @param weight the weight of the objective function.
  !! @param name the name of the objective function.
  !! @param mask_name the name of the mask.
  !! @param phi_ref target concentration used in the objective function.
  !! @param scalar_name name of the scalar field.
  subroutine scalar_mixing_init_attributes(this, design, simulation, weight, &
       name, mask_name, phi_ref, scalar_name)
    class(scalar_mixing_objective_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(simulation_t), target, intent(inout) :: simulation
    real(kind=rp), intent(in) :: weight
    real(kind=rp), intent(in) :: phi_ref
    character(len=*), intent(in) :: mask_name
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: scalar_name

    type(adjoint_mixing_scalar_source_term_t) :: adjoint_forcing
    integer :: i_scalar, i_adjoint_scalar

    ! Start by checking if the adjoint scalar has been initialized
    if (.not.allocated(simulation%adjoint_case%adjoint_scalars)) then
       call neko_error("adjoint passive scalar not initialized")
    end if

    ! Call the base initializer
    call this%init_base(name, design%size(), weight, mask_name)

    ! Associate the integration weights
    this%coef => simulation%fluid%c_Xh

    ! Compute the masked volume
    if (this%has_mask) then
       this%domain_volume = compute_masked_volume(this%mask, this%coef)
    else
       this%domain_volume = this%coef%volume
    end if

    this%scalar_name = trim(scalar_name)

    ! figure out the index associated with the scalar and adjoint scalar.
    call get_scalar_indicies(i_scalar, i_adjoint_scalar, simulation%scalars, &
         simulation%adjoint_scalars, this%scalar_name)

    !> Associate the RHS of the passive scalar equation
    !! \f$ f_{\phi^\dagger} \f$
    associate(f_phi_adj => &
         simulation%adjoint_scalars%adjoint_scalar_fields(i_adjoint_scalar)%f_Xh)

      ! Associate json parameters
      this%phi_ref = phi_ref

      ! Associate forward passive scalar
      this%phi => simulation%scalars%scalar_fields(i_scalar)%scalar%s

      ! Initialize the scalar mixing adjoint source term
      call adjoint_forcing%init_from_components(f_phi_adj, this%phi, &
           this%get_weight(), this%phi_ref, this%mask, this%has_mask, this%coef)

    end associate

    ! append adjoint source term to the adjoint passive scalar equation
    call simulation%adjoint_scalars%adjoint_scalar_fields(i_adjoint_scalar) &
         %source_term%add_source_term(adjoint_forcing)

    !--------------------------------------------------------------------------
    ! THIS SHOULD BE REPLACED WHEN THE DESIGN UPDATE OCCURS
    this%u => simulation%fluid%u
    this%v => simulation%fluid%v
    this%w => simulation%fluid%w
    this%u_adj => simulation%adjoint_fluid%u_adj
    this%v_adj => simulation%adjoint_fluid%v_adj
    this%w_adj => simulation%adjoint_fluid%w_adj
  end subroutine scalar_mixing_init_attributes

  !> Destructor.
  subroutine scalar_mixing_free(this)
    class(scalar_mixing_objective_t), intent(inout) :: this
    call this%free_base()

    this%u => null()
    this%v => null()
    this%w => null()
    this%phi => null()
    this%u_adj => null()
    this%v_adj => null()
    this%w_adj => null()

  end subroutine scalar_mixing_free

  !> Compute the objective function.
  !! @param this the object.
  !! @param design the design.
  subroutine scalar_mixing_update_value(this, design)
    class(scalar_mixing_objective_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(field_t), pointer :: work
    integer :: temp_indices(1), n


    ! The objective being computed is
    !
    ! \f$1/2 \int_\Omega_{obj} (\phi - \phi_ref)^2 d\Omega\f$

    ! get a working array
    call neko_scratch_registry%request_field(work, temp_indices(1), .false.)
    n = work%size()

    ! \f$ \phi \f$
    call field_copy(work, this%phi)
    ! \f$ \phi - \phi_ref \f$
    call field_cadd(work, -this%phi_ref)
    ! \f$ (\phi - \phi_ref)^2 \f$
    call field_col2(work, work)
    ! \f$ mask to \Omega_{obj} \f$
    if (this%has_mask) then
       call mask_exterior_const(work, this%mask, 0.0_rp)
    end if
    ! integrate
    if (NEKO_BCKND_DEVICE .eq. 1) then
       this%value = device_glsc2(work%x_d, this%coef%B_d, n)
    else
       this%value = glsc2(work%x, this%coef%B, n)
    end if

    ! \f$ scale by 1/2 and 1/|\Omega_{obj}| \f$
    this%value = 0.5_rp * this%value / this%domain_volume

    ! relingush the scratch field
    call neko_scratch_registry%relinquish_field(temp_indices)
  end subroutine scalar_mixing_update_value

  !> Compute the sensitivity of the objective function with respect to the
  !! design
  !! @param this The object.
  !! @param design the design.
  subroutine scalar_mixing_update_sensitivity(this, design)
    class(scalar_mixing_objective_t), intent(inout) :: this
    class(design_t), intent(in) :: design

  end subroutine scalar_mixing_update_sensitivity

end module scalar_mixing_objective
