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
!> Implements the `minimum_dissipation_objective_t` type.
!
! I promise I'll write this document properly in the future...
!
! But the Borval Peterson (I think) paper had an objective function
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
!      | dissipation |     |"lube term"|
!
! I say "lube term" because they said it came from lubrication theory...
! Anyway, we can change all this later (especially the names!)

! If the objective function \int |\nabla u|^2,
! the corresponding adjoint forcing is \int \nabla v \cdot \nabla u
!
! for the lube term, the adjoint forcing is \chi u
!
! This has always annoyed me...
! because now I see one objective and one constraint
!
module minimum_dissipation_objective
  use num_types, only: rp
  use field, only: field_t
  use field_math, only: field_col3, field_addcol3, field_cmult, field_add2s2, &
       field_copy
  use operators, only: grad
  use scratch_registry, only: neko_scratch_registry
  use adjoint_minimum_dissipation_source_term, only: &
       adjoint_minimum_dissipation_source_term_t
  use objective, only: objective_t
  use simulation, only: simulation_t
  use fluid_scheme_incompressible, only: fluid_scheme_incompressible_t
  use adjoint_scheme, only: adjoint_scheme_t
  use neko_config, only: NEKO_BCKND_DEVICE
  use math, only: glsc2, copy
  use device_math, only: device_copy, device_glsc2
  use design, only: design_t
  use topopt_design, only: topopt_design_t
  use adjoint_lube_source_term, only: adjoint_lube_source_term_t
  use point_zone, only: point_zone_t
  use mask_ops, only: mask_exterior_const
  use math_ext, only: glsc2_mask
  use utils, only: neko_error
  use json_module, only: json_file
  use json_utils, only: json_get_or_default
  implicit none
  private

  !> An objective function corresponding to minimum dissipation
  ! $ F =  \int_\Omega |\nabla u|^2 d \Omega + K \int_Omega \frac{1}{2} \chi
  ! |\mathbf{u}|^2 d \Omega $
  type, public, extends(objective_t) :: minimum_dissipation_objective_t
     private

     class(fluid_scheme_incompressible_t), pointer :: fluid
     class(adjoint_scheme_t), pointer :: adjoint

   contains
     !> The common constructor using a JSON object.
     procedure, public, pass(this) :: init_json => minimum_dissipation_init_json
     !> The direct initializer from attributes.
     procedure, public, pass(this) :: init_from_attributes => &
          minimum_dissipation_init_attributes
     !> Destructor.
     procedure, public, pass(this) :: free => minimum_dissipation_free
     !> Computes the value of the objective function.
     procedure, public, pass(this) :: update_value => &
          minimum_dissipation_update_value
     !> Computes the sensitivity with respect to the coefficient $\chi$.
     procedure, public, pass(this) :: update_sensitivity => &
          minimum_dissipation_update_sensitivity

  end type minimum_dissipation_objective_t

contains

  !> The common constructor using a JSON object.
  !! @param design the design.
  !! @param fluid the fluid scheme.
  !! @param adjoint the fluid adjoint.
  subroutine minimum_dissipation_init_json(this, json, design, simulation)
    class(minimum_dissipation_objective_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(design_t), intent(in) :: design
    type(simulation_t), target, intent(inout) :: simulation
    character(len=:), allocatable :: name
    character(len=:), allocatable :: mask_name
    real(kind=rp) :: weight

    call json_get_or_default(json, "weight", weight, 1.0_rp)
    call json_get_or_default(json, "mask_name", mask_name, "")
    call json_get_or_default(json, "name", name, "Dissipation")

    call this%init_from_attributes(design, simulation, weight, name, mask_name)
  end subroutine minimum_dissipation_init_json

!> The actual constructor.
!! @param design the design.
!! @param simulation the simulation.
!! @param weight the weight of the objective function.
!! @param mask_name the name of the mask.
  subroutine minimum_dissipation_init_attributes(this, design, simulation, &
       weight, name, mask_name)
    class(minimum_dissipation_objective_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(simulation_t), target, intent(inout) :: simulation
    real(kind=rp), intent(in) :: weight
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: mask_name

    type(adjoint_minimum_dissipation_source_term_t) :: adjoint_forcing

    call this%init_base(name, design%size(), weight, mask_name)

    ! Save the simulation and design
    this%fluid => simulation%fluid_scheme
    this%adjoint => simulation%adjoint_case%scheme

    ! you will need to init this!
    ! append a source term based on the minimum dissipation
    ! init the adjoint forcing term for the adjoint
    call adjoint_forcing%init_from_components( &
         this%adjoint%f_adj_x, this%adjoint%f_adj_y, this%adjoint%f_adj_z, &
         this%fluid%u, this%fluid%v, this%fluid%w, this%weight, &
         this%mask, this%has_mask, &
         this%adjoint%c_Xh)
    ! append adjoint forcing term based on objective function
    call this%adjoint%source_term%add_source_term(adjoint_forcing)

  end subroutine minimum_dissipation_init_attributes

  !> Destructor.
  subroutine minimum_dissipation_free(this)
    class(minimum_dissipation_objective_t), intent(inout) :: this
    call this%free_base()

    if (associated(this%fluid)) nullify(this%fluid)
    if (associated(this%adjoint)) nullify(this%adjoint)

  end subroutine minimum_dissipation_free

  !> Compute the objective function.
  !! @param design the design.
  !! @param fluid the fluid scheme.
  !! @param adjoint the fluid adjoint.
  subroutine minimum_dissipation_update_value(this, design)
    class(minimum_dissipation_objective_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(field_t), pointer :: wo1, wo2, wo3, work
    type(field_t), pointer :: objective_field
    integer :: temp_indices(5)
    integer n

    call neko_scratch_registry%request_field(wo1, temp_indices(1))
    call neko_scratch_registry%request_field(wo2, temp_indices(2))
    call neko_scratch_registry%request_field(wo3, temp_indices(3))
    call neko_scratch_registry%request_field(objective_field, temp_indices(4))
    call neko_scratch_registry%request_field(work, temp_indices(5))

    ! update_value the objective function.
    call grad(wo1%x, wo2%x, wo3%x, this%fluid%u%x, this%fluid%c_Xh)
    call field_col3(objective_field, wo1, wo1)
    call field_addcol3(objective_field, wo2, wo2)
    call field_addcol3(objective_field, wo3, wo3)

    call grad(wo1%x, wo2%x, wo3%x, this%fluid%v%x, this%fluid%c_Xh)
    call field_addcol3(objective_field, wo1, wo1)
    call field_addcol3(objective_field, wo2, wo2)
    call field_addcol3(objective_field, wo3, wo3)

    call grad(wo1%x, wo2%x, wo3%x, this%fluid%w%x, this%fluid%c_Xh)
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
          this%value = device_glsc2(work%x_d, this%fluid%c_xh%B_d, n)
       else
          this%value = glsc2_mask(objective_field%x, this%fluid%C_Xh%b, &
               n, this%mask%mask, this%mask%size)
       end if
    else
       if (neko_bcknd_device .eq. 1) then
          this%value = device_glsc2(objective_field%x_d, &
               this%fluid%C_Xh%b_d, n)
       else
          this%value = glsc2(objective_field%x, this%fluid%C_Xh%b, n)
       end if
    end if

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine minimum_dissipation_update_value

  !> update_value the sensitivity of the objective function with respect to $\chi$
  !! @param design the design.
  subroutine minimum_dissipation_update_sensitivity(this, design)
    class(minimum_dissipation_objective_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(field_t), pointer :: work
    integer :: temp_indices(1)

    call neko_scratch_registry%request_field(work, temp_indices(1))

    ! here it should just be an inner product between the forward and adjoint
    call field_col3(work, this%fluid%u, this%adjoint%u_adj)
    call field_addcol3(work, this%fluid%v, this%adjoint%v_adj)
    call field_addcol3(work, this%fluid%w, this%adjoint%w_adj)
    ! but negative
    call field_cmult(work, -1.0_rp)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(this%sensitivity%x_d, work%x_d, this%sensitivity%size())
    else
       call copy(this%sensitivity%x, work%x, this%sensitivity%size())
    end if

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine minimum_dissipation_update_sensitivity

end module minimum_dissipation_objective
