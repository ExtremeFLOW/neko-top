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
!> Implements the `augmented_lagrangian_objective_t` type.
module augmented_lagrangian_objective
  use num_types, only: rp
  use field, only: field_t
  use field_math, only: field_col3, field_addcol3, field_cmult
  use scratch_registry, only: neko_scratch_registry
  use objective, only: objective_t
  use simulation_m, only: simulation_t
  use neko_config, only: NEKO_BCKND_DEVICE
  use math, only: copy
  use device_math, only: device_copy
  use design, only: design_t
  use json_module, only: json_file
  use json_utils, only: json_get_or_default
  implicit none
  private

  !> An objective function implementing our augmented lagrangian sensitivity
  !! contribution.
  type, public, extends(objective_t) :: augmented_lagrangian_objective_t
     private

     !> Pointer to the u field.
     type(field_t), pointer :: u => null()
     !> Pointer to the v field.
     type(field_t), pointer :: v => null()
     !> Pointer to the w field.
     type(field_t), pointer :: w => null()

     !> Pointer to adjoint u field.
     type(field_t), pointer :: adjoint_u => null()
     !> Pointer to adjoint v field.
     type(field_t), pointer :: adjoint_v => null()
     !> Pointer to adjoint w field.
     type(field_t), pointer :: adjoint_w => null()

   contains
     !> The common constructor using a JSON object.
     procedure, public, pass(this) :: init_json_sim => &
          augmented_lagrangian_init_json_sim
     !> The direct initializer from attributes.
     procedure, public, pass(this) :: init_from_attributes => &
          augmented_lagrangian_init_attributes
     !> Destructor.
     procedure, public, pass(this) :: free => augmented_lagrangian_free
     !> Computes the value of the objective function.
     procedure, public, pass(this) :: update_value => &
          augmented_lagrangian_update_value
     !> Computes the sensitivity with respect to the coefficient \f$\chi\f$.
     procedure, public, pass(this) :: update_sensitivity => &
          augmented_lagrangian_update_sensitivity

  end type augmented_lagrangian_objective_t

contains

  !> The common constructor using a JSON object.
  !! @param this the objective.
  !! @param json the JSON object.
  !! @param design the design.
  !! @param simulation the simulation.
  subroutine augmented_lagrangian_init_json_sim(this, json, design, simulation)
    class(augmented_lagrangian_objective_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(design_t), intent(in) :: design
    type(simulation_t), target, intent(inout) :: simulation

    character(len=:), allocatable :: name
    character(len=:), allocatable :: mask_name
    real(kind=rp) :: weight

    call json_get_or_default(json, "weight", weight, 1.0_rp)
    call json_get_or_default(json, "mask_name", mask_name, "")
    call json_get_or_default(json, "name", name, "Augmented Lagrangian")

    call this%init_from_attributes(design, simulation, weight, name, mask_name)
  end subroutine augmented_lagrangian_init_json_sim

  !> The actual constructor.
  !! @param this the objective.
  !! @param design the design.
  !! @param simulation the simulation.
  !! @param weight the weight of the objective function.
  !! @param name the name of the objective.
  !! @param mask_name the name of the mask.
  subroutine augmented_lagrangian_init_attributes(this, design, simulation, &
       weight, name, mask_name)
    class(augmented_lagrangian_objective_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(simulation_t), target, intent(inout) :: simulation
    real(kind=rp), intent(in) :: weight
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: mask_name

    call this%init_base(name, design%size(), weight, mask_name)

    ! Save the simulation and design
    this%u => simulation%neko_case%fluid%u
    this%v => simulation%neko_case%fluid%v
    this%w => simulation%neko_case%fluid%w
    this%adjoint_u => simulation%adjoint_case%fluid_adj%u_adj
    this%adjoint_v => simulation%adjoint_case%fluid_adj%v_adj
    this%adjoint_w => simulation%adjoint_case%fluid_adj%w_adj

  end subroutine augmented_lagrangian_init_attributes

  !> Destructor.
  subroutine augmented_lagrangian_free(this)
    class(augmented_lagrangian_objective_t), intent(inout) :: this
    call this%free_base()

    if (associated(this%u)) nullify(this%u)
    if (associated(this%v)) nullify(this%v)
    if (associated(this%w)) nullify(this%w)

    if (associated(this%adjoint_u)) nullify(this%adjoint_u)
    if (associated(this%adjoint_v)) nullify(this%adjoint_v)
    if (associated(this%adjoint_w)) nullify(this%adjoint_w)

  end subroutine augmented_lagrangian_free

  !> Compute the objective function.
  !! @param this the objective.
  !! @param design the design.
  subroutine augmented_lagrangian_update_value(this, design)
    class(augmented_lagrangian_objective_t), intent(inout) :: this
    class(design_t), intent(in) :: design

  end subroutine augmented_lagrangian_update_value

  !> update_value the sensitivity of the objective function with respect to \f\f$\chi\f\f$
  !! @param this the objective.
  !! @param design the design.
  subroutine augmented_lagrangian_update_sensitivity(this, design)
    class(augmented_lagrangian_objective_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(field_t), pointer :: work
    integer :: temp_indices(1)

    call neko_scratch_registry%request_field(work, temp_indices(1))

    ! here it should just be an inner product between the forward and adjoint
    call field_col3(work, this%u, this%adjoint_u)
    call field_addcol3(work, this%v, this%adjoint_v)
    call field_addcol3(work, this%w, this%adjoint_w)
    ! but negative
    call field_cmult(work, -1.0_rp)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(this%sensitivity%x_d, work%x_d, this%sensitivity%size())
    else
       call copy(this%sensitivity%x, work%x, this%sensitivity%size())
    end if

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine augmented_lagrangian_update_sensitivity

end module augmented_lagrangian_objective
