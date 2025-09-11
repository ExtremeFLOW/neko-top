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
  use math, only: copy, col3, addcol3, col2, invcol2
  use device_math, only: device_copy, device_col3, device_addcol3, &
      device_col2, device_invcol2
  use design, only: design_t
  use json_module, only: json_file
  use json_utils, only: json_get_or_default
  use interpolation, only: interpolator_t
  use space, only: space_t, GL
  use coefs, only: coef_t
  use, intrinsic :: iso_c_binding, only: c_ptr, C_NULL_PTR
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

     ! ---- everything GLL ----
     !> The original space used in the simulation
     type(space_t), pointer :: Xh_GLL
     !> cfs of the original space in the simulation
     type(coef_t), pointer :: c_Xh_GLL

     ! ---- everything for GL ----
     !> The additional higher-order space used in dealiasing
     type(space_t), pointer :: Xh_GL
     !> cfs of the higher-order space
     type(coef_t), pointer :: c_Xh_GL
     !> Interpolator between the original and higher-order spaces
     type(interpolator_t), pointer :: GLL_to_GL
     !> If dealiasing should be applied
     logical :: dealias

     !> work arrays
     real(kind=rp), allocatable :: accumulate(:), fld_GL(:), adjoint_fld_GL(:)
     !> Device pointer for `accumulate`
     type(c_ptr) :: accumulate_d = C_NULL_PTR
     !> Device pointer for `fld_GL`
     type(c_ptr) :: fld_GL_d = C_NULL_PTR
     !> Device pointer for `adjoint_fld_GL`
     type(c_ptr) :: adjoint_fld_GL_d = C_NULL_PTR


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
    logical :: dealias

    call json_get_or_default(json, "weight", weight, 1.0_rp)
    call json_get_or_default(json, "mask_name", mask_name, "")
    call json_get_or_default(json, "name", name, "Augmented Lagrangian")
    call json_get_or_default(json, "dealias", dealias, .true.)

    call this%init_from_attributes(design, simulation, weight, name, &
         mask_name, dealias)
  end subroutine augmented_lagrangian_init_json_sim

  !> The actual constructor.
  !! @param this the objective.
  !! @param design the design.
  !! @param simulation the simulation.
  !! @param weight the weight of the objective function.
  !! @param name the name of the objective.
  !! @param mask_name the name of the mask.
  !! @param dealias should dealiasing be applied.
  subroutine augmented_lagrangian_init_attributes(this, design, simulation, &
       weight, name, mask_name, dealias)
    class(augmented_lagrangian_objective_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(simulation_t), target, intent(inout) :: simulation
    real(kind=rp), intent(in) :: weight
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: mask_name
    logical, intent(in) :: dealias
    integer :: nel, n_GL

    call this%init_base(name, design%size(), weight, mask_name)

    ! Save the simulation and design
    this%u => simulation%neko_case%fluid%u
    this%v => simulation%neko_case%fluid%v
    this%w => simulation%neko_case%fluid%w
    this%adjoint_u => simulation%adjoint_case%fluid_adj%u_adj
    this%adjoint_v => simulation%adjoint_case%fluid_adj%v_adj
    this%adjoint_w => simulation%adjoint_case%fluid_adj%w_adj

    this%dealias = dealias
    ! GLL
    this%c_Xh_GLL => simulation%neko_case%fluid%c_Xh
    this%Xh_GLL => this%c_Xh_GLL%Xh

    ! GL
    this%c_Xh_GL => simulation%adjoint_case%fluid_adj%c_Xh_GL
    this%Xh_GL => this%c_Xh_GL%Xh

    ! GLL to GL
    this%GLL_to_GL => simulation%adjoint_case%fluid_adj%GLL_to_GL

    if (this%dealias) then
    ! allocate work arrays for dealiasing
    nel = this%c_Xh_GLL%msh%nelv
    n_GL = nel * this%Xh_GL%lxyz
    allocate(this%accumulate(n_GL))
    allocate(this%fld_GL(n_GL))
    allocate(this%adjoint_fld_GL(n_GL))
    if (NEKO_BCKND_DEVICE .eq. 1) then
    call device_map(this%accumulate, this%accumulate_d, n_GL)
    call device_map(this%fld_GL, this%fld_GL_d, n_GL)
    call device_map(this%adjoint_fld_GL, this%adjoint_fld_GL_d, n_GL)
    end if
    end if

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
    real(kind=rp), dimension(this%Xh_GL%lxyz * this%c_Xh_GLL%msh%nelv) :: &
         accumulate, fld_GL, adjoint_fld_GL
    integer :: n_GL, nel

    call neko_scratch_registry%request_field(work, temp_indices(1))

    if (this%dealias) then

       nel = this%c_Xh_GLL%msh%nelv
       n_GL = nel * this%Xh_GL%lxyz
       if (NEKO_BCKND_DEVICE .eq. 1) then

          call this%GLL_to_GL%map(this%fld_GL, this%u%x, nel, this%Xh_GL)
          call this%GLL_to_GL%map(this%adjoint_fld_GL, this%adjoint_u%x, nel, this%Xh_GL)
          call device_col3(this%accumulate_d, this%fld_G_dL, this%adjoint_fld_GL_d, n_GL)

          call this%GLL_to_GL%map(this%fld_GL, this%v%x, nel, this%Xh_GL)
          call this%GLL_to_GL%map(this%adjoint_fld_GL, this%adjoint_v%x, nel, this%Xh_GL)
          call device_addcol3(this%accumulate_d, this%fld_GL_d, this%adjoint_fld_GL_d, n_GL)

          call this%GLL_to_GL%map(this%fld_GL, this%w%x, nel, this%Xh_GL)
          call this%GLL_to_GL%map(this%adjoint_fld_GL, this%adjoint_w%x, nel, this%Xh_GL)
          call device_addcol3(this%accumulate_d, this%fld_GL_d, adjoint_fld_GL_d, n_GL)

          ! multiply by GL mass matrix
          call device_col2(this%accumulate_d, this%c_Xh_GL%B_d, n_GL)
          ! map back to GLL
          call this%GLL_to_GL%map(work%x, this%accumulate, nel, this%Xh_GLL)
          ! preempt the GLL mass matrix
          call invcol2(work%x_d, this%c_Xh_GLL%B_d, work%size())

          ! but negative
          call field_cmult(work, -1.0_rp)

       else

          call this%GLL_to_GL%map(this%fld_GL, this%u%x, nel, this%Xh_GL)
          call this%GLL_to_GL%map(this%adjoint_fld_GL, this%adjoint_u%x, nel, this%Xh_GL)
          call col3(this%accumulate, this%fld_GL, this%adjoint_fld_GL, n_GL)

          call this%GLL_to_GL%map(this%fld_GL, this%v%x, nel, this%Xh_GL)
          call this%GLL_to_GL%map(this%adjoint_fld_GL, this%adjoint_v%x, nel, this%Xh_GL)
          call addcol3(this%accumulate, this%fld_GL, this%adjoint_fld_GL, n_GL)

          call this%GLL_to_GL%map(this%fld_GL, this%w%x, nel, this%Xh_GL)
          call this%GLL_to_GL%map(this%adjoint_fld_GL, this%adjoint_w%x, nel, this%Xh_GL)
          call addcol3(this%accumulate, this%fld_GL, adjoint_fld_GL, n_GL)

          ! multiply by GL mass matrix
          call col2(this%accumulate, this%c_Xh_GL%B, n_GL)
          ! map back to GLL
          call this%GLL_to_GL%map(work%x, this%accumulate, nel, this%Xh_GLL)
          ! preempt the GLL mass matrix
          call invcol2(work%x, this%c_Xh_GLL%B, work%size())

          ! but negative
          call field_cmult(work, -1.0_rp)

       end if
    else
       call field_col3(work, this%u, this%adjoint_u)
       call field_addcol3(work, this%v, this%adjoint_v)
       call field_addcol3(work, this%w, this%adjoint_w)
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(this%sensitivity%x_d, work%x_d, this%sensitivity%size())
    else
       call copy(this%sensitivity%x, work%x, this%sensitivity%size())
    end if

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine augmented_lagrangian_update_sensitivity

end module augmented_lagrangian_objective
