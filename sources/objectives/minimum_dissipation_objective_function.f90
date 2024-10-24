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
!> Implements the `minimum_dissipation_objective_function_t` type.
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
module minimum_dissipation_objective_function
  use num_types, only : rp
  use field, only: field_t
  use field_math, only: field_col3, field_addcol3, field_cmult, field_add2s2
  use operators, only: grad
  use scratch_registry, only : neko_scratch_registry
  use adjoint_minimum_dissipation_source_term, only: &
       adjoint_minimum_dissipation_source_term_t
  use objective_function, only : objective_function_t
  use fluid_scheme, only : fluid_scheme_t
  use adjoint_scheme, only : adjoint_scheme_t
  use math, only : glsc2
  use topopt_design, only: topopt_design_t
  use adjoint_lube_source_term, only: adjoint_lube_source_term_t
  use point_zone, only: point_zone_t
  use mask_ops, only: mask_exterior_const
  use math_ext, only: glsc2_mask
  implicit none
  private

  !> An objective function corresponding to minimum dissipation
  ! $ F =  \int_\Omega |\nabla u|^2 d \Omega + K \int_Omega \frac{1}{2} \chi
  ! |\mathbf{u}|^2 d \Omega $
  type, public, extends(objective_function_t) :: &
       minimum_dissipation_objective_function_t
     real(kind=rp) :: K, dissipation, lube_value
     logical :: if_lube

     ! TODO
     ! this is just for testing!
     ! actually rescaling the adjoint is a bit more involved,
     ! and we have to be careful of the brinkman term
     !> A scaling factor
     real(kind=rp) :: obj_scale

   contains
     !> The common constructor using a JSON object.
     procedure, pass(this) :: init => &
          minimum_dissipation_objective_function_init
     !> Destructor.
     procedure, pass(this) :: free => &
          minimum_dissipation_objective_function_free
     !> Computes the value of the objective function.
     procedure, pass(this) :: compute => &
          minimum_dissipation_objective_function_compute
     !> Computes the sensitivity with respect to the coefficent $\chi$.
     procedure, pass(this) :: compute_sensitivity => &
          minimum_dissipation_objective_function_compute_sensitivity
  end type minimum_dissipation_objective_function_t

contains
  !> The common constructor using a JSON object.
  !! @param design the design.
  !! @param fluid the fluid scheme.
  !! @param adjoint the fluid adjoint.
  subroutine minimum_dissipation_objective_function_init(this, &
       design, fluid, adjoint)
    class(minimum_dissipation_objective_function_t), intent(inout) :: this
    class(fluid_scheme_t), intent(inout) :: fluid
    class(adjoint_scheme_t), intent(inout) :: adjoint
    type(topopt_design_t), intent(inout) :: design
    type(adjoint_minimum_dissipation_source_term_t) :: adjoint_forcing
    type(adjoint_lube_source_term_t) :: lube_term
    character(len=:), allocatable :: objective_location_zone_name
    logical :: if_mask

    ! here we would read from the JSON (or have something passed in)
    ! about the lube term
    this%if_lube = .true.
    this%K = 1.0_rp
    this%obj_scale = 1.0_rp
    !this%obj_scale = 0.00000001_rp


    ! mask would also be read from JSON... I'm hard coding
    if_mask = .false.
    objective_location_zone_name = "objective_location"

    call this%init_base(fluid%dm_Xh, if_mask, objective_location_zone_name)

    ! you will need to init this!
    ! append a source term based on the minimum dissipation
    ! init the adjoint forcing term for the adjoint
    call adjoint_forcing%init_from_components( &
         adjoint%f_adj_x, adjoint%f_adj_y, adjoint%f_adj_z, &
         fluid%u, fluid%v, fluid%w, this%obj_scale, &
         this%mask, this%if_mask, &
         adjoint%c_Xh)
    ! append adjoint forcing term based on objective function
    call adjoint%source_term%add_source_term(adjoint_forcing)


    ! if we have the lube term we need to initialize and append that too
    if (this%if_lube) then
       ! TODO
       ! make this allocatable and only allocate it if needed!
       ! or is that allready what's happening? Tim, y/n?
       call lube_term%init_from_components(&
            adjoint%f_adj_x, adjoint%f_adj_y, adjoint%f_adj_z, design, &
            this%k*this%obj_scale, &
            fluid%u, fluid%v, fluid%w, &
            this%mask, this%if_mask, &
            adjoint%c_Xh)
       ! append adjoint forcing term based on objective function
       call adjoint%source_term%add_source_term(lube_term)
    end if

  end subroutine minimum_dissipation_objective_function_init


  !> Destructor.
  subroutine minimum_dissipation_objective_function_free(this)
    class(minimum_dissipation_objective_function_t), intent(inout) :: this
    ! TODO
    ! you probably need to deallocate the source term!

    call this%free_base()
  end subroutine minimum_dissipation_objective_function_free

  !> Compute the objective function.
  !! @param design the design.
  !! @param fluid the fluid scheme.
  !! @param adjoint the fluid adjoint.
  subroutine minimum_dissipation_objective_function_compute(this, design, fluid)
    class(minimum_dissipation_objective_function_t), intent(inout) :: this
    class(fluid_scheme_t), intent(in) :: fluid
    type(topopt_design_t), intent(inout) :: design
    type(field_t), pointer :: wo1, wo2, wo3
    type(field_t), pointer :: objective_field
    integer :: temp_indices(4)
    integer n

    call neko_scratch_registry%request_field(wo1, temp_indices(1))
    call neko_scratch_registry%request_field(wo2, temp_indices(2))
    call neko_scratch_registry%request_field(wo3, temp_indices(3))
    call neko_scratch_registry%request_field(objective_field, temp_indices(4))

    ! compute the objective function.
    call grad(wo1%x, wo2%x, wo3%x, fluid%u%x, fluid%C_Xh)
    call field_col3(objective_field, wo1, wo1)
    call field_addcol3(objective_field, wo2, wo2)
    call field_addcol3(objective_field, wo3, wo3)

    call grad(wo1%x, wo2%x, wo3%x, fluid%v%x, fluid%C_Xh)
    call field_addcol3(objective_field, wo1, wo1)
    call field_addcol3(objective_field, wo2, wo2)
    call field_addcol3(objective_field, wo3, wo3)

    call grad(wo1%x, wo2%x, wo3%x, fluid%w%x, fluid%C_Xh)
    call field_addcol3(objective_field, wo1, wo1)
    call field_addcol3(objective_field, wo2, wo2)
    call field_addcol3(objective_field, wo3, wo3)

    ! integrate the field
    n = wo1%size()
    if (this%if_mask) then
       this%dissipation = glsc2_mask(objective_field%x, fluid%C_Xh%b, &
            n, this%mask%mask, this%mask%size)
    else
       this%dissipation = glsc2(objective_field%x, fluid%C_Xh%b, n)
    end if

    if (this%if_lube) then
       ! it's becoming so stupid to pass the whole fluid and adjoint and
       ! design through
       ! I feel like every objective function should have internal pointers to
       ! u,v,w and u_adj, v_adj, w_adj and perhaps the design
       ! (the whole design, so we get all the coeffients)
       call field_col3(objective_field, fluid%u, design%brinkman_amplitude)
       call field_addcol3(objective_field, fluid%v, design%brinkman_amplitude)
       call field_addcol3(objective_field, fluid%w, design%brinkman_amplitude)
       if (this%if_mask) then
          this%lube_value = glsc2_mask(objective_field%x, fluid%C_Xh%b, &
               n, this%mask%mask, this%mask%size)
       else
          this%lube_value = glsc2(objective_field%x, fluid%C_Xh%b, n)
       end if
       this%objective_function_value = this%dissipation &
            + 0.5*this%K*this%lube_value
    else
       this%objective_function_value = this%dissipation
    end if

    ! scale everything
    this%objective_function_value = this%objective_function_value*this%obj_scale

    !TODO
    ! GPUS

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine minimum_dissipation_objective_function_compute

  !> compute the sensitivity of the objective function with respect to $\chi$
  !! @param design the design.
  !! @param fluid the fluid scheme.
  !! @param adjoint the fluid adjoint.
  subroutine minimum_dissipation_objective_function_compute_sensitivity(this, &
       design, fluid, adjoint)
    class(minimum_dissipation_objective_function_t), intent(inout) :: this
    type(topopt_design_t), intent(inout) :: design
    class(fluid_scheme_t), intent(in) :: fluid
    class(adjoint_scheme_t), intent(in) :: adjoint
    type(field_t), pointer :: lube_contribution
    integer :: temp_indices(1)


    ! here it should just be an inner product between the forward and adjoint
    call field_col3(this%sensitivity_to_coefficient, fluid%u, adjoint%u_adj)
    call field_addcol3(this%sensitivity_to_coefficient, fluid%v, adjoint%v_adj)
    call field_addcol3(this%sensitivity_to_coefficient, fluid%w, adjoint%w_adj)
    ! but negative
    call field_cmult(this%sensitivity_to_coefficient, -1.0_rp)

    ! if we have the lube term we also get an extra term in the sensitivity
    ! K*u^2
    ! TODO
    ! omfg be so careful with non-dimensionalization etc
    ! I bet this is scaled a smidge wrong (ie, track if it's 1/2 or not etc)
    ! do this later

    if (this%if_lube) then
       call neko_scratch_registry%request_field(lube_contribution, &
            temp_indices(1))
       call field_col3(lube_contribution, fluid%u, fluid%u)
       call field_addcol3(lube_contribution, fluid%v, fluid%v)
       call field_addcol3(lube_contribution, fluid%w, fluid%w)
       ! fuck be careful with these scalaing!
       call field_add2s2(this%sensitivity_to_coefficient, lube_contribution, &
            this%K*this%obj_scale)
       call neko_scratch_registry%relinquish_field(temp_indices)
    end if

    ! I don't actually think you scale the sensitivity...
    ! because the adjoint field is already scaled
    !call field_cmult(this%sensitivity_to_coefficient, this%obj_scale)


  end subroutine minimum_dissipation_objective_function_compute_sensitivity

end module minimum_dissipation_objective_function
