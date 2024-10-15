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
!> Implements the `scalar_mixing_objective_function_t` type.
!
! This is an objective function akin to that used by Andreasen et al. 2008
!  https://doi.org/10.1002/fld.1964
!
! with some modification. In their paper this objective function was based on
! a surface integral at the outlet
! $\int_{\Gamma_{out}} (\phi - \bar{\phi})^2 d\Gamma$ (with some normalization)
!
! which results in an adjoint boundary condition on $\Gamma_{out}$
!
! Since, in the current state, adjoint boundary conditions are not handled so
! well, we can replace the surface integral with a volume integral close to 
! the outlet, resulting in an adjoint source term of the form:
!
! $ \phi - \bar{\phi} $
!
module scalar_mixing_objective_function
  use num_types, only : rp
  use field_list, only : field_list_t
  use json_module, only : json_file
  use json_utils, only: json_get, json_get_or_default
  use source_term, only : source_term_t
  use coefs, only : coef_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use utils, only : neko_error
  use field, only: field_t
  use topopt_design, only: topopt_design_t
  use field_math, only: field_col3, field_addcol3, field_cmult, field_add2s2, &
  field_cadd, field_copy, field_col2
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
  use adjoint_mixing_scalar_source_term, &
       only : adjoint_mixing_scalar_source_term_t
  use objective_function, only : objective_function_t
  use fluid_scheme, only : fluid_scheme_t
  use adjoint_scheme, only : adjoint_scheme_t
  use fluid_source_term, only: fluid_source_term_t
  use math, only : glsc2
  use topopt_design, only: topopt_design_t
  use adjoint_lube_source_term, only: adjoint_lube_source_term_t
  use case, only: case_t
  use adjoint_case, only: adjoint_case_t
  implicit none
  private

  !> An objective function corresponding to enhanced mixing of a passive scalar
  ! $F = \int_{\Omega_{out}} (\phi - \bar{\phi})^2 d\Omega$ 
  type, public, extends(objective_function_t) :: &
       scalar_mixing_objective_function_t

     ! TODO
     ! this is just for testing!
     ! actually rescaling the adjoint is a bit more involved,
     ! and we have to be careful of the brinkman term
     !> A scaling factor
     real(kind=rp) :: obj_scale

   contains
     !> The common constructor using a JSON object.
     procedure, pass(this) :: init => &
          scalar_mixing_objective_function_init
     !> Destructor.
     procedure, pass(this) :: free => &
          scalar_mixing_objective_function_free
     !> Computes the value of the objective function.
     procedure, pass(this) :: compute => &
          scalar_mixing_objective_function_compute
     !> Computes the sensitivity with respect to the coefficent $\chi$.
     procedure, pass(this) :: compute_sensitivity => &
          scalar_mixing_objective_function_compute_sensitivity
  end type scalar_mixing_objective_function_t

contains
  !> The common constructor using a JSON object.
  !! @param design the design.
  !! @param primal the forward case.
  !! @param adjoint the adjoint case.
  subroutine scalar_mixing_objective_function_init(this, &
       design, primal, adjoint)
    class(scalar_mixing_objective_function_t), intent(inout) :: this
    class(case_t), intent(inout) :: primal
    class(adjoint_case_t), intent(inout) :: adjoint
    type(topopt_design_t), intent(inout) :: design
    type(adjoint_mixing_scalar_source_term_t) :: adjoint_forcing
    type(adjoint_lube_source_term_t) :: lube_term

    ! here we would read from the JSON (or have something passed in)
    ! about the lube term
    this%obj_scale = 1.0_rp


    call this%init_base(primal%fluid%dm_Xh)

    ! the adjoint source term is applied to the scalar
    ! init the adjoint forcing term for the adjoint
    call adjoint_forcing%init_from_components( &
         adjoint%scalar%f_Xh, primal%scalar%s, & 
         this%obj_scale, adjoint%scheme%c_Xh)
    ! append adjoint forcing term based on objective function
    print *, "apending"
    call adjoint%scalar%source_term%add_source_term(adjoint_forcing)
    print *, "apended", allocated(adjoint%scalar%source_term%source_terms), &
    size(adjoint%scalar%source_term%source_terms)

  end subroutine scalar_mixing_objective_function_init


  !> Destructor.
  subroutine scalar_mixing_objective_function_free(this)
    class(scalar_mixing_objective_function_t), intent(inout) :: this
    ! TODO
    ! you probably need to deallocate the source term!

    call this%free_base()
  end subroutine scalar_mixing_objective_function_free

  !> Compute the objective function.
  !! @param design the design.
  !! @param primal the forward case.
  !! @param adjoint the adjoint case.
  subroutine scalar_mixing_objective_function_compute(this, design, primal)
    class(scalar_mixing_objective_function_t), intent(inout) :: this
    class(case_t), intent(in) :: primal
    type(topopt_design_t), intent(inout) :: design
    integer :: i
    type(field_t), pointer :: wo1, wo2, wo3
    type(field_t), pointer :: objective_field
    integer :: temp_indices(1)
    integer n
    real (kind=rp) :: average_upstream

    call neko_scratch_registry%request_field(wo1, temp_indices(1))
    n = wo1%size()
    ! ok remember, in the actual objective function we're going to need to 
    ! take an average value upstream with a surface integral, and use that
    ! as $\bar{\phi}$
    !
    ! $\int_\Omega (\phi - 0.5)^2 d\Omega$
    !
    ! For now we know it's 0.5
    average_upstream = 0.5_rp
    call field_copy(wo1, primal%scalar%s)
    call field_cadd(wo1, -average_upstream)
    ! square it
    call field_col2(wo1,wo1)
    ! integrate
    this%objective_function_value = glsc2(wo1%x,primal%fluid%C_Xh%b, n)
    ! TODO
    ! so much more, specifically
    ! - masks
    ! - calcuate upstream etc
    ! - ultimately, replace with a surface integral and adjoint BCs

    ! scale everything
    this%objective_function_value = this%objective_function_value*this%obj_scale

    !TODO
    ! GPUS (really only glsc2)

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine scalar_mixing_objective_function_compute

  !> compute the sensitivity of the objective function with respect to $\chi$
  !! @param design the design.
  !! @param primal the forward case.
  !! @param adjoint the adjoint case.
  subroutine scalar_mixing_objective_function_compute_sensitivity(this, &
       design, primal, adjoint)
    class(scalar_mixing_objective_function_t), intent(inout) :: this
    type(topopt_design_t), intent(inout) :: design
    class(case_t), intent(in) :: primal
    class(adjoint_case_t), intent(in) :: adjoint
    type(field_t), pointer :: lube_contribution
    integer :: temp_indices(1)


    ! here it should just be an inner product between the forward and adjoint
    call field_col3(this%sensitivity_to_coefficient, primal%fluid%u, &
    adjoint%scheme%u_adj)
    call field_addcol3(this%sensitivity_to_coefficient, primal%fluid%v, &
    adjoint%scheme%v_adj)
    call field_addcol3(this%sensitivity_to_coefficient, primal%fluid%w, &
    adjoint%scheme%w_adj)
    ! but negative
    call field_cmult(this%sensitivity_to_coefficient, -1.0_rp)

    ! NOTE
    ! this will be different if we map to more coefficients

  end subroutine scalar_mixing_objective_function_compute_sensitivity

end module scalar_mixing_objective_function
