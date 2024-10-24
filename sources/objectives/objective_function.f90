! Copyright (c) 2024, The Neko Authors
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

 ! Implements the `objective_function_t` type.
module objective_function
  use num_types, only: rp, dp
  use field, only: field_t
  use fluid_scheme, only: fluid_scheme_t
  use adjoint_scheme, only: adjoint_scheme_t
  use topopt_design, only: topopt_design_t
  use case, only: case_t
  use adjoint_case, only: adjoint_case_t
  use point_zone, only: point_zone_t
  use dofmap, only: dofmap_t
  use point_zone_registry, only: neko_point_zone_registry
  implicit none
  private

  type, abstract, public :: objective_function_t
     !> objective function value
     real(kind=rp), public :: objective_function_value
     ! it may also be nice to have a list of the objective function value
     ! to append to instead?
     !real(kind=rp), public :: objective_function_value(:)

     !> A mask for where the objective function is evaluated
     class(point_zone_t), pointer :: mask
     !> containing a mask?
     logical :: if_mask

     ! TODO
     ! This is a bit strange,
     ! In my mind, an objective and a constraint are derrived from the
     ! same type and we would have separate instances depending on the problem
     !
     ! not all would need an adjoint forcing... eg a volume constraint.
     ! So this might be excessive..
     !
     ! TIM, if you can do this more neatly that would be great...
     ! So we have an objective_function class
     ! and then two derrived types,
     ! (1) things that involve fluid/passive scalars, ie typical objective
     ! (2) things that only require the design, eg, volume constraints
     !
     ! I'm also not so sure about the future, if it is necassarily
     ! a single "adjoint_forcing".
     ! It could maybe be forcing terms in both the velocity and passive
     ! scalar equations. Not sure yet.
     ! In princpal... it could also be strange BC's.
     !
     ! but I guess this all depends on the type of objective function,
     ! so it makes sense to derive more types.


     ! this is a field that holds the sensitivity to the coefficient.
     ! (again, we only have one for the Brinkman term, but in principle this
     ! should be a field list

     ! so this would hold dF/d\chi
     type(field_t), public :: sensitivity_to_coefficient


   contains
     !> init
     procedure, pass(this) :: init_base => objective_function_init_base
     procedure, pass(this) :: free_base => objective_function_free_base


     !> Compute the objective function
     ! note for TIM,
     ! this will REALLY need to be modified in the future...
     ! in a steady case, we just need to compute it on the last step
     ! in an unsteady case this will be a time integral
     ! so for now it's a proceedure.
     ! But it the future it may be a simulation component that will be
     ! appended to the fluid.
     !
     ! TODO
     ! this will need to be deffered in some way
     ! the init reads the JSON
     ! (or maybe we pass what objective function we have externally)
     ! based on the objective we have, we
     procedure(objective_function_compute), pass(this), deferred :: compute

     !> Compute the sensitivity
     ! again, right now this is just a procedure, but for unsteady cases this
     ! may require simulation components for time integrals.
     procedure(sensitivity_compute), pass(this), &
          deferred :: compute_sensitivity

     procedure(objective_function_init), pass(this), deferred :: init
     !> Destructor
     procedure(objective_function_free), pass(this), deferred :: free
  end type objective_function_t

  abstract interface
     subroutine objective_function_init(this, design, primal, adjoint)
       import objective_function_t, case_t, adjoint_case_t, &
            topopt_design_t
       class(objective_function_t), intent(inout) :: this
       ! these ones are inout because we may need to append source terms etc
       ! TODO
       ! At first I thought we could just pass the fluid, but if we have 
       ! passive scalars it's better to pass the whole case down.
       class(case_t), intent(inout) :: primal
       class(adjoint_case_t), intent(inout) :: adjoint
       ! TODO
       ! these should all be `class(design_variable)` in the future
       type(topopt_design_t), intent(inout) :: design
     end subroutine objective_function_init
  end interface

  abstract interface
     subroutine objective_function_compute(this, design, primal)
       import objective_function_t, case_t, topopt_design_t
       class(objective_function_t), intent(inout) :: this
       class(case_t), intent(in) :: primal
       type(topopt_design_t), intent(inout) :: design
     end subroutine objective_function_compute
  end interface

  abstract interface
     subroutine sensitivity_compute(this, design, primal, adjoint)
       import objective_function_t, case_t, adjoint_case_t, &
            topopt_design_t
       class(objective_function_t), intent(inout) :: this
       class(case_t), intent(in) :: primal
       class(adjoint_case_t), intent(in) :: adjoint
       type(topopt_design_t), intent(inout) :: design
     end subroutine sensitivity_compute
  end interface

  abstract interface
     subroutine objective_function_free(this)
       import objective_function_t
       class(objective_function_t), intent(inout) :: this
     end subroutine objective_function_free
  end interface

contains


  subroutine objective_function_init_base(this, dm_Xh, if_mask, mask_name)
    class(objective_function_t), target, intent(inout) :: this
    type(dofmap_t) :: dm_Xh !< Dofmap associated with \f$ X_h \f$
    character(len=*), intent(in), optional :: mask_name
    logical, intent(in) :: if_mask
    ! not sure what we need here yet

    ! initialize sensitivity field
    ! TODO
    ! Think about field lists in the future!
    call this%sensitivity_to_coefficient%init(dm_Xh, "design_indicator")

    this%if_mask = if_mask
    if (this%if_mask) then
       this%mask => &
       neko_point_zone_registry%get_point_zone(mask_name)
    end if

  end subroutine objective_function_init_base



  subroutine objective_function_free_base(this)
    class(objective_function_t), target, intent(inout) :: this
    ! not sure what we need here yet

  end subroutine objective_function_free_base
end module objective_function

