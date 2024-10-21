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

 ! Implements the `topopt_design_t` type.
module topopt_design
  use num_types, only: rp, sp
  use field, only: field_t
  use json_module, only: json_file
  use mapping, only: mapping_t
  use PDE_filter, only: PDE_filter_t
  use RAMP_mapping, only: RAMP_mapping_t
  use coefs, only: coef_t
  use scratch_registry, only: neko_scratch_registry
  use fld_file_output, only: fld_file_output_t

  implicit none
  private

  !> A topology optimization design variable
  type, extends(design_t), public :: topopt_design_t

     ! TODO
     ! in the future make this a derived type of a `design_variable`
     ! type, public, extends(design_variable_t) :: topopt_design_t
     !> the unfilitered design
     type(field_t), public :: design_indicator

     !> the mapped coefficient (Brinkman term)
     ! TODO
     ! NOTE: Tim, right now we only have the brinkman term
     ! in the future we may map to other coeeficients for other equations...
     ! in which case this should be a field list
     !
     ! or as I describe below, we also have multiple constraints,
     ! so a list-of-lists may be the correct way forward
     type(field_t), public :: brinkman_amplitude

     ! NOTE:
     ! again, we have to be so clear with nomenclature.
     ! If we have an objective function F.
     ! I also like to use \rho to describe the design_indicator
     !
     ! and let's say, for completness, we're doing conjugate heat transfer,
     ! so we have to map
     ! to 3 coefficients, \chi, C and \kappa.
     !
     ! then we will perform 3 mapping,
     ! \rho -> \chi
     ! \rho -> C
     ! \rho -> \kappa
     !
     ! Then during the forward/adjoint looping there will be an additional
     ! object, the "objective_function" object that will be responsible for
     ! computing the sensitivity of the objective function with respect to the
     ! coefficients
     ! ie,
     ! dF/d\chi, dF/dC and dF/d\kappa
     !
     ! What I'm calling "sensitivity" here, is the sensitivity with respect to
     ! the design indicator
     ! so dF/d\rho
     !
     ! so the proceedure "map_backwards" will take in the field list
     ! dF/d\chi, dF/dC and dF/d\kappa
     !and chain rule it's way back to
     ! dF/d\rho
     ! and store it here          v
     type(field_t), public :: sensitivity
     ! have a field list here
     ! type(filed_list_t), public :: constraint_sensitivity
     ! HOWEVER !
     ! What is a bit confusing to me... is how we'll deal with constraints.
     !
     ! Because in principle we could have constraints C1, C2, C3 etc
     ! Implying we also want dC1/d\rho, dC2/d\rho etc
     ! So this "sensitivity" should also be a list...
     !
     ! So I bet you'll have a nice abstract way of constructing this in the
     ! future but I think a sort of list-of-lists would be nice.
     !
     ! For F:
     ! \rho -> \tild{\rho} -> \chi
     ! \rho -> \tild{\rho} -> C
     ! \rho -> \tild{\rho} -> \kappa
     !
     ! For C1:
     ! \rho -> \tild{\rho}  (eg for a volume constraint or so)
     !
     ! For C2:
     ! \rho -> \tild{\rho}  (for something else)
     !
     ! etc..
     !
     ! So now we have multiple instances of the "objective" type,
     ! each for F, C1, C2 etc
     ! each of them can pass their dF\d_coefficents to the design
     ! and we continue from there.
     !
     ! perhaps even the "objective" type is defined as a list of objectives.

     ! Let's say way have a chain of two mappings
     type(PDE_filter_t) :: filter
     type(RAMP_mapping_t) :: mapping

     ! and we need to hold onto a field for the chain of mappings
     type(field_t) :: filtered_design


     !
     ! TODO
     ! you also had lots of masks etc that was a nice idea,
     ! but we'll cross that bridge later

     ! TODO
     ! you also had logicals for convergence etc,
     ! I feel they should live in "problem"

     ! afield writer would be nice too
     type(fld_file_output_t), private :: output


   contains
     !> init (will make this legit at some point)
     procedure, pass(this) :: init => topopt_design_init
     !> map (this will include everything from mapping
     ! design_indicator -> filtering -> chi
     ! and ultimately handle mapping different coeficients!
     procedure, pass(this) :: map_forward => topopt_design_map_forward
     !> this will contain chain rule for going backwards
     ! d_design_indicator <- d_filtering <- d_chi
     ! and ultimately handle mapping different coeficients!
     procedure, pass(this) :: map_backward => topopt_design_map_backward
     ! TODO
     ! maybe it would have been smarter to have a "coeficient" type,
     ! which is just a scalar field and set of mappings going from
     ! design_indicator -> coeficient and their corresponding chain rules
     ! maybe also some information about what equation they live in...

     ! a sampler being called from outside would be nice
     procedure, pass(this) :: sample => topopt_design_sample

     !> Destructor
     procedure, pass(this) :: free => topopt_design_free
     ! TODO
     ! I'm not sure who owns the optimizer...
     ! but it would make sense to have it in here so you provide it
     ! with dF/d_design_indicator and it updates itself.
     ! procedure, pass(this) :: update => topopt_design_update
  end type topopt_design_t


contains

  subroutine topopt_design_init(this, json, coef)
    class(topopt_design_t), target, intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(coef_t), intent(inout) :: coef
    integer :: n, i
    ! init the fields
    call this%design_indicator%init(coef%dof, "design_indicator")
    call this%brinkman_amplitude%init(coef%dof, "brinkman_amplitude")
    call this%sensitivity%init(coef%dof, "sensitivity")
    call this%filtered_design%init(coef%dof, "filtered_design")

    ! TODO
    ! this is where we steal basically everything in
    ! brinkman_source_term regarding loading initial fields
    ! for now, make it a cylinder by hand
    this%design_indicator = 0.0_rp
    this%brinkman_amplitude = 0.0_rp

    n = this%design_indicator%dof%size()
    do i = 1, n
       if (sqrt((this%design_indicator%dof%x(i,1,1,1) - 0.5_rp)**2 + &
            (this%design_indicator%dof%y(i,1,1,1) &
            - 0.5_rp)**2) .lt. 0.25_rp) then
          this%design_indicator%x(i,1,1,1) = 1.0_rp
       end if
    end do

    ! TODO
    ! we would also need to make a mapping type that reads in
    ! parameters etc about filtering and mapping
    ! ie,
    ! call mapper%init(this woud be from JSON)
    call this%filter%init(json, coef)
    call this%mapping%init(json, coef)

    ! and then we would map for the first one
    call this%map_forward()


    ! a field writer would be nice to output
    ! - design indicator (\rho)
    ! - mapped design (\chi)
    ! - sensitivity (dF/d\chi)
    ! TODO
    ! do this properly with JSON
    ! TODO
    ! obviously when we do the mappings properly, to many coeficients,
    ! we'll also have to modify this
    call this%output%init(sp, 'design', 3)
    call this%output%fields%assign(1, this%design_indicator)
    call this%output%fields%assign(2, this%brinkman_amplitude)
    call this%output%fields%assign(3, this%sensitivity)


  end subroutine topopt_design_init


  subroutine topopt_design_map_forward(this)
    class(topopt_design_t), target, intent(inout) :: this


    ! TODO
    ! this should be somehow deffered so we can pick different mappings!!!
    ! so this would be:
    ! call mapper%forward(fld_out, fld_in)

    call this%filter%apply_forward(this%filtered_design, &
         this%design_indicator)

    call this%mapping%apply_forward(this%brinkman_amplitude, &
         this%filtered_design)


  end subroutine topopt_design_map_forward

  subroutine topopt_design_map_backward(this, df_dchi)
    class(topopt_design_t), target, intent(inout) :: this
    type(field_t), intent(in) :: df_dchi
    type(field_t), pointer :: dF_dfiltered_design
    integer :: temp_indices(1)

    ! TODO
    ! again..
    ! so this would be:
    ! call mapper%backward(fld_out, fld_in)
    call neko_scratch_registry%request_field(dF_dfiltered_design, &
         temp_indices(1))

    call this%mapping%apply_backward(dF_dfiltered_design, df_dchi, &
         this%filtered_design)

    call this%filter%apply_backward(this%sensitivity, dF_dfiltered_design, &
         this%filtered_design)

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine topopt_design_map_backward

  subroutine topopt_design_free(this)
    class(topopt_design_t), target, intent(inout) :: this
    call this%brinkman_amplitude%free()
    call this%design_indicator%free()

  end subroutine topopt_design_free

  subroutine topopt_design_sample(this,t)
    class(topopt_design_t), target, intent(inout) :: this
    real(kind=rp), intent(in) :: t

    call this%output%sample(t)

  end subroutine topopt_design_sample
end module topopt_design
