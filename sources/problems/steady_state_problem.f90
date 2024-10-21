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
!> Implements the `steady_problem_t` type.
! Here, we simply march forward to steady state solutions
module steady_state_problem
  use num_types, only: rp, sp
  use problem, only: problem_t
  use adjoint_case, only: adjoint_case_t, adjoint_init, adjoint_free
  use simulation_adjoint, only: solve_adjoint
  use neko, only: neko_init, neko_finalize, neko_solve
  use case, only: case_t
  use topopt_design, only: topopt_design_t
  use json_file_module, only: json_file
  use minimum_dissipation_objective_function, only: &
       minimum_dissipation_objective_function_t
  use volume_constraint, only: volume_constraint_t
  use fld_file_output, only: fld_file_output_t
  use steady_simcomp, only: steady_simcomp_t
  use simple_brinkman_source_term, only: simple_brinkman_source_term_t
  use json_utils_ext, only: json_key_fallback, json_get_subdict
  use json_value_module, only: json_value
  use utils, only: neko_error
  use user_intf, only: simulation_component_user_settings
  use simcomp_executor, only: neko_simcomps
  implicit none
  private

  !> To compute a steady state problem
  type, public, extends(problem_t) :: steady_state_problem_t

     !> and primal case
     type(case_t), public :: C
     !> and adjoint case
     type(adjoint_case_t), public :: adj

     type(topopt_design_t), pointer :: design

     !> TODO
     ! we need a `objective_list` which is allocatable and contains a factory
     ! to fill itself up with from the JSON
     ! for now, I'm hardcoding these two
     type(minimum_dissipation_objective_function_t) :: objective_function
     type(volume_constraint_t) :: volume_constraint

     !> a steady simulation component to append to the forward
     type(steady_simcomp_t) :: steady_comp

   contains
     !> The common constructor using a JSON object.
     ! TODO
     ! we need to sort out the case file
     procedure, pass(this) :: init_base => steady_state_problem_init_base
     procedure, pass(this) :: init_design => &
          steady_state_problem_init_design_topopt
     ! but we could point to more depending on what design comes in
     !> Destructor.
     procedure, pass(this) :: free => steady_state_problem_free
     !> Computes the value of the objective and all constraints.
     !> ie, a forward simulation
     procedure, pass(this) :: compute_topopt => steady_state_problem_compute
     !> Computes the first order gradient of the objective function and
     ! all the constraints, and stores them in the design.
     procedure, pass(this) :: compute_sensitivity => &
          steady_state_problem_compute_sensitivity_topopt
     ! but we could point to more depending on what design is coming in
  end type steady_state_problem_t

contains
  !> The constructor for the base problem.
  subroutine steady_state_problem_init_base(this)
    class(steady_state_problem_t), intent(inout) :: this
    type(json_file) :: simcomp_settings

    ! append a steady state simcomp
    this%C%usr%init_user_simcomp => steady_state_simcomp
    ! call user_setup(this%C%usr)


    ! initialize the primal
    call neko_init(this%C)
    ! initialize the adjoint
    call adjoint_init(this%adj, this%C)

    ! TODO
    ! here we would read through our JSON to find out all of our constraints
    ! and objectives. NOTE, perhaps we'll just populate the list but not
    ! initialize them yet! As they may depend on the design.


  end subroutine steady_state_problem_init_base

  !> The constructor if a `topopt_design_t` is passed
  ! again, this is the only type of design we have so far...
  ! but in the future we may add other types of `design_variable_t`
  subroutine steady_state_problem_init_design_topopt(this, design)
    class(steady_state_problem_t), intent(inout) :: this
    type(topopt_design_t), target, intent(inout) :: design

    type(simple_brinkman_source_term_t) :: forward_brinkman, adjoint_brinkman

    ! Point the design to your own design
    this%design => design

    ! init the simple brinkman term for the forward problem
    call forward_brinkman%init_from_components( &
         this%C%fluid%f_x, this%C%fluid%f_y, this%C%fluid%f_z, &
         design, &
         this%C%fluid%u, this%C%fluid%v, this%C%fluid%w, &
         this%C%fluid%c_Xh)
    ! append brinkman source term to the forward problem
    call this%C%fluid%source_term%add(forward_brinkman)

    ! init the simple brinkman term for the adjoint
    call adjoint_brinkman%init_from_components( &
         this%adj%scheme%f_adj_x, this%adj%scheme%f_adj_y, &
         this%adj%scheme%f_adj_z, &
         design, &
         this%adj%scheme%u_adj, this%adj%scheme%v_adj, this%adj%scheme%w_adj, &
         this%adj%scheme%c_Xh)
    ! append brinkman source term based on design
    call this%adj%scheme%source_term%add(adjoint_brinkman)

    ! TODO
    ! Note, Tim, while you're reading this I'm sure you can already see we need
    ! to unmangle a lot of this.
    ! for instance,
    ! FIRST we have to read what our objective is and all of our constraints
    ! now if an objective involve the fluid (which it will) THIS will tell us
    ! we need to init a fluid and an adjoint
    ! SECOND now we know how many coefficients we need to map in the design
    ! THIRD we can start adding adjoint forcing etc...
    ! so the order is a little fucked up here...
    ! technically


    ! init the objective function
!---------------------------------------------------------
    ! - somehow append a user_check
    ! TODO:
    ! Tim, I loved what you did with with the source term handler. I'm hoping
    ! when you get a chance you can do something
    ! similar with simulation components?
    ! as in, this kind of post processing isn't just one function,
    ! but a list of post processing modules that can be appended
    ! to a simulation (and then we could append others to our adjoint!)
    !
    ! The thing is, because right now we're doing steady calculations,
    ! so the computations of:
    ! - The objective function: performed at the end of the steady run
    ! - The sensitivity:        performed at the end of the adjoint run
    !
    ! but when we move to unsteady calculations we'll have:
    ! - The objective function: accumulated DURING the forward run
    ! - The sensitivity:        accumulated DURING the adjoint run
    !
    ! So they'll have to be simcomps that get appended to C and adj.
    ! I trust you can whip that up lickity split!
    !
    ! in the future, the "problem" init will have already read out all the
    ! types of objective functions & constraints
    ! so in this place, we would already know from the JSON what
    ! objectives/constraints we need to init
    !
    ! So we need something akin to the `source_term_t`
    ! where we can have a list of them.
    !
    ! for this test we'll have 2
    ! minimum dissipation objective function
    call this%objective_function%init(design, this%C%fluid, this%adj%scheme)
    ! volume constraint
    call this%volume_constraint%init(design, this%C%fluid, this%adj%scheme)

    ! init the sampler
    !---------------------------------------------------------
    ! TODO
    ! obviously when we do the mappings properly, to many coeficients, we'll
    ! also have to modify this
    ! for now:
    ! - forward (p,u,v,w)                      1,2,3,4           p,vx,vy,vz
    ! - design (\rho)                          5                 temperature
    ! - adjoint (u,v,w,p)                      6,7,8,9           s1,s2,s3,s4
    ! - mapped (\chi)                          10                s5
    ! - sensitivity (dF/d\chi and dC/d\chi)    11, 12            s6,s7
    ! - sensitivity (dF/d\rho and dC/d\rho)    13, 14            s8,s9

    call this%output%init(sp, 'optimization', 13)
    call this%output%fields%assign(1, this%C%fluid%p)
    call this%output%fields%assign(2, this%C%fluid%u)
    call this%output%fields%assign(3, this%C%fluid%v)
    call this%output%fields%assign(4, this%C%fluid%w)
    ! I don't know why these ones need assign_to_field?
    call this%output%fields%assign_to_field(5, design%design_indicator)
    call this%output%fields%assign(6, this%adj%scheme%u_adj)
    call this%output%fields%assign(7, this%adj%scheme%v_adj)
    call this%output%fields%assign(8, this%adj%scheme%w_adj)
    call this%output%fields%assign(9, this%adj%scheme%p_adj)
    call this%output%fields%assign_to_field(10, design%brinkman_amplitude)
    call this%output%fields%assign_to_field(11, &
         this%objective_function%sensitivity_to_coefficient)
    call this%output%fields%assign_to_field(12, &
         this%volume_constraint%sensitivity_to_coefficient)
    call this%output%fields%assign_to_field(13, design%sensitivity)
    ! TODO
    ! I still haven't done the design%sensitivity as a field list!
    ! so it will eventually be
    ! call this%output%fields%assign(13, design%sensitivity(1))
    ! call this%output%fields%assign(14, design%sensitivity(2))
    ! or something to this effect

!------------------------------------------------------------------------------
! TODO
! the proceedure `steady_state_problem_compute_sensitivity_topopt` is currently
! the only one we have...
! but if we have a more abstract `design_variable_t` then we will need to
! including something here in the init that assigns the correct way of computing
! sensitivity given the `design_variable_t`

  end subroutine steady_state_problem_init_design_topopt

  !> Destructor.
  subroutine steady_state_problem_free(this)
    class(steady_state_problem_t), intent(inout) :: this
    call adjoint_free(this%adj)
    call neko_finalize(this%C)
    ! TODO
    ! probably also objective functions etc

  end subroutine steady_state_problem_free

  !> Here we compute all the objectives and constraints
  subroutine steady_state_problem_compute(this, design)
    class(steady_state_problem_t), intent(inout) :: this
    type(topopt_design_t), intent(in) :: design

    call neko_solve(this%C)
    ! TODO
    ! In the future, the objective_function_t will potentially include
    ! simulation components so that we can
    ! accumulate the objective function during the run...
    ! here, we just compute it on the last step
    ! TODO
    ! We would presumable have a list that holds all of objective functions
    ! and constraints, such that this would be a
    ! objectives%compute()
    call this%objective_function%compute(this%design, this%C%fluid)
    call this%volume_constraint%compute(this%design, this%C%fluid)
    print *, 'OBJECTIVE FUNCTION', &
         this%objective_function%objective_function_value
    print *, 'VOLUME CONSTRAINT', &
         this%volume_constraint%objective_function_value, &
         this%volume_constraint%volume


    ! TODO
    ! the steady simcomp only works for the forward, we either hardcode
    ! another one, or we have a forward adjoint flag
    ! or probably the smartest thing to do would be to accept a list of
    ! fields in the registry that will be checked...
    ! anyhow, I don't like the termination condition based on simulation
    ! time anyway...
    !



  end subroutine steady_state_problem_compute

  !> The computation of the sensitivity if we have a `topopt_design_t`.
  subroutine steady_state_problem_compute_sensitivity_topopt(this)
    class(steady_state_problem_t), intent(inout) :: this
    call solve_adjoint(this%adj)

    ! again, in the future, the objective_function_t will potentially include
    ! simulation components so that we can
    ! accumulate the sensitivity during the run...
    ! here, we just compute it on the last step

    ! TODO
    ! now that I look at this, we could have easily had design, fluid and
    ! adjoint passed in on init, and then have
    ! pointers pointing to what we're reffering to.
    ! So then compute would take in t and t-step.
    ! this would probably be smarter than passing the whole thing down!!!
    ! TODO
    ! We would presumable have a list that holds all of objective functions
    ! and constraints, such that this would be a
    ! objectives%compute_sensitivity()
    ! and it would cycled through the list.
    call this%objective_function%compute_sensitivity(&
         this%design, this%C%fluid, this%adj%scheme)
    call this%volume_constraint%compute_sensitivity(&
         this%design, this%C%fluid, this%adj%scheme)
    ! it would be nice to visualize this

    ! do the adjoint mapping
    call this%design%map_backward(&
         this%objective_function%sensitivity_to_coefficient)
    ! ok now you've fucked up the whole "list of sensitivity fields" aspect...
    ! we somehow need to populate the list


    ! TODO
    ! you've done this very incorrectly for the future,
    ! to do this properly you need:
    ! - a field list in the design holding sensitivity
    ! - a list of mappings from rho -> whatever enters the PDE
    ! - a list of adjoint mappings
    ! (the only reason this works is because we're considering the volume
    ! constraint on the unmapped field)
    ! - a subroutine in design%update that compiles all this stuff into a
    ! way readable for MMA

  end subroutine steady_state_problem_compute_sensitivity_topopt


  subroutine steady_state_simcomp(params)
    type(json_file), intent(inout) :: params
    type(steady_simcomp_t), allocatable :: steady_comp
    type(json_file) :: simcomp_settings

    ! Allocate a simulation component
    allocate(steady_comp)
    simcomp_settings = simulation_component_user_settings("steady", params)

    call neko_simcomps%add_user_simcomp(steady_comp, simcomp_settings)

  end subroutine steady_state_simcomp
end module steady_state_problem
