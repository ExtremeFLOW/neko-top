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
  use design, only: design_t
  use topopt_design, only: topopt_design_t
  use json_file_module, only: json_file
  use objective, only: objective_t, objective_factory
  use constraint, only: constraint_t, constraint_factory
  use fld_file_output, only: fld_file_output_t
  use steady_simcomp, only: steady_simcomp_t
  use simple_brinkman_source_term, only: simple_brinkman_source_term_t
  use json_utils_ext, only: json_key_fallback, json_get_subdict
  use json_value_module, only: json_value
  use utils, only: neko_error
  use user_intf, only: simulation_component_user_settings
  use simcomp_executor, only: neko_simcomps
  use neko_ext, only: reset
  use field_math, only: field_rzero
  use simulation, only: simulation_t
  use field, only: field_t
  use scratch_registry, only: neko_scratch_registry
  use math, only: copy
  use vector, only: vector_t
  use matrix, only: matrix_t
  implicit none
  private

  !> To compute a steady state problem
  type, public, extends(problem_t) :: steady_state_problem_t

     !> a steady simulation component to append to the forward
     type(steady_simcomp_t) :: steady_comp

   contains
     !> The common constructor using a JSON object.
     procedure, pass(this) :: init => steady_state_problem_init

     !> Constructor of a generic design problem.
     procedure, pass(this) :: init_design => steady_state_problem_init_design
     procedure, pass(this) :: init_design_topopt => &
          steady_state_problem_init_design_topopt
     ! but we could point to more depending on what design comes in
     !> Destructor.
     procedure, pass(this) :: free => steady_state_problem_free

     !> Generic compute function.
     procedure, pass(this) :: compute => steady_state_problem_compute
     !> Generic compute function for sensitivity.
     procedure, pass(this) :: compute_sensitivity => &
          steady_state_problem_compute_sensitivity

  end type steady_state_problem_t

contains

  !> The constructor for the base problem.
  subroutine steady_state_problem_init(this)
    class(steady_state_problem_t), intent(inout) :: this

    call this%simulation%init()
    call this%init_base(this%simulation%neko_case%fluid%dm_Xh%size())

    ! TODO
    ! here we would read through our JSON to find out all of our constraints
    ! and objectives. NOTE, perhaps we'll just populate the list but not
    ! initialize them yet! As they may depend on the design.


  end subroutine steady_state_problem_init

  !> The constructor of a generic design problem.
  !! @params this: The problem to initialize.
  !! @params design: The design to initialize the problem with.
  subroutine steady_state_problem_init_design(this, design)
    class(steady_state_problem_t), intent(inout) :: this
    class(design_t), target, intent(inout) :: design

    select type(design)
      type is(topopt_design_t)
       call this%init_design_topopt(design)
      class default
       call neko_error('Only topopt_design_t is supported for now')
    end select

  end subroutine steady_state_problem_init_design

  !> The constructor if a `topopt_design_t` is passed
  ! again, this is the only type of design we have so far...
  ! but in the future we may add other types of `design_variable_t`
  subroutine steady_state_problem_init_design_topopt(this, design)
    class(steady_state_problem_t), intent(inout) :: this
    type(topopt_design_t), target, intent(inout) :: design

    type(simple_brinkman_source_term_t) :: forward_brinkman, adjoint_brinkman

    !> TODO
    ! we need a `objective_list` which is allocatable and contains a factory
    ! to fill itself up with from the JSON
    ! for now, I'm hardcoding these two
    ! class(objective_t), allocatable :: objective_function
    ! class(constraint_t), allocatable :: volume_constraint

    ! init the design
    call design%init(this%simulation%neko_case%params, this%simulation%neko_case%fluid%c_Xh)

    ! init the simple brinkman term for the forward problem
    call forward_brinkman%init_from_components( &
         this%simulation%neko_case%fluid%f_x, &
         this%simulation%neko_case%fluid%f_y, &
         this%simulation%neko_case%fluid%f_z, &
         design, &
         this%simulation%neko_case%fluid%u, &
         this%simulation%neko_case%fluid%v, &
         this%simulation%neko_case%fluid%w, &
         this%simulation%neko_case%fluid%c_Xh)
    ! append brinkman source term to the forward problem
    call this%simulation%neko_case%fluid%source_term%add(forward_brinkman)

    ! init the simple brinkman term for the adjoint
    call adjoint_brinkman%init_from_components( &
         this%simulation%adjoint_case%scheme%f_adj_x, &
         this%simulation%adjoint_case%scheme%f_adj_y, &
         this%simulation%adjoint_case%scheme%f_adj_z, &
         design, &
         this%simulation%adjoint_case%scheme%u_adj, &
         this%simulation%adjoint_case%scheme%v_adj, &
         this%simulation%adjoint_case%scheme%w_adj, &
         this%simulation%adjoint_case%scheme%c_Xh)
    ! append brinkman source term based on design
    call this%simulation%adjoint_case%scheme%source_term%add(adjoint_brinkman)

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
    call this%read_objectives(this%simulation%neko_case%params, design)

    ! volume constraint
    call this%read_constraints(this%simulation%neko_case%params, design)

    ! init the sampler
    !---------------------------------------------------------
    ! TODO
    ! obviously when we do the mappings properly, to many coefficients, we'll
    ! also have to modify this
    ! for now:
    ! - forward (p,u,v,w)                      1,2,3,4           p,vx,vy,vz
    ! - design (\rho)                          5                 temperature
    ! - adjoint (u,v,w,p)                      6,7,8,9           s1,s2,s3,s4
    ! - mapped (\chi)                          10                s5
    ! - sensitivity (dF/d\chi and dC/d\chi)    11, 12            s6,s7
    ! - sensitivity (dF/d\rho and dC/d\rho)    13, 14            s8,s9

    ! Allocate the output type
    call this%output%init(sp, 'optimization', 10)
    call this%output%fields%assign(1, this%simulation%neko_case%fluid%p)
    call this%output%fields%assign(2, this%simulation%neko_case%fluid%u)
    call this%output%fields%assign(3, this%simulation%neko_case%fluid%v)
    call this%output%fields%assign(4, this%simulation%neko_case%fluid%w)
    ! I don't know why these ones need assign_to_field?
    call this%output%fields%assign(5, design%design_indicator)
    call this%output%fields%assign(6, this%simulation%adjoint_case%scheme%u_adj)
    call this%output%fields%assign(7, this%simulation%adjoint_case%scheme%v_adj)
    call this%output%fields%assign(8, this%simulation%adjoint_case%scheme%w_adj)
    call this%output%fields%assign(9, this%simulation%adjoint_case%scheme%p_adj)
    call this%output%fields%assign(10, design%brinkman_amplitude)

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

    call this%free_base()
    call this%simulation%free()

  end subroutine steady_state_problem_free

  !> The computation of the objective function and constraints.
  subroutine steady_state_problem_compute(this, design)
    class(steady_state_problem_t), intent(inout) :: this
    class(design_t), intent(inout) :: design

    call this%simulation%run_forward()

    call this%update_objectives(design)
    call this%update_constraints(design)

  end subroutine steady_state_problem_compute

  !> The computation of the objective function and constraints.
  subroutine steady_state_problem_compute_sensitivity(this, design)
    class(steady_state_problem_t), intent(inout) :: this
    class(design_t), intent(inout) :: design

    select type (design)
      type is (topopt_design_t)
       call steady_state_problem_compute_sensitivity_topopt(this, design)
      class default
       call neko_error('Only topopt_design_t is supported for now')
    end select

  end subroutine steady_state_problem_compute_sensitivity

  !> The computation of the sensitivity if we have a `topopt_design_t`.
  subroutine steady_state_problem_compute_sensitivity_topopt(this, design)
    class(steady_state_problem_t), intent(inout) :: this
    type(topopt_design_t), intent(inout) :: design

    type(vector_t) :: objective_sensitivity

    type(field_t), pointer :: objective_sensitivity_field
    integer, dimension(1) :: temp_indices

    call this%simulation%run_backward()

    call this%update_objective_sensitivities(design)
    call this%update_constraint_sensitivities(design)

    call this%get_objective_sensitivities(objective_sensitivity)

    ! it would be nice to visualize this

    call neko_scratch_registry%request_field(objective_sensitivity_field, &
         temp_indices(1))
    call copy(objective_sensitivity_field%x, objective_sensitivity%x, &
         this%get_n_design())

    ! do the adjoint mapping
    call design%map_backward(objective_sensitivity_field)
    ! ok now you've fucked up the whole "list of sensitivity fields" aspect...
    ! we somehow need to populate the list

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine steady_state_problem_compute_sensitivity_topopt

end module steady_state_problem
