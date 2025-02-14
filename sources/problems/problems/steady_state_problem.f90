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
  subroutine steady_state_problem_init(this, parameters, design, simulation)
    class(steady_state_problem_t), intent(inout) :: this
    type(json_file), intent(inout) :: parameters
    class(design_t), intent(in) :: design
    type(simulation_t), intent(inout) :: simulation

    call this%init_base(design%size())

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
    call this%output%fields%assign(1, simulation%fluid_scheme%p)
    call this%output%fields%assign(2, simulation%fluid_scheme%u)
    call this%output%fields%assign(3, simulation%fluid_scheme%v)
    call this%output%fields%assign(4, simulation%fluid_scheme%w)
    ! I don't know why these ones need assign_to_field?
    call this%output%fields%assign(6, simulation%adjoint_case%scheme%u_adj)
    call this%output%fields%assign(7, simulation%adjoint_case%scheme%v_adj)
    call this%output%fields%assign(8, simulation%adjoint_case%scheme%w_adj)
    call this%output%fields%assign(9, simulation%adjoint_case%scheme%p_adj)

    ! TODO
    ! here we would read through our JSON to find out all of our constraints
    ! and objectives. NOTE, perhaps we'll just populate the list but not
    ! initialize them yet! As they may depend on the design.

    select type(d => design)
    type is(topopt_design_t)

       call this%output%fields%assign_to_field(5, d%design_indicator)
       call this%output%fields%assign_to_field(10, d%brinkman_amplitude)

    class default
       call neko_error('Only topopt_design_t is supported for now')
    end select

    ! minimum dissipation objective function
    call this%read_objectives(parameters, simulation, design)

    ! volume constraint
    call this%read_constraints(parameters, simulation, design)

  end subroutine steady_state_problem_init

  !> Destructor.
  subroutine steady_state_problem_free(this)
    class(steady_state_problem_t), intent(inout) :: this

    call this%free_base()

  end subroutine steady_state_problem_free

  !> The computation of the objective function and constraints.
  subroutine steady_state_problem_compute(this, design)
    class(steady_state_problem_t), intent(inout) :: this
    class(design_t), intent(inout) :: design

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

    call this%update_objective_sensitivities(design)
    call this%update_constraint_sensitivities(design)

    call this%get_objective_sensitivities(objective_sensitivity)

    call design%map_backward(objective_sensitivity)

  end subroutine steady_state_problem_compute_sensitivity_topopt

end module steady_state_problem
