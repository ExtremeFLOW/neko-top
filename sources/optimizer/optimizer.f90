!> @file optimizer.f90
!! @copyright
!! Copyright (c) 2025, The Neko-TOP Authors
!! All rights reserved.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions
!! are met:
!!
!!   * Redistributions of source code must retain the above copyright
!!     notice, this list of conditions and the following disclaimer.
!!
!!   * Redistributions in binary form must reproduce the above
!!     copyright notice, this list of conditions and the following
!!     disclaimer in the documentation and/or other materials provided
!!     with the distribution.
!!
!!   * Neither the name of the authors nor the names of its
!!     contributors may be used to endorse or promote products derived
!!     from this software without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
!! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
!! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
!! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
!! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
!! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
!! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!! POSSIBILITY OF SUCH DAMAGE.

module optimizer
  !-----------------------------------------------------------!
  ! An abstract type optimizer is defined to solve the design !
  ! optimization problem using a specific type of optimizer   !
  ! algorithm, e.g. MMA.                                      !
  !-----------------------------------------------------------!

  use json_module, only: json_file
  use simulation_m, only: simulation_t
  use problem, only: problem_t
  use design, only: design_t
  use num_types, only: rp
  use logger, only: neko_log
  use profiler, only: profiler_start_region, profiler_end_region

  implicit none
  private

  !> Abstract optimizer class.
  type, abstract :: optimizer_t
     private

     !> The maximum number of iterations
     integer, public :: max_iterations
     !> The tolerance for the optimization loop
     real(kind=rp), public :: tolerance

   contains
     !> Initialize the optimizer, associate it with a specific problem
     procedure(optimizer_init_from_json), pass(this), public, deferred :: &
          init_from_json
     !> Run the optimization loop
     procedure, pass(this), public :: run => optimizer_run
     !> Prepare the optimizer before starting the optimization loop
     procedure(optimizer_initialize), pass(this), public, deferred :: &
          initialize
     !> Perform a single optimization step
     procedure(optimizer_step), pass(this), public, deferred :: step
     !> Free resources.
     procedure(optimizer_free), pass(this), public, deferred :: free

     !> Validate the solution
     procedure(optimizer_validate), pass(this), public, deferred :: validate
     !> Write the progress of the optimizer to the log file
     procedure(optimizer_write), pass(this), public, deferred :: write

     !> The base initializer
     procedure, pass(this) :: init_base => optimizer_init_base
     !> Free base resources.
     procedure, pass(this) :: free_base => optimizer_free_base

  end type optimizer_t

  ! -------------------------------------------------------------------------- !
  ! Interface for the optimizer module.

  abstract interface
     !> Interface for optimizer initialization
     subroutine optimizer_init_from_json(this, parameters, problem, design, &
          simulation)
       import optimizer_t, json_file, simulation_t, problem_t, design_t, rp
       class(optimizer_t), intent(inout) :: this
       type(json_file), intent(inout) :: parameters
       class(problem_t), intent(inout) :: problem
       class(design_t), intent(in) :: design
       type(simulation_t), optional, intent(in) :: simulation
     end subroutine optimizer_init_from_json

     !> Interface for running an optimization initialization
     !! This subroutine initializes the optimizer before starting the
     !! optimization loop.
     subroutine optimizer_initialize(this, problem, design, simulation)
       import optimizer_t, simulation_t, problem_t, design_t
       class(optimizer_t), intent(inout) :: this
       class(problem_t), intent(inout) :: problem
       class(design_t), intent(inout) :: design
       type(simulation_t), optional, intent(inout) :: simulation
     end subroutine optimizer_initialize

     !> Interface for running an optimization step
     logical function optimizer_step(this, iter, problem, design, simulation)
       import optimizer_t, simulation_t, problem_t, design_t
       integer, intent(in) :: iter
       class(optimizer_t), intent(inout) :: this
       class(problem_t), intent(inout) :: problem
       class(design_t), intent(inout) :: design
       type(simulation_t), optional, intent(inout) :: simulation
     end function optimizer_step

     !> Interface for writing the optimizer progress
     subroutine optimizer_write(this, iter, problem)
       import optimizer_t, simulation_t, problem_t, design_t
       class(optimizer_t), intent(inout) :: this
       integer, intent(in) :: iter
       class(problem_t), intent(in) :: problem
     end subroutine optimizer_write

     !> Interface for freeing resources
     subroutine optimizer_free(this)
       import optimizer_t
       class(optimizer_t), intent(inout) :: this
     end subroutine optimizer_free

     !> Interface for validating the solution
     subroutine optimizer_validate(this, problem, design)
       import optimizer_t, problem_t, design_t
       class(optimizer_t), intent(inout) :: this
       class(problem_t), intent(in) :: problem
       class(design_t), intent(in) :: design
     end subroutine optimizer_validate
  end interface

  ! -------------------------------------------------------------------------- !
  ! Interfaces for the factory functions

  !> Factory function for the optimizer
  !! @param object The optimizer object to be created.
  !! @param parameters The JSON file containing the optimizer parameters.
  !! @param problem The problem object.
  !! @param design The design object.
  !! @param simulation The simulation object.
  interface optimizer_factory
     module subroutine optimizer_factory(object, parameters, problem, design, &
          simulation)
       class(optimizer_t), allocatable, intent(inout) :: object
       type(json_file), intent(inout) :: parameters
       class(problem_t), intent(inout) :: problem
       class(design_t), intent(in) :: design
       type(simulation_t), optional, intent(in) :: simulation
     end subroutine optimizer_factory
  end interface optimizer_factory

  public :: optimizer_t, optimizer_factory

contains

  ! -------------------------------------------------------------------------- !
  ! Base initializer and free routines

  !> Base initializer for the optimizer
  !! @param this The optimizer object.
  !! @param max_iterations The maximum number of iterations.
  !! @param tolerance The tolerance for the optimization loop.
  subroutine optimizer_init_base(this, max_iterations, tolerance)
    class(optimizer_t), intent(inout) :: this
    integer, intent(in) :: max_iterations
    real(kind=rp), intent(in) :: tolerance

    this%max_iterations = max_iterations
    this%tolerance = tolerance

  end subroutine optimizer_init_base

  !> Base free routine for the optimizer
  !! @param this The optimizer object.
  subroutine optimizer_free_base(this)
    class(optimizer_t), intent(inout) :: this
  end subroutine optimizer_free_base

  ! -------------------------------------------------------------------------- !
  ! Optimization loop routine

  !> Define the optimization loop
  !! This subroutine runs the optimization loop until convergence
  !! or the maximum number of iterations is reached.
  !!
  !! @param this The optimizer object.
  !! @param problem The problem object.
  !! @param design The design object.
  !! @param simulation The simulation object.
  subroutine optimizer_run(this, problem, design, simulation)
    class(optimizer_t), intent(inout) :: this
    class(problem_t), intent(inout) :: problem
    class(design_t), intent(inout) :: design
    type(simulation_t), optional, intent(inout) :: simulation
    logical :: converged
    character(len=256) :: msg
    integer :: iter

    ! Prepare the problem state before starting the optimization
    call this%initialize(problem, design, simulation)
    call this%write(0, problem)
    call design%write(0)

    call neko_log%section('Optimization Loop')

    iter = 1
    converged = .false.
    do while (iter .le. this%max_iterations .and. .not. converged)

       call profiler_start_region('Optimizer iteration')
       converged = this%step(iter, problem, design, simulation)
       call profiler_end_region('Optimizer iteration')

       ! Log the progress and outputs
       call this%write(iter, problem)
       call design%write(iter)

       iter = iter + 1
    end do

    call neko_log%end_section()

    ! Check that the final design is valid
    call this%validate(problem, design)

    if (.not. converged) then
       write(msg, '(A,I0,A)') 'Optimizer did not converge in ', &
            this%max_iterations, ' iterations.'
       call neko_log%warning(msg)
    else
       write(msg, '(A,I0,A)') 'Optimizer converged after ', iter, &
            ' iterations.'
       call neko_log%message(msg)
    end if

  end subroutine optimizer_run


end module optimizer
