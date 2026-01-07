!> @file optimizer.f90
!! @copyright
!! Copyright (c) 2025-2026, The Neko-TOP Authors
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

!> Defines the abstract type `optimizer`
!! The optimizer type is defined to provide a generic interface to underlying
!! optimization methods. Specific optimizers should extend this type and
!! implement the deferred methods.
module optimizer

  use json_module, only: json_file
  use simulation_m, only: simulation_t
  use problem, only: problem_t
  use design, only: design_t
  use num_types, only: rp
  use logger, only: neko_log
  use profiler, only: profiler_start_region, profiler_end_region
  use mpi_f08, only: MPI_Wtime
  use utils, only: neko_error

  implicit none
  private

  !> Abstract optimizer class.
  type, abstract, public :: optimizer_t
     private

     !> The maximum number of iterations
     integer, public :: max_iterations = 0
     !> The tolerance for the optimization loop
     real(kind=rp), public :: tolerance = 0.0_rp
     !> Maximum runtime in seconds
     real(kind=rp), public :: max_runtime = -1.0_rp

   contains

     !  ---------------------------------------------------------------------- !
     ! Deferred procedures for specific optimizers

     !> Initialize the optimizer, associate it with a specific problem
     procedure(optimizer_init_from_json), pass(this), public, deferred :: &
          init_from_json
     !> Free resources.
     procedure(optimizer_free), pass(this), public, deferred :: free

     !> Prepare the optimizer before starting the optimization loop
     procedure(optimizer_initialize), pass(this), public, deferred :: initialize
     !> Perform a single optimization step
     procedure(optimizer_step), pass(this), public, deferred :: step
     !> Validate the solution
     procedure(optimizer_validate), pass(this), public, deferred :: validate

     !> Write the progress of the optimizer to the log file
     procedure(optimizer_write), pass(this), public, deferred :: write

     ! ----------------------------------------------------------------------- !
     ! Public procedures

     !> Run the optimization loop
     procedure, pass(this), public :: run => optimizer_run

     ! ----------------------------------------------------------------------- !
     ! Private procedures

     !> The base initializer
     procedure, pass(this) :: init_base => optimizer_init_base
     !> Free base resources.
     procedure, pass(this) :: free_base => optimizer_free_base
     !> Print status message
     procedure, pass(this) :: print_status => optimizer_print_status

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

  public :: optimizer_factory

contains

  ! -------------------------------------------------------------------------- !
  ! Base initializer and free routines

  !> Base initializer for the optimizer
  !! @param this The optimizer object.
  !! @param max_iterations The maximum number of iterations.
  !! @param tolerance The tolerance for the optimization loop.
  !! @param max_runtime The maximum runtime in seconds.
  subroutine optimizer_init_base(this, max_iterations, tolerance, max_runtime)
    class(optimizer_t), intent(inout) :: this
    integer, intent(in) :: max_iterations
    real(kind=rp), intent(in) :: tolerance
    real(kind=rp), intent(in), optional :: max_runtime

    this%max_iterations = max_iterations
    this%max_runtime = max_runtime
    this%tolerance = tolerance

  end subroutine optimizer_init_base

  !> Base free routine for the optimizer
  !! @param this The optimizer object.
  subroutine optimizer_free_base(this)
    class(optimizer_t), intent(inout) :: this

    this%max_iterations = 0
    this%tolerance = 0.0_rp
    this%max_runtime = -1.0_rp

  end subroutine optimizer_free_base

  ! -------------------------------------------------------------------------- !
  ! Optimization loop routine

  !> Define the optimization loop
  !! This subroutine runs the optimization loop until convergence
  !! or the maximum number of iterations is reached.
  !!
  !! The optimization loop can be terminated based on a maximum runtime. In this
  !! case, a cumulative average is used to determine if the next iteration would
  !! exceed the maximum runtime.
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
    real(kind=rp) :: start_time, elapsed_time
    real(kind=rp) :: iteration_start_time, iteration_end_time
    real(kind=rp) :: iteration_average_time
    logical :: converged
    integer :: iter, stop_flag

    ! Prepare the problem state before starting the optimization
    call this%initialize(problem, design, simulation)
    call this%write(0, problem)
    call design%write(0)

    call neko_log%section('Optimization Loop')

    stop_flag = 1
    converged = .false.
    start_time = MPI_Wtime()
    iteration_average_time = 0.0_rp
    do iter = 1, this%max_iterations

       call profiler_start_region('Optimizer iteration')
       iteration_start_time = MPI_Wtime()

       converged = this%step(iter, problem, design, simulation)

       iteration_end_time = MPI_Wtime()
       call profiler_end_region('Optimizer iteration')

       ! Log the progress and outputs
       call this%write(iter, problem)
       call design%write(iter)

       ! --------------------------------------------------------------------- !
       ! Check stopping criteria

       if (converged) then
          stop_flag = 0
          exit
       else if (this%max_runtime .gt. 0.0_rp) then
          elapsed_time = MPI_Wtime() - start_time

          ! Estimate Cumulative Average iteration time
          iteration_average_time = iteration_average_time * &
               (real(iter, kind=rp) / real(iter + 1, kind=rp)) + &
               (iteration_end_time - iteration_start_time) / &
               real(iter + 1, kind=rp)

          if (elapsed_time + iteration_average_time .gt. this%max_runtime) then
             stop_flag = 2
             exit
          end if
       end if
    end do

    ! Check that the final design is valid
    call this%validate(problem, design)
    call this%print_status(stop_flag, iter)

    call neko_log%end_section()

  end subroutine optimizer_run

  !> Print status message
  !! Supported flags:
  !! 0: Converged successfully             (SUCCESS)
  !! 1: Did not converge in max iterations (WARNING)
  !! 2: Stopped after reaching max runtime (ERROR)
  !!
  !! @param this The optimizer object.
  !! @param stop_flag The stopping flag.
  !! @param iter The number of iterations performed.
  subroutine optimizer_print_status(this, stop_flag, iter)
    class(optimizer_t), intent(in) :: this
    integer, intent(in) :: stop_flag
    integer, intent(in) :: iter
    character(len=256) :: msg

    select case (stop_flag)
    case (0)
       write(msg, '(A,I0,A)') 'Optimizer converged successfully after ', &
            iter, ' iterations.'
       call neko_log%message(msg)
    case (1)
       write(msg, '(A,I0,A)') 'Optimizer did not converge in ', &
            this%max_iterations, ' iterations.'
       call neko_log%warning(msg)
    case (2)
       write(msg, '(A,A,F8.2,A)') 'Optimizer stopped after reaching the ', &
            'maximum runtime of ', this%max_runtime, ' seconds.'
       call neko_error(trim(msg))

    case default
       write(msg, '(A)') 'Optimizer stopped for an unknown reason.'
       call neko_error(msg)
    end select
  end subroutine optimizer_print_status

end module optimizer
