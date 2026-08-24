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

!> Defines the abstract type `optimizer`.
!! @details
!! The optimizer type is defined to provide a generic interface to underlying
!! optimization methods. Specific optimizers should extend this type and
!! implement the deferred methods.
module optimizer
  use json_module, only: json_file
  use simulation_m, only: simulation_t
  use problem, only: problem_t
  use design, only: design_t
  use num_types, only: rp
  use logger, only: neko_log, LOG_SIZE
  use profiler, only: profiler_start_region, profiler_end_region
  use mpi_f08, only: MPI_Wtime, MPI_Allreduce, MPI_MAX
  use utils, only: neko_error, filename_suffix, read_duration
  use csv_file, only: csv_file_t
  use vector, only: vector_t
  use json_utils, only: json_get_or_default
  use comm, only: pe_rank, MPI_REAL_PRECISION, NEKO_COMM
  use continuation_scheduler, only: nekotop_continuation
  use math, only: Neko_EPS

  implicit none
  private

  !> Abstract optimizer class.
  type, abstract, public :: optimizer_t

     !> The type of the optimizer
     character(len=64), private :: optimizer_type = ''
     !> The maximum number of iterations
     integer, private :: max_iterations = 0
     !> The current iteration number
     integer, private :: current_iteration = 0
     !> The smallest change in design variables.
     real(kind=rp), private :: stop_design_change = NEKO_EPS
     !> The maximum observed design change.
     real(kind=rp), public :: max_design_change = 0.0_rp

     ! ----------------------------------------------------------------------- !
     ! Restart related members

     !> Checkpoint file to be restarted from.
     character(len=256), private :: checkpoint_file = ''

     ! Checkpoint related information
     character(len=256), private :: checkpoint_path = './checkpoints/'
     character(len=256), private :: checkpoint_base = 'optimizer_checkpoint'
     character(len=256), private :: checkpoint_format = 'hdf5'
     integer, private :: checkpoint_interval = -1

     ! Variables for the runtime-based stopping criteria
     real(kind=rp), private :: max_runtime = -1.0_rp
     real(kind=rp), private :: start_time = 0.0_rp
     real(kind=rp), private :: average_time = 0.0_rp
     real(kind=rp), private :: step_count = 0.0_rp

     ! Logging state
     logical, private :: log_initialized = .false.
     logical, private :: log_include_constraints = .true.
     integer, private :: log_extra_size = 0
     type(csv_file_t), private :: log_file
     type(vector_t), private :: log_data

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
     !> Save optimizer-specific components to checkpoint
     procedure(optimizer_save_checkpoint_components), pass(this), deferred :: &
          save_checkpoint_components
     !> Load optimizer-specific components from checkpoint
     procedure(optimizer_load_checkpoint_components), pass(this), deferred :: &
          load_checkpoint_components

     !> Save a checkpoint of the optimizer state
     procedure, pass(this) :: save_checkpoint => optimizer_save_checkpoint
     !> Restore the optimizer state from a checkpoint
     procedure, pass(this) :: load_checkpoint => optimizer_load_checkpoint

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
     !> Read settings from JSON parameters file
     procedure, pass(this) :: read_base_settings => optimizer_read_base_settings
     !> Print status message
     procedure, pass(this) :: print_status => optimizer_print_status
     !> Estimate if we are out of time
     procedure, pass(this) :: out_of_time => optimizer_out_of_time
     !> Initialize optimization log
     procedure, pass(this) :: init_log => optimizer_init_log
     !> Write optimization log entry
     procedure, pass(this) :: write_log => optimizer_write_log
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

     !> Interface for freeing resources
     subroutine optimizer_free(this)
       import optimizer_t
       class(optimizer_t), intent(inout) :: this
     end subroutine optimizer_free

     !> Interface for running an optimization step
     logical function optimizer_step(this, iter, problem, design, simulation)
       import optimizer_t, simulation_t, problem_t, design_t
       class(optimizer_t), intent(inout) :: this
       integer, intent(in) :: iter
       class(problem_t), intent(inout) :: problem
       class(design_t), intent(inout) :: design
       type(simulation_t), optional, intent(inout) :: simulation
     end function optimizer_step

     !> Interface for validating the solution
     subroutine optimizer_validate(this, problem, design)
       import optimizer_t, problem_t, design_t
       class(optimizer_t), intent(inout) :: this
       class(problem_t), intent(in) :: problem
       class(design_t), intent(in) :: design
     end subroutine optimizer_validate

     !> Interface for writing the optimizer progress
     subroutine optimizer_write(this, iter, problem)
       import optimizer_t, simulation_t, problem_t, design_t
       class(optimizer_t), intent(inout) :: this
       integer, intent(in) :: iter
       class(problem_t), intent(inout) :: problem
     end subroutine optimizer_write

     !> Interface for saving optimizer-specific components to checkpoint
     subroutine optimizer_save_checkpoint_components(this, filename, overwrite)
       import optimizer_t
       class(optimizer_t), intent(inout) :: this
       character(len=*), intent(in) :: filename
       logical, intent(in), optional :: overwrite
     end subroutine optimizer_save_checkpoint_components

     !> Interface for loading optimizer-specific components from checkpoint
     subroutine optimizer_load_checkpoint_components(this, filename)
       import optimizer_t
       class(optimizer_t), intent(inout) :: this
       character(len=*), intent(in) :: filename
     end subroutine optimizer_load_checkpoint_components

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

  ! -------------------------------------------------------------------------- !
  ! IO routines for HDF5 checkpoints

  interface
     !> Interface for writing a checkpoint
     module subroutine optimizer_save_checkpoint_hdf5(object, filename, iter, &
          overwrite)
       class(optimizer_t), intent(inout) :: object
       character(len=*), intent(in) :: filename
       integer, intent(in) :: iter
       logical, intent(in), optional :: overwrite
     end subroutine optimizer_save_checkpoint_hdf5

     !> Interface for reading a checkpoint
     module subroutine optimizer_load_checkpoint_hdf5(object, filename, iter)
       class(optimizer_t), intent(inout) :: object
       character(len=*), intent(in) :: filename
       integer, intent(out) :: iter
     end subroutine optimizer_load_checkpoint_hdf5
  end interface

  public :: optimizer_factory

contains

  ! -------------------------------------------------------------------------- !
  ! Base initializer and free routines

  !> Base initializer for the optimizer
  !! @param this The optimizer object.
  !! @param optimizer_type The type of the optimizer.
  !! @param max_iterations The maximum number of iterations.
  !! @param max_runtime The maximum runtime in seconds.
  !! @param stop_design_change The stopping criteria for design change.
  !! @param checkpoint_file The checkpoint file to restart from.
  !! @param checkpoint_path The path for saving checkpoint files.
  !! @param checkpoint_base The base name for checkpoint files.
  !! @param checkpoint_format The file format for checkpoint files.
  !! @param checkpoint_interval The interval for saving checkpoints in
  !!        iterations.
  subroutine optimizer_init_base(this, optimizer_type, max_iterations, &
       max_runtime, stop_design_change, checkpoint_file, checkpoint_path, &
       checkpoint_base, checkpoint_format, checkpoint_interval)
    class(optimizer_t), intent(inout) :: this
    character(len=*), intent(in) :: optimizer_type
    integer, intent(in) :: max_iterations
    real(kind=rp), intent(in), optional :: max_runtime
    real(kind=rp), intent(in), optional :: stop_design_change
    character(len=*), intent(in), optional :: checkpoint_file
    character(len=*), intent(in), optional :: checkpoint_path
    character(len=*), intent(in), optional :: checkpoint_base
    character(len=*), intent(in), optional :: checkpoint_format
    integer, intent(in), optional :: checkpoint_interval

    ! Mandatory settings
    this%optimizer_type = optimizer_type
    this%max_iterations = max_iterations

    ! Optional settings
    if (present(max_runtime)) this%max_runtime = max_runtime
    if (present(stop_design_change)) this%stop_design_change = stop_design_change
    if (present(checkpoint_file)) this%checkpoint_file = checkpoint_file
    if (present(checkpoint_path)) this%checkpoint_path = checkpoint_path
    if (present(checkpoint_base)) this%checkpoint_base = checkpoint_base
    if (present(checkpoint_format)) this%checkpoint_format = checkpoint_format
    if (present(checkpoint_interval)) then
       this%checkpoint_interval = checkpoint_interval
    end if

    ! Initialize internals
    this%start_time = MPI_Wtime()

  end subroutine optimizer_init_base

  !> Base free routine for the optimizer
  !! @param this The optimizer object.
  subroutine optimizer_free_base(this)
    class(optimizer_t), intent(inout) :: this

    this%optimizer_type = ''
    this%max_iterations = 0
    this%max_runtime = -1.0_rp
    this%stop_design_change = -1.0_rp
    this%max_design_change = 0.0_rp
    this%checkpoint_file = ''
    this%checkpoint_path = './checkpoints/'
    this%checkpoint_base = 'optimizer_checkpoint'
    this%checkpoint_format = 'hdf5'
    this%checkpoint_interval = -1

    this%start_time = 0.0_rp
    this%current_iteration = 0
    call this%log_data%free()
    this%log_initialized = .false.
    this%log_extra_size = 0
    this%log_include_constraints = .true.

  end subroutine optimizer_free_base

  !> Read settings from JSON parameters file.
  !! @param this The optimizer object.
  !! @param solver_params The JSON file containing the optimizer parameters.
  subroutine optimizer_read_base_settings(this, solver_params)
    class(optimizer_t), intent(inout) :: this
    type(json_file), intent(inout) :: solver_params
    integer :: read_int
    character(len=:), allocatable :: read_str

    call json_get_or_default(solver_params, 'max_runtime', read_str, "")
    call read_duration(read_str, this%max_runtime)
    call json_get_or_default(solver_params, 'stop_design_change', &
         this%stop_design_change, NEKO_EPS)

    call json_get_or_default(solver_params, 'restart_file', read_str, &
         this%checkpoint_file)
    this%checkpoint_file = read_str

    call json_get_or_default(solver_params, 'checkpoint.path', read_str, &
         this%checkpoint_path)
    this%checkpoint_path = read_str
    call json_get_or_default(solver_params, 'checkpoint.base', read_str, &
         this%checkpoint_base)
    this%checkpoint_base = read_str
    call json_get_or_default(solver_params, 'checkpoint.format', read_str, &
         this%checkpoint_format)
    this%checkpoint_format = read_str
    call json_get_or_default(solver_params, 'checkpoint.interval', read_int, &
         this%checkpoint_interval)
    this%checkpoint_interval = read_int

  end subroutine optimizer_read_base_settings

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
    real(kind=rp) :: iteration_time
    character(len=1024) :: checkpoint_file
    logical :: converged, file_exists
    integer :: stop_flag

    ! Initialize variables
    stop_flag = 1
    converged = .false.

    ! Restart from checkpoint if available
    if (trim(this%checkpoint_file) .ne. '') then
       checkpoint_file = trim(this%checkpoint_file)
    else
       checkpoint_file = optimizer_checkpoint_filename(this, &
            basename = 'optimizer_rt_checkpoint')
    end if

    inquire(file = checkpoint_file, exist = file_exists)
    if (file_exists) then
       call this%load_checkpoint(checkpoint_file, this%current_iteration, &
            design)
    end if

    ! compute potential internals for for a potential restart
    if (this%current_iteration .ne. 0) then
       call design%set_output_counter(this%current_iteration - 1)
       if (present(simulation)) then
          call simulation%set_output_counter(this%current_iteration - 1)
       end if
    end if

    ! Prepare the problem state before starting the optimization
    call this%initialize(problem, design, simulation)

    call this%write(this%current_iteration, problem)
    call design%write(this%current_iteration)

    call neko_log%section('Optimization Loop')

    do while (this%current_iteration .lt. this%max_iterations)
       this%current_iteration = this%current_iteration + 1
       if (pe_rank .eq. 0) then
          write(*,*) 'Starting iteration ', this%current_iteration
       end if

       call profiler_start_region('Optimizer iteration')
       iteration_time = MPI_Wtime()

       ! Update the parameters in continuation scheduler
       call nekotop_continuation%update(this%current_iteration)

       converged = this%step(this%current_iteration, problem, design, &
            simulation)

       iteration_time = MPI_Wtime() - iteration_time
       call profiler_end_region('Optimizer iteration')

       ! Log the progress and outputs
       call this%write(this%current_iteration, problem)
       call design%write(this%current_iteration)

       ! Save checkpoint if enabled
       if (this%checkpoint_interval .gt. 0 .and. &
            mod(this%current_iteration, this%checkpoint_interval) == 0) then
          call this%save_checkpoint(this%current_iteration, design, .false.)
       end if

       ! --------------------------------------------------------------------- !
       ! Check stopping criteria

       if (converged) then
          stop_flag = 0
          exit
       else if (this%current_iteration .ge. this%max_iterations) then
          stop_flag = 1
          exit
       else if (this%max_design_change .lt. this%stop_design_change) then
          stop_flag = 2
          exit
       else if (this%out_of_time(iteration_time)) then
          call this%save_checkpoint(this%current_iteration, design, .true., &
               basename = 'optimizer_rt_checkpoint')
          stop_flag = 3
          exit
       end if
    end do

    ! Check that the final design is valid
    call this%validate(problem, design)
    call this%print_status(stop_flag, this%current_iteration)

    call neko_log%end_section()

  end subroutine optimizer_run

  ! ========================================================================== !
  ! Helper routines

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
       write(msg, '(A,A,F8.2,A)') 'Optimizer stopped, design change ', &
            'below threshold of: ', this%stop_design_change, '.'
       call neko_log%warning(msg)
    case (3)
       write(msg, '(A)') 'Optimizer stopped due to runtime limit.'
       call neko_error(msg)

    case default
       write(msg, '(A)') 'Optimizer stopped for an unknown reason.'
       call neko_error(msg)
    end select
  end subroutine optimizer_print_status

  !> Estimate if we are out of time.
  !! This function uses a cumulative average of iteration times to
  !! estimate if the next iteration would exceed the maximum runtime.
  !! @param this The optimizer object.
  !! @param step_time The time taken for the latest iteration.
  !! @return out_of_time Logical indicating if we are out of time.
  function optimizer_out_of_time(this, step_time) result(out_of_time)
    class(optimizer_t), intent(inout) :: this
    real(kind=rp), intent(in) :: step_time
    logical :: out_of_time
    real(kind=rp) :: elapsed_time, time, old_avg_weight

    out_of_time = .false.

    if (this%max_runtime .lt. 0.0_rp) then
       return
    end if

    call MPI_Allreduce(step_time, time, 1, MPI_REAL_PRECISION, MPI_MAX, &
         NEKO_COMM)

    elapsed_time = MPI_Wtime() - this%start_time
    this%step_count = this%step_count + 1.0_rp
    old_avg_weight = (this%step_count - 1) / this%step_count

    ! Estimate Cumulative Average iteration time
    this%average_time = time / this%step_count + &
         this%average_time * old_avg_weight

    ! Determine if next iteration would exceed max runtime
    out_of_time = (elapsed_time + this%average_time) .gt. this%max_runtime

  end function optimizer_out_of_time

  ! ========================================================================== !
  ! Logging helpers

  !> Initialize optimization log.
  !! @param[inout] this The optimizer object.
  !! @param[in] problem The problem object.
  !! @param[in] extra_headers Header labels for extra log entries.
  !! @param[in] include_constraints Include constraints in the log.
  !! @param[in] filename Output filename for the log.
  subroutine optimizer_init_log(this, problem, extra_headers, &
       include_constraints, filename)
    class(optimizer_t), intent(inout) :: this
    class(problem_t), intent(in) :: problem
    character(len=*), intent(in), optional :: extra_headers(:)
    logical, intent(in), optional :: include_constraints
    character(len=*), intent(in), optional :: filename

    character(len=4096) :: header
    integer :: total_size, base_size, i, n_cont
    character(len=256) :: log_name

    if (present(include_constraints)) then
       this%log_include_constraints = include_constraints
    end if

    n_cont = nekotop_continuation%get_n_params()
    base_size = problem%get_log_size(this%log_include_constraints)
    this%log_extra_size = 0
    if (present(extra_headers)) this%log_extra_size = size(extra_headers)

    total_size = 1 + base_size + this%log_extra_size + 1 + n_cont
    call this%log_data%init(total_size)

    if (present(filename)) then
       log_name = trim(filename)
    else
       log_name = 'optimization_data.csv'
    end if

    call this%log_file%init(trim(log_name))

    header = 'iter, ' // &
         trim(problem%get_log_header(this%log_include_constraints))
    if (present(extra_headers)) then
       do i = 1, size(extra_headers)
          if (trim(extra_headers(i)) .eq. '') then
             call neko_error('some headers are empty')
          end if
          header = trim(header) // ', ' // trim(extra_headers(i))
       end do
    end if

    header = trim(header) // ', max_design_change'

    ! continuation parameters
    do i = 1, n_cont
       header = trim(header) // ', ' // &
            trim(nekotop_continuation%get_param_name(i))
    end do
    call this%log_file%set_header(trim(header))

    this%log_initialized = .true.
  end subroutine optimizer_init_log

  !> Write optimization log entry.
  !! @param[inout] this The optimizer object.
  !! @param[in] iter Current iteration number.
  !! @param[in] problem The problem object.
  !! @param[in] extra_values Extra log values appended after problem entries.
  subroutine optimizer_write_log(this, iter, problem, extra_values)
    class(optimizer_t), intent(inout) :: this
    integer, intent(in) :: iter
    class(problem_t), intent(in) :: problem
    real(kind=rp), intent(in), optional :: extra_values(:)
    integer :: base_size, offset, n_cont, i

    if (.not. this%log_initialized) return

    ! Number of continuation parameters
    n_cont = nekotop_continuation%get_n_params()

    base_size = problem%get_log_size(this%log_include_constraints)
    this%log_data%x(1) = real(iter, kind=rp)

    call problem%get_log_values( &
         this%log_data%x(2:1 + base_size), &
         this%log_include_constraints)

    offset = 2 + base_size
    if (present(extra_values)) then
       if (this%log_extra_size .eq. 0) then
          call neko_error('got extra values but no headers')
       end if
       this%log_data%x(offset:offset + size(extra_values) - 1) = extra_values
    end if

    ! Save maximum design change
    offset = offset + this%log_extra_size
    this%log_data%x(offset:offset + 1 - 1) = this%max_design_change

    ! Continuation parameter values
    offset = offset + 1
    do i = 1, n_cont
       this%log_data%x(offset + i) = &
            nekotop_continuation%params(i)%target
    end do

    call this%log_file%write(this%log_data)

  end subroutine optimizer_write_log

  ! ========================================================================== !
  ! IO Functions

  !> Construct a checkpoint filename from path/basename/format, applying the
  !! optimizer's stored `checkpoint_path`/`checkpoint_base`/`checkpoint_format`
  !! defaults for any component not given. Shared by `optimizer_load_checkpoint`
  !! and `optimizer_run`'s restart logic so the default path is built in one
  !! place.
  !! @param this The optimizer object.
  !! @param path The path where the checkpoint file is located.
  !! @param basename The base name of the checkpoint file.
  !! @param format The file format of the checkpoint file.
  function optimizer_checkpoint_filename(this, path, basename, format) &
       result(file_full)
    class(optimizer_t), intent(in) :: this
    character(len=*), intent(in), optional :: path
    character(len=*), intent(in), optional :: basename
    character(len=*), intent(in), optional :: format
    character(len=256) :: file_full
    character(len=:), allocatable :: checkpoint_format
    character(len=256) :: file_path, file_base, file_ext

    file_path = trim(this%checkpoint_path)
    file_base = trim(this%checkpoint_base)
    checkpoint_format = trim(this%checkpoint_format)

    if (present(path)) file_path = trim(path)
    if (present(basename)) file_base = trim(basename)
    if (present(format)) checkpoint_format = trim(format)

    if (len_trim(file_path) .eq. 0) then
       file_path = './'
    else if (file_path(len_trim(file_path):len_trim(file_path)) &
         .ne. '/') then
       file_path = trim(file_path) // '/'
    end if

    select case (trim(checkpoint_format))
    case ('h5', 'hdf5', 'hf5', 'hdf')
       file_ext = 'h5'
    case default
       call neko_error('optimizer: Unsupported checkpoint format: "' // &
            trim(checkpoint_format) // '"')
    end select

    write(file_full, '(4A)') &
         trim(file_path), trim(file_base), ".", trim(file_ext)

  end function optimizer_checkpoint_filename

  !> Save the optimizer checkpoint to a file.
  !! @param this The optimizer object.
  !! @param iter The current iteration number.
  !! @param design The design object.
  !! @param overwrite Whether to overwrite the file if it exists.
  !! @param path The path where the checkpoint file will be saved.
  !! @param basename The base name of the file to save the checkpoint.
  !! @param format The file format to use for the checkpoint file.
  subroutine optimizer_save_checkpoint(this, iter, design, overwrite, &
       path, basename, format)
    class(optimizer_t), intent(inout) :: this
    integer, intent(in) :: iter
    class(design_t), intent(inout) :: design
    logical, intent(in) :: overwrite
    character(len=*), intent(in), optional :: path
    character(len=*), intent(in), optional :: basename
    character(len=*), intent(in), optional :: format
    character(len=:), allocatable :: checkpoint_format
    character(len=256) :: file_path, file_base, file_ext, file_full
    character(len=LOG_SIZE) :: msg
    real(kind=rp) :: t_start, t_total
    logical :: exist

    call neko_log%section('Optimizer checkpoint')
    t_start = MPI_Wtime()

    ! Set default behaviour
    file_path = trim(this%checkpoint_path)
    file_base = trim(this%checkpoint_base)
    checkpoint_format = trim(this%checkpoint_format)

    ! Overwrite any user supplied components
    if (present(path)) file_path = trim(path)
    if (present(basename)) file_base = trim(basename)
    if (present(format)) checkpoint_format = trim(format)

    ! Make sure path is valid and exists
    if (len_trim(file_path) .eq. 0) then
       file_path = './'
    else if (file_path(len_trim(file_path):len_trim(file_path)) .ne. '/') then
       file_path = trim(file_path) // '/'
    end if

    inquire(file=file_path, exist=exist)
    if (.not. exist) then
       call execute_command_line('mkdir -p "' // trim(file_path) // '"')
    end if

    select case (trim(checkpoint_format))
    case ('h5', 'hdf5', 'hf5', 'hdf')
       file_ext = 'h5'
    case default
       call neko_error('optimizer: Unsupported checkpoint format: "' // &
            trim(checkpoint_format) // '"')
    end select

    ! Construct the full filename based on overwrite flag
    if (overwrite) then
       write(file_full, '(4A)') &
            trim(file_path), trim(file_base), ".", trim(file_ext)
    else
       write(file_full, '(3A,I5.5,2A)') &
            trim(file_path), trim(file_base), "_", iter, ".", trim(file_ext)
    end if

    call neko_log%message('Save general optimizer components')
    select case (trim(file_ext))
    case ('h5', 'hdf5', 'hf5')
       call optimizer_save_checkpoint_hdf5(this, file_full, iter, overwrite)
    case default
       call neko_error('optimizer: Unsupported checkpoint format: "' // &
            trim(file_ext) // '"')
    end select

    call neko_log%message('Saving components of ' // this%optimizer_type)
    call this%save_checkpoint_components(file_full, overwrite)

    call neko_log%message('Save design checkpoint')
    call design%save_checkpoint(file_full, overwrite)

    t_total = MPI_Wtime() - t_start
    write(msg, '(A,F6.2)') "Checkpoint time: ", t_total
    call neko_log%end_section(msg)

  end subroutine optimizer_save_checkpoint

  !> Load the optimizer checkpoint from a file based on file suffix.
  !! If `filename` is not given, the file is located from `path`/`basename`/
  !! `format`, defaulting to `this%checkpoint_path`/`checkpoint_base`/
  !! `checkpoint_format` exactly as `optimizer_save_checkpoint` does when
  !! saving in overwrite mode.
  !! @param this The optimizer object.
  !! @param filename The name of the file to load the checkpoint from.
  !! @param iter The iteration number read from the checkpoint.
  !! @param design The design object.
  !! @param path The path where the checkpoint file is located.
  !! @param basename The base name of the file to load the checkpoint from.
  !! @param format The file format to use for the checkpoint file.
  subroutine optimizer_load_checkpoint(this, filename, iter, design, &
       path, basename, format)
    class(optimizer_t), intent(inout) :: this
    character(len=*), intent(in), optional :: filename
    integer, intent(out) :: iter
    class(design_t), intent(inout) :: design
    character(len=*), intent(in), optional :: path
    character(len=*), intent(in), optional :: basename
    character(len=*), intent(in), optional :: format
    character(len=256) :: file_full
    character(len=12) :: suffix

    if (present(filename)) then
       file_full = trim(filename)
    else
       file_full = optimizer_checkpoint_filename(this, path, basename, format)
    end if

    ! Get the file extension
    call filename_suffix(trim(file_full), suffix)

    select case (trim(suffix))
    case ('h5', 'hdf5', 'hf5')
       call optimizer_load_checkpoint_hdf5(this, trim(file_full), iter)
    case default
       call neko_error('optimizer: Unsupported checkpoint format: "' // &
            trim(suffix) // '"')
    end select

    call this%load_checkpoint_components(trim(file_full))
    call design%load_checkpoint(trim(file_full))

    ! Set the current iteration to the loaded iteration
    this%current_iteration = iter

    if (pe_rank .eq. 0) then
       write(*,*) 'Restarted simulation from checkpoint.'
       write(*,*) '    Checkpoint file: "', trim(file_full), '"'
       write(*,*) '    Iteration      : ', this%current_iteration
    end if

  end subroutine optimizer_load_checkpoint

  ! ========================================================================== !
  ! Dummy implementations for module procedures

#if !HAVE_HDF5
  module subroutine optimizer_save_checkpoint_hdf5(object, filename, iter, &
       overwrite)
    class(optimizer_t), intent(inout) :: object
    character(len=*), intent(in) :: filename
    integer, intent(in) :: iter
    logical, intent(in), optional :: overwrite
    call neko_error('optimizer: HDF5 support not enabled rebuild with ' // &
         'HAVE_HDF5')
  end subroutine optimizer_save_checkpoint_hdf5

  module subroutine optimizer_load_checkpoint_hdf5(object, filename, iter)
    class(optimizer_t), intent(inout) :: object
    character(len=*), intent(in) :: filename
    integer, intent(out) :: iter
    call neko_error('optimizer: HDF5 support not enabled rebuild with ' // &
         'HAVE_HDF5')
  end subroutine optimizer_load_checkpoint_hdf5
#endif

end module optimizer
