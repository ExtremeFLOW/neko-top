!> @file step_count_tester.f90
!! Driver asserting that a run stops at its `end_time`, on both the steady and
!! the unsteady path.
!!
!! The case is run once through each branch of `problem_t%compute`, and each
!! run must take exactly `(end_time - start_time) / dt` timesteps and finish
!! at `end_time`.
!!
!! The case file must place `end_time` on a step boundary, because that is
!! where the check bites. The time is built by repeated addition of `dt`, so
!! on a boundary it lands a few ulps either side of `end_time`; when it lands
!! short, a time loop comparing the two without a tolerance takes one step
!! too many. With `dt = 0.005` and `end_time = 0.05` the unsteady path used
!! to run eleven steps and stop at `t = 0.055`.
!!
!! The adjoint of the unsteady run is then solved once, and must not write a
!! job-limit checkpoint. With `output_at_end` off, the adjoint's finalisation
!! writes one when it judges the run to have stopped short of `end_time`. The
!! adjoint replays the forward run's step count, so its clock ends where the
!! forward run's did, and needs the same tolerance to recognise the run as
!! complete. The case file must therefore turn `output_at_end` off.
program step_count_tester

  use simulation_m, only: simulation_t
  use brinkman_design, only: brinkman_design_t
  use problem, only: problem_t

  ! Standard modules shared by most of our tests
  use neko, only: neko_init, neko_finalize
  use json_module, only: json_file
  use json_utils, only: json_get, json_get_or_default
  use json_utils_ext, only: json_read_file
  use utils, only: neko_error
  use neko_top, only: neko_top_register_types
  use logger, only: neko_log

  ! Modules specific to this test
  use num_types, only: dp
  use file, only: file_t
  use objectives_user, only: objectives_user_setup
  implicit none

  ! JSON related arguments
  integer :: argc
  character(len=256) :: parameter_file
  type(json_file) :: parameters, design_parameters

  !> The simulation we are working with
  type(simulation_t) :: sim
  !> The design type
  type(brinkman_design_t) :: des
  !> The problem type
  type(problem_t) :: prob

  !> How far, as a fraction of a timestep, a time may sit from a step
  !! boundary and still count as being on it.
  real(kind=dp), parameter :: step_tolerance = 1.0e-6_dp

  real(kind=dp) :: start_time, end_time, dt, n_steps_exact, final_time
  integer :: n_expected, n_taken, i
  logical :: failed, output_at_end, joblimit_written
  !> The job-limit checkpoint the adjoint would write, named as it names it.
  type(file_t) :: joblimit_file
  character(len=1024) :: joblimit_fname
  character(len=:), allocatable :: checkpoint_format
  ! Wider than LOG_SIZE: these lines carry full-precision values.
  character(len=256) :: log_buf

  ! -------------------------------------------------------------------------- !
  ! Initialize the Neko environment

  call neko_init()
  call neko_top_register_types()

  ! -------------------------------------------------------------------------- !
  ! Read the parameters file as the first terminal argument

  argc = command_argument_count()
  if (argc .lt. 1) call neko_error('Missing parameter file')
  call get_command_argument(1, parameter_file)

  parameters = json_read_file(trim(parameter_file))
  call json_get(parameters, 'optimization.design', design_parameters)

  ! -------------------------------------------------------------------------- !
  ! Initialization of the components

  ! Register the user-defined inflow before the case is built.
  call objectives_user_setup(sim%neko_case%user)

  call sim%init(parameters)
  call des%init(design_parameters, sim)
  call prob%init(parameters, des, sim)

  ! The unsteady path needs the state recovery that only an unsteady problem
  ! sets up; the steady path is reached by clearing the flag afterwards.
  if (.not. sim%unsteady) then
     call neko_error('The unsteady path needs state recovery; ' // &
          'set "unsteady": true in the case file.')
  end if

  start_time = sim%neko_case%time%start_time
  end_time = sim%neko_case%time%end_time
  dt = sim%neko_case%time%dt

  n_steps_exact = (end_time - start_time) / dt
  n_expected = nint(n_steps_exact)
  if (abs(n_steps_exact - real(n_expected, kind=dp)) .gt. step_tolerance) then
     call neko_error('This test needs an end_time on a step boundary.')
  end if

  call json_get_or_default(parameters, 'case.output_at_end', output_at_end, &
       .true.)
  if (output_at_end) then
     call neko_error('This test needs "output_at_end": false.')
  end if

  write(log_buf, '(A,E22.15,A,E22.15,A,I0,A)') 'Expecting a run to t = ', &
       end_time, ' with dt = ', dt, ' to take ', n_expected, ' steps'
  call neko_log%message(trim(log_buf))

  ! -------------------------------------------------------------------------- !
  ! Run the case through each path and check where it stopped

  ! The steady path goes first, as in `time_window_tester`. Both paths zero
  ! the step counter before their loop, so a path that silently took no steps
  ! reports zero rather than the other path's count.
  failed = .false.
  do i = 1, 2
     sim%unsteady = i .eq. 2

     write(log_buf, '(A,A,A)') 'Step count test, ', &
          trim(path_name(sim%unsteady)), ' run'
     call neko_log%section(trim(log_buf))

     call prob%compute(des, sim)

     n_taken = sim%n_timesteps
     final_time = sim%neko_case%time%t

     call neko_log%end_section()

     write(log_buf, '(A,A,A,I0,A,E22.15)') 'The ', &
          trim(path_name(sim%unsteady)), ' run took ', n_taken, &
          ' steps and stopped at t = ', final_time
     call neko_log%message(trim(log_buf))

     if (n_taken .eq. n_expected .and. &
          abs(final_time - end_time) .le. step_tolerance * dt) cycle

     ! A warning rather than `neko_log%error`, which would abort here and
     ! leave the other path unreported.
     failed = .true.
     write(log_buf, '(A,A,A,I0,A,E22.15,A,I0,A,E22.15)') 'The ', &
          trim(path_name(sim%unsteady)), ' run took ', n_taken, &
          ' steps to t = ', final_time, ' instead of ', n_expected, &
          ' to t = ', end_time
     call neko_log%warning(trim(log_buf))
  end do

  ! Put the flag back the way the case file declared it.
  sim%unsteady = .true.

  ! -------------------------------------------------------------------------- !
  ! Solve the adjoint of the unsteady run, the last one made, and check it
  ! treated the run as complete

  ! Build the name the way `simulation_adjoint_joblimit_chkp` does, so that
  ! this looks for the file it would actually write.
  call json_get_or_default(parameters, 'case.checkpoint_format', &
       checkpoint_format, 'chkp')
  if (checkpoint_format .eq. 'hdf5') then
     call joblimit_file%init(sim%neko_case%output_directory // 'joblimit.h5')
  else
     call joblimit_file%init(sim%neko_case%output_directory // 'joblimit.chkp')
  end if
  joblimit_fname = joblimit_file%file_type%get_next_output_fname()
  call joblimit_file%free()

  ! Clear any left over from an earlier run of this test.
  call delete_file(joblimit_fname)

  call neko_log%section('Step count test, adjoint run')
  call prob%compute_sensitivity(des, sim)
  call neko_log%end_section()

  inquire(file = trim(joblimit_fname), exist = joblimit_written)
  if (joblimit_written) then
     failed = .true.
     call neko_log%warning('The adjoint wrote ' // trim(joblimit_fname) // &
          ', so it judged the run to have stopped short of end_time')
  else
     call neko_log%message('The adjoint wrote no job-limit checkpoint')
  end if

  if (failed) then
     call neko_error('A run or its adjoint mishandled end_time; see the log')
  end if

  call neko_log%message('Step count check passed')

  ! -------------------------------------------------------------------------- !
  ! Clean up

  call prob%free()
  call des%free()
  call sim%free()

  call neko_finalize()

contains

  !> Name of the path a run is taken through, for the log.
  !! @param unsteady Whether the run uses the unsteady path.
  !! @return `'unsteady'` or `'steady'`.
  function path_name(unsteady) result(name)
    logical, intent(in) :: unsteady
    character(len=8) :: name

    if (unsteady) then
       name = 'unsteady'
    else
       name = 'steady'
    end if
  end function path_name

  !> Delete a file, if it exists.
  !! @param fname Name of the file.
  subroutine delete_file(fname)
    character(len=*), intent(in) :: fname
    logical :: exists
    integer :: unit, stat

    inquire(file = trim(fname), exist = exists)
    if (.not. exists) return

    open(newunit = unit, file = trim(fname), iostat = stat, status = 'old')
    if (stat .eq. 0) close(unit, status = 'delete')
  end subroutine delete_file

end program step_count_tester
