!> @file time_window_tester.f90
!! Driver asserting that objective time windows behave correctly when the
!! window does not span the whole simulation.
!!
!! Two independent properties are checked, selected from the case file's
!! `optimization.time_window_test` object:
!!
!!  - **Run-length invariance.** With `end_times` listing more than one time,
!!    the same case is run to each of them and every objective must report the
!!    same value. An objective whose window is closed and lies inside all the
!!    runs measures the same interval in each, so its value cannot depend on
!!    how long the run continued afterwards.
!!
!!  - **Window-form equivalence.** With `pairwise` set, the objectives are read
!!    as consecutive pairs which are set up to cover the same samples by
!!    different means -- an open-ended window against the equivalent closed
!!    one, for instance -- and each pair must agree.
program time_window_tester

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
  use num_types, only: rp, dp
  use vector, only: vector_t
  use objectives_user, only: objectives_user_setup
  implicit none

  ! JSON related arguments
  integer :: argc, stem_end
  character(len=256) :: parameter_file
  !> Parameter file name without its '.case' suffix.
  character(len=:), allocatable :: case_stem
  type(json_file) :: parameters, design_parameters

  !> The simulation we are working with
  type(simulation_t) :: sim
  !> The design type
  type(brinkman_design_t) :: des
  !> The problem type
  type(problem_t) :: prob

  !> End times requested for the unsteady runs, one run each.
  real(kind=dp), allocatable :: end_times(:)
  !> Objective values, one row per run.
  real(kind=rp), allocatable :: results(:, :)
  type(vector_t) :: values

  real(kind=rp) :: tolerance, reference, difference
  logical :: pairwise, failed, found
  integer :: n_runs
  integer :: n_objectives, n_declared, i, j
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

  ! The tests here share a working directory and may run concurrently, so each
  ! case keeps its state recovery files in a directory named after it.
  stem_end = index(parameter_file, '.case', back = .true.)
  if (stem_end .eq. 0) then
     case_stem = trim(parameter_file)
  else
     case_stem = parameter_file(1:stem_end - 1)
  end if
  call parameters%update('state_recovery.path', case_stem // '_checkpoints/', &
       found)

  call json_get(parameters, 'optimization.design', design_parameters)

  call json_get_or_default(parameters, &
       'optimization.time_window_test.tolerance', tolerance, 1.0e-12_rp)
  call json_get_or_default(parameters, &
       'optimization.time_window_test.pairwise', pairwise, .false.)

  ! -------------------------------------------------------------------------- !
  ! Initialization of the components

  ! Register the user-defined scalar inflow before the case is built. Only
  ! the cases that ask for it in their boundary conditions will call it.
  call objectives_user_setup(sim%neko_case%user)

  call sim%init(parameters)
  call des%init(design_parameters, sim)
  call prob%init(parameters, des, sim)

  if (.not. sim%unsteady) then
     call neko_error('The time window only applies to unsteady problems; ' // &
          'set "unsteady": true in the case file.')
  end if

  ! `problem_t` appends an internal augmented-Lagrangian objective of its
  ! own, so only the objectives declared in the case file are checked.
  n_objectives = prob%get_n_objectives()
  call parameters%info('optimization.objectives', n_children = n_declared)
  if (n_declared .lt. 1) call neko_error('No objectives to check')
  if (n_declared .gt. n_objectives) then
     call neko_error('Fewer objectives were built than the case declares')
  end if

  write(log_buf, '(A,I0,A,I0,A)') 'Objectives under test: ', n_declared, &
       ' (of ', n_objectives, ' built)'
  call neko_log%message(trim(log_buf))

  if (parameters%valid_path('optimization.time_window_test.end_times')) then
     call json_get(parameters, 'optimization.time_window_test.end_times', &
          end_times)
  else
     allocate(end_times(1))
     end_times(1) = sim%neko_case%time%end_time
  end if
  n_runs = size(end_times)

  ! -------------------------------------------------------------------------- !
  ! Run the case once per requested end time

  allocate(results(n_runs, n_objectives))
  call values%init(n_objectives)

  do i = 1, n_runs
     sim%neko_case%time%end_time = end_times(i)

     write(log_buf, '(A,E13.6)') 'Time window test, run to t = ', &
          end_times(i)
     call neko_log%section(trim(log_buf))

     call prob%compute(des, sim)
     call prob%get_all_objective_values(values)
     results(i, :) = values%x

     call neko_log%end_section()
  end do

  ! -------------------------------------------------------------------------- !
  ! Report and check

  failed = .false.

  call neko_log%section('Objective values')
  do i = 1, n_runs
     do j = 1, n_declared
        write(log_buf, '(A,I0,A,E13.6,A,I0,A,E22.15)') 'run ', i, &
             ' to t = ', end_times(i), ', objective ', j, ' = ', results(i, j)
        call neko_log%message(trim(log_buf))
     end do
  end do
  call neko_log%end_section()

  ! Run-length invariance: every run must agree with the first.
  do i = 2, n_runs
     do j = 1, n_declared
        reference = results(1, j)
        difference = relative_difference(results(i, j), reference)
        if (difference .le. tolerance) cycle

        failed = .true.
        write(log_buf, '(A,I0,A,E22.15,A,E22.15,A,E13.6)') &
             'Objective ', j, ' moved with the run length: ', reference, &
             ' vs ', results(i, j), ', rel. diff ', difference
        call neko_log%error(trim(log_buf))
     end do
  end do

  ! Window-form equivalence: consecutive objectives must agree in pairs.
  if (pairwise) then
     if (mod(n_declared, 2) .ne. 0) then
        call neko_error('A pairwise check needs an even number of objectives')
     end if

     do j = 1, n_declared, 2
        reference = results(1, j)
        difference = relative_difference(results(1, j + 1), reference)
        if (difference .le. tolerance) cycle

        failed = .true.
        write(log_buf, '(A,I0,A,I0,A,E22.15,A,E22.15,A,E13.6)') &
             'Objectives ', j, ' and ', j + 1, ' disagree: ', reference, &
             ' vs ', results(1, j + 1), ', rel. diff ', difference
        call neko_log%error(trim(log_buf))
     end do
  end if

  if (failed) then
     call neko_error('Objective values are not invariant to the time window')
  end if

  call neko_log%message('Time window check passed')

  ! -------------------------------------------------------------------------- !
  ! Clean up

  call values%free()
  deallocate(results)
  deallocate(end_times)
  call prob%free()
  call des%free()
  call sim%free()

  call neko_finalize()

contains

  !> Difference of two values relative to the larger magnitude, falling back to
  !! the absolute difference when both are effectively zero.
  !! @param a First value.
  !! @param b Second value.
  !! @return The relative difference.
  function relative_difference(a, b) result(difference)
    real(kind=rp), intent(in) :: a
    real(kind=rp), intent(in) :: b
    real(kind=rp) :: difference
    real(kind=rp) :: scale

    scale = max(abs(a), abs(b))
    if (scale .le. tiny(0.0_rp)) then
       difference = abs(a - b)
    else
       difference = abs(a - b) / scale
    end if
  end function relative_difference

end program time_window_tester
