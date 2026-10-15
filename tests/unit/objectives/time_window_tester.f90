!> @file time_window_tester.f90
!! Driver asserting that objective time windows behave correctly when the
!! window does not span the whole simulation.
!!
!! Three independent properties are checked, selected from the case file's
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
!!
!!  - **Steady/unsteady equivalence.** With `compare_steady` set, the case is
!!    also run through the steady path -- `problem_t%compute` evaluating each
!!    objective on the final field rather than accumulating it over the run
!!    -- and must agree with the unsteady runs. The case file is responsible
!!    for choosing windows that make the two comparable; see the readme in
!!    this directory. `require_steady_state` additionally asserts that each
!!    run actually reached a steady state, which is what makes a window
!!    average comparable to a single converged field in the first place.
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
  use num_types, only: rp
  use vector, only: vector_t
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

  !> End times requested for the unsteady runs, one run each.
  real(kind=rp), allocatable :: end_times(:)
  !> End time of every run, in the order they are run.
  real(kind=rp), allocatable :: run_end_times(:)
  !> Which path each run is taken through.
  logical, allocatable :: run_unsteady(:)
  !> Objective values, one row per run.
  real(kind=rp), allocatable :: results(:, :)
  type(vector_t) :: values

  real(kind=rp) :: tolerance, reference, difference
  logical :: pairwise, compare_steady, require_steady_state, failed
  integer :: n_runs, n_unsteady, first_unsteady
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
  call json_get(parameters, 'optimization.design', design_parameters)

  call json_get_or_default(parameters, &
       'optimization.time_window_test.tolerance', tolerance, 1.0e-12_rp)
  call json_get_or_default(parameters, &
       'optimization.time_window_test.pairwise', pairwise, .false.)
  call json_get_or_default(parameters, &
       'optimization.time_window_test.compare_steady', compare_steady, .false.)
  call json_get_or_default(parameters, &
       'optimization.time_window_test.require_steady_state', &
       require_steady_state, .false.)

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
  n_unsteady = size(end_times)

  ! The steady run, when asked for, is one more run of the same case at the
  ! first end time, with `problem_t%compute` sent down its steady branch.
  !
  ! It goes first on purpose. Every run leaves its values in the same
  ! objectives, so a path that silently failed to evaluate them would inherit
  ! the previous run's numbers and agree with it. Going first, such a path
  ! reports the objectives' initial zero instead, and the comparison catches
  ! it. The unsteady path zeroes them itself before accumulating, so it is
  ! safe anywhere in the order.
  n_runs = n_unsteady
  if (compare_steady) n_runs = n_runs + 1

  allocate(run_end_times(n_runs))
  allocate(run_unsteady(n_runs))
  first_unsteady = n_runs - n_unsteady + 1
  run_end_times(first_unsteady:n_runs) = end_times
  run_unsteady(first_unsteady:n_runs) = .true.
  if (compare_steady) then
     run_end_times(1) = end_times(1)
     run_unsteady(1) = .false.
  end if

  ! -------------------------------------------------------------------------- !
  ! Run the case once per requested end time

  allocate(results(n_runs, n_objectives))
  call values%init(n_objectives)

  do i = 1, n_runs
     sim%unsteady = run_unsteady(i)
     sim%neko_case%time%end_time = run_end_times(i)

     write(log_buf, '(A,A,A,E13.6)') 'Time window test, ', &
          trim(path_name(run_unsteady(i))), ' run to t=', run_end_times(i)
     call neko_log%section(trim(log_buf))

     call prob%compute(des, sim)
     call prob%get_all_objective_values(values)
     results(i, :) = values%x

     ! A comparison across the two paths only means anything once the run has
     ! settled. Say so directly rather than leaving it to be inferred from a
     ! value mismatch: the `steady` simulation component freezes the fluid on
     ! convergence, so an unfrozen fluid here is a run that never got there.
     if (require_steady_state .and. .not. sim%neko_case%fluid%freeze) then
        call neko_error('The run finished without reaching a steady ' // &
             'state; lengthen it, or loosen the "steady" component''s tol.')
     end if

     call neko_log%end_section()
  end do

  ! Put the flag back the way the case file declared it; the driver refuses to
  ! start unless that was `.true.`.
  sim%unsteady = .true.

  ! -------------------------------------------------------------------------- !
  ! Report and check

  failed = .false.

  call neko_log%section('Objective values')
  do i = 1, n_runs
     do j = 1, n_declared
        write(log_buf, '(A,I0,A,A,A,E13.6,A,I0,A,E22.15)') 'run ', i, &
             ' (', trim(path_name(run_unsteady(i))), ') to t=', &
             run_end_times(i), ', objective ', j, ' = ', results(i, j)
        call neko_log%message(trim(log_buf))
     end do
  end do
  call neko_log%end_section()

  ! Run-length and steady/unsteady invariance: every run must agree with the
  ! first.
  do i = 2, n_runs
     do j = 1, n_declared
        reference = results(1, j)
        difference = relative_difference(results(i, j), reference)
        if (difference .le. tolerance) cycle

        failed = .true.
        if (run_unsteady(i) .eqv. run_unsteady(1)) then
           write(log_buf, '(A,I0,A,E22.15,A,E22.15,A,E13.6)') &
                'Objective ', j, ' moved with the run length: ', reference, &
                ' vs ', results(i, j), ', rel. diff ', difference
        else
           write(log_buf, '(A,I0,A,E22.15,A,E22.15,A,E13.6)') &
                'Objective ', j, ' differs between the steady and unsteady ' &
                // 'paths: ', reference, ' vs ', results(i, j), &
                ', rel. diff ', difference
        end if
        call neko_log%error(trim(log_buf))
     end do
  end do

  ! Window-form equivalence: consecutive objectives must agree in pairs.
  if (pairwise) then
     if (mod(n_declared, 2) .ne. 0) then
        call neko_error('A pairwise check needs an even number of objectives')
     end if

     do j = 1, n_declared, 2
        reference = results(first_unsteady, j)
        difference = relative_difference( &
             results(first_unsteady, j + 1), reference)
        if (difference .le. tolerance) cycle

        failed = .true.
        write(log_buf, '(A,I0,A,I0,A,E22.15,A,E22.15,A,E13.6)') &
             'Objectives ', j, ' and ', j + 1, ' disagree: ', reference, &
             ' vs ', results(first_unsteady, j + 1), ', rel. diff ', &
             difference
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
  deallocate(run_unsteady)
  deallocate(run_end_times)
  deallocate(end_times)
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
