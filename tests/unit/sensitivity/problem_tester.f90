program problem_tester

  use simulation_m, only: simulation_t
  use brinkman_design, only: brinkman_design_t
  use problem, only : problem_t

  ! Standard modules shared by most of our tests
  use neko, only: neko_init, neko_finalize
  use json_module, only: json_file
  use json_utils, only: json_get, json_get_or_default
  use json_utils_ext, only: json_read_file
  use utils, only: neko_error
  use neko_top, only: neko_top_register_types
  use mask_ops, only: mask_exterior_const
  use neko_config, only: NEKO_BCKND_DEVICE
  use device, only: device_memcpy, DEVICE_TO_HOST, HOST_TO_DEVICE

  ! Modules specific to this test
  use num_types, only: rp
  use vector, only: vector_t
  use matrix, only: matrix_t
  use math, only: abscmp
  use sensitivity, only: compute_sensitivity, fd_read_strict_options, &
       fd_read_perturbations, fd_assertion_skipped, FD_SKIP_EXIT_CODE, &
       fd_read_targets, fd_target_t
  use fd_criterion, only: fd_strict_options_t
  use logger, only: neko_log, LOG_SIZE
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

  ! Test specific variables. The tolerance may be overridden per case via the
  ! optional JSON key `optimization.fd_test_tolerance`; it defaults to a tight
  ! round-off floor suited to linear functionals (e.g. the volume constraint).
  ! PDE-coupled objectives, whose finite-difference floor is set by the
  ! discretisation/steady-state, set a looser value in their case file.
  real(kind=rp) :: tolerance
  !> The sweep this tier has always run. Kept as the default so that the
  !! logged CSVs stay comparable run to run; the strict criterion asks for a
  !! geometric sweep instead, because its order estimate is a statement about
  !! successive ratios and these alternate 5, 2, 5, 2.
  real(kind=rp), parameter :: default_perturbations(8) = [ &
       5e-1_rp, 1e-1_rp, 5e-2_rp, 1e-2_rp, 5e-3_rp, 1e-3_rp, 5e-4_rp, 1e-4_rp]
  real(kind=rp), allocatable :: perturbations(:)
  !> Settings of the strict assertion criterion. Disabled unless the case
  !! file asks for it, so the historical assertion stays the default.
  type(fd_strict_options_t) :: strict_options

  !> The analytic sensitivity of each target, read out before any sweep
  !! perturbs the design.
  type(vector_t), allocatable :: target_sensitivities(:)
  type(matrix_t) :: constraint_sensitivity

  integer :: i_max

  !> What this run checks: the weighted objective total, individual
  !! constraints, or any mixture of the two. Defaults to the historical
  !! single target when the case file names none.
  type(fd_target_t), allocatable :: targets(:)
  integer :: it, j, n_targets
  !> True once the constraint sensitivity matrix has been filled, so that a
  !! run checking several constraints reads it once rather than per target.
  logical :: constraints_read

  ! -------------------------------------------------------------------------- !
  ! Initialize the Neko environment

  call neko_init()
  call neko_top_register_types()

  ! -------------------------------------------------------------------------- !
  ! Read the parameters file as the first terminal argument

  argc = command_argument_count()
  if (argc .lt. 1) call neko_error('Missing parameter file')
  call get_command_argument(1, parameter_file)

  ! Read the parameters file
  parameters = json_read_file(trim(parameter_file))
  call json_get(parameters, 'optimization.design', design_parameters)
  call json_get_or_default(parameters, 'optimization.fd_test_tolerance', &
       tolerance, 1e-5_rp)

  ! This tier always takes a one-sided difference, hence the .false.
  call fd_read_strict_options(parameters, .false., strict_options)
  if (strict_options%enabled) then
     call fd_read_perturbations(parameters, .true., perturbations)
  else
     allocate(perturbations(size(default_perturbations)))
     perturbations = default_perturbations
  end if

  ! -------------------------------------------------------------------------- !
  ! Initialization of the components

  call sim%init(parameters)
  call des%init(design_parameters, sim)
  call prob%init(parameters, des, sim)

  ! -------------------------------------------------------------------------- !
  ! Determine what to check. A case may name any mixture of the weighted
  ! objective total and individual constraints; a case that names nothing
  ! falls back to the historical dispatch.

  call fd_read_targets(parameters, prob, targets)
  n_targets = size(targets)

  ! -------------------------------------------------------------------------- !
  ! Compute the sensitivity with our method

  call prob%compute(des, sim)

  ! Diagnostic only: print each individual objective's own value (the
  ! weighted-sum total the FD check gates on hides this). Cheap and
  ! harmless for every existing single-objective case; exists so a
  ! multi-objective case (e.g. a Pe study) can be inspected from the log
  ! alone instead of new throwaway Fortran per investigation.
  block
    type(vector_t) :: obj_values
    character(len=LOG_SIZE) :: log_buf
    integer :: obj_i
    if (prob%get_n_objectives() .gt. 1) then
       call obj_values%init(prob%get_n_objectives())
       call prob%get_all_objective_values(obj_values)
       do obj_i = 1, prob%get_n_objectives()
          write(log_buf, '(A,I0,A,ES16.9)') ' Objective ', obj_i, &
               ' value = ', obj_values%x(obj_i)
          call neko_log%message(log_buf)
       end do
       call obj_values%free()
    end if
  end block

  call prob%compute_sensitivity(des, sim)

  ! Read every target's analytic sensitivity out *before* running any sweep:
  ! a sweep perturbs the design and re-evaluates the problem, and only
  ! `compute_sensitivity` refills the sensitivity fields it would read.
  allocate(target_sensitivities(n_targets))
  constraints_read = .false.
  do it = 1, n_targets
     call target_sensitivities(it)%init(des%size())

     if (targets(it)%is_objective) then
        call des%get_sensitivity(target_sensitivities(it))
     else
        if (.not. constraints_read) then
           call constraint_sensitivity%init(prob%get_n_constraints(), &
                des%size())
           call prob%get_constraint_sensitivities(constraint_sensitivity)

           ! The layout is (n_constraints, n_design), so constraint i is a
           ! *row*. Asserted at the call site as well as inside the getter
           ! because the transposed allocation this replaces was silent at
           ! one constraint and corrupted memory at two.
           if (constraint_sensitivity%get_nrows() .ne. &
                prob%get_n_constraints() .or. &
                constraint_sensitivity%get_ncols() .ne. des%size()) then
              call neko_error('The constraint sensitivity matrix is not ' &
                   // 'shaped (n_constraints, n_design)')
           end if
           constraints_read = .true.
        end if

        do j = 1, des%size()
           target_sensitivities(it)%x(j) = &
                constraint_sensitivity%x(targets(it)%constraint_index, j)
        end do
        call target_sensitivities(it)%copy_from(HOST_TO_DEVICE, &
             sync = .true.)
     end if

     call des%convert_to_directional_derivative(target_sensitivities(it))
  end do

  call des%write(1)
  ! --------------------------------------
  ! Reset the simulation
  call sim%reset()

  ! -------------------------------------------------------------------------- !
  ! Loop over the perturbations and compare the finite difference estimate with
  ! the sensitivity computed by our method, once per target.

  do it = 1, n_targets
     if (NEKO_BCKND_DEVICE .eq. 1) then
        call device_memcpy(target_sensitivities(it)%x, &
             target_sensitivities(it)%x_d, target_sensitivities(it)%size(), &
             DEVICE_TO_HOST, .true.)
     end if

     i_max = maxloc(abs(target_sensitivities(it)%x), dim=1)

     call compute_sensitivity(prob, sim, des, target_sensitivities(it), &
          i_max, perturbations, tolerance, trim(parameter_file), &
          targets(it)%is_objective, sim%fluid%gs_Xh, &
          strict_options = strict_options, &
          constraint_index = targets(it)%constraint_index, &
          name_suffix = trim(targets(it)%suffix))
  end do

  ! -------------------------------------------------------------------------- !
  ! Clean up the components

  if (allocated(perturbations)) deallocate(perturbations)

  do it = 1, n_targets
     call target_sensitivities(it)%free()
  end do
  deallocate(target_sensitivities)
  deallocate(targets)
  call constraint_sensitivity%free()

  call prob%free()
  call des%free()
  call sim%free()

  ! Finalize the Neko environment
  call neko_finalize()

  ! A build that skipped the assertion made no gradient check at all, and a
  ! zero exit would be reported by CTest as a pass. Exit with the skip code
  ! instead -- `SKIP_RETURN_CODE` in this test's CMakeLists turns it into a
  ! reported SKIP.
  if (fd_assertion_skipped()) stop FD_SKIP_EXIT_CODE

end program problem_tester
