program problem_tester

  use simulation_m, only: simulation_t
  use brinkman_design, only: brinkman_design_t
  use problem, only : problem_t

  ! Standard modules shared by most of our tests
  use json_module, only: json_file
  use json_utils, only: json_get, json_get_or_default
  use json_utils_ext, only: json_read_file
  use utils, only: neko_error
  use neko_top, only: neko_top_register_types
  use neko, only: neko_init, neko_finalize
  use mask_ops, only: mask_exterior_const
  use neko_config, only: NEKO_BCKND_DEVICE
  use device, only: device_memcpy, DEVICE_TO_HOST, HOST_TO_DEVICE

  ! Modules specific to this test
  use num_types, only: rp, i8
  use vector, only: vector_t
  use matrix, only: matrix_t
  use math, only: abscmp, glmax
  use comm, only: pe_rank
  use sensitivity, only: compute_sensitivity, &
       compute_sensitivity_directional, fd_read_perturbations, &
       fd_read_central_difference, fd_read_mode, fd_read_probe_index, &
       fd_resolve_probe_index, fd_read_strict_options, &
       fd_assertion_skipped, FD_SKIP_EXIT_CODE, fd_read_targets, fd_target_t
  use fd_criterion, only: fd_strict_options_t
  use user, only: user_setup
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
  ! optional JSON key `optimization.fd_test_tolerance`; it defaults to a value
  ! appropriate for the fully steady-state-converged regression cases.
  !
  ! The perturbation sweep, the one-sided/central choice and the dof/directional
  ! mode are likewise optional per case (`optimization.fd_test_perturbations`,
  ! `optimization.fd_test_central_difference` and `optimization.fd_test_mode`,
  ! each also overridable from the environment for a sweep driven across
  ! several cases at once), and default to the historical four-point one-sided
  ! single-dof sweep. See `fd_read_perturbations`, `fd_read_central_difference`
  ! and `fd_read_mode` in tests/shared/sensitivity.f90.
  real(kind=rp) :: tolerance
  real(kind=rp), allocatable :: perturbations(:)
  logical :: use_central
  !> Settings of the strict assertion criterion. Disabled unless the case
  !! file asks for it, so the historical assertion stays the default.
  type(fd_strict_options_t) :: strict_options
  !> True to perturb the whole design along the normalised sensitivity
  !! direction (the Taylor test) rather than probing a single dof.
  logical :: fd_directional
  !> A fixed global design index to probe, from NEKO_TOP_FD_PROBE_INDEX. The
  !! default probe is the argmax of |sensitivity|, which *moves* whenever the
  !! sensitivity field moves, so two runs of the same case can otherwise
  !! silently differentiate with respect to two different design variables.
  integer(kind=i8) :: probe_dof
  logical :: probe_dof_set

  !> The analytic sensitivity of each target, read out before any sweep
  !! perturbs the design.
  type(vector_t), allocatable :: target_sensitivities(:)
  type(matrix_t) :: constraint_sensitivity

  integer :: i_max
  integer :: i_local
  real(kind=rp) :: local_abs_max, global_abs_max
  real(kind=rp) :: key_local, key_global, key_arr(1)

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
       tolerance, 1e-3_rp)
  ! Order matters: the expected truncation order defaults from the
  ! one-sided/central choice, and the default sweep from whether the strict
  ! criterion is in use.
  call fd_read_central_difference(parameters, use_central)
  call fd_read_strict_options(parameters, use_central, strict_options)
  call fd_read_perturbations(parameters, strict_options%enabled, &
       perturbations)
  call fd_read_mode(parameters, fd_directional)
  call fd_read_probe_index(probe_dof, probe_dof_set)

  ! -------------------------------------------------------------------------- !
  ! Initialization of the components

  call user_setup(sim%neko_case%user)
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
  ! weighted-sum total the FD check gates on hides this). Ported verbatim
  ! from the unit-tier driver (tests/unit/sensitivity/problem_tester.f90),
  ! which grew it for the scalar_mixing case before that case moved into
  ! this tier. Cheap and harmless for every existing single-objective
  ! case; exists so a multi-objective case can be inspected from the log
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

     i_local = maxloc(abs(target_sensitivities(it)%x), dim=1)
     local_abs_max = abs(target_sensitivities(it)%x(i_local))
     global_abs_max = glmax(abs(target_sensitivities(it)%x), &
          target_sensitivities(it)%size()) ! DEVICE?

     ! rank-based tie-breaker: only those within eps of the global max get
     ! key=1
     if (abscmp(local_abs_max, global_abs_max)) then
        key_local = 1.0_rp + 1.0e-12_rp*real(pe_rank, rp)
     else
        key_local = 0.0_rp
     end if

     ! reduce the key to pick a single owner
     key_arr(1) = key_local
     key_global = glmax(key_arr, 1)
     if (abscmp(key_local, key_global)) then
        i_max = i_local
     else
        i_max = -1 ! to indicate that this proc doesn't participate
     end if

     ! Override the argmax probe with an explicitly named design dof, if one
     ! was requested. This is what makes two runs comparable: the argmax
     ! above moves whenever the sensitivity field moves, so two runs of the
     ! same case can otherwise silently differentiate with respect to two
     ! different design variables. The global design index to use is printed
     ! by every run (the 'FD probe' line), so pinning one run to another is a
     ! copy-paste.
     if (probe_dof_set) then
        if (fd_directional) then
           if (pe_rank .eq. 0 .and. it .eq. 1) then
              write(*, '(A)') ' FD probe: WARNING -- ' // &
                   'NEKO_TOP_FD_PROBE_INDEX is set but the mode is ' // &
                   'directional, which perturbs'
              write(*, '(A)') ' FD probe: WARNING -- every design dof at ' &
                   // 'once. The requested dof is ignored.'
           end if
        else
           call fd_resolve_probe_index(des%size(), probe_dof, i_max)
        end if
     end if

     if (fd_directional) then
        ! Perturb the whole design along s = g/||g|| and compare against
        ! <g, s>. Validates the entire gradient field in one sweep rather
        ! than one argmax-selected component of it.
        call compute_sensitivity_directional(prob, sim, des, &
             target_sensitivities(it), perturbations, tolerance, &
             trim(parameter_file), targets(it)%is_objective, &
             sim%fluid%gs_Xh, use_central, strict_options, &
             targets(it)%constraint_index, trim(targets(it)%suffix))
     else
        call compute_sensitivity(prob, sim, des, target_sensitivities(it), &
             i_max, perturbations, tolerance, trim(parameter_file), &
             targets(it)%is_objective, sim%fluid%gs_Xh, use_central, &
             strict_options, targets(it)%constraint_index, &
             trim(targets(it)%suffix))
     end if
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
  ! instead; this lane's CMakeLists turns it into a reported SKIP through
  ! `SKIP_REGULAR_EXPRESSION`, its `SKIP_RETURN_CODE` being already spoken
  ! for by the opt-in gate.
  if (fd_assertion_skipped()) stop FD_SKIP_EXIT_CODE

end program problem_tester
