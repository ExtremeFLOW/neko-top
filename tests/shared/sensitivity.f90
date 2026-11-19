!> Shared finite-difference sensitivity checker.
!!
!! This single module is compiled into both the unit driver
!! (`tests/unit/sensitivity/problem_tester.f90`) and the regression driver
!! (`tests/regression/sensitivity/problem_tester.f90`). It contains:
!!  * MPI-safe perturbation/reduction handling (`glsum`, `pe_rank` guards) so it
!!    is correct both single-rank (unit) and multi-rank (regression, `mpirun`).
!!  * A real tolerance assertion on the analytic-vs-finite-difference agreement,
!!    so both lanes actually gate on correctness rather than only logging a CSV.
!!
!! The assertion is deliberately made on the error at the *smallest*
!! perturbation. For a one-sided forward difference the error at large
!! perturbations is dominated by truncation (\f$O(\epsilon)\f$) and must not be
!! asserted against; only as the perturbation shrinks does the estimate approach
!! the analytic sensitivity, exposing any genuine bias (missing coupling term,
!! wrong sign, wrong weighting) as a non-vanishing floor. Linear functionals
!! (e.g. the volume constraint) hit that floor at round-off; PDE-coupled
!! objectives hit a discretisation/steady-state floor set by the case.
!!
!! That assertion is a false green whenever the signed error crosses zero
!! inside the sweep -- the point nearest the crossing can be arbitrarily
!! small while the bias is not. Setting `optimization.fd_test_strict` swaps it
!! for the criterion in `fd_criterion`, which separates the bias from
!! truncation before asserting on it, distinguishes an inadequate sweep from a
!! wrong gradient, and writes its verdict to `FD_verdict_<case>.csv`. It is
!! off by default, so no existing case changes verdict; see
!! `fd_read_strict_options`.
!!
!! `report_sweep` additionally *reports* the minimum \f$|error|\f$ over the
!! sweep, and whether that minimum is interior to the sweep (a genuine
!! finite-difference floor) or sits at one of its ends (only a bound on one).
!! That is diagnostic output only: nothing is asserted on it. Which point of a
!! sweep a gating assertion ought to use is a separate question, deliberately
!! left open here rather than answered in passing.
!!
!! Both the perturbation list and the choice of one-sided or central differences
!! are configurable per case -- see `fd_read_perturbations` and
!! `fd_read_central_difference` -- and default to the historical behaviour.
!!
!! The sweep can differentiate along either of two directions, selected by
!! `fd_read_mode`:
!!  * `dof` (the default, historical) perturbs one design degree of freedom and
!!    compares against the assembled derivative at that dof. Which dof that is
!!    is now always reported, because the caller's usual choice -- the
!!    largest-|sensitivity| dof -- *moves* whenever the sensitivity field
!!    moves, so two runs of the same case can otherwise silently differentiate
!!    with respect to two different design variables;
!!  * `directional` perturbs the *whole* design along the normalised
!!    sensitivity direction and compares against the projection of the
!!    gradient onto it. That validates the entire gradient field at once rather
!!    than one component of it, and raises the finite-difference signal by
!!    roughly \f$\|g\|/|g_i|\f$, which buys back round-off headroom.
!!
!! Both go through one perturbation code path (`fd_run_sweep` and
!! `evaluate_perturbed`, which take a direction vector); the single-dof probe
!! is the one-hot special case of it, so the two cannot drift apart.
!!
!! **Inner product.** The design-space inner product used for both the
!! normalisation and the projection is the Euclidean one on the *assembled
!! nodal design coefficients*. It is not additionally weighted by the mass
!! matrix \f$B\f$, because the \f$B\f$ weighting is already inside the
!! sensitivities the harness is handed: the drivers call
!! `design%convert_to_directional_derivative`, which post-multiplies the
!! adjoint's \f$L^2\f$ Riesz representative \f$g_{L^2}\f$ by \f$B\f$, so what
!! arrives here is \f$g = B g_{L^2}\f$, the vector of partial derivatives with
!! respect to the nodal coefficients. Those are the variables the
!! finite difference actually perturbs, so
!! \f$ \langle g, s\rangle_2 = s^\mathsf{T} B g_{L^2}
!!     = \int_\Omega g_{L^2}\, s \,\mathrm{d}\Omega \f$:
!! the Euclidean pairing on coefficients *is* the \f$B\f$-weighted \f$L^2\f$
!! pairing of the underlying functions, expressed in the right variables.
!! Applying \f$B\f$ a second time would count it twice and is exactly the
!! inconsistency that would silently invalidate the test.
module sensitivity
  use simulation_m, only: simulation_t
  use design, only: design_t
  use utils, only: neko_error
  use num_types, only: rp, sp, i8
  use math, only: abscmp, NEKO_EPS, glsum, glmin, glmax, glimax
  use vector, only: vector_t
  use neko_config, only: NEKO_BCKND_DEVICE
  use problem, only : problem_t
  use device, only: device_memcpy, DEVICE_TO_HOST, HOST_TO_DEVICE
  use csv_file, only : csv_file_t
  use comm, only: pe_rank, NEKO_COMM
  use mpi_f08, only: MPI_Allreduce, MPI_Exscan, MPI_SUM, MPI_INTEGER8
  use gather_scatter, only : gs_t, GS_OP_ADD
  use json_module, only: json_file
  use json_utils, only: json_get, json_get_or_default, json_extract_item
  use fd_criterion, only: fd_strict_options_t, fd_verdict_t, fd_evaluate, &
       fd_verdict_print, fd_branch_name, fd_status_name, fd_c_hat_kind, &
       FD_STATUS_OK, FD_STATUS_FLOOR_EXCEEDED, FD_STATUS_INCONCLUSIVE, &
       FD_STATUS_UNREACHABLE
  implicit none

  ! Everything the drivers do not call. The module parameters belong here
  ! too: they are implementation constants of the harness, and the list named
  ! only one of them, silently leaving the rest part of the module's
  ! interface.
  private :: count_tokens, fd_lowercase, fd_sync_to_host, fd_project, &
       fd_step_headroom, fd_report_probe, fd_run_sweep, fd_global_offset, &
       fd_write_verdict, fd_sensitivity_floor, fd_sensitivity_scale, &
       fd_assert_verdict, fd_skipped_assertion, fd_env_truncated, &
       fd_name_targets, &
       fd_default_perturbations, fd_strict_perturbations, fd_design_lower, &
       fd_design_upper, fd_no_limit, fd_direction_floor

  !> One thing a finite-difference driver checks: either the weighted total
  !! of every objective, or one individual constraint.
  !!
  !! A run checks a *list* of these, so that a case can gate on its objective
  !! and its constraints in the same invocation rather than being restricted
  !! to one or the other.
  type, public :: fd_target_t
     !> True when the target is the weighted total of all objectives, false
     !! when it is the constraint at `constraint_index`.
     logical :: is_objective = .true.
     !> 1-based index of the constraint to check. Meaningless, and ignored,
     !! when `is_objective` is true.
     integer :: constraint_index = 1
     !> Suffix appended to this target's CSV file names. Deliberately empty
     !! whenever the run has a single target, so that a one-target run keeps
     !! writing exactly `FD_check_<case>.csv` and the reference data it is
     !! compared against keeps matching.
     character(len=32) :: suffix = ''
  end type fd_target_t

  !> Status returned by `get_environment_variable` when the value did not fit
  !! in the buffer supplied. Silently accepting a truncated list would run a
  !! different sweep from the one that was asked for.
  integer, parameter :: fd_env_truncated = -1

  !> The historical perturbation sweep: four points spanning three decades.
  !! Used whenever a case supplies no sweep of its own, so that every existing
  !! case file and the CI lane behave exactly as before.
  real(kind=rp), parameter :: fd_default_perturbations(4) = [ &
       1e-1_rp, 1e-2_rp, 1e-3_rp, 1e-4_rp]

  !> The sweep the strict criterion defaults to: nine points, geometric with
  !! ratio \f$\sqrt{10}\f$, from 1e-1 down to 1e-5. Geometric because the
  !! order estimate is a statement about consecutive *ratios*; nine points so
  !! that a truncation run and a plateau both have room to appear; and one
  !! decade deeper than the historical sweep because the four-point sweep
  !! stops before the round-off upturn on every case measured.
  real(kind=rp), parameter :: fd_strict_perturbations(9) = [ &
       1.0e-1_rp, 3.1622776601683794e-2_rp, 1.0e-2_rp, &
       3.1622776601683794e-3_rp, 1.0e-3_rp, 3.1622776601683794e-4_rp, &
       1.0e-4_rp, 3.1622776601683794e-5_rp, 1.0e-5_rp]

  !> Fraction of the largest sensitivity in the field below which an
  !! individual analytic sensitivity is treated as degenerate. The historical
  !! `max(|b|, NEKO_EPS)` guard is far too permissive: a sensitivity of 1e-12
  !! passes it and the quotient is then reported as a relative error, which
  !! it is not. Anything below this floor is asserted on as an absolute
  !! difference instead.
  real(kind=rp), parameter :: fd_sensitivity_floor = 1e-10_rp

  !> Bounds assumed for a design variable. These match the clamping heuristic
  !! the one-sided sweep has always used (perturb downwards once the design
  !! value is past the midpoint) and are what the central-difference sweep
  !! clamps its symmetric half-step against.
  real(kind=rp), parameter :: fd_design_lower = 0.0_rp
  real(kind=rp), parameter :: fd_design_upper = 1.0_rp

  !> Fraction of the largest entry of the perturbation direction below which
  !! an entry is ignored when the step headroom is measured. A normalised
  !! gradient can carry entries many orders below its largest, down to
  !! denormals, and the headroom at such an entry is `(1 - x)/|s_j|`, which
  !! overflows -- and, under the Testing build's `-ffpe-trap=overflow`,
  !! SIGFPEs -- long before the entry constrains anything: at this floor the
  !! entry moves the design by at most 1e-12 of the step, which is below the
  !! round-off of the design value it is added to.
  real(kind=rp), parameter :: fd_direction_floor = 1e-12_rp

  !> Sentinel step limit used where a design coordinate constrains nothing --
  !! larger than any perturbation a sweep could sensibly request, so `min`
  !! against it is a no-op.
  real(kind=rp), parameter :: fd_no_limit = huge(1.0_rp)

  !> Exit code a driver returns when the finite-difference assertion was
  !! *skipped* rather than made. Deliberately not 77, which the regression
  !! lane's opt-in gate already uses, so that the two reasons a sweep did not
  !! run stay distinguishable in the CTest output. Wired to `SKIP_RETURN_CODE`
  !! in `tests/unit/sensitivity/CMakeLists.txt`.
  integer, parameter :: FD_SKIP_EXIT_CODE = 78

  !> True once a sweep has declined to make its assertion. A skipped
  !! assertion that still exits zero is reported by CTest as a pass, which is
  !! the one outcome a gradient check must never produce, so the drivers read
  !! this and exit with `FD_SKIP_EXIT_CODE` instead.
  logical :: fd_skipped_assertion = .false.

  interface compute_sensitivity
     module procedure compute_sensitivity_list, &
          compute_sensitivity_i
  end interface compute_sensitivity

contains

  !> Read the finite-difference perturbation sweep to use for a case.
  !!
  !! Sources, highest precedence first:
  !!  1. the environment variable `NEKO_TOP_FD_PERTURBATIONS`, a comma- or
  !!     space-separated list such as `1e-1,5e-2,1e-2`. This exists because a
  !!     sweep is a property of the *investigation*, not of the case: it lets
  !!     one wide sweep be driven across every case file without editing any
  !!     of them;
  !!  2. the case-file key `optimization.fd_test_perturbations`, a JSON array
  !!     of numbers, sitting beside the existing `fd_test_tolerance`;
  !!  3. `fd_strict_perturbations` when the strict criterion is in use, and
  !!     `fd_default_perturbations` -- the historical four-point sweep --
  !!     otherwise.
  !!
  !! Entries are perturbation *magnitudes* and must be strictly positive; the
  !! sign is chosen by the sweep itself so that the perturbed design stays
  !! within its bounds. The list need not be ordered or evenly spaced.
  !!
  !! @param params The case file to read the optional key from.
  !! @param strict True when the strict criterion is in use, which defaults
  !!        the sweep to `fd_strict_perturbations` rather than the historical
  !!        four points. An explicit sweep, from either source, still wins.
  !! @param perturbations Allocated and filled with the sweep to use.
  subroutine fd_read_perturbations(params, strict, perturbations)
    type(json_file), intent(inout) :: params
    logical, intent(in) :: strict
    real(kind=rp), allocatable, intent(out) :: perturbations(:)

    character(len=*), parameter :: env_name = 'NEKO_TOP_FD_PERTURBATIONS'
    !> Characters accepted between values. A blank is included so that
    !! replacing one below is a harmless no-op.
    character(len=*), parameter :: separators = ', ;'
    character(len=1024) :: env_value
    integer :: env_status, n_values, ios, ic

    call get_environment_variable(env_name, env_value, status = env_status)
    if (env_status .eq. fd_env_truncated) then
       call neko_error(env_name // ' is longer than the buffer available ' // &
            'to read it; shorten the sweep or use the case file')
    end if

    if (env_status .eq. 0 .and. len_trim(env_value) .gt. 0) then
       ! Accept commas and semicolons as separators by turning them into the
       ! blanks that list-directed input already understands.
       do ic = 1, len_trim(env_value)
          if (scan(env_value(ic:ic), separators) .gt. 0) then
             env_value(ic:ic) = ' '
          end if
       end do

       n_values = count_tokens(env_value)
       if (n_values .lt. 1) then
          call neko_error(env_name // ' is set but contains no values')
       end if

       allocate(perturbations(n_values))
       ! List-directed input stops quietly at a '/' and expands repeat
       ! counts such as 2*1e-2, both with iostat still zero, leaving any
       ! trailing elements unassigned. Pre-fill with a negative sentinel so
       ! every element the read skips fails the positivity check below
       ! instead of carrying whatever the allocation happened to hold.
       perturbations = -1.0_rp
       read(env_value, *, iostat = ios) perturbations
       if (ios .ne. 0) then
          call neko_error(env_name // ' could not be parsed as a list of ' // &
               'numbers')
       end if

    else if (params%valid_path('optimization.fd_test_perturbations')) then
       call json_get(params, 'optimization.fd_test_perturbations', &
            perturbations)
       ! json-fortran deallocates its result on an exception and Neko's
       ! json_get never checks failed(), so a key of the wrong type comes
       ! back unallocated rather than raising an error.
       if (.not. allocated(perturbations)) then
          call neko_error('optimization.fd_test_perturbations could not ' // &
               'be read as an array of numbers')
       end if
       if (size(perturbations) .lt. 1) then
          call neko_error('optimization.fd_test_perturbations is empty')
       end if

    else if (strict) then
       allocate(perturbations(size(fd_strict_perturbations)))
       perturbations = fd_strict_perturbations

    else
       allocate(perturbations(size(fd_default_perturbations)))
       perturbations = fd_default_perturbations
    end if

    if (any(perturbations .le. 0.0_rp)) then
       call neko_error('Finite-difference perturbations must be strictly ' // &
            'positive magnitudes')
    end if

  end subroutine fd_read_perturbations

  !> Read whether the sweep uses a central rather than a one-sided difference.
  !!
  !! Sources, highest precedence first: the environment variable
  !! `NEKO_TOP_FD_CENTRAL` (true: `1`, `true`, `yes`, `on`; false: `0`,
  !! `false`, `no`, `off`; all case-insensitive, anything else is an error
  !! rather than a silent false that would also silently override the case
  !! file), then the case-file key `optimization.fd_test_central_difference`,
  !! then `.false.` -- the historical one-sided forward difference.
  !!
  !! A central difference has \f$O(\epsilon^2)\f$ truncation error instead of
  !! \f$O(\epsilon)\f$, so it reaches the finite-difference floor at a far
  !! larger perturbation and separates floor from truncation much more
  !! cleanly. It costs two forward solves per perturbation instead of one.
  !!
  !! @param params The case file to read the optional key from.
  !! @param central True to use a central difference.
  subroutine fd_read_central_difference(params, central)
    type(json_file), intent(inout) :: params
    logical, intent(out) :: central

    character(len=*), parameter :: env_name = 'NEKO_TOP_FD_CENTRAL'
    character(len=32) :: env_value
    integer :: env_status

    call get_environment_variable(env_name, env_value, status = env_status)
    if (env_status .eq. fd_env_truncated) then
       call neko_error(env_name // ' is longer than the buffer available ' // &
            'to read it')
    end if

    if (env_status .eq. 0 .and. len_trim(env_value) .gt. 0) then
       call fd_lowercase(env_value)
       ! Strict: a value that is not recognisably true or false must not
       ! quietly mean false -- it would also quietly override the case file.
       select case (trim(adjustl(env_value)))
       case ('1', 'true', 'yes', 'on')
          central = .true.
       case ('0', 'false', 'no', 'off')
          central = .false.
       case default
          call neko_error(env_name // ' must be one of 1/true/yes/on or ' // &
               '0/false/no/off, not "' // trim(adjustl(env_value)) // '"')
       end select
    else if (params%valid_path('optimization.fd_test_central_difference')) then
       call json_get(params, 'optimization.fd_test_central_difference', &
            central)
    else
       central = .false.
    end if

  end subroutine fd_read_central_difference

  !> Read the strict finite-difference criterion's settings for a case.
  !!
  !! All four keys are optional and every default reproduces today's
  !! behaviour, `optimization.fd_test_strict` most of all: with it false the
  !! harness keeps asserting on the error at the smallest perturbation, so
  !! adding this criterion changes no existing case's verdict.
  !!
  !!  * `optimization.fd_test_strict` -- apply the strict criterion.
  !!  * `optimization.fd_test_order` -- the truncation order to expect.
  !!    Defaults to 1 for a one-sided difference and 2 for a central one,
  !!    which is why `central` has to be read first.
  !!  * `optimization.fd_test_order_tolerance` -- half-width of the accepted
  !!    band around that order.
  !!  * `optimization.fd_test_plateau_fraction` -- how flat, as a fraction of
  !!    the tolerance, a window must be to bound the bias.
  !!
  !! The last two are the criterion's calibration constants and are exposed
  !! as case-file keys precisely so that re-tuning them does not need a
  !! recompile.
  !!
  !! @param params The case file to read the optional keys from.
  !! @param central True if the sweep uses a central difference.
  !! @param opts The settings to use.
  subroutine fd_read_strict_options(params, central, opts)
    type(json_file), intent(inout) :: params
    logical, intent(in) :: central
    type(fd_strict_options_t), intent(out) :: opts

    type(fd_strict_options_t) :: defaults
    real(kind=rp) :: default_order

    call json_get_or_default(params, 'optimization.fd_test_strict', &
         opts%enabled, defaults%enabled)

    default_order = defaults%p_expected
    if (central) default_order = 2.0_rp
    call json_get_or_default(params, 'optimization.fd_test_order', &
         opts%p_expected, default_order)
    call json_get_or_default(params, &
         'optimization.fd_test_order_tolerance', opts%order_tolerance, &
         defaults%order_tolerance)
    call json_get_or_default(params, &
         'optimization.fd_test_plateau_fraction', opts%plateau_fraction, &
         defaults%plateau_fraction)

    ! A non-positive setting here does not merely misbehave: it makes the
    ! criterion inert -- no ratio can match a band of zero width, and no
    ! window can be flatter than a zero spread -- while still reporting a
    ! verdict, so it must be refused rather than accepted.
    if (opts%p_expected .le. 0.0_rp) then
       call neko_error('optimization.fd_test_order must be positive')
    end if
    if (opts%order_tolerance .le. 0.0_rp) then
       call neko_error('optimization.fd_test_order_tolerance must be ' // &
            'positive')
    end if
    if (opts%plateau_fraction .le. 0.0_rp) then
       call neko_error('optimization.fd_test_plateau_fraction must be ' // &
            'positive')
    end if

    ! A band at least as wide as the order it is centred on is meaningless --
    ! it accepts every measured order down to zero as "the expected one" --
    ! and it is not merely useless: it admits p_hat ~ 0, where the Richardson
    ! extrapolation is singular. `fd_solve_order` bisects [0, 6] sixty times,
    ! so its smallest non-zero root is 3/2**60 = 2.6e-18, at which
    ! `ratio**p - 1.0_rp` is exactly 0.0 in double precision. Removing the
    ! invalid state here is preferable to guarding the division downstream.
    if (opts%order_tolerance .ge. opts%p_expected) then
       call neko_error('optimization.fd_test_order_tolerance must be ' // &
            'smaller than optimization.fd_test_order: a band at least as ' // &
            'wide as the expected order accepts a measured order of zero, ' &
            // 'which is not a truncation branch at all and makes the ' // &
            'Richardson extrapolation singular.')
    end if

  end subroutine fd_read_strict_options

  !> Read whether the sweep perturbs a single design degree of freedom or the
  !! whole design along the sensitivity direction.
  !!
  !! Sources, highest precedence first: the environment variable
  !! `NEKO_TOP_FD_MODE` (`dof` or `directional`, case-insensitive), then the
  !! case-file key `optimization.fd_test_mode`, then `dof` -- the historical
  !! single-degree-of-freedom probe. Anything else is an error rather than a
  !! silent fall back to `dof`, which would also silently override the case
  !! file, matching how `fd_read_central_difference` treats its value.
  !!
  !! `directional` perturbs every design dof at once along
  !! \f$ s = g/\|g\|_2 \f$ and compares the finite difference against
  !! \f$ \langle g, s \rangle = \|g\|_2 \f$, so it tests the whole gradient
  !! field in one sweep instead of one lottery-selected component of it, and
  !! its finite-difference signal is larger by roughly \f$\|g\|/|g_i|\f$.
  !!
  !! @param params The case file to read the optional key from.
  !! @param directional True to perturb along the sensitivity direction.
  subroutine fd_read_mode(params, directional)
    type(json_file), intent(inout) :: params
    logical, intent(out) :: directional

    character(len=*), parameter :: env_name = 'NEKO_TOP_FD_MODE'
    character(len=32) :: env_value
    character(len=:), allocatable :: mode_string
    integer :: env_status

    call get_environment_variable(env_name, env_value, status = env_status)
    if (env_status .eq. fd_env_truncated) then
       call neko_error(env_name // ' is longer than the buffer available ' // &
            'to read it')
    end if

    if (env_status .eq. 0 .and. len_trim(env_value) .gt. 0) then
       mode_string = trim(adjustl(env_value))
    else if (params%valid_path('optimization.fd_test_mode')) then
       call json_get(params, 'optimization.fd_test_mode', mode_string)
       ! json-fortran deallocates its result on an exception and Neko's
       ! json_get never checks failed(), so a key of the wrong type comes
       ! back unallocated rather than raising an error.
       if (.not. allocated(mode_string)) then
          call neko_error('optimization.fd_test_mode could not be read as ' // &
               'a string')
       end if
    else
       mode_string = 'dof'
    end if

    call fd_lowercase(mode_string)

    select case (trim(mode_string))
    case ('dof')
       directional = .false.
    case ('directional')
       directional = .true.
    case default
       call neko_error('The finite-difference mode must be "dof" or ' // &
            '"directional", not "' // trim(mode_string) // '"')
    end select

    deallocate(mode_string)

  end subroutine fd_read_mode

  !> Read the fixed design degree of freedom to probe, if one was requested.
  !!
  !! `NEKO_TOP_FD_PROBE_INDEX` names a **global design index**: a 1-based index
  !! into the globally concatenated design vector, i.e. the owning rank's
  !! offset plus the local index. It is what `fd_report_probe` prints for
  !! whichever dof a run actually probed, so pinning it is a copy-paste out of
  !! a previous log -- and pinning it is what makes two runs comparable, since
  !! the default probe is the largest-|sensitivity| dof and that argmax moves
  !! whenever the sensitivity field moves. Without it two runs of the same case
  !! can silently differentiate with respect to two different design variables.
  !!
  !! The index is unique and reproducible for a given mesh **and rank count**,
  !! which is the only scope in which two runs are comparable anyway. Run the
  !! same number on a different number of ranks and it names a different dof;
  !! the coordinates printed beside it make that visible immediately, so check
  !! them rather than assuming.
  !!
  !! Only meaningful in `dof` mode; the directional sweep probes no single dof.
  !!
  !! @param probe_index The global design index requested, unset if `is_set` is
  !!        false.
  !! @param is_set True if the environment variable was set and non-blank.
  subroutine fd_read_probe_index(probe_index, is_set)
    integer(kind=i8), intent(out) :: probe_index
    logical, intent(out) :: is_set

    character(len=*), parameter :: env_name = 'NEKO_TOP_FD_PROBE_INDEX'
    character(len=64) :: env_value
    integer :: env_status, ios

    probe_index = -1_i8
    is_set = .false.

    call get_environment_variable(env_name, env_value, status = env_status)
    if (env_status .eq. fd_env_truncated) then
       call neko_error(env_name // ' is longer than the buffer available ' // &
            'to read it')
    end if
    if (env_status .ne. 0 .or. len_trim(env_value) .eq. 0) return

    read(env_value, *, iostat = ios) probe_index
    if (ios .ne. 0) then
       call neko_error(env_name // ' could not be parsed as an integer ' // &
            'global design index')
    end if
    if (probe_index .lt. 1_i8) then
       call neko_error(env_name // ' must be a positive global design index')
    end if
    is_set = .true.

  end subroutine fd_read_probe_index

  !> Read the list of finite-difference targets to check for a case.
  !!
  !! Optional case-file key `optimization.fd_test_targets`, an array of
  !! objects:
  !!
  !! ```json
  !! "fd_test_targets": [
  !!     { "kind": "objective" },
  !!     { "kind": "constraint", "index": 1 },
  !!     { "kind": "constraint", "name": "volume" }
  !! ]
  !! ```
  !!
  !! A constraint is named by `index` (1-based, in case-file order) or by
  !! `name` (the constraint's own name), exactly one of the two.
  !!
  !! `kind: "objective"` means the **weighted total of every objective**, and
  !! only that; supplying an `index` beside it is an error rather than a
  !! per-objective check. See the error text raised below for why that is a
  !! structural property of the adjoint rather than a missing feature.
  !!
  !! **When the key is absent the historical dispatch is reproduced exactly**:
  !! the objective total if the case has no constraints, and constraint 1
  !! otherwise. Every existing case file therefore behaves as before.
  !!
  !! @param params The case file to read the optional key from.
  !! @param problem The problem, used to bound the constraint indices and to
  !!        resolve constraint names.
  !! @param targets The targets to check, in the order they were listed.
  subroutine fd_read_targets(params, problem, targets)
    type(json_file), intent(inout) :: params
    class(problem_t), intent(in) :: problem
    type(fd_target_t), allocatable, intent(out) :: targets(:)

    character(len=*), parameter :: path = 'optimization.fd_test_targets'
    type(json_file) :: target_json
    character(len=:), allocatable :: kind, name
    character(len=32) :: index_str, count_str
    integer :: n_targets, n_constraints, i, j, constraint_index
    logical :: has_index, has_name, found

    n_constraints = problem%get_n_constraints()

    if (.not. params%valid_path(path)) then
       allocate(targets(1))
       if (n_constraints .eq. 0) then
          targets(1)%is_objective = .true.
       else
          targets(1)%is_objective = .false.
          targets(1)%constraint_index = 1

          ! Not an error, because the historical dispatch is what is being
          ! reproduced here, but not silent either: a case with several
          ! constraints and no target list gets exactly one of them checked,
          ! and that is worth saying out loud.
          if (n_constraints .gt. 1 .and. pe_rank .eq. 0) then
             write(*, '(A,I0,A)') ' FD targets: WARNING -- this case has ', &
                  n_constraints, ' constraints but no'
             write(*, '(A)') ' FD targets: WARNING -- ' // path // &
                  ', so only constraint 1 is checked. List the others'
             write(*, '(A)') ' FD targets: WARNING -- explicitly to check' &
                  // ' them too.'
          end if
       end if
       call fd_name_targets(targets)
       return
    end if

    call params%info(path, n_children = n_targets)
    if (n_targets .lt. 1) then
       call neko_error(path // ' is present but lists no targets. Remove ' &
            // 'the key to take the default target, or list at least one.')
    end if

    allocate(targets(n_targets))
    do i = 1, n_targets
       call json_extract_item(params, path, i, target_json)
       call json_get(target_json, 'kind', kind)
       call fd_lowercase(kind)

       has_index = target_json%valid_path('index')
       has_name = target_json%valid_path('name')

       select case (trim(kind))
       case ('objective')

          ! The one place this restriction is stated. An objective gradient
          ! does not decompose per objective, so an index here would assert
          ! something meaningless rather than something narrower.
          if (has_index .or. has_name) then
             call neko_error('fd_test_targets: an "objective" target ' // &
                  'takes neither "index" nor "name" -- it is always the ' // &
                  'weighted total of every objective. The objective ' // &
                  'gradient does not decompose per objective: ' // &
                  'problem_read_objectives injects an augmented ' // &
                  'Lagrangian objective into every case, and that single ' // &
                  'object carries all of the PDE-mediated d/dchi, ' // &
                  'computed from one adjoint solve forced by the sum of ' // &
                  'every objective''s forcing. The individual objectives ' // &
                  '(viscous_dissipation, scalar_mixing) have empty ' // &
                  'update_sensitivity bodies and contribute only that ' // &
                  'forcing. Objective values decompose per objective; ' // &
                  'objective gradients do not, so there is no "objective ' // &
                  'i''s share of u_adj" to check against.')
          end if

          targets(i)%is_objective = .true.

       case ('constraint')

          if (has_index .and. has_name) then
             call neko_error('fd_test_targets: a "constraint" target ' // &
                  'takes either "index" or "name", not both.')
          else if (.not. (has_index .or. has_name)) then
             call neko_error('fd_test_targets: a "constraint" target ' // &
                  'needs an "index" (1-based, in case-file order) or a ' // &
                  '"name".')
          end if

          if (has_index) then
             call json_get(target_json, 'index', constraint_index)
          else
             call json_get(target_json, 'name', name)
             found = .false.
             constraint_index = 0
             do j = 1, n_constraints
                if (trim(problem%get_constraint_name(j)) .eq. trim(name)) &
                     then
                   if (found) then
                      call neko_error('fd_test_targets: the constraint ' // &
                           'name "' // trim(name) // '" is carried by ' // &
                           'more than one constraint in this case; ' // &
                           'select it by "index" instead.')
                   end if
                   found = .true.
                   constraint_index = j
                end if
             end do
             if (.not. found) then
                call neko_error('fd_test_targets: no constraint in this ' // &
                     'case is named "' // trim(name) // '".')
             end if
          end if

          if (constraint_index .lt. 1 .or. &
               constraint_index .gt. n_constraints) then
             write(index_str, '(I0)') constraint_index
             write(count_str, '(I0)') n_constraints
             call neko_error('fd_test_targets: constraint index ' // &
                  trim(index_str) // ' is out of range; this case has ' // &
                  trim(count_str) // ' constraints.')
          end if

          targets(i)%is_objective = .false.
          targets(i)%constraint_index = constraint_index

       case default
          call neko_error('fd_test_targets: unknown target kind "' // &
               trim(kind) // '". The kinds are "objective" and ' // &
               '"constraint".')
       end select
    end do

    call fd_name_targets(targets)

  end subroutine fd_read_targets

  !> Give each target the suffix its CSV files carry.
  !!
  !! A single-target run gets an **empty** suffix, so that it keeps writing
  !! `FD_check_<case>.csv` and `FD_verdict_<case>.csv` under exactly the
  !! names it has always used -- the reference data those are compared
  !! against is matched by name. Only a multi-target run, which would
  !! otherwise interleave several sweeps into one file, suffixes them.
  !!
  !! @param targets The targets to name, updated in place.
  subroutine fd_name_targets(targets)
    type(fd_target_t), intent(inout) :: targets(:)

    character(len=16) :: index_str
    integer :: i

    if (size(targets) .eq. 1) then
       targets(1)%suffix = ''
       return
    end if

    do i = 1, size(targets)
       if (targets(i)%is_objective) then
          targets(i)%suffix = '__objective'
       else
          write(index_str, '(I0)') targets(i)%constraint_index
          targets(i)%suffix = '__constraint_' // trim(index_str)
       end if
    end do

  end subroutine fd_name_targets

  !> Offset of this rank's design degrees of freedom within the globally
  !! concatenated design vector.
  !!
  !! Collective: every rank must call this.
  !!
  !! @param n Number of design degrees of freedom held locally.
  !! @return The number held by all lower-numbered ranks together.
  function fd_global_offset(n) result(offset)
    integer, intent(in) :: n
    integer(kind=i8) :: offset

    integer(kind=i8) :: local(1), scanned(1)
    integer :: ierr

    local(1) = int(n, i8)
    scanned(1) = 0_i8
    call MPI_Exscan(local, scanned, 1, MPI_INTEGER8, MPI_SUM, NEKO_COMM, ierr)
    ! The result on rank 0 is undefined by the standard, so do not read it.
    if (pe_rank .eq. 0) scanned(1) = 0_i8
    offset = scanned(1)

  end function fd_global_offset

  !> Resolve a global design index to a local index on exactly one rank.
  !!
  !! The per-rank index ranges are disjoint by construction, so at most one
  !! rank can claim the index; a count of claims is reduced so that an index
  !! past the end of the global design vector is an error rather than a silent
  !! fall back to whatever the argmax happened to pick.
  !!
  !! Collective: every rank must call this.
  !!
  !! @param n Number of design degrees of freedom held locally.
  !! @param probe_index The requested global design index.
  !! @param i Local index on the owning rank, -1 on every other rank.
  subroutine fd_resolve_probe_index(n, probe_index, i)
    integer, intent(in) :: n
    integer(kind=i8), intent(in) :: probe_index
    integer, intent(out) :: i

    integer(kind=i8) :: offset, local, claims(1), total_claims(1)
    integer :: ierr
    character(len=32) :: index_str

    offset = fd_global_offset(n)
    local = probe_index - offset

    if (local .ge. 1_i8 .and. local .le. int(n, i8)) then
       i = int(local)
       claims(1) = 1_i8
    else
       i = -1
       claims(1) = 0_i8
    end if

    call MPI_Allreduce(claims, total_claims, 1, MPI_INTEGER8, MPI_SUM, &
         NEKO_COMM, ierr)

    if (total_claims(1) .ne. 1_i8) then
       write(index_str, '(I0)') probe_index
       call neko_error('The requested probe index (' // trim(index_str) // &
            ') is outside the global design vector; take the number from ' // &
            'the "FD probe" line of a previous run of the same case, mesh ' // &
            'and rank count')
    end if

  end subroutine fd_resolve_probe_index

  !> Fold a string to lower case, in place.
  !! @param string The string to fold.
  subroutine fd_lowercase(string)
    character(len=*), intent(inout) :: string
    integer :: ic

    do ic = 1, len(string)
       if (string(ic:ic) .ge. 'A' .and. string(ic:ic) .le. 'Z') then
          string(ic:ic) = achar(iachar(string(ic:ic)) + 32)
       end if
    end do

  end subroutine fd_lowercase

  !> Copy a vector's device buffer back to the host, if there is one.
  !!
  !! Every reduction and loop in this module works on the host array, so a
  !! device build has to synchronise first.
  !!
  !! @param vec The vector to synchronise.
  subroutine fd_sync_to_host(vec)
    type(vector_t), intent(inout) :: vec

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(vec%x, vec%x_d, vec%size(), DEVICE_TO_HOST, .true.)
    end if

  end subroutine fd_sync_to_host

  !> Count the blank-separated tokens in a string.
  !! @param string The string to inspect.
  !! @return The number of tokens found.
  function count_tokens(string) result(n_tokens)
    character(len=*), intent(in) :: string
    integer :: n_tokens
    integer :: ic, n
    logical :: in_token

    n = len_trim(string)
    n_tokens = 0
    in_token = .false.
    do ic = 1, n
       if (string(ic:ic) .eq. ' ' .or. string(ic:ic) .eq. achar(9)) then
          in_token = .false.
       else if (.not. in_token) then
          in_token = .true.
          n_tokens = n_tokens + 1
       end if
    end do

  end function count_tokens

  !> Pair the per-copy analytic sensitivities with a perturbation direction.
  !!
  !! `sensitivities` holds, at each element-local index, the derivative of the
  !! functional with respect to *that copy* of the design dof; the copies are
  !! not assembled. `direction` is single-valued (equal on every copy of a
  !! shared dof), so
  !! \f$ \sum_j g_j s_j = \sum_I s_I \sum_{j\in\mathrm{copies}(I)} g_j
  !!     = \sum_I s_I \, \partial F/\partial x_I \f$:
  !! the Euclidean inner product of the *assembled* gradient with the
  !! direction, obtained without ever dividing by a multiplicity. For a
  !! one-hot direction this reduces exactly to the assembled derivative at the
  !! single dof, which is what the single-dof sweep has always compared
  !! against.
  !!
  !! @param sensitivities Per-copy analytic derivatives.
  !! @param direction Single-valued perturbation direction.
  !! @param n Number of design degrees of freedom held locally.
  !! @return The globally reduced inner product.
  function fd_project(sensitivities, direction, n) result(projection)
    real(kind=rp), intent(in) :: sensitivities(:)
    real(kind=rp), intent(in) :: direction(:)
    integer, intent(in) :: n
    real(kind=rp) :: projection

    real(kind=rp) :: work_arr(1)
    integer :: j

    work_arr(1) = 0.0_rp
    do j = 1, n
       if (abs(direction(j)) .gt. 0.0_rp) then
          work_arr(1) = work_arr(1) + sensitivities(j) * direction(j)
       end if
    end do
    projection = glsum(work_arr, 1)

  end function fd_project

  !> The floor below which an individual analytic sensitivity is degenerate.
  !!
  !! Relative to the *field*, not to itself: what makes a relative error
  !! meaningful is that the quantity it is divided by is a real number on the
  !! scale of the problem. The historical `max(|b|, NEKO_EPS)` guard only
  !! rules out an exact zero, so a sensitivity of 1e-12 in a field whose
  !! largest entry is 1 passes it and the quotient is then reported, and
  !! asserted on, as though it were a relative error.
  !!
  !! Collective: every rank must call this, and every rank gets the same
  !! answer, so nothing downstream of it needs to communicate.
  !!
  !! @param target_sensitivities The analytic sensitivities, already
  !!        synchronised to the host.
  !! @param n Number of design degrees of freedom held locally.
  !! @return The degeneracy floor for this field.
  function fd_sensitivity_scale(target_sensitivities, n) result(eps_sens)
    type(vector_t), intent(in) :: target_sensitivities
    integer, intent(in) :: n
    real(kind=rp) :: eps_sens

    real(kind=rp) :: work_arr(1)

    work_arr(1) = 0.0_rp
    if (n .gt. 0) then
       work_arr(1) = maxval(abs(target_sensitivities%x(1:n)))
    end if
    eps_sens = max(NEKO_EPS, fd_sensitivity_floor * glmax(work_arr, 1))

  end function fd_sensitivity_scale

  !> Largest step magnitudes that keep the perturbed design inside its bounds.
  !!
  !! Returns the largest \f$t\f$ for which \f$x + t\,s\f$ (and, separately,
  !! \f$x - t\,s\f$) is within \f$[0,1]\f$ at *every* design coordinate. Both
  !! limits are globally reduced, so every rank agrees: a shared dof is held by
  !! several ranks and each must clamp identically or the design stops being
  !! single-valued. Coordinates whose direction entry is below
  !! `fd_direction_floor` times the largest entry of the direction barely
  !! move and place no limit; they are skipped rather than divided by, since
  !! dividing by a denormal overflows.
  !!
  !! For the one-hot direction of the single-dof sweep this reduces to
  !! \f$x_i\f$ of headroom downwards and \f$1-x_i\f$ upwards, exactly the
  !! clamp that sweep has always used.
  !!
  !! @param design_values The unperturbed design.
  !! @param direction The single-valued perturbation direction.
  !! @param n Number of design degrees of freedom held locally.
  !! @param limit_plus Largest step magnitude along `+direction`.
  !! @param limit_minus Largest step magnitude along `-direction`.
  subroutine fd_step_headroom(design_values, direction, n, limit_plus, &
       limit_minus)
    real(kind=rp), intent(in) :: design_values(:)
    real(kind=rp), intent(in) :: direction(:)
    integer, intent(in) :: n
    real(kind=rp), intent(out) :: limit_plus, limit_minus

    real(kind=rp) :: local_plus(1), local_minus(1), up, down
    real(kind=rp) :: local_scale(1), floor_value
    integer :: j

    ! The floor is relative to the whole direction, so it must be measured
    ! over the whole direction: the reduction is collective and sits outside
    ! any ownership test, and it also makes every rank skip the same entries
    ! of a shared dof.
    local_scale(1) = 0.0_rp
    do j = 1, n
       local_scale(1) = max(local_scale(1), abs(direction(j)))
    end do
    floor_value = fd_direction_floor * glmax(local_scale, 1)

    local_plus(1) = fd_no_limit
    local_minus(1) = fd_no_limit
    do j = 1, n
       if (abs(direction(j)) .le. floor_value) cycle
       up = (fd_design_upper - design_values(j)) / abs(direction(j))
       down = (design_values(j) - fd_design_lower) / abs(direction(j))
       if (direction(j) .gt. 0.0_rp) then
          local_plus(1) = min(local_plus(1), up)
          local_minus(1) = min(local_minus(1), down)
       else
          local_plus(1) = min(local_plus(1), down)
          local_minus(1) = min(local_minus(1), up)
       end if
    end do

    limit_plus = glmin(local_plus, 1)
    limit_minus = glmin(local_minus, 1)

  end subroutine fd_step_headroom

  !> Report which design variable (or which direction) the sweep differentiates
  !! with respect to.
  !!
  !! Printed in **every** mode, unconditionally. The default single-dof probe
  !! is the largest-|sensitivity| dof, and that argmax *moves* whenever the
  !! sensitivity field moves: without this line two runs of the same case can
  !! differentiate with respect to two different design variables and nothing
  !! in the output says so. Measured on this suite, that alone spread the
  !! relative error of one case at one timestep from 0.42% to 23%.
  !!
  !! The number printed is the **global design index**: this rank's offset in
  !! the globally concatenated design vector plus the local index. That is
  !! exactly what `NEKO_TOP_FD_PROBE_INDEX` consumes, so a later run can be
  !! pinned to this run's dof by copying the number out of the log. It is
  !! reproducible for a given mesh and rank count, which is the only scope in
  !! which two runs compare anyway; the coordinates are printed beside it so
  !! that a mismatch is obvious rather than assumed away.
  !!
  !! Collective: every rank must call this.
  !!
  !! @param des The design, for the probed dof's coordinates.
  !! @param n Number of design degrees of freedom held locally.
  !! @param probe_index Local index of the probed dof, negative on ranks that
  !!        do not own it, and negative everywhere in directional mode.
  !! @param design_value The design value at the probed dof.
  !! @param directional True when the whole design is perturbed along a
  !!        direction rather than one dof being probed.
  subroutine fd_report_probe(des, n, probe_index, design_value, directional)
    class(design_t), intent(inout) :: des
    integer, intent(in) :: n
    integer, intent(in) :: probe_index
    real(kind=rp), intent(in) :: design_value
    logical, intent(in) :: directional

    integer(kind=i8) :: send(4), recv(4), offset
    real(kind=rp) :: coords(3), work_arr(1)
    integer :: ierr, k, owner, local_index
    integer(kind=i8) :: global_index

    if (directional) then
       if (pe_rank .eq. 0) then
          write(*, '(A)') ' FD probe: mode = directional -- every design ' // &
               'degree of freedom is perturbed together'
          write(*, '(A)') ' FD probe: along s = g/||g||_2, so no single dof' &
               // ' is probed and NEKO_TOP_FD_PROBE_INDEX'
          write(*, '(A)') ' FD probe: does not apply.'
       end if
       return
    end if

    ! Collective on every rank, so it must sit outside the ownership test.
    offset = fd_global_offset(n)

    ! Reduce the probed dof's identity off the owning rank so that rank 0 --
    ! which owns the log -- can report it. Reduced as integers rather than
    ! through `glsum` so that a large global index stays exact.
    send = 0_i8
    coords = 0.0_rp
    if (probe_index .gt. 0) then
       send(1) = offset + int(probe_index, i8)
       send(2) = int(pe_rank, i8)
       send(3) = int(probe_index, i8)
       send(4) = 1_i8
       coords(1) = des%x(probe_index)
       coords(2) = des%y(probe_index)
       coords(3) = des%z(probe_index)
    end if
    call MPI_Allreduce(send, recv, 4, MPI_INTEGER8, MPI_SUM, NEKO_COMM, ierr)
    do k = 1, 3
       work_arr(1) = coords(k)
       coords(k) = glsum(work_arr, 1)
    end do

    ! The sweep perturbs one dof and compares against the derivative summed
    ! over its copies, which is only the right pairing if exactly one rank
    ! nominated it. Two owners would perturb two design variables while
    ! comparing against a one-dof derivative, and report a spurious error.
    if (recv(4) .ne. 1_i8) then
       call neko_error('The finite-difference probe must be owned by ' // &
            'exactly one rank; the caller nominated it on none or on ' // &
            'several. Use the rank tie-break the drivers already apply.')
    end if

    global_index = recv(1)
    owner = int(recv(2))
    local_index = int(recv(3))

    if (pe_rank .eq. 0) then
       write(*, '(A,I0,A,I0,A,I0,A)') ' FD probe: mode = dof -- global ' &
            // 'design index ', global_index, ' (rank ', owner, &
            ', local index ', local_index, ')'
       write(*, '(A,F12.6,A,3F11.6,A)') ' FD probe: design value', &
            design_value, ', location [', coords(1), coords(2), coords(3), ']'
       write(*, '(A)') ' FD probe: reproduce this exact dof in another run ' &
            // 'by setting NEKO_TOP_FD_PROBE_INDEX to the number above.'
    end if

  end subroutine fd_report_probe

  !> Sweep a set of perturbations at one design degree of freedom and assert
  !! that the finite-difference estimate agrees with the analytic sensitivity.
  !!
  !! A thin wrapper over `fd_run_sweep`: it builds the one-hot perturbation
  !! direction for the requested dof and hands it to the single perturbation
  !! path the directional sweep also uses.
  !!
  !! @param problem The problem supplying the objective/constraint value.
  !! @param sim The simulation driving each forward solve.
  !! @param des The design being perturbed.
  !! @param target_sensitivities The analytic sensitivities to check against.
  !! @param i Local index of the dof to perturb, negative on ranks that do not
  !!          own it.
  !! @param perturbations The sweep of perturbation magnitudes, in any order.
  !! @param tolerance The largest acceptable relative error at the smallest
  !!        perturbation of the sweep.
  !! @param file_name The case file name, used to name the CSV log.
  !! @param is_objective True to test an objective, false a constraint.
  !! @param gs_h Gather-scatter handle for the design's dofmap. Used to keep a
  !!             perturbed *shared* design dof consistent across every element
  !!             that owns a copy of it -- see `evaluate_perturbed`.
  !! @param central_difference (Optional) True to use a central difference,
  !!        which halves the truncation order at the cost of a second forward
  !!        solve per perturbation. Defaults to the one-sided forward
  !!        difference used historically.
  !! @param strict_options (Optional) Settings of the strict criterion.
  !!        Defaults to disabled, i.e. the historical assertion.
  !! @param constraint_index (Optional) 1-based index of the constraint to
  !!        read when `is_objective` is false. Defaults to 1.
  !! @param name_suffix (Optional) Suffix for the CSV file names, so that
  !!        several targets checked in one run do not share a file. Defaults
  !!        to none, which keeps the historical `FD_check_<case>.csv`.
  subroutine compute_sensitivity_i(problem, sim, des, target_sensitivities, i, &
       perturbations, tolerance, file_name, is_objective, gs_h, &
       central_difference, strict_options, constraint_index, name_suffix)
    class(problem_t), intent(inout) :: problem
    type(simulation_t), intent(inout) :: sim
    class(design_t), intent(inout) :: des
    type(vector_t), intent(inout) :: target_sensitivities
    integer, intent(in) :: i
    real(kind=rp), intent(in) :: perturbations(:)
    real(kind=rp), intent(in) :: tolerance
    character(len=*), intent(in) :: file_name
    logical, intent(in) :: is_objective
    type(gs_t), intent(inout) :: gs_h
    logical, intent(in), optional :: central_difference
    type(fd_strict_options_t), intent(in), optional :: strict_options
    integer, intent(in), optional :: constraint_index
    character(len=*), intent(in), optional :: name_suffix

    real(kind=rp), allocatable :: direction(:)
    real(kind=rp) :: target_sensitivity_i, eps_sens
    type(fd_strict_options_t) :: strict
    character(len=64) :: suffix
    logical :: central
    integer :: n, i_constraint

    central = .false.
    if (present(central_difference)) central = central_difference
    if (present(strict_options)) strict = strict_options
    i_constraint = 1
    if (present(constraint_index)) i_constraint = constraint_index
    suffix = ''
    if (present(name_suffix)) suffix = name_suffix

    n = des%size()
    call fd_sync_to_host(target_sensitivities)
    eps_sens = fd_sensitivity_scale(target_sensitivities, n)

    ! Build the one-hot direction: set the owning entry to 1 and
    ! gather-scatter with GS_OP_ADD -- the same mechanism Neko uses to
    ! propagate values between elements -- so that it is 1 on *every*
    ! element-local copy of the dof, on every rank holding one. Perturbing only
    ! the owning index would leave the design multi-valued at that point, which
    ! is not a perturbation of any real design variable.
    allocate(direction(n))
    direction = 0.0_rp
    ! `.gt. 0`, not `.ge. 0`: the index is 1-based, so 0 is not an owned
    ! entry -- it is what `maxloc` returns on a zero-size local design, and
    ! `direction(0)` is a write past the start of the array.
    if (i .gt. 0) direction(i) = 1.0_rp
    call gs_h%op(direction, n, GS_OP_ADD)

    ! The target is the derivative summed over every copy of the dof, which is
    ! what moving all the copies together responds to -- see `fd_project`.
    ! Measured directly: at a dof shared by two elements, moving both copies
    ! doubles the finite-difference response, and comparing that against a
    ! single copy's derivative reports a spurious ~-95% error while comparing
    ! against the assembled sum reports the same ~2.5% as the single-copy
    ! formulation. For an element-interior dof the sum is over one entry.
    target_sensitivity_i = fd_project(target_sensitivities%x, direction, n)

    call fd_run_sweep(problem, sim, des, direction, target_sensitivity_i, i, &
         perturbations, tolerance, file_name, is_objective, central, &
         .false., eps_sens, strict, i_constraint, trim(suffix))

    deallocate(direction)

  end subroutine compute_sensitivity_i

  !> Sweep a set of perturbations along the normalised sensitivity direction
  !! and assert that the finite-difference estimate agrees with the projection
  !! of the analytic gradient onto it.
  !!
  !! This is the Taylor test: instead of one lottery-selected design variable
  !! it validates the *whole* gradient field at once, against
  !! \f$ \langle g, s\rangle \f$ with \f$ s = g/\|g\|_2 \f$. It also raises the
  !! finite-difference signal by roughly \f$\|g\|/|g_i|\f$ relative to the
  !! single-dof probe, which buys back round-off headroom.
  !!
  !! **Inner product** (see also the module header): Euclidean on the assembled
  !! nodal design coefficients, for both the normalisation and the projection.
  !! The drivers have already applied `convert_to_directional_derivative`, so
  !! the sensitivities handed in are \f$g = B g_{L^2}\f$ -- derivatives with
  !! respect to the nodal coefficients, which are the variables perturbed here.
  !! Weighting again by \f$B\f$ would count the mass matrix twice.
  !!
  !! The direction is built from the **assembled** gradient, so it is
  !! single-valued on shared dofs; that is what makes the projection exact and
  !! makes the perturbation a perturbation of real design variables. The
  !! normalisation is itself a free scaling -- what is load-bearing is that
  !! \f$\langle g,s\rangle\f$ uses exactly the same \f$s\f$ that is added to
  !! the design, which it does by construction.
  !!
  !! @param problem The problem supplying the objective/constraint value.
  !! @param sim The simulation driving each forward solve.
  !! @param des The design being perturbed.
  !! @param target_sensitivities The analytic sensitivities to check against.
  !! @param perturbations The sweep of perturbation magnitudes, in any order.
  !! @param tolerance The largest acceptable relative error at the smallest
  !!        perturbation of the sweep.
  !! @param file_name The case file name, used to name the CSV log.
  !! @param is_objective True to test an objective, false a constraint.
  !! @param gs_h Gather-scatter handle for the design's dofmap, used both to
  !!             assemble the gradient and to count the copies of each dof.
  !! @param central_difference (Optional) True to use a central difference.
  !! @param strict_options (Optional) Settings of the strict criterion.
  !!        Defaults to disabled, i.e. the historical assertion.
  !! @param constraint_index (Optional) 1-based index of the constraint to
  !!        read when `is_objective` is false. Defaults to 1.
  !! @param name_suffix (Optional) Suffix for the CSV file names. Defaults to
  !!        none, which keeps the historical `FD_check_<case>.csv`.
  subroutine compute_sensitivity_directional(problem, sim, des, &
       target_sensitivities, perturbations, tolerance, file_name, &
       is_objective, gs_h, central_difference, strict_options, &
       constraint_index, name_suffix)
    class(problem_t), intent(inout) :: problem
    type(simulation_t), intent(inout) :: sim
    class(design_t), intent(inout) :: des
    type(vector_t), intent(inout) :: target_sensitivities
    real(kind=rp), intent(in) :: perturbations(:)
    real(kind=rp), intent(in) :: tolerance
    character(len=*), intent(in) :: file_name
    logical, intent(in) :: is_objective
    type(gs_t), intent(inout) :: gs_h
    logical, intent(in), optional :: central_difference
    type(fd_strict_options_t), intent(in), optional :: strict_options
    integer, intent(in), optional :: constraint_index
    character(len=*), intent(in), optional :: name_suffix

    real(kind=rp), allocatable :: direction(:), copies(:)
    real(kind=rp) :: work_arr(1), gradient_norm, projection, max_component
    real(kind=rp) :: n_design, eps_sens
    type(fd_strict_options_t) :: strict
    character(len=64) :: suffix
    logical :: central
    integer :: j, n, i_constraint

    central = .false.
    if (present(central_difference)) central = central_difference
    if (present(strict_options)) strict = strict_options
    i_constraint = 1
    if (present(constraint_index)) i_constraint = constraint_index
    suffix = ''
    if (present(name_suffix)) suffix = name_suffix

    n = des%size()
    if (target_sensitivities%size() .lt. n) then
       call neko_error('compute_sensitivity_directional: the sensitivity ' // &
            'vector is shorter than the design')
    end if
    call fd_sync_to_host(target_sensitivities)
    eps_sens = fd_sensitivity_scale(target_sensitivities, n)

    allocate(direction(n))
    allocate(copies(n))

    ! Assemble the gradient. `target_sensitivities` holds one partial
    ! derivative per *copy* of a shared dof; GS_OP_ADD leaves their sum -- the
    ! derivative with respect to the single design variable -- at every copy.
    ! The direction is therefore single-valued by construction, which is both
    ! what makes it a perturbation of real design variables and what makes
    ! `fd_project` exact.
    do j = 1, n
       direction(j) = target_sensitivities%x(j)
    end do
    call gs_h%op(direction, n, GS_OP_ADD)

    ! Number of copies of each dof, built the way Neko builds multiplicity, so
    ! that the norm below counts each design variable exactly once instead of
    ! once per element touching it.
    copies = 1.0_rp
    call gs_h%op(copies, n, GS_OP_ADD)

    work_arr(1) = 0.0_rp
    do j = 1, n
       work_arr(1) = work_arr(1) + direction(j)*direction(j) / copies(j)
    end do
    gradient_norm = sqrt(glsum(work_arr, 1))

    if (.not. (gradient_norm .gt. 0.0_rp)) then
       call neko_error('The assembled sensitivity is identically zero, so ' // &
            'there is no direction to differentiate along and the ' // &
            'directional finite-difference test would compare zero with ' // &
            'zero. Check that the adjoint ran and that the design domain ' // &
            'influences the functional.')
    end if

    do j = 1, n
       direction(j) = direction(j) / gradient_norm
    end do

    projection = fd_project(target_sensitivities%x, direction, n)

    ! Reported so that an inconsistency between the normalisation and the
    ! projection is visible rather than silent: with a consistent inner
    ! product these two numbers are the same, because <g, g/||g||> = ||g||.
    work_arr(1) = 0.0_rp
    do j = 1, n
       work_arr(1) = work_arr(1) + 1.0_rp / copies(j)
    end do
    n_design = glsum(work_arr, 1)
    copies = abs(direction)
    max_component = glmax(copies, n)

    if (pe_rank .eq. 0) then
       write(*, '(A,E15.6E3,A,E15.6E3)') &
            ' FD directional: gradient 2-norm = ', gradient_norm, &
            '   projection onto s = ', projection
       write(*, '(A,E15.6E3,A,I0)') &
            ' FD directional: largest direction entry = ', max_component, &
            '   design variables = ', nint(n_design)
       write(*, '(A)') ' FD directional: inner product = Euclidean on the ' &
            // 'assembled nodal design'
       write(*, '(A)') ' FD directional: coefficients (the mass-matrix ' &
            // 'weighting is already inside g).'
    end if

    call fd_run_sweep(problem, sim, des, direction, projection, -1, &
         perturbations, tolerance, file_name, is_objective, central, .true., &
         eps_sens, strict, i_constraint, trim(suffix))

    deallocate(direction)
    deallocate(copies)

  end subroutine compute_sensitivity_directional

  !> Run a perturbation sweep along a given direction and assert that the
  !! finite-difference estimate agrees with the analytic directional
  !! derivative.
  !!
  !! The single perturbation path for both modes: the single-dof probe is the
  !! one-hot special case of it, so the two cannot drift apart.
  !!
  !! @param problem The problem supplying the objective/constraint value.
  !! @param sim The simulation driving each forward solve.
  !! @param des The design being perturbed.
  !! @param direction Single-valued perturbation direction, one entry per local
  !!        design degree of freedom.
  !! @param target_derivative The analytic derivative along `direction`.
  !! @param probe_index Local index of the probed dof in single-dof mode,
  !!        negative on ranks that do not own it and negative everywhere in
  !!        directional mode.
  !! @param perturbations The sweep of perturbation magnitudes, in any order.
  !! @param tolerance The largest acceptable relative error at the smallest
  !!        perturbation of the sweep.
  !! @param file_name The case file name, used to name the CSV log.
  !! @param is_objective True to test an objective, false a constraint.
  !! @param central True to use a central difference.
  !! @param directional True when the whole design is perturbed along
  !!        `direction`; false for the historical single-dof probe.
  !! @param eps_sens Floor below which the analytic derivative is degenerate
  !!        and the assertion is made on the absolute difference instead.
  !! @param strict Settings of the strict criterion. Disabled leaves the
  !!        historical assertion on the smallest perturbation in place.
  !! @param constraint_index 1-based index of the constraint to read when
  !!        `is_objective` is false.
  !! @param suffix Suffix inserted before `.csv` in both the sweep log and
  !!        the verdict log. Empty for a single-target run, which therefore
  !!        keeps the historical file names.
  subroutine fd_run_sweep(problem, sim, des, direction, target_derivative, &
       probe_index, perturbations, tolerance, file_name, is_objective, &
       central, directional, eps_sens, strict, constraint_index, suffix)
    class(problem_t), intent(inout) :: problem
    type(simulation_t), intent(inout) :: sim
    class(design_t), intent(inout) :: des
    real(kind=rp), intent(in) :: direction(:)
    real(kind=rp), intent(in) :: target_derivative
    integer, intent(in) :: probe_index
    real(kind=rp), intent(in) :: perturbations(:)
    real(kind=rp), intent(in) :: tolerance
    character(len=*), intent(in) :: file_name
    logical, intent(in) :: is_objective
    logical, intent(in) :: central
    logical, intent(in) :: directional
    real(kind=rp), intent(in) :: eps_sens
    type(fd_strict_options_t), intent(in) :: strict
    integer, intent(in) :: constraint_index
    character(len=*), intent(in) :: suffix

    character(len=*), parameter :: fmt_head = '(4X,A12,4X,A10,6X,A11,5X,A5,10X)'
    character(len=*), parameter :: fmt_data = '(4X,4E15.6E3)'

    integer :: n_perturbations, ip
    real(kind=rp) :: perturb, work_arr(1)
    type(vector_t) :: design_vector, design_perturbed, log_data, constraint_vec
    real(kind=rp) :: constraint, perturbed_constraint, minus_constraint
    real(kind=rp) :: restored_constraint
    real(kind=rp) :: fd_estimate, fd_error
    real(kind=rp) :: min_error, min_perturb, smallest_perturb, largest_perturb
    real(kind=rp) :: floor_error, design_value, limit_plus, limit_minus
    real(kind=rp), allocatable :: sweep_perturbs(:), sweep_errors(:)
    real(kind=rp), allocatable :: sweep_differences(:)
    type(fd_verdict_t) :: verdict
    logical :: degenerate
    character(len=16) :: value_str
    integer :: min_index, i_write, floor_index
    logical :: prefer_negative
    type(csv_file_t) :: logger
    integer :: n, name_len, slash
    !> Status of the rank-0 verdict-CSV write, reduced so that a failure is
    !! acted on by every rank rather than by rank 0 alone.
    integer :: write_status(1)

    n = des%size()

    ! Initialize the vectors
    call design_vector%init(n)
    call design_perturbed%init(n)
    call log_data%init(4)
    call constraint_vec%init(problem%get_n_constraints())

    ! Get the design vector for reference
    ! This is the design vector we will perturb
    call des%get_values(design_vector)

    if (is_objective) then
       call problem%get_objective_value(constraint)
    else
       call problem%get_constraint_values(constraint_vec)
       constraint = constraint_vec%x(constraint_index)
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(design_vector%x, design_vector%x_d, &
            design_vector%size(), DEVICE_TO_HOST, .true.)
    end if

    ! Broadcast the design value at the probed dof -- only the owning rank
    ! holds it, and it decides the sign of the one-sided step and is quoted in
    ! the central-difference error message, both of which every rank must
    ! compose identically.
    if (probe_index .gt. 0) then
       work_arr(1) = design_vector%x(probe_index)
    else
       work_arr(1) = 0.0_rp
    end if
    design_value = glsum(work_arr, 1)

    ! Historical rule, kept exactly: a single dof past the middle of its range
    ! is perturbed downwards, so that a one-sided step on a design initialised
    ! as a 0/1 indicator field stays inside the bounds. A whole-field direction
    ! has no such freedom -- flipping it would flip the sign of the projection
    ! it is compared against -- so it is clamped in magnitude instead.
    prefer_negative = (.not. directional) .and. (design_value .gt. &
         0.5_rp * (fd_design_lower + fd_design_upper))

    call fd_step_headroom(design_vector%x, direction, n, limit_plus, &
         limit_minus)

    call fd_report_probe(des, n, probe_index, design_value, directional)

    if (pe_rank .eq. 0) then
       write(*, fmt_head) "Perturbation", "Constraint", "FD Estimate", "Error"
       write(*, fmt_data) 0.0_rp, constraint, target_derivative, 0.0_rp
    end if

    ! Init the csv writer. Only rank 0 owns the file: under MPI every rank
    ! reaches this code (the surrounding solve is collective), and without
    ! this guard multiple ranks independently opening/writing the same path
    ! is a genuine race — depending on scheduling it can silently clobber
    ! down to one rank's data (usually harmless, since all ranks compute
    ! identical globally-reduced rows) or, observed in practice on a long
    ! (~1000-timestep) run, duplicate every row once per rank.
    !
    ! The passed file name may carry a directory prefix (a stray future
    ! caller might pass one, e.g. 'cases/<name>.case'); strip any leading
    ! directory and the '.case' extension. Today both unit and regression
    ! drivers pass a bare '<name>.case'.
    if (pe_rank .eq. 0) then
       name_len = len_trim(file_name)
       slash = index(file_name(:name_len), '/', back = .true.)
       call logger%init('FD_check_' // &
            trim(file_name(slash+1:name_len-5)) // trim(suffix) // '.csv')
       call logger%set_header('perturbation,F,dFdx,error')
    end if

    n_perturbations = size(perturbations)
    ! Collective-safe: the sweep size is the same on every rank.
    if (n_perturbations .lt. 1) then
       call neko_error('The perturbation sweep is empty')
    end if

    allocate(sweep_perturbs(n_perturbations))
    allocate(sweep_errors(n_perturbations))
    allocate(sweep_differences(n_perturbations))
    i_write = 0

    do ip = 1, n_perturbations

       ! Decide the step. Every quantity it is built from -- `design_value`
       ! and the two headroom limits -- is already globally reduced, so every
       ! rank holding a copy of a shared dof necessarily applies the *same*
       ! perturbation without a further broadcast.
       if (central) then
          ! A central difference needs room on *both* sides of the design,
          ! so the symmetric half-step is clamped to whichever bound is
          ! nearer. The alternative -- stepping outside the bounds -- would
          ! evaluate a design that does not exist.
          perturb = min(abs(perturbations(ip)), limit_plus, limit_minus)
       else
          ! Ensure the perturbation stays within the bounds
          perturb = abs(perturbations(ip))
          if (prefer_negative) perturb = -perturb
          if (directional) perturb = min(perturb, limit_plus)
       end if

       ! Clamping diagnostics, shared by every path that clamps: the
       ! central-difference half-step, and the one-sided directional step.
       ! The historical one-sided single-dof step is not clamped at all -- it
       ! keeps itself in bounds by its sign -- and so is deliberately excluded.
       if (central .or. directional) then
          if (perturb .le. 0.0_rp) then
             if (directional) then
                call neko_error('No room to step along the sensitivity ' // &
                     'direction: the design already sits on a bound that ' // &
                     'the direction pushes it through. The directional ' // &
                     'sweep needs a case whose design is initialised ' // &
                     'strictly inside its bounds, rather than as a 0/1 ' // &
                     'indicator field.')
             else
                ! Note that a design initialised as a 0/1 indicator field --
                ! as most cases in this suite are -- puts *every* dof on a
                ! bound, so this is the normal outcome there rather than a
                ! corner case.
                write(value_str, '(F12.6)') design_value
                call neko_error('Central difference impossible: the ' // &
                     'design variable being perturbed is ' // &
                     trim(adjustl(value_str)) // ', on or outside its ' // &
                     'bounds, so there is no room for a symmetric step. ' // &
                     'Either initialise the design strictly inside its ' // &
                     'bounds, or use the one-sided difference (unset ' // &
                     'NEKO_TOP_FD_CENTRAL / set ' // &
                     'fd_test_central_difference false).')
             end if
          end if

          ! A step clamped far below the requested perturbation no longer
          ! probes the sweep point that was asked for: the design is so close
          ! to a bound that the finite difference is round-off noise at a
          ! perturbation nobody requested, and its error could spuriously
          ! become the sweep minimum. Choice made here: abort rather than
          ! silently skip the point, so the sweep that ran is always the sweep
          ! that was requested. Collective-safe: every input to `perturb` is
          ! globally reduced.
          if (abs(perturb) .lt. 1e-3_rp * abs(perturbations(ip))) then
             call neko_error('Finite-difference step clamped below 1/1000 ' // &
                  'of the requested perturbation: the design is too close ' // &
                  'to a bound for this sweep point to be meaningful. ' // &
                  'Start the sweep at smaller perturbations or move the ' // &
                  'design away from its bounds.')
          end if

          ! A clamped step is a different sweep point from the one requested,
          ! so say so rather than let the CSV be read as if the requested
          ! value had been used.
          if (pe_rank .eq. 0 .and. &
               abs(perturb) .lt. (1.0_rp - 1e-12_rp)*abs(perturbations(ip))) &
               then
             write(*, '(A,E15.6E3,A,E15.6E3,A)') &
                  ' FD sweep: requested step', abs(perturbations(ip)), &
                  ' clamped to', abs(perturb), ' by the design bounds.'
          end if
       end if

       if (central) then
          i_write = i_write + 1
          call evaluate_perturbed(problem, sim, des, design_vector, &
               design_perturbed, direction, perturb, is_objective, &
               constraint_index, i_write, constraint_vec, &
               perturbed_constraint)

          i_write = i_write + 1
          call evaluate_perturbed(problem, sim, des, design_vector, &
               design_perturbed, direction, -perturb, is_objective, &
               constraint_index, i_write, constraint_vec, minus_constraint)

          fd_estimate = perturbed_constraint - minus_constraint
          if (.not. abscmp(fd_estimate, 0.0_rp)) then
             fd_estimate = fd_estimate / (2.0_rp * perturb)
          end if
       else
          i_write = i_write + 1
          call evaluate_perturbed(problem, sim, des, design_vector, &
               design_perturbed, direction, perturb, is_objective, &
               constraint_index, i_write, constraint_vec, &
               perturbed_constraint)

          fd_estimate = perturbed_constraint - constraint
          if (.not. abscmp(fd_estimate, 0.0_rp)) then
             fd_estimate = fd_estimate / perturb
          end if
       end if

       fd_error = relative_error(fd_estimate, target_derivative)

       ! The logged `F` is the forward leg in both variants, so the CSV schema
       ! is the same either way.
       if (pe_rank .eq. 0) then
          write(*, fmt_data) perturb, perturbed_constraint, fd_estimate, &
               fd_error
       end if
       if (pe_rank .eq. 0) then
          log_data%x(1) = perturb
          log_data%x(2) = perturbed_constraint
          log_data%x(3) = fd_estimate
          log_data%x(4) = fd_error
          call logger%write(log_data)
       end if

       ! Record the whole sweep so that `report_sweep` can describe it. The
       ! unnormalised difference is kept beside the relative error because a
       ! degenerate analytic derivative is asserted on in absolute terms,
       ! where the relative error means nothing.
       sweep_perturbs(ip) = perturb
       sweep_errors(ip) = fd_error
       sweep_differences(ip) = fd_estimate - target_derivative
    end do

    ! Restore the unperturbed design: the loop leaves the last perturbed
    ! design in place, and anything run after this sweep (another probe,
    ! another sweep) must start from the baseline, not from a leftover
    ! perturbation. The zero-step evaluation doubles as a free
    ! reproducibility check: it should reproduce the stored baseline value.
    ! Reported only, never asserted on.
    i_write = i_write + 1
    call evaluate_perturbed(problem, sim, des, design_vector, &
         design_perturbed, direction, 0.0_rp, is_objective, &
         constraint_index, i_write, constraint_vec, restored_constraint)
    if (pe_rank .eq. 0) then
       write(*, '(A,E15.6E3)') ' FD sweep: baseline restored; ' // &
            're-evaluated functional minus stored baseline =', &
            restored_constraint - constraint
    end if

    ! Report the minimum |error| over the sweep and whether it is interior to
    ! the sweep. Diagnostic output only -- the assertion below is the
    ! historical one, on the smallest perturbation, and is unaffected by this.
    min_index = 1
    min_error = sweep_errors(1)
    min_perturb = sweep_perturbs(1)
    do ip = 2, n_perturbations
       if (abs(sweep_errors(ip)) .lt. abs(min_error)) then
          min_error = sweep_errors(ip)
          min_perturb = sweep_perturbs(ip)
          min_index = ip
       end if
    end do

    smallest_perturb = minval(abs(sweep_perturbs))
    largest_perturb = maxval(abs(sweep_perturbs))

    call report_sweep(min_error, min_perturb, smallest_perturb, &
         largest_perturb, min_index, n_perturbations)

    ! Everything from here on is built from globally reduced quantities and
    ! is therefore bit-identical on every rank, so every rank runs the whole
    ! analysis -- outside any rank guard, with only the writes guarded -- and
    ! reaches the same verdict without a single collective. `minloc` returns
    ! the first occurrence of the minimum, matching the first-wins tie-break
    ! the historical assertion has always used.
    floor_index = minloc(abs(sweep_perturbs), dim = 1)
    floor_error = sweep_errors(floor_index)
    degenerate = abs(target_derivative) .le. eps_sens

    if (rp .eq. sp) then
       ! Recorded before the banner and outside the rank guard: every rank
       ! takes this branch (the precision is a compile-time property), and
       ! the driver turns it into a non-zero exit code so that CTest reports
       ! a skip rather than a pass.
       fd_skipped_assertion = .true.

       ! In single precision the functional itself is only reproducible to
       ! about 1e-7 relative, so the smallest perturbation at which a
       ! one-sided difference resolves anything is of order sqrt(1e-7) ~ 0.34
       ! -- larger than the whole sweep. Any verdict here would be a verdict
       ! about round-off, so refuse to pretend otherwise.
       if (pe_rank .eq. 0) then
          write(*, '(A)') ' FD sweep: SKIPPED -- this is a single-precision' &
               // ' build. The functional is reproducible to only'
          write(*, '(A)') ' FD sweep: ~1e-7 relative, which puts the' &
               // ' smallest usable perturbation above 0.3, beyond any'
          write(*, '(A)') ' FD sweep: sweep. The finite-difference' &
               // ' assertion is not enforceable and is NOT being made.'
          write(*, '(A)') ' FD sweep: Configure Neko with --enable-real' &
               // ' set to dp and rebuild to gate on this test.'
       end if
    else
       if (strict%enabled) then
          call fd_evaluate(sweep_perturbs, sweep_errors, tolerance, strict, &
               verdict, degenerate)
          write_status(1) = 0
          if (pe_rank .eq. 0) then
             call fd_verdict_print(verdict, tolerance)
             call fd_write_verdict(file_name, suffix, verdict, &
                  write_status(1))
          end if

          ! Raised outside the rank guard, on a globally reduced status:
          ! `neko_error` is an `error stop`, so raising it inside the guard
          ! would abort rank 0 while every other rank walked into the next
          ! collective -- a hang rather than a failure.
          if (glimax(write_status, 1) .ne. 0) then
             call neko_error('Could not open the strict-verdict CSV ' // &
                  'for this case for writing')
          end if
       end if

       if (degenerate) then
          ! The relative error is normalised by a guard rather than by the
          ! sensitivity, so it is not a relative error and must not be
          ! asserted on as one. The absolute difference still means
          ! something, so that is what is gated on instead.
          if (pe_rank .eq. 0) then
             write(*, '(A,E15.6E3,A,E15.6E3)') &
                  ' FD sweep: DEGENERATE SENSITIVITY -- the analytic ' // &
                  'derivative', target_derivative, ' is at or below the ' // &
                  'field floor', eps_sens
             write(*, '(A)') ' FD sweep: so the logged relative error is ' &
                  // 'normalised by that floor and is not a relative'
             write(*, '(A,E15.6E3)') ' FD sweep: error. Asserting on the ' &
                  // 'absolute difference instead:', &
                  sweep_differences(floor_index)
          end if
          if (abs(sweep_differences(floor_index)) .gt. tolerance) then
             call neko_error('Finite difference estimate does not match ' // &
                  'sensitivity: the analytic sensitivity is degenerate, ' // &
                  'and the absolute difference at the smallest ' // &
                  'perturbation exceeds the tolerance')
          end if

       else if (strict%enabled) then
          call fd_assert_verdict(verdict)

       else if (abs(floor_error) .gt. tolerance) then
          call neko_error('Finite difference estimate does not match ' // &
               'sensitivity')
       end if
    end if

    ! Free the internal vectors
    call design_vector%free()
    call design_perturbed%free()
    call log_data%free()
    call constraint_vec%free()
    deallocate(sweep_perturbs, sweep_errors, sweep_differences)

  end subroutine fd_run_sweep

  !> Turn a strict verdict that did not certify into a failure, with a
  !! message that says which kind of failure it is.
  !!
  !! The distinction is the whole point of the strict criterion: a floor that
  !! exceeds the tolerance is a statement about the *gradient*, while an
  !! inadequate sweep or an unreachable tolerance are statements about the
  !! *test*. Reporting all three as one red is how a bad sweep gets read as a
  !! bad adjoint.
  !!
  !! @param verdict The verdict to act on.
  subroutine fd_assert_verdict(verdict)
    type(fd_verdict_t), intent(in) :: verdict

    if (verdict%passed) return

    select case (verdict%status)
    case (FD_STATUS_FLOOR_EXCEEDED)
       call neko_error('Finite difference estimate does not match ' // &
            'sensitivity: the finite-difference floor of this sweep ' // &
            'exceeds the tolerance. The FD strict lines above give the ' // &
            'floor, the measured truncation order and the window it was ' // &
            'taken over.')
    case (FD_STATUS_UNREACHABLE)
       call neko_error('The finite-difference tolerance asked of this ' // &
            'case is below what the functional''s own reproducibility ' // &
            'permits, so no sweep of any depth can decide the gradient ' // &
            'at it. The FD strict lines above give the smallest reachable ' &
            // 'tolerance. This is NOT evidence that the gradient is wrong.')
    case default
       call neko_error('The finite-difference SWEEP is inadequate to ' // &
            'decide this gradient: it contains neither a plateau that ' // &
            'bounds the bias nor a truncation run of the expected order. ' // &
            'Widen or deepen the sweep ' // &
            '(optimization.fd_test_perturbations). This is NOT evidence ' // &
            'that the gradient is wrong.')
    end select

  end subroutine fd_assert_verdict

  !> Append a strict verdict to `FD_verdict_<case>.csv`.
  !!
  !! Deliberately a *separate* file from `FD_check_<case>.csv`, whose schema
  !! and contents are left exactly as they were: that file is compared
  !! against reference data, so a column added to it is a broken comparison.
  !!
  !! One row per sweep, appended, matching how the sweep log itself
  !! accumulates. Rank 0 only -- every rank holds the same verdict, and
  !! several ranks opening one path is a genuine race.
  !!
  !! The open status is *returned* rather than raised here: this runs inside
  !! a rank-0 guard, and `neko_error` is an `error stop`, so failing here
  !! would take rank 0 down alone and hang the others in the next collective.
  !! The caller reduces the status and raises the error on every rank.
  !!
  !! @param file_name The case file name, used to name the CSV.
  !! @param suffix Suffix inserted before `.csv`, empty for a single-target
  !!        run so that the historical name is kept.
  !! @param verdict The verdict to record.
  !! @param write_status Zero on success, the `iostat` of the failed open
  !!        otherwise.
  subroutine fd_write_verdict(file_name, suffix, verdict, write_status)
    character(len=*), intent(in) :: file_name
    character(len=*), intent(in) :: suffix
    type(fd_verdict_t), intent(in) :: verdict
    integer, intent(out) :: write_status

    !> Field separator, written as a named constant rather than inline so
    !! that the row and its header cannot disagree about it.
    character(len=1), parameter :: sep = ','
    character(len=*), parameter :: header = &
         'p_hat' // sep // 'C_hat' // sep // 'branch' // sep // &
         'status' // sep // 'bracketed' // sep // 'min_abs_error' // sep // &
         'min_perturbation' // sep // 'n_truncation_points' // sep // &
         'n_sign_crossings' // sep // 'C_hat_kind' // sep // 'tol_min'
    character(len=512) :: path, row
    character(len=32) :: p_str, c_str, min_str, pert_str, tol_str
    character(len=32) :: trunc_str, cross_str
    integer :: unit_id, ios, name_len, slash
    logical :: exists

    name_len = len_trim(file_name)
    slash = index(file_name(:name_len), '/', back = .true.)
    path = 'FD_verdict_' // trim(file_name(slash+1:name_len-5)) // &
         trim(suffix) // '.csv'

    ! A field the criterion could not measure is written as `nan` rather
    ! than as a zero that would read as a measurement.
    p_str = 'nan'
    if (verdict%has_p_hat) write(p_str, '(E17.10E3)') verdict%p_hat
    c_str = 'nan'
    if (verdict%has_c_hat) write(c_str, '(E17.10E3)') verdict%c_hat
    tol_str = 'nan'
    if (verdict%has_tol_min) write(tol_str, '(E17.10E3)') verdict%tol_min
    write(min_str, '(E17.10E3)') verdict%min_abs_error
    write(pert_str, '(E17.10E3)') verdict%min_perturbation
    write(trunc_str, '(I0)') verdict%n_truncation_points
    write(cross_str, '(I0)') verdict%n_sign_crossings

    row = trim(adjustl(p_str)) // sep // trim(adjustl(c_str)) // sep // &
         trim(fd_branch_name(verdict%branch)) // sep // &
         trim(fd_status_name(verdict%status)) // sep
    if (verdict%bracketed) then
       row = trim(row) // 'true' // sep
    else
       row = trim(row) // 'false' // sep
    end if
    row = trim(row) // trim(adjustl(min_str)) // sep // &
         trim(adjustl(pert_str)) // sep // trim(adjustl(trunc_str)) // &
         sep // trim(adjustl(cross_str)) // sep // &
         trim(fd_c_hat_kind(verdict)) // sep // trim(adjustl(tol_str))

    inquire(file = trim(path), exist = exists)
    open(newunit = unit_id, file = trim(path), action = 'write', &
         position = 'append', status = 'unknown', iostat = ios)
    write_status = ios
    if (ios .ne. 0) return
    if (.not. exists) write(unit_id, '(A)') header
    write(unit_id, '(A)') trim(row)
    close(unit_id)

  end subroutine fd_write_verdict

  !> Evaluate the objective (or constraint) at the design displaced by
  !! `step * direction`.
  !!
  !! The one perturbation path: a single-dof probe passes a one-hot
  !! `direction`, the Taylor test passes a whole normalised gradient field, and
  !! nothing downstream distinguishes them. Factored out of the sweep so that
  !! the central-difference variant can reuse it for both the \f$x + h s\f$ and
  !! \f$x - h s\f$ legs, and the sweep can reuse it once more with a zero step
  !! to restore the baseline design.
  !!
  !! The design field is stored per element, (lx, ly, lz, nelv), so a dof on an
  !! element interface has one copy per adjoining element, possibly on
  !! different MPI ranks. `direction` must therefore be **single-valued**:
  !! equal on every copy of a shared dof, on every rank holding one. Both
  !! callers guarantee that by gather-scattering it with GS_OP_ADD -- the same
  !! mechanism Neko uses to propagate values between elements. Displacing only
  !! one copy would leave the design multi-valued at that point, which is not a
  !! perturbation of any real design variable.
  !!
  !! @param problem The problem supplying the objective/constraint value.
  !! @param sim The simulation driving the forward solve.
  !! @param des The design, updated in place with the perturbed values.
  !! @param design_vector The unperturbed design values.
  !! @param design_perturbed Work vector receiving the perturbed design.
  !! @param direction Single-valued perturbation direction; zero entries are
  !!        left untouched.
  !! @param step The signed step scaling `direction`.
  !! @param is_objective True to read the objective, false the constraint.
  !! @param constraint_index 1-based index of the constraint to read when
  !!        `is_objective` is false.
  !! @param write_index Output index handed to `sim%write`.
  !! @param constraint_vec Work vector for the constraint values.
  !! @param functional The resulting objective or constraint value.
  subroutine evaluate_perturbed(problem, sim, des, design_vector, &
       design_perturbed, direction, step, is_objective, constraint_index, &
       write_index, constraint_vec, functional)
    class(problem_t), intent(inout) :: problem
    type(simulation_t), intent(inout) :: sim
    class(design_t), intent(inout) :: des
    type(vector_t), intent(in) :: design_vector
    type(vector_t), intent(inout) :: design_perturbed
    real(kind=rp), intent(in) :: direction(:)
    real(kind=rp), intent(in) :: step
    logical, intent(in) :: is_objective
    integer, intent(in) :: constraint_index
    integer, intent(in) :: write_index
    type(vector_t), intent(inout) :: constraint_vec
    real(kind=rp), intent(out) :: functional

    integer :: j

    ! Reset the design field
    design_perturbed%x = design_vector%x

    ! Displace along the direction, keeping the design single-valued. Entries
    ! with a zero direction are skipped rather than having zero added, so that
    ! the untouched part of the design is bit-for-bit the baseline.
    do j = 1, size(direction)
       if (abs(direction(j)) .gt. 0.0_rp) then
          design_perturbed%x(j) = design_vector%x(j) + step * direction(j)
       end if
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(design_perturbed%x, design_perturbed%x_d, &
            design_perturbed%size(), HOST_TO_DEVICE, .true.)
    end if
    call des%update_design(design_perturbed)

    ! Compute the objective value of the perturbed design
    call problem%compute(des, sim)
    if (is_objective) then
       call problem%get_objective_value(functional)
    else
       call problem%get_constraint_values(constraint_vec)
       functional = constraint_vec%x(constraint_index)
    end if
    call sim%write(write_index)
    call sim%reset()

  end subroutine evaluate_perturbed

  !> Report the outcome of a perturbation sweep on rank 0.
  !!
  !! Purely diagnostic: this routine prints, it does not assert. The gating
  !! assertion is made by the caller, on the error at the smallest
  !! perturbation.
  !!
  !! States whether the minimum error is interior to the sweep. If it is, the
  !! error fell and rose again around it and the finite-difference floor has
  !! genuinely been bracketed. If instead the minimum sits at an *end* of the
  !! sweep, nothing was bracketed: at the smallest perturbation the sweep
  !! stopped before the round-off upturn, and at the largest the error was
  !! still falling as the perturbation grew. Either way the number reported
  !! is only a bound on the floor, which must be said rather than left to be
  !! assumed.
  !!
  !! @param min_error The signed relative error at the minimum of the sweep.
  !! @param min_perturb The perturbation at which that minimum occurred.
  !! @param smallest_perturb The smallest perturbation the sweep reached.
  !! @param largest_perturb The largest perturbation the sweep reached.
  !! @param min_index Index into the sweep of the minimum.
  !! @param n_perturbations Number of points in the sweep.
  subroutine report_sweep(min_error, min_perturb, smallest_perturb, &
       largest_perturb, min_index, n_perturbations)
    real(kind=rp), intent(in) :: min_error, min_perturb, smallest_perturb
    real(kind=rp), intent(in) :: largest_perturb
    integer, intent(in) :: min_index, n_perturbations

    logical :: at_smallest, at_largest

    if (pe_rank .ne. 0) return

    write(*, '(A,E15.6E3,A,E15.6E3,A,I0,A,I0,A)') &
         ' FD sweep: minimum |error| of', abs(min_error), &
         ' at perturbation', min_perturb, &
         ' (point ', min_index, ' of ', n_perturbations, ')'

    ! Compared by magnitude rather than by index, so that a sweep whose
    ! perturbations collapse to the same effective value -- which the
    ! central-difference clamping can do -- is not misreported as bracketed.
    at_smallest = abs(min_perturb) .le. abs(smallest_perturb)
    at_largest = abs(min_perturb) .ge. abs(largest_perturb)

    if (at_smallest .and. at_largest) then
       write(*, '(A)') ' FD sweep: NOT BRACKETED -- the sweep has too few' &
            // ' distinct perturbations to bracket a floor at all; the value'
       write(*, '(A)') ' FD sweep: above is a single sample, not a floor.' &
            // ' Sweep several decades of perturbations to bracket one.'
    else if (at_smallest) then
       write(*, '(A)') ' FD sweep: NOT BRACKETED -- the minimum lies at the' &
            // ' smallest perturbation of the sweep, so no round-off upturn'
       write(*, '(A)') ' FD sweep: was observed. The value above is an upper' &
            // ' bound on the finite-difference floor, not the floor itself;'
       write(*, '(A)') ' FD sweep: extend the sweep to smaller perturbations' &
            // ' to bracket it.'
    else if (at_largest) then
       write(*, '(A)') ' FD sweep: NOT BRACKETED -- the minimum lies at the' &
            // ' largest perturbation of the sweep: the error was still'
       write(*, '(A)') ' FD sweep: falling as the perturbation grew, so no' &
            // ' floor was observed and the value above is only a bound on'
       write(*, '(A)') ' FD sweep: it. Extend the sweep to larger' &
            // ' perturbations to bracket the floor.'
    else
       write(*, '(A)') ' FD sweep: BRACKETED -- the error falls and rises' &
            // ' again around the minimum (round-off upturn observed at the'
       write(*, '(A)') ' FD sweep: small end), so the minimum above is a' &
            // ' genuine finite-difference floor.'
    end if

  end subroutine report_sweep

  !> Sweep the perturbations at each of a list of design degrees of freedom.
  !!
  !! @param problem The problem supplying the objective/constraint value.
  !! @param sim The simulation driving each forward solve.
  !! @param des The design being perturbed.
  !! @param target_sensitivities The analytic sensitivities to check against.
  !! @param list Local indices of the dofs to perturb.
  !! @param perturbations The sweep of perturbation magnitudes.
  !! @param tolerance The largest acceptable relative error at the smallest
  !!        perturbation of the sweep.
  !! @param file_name The case file name, used to name the CSV log.
  !! @param is_objective True to test an objective, false a constraint.
  !! @param gs_h Gather-scatter handle for the design's dofmap.
  !! @param central_difference (Optional) True to use a central difference.
  !! @param strict_options (Optional) Settings of the strict criterion.
  !! @param constraint_index (Optional) 1-based index of the constraint to
  !!        read when `is_objective` is false. Defaults to 1.
  !! @param name_suffix (Optional) Suffix for the CSV file names.
  subroutine compute_sensitivity_list(problem, sim, des, target_sensitivities, &
       list, perturbations, tolerance, file_name, is_objective, gs_h, &
       central_difference, strict_options, constraint_index, name_suffix)
    class(problem_t), intent(inout) :: problem
    type(simulation_t), intent(inout) :: sim
    class(design_t), intent(inout) :: des
    type(vector_t), intent(inout) :: target_sensitivities
    integer, dimension(:), intent(in) :: list
    real(kind=rp), dimension(:), intent(in) :: perturbations
    real(kind=rp), intent(in) :: tolerance
    character(len=*), intent(in) :: file_name
    logical, intent(in) :: is_objective
    type(gs_t), intent(inout) :: gs_h
    logical, intent(in), optional :: central_difference
    type(fd_strict_options_t), intent(in), optional :: strict_options
    integer, intent(in), optional :: constraint_index
    character(len=*), intent(in), optional :: name_suffix

    integer :: i, n

    n = size(list)
    do i = 1, n
       call compute_sensitivity_i(problem, sim, des, target_sensitivities, &
            list(i), perturbations, tolerance, file_name, is_objective, gs_h, &
            central_difference, strict_options, constraint_index, name_suffix)
    end do
  end subroutine compute_sensitivity_list

  !> Computes the relative difference between two numbers
  !! \f$ \frac{a - b}{|b|} \f$
  !> Whether any sweep in this run declined to make its assertion.
  !!
  !! Read by the drivers, which exit with `FD_SKIP_EXIT_CODE` when it is
  !! true: a gradient check that silently performed no check must not be
  !! reported as a pass.
  !!
  !! @return True if an assertion was skipped rather than made.
  function fd_assertion_skipped() result(skipped)
    logical :: skipped

    skipped = fd_skipped_assertion

  end function fd_assertion_skipped

  function relative_error(a, b) result(err)
    real(kind=rp), intent(in) :: a, b
    real(kind=rp) :: err

    err = (a - b) / max(abs(b), NEKO_EPS)

  end function relative_error

end module sensitivity
