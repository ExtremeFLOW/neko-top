!> Strict assertion criterion for a finite-difference sensitivity sweep.
!!
!! The historical assertion the sensitivity harness makes -- that the signed
!! relative error at the *smallest* perturbation of the sweep is within the
!! tolerance -- is a false green whenever that error crosses zero inside the
!! sweep: a gradient with a genuine bias \f$C\f$ can be sampled arbitrarily
!! close to the crossing and pass. This module implements the strict rule that
!! replaces it, and is deliberately kept free of MPI and of I/O so that it can
!! be exercised on recorded sweeps outside a simulation.
!!
!! **Model.** Over a sweep of perturbation magnitudes \f$m_k\f$ the signed
!! relative error behaves as
!! \f[ e(m) = C + A m^p + B m^{p+1} + \mathrm{noise}/m, \f]
!! with \f$p = 1\f$ for a one-sided and \f$p = 2\f$ for a central difference.
!! The quantity under test is the bias \f$C\f$: the value the error converges
!! to once truncation has vanished but before round-off takes over.
!!
!! **Two certificates, both always computed; the sweep passes if either one
!! certifies \f$|C| \le \mathrm{tolerance}\f$.**
!!
!!  * `BOUNDED`: a window of at least `fd_plateau_points` consecutive points
!!    whose signed errors span no more than `plateau_fraction * tolerance`.
!!    Truncation is negligible across such a window, so the window *bounds*
!!    the bias: \f$|C| \le |\mathrm{mean}| + k\,\mathrm{spread}\f$ with
!!    \f$k = \f$ `fd_spread_multiplier`. This branch is load-bearing rather
!!    than a convenience -- measured on the four checked-in unit cases,
!!    `volume` and `volume_filtered` are exactly linear in the design
!!    variable, have no truncation branch at all and find no order run at any
!!    order tolerance up to 2.0, so without `BOUNDED` they are two false reds
!!    on must-pass cases.
!!  * `FIT`: the differences \f$d_k = e_k - e_{k+1}\f$ cancel \f$C\f$ exactly,
!!    so the ratio \f$d_k/d_{k+1}\f$ measures the truncation order without
!!    knowing \f$C\f$. Over the longest run of consecutive ratios that are
!!    single-signed and of the expected order, \f$C\f$ is recovered by
!!    Richardson extrapolation from the two smallest points of the run and
!!    cross-checked against a three-term least-squares fit.
!!
!! `BOUNDED` yields a *bound* on the bias and `FIT` an *estimate* of it, which
!! is why they are reported under different labels; where both exist the
!! estimate is the one quoted, because it is the more accurate of the two.
!! When neither certificate can be formed at all the verdict is
!! `INCONCLUSIVE`: a loud failure saying the *sweep* is inadequate, which is
!! explicitly not the same statement as the gradient being wrong.
!!
!! **The order is solved for, not read off a single sweep ratio.** The closed
!! form \f$\hat p = \log(d_k/d_{k+1})/\log r\f$ needs a geometric sweep, and
!! the sweep the unit drivers actually run,
!! `[5e-1, 1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 5e-4, 1e-4]`, is not one: its
!! successive ratios alternate 5, 2, 5, 2. Applying the closed form there
!! returns a wrong order with no indication that anything is wrong. The
!! general relation
!! \f[ R(p) = \frac{m_0^p - m_1^p}{m_1^p - m_2^p} \f]
!! is strictly increasing in \f$p\f$, so it is inverted by bisection instead;
!! on a geometric sweep that reduces exactly, not approximately, to the
!! closed form.
!!
!! **No bracketing gate.** Whether the sweep brackets the round-off upturn is
!! reported and never gated on. Measured over the recorded sweeps, requiring
!! an upturn rejects the very sweep that exposed the false green this rule
!! exists to fix, and lifts the false-red rate from 0.012 to 0.73 with no
!! value of the bracketing factor recovering it. Sign crossings of \f$e\f$ are
!! likewise counted and reported but never used to exclude points.
!!
!! **What the sweep itself must satisfy.** Three properties are demanded of
!! the recorded sweep before any branch is evaluated, and each is refused
!! rather than worked around: it must be finite (a NaN makes every comparison
!! false and quietly shortens the analysis), single-signed (the model is the
!! truncation series of *one* one-sided difference), and its magnitudes must
!! be distinct (a repeated magnitude is one measurement counted twice, and
!! adds plateau length at zero spread).
module fd_criterion
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use num_types, only: rp
  use math, only: NEKO_EPS
  use utils, only: neko_error
  implicit none
  private

  !> Branch that produced the reported bias.
  integer, parameter, public :: FD_BRANCH_NONE = 0
  integer, parameter, public :: FD_BRANCH_BOUNDED = 1
  integer, parameter, public :: FD_BRANCH_FIT = 2

  !> Outcome of the criterion.
  integer, parameter, public :: FD_STATUS_OK = 0
  integer, parameter, public :: FD_STATUS_FLOOR_EXCEEDED = 1
  integer, parameter, public :: FD_STATUS_INCONCLUSIVE = 2
  integer, parameter, public :: FD_STATUS_DEGENERATE = 3
  integer, parameter, public :: FD_STATUS_UNREACHABLE = 4

  !> Smallest window, in sweep points, that the bounded branch will accept.
  integer, parameter, public :: fd_plateau_points = 3

  !> Smallest order run, in difference ratios, that the fit branch will
  !! accept. Two ratios span four consecutive sweep points.
  integer, parameter, public :: fd_min_order_ratios = 2

  !> Multiplier on the spread of a plateau before it is compared against the
  !! tolerance. Not 1: over a descending geometric window the mean of the
  !! window can sit \f$\mathrm{spread}/(1 - r^{-p(n-1)}) = 1.11\f$ spreads
  !! away from the bias, so a unit multiplier is optimistic by about 11%.
  !! Two is rigorous with margin and costs nothing measurable -- the four
  !! checked-in unit cases certify at 3.46e-3, 4.99e-3, 2.76e-9 and 1.23e-7
  !! against tolerances of 1e-2, 1e-2, 1e-5 and 1e-5.
  real(kind=rp), parameter, public :: fd_spread_multiplier = 2.0_rp

  !> Smallest magnitude ratio between the two points Richardson extrapolates
  !! from. Below this \f$r^p - 1\f$ is small enough that the extrapolation
  !! amplifies the difference of two nearly equal errors without limit.
  real(kind=rp), parameter, public :: fd_min_richardson_ratio = 1.3_rp

  !> Factor by which the Richardson and three-term floor estimates may differ
  !! before the disagreement is warned about. Warned about only: the assertion
  !! is made on the Richardson estimate either way.
  real(kind=rp), parameter, public :: fd_cross_check_factor = 3.0_rp

  !> Bracket of the bisection that recovers the truncation order from a
  !! difference ratio, and the number of bisections taken. Sixty halvings of
  !! [0, 6] land far below what `rp` can resolve, so the result is the exact
  !! root to working precision.
  real(kind=rp), parameter, public :: fd_order_lower = 0.0_rp
  real(kind=rp), parameter, public :: fd_order_upper = 6.0_rp
  integer, parameter, public :: fd_order_bisections = 60

  !> Fitted truncation slope below which a functional is taken to have no
  !! truncation branch at all -- exactly linear in the design variable, as
  !! the volume constraint is -- and the reachable-tolerance floor therefore
  !! does not apply.
  real(kind=rp), parameter, public :: fd_min_slope = 1e-12_rp

  !> Constant in the reachable-tolerance floor
  !! \f$\mathrm{tol}_{\min} = \sqrt{40 N |A|}\f$. It is
  !! \f$1/(\mathrm{plateau\_fraction} \times 0.1)\f$: certifying a floor at
  !! the tolerance needs the truncation term below
  !! `plateau_fraction * tolerance`, which caps the smallest useful
  !! perturbation at \f$0.25\,\mathrm{tol}/|A|\f$, while round-off needs
  !! \f$N/m \le 0.1\,\mathrm{tol}\f$, which floors it at
  !! \f$N/(0.1\,\mathrm{tol})\f$. The two windows only overlap above
  !! \f$\sqrt{40 N |A|}\f$.
  real(kind=rp), parameter, public :: fd_reachability_constant = 40.0_rp

  !> Points used to fit the truncation slope, and to measure the noise, for
  !! the reachable-tolerance floor: the smallest of the sweep, where the
  !! functional is least contaminated by its own nonlinearity.
  integer, parameter, public :: fd_slope_points = 4
  integer, parameter, public :: fd_noise_points = 6

  !> Length of one diagnostic message line.
  integer, parameter, public :: fd_message_len = 100
  !> Most diagnostic lines a verdict can carry.
  integer, parameter, public :: fd_max_messages = 16

  !> Tunables of the strict criterion, all read from the case file.
  type, public :: fd_strict_options_t
     !> Apply the strict criterion instead of the historical assertion on the
     !! smallest perturbation.
     logical :: enabled = .false.
     !> Truncation order the sweep is expected to show.
     real(kind=rp) :: p_expected = 1.0_rp
     !> Half-width of the accepted band around `p_expected`.
     real(kind=rp) :: order_tolerance = 0.10_rp
     !> Width of an accepted plateau, as a fraction of the tolerance.
     real(kind=rp) :: plateau_fraction = 0.25_rp
  end type fd_strict_options_t

  !> The verdict of the criterion, together with the diagnostics reported
  !! beside it. Every field is set by `fd_evaluate`, so a verdict never
  !! carries a value left over from an earlier sweep.
  type, public :: fd_verdict_t
     !> Branch the reported bias came from, `FD_BRANCH_*`.
     integer :: branch = FD_BRANCH_NONE
     !> Outcome, `FD_STATUS_*`.
     integer :: status = FD_STATUS_INCONCLUSIVE
     !> True only if either certificate bounded the bias within the
     !! tolerance.
     logical :: passed = .false.
     !> True if `p_hat` holds a measured order.
     logical :: has_p_hat = .false.
     !> Measured truncation order over the order run.
     real(kind=rp) :: p_hat = 0.0_rp
     !> True if `c_hat` holds an estimate of, or a bound on, the bias.
     logical :: has_c_hat = .false.
     !> The reported bias: a Richardson *estimate* on the fit branch, a
     !! *bound* on the bounded branch. `fd_c_hat_kind` names which.
     real(kind=rp) :: c_hat = 0.0_rp
     !> True if `c_cross` holds the three-term cross-check.
     logical :: has_c_cross = .false.
     !> Three-term least-squares estimate of the bias.
     real(kind=rp) :: c_cross = 0.0_rp
     !> True if a plateau certified the bias.
     logical :: bounded_certified = .false.
     !> True if the fit branch certified the bias.
     logical :: fit_certified = .false.
     !> True if a plateau window was found at all.
     logical :: has_window = .false.
     !> Bound the plateau places on the bias, \f$|mean| + k\,spread\f$.
     real(kind=rp) :: window_bound = 0.0_rp
     !> True if the minimum \f$|e|\f$ is interior to the sweep. Diagnostic.
     logical :: bracketed = .false.
     !> Smallest \f$|e|\f$ over the sweep.
     real(kind=rp) :: min_abs_error = 0.0_rp
     !> Perturbation magnitude at which that minimum occurred.
     real(kind=rp) :: min_perturbation = 0.0_rp
     !> Sweep points covered by the order run, zero if none was found.
     integer :: n_truncation_points = 0
     !> Number of sign changes of the signed error across the sweep.
     integer :: n_sign_crossings = 0
     !> Largest perturbation magnitude of the window the bias came from.
     real(kind=rp) :: window_m_max = 0.0_rp
     !> Smallest perturbation magnitude of that window.
     real(kind=rp) :: window_m_min = 0.0_rp
     !> True if the reachable-tolerance floor could be measured.
     logical :: has_tol_min = .false.
     !> Smallest tolerance this functional's reproducibility permits.
     real(kind=rp) :: tol_min = 0.0_rp
     !> Number of diagnostic lines in `message`.
     integer :: n_messages = 0
     !> Diagnostic lines, printed by `fd_verdict_print`.
     character(len=fd_message_len) :: message(fd_max_messages) = ''
  end type fd_verdict_t

  public :: fd_evaluate, fd_verdict_print, fd_branch_name, fd_status_name, &
       fd_c_hat_kind

contains

  !> Apply the strict criterion to a recorded sweep.
  !!
  !! Collective-free by construction: every input is already globally reduced
  !! and therefore bit-identical on every rank, so every rank runs the whole
  !! analysis and reaches the same verdict without communicating. Nothing here
  !! prints; `fd_verdict_print` does that, and only its caller knows which
  !! rank owns the log.
  !!
  !! @param perturbations Signed perturbations of the sweep, in any order.
  !! @param errors Signed relative errors, in the same order.
  !! @param tolerance Largest acceptable magnitude of the bias.
  !! @param opts Tunables of the criterion.
  !! @param verdict The verdict and its diagnostics.
  !! @param degenerate (Optional) True when the analytic sensitivity is below
  !!        the field-relative floor, so the relative errors are not relative
  !!        errors at all. The diagnostics are still filled in, but no branch
  !!        is evaluated and the status is `FD_STATUS_DEGENERATE`; the caller
  !!        asserts on the absolute difference instead.
  subroutine fd_evaluate(perturbations, errors, tolerance, opts, verdict, &
       degenerate)
    real(kind=rp), intent(in) :: perturbations(:)
    real(kind=rp), intent(in) :: errors(:)
    real(kind=rp), intent(in) :: tolerance
    type(fd_strict_options_t), intent(in) :: opts
    type(fd_verdict_t), intent(out) :: verdict
    logical, intent(in), optional :: degenerate

    real(kind=rp) :: m(size(perturbations)), e(size(perturbations))
    real(kind=rp) :: p_hat(max(size(perturbations) - 2, 1))
    real(kind=rp) :: level, spread, c_rich, c_three, ratio
    integer :: n, k0, k1, i0, i1
    logical :: is_degenerate, has_run, has_fit, fit_ok
    character(len=fd_message_len) :: line, reason

    is_degenerate = .false.
    if (present(degenerate)) is_degenerate = degenerate

    n = size(perturbations)
    if (size(errors) .ne. n) then
       call neko_error('fd_evaluate: the perturbation and error lists ' // &
            'have different lengths')
    end if
    if (n .lt. 1) call neko_error('fd_evaluate: the sweep is empty')

    call fd_check_finite(perturbations, errors, n)
    call fd_check_single_sign(perturbations, n)
    call fd_sort_sweep(perturbations, errors, n, m, e)
    call fd_check_distinct(m, n)

    ! Diagnostics. Reported for every sweep and never gated on -- see the
    ! module header for why bracketing in particular must not be a gate.
    call fd_diagnostics(m, e, n, verdict)

    call fd_order_run(m, e, n, opts%p_expected, opts%order_tolerance, &
         p_hat, has_run, k0, k1, verdict%p_hat, reason)
    if (has_run) then
       verdict%has_p_hat = .true.
       verdict%n_truncation_points = k1 - k0 + 3
    end if

    if (is_degenerate) then
       verdict%branch = FD_BRANCH_NONE
       verdict%status = FD_STATUS_DEGENERATE
       verdict%passed = .false.
       call fd_add_message(verdict, 'DEGENERATE_SENSITIVITY: the analytic ' &
            // 'sensitivity is below the field-relative')
       call fd_add_message(verdict, 'floor, so the logged error is ' // &
            'normalised by a guard rather than by the')
       call fd_add_message(verdict, 'sensitivity. The assertion is made ' // &
            'on the absolute difference instead.')
       return
    end if

    call fd_reachable_tolerance(m, e, n, opts%p_expected, verdict)

    ! ---- BOUNDED --------------------------------------------------------
    call fd_bounded_window(e, n, tolerance, opts%plateau_fraction, &
         verdict%has_window, i0, i1, level, spread)
    if (verdict%has_window) then
       verdict%window_bound = abs(level) + fd_spread_multiplier * spread
       verdict%bounded_certified = verdict%window_bound .le. tolerance
    end if

    ! ---- FIT ------------------------------------------------------------
    has_fit = .false.
    if (has_run) then
       call fd_richardson(m, e, k1 + 1, verdict%p_hat, c_rich, has_fit, ratio)
       if (.not. has_fit) then
          call fd_add_message(verdict, 'An order run was found, but the two' &
               // ' smallest points of it are too')
          write(line, '(A,F8.4,A)') 'close in magnitude (ratio ', ratio, &
               ') for a stable Richardson step.'
          call fd_add_message(verdict, line)
       end if
    end if
    if (has_fit) verdict%fit_certified = abs(c_rich) .le. tolerance

    ! Both certificates are now in hand. The fit estimate is quoted where it
    ! exists because it is the more accurate of the two; the plateau only
    ! ever bounds the bias.
    if (has_fit) then
       verdict%branch = FD_BRANCH_FIT
       verdict%has_c_hat = .true.
       verdict%c_hat = c_rich
       verdict%window_m_max = m(k0)
       verdict%window_m_min = m(k1 + 2)

       call fd_fit_three_term(m, e, k0, k1 + 2, verdict%p_hat, c_three, &
            fit_ok)
       if (fit_ok) then
          verdict%has_c_cross = .true.
          verdict%c_cross = c_three
       end if

       write(line, '(A,I0,A,F8.4)') 'FIT: order run over ', &
            verdict%n_truncation_points, ' consecutive points, order ', &
            verdict%p_hat
       call fd_add_message(verdict, line)
       write(line, '(A,E13.6E3,A,E13.6E3,A)') &
            'Richardson bias taken over m in [', m(k1 + 2), ', ', &
            m(k1 + 1), ']'
       call fd_add_message(verdict, line)
       call fd_cross_check(verdict, c_rich, c_three)
    else if (verdict%has_window) then
       verdict%branch = FD_BRANCH_BOUNDED
       verdict%has_c_hat = .true.
       verdict%c_hat = level
       verdict%window_m_max = m(i0)
       verdict%window_m_min = m(i1)
       write(line, '(A,I0,A,E13.6E3)') 'BOUNDED: ', i1 - i0 + 1, &
            ' consecutive points bound the bias by ', verdict%window_bound
       call fd_add_message(verdict, line)
       call fd_add_message(verdict, 'The reported bias is a BOUND, not an ' &
            // 'estimate: over this window')
       call fd_add_message(verdict, 'truncation is too small to separate ' &
            // 'from the bias.')
    else if (len_trim(reason) .gt. 0) then
       ! `fd_order_run` leaves `reason` empty when it *did* find a run, so an
       ! unconditional append prints a bare ' FD strict: ' line.
       call fd_add_message(verdict, reason)
    end if

    ! ---- the verdict ----------------------------------------------------
    verdict%passed = verdict%bounded_certified .or. verdict%fit_certified

    if (verdict%has_tol_min .and. verdict%tol_min .gt. tolerance) then
       ! Whether the gradient is right cannot be decided at this tolerance
       ! by any sweep, so say that rather than report a bare red.
       verdict%status = FD_STATUS_UNREACHABLE
       verdict%passed = .false.
       write(line, '(A,E13.6E3)') 'TOLERANCE_UNREACHABLE: the requested ' &
            // 'tolerance is ', tolerance
       call fd_add_message(verdict, line)
       write(line, '(A,E13.6E3,A)') 'but this functional''s ' // &
            'reproducibility permits at best ', verdict%tol_min, '.'
       call fd_add_message(verdict, line)
       call fd_add_message(verdict, 'Loosen fd_test_tolerance to above ' // &
            'that. This is NOT evidence that')
       call fd_add_message(verdict, 'the gradient is wrong.')
    else if (verdict%passed) then
       verdict%status = FD_STATUS_OK
       call fd_report_upturn(verdict)
    else if (verdict%has_window .or. has_fit) then
       verdict%status = FD_STATUS_FLOOR_EXCEEDED
       write(line, '(A,E13.6E3,A,E13.6E3)') 'FLOOR_EXCEEDED: bias ', &
            verdict%c_hat, ' exceeds the tolerance ', tolerance
       call fd_add_message(verdict, line)
       write(line, '(A,E13.6E3,A,E13.6E3,A)') 'over the window m in [', &
            verdict%window_m_min, ', ', verdict%window_m_max, ']'
       call fd_add_message(verdict, line)
    else
       verdict%status = FD_STATUS_INCONCLUSIVE
       call fd_add_message(verdict, 'INCONCLUSIVE: the SWEEP is ' // &
            'inadequate for this case. Widen the')
       call fd_add_message(verdict, 'perturbation range or reduce the ' // &
            'noise. This is NOT evidence that')
       call fd_add_message(verdict, 'the gradient is wrong.')
       call fd_report_sign_pattern(e, n, verdict)
       call fd_report_roundoff(e, n, verdict)
    end if

  end subroutine fd_evaluate

  !> Reject a sweep carrying a non-finite error.
  !!
  !! Not a formality. Every comparison against a NaN is false, so a NaN error
  !! silently *shortens* the order run rather than breaking the analysis, and
  !! the criterion was measured to return a clean pass with one injected. The
  !! only safe response is to refuse the sweep outright.
  !!
  !! @param perturbations The sweep perturbations.
  !! @param errors The sweep errors.
  !! @param n Number of sweep points.
  subroutine fd_check_finite(perturbations, errors, n)
    real(kind=rp), intent(in) :: perturbations(:)
    real(kind=rp), intent(in) :: errors(:)
    integer, intent(in) :: n

    integer :: k
    character(len=32) :: index_str

    do k = 1, n
       if (ieee_is_nan(errors(k)) .or. ieee_is_nan(perturbations(k))) then
          write(index_str, '(I0)') k
          call neko_error('The finite-difference sweep contains a NaN at ' &
               // 'point ' // trim(index_str) // '. Every comparison ' // &
               'against a NaN is false, so the strict criterion would ' // &
               'quietly analyse a shorter sweep and could still report a ' // &
               'pass. Fix the sweep, not the criterion.')
       end if
       if (abs(errors(k)) .gt. 0.5_rp * huge(1.0_rp) .or. &
            abs(perturbations(k)) .gt. 0.5_rp * huge(1.0_rp)) then
          write(index_str, '(I0)') k
          call neko_error('The finite-difference sweep overflows at point ' &
               // trim(index_str) // '; the strict criterion cannot be ' // &
               'applied to it.')
       end if
    end do

  end subroutine fd_check_finite

  !> Require every perturbation of the sweep to share one sign.
  !!
  !! The magnitude model the criterion fits is the truncation series of *one*
  !! one-sided difference: the harness flips the sign of the step once and for
  !! all when the probed design value sits past the middle of its range, and a
  !! sweep mixing both signs mixes two different series.
  !!
  !! @param perturbations The sweep perturbations.
  !! @param n Number of sweep points.
  subroutine fd_check_single_sign(perturbations, n)
    real(kind=rp), intent(in) :: perturbations(:)
    integer, intent(in) :: n

    integer :: k

    do k = 1, n
       if (abs(perturbations(k)) .le. 0.0_rp) then
          call neko_error('The finite-difference sweep contains a zero ' // &
               'perturbation; there is no difference to take there.')
       end if
       if (perturbations(k) * perturbations(1) .lt. 0.0_rp) then
          call neko_error('The finite-difference sweep mixes positive and ' &
               // 'negative perturbations. The strict criterion fits one ' // &
               'one-sided truncation series and cannot span both.')
       end if
    end do

  end subroutine fd_check_single_sign

  !> Order the sweep by descending perturbation magnitude.
  !!
  !! A deterministic insertion sort: the sweep is short, and two runs of the
  !! same case must order an equal pair of magnitudes the same way.
  !!
  !! @param perturbations The sweep perturbations, in any order.
  !! @param errors The matching errors.
  !! @param n Number of sweep points.
  !! @param m Perturbation magnitudes, descending.
  !! @param e Errors in the same order.
  subroutine fd_sort_sweep(perturbations, errors, n, m, e)
    real(kind=rp), intent(in) :: perturbations(:)
    real(kind=rp), intent(in) :: errors(:)
    integer, intent(in) :: n
    real(kind=rp), intent(out) :: m(:)
    real(kind=rp), intent(out) :: e(:)

    real(kind=rp) :: key_m, key_e
    integer :: j, k

    do k = 1, n
       m(k) = abs(perturbations(k))
       e(k) = errors(k)
    end do

    do k = 2, n
       key_m = m(k)
       key_e = e(k)
       j = k - 1
       do while (j .ge. 1)
          if (m(j) .ge. key_m) exit
          m(j + 1) = m(j)
          e(j + 1) = e(j)
          j = j - 1
       end do
       m(j + 1) = key_m
       e(j + 1) = key_e
    end do

  end subroutine fd_sort_sweep

  !> Require every perturbation magnitude of the sweep to be distinct.
  !!
  !! Not a tidiness check: a repeated magnitude is the same measurement
  !! written down twice, and the plateau branch prices evidence as the spread
  !! \f$\max e - \min e\f$ over a window. Duplicates therefore lengthen a
  !! window at *zero* spread, so a sweep whose points have collapsed onto one
  !! bound-limited step satisfies `fd_plateau_points` from a single
  !! measurement and certifies a bound -- which is precisely the
  !! assert-on-one-point rule this module exists to abolish, re-entering from
  !! behind. The order branch is no safer: three copies of one point make two
  !! vanishing differences and an undefined ratio.
  !!
  !! Refusing the sweep, rather than silently de-duplicating it, is the same
  !! argument the NaN guard makes: the sweep that ran is not the sweep that
  !! was asked for, and the honest response is to say so.
  !!
  !! @param m Perturbation magnitudes, descending -- `fd_sort_sweep` must
  !!        have run first, so that equal magnitudes are adjacent.
  !! @param n Number of sweep points.
  subroutine fd_check_distinct(m, n)
    real(kind=rp), intent(in) :: m(:)
    integer, intent(in) :: n

    integer :: k
    character(len=32) :: value_str

    do k = 1, n - 1
       ! Descending, so the difference cannot be negative: a non-positive
       ! one is an exact repeat, tested without comparing reals for equality.
       if (m(k) - m(k + 1) .le. 0.0_rp) then
          write(value_str, '(E13.6E3)') m(k)
          call neko_error('The finite-difference sweep repeats the ' // &
               'perturbation magnitude ' // trim(adjustl(value_str)) // &
               '. A repeated magnitude is one measurement recorded twice: ' &
               // 'it lengthens a plateau window at zero spread, so the ' // &
               'criterion would certify a bound from a single point. The ' &
               // 'likely cause is the step being clamped against a ' // &
               'design bound -- look for the "FD sweep: requested step ' // &
               '... clamped to" lines -- which collapses several ' // &
               'requested points onto the same bound-limited value. Ask ' // &
               'for perturbations that fit inside the available headroom, ' &
               // 'or move the design away from its bounds.')
       end if
    end do

  end subroutine fd_check_distinct

  !> Fill the reported-only diagnostics of a verdict.
  !!
  !! @param m Perturbation magnitudes, descending.
  !! @param e Signed relative errors.
  !! @param n Number of sweep points.
  !! @param verdict The verdict to fill.
  subroutine fd_diagnostics(m, e, n, verdict)
    real(kind=rp), intent(in) :: m(:), e(:)
    integer, intent(in) :: n
    type(fd_verdict_t), intent(inout) :: verdict

    integer :: k, k_min

    k_min = 1
    do k = 2, n
       if (abs(e(k)) .lt. abs(e(k_min))) k_min = k
    end do
    verdict%min_abs_error = abs(e(k_min))
    verdict%min_perturbation = m(k_min)
    verdict%bracketed = (k_min .gt. 1) .and. (k_min .lt. n)

    verdict%n_sign_crossings = 0
    do k = 1, n - 1
       if (e(k) * e(k + 1) .lt. 0.0_rp) then
          verdict%n_sign_crossings = verdict%n_sign_crossings + 1
       end if
    end do

  end subroutine fd_diagnostics

  !> Longest window of consecutive points over which the signed error is flat.
  !!
  !! Flat means a spread of at most `plateau_fraction * tolerance`. Ties in
  !! length are broken towards the largest perturbations, which is where
  !! round-off contaminates the error least.
  !!
  !! @param e Signed relative errors, descending in perturbation.
  !! @param n Number of sweep points.
  !! @param tolerance Largest acceptable bias.
  !! @param plateau_fraction Accepted spread, as a fraction of the tolerance.
  !! @param found True if any window qualified.
  !! @param i0 First index of the window.
  !! @param i1 Last index of the window.
  !! @param level Mean signed error over the window.
  !! @param spread Spread of the signed error over the window.
  subroutine fd_bounded_window(e, n, tolerance, plateau_fraction, found, i0, &
       i1, level, spread)
    real(kind=rp), intent(in) :: e(:)
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: tolerance, plateau_fraction
    logical, intent(out) :: found
    integer, intent(out) :: i0, i1
    real(kind=rp), intent(out) :: level, spread

    real(kind=rp) :: window_spread, window_mean
    integer :: i, j

    found = .false.
    i0 = 0
    i1 = 0
    level = 0.0_rp
    spread = 0.0_rp

    do i = 1, n - fd_plateau_points + 1
       do j = i + fd_plateau_points - 1, n
          window_spread = maxval(e(i:j)) - minval(e(i:j))
          if (window_spread .gt. plateau_fraction * tolerance) cycle
          window_mean = sum(e(i:j)) / real(j - i + 1, rp)
          if (found .and. (j - i) .le. (i1 - i0)) cycle
          found = .true.
          i0 = i
          i1 = j
          level = window_mean
          spread = window_spread
       end do
    end do

  end subroutine fd_bounded_window

  !> Longest run of consecutive difference ratios showing the expected order.
  !!
  !! \f$d_k = e_k - e_{k+1}\f$ cancels the bias exactly, so the ratio
  !! \f$d_k/d_{k+1}\f$ determines the truncation order without needing
  !! \f$C\f$. A ratio is usable only when both differences are non-zero --
  !! `volume` has two of its seven differences exactly zero, so this is a
  !! real-data requirement -- and when they share a sign, tested as the
  !! *product* \f$d_k d_{k+1}\f$ and never as a quotient, so that a
  !! difference at the bottom of the exponent range cannot manufacture an
  !! overflow or a NaN where the test only needed a sign.
  !!
  !! @param m Perturbation magnitudes, descending.
  !! @param e Signed relative errors.
  !! @param n Number of sweep points.
  !! @param p_expected Expected truncation order.
  !! @param order_tol Half-width of the accepted band around `p_expected`.
  !! @param p_hat Per-ratio order estimates, zero where unusable.
  !! @param found True if a long enough run was found.
  !! @param k0 First ratio index of the run.
  !! @param k1 Last ratio index of the run.
  !! @param p_median Median order over the run.
  !! @param reason Why no run was found, when `found` is false.
  subroutine fd_order_run(m, e, n, p_expected, order_tol, p_hat, found, k0, &
       k1, p_median, reason)
    real(kind=rp), intent(in) :: m(:), e(:)
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: p_expected, order_tol
    real(kind=rp), intent(out) :: p_hat(:)
    logical, intent(out) :: found
    integer, intent(out) :: k0, k1
    real(kind=rp), intent(out) :: p_median
    character(len=*), intent(out) :: reason

    real(kind=rp) :: d(max(n - 1, 1))
    logical :: usable(max(n - 2, 1)), solved
    integer :: k, best_length, best_start, length, start

    found = .false.
    k0 = 0
    k1 = 0
    p_median = 0.0_rp
    reason = ''
    p_hat = 0.0_rp

    if (n .lt. 4) then
       reason = 'The order check needs at least 4 sweep points; this ' // &
            'sweep has fewer.'
       return
    end if

    do k = 1, n - 1
       d(k) = e(k) - e(k + 1)
    end do

    usable = .false.
    do k = 1, n - 2
       ! A vanished difference carries no order, and the ratio is undefined.
       if (abs(d(k)) .le. 0.0_rp .or. abs(d(k + 1)) .le. 0.0_rp) cycle
       ! Sign agreement as a product: a quotient of two differences at the
       ! bottom of the exponent range can overflow where the sign test cannot.
       if (d(k) * d(k + 1) .le. 0.0_rp) cycle
       call fd_solve_order(d(k) / d(k + 1), m(k), m(k + 1), m(k + 2), &
            p_hat(k), solved)
       if (.not. solved) cycle
       usable(k) = abs(p_hat(k) - p_expected) .le. order_tol
    end do

    best_length = 0
    best_start = 0
    length = 0
    start = 0
    do k = 1, n - 2
       if (usable(k)) then
          if (length .eq. 0) start = k
          length = length + 1
          if (length .gt. best_length) then
             best_length = length
             best_start = start
          end if
       else
          length = 0
       end if
    end do

    if (best_length .lt. fd_min_order_ratios) then
       write(reason, '(A,I0,A)') 'The longest usable order run is ', &
            best_length, ' ratio(s); at least 2 are needed (4 points).'
       return
    end if

    k0 = best_start
    k1 = best_start + best_length - 1
    p_median = fd_median(p_hat(k0:k1), best_length)

    ! Implied by the run condition, since every member of the run is already
    ! within `order_tol` of `p_expected`; asserted anyway because it is the
    ! statement the criterion actually makes.
    if (abs(p_median - p_expected) .gt. order_tol) then
       write(reason, '(A,F8.4,A,F6.2)') 'The median order ', p_median, &
            ' is outside the accepted band around ', p_expected
       return
    end if

    found = .true.

  end subroutine fd_order_run

  !> The difference ratio a truncation branch of order `p` would show over
  !! three consecutive magnitudes.
  !!
  !! \f$ R(p) = (m_0^p - m_1^p)/(m_1^p - m_2^p) \f$, continued at
  !! \f$p \le 0\f$ by its limit \f$\log(m_0/m_1)/\log(m_1/m_2)\f$ so that the
  !! bisection has a finite value at the bottom of its bracket.
  !!
  !! Neither divisor can vanish: `fd_check_distinct` has already refused any
  !! sweep with a repeated magnitude, so \f$m_0 > m_1 > m_2 > 0\f$ strictly,
  !! which makes \f$\log(m_1/m_2) > 0\f$ and \f$m_1^p > m_2^p\f$ for every
  !! \f$p > 0\f$. That guarantee is a property of this module rather than of
  !! its caller, so it holds for the recorded-sweep use too.
  !!
  !! @param p The order.
  !! @param m0 Largest of the three magnitudes.
  !! @param m1 Middle magnitude.
  !! @param m2 Smallest magnitude.
  !! @return The ratio.
  function fd_ratio_of_order(p, m0, m1, m2) result(ratio)
    real(kind=rp), intent(in) :: p, m0, m1, m2
    real(kind=rp) :: ratio

    if (p .le. 0.0_rp) then
       ratio = log(m0 / m1) / log(m1 / m2)
    else
       ratio = (m0**p - m1**p) / (m1**p - m2**p)
    end if

  end function fd_ratio_of_order

  !> Recover the truncation order from a measured difference ratio.
  !!
  !! \f$R(p)\f$ is strictly increasing in \f$p\f$, so the root is found by
  !! bisection with no derivative and no failure mode beyond the ratio simply
  !! lying outside the attainable range. This is a strict generalisation of
  !! the geometric-sweep closed form \f$\log R/\log r\f$, which it reproduces
  !! exactly when \f$m_0/m_1 = m_1/m_2\f$ -- and unlike that closed form it
  !! stays correct on the non-geometric sweep the drivers actually run.
  !!
  !! @param ratio The measured difference ratio.
  !! @param m0 Largest of the three magnitudes.
  !! @param m1 Middle magnitude.
  !! @param m2 Smallest magnitude.
  !! @param p The recovered order.
  !! @param solved False if the ratio is not attainable by any order in the
  !!        bracket, in which case the three points are not a truncation
  !!        branch at all.
  subroutine fd_solve_order(ratio, m0, m1, m2, p, solved)
    real(kind=rp), intent(in) :: ratio, m0, m1, m2
    real(kind=rp), intent(out) :: p
    logical, intent(out) :: solved

    real(kind=rp) :: lo, hi, mid
    integer :: it

    p = 0.0_rp
    solved = .false.
    if (ieee_is_nan(ratio)) return
    if (ratio .le. 0.0_rp) return
    if (ratio .le. fd_ratio_of_order(fd_order_lower, m0, m1, m2)) return
    if (ratio .ge. fd_ratio_of_order(fd_order_upper, m0, m1, m2)) return

    lo = fd_order_lower
    hi = fd_order_upper
    do it = 1, fd_order_bisections
       mid = 0.5_rp * (lo + hi)
       if (fd_ratio_of_order(mid, m0, m1, m2) .lt. ratio) then
          lo = mid
       else
          hi = mid
       end if
    end do
    p = 0.5_rp * (lo + hi)
    solved = .true.

  end subroutine fd_solve_order

  !> Median of a short list, averaging the middle pair for an even count.
  !! @param v The values.
  !! @param n Number of values.
  !! @return The median.
  function fd_median(v, n) result(median)
    real(kind=rp), intent(in) :: v(:)
    integer, intent(in) :: n
    real(kind=rp) :: median

    real(kind=rp) :: sorted(n), key
    integer :: j, k

    sorted(1:n) = v(1:n)
    do k = 2, n
       key = sorted(k)
       j = k - 1
       do while (j .ge. 1)
          if (sorted(j) .le. key) exit
          sorted(j + 1) = sorted(j)
          j = j - 1
       end do
       sorted(j + 1) = key
    end do

    if (mod(n, 2) .eq. 1) then
       median = sorted((n + 1) / 2)
    else
       median = 0.5_rp * (sorted(n / 2) + sorted(n / 2 + 1))
    end if

  end function fd_median

  !> Richardson extrapolation of the bias from two neighbouring points.
  !!
  !! \f$ C = e_{k+1} - (e_k - e_{k+1})/(r^p - 1) \f$ with \f$r\f$ the *local*
  !! magnitude ratio of the two points, taken over the two smallest points of
  !! the order run: the least contaminated by the next truncation term, and
  !! so the most faithful estimate of the bias the sweep can offer.
  !!
  !! @param m Perturbation magnitudes, descending.
  !! @param e Signed relative errors.
  !! @param k Index of the larger of the two points.
  !! @param p Measured truncation order.
  !! @param c_hat The extrapolated bias.
  !! @param ok False if the two points are too close in magnitude to
  !!        extrapolate between.
  !! @param ratio The magnitude ratio of the two points.
  subroutine fd_richardson(m, e, k, p, c_hat, ok, ratio)
    real(kind=rp), intent(in) :: m(:), e(:)
    integer, intent(in) :: k
    real(kind=rp), intent(in) :: p
    real(kind=rp), intent(out) :: c_hat
    logical, intent(out) :: ok
    real(kind=rp), intent(out) :: ratio

    c_hat = 0.0_rp
    ratio = m(k) / m(k + 1)
    ok = ratio .ge. fd_min_richardson_ratio
    if (.not. ok) return

    c_hat = e(k + 1) - (e(k) - e(k + 1)) / (ratio**p - 1.0_rp)

  end subroutine fd_richardson

  !> Three-term least-squares fit \f$e = C + A m^p + B m^{p+1}\f$.
  !!
  !! Used twice: as the cross-check on the Richardson bias, where it uses
  !! every point of the truncation window rather than the two smallest and
  !! carries the next truncation term explicitly; and to recover the
  !! truncation slope \f$A\f$ and the noise level behind the
  !! reachable-tolerance floor.
  !!
  !! Both columns are rescaled by their own maximum before the normal
  !! equations are formed, and the recovered coefficients are divided back
  !! out. Measured at \f$p = 2\f$ over several decades, that drops the
  !! condition number of the normal matrix from 1.4e5 to 1.4e1. The fit is
  !! deliberately unweighted: \f$1/x\f$ weighting spans twelve orders of
  !! magnitude over the same sweep, drives the condition number to 1e16 and
  !! was measured to produce a NaN bias on 5.6% of cases.
  !!
  !! @param m Perturbation magnitudes, descending.
  !! @param e Signed relative errors.
  !! @param i0 First point of the window.
  !! @param i1 Last point of the window.
  !! @param p Truncation order to fit at.
  !! @param c_hat The fitted bias.
  !! @param ok False if the window is too short or the fit is singular.
  !! @param a_hat (Optional) The fitted truncation slope.
  !! @param b_hat (Optional) The fitted next-term coefficient.
  subroutine fd_fit_three_term(m, e, i0, i1, p, c_hat, ok, a_hat, b_hat)
    real(kind=rp), intent(in) :: m(:), e(:)
    integer, intent(in) :: i0, i1
    real(kind=rp), intent(in) :: p
    real(kind=rp), intent(out) :: c_hat
    logical, intent(out) :: ok
    real(kind=rp), intent(out), optional :: a_hat, b_hat

    real(kind=rp) :: x1(i1 - i0 + 1), x2(i1 - i0 + 1)
    real(kind=rp) :: matrix(3, 3), rhs(3), solution(3), lambda(3)
    real(kind=rp) :: scale1, scale2, condition
    integer :: k, n_window

    c_hat = 0.0_rp
    ok = .false.
    if (present(a_hat)) a_hat = 0.0_rp
    if (present(b_hat)) b_hat = 0.0_rp
    n_window = i1 - i0 + 1
    if (n_window .lt. 4) return

    do k = 1, n_window
       x1(k) = m(i0 + k - 1)**p
       x2(k) = m(i0 + k - 1)**(p + 1.0_rp)
    end do
    scale1 = maxval(x1)
    scale2 = maxval(x2)
    if (scale1 .le. 0.0_rp .or. scale2 .le. 0.0_rp) return
    x1 = x1 / scale1
    x2 = x2 / scale2

    matrix(1, 1) = real(n_window, rp)
    matrix(1, 2) = sum(x1)
    matrix(1, 3) = sum(x2)
    matrix(2, 1) = matrix(1, 2)
    matrix(2, 2) = sum(x1 * x1)
    matrix(2, 3) = sum(x1 * x2)
    matrix(3, 1) = matrix(1, 3)
    matrix(3, 2) = matrix(2, 3)
    matrix(3, 3) = sum(x2 * x2)

    rhs(1) = sum(e(i0:i1))
    rhs(2) = sum(x1 * e(i0:i1))
    rhs(3) = sum(x2 * e(i0:i1))

    ! The normal matrix is symmetric positive semi-definite by construction,
    ! so a non-positive eigenvalue means the columns are degenerate.
    call fd_symmetric_eigenvalues(matrix, lambda)
    if (minval(lambda) .le. 0.0_rp) return
    condition = maxval(lambda) / minval(lambda)
    if (condition .ge. 1.0_rp / sqrt(NEKO_EPS)) return

    call fd_solve_3x3(matrix, rhs, solution, ok)
    if (.not. ok) return
    c_hat = solution(1)
    if (present(a_hat)) a_hat = solution(2) / scale1
    if (present(b_hat)) b_hat = solution(3) / scale2

  end subroutine fd_fit_three_term

  !> Smallest tolerance this functional's own reproducibility permits.
  !!
  !! Two constraints fight over the smallest useful perturbation: certifying
  !! a floor at the tolerance needs the truncation term below the plateau
  !! threshold, \f$|A| m \le f\,\mathrm{tol}\f$, while round-off needs
  !! \f$N/m \le 0.1\,\mathrm{tol}\f$. The two only leave a window open above
  !! \f$\mathrm{tol}_{\min} = \sqrt{40 N |A|}\f$. Below it no sweep of any
  !! depth can certify the gradient, so a red verdict there says nothing
  !! about the gradient and must not be reported as if it did.
  !!
  !! \f$A\f$ is fitted over the smallest `fd_slope_points` of the sweep,
  !! where the functional is least contaminated by its own nonlinearity, and
  !! \f$N\f$ is the median of \f$|residual| \times m\f$ over the smallest
  !! `fd_noise_points`, the residual being taken against the same three-term
  !! model.
  !!
  !! @param m Perturbation magnitudes, descending.
  !! @param e Signed relative errors.
  !! @param n Number of sweep points.
  !! @param p Expected truncation order.
  !! @param verdict The verdict to record the floor in.
  subroutine fd_reachable_tolerance(m, e, n, p, verdict)
    real(kind=rp), intent(in) :: m(:), e(:)
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: p
    type(fd_verdict_t), intent(inout) :: verdict

    real(kind=rp) :: c_slope, a_slope, b_slope, c_noise, a_noise, b_noise
    real(kind=rp) :: scaled(max(n, 1)), noise
    integer :: i0, i1, k, n_noise
    logical :: ok

    if (n .lt. fd_slope_points) return

    call fd_fit_three_term(m, e, n - fd_slope_points + 1, n, p, c_slope, &
         ok, a_slope, b_slope)
    if (.not. ok) return

    if (abs(a_slope) .le. fd_min_slope) then
       ! An exactly linear functional has no truncation branch, so no depth
       ! of sweep is ruled out and every tolerance is reachable.
       verdict%has_tol_min = .true.
       verdict%tol_min = 0.0_rp
       return
    end if

    i0 = max(1, n - fd_noise_points + 1)
    i1 = n
    n_noise = i1 - i0 + 1
    call fd_fit_three_term(m, e, i0, i1, p, c_noise, ok, a_noise, b_noise)
    if (.not. ok) return

    do k = 1, n_noise
       scaled(k) = abs(e(i0 + k - 1) - (c_noise &
            + a_noise * m(i0 + k - 1)**p &
            + b_noise * m(i0 + k - 1)**(p + 1.0_rp))) * m(i0 + k - 1)
    end do
    noise = fd_median(scaled(1:n_noise), n_noise)

    verdict%has_tol_min = .true.
    verdict%tol_min = sqrt(fd_reachability_constant * noise * abs(a_slope))

  end subroutine fd_reachable_tolerance

  !> Eigenvalues of a symmetric 3x3 matrix, in closed form.
  !!
  !! Used only to condition-check the normal matrix of the three-term fit,
  !! which is why a closed form is preferred to pulling in a solver: it is
  !! allocation-free and gives the same answer on every rank.
  !!
  !! @param a The symmetric matrix.
  !! @param lambda Its three eigenvalues.
  subroutine fd_symmetric_eigenvalues(a, lambda)
    real(kind=rp), intent(in) :: a(3, 3)
    real(kind=rp), intent(out) :: lambda(3)

    real(kind=rp), parameter :: pi = 4.0_rp * atan(1.0_rp)
    real(kind=rp) :: b(3, 3), p1, p2, q, p, r, phi
    integer :: k

    p1 = a(1, 2)**2 + a(1, 3)**2 + a(2, 3)**2
    if (p1 .le. 0.0_rp) then
       lambda(1) = a(1, 1)
       lambda(2) = a(2, 2)
       lambda(3) = a(3, 3)
       return
    end if

    q = (a(1, 1) + a(2, 2) + a(3, 3)) / 3.0_rp
    p2 = (a(1, 1) - q)**2 + (a(2, 2) - q)**2 + (a(3, 3) - q)**2 + 2.0_rp * p1
    p = sqrt(p2 / 6.0_rp)

    b = a / p
    do k = 1, 3
       b(k, k) = (a(k, k) - q) / p
    end do

    r = 0.5_rp * (b(1, 1) * (b(2, 2) * b(3, 3) - b(2, 3) * b(3, 2)) &
         - b(1, 2) * (b(2, 1) * b(3, 3) - b(2, 3) * b(3, 1)) &
         + b(1, 3) * (b(2, 1) * b(3, 2) - b(2, 2) * b(3, 1)))
    r = max(-1.0_rp, min(1.0_rp, r))
    phi = acos(r) / 3.0_rp

    lambda(1) = q + 2.0_rp * p * cos(phi)
    lambda(3) = q + 2.0_rp * p * cos(phi + 2.0_rp * pi / 3.0_rp)
    lambda(2) = 3.0_rp * q - lambda(1) - lambda(3)

  end subroutine fd_symmetric_eigenvalues

  !> Solve a 3x3 system by Gaussian elimination with partial pivoting.
  !! @param a The matrix.
  !! @param b The right-hand side.
  !! @param x The solution.
  !! @param ok False if the matrix is numerically singular.
  subroutine fd_solve_3x3(a, b, x, ok)
    real(kind=rp), intent(in) :: a(3, 3)
    real(kind=rp), intent(in) :: b(3)
    real(kind=rp), intent(out) :: x(3)
    logical, intent(out) :: ok

    real(kind=rp) :: w(3, 4), factor, row(4)
    integer :: i, j, k, pivot

    ok = .false.
    x = 0.0_rp
    w(:, 1:3) = a
    w(:, 4) = b

    do k = 1, 3
       pivot = k
       do i = k + 1, 3
          if (abs(w(i, k)) .gt. abs(w(pivot, k))) pivot = i
       end do
       if (abs(w(pivot, k)) .le. 0.0_rp) return
       if (pivot .ne. k) then
          row = w(k, :)
          w(k, :) = w(pivot, :)
          w(pivot, :) = row
       end if
       do i = k + 1, 3
          factor = w(i, k) / w(k, k)
          do j = k, 4
             w(i, j) = w(i, j) - factor * w(k, j)
          end do
       end do
    end do

    do i = 3, 1, -1
       x(i) = w(i, 4)
       do j = i + 1, 3
          x(i) = x(i) - w(i, j) * x(j)
       end do
       x(i) = x(i) / w(i, i)
    end do

    ok = .true.

  end subroutine fd_solve_3x3

  !> Compare the Richardson bias against the three-term fit.
  !!
  !! A disagreement says the truncation window is contaminated, which is worth
  !! knowing but is not itself evidence about the gradient -- so it warns and
  !! the assertion is still made on the Richardson estimate.
  !!
  !! @param verdict The verdict to annotate.
  !! @param c_rich The Richardson bias.
  !! @param c_three The three-term bias.
  subroutine fd_cross_check(verdict, c_rich, c_three)
    type(fd_verdict_t), intent(inout) :: verdict
    real(kind=rp), intent(in) :: c_rich, c_three

    character(len=fd_message_len) :: line

    if (.not. verdict%has_c_cross) then
       call fd_add_message(verdict, 'The 3-term cross-check could not be ' &
            // 'formed on this window; the')
       call fd_add_message(verdict, 'Richardson estimate is uncorroborated.')
       return
    end if
    if (abs(c_rich) .le. 0.0_rp) return
    if (abs(c_three) .le. fd_cross_check_factor * abs(c_rich) .and. &
         abs(c_three) * fd_cross_check_factor .ge. abs(c_rich)) return

    write(line, '(A,E13.6E3,A,E13.6E3)') 'WARNING: Richardson ', c_rich, &
         ' and the 3-term fit ', c_three
    call fd_add_message(verdict, line)
    call fd_add_message(verdict, 'WARNING: disagree by more than the ' // &
         'cross-check factor, so the truncation')
    call fd_add_message(verdict, 'WARNING: window is probably ' // &
         'contaminated. Asserting on Richardson.')

  end subroutine fd_cross_check

  !> Note, on a passing verdict, that the sweep never reached the round-off
  !! upturn. Informational: a monotone sweep is allowed to pass, because the
  !! bias it certifies is extrapolated rather than observed at a minimum.
  !! @param verdict The verdict to annotate.
  subroutine fd_report_upturn(verdict)
    type(fd_verdict_t), intent(inout) :: verdict

    if (verdict%bracketed) return
    call fd_add_message(verdict, 'No round-off upturn was observed in this ' &
         // 'sweep; the certified floor is')
    call fd_add_message(verdict, 'extrapolated, not sampled. Reported, not ' &
         // 'asserted on.')

  end subroutine fd_report_upturn

  !> Record the sign pattern of the signed error across the sweep, largest
  !! perturbation first. Printed whenever the order run breaks, because a
  !! broken run is almost always a sign pattern that is not a single
  !! truncation branch.
  !! @param e Signed relative errors.
  !! @param n Number of sweep points.
  !! @param verdict The verdict to annotate.
  subroutine fd_report_sign_pattern(e, n, verdict)
    real(kind=rp), intent(in) :: e(:)
    integer, intent(in) :: n
    type(fd_verdict_t), intent(inout) :: verdict

    character(len=fd_message_len) :: line
    integer :: k, n_shown

    line = 'Sign pattern of the error, largest perturbation first: '
    n_shown = min(n, fd_message_len - len_trim(line) - 1)
    do k = 1, n_shown
       if (e(k) .gt. 0.0_rp) then
          line = trim(line) // '+'
       else if (e(k) .lt. 0.0_rp) then
          line = trim(line) // '-'
       else
          line = trim(line) // '0'
       end if
    end do
    call fd_add_message(verdict, line)

  end subroutine fd_report_sign_pattern

  !> Say so when the whole sweep sits in the round-off regime -- the error
  !! growing as the perturbation shrinks, everywhere -- since that is an
  !! inadequate sweep for a reason the caller can act on directly.
  !! @param e Signed relative errors.
  !! @param n Number of sweep points.
  !! @param verdict The verdict to annotate.
  subroutine fd_report_roundoff(e, n, verdict)
    real(kind=rp), intent(in) :: e(:)
    integer, intent(in) :: n
    type(fd_verdict_t), intent(inout) :: verdict

    integer :: k
    logical :: growing

    growing = n .ge. 2
    do k = 1, n - 1
       if (abs(e(k + 1)) .le. abs(e(k))) growing = .false.
    end do
    if (.not. growing) return

    call fd_add_message(verdict, 'Every point of this sweep is round-off ' &
         // 'dominated: the error grows at')
    call fd_add_message(verdict, 'every step towards smaller ' // &
         'perturbations. Start the sweep larger.')

  end subroutine fd_report_roundoff

  !> Append a diagnostic line to a verdict, silently dropping anything past
  !! the last slot rather than overwriting an earlier, more specific line.
  !! @param verdict The verdict to annotate.
  !! @param text The line to append.
  subroutine fd_add_message(verdict, text)
    type(fd_verdict_t), intent(inout) :: verdict
    character(len=*), intent(in) :: text

    if (verdict%n_messages .ge. fd_max_messages) return
    verdict%n_messages = verdict%n_messages + 1
    verdict%message(verdict%n_messages) = text

  end subroutine fd_add_message

  !> Name of a branch, for the log and the verdict CSV.
  !! @param branch The branch code.
  !! @return Its name.
  function fd_branch_name(branch) result(name)
    integer, intent(in) :: branch
    character(len=8) :: name

    select case (branch)
    case (FD_BRANCH_BOUNDED)
       name = 'BOUNDED'
    case (FD_BRANCH_FIT)
       name = 'FIT'
    case default
       name = 'NONE'
    end select

  end function fd_branch_name

  !> Name of a status, for the log and the verdict CSV.
  !! @param status The status code.
  !! @return Its name.
  function fd_status_name(status) result(name)
    integer, intent(in) :: status
    character(len=22) :: name

    select case (status)
    case (FD_STATUS_OK)
       name = 'OK'
    case (FD_STATUS_FLOOR_EXCEEDED)
       name = 'FLOOR_EXCEEDED'
    case (FD_STATUS_DEGENERATE)
       name = 'DEGENERATE_SENSITIVITY'
    case (FD_STATUS_UNREACHABLE)
       name = 'TOLERANCE_UNREACHABLE'
    case default
       name = 'INCONCLUSIVE'
    end select

  end function fd_status_name

  !> Whether the reported bias is an estimate or only a bound. The two are
  !! different claims and are never to be read as the same number.
  !! @param verdict The verdict.
  !! @return `estimate`, `bound`, or `none`.
  function fd_c_hat_kind(verdict) result(kind_name)
    type(fd_verdict_t), intent(in) :: verdict
    character(len=8) :: kind_name

    if (.not. verdict%has_c_hat) then
       kind_name = 'none'
    else if (verdict%branch .eq. FD_BRANCH_FIT) then
       kind_name = 'estimate'
    else
       kind_name = 'bound'
    end if

  end function fd_c_hat_kind

  !> Print a verdict and its diagnostics.
  !!
  !! Prints unconditionally: the caller decides which rank owns the log, as
  !! the rest of the harness does.
  !!
  !! @param verdict The verdict to print.
  !! @param tolerance The tolerance it was judged against.
  subroutine fd_verdict_print(verdict, tolerance)
    type(fd_verdict_t), intent(in) :: verdict
    real(kind=rp), intent(in) :: tolerance

    character(len=*), parameter :: tag = ' FD strict: '
    integer :: k

    write(*, '(A,A,A,A)') tag, 'branch = ', trim(fd_branch_name( &
         verdict%branch)), ', status = ' // trim(fd_status_name( &
         verdict%status))
    if (verdict%has_c_hat) then
       write(*, '(A,A,A,E15.6E3,A,E15.6E3)') tag, &
            trim(fd_c_hat_kind(verdict)), ' of the bias = ', verdict%c_hat, &
            '   tolerance = ', tolerance
    else
       write(*, '(A,A,E15.6E3)') tag, 'bias not estimated; tolerance = ', &
            tolerance
    end if
    if (verdict%has_p_hat) then
       write(*, '(A,A,F8.4,A,I0,A)') tag, 'measured order = ', &
            verdict%p_hat, ' over ', verdict%n_truncation_points, &
            ' consecutive points'
    end if
    if (verdict%has_c_cross) then
       write(*, '(A,A,E15.6E3)') tag, '3-term cross-check bias = ', &
            verdict%c_cross
    end if
    if (verdict%has_window) then
       write(*, '(A,A,E15.6E3,A,L1)') tag, 'plateau bound on the bias = ', &
            verdict%window_bound, ', certifies = ', verdict%bounded_certified
    end if
    if (verdict%has_tol_min) then
       write(*, '(A,A,E15.6E3)') tag, 'smallest reachable tolerance for ' // &
            'this functional = ', verdict%tol_min
    end if
    write(*, '(A,A,E15.6E3,A,E15.6E3)') tag, 'minimum |error| = ', &
         verdict%min_abs_error, ' at perturbation ', verdict%min_perturbation
    write(*, '(A,A,L1,A,I0)') tag, 'upturn bracketed (diagnostic only) = ', &
         verdict%bracketed, ', sign crossings = ', verdict%n_sign_crossings

    do k = 1, verdict%n_messages
       write(*, '(A,A)') tag, trim(verdict%message(k))
    end do

  end subroutine fd_verdict_print

end module fd_criterion
