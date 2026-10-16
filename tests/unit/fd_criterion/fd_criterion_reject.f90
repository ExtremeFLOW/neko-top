!> Standalone probe for the four strict-FD guards that reject a sweep or an
!! option set by calling `neko_error`, which is an `error stop` and so
!! cannot be caught in-process by pFUnit -- see fd_criterion_test.pf's
!! header and mapping_finite_difference.pf's test_convex_up_detection for
!! why.
!!
!! Each scenario constructs an input that must be refused. If the guard is
!! present the program aborts (non-zero exit, "ERROR STOP"); if it is ever
!! removed or weakened the program instead reaches the final `print` and
!! exits 0. CMakeLists.txt wires each scenario up as a plain `add_test`
!! with `WILL_FAIL TRUE`, so the ctest result is exactly backwards from a
!! normal test: exit 0 here is the failure.
!!
!! @param 1 argument: which scenario to run (1-4, see the four `case`
!!        branches below).
program fd_criterion_reject
  use num_types, only: rp
  use fd_criterion, only: fd_strict_options_t, fd_verdict_t, fd_evaluate, &
       fd_status_name
  use sensitivity, only: fd_read_strict_options
  use json_module, only: json_file
  use, intrinsic :: ieee_arithmetic, only: ieee_value, IEEE_QUIET_NAN
  implicit none

  type(fd_strict_options_t) :: opts
  type(fd_verdict_t) :: v
  type(json_file) :: params
  real(kind=rp) :: m(9), e(9)
  character(len=32) :: arg
  integer :: scenario

  call get_command_argument(1, arg)
  read(arg, '(I32)') scenario

  select case (scenario)
  case (1)
     ! A sweep repeating a perturbation magnitude -- the demonstrating
     ! input from the module header: 9 array entries, but only 3 distinct
     ! magnitudes (1e-4 repeated 7x, plus 3.16e-5 and 1e-5). Before the
     ! fix, this returned passed=T status=OK branch=BOUNDED, claiming "8
     ! consecutive points" from what are really two real measurements.
     m = [-1.0e-4_rp, -1.0e-4_rp, -1.0e-4_rp, -1.0e-4_rp, -1.0e-4_rp, &
          -1.0e-4_rp, -1.0e-4_rp, -3.16e-5_rp, -1.0e-5_rp]
     e = [1.0e-3_rp, 1.0e-3_rp, 1.0e-3_rp, 1.0e-3_rp, 1.0e-3_rp, &
          1.0e-3_rp, 1.0e-3_rp, 3.0e-3_rp, 9.0e-3_rp]
     call fd_evaluate(m, e, 1.0e-2_rp, opts, v)

  case (2)
     ! A NaN error entry. Every comparison against a NaN is false, so
     ! without an explicit ieee_is_nan check the analysis would simply run
     ! on a shorter sweep and could still report a pass.
     m(1:4) = [-1.0e-1_rp, -1.0e-2_rp, -1.0e-3_rp, -1.0e-4_rp]
     e(1:4) = [-0.11926740394713609e-1_rp, &
          ieee_value(1.0_rp, IEEE_QUIET_NAN), &
          -0.10733324016689514e-3_rp, 0.23047656070843778e-5_rp]
     call fd_evaluate(m(1:4), e(1:4), 1.0e-3_rp, opts, v)

  case (3)
     ! Mixed-sign perturbations: the model is the truncation series of one
     ! one-sided difference and cannot span both signs.
     m(1:4) = [-1.0e-1_rp, -1.0e-2_rp, 1.0e-3_rp, -1.0e-4_rp]
     e(1:4) = [-0.11926740394713609e-1_rp, -0.12022146350578529e-2_rp, &
          -0.10733324016689514e-3_rp, 0.23047656070843778e-5_rp]
     call fd_evaluate(m(1:4), e(1:4), 1.0e-3_rp, opts, v)

  case (4)
     ! optimization.fd_test_order_tolerance as wide as, or wider than,
     ! optimization.fd_test_order: a band at least as wide as the expected
     ! order accepts a measured order of zero, which is not a truncation
     ! branch at all and makes the Richardson extrapolation singular.
     call params%initialize()
     call params%load_from_string('{"optimization":{"fd_test_order":' // &
          '1.0,"fd_test_order_tolerance":1.0}}')
     call fd_read_strict_options(params, .false., opts)

  case default
     print *, 'usage: fd_criterion_reject <1-4>'
     call exit(2)
  end select

  print *, 'BUG: scenario ', scenario, ' was not rejected; verdict status = ' &
       // trim(fd_status_name(v%status))

end program fd_criterion_reject
