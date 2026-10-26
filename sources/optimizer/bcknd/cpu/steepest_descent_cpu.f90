!> @file steepest_descent_cpu.f90
!! @copyright
!! Copyright (c) 2026, The Neko-TOP Authors
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

!> CPU backend for the steepest-descent optimizer kernels.
!! @details
!! Elementwise kernels of the steepest descent method with move limits,
!! Appendix D of da Silva, Aage, Beck & Sigmund (2021), Eqs. D22-D24. The
!! device backends in `bcknd/device` implement the same formulas.
module steepest_descent_cpu
  use num_types, only: rp
  implicit none
  private

  public :: sd_direction_cpu, sd_update_move_limits_cpu, sd_update_design_cpu

contains

  !> @brief Bound-projected steepest-descent direction on CPU (Eq. D22).
  !! @details On entry `s` holds the gradient G. On exit it holds S = -G,
  !! with every component that points out of the box at an active bound set
  !! to zero.
  !! @param[inout] s Gradient on entry, projected direction on exit.
  !! @param[in] x Design variables.
  !! @param[in] xmin Lower bound of the design variables.
  !! @param[in] xmax Upper bound of the design variables.
  !! @param[in] n Number of design variables.
  subroutine sd_direction_cpu(s, x, xmin, xmax, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: s
    real(kind=rp), dimension(n), intent(in) :: x
    real(kind=rp), intent(in) :: xmin, xmax
    integer :: i

    do i = 1, n
       s(i) = -s(i)
       if ((x(i) .ge. xmax .and. s(i) .gt. 0.0_rp) .or. &
            (x(i) .le. xmin .and. s(i) .lt. 0.0_rp)) then
          s(i) = 0.0_rp
       end if
    end do

  end subroutine sd_direction_cpu

  !> @brief Adapt the move limits on CPU (Eq. D23).
  !! @details The move limit of a variable shrinks when its last two changes
  !! have opposite signs, grows when they have the same sign, and is then
  !! clamped to [dmin, dmax].
  !! @param[inout] delta Move limits.
  !! @param[in] x Design variables at the current iteration.
  !! @param[in] x1 Design variables one iteration back.
  !! @param[in] x2 Design variables two iterations back.
  !! @param[in] dec Factor applied on oscillation.
  !! @param[in] inc Factor applied on monotone progress.
  !! @param[in] dmin Smallest allowed move limit.
  !! @param[in] dmax Largest allowed move limit.
  !! @param[in] n Number of design variables.
  subroutine sd_update_move_limits_cpu(delta, x, x1, x2, dec, inc, dmin, &
       dmax, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: delta
    real(kind=rp), dimension(n), intent(in) :: x, x1, x2
    real(kind=rp), intent(in) :: dec, inc, dmin, dmax
    real(kind=rp) :: d
    integer :: i

    do i = 1, n
       d = (x(i) - x1(i)) * (x1(i) - x2(i))
       if (d .lt. 0.0_rp) then
          delta(i) = dec * delta(i)
       else if (d .gt. 0.0_rp) then
          delta(i) = inc * delta(i)
       end if
       delta(i) = min(dmax, max(dmin, delta(i)))
    end do

  end subroutine sd_update_move_limits_cpu

  !> @brief Move-limited, box-constrained design update on CPU (Eq. D24).
  !! @details Computes
  !! x_new = min(xmax, max(xmin, min(x + delta, max(x - delta,
  !! x + psi * (s / smax))))), i.e. the move limit first and the box second.
  !! This equals Eq. D24 whenever x lies in [xmin, xmax], and stays well
  !! defined for infeasible x.
  !! @param[out] x_new Updated design variables.
  !! @param[in] x Current design variables.
  !! @param[in] s Bound-projected steepest-descent direction.
  !! @param[in] delta Move limits.
  !! @param[in] psi Step length.
  !! @param[in] smax Global maximum of the absolute value of `s`.
  !! @param[in] xmin Lower bound of the design variables.
  !! @param[in] xmax Upper bound of the design variables.
  !! @param[in] n Number of design variables.
  subroutine sd_update_design_cpu(x_new, x, s, delta, psi, smax, xmin, xmax, &
       n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(out) :: x_new
    real(kind=rp), dimension(n), intent(in) :: x, s, delta
    real(kind=rp), intent(in) :: psi, smax, xmin, xmax
    integer :: i

    do i = 1, n
       x_new(i) = min(xmax, max(xmin, min(x(i) + delta(i), &
            max(x(i) - delta(i), x(i) + psi * (s(i) / smax)))))
    end do

  end subroutine sd_update_design_cpu

end module steepest_descent_cpu
