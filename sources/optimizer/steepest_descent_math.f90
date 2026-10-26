!> @file steepest_descent_math.f90
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

!> Vector-level kernels of the steepest-descent optimizer.
!! @details
!! Elementwise kernels of the steepest descent method with move limits,
!! Appendix D of da Silva, Aage, Beck & Sigmund (2021), Eqs. D22-D24. Each
!! routine operates on `vector_t` objects and runs on the device backend
!! when one is configured, and on the CPU otherwise. All three are local
!! operations; the global reduction of the direction is left to the caller.
module steepest_descent_math
  use num_types, only: rp
  use vector, only: vector_t
  use neko_config, only: NEKO_BCKND_DEVICE
  use steepest_descent_cpu, only: sd_direction_cpu, &
       sd_update_move_limits_cpu, sd_update_design_cpu
  use device_steepest_descent, only: device_sd_direction, &
       device_sd_update_move_limits, device_sd_update_design
  implicit none
  private

  public :: sd_direction, sd_update_move_limits, sd_update_design

contains

  !> @brief Bound-projected steepest-descent direction (Eq. D22).
  !! @details On entry `s` holds the gradient G. On exit it holds S = -G,
  !! with S_e = 0 where x_e >= xmax and S_e > 0, or x_e <= xmin and S_e < 0.
  !! @param[inout] s Gradient on entry, projected direction on exit.
  !! @param[in] x Design variables.
  !! @param[in] xmin Lower bound of the design variables.
  !! @param[in] xmax Upper bound of the design variables.
  !! @param[in] n Number of entries to process, the size of `s` by default.
  subroutine sd_direction(s, x, xmin, xmax, n)
    type(vector_t), intent(inout) :: s
    type(vector_t), intent(in) :: x
    real(kind=rp), intent(in) :: xmin, xmax
    integer, intent(in), optional :: n
    integer :: n_local

    if (present(n)) then
       n_local = n
    else
       n_local = s%size()
    end if

    ! An empty range has nothing to do, and its host array may be unallocated
    if (n_local .lt. 1) return

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_sd_direction(s%x_d, x%x_d, xmin, xmax, n_local)
    else
       call sd_direction_cpu(s%x, x%x, xmin, xmax, n_local)
    end if

  end subroutine sd_direction

  !> @brief Adapt the move limits (Eq. D23).
  !! @details With d_e = (x_e - x1_e) * (x1_e - x2_e), the move limit is
  !! scaled by `dec` where d_e < 0, by `inc` where d_e > 0, left unchanged
  !! where d_e = 0, and then clamped to [dmin, dmax].
  !! @param[inout] delta Move limits.
  !! @param[in] x Design variables at the current iteration.
  !! @param[in] x1 Design variables one iteration back.
  !! @param[in] x2 Design variables two iterations back.
  !! @param[in] dec Factor applied on oscillation.
  !! @param[in] inc Factor applied on monotone progress.
  !! @param[in] dmin Smallest allowed move limit.
  !! @param[in] dmax Largest allowed move limit.
  !! @param[in] n Number of entries to process, the size of `delta` by
  !! default.
  subroutine sd_update_move_limits(delta, x, x1, x2, dec, inc, dmin, dmax, n)
    type(vector_t), intent(inout) :: delta
    type(vector_t), intent(in) :: x, x1, x2
    real(kind=rp), intent(in) :: dec, inc, dmin, dmax
    integer, intent(in), optional :: n
    integer :: n_local

    if (present(n)) then
       n_local = n
    else
       n_local = delta%size()
    end if

    ! An empty range has nothing to do, and its host array may be unallocated
    if (n_local .lt. 1) return

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_sd_update_move_limits(delta%x_d, x%x_d, x1%x_d, x2%x_d, &
            dec, inc, dmin, dmax, n_local)
    else
       call sd_update_move_limits_cpu(delta%x, x%x, x1%x, x2%x, &
            dec, inc, dmin, dmax, n_local)
    end if

  end subroutine sd_update_move_limits

  !> @brief Move-limited, box-constrained design update (Eq. D24).
  !! @details Computes
  !! x_new = min(xmax, max(xmin, min(x + delta, max(x - delta,
  !! x + psi * (s / smax))))). This equals Eq. D24 whenever x lies in
  !! [xmin, xmax], and stays well defined for infeasible x.
  !! @param[inout] x_new Updated design variables.
  !! @param[in] x Current design variables.
  !! @param[in] s Bound-projected steepest-descent direction.
  !! @param[in] delta Move limits.
  !! @param[in] psi Step length.
  !! @param[in] smax Global maximum of the absolute value of `s`, which must
  !! be positive.
  !! @param[in] xmin Lower bound of the design variables.
  !! @param[in] xmax Upper bound of the design variables.
  !! @param[in] n Number of entries to process, the size of `x_new` by
  !! default.
  subroutine sd_update_design(x_new, x, s, delta, psi, smax, xmin, xmax, n)
    type(vector_t), intent(inout) :: x_new
    type(vector_t), intent(in) :: x, s, delta
    real(kind=rp), intent(in) :: psi, smax, xmin, xmax
    integer, intent(in), optional :: n
    integer :: n_local

    if (present(n)) then
       n_local = n
    else
       n_local = x_new%size()
    end if

    ! An empty range has nothing to do, and its host array may be unallocated
    if (n_local .lt. 1) return

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_sd_update_design(x_new%x_d, x%x_d, s%x_d, delta%x_d, &
            psi, smax, xmin, xmax, n_local)
    else
       call sd_update_design_cpu(x_new%x, x%x, s%x, delta%x, &
            psi, smax, xmin, xmax, n_local)
    end if

  end subroutine sd_update_design

end module steepest_descent_math
