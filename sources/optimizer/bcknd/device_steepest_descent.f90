!> @file device_steepest_descent.f90
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

!> Device backend dispatch for the steepest-descent optimizer kernels.
!! @details
!! Binds the CUDA and HIP launch wrappers in `bcknd/device/{cuda,hip}` and
!! exposes one backend-independent wrapper per kernel. The formulas are those
!! of `steepest_descent_cpu`.
module device_steepest_descent
  use utils, only: neko_error
  use num_types, only: rp, c_rp
  use, intrinsic :: iso_c_binding, only: c_ptr, c_int
  implicit none
  private

  public :: device_sd_direction, device_sd_update_move_limits, &
       device_sd_update_design

#if HAVE_HIP

  interface
     subroutine hip_sd_direction(s_d, x_d, xmin, xmax, n) &
          bind(c, name = 'hip_sd_direction')
       import c_rp, c_ptr, c_int
       type(c_ptr), value :: s_d
       type(c_ptr), value :: x_d
       real(c_rp) :: xmin
       real(c_rp) :: xmax
       integer(c_int) :: n
     end subroutine hip_sd_direction
  end interface

  interface
     subroutine hip_sd_update_move_limits(delta_d, x_d, x1_d, x2_d, &
          dec, inc, dmin, dmax, n) &
          bind(c, name = 'hip_sd_update_move_limits')
       import c_rp, c_ptr, c_int
       type(c_ptr), value :: delta_d
       type(c_ptr), value :: x_d
       type(c_ptr), value :: x1_d
       type(c_ptr), value :: x2_d
       real(c_rp) :: dec
       real(c_rp) :: inc
       real(c_rp) :: dmin
       real(c_rp) :: dmax
       integer(c_int) :: n
     end subroutine hip_sd_update_move_limits
  end interface

  interface
     subroutine hip_sd_update_design(x_new_d, x_d, s_d, delta_d, &
          psi, smax, xmin, xmax, n) &
          bind(c, name = 'hip_sd_update_design')
       import c_rp, c_ptr, c_int
       type(c_ptr), value :: x_new_d
       type(c_ptr), value :: x_d
       type(c_ptr), value :: s_d
       type(c_ptr), value :: delta_d
       real(c_rp) :: psi
       real(c_rp) :: smax
       real(c_rp) :: xmin
       real(c_rp) :: xmax
       integer(c_int) :: n
     end subroutine hip_sd_update_design
  end interface

#elif HAVE_CUDA

  interface
     subroutine cuda_sd_direction(s_d, x_d, xmin, xmax, n) &
          bind(c, name = 'cuda_sd_direction')
       import c_rp, c_ptr, c_int
       type(c_ptr), value :: s_d
       type(c_ptr), value :: x_d
       real(c_rp) :: xmin
       real(c_rp) :: xmax
       integer(c_int) :: n
     end subroutine cuda_sd_direction
  end interface

  interface
     subroutine cuda_sd_update_move_limits(delta_d, x_d, x1_d, x2_d, &
          dec, inc, dmin, dmax, n) &
          bind(c, name = 'cuda_sd_update_move_limits')
       import c_rp, c_ptr, c_int
       type(c_ptr), value :: delta_d
       type(c_ptr), value :: x_d
       type(c_ptr), value :: x1_d
       type(c_ptr), value :: x2_d
       real(c_rp) :: dec
       real(c_rp) :: inc
       real(c_rp) :: dmin
       real(c_rp) :: dmax
       integer(c_int) :: n
     end subroutine cuda_sd_update_move_limits
  end interface

  interface
     subroutine cuda_sd_update_design(x_new_d, x_d, s_d, delta_d, &
          psi, smax, xmin, xmax, n) &
          bind(c, name = 'cuda_sd_update_design')
       import c_rp, c_ptr, c_int
       type(c_ptr), value :: x_new_d
       type(c_ptr), value :: x_d
       type(c_ptr), value :: s_d
       type(c_ptr), value :: delta_d
       real(c_rp) :: psi
       real(c_rp) :: smax
       real(c_rp) :: xmin
       real(c_rp) :: xmax
       integer(c_int) :: n
     end subroutine cuda_sd_update_design
  end interface

#endif

contains

  !> @brief Bound-projected steepest-descent direction on the device
  !! (Eq. D22).
  !! @param s_d Device pointer to the gradient on entry, and to the projected
  !! direction on exit.
  !! @param x_d Device pointer to the design variables.
  !! @param[in] xmin Lower bound of the design variables.
  !! @param[in] xmax Upper bound of the design variables.
  !! @param[in] n Number of design variables.
  subroutine device_sd_direction(s_d, x_d, xmin, xmax, n)
    type(c_ptr) :: s_d
    type(c_ptr) :: x_d
    real(kind=rp), intent(in) :: xmin
    real(kind=rp), intent(in) :: xmax
    integer, intent(in) :: n
#if HAVE_HIP
    call hip_sd_direction(s_d, x_d, xmin, xmax, n)
#elif HAVE_CUDA
    call cuda_sd_direction(s_d, x_d, xmin, xmax, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_sd_direction

  !> @brief Adapt the move limits on the device (Eq. D23).
  !! @param delta_d Device pointer to the move limits, updated in place.
  !! @param x_d Device pointer to the design variables at the current
  !! iteration.
  !! @param x1_d Device pointer to the design variables one iteration back.
  !! @param x2_d Device pointer to the design variables two iterations back.
  !! @param[in] dec Factor applied on oscillation.
  !! @param[in] inc Factor applied on monotone progress.
  !! @param[in] dmin Smallest allowed move limit.
  !! @param[in] dmax Largest allowed move limit.
  !! @param[in] n Number of design variables.
  subroutine device_sd_update_move_limits(delta_d, x_d, x1_d, x2_d, &
       dec, inc, dmin, dmax, n)
    type(c_ptr) :: delta_d
    type(c_ptr) :: x_d
    type(c_ptr) :: x1_d
    type(c_ptr) :: x2_d
    real(kind=rp), intent(in) :: dec
    real(kind=rp), intent(in) :: inc
    real(kind=rp), intent(in) :: dmin
    real(kind=rp), intent(in) :: dmax
    integer, intent(in) :: n
#if HAVE_HIP
    call hip_sd_update_move_limits(delta_d, x_d, x1_d, x2_d, &
         dec, inc, dmin, dmax, n)
#elif HAVE_CUDA
    call cuda_sd_update_move_limits(delta_d, x_d, x1_d, x2_d, &
         dec, inc, dmin, dmax, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_sd_update_move_limits

  !> @brief Move-limited, box-constrained design update on the device
  !! (Eq. D24).
  !! @param x_new_d Device pointer to the updated design variables.
  !! @param x_d Device pointer to the current design variables.
  !! @param s_d Device pointer to the bound-projected direction.
  !! @param delta_d Device pointer to the move limits.
  !! @param[in] psi Step length.
  !! @param[in] smax Global maximum of the absolute value of the direction.
  !! @param[in] xmin Lower bound of the design variables.
  !! @param[in] xmax Upper bound of the design variables.
  !! @param[in] n Number of design variables.
  subroutine device_sd_update_design(x_new_d, x_d, s_d, delta_d, &
       psi, smax, xmin, xmax, n)
    type(c_ptr) :: x_new_d
    type(c_ptr) :: x_d
    type(c_ptr) :: s_d
    type(c_ptr) :: delta_d
    real(kind=rp), intent(in) :: psi
    real(kind=rp), intent(in) :: smax
    real(kind=rp), intent(in) :: xmin
    real(kind=rp), intent(in) :: xmax
    integer, intent(in) :: n
#if HAVE_HIP
    call hip_sd_update_design(x_new_d, x_d, s_d, delta_d, &
         psi, smax, xmin, xmax, n)
#elif HAVE_CUDA
    call cuda_sd_update_design(x_new_d, x_d, s_d, delta_d, &
         psi, smax, xmin, xmax, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_sd_update_design

end module device_steepest_descent
