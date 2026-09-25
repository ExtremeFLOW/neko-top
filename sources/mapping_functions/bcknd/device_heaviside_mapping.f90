!> @file device_heaviside_mapping.f90
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
!
module device_heaviside_mapping
  use utils, only: neko_error
  use num_types, only: rp, c_rp
  use, intrinsic :: iso_c_binding, only: c_ptr, c_int
  implicit none
  private

  public :: device_heaviside_mapping_apply, &
       device_heaviside_mapping_apply_backward

#if HAVE_HIP
  interface
     !> HIP kernel launcher for forward Heaviside mapping.
     subroutine hip_heaviside_mapping_apply(beta, eta, X_out_d, &
          X_in_d, n) bind(c, name = 'hip_heaviside_mapping_apply')
       import c_rp, c_ptr, c_int
       real(c_rp) :: beta
       real(c_rp) :: eta
       type(c_ptr), value :: X_out_d
       type(c_ptr), value :: X_in_d
       integer(c_int) :: n
     end subroutine hip_heaviside_mapping_apply
  end interface

  interface
     !> HIP kernel launcher for Heaviside mapping chain rule.
     subroutine hip_heaviside_mapping_apply_backward(beta, eta, &
          sens_out_d, sens_in_d, X_in_d, n) &
          bind(c, name = 'hip_heaviside_mapping_apply_backward')
       import c_rp, c_ptr, c_int
       real(c_rp) :: beta
       real(c_rp) :: eta
       type(c_ptr), value :: sens_out_d
       type(c_ptr), value :: sens_in_d
       type(c_ptr), value :: X_in_d
       integer(c_int) :: n
     end subroutine hip_heaviside_mapping_apply_backward
  end interface
#elif HAVE_CUDA
  interface
     !> CUDA kernel launcher for forward Heaviside mapping.
     subroutine cuda_heaviside_mapping_apply(beta, eta, X_out_d, &
          X_in_d, n) bind(c, name = 'cuda_heaviside_mapping_apply')
       import c_rp, c_ptr, c_int
       real(c_rp) :: beta
       real(c_rp) :: eta
       type(c_ptr), value :: X_out_d
       type(c_ptr), value :: X_in_d
       integer(c_int) :: n
     end subroutine cuda_heaviside_mapping_apply
  end interface

  interface
     !> CUDA kernel launcher for Heaviside mapping chain rule.
     subroutine cuda_heaviside_mapping_apply_backward(beta, eta, &
          sens_out_d, sens_in_d, X_in_d, n) &
          bind(c, name = 'cuda_heaviside_mapping_apply_backward')
       import c_rp, c_ptr, c_int
       real(c_rp) :: beta
       real(c_rp) :: eta
       type(c_ptr), value :: sens_out_d
       type(c_ptr), value :: sens_in_d
       type(c_ptr), value :: X_in_d
       integer(c_int) :: n
     end subroutine cuda_heaviside_mapping_apply_backward
  end interface
#elif HAVE_OPENCL
#endif

contains

  !> Dispatch forward Heaviside mapping to active device backend.
  !! @param beta Projection sharpness parameter.
  !! @param eta Projection threshold parameter.
  !! @param X_out_d Output field on device.
  !! @param X_in_d Input field on device.
  !! @param n Number of dofs.
  subroutine device_heaviside_mapping_apply(beta, eta, X_out_d, &
       X_in_d, n)
    real(kind=rp), intent(in) :: beta
    real(kind=rp), intent(in) :: eta
    type(c_ptr) :: X_out_d
    type(c_ptr) :: X_in_d
    integer :: n
#if HAVE_HIP
    call hip_heaviside_mapping_apply(beta, eta, X_out_d, X_in_d, n)
#elif HAVE_CUDA
    call cuda_heaviside_mapping_apply(beta, eta, X_out_d, X_in_d, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_heaviside_mapping_apply

  !> Dispatch Heaviside mapping chain rule to active device backend.
  !! @param beta Projection sharpness parameter.
  !! @param eta Projection threshold parameter.
  !! @param sens_out_d Output sensitivity on device.
  !! @param sens_in_d Input sensitivity on device.
  !! @param X_in_d Input field on device.
  !! @param n Number of dofs.
  subroutine device_heaviside_mapping_apply_backward(beta, eta, &
       sens_out_d, sens_in_d, X_in_d, n)
    real(kind=rp), intent(in) :: beta
    real(kind=rp), intent(in) :: eta
    type(c_ptr) :: sens_out_d
    type(c_ptr) :: sens_in_d
    type(c_ptr) :: X_in_d
    integer :: n
#if HAVE_HIP
    call hip_heaviside_mapping_apply_backward(beta, eta, sens_out_d, &
         sens_in_d, X_in_d, n)
#elif HAVE_CUDA
    call cuda_heaviside_mapping_apply_backward(beta, eta, sens_out_d, &
         sens_in_d, X_in_d, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_heaviside_mapping_apply_backward

end module device_heaviside_mapping
