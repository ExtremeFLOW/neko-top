! Copyright (c) 2021-2023, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
module device_RAMP_mapping
  use utils, only : neko_error
  use num_types, only : rp, c_rp
  use, intrinsic :: iso_c_binding
  implicit none

#if HAVE_HIP

#elif HAVE_CUDA

  interface
     subroutine cuda_convex_down_RAMP_mapping_apply(f_min, f_max, q, &           
          X_out_d, X_in_d, n) &
          bind(c, name = 'cuda_convex_down_RAMP_mapping_apply')
       use, intrinsic :: iso_c_binding
       import c_rp
       real(c_rp) :: f_min
       real(c_rp) :: f_max
       real(c_rp) :: q
       type(c_ptr), value :: X_out_d
       type(c_ptr), value :: X_in_d
       integer(c_int) :: n
     end subroutine cuda_convex_down_RAMP_mapping_apply
  end interface

  interface
     subroutine cuda_convex_down_RAMP_mapping_apply_backward(f_min, f_max, q, &           
          dF_dX_in_d, dF_dX_out_d, X_in_d, n) &
          bind(c, name = 'cuda_convex_down_RAMP_mapping_apply_backward')
       use, intrinsic :: iso_c_binding
       import c_rp
       real(c_rp) :: f_min
       real(c_rp) :: f_max
       real(c_rp) :: q
       type(c_ptr), value :: dF_dX_in_d
       type(c_ptr), value :: dF_dX_out_d
       type(c_ptr), value :: X_in_d
       integer(c_int) :: n
     end subroutine cuda_convex_down_RAMP_mapping_apply_backward
  end interface

  interface
     subroutine cuda_convex_up_RAMP_mapping_apply(f_min, f_max, q, &           
          X_out_d, X_in_d, n) &
          bind(c, name = 'cuda_convex_up_RAMP_mapping_apply')
       use, intrinsic :: iso_c_binding
       import c_rp
       real(c_rp) :: f_min
       real(c_rp) :: f_max
       real(c_rp) :: q
       type(c_ptr), value :: X_out_d
       type(c_ptr), value :: X_in_d
       integer(c_int) :: n
     end subroutine cuda_convex_up_RAMP_mapping_apply
  end interface

  interface
     subroutine cuda_convex_up_RAMP_mapping_apply_backward(f_min, f_max, q, &           
          dF_dX_in_d, dF_dX_out_d, X_in_d, n) &
          bind(c, name = 'cuda_convex_up_RAMP_mapping_apply_backward')
       use, intrinsic :: iso_c_binding
       import c_rp
       real(c_rp) :: f_min
       real(c_rp) :: f_max
       real(c_rp) :: q
       type(c_ptr), value :: dF_dX_in_d
       type(c_ptr), value :: dF_dX_out_d
       type(c_ptr), value :: X_in_d
       integer(c_int) :: n
     end subroutine cuda_convex_up_RAMP_mapping_apply_backward
  end interface
#elif HAVE_OPENCL

#endif

contains

  subroutine device_convex_down_RAMP_mapping_apply(f_min, f_max, q, &          
          X_out_d, X_in_d, n)
    real(kind=rp), intent(in) :: f_min
    real(kind=rp), intent(in) :: f_max
    real(kind=rp), intent(in) :: q
    type(c_ptr) :: X_out_d
    type(c_ptr) :: X_in_d
    integer :: n
#if HAVE_CUDA
    call cuda_convex_down_RAMP_mapping_apply(f_min, f_max, q, &          
          X_out_d, X_in_d, n)
#else
    call neko_error('No device backend configured for device_convex_down_RAMP_mapping_apply')
#endif
  end subroutine device_convex_down_RAMP_mapping_apply

  subroutine device_convex_down_RAMP_mapping_apply_backward(f_min, f_max, q, &          
          dF_dX_in_d, dF_dX_out_d, X_in_d, n)
    real(kind=rp), intent(in) :: f_min
    real(kind=rp), intent(in) :: f_max
    real(kind=rp), intent(in) :: q
    type(c_ptr) :: dF_dX_in_d
    type(c_ptr) :: dF_dX_out_d
    type(c_ptr) :: X_in_d
    integer :: n
#if HAVE_CUDA
    call cuda_convex_down_RAMP_mapping_apply_backward(f_min, f_max, q, &          
          dF_dX_in_d, dF_dX_out_d, X_in_d, n)
#else
    call neko_error('No device backend configured for device_convex_down_RAMP_mapping_apply_backward')
#endif
  end subroutine device_convex_down_RAMP_mapping_apply_backward

  subroutine device_convex_up_RAMP_mapping_apply(f_min, f_max, q, &          
          X_out_d, X_in_d, n)
    real(kind=rp), intent(in) :: f_min
    real(kind=rp), intent(in) :: f_max
    real(kind=rp), intent(in) :: q
    type(c_ptr) :: X_out_d
    type(c_ptr) :: X_in_d
    integer :: n
#if HAVE_CUDA
    call cuda_convex_up_RAMP_mapping_apply(f_min, f_max, q, &          
          X_out_d, X_in_d, n)
#else
    call neko_error('No device backend configured for device_convex_up_RAMP_mapping_apply')
#endif
  end subroutine device_convex_up_RAMP_mapping_apply

  subroutine device_convex_up_RAMP_mapping_apply_backward(f_min, f_max, q, &          
          dF_dX_in_d, dF_dX_out_d, X_in_d, n)
    real(kind=rp), intent(in) :: f_min
    real(kind=rp), intent(in) :: f_max
    real(kind=rp), intent(in) :: q
    type(c_ptr) :: dF_dX_in_d
    type(c_ptr) :: dF_dX_out_d
    type(c_ptr) :: X_in_d
    integer :: n
#if HAVE_CUDA
    call cuda_convex_up_RAMP_mapping_apply_backward(f_min, f_max, q, &          
          dF_dX_in_d, dF_dX_out_d, X_in_d, n)
#else
    call neko_error('No device backend configured for device_convex_up_RAMP_mapping_apply_backward')
#endif
  end subroutine device_convex_up_RAMP_mapping_apply_backward
end module device_RAMP_mapping
