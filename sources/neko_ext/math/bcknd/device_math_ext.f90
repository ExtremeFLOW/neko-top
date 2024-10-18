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
module device_math_ext
  use utils, only : neko_error
  use num_types, only : rp, c_rp
  use, intrinsic :: iso_c_binding
  implicit none

#ifdef HAVE_HIP

#elif HAVE_CUDA

  interface
     subroutine cuda_cadd_mask(a_d, c, size, mask_d, mask_size) &
          bind(c, name='cuda_cadd_mask')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: a_d
       real(c_rp) :: c
       integer(c_int) :: size
       type(c_ptr), value :: mask_d
       integer(c_int) :: mask_size
     end subroutine cuda_cadd_mask
  end interface
  interface
     subroutine cuda_invcol1_mask(a_d, size, mask_d, mask_size) &
          bind(c, name='cuda_invcol1_mask')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: a_d
       integer(c_int) :: size
       type(c_ptr), value :: mask_d
       integer(c_int) :: mask_size
     end subroutine cuda_invcol1_mask
  end interface
  interface
     subroutine cuda_col2_mask(a_d, b_d, size, mask_d, mask_size) &
          bind(c, name='cuda_col2_mask')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: a_d
       type(c_ptr), value :: b_d
       integer(c_int) :: size
       type(c_ptr), value :: mask_d
       integer(c_int) :: mask_size
     end subroutine cuda_col2_mask
  end interface
  interface
     subroutine cuda_col3_mask(a_d, b_d, c_d, size, mask_d, mask_size) &
          bind(c, name='cuda_col3_mask')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: a_d
       type(c_ptr), value :: b_d
       type(c_ptr), value :: c_d
       integer(c_int) :: size
       type(c_ptr), value :: mask_d
       integer(c_int) :: mask_size
     end subroutine cuda_col3_mask
  end interface
  interface
     subroutine cuda_sub3_mask(a_d, b_d, c_d, size, mask_d, mask_size) &
          bind(c, name='cuda_sub3_mask')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: a_d
       type(c_ptr), value :: b_d
       type(c_ptr), value :: c_d
       integer(c_int) :: size
       type(c_ptr), value :: mask_d
       integer(c_int) :: mask_size
     end subroutine cuda_sub3_mask
  end interface

#elif HAVE_OPENCL

#endif

contains

  subroutine device_cadd_mask(a_d, c, size, mask_d, mask_size)
    type(c_ptr) :: a_d
    real(kind=rp), intent(in) :: c
    integer :: size
    type(c_ptr) :: mask_d
    integer :: mask_size
#ifdef HAVE_CUDA
    call cuda_cadd_mask(a_d, c, size, mask_d, mask_size)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_cadd_mask

  subroutine device_invcol1_mask(a_d, size, mask_d, mask_size)
    type(c_ptr) :: a_d
    integer :: size
    type(c_ptr) :: mask_d
    integer :: mask_size
#ifdef HAVE_CUDA
    call cuda_invcol1_mask(a_d, size, mask_d, mask_size)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_invcol1_mask

  subroutine device_col2_mask(a_d, b_d, size, mask_d, mask_size)
    type(c_ptr) :: a_d
    type(c_ptr) :: b_d
    integer :: size
    type(c_ptr) :: mask_d
    integer :: mask_size
#ifdef HAVE_CUDA
    call cuda_col2_mask(a_d, b_d, size, mask_d, mask_size)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_col2_mask

  subroutine device_col3_mask(a_d, b_d, c_d, size, mask_d, mask_size)
    type(c_ptr) :: a_d
    type(c_ptr) :: b_d
    type(c_ptr) :: c_d
    integer :: size
    type(c_ptr) :: mask_d
    integer :: mask_size
#ifdef HAVE_CUDA
    call cuda_col3_mask(a_d, b_d, c_d, size, mask_d, mask_size)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_col3_mask

  subroutine device_sub3_mask(a_d, b_d, c_d, size, mask_d, mask_size)
    type(c_ptr) :: a_d
    type(c_ptr) :: b_d
    type(c_ptr) :: c_d
    integer :: size
    type(c_ptr) :: mask_d
    integer :: mask_size
#ifdef HAVE_CUDA
    call cuda_sub3_mask(a_d, b_d, c_d, size, mask_d, mask_size)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_sub3_mask


end module device_math_ext
