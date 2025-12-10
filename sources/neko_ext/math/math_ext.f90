!> @file math_ext.f90
!! @copyright
!! Copyright (c) 2024-2025, The Neko-TOP Authors
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

module math_ext
  use num_types, only: rp, xp
  use comm, only: NEKO_COMM, MPI_EXTRA_PRECISION
  use mpi_f08, only: MPI_Allreduce, MPI_SUM, MPI_IN_PLACE
  implicit none

contains

  !> @brief Copy within a mask. NOTE, this differs from `masked_copy` in the
  !! indexing.
  !! \f$ a_i = b_i, for i in mask \f$
  subroutine copy_mask(a, b, size, mask, mask_size)
    integer, intent(in) :: size, mask_size
    real(kind=rp), dimension(size), intent(inout) :: a
    real(kind=rp), dimension(size), intent(in) :: b
    integer, dimension(mask_size), intent(in) :: mask
    integer :: i

    do i = 1, mask_size
       a(mask(i)) = b(mask(i))
    end do

  end subroutine copy_mask

  !> @brief Add a constant to a masked vector.
  !! \f$ a_i = a_i + c, for i in mask \f$
  subroutine cadd_mask(a, c, size, mask, mask_size)
    integer, intent(in) :: size, mask_size
    real(kind=rp), dimension(size), intent(inout) :: a
    real(kind=rp), intent(in) :: c
    integer, dimension(mask_size), intent(in) :: mask
    integer :: i

    do i = 1, mask_size
       a(mask(i)) = a(mask(i)) + c
    end do

  end subroutine cadd_mask

  !> @brief Invert a masked vector.
  !! \f$ a_i = 1/a_i, for i in mask \f$
  subroutine invcol1_mask(a, size, mask, mask_size)
    integer, intent(in) :: size, mask_size
    real(kind=rp), dimension(size), intent(inout) :: a
    integer, dimension(mask_size), intent(in) :: mask
    integer :: i

    do i = 1, mask_size
       a(mask(i)) = 1.0_rp / a(mask(i))
    end do

  end subroutine invcol1_mask

  !> @brief Multiply a masked vector by a constant.
  !! \f$ a_i = c * a_i, for i in mask \f$
  subroutine cmult_mask(a, c, size, mask, mask_size)
    integer, intent(in) :: size, mask_size
    real(kind=rp), dimension(size), intent(inout) :: a
    real(kind=rp), intent(in) :: c
    integer, dimension(mask_size), intent(in) :: mask
    integer :: i

    do i = 1, mask_size
       a(mask(i)) = c * a(mask(i))
    end do

  end subroutine cmult_mask

  !> @brief Multiply 2 masked vectors. Save the result in a new vector.
  !! \f$ a_i = b_i * c_i, for i in mask \f$
  subroutine col2_mask(a, b, size, mask, mask_size)
    integer, intent(in) :: size, mask_size
    real(kind=rp), dimension(size), intent(inout) :: a
    real(kind=rp), dimension(size), intent(in) :: b
    integer, dimension(mask_size), intent(in) :: mask
    integer :: i

    do i = 1, mask_size
       a(mask(i)) = a(mask(i)) * b(mask(i))
    end do

  end subroutine col2_mask

  !> @brief Multiply 2 masked vectors. Save the result in a new vector.
  !! \f$ a_i = b_i * c_i, for i in mask \f$
  subroutine col3_mask(a, b, c, size, mask, mask_size)
    integer, intent(in) :: size, mask_size
    real(kind=rp), dimension(size), intent(inout) :: a
    real(kind=rp), dimension(size), intent(in) :: b, c
    integer, dimension(mask_size), intent(in) :: mask
    integer :: i

    do i = 1, mask_size
       a(mask(i)) = b(mask(i)) * c(mask(i))
    end do

  end subroutine col3_mask

  !> @brief Subtract 2 masked vectors. Save the result in a new vector.
  !! \f$ a_i = b_i - c_i, for i in mask \f$
  subroutine sub3_mask(a, b, c, size, mask, mask_size)
    integer, intent(in) :: size, mask_size
    real(kind=rp), dimension(size), intent(inout) :: a
    real(kind=rp), dimension(size), intent(in) :: b, c
    integer, dimension(mask_size), intent(in) :: mask
    integer :: i

    do i = 1, mask_size
       a(mask(i)) = b(mask(i)) - c(mask(i))
    end do

  end subroutine sub3_mask

  !> @brief Weighted inner product
  !! \f$ a^T b \f$ for indices in the mask
  function glsc2_mask(a, b, size, mask, mask_size)
    integer, intent(in) :: size, mask_size
    real(kind=rp), dimension(size), intent(in) :: a
    real(kind=rp), dimension(size), intent(in) :: b
    integer, dimension(mask_size), intent(in) :: mask
    real(kind=rp) :: glsc2_mask
    real(kind=xp) :: tmp
    integer :: i, ierr

    tmp = 0.0_xp
    do i = 1, mask_size
       tmp = tmp + a(mask(i)) * b(mask(i))
    end do

    call MPI_Allreduce(MPI_IN_PLACE, tmp, 1, &
         MPI_EXTRA_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    glsc2_mask = real(tmp, kind=rp)
  end function glsc2_mask
end module math_ext
