! Copyright (c) 2024, The Neko Authors
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
!> Some common Masking operations we may need
!! TODO
!! I know I got outvoted on defining a mask as a field of booleans...
!! You were right, `point_zones` are convenient, but I don't like that the
!! `invert` is defined on initialization (maybe I'm misunderstanding)
!! So we can chose either the interior or the exterior at the beginning.
!!
!! There's probably a smarter way of doing this, but if you want to zero the
!! interior vs the exterior, one will involve some awkward memory shifting.
module mask_ops
  use field, only: field_t
  use neko_config, only: NEKO_BCKND_DEVICE
  use num_types, only: rp, xp
  use utils, only: neko_error
  use point_zone, only: point_zone_t
  use scratch_registry, only : neko_scratch_registry
  use field_math, only: field_cfill, field_copy
  use comm
  implicit none

  private
  public :: mask_exterior_const, mask_interior_const, masked_glsc2 

contains

  !> @brief Force everything outside the mask to be a constant value
  !! @param[in,out] fld The field being masked
  !! @param[in,out] The mask being applied.
  !! @param[in] The value to be filled
  !! NOTE! we always assume the incoming point zone describes the interior of
  !! the mask, ie, `invert = .false.`
  subroutine mask_exterior_const(fld, mask, const)
    type(field_t), intent(inout) :: fld
    class(point_zone_t), intent(inout) :: mask
    real(kind=rp), intent(in) :: const
    type(field_t), pointer :: work
    integer :: temp_indices(1)
    integer :: i

    call neko_scratch_registry%request_field(work , temp_indices(1))

    ! fill background fld
    call field_cfill(work, const)

    ! copy the fld in the masked region
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_error('GPU not supported for masks yet')
    else
       do i = 1, mask%size
          work%x(mask%mask(i), 1, 1, 1) = fld%x(mask%mask(i), 1, 1, 1)
       end do
    end if

    ! copy over
    call field_copy(fld, work)

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine mask_exterior_const

  !> @brief Force everything inside the mask to be a constant value
  !! @param[in,out] fld The field being masked
  !! @param[in,out] The mask being applied.
  !! @param[in] The value to be filled
  !! NOTE! we always assume the incoming point zone describes the interior of
  !! the mask, ie, `invert = .false.`
  subroutine mask_interior_const(fld, mask, const)
    type(field_t), intent(inout) :: fld
    class(point_zone_t), intent(inout) :: mask
    real(kind=rp), intent(in) :: const
    integer :: i

    ! copy the fld in the masked region
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_error('GPU not supported for masks yet')
    else
       do i = 1, mask%size
          fld%x(mask%mask(i), 1, 1, 1) = const 
       end do
    end if

  end subroutine mask_interior_const

  !> Weighted inner product \f$ a^T b \f$ masked
  !! we probably don't need the n as an input..
  function masked_glsc2(a, b, mask, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    class(point_zone_t), intent(inout) :: mask
    real(kind=rp) :: masked_glsc2
    real(kind=xp) :: tmp
    integer :: i, ierr

    tmp = 0.0_xp
    do i = 1, mask%size
       tmp = tmp + a(mask%mask(i)) * b(mask%mask(i))
    end do

    call MPI_Allreduce(MPI_IN_PLACE, tmp, 1, &
         MPI_EXTRA_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    masked_glsc2 = tmp
  end function masked_glsc2
end module mask_ops