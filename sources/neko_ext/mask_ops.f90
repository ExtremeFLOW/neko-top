!> @file mask_ops.f90
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
!
!> Some common Masking operations we may need
module mask_ops
  use field, only: field_t
  use neko_config, only: NEKO_BCKND_DEVICE
  use num_types, only: rp, xp
  use utils, only: neko_error
  use point_zone, only: point_zone_t
  use scratch_registry, only: neko_scratch_registry
  use field_math, only: field_cfill, field_copy, field_rone
  use device_math, only: device_copy, device_glsc2, device_masked_copy_aligned
  use math_ext, only: copy_mask
  use math, only: copy, glsc2
  use vector, only: vector_t
  use coefs, only: coef_t
  implicit none

  private
  public :: mask_exterior_const, mask_exterior_fld, compute_masked_volume

  interface mask_exterior_const
     module procedure mask_exterior_const_vec
     module procedure mask_exterior_const_fld
  end interface mask_exterior_const

contains

  !> @brief Force everything outside the mask to be a constant value
  !! @param[in,out] vec The field being masked
  !! @param[in,out] zone The zone being applied.
  !! @param[in] const The value to be filled
  subroutine mask_exterior_const_vec(vec, zone, const)
    type(vector_t), intent(inout) :: vec
    class(point_zone_t), intent(inout) :: zone
    real(kind=rp), intent(in) :: const
    type(field_t), pointer :: work
    integer :: temp_indices(1), i

    ! To be discussed
    ! From what I understand with this vector/field distinction is that
    ! ultimately the design will only contain the GLL pts inside the
    ! the optimization domain, correct?
    !
    ! If this is the case, then this function makes no sense since it forces
    ! the vector and field to be the same size.
    !
    ! Alternatively, this vector/field distinction is just to make the types
    ! compatible with MMA, in which case we can continue how things are here.
    !
    ! In any case, it's a bit confusing and we should throw an error if the
    ! sizes are different
    call neko_scratch_registry%request_field(work, temp_indices(1), .false.)

    if (vec%size() .ne. work%size()) then
       call neko_error('vector and field are of incompatible dimension')
    end if

    ! fill background fld
    call field_cfill(work, const)

    ! copy the fld in the masked region
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_masked_copy_aligned(work%x_d, vec%x_d, &
            zone%mask%get_d(), work%size(), zone%size)
    else
       do i = 1, zone%size
          work%x(zone%mask%get(i), 1, 1, 1) = vec%x(zone%mask%get(i))
       end do
    end if

    ! copy over
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(vec%x_d, work%x_d, work%size())
    else
       call copy(vec%x, work%x, work%size())
    end if

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine mask_exterior_const_vec

  !> @brief Force everything outside the mask to be a constant value
  !! @param[in,out] fld The field being masked
  !! @param[in,out] zone The zone being applied.
  !! @param[in] const The value to be filled
  subroutine mask_exterior_const_fld(fld, zone, const)
    type(field_t), intent(inout) :: fld
    class(point_zone_t), intent(inout) :: zone
    real(kind=rp), intent(in) :: const
    type(field_t), pointer :: work
    integer :: temp_indices(1)

    call neko_scratch_registry%request_field(work, temp_indices(1), .false.)

    ! fill background fld
    call field_cfill(work, const)

    ! copy the fld in the masked region
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_masked_copy_aligned(work%x_d, fld%x_d, &
            zone%mask%get_d(), fld%size(), zone%size)
    else
       call copy_mask(work%x, fld%x, fld%size(), zone%mask%get(), zone%size)
    end if

    ! copy over
    call field_copy(fld, work)

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine mask_exterior_const_fld

  !> @brief Force everything outside the mask to be a background field
  !! @param[in,out] fld The field being masked
  !! @param[in,out] zone The zone being applied.
  !! @param[in] background The background field
  subroutine mask_exterior_fld(fld, zone, background)
    type(field_t), intent(inout) :: fld
    class(point_zone_t), intent(inout) :: zone
    type(field_t), intent(inout) :: background
    type(field_t), pointer :: work
    integer :: temp_indices(1)

    call neko_scratch_registry%request_field(work, temp_indices(1), .false.)

    ! fill background fld
    call field_copy(work, background)

    ! copy the fld in the masked region
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_masked_copy_aligned(work%x_d, fld%x_d, &
            zone%mask%get_d(), fld%size(), zone%size)
    else
       call copy_mask(work%x, fld%x, fld%size(), zone%mask%get(), zone%size)
    end if

    ! copy over
    call field_copy(fld, work)

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine mask_exterior_fld

  !> @brief Compute the volume of the domain contained within the mask
  !! @param[in] mask The mask considered.
  !! @param[in] coef Coefficients defined on a given mesh.
  function compute_masked_volume(mask, coef)
    class(point_zone_t), intent(inout) :: mask !this should be (in)
    class(coef_t), intent(in) :: coef
    real(kind=rp) :: compute_masked_volume
    type(field_t), pointer :: work
    integer :: temp_indices(1)
    real(kind=rp) :: tmp
    integer :: n

    ! This would be much smarter with a kernel similar to masked_glsc2
    ! When that kernel get written, we can update this function.

    call neko_scratch_registry%request_field(work, temp_indices(1), .false.)
    n = work%size()
    call field_rone(work)
    call mask_exterior_const_fld(work, mask, 0.0_rp)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       tmp = device_glsc2(work%x_d, coef%B_d, n)
    else
       tmp = glsc2(work%x, coef%B, n)
    end if
    call neko_scratch_registry%relinquish_field(temp_indices)
    compute_masked_volume = tmp
  end function compute_masked_volume

end module mask_ops
