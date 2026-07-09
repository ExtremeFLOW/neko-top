!> @file Borrvall_Petersson_mapping_cpu.f90
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
!> CPU backend for Borrvall & Petersson mapping operations.
module Borrvall_Petersson_mapping_cpu
  use num_types, only: rp
  implicit none
  private

  public :: Borrvall_Petersson_mapping_apply_cpu
  public :: Borrvall_Petersson_mapping_apply_backward_cpu

contains

  !> @brief Apply Borrvall & Petersson forward mapping on CPU.
  !! @param[in] f_min Minimum mapped value.
  !! @param[in] f_max Maximum mapped value.
  !! @param[in] q penalty parameter.
  !! @param[out] X_out Mapped field values.
  !! @param[in] X_in Unmapped field values.
  !! @param[in] n Number of degrees of freedom.
  subroutine Borrvall_Petersson_mapping_apply_cpu(f_min, f_max, q, X_out, &
       X_in, n)
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: f_min, f_max, q
    real(kind=rp), dimension(n), intent(out) :: X_out
    real(kind=rp), dimension(n), intent(in) :: X_in

    X_out = Borrvall_Petersson_mapping_kernel(f_min, f_max, q, X_in)

  end subroutine Borrvall_Petersson_mapping_apply_cpu

  !> @brief Apply Borrvall & Petersson chain rule on CPU.
  !! @param[in] f_min Minimum mapped value.
  !! @param[in] f_max Maximum mapped value.
  !! @param[in] q penalty parameter.
  !! @param[out] sens_out Sensitivity with respect to the unmapped field.
  !! @param[in] sens_in Sensitivity with respect to the mapped field.
  !! @param[in] X_in Unmapped field values.
  !! @param[in] n Number of degrees of freedom.
  subroutine Borrvall_Petersson_mapping_apply_backward_cpu(f_min, f_max, q, &
       sens_out, sens_in, X_in, n)
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: f_min, f_max, q
    real(kind=rp), dimension(n), intent(out) :: sens_out
    real(kind=rp), dimension(n), intent(in) :: sens_in
    real(kind=rp), dimension(n), intent(in) :: X_in

    sens_out = Borrvall_Petersson_mapping_backward_kernel(f_min, f_max, q, &
         sens_in, X_in)

  end subroutine Borrvall_Petersson_mapping_apply_backward_cpu

  !> @brief Elemental kernel for Borrvall & Petersson forward mapping.
  !! @param[in] f_min Minimum mapped value.
  !! @param[in] f_max Maximum mapped value.
  !! @param[in] q penalty parameter.
  !! @param[in] X_in Unmapped scalar value.
  !! @return Mapped scalar value.
  elemental function Borrvall_Petersson_mapping_kernel(f_min, f_max, q, X_in) &
       result(X_out)
    real(kind=rp), intent(in) :: f_min, f_max, q
    real(kind=rp), intent(in) :: X_in
    real(kind=rp) :: X_out

    X_out = f_min + (f_max - f_min) * X_in * (1.0_rp + q) / (X_in + q)

  end function Borrvall_Petersson_mapping_kernel

  !> @brief Elemental kernel for Borrvall & Petersson chain rule.
  !! @param[in] f_min Minimum mapped value.
  !! @param[in] f_max Maximum mapped value.
  !! @param[in] q penalty parameter.
  !! @param[in] sens_in Sensitivity with respect to mapped scalar value.
  !! @param[in] X_in Unmapped scalar value.
  !! @return Sensitivity with respect to unmapped scalar value.
  elemental function Borrvall_Petersson_mapping_backward_kernel(f_min, f_max, &
       q, sens_in, X_in) result(sens_out)
    real(kind=rp), intent(in) :: f_min, f_max, q
    real(kind=rp), intent(in) :: sens_in, X_in
    real(kind=rp) :: sens_out

    sens_out = sens_in * (f_max - f_min) * (q + 1.0_rp) * q / ((X_in + q)**2)

  end function Borrvall_Petersson_mapping_backward_kernel

end module Borrvall_Petersson_mapping_cpu
