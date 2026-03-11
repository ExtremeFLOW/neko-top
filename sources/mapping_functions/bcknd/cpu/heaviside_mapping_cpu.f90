!> @file heaviside_mapping_cpu.f90
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
!> CPU backend for smooth Heaviside mapping operations.
module heaviside_mapping_cpu
  use num_types, only: rp
  implicit none
  private

  public :: heaviside_mapping_apply_cpu, &
       heaviside_mapping_apply_backward_cpu

contains

  !> @brief Apply smooth Heaviside mapping on CPU.
  !! @param[in] beta Projection sharpness parameter.
  !! @param[in] eta Projection threshold parameter.
  !! @param[out] X_out Mapped field values.
  !! @param[in] X_in Unmapped field values.
  !! @param[in] n Number of degrees of freedom.
  subroutine heaviside_mapping_apply_cpu(beta, eta, X_out, X_in, n)
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: beta, eta
    real(kind=rp), dimension(n), intent(out) :: X_out
    real(kind=rp), dimension(n), intent(in) :: X_in
    real(kind=rp) :: den, tanh_beta_eta

    tanh_beta_eta = tanh(beta * eta)
    den = tanh_beta_eta + tanh(beta * (1.0_rp - eta))

    X_out = heaviside_mapping_kernel(beta, eta, den, tanh_beta_eta, X_in)

  end subroutine heaviside_mapping_apply_cpu

  !> @brief Apply smooth Heaviside mapping chain rule on CPU.
  !! @param[in] beta Projection sharpness parameter.
  !! @param[in] eta Projection threshold parameter.
  !! @param[out] sens_out Sensitivity with respect to unprojected field.
  !! @param[in] sens_in Sensitivity with respect to projected field.
  !! @param[in] X_in Unprojected field values.
  !! @param[in] n Number of degrees of freedom.
  subroutine heaviside_mapping_apply_backward_cpu(beta, eta, &
       sens_out, sens_in, X_in, n)
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: beta, eta
    real(kind=rp), dimension(n), intent(out) :: sens_out
    real(kind=rp), dimension(n), intent(in) :: sens_in
    real(kind=rp), dimension(n), intent(in) :: X_in
    real(kind=rp) :: den

    den = tanh(beta * eta) + tanh(beta * (1.0_rp - eta))

    sens_out = heaviside_mapping_backward_kernel(beta, eta, den, &
         sens_in, X_in)

  end subroutine heaviside_mapping_apply_backward_cpu

  !> @brief Elemental kernel for smooth Heaviside mapping.
  !! @param[in] beta Projection sharpness parameter.
  !! @param[in] eta Projection threshold parameter.
  !! @param[in] den Projection denominator.
  !! @param[in] tanh_beta_eta Precomputed \f$\tanh(\beta\eta)\f$.
  !! @param[in] X_in Unmapped scalar value.
  !! @return Mapped scalar value.
  elemental function heaviside_mapping_kernel(beta, eta, den, &
       tanh_beta_eta, X_in) result(X_out)
    real(kind=rp), intent(in) :: beta, eta, den, tanh_beta_eta
    real(kind=rp), intent(in) :: X_in
    real(kind=rp) :: X_out

    X_out = (tanh_beta_eta + tanh(beta * (X_in - eta))) / den

  end function heaviside_mapping_kernel

  !> @brief Elemental kernel for smooth Heaviside mapping chain rule.
  !! @param[in] beta Projection sharpness parameter.
  !! @param[in] eta Projection threshold parameter.
  !! @param[in] den Projection denominator.
  !! @param[in] sens_in Sensitivity with respect to projected scalar value.
  !! @param[in] X_in Unprojected scalar value.
  !! @return Sensitivity with respect to unprojected scalar value.
  elemental function heaviside_mapping_backward_kernel(beta, eta, den, &
       sens_in, X_in) result(sens_out)
    real(kind=rp), intent(in) :: beta, eta, den
    real(kind=rp), intent(in) :: sens_in, X_in
    real(kind=rp) :: sens_out, arg, tanh_arg

    arg = beta * (X_in - eta)
    tanh_arg = tanh(arg)
    sens_out = beta * (1.0_rp - tanh_arg * tanh_arg) / den * sens_in

  end function heaviside_mapping_backward_kernel

end module heaviside_mapping_cpu
