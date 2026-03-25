!> @file heaviside_mapping.f90
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
!> Smooth Heaviside mapping
module heaviside_mapping
  use num_types, only: rp
  use mapping, only: mapping_t
  use json_module, only: json_file
  use field, only: field_t
  use coefs, only: coef_t
  use neko_config, only: NEKO_BCKND_DEVICE
  use device_heaviside_mapping, only: &
       device_heaviside_mapping_apply, &
       device_heaviside_mapping_apply_backward
  use heaviside_mapping_cpu, only: &
       heaviside_mapping_apply_cpu, &
       heaviside_mapping_apply_backward_cpu
  use json_utils, only: json_get, json_get_or_default
  use utils, only: neko_error
  implicit none
  private

  !> Smooth Heaviside mapping
  !!
  !! \f[
  !! X_\mathrm{out} =
  !! \frac{\tanh(\beta \eta) + \tanh(\beta (X_\mathrm{in} - \eta))}
  !!      {\tanh(\beta \eta) + \tanh(\beta (1-\eta))}
  !! \f]
  !!
  !! with \f$\beta > 0\f$ and \f$\eta \in [0,1]\f$.
  type, public, extends(mapping_t) :: heaviside_mapping_t
     !> Projection sharpness parameter
     real(kind=rp) :: beta
     !> Threshold parameter
     real(kind=rp) :: eta

   contains
     !> Constructor from json.
     procedure, pass(this) :: init => heaviside_mapping_init_from_json
     !> Constructor from attributes.
     procedure, pass(this) :: init_from_attributes => &
          heaviside_mapping_init_from_attributes
     !> Destructor.
     procedure, pass(this) :: free => heaviside_mapping_free
     !> Apply the forward mapping.
     procedure, pass(this) :: forward_mapping => heaviside_mapping_forward
     !> Apply the chain-rule mapping.
     procedure, pass(this) :: backward_mapping => heaviside_mapping_backward
  end type heaviside_mapping_t

contains

  !> Constructor from json.
  subroutine heaviside_mapping_init_from_json(this, json, coef)
    class(heaviside_mapping_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(coef_t), intent(inout) :: coef
    real(kind=rp) :: beta, eta

    call json_get(json, 'beta', beta)
    call json_get_or_default(json, 'eta', eta, 0.5_rp)

    call this%init_base(json, coef, "heaviside_mapping")
    call this%init_from_attributes(coef, beta, eta)
  end subroutine heaviside_mapping_init_from_json

  !> Constructor from attributes.
  subroutine heaviside_mapping_init_from_attributes(this, coef, beta, eta)
    class(heaviside_mapping_t), intent(inout) :: this
    type(coef_t), intent(inout) :: coef
    real(kind=rp), intent(in) :: beta
    real(kind=rp), intent(in) :: eta

    if (beta .le. 0.0_rp) then
       call neko_error('"beta" must be > 0 in heaviside mapping')
    end if

    if (eta .lt. 0.0_rp .or. eta .gt. 1.0_rp) then
       call neko_error('"eta" must be in [0, 1] in heaviside mapping')
    end if

    this%beta = beta
    this%eta = eta
  end subroutine heaviside_mapping_init_from_attributes

  !> Destructor.
  subroutine heaviside_mapping_free(this)
    class(heaviside_mapping_t), intent(inout) :: this

    call this%free_base()
  end subroutine heaviside_mapping_free

  !> Apply the forward mapping.
  subroutine heaviside_mapping_forward(this, X_out, X_in)
    class(heaviside_mapping_t), intent(inout) :: this
    type(field_t), intent(in) :: X_in
    type(field_t), intent(inout) :: X_out
    call heaviside_mapping_apply(this%beta, this%eta, X_out, X_in)
  end subroutine heaviside_mapping_forward

  !> Apply the chain rule.
  subroutine heaviside_mapping_backward(this, sens_out, sens_in, X_in)
    class(heaviside_mapping_t), intent(inout) :: this
    type(field_t), intent(in) :: X_in
    type(field_t), intent(in) :: sens_in
    type(field_t), intent(inout) :: sens_out
    call heaviside_mapping_apply_backward(this%beta, this%eta, &
         sens_out, sens_in, X_in)
  end subroutine heaviside_mapping_backward

  !> Apply smooth Heaviside mapping.
  !! @param beta Projection sharpness parameter.
  !! @param eta Projection threshold parameter.
  !! @param X_out Mapped field.
  !! @param X_in Unmapped field.
  subroutine heaviside_mapping_apply(beta, eta, X_out, X_in)
    real(kind=rp), intent(in) :: beta, eta
    type(field_t), intent(in) :: X_in
    type(field_t), intent(inout) :: X_out
    integer :: n
    real(kind=rp) :: den, beta_eta

    beta_eta = beta * eta
    den = tanh(beta_eta) + tanh(beta * (1.0_rp - eta))

    if (abs(den) .le. tiny(den)) then
       call neko_error('invalid denominator in heaviside mapping')
    end if

    n = X_in%dof%size()
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_heaviside_mapping_apply(beta, eta, &
            X_out%x_d, X_in%x_d, n)
    else
       call heaviside_mapping_apply_cpu(beta, eta, &
            X_out%x, X_in%x, n)
    end if
  end subroutine heaviside_mapping_apply

  !> Apply chain rule for smooth Heaviside mapping.
  !! @param beta Projection sharpness parameter.
  !! @param eta Projection threshold parameter.
  !! @param sens_out Sensitivity with respect to unprojected field.
  !! @param sens_in Sensitivity with respect to projected field.
  !! @param X_in Unprojected field.
  subroutine heaviside_mapping_apply_backward(beta, eta, sens_out, sens_in, &
       X_in)
    real(kind=rp), intent(in) :: beta, eta
    type(field_t), intent(in) :: X_in
    type(field_t), intent(in) :: sens_in
    type(field_t), intent(inout) :: sens_out
    integer :: n
    real(kind=rp) :: den, beta_eta

    beta_eta = beta * eta
    den = tanh(beta_eta) + tanh(beta * (1.0_rp - eta))

    if (abs(den) .le. tiny(den)) then
       call neko_error('invalid denominator in heaviside mapping')
    end if

    n = X_in%dof%size()
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_heaviside_mapping_apply_backward(beta, eta, &
            sens_out%x_d, sens_in%x_d, X_in%x_d, n)
    else
       call heaviside_mapping_apply_backward_cpu(beta, eta, &
            sens_out%x, sens_in%x, X_in%x, n)
    end if
  end subroutine heaviside_mapping_apply_backward

end module heaviside_mapping
