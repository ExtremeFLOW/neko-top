! Copyright (c) 2023, The Neko Authors
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
!
!> A SIMP mapping of coefficients
module SIMP_mapping
  use num_types, only: rp
  use mapping, only: mapping_t
  use json_module, only: json_file
  use field, only: field_t
  use coefs, only: coef_t
  use neko_config, only: NEKO_BCKND_DEVICE
  use device_SIMP_mapping, only: device_SIMP_mapping_apply, &
       device_SIMP_mapping_apply_backward
  use SIMP_mapping_cpu, only: SIMP_mapping_apply_cpu, &
       SIMP_mapping_apply_backward_cpu
  use json_utils, only: json_get, json_get_or_default
  implicit none
  private

  !> A SIMP mapping of coefficients
  !!
  !! \f$f(x) = f_{min} + (f_{max} - f_{min}) (x)^p\f$

  type, public, extends(mapping_t) :: SIMP_mapping_t
     !> minimum value
     real(kind=rp) :: f_min
     !> maximum value
     real(kind=rp) :: f_max
     !> penalty parameter
     real(kind=rp) :: p

   contains
     !> Constructor from json.
     procedure, pass(this) :: init => SIMP_mapping_init_from_json
     !> Actual constructor.
     procedure, pass(this) :: init_from_attributes => &
          SIMP_mapping_init_from_attributes
     !> Destructor.
     procedure, pass(this) :: free => SIMP_mapping_free
     !> Apply the forward mapping
     procedure, pass(this) :: forward_mapping => SIMP_forward_mapping
     !> Apply the adjoint mapping
     procedure, pass(this) :: backward_mapping => SIMP_backward_mapping
  end type SIMP_mapping_t

contains

  !> Constructor from json.
  subroutine SIMP_mapping_init_from_json(this, json, coef)
    class(SIMP_mapping_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(coef_t), intent(inout) :: coef
    real(kind=rp) :: f_min, f_max, p

    call json_get_or_default(json, 'f_min', f_min, 0.0_rp)
    call json_get(json, 'f_max', f_max)
    call json_get_or_default(json, 'p', p, 1.0_rp)

    call this%init_base(json, coef)
    call this%init_from_attributes(coef, f_min, f_max, p)

  end subroutine SIMP_mapping_init_from_json

  !> Actual constructor.
  subroutine SIMP_mapping_init_from_attributes(this, coef, f_min, f_max, p)
    class(SIMP_mapping_t), intent(inout) :: this
    type(coef_t), intent(inout) :: coef
    real(kind=rp), intent(in) :: f_min, f_max, p

    this%f_min = f_min
    this%f_max = f_max
    this%p = p

  end subroutine SIMP_mapping_init_from_attributes

  !> Destructor.
  subroutine SIMP_mapping_free(this)
    class(SIMP_mapping_t), intent(inout) :: this

    call this%free_base()

  end subroutine SIMP_mapping_free

  !> Apply the mapping
  !! @param this mapping
  !! @param X_out mapped field
  !! @param X_in unmapped field
  subroutine SIMP_forward_mapping(this, X_out, X_in)
    class(SIMP_mapping_t), intent(inout) :: this
    type(field_t), intent(in) :: X_in
    type(field_t), intent(inout) :: X_out
    integer :: n

    n = X_in%dof%size()
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_SIMP_mapping_apply(this%f_min, this%f_max, this%p, &
            X_out%x_d, X_in%x_d, n)
    else
       call SIMP_mapping_apply_cpu(this%f_min, this%f_max, this%p, &
            X_out%x, X_in%x, n)
    end if

  end subroutine SIMP_forward_mapping


  !> Apply the  chain rule
  !! @param this mapping
  !! @param sens_out is the sensitivity with respect to the unfiltered design
  !! @param sens_in is the sensitivity with respect to the filtered design
  !! @param X_in unmapped field
  subroutine SIMP_backward_mapping(this, sens_out, sens_in, X_in)
    class(SIMP_mapping_t), intent(inout) :: this
    type(field_t), intent(in) :: X_in
    type(field_t), intent(in) :: sens_in
    type(field_t), intent(inout) :: sens_out
    integer :: n

    n = X_in%dof%size()
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_SIMP_mapping_apply_backward(this%f_min, this%f_max, this%p, &
            sens_out%x_d, sens_in%x_d, X_in%x_d, n)
    else
       call SIMP_mapping_apply_backward_cpu(this%f_min, this%f_max, this%p, &
            sens_out%x, sens_in%x, X_in%x, n)
    end if

  end subroutine SIMP_backward_mapping

end module SIMP_mapping
