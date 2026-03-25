!> @file linear_mapping.f90
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
!> A linear mapping of coefficients
module linear_mapping
  use num_types, only: rp
  use field_math, only: field_copy, field_cmult, field_cadd
  use mapping, only: mapping_t
  use num_types, only: rp
  use json_module, only: json_file
  use field, only: field_t
  use coefs, only: coef_t
  use json_utils, only: json_get, json_get_or_default
  implicit none
  private

  !> A linear mapping of coefficients $f(x) = f_{min} + (f_{max} - f_{min}) x$
  type, public, extends(mapping_t) :: linear_mapping_t
     !> minimum value
     real(kind=rp) :: f_min
     !> maximum value
     real(kind=rp) :: f_max

   contains
     !> Constructor from json.
     procedure, pass(this) :: init => linear_mapping_init_from_json
     !> Actual constructor.
     procedure, pass(this) :: init_from_attributes => &
          linear_mapping_init_from_attributes
     !> Destructor.
     procedure, pass(this) :: free => linear_mapping_free
     !> Apply the forward mapping
     procedure, pass(this) :: forward_mapping => linear_forward_mapping
     !> Apply the adjoint mapping
     procedure, pass(this) :: backward_mapping => &
          linear_backward_mapping
  end type linear_mapping_t

contains

  !> Constructor from json.
  subroutine linear_mapping_init_from_json(this, json, coef)
    class(linear_mapping_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(coef_t), intent(inout) :: coef
    real(kind=rp) :: f_min, f_max

    call json_get_or_default(json, 'f_min', f_min, 0.0_rp)
    call json_get(json, 'f_max', f_max)

    call this%init_base(json, coef, "linear_mapping")
    call linear_mapping_init_from_attributes(this, coef, f_min, f_max)

  end subroutine linear_mapping_init_from_json

  !> Actual constructor.
  subroutine linear_mapping_init_from_attributes(this, coef, f_min, f_max)
    class(linear_mapping_t), intent(inout) :: this
    type(coef_t), intent(inout) :: coef
    real(kind=rp), intent(in) :: f_min, f_max

    this%f_min = f_min
    this%f_max = f_max

  end subroutine linear_mapping_init_from_attributes

  !> Destructor.
  subroutine linear_mapping_free(this)
    class(linear_mapping_t), intent(inout) :: this

    call this%free_base()

  end subroutine linear_mapping_free

  !> Apply the mapping
  !! @param this the mapping
  !! @param X_out mapped field
  !! @param X_in unmapped field
  subroutine linear_forward_mapping(this, X_out, X_in)
    class(linear_mapping_t), intent(inout) :: this
    type(field_t), intent(in) :: X_in
    type(field_t), intent(inout) :: X_out

    ! x_out = f_min + (f_max - f_min) * x_in
    call field_copy(X_out, X_in)
    call field_cmult(X_out, this%f_max - this%f_min)
    call field_cadd(X_out, this%f_min)

  end subroutine linear_forward_mapping


  !> Apply the  chain rule
  !! @param this the mapping
  !! @param X_in unmapped field
  !! @param sens_out is the sensitivity with respect to the unfiltered design
  !! @param sens_in is the sensitivity with respect to the filtered design
  subroutine linear_backward_mapping(this, sens_out, sens_in, X_in)
    class(linear_mapping_t), intent(inout) :: this
    type(field_t), intent(in) :: X_in
    type(field_t), intent(in) :: sens_in
    type(field_t), intent(inout) :: sens_out

    ! df/dx_in = df/dx_out * dx_out/dx_in

    ! dx_out/dx_in = (f_max - f_min)

    call field_copy(sens_out, sens_in)
    call field_cmult(sens_out, this%f_max - this%f_min)

  end subroutine linear_backward_mapping

end module linear_mapping
