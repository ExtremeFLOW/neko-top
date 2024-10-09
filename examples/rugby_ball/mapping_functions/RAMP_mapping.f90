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
!> A RAMP mapping of coefficients
module RAMP_mapping
  use num_types, only : rp
  use field_math
  use mapping, only: mapping_t
  use num_types, only : rp
  use json_module, only : json_file
  use field_registry, only : neko_field_registry
  use field, only : field_t
  	 use coefs, only: coef_t
    implicit none
  private

	!> A RAMP mapping of coefficients 
	!! $f(x) = f_{min} + (f_{max} - f_{min}) \frac{q + 1}{q + x}$ 
  type, public, extends(mapping_t) :: RAMP_mapping_t
  !> minimum value
  real(kind=rp) :: f_min
  !> maximum value
  real(kind=rp) :: f_max
  !> penalty parameter
  real(kind=rp) :: q



   contains
     !> Constructor from json.
     procedure, pass(this) :: init => RAMP_mapping_init_from_json
     !> Actual constructor.
     procedure, pass(this) :: init_from_attributes => &
          RAMP_mapping_init_from_attributes
     !> Destructor.
     procedure, pass(this) :: free => RAMP_mapping_free
     !> Apply the forward mapping
     procedure, pass(this) :: apply_forward => RAMP_mapping_apply
     !> Apply the adjoint mapping
     procedure, pass(this) :: apply_backward => RAMP_mapping_apply_backward
  end type RAMP_mapping_t

contains

  !> Constructor from json.
  subroutine RAMP_mapping_init_from_json(this, json, coef)
    class(RAMP_mapping_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(coef_t), intent(inout) :: coef

    ! do the JSON stuff later
    this%f_min = 0.0_rp
    this%f_max = 1000.0_rp
    this%q = 1.0_rp

    call this%init_base(json, coef)
    call RAMP_mapping_init_from_attributes(this, coef)
   
  end subroutine RAMP_mapping_init_from_json

  !> Actual constructor.
  subroutine RAMP_mapping_init_from_attributes(this, coef)
    class(RAMP_mapping_t), intent(inout) :: this
    type(coef_t), intent(inout) :: coef

    ! there's actually nothing to do here.

  end subroutine RAMP_mapping_init_from_attributes

  !> Destructor.
  subroutine RAMP_mapping_free(this)
    class(RAMP_mapping_t), intent(inout) :: this

    call this%free_base()

  end subroutine RAMP_mapping_free

  !> Apply the mapping
  !! @param X_out mapped field
  !! @param X_in unmapped field
  subroutine RAMP_mapping_apply(this, X_out, X_in)
    class(RAMP_mapping_t), intent(inout) :: this
    type(field_t), intent(in) ::  X_in
    type(field_t), intent(inout) ::  X_out

	 ! x_out = f_min + (f_max - f_min) * (q + 1) / (x_in + q) 

	 ! TODO
	 ! Here we should make a descision, either use field_math for everything
	 ! or write individual backends.
	 !
	 ! I'm pro-field_math simply because these functions don't get executed too
	 ! often...
	 !
	 ! However! If we do descide to write backends for them, I (Harry) call
	 ! shotgun writing them, because they look easy and would be a good 
	 ! introduction to writing a kernel from scratch.
    call field_copy(X_out,X_in)
    call field_cadd(X_out, this%q)
    call field_invcol1(X_out)
    call field_cmult(X_out, (this%f_max - this%f_min)*(this%q + 1.0_rp))
    call field_cadd(X_out, this%f_min)

  end subroutine RAMP_mapping_apply


  !> Apply the  chain rule
  !! @param X_in unmapped field
  !! @param dF_dX_in is the sensitivity with respect to the unfiltered design
  !! @param dF_dX_out is the sensitivity with respect to the filtered design
  subroutine RAMP_mapping_apply_backward(this, dF_dX_in, dF_dX_out, X_in)
    class(RAMP_mapping_t), intent(inout) :: this
    type(field_t), intent(in) ::  X_in
    type(field_t), intent(in) ::  dF_dX_out
    type(field_t), intent(inout) ::  dF_dX_in

	 ! df/dx_in = df/dx_out * dx_out/dx_in 

	 ! dx_out/dx_in = (f_min - f_max) * (q + 1) / (q + x)**2 

    call field_copy(dF_dX_in, X_in)
    call field_cadd(dF_dX_in, this%q)
    call field_invcol1(dF_dX_in)
    call field_col2(dF_dX_in, dF_dX_in)
    call field_cmult(dF_dX_in, (this%f_max - this%f_min)*(this%q + 1.0_rp))
    call field_col2(dF_dX_in, dF_dX_out)

    ! hmmmmmm.... maybe I'm a getting a little anti-field_math hahahaha
    ! That's really inefficient and reads poorly.

  end subroutine RAMP_mapping_apply_backward

end module RAMP_mapping
