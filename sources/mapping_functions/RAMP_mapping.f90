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
  use mapping, only: mapping_t
  use num_types, only : rp
  use json_module, only : json_file
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use coefs, only: coef_t
  implicit none
  private

  !> A RAMP mapping of coefficients
  !! This is the standard RAMP described in
  !! https://doi.org/10.1007/s001580100129
  !!
  !! $f(x) = f_{min} + (f_{max} - f_{min}) \frac{x}{1 + q(1 - x)}$
  !!
  !!
  !!  |        .
  !!  |        .
  !!  |       .
  !!  |     ..
  !!  |  ...
  !!  |_________
  !!
  !! or a convex up equivelent used by Borrvall & Peterson
  !! https://doi.org/10.1002/fld.1964
  !!
  !! $f(x) = f_{min} + (f_{max} - f_{min}) x \frac{q + 1}{q + x}$
  !!
  !! It seems very similar to RAMP but with the convexity the other way
  !!
  !!  |       ...
  !!  |    ..
  !!  |  .
  !!  | .
  !!  |.
  !!  |_________
  
  type, public, extends(mapping_t) :: RAMP_mapping_t
  !> minimum value
  real(kind=rp) :: f_min
  !> maximum value
  real(kind=rp) :: f_max
  !> penalty parameter
  real(kind=rp) :: q
  !> Convexity of the mapping (with lower being the standard RAMP and 
  !! upper being that used by Borrvall & Peterson)
  logical :: convex_up

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
     procedure, pass(this) :: apply_backward => &
     RAMP_mapping_apply_backward
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
    this%convex_up = .true.

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

    if (this%convex_up .eqv. .true.) then
       call convex_up_RAMP_mapping_apply(this%f_min, this%f_max, &
       this%q, X_out, X_in)
    else
       call convex_down_RAMP_mapping_apply(this%f_min, this%f_max, &
       this%q, X_out, X_in)
    end if

    ! TODO
    ! memcopy or GPU backend

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

    if (this%convex_up .eqv. .true.) then
       call convex_up_RAMP_mapping_apply_backward(this%f_min, this%f_max, &
            this%q, dF_dX_in, dF_dX_out, X_in)
    else
       call convex_down_RAMP_mapping_apply_backward(this%f_min, this%f_max, &
            this%q, dF_dX_in, dF_dX_out, X_in)
    end if

    ! TODO
    ! memcopy or GPU backend

  end subroutine RAMP_mapping_apply_backward

  !> Apply the mapping
  !! @param X_out mapped field
  !! @param X_in unmapped field
  subroutine convex_down_RAMP_mapping_apply(f_min, f_max, q, X_out, X_in)
    real(kind=rp), intent(in) :: q, f_min, f_max
    type(field_t), intent(in) ::  X_in
    type(field_t), intent(inout) ::  X_out
    integer :: n, i

    ! x_out = f_min + (f_max - f_min) * x_in / (1 + q * (1 - x_in) )
    
    n = X_in%dof%size()
    do i = 1, n
       X_out%x(i,1,1,1) = f_min + (f_max - f_min) * &
       X_in%x(i,1,1,1) / (1.0_rp + q * (1.0_rp - X_in%x(i,1,1,1) ) )
    end do

  end subroutine convex_down_RAMP_mapping_apply


  !> Apply the  chain rule
  !! @param X_in unmapped field
  !! @param dF_dX_in is the sensitivity with respect to the unfiltered design
  !! @param dF_dX_out is the sensitivity with respect to the filtered design
  subroutine convex_down_RAMP_mapping_apply_backward(f_min, f_max, q, &
  dF_dX_in, dF_dX_out, X_in)
    real(kind=rp), intent(in) :: f_min, f_max, q
    type(field_t), intent(in) ::  X_in
    type(field_t), intent(in) ::  dF_dX_out
    type(field_t), intent(inout) ::  dF_dX_in
    integer :: n, i

    ! df/dx_in = df/dx_out * dx_out/dx_in
    
    ! dx_out/dx_in = (f_min - f_max) * (q + 1) / (1 - q*(x - 1))**2
    
    n = X_in%dof%size()
    do i = 1, n
       dF_dX_in%x(i,1,1,1) = (f_max - f_min) * (q + 1.0_rp) / &
       ((1.0_rp - q * (X_in%x(i,1,1,1) - 1.0_rp))**2) * &
       dF_dX_out%x(i,1,1,1)
    end do

  end subroutine convex_down_RAMP_mapping_apply_backward

  !> Apply the mapping
  !! @param X_out mapped field
  !! @param X_in unmapped field
  subroutine convex_up_RAMP_mapping_apply(f_min, f_max, q, X_out, X_in)
    real(kind=rp), intent(in) :: f_min, f_max, q
    type(field_t), intent(in) ::  X_in
    type(field_t), intent(inout) ::  X_out
    integer :: n, i

    ! x_out = f_min + (f_max - f_min) * x_in * (q + 1) / (x_in + q)
    n = X_in%dof%size()
    do i = 1, n
       X_out%x(i,1,1,1) = f_min + (f_max - f_min) * &
       X_in%x(i,1,1,1) * (1.0_rp + q ) / (X_in%x(i,1,1,1) + q) 
    end do
    

  end subroutine convex_up_RAMP_mapping_apply


  !> Apply the  chain rule
  !! @param X_in unmapped field
  !! @param dF_dX_in is the sensitivity with respect to the unfiltered design
  !! @param dF_dX_out is the sensitivity with respect to the filtered design
  subroutine convex_up_RAMP_mapping_apply_backward(f_min, f_max, q, &
  dF_dX_in, dF_dX_out, X_in)
    real(kind=rp), intent(in) :: f_min, f_max, q
    type(field_t), intent(in) ::  X_in
    type(field_t), intent(in) ::  dF_dX_out
    type(field_t), intent(inout) ::  dF_dX_in
    integer :: n, i

    ! df/dx_in = df/dx_out * dx_out/dx_in
    
    ! dx_out/dx_in = (f_min - f_max) * (q + 1) / (q + x)**2

    n = X_in%dof%size()
    do i = 1, n
       dF_dX_in%x(i,1,1,1) = (f_max - f_min) * (q + 1.0_rp) / &
       ( (X_in%x(i,1,1,1) + q)**2) * &
       dF_dX_out%x(i,1,1,1)
    end do

  end subroutine convex_up_RAMP_mapping_apply_backward
end module RAMP_mapping
