!> @file Borrvall_Peterson_mapping.f90
!! @copyright
!! Copyright (c) 2024-2026, The Neko-TOP Authors
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
!> A Borrvall & Peterson mapping of coefficients
module Borrvall_Peterson_mapping
  use num_types, only: rp
  use mapping, only: mapping_t
  use json_module, only: json_file
  use field, only: field_t
  use coefs, only: coef_t
  use neko_config, only: NEKO_BCKND_DEVICE
  use device_Borrvall_Peterson_mapping, only: &
       device_Borrvall_Peterson_mapping_apply, &
       device_Borrvall_Peterson_mapping_apply_backward
  use Borrvall_Peterson_mapping_cpu, only: &
       Borrvall_Peterson_mapping_apply_cpu, &
       Borrvall_Peterson_mapping_apply_backward_cpu
  use json_utils, only: json_get, json_get_or_default
  use logger, only: neko_log
  use continuation_scheduler, only: nekotop_continuation
  implicit none
  private

  !> A Borrvall & Peterson mapping of coefficients
  !! This is the material interpolation described by Borrvall & Peterson
  !! https://doi.org/10.1002/fld.1964
  !!
  !! It reuses the shape function
  !! \f$S(x) = x \frac{q + 1}{x + q}\f$ (with \f$S(0)=0\f$, \f$S(1)=1\f$),
  !! assembled so that the mapping is decreasing:
  !!
  !! \f$f(x) = f_{max} + (f_{min} - f_{max}) S(x)\f$
  !!
  !! with \f$f(0) = f_{max}\f$ (solid) and \f$f(1) = f_{min}\f$ (fluid),
  !! steepest near \f$x \to 0\f$, for the convention x=0: solid, x=1: fluid.
  !! This is the converse of the ordinary RAMP mapping.
  !!
  !!  |.
  !!  | .
  !!  |  .
  !!  |    ..
  !!  |       ...
  !!  |_________
  !!
  !! The shape function is identical to the historic RAMP "convex up" form;
  !! only the assembly of f_min/f_max is corrected. This is realised by
  !! calling the (unmodified) kernel entry points with f_min and f_max
  !! swapped at the call site (see forward/backward mapping below).

  type, public, extends(mapping_t) :: Borrvall_Peterson_mapping_t
     !> minimum value (the fluid-side value)
     real(kind=rp) :: f_min
     !> maximum value (the solid-side value)
     real(kind=rp) :: f_max
     !> penalty parameter
     real(kind=rp) :: q

   contains
     !> Constructor from json.
     procedure, pass(this) :: init => Borrvall_Peterson_mapping_init_from_json
     !> Actual constructor.
     procedure, pass(this) :: init_from_attributes => &
          Borrvall_Peterson_mapping_init_from_attributes
     !> Destructor.
     procedure, pass(this) :: free => Borrvall_Peterson_mapping_free
     !> Apply the forward mapping
     procedure, pass(this) :: forward_mapping => &
          Borrvall_Peterson_forward_mapping
     !> Apply the adjoint mapping
     procedure, pass(this) :: backward_mapping => &
          Borrvall_Peterson_backward_mapping
  end type Borrvall_Peterson_mapping_t

contains

  !> Constructor from json.
  subroutine Borrvall_Peterson_mapping_init_from_json(this, json, coef)
    class(Borrvall_Peterson_mapping_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(coef_t), intent(inout) :: coef
    real(kind=rp) :: f_min, f_max, q

    call json_get_or_default(json, 'f_min', f_min, 0.0_rp)
    call nekotop_continuation%json_get_or_register(json, 'f_max', this%f_max, &
         f_max)
    call nekotop_continuation%json_get_or_register(json, 'q', this%q, q, 1.0_rp)

    call this%init_base(json, coef, "Borrvall_Peterson_mapping")
    call this%init_from_attributes(coef, f_min, f_max, q)

  end subroutine Borrvall_Peterson_mapping_init_from_json

  !> Actual constructor.
  subroutine Borrvall_Peterson_mapping_init_from_attributes(this, coef, &
       f_min, f_max, q)
    class(Borrvall_Peterson_mapping_t), intent(inout) :: this
    type(coef_t), intent(inout) :: coef
    real(kind=rp), intent(in) :: f_min, f_max, q
    character(len=256) :: msg

    this%f_min = f_min
    this%f_max = f_max
    this%q = q

    call neko_log%section('Borrvall & Peterson Mapping')
    write(msg, '(A,F8.4)') '  f_min: ', this%f_min
    call neko_log%message(msg)
    write(msg, '(A,F8.4)') '  f_max: ', this%f_max
    call neko_log%message(msg)
    write(msg, '(A,F8.4)') '  q:     ', this%q
    call neko_log%message(msg)
    call neko_log%end_section()

  end subroutine Borrvall_Peterson_mapping_init_from_attributes

  !> Destructor.
  subroutine Borrvall_Peterson_mapping_free(this)
    class(Borrvall_Peterson_mapping_t), intent(inout) :: this

    call this%free_base()

  end subroutine Borrvall_Peterson_mapping_free

  !> Apply the mapping
  !! @param this mapping
  !! @param X_out mapped field
  !! @param X_in unmapped field
  subroutine Borrvall_Peterson_forward_mapping(this, X_out, X_in)
    class(Borrvall_Peterson_mapping_t), intent(inout) :: this
    type(field_t), intent(in) :: X_in
    type(field_t), intent(inout) :: X_out
    integer :: n

    ! The moved (unmodified) kernel computes
    !   arg1 + (arg2 - arg1) * S(x)
    ! with S(x) = x (q+1) / (x+q). Passing f_max as arg1 and f_min as arg2
    ! (i.e. SWAPPED with respect to the historic RAMP call) yields the correct
    ! Borrvall & Peterson assembly
    !   f(x) = f_max + (f_min - f_max) * S(x).
    n = X_in%dof%size()
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_Borrvall_Peterson_mapping_apply(this%f_max, this%f_min, &
            this%q, X_out%x_d, X_in%x_d, n)
    else
       call Borrvall_Peterson_mapping_apply_cpu(this%f_max, this%f_min, &
            this%q, X_out%x, X_in%x, n)
    end if

  end subroutine Borrvall_Peterson_forward_mapping


  !> Apply the  chain rule
  !! @param this mapping
  !! @param sens_out is the sensitivity with respect to the unfiltered design
  !! @param sens_in is the sensitivity with respect to the filtered design
  !! @param X_in unmapped field
  subroutine Borrvall_Peterson_backward_mapping(this, sens_out, sens_in, X_in)
    class(Borrvall_Peterson_mapping_t), intent(inout) :: this
    type(field_t), intent(in) :: X_in
    type(field_t), intent(in) :: sens_in
    type(field_t), intent(inout) :: sens_out
    integer :: n

    ! Same SWAP of f_min/f_max as the forward mapping. The moved kernel gives
    !   dx_out/dx_in = (arg2 - arg1) * (q+1) * q / (x+q)**2
    ! which with arg1 = f_max, arg2 = f_min becomes
    !   f'(x) = q (q+1) (f_min - f_max) / (x+q)**2.
    n = X_in%dof%size()
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_Borrvall_Peterson_mapping_apply_backward(this%f_max, &
            this%f_min, this%q, sens_out%x_d, sens_in%x_d, X_in%x_d, n)
    else
       call Borrvall_Peterson_mapping_apply_backward_cpu(this%f_max, &
            this%f_min, this%q, sens_out%x, sens_in%x, X_in%x, n)
    end if

  end subroutine Borrvall_Peterson_backward_mapping

end module Borrvall_Peterson_mapping
