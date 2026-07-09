!> @file RAMP_mapping.f90
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
!> A RAMP mapping of coefficients
module RAMP_mapping
  use num_types, only: rp
  use mapping, only: mapping_t
  use json_module, only: json_file
  use field, only: field_t
  use coefs, only: coef_t
  use neko_config, only: NEKO_BCKND_DEVICE
  use device_RAMP_mapping, only: device_convex_down_RAMP_mapping_apply, &
       device_convex_down_RAMP_mapping_apply_backward
  use RAMP_mapping_cpu, only: convex_down_RAMP_mapping_apply_cpu, &
       convex_down_RAMP_mapping_apply_backward_cpu
  use json_utils, only: json_get, json_get_or_default
  use logger, only: neko_log
  use utils, only: neko_error
  use continuation_scheduler, only: nekotop_continuation
  implicit none
  private

  !> A RAMP mapping of coefficients
  !! This is the standard RAMP described in
  !! https://doi.org/10.1007/s001580100129
  !!
  !! \f$f(x) = f_{min} + (f_{max} - f_{min}) \frac{x}{1 + q(1 - x)}\f$
  !!
  !! This is increasing, with \f$f(0) = f_{min}\f$ (fluid) and
  !! \f$f(1) = f_{max}\f$ (solid), for the convention x=0: fluid, x=1: solid.
  !!
  !!  |        .
  !!  |        .
  !!  |       .
  !!  |     ..
  !!  |  ...
  !!  |_________
  !!
  !! For the converse convention (x=0: solid, x=1: fluid) use the
  !! Borrvall & Petersson mapping instead.

  type, public, extends(mapping_t) :: RAMP_mapping_t
     !> minimum value (the fluid-side value)
     real(kind=rp) :: f_min
     !> maximum value (the solid-side value)
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
     procedure, pass(this) :: forward_mapping => RAMP_forward_mapping
     !> Apply the adjoint mapping
     procedure, pass(this) :: backward_mapping => RAMP_backward_mapping
  end type RAMP_mapping_t

contains

  !> Constructor from json.
  subroutine RAMP_mapping_init_from_json(this, json, coef)
    class(RAMP_mapping_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(coef_t), intent(inout) :: coef
    real(kind=rp) :: f_min, f_max, q

    ! The 'convex_up' option has been removed. Fail loud rather than silently
    ! reinterpreting old case files.
    if (json%valid_path('convex_up')) then
       call neko_error("RAMP mapping's 'convex_up' option has been " // &
            "removed -- plain RAMP is now always the down/increasing form " // &
            "(x=0:fluid, x=1:solid). For the x=0:solid/x=1:fluid form, " // &
            "use type 'Borrvall_Petersson' instead.")
    end if

    call json_get_or_default(json, 'f_min', f_min, 0.0_rp)
    call nekotop_continuation%json_get_or_register(json, 'f_max', this%f_max, &
         f_max)
    call nekotop_continuation%json_get_or_register(json, 'q', this%q, q, 1.0_rp)

    call this%init_base(json, coef, "RAMP_mapping")
    call this%init_from_attributes(coef, f_min, f_max, q)

  end subroutine RAMP_mapping_init_from_json

  !> Actual constructor.
  subroutine RAMP_mapping_init_from_attributes(this, coef, f_min, f_max, q)
    class(RAMP_mapping_t), intent(inout) :: this
    type(coef_t), intent(inout) :: coef
    real(kind=rp), intent(in) :: f_min, f_max, q
    character(len=256) :: msg

    this%f_min = f_min
    this%f_max = f_max
    this%q = q

    call neko_log%section('RAMP Mapping')
    write(msg, '(A,F8.4)') '  f_min: ', this%f_min
    call neko_log%message(msg)
    write(msg, '(A,F8.4)') '  f_max: ', this%f_max
    call neko_log%message(msg)
    write(msg, '(A,F8.4)') '  q:     ', this%q
    call neko_log%message(msg)
    call neko_log%end_section()

  end subroutine RAMP_mapping_init_from_attributes

  !> Destructor.
  subroutine RAMP_mapping_free(this)
    class(RAMP_mapping_t), intent(inout) :: this

    call this%free_base()

  end subroutine RAMP_mapping_free

  !> Apply the mapping
  !! @param this mapping
  !! @param X_out mapped field
  !! @param X_in unmapped field
  subroutine RAMP_forward_mapping(this, X_out, X_in)
    class(RAMP_mapping_t), intent(inout) :: this
    type(field_t), intent(in) :: X_in
    type(field_t), intent(inout) :: X_out

    call convex_down_RAMP_mapping_apply(this%f_min, this%f_max, &
         this%q, X_out, X_in)

  end subroutine RAMP_forward_mapping


  !> Apply the  chain rule
  !! @param this mapping
  !! @param sens_out is the sensitivity with respect to the unfiltered design
  !! @param sens_in is the sensitivity with respect to the filtered design
  !! @param X_in unmapped field
  subroutine RAMP_backward_mapping(this, sens_out, sens_in, X_in)
    class(RAMP_mapping_t), intent(inout) :: this
    type(field_t), intent(in) :: X_in
    type(field_t), intent(in) :: sens_in
    type(field_t), intent(inout) :: sens_out

    call convex_down_RAMP_mapping_apply_backward(this%f_min, this%f_max, &
         this%q, sens_out, sens_in, X_in)

  end subroutine RAMP_backward_mapping

  !> Apply the mapping
  !! @param f_min minimum value
  !! @param f_max maximum value
  !! @param q penalty parameter
  !! @param X_out mapped field
  !! @param X_in unmapped field
  subroutine convex_down_RAMP_mapping_apply(f_min, f_max, q, X_out, X_in)
    real(kind=rp), intent(in) :: q, f_min, f_max
    type(field_t), intent(in) :: X_in
    type(field_t), intent(inout) :: X_out
    integer :: n

    ! x_out = f_min + (f_max - f_min) * x_in / (1 + q * (1 - x_in) )

    n = X_in%dof%size()
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_convex_down_RAMP_mapping_apply(f_min, f_max, q, &
            X_out%x_d, X_in%x_d, n)
    else
       call convex_down_RAMP_mapping_apply_cpu(f_min, f_max, q, &
            X_out%x, X_in%x, n)
    end if

  end subroutine convex_down_RAMP_mapping_apply


  !> Apply the  chain rule
  !! @param f_min minimum value
  !! @param f_max maximum value
  !! @param q penalty parameter
  !! @param X_in unmapped field
  !! @param sens_out is the sensitivity with respect to the unfiltered design
  !! @param sens_in is the sensitivity with respect to the filtered design
  subroutine convex_down_RAMP_mapping_apply_backward(f_min, f_max, q, &
       sens_out, sens_in, X_in)
    real(kind=rp), intent(in) :: f_min, f_max, q
    type(field_t), intent(in) :: X_in
    type(field_t), intent(in) :: sens_in
    type(field_t), intent(inout) :: sens_out
    integer :: n

    ! df/dx_in = df/dx_out * dx_out/dx_in

    ! dx_out/dx_in = (f_max - f_min) * (q + 1) / (1 - q*(x - 1))**2

    n = X_in%dof%size()

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_convex_down_RAMP_mapping_apply_backward(f_min, f_max, q, &
            sens_out%x_d, sens_in%x_d, X_in%x_d, n)
    else
       call convex_down_RAMP_mapping_apply_backward_cpu(f_min, f_max, q, &
            sens_out%x, sens_in%x, X_in%x, n)
    end if

  end subroutine convex_down_RAMP_mapping_apply_backward

end module RAMP_mapping
