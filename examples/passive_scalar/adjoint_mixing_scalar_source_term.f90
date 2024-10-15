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
!> Implements the `adjoint_mixing_scalar_source_term` type.
! this is a such a dumb name
module adjoint_mixing_scalar_source_term
  use num_types, only : rp
  use field_list, only : field_list_t
  use field, only: field_t
  use field_registry, only: neko_field_registry
  use scratch_registry, only: neko_scratch_registry
  use json_module, only : json_file
  use json_utils, only: json_get, json_get_or_default
  use source_term, only : source_term_t
  use coefs, only : coef_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use utils, only : neko_error
  use field_math, only: field_addcol3, field_add2, field_cadd, field_col2, &
  field_add2s2, field_copy
  use operators, only: opgrad
  implicit none
  private

  ! this will be the adjoint forcing from Casper's objective function
  ! TODO
  ! it actually isn't this. Infact, his is a surface integral on the outlet
  ! Which means we get strange BC's in our adjoint problem.
  ! This source term would be if we had a certain volume that we wanted more 
  ! mixed
  ! the forcing is of the form:
  ! $ s - \bar{s} $
  ! ie, difference between it and the average.
  ! I'm going to hardcode 0.5...
  ! because we'll never actually use this source term I don't think
  type, public, extends(source_term_t) :: &
  adjoint_mixing_scalar_source_term_t
     !> here we have s (coming from the primal)
     type(field_t), pointer :: s
     !> A scalaing factor
     real(kind=rp) :: obj_scale
   contains
     !> The common constructor using a JSON object.
     procedure, pass(this) :: init => &
     adjoint_mixing_scalar_source_term_init_from_json
     !> The constructor from type components.
     procedure, pass(this) :: init_from_components => &
       adjoint_mixing_scalar_source_term_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => adjoint_mixing_scalar_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute_ => &
     adjoint_mixing_scalar_source_term_compute
  end type adjoint_mixing_scalar_source_term_t

contains
  !> The common constructor using a JSON object.
  !! @param json The JSON object for the source.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  subroutine adjoint_mixing_scalar_source_term_init_from_json(this, &
  json, fields, coef)
    class(adjoint_mixing_scalar_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(inout), target :: fields
    type(coef_t), intent(inout), target :: coef
    real(kind=rp), allocatable :: values(:)
    real(kind=rp) :: start_time, end_time


  end subroutine adjoint_mixing_scalar_source_term_init_from_json


  subroutine adjoint_mixing_scalar_source_term_init_from_components(this,&
       f_s, s, obj_scale, coef)
    class(adjoint_mixing_scalar_source_term_t), intent(inout) :: this
    type(field_t), pointer, intent(in) :: f_s
    type(field_list_t) :: fields
    type(coef_t) :: coef
    real(kind=rp) :: start_time
    real(kind=rp) :: end_time
    real(kind=rp) :: obj_scale
    type(field_t), intent(in), target :: s
    ! TODo
    ! do masks later
    !type(field_t), intent(in), target :: mask

    ! I wish you didn't need a start time and end time...
    ! but I'm just going to set a super big number...
    start_time = 0.0_rp
    end_time = 100000000.0_rp

    call this%free()

    ! this is copying the fluid source term init
    ! We package the fields for the source term to operate on in a field list.
    call fields%init(1)
    call fields%assign(1, f_s)

    call this%init_base(fields, coef, start_time, end_time)

    ! point everything in the correct places
    this%s => s
    this%obj_scale = obj_scale

    print *, "source term initialized"


    ! TODO
    !this%mask => mask
  end subroutine adjoint_mixing_scalar_source_term_init_from_components

  !> Destructor.
  subroutine adjoint_mixing_scalar_source_term_free(this)
    class(adjoint_mixing_scalar_source_term_t), intent(inout) :: this

    call this%free_base()
  end subroutine adjoint_mixing_scalar_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine adjoint_mixing_scalar_source_term_compute(this, t, tstep)
    class(adjoint_mixing_scalar_source_term_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    integer :: n_fields, i, n
    type(field_t), pointer :: fs

    
    fs => this%fields%get(1)

    ! TODO
    ! double check if add or subtract

    ! TODO
    ! Normalize

    ! TODO
    ! Masks

    call field_add2s2(fs,this%s, this%obj_scale)

    ! TODO
    ! I'm just subtracting 0.5, but it should be the average
    call field_cadd(fs, -0.5_rp*this%obj_scale)
    
  end subroutine adjoint_mixing_scalar_source_term_compute

end module adjoint_mixing_scalar_source_term
