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
!> This provides a template to be modified in constructing new source terms.
!> For a detailed description see ///
!
!------------------------------------------------------------------------------
! TO BE FILLED: a description of the source term
!------------------------------------------------------------------------------
!
module adjoint_TEMPLATE_source_term
  use num_types, only : rp
  use field_list, only : field_list_t
  use json_module, only : json_file
  use source_term, only : source_term_t
  use coefs, only : coef_t
  use field, only: field_t
  !----------------------------------------------------------------------------
  ! TO BE FILLED: Add additional modules
  !----------------------------------------------------------------------------
  implicit none
  private

  !----------------------------------------------------------------------------
  !> TO BE FILLED: a description of the source term
  !----------------------------------------------------------------------------
  type, public, extends(source_term_t) :: adjoint_TEMPLATE_source_term_t

     !> A mask for where the source term is evaluated
     class(point_zone_t), pointer :: mask
     !> containing a mask?
     logical :: if_mask
     !-------------------------------------------------------------------------
     ! TO BE FILLED: additional private variables used by you source term
     !-------------------------------------------------------------------------
   contains
     !> The common constructor using a JSON object.
     procedure, pass(this) :: init => &
          adjoint_TEMPLATE_source_term_init_from_json
     !> The constructor from type components.
     procedure, pass(this) :: init_from_components => &
          adjoint_TEMPLATE_source_term_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => adjoint_TEMPLATE_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute_ => adjoint_TEMPLATE_source_term_compute
  end type adjoint_TEMPLATE_source_term_t

contains
  !> The common constructor using a JSON object.
  !! @param json The JSON object for the source.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  subroutine adjoint_TEMPLATE_source_term_init_from_json(this, json, fields, &
       coef)
    class(adjoint_TEMPLATE_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef



  end subroutine adjoint_TEMPLATE_source_term_init_from_json

  !> The constructor from type components.
  !----------------------------------------------------------------------------
  ! TO BE FILLED: Document your additional parameters
  !----------------------------------------------------------------------------
  subroutine adjoint_TEMPLATE_source_term_init_from_components(this, mask, &
       if_mask, coef)
    class(adjoint_TEMPLATE_source_term_t), intent(inout) :: this
    type(field_list_t) :: fields
    type(coef_t) :: coef
    real(kind=rp) :: start_time
    real(kind=rp) :: end_time
    type(topopt_design_t), intent(in), target :: design
    class(point_zone_t), intent(in), target :: mask
    logical :: if_mask
    real(kind=rp) :: K
    type(field_t), intent(in), target :: u, v, w

    !--------------------------------------------------------------------------
    ! TO BE FILLED: Associate the RHS fields you aim to influence
    !--------------------------------------------------------------------------
    ! eg,

    ! call fields%init(1)
    ! call fields%assign(1, f_x)

    ! assume adjoint source term acts at all times
    start_time = 0.0_rp
    end_time = 100000000.0_rp

    ! initialize base
    call this%free()
    call this%init_base(fields, coef, start_time, end_time)


    !--------------------------------------------------------------------------
    ! TO BE FILLED: Assign pointers to internal variables
    !--------------------------------------------------------------------------

    ! Associate the mask
    this%if_mask = if_mask
    if (this%if_mask) then
       this%mask => mask
    end if

  end subroutine adjoint_TEMPLATE_source_term_init_from_components

  !> Computes the source term and adds the result to `fields`.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine adjoint_TEMPLATE_source_term_compute(this, t, tstep)
    class(adjoint_TEMPLATE_source_term_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    !--------------------------------------------------------------------------
    ! TO BE FILLED: Compute adjoint source terms
    !--------------------------------------------------------------------------

  end subroutine adjoint_TEMPLATE_source_term_compute

  !> Destructor.
  subroutine adjoint_TEMPLATE_source_term_free(this)
    class(adjoint_TEMPLATE_source_term_t), intent(inout) :: this

    !--------------------------------------------------------------------------
    ! TO BE FILLED: Free additional objects
    !--------------------------------------------------------------------------

    call this%free_base()

  end subroutine adjoint_TEMPLATE_source_term_free
end module adjoint_TEMPLATE_source_term
