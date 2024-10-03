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
!> Implements the `simple_brinkman_source_term_t` type.
module simple_brinkman_source_term
  use num_types, only : rp
  use field_list, only : field_list_t
  use json_module, only : json_file
  use json_utils, only: json_get, json_get_or_default
  use source_term, only : source_term_t
  use coefs, only : coef_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use utils, only : neko_error
  use field, only: field_t
  use new_design, only: new_design_t
  use field_math, only: field_subcol3
  implicit none
  private

  !> A constant source term.
  !! The strength is specified with the `values` keyword, which should be an
  !! array, with a value for each component of the source.
  type, public, extends(source_term_t) :: simple_brinkman_source_term_t
     ! We have a source term of the form $\chi u$
     !> the fields corresponding to \chi, u, v and w
     type(field_t), pointer :: chi, u, v, w

   contains
     !> The common constructor using a JSON object.
     procedure, pass(this) :: init => &
          simple_brinkman_source_term_init_from_json
     !> The constructor from type components.
     procedure, pass(this) :: init_from_components => &
          simple_brinkman_source_term_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => simple_brinkman_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute_ => simple_brinkman_source_term_compute
  end type simple_brinkman_source_term_t

contains
  !> The common constructor using a JSON object.
  !! @param json The JSON object for the source.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  subroutine simple_brinkman_source_term_init_from_json(this, json, fields, &
       coef)
    class(simple_brinkman_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(inout), target :: fields
    type(coef_t), intent(inout), target :: coef
    real(kind=rp), allocatable :: values(:)
    real(kind=rp) :: start_time, end_time


    ! we shouldn't be initializing this from JSON
    ! maybe throw an error?


  end subroutine simple_brinkman_source_term_init_from_json

  !> The constructor from type components.
  !! @param fields A list of fields for adding the source values.
  !! @param values The array of values, one for each field.
  !! @param coef The SEM coeffs.
  !! @param start_time When to start adding the source term.
  !! @param end_time When to stop adding the source term.
  subroutine simple_brinkman_source_term_init_from_components(this, &
       f_x, f_y, f_z, design, u, v, w, coef)
    class(simple_brinkman_source_term_t), intent(inout) :: this
    type(field_t), pointer, intent(in) :: f_x, f_y, f_z
    type(field_list_t) :: fields
    type(coef_t) :: coef
    real(kind=rp) :: start_time
    real(kind=rp) :: end_time
    type(field_t), intent(in), target :: u, v, w
    type(new_design_t), intent(in), target :: design

    ! I wish you didn't need a start time and end time...
    ! Tim you're going to hate this...
    ! but I'm just going to set a super big number...
    start_time = 0.0_rp
    end_time = 100000000.0_rp

    call this%free()

    ! this is copying the fluid source term init
    ! We package the fields for the source term to operate on in a field list.
    call fields%init(3)
    call fields%assign(1, f_x)
    call fields%assign(2, f_y)
    call fields%assign(3, f_z)



    call this%init_base(fields, coef, start_time, end_time)

    ! Real stuff

    ! point everything in the correct places
    this%u => u
    this%v => v
    this%w => w
    ! and get chi out of the design
    this%chi => design%brinkman_amplitude





  end subroutine simple_brinkman_source_term_init_from_components

  !> Destructor.
  subroutine simple_brinkman_source_term_free(this)
    class(simple_brinkman_source_term_t), intent(inout) :: this

    call this%free_base()
  end subroutine simple_brinkman_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine simple_brinkman_source_term_compute(this, t, tstep)
    class(simple_brinkman_source_term_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(field_t), pointer :: fu, fv, fw
    integer :: n

    fu => this%fields%get_by_index(1)
    fv => this%fields%get_by_index(2)
    fw => this%fields%get_by_index(3)

    call field_subcol3(fu, this%u, this%chi)
    call field_subcol3(fv, this%v, this%chi)
    call field_subcol3(fw, this%w, this%chi)


  end subroutine simple_brinkman_source_term_compute

end module simple_brinkman_source_term
