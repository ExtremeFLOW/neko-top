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
!> Implements the `adjoint_lube_source_term_t` type.
!
!
! I know this is a stupid naming convention...
! The `lube` aspect came from a paper that attributed this term to out of plane
! stresses based on lubrication theory.
!
! I preffer to think of it as a constraint that penalizes non-binary designs
!
! The term is $K \int_\Omega \frac{1}{2}\chi|\mathbf{u}|^2$
!
! the corresponding adjoint forcing is $K \chi \mathbf{u}$
module adjoint_lube_source_term
  use num_types, only : rp
  use field_list, only : field_list_t
  use json_module, only : json_file
  use source_term, only : source_term_t
  use coefs, only : coef_t
  use field, only: field_t
  use topopt_design, only: topopt_design_t
  use field_math, only: field_addcol3, field_copy, field_cmult
  use scratch_registry, only : neko_scratch_registry
  implicit none
  private

  !> A adjoint source term corresponding to an objective of
  ! $K \int_\Omega \frac{1}{2}\chi|\mathbf{u}|^2$.
  type, public, extends(source_term_t) :: adjoint_lube_source_term_t

     !> $u,v,w$ corresponding to the baseflow
     type(field_t), pointer :: u,v,w
     !> $\chi$ the Brinkman amplitude
     type(field_t), pointer :: chi
     ! TODO
     ! as mask
     ! type(field_t), pointer :: chi
     !> a scale for this term
     real(kind=rp) :: K

   contains
     !> The common constructor using a JSON object.
     procedure, pass(this) :: init => adjoint_lube_source_term_init_from_json
     !> The constructor from type components.
     procedure, pass(this) :: init_from_components => &
          adjoint_lube_source_term_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => adjoint_lube_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute_ => adjoint_lube_source_term_compute
  end type adjoint_lube_source_term_t

contains
  !> The common constructor using a JSON object.
  !! @param json The JSON object for the source.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  subroutine adjoint_lube_source_term_init_from_json(this, json, fields, coef)
    class(adjoint_lube_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(inout), target :: fields
    type(coef_t), intent(inout), target :: coef
    ! real(kind=rp), allocatable :: values(:)
    ! real(kind=rp) :: start_time, end_time


    ! we shouldn't be initializing this from JSON
    ! maybe throw an error?


  end subroutine adjoint_lube_source_term_init_from_json

  !> The constructor from type components.
  !! @param f_x, f_y, f_z the RHS of the adjoint
  !! @param design the design
  !! @param K a scale
  !! @param u, v, w the velocity fields of the primal
  ! $u,v,w$ reffer to the primal, not the adjoint
  subroutine adjoint_lube_source_term_init_from_components(this, &
       f_x, f_y, f_z, design, K, &
       u, v, w, coef)
    class(adjoint_lube_source_term_t), intent(inout) :: this
    type(field_t), pointer, intent(in) :: f_x, f_y, f_z
    type(field_list_t) :: fields
    type(coef_t) :: coef
    real(kind=rp) :: start_time
    real(kind=rp) :: end_time
    type(topopt_design_t), intent(in), target :: design
    real(kind=rp) :: K
    type(field_t), intent(in), target :: u, v, w
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
    call fields%init(3)
    call fields%assign(1, f_x)
    call fields%assign(2, f_y)
    call fields%assign(3, f_z)

    call this%init_base(fields, coef, start_time, end_time)

    ! point everything in the correct places
    ! NOTE!!!
    ! this is the primal!
    this%u => u
    this%v => v
    this%w => w

    this%chi => design%brinkman_amplitude
    this%K = K

    ! TODO
    !this%mask => mask

  end subroutine adjoint_lube_source_term_init_from_components

  !> Destructor.
  subroutine adjoint_lube_source_term_free(this)
    class(adjoint_lube_source_term_t), intent(inout) :: this

    call this%free_base()
  end subroutine adjoint_lube_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine adjoint_lube_source_term_compute(this, t, tstep)
    class(adjoint_lube_source_term_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(field_t), pointer :: fu, fv, fw
    type(field_t), pointer :: work
    integer :: temp_indices(1)

    fu => this%fields%get_by_index(1)
    fv => this%fields%get_by_index(2)
    fw => this%fields%get_by_index(3)

    ! BE SO CAREFUL HERE
    ! It make look the same as the Brinkman term, but it's assumed that
    ! this source term acts on the adjoint, and the u,v,w here come from
    ! the primal
    call neko_scratch_registry%request_field(work, temp_indices(1))
    call field_copy(work, this%chi)

    ! scale by K
    call field_cmult(work, this%K)

    ! multiple and add the RHS
    call field_addcol3(fu, this%u, work)
    call field_addcol3(fv, this%v, work)
    call field_addcol3(fw, this%w, work)
    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine adjoint_lube_source_term_compute

end module adjoint_lube_source_term
