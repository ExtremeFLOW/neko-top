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
!> Implements the `adjoint_scalar_convection_source_term` type.
! this is a such a dumb name
module adjoint_scalar_convection_source_term
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
  use field_math, only: field_addcol3
  use operators, only: opgrad
  implicit none
  private

  ! I don't know how to name this term, but when you have a passive
  ! scalar you get an extra term in the adjoint velocity equation, which comes
  ! from the convective term in the passive scalar equation.
  ! Maybe this should be called just `adjoint_scalar_convection`? 
  ! but it does come in a source term...

  ! In any case,
  ! it's a source term acting on the adjoint velocity equations, of the form:
  ! $\nabla s s_adj$
  type, public, extends(source_term_t) :: &
  adjoint_scalar_convection_source_term_t
     !> here we have s and s_adj
     type(field_t), pointer :: s, s_adj
   contains
     !> The common constructor using a JSON object.
     procedure, pass(this) :: init => &
     adjoint_scalar_convection_source_term_init_from_json
     !> The constructor from type components.
     procedure, pass(this) :: init_from_components => &
       adjoint_scalar_convection_source_term_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => adjoint_scalar_convection_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute_ => &
     adjoint_scalar_convection_source_term_compute
  end type adjoint_scalar_convection_source_term_t

contains
  !> The common constructor using a JSON object.
  !! @param json The JSON object for the source.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  subroutine adjoint_scalar_convection_source_term_init_from_json(this, &
  json, fields, coef)
    class(adjoint_scalar_convection_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(inout), target :: fields
    type(coef_t), intent(inout), target :: coef
    real(kind=rp), allocatable :: values(:)
    real(kind=rp) :: start_time, end_time

    ! this is a bit weird... because I don't think this should come from the 
    ! JSON. 
    ! Maybe we should think of all these source terms as only "appendable"
    !
    ! Because we'll never have the whole case here, so we'll never be able
    ! init from components anyway...


  end subroutine adjoint_scalar_convection_source_term_init_from_json


  subroutine adjoint_scalar_convection_source_term_init_from_components(this,&
       f_x, f_y, f_z, &
       s, s_adj, coef)
    class(adjoint_scalar_convection_source_term_t), intent(inout) :: this
    type(field_t), pointer, intent(in) :: f_x, f_y, f_z
    type(field_list_t) :: fields
    type(coef_t) :: coef
    real(kind=rp) :: start_time
    real(kind=rp) :: end_time
    real(kind=rp) :: obj_scale
    type(field_t), intent(in), target :: s, s_adj
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
    this%s => s
    this%s_adj => s_adj


    ! TODO
    !this%mask => mask
  end subroutine adjoint_scalar_convection_source_term_init_from_components

  !> Destructor.
  subroutine adjoint_scalar_convection_source_term_free(this)
    class(adjoint_scalar_convection_source_term_t), intent(inout) :: this

    call this%free_base()
  end subroutine adjoint_scalar_convection_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine adjoint_scalar_convection_source_term_compute(this, t, tstep)
    class(adjoint_scalar_convection_source_term_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    integer :: n_fields, i, n
    type(field_t), pointer :: fu, fv, fw
    integer :: temp_indices(3)
    type(field_t), pointer :: dsdx, dsdy, dsdz

    
    call neko_scratch_registry%request_field(dsdx, temp_indices(1))
    call neko_scratch_registry%request_field(dsdy, temp_indices(2))
    call neko_scratch_registry%request_field(dsdz, temp_indices(3))
    
    fu => this%fields%get(1)
    fv => this%fields%get(2)
    fw => this%fields%get(3)
    
    
    
    ! we basically just need the term 
    ! $\nabla s s_adj$
    ! TODO
    ! be super careful with opgrad vs grad.
    ! I'm not sure which is correct here
    ! From memory they only differ by a jac_inv or a B or something like that
    ! TODO
    ! I think this actually works on GPU but I haven't checked..
    call opgrad(dsdx%x,dsdy%x,dsdz%x,this%s%x, this%coef)
    ! TODO
    ! double check if add or subtract
    call field_addcol3(fu,this%s_adj,dsdx)
    call field_addcol3(fv,this%s_adj,dsdy)
    call field_addcol3(fw,this%s_adj,dsdz)
    
    ! free the scratch
    call neko_scratch_registry%relinquish_field(temp_indices)
  end subroutine adjoint_scalar_convection_source_term_compute

end module adjoint_scalar_convection_source_term
