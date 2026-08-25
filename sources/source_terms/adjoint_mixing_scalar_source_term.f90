!> @file adjoint_mixing_scalar_source_term.f90
!! @copyright
!! Copyright (c) 2025-2026, The Neko-TOP Authors
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
!> Implements the `adjoint_mixing_scalar_source_term` type.
! this is a such a dumb name
module adjoint_mixing_scalar_source_term
  use num_types, only : rp
  use field_list, only : field_list_t
  use field, only: field_t
  use registry, only: neko_registry
  use scratch_registry, only: neko_scratch_registry
  use json_module, only : json_file
  use time_state, only: time_state_t
  use json_utils, only: json_get, json_get_or_default
  use source_term, only : source_term_t
  use coefs, only : coef_t
  use field_math, only: field_add2s2, field_copy, field_cadd
  use mask_ops, only: mask_exterior_const, compute_masked_volume
  use point_zone, only: point_zone_t
  use comm, only: pe_rank
  use math, only: glsc2
  use iso_fortran_env, only: error_unit
  implicit none
  private
  public :: adjoint_mixing_scalar_source_term_allocate

  ! this will be the adjoint forcing from Casper's objective function
  ! TODO
  ! it actually isn't this. Infact, his is a surface integral on the outlet
  ! Which means we get strange BC's in our adjoint problem.
  ! This source term would be if we had a certain volume that we wanted more
  ! mixed
  ! the forcing is of the form:
  ! \f$ \phi - \phi_ref \f$
  ! ie, difference between it and the average.
  type, public, extends(source_term_t) :: adjoint_mixing_scalar_source_term_t
     !> The forward scalar field
     type(field_t), pointer :: s => null()
     !> A scalaing factor
     real(kind=rp) :: obj_scale
     !> Reference concentration
     real(kind=rp) :: phi_ref
     !> A mask for where the source term is evaluated
     class(point_zone_t), pointer :: mask => null()
     !> containing a mask?
     logical :: if_mask
     !> The volume of the masked region (or whole domain)
     real(kind=rp) :: mask_volume
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
  !> Allocator for the adjoint mixing scalar source term.
  subroutine adjoint_mixing_scalar_source_term_allocate(obj)
    class(source_term_t), allocatable, intent(inout) :: obj
    allocate(adjoint_mixing_scalar_source_term_t::obj)
  end subroutine adjoint_mixing_scalar_source_term_allocate

  !> The common constructor using a JSON object.
  !! @param this The object.
  !! @param json The JSON object for the source.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  !! @param variable_name The name of the variable where the source term acts.
  subroutine adjoint_mixing_scalar_source_term_init_from_json(this, &
       json, fields, coef, variable_name)
    class(adjoint_mixing_scalar_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef
    character(len=*), intent(in) :: variable_name


  end subroutine adjoint_mixing_scalar_source_term_init_from_json

  !> The constructor from type components.
  !! @param this The source term.
  !! @param f_s RHS of the adjoint scalar.
  !! @param s the forward scalar field.
  !! @param obj_scale a scaling factor.
  !! @param phi_ref target concentration.
  !! @param mask the mask for the source term.
  !! @param if_mask whether to use the mask.
  !! @param coef The SEM coeffs.
  !! @param start_time start of the integration window.
  !! @param end_time end of the integration window.
  subroutine adjoint_mixing_scalar_source_term_init_from_components(this, &
       f_s, s, obj_scale, phi_ref, mask, if_mask, coef, start_time, end_time)
    class(adjoint_mixing_scalar_source_term_t), intent(inout) :: this
    type(field_t), pointer, intent(in) :: f_s
    type(field_t), intent(in), target :: s
    real(kind=rp), intent(in) :: obj_scale
    real(kind=rp), intent(in) :: phi_ref
    class(point_zone_t), intent(in), target :: mask
    logical, intent(in) :: if_mask
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(in) :: start_time
    real(kind=rp), intent(in) :: end_time

    type(field_list_t) :: fields

    call this%free()

    ! this is copying the fluid source term init
    ! We package the fields for the source term to operate on in a field list.
    call fields%init(1)
    call fields%assign(1, f_s)

    call this%init_base(fields, coef, start_time, end_time)
    call fields%free()

    ! point everything in the correct places
    this%s => s
    this%obj_scale = obj_scale
    this%phi_ref = phi_ref
    this%if_mask = if_mask

    if (this%if_mask) then
       this%mask => mask
       this%mask_volume = compute_masked_volume(this%mask, coef)
    else
       this%mask_volume = coef%volume
    end if

  end subroutine adjoint_mixing_scalar_source_term_init_from_components

  !> Destructor.
  subroutine adjoint_mixing_scalar_source_term_free(this)
    class(adjoint_mixing_scalar_source_term_t), intent(inout) :: this

    call this%free_base()
    nullify(this%s)
    nullify(this%mask)

  end subroutine adjoint_mixing_scalar_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param this The object.
  !! @param time The time state.
  subroutine adjoint_mixing_scalar_source_term_compute(this, time)
    class(adjoint_mixing_scalar_source_term_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: fs
    type(field_t), pointer :: work
    integer :: temp_indices(1)
    real(kind=rp) :: dbg_work, dbg_fs


    fs => this%fields%get(1)

    call neko_scratch_registry%request_field(work, temp_indices(1), .false.)
    ! \phi
    call field_copy(work, this%s)
    ! \phi - \phi_ref
    call field_cadd(work, -this%phi_ref)
    ! mask
    if (this%if_mask) then
       call mask_exterior_const(work, this%mask, 0.0_rp)
    end if

    ! append to RHS with scaling and mask volume
    call field_add2s2(fs, work, this%obj_scale / this%mask_volume)

    ! DEBUG (temporary, for investigate-passive-scalar-dealias diagnosis)
    ! NOTE: glsc2 is collective (MPI_Allreduce) -- it must be called on EVERY
    ! rank, outside the pe_rank guard, or rank 0 deadlocks waiting for the
    ! others to join the reduction. Only the write() is rank-guarded.
    if (mod(time%tstep, 200) .eq. 0) then
       dbg_work = sqrt(glsc2(work%x, work%x, work%size()))
       dbg_fs = sqrt(glsc2(fs%x, fs%x, fs%size()))
       if (pe_rank .eq. 0) then
          write(error_unit, '(A,I0,A,E15.6,A,E15.6,A,E15.6)') &
               'DEBUG mixing_scalar_src tstep=', time%tstep, &
               ' |work|=', dbg_work, &
               ' obj_scale/mask_volume=', this%obj_scale / this%mask_volume, &
               ' |fs_after|=', dbg_fs
          flush(error_unit)
       end if
    end if

    call neko_scratch_registry%relinquish_field(temp_indices)


  end subroutine adjoint_mixing_scalar_source_term_compute

end module adjoint_mixing_scalar_source_term
