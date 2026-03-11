!> @file adjoint_lube_source_term.f90
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
  use num_types, only: rp
  use field_list, only: field_list_t
  use json_module, only: json_file
  use source_term, only: source_term_t
  use coefs, only: coef_t
  use interpolation, only: interpolator_t
  use space, only: space_t, GL
  use field, only: field_t
  use time_state, only: time_state_t
  use design, only: design_t
  use brinkman_design, only: brinkman_design_t
  use field_math, only: field_addcol3, field_copy, field_cmult, field_add2, &
       field_col3
  use mask_ops, only: mask_exterior_const
  use point_zone, only: point_zone_t
  use utils, only: neko_error
  use registry, only : neko_registry
  use scratch_registry, only: neko_scratch_registry, scratch_registry_t
  use neko_config, only: NEKO_BCKND_DEVICE
  use math, only: col2, invcol2
  use device_math, only: device_col2, device_invcol2
  implicit none
  private
  public :: adjoint_lube_source_term_allocate

  !> A adjoint source term corresponding to an objective of
  ! $K \int_\Omega \frac{1}{2}\chi|\mathbf{u}|^2$.
  type, public, extends(source_term_t) :: adjoint_lube_source_term_t

     !> u of the primal
     type(field_t), pointer :: u => null()
     !> v of the primal
     type(field_t), pointer :: v => null()
     !> w of the primal
     type(field_t), pointer :: w => null()
     !> \f$\chi\f$ the Brinkman amplitude
     type(field_t), pointer :: chi => null()
     !> a scale for this term
     real(kind=rp) :: K
     !> A mask for where the source term is evaluated
     class(point_zone_t), pointer :: mask => null()
     !> containing a mask?
     logical :: if_mask
     !> should dealiasing be used?
     logical :: dealias
     !> The original space used in the simulation
     type(space_t), pointer :: Xh_GLL
     !> The additional higher-order space used in dealiasing
     type(space_t), pointer :: Xh_GL
     !> cfs of the higher-order space
     type(coef_t), pointer :: c_Xh_GL
     !> Interpolator between the original and higher-order spaces
     type(interpolator_t), pointer :: GLL_to_GL
     !> GL scratch registry
     type(scratch_registry_t), pointer :: scratch_GL
     !> Volume of the objective domain
     real(kind=rp) :: volume
     !> Physical dimension
     integer :: gdim

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

  !> Allocator for the adjoint lube source term.
  subroutine adjoint_lube_source_term_allocate(obj)
    class(source_term_t), allocatable, intent(inout) :: obj
    allocate(adjoint_lube_source_term_t::obj)
  end subroutine adjoint_lube_source_term_allocate

  !> The common constructor using a JSON object.
  !! @param this The source term.
  !! @param json The JSON object for the source.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  !! @param variable_name The name of the variable where the source term acts.
  subroutine adjoint_lube_source_term_init_from_json(this, json, fields, coef, &
       variable_name)
    class(adjoint_lube_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef
    character(len=*), intent(in) :: variable_name
    ! real(kind=rp), allocatable :: values(:)
    ! real(kind=rp) :: start_time, end_time


    ! we shouldn't be initializing this from JSON
    ! maybe throw an error?


  end subroutine adjoint_lube_source_term_init_from_json

  !> The constructor from type components.
  !! @param this The source term.
  !! @param f_x, f_y, f_z the RHS of the adjoint
  !! @param design the design
  !! @param K a scale
  !! @param u, v, w the velocity fields of the primal
  !! @param mask the mask for the source term
  !! @param if_mask whether to use the mask
  !! @param coef The SEM coeffs.
  !! @param c_Xh_GL The SEM coeffs on the over integration mesh.
  !! @param GLL_to_GL Interpolator between GLL and GL.
  !! @param dealias weather this term should be overintegrated.
  !! @param volume volume of the objective domain.
  !! @param scratch_GL A scratch registry on the GL space.
  !! @param gdim physical dimension.
  subroutine adjoint_lube_source_term_init_from_components(this, &
       f_x, f_y, f_z, design, K, &
       u, v, w, &
       mask, if_mask, &
       coef, c_Xh_GL, GLL_to_GL, dealias, volume, scratch_GL, gdim)
    class(adjoint_lube_source_term_t), intent(inout) :: this
    type(field_t), pointer, intent(in) :: f_x, f_y, f_z
    class(design_t), intent(in), target :: design
    real(kind=rp), intent(in) :: K
    type(field_t), intent(in), target :: u, v, w
    class(point_zone_t), intent(in), target :: mask
    logical :: if_mask
    type(coef_t), intent(in) :: coef
    type(coef_t), intent(in), target :: c_Xh_GL
    type(interpolator_t), intent(in), target :: GLL_to_GL
    logical, intent(in) :: dealias
    real(kind=rp), intent(in) :: volume
    type(scratch_registry_t), intent(in), target :: scratch_GL
    integer, intent(in) :: gdim
    real(kind=rp) :: start_time
    real(kind=rp) :: end_time
    type(field_list_t) :: fields

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
    this%c_Xh_GL => c_Xh_GL
    this%Xh_GL => this%c_Xh_GL%Xh
    this%Xh_GLL => this%coef%Xh
    this%GLL_to_GL => GLL_to_GL
    this%dealias = dealias
    this%scratch_GL => scratch_GL
    this%u => u
    this%v => v
    this%w => w
    this%gdim = gdim

    this%volume = volume

    select type (design)
    type is (brinkman_design_t)
       this%chi => neko_registry%get_field("brinkman_amplitude")
    class default
       call neko_error('Unknown design type')
    end select

    this%K = K
    this%if_mask = if_mask
    if (this%if_mask) then
       this%mask => mask
    end if

  end subroutine adjoint_lube_source_term_init_from_components

  !> Destructor.
  subroutine adjoint_lube_source_term_free(this)
    class(adjoint_lube_source_term_t), intent(inout) :: this

    call this%free_base()
    nullify(this%u)
    nullify(this%v)
    nullify(this%w)
    nullify(this%c_Xh_GL)
    nullify(this%Xh_GL)
    nullify(this%Xh_GLL)
    nullify(this%GLL_to_GL)
    nullify(this%mask)
    nullify(this%scratch_GL)

  end subroutine adjoint_lube_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param this The source term.
  !! @param time The time state.
  subroutine adjoint_lube_source_term_compute(this, time)
    class(adjoint_lube_source_term_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: fu, fv, fw
    type(field_t), pointer :: work
    integer :: temp_indices(1)
    type(field_t), pointer :: accumulate, fld_GL, chi_GL
    integer :: temp_indices_GL(3)
    integer :: n_GL, nel

    fu => this%fields%get_by_index(1)
    fv => this%fields%get_by_index(2)
    fw => this%fields%get_by_index(3)

    ! BE SO CAREFUL HERE
    ! It make look the same as the Brinkman term, but it's assumed that
    ! this source term acts on the adjoint, and the u,v,w here come from
    ! the primal
    call neko_scratch_registry%request_field(work, temp_indices(1), .false.)
    call field_copy(work, this%chi)

    ! scale by K and volume
    call field_cmult(work, this%K / this%volume)

    ! mask
    if (this%if_mask) then
       call mask_exterior_const(work, this%mask, 0.0_rp)
    end if

    if (this%dealias) then
       nel = this%coef%msh%nelv
       n_GL = nel * this%Xh_GL%lxyz
       call this%scratch_GL%request_field(accumulate, temp_indices_GL(1), &
            .false.)
       call this%scratch_GL%request_field(fld_GL, temp_indices_GL(2), .false.)
       call this%scratch_GL%request_field(chi_GL, temp_indices_GL(3), .false.)

       call this%GLL_to_GL%map(chi_GL%x, work%x, nel, this%Xh_GL)

       ! u
       call this%GLL_to_GL%map(fld_GL%x, this%u%x, nel, this%Xh_GL)
       call field_col3(accumulate, chi_GL, fld_GL)
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_col2(accumulate%x_d, this%c_Xh_GL%B_d, n_GL)
          call this%GLL_to_GL%map(work%x, accumulate%x, nel, this%Xh_GLL)
          call device_invcol2(work%x_d, this%coef%B_d, work%size())
       else
          call col2(accumulate%x, this%c_Xh_GL%B, n_GL)
          call this%GLL_to_GL%map(work%x, accumulate%x, nel, this%Xh_GLL)
          call invcol2(work%x, this%coef%B, work%size())
       end if
       call field_add2(fu, work)

       ! v
       call this%GLL_to_GL%map(fld_GL%x, this%v%x, nel, this%Xh_GL)
       call field_col3(accumulate, chi_GL, fld_GL)
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_col2(accumulate%x_d, this%c_Xh_GL%B_d, n_GL)
          call this%GLL_to_GL%map(work%x, accumulate%x, nel, this%Xh_GLL)
          call device_invcol2(work%x_d, this%coef%B_d, work%size())
       else
          call col2(accumulate%x, this%c_Xh_GL%B, n_GL)
          call this%GLL_to_GL%map(work%x, accumulate%x, nel, this%Xh_GLL)
          call invcol2(work%x, this%coef%B, work%size())
       end if
       call field_add2(fv, work)

       ! w
       if (this%gdim .eq. 3) then
          call this%GLL_to_GL%map(fld_GL%x, this%w%x, nel, this%Xh_GL)
          call field_col3(accumulate, chi_GL, fld_GL)
          if (NEKO_BCKND_DEVICE .eq. 1) then
             call device_col2(accumulate%x_d, this%c_Xh_GL%B_d, n_GL)
             call this%GLL_to_GL%map(work%x, accumulate%x, nel, this%Xh_GLL)
             call device_invcol2(work%x_d, this%coef%B_d, work%size())
          else
             call col2(accumulate%x, this%c_Xh_GL%B, n_GL)
             call this%GLL_to_GL%map(work%x, accumulate%x, nel, this%Xh_GLL)
             call invcol2(work%x, this%coef%B, work%size())
          end if
          call field_add2(fw, work)
       end if

       call this%scratch_GL%relinquish_field(temp_indices_GL)

    else
       ! multiple and add the RHS
       call field_addcol3(fu, this%u, work)
       call field_addcol3(fv, this%v, work)
       if (this%gdim .eq. 3) then
          call field_addcol3(fw, this%w, work)
       end if

    end if

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine adjoint_lube_source_term_compute

end module adjoint_lube_source_term
