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
!> Implements the `adjoint_minimum_dissipation_source_term_t` type.
!
!
! If the objective function $\int |\nabla u|^2$,
! the corresponding adjoint forcing is $ \int \nabla v \cdot \nabla u $ in weak
! form.
module adjoint_minimum_dissipation_source_term
  use num_types, only: rp
  use field_list, only: field_list_t
  use json_module, only: json_file
  use json_utils, only: json_get, json_get_or_default
  use source_term, only: source_term_t
  use coefs, only: coef_t
  use neko_config, only: NEKO_BCKND_DEVICE
  use utils, only: neko_error
  use field, only: field_t
  use field_math, only: field_subcol3, field_add2, field_add2s2, field_rzero
  use json_module, only: json_file
  use steady_simcomp, only: steady_simcomp_t
  use simcomp_executor, only: neko_simcomps
  use fluid_user_source_term, only: fluid_user_source_term_t
  use num_types, only: rp
  use field, only: field_t
  use field_registry, only: neko_field_registry
  use math, only: rzero, copy, chsign, invcol2
  use device_math, only: device_copy, device_cmult, device_invcol2
  use neko_config, only: NEKO_BCKND_DEVICE
  use operators, only: opgrad, cdtp
  use scratch_registry, only: neko_scratch_registry
  use mask_ops, only: mask_exterior_const
  use point_zone, only: point_zone_t
  implicit none
  private

  !> An adjoint source term for objectives of minimum dissipation
  ! $\int \nabla v \cdot \nabla u $
  type, public, extends(source_term_t) :: &
       adjoint_minimum_dissipation_source_term_t
     !> u of the primal
     type(field_t), pointer :: u => null()
     !> v of the primal
     type(field_t), pointer :: v => null()
     !> w of the primal
     type(field_t), pointer :: w => null()
     !> a scale for the source term
     real(kind=rp) :: obj_scale
     !> A mask for where the source term is evaluated
     class(point_zone_t), pointer :: mask => null()
     !> containing a mask?
     logical :: if_mask

   contains
     !> The common constructor using a JSON object.
     procedure, pass(this) :: init => &
          adjoint_minimum_dissipation_source_term_init_from_json
     !> The constructor from type components.
     procedure, pass(this) :: init_from_components => &
          adjoint_minimum_dissipation_source_term_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => &
          adjoint_minimum_dissipation_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute_ => &
          adjoint_minimum_dissipation_source_term_compute
  end type adjoint_minimum_dissipation_source_term_t

contains
  !> The common constructor using a JSON object.
  !! @param json The JSON object for the source.
  !! @param this The source term.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  !! @param variable_name The name of the variable where the source term acts.
  subroutine adjoint_minimum_dissipation_source_term_init_from_json(this, &
       json, fields, coef, variable_name)
    class(adjoint_minimum_dissipation_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef
    character(len=*), intent(in) :: variable_name

    ! we shouldn't be initializing this from JSON
    ! maybe throw an error?


  end subroutine adjoint_minimum_dissipation_source_term_init_from_json

  !> The constructor from type components.
  !! @param this The source term.
  !! @param f_x, f_y, f_z the RHS of the adjoint equations
  !! @param u, v, w the flow fields of the primal
  !! @param obj_scale a scaling factor
  !! @param mask the mask for the source term
  !! @param if_mask whether to use the mask
  !! @param coef The SEM coeffs.
  subroutine adjoint_minimum_dissipation_source_term_init_from_components(this,&
       f_x, f_y, f_z, &
       u, v, w, obj_scale, &
       mask, if_mask, &
       coef)
    class(adjoint_minimum_dissipation_source_term_t), intent(inout) :: this
    type(field_t), pointer, intent(in) :: f_x, f_y, f_z
    type(field_list_t) :: fields
    type(coef_t) :: coef
    real(kind=rp) :: start_time
    real(kind=rp) :: end_time
    real(kind=rp) :: obj_scale
    type(field_t), intent(in), target :: u, v, w
    class(point_zone_t), intent(in), target :: mask
    logical :: if_mask

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
    this%u => u
    this%v => v
    this%w => w

    this%obj_scale = obj_scale

    this%if_mask = if_mask
    if (this%if_mask) then
       this%mask => mask
    end if

  end subroutine adjoint_minimum_dissipation_source_term_init_from_components

  !> Destructor.
  subroutine adjoint_minimum_dissipation_source_term_free(this)
    class(adjoint_minimum_dissipation_source_term_t), intent(inout) :: this

    call this%free_base()
  end subroutine adjoint_minimum_dissipation_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param this The source term.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine adjoint_minimum_dissipation_source_term_compute(this, t, tstep)
    class(adjoint_minimum_dissipation_source_term_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(field_t), pointer :: u, v, w
    type(field_t), pointer :: fu, fv, fw
    type(field_t), pointer :: dudx, dudy, dudz
    type(field_t), pointer :: dvdx, dvdy, dvdz
    type(field_t), pointer :: dwdx, dwdy, dwdz
    type(field_t), pointer :: work, result
    type(field_t), pointer :: t1 , t2
    integer :: temp_indices(11)
    integer n


    fu => this%fields%get_by_index(1)
    fv => this%fields%get_by_index(2)
    fw => this%fields%get_by_index(3)

    n = fu%size()

    ! fuck I'm not sure about this... I need a pen and paper
    ! also there should be a way to pre-process this forcing term...
    ! instead of recalculating it every time
    u => this%u
    v => this%v
    w => this%w
    call neko_scratch_registry%request_field(dudx, temp_indices(1))
    call neko_scratch_registry%request_field(dudy, temp_indices(2))
    call neko_scratch_registry%request_field(dudz, temp_indices(3))
    call neko_scratch_registry%request_field(dvdx, temp_indices(4))
    call neko_scratch_registry%request_field(dvdy, temp_indices(5))
    call neko_scratch_registry%request_field(dvdz, temp_indices(6))
    call neko_scratch_registry%request_field(dwdx , temp_indices(7))
    call neko_scratch_registry%request_field(dwdy , temp_indices(8))
    call neko_scratch_registry%request_field(dwdz , temp_indices(9))
    call neko_scratch_registry%request_field(work , temp_indices(10))
    call neko_scratch_registry%request_field(result , temp_indices(11))


    ! The forcing term, in weak form is \nabla u . \nabla v 
    ! Note that in neko all forces are assumed to be in strong form and are
    ! hence multiplied with the mass matrix before adding to the RHS, which
    ! we need to undo in preparation.
    !
    ! Discretely this is D^T G D u (see Deville pg)

    associate (coef => this%coef)

    ! G D u
    call opgrad(dudx%x, dudy%x, dudz%x, u%x, coef)
    call opgrad(dvdx%x, dvdy%x, dvdz%x, v%x, coef)
    call opgrad(dwdx%x, dwdy%x, dwdz%x, w%x, coef)

    if (NEKO_BCKND_DEVICE .eq. 1) then
        call device_col2(dudx%x_d, coef%g_11_d, result%size())
        call device_col2(dudx%x_d, coef%g_11_d, result%size())
    else
        call invcol2(result%x, coef%B, result%size())
    end if
    !--------------------------
    ! u component of D^T G D u
    call field_rzero(result)
    call cdtp(work%x, dudx%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call field_add2(result, work)
    call cdtp(work%x, dvdx%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call field_add2(result, work)
    call cdtp(work%x, dwdx%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
    call field_add2(result, work)

    ! pre-divide out the mass matrix to counteract it's multiplication
    if (NEKO_BCKND_DEVICE .eq. 1) then
        call device_invcol2(result%x_d, coef%B_d, result%size())
    else
        call invcol2(result%x, coef%B, result%size())
    end if

    ! mask
    if (this%if_mask) then
       call mask_exterior_const(result, this%mask, 0.0_rp)
    end if

    ! add to RHS
    call field_add2s2(fu, result, this%obj_scale)

    !--------------------------
    ! v component of D^T G D u
    call field_rzero(result)
    call cdtp(work%x, dudy%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call field_add2(result, work)
    call cdtp(work%x, dvdy%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call field_add2(result, work)
    call cdtp(work%x, dwdy%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
    call field_add2(result, work)

    ! pre-divide out the mass matrix to counteract it's multiplication
    if (NEKO_BCKND_DEVICE .eq. 1) then
        call device_invcol2(result%x_d, coef%B_d, result%size())
    else
        call invcol2(result%x, coef%B, result%size())
    end if

    ! mask
    if (this%if_mask) then
       call mask_exterior_const(result, this%mask, 0.0_rp)
    end if

    ! add to RHS
    call field_add2s2(fv, result, this%obj_scale)

    !--------------------------
    ! w component of D^T G D u
    call field_rzero(result)
    call cdtp(work%x, dudz%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call field_add2(result, work)
    call cdtp(work%x, dvdz%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call field_add2(result, work)
    call cdtp(work%x, dwdz%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
    call field_add2(result, work)

    ! pre-divide out the mass matrix to counteract it's multiplication
    if (NEKO_BCKND_DEVICE .eq. 1) then
        call device_invcol2(result%x_d, coef%B_d, result%size())
    else
        call invcol2(result%x, coef%B, result%size())
    end if

    ! mask
    if (this%if_mask) then
       call mask_exterior_const(result, this%mask, 0.0_rp)
    end if

    ! add to RHS
    call field_add2s2(fw, result, this%obj_scale)

    end associate

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine adjoint_minimum_dissipation_source_term_compute

end module adjoint_minimum_dissipation_source_term
