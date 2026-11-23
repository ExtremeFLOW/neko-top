!> @file adjoint_viscous_dissipation_source_term.f90
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
!> Implements the `adjoint_viscous_dissipation_source_term_t` type.
!
!
! If the objective function is $\frac{\mu}{2} \int_M |\nabla u|^2$, the
! corresponding adjoint forcing is
! $ \mu \int_M \nabla v \cdot \nabla u $ in weak form, where $M$ is the
! masked region, or the whole domain if there is no mask.
module adjoint_viscous_dissipation_source_term
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
  use time_state, only: time_state_t
  use steady_simcomp, only: steady_simcomp_t
  use simcomp_executor, only: neko_simcomps
  use user_source_term, only: user_source_term_t
  use num_types, only: rp
  use field, only: field_t
  use registry, only: neko_registry
  use math, only: rzero, copy, chsign, cfill, invcol2, abscmp
  use device_math, only: device_copy, device_cmult, device_cfill, device_invcol2
  use neko_config, only: NEKO_BCKND_DEVICE
  use scratch_registry, only: neko_scratch_registry
  use mask_ops, only: mask_exterior_const
  use point_zone, only: point_zone_t
  use ax_product, only : ax_t, ax_helm_allocator
  use operators, only: grad, cdtp

  implicit none
  private
  public :: adjoint_viscous_dissipation_source_term_allocate

  !> An adjoint source term for objectives of viscous dissipation
  ! $\int_M \nabla v \cdot \nabla u $
  type, public, extends(source_term_t) :: &
       adjoint_viscous_dissipation_source_term_t
     !> u of the primal
     type(field_t), pointer :: u => null()
     !> v of the primal
     type(field_t), pointer :: v => null()
     !> w of the primal
     type(field_t), pointer :: w => null()
     !> the weight of the objective. Pointed to, so that continuation
     !! updates of the weight reach the forcing.
     real(kind=rp), pointer :: weight => null()
     !> a further, constant scale for the source term
     real(kind=rp) :: obj_scale
     !> Scale of the forcing, owned by the objective. Zero switches the
     !! forcing off.
     real(kind=rp), pointer :: scale => null()
     !> A mask for the region the objective is integrated over
     class(point_zone_t), pointer :: mask => null()
     !> containing a mask?
     logical :: if_mask
     !> an ax_helm type to compute weak laplacian
     class(ax_t), allocatable :: Ax
     !> volume of the objective domain.
     real(kind=rp) :: volume

   contains
     !> The common constructor using a JSON object.
     procedure, pass(this) :: init => &
          adjoint_viscous_dissipation_source_term_init_from_json
     !> The constructor from type components.
     procedure, pass(this) :: init_from_components => &
          adjoint_viscous_dissipation_source_term_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => &
          adjoint_viscous_dissipation_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute_ => &
          adjoint_viscous_dissipation_source_term_compute
  end type adjoint_viscous_dissipation_source_term_t

contains

  !> Allocator for the adjoint viscous dissipation source term.
  subroutine adjoint_viscous_dissipation_source_term_allocate(obj)
    class(source_term_t), allocatable, intent(inout) :: obj
    allocate(adjoint_viscous_dissipation_source_term_t::obj)
  end subroutine adjoint_viscous_dissipation_source_term_allocate

  !> The common constructor using a JSON object.
  !! @param json The JSON object for the source.
  !! @param this The source term.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  !! @param variable_name The name of the variable where the source term acts.
  subroutine adjoint_viscous_dissipation_source_term_init_from_json(this, &
       json, fields, coef, variable_name)
    class(adjoint_viscous_dissipation_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef
    character(len=*), intent(in) :: variable_name

    ! we shouldn't be initializing this from JSON
    ! maybe throw an error?


  end subroutine adjoint_viscous_dissipation_source_term_init_from_json

  !> The constructor from type components.
  !! @param this The source term.
  !! @param f_x, f_y, f_z the RHS of the adjoint equations
  !! @param u, v, w the flow fields of the primal
  !! @param weight the weight of the objective. Must stay valid while the
  !! term is computed.
  !! @param obj_scale a further, constant scaling factor
  !! @param mask the mask for the source term
  !! @param if_mask whether to use the mask
  !! @param coef The SEM coeffs.
  !! @param volume volume of the objective domain.
  !! @param start_time when to start applying the source term.
  !! @param end_time when to stop applying the source term.
  !! @param scale scale of the forcing, owned by the objective.
  subroutine adjoint_viscous_dissipation_source_term_init_from_components( &
       this, f_x, f_y, f_z, u, v, w, weight, obj_scale, mask, if_mask, &
       coef, volume, start_time, end_time, scale)
    class(adjoint_viscous_dissipation_source_term_t), intent(inout) :: this
    type(field_t), pointer, intent(in) :: f_x, f_y, f_z
    type(field_t), intent(in), target :: u, v, w
    real(kind=rp), intent(in), target :: weight
    real(kind=rp), intent(in) :: obj_scale
    class(point_zone_t), intent(in), pointer :: mask
    logical, intent(in) :: if_mask
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(in) :: volume
    real(kind=rp), intent(in) :: start_time
    real(kind=rp), intent(in) :: end_time
    real(kind=rp), pointer, intent(in) :: scale

    type(field_list_t) :: fields

    call this%free()

    ! this is copying the fluid source term init
    ! We package the fields for the source term to operate on in a field list.
    call fields%init(3)
    call fields%assign(1, f_x)
    call fields%assign(2, f_y)
    call fields%assign(3, f_z)

    call this%init_base(fields, coef, start_time, end_time)
    call fields%free()

    ! point everything in the correct places
    this%u => u
    this%v => v
    this%w => w

    this%weight => weight
    this%obj_scale = obj_scale
    this%scale => scale
    this%volume = volume

    this%if_mask = if_mask
    if (this%if_mask) then
       this%mask => mask
    end if

    ! Initialize the ax_helm object
    call ax_helm_allocator(this%Ax, type_name = "standard")

  end subroutine adjoint_viscous_dissipation_source_term_init_from_components

  !> Destructor.
  subroutine adjoint_viscous_dissipation_source_term_free(this)
    class(adjoint_viscous_dissipation_source_term_t), intent(inout) :: this

    call this%free_base()
    nullify(this%u)
    nullify(this%v)
    nullify(this%w)
    nullify(this%mask)
    nullify(this%weight)
    nullify(this%scale)
    if (allocated(this%Ax)) then
       deallocate(this%Ax)
    end if

  end subroutine adjoint_viscous_dissipation_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param this The source term.
  !! @param time The time state.
  subroutine adjoint_viscous_dissipation_source_term_compute(this, time)
    class(adjoint_viscous_dissipation_source_term_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: u, v, w
    type(field_t), pointer :: fu, fv, fw
    type(field_t), pointer :: work
    integer :: temp_indices(1)
    integer n

    ! Off outside the objective's time window
    if (abscmp(this%scale, 0.0_rp)) return

    fu => this%fields%get_by_index(1)
    fv => this%fields%get_by_index(2)
    fw => this%fields%get_by_index(3)

    n = fu%size()

    u => this%u
    v => this%v
    w => this%w

    call neko_scratch_registry%request_field(work, temp_indices(1), .false.)

    associate(coef => this%coef)

      ! Note that axhelm computes h1 * lap u + h2 * u
      ! we set h1 = 1 and h2 = 0 to compute the weak laplacian.
      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_cfill(coef%h1_d, 1.0_rp, n)
         call device_cfill(coef%h2_d, 0.0_rp, n)
      else
         call cfill(coef%h1, 1.0_rp, n)
         call cfill(coef%h2, 0.0_rp, n)
      end if
      coef%ifh2 = .false.

      ! With a mask, the weak Laplacian must be integrated over the masked
      ! region only. Masking the full weak Laplacian instead would drop the
      ! flux through the mask boundary on the elements that boundary cuts.

      ! ------------------------------------------------------------------------
      ! u

      if (this%if_mask) then
         call masked_weak_laplacian(work, u, this%mask, coef)
      else
         call this%Ax%compute(work%x, u%x, coef, coef%msh, coef%xh)
      end if

      ! pre-divide out the mass matrix to counteract it's multiplication
      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_invcol2(work%x_d, coef%B_d, work%size())
      else
         call invcol2(work%x, coef%B, work%size())
      end if

      ! add to RHS
      call field_add2s2(fu, work, &
           this%scale * (this%weight * this%obj_scale) / this%volume)

      ! ------------------------------------------------------------------------
      ! v

      if (this%if_mask) then
         call masked_weak_laplacian(work, v, this%mask, coef)
      else
         call this%Ax%compute(work%x, v%x, coef, coef%msh, coef%xh)
      end if

      ! pre-divide out the mass matrix to counteract it's multiplication
      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_invcol2(work%x_d, coef%B_d, work%size())
      else
         call invcol2(work%x, coef%B, work%size())
      end if

      ! add to RHS
      call field_add2s2(fv, work, &
           this%scale * (this%weight * this%obj_scale) / this%volume)

      ! ------------------------------------------------------------------------
      ! w

      if (this%if_mask) then
         call masked_weak_laplacian(work, w, this%mask, coef)
      else
         call this%Ax%compute(work%x, w%x, coef, coef%msh, coef%xh)
      end if

      ! pre-divide out the mass matrix to counteract it's multiplication
      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_invcol2(work%x_d, coef%B_d, work%size())
      else
         call invcol2(work%x, coef%B, work%size())
      end if

      ! add to RHS
      call field_add2s2(fw, work, &
           this%scale * (this%weight * this%obj_scale) / this%volume)

      call neko_scratch_registry%relinquish_field(temp_indices)

    end associate

  end subroutine adjoint_viscous_dissipation_source_term_compute

  !> Computes the weak Laplacian integrated over a masked region,
  !! \f$ \sum_i D_i^T I_M B D_i c \f$, element by element and unassembled.
  !! @details \f$ D_i \f$ is the derivative in direction \f$ i \f$, \f$ B \f$
  !! the mass matrix and \f$ I_M \f$ zeroes the nodes outside the mask. This
  !! is the gradient of \f$ \frac{1}{2} \sum_{k \in M} B_k |\nabla c|_k^2 \f$,
  !! the masked integral the objective evaluates. On the elements the mask
  !! boundary cuts, it differs from the masked full weak Laplacian
  !! \f$ I_M \sum_i D_i^T B D_i c \f$ by the flux through that boundary.
  !! @param res The result.
  !! @param c The field the operator is applied to.
  !! @param mask The masked region.
  !! @param coef The SEM coefficients.
  subroutine masked_weak_laplacian(res, c, mask, coef)
    type(field_t), intent(inout) :: res
    type(field_t), intent(in) :: c
    class(point_zone_t), intent(inout) :: mask
    type(coef_t), intent(in) :: coef
    type(field_t), pointer :: dcdx, dcdy, dcdz, work
    integer :: temp_indices(4)

    call neko_scratch_registry%request_field(dcdx, temp_indices(1), .false.)
    call neko_scratch_registry%request_field(dcdy, temp_indices(2), .false.)
    call neko_scratch_registry%request_field(dcdz, temp_indices(3), .false.)
    call neko_scratch_registry%request_field(work, temp_indices(4), .false.)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call grad(dcdx%x_d, dcdy%x_d, dcdz%x_d, c%x_d, coef)
    else
       call grad(dcdx%x, dcdy%x, dcdz%x, c%x, coef)
    end if

    call mask_exterior_const(dcdx, mask, 0.0_rp)
    call mask_exterior_const(dcdy, mask, 0.0_rp)
    call mask_exterior_const(dcdz, mask, 0.0_rp)

    call cdtp(res%x, dcdx%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call cdtp(work%x, dcdy%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call field_add2(res, work)
    call cdtp(work%x, dcdz%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
    call field_add2(res, work)

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine masked_weak_laplacian

end module adjoint_viscous_dissipation_source_term
