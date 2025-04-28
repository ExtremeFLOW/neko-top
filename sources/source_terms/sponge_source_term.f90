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
!> Implements the `sponge_source_term_t` type.
!
!
! This should be implemented into neko. Essentially just a sponge.
! There could be many more json params but I need to hard code and test this
! quickly.
module sponge_source_term
  use num_types, only: rp
  use field_list, only: field_list_t
  use json_module, only: json_file
  use source_term, only: source_term_t
  use coefs, only: coef_t
  use field, only: field_t
  use design, only: design_t
  use brinkman_design, only: brinkman_design_t
  use field_math, only: field_add2s2, field_sub3, field_rzero, field_col2, field_copy, field_sub2
  use scratch_registry, only: neko_scratch_registry
  use mask_ops, only: mask_exterior_const
  use point_zone, only: point_zone_t
  use utils, only: neko_error
  use field_registry, only: neko_field_registry
  use device, only : device_memcpy, HOST_TO_DEVICE
  use neko_config, only : NEKO_BCKND_DEVICE
  implicit none
  private

  !> A adjoint source term corresponding to an objective of
  ! \f$K \int_\Omega \frac{1}{2}\chi|\mathbf{u}|^2\f$.
  type, public, extends(source_term_t) :: sponge_source_term_t

     !> u,v,w corresponding to field on which the sponge acts
     type(field_t), pointer :: u,v,w
     !> sponge zone
     type(field_t), pointer :: sponge => null()
     !> sponge velocity
     type(field_t), pointer :: sponge_u => null()
     type(field_t), pointer :: sponge_v => null()
     type(field_t), pointer :: sponge_w => null()
     !> the strength
     real(kind=rp) :: A

   contains
     !> The common constructor using a JSON object.
     procedure, pass(this) :: init => sponge_source_term_init_from_json
     !> The constructor from type components.
     procedure, pass(this) :: init_from_components => &
          sponge_source_term_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => sponge_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute_ => sponge_source_term_compute
  end type sponge_source_term_t

contains
  !> The common constructor using a JSON object.
  !! @param this The source term.
  !! @param json The JSON object for the source.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  subroutine sponge_source_term_init_from_json(this, json, fields, coef)
    class(sponge_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef
    ! real(kind=rp), allocatable :: values(:)
    ! real(kind=rp) :: start_time, end_time

    


  end subroutine sponge_source_term_init_from_json

  !> The constructor from type components.
  subroutine sponge_source_term_init_from_components(this, &
       f_x, f_y, f_z, u, v, w, coef)
    class(sponge_source_term_t), intent(inout) :: this
    type(field_t), pointer, intent(in) :: f_x, f_y, f_z
    type(field_t), intent(in), target :: u, v, w
    type(coef_t), intent(in) :: coef
    real(kind=rp) :: start_time
    real(kind=rp) :: end_time
    type(field_list_t) :: fields
    integer :: i, n
    real :: x_end, x_start, x_drop, x_current, tmp_real, squashed, &
       x_width, x_fin_drop
    


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
    ! this is the perturbation!
    this%u => u
    this%v => v
    this%w => w

    !---------------------------------------------------------------------------
    ! Hard code this
    !---------------------------------------------------------------------------
    x_end = 35.0_rp

    x_width = 12.0_rp
    x_drop = 8.0_rp
    this%A = 2.0_rp

    x_start = x_end - x_width
    x_fin_drop = x_start + x_drop
    !---------------------------------------------------------------------------

    ! init
    allocate(this%sponge)
    allocate(this%sponge_u)
    allocate(this%sponge_v)
    allocate(this%sponge_w)
    call this%sponge%init(this%coef%dof, fld_name = "sponge")
    call this%sponge_u%init(this%coef%dof, fld_name = "sponge_u")
    call this%sponge_v%init(this%coef%dof, fld_name = "sponge_v")
    call this%sponge_w%init(this%coef%dof, fld_name = "sponge_w")

    ! zero
    call field_rzero(this%sponge_u)
    call field_rzero(this%sponge_v)
    call field_rzero(this%sponge_w)

    ! hard coded sponge
    n = this%sponge%size()
    do i = 1, n
        x_current = this%coef%dof%x(i,1,1,1)
        if (x_current .le. x_start) then
           tmp_real = 0.0_rp
        elseif (x_current .le. x_fin_drop .and. x_current .gt. x_start) then
           squashed = (x_current - x_start)/(x_drop)
           tmp_real = math_stepf(squashed)
        else
           tmp_real = 1.0_rp
        end if
        this%sponge%x(i,1,1,1) = tmp_real
    end do

    ! put on device if needed
    if (NEKO_BCKND_DEVICE .eq. 1) then
      call device_memcpy(this%sponge%x, this%sponge%x_d, n, host_to_device, .true.)
    end if


  end subroutine sponge_source_term_init_from_components

  !> Destructor.
  subroutine sponge_source_term_free(this)
    class(sponge_source_term_t), intent(inout) :: this

    call this%free_base()
  end subroutine sponge_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param this The source term.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine sponge_source_term_compute(this, t, tstep)
    class(sponge_source_term_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(field_t), pointer :: fu, fv, fw
    type(field_t), pointer :: work
    integer :: temp_indices(1)

    fu => this%fields%get_by_index(1)
    fv => this%fields%get_by_index(2)
    fw => this%fields%get_by_index(3)


    call neko_scratch_registry%request_field(work, temp_indices(1))

    ! u
    call field_copy(work, this%sponge_u)
    call field_sub2(work, this%u)
    ! call field_sub3(work, this%sponge_u, this%u)
    call field_col2(work, this%sponge)
    call field_add2s2 (fu, work, this%A)
    ! v
    call field_copy(work, this%sponge_v)
    call field_sub2(work, this%v)
    !call field_sub3(work, this%sponge_v, this%v)
    call field_col2(work, this%sponge)
    call field_add2s2 (fv, work, this%A)
     ! w
    !call field_sub3(work, this%sponge_w, this%w)
    call field_copy(work, this%sponge_w)
    call field_sub2(work, this%w)
    call field_col2(work, this%sponge)
    call field_add2s2 (fw, work, this%A)

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine sponge_source_term_compute

      real function math_stepf(x)
      implicit none

      ! argument list
      real x

      ! local variables
      real xdmin, xdmax
      parameter (xdmin = 0.0001, xdmax = 0.9999)
!-----------------------------------------------------------------------
      ! get function vale
      if (x.le.xdmin) then
         math_stepf = 0.0
      else if (x.le.xdmax) then
         math_stepf = 1./( 1. + exp(1./(x - 1.) + 1./x) )
      else
         math_stepf = 1.
      end if

      return
      end function math_stepf

end module sponge_source_term
