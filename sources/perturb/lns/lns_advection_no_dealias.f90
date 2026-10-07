!> @file lns_advection_no_dealias.f90
!! @copyright
!! Copyright (c) 2024-2025, The Neko-TOP Authors
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
!> Subroutines to add perturbed advection terms to the RHS of a transport
!! equation.
module lns_advection_no_dealias
  use perturb_advection, only: perturb_advection_t
  use num_types, only: rp
  use math, only: subcol3, rzero
  use space, only: space_t
  use field, only: field_t
  use coefs, only: coef_t
  use scratch_registry, only: neko_scratch_registry
  use neko_config, only: NEKO_BCKND_DEVICE, NEKO_BCKND_SX, NEKO_BCKND_XSMM, &
       NEKO_BCKND_OPENCL, NEKO_BCKND_CUDA, NEKO_BCKND_HIP
  use operators, only: opgrad, conv1, cdtp
  use utils, only: neko_error
  use device_math, only: device_vdot3, device_sub2, device_add4, &
       device_col3, device_subcol3, device_rzero
  use device, only: device_free, device_map, device_get_ptr
  use, intrinsic :: iso_c_binding, only: c_ptr, C_NULL_PTR, &
       c_associated
  implicit none
  private

  !> Type encapsulating advection routines with no dealiasing applied
  type, public, extends(perturb_advection_t) :: lns_advection_no_dealias_t
     real(kind=rp), allocatable :: temp(:)
     type(c_ptr) :: temp_d = C_NULL_PTR
     logical :: if_adjoint
   contains
     !> Add the linearized advection term for the fluid, i.e.
     !! \f$u' \cdot \nabla \bar{U} + \bar{U} \cdot \nabla u' \f$, to
     !! the RHS.
     procedure, pass(this) :: compute_linear => linear_advection_no_dealias
     !> Add the adjoint advection term for the fluid in weak form, i.e.
     !! \f$ \int_\Omega v \cdot u' (\nabla \bar{U})^T u^\dagger d\Omega
     !! + \int_\Omega \nabla v \cdot (\bar{U} \otimes u^\dagger) d \Omega \f$
     !! , to
     !! the RHS.
     procedure, pass(this) :: compute_adjoint => lns_advection_no_dealias_no_adjoint
     
     !> Compute the adjoint passive scalar.
     ! If one integrates by parts, this essentially switches sign and adds some
     ! boundary terms.
     ! We keep the differential operator on the test function
     procedure, pass(this) :: compute_adjoint_scalar => &
          lns_advection_no_dealias_no_adjoint_scalar
     ! NOTE
     ! This linearized advection term is the same as a normal advection term
     ! so not sure what to do here...

     !> Constructor
     procedure, pass(this) :: init => init_no_dealias
     !> Destructor
     procedure, pass(this) :: free => free_no_dealias
  end type lns_advection_no_dealias_t

contains

  !> Constructor
  !! @param this The object.
  !! @param coef The coefficients of the (space, mesh) pair.
  subroutine init_no_dealias(this, coef, if_adjoint)
    class(lns_advection_no_dealias_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    logical, intent(in) :: if_adjoint

    this%if_adjoint = if_adjoint

    allocate(this%temp(coef%dof%size()))

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%temp, this%temp_d, coef%dof%size())
    end if

  end subroutine init_no_dealias

  !> Destructor
  subroutine free_no_dealias(this)
    class(lns_advection_no_dealias_t), intent(inout) :: this

    if (allocated(this%temp)) then
       deallocate(this%temp)
    end if
    if (c_associated(this%temp_d)) then
       call device_free(this%temp_d)
    end if
  end subroutine free_no_dealias
  subroutine linear_advection_no_dealias(this, vx, vy, vz, vxb, vyb, vzb, fx, &
       fy, fz, Xh, coef, n)
    implicit none
    class(lns_advection_no_dealias_t), intent(inout) :: this
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: vx, vy, vz
    type(field_t), intent(inout) :: vxb, vyb, vzb
    integer, intent(in) :: n
    type(field_t), intent(inout) :: fx, fy, fz
    type(c_ptr) :: fx_d, fy_d, fz_d

    if (NEKO_BCKND_DEVICE .eq. 1) then
       fx_d = fx%x_d
       fy_d = fy%x_d
       fz_d = fz%x_d

       ! (x)
       call conv1(this%temp, vx%x, vxb%x, vyb%x, vzb%x, Xh, coef)
       call device_subcol3 (fx_d, coef%B_d, this%temp_d, n)
       call conv1(this%temp, vxb%x, vx%x, vy%x, vz%x, Xh, coef)
       call device_subcol3 (fx_d, coef%B_d, this%temp_d, n)

       ! (y)
       call conv1(this%temp, vy%x, vxb%x, vyb%x, vzb%x, Xh, coef)
       call device_subcol3 (fy_d, coef%B_d, this%temp_d, n)
       call conv1(this%temp, vyb%x, vx%x, vy%x, vz%x, Xh, coef)
       call device_subcol3 (fy_d, coef%B_d, this%temp_d, n)

       ! (z)
       if (coef%Xh%lz .eq. 1) then
          call device_rzero (this%temp_d, n)
       else
          call conv1(this%temp, vz%x, vxb%x, vyb%x, vzb%x, Xh, coef)
          call device_subcol3(fz_d, coef%B_d, this%temp_d, n)
          call conv1(this%temp, vzb%x, vx%x, vy%x, vz%x, Xh, coef)
          call device_subcol3(fz_d, coef%B_d, this%temp_d, n)
       end if
    else
       ! (x)
       call conv1(this%temp, vx%x, vxb%x, vyb%x, vzb%x, Xh, coef)
       call subcol3 (fx%x, coef%B, this%temp, n)
       call conv1(this%temp, vxb%x, vx%x, vy%x, vz%x, Xh, coef)
       call subcol3 (fx%x, coef%B, this%temp, n)

       ! (y)
       call conv1(this%temp, vy%x, vxb%x, vyb%x, vzb%x, Xh, coef)
       call subcol3 (fy%x, coef%B, this%temp, n)
       call conv1(this%temp, vyb%x, vx%x, vy%x, vz%x, Xh, coef)
       call subcol3 (fy%x, coef%B, this%temp, n)

       ! (z)
       if (coef%Xh%lz .eq. 1) then
          call rzero (this%temp, n)
       else
          call conv1(this%temp, vz%x, vxb%x, vyb%x, vzb%x, Xh, coef)
          call subcol3(fz%x, coef%B, this%temp, n)
          call conv1(this%temp, vzb%x, vx%x, vy%x, vz%x, Xh, coef)
          call subcol3(fz%x, coef%B, this%temp, n)
       end if
    end if

  end subroutine linear_advection_no_dealias

  subroutine lns_advection_no_dealias_no_adjoint(this, vx, vy, vz, vxb, vyb, vzb, fx, &
       fy, fz, Xh, coef, n)
    class(lns_advection_no_dealias_t), intent(inout) :: this
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: vx, vy, vz, vxb, vyb, vzb, fx, fy, fz
    integer, intent(in) :: n
    call neko_error("Adjoint advection called on an LNS advection operator.")
  end subroutine lns_advection_no_dealias_no_adjoint

  subroutine lns_advection_no_dealias_no_adjoint_scalar(this, vxb, vyb, vzb, s, &
       fs, Xh, coef, n, dt)
    class(lns_advection_no_dealias_t), intent(inout) :: this
    type(field_t), intent(inout) :: vxb, vyb, vzb, s, fs
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    integer, intent(in) :: n
    real(kind=rp), intent(in), optional :: dt
    call neko_error("Adjoint scalar advection called on an LNS advection operator.")
  end subroutine lns_advection_no_dealias_no_adjoint_scalar

end module lns_advection_no_dealias
