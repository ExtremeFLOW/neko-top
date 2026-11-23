!> @file adv_adjoint_dummy.f90
!! @copyright
!! Copyright (c) 2026, The Neko-TOP Authors
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
!> Implements `adv_adjoint_dummy_t`, a zero-valued adjoint advection.
module adv_adjoint_dummy
  use advection_adjoint, only: advection_adjoint_t
  use num_types, only: rp
  use space, only: space_t
  use field, only: field_t
  use coefs, only: coef_t
  implicit none
  private

  !> A zero-valued adjoint advection, used when the primal has no advection
  !! term (Stokes flow, `case.fluid.advection = false`). The linearised and
  !! adjoint convection terms are then identically zero, so every method is a
  !! no-op. Owns no memory, on the host or on the device.
  type, public, extends(advection_adjoint_t) :: adv_adjoint_dummy_t
   contains
     !> Constructor.
     procedure, pass(this) :: init => init_adv_adjoint_dummy
     !> Destructor.
     procedure, pass(this) :: free => free_adv_adjoint_dummy
     !> Add the linearised advection term for the fluid. Does nothing.
     procedure, pass(this) :: compute_linear => compute_linear_adv_dummy
     !> Add the adjoint advection term for the fluid. Does nothing.
     procedure, pass(this) :: compute_adjoint => compute_adjoint_adv_dummy
     !> Add the adjoint advection term for a scalar. Does nothing.
     procedure, pass(this) :: compute_adjoint_scalar => &
          compute_adjoint_scalar_adv_dummy
  end type adv_adjoint_dummy_t

contains

  !> Constructor.
  !! @param this The object.
  !! @param coef The coefficients of the (space, mesh) pair.
  subroutine init_adv_adjoint_dummy(this, coef)
    class(adv_adjoint_dummy_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef

  end subroutine init_adv_adjoint_dummy

  !> Destructor.
  !! @param this The object.
  subroutine free_adv_adjoint_dummy(this)
    class(adv_adjoint_dummy_t), intent(inout) :: this

  end subroutine free_adv_adjoint_dummy

  !> Add the linearised advection term for the fluid to the RHS.
  !! @details Here, does nothing.
  !! @param this The object.
  !! @param vx The x component of the perturbation velocity.
  !! @param vy The y component of the perturbation velocity.
  !! @param vz The z component of the perturbation velocity.
  !! @param vxb The x component of the baseflow.
  !! @param vyb The y component of the baseflow.
  !! @param vzb The z component of the baseflow.
  !! @param fx The x component of the source term.
  !! @param fy The y component of the source term.
  !! @param fz The z component of the source term.
  !! @param Xh The function space.
  !! @param coef The coefficients of the (Xh, mesh) pair.
  !! @param n Typically the size of the mesh.
  subroutine compute_linear_adv_dummy(this, vx, vy, vz, vxb, vyb, vzb, &
       fx, fy, fz, Xh, coef, n)
    class(adv_adjoint_dummy_t), intent(inout) :: this
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: vx, vy, vz
    type(field_t), intent(inout) :: vxb, vyb, vzb
    type(field_t), intent(inout) :: fx, fy, fz
    integer, intent(in) :: n

  end subroutine compute_linear_adv_dummy

  !> Add the adjoint advection term for the fluid to the RHS.
  !! @details Here, does nothing.
  !! @param this The object.
  !! @param vx The x component of the adjoint velocity.
  !! @param vy The y component of the adjoint velocity.
  !! @param vz The z component of the adjoint velocity.
  !! @param vxb The x component of the baseflow.
  !! @param vyb The y component of the baseflow.
  !! @param vzb The z component of the baseflow.
  !! @param fx The x component of the source term.
  !! @param fy The y component of the source term.
  !! @param fz The z component of the source term.
  !! @param Xh The function space.
  !! @param coef The coefficients of the (Xh, mesh) pair.
  !! @param n Typically the size of the mesh.
  subroutine compute_adjoint_adv_dummy(this, vx, vy, vz, vxb, vyb, vzb, &
       fx, fy, fz, Xh, coef, n)
    class(adv_adjoint_dummy_t), intent(inout) :: this
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: vx, vy, vz
    type(field_t), intent(inout) :: vxb, vyb, vzb
    type(field_t), intent(inout) :: fx, fy, fz
    integer, intent(in) :: n

  end subroutine compute_adjoint_adv_dummy

  !> Add the adjoint advection term for a scalar to the RHS.
  !! @details Here, does nothing.
  !! @param this The object.
  !! @param vxb The x component of the baseflow.
  !! @param vyb The y component of the baseflow.
  !! @param vzb The z component of the baseflow.
  !! @param s The adjoint scalar.
  !! @param fs The source term.
  !! @param Xh The function space.
  !! @param coef The coefficients of the (Xh, mesh) pair.
  !! @param n Typically the size of the mesh.
  !! @param dt Current time-step, not required for this method.
  subroutine compute_adjoint_scalar_adv_dummy(this, vxb, vyb, vzb, s, fs, &
       Xh, coef, n, dt)
    class(adv_adjoint_dummy_t), intent(inout) :: this
    type(field_t), intent(inout) :: vxb, vyb, vzb
    type(field_t), intent(inout) :: s
    type(field_t), intent(inout) :: fs
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    integer, intent(in) :: n
    real(kind=rp), intent(in), optional :: dt

  end subroutine compute_adjoint_scalar_adv_dummy

end module adv_adjoint_dummy
