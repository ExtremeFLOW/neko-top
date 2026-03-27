!> @file design_sem.f90
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

!> Intermediate base class for SEM-backed designs.
module design_sem
  use design, only: design_t
  use coefs, only: coef_t
  use vector, only: vector_t
  use num_types, only: rp
  use math, only: copy, glsum, invcol1
  use device_math, only: device_copy, device_glsum, device_invcol1
  use device_math_ext, only: device_sqrt_inplace
  use neko_config, only: NEKO_BCKND_DEVICE
  use utils, only: neko_error
  use vector_math, only: vector_cmult, vector_rone, vector_copy
  implicit none
  private

  type, abstract, extends(design_t), public :: design_sem_t
     private
     !> MMA variable map (A in y = A x), built from SEM coefficients.
     type(vector_t) :: mma_map
     !> Gradient scaling in MMA space.
     type(vector_t) :: gradient_scale
     !> Whether to match norms of the gradient
     logical :: match_gradient_norm
     !> SEM MMA scaling option.
     integer :: sem_map_option = 2
   contains
     !> Initialize base design data and SEM-derived MMA map.
     procedure, pass(this) :: init_sem_base
     !> Free SEM base resources.
     procedure, pass(this) :: free_sem_base
     !> Set SEM MMA scaling option.
     procedure, pass(this) :: set_sem_map_option
     !> Provide MMA variable map.
     procedure, pass(this) :: get_mma_variable_map
     !> Build SEM map from coefficients.
     procedure, pass(this), private :: build_sem_map
     !> SEM designs provide objective sensitivity directly.
     procedure(design_sem_get_sensitivity), public, pass(this), deferred :: &
          get_sensitivity
  end type design_sem_t

  abstract interface
     subroutine design_sem_get_sensitivity(this, values)
       import design_sem_t, vector_t
       class(design_sem_t), intent(in) :: this
       type(vector_t), intent(inout) :: values
     end subroutine design_sem_get_sensitivity
  end interface

contains

  !> Initialize common SEM design state.
  subroutine init_sem_base(this, name, n, coef)
    class(design_sem_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    type(coef_t), intent(in) :: coef

    call this%init_base(name, n)
    call this%build_sem_map(coef)
  end subroutine init_sem_base

  !> Free SEM design base resources.
  subroutine free_sem_base(this)
    class(design_sem_t), intent(inout) :: this

    call this%mma_map%free()
    call this%gradient_scale%free()
    call this%free_base()
  end subroutine free_sem_base

  !> Set SEM map/scaling strategy.
  subroutine set_sem_map_option(this, option)
    class(design_sem_t), intent(inout) :: this
    integer, intent(in) :: option

    if (option .lt. 0 .or. option .gt. 7) then
      call neko_error('design_sem_t: sem_map_option must be in [0, 7]')
    end if

    this%sem_map_option = option
  end subroutine set_sem_map_option

  !> Return the SEM MMA map.
  subroutine get_mma_variable_map(this, map, gradient_scale, &
       match_gradient_norm)
    class(design_sem_t), intent(in) :: this
    type(vector_t), intent(inout) :: map
    type(vector_t), intent(inout) :: gradient_scale
    logical, intent(out) :: match_gradient_norm

    if (this%mma_map%size() .le. 0) then
      call neko_error('design_sem_t: MMA map is not initialized')
    end if
    if (this%gradient_scale%size() .le. 0) then
      call neko_error('design_sem_t: gradient scaling is not initialized')
    end if

    map = this%mma_map
    gradient_scale = this%gradient_scale
    match_gradient_norm = this%match_gradient_norm
  end subroutine get_mma_variable_map

  !> Build map A and gradient scaling.
  subroutine build_sem_map(this, coef)
    class(design_sem_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    integer :: n, i
    real(kind=rp) :: sum_B_half, global_count
    real(kind=rp) :: local_count(1)
    type(vector_t) :: sqrtB, B_inv
    integer :: option

    n = this%size()
    if (coef%dof%size() .ne. n) then
      call neko_error('coefficient/dof size mismatch')
    end if

    call this%mma_map%free()
    call this%mma_map%init(n)
    call this%gradient_scale%free()
    call this%gradient_scale%init(n)

    call sqrtB%init(n)
    ! note Binv in coef differs to this by the multiplicity
    call B_inv%init(n)


    ! Fill in some vectors
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(B_inv%x_d, coef%B_d, n)
       call device_invcol1(B_inv%x_d, n)
       call device_copy(sqrtB%x_d, coef%B_d, n)
       call device_sqrt_inplace(sqrtB%x_d, n)
    else
       call copy(B_inv%x, coef%B, n)
       call invcol1(B_inv%x, n)
       call copy(sqrtB%x, coef%B, n)
       do i = 1, n
          sqrtB%x(i) = sqrt(sqrtB%x(i))
       end do
    end if

    ! Some totals
    if (NEKO_BCKND_DEVICE .eq. 1) then
       sum_B_half = device_glsum(sqrtB%x_d, n)
    else
       sum_B_half = glsum(sqrtB%x, n)
    end if
    local_count(1) = real(n, rp)
    global_count = glsum(local_count, 1)

    ! Decision for scaling strategy.
    option = this%sem_map_option

    ! in the document A = mma_map, B = gradient_scale

    if (option .eq. 0) then
       ! do nothing (what we have originally)
       call vector_rone(this%gradient_scale)
       call vector_rone(this%mma_map)
       this%match_gradient_norm = .false.

    else if (option .eq. 1) then
       ! remove mass matrix, no change to MMA
       call vector_copy(this%gradient_scale, B_inv)
       call vector_rone(this%mma_map)
       this%match_gradient_norm = .false.

    else if (option .eq. 2) then
       ! remove mass matrix, rescale gradient, no change to MMA
       call vector_copy(this%gradient_scale, B_inv)
       call vector_rone(this%mma_map)
       this%match_gradient_norm = .true.

    else if (option .eq. 3) then
       ! remove mass matrix, map MMA (this is what's written in Martin's notes)
       call vector_copy(this%gradient_scale, B_inv)
       call vector_copy(this%mma_map, sqrtB)
       this%match_gradient_norm = .false.

    else if (option .eq. 4) then
       ! remove mass matrix, rescale gradient, map MMA
       call vector_copy(this%gradient_scale, B_inv)
       call vector_copy(this%mma_map, sqrtB)
       this%match_gradient_norm = .true.

    else if (option .eq. 5) then
       ! remove mass matrix, map MMA, constant scale of
       ! 1 / avg(sqrt(B)) (in MMA only)
       call vector_copy(this%gradient_scale, B_inv)
       call vector_copy(this%mma_map, sqrtB)
       call vector_cmult(this%mma_map, 1.0_rp / (sum_B_half / global_count))
       this%match_gradient_norm = .false.

    else if (option .eq. 6) then
       ! remove mass matrix, map MMA, constant scale of
       ! 1 / avg(sqrt(B)) (in MMA and squared in gradient)
       call vector_copy(this%gradient_scale, B_inv)
       call vector_copy(this%mma_map, sqrtB)
       call vector_cmult(this%mma_map, 1.0_rp / (sum_B_half / global_count))
       call vector_cmult(this%gradient_scale, (sum_B_half / global_count) ** 2)
       this%match_gradient_norm = .false.
    else if (option .eq. 7) then
       ! remove mass matrix, map MMA (this is what's written in Martin's notes)
       call vector_rone(this%gradient_scale)
       call vector_copy(this%mma_map, sqrtB)
       this%match_gradient_norm = .false.

    else
       call neko_error('design_sem_t: invalid map scaling option')
    end if

    call sqrtB%free()
    call B_inv%free()

  end subroutine build_sem_map

end module design_sem
