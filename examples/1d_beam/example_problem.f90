! Copyright (c) 2024, The Neko Authors
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

module example_problem
  use num_types, only: rp

  use objective, only: objective_t
  use constraint, only: constraint_t

  use design, only: design_t
  use math, only: glsum
  use json_module, only: json_file
  use vector, only: vector_t

  use device, only: device_memcpy, DEVICE_TO_HOST
  use neko_config, only: NEKO_BCKND_DEVICE

  use mpi_f08, only: mpi_exscan, mpi_sum, MPI_INTEGER, MPI_Allreduce
  use comm, only: pe_rank, pe_size, neko_comm, mpi_real_precision
  use vector_math, only:  vector_glsum, vector_cmult, vector_cadd
  implicit none
  private

  ! ========================================================================== !
  ! Global beam parameters
  real(rp), public :: L_total   = 2.0_rp      ! total length (m)
  real(rp), public :: b         = 0.02_rp     ! width (m)
  real(rp), public :: rho       = 7800.0_rp   ! density (kg/m3)
  real(rp), public :: E         = 210.0e9_rp  ! Young's modulus (Pa)
  real(rp), public :: P         = 1000.0_rp   ! tip load (N)
  real(rp), public :: h_min     = 0.005_rp    ! min height (m)
  real(rp), public :: h_max     = 0.05_rp     ! max height (m)

  ! ========================================================================== !
  ! Objective: tip deflection
  type, public, extends(objective_t) :: mma_obj
   contains
     procedure, public, pass(this) :: init_json => mma_obj_init_from_json
     procedure, public, pass(this) :: init_from_components => &
          mma_obj_init_from_components
     procedure, public, pass(this) :: free => mma_obj_free
     procedure, public, pass(this) :: update_value => mma_obj_update_value
     procedure, public, pass(this) :: update_sensitivity => &
          mma_obj_update_sensitivity
  end type mma_obj
  ! ========================================================================== !
  ! Objective: beam weight
  type, public, extends(objective_t) :: beamweight_obj
   contains
     procedure, public, pass(this) :: init_json => beamweight_init_from_json
     procedure, public, pass(this) :: init_from_components => &
          beamweight_init_from_components
     procedure, public, pass(this) :: free => beamweight_free
     procedure, public, pass(this) :: update_value => beamweight_update_value
     procedure, public, pass(this) :: update_sensitivity => &
          beamweight_update_sensitivity
  end type beamweight_obj

  ! ========================================================================== !
  ! stress constraints
  type, public, extends(constraint_t) :: mma_con
     real(kind=rp) :: sign = 1.0_rp
   contains
     procedure, public, pass(this) :: init_json => mma_con_init_from_json
     procedure, public, pass(this) :: init_from_components => &
          mma_con_init_from_components
     procedure, public, pass(this) :: free => mma_con_free
     procedure, public, pass(this) :: update_value => mma_con_update_value
     procedure, public, pass(this) :: update_sensitivity => &
          mma_con_update_sensitivity
  end type mma_con

contains

  ! ========================================================================== !
  ! Methods for the Objective Function

  subroutine mma_obj_init_from_json(this, json, design)
    class(mma_obj), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(design_t), intent(in) :: design
    character(len=256), parameter :: name = 'mma_obj'
    real(kind=rp) :: weight

    call this%init_from_components(name, design, weight)

  end subroutine mma_obj_init_from_json

  subroutine mma_obj_init_from_components(this, name, design, weight)
    class(mma_obj), intent(inout) :: this
    character(len=*), intent(in) :: name
    class(design_t), intent(in) :: design
    real(kind=rp), intent(in) :: weight

    call this%init_base(name, design%size(), 1.0_rp)
  end subroutine mma_obj_init_from_components

  subroutine mma_obj_free(this)
    class(mma_obj), intent(inout) :: this
    call this%free_base()
  end subroutine mma_obj_free

  subroutine mma_obj_update_value(this, design)
    class(mma_obj), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(vector_t) :: design_values !h = (h_max-h_min) x + h_min
    integer :: offset, ierr, n, k
    real(rp) :: u_local, u_global, I_k, Delta_k, x_k, x_km1, Le

    ! Get local design values
    design_values = design%get_values()
    n = design%size()

    ! Project design variables to physical height
    call vector_cmult(design_values, (h_max - h_min) )
    call vector_cadd(design_values, h_min)
    

    ! Exclusive prefix sum for global indexing
    call MPI_Exscan(n, offset, 1, MPI_INTEGER, MPI_SUM, neko_comm, ierr)
    if (pe_rank == 0) offset = 0

    ! Global element length
    Le = L_total / real(design%size_global(), kind=rp)

    ! Local contribution to tip deflection
    u_local = 0.0_rp
    do k = 1, n
      x_km1 = Le * real(offset + k - 1, kind=rp)
      x_k   = Le * real(offset + k, kind=rp)
      Delta_k = ((L_total - x_km1)**3 - (L_total - x_k)**3)/3.0_rp
      I_k = b * design_values%x(k)**3 / 12.0_rp
      u_local = u_local + Delta_k / (E * I_k)
    end do

    ! Reduce to global tip deflection
    call MPI_Allreduce(u_local, u_global, 1, mpi_real_precision, MPI_SUM, &
         neko_comm, ierr)
    this%value = u_global

    if (pe_rank == 0) print *, "Global tip deflection =", this%value
  end subroutine mma_obj_update_value

  subroutine mma_obj_update_sensitivity(this, design)
    class(mma_obj), intent(inout) :: this
    class(design_t), intent(in) :: design

    this%sensitivity = 1.0_rp / real(design%size_global(), kind=rp)

  end subroutine mma_obj_update_sensitivity

  ! ========================================================================== !
  ! Methods for the Beam Weight Objective

  subroutine beamweight_init_from_json(this, json, design)
    class(beamweight_obj), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(design_t), intent(in) :: design
    character(len=256), parameter :: name = 'beamweight_obj'
    real(kind=rp) :: weight

    call this%init_from_components(name, design, weight)
  end subroutine beamweight_init_from_json

  subroutine beamweight_init_from_components(this, name, design, weight)
    class(beamweight_obj), intent(inout) :: this
    character(len=*), intent(in) :: name
    class(design_t), intent(in) :: design
    real(kind=rp), intent(in) :: weight

    call this%init_base(name, design%size(), 1.0_rp)
  end subroutine beamweight_init_from_components

  subroutine beamweight_free(this)
    class(beamweight_obj), intent(inout) :: this
    call this%free_base()
  end subroutine beamweight_free

  subroutine beamweight_update_value(this, design)
    class(beamweight_obj), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(vector_t) :: design_values     !h = (h_max-h_min) x + h_min
    integer :: ierr, n, k, offset
    real(rp) :: Le, local_mass, global_mass

    ! Get local design values
    design_values = design%get_values()
    n = design%size()

    ! Project design variables to physical height
    call vector_cmult(design_values, (h_max - h_min) )
    call vector_cadd(design_values, h_min)


    ! Exclusive prefix sum for global indexing
    call MPI_Exscan(n, offset, 1, MPI_INTEGER, MPI_SUM, neko_comm, ierr)
    if (pe_rank == 0) offset = 0

    ! Global element length
    Le = L_total / real(design%size_global(), kind=rp)

    global_mass =  rho * b * Le * vector_glsum(design_values, n)

    this%value = global_mass

    if (pe_rank == 0) print *, "Global beam weight =", this%value

  end subroutine beamweight_update_value

  subroutine beamweight_update_sensitivity(this, design)
    class(beamweight_obj), intent(inout) :: this
    class(design_t), intent(in) :: design
    integer :: n
    real(rp) :: Le

    n = design%size()
    Le = L_total / real(design%size_global(), kind=rp)

    ! Sensitivity: dm/dh_k = rho * b * Le
    this%sensitivity = rho * b * Le
  end subroutine beamweight_update_sensitivity

  ! ========================================================================== !
  ! Methods for the Constraint Function

  subroutine mma_con_init_from_json(this, json, design)
    class(mma_con), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(design_t), intent(in) :: design
    character(len=256), parameter :: name = 'mma_con'
    integer :: sign

    call this%init_from_components(name, design, sign)

  end subroutine mma_con_init_from_json

  subroutine mma_con_init_from_components(this, name, design, sign)
    class(mma_con), intent(inout) :: this
    character(len=*), intent(in) :: name
    class(design_t), intent(in) :: design
    integer, intent(in), optional :: sign
    integer :: sign_ = 1

    call this%init_base(name, design%size())

    if (present(sign)) sign_ = sign
    if (sign_ .lt. 0) this%sign = -1.0_rp
    if (sign_ .ge. 0) this%sign = 1.0_rp

  end subroutine mma_con_init_from_components

  subroutine mma_con_free(this)
    class(mma_con), intent(inout) :: this
    call this%free_base()
  end subroutine mma_con_free

  subroutine mma_con_update_value(this, design)
    class(mma_con), intent(inout) :: this
    class(design_t), intent(in) :: design


  end subroutine mma_con_update_value

  subroutine mma_con_update_sensitivity(this, design)
    class(mma_con), intent(inout) :: this
    class(design_t), intent(in) :: design

  end subroutine mma_con_update_sensitivity
end module example_problem