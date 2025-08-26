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

  use device, only: device_memcpy, HOST_TO_DEVICE, DEVICE_TO_HOST
  use neko_config, only: NEKO_BCKND_DEVICE

  use mpi_f08, only: mpi_exscan, mpi_sum, MPI_INTEGER, MPI_Allreduce
  use comm, only: pe_rank, pe_size, neko_comm, mpi_real_precision
  use vector_math, only:  vector_glsum, vector_cmult, vector_cadd, &
       vector_col2, vector_invcol2, vector_copy

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
  real(rp), public :: u_tip_max = 0.25_rp     ! max tip deflection (m)
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
  ! Stress constraint for specific global elements
  type, public, extends(constraint_t) :: stress_con
    integer :: global_element_index ! Global element index this constraint 
    real(rp) :: sigma_max           ! Max stress for this element
    logical :: is_local             ! If this element is on current MPI-rank
    integer :: local_index          ! Local index if is_local
   contains
     procedure, public, pass(this) :: init_stress_con
     procedure, public, pass(this) :: free => stress_con_free
     procedure, public, pass(this) :: update_value => stress_con_update_value
     procedure, public, pass(this) :: update_sensitivity => &
          stress_con_update_sensitivity
  end type stress_con

  ! ========================================================================== !
  ! stress constraints
  ! type, public, extends(constraint_t) :: mma_con
  !    real(kind=rp) :: sign = 1.0_rp
  !  contains
  !    procedure, public, pass(this) :: init_json => mma_con_init_from_json
  !    procedure, public, pass(this) :: init_from_components => &
  !         mma_con_init_from_components
  !    procedure, public, pass(this) :: free => mma_con_free
  !    procedure, public, pass(this) :: update_value => mma_con_update_value
  !    procedure, public, pass(this) :: update_sensitivity => &
  !         mma_con_update_sensitivity
  ! end type mma_con

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

    call this%init_base(name, design%size(), weight)
  end subroutine mma_obj_init_from_components

  subroutine mma_obj_free(this)
    class(mma_obj), intent(inout) :: this
    call this%free_base()
  end subroutine mma_obj_free

  subroutine mma_obj_update_value(this, design)
    class(mma_obj), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(vector_t) :: h, I, contrib
    integer :: ierr, n, offset, k
    real(rp) :: Le, u_global
    real(rp), allocatable :: Delta(:)

    ! Local design values
    h = design%get_values()
    n = design%size()

    ! Project design variables to physical height h = (h_max - h_min)*x + h_min
    call vector_cmult(h, (h_max - h_min), n)
    call vector_cadd(h, h_min, n)

    ! Global element length
    Le = L_total / real(design%size_global(), kind=rp)

    ! Exclusive prefix sum for global offset
    call MPI_Exscan(n, offset, 1, MPI_INTEGER, MPI_SUM, neko_comm, ierr)
    if (pe_rank == 0) offset = 0

    ! Build elementwise Δ_k on host (geometry-related, not design-related)
    allocate(Delta(n))
    do k = 1, n
      Delta(k) = ((L_total - Le*real(offset+k-1,rp))**3 - &
                  (L_total - Le*real(offset+k,rp))**3) / 3.0_rp
    end do

    call I%init(n)

    ! Now compute I_k = b * h^3 / 12
    call vector_copy(I, h, n)             ! I = h
    call vector_col2(I, h, n)             ! I = h^2
    call vector_col2(I, h, n)             ! I = h^3
    call vector_cmult(I, b/12.0_rp, n)    ! I = b * h^3 / 12

    ! contrib = Δ / (E * I)
    call contrib%init(n)
    
    contrib%x = Delta                    ! copy Δ into contrib
    if (neko_bcknd_device .eq. 1) then
       call device_memcpy(contrib%x,contrib%x_d, n, &
            HOST_TO_DEVICE, sync = .true.)
    end if
    call vector_invcol2(contrib, I, n)   ! contrib = contrib / I
    call vector_cmult(contrib, P/E, n)   ! contrib = contrib * P / E

    ! Global sum
    u_global = vector_glsum(contrib, n)

    ! if (pe_rank == 0) print *, "Global tip deflection =", u_global

    ! Normalizing with u_tip_max
    u_global = u_global/u_tip_max - 1

    this%value = u_global

    ! if (pe_rank == 0) print *, "Normalized Global tip deflection =", this%value

    deallocate(Delta)
  end subroutine mma_obj_update_value
 

  subroutine mma_obj_update_sensitivity(this, design)
    class(mma_obj), intent(inout) :: this
    class(design_t), intent(in) :: design

    real(rp) :: Le
    type(vector_t) :: h, sensitivity
    integer :: ierr, n, offset, k
    real(rp), allocatable :: Delta(:)

    ! Local design values
    h = design%get_values()
    n = design%size()
    call sensitivity%init(n)
    
    ! Project design variables to physical height
    call vector_cmult(h, (h_max - h_min), n)
    call vector_cadd(h, h_min, n)

    ! Global element length
    Le = L_total / real(design%size_global(), kind=rp)

    ! Exclusive prefix sum for global offset
    call MPI_Exscan(n, offset, 1, MPI_INTEGER, MPI_SUM, neko_comm, ierr)
    if (pe_rank == 0) offset = 0

    ! Build elementwise Δ_k
    allocate(Delta(n))
    do k = 1, n
      Delta(k) = ((L_total - Le*real(offset+k-1,rp))**3 - &
                  (L_total - Le*real(offset+k,rp))**3) / 3.0_rp
    end do

    ! Compute sensitivity: 
    ! dg/dx = P * Δ * (-36.0 / (E * b * h^4)) * (h_max - h_min)
    do k = 1, n
      sensitivity%x(k) = P * Delta(k) * (-36.0_rp) / &
           (E * b * h%x(k)**4) * (h_max - h_min)
    end do
    if (neko_bcknd_device .eq. 1) then
       call device_memcpy(sensitivity%x,sensitivity%x_d, n, &
            HOST_TO_DEVICE, sync = .true.)
    end if

    ! Normalize by u_tip_max
    call vector_cmult(sensitivity, 1.0_rp/u_tip_max, n)
    
    call vector_copy(this%sensitivity, sensitivity, n)

    deallocate(Delta)
    call sensitivity%free()
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

    call this%init_base(name, design%size(), weight)
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

    ! Global element length
    Le = L_total / real(design%size_global(), kind=rp)

    global_mass =  rho * b * Le * vector_glsum(design_values, n)

    this%value = global_mass

    ! if (pe_rank == 0) print *, "Global beam weight =", this%value

  end subroutine beamweight_update_value

  subroutine beamweight_update_sensitivity(this, design)
    class(beamweight_obj), intent(inout) :: this
    class(design_t), intent(in) :: design
    integer :: n
    real(rp) :: Le

    n = design%size()
    Le = L_total / real(design%size_global(), kind=rp)

    ! Sensitivity: dm/dh_k = rho * b * Le
    this%sensitivity%x = rho * b * Le * (h_max - h_min)
    if (neko_bcknd_device .eq. 1) then
       call device_memcpy(this%sensitivity%x,this%sensitivity%x_d, n, &
            HOST_TO_DEVICE, sync = .false.)
    end if
  end subroutine beamweight_update_sensitivity

  ! ========================================================================== !
  ! Methods for the Constraint Function
  subroutine init_stress_con(this, name, design, global_index, sigma_max_value)
    class(stress_con), intent(inout) :: this
    character(len=*), intent(in) :: name
    class(design_t), intent(in) :: design
    integer, intent(in) :: global_index
    real(rp), intent(in) :: sigma_max_value
    
    integer :: ierr, n, offset, local_offset
    
    call this%init_base(name, design%size())  ! 1 constraint per object
    
    this%global_element_index = global_index
    this%sigma_max = sigma_max_value
    
    ! Get MPI distribution info
    n = design%size()
    call MPI_Exscan(n, offset, 1, MPI_INTEGER, MPI_SUM, neko_comm, ierr)
    if (pe_rank == 0) offset = 0
    
    ! Check if this global element is on our rank
    this%is_local = (global_index > offset .and. global_index <= offset + n)
    if (this%is_local) then
      this%local_index = global_index - offset
    else
      this%local_index = -1
    endif
    ! print *, "rank=", pe_rank, "global_index=", global_index, "this%is_local=", this%is_local, &
    !      "this%local_index =", this%local_index 
  end subroutine init_stress_con

  subroutine stress_con_free(this)
    class(stress_con), intent(inout) :: this
    call this%free_base()
  end subroutine stress_con_free

  subroutine stress_con_update_value(this, design)
    class(stress_con), intent(inout) :: this
    class(design_t), intent(in) :: design
    
    type(vector_t) :: h
    integer :: ierr, n, offset
    real(rp) :: Le, x_e, I_e, c_e, M_e, sigma_e
    real(rp) :: global_value
    
    ! Initialize to safe value
    this%value = 0.0_rp
    if (this%is_local) then
       ! This element is  on our rank
       ! Local design values
       h = design%get_values()
       n = design%size()

       ! Project design variables to physical height
       call vector_cmult(h, (h_max - h_min), n)
       call vector_cadd(h, h_min, n)
       
       ! Global element length
       Le = L_total / real(design%size_global(), kind=rp)
       
       ! Start coordinate of this element
       x_e = Le * real(this%global_element_index - 1, rp)
       
       ! Section properties
       I_e = b * h%x(this%local_index)**3 / 12.0_rp
       c_e = h%x(this%local_index) / 2.0_rp
       
       ! Bending moment at this element
       M_e = P * (L_total - x_e)
       
       ! Stress in this element
       sigma_e = M_e * c_e / I_e
       
       ! Constraint value: stress <= sigma_max
       this%value = sigma_e / this%sigma_max - 1.0_rp
    endif

    ! Sum across all MPI ranks (only one rank will have non-zero value)
    call MPI_Allreduce(this%value, global_value, 1, mpi_real_precision, &
                      MPI_SUM, neko_comm, ierr)
    ! print *, "pe_rank=", pe_rank, this%name, "value=", this%value, "global_value=", global_value
    
    this%value = global_value
  end subroutine stress_con_update_value

  subroutine stress_con_update_sensitivity(this, design)
    class(stress_con), intent(inout) :: this
    class(design_t), intent(in) :: design
    
    type(vector_t) :: h
    integer :: ierr, n
    real(rp) :: Le, x_e, I_e, c_e, M_e, dsigma_dh
    real(rp), allocatable :: local_sensitivity(:)
    
    ! Initialize local sensitivity to zero
    if (allocated(local_sensitivity)) deallocate(local_sensitivity)
    allocate(local_sensitivity(design%size_global()))
    local_sensitivity = 0.0_rp
    
    if (this%is_local) then
      ! This element is on our rank
      ! Local design values
      h = design%get_values()
      n = design%size()
      
      ! Project design variables to physical height
      call vector_cmult(h, (h_max - h_min), n)
      call vector_cadd(h, h_min, n)
      
      ! Global element length
      Le = L_total / real(design%size_global(), kind=rp)
      
      ! Start coordinate of this element
      x_e = Le * real(this%global_element_index - 1, rp)
      
      ! Section properties
      I_e = b * h%x(this%local_index)**3 / 12.0_rp
      c_e = h%x(this%local_index) / 2.0_rp
      
      ! Bending moment at this element
      M_e = P * (L_total - x_e)
      
      ! Sensitivity wrt h
      dsigma_dh = M_e * ( (1.0_rp / (2.0_rp * I_e)) - &
          (c_e * 3.0_rp * b * h%x(this%local_index)**2 / 12.0_rp) / (I_e**2) )
      
      ! Chain rule for normalized variable
      local_sensitivity(this%local_index) = dsigma_dh * (h_max - h_min) / this%sigma_max
    endif
    
    ! If the constraint is not local we dont need its sensitivity on the current node
    ! call MPI_Allreduce(local_sensitivity, this%sensitivity%x, design%size_global(), &
    !                   mpi_real_precision, MPI_SUM, neko_comm, ierr)
    this%sensitivity%x = local_sensitivity

    ! Update device memory if needed
    if (neko_bcknd_device .eq. 1) then
      call device_memcpy(this%sensitivity%x, this%sensitivity%x_d, design%size_global(), &
                        HOST_TO_DEVICE, sync = .false.)
    end if
    ! print *, "pe_rank=", pe_rank, this%name, "maxval(this%sensitivity%x)=", maxval(this%sensitivity%x),&
    !      "sum(this%sensitivity%x)=", sum(this%sensitivity%x), "minval(this%sensitivity%x)=", minval(this%sensitivity%x)

    deallocate(local_sensitivity)
  end subroutine stress_con_update_sensitivity

  ! subroutine mma_con_init_from_json(this, json, design)
  !   class(mma_con), intent(inout) :: this
  !   type(json_file), intent(inout) :: json
  !   class(design_t), intent(in) :: design
  !   character(len=256), parameter :: name = 'mma_con'
  !   integer :: sign

  !   call this%init_from_components(name, design, sign)

  ! end subroutine mma_con_init_from_json

  ! subroutine mma_con_init_from_components(this, name, design, sign)
  !   class(mma_con), intent(inout) :: this
  !   character(len=*), intent(in) :: name
  !   class(design_t), intent(in) :: design
  !   integer, intent(in), optional :: sign
  !   integer :: sign_ = 1

  !   call this%init_base(name, design%size())

  !   if (present(sign)) sign_ = sign
  !   if (sign_ .lt. 0) this%sign = -1.0_rp
  !   if (sign_ .ge. 0) this%sign = 1.0_rp

  ! end subroutine mma_con_init_from_components

  ! subroutine mma_con_free(this)
  !   class(mma_con), intent(inout) :: this
  !   call this%free_base()
  ! end subroutine mma_con_free

  ! subroutine mma_con_update_value(this, design)
  !   class(mma_con), intent(inout) :: this
  !   class(design_t), intent(in) :: design


  ! end subroutine mma_con_update_value

  ! subroutine mma_con_update_sensitivity(this, design)
  !   class(mma_con), intent(inout) :: this
  !   class(design_t), intent(in) :: design

  ! end subroutine mma_con_update_sensitivity
end module example_problem