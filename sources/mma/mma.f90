!> @file mma.f90
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

!===========================================================================!
!                       Method of Moving Asymptotes                         !
! This implementation is based on the following documents:                  !
!        1-https://people.kth.se/~krille/mmagcmma.pdf                       !
!        2-https://people.kth.se/~krille/originalmma.pdf                    !
!        2-https://comsolyar.com/wp-content/uploads/2020/03/gcmma.pdf       !
! ------------------------------------------------------------------------- !
!                                                                           !
! This module solves the following original optimization problem:           !
!                                                                           !
!      Minimize  f_0(x) + a_0*z + sum( c_i*y_i + 0.5*d_i*(y_i)^2 )          !
!    subject to  f_i(x) - a_i*z - y_i <= 0,  i = 1,...,m                    !
!                xmin_j <= x_j <= xmax_j,    j = 1,...,n                    !
!                z >= 0,   y_i >= 0,         i = 1,...,m                    !
!                                                                           !
! by first creating the following convex approximation of the original      !
! problem:                                                                  !
!                                                                           !
!      Minimize sum_{j = 1,...,n} (p0j / (upp_j-x_j) + q0j / (x_j-low_j)) + !
!                        a0*z + sum_i = 1,...,m(c_i*y_i + 0.5*d_i*y_i^2)    !
!    subject to sum_{j = 1,...,n} (pij / (upp_j-x_j) + qij / (x_j-low_j)) + !
!                    a_i*z + y_i <= b_i,                       i = 1,...,m  !
!               xmin_j <= alpha_j <= x_j <= beta_j <= xmax_j   j = 1,...,n  !
!               y_i>=0                                         i = 1,...,m  !
!               z>=0                                                        !
!                                                                           !
! note that based on eq(3.5) there should be r0 in the approximated problem !
! however since it is just a constant added to a minimization problem, it   !
! is ignored.                                                               !
! A primal-dual algorithm is then employed to solve the approximated        !
! problem using interior point method.                                      !
!===========================================================================!

!> MMA module
module mma
  ! Inclusions from Neko
  use num_types, only: rp, dp, sp
  use json_module, only: json_file
  use json_utils, only: json_get_or_default
  use vector, only: vector_t
  use matrix, only: matrix_t
  use comm, only: pe_rank, NEKO_COMM, pe_size, MPI_REAL_PRECISION
  use utils, only: neko_error, filename_suffix
  use neko_config, only: NEKO_BCKND_DEVICE, NEKO_BCKND_CUDA, NEKO_BCKND_HIP, &
       NEKO_BCKND_OPENCL
  use device, only: device_memcpy, HOST_TO_DEVICE, DEVICE_TO_HOST
  use, intrinsic :: iso_c_binding, only: c_ptr
  use logger, only: neko_log
  use mpi_f08, only: MPI_SUM, MPI_Allreduce, MPI_INTEGER
  use scratch_registry, only: scratch_registry_t

  implicit none
  private

  !> MMA type
  type, public :: mma_t
     private
     integer :: n, m, n_global, max_iter
     real(kind=rp) :: a0, asyinit, asyincr, asydecr, epsimin, &
          residumax, residunorm
     real(kind=rp) :: move_limit = 0.2_rp
     type(vector_t) :: xold1, xold2, low, upp, alpha, beta, a, c, d, xmax, xmin
     logical :: is_initialized = .false.
     logical :: is_updated = .false.
     logical :: unconstrained_problem = .false.
     type(scratch_registry_t) :: scratch
     character(len=:), allocatable :: subsolver, bcknd

     ! Internal dummy variables for MMA
     type(vector_t) :: p0j, q0j
     type(matrix_t) :: pij, qij
     type(vector_t) :: bi

     !---nesessary for KKT check after updating df0dx, fval, dfdx --------
     real(kind=rp) :: z, zeta
     type(vector_t) :: y, lambda, s, mu
     type(vector_t) :: xsi, eta
   contains
     !> Interface for initializing the MMA object
     generic, public :: init => init_from_json, init_from_components
     procedure, public, pass(this) :: init_from_json => mma_init_from_json
     procedure, public, pass(this) :: init_from_components => &
          mma_init_from_components
     procedure, public, pass(this) :: free => mma_free
     procedure, public, pass(this) :: get_n => mma_get_n
     procedure, public, pass(this) :: get_m => mma_get_m
     procedure, public, pass(this) :: get_residumax => mma_get_residumax
     procedure, public, pass(this) :: get_residunorm => mma_get_residunorm
     procedure, public, pass(this) :: get_max_iter => mma_get_max_iter
     procedure, public, pass(this) :: get_backend_and_subsolver => &
          mma_get_backend_and_subsolver

     generic, public :: update => update_vector, update_cpu, update_device
     procedure, pass(this) :: update_vector => mma_update_vector
     procedure, pass(this) :: update_cpu => mma_update_cpu
     procedure, pass(this) :: update_device => mma_update_device

     generic, public :: kkt => KKT_vector, KKT_cpu, KKT_device
     procedure, pass(this) :: KKT_vector => mma_KKT_vector
     procedure, pass(this) :: KKT_cpu => mma_KKT_cpu
     procedure, pass(this) :: KKT_device => mma_KKT_device

     procedure, pass(this) :: save_checkpoint => mma_save_checkpoint
     procedure, pass(this) :: load_checkpoint => mma_load_checkpoint

     ! Private utilities
     procedure, pass(this) :: copy_from => mma_copy_from

  end type mma_t

  ! ========================================================================== !
  ! Default parameters

  real(kind=rp), parameter :: a0_default = 1.0_rp
  real(kind=rp), parameter :: a_default = 0.0_rp
  real(kind=rp), parameter :: c_default = 100.0_rp
  real(kind=rp), parameter :: d_default = 0.0_rp
  real(kind=rp), parameter :: xmin_default = 0.0_rp
  real(kind=rp), parameter :: xmax_default = 1.0_rp

  real(kind=rp), parameter :: asyinit_default = 0.2_rp
  real(kind=rp), parameter :: asyincr_default = 1.05_rp
  real(kind=rp), parameter :: asydecr_default = 0.65_rp
  real(kind=rp), parameter :: move_limit_default = 0.2_rp

  integer, parameter :: max_iter_default = 100
  character(len=*), parameter :: subsolver_default = "dip"
  real(kind=rp), parameter :: scale_default = 1.0_rp
  logical, parameter :: auto_scale_default = .false.

  ! ========================================================================== !
  ! interface for cpu backend module subroutines

  interface
     ! CPU update function, runs one iteration of MMA
     module subroutine mma_update_cpu(this, iter, x, df0dx, fval, dfdx)
       class(mma_t), intent(inout) :: this
       integer, intent(in) :: iter
       real(kind=rp), dimension(this%n), intent(inout) :: x
       real(kind=rp), dimension(this%n), intent(in) :: df0dx
       real(kind=rp), dimension(this%m), intent(in) :: fval
       real(kind=rp), dimension(this%m, this%n), intent(in) :: dfdx
     end subroutine mma_update_cpu

     ! CPU KKT check for convergence
     module subroutine mma_KKT_cpu(this, x, df0dx, fval, dfdx)
       class(mma_t), intent(inout) :: this
       real(kind=rp), dimension(this%n), intent(in) :: x
       real(kind=rp), dimension(this%n), intent(in) :: df0dx
       real(kind=rp), dimension(this%m), intent(in) :: fval
       real(kind=rp), dimension(this%m, this%n), intent(in) :: dfdx
     end subroutine mma_KKT_cpu

     ! ========================================================================== !
     ! interface for device backend module subroutines

     ! Device update function, runs one iteration of MMA
     module subroutine mma_update_device(this, iter, x, df0dx, fval, dfdx)
       class(mma_t), intent(inout) :: this
       integer, intent(in) :: iter
       type(c_ptr), intent(inout) :: x
       type(c_ptr), intent(in) :: df0dx, fval, dfdx
     end subroutine mma_update_device

     !> Device KKT check for convergence
     module subroutine mma_KKT_device(this, x, df0dx, fval, dfdx)
       class(mma_t), intent(inout) :: this
       type(c_ptr), intent(in) :: x, df0dx, fval, dfdx
     end subroutine mma_KKT_device

  end interface

  ! ========================================================================== !
  ! Interface for IO routines

  interface
     module subroutine mma_save_checkpoint_hdf5(object, filename, overwrite)
       class(mma_t), intent(inout) :: object
       character(len=*), intent(in) :: filename
       logical, intent(in), optional :: overwrite
     end subroutine mma_save_checkpoint_hdf5

     module subroutine mma_load_checkpoint_hdf5(object, filename)
       class(mma_t), intent(inout) :: object
       character(len=*), intent(in) :: filename
     end subroutine mma_load_checkpoint_hdf5
  end interface

contains

  ! ========================================================================== !
  ! Initializers and destructors

  !> Read attributes from the case file, and calling the init function
  subroutine mma_init_from_json(this, x, n, m, json, scale, auto_scale, &
       unconstrained_problem)
    ! ----------------------------------------------------- !
    ! Initializing the mma object and all the parameters    !
    ! required for MMA method. (a_i, c_i, d_i, ...)         !
    ! x: the design varaibles(DV), n: number of DV,         !
    ! m: number of constraints                              !
    !                                                       !
    ! Note that residumax & residunorm of the KKT conditions!
    ! are initialized with 10^5. This is done to avoid      !
    ! unnecessary extera computation of KKT norms for the   !
    ! initial design.                                       !
    ! ----------------------------------------------------- !
    class(mma_t), intent(inout) :: this
    integer, intent(in) :: n, m
    type(vector_t), intent(in) :: x

    type(json_file), intent(inout) :: json

    ! Read the scaling info for fval and dfdx from json
    real(kind=rp), intent(out) :: scale
    logical, intent(out) :: auto_scale
    logical, intent(in) :: unconstrained_problem
    ! -------------------------------------------------------------------!
    !      Internal parameters for MMA                                   !
    !      Minimize  f_0(x) + a_0*z + sum( c_i*y_i + 0.5*d_i*(y_i)^2 )   !
    !    subject to  f_i(x) - a_i*z - y_i <= 0,  i = 1,...,m             !
    !                xmin_j <= x_j <= xmax_j,    j = 1,...,n             !
    !                z >= 0,   y_i >= 0,         i = 1,...,m             !
    ! -------------------------------------------------------------------!
    real(kind=rp), dimension(n) :: xmax, xmin
    real(kind=rp), dimension(m) :: a, c, d
    character(len=:), allocatable :: subsolver, bcknd, bcknd_default

    ! For reading the values from json and then set the value for the arrays
    real(kind=rp) :: a0 , xmax_const, xmin_const, a_const, c_const, d_const
    real(kind=rp) :: move_limit

    integer :: max_iter, n_global, ierr
    real(kind=rp) :: epsimin, asyinit, asyincr, asydecr

    call MPI_Allreduce(n, n_global, 1, MPI_INTEGER, &
         MPI_SUM, NEKO_COMM, ierr)

    ! Assign default values for the backend based on the NEKO_BCKND_DEVICE
    if (NEKO_BCKND_DEVICE .eq. 1) then
       bcknd_default = "device"
    else
       bcknd_default = "cpu"
    end if

    ! ------------------------------------------------------------------------ !
    ! Assign defaults if nothing is parsed
    ! based on the Cpp Code by Niels
    call json_get_or_default(json, 'mma.epsimin', epsimin, &
         1.0e-9_rp * sqrt(real(m + n_global, rp)))
    call json_get_or_default(json, 'mma.max_iter', max_iter, max_iter_default)

    ! Following parameters are set based on eq.3.8:--------
    call json_get_or_default(json, 'mma.asyinit', asyinit, asyinit_default)
    call json_get_or_default(json, 'mma.asyincr', asyincr, asyincr_default)
    call json_get_or_default(json, 'mma.asydecr', asydecr, asydecr_default)

    call json_get_or_default(json, 'mma.backend', bcknd, bcknd_default)
    call json_get_or_default(json, 'mma.subsolver', subsolver, subsolver_default)

    call json_get_or_default(json, 'mma.xmin', xmin_const, xmin_default)
    call json_get_or_default(json, 'mma.xmax', xmax_const, xmax_default)
    call json_get_or_default(json, 'mma.a0', a0, a0_default)
    call json_get_or_default(json, 'mma.a', a_const, a_default)
    call json_get_or_default(json, 'mma.c', c_const, c_default)
    call json_get_or_default(json, 'mma.d', d_const, d_default)
    call json_get_or_default(json, 'mma.move_limit', move_limit, move_limit_default)

    call json_get_or_default(json, 'mma.scale', scale, scale_default)
    call json_get_or_default(json, 'mma.auto_scale', auto_scale, auto_scale_default)

    ! Initialize the MMA object with the parsed parameters
    a = a_const
    c = c_const
    d = d_const
    xmin = xmin_const
    xmax = xmax_const

    ! ------------------------------------------------------------------------ !
    ! Initialize the MMA object with the parameters read from json

    call this%init(x, n, m, a0, a, c, d, xmin, xmax, &
         max_iter, epsimin, asyinit, asyincr, asydecr, bcknd, subsolver, &
         move_limit, unconstrained_problem)

  end subroutine mma_init_from_json

  !> Deallocate mma object
  subroutine mma_free(this)
    class(mma_t), intent(inout) :: this
    ! Deallocate the internal vectors
    call this%xold1%free()
    call this%xold2%free()
    call this%alpha%free()
    call this%beta%free()
    call this%a%free()
    call this%c%free()
    call this%d%free()
    call this%low%free()
    call this%upp%free()
    call this%xmax%free()
    call this%xmin%free()
    call this%p0j%free()
    call this%q0j%free()
    call this%bi%free()
    call this%y%free()
    call this%lambda%free()
    call this%s%free()
    call this%mu%free()
    call this%xsi%free()
    call this%eta%free()

    ! Deallocate the internal dummy matrices
    call this%pij%free()
    call this%qij%free()
    call this%scratch%free()

    this%is_initialized = .false.
    this%is_updated = .false.
  end subroutine mma_free

  !> Initialize the mma object based on the attributes from the json file
  subroutine mma_init_from_components(this, x, n, m, a0, a, c, d, xmin, xmax, &
       max_iter, epsimin, asyinit, asyincr, asydecr, bcknd, subsolver, &
       move_limit, unconstrained_problem)
    ! ----------------------------------------------------- !
    ! Initializing the mma object and all the parameters    !
    ! required for MMA method. (a_i, c_i, d_i, ...)         !
    ! x: the design varaibles(DV), n: number of DV,         !
    ! m: number of constraints                              !
    !                                                       !
    ! Note that residumax & residunorm of the KKT conditions!
    ! are initialized with 10^5. This is done to avoid      !
    ! unnecessary extera computation of KKT norms for the   !
    ! initial design.                                       !
    ! ----------------------------------------------------- !
    class(mma_t), intent(inout) :: this
    integer, intent(in) :: n, m
    type(vector_t), intent(in) :: x
    ! -------------------------------------------------------------------!
    !      Internal parameters for MMA                                   !
    !      Minimize  f_0(x) + a_0*z + sum( c_i*y_i + 0.5*d_i*(y_i)^2 )   !
    !    subject to  f_i(x) - a_i*z - y_i <= 0,  i = 1,...,m             !
    !                xmin_j <= x_j <= xmax_j,    j = 1,...,n             !
    !                z >= 0,   y_i >= 0,         i = 1,...,m             !
    ! -------------------------------------------------------------------!
    real(kind=rp), intent(in), dimension(n) :: xmax, xmin
    real(kind=rp), intent(in), dimension(m) :: a, c, d
    real(kind=rp), intent(in) :: a0
    integer, intent(in), optional :: max_iter
    real(kind=rp), intent(in), optional :: epsimin, asyinit, asyincr, asydecr
    real(kind=rp), intent(in), optional :: move_limit
    character(len=*), intent(in), optional :: bcknd, subsolver
    logical, intent(in) :: unconstrained_problem
    character(len=256) :: log_msg
    integer :: i, ierr

    call this%free()
    call this%scratch%init()

    this%unconstrained_problem = unconstrained_problem
    this%n = n
    this%m = m

    call this%xold1%init(n)
    call this%xold2%init(n)
    this%xold1 = x
    this%xold2 = x

    call this%alpha%init(n)
    call this%beta%init(n)

    call this%a%init(m)
    call this%c%init(m)
    call this%d%init(m)
    call this%low%init(n)
    call this%upp%init(n)
    call this%xmax%init(n)
    call this%xmin%init(n)

    !internal dummy variables for MMA
    call this%p0j%init(n)
    call this%q0j%init(n)
    call this%pij%init(m, n)
    call this%qij%init(m, n)
    call this%bi%init(m)

    ! Necessary for KKT check after updating df0dx, fval, dfdx
    call this%y%init(m)
    call this%lambda%init(m)
    call this%s%init(m)
    call this%mu%init(m)
    call this%xsi%init(n)
    call this%eta%init(n)

    this%a0 = a0
    this%a%x = a
    this%c%x = c
    this%d%x = d

    ! Set the bounds for the design variable based on the problem
    this%xmax%x = xmax
    this%xmin%x = xmin

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call this%a%copy_from(HOST_TO_DEVICE, sync = .false.)
       call this%c%copy_from(HOST_TO_DEVICE, sync = .false.)
       call this%d%copy_from(HOST_TO_DEVICE, sync = .false.)
       call this%xmax%copy_from(HOST_TO_DEVICE, sync = .false.)
       call this%xmin%copy_from(HOST_TO_DEVICE, sync = .true.)
    end if

    ! Set KKT norms to a large number for the initial design
    this%residumax = huge(0.0_rp)
    this%residunorm = huge(0.0_rp)

    ! Sync parameters across MPI
    call MPI_Allreduce(n, this%n_global, 1, MPI_INTEGER, MPI_SUM, neko_comm, &
         ierr)

    ! ------------------------------------------------------------------------ !
    ! Assign defaults if nothing is parsed

    ! Based on the Cpp Code by Niels
    if (.not. present(max_iter)) this%max_iter = max_iter_default
    if (.not. present(epsimin)) then
       this%epsimin = 1.0e-9_rp * sqrt(real(this%m + this%n_global, rp))
    end if

    ! Following parameters are set based on eq.3.8
    if (.not. present(asyinit)) this%asyinit = asyinit_default
    if (.not. present(asyincr)) this%asyincr = asyincr_default
    if (.not. present(asydecr)) this%asydecr = asydecr_default
    if (.not. present(move_limit)) this%move_limit = move_limit_default

    ! Set default backend based on NEKO_BCKND_DEVICE
    if (.not. present(bcknd) .and. NEKO_BCKND_DEVICE .eq. 0) then
       this%bcknd = "cpu"
    else if (.not. present(bcknd)) then
       this%bcknd = "device"
    end if

    ! Set default subsolver
    if (.not. present(subsolver)) this%subsolver = subsolver_default

    ! Assign values from inputs when present
    if (present(max_iter)) this%max_iter = max_iter
    if (present(epsimin)) this%epsimin = epsimin
    if (present(asyinit)) this%asyinit = asyinit
    if (present(asyincr)) this%asyincr = asyincr
    if (present(asydecr)) this%asydecr = asydecr
    if (present(move_limit)) this%move_limit = move_limit
    if (present(bcknd)) this%bcknd = bcknd
    if (present(subsolver)) this%subsolver = subsolver

    call neko_log%section('MMA Parameters')

    write(log_msg, '(A10,1X,A)') 'backend   ', trim(this%bcknd)
    call neko_log%message(log_msg)
    write(log_msg, '(A10,1X,A)') 'subsolver ', trim(this%subsolver)
    call neko_log%message(log_msg)

    write(log_msg, '(A10,1X,I0)') 'n         ', this%n_global
    call neko_log%message(log_msg)
    write(log_msg, '(A10,1X,I0)') 'm         ', this%m
    call neko_log%message(log_msg)
    write(log_msg, '(A10,1X,I0)') 'max_iter  ', this%max_iter
    call neko_log%message(log_msg)

    write(log_msg, '(A10,1X,E11.5)') 'epsimin   ', this%epsimin
    call neko_log%message(log_msg)

    write(log_msg, '(A10,1X,E11.5)') 'asyinit   ', this%asyinit
    call neko_log%message(log_msg)
    write(log_msg, '(A10,1X,E11.5)') 'asyincr   ', this%asyincr
    call neko_log%message(log_msg)
    write(log_msg, '(A10,1X,E11.5)') 'asydecr   ', this%asydecr
    call neko_log%message(log_msg)
    write(log_msg, '(A10,1X,E11.5)') 'a0        ', this%a0
    call neko_log%message(log_msg)
    write(log_msg, '(A10,1X,E11.5)') 'movelimit ', this%move_limit
    call neko_log%message(log_msg)

    call neko_log%message('Parameters a')
    do i = 1, this%m
       write(log_msg, '(3X,A,I2,A,E11.5)') 'a(', i, ') = ', this%a%x(i)
       call neko_log%message(log_msg)
    end do
    call neko_log%message('Parameters c')
    do i = 1, this%m
       write(log_msg, '(3X,A,I2,A,E11.5)') 'c(', i, ') = ', this%c%x(i)
       call neko_log%message(log_msg)
    end do
    call neko_log%message('Parameters d')
    do i = 1, this%m
       write(log_msg, '(3X,A,I2,A,E11.5)') 'd(', i, ') = ', this%d%x(i)
       call neko_log%message(log_msg)
    end do

    call neko_log%end_section()

    ! The object is correctly initialized
    this%is_initialized = .true.
  end subroutine mma_init_from_components

  ! ========================================================================== !
  ! Updator and KKT checker

  !> Call the update function based on the backend
  subroutine mma_update_vector(this, iter, x, df0dx, fval, dfdx)
    class(mma_t), intent(inout) :: this
    integer, intent(in) :: iter
    type(vector_t), intent(inout) :: x
    type(vector_t), intent(inout) :: df0dx, fval
    type(matrix_t), intent(inout) :: dfdx

    ! Select backend type
    select case (this%bcknd)
    case ("cpu")
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call x%copy_from(DEVICE_TO_HOST, sync = .false.)
          call df0dx%copy_from(DEVICE_TO_HOST, sync = .false.)
          call fval%copy_from(DEVICE_TO_HOST, sync = .false.)
          call dfdx%copy_from(DEVICE_TO_HOST, sync = .true.)
       end if

       call mma_update_cpu(this, iter, x%x, df0dx%x, fval%x, dfdx%x)

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call x%copy_from(HOST_TO_DEVICE, sync = .true.)
       end if

    case ("device")
       call mma_update_device(this, iter, x%x_d, df0dx%x_d, fval%x_d, dfdx%x_d)
    end select

  end subroutine mma_update_vector

  !> Call the KKT ckeck function based on the backend
  subroutine mma_KKT_vector(this, x, df0dx, fval, dfdx)
    class(mma_t), intent(inout) :: this
    type(vector_t), intent(inout) :: x, df0dx, fval
    type(matrix_t), intent(inout) :: dfdx

    ! Select backend type
    select case (this%bcknd )
    case ("cpu")
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_memcpy(x%x, x%x_d, this%n, DEVICE_TO_HOST, &
               sync = .false.)
          call device_memcpy(df0dx%x, df0dx%x_d, this%n, DEVICE_TO_HOST, &
               sync = .false.)
          call device_memcpy(fval%x, fval%x_d, this%m, DEVICE_TO_HOST, &
               sync = .false.)
          call device_memcpy(dfdx%x, dfdx%x_d, this%m * this%n, DEVICE_TO_HOST,&
               sync = .true.)
       end if

       call mma_KKT_cpu(this, x%x, df0dx%x, fval%x, dfdx%x)
    case ("device")
       call mma_KKT_device(this, x%x_d, df0dx%x_d, fval%x_d, dfdx%x_d)
    end select
  end subroutine mma_KKT_vector

  ! ========================================================================== !
  ! IO Functions

  !> Save the mma checkpoint to a file based on file suffix.
  !! @param this The mma object.
  !! @param filename The name of the file to save the checkpoint.
  !! @param overwrite Whether to overwrite the file if it exists.
  subroutine mma_save_checkpoint(this, filename, overwrite)
    class(mma_t), intent(inout) :: this
    character(len=*), intent(in) :: filename
    logical, intent(in), optional :: overwrite
    character(len=12) :: file_ext

    ! Get the file extension
    call filename_suffix(filename, file_ext)

    select case (trim(file_ext))
    case ('h5', 'hdf5', 'hf5')
       call mma_save_checkpoint_hdf5(this, filename, overwrite)
    case default
       call neko_error('mma_save_checkpoint: Unsupported file format: ' // &
            trim(file_ext))
    end select

  end subroutine mma_save_checkpoint

  !> Load the mma checkpoint from a file based on file suffix.
  !! @param this The mma object.
  !! @param filename The name of the file to load the checkpoint from.
  subroutine mma_load_checkpoint(this, filename)
    class(mma_t), intent(inout) :: this
    character(len=*), intent(in) :: filename
    character(len=12) :: file_ext

    ! Get the file extension
    call filename_suffix(filename, file_ext)

    select case (trim(file_ext))
    case ('h5', 'hdf5', 'hf5')
       call mma_load_checkpoint_hdf5(this, filename)
    case default
       call neko_error('mma_load_checkpoint: Unsupported file format: ' // &
            trim(file_ext))
    end select
  end subroutine mma_load_checkpoint

  ! ========================================================================== !
  ! Getters and setters

  !> Get the number of design variables (nloc)
  pure function mma_get_n(this) result(n)
    class(mma_t), intent(in) :: this
    integer :: n
    n = this%n
  end function mma_get_n

  !> Get the number of constriants
  pure function mma_get_m(this) result(m)
    class(mma_t), intent(in) :: this
    integer :: m
    m = this%m
  end function mma_get_m

  !> Get L^{inf} norm (Max Norm) of the KKT conditions
  pure function mma_get_residumax(this) result(residumax)
    class(mma_t), intent(in) :: this
    real(kind=rp) :: residumax
    residumax = this%residumax
  end function mma_get_residumax

  !> Get L^{2} norm (Euclidean Norm) of the KKT conditions
  pure function mma_get_residunorm(this) result(residunorm)
    class(mma_t), intent(in) :: this
    real(kind=rp) :: residunorm
    residunorm = this%residunorm
  end function mma_get_residunorm

  !> Get the maximum number of iterations for the mma_subsolve inner loop
  pure function mma_get_max_iter(this) result(max_iter_value)
    class(mma_t), intent(in) :: this
    integer :: max_iter_value
    max_iter_value = this%max_iter
  end function mma_get_max_iter

  !> Get the maximum number of iterations for the mma_subsolve inner loop
  pure function mma_get_backend_and_subsolver(this) result(backend_subsolver)
    class(mma_t), intent(in) :: this
    character(len=:), allocatable :: backend_subsolver
    character(len=:), allocatable :: backend

    if (NEKO_BCKND_CUDA .eq. 1) then
       backend = "cuda"
    else if (NEKO_BCKND_HIP .eq. 1) then
       backend = "hip"
    else if (NEKO_BCKND_OPENCL .eq. 1) then
       backend = "opencl"
    else
       backend = "cpu"
    end if

    backend_subsolver = 'backend:' // trim(backend) // ', subsolver:' // &
         trim(this%subsolver)
  end function mma_get_backend_and_subsolver

  ! ========================================================================== !
  ! Private utilities

  !> Sync device memory to host for all internal vectors/matrices
  subroutine mma_copy_from(this, direction, sync)
    class(mma_t), intent(inout) :: this
    integer, intent(in) :: direction
    logical, intent(in) :: sync

    call this%xold1%copy_from(direction, sync = .false.)
    call this%xold2%copy_from(direction, sync = .false.)
    call this%xmax%copy_from(direction, sync = .false.)
    call this%xmin%copy_from(direction, sync = .false.)

    call this%low%copy_from(direction, sync = .false.)
    call this%upp%copy_from(direction, sync = .false.)

    call this%a%copy_from(direction, sync = .false.)
    call this%c%copy_from(direction, sync = .false.)
    call this%d%copy_from(direction, sync = .false.)
    call this%y%copy_from(direction, sync = .false.)
    call this%s%copy_from(direction, sync = .false.)

    call this%p0j%copy_from(direction, sync = .false.)
    call this%q0j%copy_from(direction, sync = .false.)
    call this%pij%copy_from(direction, sync = .false.)
    call this%qij%copy_from(direction, sync = .false.)
    call this%bi%copy_from(direction, sync = .false.)

    call this%alpha%copy_from(direction, sync = .false.)
    call this%beta%copy_from(direction, sync = .false.)
    call this%lambda%copy_from(direction, sync = .false.)
    call this%mu%copy_from(direction, sync = .false.)
    call this%xsi%copy_from(direction, sync = .false.)
    call this%eta%copy_from(direction, sync = sync)

  end subroutine mma_copy_from

  ! ========================================================================== !
  ! Dummy implementations for module procedures

#if !HAVE_HDF5
  module subroutine mma_save_checkpoint_hdf5(object, filename, overwrite)
    class(mma_t), intent(inout) :: object
    character(len=*), intent(in) :: filename
    logical, intent(in), optional :: overwrite
    call neko_error('mma: HDF5 support not enabled rebuild with HAVE_HDF5')
  end subroutine mma_save_checkpoint_hdf5

  module subroutine mma_load_checkpoint_hdf5(object, filename)
    class(mma_t), intent(inout) :: object
    character(len=*), intent(in) :: filename
    call neko_error('mma: HDF5 support not enabled rebuild with HAVE_HDF5')
  end subroutine mma_load_checkpoint_hdf5
#endif

end module mma
