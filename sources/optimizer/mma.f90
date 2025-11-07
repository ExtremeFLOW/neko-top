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
! A primal-dual algorithm is then employed to solve the aproximated problem !
! using interior point method.                                              !
!===========================================================================!

module mma
  ! Inclusions from Neko
  use num_types, only: rp
  use json_module, only: json_file
  use json_utils, only: json_get_or_default
  use vector, only: vector_t
  use matrix, only: matrix_t
  use mpi_f08, only: MPI_Allreduce, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD
  use comm, only: pe_rank
  use utils, only: neko_error
  use neko_config, only: NEKO_BCKND_DEVICE, NEKO_BCKND_CUDA, NEKO_BCKND_HIP, &
       NEKO_BCKND_OPENCL
  use device, only: device_memcpy, HOST_TO_DEVICE, DEVICE_TO_HOST
  use, intrinsic :: iso_c_binding, only: c_ptr

  implicit none
  private

  type, public :: mma_t
     private
     integer :: n, m, max_iter
     real(kind=rp) :: a0, f0val, asyinit, asyincr, asydecr, epsimin, &
          residumax, residunorm
     type(vector_t) :: xold1, xold2, low, upp, alpha, beta, a, c, d, xmax, xmin
     logical :: is_initialized = .false.
     logical :: is_updated = .false.
     character(len=:), allocatable :: bcknd, subsolver

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

     generic, public :: KKT => KKT_vector, KKT_cpu, KKT_device
     procedure, pass(this) :: KKT_vector => mma_KKT_vector
     procedure, pass(this) :: KKT_cpu => mma_KKT_cpu
     procedure, pass(this) :: KKT_device => mma_KKT_device

  end type mma_t

  interface
     ! ======================================================================= !
     ! interface for cpu backend module subroutines

     !> CPU update function, runs one iteration of MMA
     module subroutine mma_update_cpu(this, iter, x, df0dx, fval, dfdx)
       class(mma_t), intent(inout) :: this
       integer, intent(in) :: iter
       real(kind=rp), dimension(this%n), intent(inout) :: x
       real(kind=rp), dimension(this%n), intent(in) :: df0dx
       real(kind=rp), dimension(this%m), intent(in) :: fval
       real(kind=rp), dimension(this%m, this%n), intent(in) :: dfdx
     end subroutine mma_update_cpu

     !> CPU KKT check for convergence
     module subroutine mma_KKT_cpu(this, x, df0dx, fval, dfdx)
       class(mma_t), intent(inout) :: this
       real(kind=rp), dimension(this%n), intent(in) :: x
       real(kind=rp), dimension(this%n), intent(in) :: df0dx
       real(kind=rp), dimension(this%m), intent(in) :: fval
       real(kind=rp), dimension(this%m, this%n), intent(in) :: dfdx
     end subroutine mma_KKT_cpu

     ! ======================================================================= !
     ! interface for device backend module subroutines

     !> Device update function, runs one iteration of MMA
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

contains

  !> Read attributes from the case file, and calling the init function
  subroutine mma_init_from_json(this, x, n, m, json, scale, auto_scale)
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

    integer :: max_iter, n_global, ierr
    real(kind=rp) :: epsimin, asyinit, asyincr, asydecr

    call MPI_Allreduce(n, n_global, 1, MPI_INTEGER, &
         MPI_SUM, MPI_COMM_WORLD, ierr)

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
    call json_get_or_default(json, 'mma.max_iter', max_iter, 100)

    ! Following parameters are set based on eq.3.8:--------
    call json_get_or_default(json, 'mma.asyinit', asyinit, 0.5_rp)
    call json_get_or_default(json, 'mma.asyincr', asyincr, 1.2_rp)
    call json_get_or_default(json, 'mma.asydecr', asydecr, 0.7_rp)

    call json_get_or_default(json, 'mma.backend', bcknd, bcknd_default)
    call json_get_or_default(json, 'mma.subsolver', subsolver, 'dip')

    call json_get_or_default(json, 'mma.xmin', xmin_const, 0.0_rp)
    call json_get_or_default(json, 'mma.xmax', xmax_const, 1.0_rp)
    call json_get_or_default(json, 'mma.a0', a0, 1.0_rp)
    call json_get_or_default(json, 'mma.a', a_const, 0.0_rp)
    call json_get_or_default(json, 'mma.c', c_const, 100.0_rp)
    call json_get_or_default(json, 'mma.d', d_const, 0.0_rp)

    call json_get_or_default(json, 'mma.scale', scale, 10.0_rp)
    call json_get_or_default(json, 'mma.auto_scale', auto_scale, .false.)

    ! Initialize the MMA object with the parsed parameters
    a = a_const
    c = c_const
    d = d_const
    xmin = xmin_const
    xmax = xmax_const
    ! initializing the mma concrete type (mma_cpu_t or mma_device_t)
    if (pe_rank .eq. 0) then
       print *, "Initializing MMA backend to >>> ", bcknd
    end if

    ! ------------------------------------------------------------------------ !
    ! Initialize the MMA object with the parameters read from json

    call this%init(x, n, m, a0, a, c, d, xmin, xmax, &
         max_iter, epsimin, asyinit, asyincr, asydecr, bcknd, subsolver)

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

    this%is_initialized = .false.
    this%is_updated = .false.
  end subroutine mma_free

  !> Initialize the mma object based on the attributes from the json file
  subroutine mma_init_from_components(this, x, n, m, a0, a, c, d, xmin, xmax, &
       max_iter, epsimin, asyinit, asyincr, asydecr, bcknd, subsolver)
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
    character(len=:), intent(in), allocatable :: bcknd, subsolver

    call this%free()

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
       call device_memcpy(this%a%x, this%a%x_d, m, HOST_TO_DEVICE, &
            sync = .false.)
       call device_memcpy(this%c%x, this%c%x_d, m, HOST_TO_DEVICE, &
            sync = .false.)
       call device_memcpy(this%d%x, this%d%x_d, m, HOST_TO_DEVICE, &
            sync = .false.)
       call device_memcpy(this%xmax%x, this%xmax%x_d, n, HOST_TO_DEVICE, &
            sync = .false.)
       call device_memcpy(this%xmin%x, this%xmin%x_d, n, HOST_TO_DEVICE, &
            sync = .true.)
    end if


    ! Set KKT norms to a large number for the initial design
    this%residumax = huge(0.0_rp)
    this%residunorm = huge(0.0_rp)

    ! Select backend type
    select case (bcknd)
    case ("cpu")
       if (pe_rank == 0) then
          print *, "MMA initialized with CPU backend!"
       end if
    case ("device")
       if (pe_rank == 0) then
          if (NEKO_BCKND_CUDA .eq. 1) then
             print *, "MMA initialized with CUDA backend!"
          else if (NEKO_BCKND_HIP .eq. 1) then
             print *, "MMA initialized with HIP backend!"
          else if (NEKO_BCKND_OPENCL .eq. 1) then
             print *, "MMA initialized with OPENCL backend!"
          else
             call neko_error('Unknown backend device in mma_init_components')
          end if
       end if
    case default
       call neko_error('Unknown backend in mma_init_components')
    end select

    ! ------------------------------------------------------------------------ !
    ! Assign defaults if nothing is parsed

    ! Based on the Cpp Code by Niels
    if (.not. present(epsimin)) this%epsimin = 1.0e-9_rp * sqrt(real(m + n, rp))
    if (.not. present(max_iter)) this%max_iter = 100

    ! Following parameters are set based on eq.3.8
    if (.not. present(asyinit)) this%asyinit = 0.5_rp
    if (.not. present(asyincr)) this%asyincr = 1.2_rp
    if (.not. present(asydecr)) this%asydecr = 0.7_rp

    ! Assign values from inputs when present
    if (present(max_iter)) this%max_iter = max_iter
    if (present(epsimin)) this%epsimin = epsimin
    if (present(asyinit)) this%asyinit = asyinit
    if (present(asyincr)) this%asyincr = asyincr
    if (present(asydecr)) this%asydecr = asydecr
    this%bcknd = bcknd
    this%subsolver = subsolver
    if (pe_rank == 0) then
       if (this%subsolver .eq. "dip") then
          print *, "Using dual solver for MMA subsolve."
       elseif (this%subsolver .eq. "dpip") then
          print *, "Using dual-primal solver for MMA subsolve."
       else
          call neko_error('Unknown subsolver for MMA, mma_init_from_components')
       end if
    end if

    if (pe_rank .eq. 0) then
       print *, "MMA is initialized with a0 = ", a0, ", a = ", a, ", c = ", c, &
            ", d = ", d, "epsimin = ", this%epsimin
    end if
    ! The object is correctly initialized
    this%is_initialized = .true.
  end subroutine mma_init_from_components

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
          call device_memcpy(x%x, x%x_d, this%n, DEVICE_TO_HOST, &
               sync = .false.)
          call device_memcpy(df0dx%x, df0dx%x_d, this%n, DEVICE_TO_HOST, &
               sync = .false.)
          call device_memcpy(fval%x, fval%x_d, this%m, DEVICE_TO_HOST, &
               sync = .false.)
          call device_memcpy(dfdx%x, dfdx%x_d, this%m * this%n, DEVICE_TO_HOST,&
               sync = .true.)
       end if

       call mma_update_cpu(this, iter, x%x, df0dx%x, fval%x, dfdx%x)

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_memcpy(x%x, x%x_d, this%n, HOST_TO_DEVICE, sync = .true.)
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
    backend_subsolver = ', backend:' // this%bcknd // ',subsolver:' // &
         this%subsolver
  end function mma_get_backend_and_subsolver
end module mma
