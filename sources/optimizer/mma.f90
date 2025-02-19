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
  use mpi_f08, only: MPI_Allreduce, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, &
       mpi_min, mpi_max
  use comm, only: pe_rank, neko_comm, mpi_real_precision
  use, intrinsic :: iso_fortran_env, only: stderr => error_unit
  use utils, only: neko_error

  implicit none
  private

  type, public :: mma_t
     private
     integer :: n, m, max_iter
     real(kind=rp) :: a0, f0val, asyinit, asyincr, asydecr, epsimin, residumax, &
          residunorm
     type(vector_t) :: xold1, xold2, low, upp, alpha, beta, a, c, d, xmax, xmin
     logical :: is_initialized = .false.
     logical :: is_updated = .false.
     character(len=:), allocatable :: backnd

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
     procedure, public, pass(this) :: init_json => mma_init_json
     procedure, public, pass(this) :: free => mma_free
     procedure, public, pass(this) :: get_n => mma_get_n
     procedure, public, pass(this) :: get_m => mma_get_m
     procedure, public, pass(this) :: get_residumax => mma_get_residumax
     procedure, public, pass(this) :: get_residunorm => mma_get_residunorm
     procedure, public, pass(this) :: get_max_iter => mma_get_max_iter

     procedure, public, pass(this) :: init => mma_init_attributes
     procedure, public, pass(this) :: KKT => mma_KKT
     procedure, public, pass(this) :: update => mma_update

  end type mma_t

  interface
     ! ======================================================================== !
     !! interface for cpu backend module subroutines
     module subroutine mma_init_attributes_cpu(this, x, n, m, a0, a, c, d, &
          xmin, xmax, max_iter, epsimin, asyinit, asyincr, asydecr)
       class(mma_t), intent(inout) :: this
       integer, intent(in) :: n, m
       real(kind=rp), intent(in), dimension(n) :: x
       real(kind=rp), intent(in), dimension(n) :: xmax, xmin
       real(kind=rp), intent(in), dimension(m) :: a, c, d
       real(kind=rp), intent(in) :: a0
       integer, intent(in), optional :: max_iter
       real(kind=rp), intent(in), optional :: epsimin, asyinit, asyincr, &
            asydecr
     end subroutine mma_init_attributes_cpu

     module subroutine mma_update_cpu(this, iter, x, df0dx, fval, dfdx)
       class(mma_t), intent(inout) :: this
       integer, intent(in) :: iter
       real(kind=rp), dimension(this%n), intent(inout) :: x
       type(vector_t) :: df0dx, fval
       type(matrix_t) :: dfdx
     end subroutine mma_update_cpu

     module subroutine mma_KKT_cpu(this, x, df0dx, fval, dfdx)
       class(mma_t), intent(inout) :: this
       real(kind=rp), dimension(this%n), intent(in) :: x
       type(vector_t), intent(in) :: df0dx, fval
       type(matrix_t), intent(in) :: dfdx
     end subroutine mma_KKT_cpu


     !! interface for device backend module subroutines
     module subroutine mma_init_attributes_device(this, x, n, m, a0, a, c, d, &
          xmin, xmax, max_iter, epsimin, asyinit, asyincr, asydecr)
       class(mma_t), intent(inout) :: this
       integer, intent(in) :: n, m
       real(kind=rp), intent(in), dimension(n) :: x
       real(kind=rp), intent(in), dimension(n) :: xmax, xmin
       real(kind=rp), intent(in), dimension(m) :: a, c, d
       real(kind=rp), intent(in) :: a0
       integer, intent(in), optional :: max_iter
       real(kind=rp), intent(in), optional :: epsimin, asyinit, asyincr, asydecr
     end subroutine mma_init_attributes_device

     module subroutine mma_update_device(this, iter, x, df0dx, fval, dfdx)
       class(mma_t), intent(inout) :: this
       integer, intent(in) :: iter
       real(kind=rp), dimension(this%n), intent(inout) :: x
       type(vector_t) :: df0dx, fval
       type(matrix_t) :: dfdx
     end subroutine mma_update_device

     module subroutine mma_KKT_device(this, x, df0dx, fval, dfdx)
       class(mma_t), intent(inout) :: this
       real(kind=rp), dimension(this%n), intent(in) :: x
       type(vector_t), intent(in) :: df0dx, fval
       type(matrix_t), intent(in) :: dfdx
     end subroutine mma_KKT_device

  end interface

contains

  subroutine mma_init_json( this, x, n, m, json, scale, auto_scale)
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
    !class(mma_t), allocatable :: mma
    ! class(mma_t), allocatable, intent(inout) :: this
    class(mma_t), intent(inout) :: this
    integer, intent(in) :: n, m
    real(kind=rp), intent(in), dimension(n) :: x

    type(json_file), intent(inout) :: json
    ! -------------------------------------------------------------------!
    !      Internal parameters for MMA                                   !
    !      Minimize  f_0(x) + a_0*z + sum( c_i*y_i + 0.5*d_i*(y_i)^2 )   !
    !    subject to  f_i(x) - a_i*z - y_i <= 0,  i = 1,...,m             !
    !                xmin_j <= x_j <= xmax_j,    j = 1,...,n             !
    !                z >= 0,   y_i >= 0,         i = 1,...,m             !
    ! -------------------------------------------------------------------!
    real(kind=rp), dimension(n) :: xmax, xmin
    real(kind=rp), dimension(m) :: a, c, d

    !! For reading the values from json and then set the value for the arrays
    real(kind=rp) :: a0 , xmax_const, xmin_const, a_const, c_const, d_const

    integer :: max_iter, n_global, ierr
    real(kind=rp) :: epsimin, asyinit, asyincr, asydecr

    !! Read the scaling info for fval and dfdx from json
    real(kind=rp), intent(out) :: scale
    logical, intent(out) :: auto_scale

    call MPI_Allreduce(n, n_global, 1, MPI_INTEGER, &
         MPI_SUM, MPI_COMM_WORLD, ierr)

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

    call json_get_or_default(json, 'mma.backend', this%backnd, 'cpu')

    call json_get_or_default(json, 'mma.xmin', xmin_const, 0.0_rp)
    call json_get_or_default(json, 'mma.xmax', xmax_const, 1.0_rp)
    call json_get_or_default(json, 'mma.a0', a0, 1.0_rp)
    call json_get_or_default(json, 'mma.a', a_const, 0.0_rp)
    call json_get_or_default(json, 'mma.c', c_const, 100.0_rp)
    call json_get_or_default(json, 'mma.d', d_const, 0.0_rp)

    call json_get_or_default(json, 'mma.scale', scale, 10.0_rp)
    call json_get_or_default(json, 'mma.auto_scale', auto_scale, .true.)

    call json_get_or_default(json, 'mma.epsimin', epsimin, 1.0e-9_rp)
    ! Initialize the MMA object with the parsed parameters
    a = a_const
    c = c_const
    d = d_const
    xmin = xmin_const
    xmax = xmax_const
    ! initializing the mma concrete type (mma_cpu_t or mma_device_t)
    if (pe_rank .eq. 0) then
       print *, "Initializing MMA backend to >>> ", this%backnd
    end if

    ! ------------------------------------------------------------------------ !
    ! Initialize the MMA object with the parameters read from json
    ! call this%init(x, n, m, a0, a, c, d, xmin, xmax, &
    !      max_iter, epsimin, asyinit, asyincr, asydecr, backnd)
    call this%init(x, n, m, a0, a, c, d, xmin, xmax, &
         max_iter, epsimin, asyinit, asyincr, asydecr)

  end subroutine mma_init_json

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

  subroutine mma_init_attributes(this, x, n, m, a0, a, c, d, xmin, xmax, &
       max_iter, epsimin, asyinit, asyincr, asydecr)
    class(mma_t), intent(inout) :: this
    integer, intent(in) :: n, m
    real(kind=rp), intent(in), dimension(n) :: x

    real(kind=rp), intent(in), dimension(n) :: xmax, xmin
    real(kind=rp), intent(in), dimension(m) :: a, c, d
    real(kind=rp), intent(in) :: a0
    integer, intent(in), optional :: max_iter
    real(kind=rp), intent(in), optional :: epsimin, asyinit, asyincr, asydecr


    ! Select backend type
    select case (this%backnd)
    case ("cpu")
       call mma_init_attributes_cpu(this, x, n, m, a0, a, c, d, xmin, xmax, &
            max_iter, epsimin, asyinit, asyincr, asydecr)
       if (pe_rank == 0) then
          print *, "MMA initialized with CPU backend!"
       end if
    case ("cuda")
       call mma_init_attributes_device(this, x, n, m, a0, a, c, d, xmin, xmax, &
            max_iter, epsimin, asyinit, asyincr, asydecr)
       if (pe_rank == 0) then
          print *, "MMA initialized with CUDA backend!"
       end if
    case default
       call neko_error('Unknown backend in mma_init_attributes')
    end select

  end subroutine mma_init_attributes

  subroutine mma_update(this, iter, x, df0dx, fval, dfdx)
    class(mma_t), intent(inout) :: this
    integer, intent(in) :: iter
    real(kind=rp), dimension(this%n), intent(inout) :: x
    type(vector_t) :: df0dx, fval
    type(matrix_t) :: dfdx

    ! Select backend type
    select case (this%backnd )
    case ("cpu")
       call mma_update_cpu(this, iter, x, df0dx, fval, dfdx)
    case ("cuda")
       call mma_update_device(this, iter, x, df0dx, fval, dfdx)
    case default
       call mma_update_cpu(this, iter, x, df0dx, fval, dfdx)
    end select
  end subroutine mma_update

  subroutine mma_KKT(this, x, df0dx, fval, dfdx)
    class(mma_t), intent(inout) :: this
    real(kind=rp), dimension(this%n), intent(in) :: x
    type(vector_t), intent(in) :: df0dx, fval
    type(matrix_t), intent(in) :: dfdx

    ! Select backend type
    select case (this%backnd )
    case ("cpu")
       call mma_KKT_cpu(this, x, df0dx, fval, dfdx)
    case ("cuda")
       call mma_KKT_device(this, x, df0dx, fval, dfdx)
    case default
       call mma_KKT_cpu(this,x, df0dx, fval, dfdx)
    end select
  end subroutine mma_KKT

  ! ========================================================================== !
  ! Getters and setters
  pure function mma_get_n(this) result(n)
    class(mma_t), intent(in) :: this
    integer :: n
    n = this%n
  end function mma_get_n

  pure function mma_get_m(this) result(m)
    class(mma_t), intent(in) :: this
    integer :: m
    m = this%m
  end function mma_get_m

  pure function mma_get_residumax(this) result(residumax)
    class(mma_t), intent(in) :: this
    real(kind=rp) :: residumax
    residumax = this%residumax
  end function mma_get_residumax

  pure function mma_get_residunorm(this) result(residunorm)
    class(mma_t), intent(in) :: this
    real(kind=rp) :: residunorm
    residunorm = this%residunorm
  end function mma_get_residunorm

  pure function mma_get_max_iter(this) result(max_iter_value)
    class(mma_t), intent(in) :: this
    integer :: max_iter_value
    max_iter_value = this%max_iter
  end function mma_get_max_iter

end module mma

