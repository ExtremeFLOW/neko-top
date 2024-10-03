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
  use neko_config, only: NEKO_BCKND_DEVICE
  use vector, only: vector_t
  use matrix, only: matrix_t
  use json_module, only: json_file
  use json_utils, only: json_get_or_default

  ! Inclusions from external dependencies and standard libraries
  use, intrinsic :: iso_fortran_env, only: stderr => error_unit
  implicit none
  private

  type, public :: mma_t
     private

     real(kind=rp) :: a0, f0val, asyinit, asyincr, asydecr, epsimin, &
          residumax
     ! we need to have the KKT residual public
     ! or alternatively, have KKT a function not a subroutine, 
     ! that return residunrom
     real(kind=rp), public :: residunorm
     integer :: n, m, max_iter
     type(vector_t) :: xold1, xold2, low, upp, alpha, beta, a, c, d, xmax, xmin

     logical :: is_initialized = .false.
     logical :: is_updated = .false.
     character(len=:), allocatable :: backend

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
     procedure, public, pass(this) :: init => mma_init_attributes
     procedure, public, pass(this) :: init_json => mma_init_json

     procedure, public, pass(this) :: free => mma_free
     procedure, public, pass(this) :: KKT => mma_KKT

     !> Interface for updating the MMA
     generic, public :: update => mma_update_cpu, mma_update_vector
     ! for debugging I'm making this public
     procedure, public, pass(this) :: mma_update_cpu
     procedure, pass(this) :: mma_update_vector

     ! Getters for the MMA object
     procedure, public, pass(this) :: get_n => mma_get_n
     procedure, public, pass(this) :: get_m => mma_get_m
     procedure, public, pass(this) :: get_residumax => mma_get_residumax
     procedure, public, pass(this) :: get_residunorm => mma_get_residunorm

     !> Interface for generating the approximation sub problem
     generic :: gensub => mma_gensub_cpu
     procedure, pass(this) :: mma_gensub_cpu

     !> Interface for solving the dual with an interior point method
     generic :: subsolve => mma_subsolve_dpip_cpu
     procedure, pass(this) :: mma_subsolve_dpip_cpu

  end type mma_t

  ! ========================================================================== !
  ! Interface for the CPU backend

  interface
     !> Generate the approximation sub problem on the CPU.
     module subroutine mma_gensub_cpu(this, iter, x, df0dx, fval, dfdx)
       class(mma_t), intent(inout) :: this
       real(kind=rp), dimension(this%n), intent(in) :: x
       real(kind=rp), dimension(this%n), intent(in) :: df0dx
       real(kind=rp), dimension(this%m), intent(in) :: fval
       real(kind=rp), dimension(this%m, this%n), intent(in) :: dfdx
       integer, intent(in) :: iter
     end subroutine mma_gensub_cpu

     !> Solve the dual with an interior point method on the CPU.
     module subroutine mma_subsolve_dpip_cpu(this, designx)
       class(mma_t), intent(inout) :: this
       real(kind=rp), dimension(this%n), intent(inout) :: designx
     end subroutine mma_subsolve_dpip_cpu

     !> Compute the KKT condition for a given design x on the CPU.
     module subroutine mma_KKT_cpu(this, x, df0dx, fval, dfdx)
       class(mma_t), intent(inout) :: this
       real(kind=rp), dimension(this%n), intent(in) :: x
       real(kind=rp), dimension(this%m), intent(in) :: fval
       real(kind=rp), dimension(this%n), intent(in) :: df0dx
       real(kind=rp), dimension(this%m, this%n), intent(in) :: dfdx
     end subroutine mma_KKT_cpu
  end interface

contains

  subroutine mma_init_json(this, x, n, m, a0, a, c, d, xmin, xmax, json)
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
    real(kind=rp), intent(in), dimension(n) :: x
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
    type(json_file), intent(inout) :: json

    integer :: max_iter
    real(kind=rp) :: epsimin, asyinit, asyincr, asydecr
    character(len=:), allocatable :: backend

    ! ------------------------------------------------------------------------ !
    ! Assign defaults if nothing is parsed

    ! based on the Cpp Code by Niels
    call json_get_or_default(json, 'mma.epsimin', epsimin, &
         1.0e-9_rp * sqrt(real(m + n, rp)))
    call json_get_or_default(json, 'mma.max_iter', max_iter, 100)

    ! Following parameters are set based on eq.3.8:--------
    call json_get_or_default(json, 'mma.asyinit', asyinit, 0.5_rp)
    call json_get_or_default(json, 'mma.asyincr', asyincr, 1.2_rp)
    call json_get_or_default(json, 'mma.asydecr', asydecr, 0.7_rp)

    call json_get_or_default(json, 'mma.backend', backend, 'cpu')

    ! ------------------------------------------------------------------------ !
    ! Initialize the MMA object with the parsed parameters
    call this%init(x, n, m, a0, a, c, d, xmin, xmax, &
         max_iter, epsimin, asyinit, asyincr, asydecr, backend)

  end subroutine mma_init_json

  subroutine mma_init_attributes(this, x, n, m, a0, a, c, d, xmin, xmax, &
       max_iter, epsimin, asyinit, asyincr, asydecr, backend)
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
    real(kind=rp), intent(in), dimension(n) :: x
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
    character(len=:), allocatable, intent(in), optional :: backend

    call this%free()

    this%n = n
    this%m = m

    ! allocate(this%x(n))
    ! this%x = x
    call this%xold1%init(n)
    call this%xold2%init(n)
    this%xold1%x = x
    this%xold2%x = x

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
    call this%pij%init(m,n)
    call this%qij%init(m,n)
    call this%bi%init(m)

    !---nesessary for KKT check after updating df0dx, fval, dfdx --------
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

    !setting the bounds for the design variable based on the problem
    this%xmax%x = xmax
    this%xmin%x = xmin

    this%low%x(:) = minval(x)
    this%upp%x(:) = maxval(x)

    !setting KKT norms to a large number for the initial design
    this%residumax = huge(0.0_rp)
    this%residunorm = huge(0.0_rp)

    ! ------------------------------------------------------------------------ !
    ! Assign defaults if nothing is parsed

    ! based on the Cpp Code by Niels
    if (.not. present(epsimin)) this%epsimin = 1.0e-9_rp * sqrt(real(m + n, rp))
    if (.not. present(max_iter)) this%max_iter = 100

    ! Following parameters are set based on eq.3.8:--------
    if (.not. present(asyinit)) this%asyinit = 0.5_rp
    if (.not. present(asyincr)) this%asyincr = 1.2_rp
    if (.not. present(asydecr)) this%asydecr = 0.7_rp

    if (.not. present(backend)) this%backend = 'cpu'

    ! Assign values from inputs when present
    if (present(max_iter)) this%max_iter = max_iter
    if (present(epsimin)) this%epsimin = epsimin
    if (present(asyinit)) this%asyinit = asyinit
    if (present(asyincr)) this%asyincr = asyincr
    if (present(asydecr)) this%asydecr = asydecr
    if (present(backend)) this%backend = backend

    !the object is correctly initialized
    this%is_initialized = .true.
  end subroutine mma_init_attributes

  subroutine mma_update_cpu(this, iter, x, df0dx, fval, dfdx)
    ! ----------------------------------------------------- !
    ! Update the design variable x by solving the convex    !
    ! approximation of the problem.                         !
    !                                                       !
    ! This subroutine is called in each iteration of the    !
    ! optimization loop                                     !
    ! ----------------------------------------------------- !
    class(mma_t), intent(inout) :: this
    integer, intent(in) :: iter
    real(kind=rp), dimension(this%n), intent(inout) :: x
    real(kind=rp), dimension(this%n), intent(in) :: df0dx
    real(kind=rp), dimension(this%m), intent(in) :: fval
    real(kind=rp), dimension(this%m, this%n), intent(in) :: dfdx

    if (.not. this%is_initialized) then
       write(stderr, *) "The MMA object is not initialized."
       error stop
    end if

    ! generate a convex approximation of the problem
    call this%gensub(iter, x, df0dx, fval, dfdx)

    !solve the approximation problem using interior point method
    call this%subsolve(x)

    this%is_updated = .true.
  end subroutine mma_update_cpu

  subroutine mma_update_vector(this, iter, x, df0dx, fval, dfdx)
    ! ----------------------------------------------------- !
    ! Update the design variable x by solving the convex    !
    ! approximation of the problem.                         !
    !                                                       !
    ! This subroutine is called in each iteration of the    !
    ! optimization loop                                     !
    ! ----------------------------------------------------- !
    class(mma_t), intent(inout) :: this
    integer, intent(in) :: iter
    type(vector_t), intent(inout) :: x
    type(vector_t), intent(in) :: df0dx
    type(vector_t), intent(in) :: fval
    type(matrix_t), intent(in) :: dfdx

    if (.not. this%is_initialized) then
       write(stderr, *) "The MMA object is not initialized."
       error stop
    end if

    if (this%backend == 'cpu') then
       call this%mma_update_cpu(iter, x%x, df0dx%x, fval%x, dfdx%x)
    else
       write(stderr, *) "Device not supported for MMA."
       error stop
    end if

  end subroutine mma_update_vector

  subroutine mma_KKT(this, x, df0dx, fval, dfdx)
    ! ----------------------------------------------------- !
    ! Compute the KKT condition right hand side for a given !
    ! design x and set the max and norm values of the       !
    ! residue of KKT system to this%residumax and           !
    ! this%residunorm.                                      !
    !                                                       !
    ! The left hand sides of the KKT conditions are computed!
    ! for the following nonlinear programming problem:      !
    ! Minimize  f_0(x) + a_0*z +                            !
    !                       sum( c_i*y_i + 0.5*d_i*(y_i)^2 )!
    !   subject to  f_i(x) - a_i*z - y_i <= 0,  i = 1,...,m !
    !         xmax_j <= x_j <= xmin_j,    j = 1,...,n       !
    !        z >= 0,   y_i >= 0,         i = 1,...,m        !
    !                                                       !
    !                                                       !
    ! Note that before calling this function, the function  !
    ! values (f0val, fval, dfdx, ...) should be updated     !
    ! using the new x values.                               !
    ! ----------------------------------------------------- !
    class(mma_t), intent(inout) :: this
    real(kind=rp), dimension(this%n), intent(in) :: x
    real(kind=rp), dimension(this%m), intent(in) :: fval
    real(kind=rp), dimension(this%n), intent(in) :: df0dx
    real(kind=rp), dimension(this%m, this%n), intent(in) :: dfdx

    if (.not. this%is_initialized) then
       write(stderr, *) "The MMA object is not initialized."
       error stop
    end if

    if (NEKO_BCKND_DEVICE .eq. 0) then
       call mma_KKT_cpu(this, x, df0dx, fval, dfdx)
    else
       write(stderr, *) "Device not supported for MMA."
       error stop
    end if

  end subroutine mma_KKT

  !> Deallocate the MMA object.
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

end module mma

