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
  implicit none
  private

  !> Abstract type to compute pressure residual
  type, public, abstract :: mma_t
    integer :: n, m
   contains
    !> Interface for initializing the MMA object
    procedure, public, nopass :: init_json => mma_init_json

    ! Add interfaces to the abstract type procedure
    procedure(mma_init), pass(this), deferred :: init
    procedure(mma_update), pass(this), deferred :: update

    ! Deferred methods for accessing properties of concrete types of mma_t
    procedure(get_n), public, deferred :: get_n
    procedure(get_m), public, deferred :: get_m
    procedure(get_residumax), public, deferred :: get_residumax
    procedure(get_residunorm), public, deferred :: get_residunorm
    procedure(get_max_iter), public, deferred :: get_max_iter

  end type mma_t

  abstract interface
    subroutine mma_init(this, x, n, m, a0, a, c, d, xmin, xmax, &
       max_iter, epsimin, asyinit, asyincr, asydecr, backend)
      import mma_t, rp
      class(mma_t), intent(inout) :: this
      integer, intent(in) :: n, m
      real(kind=rp), intent(in), dimension(n) :: x
      real(kind=rp), intent(in), dimension(n) :: xmax, xmin
      real(kind=rp), intent(in), dimension(m) :: a, c, d
      real(kind=rp), intent(in) :: a0
      integer, intent(in), optional :: max_iter
      real(kind=rp), intent(in), optional :: epsimin, asyinit, asyincr, asydecr
      character(len=:), allocatable, intent(in), optional :: backend
    end subroutine mma_init

    subroutine mma_update(this, iter, x, df0dx, fval, dfdx)
      import mma_t, rp
      class(mma_t), intent(inout) :: this
      integer, intent(in) :: iter
      real(kind=rp), dimension(this%n), intent(inout) :: x
      real(kind=rp), dimension(this%n), intent(in) :: df0dx
      real(kind=rp), dimension(this%m), intent(in) :: fval
      real(kind=rp), dimension(this%m, this%n), intent(in) :: dfdx
    end subroutine mma_update
  end interface

  abstract interface
    pure function get_n(this) result(n)
      import mma_t
      class(mma_t), intent(in) :: this
      integer :: n
    end function get_n

    pure function get_m(this) result(m)
      import mma_t
      class(mma_t), intent(in) :: this
      integer :: m
    end function get_m

    pure function get_residumax(this) result(residumax)
      import mma_t, rp
      class(mma_t), intent(in) :: this
      real(kind=rp) :: residumax
    end function get_residumax

    pure function get_residunorm(this) result(residunorm)
      import mma_t, rp
      class(mma_t), intent(in) :: this
      real(kind=rp) :: residunorm
    end function get_residunorm

    pure function get_max_iter(this) result(max_iter_value)
      import mma_t
      class(mma_t), intent(in) :: this
      integer :: max_iter_value
    end function get_max_iter
  end interface

  interface

     !> Factory for the mma_t
     !! @details Only selects the compute backend.
     !! @param object The object to be allocated by the factory.
     module subroutine mma_factory(object)
       class(mma_t), allocatable, intent(inout) :: object
     end subroutine mma_factory
  end interface
  public :: mma_factory
  
  contains 

  subroutine mma_init_json( x, n, json, scale, auto_scale)
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
    class(mma_t), allocatable :: mma
    ! class(pnpn_prs_res_t), allocatable :: prs_res
    integer, intent(in) :: n
    integer, intent(in) :: m
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
    real(kind=rp), allocatable :: a(:), c(:), d(:)

    !! For reading the values from json and then set the value for the arrays
    real(kind=rp) :: a0 , xmax_const, xmin_const, a_const, c_const, d_const

    integer :: max_iter
    real(kind=rp) :: epsimin, asyinit, asyincr, asydecr
    character(len=:), allocatable :: backend

    !! Read the scaling info for fval and dfdx from json
    real(kind=rp), intent(out) :: scale
    logical, intent(out) :: auto_scale
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

    call json_get_or_default(json, 'mma.xmin', xmin_const, 0.0_rp)
    call json_get_or_default(json, 'mma.xmax', xmax_const, 1.0_rp)
    call json_get_or_default(json, 'mma.a0', a0, 1.0_rp)
    call json_get_or_default(json, 'mma.a', a_const, 0.0_rp)
    call json_get_or_default(json, 'mma.c', c_const, 100.0_rp)
    call json_get_or_default(json, 'mma.d', d_const, 0.0_rp)

    call json_get_or_default(json, 'mma.scale', scale, 10.0_rp)
    call json_get_or_default(json, 'mma.auto_scale', auto_scale, .true.)


    allocate(a(m))
    allocate(c(m))
    allocate(d(m))
    a = a_const
    c = c_const
    d = d_const
    xmin = xmin_const
    xmax = xmax_const


    ! calling the mma_factory(mma) to set the drived type for mma
    call mma_factory(mma)
    ! ------------------------------------------------------------------------ !
    ! Initialize the MMA object with the parsed parameters
    call mma%init(x, n, m, a0, a, c, d, xmin, xmax, &
         max_iter, epsimin, asyinit, asyincr, asydecr, backend)

  end subroutine mma_init_json


end module mma

