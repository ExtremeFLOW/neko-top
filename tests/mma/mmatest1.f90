program mmatest
!===========================================================================!
!       Test_Case_I for the MMA solver implementation                       !
! Dependencies: Neko-top and Lapack (lapack function DGESV is used to solve !
! the linear system )                                                       !
! ------------------------------------------------------------------------- !
!                                                                           !
!        TEST_CASE_1:                                                       !
!    Problem 1 in "https://people.kth.se/~krille/originalmma.pdf"           !
!                                                                           !
!===========================================================================!


    use num_types, only: rp
    use mma
    implicit none

    !!!!!!! TEST_CASE_1
    integer, parameter :: n =5
    integer, parameter :: m =1 

    integer :: iter, i, j
    real(kind=rp), dimension(n) :: x
    real(kind=rp), dimension(m) :: fval
    real(kind=rp), dimension(m,n) :: dfdx
    real(kind=rp) :: f0val, tol
    real(kind=rp), dimension(n) :: df0dx
    integer :: maxiter
    logical :: mmadone
    type(mma_t) :: optprob
    ! ----------------------------------------------------------------------!
    !      Internal parameters for MMA                                      !
    !      Minimize  f_0(x) + a_0*z + sum( c_i*y_i + 0.5*d_i*(y_i)^2 )      !
    !    subject to  f_i(x) - a_i*z - y_i <= 0,  i = 1,...,m                !
    !                xmin_j <= x_j <= xmax_j,    j = 1,...,n                !
    !                z >= 0,   y_i >= 0,         i = 1,...,m                !
    ! ----------------------------------------------------------------------!
    real(kind=rp), dimension(n) ::  xmax, xmin
    real(kind=rp), dimension(m) ::a, c, d
    real(kind=rp) :: a0

    a0=1
    a=0
    c=1000
    d=0
    !setting the bounds for the design variable based on the problem
    xmax = 10
    xmin = 0

    !!!!!! intial value
    x = [5.0, 5.0, 5.0, 5.0, 5.0]
    !!!!!! opt value
    !!!!!! x = (/6.016, 5.309, 4.494, 3.502, 2.153/)

    !!!!!! Compute the objective function(f0val), constraints(fval) and 
    !!!!!! thier derivatives(df0dx, dfdx) for the initial value x
    call func1 (n, m, x, f0val , df0dx, fval , dfdx)
    print *, "f0val=", f0val, "fval=", fval
    print '(A,I3,A,5F10.4,A,F10.7,A,F10.7)', 'iter=', 0, ' , x= ', x, &
        '-------,f0val= ', f0val, ',   fval= ', fval

    !!!!!! Maximum number of iterations for the optimization and the tolerance 
    !!!!!! for the norm of KKT conditions for convergence
    maxiter = 40
    mmadone = .false.
    tol = 1.0e-8_rp

    !!!!!! Initializing the MMA object
    call optprob%init(x, n, m , a0, a, c, d, xmin, xmax)


    !!!!!! Optimization loop
    do iter = 1, maxiter
        
        call optprob%update(iter, x, df0dx, fval, dfdx)
        
        !Updating the design variables
        ! x= optprob%x(:)


        !update the function value and derivatives
        
        !!!!!!! TEST_CASE_1
        call func1 (n, m, x, f0val , df0dx, fval , dfdx)
        call optprob%KKT(x,df0dx,fval,dfdx)
        print '(A,I3,A,5F10.4,A,F10.7,A,F10.7,A,E10.4,A,E10.4)', &
            'iter=', iter, ' , x= ', x, &
            '-------,f0val= ', f0val, ',   fval= ', fval, &
            ',  KKTmax=', optprob%residumax, ', KKTnorm2=', optprob%residunorm

        if (optprob%residunorm .lt. tol) exit
    end do

    print *, "f0val=", f0val, "fval=", fval
contains

    subroutine func1 (n, m, x, f0val, df0dx, fval , dfdx)
        implicit none

        integer :: n
        integer :: m
        real(kind=rp), dimension(n), intent(in) :: x
        real(kind=rp), intent(inout) :: f0val
        real(kind=rp), dimension(n), intent(inout) :: df0dx
        real(kind=rp), dimension(m), intent(inout) :: fval
        real(kind=rp), dimension(m,n), intent(inout) :: dfdx
        ! ----------------------------------------------------------- !
        !  This file calculates function values and gradients         !
        !  for "toy problem 1":                                       !
        !                                                             !
        !    minimize 0.0624*(x(1) + x(2) + x(3) + x(4) + x(5))       !
        !  subject to 61/x(1)^3 + 37/x(2)^3 + 19/x(3)^3 + 7/x(4)^3 +  !
        !                    1/x(5)^3 - 1 =< 0                        ! 
        !              0 =< x(j) =< 10, for j=1,2,3,4,5.              !
        ! ----------------------------------------------------------- !

        f0val = 0.0624*sum(x)

        df0dx(1) = 0.0624
        df0dx(2) = 0.0624
        df0dx(3) = 0.0624
        df0dx(4) = 0.0624
        df0dx(5) = 0.0624

        fval(1) = 61/x(1)**3 + 37/x(2)**3 + 19/x(3)**3 + 7/x(4)**3+1/x(5)**3-1

        dfdx(1,1) = -3*61/x(1)**4
        dfdx(1,2) = -3*37/x(2)**4 
        dfdx(1,3) = -3*19/x(3)**4 
        dfdx(1,4) = -3*7/x(4)**4
        dfdx(1,5) = -3*1/x(5)**4
    end subroutine func1
end program mmatest
