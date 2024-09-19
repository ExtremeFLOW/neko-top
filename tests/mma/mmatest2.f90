program mmatest
!===========================================================================!
!       Test_Case_II for the MMA solver implementation                      !
! Dependencies: Neko-top and Lapack (lapack function DGESV is used to solve !
! the linear system )                                                       !
! ------------------------------------------------------------------------- !
!                                                                           !
!        TEST_CASE_2:                                                       !
!    Problem 1 in "https://people.kth.se/~krille/originalmma.pdf"           !
!                                                                           !
!===========================================================================!


    use num_types, only: rp
    use mma
    implicit none


    !!!!!!! TEST_CASE_2
    integer, parameter :: n =3
    integer, parameter :: m =2

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

    a0=1.0
    a=0
    c=1000
    d=1
    !setting the bounds for the design variable based on the problem
    xmax(:) = 5
    xmin(:) = 0

    !!!!!! intial value
    x = [4.0,3.0,2.0]
    !!!!!!! opt value
    !!!!!!! x = [2.017526, 1.779968, 1.237558]

    !!!!!! Compute the objective function(f0val), constraints(fval) and 
    !!!!!! thier derivatives(df0dx, dfdx) for the initial value x
    call func2 (n, m, x, f0val , df0dx, fval , dfdx)
    print *, "f0val=", f0val, "fval=", fval
    print '(A,I3,A,3F10.4,A,F10.7,A,F10.7,A,F10.7)', 'iter=', 0, ' , x= ', &
       x, '-------,f0val= ', f0val, ',   fval1= ', fval(1), ',   fval2= ', &
       fval(1)





    maxiter = 40
    mmadone = .false.
    tol = 1.0e-8_rp
    
    call optprob%init(x, n, m , a0, a, c, d, xmin, xmax)

    do iter = 1, maxiter
        
        call optprob%update(iter, x, df0dx, fval, dfdx)
        
        !Updating the design variables
        ! x= optprob%x(:)


        !update the function value and derivatives
        
        !!!!!! TEST_CASE_2
        call func2 (n, m, x, f0val , df0dx, fval , dfdx)
        call optprob%KKT(x,df0dx,fval,dfdx)
            print '(A,I3,A,3F10.4,A,F10.7,A,F10.7,A,F10.7,A,E10.4,A,E10.4)', &
                'iter=', iter, &
                ' , x= ', x, '-------,f0val= ', f0val, ',   fval1= ',  &
                fval(1), ',   fval2= ', fval(1), &
           ',  KKTmax=', optprob%residumax, ', KKTnorm2=', optprob%residunorm

        if (optprob%residunorm .lt. tol) exit
    end do

    print *, "f0val=", f0val, "fval=", fval


contains

    subroutine func2 (n, m, x, f0val, df0dx, fval , dfdx)
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
        !  for "toy problem 2":                                       !
        !                                                             !
        !    minimize x(1)^2 + x(2)^2 + x(3)^2                        !
        !  subject to (x(1)-5)^2 + (x(2)-2)^2 + (x(3)-1)^2 - 9 =< 0   !
        !             (x(1)-3)^2 + (x(2)-4)^2 + (x(3)-3)^2 - 9 =< 0   !
        !              0 =< x(j) =< 5, for j=1,2,3.                   !
        ! ----------------------------------------------------------- !

        f0val = sum(x**2)
            
        df0dx = 2*x

        !fval  = [(x(1)-5)^2+(x(2)-2)^2+(x(3)-1)^2-9
        !         (x(1)-3)^2+(x(2)-4)^2+(x(3)-3)^2-9];
        fval(1) = (x(1)-5)**2 + (x(2)-2)**2 + (x(3)-1)**2 - 9.0
        fval(2) = (x(1)-3)**2 + (x(2)-4)**2 + (x(3)-3)**2 - 9.0
        
        !dfdx  = 2*[x(1)-5  x(2)-2  x(3)-1
        !           x(1)-3  x(2)-4  x(3)-3];
        dfdx(1,:) = 2.0*[x(1) - 5.0,x(2) - 2.0,x(3) - 1.0]
        dfdx(2,:) = 2.0*[x(1) - 3.0,x(2) - 4.0,x(3) - 3.0]  

    end subroutine func2

end program mmatest
