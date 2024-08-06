program mmatest
!===========================================================================!
!       This program is used to test the MMA solver implementation          !
! This solves two different test problems using the MMA implementation in   !
! the file mma.f90. Note that lapack function DGESV is used to solve linear !
! system in mma.f90. Also num_types.f90 is used for data type configuration.!
! Run by: (gfortran mmatest.f90 num_types.f90 mma.f90 -llapack)             !
! ------------------------------------------------------------------------- !
!                                                                           !
!        TEST_CASE_1:                                                       !
!    Problem 1 in "https://people.kth.se/~krille/originalmma.pdf"           !
!                                                                           !
!                                                                           !
!        TEST_CASE_2:                                                       !
!    Problem 1 in "https://people.kth.se/~krille/originalmma.pdf"           !
!                                                                           !
!                                                                           !
! ------------------------------------------------------------------------- !
!                                                                           !
! To switch between the above two test cases:                               !
! uncomment/comment the lines 36-37 and 40-41.                              !
! uncomment/comment the lines 59-66 and 76-83.                              !
! uncomment/comment the lines 106-111 and 115-120.                          !
! Also, in file mma.f90 in "subroutine mma_init(this, x, n, m)",            ! 
! uncomment/comment lines 166-172 and 175-181.                              !
!===========================================================================!


    use num_types
    use mma
    implicit none



    !!!!!!! TEST_CASE_1
    integer, parameter :: n =5
    integer, parameter :: m =1 

    !!!!!!! TEST_CASE_2
    ! integer, parameter :: n =3
    ! integer, parameter :: m =2

    integer :: iter, i, j
    real(kind=rp), dimension(n) :: x
    real(kind=rp), dimension(m) :: fval
    real(kind=rp), dimension(m,n) :: dfdx
    real(kind=rp) :: f0val
    real(kind=rp), dimension(n) :: df0dx
    integer :: maxiter
    logical :: mmadone
    type(mma_t) :: optprob

    !!!!!!! TEST_CASE_1
    !!!!!!! for test problem 1 in "https://people.kth.se/~krille/originalmma.pdf" 
    !!!!!!! change the values of n and m at the begining to 5 and 1, respectively.
    !!!!!!! integer, parameter :: n =5
    !!!!!!! integer, parameter :: m =1   
    ! ! x = [6.0, 5.0, 4.0, 3.0, 2.0]
    x = [5.0, 5.0, 5.0, 5.0, 5.0]
    !x =[6.0093,    5.3217 ,   4.5247  ,  3.4820  ,  2.1137]
    !!!!!! opt value
    !!!!!! x = (/6.016, 5.309, 4.494, 3.502, 2.153/)
    call func1 (n, m, x, f0val , df0dx, fval , dfdx)
    print *, "f0val=", f0val, "fval=", fval
    print '(A,I3,A,5F10.4,A,F10.7,A,F10.7)', 'iter=', 0, ' , x= ', x, &
        '-------,f0val= ', f0val, ',   fval= ', fval



    !!!!!!! TEST_CASE_2
    !!!!!!! For test problem 1 in:
    !!!!!!! "https://comsolyar.com/wp-content/uploads/2020/03/gcmma.pdf"
    !!!!!!! change the values of n and m at the begining to 3 and 2.
    !!!!!!! integer, parameter :: n =3
    !!!!!!! integer, parameter :: m =2  
    ! x = [4.0,3.0,2.0]
    ! !!!!!!! opt value
    ! !!!!!!! x = [2.017526, 1.779968, 1.237558]
    ! call func2 (n, m, x, f0val , df0dx, fval , dfdx)
    ! print *, "f0val=", f0val, "fval=", fval
    ! print '(A,I3,A,3F10.4,A,F10.7,A,F10.7,A,F10.7)', 'iter=', 0, ' , x= ', &
    !    x, '-------,f0val= ', f0val, ',   fval1= ', fval(1), ',   fval2= ', &
    !    fval(1)





    maxiter = 40
    mmadone = .false.
    

    call optprob%init(x,n,m)

    do iter = 1, maxiter
        
        call optprob%update(iter, x, df0dx, fval, dfdx)
        
        !Updating the design variables
        x= optprob%x(:)


        !update the function value and derivatives
        
        !!!!!!! TEST_CASE_1
        call func1 (n, m, x, f0val , df0dx, fval , dfdx)
        call optprob%KKT(df0dx,fval,dfdx)
        print '(A,I3,A,5F10.4,A,F10.7,A,F10.7,A,E10.4,A,E10.4)', &
            'iter=', iter, ' , x= ', x, &
            '-------,f0val= ', f0val, ',   fval= ', fval, &
            ',  KKTmax=', optprob%residumax, ', KKTnorm2=', optprob%residunorm


        !!!!!!! TEST_CASE_2
        ! call func2 (n, m, x, f0val , df0dx, fval , dfdx)
        ! call optprob%KKT(df0dx,fval,dfdx)
        !     print '(A,I3,A,3F10.4,A,F10.7,A,F10.7,A,F10.7,A,E10.4,A,E10.4)', 'iter=', iter, &
        !         ' , x= ', x, '-------,f0val= ', f0val, ',   fval1= ',  &
        !         fval(1), ',   fval2= ', fval(1), &
        !    ',  KKTmax=', optprob%residumax, ', KKTnorm2=', optprob%residunorm

        if (optprob%residunorm .lt. 1.0e-8) exit
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
