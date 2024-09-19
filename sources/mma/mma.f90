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
  use comm, only: neko_comm, mpi_real_precision, pe_rank, pe_size
  use neko_config, only: NEKO_BCKND_DEVICE

  ! Inclusions from external dependencies and standard libraries
  use, intrinsic :: iso_fortran_env, only: stderr => error_unit
  use mpi_f08, only: mpi_sum, MPI_Allreduce, mpi_max, mpi_min, mpi_sum, &
       mpi_integer
  implicit none
  private

  type, public :: mma_t
     private

     real(kind=rp) :: a0, f0val, asyinit, asyincr, asydecr, epsimin, &
          residumax, residunorm
     integer :: n, m, max_iter
     real(kind=rp), allocatable :: xold1(:), xold2(:), low(:), &
          upp(:), alpha(:), beta(:), a(:), c(:), d(:), xmax(:), xmin(:)

     logical :: is_initialized = .false.
     logical :: is_updated = .false.

     ! Internal dummy variables for MMA
     real(kind=rp), allocatable :: p0j(:), q0j(:)
     real(kind=rp), allocatable :: pij(:,:), qij(:,:)
     real(kind=rp), allocatable :: bi(:)

     !---nesessary for KKT check after updating df0dx, fval, dfdx --------
     real(kind=rp) :: z, zeta
     real(kind=rp), allocatable :: y(:), lambda(:), s(:), mu(:)
     real(kind=rp), allocatable :: xsi(:), eta(:)

   contains
     procedure, public, pass(this) :: init => mma_init
     procedure, public, pass(this) :: free => mma_free
     procedure, public, pass(this) :: update => mma_update
     procedure, public, pass(this) :: KKT => mma_KKT

     ! Getters for the MMA object
     procedure, public, pass(this) :: get_n => mma_get_n
     procedure, public, pass(this) :: get_m => mma_get_m
     procedure, public, pass(this) :: get_residumax => mma_get_residumax
     procedure, public, pass(this) :: get_residunorm => mma_get_residunorm

     !Generates the sub problem--the MMA convex approximation
     procedure, pass(this) :: mma_gensub
     !Solve the dual with an interior point method
     procedure, pass(this) :: mma_subsolve_dpip

  end type mma_t

  ! ========================================================================== !
  ! Interface for the CPU backend

  interface
     module subroutine mma_gensub_cpu(this, iter, x, df0dx, fval, dfdx)
       ! ----------------------------------------------------- !
       ! Generate the approximation sub problem by computing   !
       ! the lower and upper asymtotes and the other necessary !
       ! parameters (alpha, beta, p0j, q0j, pij, qij, ...).    !
       ! ----------------------------------------------------- !
       class(mma_t), intent(inout) :: this
       real(kind=rp), dimension(this%n), intent(in) :: x
       real(kind=rp), dimension(this%n), intent(in) :: df0dx(:)
       real(kind=rp), dimension(this%m), intent(in) :: fval(:)
       real(kind=rp), dimension(this%m, this%n), intent(in) :: dfdx(:,:)
       integer, intent(in) :: iter
     end subroutine mma_gensub_cpu
  end interface

contains

  subroutine mma_init(this, x, n, m, a0, a, c, d, xmin, xmax)
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
    real(kind=rp), dimension(n) :: xmax, xmin
    real(kind=rp), dimension(m) :: a, c, d
    real(kind=rp) :: a0

    this%n = n
    this%m = m

    call this%free()

    ! allocate(this%x(n))
    ! this%x = x
    allocate(this%xold1(n))
    allocate(this%xold2(n))
    this%xold1 = x
    this%xold2 = x

    allocate(this%alpha(n))
    allocate(this%beta(n))

    allocate(this%a(m))
    allocate(this%c(m))
    allocate(this%d(m))
    allocate(this%low(n))
    allocate(this%upp(n))
    allocate(this%xmax(n))
    allocate(this%xmin(n))

    !internal dummy variables for MMA
    allocate(this%p0j(n))
    allocate(this%q0j(n))
    allocate(this%pij(m,n))
    allocate(this%qij(m,n))
    allocate(this%bi(m))

    !---nesessary for KKT check after updating df0dx, fval, dfdx --------
    allocate(this%y(m))
    allocate(this%lambda(m))
    allocate(this%s(m))
    allocate(this%mu(m))
    allocate(this%xsi(n))
    allocate(this%eta(n))


    ! this%epsimin =  1.0e-10_rp
    ! based on the Cpp Code by Neils
    this%epsimin = 1.0e-9_rp * sqrt(1.0*m + 1.0*n)

    this%max_iter = 100

    this%a0 = a0
    this%a = a
    this%c = c
    this%d = d
    !setting the bounds for the design variable based on the problem
    this%xmax = xmax
    this%xmin = xmin




    this%low(:) = minval(x)
    this%upp(:) = maxval(x)

    !following parameters are set based on eq.3.8:--------
    this%asyinit = 0.5_rp !
    this%asyincr = 1.2_rp ! 1.1
    this%asydecr = 0.7_rp !0.65

    !setting KKT norms to a large number for the initial design
    this%residumax = 10**5_rp
    this%residunorm = 10**5_rp

    !the object is correctly initialized
    this%is_initialized = .true.
  end subroutine mma_init

  subroutine mma_update(this, iter, x, df0dx, fval, dfdx)
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
    real(kind=rp), dimension(this%n), intent(in) :: df0dx(:)
    real(kind=rp), dimension(this%m), intent(in) :: fval(:)
    real(kind=rp), dimension(this%m, this%n), intent(in) :: dfdx(:,:)

    if (.not. this%is_initialized) then
       write(stderr, *) "The MMA object is not initialized."
       write(stderr, *) "Please call MMA_obj.init() first and then ", &
            "call MMA_obj.update()."
       error stop
    end if

    ! generate a convex approximation of the problem
    call this%mma_gensub(iter, x, df0dx, fval, dfdx)
    this%xold2 = this%xold1
    this%xold1 = x

    !solve the approximation problem using interior point method
    call this%mma_subsolve_dpip(x)

    this%is_updated = .true.
  end subroutine mma_update

  subroutine mma_gensub(this, iter, x, df0dx, fval, dfdx)
    ! ----------------------------------------------------- !
    ! Generate the approximation sub problem by computing   !
    ! the lower and upper asymtotes and the other necessary !
    ! parameters (alpha, beta, p0j, q0j, pij, qij, ...).    !
    ! ----------------------------------------------------- !
    class(mma_t), intent(inout) :: this
    real(kind=rp), dimension(this%n), intent(in) :: x
    real(kind=rp), dimension(this%n), intent(in) :: df0dx(:)
    real(kind=rp), dimension(this%m), intent(in) :: fval(:)
    real(kind=rp), dimension(this%m, this%n), intent(in) :: dfdx(:,:)
    integer, intent(in) :: iter

    if (NEKO_BCKND_DEVICE .eq. 0) then
       call mma_gensub_cpu(this, iter, x, df0dx, fval, dfdx)
    else
       write(stderr, *) "Device not supported for MMA."
       error stop
    end if

  end subroutine mma_gensub

  subroutine mma_subsolve_dpip(this, designx)
    ! ----------------------------------------------------- !
    ! Dual-primal interior point method using Newton's step !
    ! to solve MMA sub problem.                             !
    ! A Backtracking Line Search approach is used to compute!
    ! the step size; starting with the full Newton's step   !
    ! (delta = 1) and deviding by 2 until we have a step size !
    ! that leads to a feasible point while ensuring a       !
    ! decrease in the residue.                              !
    ! ----------------------------------------------------- !
    implicit none
    class(mma_t), intent(inout) :: this
    real(kind=rp), dimension(this%n), intent(inout) :: designx
    !Note that there is a local dummy "x" in this subroutine, thus, we call
    !the current design "designx" instead of just "x"
    integer :: i, j, k, iter, ggdumiter, itto, ierr
    real(kind=rp) :: epsi, residumax, residunorm, &
         z, zeta, rez, rezeta, &
         delz, dz, dzeta, &
         steg, dummy_one, zold, zetaold, newresidu
    real(kind=rp), dimension(this%m) :: y, lambda, s, mu, &
         rey, relambda, remu, res, &
         dely, dellambda, dummyGDinv, &
         dy, dlambda, ds, dmu, &
         yold, lambdaold, sold, muold, &
         globaltmp_m
    real(kind=rp), dimension(this%n) :: x, xsi, eta, &
         rex, rexsi, reeta, &
         delx, diagx, dx, dxsi, deta, &
         xold, xsiold, etaold
    real(kind=rp), dimension(3*this%n+4*this%m+2) :: residu
    real(kind=rp), dimension(4*this%m+2) :: residu_small
    real(kind=rp), dimension(2*this%n+4*this%m+2) :: xx, dxx

    real(kind=rp), dimension(this%m, this%n) :: GG
    real(kind=rp), dimension(this%m+1) :: bb
    real(kind=rp), dimension(this%m+1, this%m+1) :: AA
    real(kind=rp), dimension(this%m, this%m) :: globaltmp_mm
    ! using DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO ) in lapack to solve
    ! the linear system which needs the following parameters
    integer :: info
    integer, dimension(this%m+1) :: ipiv

    real(kind=rp) :: re_xstuff_squ_global

    integer :: nglobal

    !intial value for the parameters in the subsolve based on
    !page 15 of "https://people.kth.se/~krille/mmagcmma.pdf"
    dummy_one = 1
    epsi = 1 !100
    x(:) = 0.5*(this%alpha(:)+this%beta(:))
    y(:) = 1
    z = 1
    zeta = 1
    lambda(:) = 1
    s(:) = 1
    xsi(:) = max(1.0, 1.0/(x(:) - this%alpha(:)))
    eta(:) = max(1.0, 1.0/(this%beta(:) - x(:)))
    mu(:) = max(1.0, 0.5*this%c(:))

    do while (epsi .gt. 0.9*this%epsimin)
       ! calculating residuals based on
       ! "https://people.kth.se/~krille/mmagcmma.pdf" for the variables
       ! x, y, z, lambda residuals based on eq(5.9a)-(5.9d), respectively.
       rex(:) = ((this%p0j(:) + matmul(transpose(this%pij(:,:)), &
            lambda(:)))/(this%upp(:) - x(:))**2 - &
            (this%q0j(:) + matmul(transpose(this%qij(:,:)), &
            lambda(:)))/(x(:) - this%low(:))**2 ) - &
            xsi(:) + eta(:)

       call MPI_Allreduce(this%n, nglobal, 1, &
            MPI_INTEGER, mpi_sum, neko_comm, ierr)

       !!!! computing without matmul and transpose
       ! rex = 0.0_rp
       ! do j = 1, this%n
       !     do i = 1, this%m
       !         rex(j) = rex(j) + this%pij(i,j) * &
       !             lambda(i)/(this%upp(j) - x(j))**2 - &
       !             this%qij(i,j) * lambda(i)/(x(j) - this%low(j))**2
       !     end do
       !     rex(j) = rex(j) + this%p0j(j)/(this%upp(j) - x(j))**2 &
       !                     - this%q0j(j)/(x(j) - this%low(j))**2 &
       !                     - xsi(j)  + eta(j)
       ! end do


       rey(:) = this%c(:) + this%d(:)*y(:) - lambda(:) - mu(:)
       rez = this%a0 - zeta - dot_product(lambda(:), this%a(:))

       ! relambda(:) = matmul(this%pij(:,:),1.0/(this%upp(:) - x(:))) + &
       !         matmul(this%qij(:,:), 1.0/(x(:) - this%low(:))) - &
       !         this%a(:)*z - y(:) + s(:) - this%bi(:)
       relambda = 0.0_rp
       do i = 1, this%m
          do j = 1, this%n !this n is global
             ! Accumulate sums for relambda (the term gi(x))
             relambda(i) = relambda(i) + &
                  this%pij(i,j)/(this%upp(j) - x(j)) &
                  + this%qij(i,j)/(x(j) - this%low(j))
          end do
       end do


       globaltmp_m = 0.0_rp
       call MPI_Allreduce(relambda, globaltmp_m, this%m, &
            mpi_real_precision, mpi_sum, neko_comm, ierr)
       relambda = globaltmp_m - this%a(:)*z - y(:) + s(:) - this%bi(:)


       rexsi(:) = xsi(:)*(x(:) - this%alpha(:)) - epsi
       reeta(:) = eta(:)*(this%beta(:) - x(:)) - epsi
       remu(:) = mu(:)*y(:) - epsi
       rezeta = zeta*z - epsi
       res(:) = lambda(:)*s(:) - epsi

       residu = [rex, rey, rez, relambda, rexsi, reeta, remu, rezeta, res]
       residumax = 0_rp

       call MPI_Allreduce(maxval(abs(residu)), residumax, 1, &
            mpi_real_precision, mpi_max, neko_comm, ierr)

       re_xstuff_squ_global = 0_rp
       call MPI_Allreduce(norm2(rex)**2+norm2(rexsi)**2+norm2(reeta)**2,&
            re_xstuff_squ_global, 1, mpi_real_precision, mpi_sum,&
            neko_comm, ierr)
       residu_small = [rey, rez, relambda, &
            remu, rezeta, res]
       residunorm = sqrt(norm2(residu_small)**2 + re_xstuff_squ_global)


       do iter = 1, this%max_iter !ittt
          if (iter .gt. (this%max_iter -2)) then
             ! print *, "The mma inner loop seems not to converge"
             ! print *, "residumax = ", residumax, "for epsi = ", epsi, &
             !         ", ittt  = ", iter, "out of ", this%max_iter
          end if
          !Check the condition
          if (residumax .lt. epsi) exit

          delx = 0_rp
          do j = 1, this%n
             do i = 1, this%m
                delx(j) = delx(j) + this%pij(i,j) * &
                     lambda(i)/(this%upp(j) - x(j))**2 &
                     - this%qij(i,j) * lambda(i)/(x(j) - this%low(j))**2
             end do
             delx(j) = delx(j) + this%p0j(j)/(this%upp(j) - x(j))**2 &
                  - this%q0j(j)/(x(j) - this%low(j))**2 &
                  - epsi/(x(j) - this%alpha(j)) &
                  + epsi/(this%beta(j) - x(j))
          end do
          dely = this%c + this%d*y - lambda - epsi/y
          delz = this%a0 - dot_product(lambda(:), this%a(:)) - epsi/z

          dellambda(:) = 0.0_rp
          do i = 1, this%m
             do j = 1, this%n !this n is global
                ! Accumulate sums for dellambda (the term gi(x))
                dellambda(i) = dellambda(i) + &
                     this%pij(i,j)/(this%upp(j) - x(j)) &
                     + this%qij(i,j)/(x(j) - this%low(j))
             end do
          end do

          globaltmp_m = 0.0_rp
          call MPI_Allreduce(dellambda, globaltmp_m, this%m, &
               mpi_real_precision, mpi_sum, neko_comm, ierr)

          dellambda = globaltmp_m - this%a*z - y - this%bi + epsi/lambda

          ! delx(:) = ((this%p0j(:) + matmul(transpose(this%pij(:,:)), &
          !     lambda(:)))/(this%upp(:) - x(:))**2 - &
          !     (this%q0j(:) + matmul(transpose(this%qij(:,:)), &
          !     lambda(:)))/(x(:) - this%low(:))**2 ) - &
          !     epsi/(x(:) - this%alpha(:)) + epsi/(this%beta(:) - x(:))

          ! dely(:) =  this%c(:) + this%d(:)*y(:) - lambda(:) - epsi/y(:)
          ! delz = this%a0 - dot_product(lambda(:), this%a(:)) - epsi/z
          ! dellambda(:) = matmul(this%pij(:,:),1.0/(this%upp(:) - x(:)))+&
          !     matmul(this%qij(:,:), 1.0/(x(:) - this%low(:))) - &
          !     this%a(:)*z - y(:) - this%bi(:) + epsi/lambda(:)

          do ggdumiter = 1, this%m
             GG(ggdumiter, :) = this%pij(ggdumiter,:)/ &
                  (this%upp(:) - x(:))**2 - &
                  this%qij(ggdumiter,:)/(x(:) - this%low(:))**2
          end do

          diagx(:) = ((this%p0j(:) + matmul(transpose(this%pij(:,:)), &
               lambda(:)))/(this%upp(:) - x(:))**3 + &
               (this%q0j(:) + matmul(transpose(this%qij(:,:)), &
               lambda(:)))/(x(:) - this%low(:))**3 )
          diagx(:) = 2*diagx(:) + xsi(:)/(x(:) - this%alpha(:)) + &
               eta(:)/(this%beta(:)- x(:))


          !Here we only consider the case m<n in the matlab code
          !assembling the right hand side matrix based on eq(5.20)
          ! bb = [dellambda + dely(:)/(this%d(:) + &
          !         (mu(:)/y(:))) - matmul(GG,delx/diagx), delz ]
          !!!!!!!!!!!!!!for MPI computation of bb!!!!!!!!!!!!!!!!!!!!!!!!!
          bb = 0.0_rp
          do i = 1, this%m
             do j = 1, this%n ! this n is global
                bb(i) = bb(i) + GG(i, j) * (delx(j) / diagx(j))
             end do
          end do
          globaltmp_m = 0.0_rp
          call MPI_Allreduce(bb(1:this%m), globaltmp_m, this%m, &
               mpi_real_precision, mpi_sum, neko_comm, ierr)
          bb(1:this%m) = globaltmp_m

          bb(1:this%m) = dellambda + dely/(this%d + (mu/y)) - bb(1:this%m)
          bb(this%m +1) = delz
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! !assembling the coefficients matrix AA based on eq(5.20)
          ! AA(1:this%m,1:this%m) =  &
          ! matmul(matmul(GG,mma_diag(1/diagx)), transpose(GG))
          ! !update diag(AA)
          ! AA(1:this%m,1:this%m) = AA(1:this%m,1:this%m) + &
          !     mma_diag(s(:)/lambda(:) + 1.0/(this%d(:) + (mu(:)/y(:))))

          AA = 0.0_rp
          !Direct computation of the matrix multiplication
          !(for better performance)
          do i = 1, this%m
             do j = 1, this%m
                ! Compute the (i, j) element of AA
                do k = 1, this%n !this n is global
                   AA(i, j) = AA(i, j) + GG(i, k) * &
                        (1.0_rp / diagx(k)) * GG(j, k)
                end do
             end do
          end do

          globaltmp_mm = 0.0_rp
          call MPI_Allreduce(AA(1:this%m, 1:this%m), globaltmp_mm, &
               this%m*this%m, mpi_real_precision, mpi_sum, neko_comm, ierr)
          AA(1:this%m,1:this%m) = globaltmp_mm
          do i = 1, this%m
             !update the diag AA
             AA(i, i) = AA(i, i) + (s(i) / lambda(i) + &
                  1.0_rp / (this%d(i) + mu(i) / y(i)))
          end do

          AA(1:this%m, this%m+1) = this%a(:)
          AA(this%m+1, 1:this%m) = this%a(:)
          AA(this%m+1, this%m+1) = -zeta/z




          call DGESV(this%m+1, 1, AA, this%m+1, ipiv, bb, this%m+1, info)
          ! if info! = 0 then DGESV is failed.
          if (info .ne. 0) then
             write(stderr, *) "DGESV failed to solve the linear system in MMA."
             write(stderr, *) "Please check mma_subsolve_dpip in mma.f90"
             error stop
          end if

          dlambda = bb(1:this%m)
          dz = bb(this%m + 1)
          ! based on eq(5.19)
          dx = -delx/diagx - matmul(transpose(GG), dlambda)/diagx

          dy = (-dely+dlambda)/(this%d(:) + (mu(:)/y(:)))
          dxsi = -xsi + (epsi-dx*xsi(:))/(x(:) - this%alpha(:))
          deta = -eta + (epsi+dx*eta(:))/(this%beta(:) - x(:))
          dmu = -mu + (epsi-mu*dy(:))/y(:)
          dzeta = -zeta + (epsi-zeta*dz)/z
          ds = -s + (epsi-dlambda*s(:))/lambda(:)

          !2*this%n+4*this%m+2
          dxx = [dy, dz, dlambda, dxsi, deta, dmu, dzeta, ds]
          xx = [y, z, lambda, xsi, eta, mu, zeta, s]
          steg = maxval([dummy_one, -1.01*dxx/xx, -1.01*dx/ &
               (x(:) - this%alpha(:)), 1.01*dx/(this%beta(:) - x(:))])
          steg = 1.0/steg

          call MPI_Allreduce(steg, steg, 1, &
               mpi_real_precision, mpi_min, neko_comm, ierr)

          xold = x
          yold = y
          zold = z
          lambdaold = lambda
          xsiold = xsi
          etaold = eta
          muold = mu
          zetaold = zeta
          sold = s

          !The innermost loop to determine the suitable step length
          !using the Backtracking Line Search approach
          newresidu = 2*residunorm
          itto = 0
          do while ((newresidu .gt. residunorm) .and. (itto .lt. 50))
             itto = itto + 1
             !update the variables
             x = xold + steg*dx
             y = yold + steg*dy
             z = zold + steg*dz
             lambda = lambdaold + steg*dlambda
             xsi = xsiold + steg*dxsi
             eta = etaold + steg*deta
             mu = muold + steg*dmu
             zeta = zetaold + steg*dzeta
             s = sold + steg*ds
             !recompute the newresidu to see if this stepsize improves
             !the residue
             rex(:) = ((this%p0j(:) + matmul(transpose(this%pij(:,:)), &
                  lambda(:)))/(this%upp(:) - x(:))**2 - &
                  (this%q0j(:) + matmul(transpose(this%qij(:,:)), &
                  lambda(:)))/(x(:) - this%low(:))**2 ) - &
                  xsi(:) + eta(:)
             rey(:) = this%c(:) + this%d(:)*y(:) - lambda(:) - mu(:)
             rez = this%a0 - zeta - dot_product(lambda(:), this%a(:))
             ! relambda(:) = matmul(this%pij(:,:),1.0/&
             !         (this%upp(:) - x(:))) + matmul(this%qij(:,:), &
             !         1.0/(x(:) - this%low(:))) - this%a(:)*z - &
             !         y(:) + s(:) - this%bi(:)
             relambda = 0.0_rp
             do i = 1, this%m
                do j = 1, this%n !this n is global
                   ! Accumulate sums for relambda (the term gi(x))
                   relambda(i) = relambda(i) + &
                        this%pij(i,j)/(this%upp(j) - x(j)) &
                        + this%qij(i,j)/(x(j) - this%low(j))
                end do
             end do
             globaltmp_m = 0.0_rp
             call MPI_Allreduce(relambda, globaltmp_m, this%m, &
                  mpi_real_precision, mpi_sum, neko_comm, ierr)
             relambda = globaltmp_m


             relambda = relambda - this%a(:)*z - y(:) + s(:) - this%bi(:)

             rexsi(:) = xsi(:)*(x(:) - this%alpha(:)) - epsi
             reeta(:) = eta(:)*(this%beta(:) - x(:)) - epsi
             remu(:) = mu(:)*y(:) - epsi
             rezeta = zeta*z - epsi
             res(:) = lambda(:)*s(:) - epsi

             residu = [rex, rey, rez, relambda, &
                  rexsi, reeta, remu, rezeta, res]

             re_xstuff_squ_global = 0_rp
             call MPI_Allreduce(norm2(rex)**2 + &
                  norm2(rexsi)**2+norm2(reeta)**2, re_xstuff_squ_global, &
                  1, mpi_real_precision, mpi_sum, neko_comm, ierr)
             residu_small = [rey, rez, relambda, &
                  remu, rezeta, res]
             newresidu = sqrt(norm2(residu_small)**2 + &
                  re_xstuff_squ_global)

             steg = steg/2
          end do

          residunorm = newresidu
          residumax = 0_rp
          call MPI_Allreduce(maxval(abs(residu)), residumax, 1, &
               mpi_real_precision, mpi_max, neko_comm, ierr)

          !correct the step size for the extra devision by 2 in the final
          !loop
          steg = 2*steg

          ! print *,"Processor ",pe_rank, "iter = ", iter, "epsi = ", epsi, &
          !     "steg = ", steg, "residunorm = ",residunorm, &
          !       "residumax = ",residumax
       end do
       epsi = 0.1*epsi

    end do

    designx = x
    !update the parameters of the MMA object nesessary to compute KKT residu
    this%y = y
    this%z = z
    this%lambda = lambda
    this%zeta = zeta
    this%xsi = xsi
    this%eta = eta
    this%mu = mu
    this%s = s

  end subroutine mma_subsolve_dpip

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

    real(kind=rp), dimension(this%m), intent(in) :: fval(:)
    real(kind=rp), dimension(this%n), intent(in) :: df0dx(:)
    real(kind=rp), dimension(this%m, this%n), intent(in) :: dfdx(:,:)

    real(kind=rp) :: residumax, residunorm, rez, rezeta
    real(kind=rp), dimension(this%m) :: rey, relambda, remu, res
    real(kind=rp), dimension(this%n) :: rex, rexsi, reeta
    real(kind=rp), dimension(3*this%n+4*this%m+2) :: residu

    real(kind=rp), dimension(4*this%m+2) :: residu_small
    integer :: ierr, pe_rank
    real(kind=rp) :: re_xstuff_squ_global

    if (.not.(this%is_initialized .and. this%is_updated)) then
       write(stderr, *) &
            'The MMA object is either not initialized or not updated.'
       write(stderr, *) &
            'call mmaobj%init and mmaobj%update and then mmaobj%KKT.'
       error stop
    end if

    rex(:) = df0dx + matmul(transpose(dfdx), this%lambda) - this%xsi(:) + &
         this%eta(:)
    rey(:) = this%c(:) + this%d(:)*this%y(:) - this%lambda(:) - this%mu(:)
    rez = this%a0 - this%zeta - dot_product(this%lambda(:), this%a(:))

    relambda(:) = fval - this%a(:)*this%z - this%y(:) + this%s(:)
    rexsi(:) = this%xsi(:)*(x(:) - this%xmin(:))
    reeta(:) = this%eta(:)*(this%xmax(:) - x(:))
    remu(:) = this%mu(:)*this%y(:)
    rezeta = this%zeta*this%z
    res(:) = this%lambda(:)*this%s(:)


    residu = [rex, rey, rez, relambda, rexsi, reeta, remu, rezeta, res]

    this%residumax = 0_rp
    call MPI_Allreduce(maxval(abs(residu)), this%residumax, 1, &
         mpi_real_precision, mpi_max, neko_comm, ierr)



    re_xstuff_squ_global = 0_rp
    call MPI_Allreduce(norm2(rex)**2 + norm2(rexsi)**2 + norm2(reeta)**2, &
         re_xstuff_squ_global, 1, mpi_real_precision, mpi_sum, neko_comm, ierr)
    residu_small = [rey, rez, relambda, &
         remu, rezeta, res]
    residunorm = sqrt(norm2(residu_small)**2 + re_xstuff_squ_global)


    this%residunorm = residunorm

  end subroutine mma_KKT

  !> Deallocate the MMA object.
  subroutine mma_free(this)

    class(mma_t), intent(inout) :: this

    if (allocated(this%xold1)) deallocate(this%xold1)
    if (allocated(this%xold2)) deallocate(this%xold2)
    if (allocated(this%low)) deallocate(this%low)
    if (allocated(this%upp)) deallocate(this%upp)
    if (allocated(this%alpha)) deallocate(this%alpha)
    if (allocated(this%beta)) deallocate(this%beta)
    if (allocated(this%a)) deallocate(this%a)
    if (allocated(this%c)) deallocate(this%c)
    if (allocated(this%d)) deallocate(this%d)
    if (allocated(this%xmax)) deallocate(this%xmax)
    if (allocated(this%xmin)) deallocate(this%xmin)
    if (allocated(this%p0j)) deallocate(this%p0j)
    if (allocated(this%q0j)) deallocate(this%q0j)
    if (allocated(this%pij)) deallocate(this%pij)
    if (allocated(this%qij)) deallocate(this%qij)
    if (allocated(this%bi)) deallocate(this%bi)
    if (allocated(this%y)) deallocate(this%y)
    if (allocated(this%lambda)) deallocate(this%lambda)
    if (allocated(this%s)) deallocate(this%s)
    if (allocated(this%mu)) deallocate(this%mu)
    if (allocated(this%xsi)) deallocate(this%xsi)
    if (allocated(this%eta)) deallocate(this%eta)

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

