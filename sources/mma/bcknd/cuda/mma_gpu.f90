submodule (mma) mma_gpu
use mpi_f08, only: MPI_INTEGER, MPI_REAL, mpi_sum, mpi_min, mpi_max, &
MPI_Allreduce
use utils, only: neko_error
use comm, only: neko_comm, mpi_real_precision

contains
module subroutine mma_gensub_gpu(this, iter, x_d, df0dx_d, fval_d, dfdx_d)
    ! ----------------------------------------------------- !
    ! Generate the approximation sub problem by computing   !
    ! the lower and upper asymtotes and the other necessary !
    ! parameters (alpha, beta, p0j, q0j, pij, qij, ...).    !
    ! ----------------------------------------------------- !
    class(mma_t), intent(inout) :: this
    type(vector_t), intent(in) :: x_d
    type(vector_t), intent(in) :: df0dx_d
    type(vector_t), intent(in) :: fval_d
    type(matrix_t), intent(in) :: dfdx_d

    integer, intent(in) :: iter
    integer :: i, j, ierr
    real(kind=rp), dimension(this%m) :: globaltmp_m
    if (iter .lt. 3) then
      call mma_gensub1_gpu(this%low%x_d, this%upp%x_d,x%x_d, this%xmin%x_d, this%xmax%x_d, this%asyinit, this%n)
  else
      call mma_gensub2_gpu(this%low%x_d, this%upp%x_d, x%x_d, this%xold1%x_d, this%xold2%x_d,this%xmin%x_d, this%xmax%x_d, &
        this%asydecr, this%asyincr, this%n)
  end if
  call mma_gensub3_gpu(x%x_d, df0dx%x_d, dfdx%x_d,this%low%x_d, this%upp%x_d, this%xmin%x_d, this%xmax%x_d,this%alpha%x_d, &
    this%beta%x_d, this%p0j%x_d, this%q0j%x_d, this%pij%x_d, this%qij%x_d, this%n, this%m) 
  call mma_gensub4_gpu(x%x_d, this%low%x_d, this%upp%x_d, this%pij%x_d, this%qij%x_d, this%n, this%m, this%bi%x_d)
  call device_memcpy(this%bi%x, this%bi%x_d, this%m, DEVICE_TO_HOST, sync=.false.)
  !call mma_gensub_gpu( x%x_d,  this%xold1%x_d, this%xold2%x_d, df0dx%x_d, dfdx%x_d, this%xlow%x_d, this%xupp%x_d, this%xmin%x_d, this%xmax%x_d, &
  !    this%alpha%x_d, this%beta%x_d, this%p0j%x_d, this%q0j%x_d, this%pij%x_d, this%qij%x_d, this%bi%x_d, this%asyinit, this%asydecr, this%asyincr,&
  !    this%n, this%m, iter) 
  !call device_memcpy(this%bi%x, this%bi%x_d, this%m, DEVICE_TO_HOST, sync=.false.)
  globaltmp_m = 0.0_rp
  call MPI_Allreduce(this%bi%x, globaltmp_m, this%m, &
    mpi_real_precision, mpi_sum, neko_comm, ierr)
  call device_memcpy(fval%x, fval%x_d, this%m, DEVICE_TO_HOST, sync=.false.)
  this%bi%x = globaltmp_m - fval%x
end subroutine mma_gensub_cpu

subroutine mma_subsolve_dpip_cpu(this, designx)
    ! ------------------------------------------------------- !
    ! Dual-primal interior point method using Newton's step   !
    ! to solve MMA sub problem.                               !
    ! A Backtracking Line Search approach is used to compute  !
    ! the step size; starting with the full Newton's step     !
    ! (delta = 1) and deviding by 2 until we have a step size !
    ! that leads to a feasible point while ensuring a         !
    ! decrease in the residue.                                !
    ! ------------------------------------------------------- !
    class(mma_t), intent(inout) :: this
    type(vector_t), intent(in) :: designx
    !Note that there is a local dummy "x" in this subroutine, thus, we call
    !the current design "designx" instead of just "x"
    integer :: i, j, k, iter, ggdumiter, itto, ierr
    real(kind=rp) :: epsi, residumax, residunorm, &
    z, zeta, rez, rezeta, &
    delz, dz, dzeta, &
    steg, dummy_one, zold, zetaold, newresidu
    type(vector_t), intent(in) ::  y, lambda, s, mu, &   !!!m
    rey, relambda, remu, res, &
    dely, dellambda, &
    dy, dlambda, ds, dmu, &
    yold, lambdaold, sold, muold, &
    globaltmp_m
    type(vector_t), intent(in) :: x, xsi, eta, & !!!!!n
    rex, rexsi, reeta, &
    delx, diagx, dx, dxsi, deta, &
    xold, xsiold, etaold
    type(vector_t), intent(in) ::residu !!!dimension(3*this%n+4*this%m+2)
    type(vector_t), intent(in) ::residu_small !!!!dimension(4*this%m+2) 
    type(vector_t), intent(in) :: xx, dxx !!!dimension(2*this%n+4*this%m+2)

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

    ! intial value for the parameters in the subsolve based on
    ! page 15 of "https://people.kth.se/~krille/mmagcmma.pdf"
    dummy_one = 1
    epsi = 1 !100
    real(kind=rp) :: cons
    cons=0.5
    call cuda_add3s2(x%x_d,this%alpha%x_d,this%beta%x_d,cons,cons,this%n)
    cons=1
    call cuda_cfill(y%x_d,cons,this%m)
    z = 1
    zeta = 1
    call cuda_cfill(lambda%x_d,cons,this%m)
    call cuda_cfill(s%x_d,cons,this%m)

    call cuda_max(x%x_d,this%alpha%x_d, this%beta%x_d, this%xsi%x_d,this%eta%x_d, this%mu%x_d, this%c%x_d, this%n%x_d) 

    do while (epsi .gt. 0.9*this%epsimin)
       ! calculating residuals based on
       ! "https://people.kth.se/~krille/mmagcmma.pdf" for the variables
       ! x, y, z, lambda residuals based on eq(5.9a)-(5.9d), respectively.
       call cuda_rex(rex%x_d,  x%x_d,  this%low%x_d, this%upp%x_d,  this%pij%x_d, this%p0j%x_d,this%qij%x_d, this%q0j%x_d, &
        lambda%x_d, xsi%x_d, eta%x_d, this%n, this%m) 

       call MPI_Allreduce(this%n, nglobal, 1, &
        MPI_INTEGER, mpi_sum, neko_comm, ierr)


       call cuda_col3(rey%x_d,this%d%x_d*y%x_d,this%m)
       call cuda_add2s1(rey%z_d,this%c%x_d,this%m)
       call cuda_sub2(rey%x_d,lambda%x_d,mu%x_d,this%m)
       rez = this%a0 - zeta - cuda_glsc2(lambda%x_d, this%a%x_d,this%m)

       ! relambda(:) = matmul(this%pij%x(:,:),1.0/(this%upp%x(:) - x(:))) + &
       !         matmul(this%qij%x(:,:), 1.0/(x(:) - this%low%x(:))) - &
       !         this%a%x(:)*z - y(:) + s(:) - this%bi%x(:)
       cons=0
       call cuda_cfill(relambda%x_d,cons, this%m)
       call cuda_relambda(relambda%x_d, x%x_d,  this%upp%x_d, this%low%x_d, this%pij%x_d, this%qij%x_d,  this%n, this%m)

       call device_memcpy(relambda%x, relambda%x_d, this%m, DEVICE_TO_HOST, sync=.false.)
       globaltmp_m = 0.0_rp
       call MPI_Allreduce(relambda%x, globaltmp_m, this%m, &
        mpi_real_precision, mpi_sum, neko_comm, ierr)
       call device_memcpy(this%a%x, this%a%x_d, this%m, DEVICE_TO_HOST, sync=.false.)
       call device_memcpy(y%x, y%x_d, this%m, DEVICE_TO_HOST, sync=.false.)
       call device_memcpy(this%bi%x, this%bi%x_d, this%m, DEVICE_TO_HOST, sync=.false.)
       relambda%x = globaltmp_m - this%a%x*z - y%x + s%x - this%bi%x
       call device_memcpy(relambda%x, relambda%x_d, this%m, HOST_TO_DEVICE, sync=.false.)



       call cuda_sub2cons2(rexsi%x_d,xsi%x_d,x%x_d,this%alpha%x_d,this%epsi)
       call cuda_sub2cons2(reeta%x_d,eta%x_d,this%beta%x_d,x%x_d,this%epsi)
       call cuda_sub2cons(remu%x_d,mu%x_d,y%x_d,this%epsi)
       rezeta = zeta*z - epsi
       call cuda_sub2cons(res%x_d,lambda%x_d,s%x_d,this%epsi)
       residu=maxval(cuda_maxval(rex%x_d,this%n), cuda_maxval(rey%x_d,this%n), rez, cuda_maxval(relambda%x_d,this%n), &
        cuda_maxval(rexsi%x_d,this%n), cuda_maxval(reeta%x_d,this%n), cuda_maxval(remu%x_d,this%n), rezeta, cuda_maxval(res%x_d,this%n))
       residumax = 0_rp
       call MPI_Allreduce(maxval(abs(residu)), residumax, 1, &
         mpi_real_precision, mpi_max, neko_comm, ierr)

       re_xstuff_squ_global = 0_rp

       globaltemp_norm=cuda_norm(rex%x_d,this%n)+cuda_norm(rexsi%x_d,this%n)+cuda_norm(reeta%x_d,this%n);
       call MPI_Allreduce(globaltemp_norm,&
         re_xstuff_squ_global, 1, mpi_real_precision, mpi_sum,&
         neko_comm, ierr)
       residu_small_norm=cuda_norm(rey%x_d,this%n)+norm2(rez)**2+cuda_norm(relambda%x_d,this%n)+cuda_norm(remu%x_d,this%n)&
       +norm2(rezeta)**2+cuda_norm(res%x_d,this%n)
       residunorm = sqrt(residu_small_norm + re_xstuff_squ_global)


       do iter = 1, this%max_iter !ittt
           if (iter .gt. (this%max_iter -2)) then
            ! print *, "The mma inner loop seems not to converge"
            ! print *, "residumax = ", residumax, "for epsi = ", epsi, &
            !         ", ittt  = ", iter, "out of ", this%max_iter
        end if
        !Check the condition
        if (residumax .lt. epsi) exit
        call cuda_delx(delx_d, x_d, xlow_d, xupp_d,  pij_d,  qij_d,  p0j_d, q0j_d, alpha_d,  beta_d, lambda_d, epsi, n, m)
        call cuda_dely(dely_d,  c_d, d_d, y_d, lambda_d, epsi, n)
        delz = this%a0 - cuda_glsc2(lambda_d, this%a%x_d, m) - epsi/z

        dellambda_d = 0.0_rp
        call cuda_dellambda( dellambda_d, x_d, xlow_d, xupp_d, pij_d, qij_d, n)

        globaltmp_m = 0.0_rp
        call MPI_Allreduce(dellambda, globaltmp_m, this%m, &
           mpi_real_precision, mpi_sum, neko_comm, ierr)

        dellambda = globaltmp_m - this%a%x*z - y - this%bi%x + epsi/lambda

        ! delx(:) = ((this%p0j%x(:) + matmul(transpose(this%pij%x(:,:)), &
        !     lambda(:)))/(this%upp%x(:) - x(:))**2 - &
        !     (this%q0j%x(:) + matmul(transpose(this%qij%x(:,:)), &
        !     lambda(:)))/(x(:) - this%low%x(:))**2 ) - &
        !     epsi/(x(:) - this%alpha%x(:)) + epsi/(this%beta%x(:) - x(:))

        ! dely(:) =  this%c%x(:) + this%d%x(:)*y(:) - lambda(:) - epsi/y(:)
        ! delz = this%a0 - dot_product(lambda(:), this%a%x(:)) - epsi/z
        ! dellambda(:) = matmul(this%pij%x(:,:),1.0/(this%upp%x(:) - x(:)))+&
        !     matmul(this%qij%x(:,:), 1.0/(x(:) - this%low%x(:))) - &
        !     this%a%x(:)*z - y(:) - this%bi%x(:) + epsi/lambda(:)
        call cuda_GG(GG_d,  x_d,  xlow_d,  xupp_d, pij_d, qij_d, n, m)

        call cuda_diagx(diagx_d, x_d, xsi_d, xlow_d, xupp_d, p0j_d,  q0j_d,  pij_d,  qij_d,  alpha_d, beta_d,  eta_d, lambda_d, n, m)


        !Here we only consider the case m<n in the matlab code
        !assembling the right hand side matrix based on eq(5.20)
        ! bb = [dellambda + dely(:)/(this%d%x(:) + &
        !         (mu(:)/y(:))) - matmul(GG,delx/diagx), delz ]
        !!!!!!!!!!!!!!for MPI computation of bb!!!!!!!!!!!!!!!!!!!!!!!!!
        call cuda_bb(bb_d, GG_d, delx_d,diagx_d,n,m)
        globaltmp_m = 0.0_rp
        call MPI_Allreduce(bb(1:this%m), globaltmp_m, this%m, &
           mpi_real_precision, mpi_sum, neko_comm, ierr)
        bb(1:this%m) = globaltmp_m

        bb(1:this%m) = dellambda + dely/(this%d%x + (mu/y)) - bb(1:this%m)
        bb(this%m +1) = delz
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! !assembling the coefficients matrix AA based on eq(5.20)
        ! AA(1:this%m,1:this%m) =  &
        ! matmul(matmul(GG,mma_diag(1/diagx)), transpose(GG))
        ! !update diag(AA)
        ! AA(1:this%m,1:this%m) = AA(1:this%m,1:this%m) + &
        !     mma_diag(s(:)/lambda(:) + 1.0/(this%d%x(:) + (mu(:)/y(:))))

        AA = 0.0_rp
        !Direct computation of the matrix multiplication
        !(for better performance)
        call cuda_AA(AA_d, GG_d,  diagx_d, n, m) 

        globaltmp_mm = 0.0_rp
        call MPI_Allreduce(AA(1:this%m, 1:this%m), globaltmp_mm, &
           this%m*this%m, mpi_real_precision, mpi_sum, neko_comm, ierr)
        AA(1:this%m,1:this%m) = globaltmp_mm
        do i = 1, this%m
         !update the diag AA
         AA(i, i) = AA(i, i) + (s(i) / lambda(i) + &
            1.0_rp / (this%d%x(i) + mu(i) / y(i)))
     end do

     AA(1:this%m, this%m+1) = this%a%x(:)
     AA(this%m+1, 1:this%m) = this%a%x(:)
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

     dy = (-dely+dlambda)/(this%d%x(:) + (mu(:)/y(:)))
     dxsi = -xsi + (epsi-dx*xsi(:))/(x(:) - this%alpha%x(:))
     deta = -eta + (epsi+dx*eta(:))/(this%beta%x(:) - x(:))
     dmu = -mu + (epsi-mu*dy(:))/y(:)
     dzeta = -zeta + (epsi-zeta*dz)/z
     ds = -s + (epsi-dlambda*s(:))/lambda(:)

     !2*this%n+4*this%m+2
     dxx = [dy, dz, dlambda, dxsi, deta, dmu, dzeta, ds]
     xx = [y, z, lambda, xsi, eta, mu, zeta, s]
     steg = maxval([dummy_one, -1.01*dxx/xx, -1.01*dx/ &
       (x(:) - this%alpha%x(:)), 1.01*dx/(this%beta%x(:) - x(:))])
     steg = 1.0/steg

     call MPI_Allreduce(steg, steg, 1, &
       mpi_real_precision, mpi_min, neko_comm, ierr)

     call cuda_copy   (  xold_d,x_d,n)
     y_old=y
     zold = z
     lambdaold = lambda
     call cuda_copy   (  xsiold_d,xsi_d,n)
     call cuda_copy   (  etaold_d,eta_d,n)
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
         rex(:) = ((this%p0j%x(:) + matmul(transpose(this%pij%x(:,:)), &
            lambda(:)))/(this%upp%x(:) - x(:))**2 - &
         (this%q0j%x(:) + matmul(transpose(this%qij%x(:,:)), &
            lambda(:)))/(x(:) - this%low%x(:))**2 ) - &
         xsi(:) + eta(:)
         rey(:) = this%c%x(:) + this%d%x(:)*y(:) - lambda(:) - mu(:)
         rez = this%a0 - zeta - dot_product(lambda(:), this%a%x(:))
         ! relambda(:) = matmul(this%pij%x(:,:),1.0/&
         !         (this%upp%x(:) - x(:))) + matmul(this%qij%x(:,:), &
         !         1.0/(x(:) - this%low%x(:))) - this%a%x(:)*z - &
         !         y(:) + s(:) - this%bi%x(:)
         relambda = 0.0_rp
         do i = 1, this%m
          do j = 1, this%n !this n is global
           ! Accumulate sums for relambda (the term gi(x))
           relambda(i) = relambda(i) + &
           this%pij%x(i,j)/(this%upp%x(j) - x(j)) &
           + this%qij%x(i,j)/(x(j) - this%low%x(j))
       end do
   end do
   globaltmp_m = 0.0_rp
   call MPI_Allreduce(relambda, globaltmp_m, this%m, &
     mpi_real_precision, mpi_sum, neko_comm, ierr)
   relambda = globaltmp_m


   relambda = relambda - this%a%x(:)*z - y(:) + s(:) - this%bi%x(:)

   rexsi(:) = xsi(:)*(x(:) - this%alpha%x(:)) - epsi
   reeta(:) = eta(:)*(this%beta%x(:) - x(:)) - epsi
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

! Save the new design
this%xold2 = this%xold1
this%xold1%x = designx
designx = x

    !update the parameters of the MMA object nesessary to compute KKT residu
this%y%x = y
this%z = z
this%lambda%x = lambda
this%zeta = zeta
this%xsi%x = xsi
this%eta%x = eta
this%mu%x = mu
this%s%x = s

end subroutine mma_subsolve_dpip_cpu

subroutine mma_KKT_cpu(this, x, df0dx, fval, dfdx)
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

    real(kind=rp) :: rez, rezeta
    real(kind=rp), dimension(this%m) :: rey, relambda, remu, res
    real(kind=rp), dimension(this%n) :: rex, rexsi, reeta
    real(kind=rp), dimension(3*this%n+4*this%m+2) :: residu

    real(kind=rp), dimension(4*this%m+2) :: residu_small
    integer :: ierr
    real(kind=rp) :: re_xstuff_squ_global

    rex(:) = df0dx + matmul(transpose(dfdx), this%lambda%x(:)) - this%xsi%x(:) + &
    this%eta%x(:)
    rey(:) = this%c%x(:) + this%d%x(:)*this%y%x(:) - this%lambda%x(:) - this%mu%x(:)
    rez = this%a0 - this%zeta - dot_product(this%lambda%x(:), this%a%x(:))

    relambda(:) = fval - this%a%x(:)*this%z - this%y%x(:) + this%s%x(:)
    rexsi(:) = this%xsi%x(:)*(x(:) - this%xmin%x(:))
    reeta(:) = this%eta%x(:)*(this%xmax%x(:) - x(:))
    remu(:) = this%mu%x(:)*this%y%x(:)
    rezeta = this%zeta*this%z
    res(:) = this%lambda%x(:)*this%s%x(:)

    residu = [rex, rey, rez, relambda, rexsi, reeta, remu, rezeta, res]

    call MPI_Allreduce(maxval(abs(residu)), this%residumax, 1, &
       mpi_real_precision, mpi_max, neko_comm, ierr)

    call MPI_Allreduce(norm2(rex)**2 + norm2(rexsi)**2 + norm2(reeta)**2, &
       re_xstuff_squ_global, 1, mpi_real_precision, mpi_sum, neko_comm, ierr)

    residu_small = [rey, rez, relambda, remu, rezeta, res]

    this%residunorm = sqrt(norm2(residu_small)**2 + re_xstuff_squ_global)

end subroutine mma_KKT_cpu
end submodule mma_cpu
