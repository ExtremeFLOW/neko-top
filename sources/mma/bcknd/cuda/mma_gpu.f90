submodule (mma) mma_gpu
use mpi_f08, only: MPI_INTEGER, MPI_REAL, mpi_sum, mpi_min, mpi_max, &
MPI_Allreduce
use utils, only: neko_error
use comm, only: neko_comm, mpi_real_precision

contains
subroutine mma_gensub_gpu(this, iter, x, df0dx, fval, dfdx)
    ! ----------------------------------------------------- !
    ! Generate the approximation sub problem by computing   !
    ! the lower and upper asymtotes and the other necessary !
    ! parameters (alpha, beta, p0j, q0j, pij, qij, ...).    !
    ! ----------------------------------------------------- !
    class(mma_t), intent(inout) :: this
    type(vector_t), intent(in) :: x
    type(vector_t), intent(in) :: df0dx
    type(vector_t), intent(in) :: fval
    type(matrix_t), intent(in) :: dfdx

    integer, intent(in) :: iter
    integer :: i, j, ierr
    type(vector_t) :: globaltmp_m
    real(kind=rp) :: cons



    call globaltmp_m%init(this%m)
    if (iter .lt. 3) then
        call device_add3s2(this%low%x_d,this%xmax%x_d,this%xmin%x_d,-this%asyinit,this%asyinit,this%n)
        call device_add2(this%low%x_d,x%x_d,this%n)

        call device_add3s2( this%upp%x_d,this%xmax%x_d,this%xmin%x_d,this%asyinit,- this%asyinit,this%n)
        call device_add2(this%upp%x_d,x%x_d,this%n)
        !call mma_gensub1_gpu(this%low%x_d, this%upp%x_d,x%x_d, this%xmin%x_d, this%xmax%x_d, this%asyinit, this%n)
    else
      call mma_gensub2_gpu(this%low%x_d, this%upp%x_d, x%x_d, this%xold1%x_d, this%xold2%x_d,this%xmin%x_d, this%xmax%x_d, &
        this%asydecr, this%asyincr, this%n)
  end if
  call mma_gensub3_gpu(x%x_d, df0dx%x_d, dfdx%x_d,this%low%x_d, this%upp%x_d, this%xmin%x_d, this%xmax%x_d,this%alpha%x_d, &
    this%beta%x_d, this%p0j%x_d, this%q0j%x_d, this%pij%x_d, this%qij%x_d, this%n, this%m) 
  call mma_gensub4_gpu(x%x_d, this%low%x_d, this%upp%x_d, this%pij%x_d, this%qij%x_d, this%n, this%m, this%bi%x_d)



  !!!!!cpu gpu transfer part
  ! globaltmp_m%x=0.0_rp
  ! call device_memcpy(this%bi%x, this%bi%x_d, this%m, DEVICE_TO_HOST, sync=.false.)
  ! call MPI_Allreduce(this%bi%x, globaltmp_m%x, this%m, &
  !   mpi_real_precision, mpi_sum, neko_comm, ierr)
  ! call device_memcpy(fval%x, fval%x_d, this%m, DEVICE_TO_HOST, sync=.false.)
  ! this%bi%x = globaltmp_m%x - fval%x

  call cuda_mpisum(this%bi%x_d, this%m)
  call device_sub(this%bi%x_d, fval%x_d, this%m)
end subroutine mma_gensub_gpu

subroutine mma_subsolve_dpip_gpu(this, designx)
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
    type(vector_t) :: globaltmp_m

    type(vector_t), intent(in) :: x, xsi, eta, & !!!!!n
    rex, rexsi, reeta, &
    delx, diagx, dx, dxsi, deta, &
    xold, xsiold, etaold

    real(kind=rp) :: :: residu
    real(kind=rp) :: :: residu_small
    real(kind=rp) :: :: ratio_xx_dxx

    type(vector_t) :: bb
    type(matrix_t) :: GG
    type(matrix_t) :: AA
    real(kind=rp), dimension(this%m, this%m) :: globaltmp_mm
    integer :: info
    integer, dimension(this%m+1) :: ipiv
    real(kind=rp) :: re_xstuff_squ_global
    integer :: nglobal
    real(kind=rp) :: cons, cons1
    real(kind=rp) ::residu_temp
    real(kind=rp) ::globaltemp_norm

    call globaltmp_m%init(this%m)
    
    call y%init(this%m)
    call lambda%init(this%m)
    call s%init(this%m)
    call mu%init(this%m)
    call rey%init(this%m)
    call relambda%init(this%m)
    call remu%init(this%m)
    call res%init(this%m)
    call dely%init(this%m)
    call dellambda%init(this%m)
    call dy%init(this%m)
    call dlambda%init(this%m)
    call ds%init(this%m)
    call dmu%init(this%m)
    call yold%init(this%m)
    call lambdaold%init(this%m)
    call sold%init(this%m)
    call muold%init(this%m)




    call x%init(this%n)
    call xsi%init(this%n)
    call eta%init(this%n)
    call rex%init(this%n)
    call rexsi%init(this%n)
    call reeta%init(this%n)
    call delx%init(this%n)
    call diagx%init(this%n)
    call dx%init(this%n)
    call dxsi%init(this%n)
    call deta%init(this%n)
    call xold%init(this%n)
    call xsiold%init(this%n)
    call etaold%init(this%n)





    !call residu%init(3*this%n+4*this%m+2)
    !call residu_small%init(4*this%m+2)
    !call xx%init(2*this%n+4*this%m+2)
    !call dxx%init(2*this%n+4*this%m+2)

    call bb%init(this%m+1)
    call GG%init(this%m, this%n)
    call AA%init(this%m+1, this%m+1)
    ! using DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO ) in lapack to solve
    ! the linear system which needs the following parameters
    

    ! intial value for the parameters in the subsolve based on
    ! page 15 of "https://people.kth.se/~krille/mmagcmma.pdf"
    dummy_one = 1
    epsi = 1 !100
    call device_add3s2(x%x_d,this%alpha%x_d,this%beta%x_d,0.5_rp,0.5_rp,this%n)
    call device_cfill(y%x_d,1.0_rp,this%m)
    z = 1
    zeta = 1
    call device_cfill(lambda%x_d,1.0_rp,this%m)
    call device_fill(s%x_d,1.0_rp,this%m)

    call cuda_mma_max(xsi%x_d,x%x_d,this%alpha%x_d,this%n)
    call cuda_mma_max(eta%x_d,this%beta%x_d,x%x_d,this%n)



    !!!!cpu 
    ! mu%x = max(1.0, 0.5*this%c%x)

    cons=1.0
    cons1=0.5
    call cuda_maxcons(mu%x_d,cons,cons1, this%c%x_d, this%n) 

    do while (epsi .gt. 0.9*this%epsimin)
       ! calculating residuals based on
       ! "https://people.kth.se/~krille/mmagcmma.pdf" for the variables
       ! x, y, z, lambda residuals based on eq(5.9a)-(5.9d), respectively.
       call cuda_rex(rex%x_d,  x%x_d,  this%low%x_d, this%upp%x_d,  this%pij%x_d, this%p0j%x_d,this%qij%x_d, this%q0j%x_d, &
        lambda%x_d, xsi%x_d, eta%x_d, this%n, this%m) 

       call MPI_Allreduce(this%n, nglobal, 1, MPI_INTEGER, mpi_sum, neko_comm, ierr)

       !!!cpu 
       !rey%x = this%c%x + this%d%x*y%x - lambda%x - mu%x

       call cuda_sub3(rey%x_d, this%c%x_d, lambda%x_d, this%m)
       call cuda_sub2(rey%x_d, mu%x_d, this%m)
       call cuda_addcol3(rey%x_d, this%d%x_d,y%x_d)


       rez = this%a0 - zeta - cuda_lcsc2(lambda%x_d, this%a%x_d,this%m)

       ! relambda(:) = matmul(this%pij%x(:,:),1.0/(this%upp%x(:) - x(:))) + &
       !         matmul(this%qij%x(:,:), 1.0/(x(:) - this%low%x(:))) - &
       !         this%a%x(:)*z - y(:) + s(:) - this%bi%x(:)
       call device_cfill(relambda%x_d,0.0_rp, this%m)
       call cuda_relambda(relambda%x_d, x%x_d,  this%upp%x_d, this%low%x_d, this%pij%x_d, this%qij%x_d,  this%n, this%m)

       call cuda_mpisum(relambda%x_d, this%m)
       call device_add2s2(relambda%x_d, this%a%x_d,-z, this%m)
       call device_sub2(relambda%x_d, y%x_d, this%m)
       call device_add2(relambda%x_d, s%x_d, this%m)
       call device_sub2(relambda%x_d, this%bi%x_d, this%m)

       !!!!cpu gpu tranfer
       ! call device_memcpy(relambda%x, relambda%x_d, this%m, DEVICE_TO_HOST, sync=.false.)
       ! globaltmp_m%x = 0.0_rp
       ! call MPI_Allreduce(relambda%x, globaltmp_m%x, this%m, &
       !  mpi_real_precision, mpi_sum, neko_comm, ierr)
       ! relambda%x = globaltmp_m%x - this%a%x*z - y%x + s%x - this%bi%x
       ! call device_memcpy(relambda%x, relambda%x_d, this%m, HOST_TO_DEVICE, sync=.false.)



       call cuda_sub2cons2(rexsi%x_d,xsi%x_d,x%x_d,this%alpha%x_d,this%epsi,this%n)
       call cuda_sub2cons2(reeta%x_d,eta%x_d,this%beta%x_d,x%x_d,this%epsi,this%n)


       !cpu
       !!remu%x = mu%x*y%x- epsi

       call device_cmult2(remu%x_d, mu%x_d, y%x_d)
       call device_cadd(remu%x_d, -epsi)


       rezeta = zeta*z - epsi


       !!cpu
       !res%x = lambda%x*s%x - epsi

       call device_cmult2(res%x_d, lambda%x_d, s%x_d)
       call device_cadd(res%x_d, -epsi)


       residumax%x = 0.0_rp
       call MPI_Allreduce(maxval(abs([cuda_maxval(rex%x_d,this%n), cuda_maxval(rey%x_d,this%m), rez, cuda_maxval(relambda%x_d,this%m), &
        cuda_maxval(rexsi%x_d,this%n), cuda_maxval(reeta%x_d,this%n), cuda_maxval(remu%x_d,this%m), &
        rezeta, cuda_maxval(res%x_d,this%m)])), residumax%x, 1, &
       mpi_real_precision, mpi_max, neko_comm, ierr)

       re_xstuff_squ_global = 0.0_rp

       globaltemp_norm = cuda_norm(rex%x_d,this%n) + cuda_norm(rexsi%x_d,this%n)+cuda_norm(reeta%x_d,this%n)
       call MPI_Allreduce(globaltemp_norm,&
         re_xstuff_squ_global, 1, mpi_real_precision, mpi_sum,&
         neko_comm, ierr)

       this%residunorm = sqrt(cuda_norm(rey%x_d,this%m)+norm2(rez)**2+cuda_norm(relambda%x_d,this%m)+cuda_norm(remu%x_d,this%m)&
        +norm2(rezeta)**2+cuda_norm(res%x_d,this%m + re_xstuff_squ_global))

       do iter = 1, this%max_iter !ittt
           if (iter .gt. (this%max_iter -2)) then
            ! print *, "The mma inner loop seems not to converge"
            ! print *, "residumax = ", residumax, "for epsi = ", epsi, &
            !         ", ittt  = ", iter, "out of ", this%max_iter
        end if
        !Check the condition
        if (residumax .lt. epsi) exit
        call cuda_delx(delx%x_d, x%x_d, this%low%x_d, this%upp%x_d,  this%pij%x_d,  this%qij%x_d,  this%p0j%x_d, this%q0j%x_d, &
            this%alpha%x_d,  this%beta%x_d, lambda%x_d, epsi, this%n, this%m)



        !dely%x = this%c%x + this%d%x*y%x - lambda%x - epsi/y%x
        !call cuda_dely(dely%x_d,  this%c%x_d, this%d%x_d, y%x_d, lambda%x_d, epsi, this%n)

        call device_copy(dely%x_d,y%x_d, this%m)
        call device_invcol1(dely%x_d, this%m)
        call device_cmult(dely%x_d, -epsi, this%m)
        call device_add2(dey%x_d, this%c%x_d, this%m)
        call device_addcol3(dely%x_d, this%d%x_d,this%y%x_d, this%m)
        call device_sub2(dely%x%x_d, lambda%x_d, this%m)


        delz = this%a0 - cuda_lcsc2(lambda%x_d, this%a%x_d, this%m) - epsi/z

        call device_cfill(relambda%x_d,0.0_rp, this%m)
        call cuda_dellambda( dellambda%x_d, x%x_d, this%low%x_d, this%upp%x_d, this%pij%x_d, this%qij%x_d, this%n)


!!!cpu gpu transfer
        ! call device_memcpy(dellambda%x, dellambda%x_d, this%m, DEVICE_TO_HOST, sync=.false.)
        ! globaltmp_m%x = 0.0_rp
        ! call MPI_Allreduce(dellambda%x, globaltmp_m%x, this%m, &
        !    mpi_real_precision, mpi_sum, neko_comm, ierr)
        ! dellambda%x = globaltmp_m%x - this%a%x*z - y%x - this%bi%x + epsi/lambda%x
        ! call device_memcpy(dellambda%x, dellambda%x_d, this%m, HOST_TO_DEVICE, sync=.false.)

        call cuda_mpisum(dellambda%x_d, this%m)
        call device_add2s2(dellambda%x_d, this%a%x_d, -z, this%m)
        call device_sub2(dellambda%x_d, y%x_d, this%m)
        call device_sub2(dellambda%x_d, this%bi%x_d, this%m)
        call cuda_add2inv2(dellambda%x_d, lambda%x_d, epsi, this%m)

        call cuda_GG(GG%x_d,  x%x_d,  this%low%x_d,  this%upp%x_d, this%pij%x_d, this%qij%x_d, this%n, this%m)

        call cuda_diagx(diagx%x_d, x%x_d, xsi%x_d, this%low%x_d, this%upp%x_d, this%p0j%x_d, this%q0j%x_d,  this%pij%x_d,&
          this%qij%x_d,  this%alpha%x_d, this%beta%x_d,  eta%x_d, lambda%x_d, this%n, this%m)


        !Here we only consider the case m<n in the matlab code
        !assembling the right hand side matrix based on eq(5.20)
        ! bb = [dellambda + dely(:)/(this%d%x(:) + &
        !         (mu(:)/y(:))) - matmul(GG,delx/diagx), delz ]
        !!!!!!!!!!!!!!for MPI computation of bb!!!!!!!!!!!!!!!!!!!!!!!!!
        call device_cfill(bb%x_d,0.0_rp,this%m+1)
        call cuda_bb(bb%x_d, GG%x_d, delx%x_d,diagx%x_d,this%n,this%m)
        call device_memcpy(bb%x, bb%x_d, this%m, DEVICE_TO_HOST, sync=.false.)

        globaltmp_m%x = 0.0_rp
        call MPI_Allreduce(bb%x(1:this%m), globaltmp_m%x, this%m, &
           mpi_real_precision, mpi_sum, neko_comm, ierr)
        bb%x(1:this%m) = globaltmp_m%x
        bb%x(1:this%m) = dellambda%x + dely%x/(this%d%x + (mu%x/y%x)) - bb%x(1:this%m)
        bb%x(this%m +1) = delz
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! !assembling the coefficients matrix AA based on eq(5.20)
        ! AA(1:this%m,1:this%m) =  &
        ! matmul(matmul(GG,mma_diag(1/diagx)), transpose(GG))
        ! !update diag(AA)
        ! AA(1:this%m,1:this%m) = AA(1:this%m,1:this%m) + &
        !     mma_diag(s(:)/lambda(:) + 1.0/(this%d%x(:) + (mu(:)/y(:))))
        call device_cfill(AA%x_d,0.0_rp, (this%m+1)*(this%m+1))
        !Direct computation of the matrix multiplication
        !(for better performance)
        call cuda_AA(AA%x_d, GG%x_d,  diagx%x_d, this%n, this%m) 
        call device_memcpy(AA%x, AA%x_d, this%m, DEVICE_TO_HOST, sync=.false.)

        globaltmp_mm = 0.0_rp
        call MPI_Allreduce(AA%x(1:this%m, 1:this%m), globaltmp_mm, &
           this%m*this%m, mpi_real_precision, mpi_sum, neko_comm, ierr)
        AA%x(1:this%m,1:this%m) = globaltmp_mm
        do i = 1, this%m
         !update the diag AA
         AA%x(i, i) = AA%x(i, i) + (s%x(i) / lambda%x(i) + 1.0_rp / (this%d%x(i) + mu%x(i) / y%x(i)))
     end do

     AA%x(1:this%m, this%m+1) = this%a%x(:)
     AA%x(this%m+1, 1:this%m) = this%a%x(:)
     AA%x(this%m+1, this%m+1) = -zeta/z
     call device_memcpy(AA%x, AA%x_d, (this%m+1)*(this%m+1), HOST_TO_DEVICE, sync=.false.)


     call DGESV(this%m+1, 1, AA%x, this%m+1, ipiv, bb%x, this%m+1, info)
     ! if info! = 0 then DGESV is failed.
     if (info .ne. 0) then
         write(stderr, *) "DGESV failed to solve the linear system in MMA."
         write(stderr, *) "Please check mma_subsolve_dpip in mma.f90"
         error stop
     end if

     dlambda%x = bb%x(1:this%m)
     call cuda_copy(dlambda%x_d, bb%x_d, this%m)

     dz = bb%x(this%m + 1)
     ! based on eq(5.19)

     call device_memcpy(dlambda%x, dlambda%x_d, this%m, HOST_TO_DEVICE, sync=.false.)
     call cuda_dx(dx%x_d, delx%x_d, diagx%x_d, GG%x_d, dlambda%x_d, this%n, this%m)

     dy%x = (-dely%x+dlambda%x)/(this%d%x(:) + (mu%x(:)/y%x(:)))
     call device_memcpy(dy%x, dy%x_d, this%m, HOST_TO_DEVICE, sync=.false.)


     call cuda_dxsi(dxsi%x_d, xsi%x_d, dx%x_d,x%x_d,this%alpha%x_d, epsi, this%n) 
     call cuda_deta(deta%x_d, eta%x_d, dx%x_d,  x%x_d, beta%x_d, epsi,this%n)

     dmu%x = -mu%x + (epsi-mu%x*dy%x(:))/y%x(:)
     call device_memcpy(dmu%x, dmu%x_d, this%m, HOST_TO_DEVICE, sync=.false.)


     dzeta = -zeta + (epsi-zeta*dz)/z
     ds%x = -s%x + (epsi-dlambda%x*s%x(:))/lambda%x(:)
     call device_memcpy(ds%x, ds%x_d, this%m, HOST_TO_DEVICE, sync=.false.)


     !2*this%n+4*this%m+2
     cons =-1.01
     steg=maxval(cons*[dummy_one/cons,dy%x/dy%x,dz/dz,dlambda%x/lambda,dmu%x/mu%x,dzeta/zeta,ds%x/s%x])
     steg=maxval([steg,cuda_maxval2(dxx%x_d, xx%x_d, cons, this%n)])
     steg=maxval([steg,cuda_maxval3(dx%x_d, x%x_d, this%alpha%x_d,cons, this%n)])
     cons=1.01
     steg=maxval([steg,cuda_maxval3(dx%x_d, this%beta%x_d,x%x_d, cons, this%n)])

     !dxx = [dy, dz, dlambda, dxsi, deta, dmu, dzeta, ds]
     !xx = [y, z, lambda, xsi, eta, mu, zeta, s]
     !steg = maxval([dummy_one, -1.01*dxx/xx, -1.01*dx/ &
     !  (x(:) - this%alpha%x(:)), 1.01*dx/(this%beta%x(:) - x(:))])
     steg = 1.0/steg

     call MPI_Allreduce(steg, steg, 1, &
       mpi_real_precision, mpi_min, neko_comm, ierr)

     call cuda_copy(xold%x_d,x%x_d,this%n)

     call cuda_copy(yold%x_d,y%x_d,this%m)
     yold%x=y%x

     zold = z

     call cuda_copy(lambdaold%x_d,lambda%x_d,this%m)
     lambdaold%x=lambda%x

     call cuda_copy(xsiold%x_d,xsi%x_d,this%n)
     call cuda_copy(etaold%x_d,eta%x_d,this%n)

     call cuda_copy(muold%x_d,mu%x_d,this%m)
     muold%x=mu%x

     zetaold = zeta

     sold%s = s%x
     call cuda_copy(sold%x_d,s%x_d,this%m)

     !The innermost loop to determine the suitable step length
     !using the Backtracking Line Search approach
     newresidu = 2*residunorm
     itto = 0
     do while ((newresidu .gt. residunorm) .and. (itto .lt. 50))
         itto = itto + 1
         !update the variables
         cons=1
         call add3s2(x%x_d,xold%x_d,dx%x_d,cons,steg,this%n)

         call add3s2(y%x_d,yold%x_d,dy%x_d,cons,steg,this%m)
         y%x = yold%x + steg*dy%x

         z = zold + steg*dz

         call add3s2(lambda%x_d,lambdaold%x_d,dy%x_d,cons,steg,this%m)
         lambda%x = lambdaold%x + steg*dlambda%x

         call add3s2(xsi%x_d,xsiold%x_d,dxsi%x_d,cons,steg)
         call add3s2(eta%x_d,etaold%x_d,deta%x_d,cons,steg)

         call add3s2(mu%x_d,muold%x_d,dy%x_d,cons,steg,this%m)
         mu%x = muold%x + steg*dmu%x

         zeta = zetaold + steg*dzeta

         call add3s2(s%x_d,sold%x_d,dy%x_d,cons,steg,this%m)
         s%x = sold%x + steg*ds%x

         !recompute the newresidu to see if this stepsize improves
         !the residue
         call cuda_rex(rex%x_d,  x%x_d,  this%low%x_d, this%upp%x_d,  this%pij%x_d, this%p0j%x_d,this%qij%x_d, this%q0j%x_d, &
            lambda%x_d, xsi%x_d, eta%x_d, this%n, this%m) 

         rey%x = this%c%x + this%d%x*y%x - lambda%x - mu%x
         call cuda_sub3(rey%x_d, this%c%x_d, lambda%x_d, this%m)
         call cuda_sub2(rey%x_d, mu%x_d, this%m)
         call cuda_addcol3(rey%x_d, this%d%x_d,y%x_d)

         rez = this%a0 - zeta - cuda_lcsc2(lambda%x_d, this%a%x_d, this%m)

         call device_cfill(relambda%x_d,0.0_rp, this%m)
         call cuda_relambda(relambda%x_d, x%x_d,  this%upp%x_d, this%low%x_d, this%pij%x_d, this%qij%x_d,  this%n, this%m)
         call device_memcpy(relambda%x, relambda%x_d, this%m, DEVICE_TO_HOST, sync=.false.)

         globaltmp_m%x= 0.0_rp
         call MPI_Allreduce(relambda%x, globaltmp_m%x, this%m, &
             mpi_real_precision, mpi_sum, neko_comm, ierr)
         relambda%x = globaltmp_m%x


         relambda%x = relambda%x - this%a%x*z - y%x + s%x - this%bi%x
         call device_memcpy(relambda%x, relambda%x_d, this%m, HOST_TO_DEVICE, sync=.false.)


         call cuda_sub2cons2(rexsi%x_d,xsi%x_d,x%x_d,this%alpha%x_d,this%epsi,this%n)
         call cuda_sub2cons2(reeta%x_d,eta%x_d,this%beta%x_d,x%x_d,this%epsi,this%n)

         call cuda_sub2cons(remu%x_d,mu%x_d,y%x_d,this%epsi,this%m)
         remu%x = mu%x*y%x - epsi

         rezeta = zeta*z - epsi


         res%x = lambda%x*s%x - epsi
         call cuda_sub2cons(res%x_d,lambda%x_d,s%x_d,this%epsi,this%m)

         residu_temp=0
         residu_temp=maxval(cuda_maxval(rex%x_d,this%n), rey%x, rez, relambda%x, &
            cuda_maxval(rexsi%x_d,this%n), cuda_maxval(reeta%x_d,this%n), remu%x, rezeta, res%x)


         re_xstuff_squ_global = 0_rp    
         globaltemp_norm=cuda_norm(rex%x_d,this%n)+cuda_norm(rexsi%x_d,this%n)+cuda_norm(reeta%x_d,this%n)

         call MPI_Allreduce(globaltemp_norm,re_xstuff_squ_global, 1, mpi_real_precision, mpi_sum,&
             neko_comm, ierr)

         !residu_small = [rey%x, rez%x, relambda%x, remu%x, rezeta%x, res%x]

         newresidu=sqrt(cuda_norm(rey%x_d,this%m)+norm2(rez)**2+cuda_norm(relambda%x_d,this%m)+cuda_norm(remu%x_d,this%m)&
            +norm2(rezeta)**2+cuda_norm(res%x_d,this%m + re_xstuff_squ_global))
         steg = steg/2
     end do

     residunorm = newresidu
     residumax = 0_rp
     call MPI_Allreduce(maxval(abs(residu_temp)), residumax, 1, &
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
call cuda_copy(this%xold2%x_d,this%xold1%x_d,this%n)
call cuda_copy(this%xold1%x_d,designx%x_d,this%n)
call cuda_copy(designx%x_d,x%x_d,this%n)
    !update the parameters of the MMA object nesessary to compute KKT residu

call cuda_copy(this%y%x_d,y%x_d,this%m)

this%z = z

call cuda_copy(this%lambda%x_d,lambda%x_d,this%m)

this%zeta = zeta

call cuda_copy(this%xsi%x_d,xsi%x_d,this%n)
call cuda_copy(this%eta%x_d,eta%x_d,this%n)

call cuda_copy(this%mu%x_d,mu%x_d,this%m)
call cuda_copy(this%s%x_d,s%x_d,this%m)


!!!!!!!!!!remove cpu later
this%y%x = y%x
this%lambda%x = lambda%x
this%mu%x = mu%x
this%s%x = s%x

end subroutine mma_subsolve_dpip_gpu

subroutine mma_KKT_gpu(this, x, df0dx, fval, dfdx)
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
    type(vector_t), intent(in) :: x
    type(vector_t), intent(in) :: fval
    type(vector_t), intent(in) :: df0dx
    type(matrix_t), intent(in) :: dfdx
    real(kind=rp) :: rez, rezeta
    type(vector_t) :: rey, relambda, remu, res
    type(vector_t) :: rex, rexsi, reeta
    real(kind=rp) ::residu_val !!!(3*this%n+4*this%m+2)
    real(kind=rp), dimension(4*this%m+2) ::residu_small !!!(4*this%m+2)
    integer :: ierr
    real(kind=rp) :: re_xstuff_squ_global
    real(kind=rp) :: cons
    real(kind=rp) :: globaltemp_norm


    call rey%init(this%m)
    call relambda%init(this%m)
    call remu%init(this%m)
    call res%init(this%m)

    call rex%init(this%n)
    call rexsi%init(this%n)
    call reeta%init(this%n)




    call cuda_kkt_rex(rex%x_d,  df0dx%x_d,  dfdx%x_d, this%xsi%x_d, this%eta%x_d, this%lambda%x_d, this%n, this%m)

    call cuda_sub3(rey%x_d, this%c%x_d, this%lambda%x_d, this%m)
    call cuda_sub2(rey%x_d, this%mu%x_d, this%m)
    call cuda_addcol3(rey%x_d, this%d%x_d,this%y%x_d)

    rez = this%a0 - this%zeta - cuda_lcsc2(this%lambda%x_d, this%a%x_d,this%m)

    cons=-this%z
    call cuda_sub3(relambda%x_d,fval%x_d,this%y%x_d,this%m)
    call cuda_add2s2(relambda%x_d, this%a%x_d, cons,this%m)
    call cuda_add2s2(relambda%x_d, this%s%x, cons, this%m)
    !relambda%x = fval%x - this%a%x*this%z - this%y%x + this%s%x

    call cuda_sub3(rexsi%x_d,x%x_d,this%xmin%x_d,this%n)
    call cuda_col2(rexsi%x_d, this%xsi%x_d,this%n)

    call cuda_sub3(reeta%x_d,this%xmax%x_d,x%x_d,this%n)
    call cuda_col2(reeta%x_d, this%eta%x_d,this%n)

    call cuda_col3(remu%x_d,this%mu%x_d,this%y%x_d,this%m)

    rezeta = this%zeta*this%z

    call cuda_col3(res%x_d,this%lambda%x_d,this%s%x_d,this%m)




    !!!!!!!cpu, remove later!!!!!!!!!!
    rey%x = this%c%x + this%d%x*this%y%x - this%lambda%x - this%mu%x

    call device_memcpy(fval%x, fval%x_d, this%m, DEVICE_TO_HOST, sync=.false.)
    relambda%x = fval%x - this%a%x*this%z - this%y%x + this%s%x
    remu%x = this%mu%x*this%y%x
    res%x = this%lambda%x*this%s%x
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    residu_val=maxval(cuda_maxval(rex%x_d,this%n), cuda_maxval(rey%x_d,this%m), rez, cuda_maxval(relambda%x_d,this%m), &
        cuda_maxval(rexsi%x_d,this%n), cuda_maxval(reeta%x_d,this%n), cuda_maxval(remu%x_d,this%m), &
        rezeta, cuda_maxval(res%x_d,this%m)
        !residu = [rex, rey, rez, relambda, rexsi, reeta, remu, rezeta, res]

        call MPI_Allreduce(maxval(abs(residu_val)), this%residumax, 1, &
           mpi_real_precision, mpi_max, neko_comm, ierr)

        globaltemp_norm=cuda_norm(rex%x_d,this%n)+cuda_norm(rexsi%x_d,this%n)+cuda_norm(reeta%x_d,this%n)
        call MPI_Allreduce(globaltemp_norm, &
           re_xstuff_squ_global, 1, mpi_real_precision, mpi_sum, neko_comm, ierr)



        !residu_small = [rey%x, rez%x, relambda%x, remu%x, rezeta%x, res%x]

        this%residunorm = sqrt(cuda_norm(rey%x_d,this%m)+norm2(rez)**2+cuda_norm(relambda%x_d,this%m)+cuda_norm(remu%x_d,this%m)&
           +norm2(rezeta)**2+cuda_norm(res%x_d,this%m + re_xstuff_squ_global)

       end subroutine mma_KKT_gpu
       end submodule mma_gpu
