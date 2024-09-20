submodule (mma) mma_cpu

contains
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
    integer :: i, j, ierr
    real(kind=rp), dimension(this%m) :: globaltmp_m

    if (iter .lt. 3) then
       do j = 1, this%n
          this%low(j) = x(j) - this%asyinit * (this%xmax(j) - &
               this%xmin(j))
          this%upp(j) = x(j) + this%asyinit * (this%xmax(j) - &
               this%xmin(j))
       end do
    else
       !Move asymptotes low and upp
       do j = 1, this%n
          if ((x(j) - this%xold1(j))*(this%xold1(j) - this%xold2(j)) &
               .lt. 0) then
             this%low(j) = x(j) - &
                  this%asydecr * (this%xold1(j) - this%low(j))
             this%upp(j) = x(j) + &
                  this%asydecr * (this%upp(j) - this%xold1(j))

          else if ((x(j) - this%xold1(j))* &
               (this%xold1(j) - this%xold2(j)) .gt. 0) then
             this%low(j) = x(j) - &
                  this%asyincr * (this%xold1(j) - this%low(j))
             this%upp(j) = x(j) + &
                  this%asyincr * (this%upp(j) - this%xold1(j))
          else
             this%low(j) = x(j) - (this%xold1(j) - this%low(j))
             this%upp(j) = x(j) + (this%upp(j) - this%xold1(j))
          end if

          ! setting a minimum and maximum for the low and upp
          ! asymptotes (eq3.9)
          this%low(j) = max(this%low(j), &
               x(j) - 10*(this%xmax(j) - this%xmin(j)))
          this%low(j) = min(this%low(j), &
               x(j) - 0.01*(this%xmax(j) - this%xmin(j)))

          this%upp(j) = min(this%upp(j), &
               x(j) + 10*(this%xmax(j) - this%xmin(j)))
          this%upp(j) = max(this%upp(j), &
               x(j) + 0.01*(this%xmax(j) - this%xmin(j)))
       end do
    end if
    ! we can move alpha and beta out of the following loop if needed as:
    ! this%alpha = max(this%xmin, this%low + &
    !     0.1*(this%x- this%low), this%x - 0.5*(this%xmax - this%xmin))
    ! this%beta = min(this%xmax, this%upp -  &
    !     0.1*(this%upp - this%x), this%x + 0.5*(this%xmax - this%xmin))
    do j = 1, this%n
       ! set the the bounds and coefficients for the approximation
       ! the move bounds (alpha and beta )are slightly more restrictive
       ! than low and upp. This is done based on eq(3.6)--eq(3.10).
       ! also check
       ! https://comsolyar.com/wp-content/uploads/2020/03/gcmma.pdf
       ! eq (2.8) and (2.9)
       this%alpha(j) = max(this%xmin(j), this%low(j) + &
            0.1*(x(j)- this%low(j)), &
            x(j) - 0.5*(this%xmax(j) - this%xmin(j)))
       this%beta(j) = min(this%xmax(j), this%upp(j) - &
            0.1*(this%upp(j) - x(j)), &
            x(j) + 0.5*(this%xmax(j) - this%xmin(j)))

       !Calculate p0j, q0j, pij, qij
       !where j = 1,2,...,n and i = 1,2,...,m  (eq(2.3)-eq(2.5))
       this%p0j(j) = (this%upp(j) - x(j))**2 * &
            (1.001*max(df0dx(j),0.0) + &
            0.001*max(-df0dx(j),0.0) + &
            (0.00001/(max(0.00001, &
            (this%xmax(j) - this%xmin(j))))))

       this%q0j(j) = (x(j) - this%low(j))**2 * &
            (0.001*max(df0dx(j),0.0) + &
            1.001*max(-df0dx(j),0.0) + &
            (0.00001/(max(0.00001, &
            (this%xmax(j) - this%xmin(j))))))

       do i = 1, this%m
          this%pij(i,j) = (this%upp(j) - x(j))**2 * &
               (1.001*max(dfdx(i,j),0.0) + &
               0.001*max(-dfdx(i,j),0.0) + &
               (0.00001/(max(0.00001, &
               (this%xmax(j) - this%xmin(j))))))
          this%qij(i,j) = (x(j) - this%low(j))**2 * &
               (0.001*max(dfdx(i, j), 0.0) + &
               1.001*max(-dfdx(i, j), 0.0) + &
               (0.00001/(max(0.00001, &
               (this%xmax(j) - this%xmin(j))))))
       end do
    end do

    !computing bi as defined in page 5
    this%bi = 0_rp
    do i = 1, this%m
       !MPI: here this%n is the global n
       do j = 1, this%n
          this%bi(i) = this%bi(i) + &
               this%pij(i,j) / (this%upp(j) - x(j)) + &
               this%qij(i,j) / (x(j) - this%low(j))
       end do
    end do


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!Showing that for double precision, bi will be different when!!!!!!!!
    !!!!!!!!!!!computed in parallel compare to sequential!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! this%bi = 0_rp
    ! longbi = 0.0
    ! do i = 1, this%m
    !     !MPI: here this%n is the global n
    !     do j = 1, this%n
    !         this%bi(i) = this%bi(i) + &
    !                     this%pij(i,j)/ (this%upp(j) - x(j)) + &
    !                     this%qij(i,j)/(x(j) - this%low(j))
    !         longbi(i) = longbi(i) + &
    !                     this%pij(i,j)/ (this%upp(j) - x(j)) + &
    !                     this%qij(i,j)/(x(j) - this%low(j))
    !     end do
    ! end do
    ! print *, "bi =  ", this%bi, "this%n = ", this%n
    ! print *, "longbi =  ", longbi
    ! ierr = 2160
    ! longbi = 0.0
    ! this%bi = 0_rp
    ! do i = 1, this%m
    !     do j = 1, ierr
    !         this%bi(i) = this%bi(i) + &
    !                     this%pij(i,j)/ (this%upp(j) - x(j)) + &
    !                     this%qij(i,j)/(x(j) - this%low(j))
    !         longbi(i) = longbi(i) + &
    !                     this%pij(i,j)/ (this%upp(j) - x(j)) + &
    !                     this%qij(i,j)/(x(j) - this%low(j))
    !     end do
    ! end do
    ! print *, "bi =  ", this%bi, "first batch(1-ierr)"
    ! print *, "longbi =  ", longbi, "first batch(1-ierr)"
    ! longbiglobal = longbi
    ! longbi = 0.0
    ! globaltmp_m = this%bi
    ! this%bi = 0_rp
    ! do i = 1, this%m
    !     do j = ierr+1, this%n
    !         this%bi(i) = this%bi(i) + &
    !                     this%pij(i,j)/ (this%upp(j) - x(j)) + &
    !                     this%qij(i,j)/(x(j) - this%low(j))
    !         longbi(i) = longbi(i) + &
    !                     this%pij(i,j)/ (this%upp(j) - x(j)) + &
    !                     this%qij(i,j)/(x(j) - this%low(j))
    !     end do
    ! end do
    ! print *, "bi =  ", this%bi, "second batch(ierr+1:end)"
    ! print *, "longbi =  ", longbi, "second batch(ierr+1:end)"
    ! print *, "bi =  ", this%bi+globaltmp_m, "first + second"
    ! print *, "longbi =  ", longbi+longbiglobal, "first + second"
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    globaltmp_m = 0.0_rp
    call MPI_Allreduce(this%bi, globaltmp_m, this%m, &
         mpi_real_precision, mpi_sum, neko_comm, ierr)
    this%bi = globaltmp_m - fval

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
         dely, dellambda, &
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

    ! intial value for the parameters in the subsolve based on
    ! page 15 of "https://people.kth.se/~krille/mmagcmma.pdf"
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

    real(kind=rp), dimension(this%m), intent(in) :: fval(:)
    real(kind=rp), dimension(this%n), intent(in) :: df0dx(:)
    real(kind=rp), dimension(this%m, this%n), intent(in) :: dfdx(:,:)

    real(kind=rp) :: rez, rezeta
    real(kind=rp), dimension(this%m) :: rey, relambda, remu, res
    real(kind=rp), dimension(this%n) :: rex, rexsi, reeta
    real(kind=rp), dimension(3*this%n+4*this%m+2) :: residu

    real(kind=rp), dimension(4*this%m+2) :: residu_small
    integer :: ierr
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

    call MPI_Allreduce(maxval(abs(residu)), this%residumax, 1, &
         mpi_real_precision, mpi_max, neko_comm, ierr)

    call MPI_Allreduce(norm2(rex)**2 + norm2(rexsi)**2 + norm2(reeta)**2, &
         re_xstuff_squ_global, 1, mpi_real_precision, mpi_sum, neko_comm, ierr)

    residu_small = [rey, rez, relambda, remu, rezeta, res]

    this%residunorm = sqrt(norm2(residu_small)**2 + re_xstuff_squ_global)

  end subroutine mma_KKT_cpu
end submodule mma_cpu
