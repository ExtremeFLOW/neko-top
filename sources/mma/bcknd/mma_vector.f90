submodule (mma) mma_vector

contains
  module subroutine mma_gensub_vector(this, iter, x, df0dx, fval, dfdx)
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

    call globaltmp_m%init(this%m)

    if (iter .lt. 3) then
       this%low = x - this%asyinit * (this%xmax - this%xmin)
       this%upp = x + this%asyinit * (this%xmax - this%xmin)
    else
       !Move asymptotes low and upp
       ! Todo: Port to vectorized operations
       do j = 1, this%n
          if ((x%x(j) - this%xold1%x(j))*(this%xold1%x(j) - this%xold2%x(j)) &
               .lt. 0) then
             this%low%x(j) = x%x(j) - &
                  this%asydecr * (this%xold1%x(j) - this%low%x(j))
             this%upp%x(j) = x%x(j) + &
                  this%asydecr * (this%upp%x(j) - this%xold1%x(j))

          else if ((x%x(j) - this%xold1%x(j))* &
               (this%xold1%x(j) - this%xold2%x(j)) .gt. 0) then
             this%low%x(j) = x%x(j) - &
                  this%asyincr * (this%xold1%x(j) - this%low%x(j))
             this%upp%x(j) = x%x(j) + &
                  this%asyincr * (this%upp%x(j) - this%xold1%x(j))
          else
             this%low%x(j) = x%x(j) - (this%xold1%x(j) - this%low%x(j))
             this%upp%x(j) = x%x(j) + (this%upp%x(j) - this%xold1%x(j))
          end if

          ! setting a minimum and maximum for the low and upp
          ! asymptotes (eq3.9)
          this%low%x(j) = max(this%low%x(j), &
               x%x(j) - 10*(this%xmax%x(j) - this%xmin%x(j)))
          this%low%x(j) = min(this%low%x(j), &
               x%x(j) - 0.01*(this%xmax%x(j) - this%xmin%x(j)))

          this%upp%x(j) = min(this%upp%x(j), &
               x%x(j) + 10*(this%xmax%x(j) - this%xmin%x(j)))
          this%upp%x(j) = max(this%upp%x(j), &
               x%x(j) + 0.01*(this%xmax%x(j) - this%xmin%x(j)))
       end do
    end if
    ! we can move alpha and beta out of the following loop if needed as:
    ! this%alpha = max(this%xmin, this%low + &
    !     0.1*(this%x- this%low), this- 0.5*(this%xmax - this%xmin))
    ! this%beta = min(this%xmax, this%upp -  &
    !     0.1*(this%upp - this), this+ 0.5*(this%xmax - this%xmin))
    ! Todo: Port to vectorized operations
    do j = 1, this%n
       ! set the the bounds and coefficients for the approximation
       ! the move bounds (alpha and beta )are slightly more restrictive
       ! than low and upp. This is done based on eq(3.6)--eq(3.10).
       ! also check
       ! https://comsolyar.com/wp-content/uploads/2020/03/gcmma.pdf
       ! eq (2.8) and (2.9)
       this%alpha%x(j) = max(this%xmin%x(j), this%low%x(j) + &
            0.1*(x%x(j)- this%low%x(j)), &
            x%x(j) - 0.5*(this%xmax%x(j) - this%xmin%x(j)))
       this%beta%x(j) = min(this%xmax%x(j), this%upp%x(j) - &
            0.1*(this%upp%x(j) - x%x(j)), &
            x%x(j) + 0.5*(this%xmax%x(j) - this%xmin%x(j)))

       !Calculate p0j, q0j, pij, qij
       !where j = 1,2,...,n and i = 1,2,...,m  (eq(2.3)-eq(2.5))
       this%p0j%x(j) = (this%upp%x(j) - x%x(j))**2 * &
            (1.001_rp*max(df0dx%x(j), 0.0_rp) + &
            0.001_rp*max(-df0dx%x(j), 0.0_rp) + &
            (0.00001_rp/(max(0.00001_rp, &
            (this%xmax%x(j) - this%xmin%x(j))))))

       this%q0j%x(j) = (x%x(j) - this%low%x(j))**2 * &
            (0.001_rp*max(df0dx%x(j),0.0_rp) + &
            1.001_rp*max(-df0dx%x(j),0.0_rp) + &
            (0.00001_rp/(max(0.00001_rp, &
            (this%xmax%x(j) - this%xmin%x(j))))))

       do i = 1, this%m
          this%pij%x(i,j) = (this%upp%x(j) - x%x(j))**2 * &
               (1.001_rp*max(dfdx%x(i, j), 0.0_rp) + &
               0.001_rp*max(-dfdx%x(i, j), 0.0_rp) + &
               (0.00001_rp/(max(0.00001_rp, &
               (this%xmax%x(j) - this%xmin%x(j))))))
          this%qij%x(i,j) = (x%x(j) - this%low%x(j))**2 * &
               (0.001_rp*max(dfdx%x(i, j), 0.0_rp) + &
               1.001_rp*max(-dfdx%x(i, j), 0.0_rp) + &
               (0.00001_rp/(max(0.00001_rp, &
               (this%xmax%x(j) - this%xmin%x(j))))))
       end do
    end do

    !computing bi as defined in page 5
    ! Todo: Port to vectorized operations
    this%bi = 0.0_rp
    do i = 1, this%m
       !MPI: here this%n is the global n
       do j = 1, this%n
          this%bi%x(i) = this%bi%x(i) + &
               this%pij%x(i,j) / (this%upp%x(j) - x%x(j)) + &
               this%qij%x(i,j) / (x%x(j) - this%low%x(j))
       end do
    end do


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!Showing that for double precision, bi will be different when!!!!!!!!
    !!!!!!!!!!!computed in parallel compare to sequential!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! this%bi= 0_rp
    ! longbi = 0.0
    ! do i = 1, this%m
    !     !MPI: here this%n is the global n
    !     do j = 1, this%n
    !         this%bi(i) = this%bi(i) + &
    !                     this%pij%x(i,j)/ (this%upp%x(j) - x(j)) + &
    !                     this%qij%x(i,j)/(x(j) - this%low%x(j))
    !         longbi(i) = longbi(i) + &
    !                     this%pij%x(i,j)/ (this%upp%x(j) - x(j)) + &
    !                     this%qij%x(i,j)/(x(j) - this%low%x(j))
    !     end do
    ! end do
    ! print *, "bi =  ", this%bi%x, "this%n = ", this%n
    ! print *, "longbi =  ", longbi
    ! ierr = 2160
    ! longbi = 0.0
    ! this%bi= 0_rp
    ! do i = 1, this%m
    !     do j = 1, ierr
    !         this%bi(i) = this%bi(i) + &
    !                     this%pij%x(i,j)/ (this%upp%x(j) - x(j)) + &
    !                     this%qij%x(i,j)/(x(j) - this%low%x(j))
    !         longbi(i) = longbi(i) + &
    !                     this%pij%x(i,j)/ (this%upp%x(j) - x(j)) + &
    !                     this%qij%x(i,j)/(x(j) - this%low%x(j))
    !     end do
    ! end do
    ! print *, "bi =  ", this%bi%x, "first batch(1-ierr)"
    ! print *, "longbi =  ", longbi, "first batch(1-ierr)"
    ! longbiglobal = longbi
    ! longbi = 0.0
    ! globaltmp_m= this%bi
    ! this%bi= 0_rp
    ! do i = 1, this%m
    !     do j = ierr+1, this%n
    !         this%bi(i) = this%bi(i) + &
    !                     this%pij%x(i,j)/ (this%upp%x(j) - x(j)) + &
    !                     this%qij%x(i,j)/(x(j) - this%low%x(j))
    !         longbi(i) = longbi(i) + &
    !                     this%pij%x(i,j)/ (this%upp%x(j) - x(j)) + &
    !                     this%qij%x(i,j)/(x(j) - this%low%x(j))
    !     end do
    ! end do
    ! print *, "bi =  ", this%bi%x, "second batch(ierr+1:end)"
    ! print *, "longbi =  ", longbi, "second batch(ierr+1:end)"
    ! print *, "bi =  ", this%bi+globaltmp_m, "first + second"
    ! print *, "longbi =  ", longbi+longbiglobal, "first + second"
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    globaltmp_m = 0.0_rp
    call MPI_Allreduce(this%bi, globaltmp_m, this%m, &
         mpi_real_precision, mpi_sum, neko_comm, ierr)
    this%bi= globaltmp_m - fval

  end subroutine mma_gensub_vector

  subroutine mma_subsolve_dpip_vector(this, designx)
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
    type(vector_t), intent(inout) :: designx
    !Note that there is a local dummy "x" in this subroutine, thus, we call
    !the current design "designx" instead of just "x"
    integer :: i, j, k, iter, ggdumiter, itto, ierr
    real(kind=rp) :: epsi, residumax, residunorm, &
         z, zeta, rez, rezeta, &
         delz, dz, dzeta, &
         steg, dummy_one, zold, zetaold, newresidu
    type(vector_t) :: y, lambda, s, mu, &
         rey, relambda, remu, res, &
         dely, dellambda, &
         dy, dlambda, ds, dmu, &
         yold, lambdaold, sold, muold, &
         globaltmp_m
    type(vector_t) :: x, xsi, eta, &
         rex, rexsi, reeta, &
         delx, diagx, dx, dxsi, deta, &
         xold, xsiold, etaold
    type(vector_t) :: residu
    type(vector_t) :: residu_small
    type(vector_t) :: xx, dxx

    type(vector_t) :: bb
    type(matrix_t) :: GG
    type(matrix_t) :: AA
    type(matrix_t) :: globaltmp_mm

    ! using DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO ) in lapack to solve
    ! the linear system which needs the following parameters
    integer :: info
    integer, dimension(this%m+1) :: ipiv

    real(kind=rp) :: re_xstuff_squ_global

    integer :: nglobal

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
    call globaltmp_m%init(this%m)

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

    call residu%init(3*this%n+4*this%m+2)
    call residu_small%init(4*this%m+2)
    call xx%init(2*this%n+4*this%m+2)
    call dxx%init(2*this%n+4*this%m+2)

    call bb%init(this%m+1)

    call GG%init(this%m, this%n)
    call AA%init(this%m+1, this%m+1)
    call globaltmp_mm%init(this%m, this%m)

    ! intial value for the parameters in the subsolve based on
    ! page 15 of "https://people.kth.se/~krille/mmagcmma.pdf"
    dummy_one = 1
    epsi = 1 !100
    x = 0.5_rp * (this%alpha + this%beta)
    y = 1.0_rp
    z = 1.0_rp
    zeta = 1.0_rp
    lambda = 1.0_rp
    s = 1.0_rp
    xsi = max(1.0_rp, 1.0_rp / (x - this%alpha))
    eta = max(1.0_rp, 1.0_rp / (this%beta - x))
    mu = max(1.0_rp, 0.5_rp * this%c)

    do while (epsi .gt. 0.9*this%epsimin)
       ! calculating residuals based on
       ! "https://people.kth.se/~krille/mmagcmma.pdf" for the variables
       ! x, y, z, lambda residuals based on eq(5.9a)-(5.9d), respectively.
       rex = ((this%p0j + matmul(transpose(this%pij), &
            lambda))/(this%upp - x)**2 - &
            (this%q0j + matmul(transpose(this%qij), &
            lambda))/(x - this%low)**2 ) - &
            xsi + eta

       call MPI_Allreduce(this%n, nglobal, 1, &
            MPI_INTEGER, mpi_sum, neko_comm, ierr)

       !!!! computing without matmul and transpose
       ! rex = 0.0_rp
       ! do j = 1, this%n
       !     do i = 1, this%m
       !         rex(j) = rex(j) + this%pij%x(i,j) * &
       !             lambda(i)/(this%upp%x(j) - x(j))**2 - &
       !             this%qij%x(i,j) * lambda(i)/(x(j) - this%low%x(j))**2
       !     end do
       !     rex(j) = rex(j) + this%p0j%x(j)/(this%upp%x(j) - x(j))**2 &
       !                     - this%q0j%x(j)/(x(j) - this%low%x(j))**2 &
       !                     - xsi(j)  + eta(j)
       ! end do


       rey = this%c + this%d*y - lambda - mu
       rez = this%a0 - zeta - dot_product(lambda, this%a)

       ! relambda = matmul(this%pij,1.0/(this%upp - x)) + &
       !         matmul(this%qij, 1.0/(x - this%low)) - &
       !         this%a*z - y + s - this%bi
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
       relambda = globaltmp_m - this%a*z - y + s - this%bi


       rexsi = xsi*(x - this%alpha) - epsi
       reeta = eta*(this%beta - x) - epsi
       remu = mu*y - epsi
       rezeta = zeta*z - epsi
       res = lambda*s - epsi

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

          delx = 0.0_rp
          do j = 1, this%n
             do i = 1, this%m
                delx(j) = delx(j) + this%pij%x(i,j) * &
                     lambda(i)/(this%upp%x(j) - x(j))**2 &
                     - this%qij%x(i,j) * lambda(i)/(x(j) - this%low%x(j))**2
             end do
             delx(j) = delx(j) + this%p0j%x(j)/(this%upp%x(j) - x(j))**2 &
                  - this%q0j%x(j)/(x(j) - this%low%x(j))**2 &
                  - epsi/(x(j) - this%alpha%x(j)) &
                  + epsi/(this%beta%x(j) - x(j))
          end do
          dely = this%c+ this%d*y - lambda - epsi/y
          delz = this%a0 - dot_product(lambda, this%a) - epsi/z

          dellambda = 0.0_rp
          do i = 1, this%m
             do j = 1, this%n !this n is global
                ! Accumulate sums for dellambda (the term gi(x))
                dellambda(i) = dellambda(i) + &
                     this%pij%x(i,j)/(this%upp%x(j) - x(j)) &
                     + this%qij%x(i,j)/(x(j) - this%low%x(j))
             end do
          end do

          globaltmp_m = 0.0_rp
          call MPI_Allreduce(dellambda, globaltmp_m, this%m, &
               mpi_real_precision, mpi_sum, neko_comm, ierr)

          dellambda = globaltmp_m - this%a*z - y - this%bi + epsi / lambda

          ! delx = ((this%p0j + matmul(transpose(this%pij), &
          !     lambda))/(this%upp - x)**2 - &
          !     (this%q0j + matmul(transpose(this%qij), &
          !     lambda))/(x - this%low)**2 ) - &
          !     epsi/(x - this%alpha) + epsi/(this%beta - x)

          ! dely =  this%c + this%d*y - lambda - epsi/y
          ! delz = this%a0 - dot_product(lambda, this%a) - epsi/z
          ! dellambda = matmul(this%pij,1.0/(this%upp - x))+&
          !     matmul(this%qij, 1.0/(x - this%low)) - &
          !     this%a*z - y - this%bi + epsi/lambda

          do ggdumiter = 1, this%m
             GG(ggdumiter, :) = this%pij%x(ggdumiter,:)/ &
                  (this%upp - x)**2 - &
                  this%qij%x(ggdumiter,:)/(x - this%low)**2
          end do

          diagx = ((this%p0j + matmul(transpose(this%pij), &
               lambda))/(this%upp - x)**3 + &
               (this%q0j + matmul(transpose(this%qij), &
               lambda))/(x - this%low)**3 )
          diagx = 2*diagx + xsi/(x - this%alpha) + &
               eta/(this%beta- x)


          !Here we only consider the case m<n in the matlab code
          !assembling the right hand side matrix based on eq(5.20)
          ! bb = [dellambda + dely/(this%d + &
          !         (mu/y)) - matmul(GG,delx/diagx), delz ]
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

          bb(1:this%m) = dellambda + dely/(this%d+ (mu/y)) - bb(1:this%m)
          bb(this%m +1) = delz
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! !assembling the coefficients matrix AA based on eq(5.20)
          ! AA(1:this%m,1:this%m) =  &
          ! matmul(matmul(GG,mma_diag(1/diagx)), transpose(GG))
          ! !update diag(AA)
          ! AA(1:this%m,1:this%m) = AA(1:this%m,1:this%m) + &
          !     mma_diag(s/lambda + 1.0/(this%d + (mu/y)))

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
                  1.0_rp / (this%d%x(i) + mu(i) / y(i)))
          end do

          AA(1:this%m, this%m+1) = this%a
          AA(this%m+1, 1:this%m) = this%a
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
          dx = -delx / diagx - matmul(transpose(GG), dlambda)/diagx

          dy = (-dely + dlambda) / (this%d + (mu / y))
          dxsi = -xsi + (epsi - dx*xsi) / (x - this%alpha)
          deta = -eta + (epsi + dx*eta) / (this%beta - x)
          dmu = -mu + (epsi - mu*dy) / y
          dzeta = -zeta + (epsi-zeta*dz)/z
          ds = -s + (epsi - dlambda*s) / lambda

          !2*this%n+4*this%m+2
          dxx = [dy, dz, dlambda, dxsi, deta, dmu, dzeta, ds]
          xx = [y, z, lambda, xsi, eta, mu, zeta, s]
          steg = maxval([dummy_one, -1.01*dxx/xx, -1.01*dx/ &
               (x - this%alpha), 1.01 * dx / (this%beta - x)])
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
             rex = ((this%p0j + matmul(transpose(this%pij), &
                  lambda))/(this%upp - x)**2 - &
                  (this%q0j + matmul(transpose(this%qij), &
                  lambda))/(x - this%low)**2 ) - &
                  xsi + eta
             rey = this%c + this%d*y - lambda - mu
             rez = this%a0 - zeta - dot_product(lambda, this%a)
             ! relambda = matmul(this%pij,1.0/&
             !         (this%upp - x)) + matmul(this%qij, &
             !         1.0/(x - this%low)) - this%a*z - &
             !         y + s - this%bi
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


             relambda = relambda - this%a*z - y + s - this%bi

             rexsi = xsi*(x - this%alpha) - epsi
             reeta = eta*(this%beta - x) - epsi
             remu = mu*y - epsi
             rezeta = zeta*z - epsi
             res = lambda*s - epsi

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
    this%xold1 = designx
    designx= x

    !update the parameters of the MMA object nesessary to compute KKT residu
    this%y= y
    this%z = z
    this%lambda= lambda
    this%zeta = zeta
    this%xsi= xsi
    this%eta= eta
    this%mu= mu
    this%s= s

  end subroutine mma_subsolve_dpip_vector

  subroutine mma_KKT_vector(this, x, df0dx, fval, dfdx)
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
    type(vector_t) :: residu

    type(vector_t) :: residu_small
    integer :: ierr
    real(kind=rp) :: re_xstuff_squ_global

    call rey%init(this%m)
    call relambda%init(this%m)
    call remu%init(this%m)
    call res%init(this%m)

    call rex%init(this%n)
    call rexsi%init(this%n)
    call reeta%init(this%n)

    call residu%init(3*this%n+4*this%m+2)
    call residu_small%init(4*this%m+2)


    rex = df0dx+ matmul(transpose(dfdx), this%lambda) - this%xsi + &
         this%eta
    rey = this%c + this%d*this%y - this%lambda - this%mu
    rez = this%a0 - this%zeta - dot_product(this%lambda, this%a)

    relambda = fval- this%a*this%z - this%y + this%s
    rexsi = this%xsi*(x - this%xmin)
    reeta = this%eta*(this%xmax - x)
    remu = this%mu*this%y
    rezeta = this%zeta*this%z
    res = this%lambda*this%s

    residu = [rex, rey, rez, relambda, rexsi, reeta, remu, rezeta, res]

    call MPI_Allreduce(maxval(abs(residu)), this%residumax, 1, &
         mpi_real_precision, mpi_max, neko_comm, ierr)

    call MPI_Allreduce(norm2(rex)**2 + norm2(rexsi)**2 + norm2(reeta)**2, &
         re_xstuff_squ_global, 1, mpi_real_precision, mpi_sum, neko_comm, ierr)

    residu_small = [rey, rez, relambda, remu, rezeta, res]

    this%residunorm = sqrt(norm2(residu_small)**2 + re_xstuff_squ_global)

  end subroutine mma_KKT_vector
end submodule mma_vector
