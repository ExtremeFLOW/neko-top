submodule (mma) mma_vector
  use math, only: glsum, vlsc2, glsc2
  use device_math, only: device_glsum
  use utils, only: neko_error
  use comm, only: mpi_real_precision, neko_comm
  use mpi_f08, only: MPI_INTEGER, MPI_IN_PLACE, MPI_Allreduce, mpi_min, mpi_sum

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
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call neko_error("Device if loop is not implemented")
       else
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
    end if

    ! we can move alpha and beta out of the following loop if needed as:
    ! this%alpha = max(this%xmin, this%low + &
    !     0.1*(this%x- this%low), this- 0.5*(this%xmax - this%xmin))
    ! this%beta = min(this%xmax, this%upp -  &
    !     0.1*(this%upp - this), this + 0.5*(this%xmax - this%xmin))
    ! Todo: Port to vectorized operations
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_error("Device if loop is not implemented")
    else
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
    end if

    !computing bi as defined in page 5
    this%bi = 0.0_rp
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_error("Device computation of bi is not implemented")
    else
       do i = 1, this%m
          !MPI: here this%n is the global n
          do j = 1, this%n
             this%bi%x(i) = this%bi%x(i) + &
                  this%pij%x(i,j) / (this%upp%x(j) - x%x(j)) + &
                  this%qij%x(i,j) / (x%x(j) - this%low%x(j))
          end do
       end do
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       this%bi = device_glsum(this%bi%x_d, this%m) - fval
    else
       this%bi = glsum(this%bi%x, this%m) - fval
    end if

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
    real(kind = rp) :: epsi, residumax, residunorm, &
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
    integer, dimension(this%m + 1) :: ipiv

    real(kind = rp) :: re_xstuff_squ_global

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

    call residu%init(3*this%n + 4*this%m + 2)
    call residu_small%init(4*this%m + 2)
    call xx%init(2*this%n + 4*this%m + 2)
    call dxx%init(2*this%n + 4*this%m + 2)

    call bb%init(this%m + 1)

    call GG%init(this%m, this%n)
    call AA%init(this%m + 1, this%m + 1)
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
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_error("Device `max` is not implemented")
    else
       xsi%x = max(1.0_rp, 1.0_rp / (x%x - this%alpha%x))
       eta%x = max(1.0_rp, 1.0_rp / (this%beta%x - x%x))
       mu%x = max(1.0_rp, 0.5_rp * this%c%x)
    end if

    do while (epsi .gt. 0.9*this%epsimin)
       ! calculating residuals based on
       ! "https://people.kth.se/~krille/mmagcmma.pdf" for the variables
       ! x, y, z, lambda residuals based on eq(5.9a)-(5.9d), respectively.
       rex = (this%p0j + this%pij%transpose()*lambda) / &
            (this%upp - x)**2 - &
            (this%q0j + this%qij%transpose()*lambda) / &
            (x - this%low)**2 - xsi + eta

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
       if (NEKO_BCKND_DEVICE .eq. 1) then
          rez = this%a0 - zeta - device_vlsc2(lambda%x_d, this%a%x_d, this%m)
       else
          rez = this%a0 - zeta - vlsc2(lambda%x, this%a%x, this%m)
       end if

       relambda = this%pij * (1.0_rp / (this%upp - x)) + &
            this%qij * (1.0_rp / (x - this%low)) - &
            this%a*z - y + s - this%bi

       if (NEKO_BCKND_DEVICE .eq. 1) then
          relambda = device_glsum(relambda%x_d, this%m) - &
               this%a*z - y + s - this%bi
       else
          relambda = glsum(relambda%x, this%m) - &
               this%a*z - y + s - this%bi
       end if

       rexsi = xsi*(x - this%alpha) - epsi
       reeta = eta*(this%beta - x) - epsi
       remu = mu*y - epsi
       rezeta = zeta*z - epsi
       res = lambda*s - epsi

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call neko_error("Device Array construction is not implemented")
       else
          residu%x = [rex%x, rey%x, rez, relambda%x, rexsi%x, reeta%x, remu%x, &
               rezeta, res%x]
       end if
       residumax = 0_rp

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call neko_error("Device `abs` is not implemented")
       else
          residumax = glmax(abs(residu%x), 3*this%n + 4*this%m + 2)
       end if

       if (NEKO_BCKND_DEVICE .eq. 1) then
          re_xstuff_squ_global = device_glsc2(rex%x_d, rex%x_d, this%n) + &
               device_glsc2(rexsi%x_d, rexsi%x_d, this%n) + &
               device_glsc2(reeta%x_d, reeta%x_d, this%n)
       else
          re_xstuff_squ_global = glsc2(rex%x, rex%x, this%n) + &
               glsc2(rexsi%x, rexsi%x, this%n) + &
               glsc2(reeta%x, reeta%x, this%n)
       end if

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call neko_error("Device Array construction is not implemented")
       else
          residu_small%x = [rey%x, rez, relambda%x, remu%x, rezeta, res%x]
       end if

       if (NEKO_BCKND_DEVICE .eq. 1) then
          residunorm = sqrt( &
               device_glsc2(residu_small%x_d, residu_small%x_d, this%m) + &
               re_xstuff_squ_global)
       else
          residunorm = sqrt(glsc2(residu_small%x, residu_small%x, this%m) + &
               re_xstuff_squ_global)
       end if


       do iter = 1, this%max_iter !ittt
          !Check the condition
          if (residumax .lt. epsi) exit

          delx = 0.0_rp
          if (NEKO_BCKND_DEVICE .eq. 1) then
             call neko_error("Device if loop is not implemented")
          else
             do j = 1, this%n
                do i = 1, this%m
                   delx%x(j) = delx%x(j) + this%pij%x(i,j) * &
                        lambda%x(i)/(this%upp%x(j) - x%x(j))**2 &
                        - this%qij%x(i,j) * lambda%x(i)/(x%x(j) - this%low%x(j))**2
                end do
                delx%x(j) = delx%x(j) + this%p0j%x(j)/(this%upp%x(j) - x%x(j))**2 &
                     - this%q0j%x(j)/(x%x(j) - this%low%x(j))**2 &
                     - epsi/(x%x(j) - this%alpha%x(j)) &
                     + epsi/(this%beta%x(j) - x%x(j))
             end do
          end if
          dely = this%c + this%d*y - lambda - epsi/y

          if (NEKO_BCKND_DEVICE .eq. 1) then
             delz = this%a0 - device_vlsc2(lambda%x_d, this%a%x_d, this%m) - epsi/z
          else
             delz = this%a0 - vlsc2(lambda%x, this%a%x, this%m) - epsi/z
          end if

          dellambda = 0.0_rp
          if (NEKO_BCKND_DEVICE .eq. 1) then
             call neko_error("Device if loop is not implemented")
          else
             do i = 1, this%m
                do j = 1, this%n !this n is global
                   ! Accumulate sums for dellambda (the term gi(x))
                   dellambda%x(i) = dellambda%x(i) + &
                        this%pij%x(i,j)/(this%upp%x(j) - x%x(j)) &
                        + this%qij%x(i,j)/(x%x(j) - this%low%x(j))
                end do
             end do
          end if

          globaltmp_m = 0.0_rp
          if (NEKO_BCKND_DEVICE .eq. 1) then
             call neko_error("Device all reduce of array is not implemented")
          else
             call MPI_Allreduce(dellambda%x, globaltmp_m%x, this%m, &
                  mpi_real_precision, mpi_sum, neko_comm, ierr)
          end if

          dellambda = globaltmp_m - this%a*z - y - this%bi + epsi / lambda

          ! delx = ((this%p0j + matmul(transpose(this%pij), &
          !     lambda))/(this%upp - x)**2 - &
          !     (this%q0j + matmul(transpose(this%qij), &
          !     lambda))/(x - this%low)**2 ) - &
          !     epsi/(x - this%alpha) + epsi/(this%beta - x)

          ! dely =  this%c + this%d*y - lambda - epsi/y
          ! delz = this%a0 - dot_product(lambda, this%a) - epsi/z
          ! dellambda = matmul(this%pij,1.0_rp/(this%upp - x)) + &
          !     matmul(this%qij, 1.0_rp/(x - this%low)) - &
          !     this%a*z - y - this%bi + epsi/lambda

          if (NEKO_BCKND_DEVICE .eq. 1) then
             call neko_error("Device GG Array construction is not implemented")
          else
             GG = 0.0_rp
             do ggdumiter = 1, this%m
                GG%x(ggdumiter, :) = this%pij%x(ggdumiter,:)/ &
                     (this%upp%x - x%x)**2 - &
                     this%qij%x(ggdumiter,:)/(x%x - this%low%x)**2
             end do
          end if

          diagx = (this%p0j + this%pij%transpose() * lambda) / &
               (this%upp - x)**3 + &
               (this%q0j + this%qij%transpose() * lambda) / &
               (x - this%low)**3
          diagx = 2.0_rp*diagx + xsi/(x - this%alpha) + &
               eta/(this%beta - x)


          !Here we only consider the case m<n in the matlab code
          !assembling the right hand side matrix based on eq(5.20)
          ! bb = [dellambda + dely/(this%d + &
          !         (mu/y)) - matmul(GG,delx/diagx), delz ]
          !!!!!!!!!!!!!!for MPI computation of bb!!!!!!!!!!!!!!!!!!!!!!!!!

          if (NEKO_BCKND_DEVICE .eq. 1) then
             call neko_error("Device BB construction is not implemented")
          else
             bb = 0.0_rp
             do i = 1, this%m
                do j = 1, this%n ! this n is global
                   bb%x(i) = bb%x(i) + GG%x(i, j) * (delx%x(j) / diagx%x(j))
                end do
             end do
             globaltmp_m = 0.0_rp
             call MPI_Allreduce(bb%x(1:this%m), globaltmp_m%x, this%m, &
                  mpi_real_precision, mpi_sum, neko_comm, ierr)
             bb%x(1:this%m) = globaltmp_m%x

             bb%x(1:this%m) = dellambda%x + dely%x/(this%d%x + (mu%x/y%x)) - bb%x(1:this%m)
             bb%x(this%m + 1) = delz
          end if
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! !assembling the coefficients matrix AA based on eq(5.20)
          ! AA(1:this%m,1:this%m) =  &
          ! matmul(matmul(GG,mma_diag(1/diagx)), transpose(GG))
          ! !update diag(AA)
          ! AA(1:this%m,1:this%m) = AA(1:this%m,1:this%m) + &
          !     mma_diag(s/lambda + 1.0_rp/(this%d + (mu/y)))

          if (NEKO_BCKND_DEVICE .eq. 1) then
             call neko_error("Device if loop is not implemented")
          else
             AA = 0.0_rp
             !Direct computation of the matrix multiplication
             !(for better performance)
             do i = 1, this%m
                do j = 1, this%m
                   ! Compute the (i, j) element of AA
                   do k = 1, this%n !this n is global
                      AA%x(i, j) = AA%x(i, j) + GG%x(i, k) * &
                           (1.0_rp / diagx%x(k)) * GG%x(j, k)
                   end do
                end do
             end do

             call MPI_Allreduce(MPI_IN_PLACE, AA%x(1:this%m, 1:this%m), &
                  this%m*this%m, mpi_real_precision, mpi_sum, neko_comm, ierr)

             do i = 1, this%m
                !update the diag AA
                AA%x(i, i) = AA%x(i, i) + (s%x(i) / lambda%x(i) + &
                     1.0_rp / (this%d%x(i) + mu%x(i) / y%x(i)))
             end do

             AA%x(1:this%m, this%m + 1) = this%a%x
             AA%x(this%m + 1, 1:this%m) = this%a%x
             AA%x(this%m + 1, this%m + 1) = -zeta/z
          end if


          if (NEKO_BCKND_DEVICE .eq. 1) then
             call neko_error("Device DGESV is not implemented")
          else
             call DGESV(this%m + 1, 1, AA%x, this%m + 1, ipiv, bb%x, this%m + 1, info)
             ! if info! = 0 then DGESV is failed.
             if (info .ne. 0) then
                write(stderr, *) "DGESV failed to solve the linear system in MMA."
                write(stderr, *) "Please check mma_subsolve_dpip in mma.f90"
                error stop
             end if
          end if

          if (NEKO_BCKND_DEVICE .eq. 1) then
             call neko_error("Device Index Array construction is not implemented")
          else
             dlambda%x = bb%x(1:this%m)
             dz = bb%x(this%m + 1)
          end if
          ! based on eq(5.19)
          dx = -delx / diagx - (GG%transpose() * dlambda)/diagx

          dy = (-dely + dlambda) / (this%d + (mu / y))
          dxsi = -xsi + (epsi - dx*xsi) / (x - this%alpha)
          deta = -eta + (epsi + dx*eta) / (this%beta - x)
          dmu = -mu + (epsi - mu*dy) / y
          dzeta = -zeta + (epsi-zeta*dz)/z
          ds = -s + (epsi - dlambda*s) / lambda

          !2*this%n + 4*this%m + 2
          if (NEKO_BCKND_DEVICE .eq. 1) then
             call neko_error("Device Array construction is not implemented")
          else
             dxx%x = [dy%x, dz, dlambda%x, dxsi%x, deta%x, dmu%x, dzeta, ds%x]
          end if

          if (NEKO_BCKND_DEVICE .eq. 1) then
             call neko_error("Device Array construction is not implemented")
          else
             xx%x = [y%x, z, lambda%x, xsi%x, eta%x, mu%x, zeta, s%x]
          end if

          if (NEKO_BCKND_DEVICE .eq. 1) then
             call neko_error("Device Array construction is not implemented")
          else
             steg = maxval([1.0_rp, -1.01_rp*dxx%x/xx%x, -1.01_rp*dx%x/ &
                  (x%x - this%alpha%x), 1.01_rp * dx%x / (this%beta%x - x%x)])
          end if
          steg = 1.0_rp/steg

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
             rex = ((this%p0j + this%pij%transpose() * lambda) / &
                  (this%upp - x)**2 - &
                  (this%q0j + this%qij%transpose()* lambda) / &
                  (x - this%low)**2 ) - &
                  xsi + eta

             rey = this%c + this%d*y - lambda - mu

             if (NEKO_BCKND_DEVICE .eq. 1) then
                rez = this%a0 - zeta - device_vlsc2(lambda%x_d, &
                     this%a%x_d, this%m)
             else
                rez = this%a0 - zeta - vlsc2(lambda%x, this%a%x, this%m)
             end if

             relambda = this%pij * (1.0_rp/ (this%upp - x)) + &
                  this%qij * ( 1.0_rp/(x - this%low)) - &
                  this%a*z - y + s - this%bi

             if (NEKO_BCKND_DEVICE .eq. 1) then
                call neko_error("Device MPI missing in mma_subsolve_dpip_vector")
             else
                call MPI_Allreduce(MPI_IN_PLACE, relambda%x, this%m, &
                     mpi_real_precision, mpi_sum, neko_comm, ierr)
             end if

             relambda = relambda - this%a*z - y + s - this%bi

             rexsi = xsi*(x - this%alpha) - epsi
             reeta = eta*(this%beta - x) - epsi
             remu = mu*y - epsi
             rezeta = zeta*z - epsi
             res = lambda*s - epsi

             if (NEKO_BCKND_DEVICE .eq. 1) then
                call neko_error("Device Array construction is not implemented")
             else
                residu%x = [rex%x, rey%x, rez, relambda%x, rexsi%x, reeta%x, &
                     remu%x, rezeta, res%x]
             end if

             if (NEKO_BCKND_DEVICE .eq. 1) then
                re_xstuff_squ_global = device_vlsc2(rex%x_d, rex%x_d, this%m) + &
                     device_vlsc2(rexsi%x_d, rexsi%x_d, this%m) + &
                     device_vlsc2(reeta%x_d, reeta%x_d, this%m)
             else
                re_xstuff_squ_global = vlsc2(rex%x, rex%x, this%m) + &
                     vlsc2(rexsi%x, rexsi%x, this%m) + &
                     vlsc2(reeta%x, reeta%x, this%m)
             endif

             call MPI_Allreduce(MPI_IN_PLACE, re_xstuff_squ_global, &
                  1, mpi_real_precision, mpi_sum, neko_comm, ierr)

             if (NEKO_BCKND_DEVICE .eq. 1) then
                call neko_error("Device Array construction is not implemented")
             else
                residu_small%x = [rey%x, rez, relambda%x, remu%x, rezeta, res%x]
             end if

             if (NEKO_BCKND_DEVICE .eq. 1) then
                newresidu = sqrt(device_vlsc2(residu_small%x_d, &
                     residu_small%x_d, this%m) + re_xstuff_squ_global)
             else
                newresidu = sqrt(vlsc2(residu_small%x, residu_small%x, this%m) + &
                     re_xstuff_squ_global)
             end if

             steg = steg/2
          end do

          residunorm = newresidu
          residumax = 0_rp
          if (NEKO_BCKND_DEVICE .eq. 1) then
             call neko_error("Device `abs` is not implemented")
          else
             residumax = glmax(abs(residu%x), 3*this%n + 4*this%m + 2)
          end if

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


    real(kind = rp) :: rez, rezeta
    type(vector_t) :: rey, relambda, remu, res
    type(vector_t) :: rex, rexsi, reeta
    type(vector_t) :: residu

    type(vector_t) :: residu_small
    integer :: ierr
    real(kind = rp) :: re_xstuff_squ_global

    call rey%init(this%m)
    call relambda%init(this%m)
    call remu%init(this%m)
    call res%init(this%m)

    call rex%init(this%n)
    call rexsi%init(this%n)
    call reeta%init(this%n)

    call residu%init(3*this%n + 4*this%m + 2)
    call residu_small%init(4*this%m + 2)


    rex = df0dx + dfdx%transpose() * this%lambda - this%xsi + &
         this%eta
    rey = this%c + this%d*this%y - this%lambda - this%mu

    if (NEKO_BCKND_DEVICE == 1)then
       rez = this%a0 - this%zeta - device_vlsc2(this%lambda%x_d, this%a%x_d, this%m)
    else
       rez = this%a0 - this%zeta - vlsc2(this%lambda%x, this%a%x, this%m)
    end if

    relambda = fval - this%a*this%z - this%y + this%s
    rexsi = this%xsi*(x - this%xmin)
    reeta = this%eta*(this%xmax - x)
    remu = this%mu*this%y
    rezeta = this%zeta*this%z
    res = this%lambda*this%s

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_error("Device Array construction is not implemented")
    else
       residu%x = [rex%x, rey%x, rez, relambda%x, rexsi%x, reeta%x, remu%x, &
            rezeta, res%x]
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call nakeo_error("Device `abs` is not implemented")
    else
       this%residumax = glmax(abs(residu%x), 3*this%n + 4*this%m + 2)
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       re_xstuff_squ_global = device_vlsc2(rex%x_d, rex%x_d, this%m) + &
            device_vlsc2(rexsi%x_d, rexsi%x_d, this%m) + &
            device_vlsc2(reeta%x_d, reeta%x_d, this%m)
    else
       re_xstuff_squ_global = vlsc2(rex%x, rex%x, this%m) + &
            vlsc2(rexsi%x, rexsi%x, this%m) + &
            vlsc2(reeta%x, reeta%x, this%m)
    endif

    call MPI_Allreduce(MPI_IN_PLACE, re_xstuff_squ_global, 1, &
         mpi_real_precision, mpi_sum, neko_comm, ierr)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_error("Device Array construction is not implemented")
    else
       residu_small%x = [rey%x, rez, relambda%x, remu%x, rezeta, res%x]
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       this%residunorm = sqrt(device_vlsc2(residu_small%x_d, &
            residu_small%x_d, this%m) + re_xstuff_squ_global)
    else
       this%residunorm = sqrt(vlsc2(residu_small%x, residu_small%x, this%m) + &
            re_xstuff_squ_global)
    end if

  end subroutine mma_KKT_vector
end submodule mma_vector
