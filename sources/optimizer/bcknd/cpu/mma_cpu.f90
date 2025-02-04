! Copyright (c) 2025, The Neko-TOP Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.

!> Submodule for the CPU implementations related to MMA
submodule (mma) mma_cpu
  use mpi_f08, only: MPI_INTEGER, MPI_REAL, mpi_sum, mpi_min, mpi_max, &
       MPI_Allreduce, MPI_IN_PLACE
  use utils, only: neko_error
  use comm, only: neko_comm, mpi_real_precision

contains

  module subroutine mma_gensub_cpu(this, iter, x, df0dx, fval, dfdx)
    ! ------------------------------------------------------ !
    ! Generate the approximation sub problem by computing    !
    ! the lower and upper asymptotes and the other necessary !
    ! parameters (alpha, beta, p0j, q0j, pij, qij, ...).     !
    ! ------------------------------------------------------ !
    class(mma_t), intent(inout) :: this
    real(kind=rp), dimension(this%n), intent(in) :: x
    real(kind=rp), dimension(this%n), intent(in) :: df0dx
    real(kind=rp), dimension(this%m), intent(in) :: fval
    real(kind=rp), dimension(this%m, this%n), intent(in) :: dfdx
    integer, intent(in) :: iter
    integer :: i, j, ierr
    real(kind=rp) :: asy_factor
    real(kind=rp), dimension(this%n) :: x_diff

    x_diff = this%xmax%x - this%xmin%x

    ! ------------------------------------------------------------------------ !
    ! Setup the current asymptotes

    associate(low => this%low%x, upp => this%upp%x, &
         x_1 => this%xold1%x, x_2 => this%xold2%x)

      if (iter .lt. 3) then
         ! Initialize the lower and upper asymptotes
         low = x - this%asyinit * x_diff
         upp = x + this%asyinit * x_diff
      else
         do j = 1, this%n
            if ((x(j) - x_1(j)) * (x_1(j) - x_2(j)) .lt. 0.0_rp) then
               asy_factor = this%asydecr
            else if ((x(j) - x_1(j)) * (x_1(j) - x_2(j)) .gt. 0.0_rp) then
               asy_factor = this%asyincr
            else
               asy_factor = 1.0_rp
            end if

            low(j) = x(j) - asy_factor * (x_1(j) - low(j))
            upp(j) = x(j) + asy_factor * (upp(j) - x_1(j))
         end do


         ! Setting a minimum and maximum for the low and upp
         ! asymptotes (eq3.9)
         low = max(low, x - 10.0_rp * x_diff)
         low = min(low, x - 0.01_rp * x_diff)

         upp = min(upp, x + 10.0_rp * x_diff)
         upp = max(upp, x + 0.01_rp * x_diff)
      end if

    end associate

    ! ------------------------------------------------------------------------ !
    ! Set the the bounds and coefficients for the approximation
    ! the move bounds (alpha and beta) are slightly more restrictive
    ! than low and upp. This is done based on eq(3.6)--eq(3.10).
    ! also check
    ! https://comsolyar.com/wp-content/uploads/2020/03/gcmma.pdf
    ! eq (2.8) and (2.9)

    associate(alpha => this%alpha%x, beta => this%beta%x, &
         xmin => this%xmin%x, xmax => this%xmax%x, &
         low => this%low%x, upp => this%upp%x)

      alpha = max(xmin, low + 0.1_rp*(x - low), x - 0.5_rp*x_diff)
      beta = min(xmax, upp - 0.1_rp*(upp - x), x + 0.5_rp*x_diff)

    end associate

    ! ------------------------------------------------------------------------ !
    ! Calculate p0j, q0j, pij, qij
    ! where j = 1,2,...,n and i = 1,2,...,m  (eq(2.3)-eq(2.5))

    associate(p0j => this%p0j%x, q0j => this%q0j%x, &
         pij => this%pij%x, qij => this%qij%x, &
         low => this%low%x, upp => this%upp%x)

      p0j = ( &
           1.001_rp * max(df0dx, 0.0_rp) &
           + 0.001_rp * max(-df0dx, 0.0_rp) &
           + 0.00001_rp / max(x_diff, 0.00001_rp) &
           ) * (upp - x)**2

      q0j = ( &
           0.001_rp * max(df0dx, 0.0_rp) &
           + 1.001_rp * max(-df0dx, 0.0_rp) &
           + 0.00001_rp / max(x_diff, 0.00001_rp)&
           ) * (x - low)**2

      do j = 1, this%n
         do i = 1, this%m
            pij(i, j) = ( &
                 1.001_rp * max(dfdx(i, j), 0.0_rp) &
                 + 0.001_rp * max(-dfdx(i, j), 0.0_rp) &
                 + 0.00001_rp / max(x_diff(j), 0.00001_rp) &
                 ) * (upp(j) - x(j))**2

            qij(i, j) = ( &
                 0.001_rp * max(dfdx(i, j), 0.0_rp) &
                 + 1.001_rp * max(-dfdx(i, j), 0.0_rp) &
                 + 0.00001_rp / max(x_diff(j), 0.00001_rp) &
                 ) * (x(j) - low(j))**2
         end do
      end do

    end associate

    ! ------------------------------------------------------------------------ !
    ! Computing bi as defined in page 5

    associate(bi => this%bi%x, &
         pij => this%pij%x, qij => this%qij%x, &
         low => this%low%x, upp => this%upp%x)

      bi = 0.0_rp
      do i = 1, this%m
         do j = 1, this%n
            bi(i) = bi(i) &
                 + pij(i, j) / (upp(j) - x(j)) &
                 + qij(i, j) / (x(j) - low(j))
         end do
      end do

      call MPI_Allreduce(MPI_IN_PLACE, bi, this%m, &
           mpi_real_precision, mpi_sum, neko_comm, ierr)
      bi = bi - fval

    end associate

  end subroutine mma_gensub_cpu

  subroutine mma_subsolve_dpip_cpu(this, designx)
    ! ------------------------------------------------------- !
    ! Dual-primal interior point method using Newton's step   !
    ! to solve MMA sub problem.                               !
    ! A Backtracking Line Search approach is used to compute  !
    ! the step size; starting with the full Newton's step     !
    ! (delta = 1) and dividing by 2 until we have a step size !
    ! that leads to a feasible point while ensuring a         !
    ! decrease in the residue.                                !
    ! ------------------------------------------------------- !
    class(mma_t), intent(inout) :: this
    real(kind=rp), dimension(this%n), intent(inout) :: designx
    !Note that there is a local dummy "x" in this subroutine, thus, we call
    !the current design "designx" instead of just "x"
    integer :: i, j, k, iter, i, itto, ierr
    real(kind=rp) :: epsi, residual_max, residual_norm, &
         z, zeta, rez, rezeta, &
         delz, dz, dzeta, &
         steg, zold, zetaold, new_residual
    real(kind=rp), dimension(this%m) :: y, lambda, s, mu, &
         rey, relambda, remu, res, &
         dely, dellambda, &
         dy, dlambda, ds, dmu, &
         yold, lambdaold, sold, muold
    real(kind=rp), dimension(this%n) :: x, xsi, eta, &
         rex, rexsi, reeta, &
         delx, diagx, dx, dxsi, deta, &
         xold, xsiold, etaold
    real(kind=rp), dimension(4*this%m + 2) :: residual_small
    real(kind=rp), dimension(3*this%n + 4*this%m + 2) :: residual
    real(kind=rp), dimension(2*this%n + 4*this%m + 2) :: xx, dxx

    real(kind=rp), dimension(this%m, this%n) :: GG
    real(kind=rp), dimension(this%m+1) :: bb
    real(kind=rp), dimension(this%m+1, this%m+1) :: AA

    ! using DGESV in lapack to solve
    ! the linear system which needs the following parameters
    integer :: info
    integer, dimension(this%m+1) :: ipiv

    ! Parameters for global communication
    real(kind=rp) :: re_sq_norm
    real(kind=rp) :: minimal_epsilon

    integer :: nglobal

    ! ------------------------------------------------------------------------ !
    ! initial value for the parameters in the subsolve based on
    ! page 15 of "https://people.kth.se/~krille/mmagcmma.pdf"

    epsi = 1.0_rp !100
    x = 0.5_rp * (this%alpha%x + this%beta%x)
    y = 1.0_rp
    z = 1.0_rp
    zeta = 1.0_rp
    lambda = 1.0_rp
    s = 1.0_rp
    xsi = max(1.0_rp, 1.0_rp / (x - this%alpha%x))
    eta = max(1.0_rp, 1.0_rp / (this%beta%x - x))
    mu = max(1.0_rp, 0.5_rp * this%c%x)

    call MPI_Allreduce(this%n, nglobal, 1, &
         MPI_INTEGER, mpi_sum, neko_comm, ierr)

    ! ------------------------------------------------------------------------ !
    ! Computing the minimal epsilon and choose the most conservative one

    minimal_epsilon = max(0.9_rp * this%epsimin, 1.0e-12_rp)
    call MPI_Allreduce(MPI_IN_PLACE, minimal_epsilon, 1, &
         mpi_real_precision, mpi_min, neko_comm, ierr)

    ! ------------------------------------------------------------------------ !
    !  The main loop of the dual-primal interior point method.

    do while (epsi .gt. minimal_epsilon)

       ! --------------------------------------------------------------------- !
       ! Calculating residuals based on
       ! "https://people.kth.se/~krille/mmagcmma.pdf" for the variables
       ! x, y, z, lambda residuals based on eq(5.9a)-(5.9d), respectively.

       associate(p0j => this%p0j%x, q0j => this%q0j%x, &
            pij => this%pij%x, qij => this%qij%x, &
            low => this%low%x, upp => this%upp%x, &
            alpha => this%alpha%x, beta => this%beta%x, &
            c => this%c%x, d => this%d%x, &
            a0 => this%a0, a => this%a%x, &
            bi => this%bi%x)

         rex = (p0j + matmul(transpose(pij), lambda)) / (upp - x)**2 &
              - (q0j + matmul(transpose(qij), lambda)) / (x - low)**2 &
              - xsi + eta

         rey = c + d * y - lambda - mu
         rez = a0 - zeta - dot_product(lambda, a)

         relambda = 0.0_rp
         do i = 1, this%m
            do j = 1, this%n
               ! Accumulate sums for relambda (the term gi(x))
               relambda(i) = relambda(i) &
                    + pij(i, j) / (upp(j) - x(j)) &
                    + qij(i, j) / (x(j) - low(j))
            end do
         end do

       end associate

       ! --------------------------------------------------------------------- !
       ! Computing the norm of the residuals

       ! Complete the computations of lambda residuals
       call MPI_Allreduce(MPI_IN_PLACE, relambda, this%m, &
            mpi_real_precision, mpi_sum, neko_comm, ierr)
       relambda = relambda - this%a%x*z - y + s - this%bi%x

       rexsi = xsi * (x - this%alpha%x) - epsi
       reeta = eta * (this%beta%x - x) - epsi
       remu = mu * y - epsi
       rezeta = zeta * z - epsi
       res = lambda * s - epsi

       ! Setup vectors of residuals and their norms
       residual = [rex, rey, rez, relambda, rexsi, reeta, remu, rezeta, res]
       residual_small = [rey, rez, relambda, remu, rezeta, res]

       residual_max = maxval(abs(residual))
       re_sq_norm = norm2(rex)**2 + norm2(rexsi)**2 + norm2(reeta)**2

       call MPI_Allreduce(MPI_IN_PLACE, residual_max, 1, &
            mpi_real_precision, mpi_max, neko_comm, ierr)

       call MPI_Allreduce(MPI_IN_PLACE, re_sq_norm, &
            1, mpi_real_precision, mpi_sum, neko_comm, ierr)

       residual_norm = sqrt(norm2(residual_small)**2 + re_sq_norm)

       ! --------------------------------------------------------------------- !
       ! Internal loop

       do iter = 1, this%max_iter

          !Check the condition
          if (residual_max .lt. epsi) exit

          delx = 0.0_rp
          do j = 1, this%n
             do i = 1, this%m
                delx(j) = delx(j) &
                     + this%pij%x(i,j) * lambda(i) / (this%upp%x(j) - x(j))**2 &
                     - this%qij%x(i,j) * lambda(i) / (x(j) - this%low%x(j))**2
             end do
          end do

          delx = delx &
               + this%p0j%x / (this%upp%x - x)**2 &
               - this%q0j%x / (x - this%low%x)**2 &
               - epsi / (x - this%alpha%x) &
               + epsi / (this%beta%x - x)

          dely = this%c%x + this%d%x * y - lambda - epsi / y
          delz = this%a0 - dot_product(lambda, this%a%x) - epsi / z

          ! Accumulate sums for dellambda (the term gi(x))
          dellambda = 0.0_rp
          do i = 1, this%m
             do j = 1, this%n
                dellambda(i) = dellambda(i) &
                     + this%pij%x(i, j) / (this%upp%x(j) - x(j)) &
                     + this%qij%x(i, j) / (x(j) - this%low%x(j))
             end do
          end do

          call MPI_Allreduce(MPI_IN_PLACE, dellambda, this%m, &
               mpi_real_precision, mpi_sum, neko_comm, ierr)

          dellambda = dellambda - this%a%x*z - y - this%bi%x + epsi / lambda

          do i = 1, this%m
             GG(i,:) = this%pij%x(i,:) / (this%upp%x - x)**2 &
                  - this%qij%x(i,:) / (x - this%low%x)**2
          end do

          diagx = &
               (this%p0j%x + matmul(transpose(this%pij%x), lambda)) &
               / (this%upp%x - x)**3 &
               + (this%q0j%x + matmul(transpose(this%qij%x), lambda)) &
               / (x - this%low%x)**3

          diagx = 2.0_rp * diagx &
               + xsi / (x - this%alpha%x) &
               + eta / (this%beta%x - x)


          !Here we only consider the case m<n in the matlab code
          !assembling the right hand side matrix based on eq(5.20)
          ! bb = [dellambda + dely/(this%d%x + &
          !         (mu/y)) - matmul(GG,delx/diagx), delz ]
          !!!!!!!!!!!!!!for MPI computation of bb!!!!!!!!!!!!!!!!!!!!!!!!!
          bb = 0.0_rp
          do i = 1, this%m
             do j = 1, this%n
                bb(i) = bb(i) + GG(i, j) * (delx(j) / diagx(j))
             end do
          end do

          call MPI_Allreduce(MPI_IN_PLACE, bb(1:this%m), this%m, &
               mpi_real_precision, mpi_sum, neko_comm, ierr)

          bb(1:this%m) = dellambda + dely / (this%d%x + mu / y) - bb(1:this%m)
          bb(this%m + 1) = delz

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! !assembling the coefficients matrix AA based on eq(5.20)
          ! AA(1:this%m,1:this%m) =  &
          ! matmul(matmul(GG,mma_diag(1/diagx)), transpose(GG))
          ! !update diag(AA)
          ! AA(1:this%m,1:this%m) = AA(1:this%m,1:this%m) + &
          !     mma_diag(s/lambda + 1.0/(this%d%x + (mu/y)))

          AA = 0.0_rp
          !Direct computation of the matrix multiplication
          !(for better performance)
          do i = 1, this%m
             do j = 1, this%m
                ! Compute the (i, j) element of AA
                do k = 1, this%n !this n is global
                   AA(i, j) = AA(i, j) &
                        + GG(i, k) * (1.0_rp / diagx(k)) * GG(j, k)
                end do
             end do
          end do

          call MPI_Allreduce(MPI_IN_PLACE, AA(1:this%m, 1:this%m), &
               this%m*this%m, mpi_real_precision, mpi_sum, neko_comm, ierr)

          do i = 1, this%m
             !update the diag AA
             AA(i, i) = AA(i, i) &
                  + s(i) / lambda(i) &
                  + 1.0_rp / (this%d%x(i) + mu(i) / y(i))
          end do

          AA(1:this%m, this%m+1) = this%a%x
          AA(this%m+1, 1:this%m) = this%a%x
          AA(this%m+1, this%m+1) = - zeta/z

          call DGESV(this%m + 1, 1, AA, this%m + 1, ipiv, bb, this%m + 1, info)

          if (info .ne. 0) then
             write(stderr, *) "DGESV failed to solve the linear system in MMA."
             write(stderr, *) "Please check mma_subsolve_dpip in mma.f90"
             error stop
          end if

          dlambda = bb(1:this%m)
          dz = bb(this%m + 1)

          ! based on eq(5.19)
          dx = - delx / diagx - matmul(transpose(GG), dlambda) / diagx
          dy = (-dely + dlambda) / (this%d%x + mu / y)

          dxsi = -xsi + (epsi - dx * xsi) / (x - this%alpha%x)
          deta = -eta + (epsi + dx * eta) / (this%beta%x - x)
          dmu = -mu + (epsi - mu * dy) / y
          dzeta = -zeta + (epsi - zeta * dz) / z
          ds = -s + (epsi - dlambda * s) / lambda

          dxx = [dy, dz, dlambda, dxsi, deta, dmu, dzeta, ds]
          xx = [y, z, lambda, xsi, eta, mu, zeta, s]

          steg = 1.0_rp / maxval([ &
               1.0_rp, &
               -1.01_rp * dxx / xx, &
               -1.01_rp * dx / (x - this%alpha%x), &
               1.01_rp * dx / (this%beta%x - x) &
               ])

          ! Save the old values
          xold = x
          yold = y
          zold = z
          lambdaold = lambda
          xsiold = xsi
          etaold = eta
          muold = mu
          zetaold = zeta
          sold = s

          ! The innermost loop to determine the suitable step length
          ! using the Backtracking Line Search approach
          new_residual = 2.0_rp * residual_norm

          ! Share the new_residual and steg values
          call MPI_Allreduce(MPI_IN_PLACE, steg, 1, &
               mpi_real_precision, mpi_min, neko_comm, ierr)
          call MPI_Allreduce(MPI_IN_PLACE, new_residual, 1, &
               mpi_real_precision, mpi_min, neko_comm, ierr)

          itto = 0
          do while ((new_residual .gt. residual_norm) .and. (itto .lt. 50))
             itto = itto + 1

             ! update the variables
             x = xold + steg*dx
             y = yold + steg*dy
             z = zold + steg*dz
             lambda = lambdaold + steg*dlambda
             xsi = xsiold + steg*dxsi
             eta = etaold + steg*deta
             mu = muold + steg*dmu
             zeta = zetaold + steg*dzeta
             s = sold + steg*ds

             ! Recompute the new_residual to see if this stepsize improves
             ! the residue
             rex = (this%p0j%x + matmul(transpose(this%pij%x), lambda)) &
                  / (this%upp%x - x)**2 &
                  - (this%q0j%x + matmul(transpose(this%qij%x), lambda)) &
                  / (x - this%low%x)**2 &
                  - xsi + eta

             rey = this%c%x + this%d%x*y - lambda - mu
             rez = this%a0 - zeta - dot_product(lambda, this%a%x)

             ! Accumulate sums for relambda (the term gi(x))
             relambda = 0.0_rp
             do i = 1, this%m
                do j = 1, this%n
                   relambda(i) = relambda(i) &
                        + this%pij%x(i, j) / (this%upp%x(j) - x(j)) &
                        + this%qij%x(i, j) / (x(j) - this%low%x(j))
                end do
             end do

             call MPI_Allreduce(MPI_IN_PLACE, relambda, this%m, &
                  mpi_real_precision, mpi_sum, neko_comm, ierr)

             relambda = relambda - this%a%x*z - y + s - this%bi%x

             rexsi = xsi * (x - this%alpha%x) - epsi
             reeta = eta * (this%beta%x - x) - epsi
             remu = mu * y - epsi
             rezeta = zeta * z - epsi
             res = lambda * s - epsi

             residual_small = [rey, rez, relambda, remu, rezeta, res]

             re_sq_norm = norm2(rex)**2 + norm2(rexsi)**2 + norm2(reeta)**2
             call MPI_Allreduce(MPI_IN_PLACE, re_sq_norm, &
                  1, mpi_real_precision, mpi_sum, neko_comm, ierr)

             new_residual = sqrt(norm2(residual_small)**2 + re_sq_norm)

             steg = steg / 2.0_rp
          end do
          steg = 2.0_rp * steg ! Correction for the final division by 2

          residual = [rex, rey, rez, relambda, rexsi, reeta, remu, rezeta, res]

          ! Update the maximum and norm of the residuals
          residual_norm = new_residual
          residual_max = maxval(abs(residual))
          call MPI_Allreduce(MPI_IN_PLACE, residual_max, 1, &
               mpi_real_precision, mpi_max, neko_comm, ierr)
       end do

       epsi = 0.1_rp * epsi
    end do

    ! Save the new designx
    this%xold2 = this%xold1
    this%xold1%x = designx
    designx = x

    !update the parameters of the MMA object nesessary to compute KKT residual
    this%y%x = y
    this%z = z
    this%lambda%x = lambda
    this%zeta = zeta
    this%xsi%x = xsi
    this%eta%x = eta
    this%mu%x = mu
    this%s%x = s

  end subroutine mma_subsolve_dpip_cpu

  !> Implementation of the KKT residual computation for the MMA algorithm.
  subroutine mma_KKT_cpu(this, x, df0dx, fval, dfdx)
    ! ----------------------------------------------------- !
    ! Compute the KKT condition right hand side for a given !
    ! designx x and set the max and norm values of the       !
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
    real(kind=rp), dimension(3*this%n+4*this%m+2) :: residual

    real(kind=rp), dimension(4*this%m+2) :: residual_small
    integer :: ierr
    real(kind=rp) :: re_sq_norm

    rex = df0dx + matmul(transpose(dfdx), this%lambda%x) &
         - this%xsi%x + this%eta%x
    rey = this%c%x + this%d%x*this%y%x - this%lambda%x - this%mu%x
    rez = this%a0 - this%zeta - dot_product(this%lambda%x, this%a%x)

    relambda = fval - this%a%x * this%z - this%y%x + this%s%x
    rexsi = this%xsi%x * (x - this%xmin%x)
    reeta = this%eta%x * (this%xmax%x - x)
    remu = this%mu%x * this%y%x
    rezeta = this%zeta * this%z
    res = this%lambda%x * this%s%x

    residual = [rex, rey, rez, relambda, rexsi, reeta, remu, rezeta, res]
    residual_small = [rey, rez, relambda, remu, rezeta, res]

    this%residumax = maxval(abs(residual))
    re_sq_norm = norm2(rex)**2 + norm2(rexsi)**2 + norm2(reeta)**2

    call MPI_Allreduce(MPI_IN_PLACE, this%residumax, 1, &
         mpi_real_precision, mpi_max, neko_comm, ierr)

    call MPI_Allreduce(MPI_IN_PLACE, re_sq_norm, 1, &
         mpi_real_precision, mpi_sum, neko_comm, ierr)

    this%residunorm = sqrt(norm2(residual_small)**2 + re_sq_norm)

  end subroutine mma_KKT_cpu
end submodule mma_cpu
