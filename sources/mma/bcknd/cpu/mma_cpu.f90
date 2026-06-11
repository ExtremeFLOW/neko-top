!> @file mma_cpu.f90
!! @copyright (c) 2024-2025, The Neko-TOP Authors
!! All rights reserved.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions
!! are met:
!!
!!   * Redistributions of source code must retain the above copyright
!!     notice, this list of conditions and the following disclaimer.
!!
!!   * Redistributions in binary form must reproduce the above
!!     copyright notice, this list of conditions and the following
!!     disclaimer in the documentation and/or other materials provided
!!     with the distribution.
!!
!!   * Neither the name of the authors nor the names of its
!!     contributors may be used to endorse or promote products derived
!!     from this software without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
!! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
!! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
!! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
!! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
!! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
!! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!! POSSIBILITY OF SUCH DAMAGE.

submodule (mma) mma_cpu
  use lapack_interfaces, only: dgesv
  use mpi_f08, only: MPI_IN_PLACE, MPI_MAX, MPI_MIN
  use comm, only: neko_comm, pe_rank, mpi_real_precision
  use math, only: NEKO_EPS
  use profiler, only: profiler_start_region, profiler_end_region

  implicit none

contains

  module subroutine mma_update_cpu(this, iter, x, df0dx, fval, dfdx)
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
    real(kind=rp), dimension(this%n), intent(in) :: df0dx
    real(kind=rp), dimension(this%m), intent(in) :: fval
    real(kind=rp), dimension(this%m, this%n), intent(in) :: dfdx

    if (.not. this%is_initialized) then
       call neko_error("The MMA object is not initialized.")
    end if

    call profiler_start_region("MMA update")

    ! generate a convex approximation of the problem
    call profiler_start_region("MMA gensub")
    call mma_gensub_cpu(this, iter, x, df0dx, fval, dfdx)
    call profiler_end_region("MMA gensub")

    !solve the approximation problem using interior point method
    call profiler_start_region("MMA subsolve")
    if (this%subsolver .eq. "dip") then
       call mma_subsolve_dip_cpu(this, x)
    else
       call mma_subsolve_dpip_cpu(this, x)
    end if
    call profiler_end_region("MMA subsolve")

    call profiler_end_region("MMA update")

    this%is_updated = .true.
  end subroutine mma_update_cpu

  !> Implementation of the KKT residual computation for the MMA algorithm.
  module subroutine mma_KKT_cpu(this, x, df0dx, fval, dfdx)
    ! ----------------------------------------------------- !
    ! Compute the KKT condition right hand side for a given !
    ! designx x and set the max and norm values of the       !
    ! residue of KKT system to this%residumax and           !
    ! this%residunorm.                                      !
    !                                                       !
    ! The left hand sides of the KKT conditions are computed!
    ! for the following nonlinear programming problem:      !
    ! Minimize  f_0(x) + a_0*z +                            !
    !                       sum(c_i*y_i + 0.5*d_i*(y_i)^2)!
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
    real(kind=rp), dimension(this%n), intent(in) :: df0dx
    real(kind=rp), dimension(this%m), intent(in) :: fval
    real(kind=rp), dimension(this%m, this%n), intent(in) :: dfdx

    call profiler_start_region("MMA KKT computation")

    if (this%subsolver .eq. "dip") then
       call mma_dip_KKT_cpu(this, x, df0dx, fval, dfdx)
    else
       call mma_dpip_KKT_cpu(this, x, df0dx, fval, dfdx)
    end if

    call profiler_end_region("MMA KKT computation")
  end subroutine mma_KKT_cpu

  !> Implementation of the KKT residual computation for dual primal interior
  ! point method (dpip) subsolve of MMA algorithm.
  module subroutine mma_dpip_KKT_cpu(this, x, df0dx, fval, dfdx)
    ! ----------------------------------------------------- !
    ! Compute the KKT condition right hand side for a given !
    ! designx x and set the max and norm values of the      !
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
    real(kind=rp), dimension(this%n), intent(in) :: df0dx
    real(kind=rp), dimension(this%m), intent(in) :: fval
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
         mpi_real_precision, MPI_MAX, neko_comm, ierr)

    call MPI_Allreduce(MPI_IN_PLACE, re_sq_norm, 1, &
         mpi_real_precision, mpi_sum, neko_comm, ierr)

    this%residunorm = sqrt(norm2(residual_small)**2 + re_sq_norm)
  end subroutine mma_dpip_KKT_cpu

  !> Implementation of the KKT residual computation for dual interior
  ! point method (dip) subsolve of MMA algorithm.
  module subroutine mma_dip_KKT_cpu(this, x, df0dx, fval, dfdx)
    ! ----------------------------------------------------- !
    ! Compute the KKT condition right hand side for a given !
    ! designx x and set the max and norm values of the      !
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
    real(kind=rp), dimension(this%n), intent(in) :: df0dx
    real(kind=rp), dimension(this%m), intent(in) :: fval
    real(kind=rp), dimension(this%m, this%n), intent(in) :: dfdx

    real(kind=rp), dimension(this%m) :: relambda, remu
    real(kind=rp), dimension(2*this%m) :: residual


    relambda = fval - this%a%x * this%z - this%y%x + this%mu%x
    ! Compute residual for mu (eta in the paper)
    remu = this%lambda%x * this%mu%x

    residual = abs([relambda, remu])
    this%residumax = maxval(residual)
    this%residunorm = norm2(residual)

  end subroutine mma_dip_KKT_cpu

  !============================================================================!
  ! private internal subroutines

  !> generat a subproblem; convex approximation of the optimization problem
  subroutine mma_gensub_cpu(this, iter, x, df0dx, fval, dfdx)
    ! ----------------------------------------------------- !
    ! Generate the approximation sub problem by computing   !
    ! the lower and upper asymtotes and the other necessary !
    ! parameters (alpha, beta, p0j, q0j, pij, qij, ...).    !
    ! ----------------------------------------------------- !
    class(mma_t), intent(inout) :: this
    real(kind=rp), dimension(this%n), intent(in) :: x
    real(kind=rp), dimension(this%n), intent(in) :: df0dx
    real(kind=rp), dimension(this%m), intent(in) :: fval
    real(kind=rp), dimension(this%m, this%n), intent(in) :: dfdx
    integer, intent(in) :: iter
    integer :: i, j, ierr
    real(kind=rp), dimension(this%n) :: xmin_eff, xmax_eff
    real(kind=rp), dimension(this%n) :: x_diff
    real(kind=rp) :: asy_factor

    xmin_eff = this%xmin%x
    xmax_eff = this%xmax%x
    if (this%move_limit .gt. 0.0_rp) then
       xmin_eff = max(xmin_eff, x - this%move_limit)
       xmax_eff = min(xmax_eff, x + this%move_limit)
    end if

    x_diff = max(xmax_eff - xmin_eff, 1.0e-5_rp)

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
         xmin => xmin_eff, xmax => xmax_eff, &
         low => this%low%x, upp => this%upp%x, x => x)

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

  !> solve the subproblem defined by this%pij, this%qij, etc using Dual-primal
  !! interior point method.
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
    ! Note that there is a local dummy "x" in this subroutine, thus, we call
    ! the current design "designx" instead of just "x"
    integer :: i, j, k, iter, itto, ierr
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
    real(kind=rp), dimension(this%m * this%m) :: AA_buffer

    ! using DGESV in lapack to solve
    ! the linear system which needs the following parameters
    integer :: info
    integer, dimension(this%m+1) :: ipiv

    ! Parameters for global communication
    real(kind=rp) :: re_sq_norm
    real(kind=rp) :: minimal_epsilon

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

    ! ------------------------------------------------------------------------ !
    ! Computing the minimal epsilon and choose the most conservative one

    minimal_epsilon = max(0.9_rp * this%epsimin, 1.0e-12_rp)
    call MPI_Allreduce(MPI_IN_PLACE, minimal_epsilon, 1, &
         mpi_real_precision, MPI_MIN, neko_comm, ierr)

    ! ------------------------------------------------------------------------ !
    ! The main loop of the dual-primal interior point method.

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
            a0 => this%a0, a => this%a%x)

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
            mpi_real_precision, MPI_MAX, neko_comm, ierr)

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

          !--------------------------------------------------------------------!
          ! for MPI computation of bb

          bb = 0.0_rp
          do i = 1, this%m
             do j = 1, this%n
                bb(i) = bb(i) + GG(i, j) * (delx(j) / diagx(j))
             end do
          end do

          call MPI_Allreduce(MPI_IN_PLACE, bb, this%m, &
               mpi_real_precision, mpi_sum, neko_comm, ierr)

          bb(1:this%m) = dellambda + dely / (this%d%x + mu / y) - bb(1:this%m)
          bb(this%m + 1) = delz

          !--------------------------------------------------------------------!
          ! assembling the coefficients matrix AA based on eq(5.20)
          ! AA(1:this%m,1:this%m) =  &
          ! matmul(matmul(GG,mma_diag(1/diagx)), transpose(GG))
          ! !update diag(AA)
          ! AA(1:this%m,1:this%m) = AA(1:this%m,1:this%m) + &
          !     mma_diag(s/lambda + 1.0/(this%d%x + (mu/y)))

          AA = 0.0_rp
          ! Direct computation of the matrix multiplication
          ! (for better performance)
          do i = 1, this%m
             do j = 1, this%m
                ! Compute the (i, j) element of AA
                do k = 1, this%n !this n is global
                   AA(i, j) = AA(i, j) &
                        + GG(i, k) * (1.0_rp / diagx(k)) * GG(j, k)
                end do
             end do
          end do

          AA_buffer = reshape(AA(1:this%m, 1:this%m), [this%m * this%m])

          call MPI_Allreduce(MPI_IN_PLACE, AA_buffer, &
               this%m*this%m, mpi_real_precision, mpi_sum, neko_comm, ierr)

          AA(1:this%m, 1:this%m) = reshape(AA_buffer, [this%m, this%m])

          do i = 1, this%m
             ! update the diag AA
             AA(i, i) = AA(i, i) &
                  + s(i) / lambda(i) &
                  + 1.0_rp / (this%d%x(i) + mu(i) / y(i))
          end do

          AA(1:this%m, this%m+1) = this%a%x
          AA(this%m+1, 1:this%m) = this%a%x
          AA(this%m+1, this%m+1) = - zeta/z

          call DGESV(this%m + 1, 1, AA, this%m + 1, ipiv, bb, this%m + 1, info)

          if (info .ne. 0) then
             call neko_error("DGESV failed to solve the linear system in " // &
                  "mma_subsolve_dpip.")
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

          new_residual = 2.0_rp * residual_norm

          ! Share the new_residual and steg values
          call MPI_Allreduce(MPI_IN_PLACE, steg, 1, &
               mpi_real_precision, MPI_MIN, neko_comm, ierr)
          call MPI_Allreduce(MPI_IN_PLACE, new_residual, 1, &
               mpi_real_precision, MPI_MIN, neko_comm, ierr)

          ! The innermost loop to determine the suitable step length
          ! using the Backtracking Line Search approach
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

             ! Compute squared norms for the residuals
             re_sq_norm = norm2(rex)**2 + norm2(rexsi)**2 + norm2(reeta)**2
             call MPI_Allreduce(MPI_IN_PLACE, re_sq_norm, &
                  1, mpi_real_precision, mpi_sum, neko_comm, ierr)

             residual_small = [rey, rez, relambda, remu, rezeta, res]
             new_residual = sqrt(norm2(residual_small)**2 + re_sq_norm)
             call MPI_Allreduce(MPI_IN_PLACE, new_residual, 1, &
                  mpi_real_precision, MPI_MIN, neko_comm, ierr)

             steg = steg / 2.0_rp
          end do
          steg = 2.0_rp * steg ! Correction for the final division by 2

          residual = [rex, rey, rez, relambda, rexsi, reeta, remu, rezeta, res]

          ! Update the maximum and norm of the residuals
          residual_norm = new_residual
          residual_max = maxval(abs(residual))
          call MPI_Allreduce(MPI_IN_PLACE, residual_max, 1, &
               mpi_real_precision, MPI_MAX, neko_comm, ierr)
       end do

       epsi = 0.1_rp * epsi
    end do

    ! Save the new designx
    this%xold2%x = this%xold1%x
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

  !> solve the subproblem defined by this%pij, this%qij, etc using a pure Dual
  !! interior point method.
  subroutine mma_subsolve_dip_cpu(this, designx)
    ! ------------------------------------------------------------------------ !
    ! -------------------------------Dual Solver------------------------------ !
    ! ------------------------------------------------------------------------ !
    ! This implementation is based on:                                         !
    ! https://doi.org/10.1007/s00158-012-0869-2                                !
    ! Definition of the Lagrangian function:                                   !
    ! (Note that the equation is slightly different with d(i)=1 and a quadratic!
    ! term for z. This is done to ensure that we have quadratic terms for both !
    ! y and z.)                                                                !
    !                                                                          !
    !     L(x, y, z, λ) =                                                      !
    !       sum_{j=1}^{n} [ (p_{0j} + sum_{i=1}^{m} λ_i * p_{ij}) / (u_j - x_j)!
    !                   + (q_{0j} + sum_{i=1}^{m} λ_i * q_{ij}) / (x_j - l_j) ]!
    !       - sum_{i=1}^{m} λ_i * b_i                                          !
    !       + sum_{i=1}^{m} [ (c_i - λ_i) * y_i + 0.5 * y_i^2 ]                !
    !       + (a_0 - sum_{i=1}^{m} λ_i * a_i) * z + 0.5 * z^2                  !
    !                                                                          !
    ! Breakdown of terms:                                                      !
    !   - Terms related to x:  L_x (the first three lines of L(x, y, z, λ))    !
    !   - Terms related to y:  L_y (the fourth line of L(x, y, z, λ))          !
    !   - Terms related to z:  L_z (the last line of L(x, y, z, λ))            !
    !                                                                          !
    ! Optimization problem if λ is given:                                      !
    !                                                                          !
    !     Minimize L(x, y, z, λ)                                               !
    !     subject to: α_j ≤ x_j ≤ β_j, z ≥ 0, and y_i ≥ 0 for all i, j.        !
    !                                                                          !
    ! Since the problem is separable:                                          !
    !     Ψ(λ) =                                                               !
    !       sum_{j=1}^{n} min_xj {L_x(x_j, λ) | α_j ≤ x_j ≤ β_j}               !
    !       + min_z {L_z(z, λ) | z ≥ 0}                                        !
    !       + sum_{i=1}^{m} min_yi {L_y(y_i, λ) | y_i ≥ 0}                     !
    !                                                                          !
    ! Maximize Ψ(λ) subject to λ_i ≥ 0 for i = 1, ..., m.                      !
    !                                                                          !
    ! ------------------------------------------------------------------------ !

    class(mma_t), intent(inout) :: this
    real(kind=rp), dimension(this%n), intent(inout) :: designx
    ! Note that there is a local dummy "x" in this subroutine, thus, we call
    ! the current design "designx" instead of just "x"
    integer :: i, j, k, iter, ierr
    real(kind=rp) :: epsi, residual_max, z, steg
    real(kind=rp), dimension(this%m) :: y, lambda, mu, &
         relambda, remu, dlambda, dmu, gradlambda
    real(kind=rp), dimension(this%n) :: x, pjlambda, qjlambda

    ! To compute the Hessian based on eq(13)
    ! https://doi.org/10.1007/s00158-012-0869-2

    ! inverse of a diag matrix:
    real(kind=rp), dimension(this%n) :: Ljjxinv ! [∇_x^2 Ljj]−1
    real(kind=rp), dimension(this%m,this%n) :: hijx ! ∇_x hij
    real(kind=rp), dimension(this%m,this%m) :: Hess
    real(kind=rp) :: Hesstrace


    ! using DGESV in lapack to solve
    ! the linear system which needs the following parameters
    integer :: info
    integer, dimension(this%m+1) :: ipiv

    ! Parameters for global communication
    real(kind=rp) :: minimal_epsilon

    ! ------------------------------------------------------------------------ !
    ! initial value for the parameters in the subsolve based on
    ! page 15 of "https://people.kth.se/~krille/mmagcmma.pdf"

    epsi = 1.0_rp !100
    ! x = 0.5_rp * (this%alpha%x + this%beta%x)
    y = 1.0_rp
    z = 0.0_rp
    lambda = max(1.0_rp, 0.5_rp * this%c%x)
    mu = 1.0_rp !this parameter is eta in Niel's paper
    ! note that mu in the paper translates to epsi in the code following the
    ! same style as the Cpp code by Neils

    ! ------------------------------------------------------------------------ !
    ! Computing the minimal epsilon and choose the most conservative one

    minimal_epsilon = max(0.9_rp * this%epsimin, 1.0e-12_rp)
    call MPI_Allreduce(MPI_IN_PLACE, minimal_epsilon, 1, &
         mpi_real_precision, MPI_MIN, neko_comm, ierr)

    ! ------------------------------------------------------------------------ !
    ! The main loop of the dual-primal interior point method.

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
         ! minimize(L_x, L_y, L_z) and compute x(λ), y(λ), z(λ) for
         ! the initial value of λ

         ! Comput the value of y that minimizes L_y for the current λ
         ! minimize (sum_{i=1}^{m} [ (c_i - λ_i) * y_i + 0.5 * y_i^2 ])
         ! dL_y/dy =0   => y= (λ_i - c_i), ensure y>=0
         y = max(0.0_rp, lambda - c)

         ! Comput the value of z that minimizes L_z for the current λ
         ! minimize ((a_0 - sum_{i=1}^{m} λ_i * a_i) * z + 0.5 * z^2)
         ! ensure z>=0
         z = max(0.0_rp,dot_product(lambda, a) - a0)

         ! Comput the value of x that minimizes L_x for the current λ
         ! minimize( sum_{j=1}^{n} [ (p_{0j} + sum_{i=1}^{m} ��_i *
         ! p_{ij}) / (u_j - x_j) + (q_{0j} + sum_{i=1}^{m} λ_i * q_{ij}) /
         ! (x_j - l_j) ] - sum_{i=1}^{m} λ_i * b_i)
         pjlambda = (p0j + matmul(transpose(pij), lambda))
         qjlambda = (q0j + matmul(transpose(qij), lambda))
         x = (sqrt(pjlambda) * low + sqrt(qjlambda) * upp) / &
              (sqrt(pjlambda) + sqrt(qjlambda))

         ! Ensure that x is feasible (alpha<=x<=beta)
         x = merge(alpha, x, x .lt. alpha)
         x = merge(beta, x, x .gt. beta)

         ! Compute the residual for the lambda and mu using eq(9) and eq(15)
         relambda = matmul(pij, 1.0_rp / (upp - x)) + &
              matmul(qij, 1.0_rp / (x - low))

         ! Global comminucation for relambda values
         call MPI_Allreduce(MPI_IN_PLACE, relambda, this%m, &
              mpi_real_precision, mpi_sum, neko_comm, ierr)
         relambda = relambda - bi - y - a * z + mu

         ! Compute residual for mu (eta in the paper)
         remu = mu * lambda - epsi

         residual_max = maxval(abs([relambda, remu]))
         call MPI_Allreduce(MPI_IN_PLACE, residual_max, 1, &
              mpi_real_precision, MPI_MAX, neko_comm, ierr)

         ! ------------------------------------------------------------------- !
         ! Internal loop
         do iter = 1, this%max_iter

            !Check the condition
            if (residual_max .lt. epsi) exit

            ! Compute dL(x, y, z, λ)/dλ for the updated x(λ), y(λ), z(λ)

            gradlambda = matmul(pij, 1.0_rp / (upp - x)) + &
                 matmul(qij, 1.0_rp/(x - low))

            ! Global comminucation for gradlambda values
            call MPI_Allreduce(MPI_IN_PLACE, gradlambda, this%m, &
                 mpi_real_precision, mpi_sum, neko_comm, ierr)
            gradlambda = gradlambda - bi - y - a * z

            ! Update gradlambda as the right hand side for Newton's method(eq10)
            gradlambda = - gradlambda - epsi / lambda

            ! Computing the Hessian as in equation (13) in
            !! https://doi.org/10.1007/s00158-012-0869-2

            !--------------contributions of x terms to Hess--------------------!
            Ljjxinv= - 1.0_rp / ( (2.0_rp * pjlambda/(upp - x)**3) + &
                 (2.0_rp * qjlambda/(x - low)**3))


            ! Remove the sensitivity for the active primal constraints
            Ljjxinv = merge(0.0_rp, Ljjxinv, x - alpha < NEKO_EPS)
            Ljjxinv = merge(0.0_rp, Ljjxinv, beta - x < NEKO_EPS)

            do i = 1, this%m
               hijx(i,:) = pij(i,:) / (upp - x)**2 &
                    - qij(i,:) / (x - low)**2
            end do

            Hess = 0.0_rp
            ! Direct computation of the matrix multiplication
            ! (for better performance)
            do i = 1, this%m
               do j = 1, this%m
                  ! Compute the (i, j) element of AA
                  do k = 1, this%n !this n is global
                     Hess(i, j) = Hess(i, j) &
                          + hijx(i, k) * (Ljjxinv(k)) * hijx(j, k)
                  end do
               end do
            end do

            call MPI_Allreduce(MPI_IN_PLACE, Hess, &
                 this%m * this%m, mpi_real_precision, mpi_sum, neko_comm, ierr)

            !---------------contributions of z terms to Hess-------------------!
            ! Only for inactive constraint, we consider contributions to Hess
            ! based on the cpp code by Niels.
            if (dot_product(lambda, a) .gt. 0.0_rp) then
               do i = 1, this%m
                  do j = 1, this%m
                     Hess(i, j) = Hess(i, j) - a(i) * a(j)
                  end do
               end do
            end if

            !---------------contributions of y terms to Hess-------------------!
            ! Only for inactive constraint, we consider contributions to Hess.
            do i = 1, this%m
               if (y(i) .gt. 0.0_rp) then
                  Hess(i, i) = Hess(i, i) - 1.0_rp
               end if
               ! Based on eq(10), note the term (-\Omega \Lambda)
               Hess(i, i) = Hess(i, i) - mu(i) / lambda(i)
            end do

            ! Improve the robustness by stablizing the Hess using
            ! Levenberg-Marquardt algorithm (heuristically)
            Hesstrace = 0.0_rp
            do i=1, this%m
               Hesstrace = Hesstrace + Hess(i, i)
            end do
            do i=1, this%m
               Hess(i,i) = Hess(i, i) - &
                    max(-1.0e-4_rp*Hesstrace/this%m, 1.0e-7_rp)
            end do

            call DGESV(this%m , 1, Hess, this%m , ipiv, &
                 gradlambda, this%m, info)

            if (info .ne. 0) then
               call neko_error("DGESV failed to solve the linear system in " // &
                    "mma_subsolve_dip.")
            end if
            dlambda = gradlambda

            ! based on eq(11) for delta eta
            dmu = -mu + epsi / lambda - dlambda * mu / lambda

            ! Compute the stepsize and update lambda and mu (eta in the paper)

            steg = 1.005_rp
            do i = 1, this%m
               steg = merge(-1.01_rp * dlambda(i) / lambda(i), steg, &
                    steg < -1.01_rp * dlambda(i) / lambda(i))
               steg = merge(-1.01_rp * dmu(i) / mu(i), &
                    steg, steg < -1.01_rp * dmu(i) / mu(i))
            end do

            steg = 1.0_rp / steg

            lambda = lambda + steg*dlambda
            mu = mu + steg*dmu

            ! minimize(L_x, L_y, L_z) and compute x(λ), y(λ), z(λ) for
            ! the updated values of λ

            ! Comput the value of y that minimizes L_y for the current λ
            ! minimize (sum_{i=1}^{m} [ (c_i - λ_i) * y_i + 0.5 * y_i^2 ])
            ! dL_y/dy =0   => y= (λ_i - c_i), ensure y>=0
            y = max(0.0_rp, lambda - c)


            ! Comput the value of z that minimizes L_z for the current λ
            ! minimize ((a_0 - sum_{i=1}^{m} λ_i * a_i) * z + 0.5 * z^2)
            ! ensure z>=0
            z = max(0.0_rp,dot_product(lambda, a) - a0)

            ! Comput the value of x that minimizes L_x for the current λ
            ! minimize( sum_{j=1}^{n} [ (p_{0j} + sum_{i=1}^{m} λ_i *
            ! p_{ij}) / (u_j - x_j) + (q_{0j} + sum_{i=1}^{m} λ_i * q_{ij}) /
            ! (x_j - l_j) ] - sum_{i=1}^{m} λ_i * b_i)
            pjlambda = (p0j + matmul(transpose(pij), lambda))
            qjlambda = (q0j + matmul(transpose(qij), lambda))
            x = (sqrt(pjlambda) * low + sqrt(qjlambda) * upp) / &
                 (sqrt(pjlambda) + sqrt(qjlambda))

            ! Ensure that x is feasible (alpha<=x<=beta)
            x = merge(alpha, x, x .lt. alpha)
            x = merge(beta, x, x .gt. beta)

            ! Compute the residual for the lambda and mu using eq(9) and eq(15)
            relambda = matmul(pij, 1.0_rp / (upp - x)) + &
                 matmul(qij, 1.0_rp / (x - low))
            ! Global comminucation for relambda values
            call MPI_Allreduce(MPI_IN_PLACE, relambda, this%m, &
                 mpi_real_precision, mpi_sum, neko_comm, ierr)
            relambda = relambda - bi - y - a * z + mu

            ! Compute residual for mu (eta in the paper)
            remu = mu * lambda - epsi

            residual_max = maxval(abs([relambda, remu]))
            call MPI_Allreduce(MPI_IN_PLACE, residual_max, 1, &
                 mpi_real_precision, MPI_MAX, neko_comm, ierr)
         end do
       end associate

       epsi = 0.1_rp * epsi
    end do

    ! Save the new designx
    this%xold2%x = this%xold1%x
    this%xold1%x = designx
    designx = x

    !update the parameters of the MMA object nesessary to compute KKT residual

    this%y%x = y
    this%z = z
    this%lambda%x = lambda
    this%mu%x = mu
  end subroutine mma_subsolve_dip_cpu

end submodule mma_cpu
