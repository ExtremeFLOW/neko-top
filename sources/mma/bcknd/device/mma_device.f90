!> @file mma_device.f90
!! @copyright
!! Copyright (c) 2025, The Neko-TOP Authors
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

submodule (mma) mma_device

  use device_math, only: device_copy, device_cmult, device_cadd, device_cfill, &
       device_add2, device_add3s2, device_invcol2, device_col2, device_col3, &
       device_sub2, device_sub3, device_add2s2, device_cadd2, device_pwmax2, &
       device_pwmin2, device_cpwmax2, device_glsum, device_cmult2
  use device_mma_math, only: device_maxval, device_norm, device_lcsc2, &
       device_maxval2, device_maxval3, device_mma_gensub3, &
       device_mma_gensub4, device_mma_max, device_max2, device_rex, &
       device_relambda, device_delx, device_add2inv2, device_gg, device_diagx, &
       device_bb, device_updatebb, device_aa, device_updateaa, device_dx, &
       device_dy, device_dxsi, device_deta, device_kkt_rex, &
       device_mma_gensub2, device_mattrans_v_mul, device_mma_dipsolvesub1, &
       device_mma_Ljjxinv, device_Hess, device_solve_linear_system, &
       device_prepare_hessian, device_prepare_aa_matrix, device_update_hessian_z, &
       device_unconstrained_kkt

  use neko_config, only: NEKO_BCKND_DEVICE, NEKO_DEVICE_MPI
  use device, only: DEVICE_TO_HOST
  use comm, only: neko_comm, pe_rank, mpi_real_precision
  use mpi_f08, only: MPI_IN_PLACE, MPI_MAX, MPI_MIN
  use profiler, only: profiler_start_region, profiler_end_region
  use math, only: NEKO_EPS

  implicit none

contains

  module subroutine mma_update_device(this, iter, x, df0dx, fval, dfdx)
    ! ----------------------------------------------------- !
    ! Update the design variable x by solving the convex    !
    ! approximation of the problem.                         !
    !                                                       !
    ! This subroutine is called in each iteration of the    !
    ! optimization loop                                     !
    ! ----------------------------------------------------- !
    class(mma_t), intent(inout) :: this
    integer, intent(in) :: iter
    type(c_ptr), intent(inout) :: x
    type(c_ptr), intent(in) :: df0dx, fval, dfdx

    if (.not. this%is_initialized) then
       call neko_error("The MMA object is not initialized.")
    end if

    call profiler_start_region("MMA gensub")
    ! generate a convex approximation of the problem
    call mma_gensub_device(this, iter, x, df0dx, fval, dfdx)
    call profiler_end_region("MMA gensub")

    !solve the approximation problem using interior point method
    call profiler_start_region("MMA subsolve")

    if (this%unconstrained_problem) then
       call mma_subsolve_unconstrained_device(this, x)
    elseif (this%subsolver .eq. "dip") then
       call mma_subsolve_dip_device(this, x)
    else if (this%subsolver .eq. "dpip") then
       call mma_subsolve_dpip_device(this, x)
    else
       call neko_error("Unrecognized subsolver for MMA in mma_device.")
    end if
    call profiler_end_region("MMA subsolve")

    this%is_updated = .true.
  end subroutine mma_update_device

  module subroutine mma_KKT_device(this, x, df0dx, fval, dfdx)
    class(mma_t), intent(inout) :: this
    type(c_ptr), intent(in) :: x, df0dx, fval, dfdx

    if (this%unconstrained_problem) then
       call mma_unconstrained_KKT_device(this, x, df0dx)
    elseif (this%subsolver .eq. "dip") then
       call mma_dip_KKT_device(this, x, df0dx, fval, dfdx)
    elseif (this%subsolver .eq. "dpip") then
       call mma_dpip_KKT_device(this, x, df0dx, fval, dfdx)
    else
       call neko_error("MMA subsolver is not recognised (NOT dip or dpip).")
    end if
  end subroutine mma_KKT_device

  module subroutine mma_unconstrained_KKT_device(this, x, df0dx)
    class(mma_t), intent(inout) :: this
    type(c_ptr), intent(in) :: x, df0dx

    type(vector_t), pointer :: rex
    integer :: ind
    real(kind=rp) :: re_sq_norm
    integer :: ierr

    ! Request temporary scratch space for rex on the device
    call this%scratch%request(rex, ind, this%n, .false.)

    ! Launch the custom element-wise KKT kernel via our abstraction layer
    call device_unconstrained_kkt(rex%x_d, x, this%xmin%x_d, this%xmax%x_d, &
         df0dx, NEKO_EPS, this%n)

    ! Compute global metrics using your existing library routines
    ! (Assuming device_maxval returns maxval(abs(a)) and device_norm returns sum(a^2))
    this%residumax = device_maxval(rex%x_d, this%n)
    re_sq_norm = device_norm(rex%x_d, this%n)

    ! Global MPI Reductions across nodes
    call MPI_Allreduce(MPI_IN_PLACE, this%residumax, 1, &
         mpi_real_precision, MPI_MAX, neko_comm, ierr)

    call MPI_Allreduce(MPI_IN_PLACE, re_sq_norm, 1, &
         mpi_real_precision, mpi_sum, neko_comm, ierr)

    this%residunorm = sqrt(re_sq_norm)

    ! Release scratch memory
    call this%scratch%relinquish(ind)
  end subroutine mma_unconstrained_KKT_device

  !> Implementation of the KKT residual computation for dual interior
  ! point method (dip) subsolve of MMA algorithm.
  module subroutine mma_dip_KKT_device(this, x, df0dx, fval, dfdx)
    class(mma_t), intent(inout) :: this
    type(c_ptr), intent(in) :: x, df0dx, fval, dfdx

    type(vector_t), pointer :: relambda, remu
    integer :: ind(2)

    call this%scratch%request(relambda, ind(1), this%m, .false.)
    call this%scratch%request(remu, ind(2), this%m, .false.)

    ! relambda = fval - this%a%x * this%z - this%y%x + this%mu%x
    call device_add3s2(relambda%x_d, fval, this%a%x_d, 1.0_rp, -this%z, &
         this%m)
    call device_sub2(relambda%x_d, this%y%x_d, this%m)
    call device_add2(relambda%x_d, this%mu%x_d, this%m)

    ! Compute residual for mu (eta in the paper)
    call device_col3(remu%x_d, this%lambda%x_d, this%mu%x_d, this%m)

    this%residumax = maxval([device_maxval(relambda%x_d, this%m), &
         device_maxval(remu%x_d, this%m)])
    this%residunorm = sqrt(device_norm(relambda%x_d, this%m)+ &
         device_norm(remu%x_d, this%m))

    call this%scratch%relinquish(ind)
  end subroutine mma_dip_KKT_device

  !> Implementation of the KKT residual computation for dual primal interior
  ! point method (dpip) subsolve of MMA algorithm.
  module subroutine mma_dpip_KKT_device(this, x, df0dx, fval, dfdx)
    class(mma_t), intent(inout) :: this
    type(c_ptr), intent(in) :: x, df0dx, fval, dfdx

    real(kind=rp) :: rez, rezeta
    type(vector_t), pointer :: rey, relambda, remu, res
    type(vector_t), pointer :: rex, rexsi, reeta
    integer :: ierr, ind(7)
    real(kind=rp) :: re_sq_norm

    call this%scratch%request(rey, ind(1), this%m, .false.)
    call this%scratch%request(relambda, ind(2), this%m, .false.)
    call this%scratch%request(remu, ind(3), this%m, .false.)
    call this%scratch%request(res, ind(4), this%m, .false.)

    call this%scratch%request(rex, ind(5), this%n, .false.)
    call this%scratch%request(rexsi, ind(6), this%n, .false.)
    call this%scratch%request(reeta, ind(7), this%n, .false.)

    call device_kkt_rex(rex%x_d, df0dx, dfdx, this%xsi%x_d, &
         this%eta%x_d, this%lambda%x_d, this%n, this%m)

    call device_col3(rey%x_d, this%d%x_d, this%y%x_d, this%m)
    call device_add2(rey%x_d, this%c%x_d, this%m)
    call device_sub2(rey%x_d, this%lambda%x_d, this%m)
    call device_sub2(rey%x_d, this%mu%x_d, this%m)

    rez = this%a0 - this%zeta - device_lcsc2(this%lambda%x_d, this%a%x_d, &
         this%m)

    call device_add3s2(relambda%x_d, fval, this%a%x_d, 1.0_rp, -this%z, &
         this%m)
    call device_sub2(relambda%x_d, this%y%x_d, this%m)
    call device_add2(relambda%x_d, this%s%x_d, this%m)

    call device_sub3(rexsi%x_d, x, this%xmin%x_d, this%n)
    call device_col2(rexsi%x_d, this%xsi%x_d, this%n)

    call device_sub3(reeta%x_d, this%xmax%x_d, x, this%n)
    call device_col2(reeta%x_d, this%eta%x_d, this%n)

    call device_col3(remu%x_d, this%mu%x_d, this%y%x_d, this%m)

    rezeta = this%zeta * this%z

    call device_col3(res%x_d, this%lambda%x_d, this%s%x_d, this%m)

    this%residumax = maxval([ &
         device_maxval(rex%x_d, this%n), &
         device_maxval(rey%x_d, this%m), &
         abs(rez), &
         device_maxval(relambda%x_d, this%m), &
         device_maxval(rexsi%x_d, this%n), &
         device_maxval(reeta%x_d, this%n), &
         device_maxval(remu%x_d, this%m), &
         abs(rezeta), &
         device_maxval(res%x_d, this%m)])

    re_sq_norm = device_norm(rex%x_d, this%n) + &
         device_norm(rexsi%x_d, this%n) + &
         device_norm(reeta%x_d, this%n)

    call MPI_Allreduce(MPI_IN_PLACE, this%residumax, 1, &
         mpi_real_precision, mpi_max, neko_comm, ierr)

    call MPI_Allreduce(MPI_IN_PLACE, re_sq_norm, 1, &
         mpi_real_precision, mpi_sum, neko_comm, ierr)

    this%residunorm = sqrt(( &
         device_norm(rey%x_d, this%m) + &
         rez**2 + &
         device_norm(relambda%x_d, this%m) + &
         device_norm(remu%x_d, this%m) + &
         rezeta**2 + &
         device_norm(res%x_d, this%m) &
         ) + re_sq_norm)

    call this%scratch%relinquish(ind)
  end subroutine mma_dpip_KKT_device

  !============================================================================!
  ! private internal subroutines

  !> generat a subproblem; convex approximation of the optimization problem
  subroutine mma_gensub_device(this, iter, x, df0dx, fval, dfdx)
    ! ----------------------------------------------------- !
    ! Generate the approximation sub problem by computing   !
    ! the lower and upper asymtotes and the other necessary !
    ! parameters (alpha, beta, p0j, q0j, pij, qij, ...).    !
    ! ----------------------------------------------------- !
    class(mma_t), intent(inout) :: this
    type(c_ptr), intent(in) :: x
    type(c_ptr), intent(in) :: df0dx
    type(c_ptr), intent(in) :: fval
    type(c_ptr), intent(in) :: dfdx

    integer, intent(in) :: iter
    integer :: ierr

    type(vector_t), pointer :: x_diff, xmin_eff, xmax_eff
    integer :: ind(3)

    call this%scratch%request(x_diff, ind(1), this%n, .false.)
    call this%scratch%request(xmin_eff, ind(2), this%n, .false.)
    call this%scratch%request(xmax_eff, ind(3), this%n, .false.)

    call device_copy(xmin_eff%x_d, this%xmin%x_d, this%n)
    call device_copy(xmax_eff%x_d, this%xmax%x_d, this%n)

    if (this%move_limit .gt. 0.0_rp) then
       call device_cadd2(xmin_eff%x_d, x, -this%move_limit, this%n)
       call device_pwmax2(xmin_eff%x_d, this%xmin%x_d, this%n)

       call device_cadd2(xmax_eff%x_d, x, this%move_limit, this%n)
       call device_pwmin2(xmax_eff%x_d, this%xmax%x_d, this%n)
    end if

    call device_sub3(x_diff%x_d, xmax_eff%x_d, xmin_eff%x_d, this%n)
    call device_cpwmax2(x_diff%x_d, 1.0e-5_rp, this%n)

    ! ------------------------------------------------------------------------ !
    ! Setup the current asymptotes

    if (iter .lt. 3) then
       call device_copy(this%low%x_d, x, this%n)
       call device_add2s2(this%low%x_d, x_diff%x_d, - this%asyinit, this%n)
       call device_copy(this%upp%x_d, x, this%n)
       call device_add2s2(this%upp%x_d, x_diff%x_d, this%asyinit, this%n)
    else
       call device_mma_gensub2(this%low%x_d, this%upp%x_d, x, &
            this%xold1%x_d, this%xold2%x_d, x_diff%x_d, &
            this%asydecr, this%asyincr, this%n)
    end if

    ! ------------------------------------------------------------------------ !
    ! Calculate p0j, q0j, pij, qij, alpha, and beta

    call device_mma_gensub3(x, df0dx, dfdx, this%low%x_d, &
         this%upp%x_d, xmin_eff%x_d, xmax_eff%x_d, this%alpha%x_d, &
         this%beta%x_d, this%p0j%x_d, this%q0j%x_d, this%pij%x_d, &
         this%qij%x_d, this%n, this%m)

    ! ------------------------------------------------------------------------ !
    ! Computing bi as defined in page 5

    call device_mma_gensub4(x, this%low%x_d, this%upp%x_d, this%pij%x_d, &
         this%qij%x_d, this%n, this%m, this%bi%x_d)

    call device_memcpy(this%bi%x, this%bi%x_d, this%m, DEVICE_TO_HOST, &
         sync = .true.)
    call MPI_Allreduce(MPI_IN_PLACE, this%bi%x, this%m, &
         mpi_real_precision, mpi_sum, neko_comm, ierr)
    call device_memcpy(this%bi%x, this%bi%x_d, this%m, HOST_TO_DEVICE, &
         sync = .true.)

    call device_sub2(this%bi%x_d, fval, this%m)

    call this%scratch%relinquish(ind)
  end subroutine mma_gensub_device

  !> solve the subproblem defined by this%pij, this%qij, etc. using dual-primal
  !! interior point method
  subroutine mma_subsolve_dpip_device(this, designx_d)
    class(mma_t), intent(inout) :: this
    type(c_ptr), intent(in) :: designx_d
    integer :: iter, itto, ierr
    real(kind=rp) :: epsi, residual_max, residual_norm, z, zeta, rez, rezeta, &
         delz, dz, dzeta, steg, zold, zetaold, new_residual
    ! vectors with size m
    type(vector_t) , pointer :: y, lambda, s, mu, rey, relambda, remu, res, &
         dely, dellambda, dy, dlambda, ds, dmu, yold, lambdaold, sold, muold

    ! vectors with size n
    type(vector_t), pointer :: x, xsi, eta, rex, rexsi, reeta, &
         delx, diagx, dx, dxsi, deta, xold, xsiold, etaold

    type(vector_t), pointer :: bb
    type(matrix_t), pointer :: GG
    type(matrix_t), pointer :: AA

    integer :: info
    real(kind=rp) :: re_sq_norm

    integer :: ind(35)

    real(kind=rp) :: minimal_epsilon

    call this%scratch%request(y, ind(1), this%m, .false.)
    call this%scratch%request(lambda, ind(2), this%m, .false.)
    call this%scratch%request(s, ind(3), this%m, .false.)
    call this%scratch%request(mu, ind(4), this%m, .false.)
    call this%scratch%request(rey, ind(5), this%m, .false.)
    call this%scratch%request(relambda, ind(6), this%m, .false.)
    call this%scratch%request(remu, ind(7), this%m, .false.)
    call this%scratch%request(res, ind(8), this%m, .false.)
    call this%scratch%request(dely, ind(9), this%m, .false.)
    call this%scratch%request(dellambda, ind(10), this%m, .false.)
    call this%scratch%request(dy, ind(11), this%m, .false.)
    call this%scratch%request(dlambda, ind(12), this%m, .false.)
    call this%scratch%request(ds, ind(13), this%m, .false.)
    call this%scratch%request(dmu, ind(14), this%m, .false.)
    call this%scratch%request(yold, ind(15), this%m, .false.)
    call this%scratch%request(lambdaold, ind(16), this%m, .false.)
    call this%scratch%request(sold, ind(17), this%m, .false.)
    call this%scratch%request(muold, ind(18), this%m, .false.)
    call this%scratch%request(x, ind(19), this%n, .false.)
    call this%scratch%request(xsi, ind(20), this%n, .false.)
    call this%scratch%request(eta, ind(21), this%n, .false.)
    call this%scratch%request(rex, ind(22), this%n, .false.)
    call this%scratch%request(rexsi, ind(23), this%n, .false.)
    call this%scratch%request(reeta, ind(24), this%n, .false.)
    call this%scratch%request(delx, ind(25), this%n, .false.)
    call this%scratch%request(diagx, ind(26), this%n, .false.)
    call this%scratch%request(dx, ind(27), this%n, .false.)
    call this%scratch%request(dxsi, ind(28), this%n, .false.)
    call this%scratch%request(deta, ind(29), this%n, .false.)
    call this%scratch%request(xold, ind(30), this%n, .false.)
    call this%scratch%request(xsiold, ind(31), this%n, .false.)
    call this%scratch%request(etaold, ind(32), this%n, .false.)
    call this%scratch%request(bb, ind(33), this%m+1, .false.)

    call this%scratch%request(GG, ind(34), this%m, this%n, .false.)
    call this%scratch%request(AA, ind(35), this%m+1, this%m+1, .false.)

    ! ------------------------------------------------------------------------ !
    ! initial value for the parameters in the subsolve based on
    ! page 15 of "https://people.kth.se/~krille/mmagcmma.pdf"

    epsi = 1.0_rp !100
    call device_add3s2(x%x_d, this%alpha%x_d, this%beta%x_d, 0.5_rp, 0.5_rp, &
         this%n)
    call device_cfill(y%x_d, 1.0_rp, this%m)
    z = 1.0_rp
    zeta = 1.0_rp
    call device_cfill(lambda%x_d, 1.0_rp, this%m)
    call device_cfill(s%x_d, 1.0_rp, this%m)
    call device_mma_max(xsi%x_d, x%x_d, this%alpha%x_d, this%n)
    call device_mma_max(eta%x_d, this%beta%x_d, x%x_d, this%n)
    call device_max2(mu%x_d, 1.0_rp, this%c%x_d, 0.5_rp, this%m)

    ! ------------------------------------------------------------------------ !
    ! Computing the minimal epsilon and choose the most conservative one

    minimal_epsilon = max(0.9_rp * this%epsimin, 1.0e-12_rp)
    call MPI_Allreduce(MPI_IN_PLACE, minimal_epsilon, 1, &
         mpi_real_precision, mpi_min, neko_comm, ierr)

    ! ------------------------------------------------------------------------ !
    ! The main loop of the dual-primal interior point method.

    do while (epsi .gt. minimal_epsilon)

       ! --------------------------------------------------------------------- !
       ! Calculating residuals based on
       ! "https://people.kth.se/~krille/mmagcmma.pdf" for the variables
       ! x, y, z, lambda residuals based on eq(5.9a)-(5.9d), respectively.

       associate(p0j => this%p0j, q0j => this%q0j, &
            pij => this%pij, qij => this%qij, &
            low => this%low, upp => this%upp, &
            alpha => this%alpha, beta => this%beta, &
            c => this%c, d => this%d, &
            a0 => this%a0, a => this%a)

         call device_rex(rex%x_d, x%x_d, low%x_d, upp%x_d, &
              pij%x_d, p0j%x_d, qij%x_d, q0j%x_d, &
              lambda%x_d, xsi%x_d, eta%x_d, this%n, this%m)

         call device_col3(rey%x_d, d%x_d, y%x_d, this%m)
         call device_add2(rey%x_d, c%x_d, this%m)
         call device_sub2(rey%x_d, lambda%x_d, this%m)
         call device_sub2(rey%x_d, mu%x_d, this%m)
         rez = a0 - zeta - device_lcsc2(lambda%x_d, a%x_d, this%m)

         call device_cfill(relambda%x_d, 0.0_rp, this%m)
         call device_relambda(relambda%x_d, x%x_d, this%upp%x_d, &
              low%x_d, pij%x_d, qij%x_d, this%n, this%m)

       end associate

       ! --------------------------------------------------------------------- !
       ! Computing the norm of the residuals

       ! Complete the computations of lambda residuals
       call device_memcpy(relambda%x, relambda%x_d, this%m, DEVICE_TO_HOST, &
            sync = .true.)
       call MPI_Allreduce(MPI_IN_PLACE, relambda%x, this%m, &
            mpi_real_precision, mpi_sum, neko_comm, ierr)
       call device_memcpy(relambda%x, relambda%x_d, this%m, HOST_TO_DEVICE, &
            sync = .true.)

       call device_add2s2(relambda%x_d, this%a%x_d, -z, this%m)
       call device_sub2(relambda%x_d, y%x_d, this%m)
       call device_add2(relambda%x_d, s%x_d, this%m)
       call device_sub2(relambda%x_d, this%bi%x_d, this%m)

       call device_sub3(rexsi%x_d, x%x_d, this%alpha%x_d, this%n)
       call device_col2(rexsi%x_d, xsi%x_d, this%n)
       call device_cadd(rexsi%x_d, - epsi, this%n)

       call device_sub3(reeta%x_d, this%beta%x_d, x%x_d, this%n)
       call device_col2(reeta%x_d, eta%x_d, this%n)
       call device_cadd(reeta%x_d, - epsi, this%n)

       call device_col3(remu%x_d, mu%x_d, y%x_d, this%m)
       call device_cadd(remu%x_d, - epsi, this%m)

       rezeta = zeta * z - epsi

       call device_col3(res%x_d, lambda%x_d, s%x_d, this%m)
       call device_cadd(res%x_d, - epsi, this%m)

       ! Setup vectors of residuals and their norms
       residual_max = maxval([device_maxval(rex%x_d, this%n), &
            device_maxval(rey%x_d, this%m), abs(rez), &
            device_maxval(relambda%x_d, this%m), &
            device_maxval(rexsi%x_d, this%n), &
            device_maxval(reeta%x_d, this%n), &
            device_maxval(remu%x_d, this%m), abs(rezeta), &
            device_maxval(res%x_d, this%m)])

       re_sq_norm = device_norm(rex%x_d, this%n) + &
            device_norm(rexsi%x_d, this%n) + device_norm(reeta%x_d, this%n)

       call MPI_Allreduce(MPI_IN_PLACE, residual_max, 1, &
            mpi_real_precision, mpi_max, neko_comm, ierr)

       call MPI_Allreduce(MPI_IN_PLACE, re_sq_norm, &
            1, mpi_real_precision, mpi_sum, neko_comm, ierr)

       residual_norm = sqrt(device_norm(rey%x_d, this%m) + &
            rez**2 + &
            device_norm(relambda%x_d, this%m) + &
            device_norm(remu%x_d, this%m)+ &
            rezeta**2 + &
            device_norm(res%x_d, this%m) &
            + re_sq_norm)

       ! --------------------------------------------------------------------- !
       ! Internal loop

       do iter = 1, this%max_iter

          if (residual_max .lt. epsi) exit

          call device_delx(delx%x_d, x%x_d, this%low%x_d, this%upp%x_d, &
               this%pij%x_d, this%qij%x_d, this%p0j%x_d, this%q0j%x_d, &
               this%alpha%x_d, this%beta%x_d, lambda%x_d, epsi, this%n, &
               this%m)

          call device_col3(dely%x_d, this%d%x_d, y%x_d, this%m)
          call device_add2(dely%x_d, this%c%x_d, this%m)
          call device_sub2(dely%x_d, lambda%x_d, this%m)
          call device_add2inv2(dely%x_d, y%x_d, - epsi, this%m)
          delz = this%a0 - device_lcsc2(lambda%x_d, this%a%x_d, this%m) - epsi/z

          ! Accumulate sums for dellambda (the term gi(x))
          call device_cfill(dellambda%x_d, 0.0_rp, this%m)
          call device_relambda(dellambda%x_d, x%x_d, this%upp%x_d, &
               this%low%x_d, this%pij%x_d, this%qij%x_d, this%n, this%m)

          call device_memcpy(dellambda%x, dellambda%x_d, this%m, &
               DEVICE_TO_HOST, sync = .true.)
          call MPI_Allreduce(MPI_IN_PLACE, dellambda%x, this%m, &
               mpi_real_precision, mpi_sum, neko_comm, ierr)
          call device_memcpy(dellambda%x, dellambda%x_d, this%m, &
               HOST_TO_DEVICE, sync = .true.)

          call device_add3s2(dellambda%x_d, dellambda%x_d, this%a%x_d, &
               1.0_rp, -z, this%m)
          call device_sub2(dellambda%x_d, y%x_d, this%m)
          call device_sub2(dellambda%x_d, this%bi%x_d, this%m)
          call device_add2inv2(dellambda%x_d, lambda%x_d, epsi, this%m)

          call device_GG(GG%x_d, x%x_d, this%low%x_d, this%upp%x_d, &
               this%pij%x_d, this%qij%x_d, this%n, this%m)

          call device_diagx(diagx%x_d, x%x_d, xsi%x_d, this%low%x_d, &
               this%upp%x_d, this%p0j%x_d, this%q0j%x_d, this%pij%x_d, &
               this%qij%x_d, this%alpha%x_d, this%beta%x_d, eta%x_d, &
               lambda%x_d, this%n, this%m)

          !Here we only consider the case m<n in the matlab code
          !assembling the right hand side matrix based on eq(5.20)
          ! bb = [dellambda + dely/(this%d%x + &
          !         (mu/y)) - matmul(GG,delx/diagx), delz ]

          !--------------------------------------------------------------------!
          ! for MPI computation of bb

          call device_bb(bb%x_d, GG%x_d, delx%x_d, diagx%x_d, this%n, &
               this%m)

          call device_memcpy(bb%x, bb%x_d, this%m + 1, DEVICE_TO_HOST, &
               sync = .true.)
          call MPI_Allreduce(MPI_IN_PLACE, bb%x, this%m + 1, &
               mpi_real_precision, mpi_sum, neko_comm, ierr)
          call device_memcpy(bb%x, bb%x_d, this%m + 1, &
               HOST_TO_DEVICE, sync = .true.)

          call device_updatebb(bb%x_d, dellambda%x_d, dely%x_d, &
               this%d%x_d, mu%x_d, y%x_d, delz, this%m)

          ! assembling the coefficients matrix AA based on eq(5.20)
          ! AA(1:this%m,1:this%m) =  &
          ! matmul(matmul(GG,mma_diag(1/diagx)), transpose(GG))
          ! !update diag(AA)
          ! AA(1:this%m,1:this%m) = AA(1:this%m,1:this%m) + &
          !     mma_diag(s/lambda + 1.0/(this%d%x + (mu/y)))

          call device_cfill(AA%x_d, 0.0_rp, (this%m+1) * (this%m+1))
          call device_AA(AA%x_d, GG%x_d, diagx%x_d, this%n, this%m)

          call device_memcpy(AA%x, AA%x_d, (this%m+1) * (this%m+1), &
               DEVICE_TO_HOST, sync = .true.)
          call MPI_Allreduce(MPI_IN_PLACE, AA%x, &
               (this%m + 1)**2, mpi_real_precision, mpi_sum, neko_comm, ierr)
          call device_memcpy(AA%x, AA%x_d, (this%m+1) * (this%m+1), &
               HOST_TO_DEVICE, sync = .true.)

          call device_prepare_aa_matrix(AA%x_d, s%x_d, lambda%x_d, &
               this%d%x_d, mu%x_d, y%x_d, this%a%x_d, zeta, z, this%m)

          ! Device solve for the linear system
          call device_solve_linear_system(AA%x_d, bb%x_d, this%m + 1, info)
          if (info .ne. 0) then
             call neko_error("Linear solver failed on the device in  " // &
                  "mma_subsolve_dpip")
          end if

          call device_copy(dlambda%x_d, bb%x_d, this%m)


          !We need to write the last element of bb to dz so this is necessary
          call device_memcpy(bb%x, bb%x_d, this%m+1, DEVICE_TO_HOST, &
               sync = .true.)
          dz = bb%x(this%m + 1)


          ! based on eq(5.19)
          call device_dx(dx%x_d, delx%x_d, diagx%x_d, GG%x_d, &
               dlambda%x_d, this%n, this%m)
          call device_dy(dy%x_d, dely%x_d, dlambda%x_d, this%d%x_d, &
               mu%x_d, y%x_d, this%m)
          call device_dxsi(dxsi%x_d, xsi%x_d, dx%x_d, x%x_d, &
               this%alpha%x_d, epsi, this%n)
          call device_deta(deta%x_d, eta%x_d, dx%x_d, x%x_d, &
               this%beta%x_d, epsi, this%n)

          call device_col3(dmu%x_d, mu%x_d, dy%x_d, this%m)
          call device_cmult(dmu%x_d, -1.0_rp, this%m)
          call device_cadd(dmu%x_d, epsi, this%m)
          call device_invcol2(dmu%x_d, y%x_d, this%m)
          call device_sub2(dmu%x_d, mu%x_d, this%m)
          dzeta = -zeta + (epsi - zeta * dz) / z
          call device_col3(ds%x_d, dlambda%x_d, s%x_d, this%m)
          call device_cmult(ds%x_d, -1.0_rp, this%m)
          call device_cadd(ds%x_d, epsi, this%m)
          call device_invcol2(ds%x_d, lambda%x_d, this%m)
          call device_sub2(ds%x_d, s%x_d, this%m)

          steg = maxval([1.0_rp, &
               device_maxval2(dy%x_d, y%x_d, -1.01_rp, this%m), &
               -1.01_rp * dz / z, &
               device_maxval2(dlambda%x_d, lambda%x_d, -1.01_rp, this%m), &
               device_maxval2(dxsi%x_d, xsi%x_d, -1.01_rp, this%n), &
               device_maxval2(deta%x_d, eta%x_d, -1.01_rp, this%n), &
               device_maxval2(dmu%x_d, mu%x_d, -1.01_rp, this%m), &
               -1.01_rp * dzeta / zeta, &
               device_maxval2(ds%x_d, s%x_d, -1.01_rp, this%m), &
               device_maxval3(dx%x_d, x%x_d, this%alpha%x_d, -1.01_rp, this%n),&
               device_maxval3(dx%x_d, this%beta%x_d, x%x_d, 1.01_rp, this%n)])

          steg = 1.0_rp / steg

          call device_copy(xold%x_d, x%x_d, this%n)
          call device_copy(yold%x_d, y%x_d, this%m)
          zold = z
          call device_copy(lambdaold%x_d, lambda%x_d, this%m)
          call device_copy(xsiold%x_d, xsi%x_d, this%n)
          call device_copy(etaold%x_d, eta%x_d, this%n)
          call device_copy(muold%x_d, mu%x_d, this%m)
          zetaold = zeta
          call device_copy(sold%x_d, s%x_d, this%m)

          new_residual = 2.0_rp * residual_norm

          ! Share the new_residual and steg values
          call MPI_Allreduce(MPI_IN_PLACE, steg, 1, &
               mpi_real_precision, mpi_min, neko_comm, ierr)
          call MPI_Allreduce(MPI_IN_PLACE, new_residual, 1, &
               mpi_real_precision, mpi_min, neko_comm, ierr)

          ! The innermost loop to determine the suitable step length
          ! using the Backtracking Line Search approach
          itto = 0
          do while ((new_residual .gt. residual_norm) .and. (itto .lt. 50))
             itto = itto + 1

             ! update the variables
             call device_add3s2(x%x_d, xold%x_d, dx%x_d, 1.0_rp, steg, this%n)
             call device_add3s2(y%x_d, yold%x_d, dy%x_d, 1.0_rp, steg, this%m)
             z = zold + steg*dz
             call device_add3s2(lambda%x_d, lambdaold%x_d, &
                  dlambda%x_d, 1.0_rp, steg, this%m)
             call device_add3s2(xsi%x_d, xsiold%x_d, dxsi%x_d, &
                  1.0_rp, steg, this%n)
             call device_add3s2(eta%x_d, etaold%x_d, deta%x_d, &
                  1.0_rp, steg, this%n)
             call device_add3s2(mu%x_d, muold%x_d, dmu%x_d, &
                  1.0_rp, steg, this%m)
             zeta = zetaold + steg*dzeta
             call device_add3s2(s%x_d, sold%x_d, ds%x_d, 1.0_rp, &
                  steg, this%m)

             ! Recompute the new_residual to see if this stepsize improves
             ! the residue
             call device_rex(rex%x_d, x%x_d, this%low%x_d, &
                  this%upp%x_d, this%pij%x_d, this%p0j%x_d, &
                  this%qij%x_d, this%q0j%x_d, lambda%x_d, xsi%x_d, &
                  eta%x_d, this%n, this%m)

             call device_col3(rey%x_d, this%d%x_d, y%x_d, this%m)
             call device_add2(rey%x_d, this%c%x_d, this%m)
             call device_sub2(rey%x_d, lambda%x_d, this%m)
             call device_sub2(rey%x_d, mu%x_d, this%m)

             rez = this%a0 - zeta - device_lcsc2(lambda%x_d, this%a%x_d, this%m)

             ! Accumulate sums for relambda (the term gi(x))
             call device_cfill(relambda%x_d, 0.0_rp, this%m)
             call device_relambda(relambda%x_d, x%x_d, this%upp%x_d, &
                  this%low%x_d, this%pij%x_d, this%qij%x_d, &
                  this%n, this%m)

             call device_memcpy(relambda%x, relambda%x_d, this%m, &
                  DEVICE_TO_HOST, sync = .true.)
             call MPI_Allreduce(MPI_IN_PLACE, relambda%x, this%m, &
                  mpi_real_precision, mpi_sum, neko_comm, ierr)
             call device_memcpy(relambda%x, relambda%x_d, &
                  this%m, HOST_TO_DEVICE, sync = .true.)

             call device_add3s2(relambda%x_d, relambda%x_d, &
                  this%a%x_d, 1.0_rp, -z, this%m)
             call device_sub2(relambda%x_d, y%x_d, this%m)
             call device_add2(relambda%x_d, s%x_d, this%m)
             call device_sub2(relambda%x_d, this%bi%x_d, this%m)

             call device_sub3(rexsi%x_d, x%x_d, this%alpha%x_d, this%n)
             call device_col2(rexsi%x_d, xsi%x_d, this%n)
             call device_cadd(rexsi%x_d, - epsi, this%n)

             call device_sub3(reeta%x_d, this%beta%x_d, x%x_d, this%n)
             call device_col2(reeta%x_d, eta%x_d, this%n)
             call device_cadd(reeta%x_d, - epsi, this%n)

             call device_col3(remu%x_d, mu%x_d, y%x_d, this%m)
             call device_cadd(remu%x_d, - epsi, this%m)

             rezeta = zeta*z - epsi

             call device_col3(res%x_d, lambda%x_d, s%x_d, this%m)
             call device_cadd(res%x_d, - epsi, this%m)

             ! Compute squared norms for the residuals
             re_sq_norm = device_norm(rex%x_d, this%n) + &
                  device_norm(rexsi%x_d, this%n) + &
                  device_norm(reeta%x_d, this%n)
             call MPI_Allreduce(MPI_IN_PLACE, re_sq_norm, 1, &
                  mpi_real_precision, mpi_sum, neko_comm, ierr)

             new_residual = sqrt(device_norm(rey%x_d, this%m) + &
                  rez**2 + &
                  device_norm(relambda%x_d, this%m) + &
                  device_norm(remu%x_d, this%m) + &
                  rezeta**2 + &
                  device_norm(res%x_d, this%m) + &
                  re_sq_norm)

             call MPI_Allreduce(MPI_IN_PLACE, new_residual, 1, &
                  mpi_real_precision, mpi_sum, neko_comm, ierr)

             steg = steg / 2.0_rp

          end do
          steg = 2.0_rp * steg ! Correction for the final division by 2

          ! Update the maximum and norm of the residuals
          residual_norm = new_residual
          residual_max = maxval([ &
               device_maxval(rex%x_d, this%n), &
               device_maxval(rey%x_d, this%m), &
               abs(rez), &
               device_maxval(relambda%x_d, this%m), &
               device_maxval(rexsi%x_d, this%n), &
               device_maxval(reeta%x_d, this%n), &
               device_maxval(remu%x_d, this%m), &
               abs(rezeta), &
               device_maxval(res%x_d, this%m)])

          call MPI_Allreduce(MPI_IN_PLACE, residual_max, 1, &
               mpi_real_precision, mpi_max, neko_comm, ierr)

       end do

       epsi = 0.1_rp * epsi
    end do

    ! Save the new designx
    call device_copy(this%xold2%x_d, this%xold1%x_d, this%n)
    call device_copy(this%xold1%x_d, designx_d, this%n)
    call device_copy(designx_d, x%x_d, this%n)

    ! update the parameters of the MMA object nesessary to compute KKT residual
    call device_copy(this%y%x_d, y%x_d, this%m)
    this%z = z
    call device_copy(this%lambda%x_d, lambda%x_d, this%m)
    this%zeta = zeta
    call device_copy(this%xsi%x_d, xsi%x_d, this%n)
    call device_copy(this%eta%x_d, eta%x_d, this%n)
    call device_copy(this%mu%x_d, mu%x_d, this%m)
    call device_copy(this%s%x_d, s%x_d, this%m)

    !free all the initiated variables in this subroutine
    call this%scratch%relinquish(ind)
  end subroutine mma_subsolve_dpip_device

  !> solve the subproblem defined by this%pij, this%qij, etc. using dual
  !! interior point method
  subroutine mma_subsolve_dip_device(this, designx_d)
    class(mma_t), intent(inout) :: this
    type(c_ptr), intent(in) :: designx_d
    integer :: iter, ierr
    real(kind=rp) :: epsi, residumax, z, steg
    ! vectors with size m
    type(vector_t), pointer :: y, lambda, mu, relambda, remu, dlambda, dmu, &
         gradlambda, zerom, dd, dummy_m
    ! vectors with size n
    type(vector_t), pointer :: x, pjlambda, qjlambda

    ! inverse of a diag matrix:
    type(vector_t), pointer :: Ljjxinv ! [∇_x^2 Ljj]−1
    type(matrix_t), pointer :: hijx ! ∇_x hij
    type(matrix_t), pointer :: Hess

    integer :: info, ind(17)

    real(kind=rp) :: minimal_epsilon

    call this%scratch%request(y, ind(1), this%m, .false.)
    call this%scratch%request(lambda, ind(2), this%m, .false.)
    call this%scratch%request(mu, ind(3), this%m, .false.)
    call this%scratch%request(relambda, ind(4), this%m, .false.)
    call this%scratch%request(remu, ind(5), this%m, .false.)
    call this%scratch%request(dlambda, ind(6), this%m, .false.)
    call this%scratch%request(dmu, ind(7), this%m, .false.)
    call this%scratch%request(gradlambda, ind(8), this%m, .false.)
    call this%scratch%request(zerom, ind(9), this%m, .false.)
    call this%scratch%request(dd, ind(10), this%m, .false.)
    call this%scratch%request(dummy_m, ind(11), this%m, .false.)

    call this%scratch%request(x, ind(12), this%n, .false.)
    call this%scratch%request(pjlambda,ind(13), this%n, .false.)
    call this%scratch%request(qjlambda, ind(14), this%n, .false.)

    call this%scratch%request(Ljjxinv, ind(15), this%n, .false.)

    call this%scratch%request(hijx, ind(16), this%m, this%n, .false.)
    call this%scratch%request(Hess, ind(17), this%m, this%m, .false.)

    ! ------------------------------------------------------------------------ !
    ! initial value for the parameters in the subsolve based on
    ! page 15 of "https://people.kth.se/~krille/mmagcmma.pdf"

    epsi = 1.0_rp !100
    call device_cfill(y%x_d, 1.0_rp, this%m)
    ! initialize lambda with an array of ones (change to this%c%x/2 if needed!)
    call device_cfill(lambda%x_d, 1.0_rp, this%m)
    call device_cmult2(dummy_m%x_d, this%c%x_d, 0.5_rp, this%m)
    call device_pwmax2(lambda%x_d, dummy_m%x_d, this%m)

    call device_cfill(mu%x_d, 1.0_rp, this%m)
    z = 0.0_rp

    ! ------------------------------------------------------------------------ !
    ! Computing the minimal epsilon and choose the most conservative one

    minimal_epsilon = max(0.9_rp * this%epsimin, 1.0e-12_rp)
    call MPI_Allreduce(MPI_IN_PLACE, minimal_epsilon, 1, &
         mpi_real_precision, mpi_min, neko_comm, ierr)

    ! ------------------------------------------------------------------------ !
    ! The main loop of the dual-primal interior point method.

    outer: do while (epsi .gt. minimal_epsilon)
       ! calculating residuals based on
       ! "https://people.kth.se/~krille/mmagcmma.pdf" for the variables
       ! x, y, z, lambda residuals based on eq(5.9a)-(5.9d), respectively.
       associate(p0j => this%p0j, q0j => this%q0j, &
            pij => this%pij, qij => this%qij, &
            low => this%low, upp => this%upp, &
            alpha => this%alpha, beta => this%beta, &
            c => this%c, a0 => this%a0, a => this%a)

         ! minimize(L_x, L_y, L_z) and compute x(λ), y(λ), z(λ) for
         ! the initial value of λ

         ! Comput the value of y that minimizes L_y for the current λ
         ! minimize (sum_{i=1}^{m} [ (c_i - λ_i) * y_i + 0.5 * y_i^2 ])
         ! dL_y/dy =0   => y= (λ_i - c_i), ensure y>=0
         call device_sub3(y%x_d, lambda%x_d, c%x_d, this%m)
         call device_pwmax2(y%x_d, zerom%x_d, this%m)

         ! Comput the value of z that minimizes L_z for the current λ
         ! minimize ((a_0 - sum_{i=1}^{m} λ_i * a_i) * z + 0.5 * z^2)
         ! ensure z>=0
         call device_col3(dummy_m%x_d, lambda%x_d, a%x_d, this%m)
         z = device_glsum(dummy_m%x_d, this%m)
         z = max(0.0_rp, z - a0)

         ! Comput the value of x that minimizes L_x for the current λ
         ! minimize( sum_{j=1}^{n} [ (p_{0j} + sum_{i=1}^{m} λ_i *
         ! p_{ij}) / (u_j - x_j) + (q_{0j} + sum_{i=1}^{m} λ_i * q_{ij}) /
         ! (x_j - l_j) ] - sum_{i=1}^{m} λ_i * b_i)
         call device_mattrans_v_mul(pjlambda%x_d, pij%x_d, lambda%x_d, this%m, this%n)
         call device_mattrans_v_mul(qjlambda%x_d, qij%x_d, lambda%x_d, this%m, this%n)
         call device_add2(pjlambda%x_d, p0j%x_d, this%n)
         call device_add2(qjlambda%x_d, q0j%x_d, this%n)

         call device_mma_dipsolvesub1(x%x_d, pjlambda%x_d, qjlambda%x_d, &
              low%x_d, upp%x_d, alpha%x_d, beta%x_d, this%n)

         call device_cfill(relambda%x_d, 0.0_rp, this%m)
         call device_relambda(relambda%x_d, x%x_d, this%upp%x_d, &
              low%x_d, pij%x_d, qij%x_d, this%n, this%m)

         ! Global comminucation for relambda values

         call device_memcpy(relambda%x, relambda%x_d, this%m, DEVICE_TO_HOST, &
              sync = .true.)
         call MPI_Allreduce(MPI_IN_PLACE, relambda%x, this%m, &
              mpi_real_precision, mpi_sum, neko_comm, ierr)
         call device_memcpy(relambda%x, relambda%x_d, this%m, &
              HOST_TO_DEVICE, sync = .true.)

         call device_add2s2(relambda%x_d, this%a%x_d, -z, this%m)
         call device_sub2(relambda%x_d, y%x_d, this%m)
         call device_add2(relambda%x_d, mu%x_d, this%m)
         call device_sub2(relambda%x_d, this%bi%x_d, this%m)

         call device_col3(remu%x_d, mu%x_d, lambda%x_d, this%m)
         call device_cadd(remu%x_d, -epsi, this%m)

         residumax = maxval([device_maxval(relambda%x_d, this%m), &
              device_maxval(remu%x_d, this%m)])

         ! ------------------------------------------------------------------- !
         ! Internal loop
         do iter = 1, this%max_iter
            !Check the condition
            if (residumax .lt. epsi) exit

            ! Compute dL(x, y, z, λ)/dλ for the updated x(λ), y(λ), z(λ)
            ! based on the implementation in the following paper by Niels
            ! https://doi.org/10.1007/s00158-012-0869-2
            ! (https://github.com/topopt/TopOpt_in_PETSc/blob/master/MMA.cc)
            ! The formula for gradlambda and relambda are basically the same:
            ! thus, we utilise gradlambda = relambda - mu for efficiency
            call device_copy(gradlambda%x_d, relambda%x_d, this%m)
            call device_sub2(gradlambda%x_d, mu%x_d, this%m)

            ! Update gradlambda as the right hand side for Newton's method(eq10)
            call device_cfill(dummy_m%x_d, epsi, this%m)
            call device_invcol2(dummy_m%x_d, lambda%x_d, this%m)
            call device_add2(gradlambda%x_d, dummy_m%x_d, this%m)
            call device_cmult(gradlambda%x_d, -1.0_rp, this%m)

            ! Computing the Hessian as in equation (13) in
            !! https://doi.org/10.1007/s00158-012-0869-2

            !--------------contributions of x terms to Hess--------------------!
            call device_mma_Ljjxinv(Ljjxinv%x_d, pjlambda%x_d, qjlambda%x_d, &
                 x%x_d, low%x_d, upp%x_d, alpha%x_d, beta%x_d, this%n)

            call device_GG(hijx%x_d, x%x_d, this%low%x_d, this%upp%x_d, &
                 this%pij%x_d, this%qij%x_d, this%n, this%m)

            call device_cfill(Hess%x_d, 0.0_rp, (this%m) * (this%m) )
            call device_Hess(Hess%x_d, hijx%x_d, Ljjxinv%x_d, this%n, this%m)

            ! download Hess to CPU, mpi reduce, upload to the device
            call device_memcpy(Hess%x, Hess%x_d, this%m*this%m, DEVICE_TO_HOST, &
                 sync = .true.)
            call MPI_Allreduce(MPI_IN_PLACE, Hess%x, &
                 this%m*this%m, mpi_real_precision, mpi_sum, neko_comm, ierr)
            call device_memcpy(Hess%x, Hess%x_d, this%m*this%m, &
                 HOST_TO_DEVICE, sync = .true.)

            !---------------contributions of z terms to Hess-------------------!
            ! Only for inactive constraint, we consider contributions to Hess
            ! based on the cpp code by Niels.
            call device_col3(dummy_m%x_d, lambda%x_d, a%x_d, this%m)
            if (device_glsum(dummy_m%x_d, this%m) .gt. 0.0_rp) then
               call device_update_hessian_z(Hess%x_d, a%x_d, this%m)
            end if

            !---------------contributions of y terms to Hess-------------------!
            ! Only for inactive constraint, we consider contributions to Hess.
            ! Note that if d(i) = 0, the y terms (just like z terms) will not
            ! contribute to the Hessian matrix.
            ! Note that since we use DGESV to solve LSE on CPU, we dont need
            ! cuda kernel for this part
            ! Also, improve the robustness by stablizing the Hess using
            ! Levenberg-Marquardt algorithm (heuristically)
            call device_prepare_hessian(Hess%x_d, y%x_d, mu%x_d, lambda%x_d, &
                 this%m)

            ! Device solve for the linear system
            call device_solve_linear_system(Hess%x_d, gradlambda%x_d, &
                 this%m, info)
            if (info .ne. 0) then
               call neko_error("Linear solver failed on the device in  " // &
                    "mma_subsolve_dip")
            end if

            call device_copy(dlambda%x_d, gradlambda%x_d, this%m)

            ! based on eq(11) for delta eta
            call device_copy(dummy_m%x_d, dlambda%x_d, this%m)
            call device_col2(dummy_m%x_d, mu%x_d, this%m)
            call device_invcol2(dummy_m%x_d, lambda%x_d, this%m)

            call device_cfill(dmu%x_d, epsi, this%m)
            call device_invcol2(dmu%x_d, lambda%x_d, this%m)
            call device_add2s2(dmu%x_d, dummy_m%x_d, -1.0_rp, this%m)
            call device_sub2(dmu%x_d, mu%x_d, this%m)

            steg = maxval([1.005_rp, device_maxval2(dlambda%x_d, lambda%x_d, &
                 -1.01_rp, this%m), device_maxval2(dmu%x_d, mu%x_d, -1.01_rp, &
                 this%m)])
            steg = 1.0_rp / steg

            call device_add2s2(lambda%x_d, dlambda%x_d, steg, this%m)
            call device_add2s2(mu%x_d, dmu%x_d, steg, this%m)

            ! minimize(L_x, L_y, L_z) and compute x(λ), y(λ), z(λ) for
            ! the updated values of λ

            ! Comput the value of y that minimizes L_y for the current λ
            ! minimize (sum_{i=1}^{m} [ (c_i - λ_i) * y_i + 0.5 * y_i^2 ])
            ! dL_y/dy =0   => y= (λ_i - c_i), ensure y>=0
            call device_sub3(y%x_d, lambda%x_d, c%x_d, this%m)
            call device_pwmax2(y%x_d, zerom%x_d, this%m)

            ! Comput the value of z that minimizes L_z for the current λ
            ! minimize ((a_0 - sum_{i=1}^{m} λ_i * a_i) * z + 0.5 * z^2)
            ! ensure z>=0
            call device_col3(dummy_m%x_d, lambda%x_d, a%x_d, this%m)
            z = device_glsum(dummy_m%x_d, this%m)
            z = max(0.0_rp, z - a0)

            ! Comput the value of x that minimizes L_x for the current λ
            ! minimize( sum_{j=1}^{n} [ (p_{0j} + sum_{i=1}^{m} λ_i *
            ! p_{ij}) / (u_j - x_j) + (q_{0j} + sum_{i=1}^{m} λ_i * q_{ij}) /
            ! (x_j - l_j) ] - sum_{i=1}^{m} λ_i * b_i)
            call device_mattrans_v_mul(pjlambda%x_d, pij%x_d, lambda%x_d, this%m, this%n)
            call device_mattrans_v_mul(qjlambda%x_d, qij%x_d, lambda%x_d, this%m, this%n)
            call device_add2(pjlambda%x_d, p0j%x_d, this%n)
            call device_add2(qjlambda%x_d, q0j%x_d, this%n)

            call device_mma_dipsolvesub1(x%x_d, pjlambda%x_d, qjlambda%x_d, &
                 low%x_d, upp%x_d, alpha%x_d, beta%x_d, this%n)

            ! Compute the residual for the lambda and mu using eq(9) and eq(15)

            call device_cfill(relambda%x_d, 0.0_rp, this%m)
            call device_relambda(relambda%x_d, x%x_d, this%upp%x_d, &
                 low%x_d, pij%x_d, qij%x_d, this%n, this%m)

            ! Global comminucation for relambda values

            call device_memcpy(relambda%x, relambda%x_d, this%m, DEVICE_TO_HOST, &
                 sync = .true.)
            call MPI_Allreduce(MPI_IN_PLACE, relambda%x, this%m, &
                 mpi_real_precision, mpi_sum, neko_comm, ierr)
            call device_memcpy(relambda%x, relambda%x_d, this%m, &
                 HOST_TO_DEVICE, sync = .true.)

            call device_add2s2(relambda%x_d, this%a%x_d, -z, this%m)
            call device_sub2(relambda%x_d, y%x_d, this%m)
            call device_add2(relambda%x_d, mu%x_d, this%m)
            call device_sub2(relambda%x_d, this%bi%x_d, this%m)

            call device_col3(remu%x_d, mu%x_d, lambda%x_d, this%m)
            call device_cadd(remu%x_d, -epsi, this%m)

            residumax = maxval([device_maxval(relambda%x_d, this%m), &
                 device_maxval(remu%x_d, this%m)])
         end do
       end associate
       epsi = 0.1_rp * epsi
    end do outer

    ! Save the new designx
    call device_copy(this%xold2%x_d, this%xold1%x_d, this%n)
    call device_copy(this%xold1%x_d, designx_d, this%n)
    call device_copy(designx_d, x%x_d, this%n)

    ! update the parameters of the MMA object nesessary to compute KKT residual
    call device_copy(this%y%x_d, y%x_d, this%m)
    this%z = z
    call device_copy(this%lambda%x_d, lambda%x_d, this%m)
    call device_copy(this%mu%x_d, mu%x_d, this%m)

    call this%scratch%relinquish(ind)
  end subroutine mma_subsolve_dip_device

  subroutine mma_subsolve_unconstrained_device(this, designx_d)
    class(mma_t), intent(inout) :: this
    type(c_ptr), intent(in) :: designx_d
    type(vector_t), pointer :: x
    integer :: ind

    call this%scratch%request(x, ind, this%n, .false.)
    associate(p0j => this%p0j, q0j => this%q0j, &
         low => this%low, upp => this%upp, &
         alpha => this%alpha, beta => this%beta)
      call device_mma_dipsolvesub1(x%x_d, p0j%x_d, q0j%x_d, &
           low%x_d, upp%x_d, alpha%x_d, beta%x_d, this%n)
      ! Closed-form primal solution for unconstrained MMA subproblem

    end associate

    ! Save the new designx
    call device_copy(this%xold2%x_d, this%xold1%x_d, this%n)
    call device_copy(this%xold1%x_d, designx_d, this%n)
    call device_copy(designx_d, x%x_d, this%n)
    call this%scratch%relinquish(ind)
  end subroutine mma_subsolve_unconstrained_device

end submodule mma_device
