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

submodule (mma) mma_device

  use device_math, only: device_copy, device_cmult, device_cadd, device_cfill, &
       device_add2, device_add3s2, device_invcol2, device_col2, device_col3, &
       device_sub2, device_sub3, device_add2s2, device_cadd2, device_pwmax, &
       device_glsum, device_cmult2, device_pwmax
  use device_mma_math, only: device_maxval, device_norm, device_lcsc2, &
       device_maxval2, device_maxval3, device_mma_gensub3, &
       device_mma_gensub4, device_mma_max, device_max2, device_rex, &
       device_relambda, device_delx, device_add2inv2, device_gg, device_diagx, &
       device_bb, device_updatebb, device_aa, device_updateaa, device_dx, &
       device_dy, device_dxsi, device_deta, device_kkt_rex, &
       device_mma_gensub2, device_mattrans_v_mul, device_mma_dipsolvesub1, &
       device_mma_Ljjxinv, device_Hess

  use neko_config, only: NEKO_BCKND_DEVICE
  use device, only: DEVICE_TO_HOST
  use comm, only: neko_comm, pe_rank, mpi_real_precision
  use mpi_f08, only: MPI_IN_PLACE, MPI_MAX, MPI_MIN

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
    real(kind=rp), dimension(this%n), intent(inout) :: x

    type(vector_t) :: df0dx, fval, xdesign
    type(matrix_t) :: dfdx

    if (.not. this%is_initialized) then
       call neko_error("The MMA object is not initialized.")
    end if

    call xdesign%init(this%n)
    xdesign%x = x
    call device_memcpy(xdesign%x, xdesign%x_d, this%n, HOST_TO_DEVICE, sync = .true.)

     ! call test_gensub_cpu(this, iter, x, df0dx, fval, dfdx)
     ! print *, "on cpu:"
     ! print *, "p0j=", sum(this%p0j%x), &
     !      "q0j=", sum(this%q0j%x), &
     !      "pij=", sum(this%pij%x), &
     !      "qij=", sum(this%qij%x), &
     !      "bi=", sum(this%bi%x), &
     !      "alpha=", sum(this%alpha%x), &
     !      "beta=", sum(this%beta%x), &
     !      "upp=", sum(this%upp%x), &
     !      "low=", sum(this%low%x)

    ! generate a convex approximation of the problem
         call mma_gensub_device(this, iter, xdesign, df0dx, fval, dfdx)
          call device_memcpy(this%p0j%x, this%p0j%x_d, this%n, &
               DEVICE_TO_HOST, sync = .true.)
          call device_memcpy(this%q0j%x, this%q0j%x_d, this%n, &
               DEVICE_TO_HOST, sync = .true.)
          call device_memcpy(this%pij%x, this%pij%x_d, this%n*this%m, &
               DEVICE_TO_HOST, sync = .true.)
          call device_memcpy(this%qij%x, this%qij%x_d, this%n*this%m, &
               DEVICE_TO_HOST, sync = .true.)
          call device_memcpy(this%bi%x, this%bi%x_d, this%m, &
               DEVICE_TO_HOST, sync = .true.)
          call device_memcpy(this%alpha%x, this%alpha%x_d, this%n, &
               DEVICE_TO_HOST, sync = .true.)
          call device_memcpy(this%beta%x, this%beta%x_d, this%n, &
               DEVICE_TO_HOST, sync = .true.)
          call device_memcpy(this%upp%x, this%upp%x_d, this%n, &
               DEVICE_TO_HOST, sync = .true.)
          call device_memcpy(this%low%x, this%low%x_d, this%n, &
               DEVICE_TO_HOST, sync = .true.)
     ! print *, "on cuda:"
     ! print *, "p0j=", sum(this%p0j%x), &
     !      "q0j=", sum(this%q0j%x), &
     !      "pij=", sum(this%pij%x), &
     !      "qij=", sum(this%qij%x), &
     !      "bi=", sum(this%bi%x), &
     !      "alpha=", sum(this%alpha%x), &
     !      "beta=", sum(this%beta%x), &
     !      "upp=", sum(this%upp%x), &
     !      "low=", sum(this%low%x)



    !solve the approximation problem using interior point method
    if (this%subsolver .eq. "dip") then
       call mma_subsolve_dip_device(this, xdesign)
    else
       call mma_subsolve_dpip_device(this, xdesign)
     !   call test_subsolve_dpip_cpu(this, x)
    end if
     ! call device_memcpy(this%xold1%x, this%xold1%x_d, this%n, &
     !      HOST_TO_DEVICE, sync = .true.)
     ! call device_memcpy(this%xold2%x, this%xold2%x_d, this%n, &
     !      HOST_TO_DEVICE, sync = .true.)

     ! call device_memcpy(this%y%x, this%y%x_d, this%m, &
     !      HOST_TO_DEVICE, sync = .true.)
     ! call device_memcpy(this%lambda%x, this%lambda%x_d, this%m, &
     !      HOST_TO_DEVICE, sync = .true.)
     ! call device_memcpy(this%xsi%x, this%xsi%x_d, this%n, &
     !      HOST_TO_DEVICE, sync = .true.)
     ! call device_memcpy(this%eta%x, this%eta%x_d, this%n, &
     !      HOST_TO_DEVICE, sync = .true.)
     ! call device_memcpy(this%mu%x, this%mu%x_d, this%m, &
     !      HOST_TO_DEVICE, sync = .true.)
     ! call device_memcpy(this%s%x, this%s%x_d, this%m, &
     !      HOST_TO_DEVICE, sync = .true.)



    !update the design vector x on the host
    call device_memcpy(x, xdesign%x_d, this%n, DEVICE_TO_HOST, sync = .false.)

    this%is_updated = .true.
    call xdesign%free()
  end subroutine mma_update_device

  module subroutine mma_KKT_device(this, x, df0dx, fval, dfdx)
    class(mma_t), intent(inout) :: this
    real(kind=rp), dimension(this%n), intent(in) :: x
    type(vector_t), intent(in) :: fval, df0dx
    type(matrix_t), intent(in) :: dfdx

    if (this%subsolver .eq. "dip") then
       call mma_dip_KKT_device(this, x, df0dx, fval, dfdx)
    else
       call mma_dpip_KKT_device(this, x, df0dx, fval, dfdx)
     !   call mma_KKT_cpu(this, x, df0dx, fval, dfdx)
    end if
  end subroutine mma_KKT_device

  !> Implementation of the KKT residual computation for dual interior
  ! point method (dip) subsolve of MMA algorithm.
  module subroutine mma_dip_KKT_device(this, x, df0dx, fval, dfdx)
    class(mma_t), intent(inout) :: this
    real(kind=rp), dimension(this%n), intent(in) :: x
    type(vector_t), intent(in) :: fval, df0dx
    type(matrix_t), intent(in) :: dfdx

    type(vector_t) :: relambda, remu

    call relambda%init(this%m)
    call remu%init(this%m)

    ! relambda = fval - this%a%x * this%z - this%y%x + this%mu%x
    call device_add3s2(relambda%x_d, fval%x_d, this%a%x_d, 1.0_rp, -this%z, &
         this%m)
    call device_sub2(relambda%x_d, this%y%x_d, this%m)
    call device_add2(relambda%x_d, this%mu%x_d, this%m)

    ! Compute residual for mu (eta in the paper)
    call device_col3 (remu%x_d, this%lambda%x_d, this%mu%x_d, this%m)


    this%residumax = maxval([device_maxval(relambda%x_d, this%m), &
         device_maxval(remu%x_d, this%m)])
    this%residunorm = sqrt(device_norm(relambda%x_d, this%m)+ &
         device_norm(remu%x_d, this%m))

    call relambda%free()
    call remu%free()
  end subroutine mma_dip_KKT_device

  !> Implementation of the KKT residual computation for dual primal interior
  ! point method (dpip) subsolve of MMA algorithm.
  module subroutine mma_dpip_KKT_device(this, x, df0dx, fval, dfdx)
    class(mma_t), intent(inout) :: this
    real(kind=rp), dimension(this%n), intent(in) :: x
    type(vector_t), intent(in) :: fval, df0dx
    type(matrix_t), intent(in) :: dfdx

    type(vector_t) :: designx
    real(kind=rp) :: rez, rezeta
    type(vector_t) :: rey, relambda, remu, res
    type(vector_t) :: rex, rexsi, reeta
    integer :: ierr
    real(kind=rp) :: re_sq_norm

    ! create a vector type x to have a c_ptr to point to the array designx
    call designx%init(this%n)
    designx%x = x
    call device_memcpy(designx%x, designx%x_d, this%n, HOST_TO_DEVICE, &
         sync = .false.)

    call rey%init(this%m)
    call relambda%init(this%m)
    call remu%init(this%m)
    call res%init(this%m)

    call rex%init(this%n)
    call rexsi%init(this%n)
    call reeta%init(this%n)

    call device_kkt_rex(rex%x_d, df0dx%x_d, dfdx%x_d, this%xsi%x_d, &
         this%eta%x_d, this%lambda%x_d, this%n, this%m)

    call device_col3(rey%x_d, this%d%x_d, this%y%x_d, this%m)
    call device_add2(rey%x_d, this%c%x_d, this%m)
    call device_sub2(rey%x_d, this%lambda%x_d, this%m)
    call device_sub2(rey%x_d, this%mu%x_d, this%m)

    rez = this%a0 - this%zeta - device_lcsc2(this%lambda%x_d, this%a%x_d, &
         this%m)

    call device_add3s2(relambda%x_d, fval%x_d, this%a%x_d, 1.0_rp, -this%z, &
         this%m)
    call device_sub2(relambda%x_d, this%y%x_d, this%m)
    call device_add2(relambda%x_d, this%s%x_d, this%m)

    call device_sub3(rexsi%x_d, designx%x_d, this%xmin%x_d, this%n)
    call device_col2(rexsi%x_d, this%xsi%x_d, this%n)

    call device_sub3(reeta%x_d, this%xmax%x_d, designx%x_d, this%n)
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

    this%residunorm = sqrt( ( &
         device_norm(rey%x_d, this%m) + &
         rez**2 + &
         device_norm(relambda%x_d, this%m) + &
         device_norm(remu%x_d, this%m) + &
         rezeta**2 + &
         device_norm(res%x_d, this%m) &
         ) + re_sq_norm)

    write (*, '(A)') "==============================================================="
    write (*, '(A)') "MMA KKT MAX residuals"
    write (*, '(A10,ES16.7E2)') "rex", device_maxval(rex%x_d, this%n)
    write (*, '(A10,ES16.7E2)') "rey", device_maxval(rey%x_d, this%m)
    write (*, '(A10,ES16.7E2)') "rez", abs(rez)
    write (*, '(A10,ES16.7E2)') "relambda", device_maxval(relambda%x_d, this%m)
    write (*, '(A10,ES16.7E2)') "rexsi", device_maxval(rexsi%x_d, this%n)
    write (*, '(A10,ES16.7E2)') "reeta", device_maxval(reeta%x_d, this%n)
    write (*, '(A10,ES16.7E2)') "remu", device_maxval(remu%x_d, this%m)
    write (*, '(A10,ES16.7E2)') "rezeta", abs(rezeta)
    write (*, '(A10,ES16.7E2)') "res", device_maxval(res%x_d, this%m)
    write (*, '(A)') "==============================================================="

    call designx%free()
    call rey%free()
    call relambda%free()
    call remu%free()
    call res%free()
    call rex%free()
    call rexsi%free()
    call reeta%free()
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
    type(vector_t), intent(in) :: x
    type(vector_t), intent(in) :: df0dx
    type(vector_t), intent(in) :: fval
    type(matrix_t), intent(in) :: dfdx

    integer, intent(in) :: iter
    integer :: ierr

    type(vector_t):: x_diff

    call x_diff%init(this%n)
    call device_sub3 (x_diff%x_d, this%xmax%x_d, this%xmin%x_d, this%n)
     call device_memcpy(x_diff%x, x_diff%x_d, this%n, &
        DEVICE_TO_HOST, sync = .true.)
     ! print *, "iter=", iter, "device x_diff=", sum(x_diff%x)
    ! ------------------------------------------------------------------------ !
    ! Setup the current asymptotes

    if (iter .lt. 3) then
       call device_copy(this%low%x_d, x%x_d, this%n)
       call device_add2s2(this%low%x_d, x_diff%x_d, - this%asyinit, this%n)
       call device_copy(this%upp%x_d, x%x_d, this%n)
       call device_add2s2(this%upp%x_d, x_diff%x_d, this%asyinit, this%n)
        
   
    else
       call device_mma_gensub2(this%low%x_d, this%upp%x_d, x%x_d, &
            this%xold1%x_d, this%xold2%x_d, x_diff%x_d, &
            this%asydecr, this%asyincr, this%n)


     call device_memcpy(this%low%x, this%low%x_d, this%n, &
        DEVICE_TO_HOST, sync = .true.)
     call device_memcpy(this%upp%x, this%upp%x_d, this%n, &
        DEVICE_TO_HOST, sync = .true.)    
     !     print *, "this%upp%x and this%low%x =", sum(this%upp%x), sum(this%low%x)
       
    end if

    ! ------------------------------------------------------------------------ !
    ! Calculate p0j, q0j, pij, qij, alpha, and beta

    call device_mma_gensub3(x%x_d, df0dx%x_d, dfdx%x_d, this%low%x_d, &
         this%upp%x_d, this%xmin%x_d, this%xmax%x_d, this%alpha%x_d, &
         this%beta%x_d, this%p0j%x_d, this%q0j%x_d, this%pij%x_d, &
         this%qij%x_d, this%n, this%m)

    ! ------------------------------------------------------------------------ !
    ! Computing bi as defined in page 5


    call device_mma_gensub4(x%x_d, this%low%x_d, this%upp%x_d, this%pij%x_d, &
         this%qij%x_d, this%n, this%m, this%bi%x_d)

    call device_memcpy(this%bi%x, this%bi%x_d, this%m, DEVICE_TO_HOST, &
         sync = .true.)
    call MPI_Allreduce(MPI_IN_PLACE, this%bi%x, this%m, &
         mpi_real_precision, mpi_sum, neko_comm, ierr)
    call device_memcpy(this%bi%x, this%bi%x_d, this%m, HOST_TO_DEVICE, &
         sync = .true.)
    call device_sub2(this%bi%x_d, fval%x_d, this%m)

  end subroutine mma_gensub_device

  !> solve the subproblem defined by this%pij, this%qij, etc. using dual-primal
  !! interior point method
  subroutine mma_subsolve_dpip_device(this, designx)
    class(mma_t), intent(inout) :: this
    type(vector_t), intent(in) :: designx
    integer :: iter, itto, ierr
    real(kind=rp) :: epsi, residual_max, residual_norm, z, zeta, rez, rezeta, &
         delz, dz, dzeta, steg, zold, zetaold, new_residual
    ! vectors with size m
    type(vector_t) :: y, lambda, s, mu, rey, relambda, remu, res, &
         dely, dellambda, dy, dlambda, ds, dmu, yold, lambdaold, sold, muold

    ! vectors with size n
    type(vector_t) :: x, xsi, eta, rex, rexsi, reeta, &
         delx, diagx, dx, dxsi, deta, xold, xsiold, etaold

    type(vector_t) :: bb
    type(matrix_t) :: GG
    type(matrix_t) :: AA

    real(kind=rp), dimension(this%m*this%m) :: AA_buffer
    integer :: info
    integer, dimension(this%m+1) :: ipiv
    real(kind=rp) :: re_sq_norm

    integer :: i

    real(kind=rp) :: minimal_epsilon

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
    call bb%init(this%m+1)

    call GG%init(this%m, this%n)
    call AA%init(this%m+1, this%m+1)

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

          !--------------------------------------------------------------------!
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

          AA_buffer = reshape(AA%x(1:this%m, 1:this%m), [this%m * this%m])
          call MPI_Allreduce(MPI_IN_PLACE, AA_buffer, &
               this%m*this%m, mpi_real_precision, mpi_sum, neko_comm, ierr)

          AA%x(1:this%m, 1:this%m) = reshape(AA_buffer, [this%m, this%m])

          call device_memcpy(lambda%x, lambda%x_d, this%m, DEVICE_TO_HOST, &
               sync = .false.)
          call device_memcpy(mu%x, mu%x_d, this%m, DEVICE_TO_HOST, &
               sync = .false.)
          call device_memcpy(y%x, y%x_d, this%m, DEVICE_TO_HOST, &
               sync = .false.)
          call device_memcpy(s%x, s%x_d, this%m, DEVICE_TO_HOST, &
               sync = .true.)
          do i = 1, this%m
             ! update the diag AA
             AA%x(i, i) = AA%x(i, i) &
                  + s%x(i) / lambda%x(i) &
                  + 1.0_rp / (this%d%x(i) + mu%x(i) / y%x(i))
          end do
          AA%x(1:this%m, this%m+1) = this%a%x
          AA%x(this%m+1, 1:this%m) = this%a%x
          AA%x(this%m+1, this%m+1) = - zeta/z

          call device_memcpy(AA%x, AA%x_d, &
               (this%m + 1) * (this%m + 1), HOST_TO_DEVICE, sync = .true.)

          call device_memcpy(bb%x, bb%x_d, this%m+1, DEVICE_TO_HOST, &
               sync = .true.)
          call DGESV(this%m + 1, 1, AA%x, this%m + 1, ipiv, bb%x, this%m + 1, &
               info)

          if (info .ne. 0) then
             call neko_error("DGESV failed to solve the linear system in " // &
                  "mma_subsolve_dpip (device).")
          end if

          call device_memcpy(bb%x, bb%x_d, this%m+1, HOST_TO_DEVICE, &
               sync = .true.)

          dlambda%x = bb%x(1:this%m)
          call device_memcpy(dlambda%x, dlambda%x_d, this%m, HOST_TO_DEVICE, &
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
    call device_copy(this%xold1%x_d, designx%x_d, this%n)
    call device_copy(designx%x_d, x%x_d, this%n)

    ! update the parameters of the MMA object nesessary to compute KKT residual
    call device_copy(this%y%x_d, y%x_d, this%m)
    this%z = z
    call device_copy(this%lambda%x_d, lambda%x_d, this%m)
    this%zeta = zeta
    call device_copy(this%xsi%x_d, xsi%x_d, this%n)
    call device_copy(this%eta%x_d, eta%x_d, this%n)
    call device_copy(this%mu%x_d, mu%x_d, this%m)
    call device_copy(this%s%x_d, s%x_d, this%m)

     !     call device_memcpy(this%y%x, this%y%x_d, this%m, &
     !          DEVICE_TO_HOST, sync = .true.)
     !     call device_memcpy(this%lambda%x, this%lambda%x_d, this%m, &
     !          DEVICE_TO_HOST, sync = .true.)
     !     call device_memcpy(this%xsi%x, this%xsi%x_d, this%m, &
     !          DEVICE_TO_HOST, sync = .true.)
     !     call device_memcpy(this%eta%x, this%eta%x_d, this%m, &
     !          DEVICE_TO_HOST, sync = .true.)
     !     call device_memcpy(this%mu%x, this%mu%x_d, this%m, &
     !          DEVICE_TO_HOST, sync = .true.)
     !     call device_memcpy(this%s%x, this%s%x_d, this%m, &
     !          DEVICE_TO_HOST, sync = .true.)

    !free all the initiated variables in this subroutine
    call y%free()
    call lambda%free()
    call s%free()
    call mu%free()
    call rey%free()
    call relambda%free()
    call remu%free()
    call res%free()
    call dely%free()
    call dellambda%free()
    call dy%free()
    call dlambda%free()
    call ds%free()
    call dmu%free()
    call yold%free()
    call lambdaold%free()
    call sold%free()
    call muold%free()
    call x%free()
    call xsi%free()
    call eta%free()
    call rex%free()
    call rexsi%free()
    call reeta%free()
    call delx%free()
    call diagx%free()
    call dx%free()
    call dxsi%free()
    call deta%free()
    call xold%free()
    call xsiold%free()
    call etaold%free()
    call bb%free()

  end subroutine mma_subsolve_dpip_device

  !> solve the subproblem defined by this%pij, this%qij, etc. using dual
  !! interior point method
  subroutine mma_subsolve_dip_device(this, designx)
    class(mma_t), intent(inout) :: this
    type(vector_t), intent(in) :: designx
    integer :: iter, ierr
    real(kind=rp) :: epsi, residumax, z, steg
    ! vectors with size m
    type(vector_t) :: y, lambda, mu, relambda, remu, dlambda, dmu, &
         gradlambda, zerom, dd, dummy_m
    ! vectors with size n
    type(vector_t) :: x, pjlambda, qjlambda

    ! inverse of a diag matrix:
    type(vector_t) :: Ljjxinv ! [∇_x^2 Ljj]−1
    type(matrix_t) :: hijx ! ∇_x hij
    type(matrix_t) :: Hess
    real(kind=rp) :: Hesstrace

    integer :: info
    integer, dimension(this%m+1) :: ipiv
    integer :: i

    real(kind=rp) :: minimal_epsilon

    call y%init(this%m)
    call lambda%init(this%m)
    call mu%init(this%m)
    call relambda%init(this%m)
    call remu%init(this%m)
    call dlambda%init(this%m)
    call dmu%init(this%m)
    call gradlambda%init(this%m)
    call zerom%init(this%m)
    call dd%init(this%m)
    call dummy_m%init(this%m)

    call x%init(this%n)
    call pjlambda%init(this%n)
    call qjlambda%init(this%n)

    call Ljjxinv%init(this%n)
    call hijx%init(this%m,this%n)
    call Hess%init(this%m,this%m)

    call device_cfill(zerom%x_d, 0.0_rp, this%m)

    ! ------------------------------------------------------------------------ !
    ! initial value for the parameters in the subsolve based on
    ! page 15 of "https://people.kth.se/~krille/mmagcmma.pdf"

    epsi = 1.0_rp !100
    call device_cfill(y%x_d, 1.0_rp, this%m)
    ! initialize lambda with an array of ones (change to this%c%x/2 if needed!)
    call device_cfill(lambda%x_d, 1.0_rp, this%m)
    call device_cmult2(dummy_m%x_d, this%c%x_d, 0.5_rp, this%m)
    call device_pwmax(lambda%x_d, dummy_m%x_d, this%m)

    call device_cfill(mu%x_d, 1.0_rp, this%m)
    z = 0.0_rp

    ! dd is defined as this%d + 1.0e-8_rp, to avoid devision by 0 in computing y
    call device_cadd2(dd%x_d, this%d%x_d, 1.0e-8_rp, this%m)

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
         ! minimize (sum_{i=1}^{m} [ (c_i - λ_i) * y_i + 0.5 * d_i * y_i^2 ])
         ! dL_y/dy =0   => y= (λ_i - c_i)/d_i, ensure y>=0
         call device_sub3(y%x_d, lambda%x_d, c%x_d, this%m)
         ! division by dd to avoid devision by 0 (in case this%d%x_d)
         call device_invcol2(y%x_d, dd%x_d, this%m)
         call device_pwmax(y%x_d, zerom%x_d, this%m)

         ! Comput the value of z that minimizes L_z for the current λ
         ! minimize ((a_0 - sum_{i=1}^{m} λ_i * a_i) * z)
         ! if (a_0-dot_product(lambda, a)>=0) z=0 else z= 1.0
         ! ensure z>=0
         call device_col3(dummy_m%x_d, lambda%x_d, a%x_d, this%m)
         z = device_glsum(dummy_m%x_d, this%m)
         z = merge(0.0_rp, 1.0_rp, a0 - z >= 0.0)

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

         ! Download the re(lambda, mu) to CPU to calculate residumax

         call device_memcpy(relambda%x, relambda%x_d, this%m, DEVICE_TO_HOST, &
              sync = .true.)
         call device_memcpy(remu%x, remu%x_d, this%m, DEVICE_TO_HOST, &
              sync = .true.)
         residumax = maxval(abs([relambda%x, remu%x]))

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

            call device_memcpy(hijx%x, hijx%x_d, this%n*this%m, DEVICE_TO_HOST, &
                 sync = .true.)

            call device_cfill(Hess%x_d, 0.0_rp, (this%m) * (this%m) )
            call device_Hess(Hess%x_d, hijx%x_d, Ljjxinv%x_d, this%n, this%m)

            ! download Hess to CPU, mpi reduce, upload to the device
            call device_memcpy(Hess%x, Hess%x_d, this%m*this%m, DEVICE_TO_HOST, &
                 sync = .true.)
            call MPI_Allreduce(MPI_IN_PLACE, Hess%x, &
                 this%m*this%m, mpi_real_precision, mpi_sum, neko_comm, ierr)
            ! No need to upload to device since we solve LSE on CPU
            ! call device_memcpy(Hess%x, Hess%x_d, this%m*this%m, &
            ! HOST_TO_DEVICE, sync = .true.)

            !---------------contributions of z terms to Hess-------------------!
            ! There is no contibution to the Hess from z terms as z terms are
            ! linear w.r.t λ


            !---------------contributions of y terms to Hess-------------------!
            ! Only for inactive constraint, we consider contributions to Hess.
            ! Note that if d(i) = 0, the y terms (just like z terms) will not
            ! contribute to the Hessian matrix.
            ! Note that since we use DGESV to solve LSE on CPU, we dont need
            ! cuda kernel for this part

            call device_memcpy(lambda%x, lambda%x_d, this%m, DEVICE_TO_HOST, &
                 sync = .true.)
            call device_memcpy(mu%x, mu%x_d, this%m, DEVICE_TO_HOST, &
                 sync = .true.)
            call device_memcpy(y%x, y%x_d, this%m, DEVICE_TO_HOST, &
                 sync = .true.)
            do i = 1, this%m
               if (y%x(i) .gt. 0.0_rp) then
                  if (abs(this%d%x(i)) < 1.0e-15_rp) then
                     ! Hess(i, i) = Hess(i, i) - 1.0_rp/1.0e-8_rp
                  else
                     Hess%x(i, i) = Hess%x(i, i) - 1.0_rp/this%d%x(i)
                  end if
               end if
               ! Based on eq(10), note the term (-\Omega \Lambda)
               Hess%x(i, i) = Hess%x(i, i) - mu%x(i) / lambda%x(i)
            end do

            ! Improve the robustness by stablizing the Hess using
            ! Levenberg-Marquardt algorithm (heuristically)
            Hesstrace = 0.0_rp
            do i=1, this%m
               Hesstrace = Hesstrace + Hess%x(i, i)
            end do
            do i=1, this%m
               Hess%x(i,i) = Hess%x(i, i) - &
                    max(-1.0e-4_rp*Hesstrace/this%m, 1.0e-7_rp)
            end do

            call device_memcpy(gradlambda%x, gradlambda%x_d, this%m, DEVICE_TO_HOST, &
                 sync = .true.)
            call DGESV(this%m , 1, Hess%x, this%m , ipiv, &
                 gradlambda%x, this%m, info)

            if (info .ne. 0) then
               call neko_error("DGESV failed to solve the linear system in " // &
                    "mma_subsolve_dip (device).")
            end if
            call device_memcpy(gradlambda%x, gradlambda%x_d, this%m, HOST_TO_DEVICE, &
                 sync = .true.)

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

            call device_memcpy(lambda%x, lambda%x_d, this%m, DEVICE_TO_HOST, &
                 sync = .true.)
            call device_memcpy(mu%x, mu%x_d, this%m, DEVICE_TO_HOST, &
                 sync = .true.)

            ! minimize(L_x, L_y, L_z) and compute x(λ), y(λ), z(λ) for
            ! the updated values of λ

            ! Comput the value of y that minimizes L_y for the current λ
            ! minimize (sum_{i=1}^{m} [ (c_i - λ_i) * y_i + 0.5 * d_i * y_i^2 ])
            ! dL_y/dy =0   => y= (λ_i - c_i)/d_i, ensure y>=0

            call device_sub3(y%x_d, lambda%x_d, c%x_d, this%m)
            ! division by dd to avoid devision by 0 (in case this%d%x_d)
            call device_invcol2(y%x_d, dd%x_d, this%m)
            call device_pwmax(y%x_d, zerom%x_d, this%m)

            ! Comput the value of z that minimizes L_z for the current λ
            ! minimize ((a_0 - sum_{i=1}^{m} λ_i * a_i) * z)
            ! if (a_0-dot_product(lambda, a)>=0) z=0 else z= 1.0
            ! ensure z>=0
            call device_col3(dummy_m%x_d, lambda%x_d, a%x_d, this%m)
            z = device_glsum(dummy_m%x_d, this%m)
            z = merge(0.0_rp, 1.0_rp, a0 - z >= 0.0)

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


            !> Download the re(lambda, mu) to CPU to calculate residumax

            call device_memcpy(relambda%x, relambda%x_d, this%m, DEVICE_TO_HOST, &
                 sync = .true.)
            call device_memcpy(remu%x, remu%x_d, this%m, DEVICE_TO_HOST, &
                 sync = .true.)
            residumax = maxval(abs([relambda%x, remu%x]))
         end do
       end associate
       epsi = 0.1_rp * epsi
    end do outer

    ! Save the new designx
    call device_copy(this%xold2%x_d, this%xold1%x_d, this%n)
    call device_copy(this%xold1%x_d, designx%x_d, this%n)
    call device_copy(designx%x_d, x%x_d, this%n)

    ! update the parameters of the MMA object nesessary to compute KKT residual
    call device_copy(this%y%x_d, y%x_d, this%m)
    this%z = z
    call device_copy(this%lambda%x_d, lambda%x_d, this%m)
    call device_copy(this%mu%x_d, mu%x_d, this%m)

    call y%free()
    call lambda%free()
    call mu%free()
    call relambda%free()
    call remu%free()
    call dlambda%free()
    call dmu%free()
    call gradlambda%free()
    call zerom%free()
    call dd%free()
    call dummy_m%free()

    call x%free()
    call pjlambda%free()
    call qjlambda%free()

    call Ljjxinv%free()
    call hijx%free()
    call Hess%free()
  end subroutine mma_subsolve_dip_device





  subroutine test_subsolve_dpip_cpu(this, designx)
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

             ! Compute squared norms for the residuals
             re_sq_norm = norm2(rex)**2 + norm2(rexsi)**2 + norm2(reeta)**2
             call MPI_Allreduce(MPI_IN_PLACE, re_sq_norm, &
                  1, mpi_real_precision, mpi_sum, neko_comm, ierr)

             residual_small = [rey, rez, relambda, remu, rezeta, res]
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

  end subroutine test_subsolve_dpip_cpu


  subroutine test_gensub_cpu(this, iter, xdesign, df0dx, fval, dfdx)
    ! ----------------------------------------------------- !
    ! Generate the approximation sub problem by computing   !
    ! the lower and upper asymtotes and the other necessary !
    ! parameters (alpha, beta, p0j, q0j, pij, qij, ...).    !
    ! ----------------------------------------------------- !
    class(mma_t), intent(inout) :: this
    real(kind=rp), dimension(this%n), intent(in) :: xdesign
    type(vector_t) :: df0dx, fval
    type(matrix_t) :: dfdx
    integer, intent(in) :: iter
    integer :: i, j, ierr
    real(kind=rp), dimension(this%n) :: x_diff
    real(kind=rp) :: asy_factor

    x_diff = this%xmax%x - this%xmin%x
     ! print *, "iter=", iter,  "cpu x_diff=", sum(x_diff)

    ! ------------------------------------------------------------------------ !
    ! Setup the current asymptotes
    associate(low => this%low%x, upp => this%upp%x, &
         x_1 => this%xold1%x, x_2 => this%xold2%x, x => xdesign)

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
     !     print *, "upp and low =", sum(upp), sum(low)

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
         low => this%low%x, upp => this%upp%x, x => xdesign)

      alpha = max(xmin, low + 0.1_rp*(x - low), x - 0.5_rp*x_diff)
      beta = min(xmax, upp - 0.1_rp*(upp - x), x + 0.5_rp*x_diff)
    end associate

    ! ------------------------------------------------------------------------ !
    ! Calculate p0j, q0j, pij, qij
    ! where j = 1,2,...,n and i = 1,2,...,m  (eq(2.3)-eq(2.5))

    associate(p0j => this%p0j%x, q0j => this%q0j%x, &
         pij => this%pij%x, qij => this%qij%x, &
         low => this%low%x, upp => this%upp%x, x => xdesign, &
         dfdx => dfdx%x, df0dx => df0dx%x)

      p0j = (&
           1.001_rp * max(df0dx, 0.0_rp) &
           + 0.001_rp * max(-df0dx, 0.0_rp) &
           + 0.00001_rp / max(x_diff, 0.00001_rp) &
           ) * (upp - x)**2

      q0j = (&
           0.001_rp * max(df0dx, 0.0_rp) &
           + 1.001_rp * max(-df0dx, 0.0_rp) &
           + 0.00001_rp / max(x_diff, 0.00001_rp)&
           ) * (x - low)**2

      do j = 1, this%n
         do i = 1, this%m
            pij(i, j) = (&
                 1.001_rp * max(dfdx(i, j), 0.0_rp) &
                 + 0.001_rp * max(-dfdx(i, j), 0.0_rp) &
                 + 0.00001_rp / max(x_diff(j), 0.00001_rp) &
                 ) * (upp(j) - x(j))**2

            qij(i, j) = (&
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
         low => this%low%x, upp => this%upp%x, x => xdesign)

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
      bi = bi - fval%x

    end associate
  end subroutine test_gensub_cpu
end submodule mma_device
