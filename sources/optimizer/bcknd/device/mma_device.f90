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

  use device_math
  use device_mma_math
  use neko_config, only: NEKO_BCKND_DEVICE
  use device_math_ext
  use device, only : device_memcpy, HOST_TO_DEVICE, DEVICE_TO_HOST
  use comm, only : pe_rank

  ! for printing the result of kernels and debugging
  use, intrinsic :: iso_c_binding

  implicit none

contains


  module subroutine mma_init_attributes_device(this, x, n, m, a0, a, c, d, &
       xmin, xmax, max_iter, epsimin, asyinit, asyincr, asydecr)
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
    real(kind=rp), intent(in), dimension(n) :: xmax, xmin
    real(kind=rp), intent(in), dimension(m) :: a, c, d
    real(kind=rp), intent(in) :: a0
    integer, intent(in), optional :: max_iter
    real(kind=rp), intent(in), optional :: epsimin, asyinit, asyincr, asydecr

    call this%free()

    this%n = n
    this%m = m

    call this%xold1%init(n)
    call this%xold2%init(n)
    this%xold1%x = x
    this%xold2%x = x

    call this%alpha%init(n)
    call this%beta%init(n)

    call this%a%init(m)
    call this%c%init(m)
    call this%d%init(m)
    call this%low%init(n)
    call this%upp%init(n)
    call this%xmax%init(n)
    call this%xmin%init(n)

    !internal dummy variables for MMA
    call this%p0j%init(n)
    call this%q0j%init(n)
    call this%pij%init(m, n)
    call this%qij%init(m, n)
    call this%bi%init(m)

    !---nesessary for KKT check after updating df0dx, fval, dfdx --------
    call this%y%init(m)
    call this%lambda%init(m)
    call this%s%init(m)
    call this%mu%init(m)
    call this%xsi%init(n)
    call this%eta%init(n)

    this%a0 = a0
    this%a%x = a
    this%c%x = c
    this%d%x = d

    !setting the bounds for the design variable based on the problem
    this%xmax%x = xmax
    this%xmin%x = xmin

    this%low%x(:) = minval(x)
    this%upp%x(:) = maxval(x)

    !setting KKT norms to a large number for the initial design
    this%residumax = huge(0.0_rp)
    this%residunorm = huge(0.0_rp)

    ! upload all init values to device pointers
    call device_memcpy(this%xold1%x, this%xold1%x_d, this%n, HOST_TO_DEVICE, &
         sync=.false.)
    call device_memcpy(this%xold1%x, this%xold2%x_d, this%n, HOST_TO_DEVICE, &
         sync=.false.)
    call device_memcpy(this%a%x, this%a%x_d, this%m, HOST_TO_DEVICE, &
         sync=.false.)
    call device_memcpy(this%c%x, this%c%x_d, this%m, HOST_TO_DEVICE, &
         sync=.false.)
    call device_memcpy(this%d%x, this%d%x_d, this%m, HOST_TO_DEVICE, &
         sync=.false.)
    call device_memcpy(this%xmax%x, this%xmax%x_d, this%n, HOST_TO_DEVICE, &
         sync=.false.)
    call device_memcpy(this%xmin%x, this%xmin%x_d, this%n, HOST_TO_DEVICE, &
         sync=.false.)
    ! ------------------------------------------------------------------------ !
    ! Assign defaults if nothing is parsed

    ! based on the Cpp Code by Niels
    if (.not. present(epsimin)) this%epsimin = 1.0e-9_rp * sqrt(real(m + n, rp))
    if (.not. present(max_iter)) this%max_iter = 100

    ! Following parameters are set based on eq.3.8:--------
    if (.not. present(asyinit)) this%asyinit = 0.5_rp
    if (.not. present(asyincr)) this%asyincr = 1.2_rp
    if (.not. present(asydecr)) this%asydecr = 0.7_rp

    ! Assign values from inputs when present
    if (present(max_iter)) this%max_iter = max_iter
    if (present(epsimin)) this%epsimin = epsimin
    if (present(asyinit)) this%asyinit = asyinit
    if (present(asyincr)) this%asyincr = asyincr
    if (present(asydecr)) this%asydecr = asydecr
!     if (present(backnd)) this%backnd = backnd
    if (pe_rank .eq. 0) then
       print *, "MMA is initialized with a0=", a0, ", a=", a, ", c=", c, &
            ", d=", d, "epsimin =", this%epsimin
    end if
    !the object is correctly initialized
    this%is_initialized = .true.
  end subroutine mma_init_attributes_device

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
       write(stderr, *) "The MMA object is not initialized."
       error stop
    end if
    

    call xdesign%init(this%n)
    call device_memcpy(x, xdesign%x_d, this%n, HOST_TO_DEVICE, sync=.false.)

    ! generate a convex approximation of the problem
    call mma_gensub_device(this, iter, xdesign, df0dx, fval, dfdx)
    !solve the approximation problem using interior point method
    call mma_subsolve_dpip_device(this, xdesign)
    !update the design vector x on the host
    call device_memcpy(x, xdesign%x_d, this%n, DEVICE_TO_HOST, sync=.false.)

    this%is_updated = .true.
  end subroutine mma_update_device

  !> Deallocate the MMA object.
  module subroutine mma_free_device(this)

    class(mma_t), intent(inout) :: this

    ! Deallocate the internal vectors
    call this%xold1%free()
    call this%xold2%free()
    call this%alpha%free()
    call this%beta%free()
    call this%a%free()
    call this%c%free()
    call this%d%free()
    call this%low%free()
    call this%upp%free()
    call this%xmax%free()
    call this%xmin%free()
    call this%p0j%free()
    call this%q0j%free()
    call this%bi%free()
    call this%y%free()
    call this%lambda%free()
    call this%s%free()
    call this%mu%free()
    call this%xsi%free()
    call this%eta%free()

    ! Deallocate the internal dummy matrices
    call this%pij%free()
    call this%qij%free()

    this%is_initialized = .false.
    this%is_updated = .false.

  end subroutine mma_free_device

  module subroutine mma_KKT_device(this, x, df0dx, fval, dfdx)
     class(mma_t), intent(inout) :: this
     real(kind=rp), dimension(this%n), intent(in) :: x
     type(vector_t), intent(in) :: fval, df0dx
     type(matrix_t), intent(in) :: dfdx
     
     type(vector_t) :: designx
     real(kind=rp) :: rez, rezeta
     type(vector_t) :: rey, relambda, remu, res
     type(vector_t) :: rex, rexsi, reeta
     real(kind=rp) :: residu_val !!!(3*this%n+4*this%m+2)
     real(kind=rp), dimension(4*this%m+2) :: residu_small !!!(4*this%m+2)
     integer :: ierr
     real(kind=rp) :: re_xstuff_squ_global
     real(kind=rp) :: globaltemp_norm

     ! create a vector type x to have a c_ptr to point to the array designx
     call designx%init(this%n)
     designx%x = x
     call device_memcpy(designx%x, designx%x_d, this%n, HOST_TO_DEVICE, &
          sync=.false.)


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

     rezeta = this%zeta*this%z

     call device_col3(res%x_d, this%lambda%x_d, this%s%x_d, this%m)

     residu_val = maxval([device_maxval(rex%x_d, this%n), &
          device_maxval(rey%x_d, this%m), rez, &
          device_maxval(relambda%x_d, this%m), &
          device_maxval(rexsi%x_d, this%n), device_maxval(reeta%x_d, this%n), &
          device_maxval(remu%x_d, this%m), rezeta, &
          device_maxval(res%x_d, this%m)])

     call MPI_Allreduce(residu_val, this%residumax, 1, &
          mpi_real_precision, mpi_max, neko_comm, ierr)

     globaltemp_norm = device_norm(rex%x_d, this%n) + &
          device_norm(rexsi%x_d, this%n) + device_norm(reeta%x_d, this%n)
     call MPI_Allreduce(globaltemp_norm, re_xstuff_squ_global, 1, &
          mpi_real_precision, mpi_sum, neko_comm, ierr)
     this%residunorm = sqrt(device_norm(rey%x_d, this%m) + rez**2 + &
          device_norm(relambda%x_d, this%m) + device_norm(remu%x_d, this%m) + &
          rezeta**2+device_norm(res%x_d, this%m) + re_xstuff_squ_global)
  end subroutine mma_KKT_device

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
     integer :: i, j, ierr
     type(vector_t) :: globaltmp_m

     call globaltmp_m%init(this%m)
     if (iter .lt. 3) then
          call device_add3s2(this%low%x_d, this%xmax%x_d, this%xmin%x_d, &
               - this%asyinit, this%asyinit, this%n)
          call device_add2(this%low%x_d, x%x_d, this%n)

          call device_add3s2( this%upp%x_d, this%xmax%x_d, this%xmin%x_d, &
               this%asyinit, - this%asyinit, this%n)
          call device_add2(this%upp%x_d, x%x_d, this%n)
     else
          call device_mma_gensub2(this%low%x_d, this%upp%x_d, x%x_d, &
               this%xold1%x_d, this%xold2%x_d, this%xmin%x_d, this%xmax%x_d, &
               this%asydecr, this%asyincr, this%n)
     end if
     call device_memcpy(this%upp%x, this%upp%x_d, this%n, DEVICE_TO_HOST, &
          sync=.true.)
     call device_memcpy(this%low%x, this%low%x_d, this%n, DEVICE_TO_HOST, &
          sync=.true.)
     call device_mma_gensub3(x%x_d, df0dx%x_d, dfdx%x_d, this%low%x_d, &
          this%upp%x_d, this%xmin%x_d, this%xmax%x_d, this%alpha%x_d, &
          this%beta%x_d, this%p0j%x_d, this%q0j%x_d, this%pij%x_d, &
          this%qij%x_d, this%n, this%m)

     call device_memcpy(this%alpha%x, this%alpha%x_d, this%n, DEVICE_TO_HOST, &
         sync=.true.)
     call device_memcpy(this%beta%x, this%beta%x_d, this%n, DEVICE_TO_HOST, &
          sync=.true.)

     call device_mma_gensub4(x%x_d, this%low%x_d, this%upp%x_d, this%pij%x_d, &
          this%qij%x_d, this%n, this%m, this%bi%x_d)
     call device_memcpy(this%pij%x, this%pij%x_d, this%n*this%m, &
          DEVICE_TO_HOST, sync=.true.)
     call device_memcpy(this%qij%x, this%qij%x_d, this%n*this%m, &
          DEVICE_TO_HOST, sync=.true.)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!cpu gpu transfer and global sum
     globaltmp_m%x=0.0_rp
     call device_memcpy(this%bi%x, this%bi%x_d, this%m, DEVICE_TO_HOST, &
          sync=.true.)
     call MPI_Allreduce(this%bi%x, globaltmp_m%x, this%m, mpi_real_precision, &
          mpi_sum, neko_comm, ierr)
     call device_memcpy(globaltmp_m%x, globaltmp_m%x_d, this%m, &
          HOST_TO_DEVICE, sync=.true.)
     call device_sub3(this%bi%x_d, globaltmp_m%x_d, fval%x_d, this%m)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call device_memcpy(this%bi%x, this%bi%x_d, this%m, DEVICE_TO_HOST, &
          sync=.true.)

     call globaltmp_m%free()
  end subroutine mma_gensub_device

  subroutine mma_subsolve_dpip_device(this, designx)
     class(mma_t), intent(inout) :: this
     type(vector_t), intent(in) :: designx
     integer :: i, j, k, iter, ggdumiter, itto, ierr
     real(kind=rp) :: epsi, residumax, residunorm, z, zeta, rez, rezeta, &
          delz, dz, dzeta, steg, dummy_one, zold, zetaold, newresidu
     !! the vectors with size m     
     type(vector_t) :: y, lambda, s, mu, rey, relambda, remu, res, &
          dely, dellambda, dy, dlambda, ds, dmu, yold, lambdaold, sold, muold
     type(vector_t) :: globaltmp_m

     !! the vectors with size n
     type(vector_t) :: x, xsi, eta, rex, rexsi, reeta, &
          delx, diagx, dx, dxsi, deta, xold, xsiold, etaold
     real(kind=rp) :: residu, residu_small, ratio_xx_dxx

     type(vector_t) :: bb
     type(matrix_t) :: GG
     type(matrix_t) :: AA
     type(matrix_t) :: globaltmp_mm

     integer :: info
     integer, dimension(this%m+1) :: ipiv
     real(kind=rp) :: re_xstuff_squ_global

     integer :: nglobal

     real(kind=rp) :: cons


     call globaltmp_m%init(this%m)
     call globaltmp_mm%init(this%m, this%m)


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
     dummy_one = 1.0_rp
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
     call device_memcpy(xsi%x, xsi%x_d, this%n, DEVICE_TO_HOST, sync=.true.)
     call device_memcpy(eta%x, eta%x_d, this%n, DEVICE_TO_HOST, sync=.true.)
     call device_memcpy(mu%x, mu%x_d, this%m, DEVICE_TO_HOST, sync=.true.)

     outer: do while (epsi .gt. 0.9*this%epsimin)
          ! calculating residuals based on
          ! "https://people.kth.se/~krille/mmagcmma.pdf" for the variables
          ! x, y, z, lambda residuals based on eq(5.9a)-(5.9d), respectively.
          call device_rex(rex%x_d, x%x_d, this%low%x_d, this%upp%x_d, &
               this%pij%x_d, this%p0j%x_d, this%qij%x_d, this%q0j%x_d, &
               lambda%x_d, xsi%x_d, eta%x_d, this%n, this%m)
          call MPI_Allreduce(this%n, nglobal, 1, MPI_INTEGER, mpi_sum, &
               neko_comm, ierr)
          call device_col3(rey%x_d, this%d%x_d, y%x_d, this%m)
          call device_add2(rey%x_d, this%c%x_d, this%m)
          call device_sub2(rey%x_d, lambda%x_d, this%m)
          call device_sub2(rey%x_d, mu%x_d, this%m)
          rez = this%a0 - zeta - device_lcsc2(lambda%x_d, this%a%x_d, this%m)
          call device_cfill(relambda%x_d, 0.0_rp, this%m)
          call device_relambda(relambda%x_d, x%x_d, this%upp%x_d, &
               this%low%x_d, this%pij%x_d, this%qij%x_d, this%n, this%m)
          call device_memcpy(relambda%x, relambda%x_d, this%m, DEVICE_TO_HOST, &
               sync=.true.)

          globaltmp_m%x = 0.0_rp
          call MPI_Allreduce(relambda%x, globaltmp_m%x, this%m, &
               mpi_real_precision, mpi_sum, neko_comm, ierr)

          call device_memcpy(globaltmp_m%x, globaltmp_m%x_d, this%m, &
               HOST_TO_DEVICE, sync=.true.)
          call device_add3s2(relambda%x_d, globaltmp_m%x_d, this%a%x_d, &
               1.0_rp, -z, this%m)
          call device_sub2(relambda%x_d, y%x_d, this%m)
          call device_add2(relambda%x_d, s%x_d, this%m)
          call device_sub2(relambda%x_d, this%bi%x_d, this%m)

          call device_sub3(rexsi%x_d, x%x_d, this%alpha%x_d, this%n)
          call device_col2(rexsi%x_d, xsi%x_d, this%n)
          call device_cadd(rexsi%x_d, -epsi, this%n)

          call device_sub3(reeta%x_d, this%beta%x_d, x%x_d, this%n)
          call device_col2(reeta%x_d, eta%x_d, this%n)
          call device_cadd(reeta%x_d, -epsi, this%n)

          call device_col3(remu%x_d, mu%x_d, y%x_d, this%m)
          call device_cadd(remu%x_d, -epsi, this%m)

          rezeta = zeta*z -epsi

          call device_col3(res%x_d, lambda%x_d, s%x_d, this%m)
          call device_cadd(res%x_d, -epsi, this%m)

          cons = 0.0_rp
          cons = maxval([device_maxval(rex%x_d, this%n), &
               device_maxval(rey%x_d, this%m), rez, &
               device_maxval(relambda%x_d, this%m), &
               device_maxval(rexsi%x_d, this%n), &
               device_maxval(reeta%x_d, this%n), &
               device_maxval(remu%x_d, this%m), rezeta, &
               device_maxval(res%x_d, this%m)])
          residumax = 0.0_rp
          call MPI_Allreduce(cons, residumax, 1, mpi_real_precision, mpi_max, &
               neko_comm, ierr)

          re_xstuff_squ_global = 0.0_rp
          cons = device_norm(rex%x_d, this%n) + &
               device_norm(rexsi%x_d, this%n)+device_norm(reeta%x_d, this%n)
          call MPI_Allreduce(cons, re_xstuff_squ_global, 1, &
               mpi_real_precision, mpi_sum, neko_comm, ierr)
          cons = device_norm(rey%x_d, this%m) + rez**2 + &
               device_norm(relambda%x_d, this%m) + &
               device_norm(remu%x_d, this%m)+ &
               rezeta**2+device_norm(res%x_d, this%m)
          residunorm = sqrt(cons + re_xstuff_squ_global)



          do iter = 1, 100 !this%max_iter !ittt
               if (iter .gt. (this%max_iter -2)) then
                  ! print *, "The mma inner loop seems not to converge"
                  ! print *, "residumax = ", residumax, "for epsi = ", epsi, &
                  !         ", ittt  = ", iter, "out of ", this%max_iter
               end if
               if (residumax .lt. epsi) exit
               call device_delx(delx%x_d, x%x_d, this%low%x_d, this%upp%x_d, &
                    this%pij%x_d, this%qij%x_d, this%p0j%x_d, this%q0j%x_d, &
                    this%alpha%x_d, this%beta%x_d, lambda%x_d, epsi, this%n, &
                    this%m)

               call device_col3(dely%x_d, this%d%x_d, y%x_d, this%m)
               call device_add2(dely%x_d, this%c%x_d, this%m)
               call device_sub2(dely%x_d, lambda%x_d, this%m)
               call device_add2inv2(dely%x_d, y%x_d, -epsi, this%m)
               delz = this%a0 - device_lcsc2(lambda%x_d, this%a%x_d, this%m) - &
                    epsi/z
               call device_cfill(dellambda%x_d, 0.0_rp, this%m)
               call device_relambda(dellambda%x_d, x%x_d, this%upp%x_d, &
                    this%low%x_d, this%pij%x_d, this%qij%x_d, this%n, this%m)
               call device_memcpy(dellambda%x, dellambda%x_d, this%m, &
                    DEVICE_TO_HOST, sync=.true.)

               globaltmp_m%x = 0.0_rp
               call MPI_Allreduce(dellambda%x, globaltmp_m%x, this%m, &
                    mpi_real_precision, mpi_sum, neko_comm, ierr)

               call device_memcpy(globaltmp_m%x, globaltmp_m%x_d, this%m, &
                    HOST_TO_DEVICE, sync=.true.)
               call device_add3s2(dellambda%x_d, globaltmp_m%x_d, this%a%x_d, &
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

               call device_bb(bb%x_d, GG%x_d, delx%x_d, diagx%x_d, this%n, &
                    this%m)
               call device_memcpy(bb%x, bb%x_d, this%m, DEVICE_TO_HOST, &
                    sync=.true.)

               globaltmp_m%x = 0.0_rp
               call MPI_Allreduce(bb%x(1:this%m), globaltmp_m%x, this%m, &
                    mpi_real_precision, mpi_sum, neko_comm, ierr)
               call device_memcpy(globaltmp_m%x, globaltmp_m%x_d, this%m, &
                    HOST_TO_DEVICE, sync=.true.)

               call device_updatebb(bb%x_d, dellambda%x_d, dely%x_d, &
                    this%d%x_d, mu%x_d, y%x_d, delz, this%m)

               call device_cfill(AA%x_d, 0.0_rp, (this%m+1) * (this%m+1) )
               call device_AA(AA%x_d, GG%x_d, diagx%x_d, this%n, this%m)
               call device_memcpy(AA%x, AA%x_d, (this%m+1) * (this%m+1), &
                    DEVICE_TO_HOST, sync=.true.)
               globaltmp_mm%x = 0.0_rp
               call MPI_Allreduce(AA%x(1:this%m, 1:this%m), globaltmp_mm%x, &
                    this%m*this%m, mpi_real_precision, mpi_sum, neko_comm, ierr)
               call device_memcpy(globaltmp_mm%x, globaltmp_mm%x_d, &
                    (this%m) * (this%m), HOST_TO_DEVICE, sync=.true.)
               call device_updateAA(AA%x_d, globaltmp_mm%x_d, s%x_d, &
                    lambda%x_d, this%d%x_d, mu%x_d, y%x_d, this%a%x_d, &
                    zeta, z, this%m)
               call device_memcpy(AA%x, AA%x_d, (this%m+1)*(this%m+1), &
                    DEVICE_TO_HOST, sync=.true.)



               call device_memcpy(bb%x, bb%x_d, this%m+1, DEVICE_TO_HOST, &
                    sync=.true.)
               call DGESV(this%m+1, 1, AA%x, this%m+1, ipiv, bb%x, this%m+1, &
                    info)
               if (info .ne. 0) then
                  write(stderr, *) "DGESV failed in mma_device.f90."
                  write(stderr, *) "Please check mma_subsolve_dpip in mma.f90"
                  error stop
               end if
               call device_memcpy(bb%x, bb%x_d, this%m+1, HOST_TO_DEVICE, &
                    sync=.true.)

               call device_copy(dlambda%x_d, bb%x_d, this%m)
               dz = bb%x(this%m + 1)

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

               dzeta = -zeta + (epsi-zeta*dz)/z
               call device_col3(ds%x_d, dlambda%x_d, s%x_d, this%m)
               call device_cmult(ds%x_d, -1.0_rp, this%m)
               call device_cadd(ds%x_d, epsi, this%m)
               call device_invcol2(ds%x_d, lambda%x_d, this%m)
               call device_sub2(ds%x_d, s%x_d, this%m)


               steg = maxval([dummy_one, device_maxval2(dy%x_d, y%x_d, &
                    -1.01_rp, this%m), -1.01_rp*dz/z, &
                    device_maxval2(dlambda%x_d, lambda%x_d, &
                    -1.01_rp, this%m), &
                    device_maxval2(dxsi%x_d, xsi%x_d, -1.01_rp, this%n), &
                    device_maxval2(deta%x_d, eta%x_d, -1.01_rp, this%n), &
                    device_maxval2(dmu%x_d, mu%x_d, -1.01_rp, this%m), &
                    device_maxval2(ds%x_d, s%x_d, -1.01_rp, this%m), &
                    device_maxval3(dx%x_d, x%x_d, this%alpha%x_d, -1.01_rp, &
                    this%n), device_maxval3(dx%x_d, this%beta%x_d, x%x_d, &
                    1.01_rp, this%n), -1.01_rp*dzeta/zeta])
               steg = 1.0_rp/steg
               call MPI_Allreduce(steg, steg, 1, &
                    mpi_real_precision, mpi_min, neko_comm, ierr)


               call device_copy(xold%x_d, x%x_d, this%n)
               call device_copy(yold%x_d, y%x_d, this%m)
               zold = z
               call device_copy(lambdaold%x_d, lambda%x_d, this%m)
               call device_copy(xsiold%x_d, xsi%x_d, this%n)
               call device_copy(etaold%x_d, eta%x_d, this%n)
               call device_copy(muold%x_d, mu%x_d, this%m)
               zetaold = zeta
               call device_copy(sold%x_d, s%x_d, this%m)
               newresidu = 2.0*residunorm
               itto = 0
               do while ((newresidu .gt. residunorm) .and. (itto .lt. 50))
                    itto = itto + 1
                    call device_add3s2(x%x_d, xold%x_d, dx%x_d, 1.0_rp, &
                         steg, this%n)
                    call device_add3s2(y%x_d, yold%x_d, dy%x_d, 1.0_rp, &
                         steg, this%m)
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

                    !recompute the newresidu to see if this stepsize improves
                    !the residue
                    call device_rex(rex%x_d, x%x_d, this%low%x_d, &
                         this%upp%x_d, this%pij%x_d, this%p0j%x_d, &
                         this%qij%x_d, this%q0j%x_d, lambda%x_d, xsi%x_d, &
                         eta%x_d, this%n, this%m)

                    call device_memcpy(rex%x, rex%x_d, this%n, DEVICE_TO_HOST, &
                         sync=.true.)
                    call device_memcpy(xsi%x, xsi%x_d, this%n, DEVICE_TO_HOST, &
                         sync=.true.)
                    call device_memcpy(eta%x, eta%x_d, this%n, DEVICE_TO_HOST, &
                         sync=.true.)
                    call device_memcpy(lambda%x, lambda%x_d, this%m, &
                         DEVICE_TO_HOST, sync=.true.)



                    call device_col3(rey%x_d, this%d%x_d, y%x_d, this%m)
                    call device_add2(rey%x_d, this%c%x_d, this%m)
                    call device_sub2(rey%x_d, lambda%x_d, this%m)
                    call device_sub2(rey%x_d, mu%x_d, this%m)

                    rez = this%a0 - zeta - device_lcsc2(lambda%x_d, &
                         this%a%x_d, this%m)

                    call device_cfill(relambda%x_d, 0.0_rp, this%m)
                    call device_relambda(relambda%x_d, x%x_d, this%upp%x_d, &
                         this%low%x_d, this%pij%x_d, this%qij%x_d, &
                         this%n, this%m)
                    call device_memcpy(relambda%x, relambda%x_d, this%m, &
                         DEVICE_TO_HOST, sync=.true.)

                    globaltmp_m%x = 0.0_rp
                    call MPI_Allreduce(relambda%x, globaltmp_m%x, this%m, &
                         mpi_real_precision, mpi_sum, neko_comm, ierr)

                    call device_memcpy(globaltmp_m%x, globaltmp_m%x_d, &
                         this%m, HOST_TO_DEVICE, sync=.true.)



                    call device_add3s2(relambda%x_d, globaltmp_m%x_d, &
                         this%a%x_d, 1.0_rp, -z, this%m)
                    call device_sub2(relambda%x_d, y%x_d, this%m)
                    call device_add2(relambda%x_d, s%x_d, this%m)
                    call device_sub2(relambda%x_d, this%bi%x_d, this%m)


                    call device_sub3(rexsi%x_d, x%x_d, this%alpha%x_d, this%n)
                    call device_col2(rexsi%x_d, xsi%x_d, this%n)
                    call device_cadd(rexsi%x_d, -epsi, this%n)

                    call device_sub3(reeta%x_d, this%beta%x_d, x%x_d, this%n)
                    call device_col2(reeta%x_d, eta%x_d, this%n)
                    call device_cadd(reeta%x_d, -epsi, this%n)

                    call device_col3(remu%x_d, mu%x_d, y%x_d, this%m)
                    call device_cadd(remu%x_d, -epsi, this%m)

                    rezeta = zeta*z - epsi


                    call device_col3(res%x_d, lambda%x_d, s%x_d, this%m)
                    call device_cadd(res%x_d, -epsi, this%m)

                    re_xstuff_squ_global = 0.0_rp
                    cons = device_norm(rex%x_d, this%n) + &
                         device_norm(rexsi%x_d, this%n) + &
                         device_norm(reeta%x_d, this%n)
                    call MPI_Allreduce(cons, re_xstuff_squ_global, 1, &
                         mpi_real_precision, mpi_sum, neko_comm, ierr)

                    cons = device_norm(rey%x_d, this%m) + rez**2 + &
                         device_norm(relambda%x_d, this%m) + &
                         device_norm(remu%x_d, this%m) + &
                         rezeta**2+device_norm(res%x_d, this%m)
                    
                    newresidu = sqrt(cons+ re_xstuff_squ_global)

                    steg = steg/2.0_rp

                    cons = 0.0_rp
                    cons = maxval([device_maxval(rex%x_d, this%n), &
                         device_maxval(rey%x_d, this%m), rez, &
                         device_maxval(relambda%x_d, this%m), &
                         device_maxval(rexsi%x_d, this%n), &
                         device_maxval(reeta%x_d, this%n), &
                         device_maxval(remu%x_d, this%m), rezeta, &
                         device_maxval(res%x_d, this%m)])
                    !exit outer
               end do
               residunorm = newresidu
               residumax = 0.0_rp
               call MPI_Allreduce(cons, residumax, 1, mpi_real_precision, &
                    mpi_max, neko_comm, ierr)
               steg = 2.0_rp*steg
               !print *, "gpu iter = ", iter, "epsi = ", epsi, &
               !   "steg = ", steg, "residunorm = ", residunorm, &
               !     "residumax = ", residumax
          end do
          epsi = 0.1_rp * epsi
     end do outer
     call device_copy(this%xold2%x_d, this%xold1%x_d, this%n)
     call device_copy(this%xold1%x_d, designx%x_d, this%n)
     call device_copy(designx%x_d, x%x_d, this%n)
     call device_copy(this%y%x_d, y%x_d, this%m)
     this%z = z
     call device_copy(this%lambda%x_d, lambda%x_d, this%m)


     this%zeta = zeta
     call device_copy(this%xsi%x_d, xsi%x_d, this%n)
     call device_copy(this%eta%x_d, eta%x_d, this%n)
     call device_copy(this%mu%x_d, mu%x_d, this%m)
     call device_copy(this%s%x_d, s%x_d, this%m)
     !print *, "I am in mma_subsolve_dpip_gpu"
  end subroutine mma_subsolve_dpip_device



end submodule mma_device
