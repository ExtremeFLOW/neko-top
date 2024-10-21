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

    print *, "I am in mma_gpu.f90"
    ! call globaltmp_m%init(this%m)
    ! if (iter .lt. 3) then
        ! call device_add3s2(this%low%x_d,this%xmax%x_d,this%xmin%x_d,-this%asyinit,this%asyinit,this%n)
        ! call device_add2(this%low%x_d,x%x_d,this%n)

        ! call device_add3s2( this%upp%x_d,this%xmax%x_d,this%xmin%x_d,this%asyinit,- this%asyinit,this%n)
        ! call device_add2(this%upp%x_d,x%x_d,this%n)
        ! call mma_gensub1_gpu(this%low%x_d, this%upp%x_d,x%x_d, this%xmin%x_d, this%xmax%x_d, this%asyinit, this%n)
    ! else
      ! call mma_gensub2_gpu(this%low%x_d, this%upp%x_d, x%x_d, this%xold1%x_d, this%xold2%x_d,this%xmin%x_d, this%xmax%x_d, &
        ! this%asydecr, this%asyincr, this%n)
  ! end if
  ! call mma_gensub3_gpu(x%x_d, df0dx%x_d, dfdx%x_d,this%low%x_d, this%upp%x_d, this%xmin%x_d, this%xmax%x_d,this%alpha%x_d, &
    ! this%beta%x_d, this%p0j%x_d, this%q0j%x_d, this%pij%x_d, this%qij%x_d, this%n, this%m) 
  ! call mma_gensub4_gpu(x%x_d, this%low%x_d, this%upp%x_d, this%pij%x_d, this%qij%x_d, this%n, this%m, this%bi%x_d)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!cpu gpu transfer part
  ! globaltmp_m%x=0.0_rp
  ! call device_memcpy(this%bi%x, this%bi%x_d, this%m, DEVICE_TO_HOST, sync=.false.)
  ! call MPI_Allreduce(this%bi%x, globaltmp_m%x, this%m, &
  !   mpi_real_precision, mpi_sum, neko_comm, ierr)
  ! call device_memcpy(fval%x, fval%x_d, this%m, DEVICE_TO_HOST, sync=.false.)
  ! this%bi%x = globaltmp_m%x - fval%x
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! call cuda_mpisum(this%bi%x_d, this%m)
  ! call device_sub(this%bi%x_d, fval%x_d, this%m)
end subroutine mma_gensub_gpu


subroutine mma_subsolve_dpip_gpu(this, designx)
    class(mma_t), intent(inout) :: this
    type(vector_t), intent(in) :: designx
end subroutine mma_subsolve_dpip_gpu

subroutine mma_KKT_gpu(this, x, df0dx, fval, dfdx)
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
end subroutine mma_KKT_gpu


end submodule mma_gpu
