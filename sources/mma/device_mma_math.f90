module device_math
    use, intrinsic :: iso_c_binding, only: c_ptr, c_int
    use num_types, only: rp, c_rp
    use utils, only: neko_error
    use comm, only: NEKO_COMM, pe_size, MPI_REAL_PRECISION
    use mpi_f08, only: MPI_SUM, MPI_IN_PLACE, MPI_Allreduce

    ! ========================================================================== !
    ! Device math interfaces

    use hip_math
    use cuda_math
    use opencl_math

    implicit none
    private

    public :: 

    contains

    !> Copy a vector \f$ a = b \f$
    subroutine device_mpisum(a_d, n)
        type(c_ptr) :: a_d
        integer :: n
#if HAVE_HIP

#elif HAVE_CUDA
        call cuda_mpisum(a_d, n) 
#elif HAVE_OPENCL

#else
        call neko_error('no device backend configured')
#endif
    end subroutine device_mpisum

end module device_math
