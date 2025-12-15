program test_mma_hdf5
  use, intrinsic :: iso_c_binding
  use mpi_f08
  use neko
  use hdf5
  use num_types, only: rp
  use mma
  use vector, only: vector_t
  implicit none

  integer :: ierr
  integer :: n, m
  type(mma_t) :: obj
  type(vector_t) :: x, read_x
  real(kind=rp), allocatable :: a(:), c(:), d(:), xmin(:), xmax(:)
  real(kind=rp) :: a0
  integer :: max_iter
  real(kind=rp) :: epsimin, asyinit, asyincr, asydecr
  character(len=:), allocatable :: bcknd_arg, subsolver_arg

  character(len=256) :: fname
  integer(hid_t) :: plist_id, file_id, attr_id, dset_id
  integer(hsize_t), dimension(1) :: ddim
  integer :: read_n, read_m
  logical :: error_flag = .false.

  call neko_init()

  n = 4
  m = 2

  call x%init(n)
  x%x = [1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp]
  x%x = x%x + real(pe_rank, kind=rp) * 0.01_rp

  allocate(a(m), c(m), d(m), xmin(n), xmax(n))
  a = 0.5_rp
  c = 1.0_rp
  d = 0.0_rp
  xmin = 0.0_rp
  xmax = 1.0_rp
  a0 = 1.0_rp
  max_iter = 10
  epsimin = 1.0e-9_rp
  asyinit = 0.5_rp
  asyincr = 1.2_rp
  asydecr = 0.7_rp

  bcknd_arg = 'cpu'
  subsolver_arg = 'dip'
  call obj%init_from_components(x, n, m, a0, a, c, d, xmin, xmax, &
       max_iter, epsimin, asyinit, asyincr, asydecr, bcknd_arg, subsolver_arg)

  write(fname, '(A,I0,A)') 'mma_test_rank', pe_rank, '.h5'

#ifdef HAVE_HDF5
  call obj%write_hdf5(trim(fname))

  ! -------------------------------------------------------------------------- !
  ! Read back the data

  ! Initialize reader values
  read_n = -1
  read_m = -1

  call read_x%init(n)

  ! -------------------------------------------------------------------------- !
  ! Read from HDF5 file

  ! Open file and prepare reading
  call h5open_f(ierr)
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
  call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD%mpi_val, MPI_INFO_NULL%mpi_val, ierr)
  call h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, ierr, access_prp = plist_id)

  ! Read attributes
  ddim(1) = 1
  call h5aopen_by_name_f(file_id, '/MMA/Parameters', 'n', attr_id, ierr)
  call h5aread_f(attr_id, H5T_NATIVE_INTEGER, read_n, ddim, ierr)
  call h5aclose_f(attr_id, ierr)

  call h5aopen_by_name_f(file_id, '/MMA/Parameters', 'm', attr_id, ierr)
  call h5aread_f(attr_id, H5T_NATIVE_INTEGER, read_m, ddim, ierr)
  call h5aclose_f(attr_id, ierr)

  ! Read datasets
  ddim(1) = n
  call h5dopen_f(file_id, '/MMA/xold1', dset_id, ierr)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, read_x%x, ddim, ierr)
  call h5dclose_f(dset_id, ierr)

  ! Close file and property list
  call h5fclose_f(file_id, ierr)
  call h5pclose_f(plist_id, ierr)
  call h5close_f(ierr)

  ! -------------------------------------------------------------------------- !
  ! Check the values read from file

  if (read_n .ne. n) then
     write(*,*) 'TEST FAILED: read n=', read_n, ' expected ', n
     error_flag = .true.
  end if
  if (read_m .ne. m) then
     write(*,*) 'TEST FAILED: read m=', read_m, ' expected ', m
     error_flag = .true.
  end if

  if (any(abs(read_x%x - x%x) > 1.0e-12_rp)) then
     write(*,*) 'TEST FAILED: read x does not match written x'
     write(*,*) ' written x: ', x%x
     write(*,*) ' read x:    ', read_x%x
     error_flag = .true.
  end if

  ! Clean up
  call read_x%free()
  call obj%free()
  deallocate(a, c, d, xmin, xmax)
#else
  write(*,*) 'HDF5 not available; test skipped'
#endif

  if (error_flag) then
     write(*,*) 'TEST FAILED'
     error stop 1
  else
     write(*,*) 'TEST PASSED'
  end if
  call neko_finalize()
end program test_mma_hdf5
