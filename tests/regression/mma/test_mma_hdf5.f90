program test_mma_hdf5
  use, intrinsic :: iso_c_binding
  use mpi_f08
  use hdf5
  use num_types, only: rp
  use mma
  use vector, only: vector_t
  implicit none

  integer :: ierr, rank
  integer :: n, m
  type(mma_t) :: obj
  type(vector_t) :: x
  real(kind=rp), allocatable :: a(:), c(:), d(:), xmin(:), xmax(:)
  real(kind=rp) :: a0
  integer :: max_iter
  real(kind=rp) :: epsimin, asyinit, asyincr, asydecr
  character(len=:), allocatable :: bcknd_arg, subsolver_arg

  character(len=256) :: fname
  integer(hid_t) :: plist_id, file_id, attr_id
  integer(hsize_t), dimension(1) :: ddim
  integer :: read_n, read_m

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

  n = 4
  m = 2

  call x%init(n)
  x%x = (/ 1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp /)

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

  write(fname, '(A,I0,A)') 'mma_test_rank', rank, '.h5'

#ifdef HAVE_HDF5
  call obj%write_hdf5(trim(fname))

  ! open file and read attributes n and m
  call h5open_f(ierr)
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
  call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD%mpi_val, MPI_INFO_NULL%mpi_val, ierr)
  call h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, ierr, access_prp = plist_id)

  ddim = 1
  call h5aopen_name_f(file_id, 'n', attr_id, ierr)
  call h5aread_f(attr_id, H5T_NATIVE_INTEGER, read_n, ddim, ierr)
  call h5aclose_f(attr_id, ierr)

  call h5aopen_name_f(file_id, 'm', attr_id, ierr)
  call h5aread_f(attr_id, H5T_NATIVE_INTEGER, read_m, ddim, ierr)
  call h5aclose_f(attr_id, ierr)

  call h5fclose_f(file_id, ierr)
  call h5pclose_f(plist_id, ierr)
  call h5close_f(ierr)

  if (read_n .ne. n) then
     write(*,*) 'TEST FAILED: read n=', read_n, ' expected ', n
  else
     write(*,*) 'n OK'
  end if
  if (read_m .ne. m) then
     write(*,*) 'TEST FAILED: read m=', read_m, ' expected ', m
  else
     write(*,*) 'm OK'
  end if
#else
  write(*,*) 'HDF5 not available; test skipped'
#endif

  call MPI_Finalize(ierr)
end program test_mma_hdf5
