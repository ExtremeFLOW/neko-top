

#ifdef HAVE_HDF5

program test_mma_hdf5
  use, intrinsic :: iso_c_binding
  use mpi_f08
  use neko
  use hdf5
  use mma
  implicit none

  integer :: ierr
  integer :: n, m
  type(mma_t) :: obj
  type(vector_t) :: x, read_xold1
  real(kind=rp), allocatable :: a(:), c(:), d(:), xmin(:), xmax(:)
  real(kind=rp) :: a0
  integer :: max_iter

  real(kind=rp) :: epsimin, asyinit, asyincr, asydecr

  character(len=:), allocatable :: bcknd, subsolver

  character(len=256) :: fname
  integer(hid_t) :: H5T_NEKO_REAL
  integer(hid_t) :: plist_id, file_id, attr_id, dset_id, atype_id, mem_type_id
  integer(hsize_t), dimension(1) :: ddim
  integer(hsize_t) :: type_size
  integer :: read_n, read_m
  integer :: read_max_iter
  real(kind=rp) :: read_asyinit, read_asyincr, read_asydecr, read_epsimin
  character(len=:), allocatable :: read_bcknd, read_subsolver

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

  bcknd = 'cpu'
  subsolver = 'dip'
  call obj%init_from_components(x, n, m, a0, a, c, d, xmin, xmax, &
       max_iter, epsimin, asyinit, asyincr, asydecr, bcknd, subsolver)

  fname = 'mma_test.h5'

  call obj%write_hdf5(trim(fname))

  ! -------------------------------------------------------------------------- !
  ! Read back the data

  ! Initialize reader values
  read_n = -1
  read_m = -1
  read_max_iter = -1
  read_asyinit = -1.0_rp
  read_asyincr = -1.0_rp
  read_asydecr = -1.0_rp
  read_epsimin = -1.0_rp

  read_bcknd = ''
  read_subsolver = ''

  call read_xold1%init(n)


  ! -------------------------------------------------------------------------- !
  ! Read from HDF5 file

  ! Open file and prepare reading
  call h5open_f(ierr)
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
  call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD%mpi_val, &
       MPI_INFO_NULL%mpi_val, ierr)
  call h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, ierr, &
       access_prp = plist_id)

  ! Assign the correct HDF5 data type based on the neko real kind
  select case (rp)
  case (dp)
     H5T_NEKO_REAL = H5T_NATIVE_DOUBLE
  case (sp)
     H5T_NEKO_REAL = H5T_NATIVE_REAL
  case default
     call neko_error('mma: unsupported real kind for HDF5')
  end select

  ! Read attributes
  ddim(1) = 1
  call h5aopen_by_name_f(file_id, '/MMA/Parameters', 'n', attr_id, ierr)
  call h5aread_f(attr_id, H5T_NATIVE_INTEGER, read_n, ddim, ierr)
  call h5aclose_f(attr_id, ierr)

  call h5aopen_by_name_f(file_id, '/MMA/Parameters', 'm', attr_id, ierr)
  call h5aread_f(attr_id, H5T_NATIVE_INTEGER, read_m, ddim, ierr)
  call h5aclose_f(attr_id, ierr)

  call h5aopen_by_name_f(file_id, '/MMA/Parameters', 'max_iter', attr_id, ierr)
  call h5aread_f(attr_id, H5T_NATIVE_INTEGER, read_max_iter, ddim, ierr)
  call h5aclose_f(attr_id, ierr)

  call h5aopen_by_name_f(file_id, '/MMA/Parameters', 'asyinit', attr_id, ierr)
  call h5aread_f(attr_id, H5T_NEKO_REAL, read_asyinit, ddim, ierr)
  call h5aclose_f(attr_id, ierr)

  call h5aopen_by_name_f(file_id, '/MMA/Parameters', 'asyincr', attr_id, ierr)
  call h5aread_f(attr_id, H5T_NEKO_REAL, read_asyincr, ddim, ierr)
  call h5aclose_f(attr_id, ierr)

  call h5aopen_by_name_f(file_id, '/MMA/Parameters', 'asydecr', attr_id, ierr)
  call h5aread_f(attr_id, H5T_NEKO_REAL, read_asydecr, ddim, ierr)
  call h5aclose_f(attr_id, ierr)

  call h5aopen_by_name_f(file_id, '/MMA/Parameters', 'epsimin', attr_id, ierr)
  call h5aread_f(attr_id, H5T_NEKO_REAL, read_epsimin, ddim, ierr)
  call h5aclose_f(attr_id, ierr)

  ! Read strings
  call h5aopen_by_name_f(file_id, '/MMA/Parameters', 'bcknd', attr_id, ierr)
  ! determine attribute string size, allocate Fortran string and read with matching memory type
  call h5aget_type_f(attr_id, atype_id, ierr)
  call h5tget_size_f(atype_id, type_size, ierr)
  call h5tclose_f(atype_id, ierr)
  if (allocated(read_bcknd)) deallocate(read_bcknd)
  allocate(character(len=int(type_size)) :: read_bcknd)
  call h5tcopy_f(H5T_C_S1, mem_type_id, ierr)
  call h5tset_size_f(mem_type_id, type_size, ierr)
  call h5aread_f(attr_id, mem_type_id, read_bcknd, ddim, ierr)
  call h5tclose_f(mem_type_id, ierr)
  call h5aclose_f(attr_id, ierr)

  call h5aopen_by_name_f(file_id, '/MMA/Parameters', 'subsolver', attr_id, ierr)
  ! determine attribute string size, allocate Fortran string and read with matching memory type
  call h5aget_type_f(attr_id, atype_id, ierr)
  call h5tget_size_f(atype_id, type_size, ierr)
  call h5tclose_f(atype_id, ierr)
  if (allocated(read_subsolver)) deallocate(read_subsolver)
  allocate(character(len=int(type_size)) :: read_subsolver)
  call h5tcopy_f(H5T_C_S1, mem_type_id, ierr)
  call h5tset_size_f(mem_type_id, type_size, ierr)
  call h5aread_f(attr_id, mem_type_id, read_subsolver, ddim, ierr)
  call h5tclose_f(mem_type_id, ierr)
  call h5aclose_f(attr_id, ierr)

  ! Read datasets
  ddim(1) = n
  call h5dopen_f(file_id, '/MMA/xold1', dset_id, ierr)
  call h5dread_f(dset_id, H5T_NEKO_REAL, read_xold1%x, ddim, ierr)
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
  if (read_max_iter .ne. max_iter) then
     write(*,*) 'TEST FAILED: read max_iter=', read_max_iter, &
          ' expected ', max_iter
     error_flag = .true.
  end if
  if (abs(read_asyinit - asyinit) .gt. 1.0e-12_rp) then
     write(*,*) 'TEST FAILED: read asyinit=', read_asyinit, &
          ' expected ', asyinit
     error_flag = .true.
  end if
  if (abs(read_asyincr - asyincr) .gt. 1.0e-12_rp) then
     write(*,*) 'TEST FAILED: read asyincr=', read_asyincr, &
          ' expected ', asyincr
     error_flag = .true.
  end if
  if (abs(read_asydecr - asydecr) .gt. 1.0e-12_rp) then
     write(*,*) 'TEST FAILED: read asydecr=', read_asydecr, &
          ' expected ', asydecr
     error_flag = .true.
  end if
  if (abs(read_epsimin - epsimin) .gt. 1.0e-12_rp) then
     write(*,*) 'TEST FAILED: read epsimin=', read_epsimin, &
          ' expected ', epsimin
     error_flag = .true.
  end if

  if (trim(read_bcknd) .ne. trim(bcknd)) then
     write(*,*) 'TEST FAILED: read bcknd=', trim(read_bcknd), &
          ' expected ', trim(bcknd)
     error_flag = .true.
  end if
  if (trim(read_subsolver) .ne. trim(subsolver)) then
     write(*,*) 'TEST FAILED: read subsolver=', trim(read_subsolver), &
          ' expected ', trim(subsolver)
     error_flag = .true.
  end if

  if (any(abs(read_xold1%x - x%x) > 1.0e-12_rp)) then
     write(*,*) 'TEST FAILED: read x does not match written x'
     write(*,*) ' written x: ', x%x
     write(*,*) ' read x:    ', read_xold1%x
     error_flag = .true.
  end if

  ! -------------------------------------------------------------------------- !
  ! Clean up

  call read_xold1%free()
  call obj%free()
  deallocate(a, c, d, xmin, xmax)

  call neko_finalize()

  if (error_flag) then
     write(*,*) 'TEST FAILED'
     error stop 1
  else
     write(*,*) 'TEST PASSED'
  end if

end program test_mma_hdf5

#else

! Dummy program when HDF5 is not available
program test_mma_hdf5
  write(*,*) 'HDF5 not available; test skipped'
  stop
end program test_mma_hdf5

#endif
