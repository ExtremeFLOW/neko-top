

#ifdef HAVE_HDF5

module test_mma_util
  use num_types, only: rp
  implicit none

  logical :: test_status = .true.

  interface test
     module procedure test_integer, test_real, test_real_array, test_string
  end interface test


  public :: test, test_status

contains

  subroutine test_integer(name, target, expected)
    character(len=*), intent(in) :: name
    integer, intent(in) :: target, expected
    logical :: is_equal

    is_equal = (target == expected)
    if (.not. is_equal) then
       write(*,'(A,A,A,I0,A,I0)') 'TEST FAILED: ', trim(name), ' = ', target, &
            ' expected ', expected
    end if

    test_status = test_status .and. is_equal
  end subroutine test_integer

  subroutine test_real(name, target, expected)
    character(len=*), intent(in) :: name
    real(kind=rp), intent(in) :: target, expected
    logical :: is_equal

    is_equal = (abs(target - expected) .le. 1.0e-12_rp)
    if (.not. is_equal) then
       write(*,'(A,A,A,F12.6,A,F12.6)') 'TEST FAILED: ', trim(name), ' = ', &
            target, ' expected ', expected
    end if

    test_status = test_status .and. is_equal
  end subroutine test_real

  subroutine test_real_array(name, target, expected)
    character(len=*), intent(in) :: name
    real(kind=rp), intent(in) :: target(:), expected(:)
    logical :: is_equal
    integer :: i

    is_equal = all(abs(target - expected) .le. 1.0e-12_rp)
    if (.not. is_equal) then
       write(*,'(A,A)') 'TEST FAILED: ', trim(name)
       do i = 1, size(target)
          if (abs(target(i) - expected(i)) .gt. 1.0e-12_rp) then
             write(*,'(A,I0,A,F12.6,A,F12.6)') '  ', i, ': ', target(i), &
                  ' expected ', expected(i)
          end if
       end do
    end if

    test_status = test_status .and. is_equal
  end subroutine test_real_array

  subroutine test_string(name, target, expected)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: target, expected
    logical :: is_equal

    is_equal = (trim(target) .eq. trim(expected))
    if (.not. is_equal) then
       write(*,'(A,A,A,A,A)') 'TEST FAILED: ', trim(name), ' = ', trim(target), &
            ' expected ', trim(expected)
    end if

    test_status = test_status .and. is_equal
  end subroutine test_string

end module test_mma_util

program test_mma_hdf5
  use, intrinsic :: iso_c_binding
  use mpi_f08
  use neko
  use hdf5
  use mma
  use test_mma_util
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
  integer(hid_t) :: H5T_NEKO_REAL, str_type
  integer(hid_t) :: plist_id, file_id, attr_id, dset_id, atype_id, mem_type_id
  integer(hsize_t), dimension(1) :: ddim
  integer(hsize_t) :: type_size
  integer :: read_n, read_m
  integer :: read_max_iter
  real(kind=rp) :: read_asyinit, read_asyincr, read_asydecr, read_epsimin
  real(kind=rp) :: read_a0
  real(kind=rp), allocatable :: read_a(:), read_c(:), read_d(:)

  character(len=12) :: read_bcknd, read_subsolver

  call neko_init()

  if (pe_rank .eq. 0) then
     n = 4
  else
     n = 2
  end if
  m = 2

  call x%init(n)
  if (pe_rank .eq. 0) then
     x%x = [3.0_rp, 4.0_rp, 5.0_rp, 6.0_rp]
  else
     x%x = [1.0_rp, 2.0_rp]
  end if
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

  allocate(read_a(m), read_c(m), read_d(m))

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
  call h5aget_type_f(attr_id, str_type, ierr)
  call h5aread_f(attr_id, str_type, read_bcknd, ddim, ierr)
  call h5aclose_f(attr_id, ierr)

  call h5aopen_by_name_f(file_id, '/MMA/Parameters', 'subsolver', attr_id, ierr)
  call h5aget_type_f(attr_id, str_type, ierr)
  call h5aread_f(attr_id, str_type, read_subsolver, ddim, ierr)
  call h5aclose_f(attr_id, ierr)

  call h5tclose_f(str_type, ierr)

  ! Read penalty parameters
  ddim(1) = 1
  call h5dopen_f(file_id, '/MMA/a0', dset_id, ierr)
  call h5dread_f(dset_id, H5T_NEKO_REAL, read_a0, ddim, ierr)
  call h5dclose_f(dset_id, ierr)

  ddim(1) = m
  call h5dopen_f(file_id, '/MMA/a', dset_id, ierr)
  call h5dread_f(dset_id, H5T_NEKO_REAL, read_a, ddim, ierr)
  call h5dclose_f(dset_id, ierr)

  call h5dopen_f(file_id, '/MMA/c', dset_id, ierr)
  call h5dread_f(dset_id, H5T_NEKO_REAL, read_c, ddim, ierr)
  call h5dclose_f(dset_id, ierr)

  call h5dopen_f(file_id, '/MMA/d', dset_id, ierr)
  call h5dread_f(dset_id, H5T_NEKO_REAL, read_d, ddim, ierr)
  call h5dclose_f(dset_id, ierr)

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

  call test('n', read_n, n)
  call test('m', read_m, m)
  call test('max_iter', read_max_iter, max_iter)

  call test('asyinit', read_asyinit, asyinit)
  call test('asyincr', read_asyincr, asyincr)
  call test('asydecr', read_asydecr, asydecr)
  call test('epsimin', read_epsimin, epsimin)

  call test('bcknd', read_bcknd, bcknd)
  call test('subsolver', read_subsolver, subsolver)

  call test('a0', read_a0, a0)
  call test('a', read_a, a)
  call test('c', read_c, c)
  call test('d', read_d, d)

  call test('xold1', read_xold1%x, x%x)

  ! -------------------------------------------------------------------------- !
  ! Clean up

  call read_xold1%free()
  call obj%free()
  deallocate(a, c, d, xmin, xmax)

  call neko_finalize()

  if (test_status) then
     write(*,*) 'TEST PASSED'
  else
     write(*,*) 'TEST FAILED'
     error stop 1
  end if
end program test_mma_hdf5


#else

! Dummy program when HDF5 is not available
program test_mma_hdf5
  write(*,*) 'HDF5 not available; test skipped'
  stop
end program test_mma_hdf5

#endif
