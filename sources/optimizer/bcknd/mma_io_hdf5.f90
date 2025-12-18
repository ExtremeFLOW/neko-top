!> @file mma.f90
!! @copyright
!! Copyright (c) 2024-2025, The Neko-TOP Authors
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

!> Submodule for handling HDF5 IO for the mma object.
submodule (mma) mma_io_hdf5
#ifdef HAVE_HDF5
  use hdf5
#endif
  use mpi_f08, only: MPI_INFO_NULL, MPI_INTEGER8

contains

#ifdef HAVE_HDF5

  !> Write the MMA object to an HDF5 file (parallel-aware).
  !! This routine will first ensure device-host synchronization for all
  !! vectors and matrices, then perform the write. Currently the low-level
  !! HDF5 write is delegated to the project's I/O layer. If no I/O layer is
  !! available at link time a runtime error will be raised.
  module subroutine mma_write_hdf5(this, filename)
    class(mma_t), intent(inout) :: this
    character(len=*), intent(in) :: filename
    integer(hid_t) :: fapl_id, xf_id, file_id, dset_id, filespace, memspace, &
         attr_id, grp_id, mma_grp_id, str_type
    integer(hid_t) :: H5T_NEKO_REAL
    integer(hsize_t), dimension(1) :: ddim, dcount, doffset
    integer :: ierr, info, drank
    integer(kind=8) :: local_n8, prefix8, total_n8


    ! Ensure device state is on host
    call this%copy_from(DEVICE_TO_HOST, sync = .true.)

    call h5open_f(ierr)

    ! Assign the correct HDF5 data type based on the neko real kind
    select case (rp)
    case (dp)
       H5T_NEKO_REAL = H5T_NATIVE_DOUBLE
    case (sp)
       H5T_NEKO_REAL = H5T_NATIVE_REAL
    case default
       call neko_error('mma: unsupported real kind for HDF5')
    end select

    ! Create file with MPIO access
    call h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
    info = MPI_INFO_NULL%mpi_val
    call h5pset_fapl_mpio_f(fapl_id, NEKO_COMM%mpi_val, info, ierr)

    call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, ierr, &
         access_prp = fapl_id)
    call h5gcreate_f(file_id, "MMA", mma_grp_id, ierr, &
         lcpl_id = h5p_default_f, gcpl_id = h5p_default_f, &
         gapl_id = h5p_default_f)

    ! ------------------------------------------------------------------------ !
    ! Write basic Parameters attributes

    call h5gcreate_f(mma_grp_id, "Parameters", grp_id, ierr, &
         lcpl_id = h5p_default_f, gcpl_id = h5p_default_f, &
         gapl_id = h5p_default_f)

    call h5screate_f(H5S_SCALAR_F, filespace, ierr)
    ddim = 1

    ! Integer-valued attributes

    call h5acreate_f(grp_id, 'n', H5T_NATIVE_INTEGER, filespace, attr_id, &
         ierr, h5p_default_f, h5p_default_f)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, this%n, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5acreate_f(grp_id, 'm', H5T_NATIVE_INTEGER, filespace, attr_id, &
         ierr, h5p_default_f, h5p_default_f)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, this%m, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5acreate_f(grp_id, 'max_iter', H5T_NATIVE_INTEGER, filespace, &
         attr_id, ierr, h5p_default_f, h5p_default_f)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, this%max_iter, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    ! Real-valued attributes

    call h5acreate_f(grp_id, 'asyinit', H5T_NEKO_REAL, filespace, attr_id, &
         ierr, h5p_default_f, h5p_default_f)
    call h5awrite_f(attr_id, H5T_NEKO_REAL, this%asyinit, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5acreate_f(grp_id, 'asyincr', H5T_NEKO_REAL, filespace, attr_id, &
         ierr, h5p_default_f, h5p_default_f)
    call h5awrite_f(attr_id, H5T_NEKO_REAL, this%asyincr, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5acreate_f(grp_id, 'asydecr', H5T_NEKO_REAL, filespace, attr_id, &
         ierr, h5p_default_f, h5p_default_f)
    call h5awrite_f(attr_id, H5T_NEKO_REAL, this%asydecr, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5acreate_f(grp_id, 'epsimin', H5T_NEKO_REAL, filespace, attr_id, &
         ierr, h5p_default_f, h5p_default_f)
    call h5awrite_f(attr_id, H5T_NEKO_REAL, this%epsimin, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    ! String-valued attributes

    ! Create the string type
    call h5tcopy_f(H5T_FORTRAN_S1, str_type, ierr)
    call h5tset_strpad_f(str_type, H5T_STR_SPACEPAD_F, ierr)

    ddim(1) = len_trim(this%bcknd)
    call h5tset_size_f(str_type, ddim(1), ierr)
    call h5acreate_f(grp_id, 'bcknd', str_type, filespace, attr_id, &
         ierr, h5p_default_f, h5p_default_f)
    call h5awrite_f(attr_id, str_type, this%bcknd, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    ddim(1) = len_trim(this%subsolver)
    call h5tset_size_f(str_type, ddim(1), ierr)
    call h5acreate_f(grp_id, 'subsolver', str_type, filespace, attr_id, &
         ierr, h5p_default_f, h5p_default_f)
    call h5awrite_f(attr_id, str_type, this%subsolver, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    ! Close the string type
    call h5tclose_f(str_type, ierr)

    ! Close the filespace and Parameters group
    call h5sclose_f(filespace, ierr)
    call h5gclose_f(grp_id, ierr)

    ! ------------------------------------------------------------------------ !
    ! Global arrays datasets

    ddim(1) = 1
    call h5screate_f(H5S_SCALAR_F, filespace, ierr)

    call h5dcreate_f(mma_grp_id, 'a0', H5T_NEKO_REAL, filespace, dset_id, ierr)
    call h5dwrite_f(dset_id, H5T_NEKO_REAL, this%a0, ddim, ierr)
    call h5dclose_f(dset_id, ierr)

    ! The next batch are vectors of size m
    ddim(1) = this%m
    drank = 1

    call h5screate_simple_f(drank, ddim, filespace, ierr)

    call h5dcreate_f(mma_grp_id, 'a', H5T_NEKO_REAL, filespace, dset_id, ierr)
    call h5dwrite_f(dset_id, H5T_NEKO_REAL, this%a%x(1), ddim, ierr)
    call h5dclose_f(dset_id, ierr)

    call h5dcreate_f(mma_grp_id, 'c', H5T_NEKO_REAL, filespace, dset_id, ierr)
    call h5dwrite_f(dset_id, H5T_NEKO_REAL, this%c%x(1), ddim, ierr)
    call h5dclose_f(dset_id, ierr)

    call h5dcreate_f(mma_grp_id, 'd', H5T_NEKO_REAL, filespace, dset_id, ierr)
    call h5dwrite_f(dset_id, H5T_NEKO_REAL, this%d%x(1), ddim, ierr)
    call h5dclose_f(dset_id, ierr)

    ! Close the filespace and group
    call h5sclose_f(filespace, ierr)

    ! ------------------------------------------------------------------------ !
    ! Write per-rank datasets

    ! Define the sizes and offsets
    local_n8 = int(this%n, 8)
    total_n8 = int(this%n_global, 8)
    call MPI_Scan(local_n8, prefix8, 1, MPI_INTEGER8, MPI_SUM, NEKO_COMM, ierr)

    drank = 1
    dcount(1) = local_n8
    doffset(1) = prefix8 - local_n8
    ddim(1) = total_n8

    ! Create a file and memory space
    call h5screate_simple_f(drank, ddim, filespace, ierr)
    call h5screate_simple_f(drank, dcount, memspace, ierr)

    ! Dataset transfer property list: collective
    call h5pcreate_f(H5P_DATASET_XFER_F, xf_id, ierr)
    call h5pset_dxpl_mpio_f(xf_id, H5FD_MPIO_COLLECTIVE_F, ierr)

    ! Select hyperslab in the file space
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, doffset, dcount, &
         ierr)

    ! Write xmin
    call h5dcreate_f(mma_grp_id, 'xmin', H5T_NEKO_REAL, &
         filespace, dset_id, ierr)
    call h5dwrite_f(dset_id, H5T_NEKO_REAL, this%xmin%x(1), ddim, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = xf_id)
    call h5dclose_f(dset_id, ierr)

    ! Write xmax
    call h5dcreate_f(mma_grp_id, 'xmax', H5T_NEKO_REAL, &
         filespace, dset_id, ierr)
    call h5dwrite_f(dset_id, H5T_NEKO_REAL, this%xmax%x(1), ddim, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = xf_id)
    call h5dclose_f(dset_id, ierr)

    ! Write xold1
    call h5dcreate_f(mma_grp_id, 'xold1', H5T_NEKO_REAL, &
         filespace, dset_id, ierr)
    call h5dwrite_f(dset_id, H5T_NEKO_REAL, this%xold1%x(1), ddim, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = xf_id)
    call h5dclose_f(dset_id, ierr)

    ! Write xold2
    call h5dcreate_f(mma_grp_id, 'xold2', H5T_NEKO_REAL, &
         filespace, dset_id, ierr)
    call h5dwrite_f(dset_id, H5T_NEKO_REAL, this%xold2%x(1), ddim, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = xf_id)
    call h5dclose_f(dset_id, ierr)

    ! Write low
    call h5dcreate_f(mma_grp_id, 'low', H5T_NEKO_REAL, &
         filespace, dset_id, ierr)
    call h5dwrite_f(dset_id, H5T_NEKO_REAL, this%low%x(1), ddim, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = xf_id)
    call h5dclose_f(dset_id, ierr)

    ! Write upp
    call h5dcreate_f(mma_grp_id, 'upp', H5T_NEKO_REAL, &
         filespace, dset_id, ierr)
    call h5dwrite_f(dset_id, H5T_NEKO_REAL, this%upp%x(1), ddim, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = xf_id)
    call h5dclose_f(dset_id, ierr)

    ! Close the spaces used
    call h5sclose_f(filespace, ierr)
    call h5sclose_f(memspace, ierr)
    call h5pclose_f(xf_id, ierr)

    ! ------------------------------------------------------------------------ !
    ! Close the group and file

    call h5gclose_f(mma_grp_id, ierr)
    call h5fclose_f(file_id, ierr)
    call h5pclose_f(fapl_id, ierr)
    call h5close_f(ierr)

  end subroutine mma_write_hdf5

  !> Read the MMA object from an HDF5 file (parallel-aware).
  !! This routine will perform the read and then ensure device-host
  !! synchronization for all vectors and matrices.
  module subroutine mma_read_hdf5(this, filename)
    class(mma_t), intent(inout) :: this
    character(len=*), intent(in) :: filename
    integer(hid_t) :: fapl_id, xf_id, file_id, dset_id, filespace, memspace, &
         attr_id, grp_id, mma_grp_id, str_type
    integer(hid_t) :: H5T_NEKO_REAL
    integer(hsize_t), dimension(1) :: ddim, dcount, doffset
    integer :: ierr, info, drank
    integer(kind=8) :: local_n8, prefix8, total_n8
    integer :: n, m, max_iter
    real(kind=rp) :: asyinit, asyincr, asydecr, epsimin
    character(len=256) :: bcknd, subsolver

    ! Open file with MPIO access
    call h5open_f(ierr)

    ! Assign the correct HDF5 data type based on the neko real kind
    select case (rp)
    case (dp)
       H5T_NEKO_REAL = H5T_NATIVE_DOUBLE
    case (sp)
       H5T_NEKO_REAL = H5T_NATIVE_REAL
    case default
       call neko_error('mma: unsupported real kind for HDF5')
    end select

    call h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
    info = MPI_INFO_NULL%mpi_val
    call h5pset_fapl_mpio_f(fapl_id, NEKO_COMM%mpi_val, info, ierr)

    call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, ierr, &
         access_prp = fapl_id)
    call h5gopen_f(file_id, "MMA", mma_grp_id, ierr, gapl_id = h5p_default_f)

    ! ------------------------------------------------------------------------ !
    ! Read basic Parameters attributes

    call h5gopen_f(mma_grp_id, "Parameters", grp_id, ierr, &
         gapl_id = h5p_default_f)

    call h5screate_f(H5S_SCALAR_F, filespace, ierr)
    ddim = 1

    ! Integer-valued attributes

    call h5aopen_f(grp_id, 'n', attr_id, ierr, aapl_id = h5p_default_f)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, n, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5aopen_f(grp_id, 'm', attr_id, ierr, aapl_id = h5p_default_f)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, m, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5aopen_f(grp_id, 'max_iter', attr_id, ierr, aapl_id = h5p_default_f)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, max_iter, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    ! Real-valued attributes
    call h5aopen_f(grp_id, 'asyinit', attr_id, ierr, aapl_id = h5p_default_f)
    call h5aread_f(attr_id, H5T_NEKO_REAL, asyinit, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5aopen_f(grp_id, 'asyincr', attr_id, ierr, aapl_id = h5p_default_f)
    call h5aread_f(attr_id, H5T_NEKO_REAL, asyincr, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5aopen_f(grp_id, 'asydecr', attr_id, ierr, aapl_id = h5p_default_f)
    call h5aread_f(attr_id, H5T_NEKO_REAL, asydecr, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5aopen_f(grp_id, 'epsimin', attr_id, ierr, aapl_id = h5p_default_f)
    call h5aread_f(attr_id, H5T_NEKO_REAL, epsimin, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    ! String-valued attributes
    call h5aopen_f(grp_id, 'bcknd', attr_id, ierr, aapl_id = h5p_default_f)
    call h5aget_type_f(attr_id, str_type, ierr)
    call h5aread_f(attr_id, str_type, bcknd, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5aopen_f(grp_id, 'subsolver', attr_id, ierr, aapl_id = h5p_default_f)
    call h5aget_type_f(attr_id, str_type, ierr)
    call h5aread_f(attr_id, str_type, subsolver, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5tclose_f(str_type, ierr)
    call h5sclose_f(filespace, ierr)
    call h5gclose_f(grp_id, ierr)

    ! Ensure the MMA object is allocated with the same configuration as the file
    if (n .ne. this%n) then
       call neko_error('mma: mismatch in n during HDF5 read')
    end if
    if (m .ne. this%m) then
       call neko_error('mma: mismatch in m during HDF5 read')
    end if
    if (max_iter .ne. this%max_iter) then
       call neko_error('mma: mismatch in max_iter during HDF5 read')
    end if
    if (abs(asyinit - this%asyinit) .gt. 1.0e-12_rp) then
       call neko_error('mma: mismatch in asyinit during HDF5 read')
    end if
    if (abs(asyincr - this%asyincr) .gt. 1.0e-12_rp) then
       call neko_error('mma: mismatch in asyincr during HDF5 read')
    end if
    if (abs(asydecr - this%asydecr) .gt. 1.0e-12_rp) then
       call neko_error('mma: mismatch in asydecr during HDF5 read')
    end if
    if (abs(epsimin - this%epsimin) .gt. 1.0e-12_rp) then
       call neko_error('mma: mismatch in epsimin during HDF5 read')
    end if
    if (trim(bcknd) .ne. trim(this%bcknd)) then
       call neko_error('mma: mismatch in bcknd during HDF5 read')
    end if
    if (trim(subsolver) .ne. trim(this%subsolver)) then
       call neko_error('mma: mismatch in subsolver during HDF5 read')
    end if

    ! ------------------------------------------------------------------------ !
    ! Global arrays datasets

    ddim(1) = 1
    call h5screate_f(H5S_SCALAR_F, filespace, ierr)

    call h5dcreate_f(mma_grp_id, 'a0', H5T_NEKO_REAL, filespace, dset_id, ierr)
    call h5dread_f(dset_id, H5T_NEKO_REAL, this%a0, ddim, ierr)
    call h5dclose_f(dset_id, ierr)

    ! The next batch are vectors of size m
    ddim(1) = this%m
    drank = 1

    call h5screate_simple_f(drank, ddim, filespace, ierr)

    call h5dcreate_f(mma_grp_id, 'a', H5T_NEKO_REAL, filespace, dset_id, ierr)
    call h5dread_f(dset_id, H5T_NEKO_REAL, this%a%x(1), ddim, ierr)
    call h5dclose_f(dset_id, ierr)

    call h5dcreate_f(mma_grp_id, 'c', H5T_NEKO_REAL, filespace, dset_id, ierr)
    call h5dread_f(dset_id, H5T_NEKO_REAL, this%c%x(1), ddim, ierr)
    call h5dclose_f(dset_id, ierr)

    call h5dcreate_f(mma_grp_id, 'd', H5T_NEKO_REAL, filespace, dset_id, ierr)
    call h5dread_f(dset_id, H5T_NEKO_REAL, this%d%x(1), ddim, ierr)
    call h5dclose_f(dset_id, ierr)

    ! Close the filespace and group
    call h5sclose_f(filespace, ierr)

    ! ------------------------------------------------------------------------ !
    ! Write per-rank datasets

    ! Define the sizes and offsets
    local_n8 = int(this%n, 8)
    total_n8 = int(this%n_global, 8)
    call MPI_Scan(local_n8, prefix8, 1, MPI_INTEGER8, MPI_SUM, NEKO_COMM, ierr)

    drank = 1
    dcount(1) = local_n8
    doffset(1) = prefix8 - local_n8
    ddim(1) = total_n8

    ! Create a file and memory space
    call h5screate_simple_f(drank, ddim, filespace, ierr)
    call h5screate_simple_f(drank, dcount, memspace, ierr)

    ! Dataset transfer property list: collective
    call h5pcreate_f(H5P_DATASET_XFER_F, xf_id, ierr)
    call h5pset_dxpl_mpio_f(xf_id, H5FD_MPIO_COLLECTIVE_F, ierr)

    ! Select hyperslab in the file space
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, doffset, dcount, &
         ierr)

    ! Write xmin
    call h5dcreate_f(mma_grp_id, 'xmin', H5T_NEKO_REAL, &
         filespace, dset_id, ierr)
    call h5dread_f(dset_id, H5T_NEKO_REAL, this%xmin%x(1), ddim, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = xf_id)
    call h5dclose_f(dset_id, ierr)

    ! Write xmax
    call h5dcreate_f(mma_grp_id, 'xmax', H5T_NEKO_REAL, &
         filespace, dset_id, ierr)
    call h5dread_f(dset_id, H5T_NEKO_REAL, this%xmax%x(1), ddim, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = xf_id)
    call h5dclose_f(dset_id, ierr)

    ! Write xold1
    call h5dcreate_f(mma_grp_id, 'xold1', H5T_NEKO_REAL, &
         filespace, dset_id, ierr)
    call h5dread_f(dset_id, H5T_NEKO_REAL, this%xold1%x(1), ddim, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = xf_id)
    call h5dclose_f(dset_id, ierr)

    ! Write xold2
    call h5dcreate_f(mma_grp_id, 'xold2', H5T_NEKO_REAL, &
         filespace, dset_id, ierr)
    call h5dread_f(dset_id, H5T_NEKO_REAL, this%xold2%x(1), ddim, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = xf_id)
    call h5dclose_f(dset_id, ierr)

    ! Write low
    call h5dcreate_f(mma_grp_id, 'low', H5T_NEKO_REAL, &
         filespace, dset_id, ierr)
    call h5dread_f(dset_id, H5T_NEKO_REAL, this%low%x(1), ddim, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = xf_id)
    call h5dclose_f(dset_id, ierr)

    ! Write upp
    call h5dcreate_f(mma_grp_id, 'upp', H5T_NEKO_REAL, &
         filespace, dset_id, ierr)
    call h5dread_f(dset_id, H5T_NEKO_REAL, this%upp%x(1), ddim, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = xf_id)
    call h5dclose_f(dset_id, ierr)

    ! Close the spaces used
    call h5sclose_f(filespace, ierr)
    call h5sclose_f(memspace, ierr)
    call h5pclose_f(xf_id, ierr)

    ! ------------------------------------------------------------------------ !
    ! Close the group and file

    call h5gclose_f(mma_grp_id, ierr)
    call h5fclose_f(file_id, ierr)
    call h5pclose_f(fapl_id, ierr)
    call h5close_f(ierr)

    ! Ensure device state is updated
    call this%copy_from(HOST_TO_DEVICE, sync = .true.)

  end subroutine mma_read_hdf5

#else

  module subroutine mma_write_hdf5(this, filename)
    class(mma_t), intent(inout) :: this
    character(len=*), intent(in) :: filename
    call neko_error('mma: HDF5 support not enabled rebuild with HAVE_HDF5')
  end subroutine mma_write_hdf5

  module subroutine mma_read_hdf5(this, filename)
    class(mma_t), intent(inout) :: this
    character(len=*), intent(in) :: filename
    call neko_error('mma: HDF5 support not enabled rebuild with HAVE_HDF5')
  end subroutine mma_read_hdf5

#endif

end submodule mma_io_hdf5
