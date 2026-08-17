!> @file mma_hdf5_checkpoint.f90
!! @brief HDF5 IO submodule for the mma object.
!! @details
!! This submodule provides routines for saving and loading the mma
!! optimization object to and from HDF5 files in a parallel-aware manner.
!! @copyright
!! Copyright (c) 2024-2026, The Neko-TOP Authors
!! All rights reserved.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions
!! are met:
!!
!!   * Redistributions of source code must retain the above copyright
!!     notice, object list of conditions and the following disclaimer.
!!
!!   * Redistributions in binary form must reproduce the above
!!     copyright notice, object list of conditions and the following
!!     disclaimer in the documentation and/or other materials provided
!!     with the distribution.
!!
!!   * Neither the name of the authors nor the names of its
!!     contributors may be used to endorse or promote products derived
!!     from object software without specific prior written permission.
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
submodule (mma) mma_hdf5_io
  use hdf5
  use mpi_f08, only: MPI_Scan, MPI_INFO_NULL, MPI_INTEGER8, MPI_MAX, MPI_MIN, &
       MPI_IN_PLACE
  use math, only: abscmp
  use comm, only: pe_size

contains

  !> Write the MMA object to an HDF5 file (parallel-aware).
  !! This routine will first ensure device-host synchronization for all
  !! vectors and matrices, then perform the write. Currently the low-level
  !! HDF5 write is delegated to the project's I/O layer. If no I/O layer is
  !! available at link time a runtime error will be raised.
  !! @param object The MMA object to save.
  !! @param filename The HDF5 file to write to.
  !! @param overwrite Logical flag to allow overwriting existing files.
  module subroutine mma_save_checkpoint_hdf5(object, filename, overwrite)
    class(mma_t), intent(inout) :: object
    character(len=*), intent(in) :: filename
    logical, intent(in), optional :: overwrite
    integer(hid_t) :: fapl_id, xf_id, file_id, dset_id, filespace, memspace, &
         attr_id, grp_id, str_type, mma_grp_id
    integer(hid_t) :: H5T_NEKO_REAL
    integer(hsize_t), dimension(1) :: ddim, dcount, doffset
    integer :: ierr, info, drank
    logical :: file_exists, mma_exists, overwrite_flag
    integer :: n_accum, n_array(pe_size)

    overwrite_flag = .false.
    if (present(overwrite)) overwrite_flag = overwrite

    ! Ensure device state is on host
    call object%copy_from(DEVICE_TO_HOST, sync = .true.)

    ! Gather the per-rank n values
    n_array = -1
    n_array(pe_rank + 1) = object%n
    call MPI_Allreduce(MPI_IN_PLACE, n_array, pe_size, MPI_INTEGER, MPI_MAX, &
         NEKO_COMM, ierr)

    ! ------------------------------------------------------------------------ !
    ! Prepare the HDF5 file, settings for MPIO access and groups

    call h5open_f(ierr)
    call h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
    info = MPI_INFO_NULL%mpi_val
    call h5pset_fapl_mpio_f(fapl_id, NEKO_COMM%mpi_val, info, ierr)

    ! Handle overwriting if the file exists
    file_exists = .false.
    inquire(file=trim(filename), exist=file_exists)
    if (.not. file_exists) then
       call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, ierr, &
            access_prp = fapl_id)

       ! Create the MMA/checkpoint group
       call h5gcreate_f(file_id, "MMA", mma_grp_id, ierr)
       call h5gcreate_f(mma_grp_id, "checkpoint", grp_id, ierr)

    else
       call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, ierr, fapl_id)

       ! Check for existing MMA group
       mma_exists = .false.
       call h5lexists_f(file_id, '/MMA', mma_exists, ierr)

       if (.not. mma_exists) then
          ! Create the MMA/checkpoint group
          call h5gcreate_f(file_id, "MMA", mma_grp_id, ierr)
          call h5gcreate_f(mma_grp_id, "checkpoint", grp_id, ierr)
       else
          ! Open the existing "MMA" group
          call h5gopen_f(file_id, 'MMA', mma_grp_id, ierr)
          call h5lexists_f(mma_grp_id, 'checkpoint', mma_exists, ierr)

          if (mma_exists .and. overwrite_flag) then
             call h5gunlink_f(mma_grp_id, "checkpoint", ierr)
             call h5gcreate_f(mma_grp_id, "checkpoint", grp_id, ierr)
          else
             call neko_error('mma: HDF5 file "' // trim(filename) // '" ' // &
                  'already contains MMA group; use overwrite option to replace')
          end if
       end if
    end if
    call h5gclose_f(mma_grp_id, ierr)

    ! Assign the correct HDF5 data type based on the neko real kind
    select case (rp)
    case (dp)
       H5T_NEKO_REAL = H5T_NATIVE_DOUBLE
    case (sp)
       H5T_NEKO_REAL = H5T_NATIVE_REAL
    case default
       call neko_error('mma: unsupported real kind for HDF5')
    end select

    ! ------------------------------------------------------------------------ !
    ! Write basic Parameters attributes

    ! Integer-valued attributes
    ddim = 1
    call h5screate_f(H5S_SCALAR_F, filespace, ierr)

    call h5acreate_f(grp_id, 'n_global', H5T_NATIVE_INTEGER, filespace, &
         attr_id, ierr, h5p_default_f, h5p_default_f)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, object%n_global, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5acreate_f(grp_id, 'm', H5T_NATIVE_INTEGER, filespace, attr_id, &
         ierr, h5p_default_f, h5p_default_f)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, object%m, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5acreate_f(grp_id, 'max_iter', H5T_NATIVE_INTEGER, filespace, &
         attr_id, ierr, h5p_default_f, h5p_default_f)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, object%max_iter, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    ! Real-valued attributes

    call h5acreate_f(grp_id, 'asyinit', H5T_NEKO_REAL, filespace, attr_id, &
         ierr, h5p_default_f, h5p_default_f)
    call h5awrite_f(attr_id, H5T_NEKO_REAL, object%asyinit, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5acreate_f(grp_id, 'asyincr', H5T_NEKO_REAL, filespace, attr_id, &
         ierr, h5p_default_f, h5p_default_f)
    call h5awrite_f(attr_id, H5T_NEKO_REAL, object%asyincr, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5acreate_f(grp_id, 'asydecr', H5T_NEKO_REAL, filespace, attr_id, &
         ierr, h5p_default_f, h5p_default_f)
    call h5awrite_f(attr_id, H5T_NEKO_REAL, object%asydecr, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5acreate_f(grp_id, 'epsimin', H5T_NEKO_REAL, filespace, attr_id, &
         ierr, h5p_default_f, h5p_default_f)
    call h5awrite_f(attr_id, H5T_NEKO_REAL, object%epsimin, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5acreate_f(grp_id, 'move_limit', H5T_NEKO_REAL, filespace, attr_id, &
         ierr, h5p_default_f, h5p_default_f)
    call h5awrite_f(attr_id, H5T_NEKO_REAL, object%move_limit, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    ! String-valued attributes

    ! Create the string type
    call h5tcopy_f(H5T_FORTRAN_S1, str_type, ierr)
    call h5tset_strpad_f(str_type, H5T_STR_SPACEPAD_F, ierr)

    ddim(1) = len_trim(object%bcknd)
    call h5tset_size_f(str_type, ddim(1), ierr)
    call h5acreate_f(grp_id, 'bcknd', str_type, filespace, attr_id, &
         ierr, h5p_default_f, h5p_default_f)
    call h5awrite_f(attr_id, str_type, object%bcknd, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    ddim(1) = len_trim(object%subsolver)
    call h5tset_size_f(str_type, ddim(1), ierr)
    call h5acreate_f(grp_id, 'subsolver', str_type, filespace, attr_id, &
         ierr, h5p_default_f, h5p_default_f)
    call h5awrite_f(attr_id, str_type, object%subsolver, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    ! Close the string type
    call h5tclose_f(str_type, ierr)

    ! Close the scalar filespace
    call h5sclose_f(filespace, ierr)

    ! Array-valued attributes
    drank = 1
    ddim = pe_size
    call h5screate_simple_f(drank, ddim, filespace, ierr)

    call h5acreate_f(grp_id, 'n', H5T_NATIVE_INTEGER, filespace, attr_id, ierr)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, n_array, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5sclose_f(filespace, ierr)

    ! ------------------------------------------------------------------------ !
    ! Global arrays datasets

    ddim(1) = 1
    call h5screate_f(H5S_SCALAR_F, filespace, ierr)

    call h5dcreate_f(grp_id, 'a0', H5T_NEKO_REAL, filespace, dset_id, ierr)
    call h5dwrite_f(dset_id, H5T_NEKO_REAL, object%a0, ddim, ierr)
    call h5dclose_f(dset_id, ierr)

    ! The next batch are vectors of size m
    ddim(1) = object%m
    drank = 1

    call h5screate_simple_f(drank, ddim, filespace, ierr)

    call h5dcreate_f(grp_id, 'a', H5T_NEKO_REAL, filespace, dset_id, ierr)
    call h5dwrite_f(dset_id, H5T_NEKO_REAL, object%a%x, ddim, ierr)
    call h5dclose_f(dset_id, ierr)

    call h5dcreate_f(grp_id, 'c', H5T_NEKO_REAL, filespace, dset_id, ierr)
    call h5dwrite_f(dset_id, H5T_NEKO_REAL, object%c%x, ddim, ierr)
    call h5dclose_f(dset_id, ierr)

    call h5dcreate_f(grp_id, 'd', H5T_NEKO_REAL, filespace, dset_id, ierr)
    call h5dwrite_f(dset_id, H5T_NEKO_REAL, object%d%x, ddim, ierr)
    call h5dclose_f(dset_id, ierr)

    ! Close the filespace and group
    call h5sclose_f(filespace, ierr)

    ! ------------------------------------------------------------------------ !
    ! Write per-rank datasets

    ! Define the sizes and offsets
    call MPI_Scan(object%n, n_accum, 1, MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)

    drank = 1
    dcount(1) = object%n
    doffset(1) = n_accum - object%n
    ddim(1) = object%n_global

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
    call h5dcreate_f(grp_id, 'xmin', H5T_NEKO_REAL, &
         filespace, dset_id, ierr)
    call h5dwrite_f(dset_id, H5T_NEKO_REAL, object%xmin%x, ddim, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = xf_id)
    call h5dclose_f(dset_id, ierr)

    ! Write xmax
    call h5dcreate_f(grp_id, 'xmax', H5T_NEKO_REAL, &
         filespace, dset_id, ierr)
    call h5dwrite_f(dset_id, H5T_NEKO_REAL, object%xmax%x, ddim, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = xf_id)
    call h5dclose_f(dset_id, ierr)

    ! Write xold1
    call h5dcreate_f(grp_id, 'xold1', H5T_NEKO_REAL, &
         filespace, dset_id, ierr)
    call h5dwrite_f(dset_id, H5T_NEKO_REAL, object%xold1%x, ddim, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = xf_id)
    call h5dclose_f(dset_id, ierr)

    ! Write xold2
    call h5dcreate_f(grp_id, 'xold2', H5T_NEKO_REAL, &
         filespace, dset_id, ierr)
    call h5dwrite_f(dset_id, H5T_NEKO_REAL, object%xold2%x, ddim, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = xf_id)
    call h5dclose_f(dset_id, ierr)

    ! Write low
    call h5dcreate_f(grp_id, 'low', H5T_NEKO_REAL, &
         filespace, dset_id, ierr)
    call h5dwrite_f(dset_id, H5T_NEKO_REAL, object%low%x, ddim, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = xf_id)
    call h5dclose_f(dset_id, ierr)

    ! Write upp
    call h5dcreate_f(grp_id, 'upp', H5T_NEKO_REAL, &
         filespace, dset_id, ierr)
    call h5dwrite_f(dset_id, H5T_NEKO_REAL, object%upp%x, ddim, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = xf_id)
    call h5dclose_f(dset_id, ierr)

    ! Close the spaces used
    call h5sclose_f(filespace, ierr)
    call h5sclose_f(memspace, ierr)
    call h5pclose_f(xf_id, ierr)

    ! ------------------------------------------------------------------------ !
    ! Close the group and file

    call h5gclose_f(grp_id, ierr)
    call h5fclose_f(file_id, ierr)
    call h5pclose_f(fapl_id, ierr)
    call h5close_f(ierr)

  end subroutine mma_save_checkpoint_hdf5

  !> Read the MMA object from an HDF5 file (parallel-aware).
  !! This routine will perform the read and then ensure device-host
  !! synchronization for all vectors and matrices.
  !! @param object The MMA object to load data into.
  !! @param filename The HDF5 file to read from.
  module subroutine mma_load_checkpoint_hdf5(object, filename)
    class(mma_t), intent(inout) :: object
    character(len=*), intent(in) :: filename
    integer(hid_t) :: fapl_id, file_id, dset_id, attr_id, grp_id
    integer(hid_t) :: str_type, filespace, memspace, xf_id
    integer(hid_t) :: H5T_NEKO_REAL
    integer(hsize_t), dimension(1) :: ddim, dcount, doffset
    integer :: ierr, info, mpi_comm
    integer :: n, n_global, m, max_iter, n_accum, n_array(pe_size)
    real(kind=rp) :: asyinit, asyincr, asydecr, epsimin, move_limit
    logical :: move_limit_exists

    character(len=*), parameter :: h5_group = '/MMA/checkpoint'
    character(len=12) :: bcknd, subsolver

    ! Initialize reader values
    n = -1
    m = -1
    max_iter = -1
    asyinit = -1.0_rp
    asyincr = -1.0_rp
    asydecr = -1.0_rp
    epsimin = -1.0_rp
    move_limit = object%move_limit

    bcknd = ''
    subsolver = ''

    ! MPI interfaces
    info = MPI_INFO_NULL%mpi_val
    mpi_comm = NEKO_COMM%mpi_val

    ! Open file and prepare reading
    call h5open_f(ierr)
    call h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
    call h5pset_fapl_mpio_f(fapl_id, mpi_comm, info, ierr)
    call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, ierr, &
         access_prp = fapl_id)
    call h5gopen_f(file_id, h5_group, grp_id, ierr)

    if (ierr .ne. 0) then
       call neko_error('mma: unable to open HDF5 file "' // trim(filename) // &
            '" or group "' // trim(h5_group) // '".')
    end if

    ! Assign the correct HDF5 data type based on the neko real kind
    select case (rp)
    case (dp)
       H5T_NEKO_REAL = H5T_NATIVE_DOUBLE
    case (sp)
       H5T_NEKO_REAL = H5T_NATIVE_REAL
    case default
       call neko_error('mma: unsupported real kind for HDF5')
    end select

    ! ------------------------------------------------------------------------ !
    ! Read basic Parameters attributes

    ! Read scalar attributes
    ddim(1) = 1
    call h5aopen_f(grp_id, 'n_global', attr_id, ierr)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, n_global, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5aopen_f(grp_id, 'm', attr_id, ierr)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, m, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5aopen_f(grp_id, 'max_iter', attr_id, ierr)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, max_iter, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5aopen_f(grp_id, 'asyinit', attr_id, ierr)
    call h5aread_f(attr_id, H5T_NEKO_REAL, asyinit, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5aopen_f(grp_id, 'asyincr', attr_id, ierr)
    call h5aread_f(attr_id, H5T_NEKO_REAL, asyincr, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5aopen_f(grp_id, 'asydecr', attr_id, ierr)
    call h5aread_f(attr_id, H5T_NEKO_REAL, asydecr, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5aopen_f(grp_id, 'epsimin', attr_id, ierr)
    call h5aread_f(attr_id, H5T_NEKO_REAL, epsimin, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    move_limit_exists = .false.
    call h5aexists_f(grp_id, 'move_limit', move_limit_exists, ierr)
    if (move_limit_exists) then
       call h5aopen_f(grp_id, 'move_limit', attr_id, ierr)
       call h5aread_f(attr_id, H5T_NEKO_REAL, move_limit, ddim, ierr)
       call h5aclose_f(attr_id, ierr)
    end if

    ! Read strings
    call h5aopen_f(grp_id, 'bcknd', attr_id, ierr)
    call h5aget_type_f(attr_id, str_type, ierr)
    call h5aread_f(attr_id, str_type, bcknd, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5aopen_f(grp_id, 'subsolver', attr_id, ierr)
    call h5aget_type_f(attr_id, str_type, ierr)
    call h5aread_f(attr_id, str_type, subsolver, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5tclose_f(str_type, ierr)

    ! Read array attribute n
    ddim(1) = pe_size
    n_array = -1

    call h5aopen_f(grp_id, 'n', attr_id, ierr)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, n_array, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    ! ------------------------------------------------------------------------ !
    ! Ensure the MMA object is allocated with the same configuration as the file

    if (n_array(pe_rank + 1) .ne. object%n) then
       call neko_error('mma: mismatch in n during HDF5 read')
    end if
    if (n_global .ne. object%n_global) then
       call neko_error('mma: mismatch in n_global during HDF5 read')
    end if
    if (m .ne. object%m) then
       call neko_error('mma: mismatch in m during HDF5 read')
    end if
    if (max_iter .ne. object%max_iter) then
       call neko_error('mma: mismatch in max_iter during HDF5 read')
    end if
    if (.not. abscmp(asyinit, object%asyinit)) then
       call neko_error('mma: mismatch in asyinit during HDF5 read')
    end if
    if (.not. abscmp(asyincr, object%asyincr)) then
       call neko_error('mma: mismatch in asyincr during HDF5 read')
    end if
    if (.not. abscmp(asydecr, object%asydecr)) then
       call neko_error('mma: mismatch in asydecr during HDF5 read')
    end if
    if (.not. abscmp(epsimin, object%epsimin)) then
       call neko_error('mma: mismatch in epsimin during HDF5 read')
    end if
    if (move_limit_exists) then
       if (.not. abscmp(move_limit, object%move_limit)) then
          call neko_error('mma: mismatch in move_limit during HDF5 read')
       end if
    end if
    if (trim(bcknd) .ne. trim(object%bcknd)) then
       call neko_error('mma: mismatch in bcknd during HDF5 read')
    end if
    if (trim(subsolver) .ne. trim(object%subsolver)) then
       call neko_error('mma: mismatch in subsolver during HDF5 read')
    end if

    ! ------------------------------------------------------------------------ !
    ! Global arrays datasets

    ! Read penalty parameters
    ddim(1) = 1
    call h5dopen_f(grp_id, 'a0', dset_id, ierr)
    call h5dread_f(dset_id, H5T_NEKO_REAL, object%a0, ddim, ierr)
    call h5dclose_f(dset_id, ierr)

    ddim(1) = m
    call h5dopen_f(grp_id, 'a', dset_id, ierr)
    call h5dread_f(dset_id, H5T_NEKO_REAL, object%a%x, ddim, ierr)
    call h5dclose_f(dset_id, ierr)

    call h5dopen_f(grp_id, 'c', dset_id, ierr)
    call h5dread_f(dset_id, H5T_NEKO_REAL, object%c%x, ddim, ierr)
    call h5dclose_f(dset_id, ierr)

    call h5dopen_f(grp_id, 'd', dset_id, ierr)
    call h5dread_f(dset_id, H5T_NEKO_REAL, object%d%x, ddim, ierr)
    call h5dclose_f(dset_id, ierr)

    ! ------------------------------------------------------------------------ !
    ! Read per-rank datasets

    ddim(1) = object%n

    ! Define the sizes and offsets (each pe_rank reads its local n from the global array)
    call MPI_Scan(object%n, n_accum, 1, MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)

    dcount(1) = object%n
    doffset(1) = n_accum - object%n
    ddim(1) = object%n_global

    ! Create file and memory space (pe_rank = 1 for 1D arrays)
    call h5screate_simple_f(1, ddim, filespace, ierr)
    call h5screate_simple_f(1, dcount, memspace, ierr)

    ! Dataset transfer property list: MPI collective I/O
    call h5pcreate_f(H5P_DATASET_XFER_F, xf_id, ierr)
    call h5pset_dxpl_mpio_f(xf_id, H5FD_MPIO_COLLECTIVE_F, ierr)

    ! Select hyperslab in the file space for object pe_rank
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, doffset, dcount, &
         ierr)

    call h5dopen_f(grp_id, 'xmin', dset_id, ierr)
    call h5dread_f(dset_id, H5T_NEKO_REAL, object%xmin%x, dcount, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = xf_id)
    call h5dclose_f(dset_id, ierr)

    call h5dopen_f(grp_id, 'xmax', dset_id, ierr)
    call h5dread_f(dset_id, H5T_NEKO_REAL, object%xmax%x, dcount, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = xf_id)
    call h5dclose_f(dset_id, ierr)

    call h5dopen_f(grp_id, 'xold1', dset_id, ierr)
    call h5dread_f(dset_id, H5T_NEKO_REAL, object%xold1%x, dcount, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = xf_id)
    call h5dclose_f(dset_id, ierr)

    call h5dopen_f(grp_id, 'xold2', dset_id, ierr)
    call h5dread_f(dset_id, H5T_NEKO_REAL, object%xold2%x, dcount, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = xf_id)
    call h5dclose_f(dset_id, ierr)

    call h5dopen_f(grp_id, 'low', dset_id, ierr)
    call h5dread_f(dset_id, H5T_NEKO_REAL, object%low%x, dcount, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = xf_id)
    call h5dclose_f(dset_id, ierr)

    call h5dopen_f(grp_id, 'upp', dset_id, ierr)
    call h5dread_f(dset_id, H5T_NEKO_REAL, object%upp%x, dcount, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = xf_id)
    call h5dclose_f(dset_id, ierr)

    ! Close transfer property list and dataspaces
    call h5pclose_f(xf_id, ierr)
    call h5sclose_f(memspace, ierr)
    call h5sclose_f(filespace, ierr)

    ! ------------------------------------------------------------------------ !
    ! Close the group and file

    call h5gclose_f(grp_id, ierr)
    call h5fclose_f(file_id, ierr)
    call h5pclose_f(fapl_id, ierr)
    call h5close_f(ierr)

    ! Ensure device state is updated
    call object%copy_from(HOST_TO_DEVICE, sync = .true.)

  end subroutine mma_load_checkpoint_hdf5

end submodule mma_hdf5_io
