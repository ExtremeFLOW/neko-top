!> @file design_hdf5_checkpoint.f90
!! @brief HDF5 IO submodule for the design object.
!! @details
!! This submodule provides routines for saving and loading the design
!! object to and from HDF5 files in a parallel-aware manner.
!! @copyright
!! Copyright (c) 2024-2026, The Neko-TOP Authors
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

!> Submodule for handling HDF5 IO for the design object.
submodule (design) design_hdf5_io
  use hdf5
  use mpi_f08, only: MPI_Scan, MPI_INFO_NULL, MPI_INTEGER8, MPI_MAX, MPI_MIN, &
       MPI_IN_PLACE
  use math, only: abscmp
  use device, only: DEVICE_TO_HOST, HOST_TO_DEVICE
  use comm, only: pe_rank, pe_size
  use num_types, only: sp, dp
  use scratch_registry, only: neko_scratch_registry

contains

  !> Write the Design object to an HDF5 file (parallel-aware).
  !! This routine will first ensure device-host synchronization for all
  !! vectors and matrices, then perform the write. Currently the low-level
  !! HDF5 write is delegated to the project's I/O layer. If no I/O layer is
  !! available at link time a runtime error will be raised.
  module subroutine design_save_checkpoint_hdf5(this, filename, overwrite)
    class(design_t), intent(in) :: this
    character(len=*), intent(in) :: filename
    logical, intent(in), optional :: overwrite
    type(vector_t), pointer :: temp_vec
    integer :: temp_index

    integer(hid_t) :: fapl_id, xf_id, file_id, dset_id, filespace, memspace, &
         attr_id, grp_id, str_type, design_grp_id
    integer(hid_t) :: H5T_NEKO_REAL
    integer(hsize_t), dimension(1) :: ddim, dcount, doffset
    integer :: ierr, info, drank
    logical :: file_exists, design_exists, overwrite_flag
    integer :: n_accum, n_array(pe_size)

    overwrite_flag = .false.
    if (present(overwrite)) overwrite_flag = overwrite

    ! Gather the per-rank n values
    n_array = -1
    n_array(pe_rank + 1) = this%n
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

       ! Create the Design/checkpoint group
       call h5gcreate_f(file_id, "Design", design_grp_id, ierr)
       call h5gcreate_f(design_grp_id, "checkpoint", grp_id, ierr)

    else
       call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, ierr, fapl_id)

       ! Check for existing Design group
       design_exists = .false.
       call h5lexists_f(file_id, '/Design', design_exists, ierr)

       if (.not. design_exists) then
          ! Create the Design/checkpoint group
          call h5gcreate_f(file_id, "Design", design_grp_id, ierr)
          call h5gcreate_f(design_grp_id, "checkpoint", grp_id, ierr)
       else
          ! Open the existing "Design" group
          call h5gopen_f(file_id, 'Design', design_grp_id, ierr)
          call h5lexists_f(design_grp_id, 'checkpoint', design_exists, ierr)

          if (design_exists .and. overwrite_flag) then
             call h5gunlink_f(design_grp_id, "checkpoint", ierr)
             call h5gcreate_f(design_grp_id, "checkpoint", grp_id, ierr)
          else
             call neko_error('design: HDF5 file "' // trim(filename) // '" ' &
                  // 'already contains Design group; ' &
                  // 'use overwrite option to replace')
          end if
       end if
    end if

    ! Assign the correct HDF5 data type based on the neko real kind
    select case (rp)
    case (dp)
       H5T_NEKO_REAL = H5T_NATIVE_DOUBLE
    case (sp)
       H5T_NEKO_REAL = H5T_NATIVE_REAL
    case default
       call neko_error('design: unsupported real kind for HDF5')
    end select

    ! ------------------------------------------------------------------------ !
    ! Write basic Parameters attributes

    ! Integer-valued attributes
    ddim = 1
    call h5screate_f(H5S_SCALAR_F, filespace, ierr)

    call h5acreate_f(grp_id, 'n_global', H5T_NATIVE_INTEGER, filespace, &
         attr_id, ierr, h5p_default_f, h5p_default_f)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, this%n_global, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    ! String-valued attributes

    ! Create the string type
    call h5tcopy_f(H5T_FORTRAN_S1, str_type, ierr)
    call h5tset_strpad_f(str_type, H5T_STR_SPACEPAD_F, ierr)

    ddim(1) = len_trim(this%name)
    call h5tset_size_f(str_type, ddim(1), ierr)
    call h5acreate_f(grp_id, 'name', str_type, filespace, attr_id, &
         ierr, h5p_default_f, h5p_default_f)
    call h5awrite_f(attr_id, str_type, this%name, ddim, ierr)
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
    ! Write per-rank datasets

    ! Define the sizes and offsets
    call MPI_Scan(this%n, n_accum, 1, MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)

    drank = 1
    dcount(1) = this%n
    doffset(1) = n_accum - this%n
    ddim(1) = this%n_global

    ! Create a file and memory space
    call h5screate_simple_f(drank, ddim, filespace, ierr)
    call h5screate_simple_f(drank, dcount, memspace, ierr)

    ! Dataset transfer property list: collective
    call h5pcreate_f(H5P_DATASET_XFER_F, xf_id, ierr)
    call h5pset_dxpl_mpio_f(xf_id, H5FD_MPIO_COLLECTIVE_F, ierr)

    ! Select hyperslab in the file space
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, doffset, dcount, &
         ierr)

    call neko_scratch_registry%request(temp_vec, temp_index, this%n, .true.)
    call this%get_values(temp_vec)
    call temp_vec%copy_from(DEVICE_TO_HOST, sync = .true.)

    ! Write values
    call h5dcreate_f(grp_id, 'values', H5T_NEKO_REAL, &
         filespace, dset_id, ierr)
    call h5dwrite_f(dset_id, H5T_NEKO_REAL, temp_vec%x, ddim, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = xf_id)
    call h5dclose_f(dset_id, ierr)

    call neko_scratch_registry%relinquish(temp_index)
    nullify(temp_vec)

    ! Close the spaces used
    call h5sclose_f(filespace, ierr)
    call h5sclose_f(memspace, ierr)
    call h5pclose_f(xf_id, ierr)

    ! ------------------------------------------------------------------------ !
    ! Close the group and file

    call h5gclose_f(grp_id, ierr)
    call h5gclose_f(design_grp_id, ierr)
    call h5fclose_f(file_id, ierr)
    call h5pclose_f(fapl_id, ierr)
    call h5close_f(ierr)

  end subroutine design_save_checkpoint_hdf5

  !> Read the Design object from an HDF5 file (parallel-aware).
  !! This routine will perform the read and then ensure device-host
  !! synchronization for all vectors and matrices.
  module subroutine design_load_checkpoint_hdf5(this, filename)
    class(design_t), intent(inout) :: this
    character(len=*), intent(in) :: filename

    type(vector_t), pointer :: temp_vec
    integer :: temp_index

    character(len=64) :: name
    integer(hid_t) :: fapl_id, file_id, dset_id, &
         attr_id, grp_id, str_type, filespace, memspace, xf_id
    integer(hid_t) :: H5T_NEKO_REAL
    integer(hsize_t), dimension(1) :: ddim, dcount, doffset
    integer :: ierr, info, mpi_comm
    integer :: n, n_global, n_accum, n_array(pe_size)

    character(len=*), parameter :: h5_group = '/Design/checkpoint'

    ! Initialize reader values
    n = -1
    n_global = -1
    name = ''

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

    ! Assign the correct HDF5 data type based on the neko real kind
    select case (rp)
    case (dp)
       H5T_NEKO_REAL = H5T_NATIVE_DOUBLE
    case (sp)
       H5T_NEKO_REAL = H5T_NATIVE_REAL
    case default
       call neko_error('design: unsupported real kind for HDF5')
    end select

    ! ------------------------------------------------------------------------ !
    ! Read basic Parameters attributes

    ! Read scalar attributes
    ddim(1) = 1
    call h5aopen_f(grp_id, 'n_global', attr_id, ierr)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, n_global, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    ! Read strings
    call h5aopen_f(grp_id, 'name', attr_id, ierr)
    call h5aget_type_f(attr_id, str_type, ierr)
    call h5aread_f(attr_id, str_type, name, ddim, ierr)
    call h5tclose_f(str_type, ierr)
    call h5aclose_f(attr_id, ierr)

    ! Read array attribute n
    ddim(1) = pe_size
    n_array = -1

    call h5aopen_f(grp_id, 'n', attr_id, ierr)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, n_array, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    ! ------------------------------------------------------------------------ !
    ! Ensure the Design object is allocated with the same configuration as the file

    if (n_array(pe_rank + 1) .ne. this%n) then
       call neko_error('design: mismatch in n during HDF5 read')
    end if
    if (n_global .ne. this%n_global) then
       call neko_error('design: mismatch in n_global during HDF5 read')
    end if
    if (trim(name) .ne. trim(this%name)) then
       call neko_error('design: mismatch in name during HDF5 read' // &
            ' (file: ' // trim(name) // ', object: ' // trim(this%name) // ')')
    end if

    ! ------------------------------------------------------------------------ !
    ! Read per-rank datasets

    ddim(1) = this%n

    ! Define the sizes and offsets (each pe_rank reads its local n from the global array)
    call MPI_Scan(this%n, n_accum, 1, MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)

    dcount(1) = this%n
    doffset(1) = n_accum - this%n
    ddim(1) = this%n_global

    ! Create file and memory space (pe_rank = 1 for 1D arrays)
    call h5screate_simple_f(1, ddim, filespace, ierr)
    call h5screate_simple_f(1, dcount, memspace, ierr)

    ! Dataset transfer property list: MPI collective I/O
    call h5pcreate_f(H5P_DATASET_XFER_F, xf_id, ierr)
    call h5pset_dxpl_mpio_f(xf_id, H5FD_MPIO_COLLECTIVE_F, ierr)

    ! Select hyperslab in the file space for this pe_rank
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, doffset, dcount, &
         ierr)

    ! Read values
    call neko_scratch_registry%request(temp_vec, temp_index, this%n, .false.)
    call h5dopen_f(grp_id, 'values', dset_id, ierr)
    call h5dread_f(dset_id, H5T_NEKO_REAL, temp_vec%x, dcount, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = xf_id)
    call h5dclose_f(dset_id, ierr)

    call temp_vec%copy_from(HOST_TO_DEVICE, sync = .true.)
    call this%update_design(temp_vec)
    call neko_scratch_registry%relinquish(temp_index)
    nullify(temp_vec)

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

  end subroutine design_load_checkpoint_hdf5

end submodule design_hdf5_io
