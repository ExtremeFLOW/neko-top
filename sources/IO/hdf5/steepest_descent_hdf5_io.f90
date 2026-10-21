!> @file steepest_descent_hdf5_io.f90
!! @brief HDF5 IO submodule for the steepest descent optimizer.
!! @details
!! This submodule provides routines for saving and loading the state of the
!! steepest descent optimizer to and from HDF5 files in a parallel-aware
!! manner.
!! @copyright
!! Copyright (c) 2026, The Neko-TOP Authors
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

!> Submodule for handling HDF5 IO for the steepest descent optimizer.
!! @details
!! The state is stored in the group `/SteepestDescent/checkpoint`: the
!! attributes `n_global` and `n` (the local size of every rank), and the
!! distributed datasets `xold1`, `xold2` and `delta`. The settings are not
!! stored; on a restart the case file defines them.
submodule (steepest_descent_optimizer) steepest_descent_hdf5_io
  use hdf5
  use num_types, only: dp, sp
  use mpi_f08, only: MPI_Allreduce, MPI_Scan, MPI_IN_PLACE, MPI_INFO_NULL, &
       MPI_INTEGER, MPI_SUM, MPI_MAX
  use comm, only: NEKO_COMM, pe_rank, pe_size
  use device, only: DEVICE_TO_HOST, HOST_TO_DEVICE

  implicit none

  !> Name of the checkpoint group.
  character(len=*), parameter :: SD_GROUP = '/SteepestDescent/checkpoint'

contains

  !> Save the steepest descent state to an HDF5 checkpoint.
  !! @details The file must already exist; the base optimizer creates it
  !! before saving the components of the specific optimizer.
  !! @param object The steepest descent optimizer.
  !! @param filename The HDF5 file to write to.
  !! @param overwrite Replace an existing steepest descent group.
  module subroutine sd_optimizer_save_checkpoint_hdf5(object, filename, &
       overwrite)
    class(steepest_descent_optimizer_t), intent(inout) :: object
    character(len=*), intent(in) :: filename
    logical, intent(in), optional :: overwrite
    integer(hid_t) :: fapl_id, file_id, sd_grp_id, grp_id, filespace, attr_id
    integer(hsize_t) :: ddim(1)
    integer :: ierr, info, n, n_global, n_accum, n_array(pe_size)
    logical :: overwrite_flag, link_exists

    overwrite_flag = .false.
    if (present(overwrite)) overwrite_flag = overwrite

    ! Ensure device state is on host
    call object%xold1%copy_from(DEVICE_TO_HOST, .true.)
    call object%xold2%copy_from(DEVICE_TO_HOST, .true.)
    call object%delta%copy_from(DEVICE_TO_HOST, .true.)

    ! Gather the per-rank sizes and the offset of this rank
    n = object%xold1%size()
    n_array = 0
    n_array(pe_rank + 1) = n
    call MPI_Allreduce(MPI_IN_PLACE, n_array, pe_size, MPI_INTEGER, MPI_SUM, &
         NEKO_COMM, ierr)
    n_global = sum(n_array)
    call MPI_Scan(n, n_accum, 1, MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)

    ! ------------------------------------------------------------------------ !
    ! Open the file and prepare the checkpoint group

    call h5open_f(ierr)
    call h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
    info = MPI_INFO_NULL%mpi_val
    call h5pset_fapl_mpio_f(fapl_id, NEKO_COMM%mpi_val, info, ierr)

    call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, ierr, &
         access_prp = fapl_id)
    if (ierr .ne. 0) then
       call neko_error('steepest_descent: unable to open HDF5 file "' // &
            trim(filename) // '"')
    end if

    link_exists = .false.
    call h5lexists_f(file_id, 'SteepestDescent', link_exists, ierr)
    if (link_exists) then
       call h5gopen_f(file_id, 'SteepestDescent', sd_grp_id, ierr)
    else
       call h5gcreate_f(file_id, 'SteepestDescent', sd_grp_id, ierr)
    end if

    link_exists = .false.
    call h5lexists_f(sd_grp_id, 'checkpoint', link_exists, ierr)
    if (link_exists) then
       if (.not. overwrite_flag) then
          call neko_error('steepest_descent: HDF5 file "' // &
               trim(filename) // '" already contains a steepest descent ' // &
               'checkpoint; use overwrite option to replace')
       end if
       call h5ldelete_f(sd_grp_id, 'checkpoint', ierr)
    end if
    call h5gcreate_f(sd_grp_id, 'checkpoint', grp_id, ierr)
    call h5gclose_f(sd_grp_id, ierr)

    ! ------------------------------------------------------------------------ !
    ! Write the size attributes

    ddim(1) = 1
    call h5screate_f(H5S_SCALAR_F, filespace, ierr)
    call h5acreate_f(grp_id, 'n_global', H5T_NATIVE_INTEGER, filespace, &
         attr_id, ierr, h5p_default_f, h5p_default_f)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, n_global, ddim, ierr)
    call h5aclose_f(attr_id, ierr)
    call h5sclose_f(filespace, ierr)

    ddim(1) = int(pe_size, hsize_t)
    call h5screate_simple_f(1, ddim, filespace, ierr)
    call h5acreate_f(grp_id, 'n', H5T_NATIVE_INTEGER, filespace, attr_id, &
         ierr)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, n_array, ddim, ierr)
    call h5aclose_f(attr_id, ierr)
    call h5sclose_f(filespace, ierr)

    ! ------------------------------------------------------------------------ !
    ! Write the distributed datasets

    call sd_write_dist(grp_id, 'xold1', object%xold1, n_accum - n, n_global)
    call sd_write_dist(grp_id, 'xold2', object%xold2, n_accum - n, n_global)
    call sd_write_dist(grp_id, 'delta', object%delta, n_accum - n, n_global)

    ! ------------------------------------------------------------------------ !
    ! Close the group and file

    call h5gclose_f(grp_id, ierr)
    call h5fclose_f(file_id, ierr)
    call h5pclose_f(fapl_id, ierr)
    call h5close_f(ierr)

  end subroutine sd_optimizer_save_checkpoint_hdf5

  !> Load the steepest descent state from an HDF5 checkpoint.
  !! @details The global and local sizes must match those of the checkpoint.
  !! @param object The steepest descent optimizer.
  !! @param filename The HDF5 file to read from.
  module subroutine sd_optimizer_load_checkpoint_hdf5(object, filename)
    class(steepest_descent_optimizer_t), intent(inout) :: object
    character(len=*), intent(in) :: filename
    integer(hid_t) :: fapl_id, file_id, grp_id, attr_id, space_id
    integer(hsize_t) :: ddim(1), n_stored
    integer :: ierr, info, n, n_global, n_accum, n_array(pe_size)
    integer :: n_global_file, mismatch

    ! Local and global sizes of the current run
    n = object%xold1%size()
    call MPI_Allreduce(n, n_global, 1, MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
    call MPI_Scan(n, n_accum, 1, MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)

    ! ------------------------------------------------------------------------ !
    ! Open the file and the checkpoint group

    call h5open_f(ierr)
    call h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
    info = MPI_INFO_NULL%mpi_val
    call h5pset_fapl_mpio_f(fapl_id, NEKO_COMM%mpi_val, info, ierr)
    call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, ierr, &
         access_prp = fapl_id)
    call h5gopen_f(file_id, SD_GROUP, grp_id, ierr)

    if (ierr .ne. 0) then
       call neko_error('steepest_descent: unable to open HDF5 file "' // &
            trim(filename) // '" or group "' // SD_GROUP // '"')
    end if

    ! ------------------------------------------------------------------------ !
    ! Read and check the size attributes

    n_global_file = -1
    ddim(1) = 1
    call h5aopen_f(grp_id, 'n_global', attr_id, ierr)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, n_global_file, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    ! The read fills the stored extent, so it must match the rank count
    call h5aopen_f(grp_id, 'n', attr_id, ierr)
    call h5aget_space_f(attr_id, space_id, ierr)
    call h5sget_simple_extent_npoints_f(space_id, n_stored, ierr)
    call h5sclose_f(space_id, ierr)
    if (n_stored .ne. int(pe_size, hsize_t)) then
       call neko_error('steepest_descent: checkpoint was written with a ' // &
            'different number of ranks')
    end if

    n_array = -1
    ddim(1) = int(pe_size, hsize_t)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, n_array, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    if (n_global_file .ne. n_global) then
       call neko_error('steepest_descent: mismatch in n_global during ' // &
            'HDF5 read')
    end if

    ! The local sizes differ per rank, so agree on the outcome before the
    ! collective reads
    mismatch = 0
    if (n_array(pe_rank + 1) .ne. n) mismatch = 1
    call MPI_Allreduce(MPI_IN_PLACE, mismatch, 1, MPI_INTEGER, MPI_MAX, &
         NEKO_COMM, ierr)
    if (mismatch .ne. 0) then
       call neko_error('steepest_descent: mismatch in n during HDF5 read')
    end if

    ! ------------------------------------------------------------------------ !
    ! Read the distributed datasets

    call sd_read_dist(grp_id, 'xold1', object%xold1, n_accum - n, n_global)
    call sd_read_dist(grp_id, 'xold2', object%xold2, n_accum - n, n_global)
    call sd_read_dist(grp_id, 'delta', object%delta, n_accum - n, n_global)

    ! ------------------------------------------------------------------------ !
    ! Close the group and file

    call h5gclose_f(grp_id, ierr)
    call h5fclose_f(file_id, ierr)
    call h5pclose_f(fapl_id, ierr)
    call h5close_f(ierr)

    ! Ensure device state is updated
    call object%xold1%copy_from(HOST_TO_DEVICE, .true.)
    call object%xold2%copy_from(HOST_TO_DEVICE, .true.)
    call object%delta%copy_from(HOST_TO_DEVICE, .true.)

  end subroutine sd_optimizer_load_checkpoint_hdf5

  ! -------------------------------------------------------------------------- !
  ! Private helpers

  !> HDF5 type matching the working precision.
  !! @return The HDF5 native type of `real(kind=rp)`.
  function sd_h5t_real() result(h5t_real)
    integer(hid_t) :: h5t_real

    select case (rp)
    case (dp)
       h5t_real = H5T_NATIVE_DOUBLE
    case (sp)
       h5t_real = H5T_NATIVE_REAL
    case default
       call neko_error('steepest_descent: unsupported real kind for HDF5')
    end select

  end function sd_h5t_real

  !> Write a distributed vector as one global dataset, collectively.
  !! @details A rank without local values still takes part in the collective
  !! write, with an empty selection and a placeholder buffer, since the host
  !! array of an empty vector may be unallocated.
  !! @param grp_id The group to create the dataset in.
  !! @param name The name of the dataset.
  !! @param vec The local part of the vector.
  !! @param offset The global index of the first local value, from 0.
  !! @param n_global The global number of values.
  subroutine sd_write_dist(grp_id, name, vec, offset, n_global)
    integer(hid_t), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    type(vector_t), intent(in) :: vec
    integer, intent(in) :: offset, n_global
    integer(hid_t) :: filespace, memspace, xf_id, dset_id
    integer(hsize_t) :: ddim(1), dcount(1), doffset(1)
    real(kind=rp) :: empty(1)
    integer :: ierr

    empty = 0.0_rp
    ddim(1) = int(n_global, hsize_t)
    dcount(1) = int(vec%size(), hsize_t)
    doffset(1) = int(offset, hsize_t)

    call h5screate_simple_f(1, ddim, filespace, ierr)
    call h5screate_simple_f(1, dcount, memspace, ierr)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, doffset, dcount, &
         ierr)

    call h5pcreate_f(H5P_DATASET_XFER_F, xf_id, ierr)
    call h5pset_dxpl_mpio_f(xf_id, H5FD_MPIO_COLLECTIVE_F, ierr)

    call h5dcreate_f(grp_id, name, sd_h5t_real(), filespace, dset_id, ierr)
    if (vec%size() .gt. 0) then
       call h5dwrite_f(dset_id, sd_h5t_real(), vec%x, dcount, ierr, &
            file_space_id = filespace, mem_space_id = memspace, &
            xfer_prp = xf_id)
    else
       call h5dwrite_f(dset_id, sd_h5t_real(), empty, dcount, ierr, &
            file_space_id = filespace, mem_space_id = memspace, &
            xfer_prp = xf_id)
    end if
    call h5dclose_f(dset_id, ierr)

    call h5pclose_f(xf_id, ierr)
    call h5sclose_f(memspace, ierr)
    call h5sclose_f(filespace, ierr)

  end subroutine sd_write_dist

  !> Read this rank's part of a global dataset, collectively.
  !! @details A rank without local values still takes part in the collective
  !! read, with an empty selection and a placeholder buffer, since the host
  !! array of an empty vector may be unallocated.
  !! @param grp_id The group holding the dataset.
  !! @param name The name of the dataset.
  !! @param vec The local part of the vector, read into its host array.
  !! @param offset The global index of the first local value, from 0.
  !! @param n_global The global number of values.
  subroutine sd_read_dist(grp_id, name, vec, offset, n_global)
    integer(hid_t), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    type(vector_t), intent(inout) :: vec
    integer, intent(in) :: offset, n_global
    integer(hid_t) :: filespace, memspace, xf_id, dset_id
    integer(hsize_t) :: ddim(1), dcount(1), doffset(1)
    real(kind=rp) :: empty(1)
    integer :: ierr

    ddim(1) = int(n_global, hsize_t)
    dcount(1) = int(vec%size(), hsize_t)
    doffset(1) = int(offset, hsize_t)

    call h5screate_simple_f(1, ddim, filespace, ierr)
    call h5screate_simple_f(1, dcount, memspace, ierr)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, doffset, dcount, &
         ierr)

    call h5pcreate_f(H5P_DATASET_XFER_F, xf_id, ierr)
    call h5pset_dxpl_mpio_f(xf_id, H5FD_MPIO_COLLECTIVE_F, ierr)

    call h5dopen_f(grp_id, name, dset_id, ierr)
    if (vec%size() .gt. 0) then
       call h5dread_f(dset_id, sd_h5t_real(), vec%x, dcount, ierr, &
            file_space_id = filespace, mem_space_id = memspace, &
            xfer_prp = xf_id)
    else
       call h5dread_f(dset_id, sd_h5t_real(), empty, dcount, ierr, &
            file_space_id = filespace, mem_space_id = memspace, &
            xfer_prp = xf_id)
    end if
    call h5dclose_f(dset_id, ierr)

    call h5pclose_f(xf_id, ierr)
    call h5sclose_f(memspace, ierr)
    call h5sclose_f(filespace, ierr)

  end subroutine sd_read_dist

end submodule steepest_descent_hdf5_io
