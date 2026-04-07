!> @file optimizer_hdf5_checkpoint.f90
!! @brief HDF5 IO submodule for the optimizer object.
!! @details
!! This submodule provides routines for saving and loading the optimizer
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

!> Submodule for handling HDF5 IO for the optimizer object.
submodule (optimizer) optimizer_hdf5_io
  use hdf5
  use mpi_f08, only: MPI_INFO_NULL
  use comm, only: NEKO_COMM

contains

  ! -------------------------------------------------------------------------- !
  ! Checkpointing methods for the optimizer optimizer

  subroutine optimizer_save_checkpoint_hdf5(object, filename, iter, overwrite)
    class(optimizer_t), intent(inout) :: object
    character(len=*), intent(in) :: filename
    integer, intent(in) :: iter
    logical, intent(in), optional :: overwrite
    logical :: overwrite_flag, file_exists

    ! HDF5 variables
    integer(hid_t) :: file_id, fapl_id, filespace, str_type
    integer(hid_t) :: grp_id, optimizer_grp_id, attr_id
    integer :: ierr, info
    integer(hsize_t) :: ddim(1)
    logical :: optimizer_exists

    overwrite_flag = .false.
    if (present(overwrite)) overwrite_flag = overwrite

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

       ! Create the Optimizer/checkpoint group
       call h5gcreate_f(file_id, "Optimizer", optimizer_grp_id, ierr)
       call h5gcreate_f(optimizer_grp_id, "checkpoint", grp_id, ierr)

    else
       call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, ierr, fapl_id)

       ! Check for existing Optimizer group
       optimizer_exists = .false.
       call h5lexists_f(file_id, 'Optimizer', optimizer_exists, ierr)

       if (.not. optimizer_exists) then
          ! Create the Optimizer/checkpoint group
          call h5gcreate_f(file_id, "Optimizer", optimizer_grp_id, ierr)
          call h5gcreate_f(optimizer_grp_id, "checkpoint", grp_id, ierr)
       else
          ! Open the existing "Optimizer" group
          call h5gopen_f(file_id, 'Optimizer', optimizer_grp_id, ierr)
          call h5lexists_f(optimizer_grp_id, 'checkpoint', optimizer_exists, ierr)

          if (optimizer_exists .and. overwrite_flag) then
             call h5gunlink_f(optimizer_grp_id, "checkpoint", ierr)
             call h5gcreate_f(optimizer_grp_id, "checkpoint", grp_id, ierr)
          else
             call neko_error('optimizer: HDF5 file "' // trim(filename) // '" ' &
                  // 'already contains Optimizer group; ' &
                  // 'use overwrite option to replace')
          end if
       end if
    end if
    call h5gclose_f(optimizer_grp_id, ierr)

    ! ------------------------------------------------------------------------ !
    ! Write the optimizer optimizer checkpoint group

    ! Create the dataspace for attributes
    call h5screate_f(H5S_SCALAR_F, filespace, ierr)

    ! Save the optimizer type
    call h5tcopy_f(H5T_FORTRAN_S1, str_type, ierr)
    call h5tset_strpad_f(str_type, H5T_STR_SPACEPAD_F, ierr)

    ddim(1) = len_trim(object%optimizer_type)
    call h5tset_size_f(str_type, ddim(1), ierr)
    call h5acreate_f(grp_id, 'type', str_type, filespace, attr_id, &
         ierr, h5p_default_f, h5p_default_f)
    call h5awrite_f(attr_id, str_type, object%optimizer_type, ddim, ierr)
    call h5aclose_f(attr_id, ierr)
    call h5tclose_f(str_type, ierr)

    ! Save the current iteration number
    ddim = 1
    call h5acreate_f(grp_id, 'iter', H5T_NATIVE_INTEGER, filespace, attr_id, &
         ierr, h5p_default_f, h5p_default_f)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, iter, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    ! ------------------------------------------------------------------------ !
    ! Close HDF5 objects
    call h5sclose_f(filespace, ierr)
    call h5gclose_f(grp_id, ierr)
    call h5fclose_f(file_id, ierr)
    call h5pclose_f(fapl_id, ierr)
    call h5close_f(ierr)

  end subroutine optimizer_save_checkpoint_hdf5

  subroutine optimizer_load_checkpoint_hdf5(object, filename, iter)
    class(optimizer_t), intent(inout) :: object
    character(len=*), intent(in) :: filename
    integer, intent(out) :: iter

    ! HDF5 variables
    character(len=64) :: type_name
    integer(hid_t) :: file_id, fapl_id, grp_id, attr_id, str_type
    integer :: ierr, info
    integer(hsize_t) :: ddim(1)
    character(len=*), parameter :: h5_group = 'Optimizer/checkpoint'

    ! Initialize reader variables
    type_name = ''
    iter = -1

    ! ------------------------------------------------------------------------ !
    ! Open the HDF5 context and read identification parameters

    call h5open_f(ierr)

    ! Prepare the HDF5 settings for MPIO access
    call h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
    info = MPI_INFO_NULL%mpi_val
    call h5pset_fapl_mpio_f(fapl_id, NEKO_COMM%mpi_val, info, ierr)

    ! Open the HDF5 file
    call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, ierr, &
         access_prp = fapl_id)

    ! Open the optimizer optimizer checkpoint group
    call h5gopen_f(file_id, trim(h5_group), grp_id, ierr)

    if (ierr .ne. 0) then
       call neko_error('optimizer: unable to open HDF5 file "' // &
            trim(filename) // '"')
    end if

    ! Read the optimizer type and verify
    call h5aopen_f(grp_id, 'type', attr_id, ierr)
    call h5aget_type_f(attr_id, str_type, ierr)
    call h5aread_f(attr_id, str_type, type_name, ddim, ierr)
    call h5tclose_f(str_type, ierr)
    call h5aclose_f(attr_id, ierr)

    ! Verify that the type is optimizer
    if (trim(type_name) .ne. trim(object%optimizer_type)) then
       call neko_error("Incompatible optimizer type in checkpoint file")
    end if

    ! Read the current optimizer iteration
    ddim(1) = 1
    call h5aopen_f(grp_id, 'iter', attr_id, ierr)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, iter, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    if (iter .lt. 0) then
       call neko_error("Incompatible optimizer iteration in checkpoint file")
    end if

    ! Read the current iteration
    call h5gclose_f(grp_id, ierr)
    call h5fclose_f(file_id, ierr)
    call h5pclose_f(fapl_id, ierr)
    call h5close_f(ierr)

  end subroutine optimizer_load_checkpoint_hdf5

end submodule optimizer_hdf5_io
