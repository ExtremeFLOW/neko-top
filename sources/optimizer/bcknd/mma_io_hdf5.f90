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
    integer(hid_t) :: fapl_id, xfer_plist_id, file_id, dset_id, filespace, memspace, attr_id,&
         grp_id, mma_grp_id
    integer(hid_t) :: H5T_NEKO_REAL
    integer(hsize_t), dimension(1) :: ddim, dcount, doffset
    integer :: ierr, info, drank
    integer :: local_n
    integer(kind=8) :: local_n8, prefix8, total_n8
    integer(hsize_t), dimension(1) :: qdims, qmaxdims


    ! Ensure device state is on host
    call this%sync_host()
    call neko_log%message('mma: device memory synced to host for HDF5 write: ' &
         // trim(filename))

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

    ! Dataset transfer property list: collective
    call h5pcreate_f(H5P_DATASET_XFER_F, xfer_plist_id, ierr)
    call h5pset_dxpl_mpio_f(xfer_plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)

    ! Write basic scalars as attributes
    call h5gcreate_f(mma_grp_id, "Parameters", grp_id, ierr, lcpl_id=h5p_default_f, &
         gcpl_id=h5p_default_f, gapl_id=h5p_default_f)

    call h5screate_f(H5S_SCALAR_F, filespace, ierr)
    ddim = 1

    call h5acreate_f(grp_id, 'n', H5T_NATIVE_INTEGER, filespace, attr_id, &
         ierr, h5p_default_f, h5p_default_f)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, this%n, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5acreate_f(grp_id, 'm', H5T_NATIVE_INTEGER, filespace, attr_id, &
         ierr, h5p_default_f, h5p_default_f)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, this%m, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5acreate_f(grp_id, 'a0', H5T_NEKO_REAL, filespace, attr_id, &
         ierr, h5p_default_f, h5p_default_f)
    call h5awrite_f(attr_id, H5T_NEKO_REAL, this%a0, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5acreate_f(grp_id, 'f0val', H5T_NEKO_REAL, filespace, attr_id, &
         ierr, h5p_default_f, h5p_default_f)
    call h5awrite_f(attr_id, H5T_NEKO_REAL, this%f0val, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5acreate_f(grp_id, 'z', H5T_NEKO_REAL, filespace, attr_id, &
         ierr, h5p_default_f, h5p_default_f)
    call h5awrite_f(attr_id, H5T_NEKO_REAL, this%z, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5acreate_f(grp_id, 'zeta', H5T_NEKO_REAL, filespace, attr_id, &
         ierr, h5p_default_f, h5p_default_f)
    call h5awrite_f(attr_id, H5T_NEKO_REAL, this%zeta, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    ! Close the filespace and Parameters group
    call h5sclose_f(filespace, ierr)
    call h5gclose_f(grp_id, ierr)

    ! ........................................................................ !
    ! Reference from Neko's hdf5_io module

#ifdef NOTHING
    ! Create group fro the fields
    call h5gcreate_f(file_id, "Fields", grp_id, ierr, lcpl_id=h5p_default_f, &
         gcpl_id=h5p_default_f, gapl_id=h5p_default_f)

    ! Prepare some parameters for the writing process
    dcount(1) = int(dof%size(), 8)
    doffset(1) = int(msh%offset_el, 8) * int((dof%Xh%lx**3),8)
    ddim = int(dof%size(), 8)
    drank = 1
    call MPI_Allreduce(MPI_IN_PLACE, ddim(1), 1, &
         MPI_INTEGER8, MPI_SUM, NEKO_COMM, ierr)

    ! Create a file and memory space
    call h5screate_simple_f(drank, ddim, filespace, ierr)
    call h5screate_simple_f(drank, dcount, memspace, ierr)

    ! Write field at id: i
    call h5dcreate_f(grp_id, fp(i)%ptr%name, H5T_NEKO_REAL, &
         filespace, dset_id, ierr)
    call h5dget_space_f(dset_id, filespace, ierr)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
         doffset, dcount, ierr)
    call h5dwrite_f(dset_id, H5T_NEKO_REAL, &
         fp(i)%ptr%x(1,1,1,1), &
         ddim, ierr, file_space_id = filespace, &
         mem_space_id = memspace, xfer_prp = plist_id)
    call h5dclose_f(dset_id, ierr)

    ! Close the spaces used
    call h5sclose_f(filespace, ierr)
    call h5sclose_f(memspace, ierr)
    call h5gclose_f(grp_id, ierr)
#endif

    ! End of reference
    ! ........................................................................ !

    ! Helper to write a 1D dataset with rank-0 data only (safe parallel)
    drank = 1

    ! list of vectors to write (name, size, data) - per-rank hyperslabs

    ! Write X_old to the file
    local_n = this%xold1%size()
    local_n8 = int(local_n, 8)
    call MPI_Allreduce(local_n8, total_n8, 1, MPI_INTEGER8, MPI_SUM, NEKO_COMM, ierr)
    call MPI_Scan(local_n8, prefix8, 1, MPI_INTEGER8, MPI_SUM, NEKO_COMM, ierr)

    dcount(1) = int(this%xold1%size(), kind = 8)
    doffset(1) = prefix8 - local_n8
    ddim(1) = total_n8

    do i = 0, pe_size
       if (i .eq. pe_rank) then
          write(*,*) '----------------------------------------'
          write(*,*) 'Writing xold1:', pe_rank
          write(*,*) '  local size=', local_n, ' total size=', total_n8
          write(*,*) '  offset=', doffset(1), ' count=', dcount(1)
          write(*,*) ' Data:'
          write(*,*) pe_rank, ": ", this%xold1%x
       end if
       call MPI_Barrier(NEKO_COMM, ierr)
    end do

    if (pe_rank .eq. 0) then
       write(*,*) '----------------------------------------'
    end if

    ! Create a file and memory space
    call h5screate_simple_f(drank, ddim, filespace, ierr)
    if (ierr .ne. 0) call neko_error('mma: h5screate_simple_f(filespace) failed')
    call h5screate_simple_f(drank, dcount, memspace, ierr)
    if (ierr .ne. 0) call neko_error('mma: h5screate_simple_f(memspace) failed')

    ! Save the array type objects
    call h5dcreate_f(mma_grp_id, 'xold1', H5T_NEKO_REAL, filespace, dset_id, ierr)
    if (ierr .ne. 0) call neko_error('mma: h5dcreate_f(xold1) failed')
    call h5dget_space_f(dset_id, filespace, ierr)
    if (ierr .ne. 0) call neko_error('mma: h5dget_space_f failed')

    ! Diagnostic: query the dataspace extents
    call h5sget_simple_extent_dims_f(filespace, qdims, qmaxdims, ierr)
    if (pe_rank .eq. 0) then
       write(*,*) 'mma: dataset filespace dims=', qdims(1), ' max=', qmaxdims(1)
    end if
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, doffset, dcount, ierr)
    if (ierr .ne. 0) call neko_error('mma: h5sselect_hyperslab_f failed')

    call h5dwrite_f(dset_id, H5T_NEKO_REAL, this%xold1%x(1), ddim, ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = xfer_plist_id)
    if (ierr .ne. 0) then
       write(*,*) 'mma: h5dwrite_f failed with ierr=', ierr
       call neko_error('mma: h5dwrite_f failed')
    end if

    call h5dclose_f(dset_id, ierr)
    if (ierr .ne. 0) call neko_error('mma: h5dclose_f failed')

    ! Close the spaces used
    call h5sclose_f(filespace, ierr)
    if (ierr .ne. 0) call neko_error('mma: h5sclose_f(filespace) failed')
    call h5sclose_f(memspace, ierr)
    if (ierr .ne. 0) call neko_error('mma: h5sclose_f(memspace) failed')

    ! Ensure data is flushed to disk
    call h5fflush_f(file_id, H5F_SCOPE_GLOBAL_F, ierr)
    if (ierr .ne. 0) call neko_error('mma: h5fflush_f failed')

    ! Close the group and file
    call h5gclose_f(mma_grp_id, ierr)
    call h5pclose_f(xfer_plist_id, ierr)
    call h5pclose_f(fapl_id, ierr)
    call h5fclose_f(file_id, ierr)
    call h5close_f(ierr)

  end subroutine mma_write_hdf5


#else
  module subroutine mma_write_hdf5(this, filename)
    class(mma_t), intent(inout) :: this
    character(len=*), intent(in) :: filename
    call neko_error('mma: HDF5 support not enabled rebuild with HAVE_HDF5')
  end subroutine mma_write_hdf5
#endif

end submodule mma_io_hdf5
