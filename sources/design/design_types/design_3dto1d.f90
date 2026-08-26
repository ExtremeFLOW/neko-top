!> @file design_3dto1d.f90
!! @copyright
!! Copyright (c) 2025, The Neko-TOP Authors
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

! Implements the `design_3dto1d_t` type.
module design_3dto1d
  use num_types, only: rp, sp
  use json_module, only: json_file
  use mapping, only: mapping_t
  use PDE_filter, only: PDE_filter_t
  use RAMP_mapping, only: RAMP_mapping_t
  use coefs, only: coef_t
  use scratch_registry, only: neko_scratch_registry
  use point_zone_registry, only: neko_point_zone_registry
  use point_zone, only: point_zone_t
  use mask_ops, only: mask_exterior_const
  use neko_config, only: NEKO_BCKND_DEVICE
  use device, only: device_memcpy, HOST_TO_DEVICE, DEVICE_TO_HOST
  use device_math, only: device_copy
  use design, only: design_t
  use math, only: rzero
  use simulation_m, only: simulation_t
  use json_module, only: json_file
  use json_utils, only: json_get
  use utils, only: neko_error

  use vector, only: vector_t
  use math, only: copy

  use mpi_f08, only: mpi_exscan, mpi_sum, MPI_INTEGER, MPI_Allreduce
  use comm, only: pe_rank, pe_size, neko_comm, mpi_real_precision

  use fld_file_output, only: fld_file_output_t

  implicit none
  private

  !> A topology optimization design variable
  type, extends(design_t), public :: design_3dto1d_t
     private

     type(vector_t) :: values

   contains

     ! ----------------------------------------------------------------------- !
     ! Initializations

     !> Initialize the design
     generic, public :: init => init_from_components
     !> Initialize the design from components
     procedure, pass(this) :: init_from_components => &
          design_3dto1d_init_from_components

     !> Retrieve the design variables
     procedure, pass(this) :: get_values => design_3dto1d_get_values

     !> Update the design
     procedure, pass(this) :: update_design => design_3dto1d_update_design

     !> Write the design
     procedure, pass(this) :: write => design_3dto1d_write

     !> Destructor
     procedure, pass(this) :: free => design_3dto1d_free

     !> map (this will include everything from mapping
     procedure, pass(this) :: map_forward => design_3dto1d_map_forward
     !> this will contain chain rule for going backwards
     procedure, pass(this) :: map_backward => design_3dto1d_map_backward

  end type design_3dto1d_t

contains


  subroutine design_3dto1d_init_from_components(this, n)
    class(design_3dto1d_t), intent(inout) :: this
    integer, intent(in) :: n

    call this%init_base('design_3dto1d', n)

    call this%values%init(n)

  end subroutine design_3dto1d_init_from_components

  !> Free the design
  subroutine design_3dto1d_free(this)
    class(design_3dto1d_t), intent(inout) :: this

    call this%free_base()
    call this%values%free()
  end subroutine design_3dto1d_free

  subroutine design_3dto1d_map_forward(this)
    class(design_3dto1d_t), intent(inout) :: this

  end subroutine design_3dto1d_map_forward

  subroutine design_3dto1d_map_backward(this, sensitivity)
    class(design_3dto1d_t), intent(inout) :: this
    type(vector_t), intent(in) :: sensitivity
  end subroutine design_3dto1d_map_backward


  subroutine design_3dto1d_get_values(this, values)
    class(design_3dto1d_t), intent(in) :: this
    type(vector_t), intent(inout) :: values

    if (this%size() .ne. values%size()) then
       call neko_error('Get design: size mismatch')
    end if

    values = this%values

  end subroutine design_3dto1d_get_values

  subroutine design_3dto1d_update_design(this, values)
    class(design_3dto1d_t), intent(inout) :: this
    type(vector_t), intent(inout) :: values
    integer :: n

    n = this%size()

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(this%values%x_d, values%x_d, n)
    else
       this%values = values
    end if


  end subroutine design_3dto1d_update_design

  subroutine design_3dto1d_write(this, idx)
    class(design_3dto1d_t), intent(inout) :: this
    integer, intent(in) :: idx

    character(len=100) :: filename
    integer :: i, iunit, ierr
    real(rp) :: Le
    real(rp), allocatable :: global_values(:)
    real(rp), allocatable :: global_x(:)
    integer :: global_size, local_size
    real(rp) :: L_total
    integer, allocatable :: recvcounts(:), displs(:)
    integer :: root_rank = 0

    L_total = 2.0_rp

    ! Get local size on all ranks
    local_size = this%size()

    ! Only rank 0 handles global arrays
    if (pe_rank == root_rank) then
       ! Get global size
       global_size = this%size_global()
       allocate(global_values(global_size))
       allocate(global_x(global_size))
       allocate(recvcounts(pe_size), displs(pe_size))
       ! Calculate element length and x positions
       Le = L_total / real(global_size, kind=rp)
       do i = 1, global_size
          global_x(i) = Le * (i - 0.5_rp) ! Center of each element
       end do
    else
       ! Non-root ranks: minimal dummy allocations
       allocate(global_values(1), recvcounts(1), displs(1))
    endif

    ! First, gather all the local sizes to rank 0
    call MPI_Gather(local_size, 1, MPI_INTEGER, &
         recvcounts, 1, MPI_INTEGER, &
         root_rank, neko_comm, ierr)

    ! Calculate displacements on rank 0
    if (pe_rank == root_rank) then
       displs(1) = 0
       do i = 2, pe_size
          displs(i) = displs(i-1) + recvcounts(i-1)
       end do
    endif

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%values%x, this%values%x_d, local_size, &
            DEVICE_TO_HOST, sync = .false.)
    end if
    ! Now gather the actual data with proper displacement handling
    call MPI_Gatherv(this%values%x, local_size, mpi_real_precision, &
         global_values, recvcounts, displs, mpi_real_precision, &
         root_rank, neko_comm, ierr)

    if (pe_rank == root_rank) then
       ! Create filename with iteration index
       write(filename, '(A,I0.6,A)') 'design_iter_', idx, '.txt'
       ! Open file for writing
       open(newunit=iunit, file=trim(filename), status='replace', action='write')
       ! Write header
       write(iunit, '(A,I0)') '# Iteration: ', idx
       write(iunit, '(A)') '# x_position height'

       ! Write data
       do i = 1, global_size
          write(iunit, '(2E16.8)') global_x(i), global_values(i)
       end do

       close(iunit)

       deallocate(global_values, global_x, recvcounts, displs)
       call nekotop_log%message("Design written to " // trim(filename))
    else
       deallocate(global_values, recvcounts, displs)
    endif
  end subroutine design_3dto1d_write

end module design_3dto1d
