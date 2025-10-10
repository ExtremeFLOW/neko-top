! Copyright (c) 2024, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.

! Implements the `simple_design_t` type.
module simple_design
  use num_types, only: rp, sp
  use field, only: field_t
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
  use design, only: design_t
  use math, only: rzero
  use simulation_m, only: simulation_t
  use comm, only: pe_size
  use json_module, only: json_file
  use json_utils, only: json_get, json_get_or_default
  use simple_brinkman_source_term, only: simple_brinkman_source_term_t
  use vector, only: vector_t
  use math, only: copy
  use field_registry, only: neko_field_registry
  use mpi_f08, only: MPI_Comm_size, MPI_COMM_WORLD
  use utils, only: neko_error
  implicit none
  private

  !> A topology optimization design variable
  type, extends(design_t), public :: simple_design_t
     private

     type(vector_t) :: values
     type(vector_t) :: x_coord
     type(vector_t) :: y_coord
     type(vector_t) :: z_coord

     ! needed to write the design into a vtk file with connectivity
     integer :: nx, ny, nz

   contains

     ! ----------------------------------------------------------------------- !
     ! Initializations

     !> Initialize the design
     generic, public :: init => init_from_json, init_from_components
     !> Initialize the design from a JSON file
     procedure, pass(this) :: init_from_json => &
          design_simple_init_from_json
     !> Initialize the design from components
     procedure, pass(this) :: init_from_components => &
          design_simple_init_from_components

     !> Add mappings to the design
     procedure, pass(this) :: add_mapping => design_simple_add_mapping

     !> Retrieve the design variables
     procedure, pass(this) :: get_values => design_simple_get_values
     !> Retrieve the x location of the design variables
     procedure, pass(this) :: design_get_x => design_simple_get_x
     !> Retrieve the y location of the design variables
     procedure, pass(this) :: design_get_y => design_simple_get_y
     !> Retrieve the z location of the design variables
     procedure, pass(this) :: design_get_z => design_simple_get_z

     !> Update the design
     procedure, pass(this) :: update_design => design_simple_update_design

     !> map (this will include everything from mapping
     procedure, pass(this) :: map_forward => design_simple_map_forward

     !> this will contain chain rule for going backwards
     procedure, pass(this) :: map_backward => design_simple_map_backward

     !> Write the design
     procedure, pass(this) :: write => design_simple_write

     !> Destructor
     procedure, pass(this) :: free => design_simple_free

  end type simple_design_t

contains

  !> Initialize the design from a JSON file
  subroutine design_simple_init_from_json(this, parameters)
    class(simple_design_t), intent(inout) :: this
    type(json_file), intent(inout) :: parameters
    character(len=:), allocatable :: type, name
    integer :: n, nx, ny, nz, i, j, k, index
    real(kind=rp), dimension(:), allocatable :: limits
    type(vector_t) :: x, y, z

    call json_get(parameters, 'domain.type', type)
    call json_get_or_default(parameters, 'name', name, 'Simple Design')

    select case (trim(type))
    case ("box")
       call json_get(parameters, 'optimization.design.domain.nx', nx)
       call json_get(parameters, 'optimization.design.domain.ny', ny)
       call json_get(parameters, 'optimization.design.domain.nz', nz)
       call json_get(parameters, 'optimization.design.domain.limits', limits)
       n = (nx + 1) * (ny + 1) * (nz + 1)
       this%nx = nx
       this%ny = ny
       this%nz = nz


       call x%init(n)
       call y%init(n)
       call z%init(n)
       index = 1
       do k = 1, nz + 1
          do j = 1, ny + 1
             do i = 1, nx + 1
                x%x(index) = limits(1) + (limits(2) - limits(1)) * &
                     real(i - 1, kind=rp) / real(nx, kind=rp)
                y%x(index) = limits(3) + (limits(4) - limits(3)) * &
                     real(j - 1, kind=rp) / real(ny, kind=rp)
                z%x(index) = limits(5) + (limits(6) - limits(5)) * &
                     real(k - 1, kind=rp) / real(nz, kind=rp)
                index = index + 1
             end do
          end do
       end do

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_memcpy(x%x, x%x_d, n, HOST_TO_DEVICE, sync = .false.)
          call device_memcpy(y%x, y%x_d, n, HOST_TO_DEVICE, sync = .false.)
          call device_memcpy(z%x, z%x_d, n, HOST_TO_DEVICE, sync = .true.)
       end if

    end select

    call this%init_from_components(name, n, x, y, z)

  end subroutine design_simple_init_from_json

  subroutine design_simple_init_from_components(this, name, n, x, y, z)
    class(simple_design_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    type(vector_t), intent(in) :: x, y, z
    integer :: nproc, ierr

    ! Throw an error if we run on multiple MPI ranks
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
    if (nproc .gt. 1) then
       call neko_error('simple_design_t can only be used on a single MPI rank')
    end if

    if (pe_size .ne. 1) then
       call neko_error("Simple design can only be used with a single MPI " // &
            "process.")
    end if

    call this%init_base(name, n)

    call this%values%init(n)
    this%x_coord = x
    this%y_coord = y
    this%z_coord = z

  end subroutine design_simple_init_from_components

  !> Free the design
  subroutine design_simple_free(this)
    class(simple_design_t), intent(inout) :: this

    call this%free_base()
    call this%values%free()
    call this%x_coord%free()
    call this%y_coord%free()
    call this%z_coord%free()

  end subroutine design_simple_free

  !> Add mappings to the design
  subroutine design_simple_add_mapping(this, parameters, simulation)
    class(simple_design_t), intent(inout) :: this
    type(json_file), intent(inout) :: parameters
    type(simulation_t), intent(inout) :: simulation

  end subroutine design_simple_add_mapping


  subroutine design_simple_map_forward(this)
    class(simple_design_t), intent(inout) :: this


  end subroutine design_simple_map_forward

  subroutine design_simple_get_values(this, values)
    class(simple_design_t), intent(in) :: this
    type(vector_t), intent(inout) :: values

    values = this%values

  end subroutine design_simple_get_values

  subroutine design_simple_get_x(this, x)
    class(simple_design_t), intent(in) :: this
    type(vector_t), intent(inout) :: x

    x = this%x_coord

  end subroutine design_simple_get_x

  subroutine design_simple_get_y(this, y)
    class(simple_design_t), intent(in) :: this
    type(vector_t), intent(inout) :: y

    y = this%y_coord

  end subroutine design_simple_get_y

  subroutine design_simple_get_z(this, z)
    class(simple_design_t), intent(in) :: this
    type(vector_t), intent(inout) :: z

    z = this%z_coord

  end subroutine design_simple_get_z

  subroutine design_simple_update_design(this, values)
    class(simple_design_t), intent(inout) :: this
    type(vector_t), intent(inout) :: values

    this%values = values

  end subroutine design_simple_update_design

  subroutine design_simple_map_backward(this, sensitivity)
    class(simple_design_t), intent(inout) :: this
    type(vector_t), intent(in) :: sensitivity

  end subroutine design_simple_map_backward

  subroutine design_simple_write(this, idx)
    class(simple_design_t), intent(inout) :: this
    integer, intent(in) :: idx

    integer :: i
    integer :: npts, nx, ny, nz
    character(len=100) :: filename
    integer :: ios, funit

    nx = this%nx + 1
    ny = this%ny + 1
    nz = this%nz + 1
    npts = nx * ny * nz

    ! Synchronize the device memory if using a GPU backend is enabled
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%x_coord%x, this%x_coord%x_d, npts, &
            DEVICE_TO_HOST, sync = .false.)
       call device_memcpy(this%y_coord%x, this%y_coord%x_d, npts, &
            DEVICE_TO_HOST, sync = .false.)
       call device_memcpy(this%z_coord%x, this%z_coord%x_d, npts, &
            DEVICE_TO_HOST, sync = .false.)
       call device_memcpy(this%values%x, this%values%x_d, npts, &
            DEVICE_TO_HOST, sync = .true.)
    end if

    write(filename, '(A,I4.4,A)') 'design_', idx, '.vtk'

    open(newunit=funit, file=filename, status='replace', action='write', iostat=ios)
    if (ios .ne. 0) then
       print *, 'Error opening file ', filename
       stop
    end if

    ! Header
    write(funit, '(A)') '# vtk DataFile Version 3.0'
    write(funit, '(A)') 'Simple Scalar Field'
    write(funit, '(A)') 'ASCII'
    write(funit, '(A)') 'DATASET STRUCTURED_GRID'
    write(funit, '(A,3(1X,I0))') 'DIMENSIONS', nx, ny, nz

    ! Points
    write(funit, '(A,1X,I0,1X,A)') 'POINTS', npts, 'float'
    do i = 1, npts
       write(funit, '(3(F20.12,1X))') this%x_coord%x(i), this%y_coord%x(i), &
            this%z_coord%x(i)
    end do

    ! Scalars
    write(funit, '(A,1X,I0)') 'POINT_DATA', npts
    write(funit, '(A)') 'SCALARS density float 1'
    write(funit, '(A)') 'LOOKUP_TABLE default'
    do i = 1, npts
       write(funit, '(F20.12)') this%values%x(i)
    end do

    close(funit)
  end subroutine design_simple_write

end module simple_design
