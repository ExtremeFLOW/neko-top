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
  use design, only: design_t
  use math, only: rzero
  use simulation_m, only: simulation_t
  use json_module, only: json_file
  use json_utils, only: json_get
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
     type(vector_t) :: x
     type(vector_t) :: y
     type(vector_t) :: z

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
     procedure, pass(this) :: get_x => design_simple_get_x
     !> Retrieve the y location of the design variables
     procedure, pass(this) :: get_y => design_simple_get_y
     !> Retrieve the z location of the design variables
     procedure, pass(this) :: get_z => design_simple_get_z

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
    character(len=:), allocatable :: type
    integer :: n, nx, ny, nz, i, j, k, index
    real(kind=rp), dimension(:), allocatable :: limits
    type(vector_t) :: x, y, z

    call json_get(parameters, 'optimization.design.domain.type', type)

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

    call this%init_from_components(n, x, y, z)

  end subroutine design_simple_init_from_json

  subroutine design_simple_init_from_components(this, n, x, y, z)
    class(simple_design_t), intent(inout) :: this
    integer, intent(in) :: n
    type(vector_t), intent(in) :: x, y, z
    integer :: nproc, ierr

    ! Throw an error if we run on multiple MPI ranks
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
    if (nproc .gt. 1) then
       call neko_error('simple_design_t can only be used on a single MPI rank')
    end if

    call this%init_base(n)

    call this%values%init(n)
    call this%x%init(n)
    call this%y%init(n)
    call this%z%init(n)

    this%x = x
    this%y = y
    this%z = z

  end subroutine design_simple_init_from_components

  !> Free the design
  subroutine design_simple_free(this)
    class(simple_design_t), intent(inout) :: this

    call this%free_base()
    call this%values%free()
    call this%x%free()
    call this%y%free()
    call this%z%free()

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

  function design_simple_get_values(this) result(values)
    class(simple_design_t), intent(in) :: this
    type(vector_t) :: values

    values = this%values

  end function design_simple_get_values

  function design_simple_get_x(this) result(x)
    class(simple_design_t), intent(in) :: this
    type(vector_t) :: x

    x = this%x

  end function design_simple_get_x

  function design_simple_get_y(this) result(y)
    class(simple_design_t), intent(in) :: this
    type(vector_t) :: y

    y = this%y

  end function design_simple_get_y

  function design_simple_get_z(this) result(z)
    class(simple_design_t), intent(in) :: this
    type(vector_t) :: z

    z = this%z

  end function design_simple_get_z

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
    integer :: ios

    nx = this%nx + 1
    ny = this%ny + 1
    nz = this%nz + 1
    npts = nx * ny * nz

    ! Synchronize the device memory if using a GPU backend is enabled
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%x%x, this%x%x_d, npts, &
            DEVICE_TO_HOST, sync = .false.)
       call device_memcpy(this%y%x, this%y%x_d, npts, &
            DEVICE_TO_HOST, sync = .false.)
       call device_memcpy(this%z%x, this%z%x_d, npts, &
            DEVICE_TO_HOST, sync = .false.)
       call device_memcpy(this%values%x, this%values%x_d, npts, &
            DEVICE_TO_HOST, sync = .true.)
    end if

    write(filename, '(A,I4.4,A)') 'design_', idx, '.vtk'

    open(unit=10, file=filename, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
       print *, 'Error opening file ', filename
       stop
    end if

    ! Header
    write(10,'(A)') '# vtk DataFile Version 3.0'
    write(10,'(A)') 'Simple Scalar Field'
    write(10,'(A)') 'ASCII'
    write(10,'(A)') 'DATASET STRUCTURED_GRID'
    write(10,'(A,3(1X,I0))') 'DIMENSIONS', nx, ny, nz

    ! Points
    write(10,'(A,1X,I0,1X,A)') 'POINTS', npts, 'float'
    do i = 1, npts
       write(10,'(3(F20.12,1X))') this%x%x(i), this%y%x(i), this%z%x(i)
    end do

    ! Scalars
    write(10,'(A,1X,I0)') 'POINT_DATA', npts
    write(10,'(A)') 'SCALARS density float 1'
    write(10,'(A)') 'LOOKUP_TABLE default'
    do i = 1, npts
       write(10,'(F20.12)') this%values%x(i)
    end do

    close(10)
  end subroutine design_simple_write

end module simple_design
