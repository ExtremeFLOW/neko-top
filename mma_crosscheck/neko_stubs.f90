! Minimal stand-ins for the Neko modules that neko-top's mma.f90 / mma_cpu.f90
! depend on. Only what is needed to compile and run the CPU path.
module num_types
  implicit none
  integer, parameter :: sp = kind(1.0), dp = kind(1.0d0), rp = dp
end module num_types

module neko_config
  implicit none
  integer, parameter :: NEKO_BCKND_DEVICE = 0, NEKO_BCKND_CUDA = 0, &
       NEKO_BCKND_HIP = 0, NEKO_BCKND_OPENCL = 0
end module neko_config

module comm
  use mpi_f08
  implicit none
  integer :: pe_rank = 0, pe_size = 1
  type(MPI_Comm) :: NEKO_COMM
  type(MPI_Datatype) :: MPI_REAL_PRECISION
end module comm

module utils
  implicit none
contains
  subroutine neko_error(msg)
    character(len=*), intent(in) :: msg
    write(*,*) 'NEKO ERROR: ', trim(msg)
    error stop 1
  end subroutine neko_error
  subroutine filename_suffix(fname, suffix)
    character(len=*), intent(in) :: fname
    character(len=*), intent(out) :: suffix
    integer :: i
    i = index(fname, '.', back = .true.)
    suffix = fname(i+1:)
  end subroutine filename_suffix
end module utils

module math
  use num_types
  implicit none
  real(kind=rp), parameter :: NEKO_EPS = epsilon(1.0_rp)
end module math

module profiler
  implicit none
contains
  subroutine profiler_start_region(name, region_id)
    character(len=*), intent(in) :: name
    integer, intent(in), optional :: region_id
  end subroutine profiler_start_region
  subroutine profiler_end_region(name, region_id)
    character(len=*), intent(in), optional :: name
    integer, intent(in), optional :: region_id
  end subroutine profiler_end_region
end module profiler

module lapack_interfaces
  implicit none
  interface
     subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
       integer :: n, nrhs, lda, ldb, info
       integer :: ipiv(*)
       double precision :: a(lda, *), b(ldb, *)
     end subroutine dgesv
  end interface
end module lapack_interfaces

module device
  use num_types
  use, intrinsic :: iso_c_binding
  implicit none
  integer, parameter :: HOST_TO_DEVICE = 1, DEVICE_TO_HOST = 2
  interface device_memcpy
     module procedure device_memcpy_r1, device_memcpy_r2
  end interface device_memcpy
contains
  subroutine device_memcpy_r1(x, x_d, n, dir, sync)
    real(rp), intent(inout) :: x(:)
    type(c_ptr), intent(inout) :: x_d
    integer, intent(in) :: n, dir
    logical, intent(in), optional :: sync
  end subroutine device_memcpy_r1
  subroutine device_memcpy_r2(x, x_d, n, dir, sync)
    real(rp), intent(inout) :: x(:,:)
    type(c_ptr), intent(inout) :: x_d
    integer, intent(in) :: n, dir
    logical, intent(in), optional :: sync
  end subroutine device_memcpy_r2
end module device

module vector
  use num_types
  use, intrinsic :: iso_c_binding
  implicit none
  type, public :: vector_t
     real(rp), allocatable :: x(:)
     type(c_ptr) :: x_d = c_null_ptr
     integer :: n = 0
   contains
     procedure, pass(v) :: init => vector_init
     procedure, pass(v) :: free => vector_free
     procedure, pass(v) :: size => vector_size
     procedure, pass(v) :: copy_from => vector_copy_from
     procedure, pass(v) :: vector_assign_vector
     procedure, pass(v) :: vector_assign_scalar
     generic :: assignment(=) => vector_assign_vector, vector_assign_scalar
  end type vector_t
contains
  subroutine vector_init(v, n)
    class(vector_t), intent(inout) :: v
    integer, intent(in) :: n
    call v%free()
    allocate(v%x(n))
    v%x = 0.0_rp
    v%n = n
  end subroutine vector_init
  subroutine vector_free(v)
    class(vector_t), intent(inout) :: v
    if (allocated(v%x)) deallocate(v%x)
    v%n = 0
  end subroutine vector_free
  pure integer function vector_size(v)
    class(vector_t), intent(in) :: v
    vector_size = v%n
  end function vector_size
  subroutine vector_copy_from(v, memdir, sync)
    class(vector_t), intent(inout) :: v
    integer, intent(in) :: memdir
    logical, intent(in) :: sync
  end subroutine vector_copy_from
  subroutine vector_assign_vector(v, w)
    class(vector_t), intent(inout) :: v
    type(vector_t), intent(in) :: w
    if (.not. allocated(v%x)) call v%init(w%n)
    v%x = w%x
  end subroutine vector_assign_vector
  subroutine vector_assign_scalar(v, s)
    class(vector_t), intent(inout) :: v
    real(rp), intent(in) :: s
    v%x = s
  end subroutine vector_assign_scalar
end module vector

module matrix
  use num_types
  use, intrinsic :: iso_c_binding
  implicit none
  type, public :: matrix_t
     real(rp), allocatable :: x(:,:)
     type(c_ptr) :: x_d = c_null_ptr
   contains
     procedure, pass(mt) :: init => matrix_init
     procedure, pass(mt) :: free => matrix_free
     procedure, pass(mt) :: copy_from => matrix_copy_from
  end type matrix_t
contains
  subroutine matrix_init(mt, nr, nc)
    class(matrix_t), intent(inout) :: mt
    integer, intent(in) :: nr, nc
    call mt%free()
    allocate(mt%x(nr, nc))
    mt%x = 0.0_rp
  end subroutine matrix_init
  subroutine matrix_free(mt)
    class(matrix_t), intent(inout) :: mt
    if (allocated(mt%x)) deallocate(mt%x)
  end subroutine matrix_free
  subroutine matrix_copy_from(mt, memdir, sync)
    class(matrix_t), intent(inout) :: mt
    integer, intent(in) :: memdir
    logical, intent(in) :: sync
  end subroutine matrix_copy_from
end module matrix

module json_module
  implicit none
  type, public :: json_file
     integer :: dummy = 0
  end type json_file
end module json_module

module json_utils
  use num_types
  use json_module
  implicit none
  interface json_get_or_default
     module procedure jgod_real, jgod_int, jgod_logical, jgod_string
  end interface json_get_or_default
contains
  subroutine jgod_real(json, name, val, default)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: name
    real(rp), intent(out) :: val
    real(rp), intent(in) :: default
    val = default
  end subroutine jgod_real
  subroutine jgod_int(json, name, val, default)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: name
    integer, intent(out) :: val
    integer, intent(in) :: default
    val = default
  end subroutine jgod_int
  subroutine jgod_logical(json, name, val, default)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: name
    logical, intent(out) :: val
    logical, intent(in) :: default
    val = default
  end subroutine jgod_logical
  subroutine jgod_string(json, name, val, default)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: name
    character(len=:), allocatable, intent(out) :: val
    character(len=*), intent(in) :: default
    val = default
  end subroutine jgod_string
end module json_utils

module logger
  implicit none
  type :: log_t
   contains
     procedure, nopass :: section => log_section
     procedure, nopass :: message => log_message
     procedure, nopass :: end_section => log_end_section
  end type log_t
  type(log_t) :: neko_log
  logical :: log_verbose = .false.
contains
  subroutine log_section(name)
    character(len=*), intent(in) :: name
    if (log_verbose) write(*,'(A)') '--- ' // trim(name)
  end subroutine log_section
  subroutine log_message(msg)
    character(len=*), intent(in) :: msg
    if (log_verbose) write(*,'(A)') '    ' // trim(msg)
  end subroutine log_message
  subroutine log_end_section()
  end subroutine log_end_section
end module logger

module scratch_registry
  implicit none
  type, public :: scratch_registry_t
     integer :: dummy = 0
   contains
     procedure, pass(this) :: init => sr_init
     procedure, pass(this) :: free => sr_free
  end type scratch_registry_t
contains
  subroutine sr_init(this)
    class(scratch_registry_t), intent(inout) :: this
  end subroutine sr_init
  subroutine sr_free(this)
    class(scratch_registry_t), intent(inout) :: this
  end subroutine sr_free
end module scratch_registry
