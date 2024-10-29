module json_utils_ext

  use json_file_module, only: json_file
  use json_value_module, only: json_value
  use utils, only: neko_error, filename_suffix

  use mpi_f08, only: MPI_COMM_WORLD, MPI_INTEGER, MPI_CHARACTER, MPI_Bcast, MPI_Comm_rank
  implicit none
  private

  public :: json_key_fallback, json_get_subdict, read_case

contains

  !> Create a json_string based on fallback logic.
  !! If the lookup key is present in the json object, return it.
  !! If the fallback key is present in the json object, return it.
  !! Otherwise, return the lookup key.
  function json_key_fallback(json, lookup, fallback) result(string)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: lookup
    character(len=*), intent(in) :: fallback
    character(len=:), allocatable :: string

    if ((lookup .in. json)) then
       string = lookup
    else if (fallback .in. json) then
       string = fallback
    else
       string = lookup
    end if

  end function json_key_fallback

  !> Extract a sub-object from a json object.
  subroutine json_get_subdict(json, key, output)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: key
    type(json_file), intent(out) :: output

    type(json_value), pointer :: child
    logical :: valid

    nullify(child)

    valid = .false.
    call json%get(key, child, valid)
    if (.not. valid) then
       call neko_error('Parameter "' // &
            trim(key) // '" missing from the case file')
    end if

    call output%initialize()
    call output%add(child)
    nullify(child)

  end subroutine json_get_subdict

  function read_case(filename) result(case_params)
    character(len=*), intent(in) :: filename
    type(json_file) :: case_params

    integer :: pe_rank, ierr, length, argc
    character(len=:), allocatable :: json_buffer
    character(len=4) :: suffix
    logical :: mpi_is_initialized

    pe_rank = 0
    call MPI_Comm_rank(MPI_COMM_WORLD, pe_rank, ierr)

    if (pe_rank .eq. 0) then
       call case_params%load_file(filename = trim(filename))
       call case_params%print_to_string(json_buffer)
       length = len(json_buffer)
    end if


    call MPI_Bcast(length, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (pe_rank .ne. 0) allocate(character(len = length) :: json_buffer)
    call MPI_Bcast(json_buffer, length, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    call case_params%load_from_string(json_buffer)
    deallocate(json_buffer)

  end function read_case


end module json_utils_ext
