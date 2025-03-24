module json_utils_ext
  use mpi_f08, only: MPI_Comm_rank, MPI_Initialized, MPI_Bcast, &
       MPI_COMM_WORLD, MPI_INTEGER, MPI_CHARACTER
  use json_file_module, only: json_file
  use json_value_module, only: json_value
  use utils, only: neko_error, filename_suffix

  implicit none
  private

  public :: json_key_fallback, json_read_file

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

    if (json%valid_path(lookup)) then
       string = lookup
    else if (json%valid_path(fallback)) then
       string = fallback
    else
       string = lookup
    end if

  end function json_key_fallback

  !> Read a json file taking mpi into account.
  !! This function reads a json file and broadcasts it to all ranks.
  !! @params[in] filename The name of the file to read.
  !! @return The json object.
  function json_read_file(filename) result(json)
    character(len=*), intent(in) :: filename
    type(json_file) :: json

    logical :: mpi_is_initialized
    integer :: rank, ierr, length
    character(len=:), allocatable :: json_buffer
    character(len=4) :: suffix

    ! Check if the file is a json or Neko case file.
    call filename_suffix(filename, suffix)

    if (trim(suffix) .ne. 'json' .and. trim(suffix) .ne. 'case' ) then
       call neko_error('Invalid case file')
    end if

    ! Initialize the mpi variables.
    rank = 0
    mpi_is_initialized = .false.

    ! Check if MPI is initialized and get the rank if it is.
    call MPI_Initialized(mpi_is_initialized, ierr)
    if (mpi_is_initialized) call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    ! Read the json file and broadcast it to all ranks.
    if (rank .eq. 0) then
       call json%load_file(filename = trim(filename))
    end if

    ! Serialize the json object to a string so it can be broadcast.
    if (mpi_is_initialized) then
       if (rank .eq. 0) call json%print_to_string(json_buffer)

       length = len(json_buffer)
       call MPI_Bcast(length, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

       if (rank .ne. 0) allocate(character(len=length) :: json_buffer)
       call MPI_Bcast(json_buffer, length, MPI_CHARACTER, 0, MPI_COMM_WORLD, &
            ierr)

       if (rank .ne. 0) then
          call json%load_from_string(json_buffer)
          deallocate(json_buffer)
       end if
    end if

  end function json_read_file

end module json_utils_ext
