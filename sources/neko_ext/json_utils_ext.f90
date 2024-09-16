module json_utils_ext

  use json_file_module, only: json_file
  use json_value_module, only: json_value, json_core
  use utils, only: neko_error
  implicit none
  private

  public :: json_key_fallback, json_get_subdict

contains

  !> Create a json_string based on fallback logic.
  function json_key_fallback(json, lookup, fallback) result(string)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: lookup
    character(len=*), intent(in) :: fallback
    character(len=:), allocatable :: string
    type(json_value), pointer :: value

    logical :: status

    logical :: lookup_valid, fallback_valid
    lookup_valid = .false.
    fallback_valid = .false.

    lookup_valid = lookup .in. json
    fallback_valid = fallback .in. json

    if (.not. lookup_valid .and. fallback_valid) then
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

    ! type(json_core) :: core
    type(json_value), pointer :: child
    ! character(len=:), allocatable :: buffer
    logical :: valid

    logical :: status = .false.
    character(len=:), allocatable :: message

    call json%get('case.fluid.initial_condition', child, valid)
    ! call json%check_for_errors(status, message)
    ! if (.not. status) then
    !    call neko_error('Error in json core: '//trim(message))
    ! end if

    if (.not. valid) then
       call neko_error('Parameter '//trim(key)//' missing from the case file')
    end if

    ! if (.not. associated(child)) then
    !    call neko_error('Parameter '//trim(key)//' is not a valid pointer')
    ! end if

    ! call child%print_to_string(buffer)

    ! ! call core%print_to_string(child, buffer)
    ! call output%initialize()
    ! call output%add(child)
    ! call output%check_for_errors(status, message)
    ! if (.not. status) then
    !    call neko_error('Error in json file: '//trim(message))
    ! end if

    ! call output%print()

    ! call neko_error('Error in json file: '//trim(buffer))



    ! call json%get(key, value, valid)

    ! call json%check_for_errors(status, message)

    ! if (.not. status) then
    !    call neko_error('Error in json file: '//trim(message))
    ! end if

    ! if (.not. valid) then
    !    call neko_error('Parameter '//trim(key)//' missing from the case file')
    ! end if



    call neko_error('End of json_get_subdict')

  end subroutine json_get_subdict



end module json_utils_ext
