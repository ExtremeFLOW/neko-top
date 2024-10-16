module json_utils_ext

  use json_file_module, only: json_file
  use json_value_module, only: json_value
  use utils, only: neko_error
  implicit none
  private

  public :: json_key_fallback, json_get_subdict

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



end module json_utils_ext
