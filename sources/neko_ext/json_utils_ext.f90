module json_utils_ext
  use json_module, only: json_file
  implicit none
  private

  public :: json_key_fallback

contains

  !> Create a json_string based on fallback logic.
  function json_key_fallback(json, lookup, fallback) result(string)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: lookup
    character(len=*), intent(in) :: fallback
    character(len=:), allocatable :: string

    logical :: lookup_valid, fallback_valid

    lookup_valid = json%valid_path(lookup)
    fallback_valid = json%valid_path(fallback)

    if (.not. lookup_valid .and. fallback_valid) then
       string = fallback
    else
       string = lookup
    end if

  end function json_key_fallback

end module json_utils_ext
