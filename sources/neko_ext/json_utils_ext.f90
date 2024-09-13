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

    string = trim(adjustl(lookup))
    if (.not. json%valid_path(string)) then
       string = trim(adjustl(fallback))
    end if
  end function json_key_fallback

end module json_utils_ext
