module user
  use neko
  use user_initial_conditions, only: scalar_z_split_ic
  implicit none

contains
  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%initial_conditions => initial_conditions
  end subroutine user_setup

  !> User initial condition
  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields

    type(field_t), pointer :: s

    ! See scalar.name in the case file, makes sure that we only
    ! run this for the scalar field.
    if (scheme_name .ne. 'temperature') return

    s => fields%items(1)%ptr

    call scalar_z_split_ic(s, 0.5_rp)
  end subroutine initial_conditions

end module user
