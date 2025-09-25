module user
  use neko
  use user_initial_conditions, only: scalar_z_split_ic
  implicit none

contains
  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%initial_conditions => initial_conditions
    user%dirichlet_conditions => dirichlet_update
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

  subroutine dirichlet_update(fields, bc, time)
    type(field_list_t), intent(inout) :: fields
    type(field_dirichlet_t), intent(in) :: bc
    type(time_state_t), intent(in) :: time

    type(field_t), pointer :: s
    integer :: i
    real(kind=rp) :: z

    if (time%tstep .ne. 1) return

    ! ------------------------------------------------------------------------ #
    ! Assign the boundary condition for the scalar field here.

    s => fields%get("temperature")

    do i = 1, bc%msk(0)
       z = bc%dof%z(bc%msk(i), 1, 1, 1)
       if (z .gt. 0.5_rp) then
          s%x(bc%msk(i), 1, 1, 1) = 1.0_rp
       else
          s%x(bc%msk(i), 1, 1, 1) = 0.0_rp
       end if
    end do

    call s%copy_from(HOST_TO_DEVICE, .true.)

  end subroutine dirichlet_update
end module user
