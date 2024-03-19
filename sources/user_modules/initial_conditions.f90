module initial_conditions
  use field, only: field_t
  use json_file_module, only: json_file
  implicit none

contains

  !> Set the initial condition for the scalar field
  !! @details This function will initialize the scalar field with a two part
  !! uniform value. Above z=z0 the scalar field will be 0.0 and below z=z0 the
  !! scalar field will be 1.0.
  !! z0 is read from the JSON file under the key
  !! 'case.scalar.initial_condition.value' or set to 0.0 if not found.
  !!
  !! @param[inout] s Scalar field
  !! @param[inout] params JSON file
  subroutine scalar_z_split_ic(s, params)
    use json_utils, only: json_get_or_default
    use num_types, only: rp
    implicit none

    type(field_t), intent(inout) :: s
    type(json_file), intent(inout) :: params

    real(kind=rp) :: z_value, split_value
    integer :: i

    call json_get_or_default(params, &
                             'case.scalar.initial_condition.value', &
                             split_value, 0.0_rp)

    do i = 1, s%dof%size()
       z_value = s%dof%z(i, 1, 1, 1)

       if (z_value .gt. 0.0_rp) then
          s%x(i, 1, 1, 1) = 0.0_rp
       else
          s%x(i, 1, 1, 1) = 1.0_rp
       end if

    end do

  end subroutine scalar_z_split_ic

end module initial_conditions
