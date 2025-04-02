!> @file neko_ext.f90
!! @brief Contains the module neko_ext
!!
!! This module contains extensions to the neko library required to run the
!! code. It is not part of the neko library itself.

!> @brief Contains extensions to the neko library required to run the topology
!! optimization code.
module neko_ext
  use case, only: case_t
  use json_utils, only: json_get, json_get_or_default, json_extract_object
  use num_types, only: rp
  use simcomp_executor, only: neko_simcomps
  use flow_ic, only: set_flow_ic
  use scalar_ic, only: set_scalar_ic
  use field, only: field_t
  use chkp_output, only: chkp_output_t
  use output_controller, only: output_controller_t
  ! for vector/field math
  use math, only: copy
  use device_math, only: device_copy
  use neko_config, only : NEKO_BCKND_DEVICE
  use vector, only: vector_t
  use field, only: field_t
  use utils, only: neko_error
  use json_module, only : json_file
  implicit none

  ! ========================================================================= !
  ! Module interface
  ! ========================================================================= !
  private
  public :: setup_iteration, reset, field_to_vector, vector_to_field

contains

  ! ========================================================================= !
  ! Public routines
  ! ========================================================================= !

  !> @brief Reset the case data structure
  !>
  !> @details This subroutine resets the case data structure. It is called at
  !> the beginning of each iteration. It is used to reset the time step counter,
  !> the lagged time step parameters, the external BDF coefficients, the fluid
  !> and scalar fields, and the simulation components.
  !>
  !> @param[inout] neko_case Case data structure.
  subroutine reset(neko_case)
    type(case_t), intent(inout) :: neko_case
    real(kind=rp) :: t
    integer :: i
    character(len=:), allocatable :: string_val
    logical :: has_scalar, freezeflow
    type(field_t), pointer :: u, v, w, p, s
    type(json_file) :: json_subdict

    t = 0.0_rp

    ! ------------------------------------------------------------------------ !
    ! Setup shorthand notation
    ! ------------------------------------------------------------------------ !

    u => neko_case%fluid%u
    v => neko_case%fluid%v
    w => neko_case%fluid%w
    p => neko_case%fluid%p
    if (allocated(neko_case%scalar)) then
       s => neko_case%scalar%s
    else
       nullify(s)
    end if

    ! ------------------------------------------------------------------------ !
    ! Reset the timing parameters
    ! ------------------------------------------------------------------------ !

    ! Setup lagged time step parameters
    neko_case%time%tlag(:) = t
    neko_case%time%dtlag(:) = neko_case%time%dt
    do i = 1, size(neko_case%time%tlag)
       neko_case%time%tlag(i) = t - i*neko_case%time%dtlag(i)
    end do

    ! Reset the time step counter
    call neko_case%output_controller%set_counter(neko_case%time)

    ! Restart the fields
    call neko_case%fluid%restart(neko_case%chkp)
    if (allocated(neko_case%scalar)) then
       call neko_case%scalar%restart(neko_case%chkp)
    end if

    ! Reset the external BDF coefficients
    do i = 1, size(neko_case%time%dtlag)
       call neko_case%fluid%ext_bdf%set_coeffs(neko_case%time%dtlag)
    end do

    ! Restart the simulation components
    call neko_simcomps%restart(neko_case%time)

    ! ------------------------------------------------------------------------ !
    ! Reset the fluid field to the initial condition
    ! ------------------------------------------------------------------------ !

    call json_get(neko_case%params, &
         'case.fluid.initial_condition.type', string_val)
    call json_extract_object(neko_case%params, 'case.fluid.initial_condition', &
         json_subdict)

    if (trim(string_val) .ne. 'user') then
       call set_flow_ic(u, v, w, p, &
            neko_case%fluid%c_Xh, neko_case%fluid%gs_Xh, &
            string_val, json_subdict)
    else
       call set_flow_ic(u, v, w, p, &
            neko_case%fluid%c_Xh, neko_case%fluid%gs_Xh, &
            neko_case%usr%fluid_user_ic, neko_case%params)
    end if

    ! ------------------------------------------------------------------------ !
    ! Reset the scalar field to the initial condition
    ! ------------------------------------------------------------------------ !

    call json_get_or_default(neko_case%params, &
         'case.scalar.enabled', has_scalar, .false.)

    if (has_scalar) then
       call json_get(neko_case%params, &
            'case.scalar.initial_condition.type', string_val)
       call json_extract_object(neko_case%params, &
            'case.scalar.initial_condition', json_subdict)

       if (trim(string_val) .ne. 'user') then
          call set_scalar_ic(s, &
               neko_case%scalar%c_Xh, neko_case%scalar%gs_Xh, &
               string_val, &
               json_subdict)
       else
          call set_scalar_ic(s, &
               neko_case%scalar%c_Xh, neko_case%scalar%gs_Xh, &
               neko_case%usr%scalar_user_ic, &
               neko_case%params)
       end if
    end if

    ! ------------------------------------------------------------------------ !
    ! Reset the "freeze" parameter of the flow
    ! ------------------------------------------------------------------------ !

    call json_get_or_default(neko_case%params, &
         'case.fluid.freeze_flow', freezeflow, .false.)

    neko_case%fluid%freeze = freezeflow

  end subroutine reset

  !> @brief Setup the iteration
  !!
  !! @details This subroutine sets up the iteration. It is called at the
  !! beginning of each iteration. It is used to save the initial configuration
  !! and to set the output file name.
  !!
  !! @param[inout] neko_case Case data structure.
  !! @param[in] iter Iteration number.
  subroutine setup_iteration(neko_case, iter)
    type(case_t), intent(inout) :: neko_case
    integer, intent(in) :: iter

    character(len=:), allocatable :: dirname
    character(len=80) :: file_name

    if (iter .ne. 1) then
       call reset(neko_case)
    end if

    call json_get_or_default(neko_case%params, &
         'case.output_directory', dirname, './')

    write (file_name, '(a,a,i5.5,a)') &
         trim(adjustl(dirname)), '/topopt_', iter, '_.fld'

    neko_case%f_out%output_t%file_%file_type%fname = trim(file_name)
    neko_case%f_out%output_t%file_%file_type%counter = 0
    neko_case%f_out%output_t%file_%file_type%start_counter = 0
    call neko_case%output_controller%execute(neko_case%time, .true.)

  end subroutine setup_iteration

  !> @brief Vector to field
  !!
  !! @details This subroutine converts a vector to a field in the special case
  !! that they have the same dimension.
  !!
  !! @param[out] field the output field.
  !! @param[in] vector the input vector.
  subroutine vector_to_field(field, vector)
    type(field_t), intent(inout) :: field
    type(vector_t), intent(in) :: vector

    ! first check they're the same size
    if (field%size() .ne. vector%size()) then
       call neko_error("vector and field are not the same size")
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(field%x_d, vector%x_d, field%size())
    else
       call copy(field%x, vector%x, field%size())
    end if

  end subroutine vector_to_field

  !> @brief Field to vector
  !!
  !! @details This subroutine converts a field to a vector in the special case
  !! that they have the same dimension.
  !!
  !! @param[out] vector the output vector.
  !! @param[in] field the input field.
  subroutine field_to_vector(vector, field)
    type(vector_t), intent(inout) :: vector
    type(field_t), intent(in) :: field

    ! first check they're the same size
    if (field%size() .ne. vector%size()) then
       call neko_error("vector and field are not the same size")
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(vector%x_d, field%x_d, field%size())
    else
       call copy(vector%x, field%x, field%size())
    end if

  end subroutine field_to_vector

end module neko_ext
