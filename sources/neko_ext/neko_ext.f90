!> @file neko_ext.f90
!! @brief Contains the module neko_ext
!!
!! This module contains extensions to the neko library required to run the
!! code. It is not part of the neko library itself.

!> @brief Contains extensions to the neko library required to run the topology
!! optimization code.
module neko_ext
  use case, only: case_t
  implicit none

  ! ========================================================================= !
  ! Module interface
  ! ========================================================================= !
  private
  public :: setup_iteration, reset

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
  !> @param[inout] C Case data structure.
  subroutine reset(neko_case)
    use json_utils, only: json_get, json_get_or_default
    use num_types, only: rp
    use simcomp_executor, only: neko_simcomps
    use flow_ic, only: set_flow_ic
    use scalar_ic, only: set_scalar_ic
    use field, only: field_t
    implicit none

    type(case_t), intent(inout) :: neko_case
    real(kind=rp) :: t
    integer :: i
    character(len=:), allocatable :: string_val
    logical :: has_scalar
    type(field_t), pointer :: u, v, w, p, s

    ! ------------------------------------------------------------------------ !
    ! Setup shorthand notation
    ! ------------------------------------------------------------------------ !

    u => neko_case%fluid%u
    v => neko_case%fluid%v
    w => neko_case%fluid%w
    p => neko_case%fluid%p
    if(allocated(neko_case%scalar)) then
       s => neko_case%scalar%s
    end if

    ! ------------------------------------------------------------------------ !
    ! Reset the timing parameters
    ! ------------------------------------------------------------------------ !

    ! Setup lagged time step parameters
    neko_case%tlag(:) = t
    neko_case%dtlag(:) = neko_case%dt
    do i = 1, size(neko_case%tlag)
       neko_case%tlag(i) = t - i*neko_case%dtlag(i)
    end do

    ! Reset the time step counter
    call neko_case%output_controller%set_counter(t)

    ! Restart the fields
    call neko_case%fluid%restart(neko_case%dtlag, neko_case%tlag)
    if (allocated(neko_case%scalar)) then
       call neko_case%scalar%restart(neko_case%dtlag, neko_case%tlag)
    end if

    ! Reset the external BDF coefficients
    do i = 1, size(neko_case%dtlag)
       call neko_case%ext_bdf%set_coeffs(neko_case%dtlag)
    end do

    ! Restart the simulation components
    call neko_simcomps%restart(t)

    ! ------------------------------------------------------------------------ !
    ! Reset the fluid field to the initial condition
    ! ------------------------------------------------------------------------ !

    call json_get(neko_case%params, &
         'case.fluid.initial_condition.type', string_val)

    if (trim(string_val) .ne. 'user') then
       call set_flow_ic(u, v, w, p, &
            neko_case%fluid%c_Xh, neko_case%fluid%gs_Xh, &
            string_val, neko_case%params)
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

    print *, 'has scalar?', has_scalar

    ! TODO
    ! come back and fix this!!
    !if (has_scalar) then
    !   call json_get(neko_case%params, &
    !        'case.scalar.initial_condition.type', string_val)
    !print *, 'in here'

    !   if (trim(string_val) .ne. 'user') then
    !      call set_scalar_ic(s, &
    !           neko_case%scalar%c_Xh, neko_case%scalar%gs_Xh, &
    !           string_val, &
    !           neko_case%params)
    !   else
    !      call set_scalar_ic(s, &
    !           neko_case%scalar%c_Xh, neko_case%scalar%gs_Xh, &
    !           neko_case%usr%scalar_user_ic, &
    !           neko_case%params)
    !   end if
    !end if

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
    use case, only: case_t
    use num_types, only: rp
    use chkp_output, only: chkp_output_t
    use json_utils, only: json_get_or_default
    use output_controller, only : output_controller_t
    implicit none

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
    call neko_case%output_controller%execute(0.0_rp, 0, .true.)

  end subroutine setup_iteration

end module neko_ext
