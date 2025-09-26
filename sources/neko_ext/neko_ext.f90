!> @file neko_ext.f90
!! @brief Contains the module neko_ext
!!
!! This module contains extensions to the neko library required to run the
!! code. It is not part of the neko library itself.

!> @brief Contains extensions to the neko library required to run the topology
!! optimization code.
module neko_ext
  use case, only: case_t
  use json_utils, only: json_get, json_get_or_default, json_extract_item
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
  use scalars, only: scalars_t
  use adjoint_scalars, only: adjoint_scalars_t

  implicit none

  ! ========================================================================= !
  ! Module interface
  ! ========================================================================= !
  private
  public :: reset, field_to_vector, vector_to_field, get_scalar_indicies

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
    integer :: i, n_scalars
    character(len=:), allocatable :: string_val
    logical :: has_scalar, freezeflow
    type(field_t), pointer :: u, v, w, p, s
    type(json_file) :: json_subdict, scalar_params

    ! ------------------------------------------------------------------------ !
    ! Setup shorthand notation
    ! ------------------------------------------------------------------------ !

    u => neko_case%fluid%u
    v => neko_case%fluid%v
    w => neko_case%fluid%w
    p => neko_case%fluid%p
    if (allocated(neko_case%scalars)) then
       s => neko_case%scalars%scalar_fields(1)%s
    else
       nullify(s)
    end if

    ! ------------------------------------------------------------------------ !
    ! Reset the timing parameters
    ! ------------------------------------------------------------------------ !

    t = 0.0_rp
    neko_case%time%t = t
    neko_case%time%tstep = 0

    ! Setup lagged time step parameters
    neko_case%time%tlag = t
    neko_case%time%dtlag = neko_case%time%dt
    do i = 1, size(neko_case%time%tlag)
       neko_case%time%tlag(i) = t - i*neko_case%time%dtlag(i)
    end do

    ! Reset the time step counter
    call neko_case%output_controller%set_counter(neko_case%time)

    ! Restart the fields
    call neko_case%fluid%restart(neko_case%chkp)
    if (allocated(neko_case%scalars)) then
       call neko_case%scalars%restart(neko_case%chkp)
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
    call json_get(neko_case%params, 'case.fluid.initial_condition', &
         json_subdict)

    if (trim(string_val) .ne. 'user') then
       call set_flow_ic(u, v, w, p, &
            neko_case%fluid%c_Xh, neko_case%fluid%gs_Xh, &
            string_val, json_subdict)
    else
       call set_flow_ic(u, v, w, p, &
            neko_case%fluid%c_Xh, neko_case%fluid%gs_Xh, &
            neko_case%user%initial_conditions, neko_case%fluid%name)
    end if

    ! ------------------------------------------------------------------------ !
    ! Reset the scalar field to the initial condition
    ! ------------------------------------------------------------------------ !

    call json_get_or_default(neko_case%params, &
         'case.scalar.enabled', has_scalar, .false.)

    if (has_scalar) then
       if (neko_case%params%valid_path('case.adjoint_scalar')) then
          ! we shouldn't fallback to the primal here.
          call json_get(neko_case%params, &
               'case.adjoint_scalar.initial_condition.type', string_val)
          call json_get(neko_case%params, &
               'case.adjoint_scalar.initial_condition', json_subdict)

          !call neko_log%section("Adjoint scalar initial condition ")

          if (trim(string_val) .ne. 'user') then
             call set_scalar_ic(neko_case%scalars%scalar_fields(1)%s, &
                  neko_case%fluid%c_Xh, neko_case%fluid%gs_Xh, string_val, &
                  json_subdict)
          else
             call neko_error("user defined ICs not implemented for " // &
                  "adjoint scalar")
             ! call set_scalar_ic(this%adjoint_scalars%s_adj, &
             !      this%adjoint_scalars%c_Xh, this%adjoint_scalars%gs_Xh, &
             !      this%usr%scalar_user_ic, neko_case%params)
          end if

          ! call neko_log%end_section()
       else

          ! Handle multiple scalars
          call neko_case%params%info('case.scalars', n_children = n_scalars)

          do i = 1, n_scalars
             call json_extract_item(neko_case%params, 'case.adjoint_scalars', &
                  i, scalar_params)
             call json_get(scalar_params, 'initial_condition.type', string_val)
             call json_get(scalar_params, 'initial_condition', &
                  json_subdict)

             if (trim(string_val) .ne. 'user') then
                call set_scalar_ic(neko_case%scalars%scalar_fields(i)%s, &
                     neko_case%scalars%scalar_fields(i)%c_Xh, &
                     neko_case%scalars%scalar_fields(i)%gs_Xh, string_val, &
                     json_subdict)
             else
                call neko_error("user defined ICs not implemented for " // &
                     "adjoint scalar")
             end if
          end do
       end if
    end if

    ! ------------------------------------------------------------------------ !
    ! Reset the "freeze" parameter of the flow
    ! ------------------------------------------------------------------------ !

    call json_get_or_default(neko_case%params, &
         'case.fluid.freeze_flow', freezeflow, .false.)

    neko_case%fluid%freeze = freezeflow

  end subroutine reset

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

  !> @brief get scalar indices
  !! @details Given a primal scalar name, return the indices in the scalars
  !! and adjoint_scalars list corresponding to this pair.
  !! @param[out] i_primal Index in the primal scalar list.
  !! @param[out] i_adjoint Index in the adjoint scalar list.
  !! @param[inout] scalars Primal scalars list.
  !! @param[inout] adjoint_scalars Adjoint scalars list.
  !! @param[in] primal_name Name of the primal scalar.
  subroutine get_scalar_indicies(i_primal, i_adjoint, scalars, &
       adjoint_scalars, primal_name)
    integer, intent(out) :: i_primal
    integer, intent(out) :: i_adjoint
    type(scalars_t), intent(inout) :: scalars
    type(adjoint_scalars_t), intent(inout) :: adjoint_scalars
    character(len=*), intent(in) :: primal_name
    integer :: i, n_primal_scalars, n_adjoint_scalars

    i_primal = -1
    i_adjoint = -1
    n_adjoint_scalars = size(adjoint_scalars%adjoint_scalar_fields)
    n_primal_scalars = size(scalars%scalar_fields)

    if ((n_adjoint_scalars .eq. 1) .and. (n_primal_scalars .eq. 1)) then
       i_primal = 1
       i_adjoint = 1
       return
    end if

    do i = 1, n_adjoint_scalars
       if (adjoint_scalars%adjoint_scalar_fields(i)%primal_name &
            .eq. primal_name) then
          i_adjoint = i
          exit
       end if
    end do

    do i = 1, n_primal_scalars
       if (scalars%scalar_fields(i)%name .eq. primal_name) then
          i_primal = i
          exit
       end if
    end do

    if (i_primal .le. 0 .or. i_adjoint .le. 0) then
       call neko_error('could not find matching primal and adjoint' // &
            ' scalar fields')
    end if

  end subroutine get_scalar_indicies

end module neko_ext
