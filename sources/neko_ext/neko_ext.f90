!> @file neko_ext.f90
!! @brief Contains the module neko_ext
!!
!! This module contains extensions to the neko library required to run the
!! code. It is not part of the neko library itself.
!!
!! @copyright
!! Copyright (c) 2024-2025, The Neko-TOP Authors
!! All rights reserved.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions
!! are met:
!!
!!   * Redistributions of source code must retain the above copyright
!!     notice, this list of conditions and the following disclaimer.
!!
!!   * Redistributions in binary form must reproduce the above
!!     copyright notice, this list of conditions and the following
!!     disclaimer in the documentation and/or other materials provided
!!     with the distribution.
!!
!!   * Neither the name of the authors nor the names of its
!!     contributors may be used to endorse or promote products derived
!!     from this software without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
!! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
!! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
!! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
!! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
!! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
!! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!! POSSIBILITY OF SUCH DAMAGE.

!> @brief Contains extensions to the neko library required to run the topology
!! optimization code.
module neko_ext
  use case, only: case_t
  use adjoint_case, only: adjoint_case_t
  use json_utils, only: json_get, json_get_or_default
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
  use field_math, only: field_rzero, field_copy
  use fluid_pnpn, only: fluid_pnpn_t
  use adjoint_fluid_pnpn, only: adjoint_fluid_pnpn_t

  implicit none

  ! ========================================================================= !
  ! Module interface
  ! ========================================================================= !
  private
  public :: reset, reset_adjoint, field_to_vector, vector_to_field, &
       get_scalar_indicies

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

    ! ------------------------------------------------------------------------ !
    ! Setup shorthand notation
    ! ------------------------------------------------------------------------ !

    u => neko_case%fluid%u
    v => neko_case%fluid%v
    w => neko_case%fluid%w
    p => neko_case%fluid%p
    if (allocated(neko_case%scalars)) then
       s => neko_case%scalars%scalar_fields(1)%scalar%s
    else
       nullify(s)
    end if

    ! ------------------------------------------------------------------------ !
    ! Reset the timing parameters
    ! ------------------------------------------------------------------------ !

    call neko_case%time%reset()
    t = neko_case%time%start_time
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

    ! set lags to IC
    call neko_case%fluid%ulag%set(u)
    call neko_case%fluid%vlag%set(v)
    call neko_case%fluid%wlag%set(w)
    ! zero out RHS etc
    select type (f => neko_case%fluid)
    type is (fluid_pnpn_t)
       call field_rzero(f%abx1)
       call field_rzero(f%aby1)
       call field_rzero(f%abz1)
       call field_rzero(f%abx2)
       call field_rzero(f%aby2)
       call field_rzero(f%abz2)
       call field_copy(f%u_e, u)
       call field_copy(f%v_e, v)
       call field_copy(f%w_e, w)
    end select
    call field_rzero(neko_case%fluid%f_x)
    call field_rzero(neko_case%fluid%f_y)
    call field_rzero(neko_case%fluid%f_z)
    ! ------------------------------------------------------------------------ !
    ! Reset the scalar field to the initial condition
    ! ------------------------------------------------------------------------ !

    ! check for a single scalar
    call json_get_or_default(neko_case%params, &
         'case.scalar.enabled', has_scalar, .false.)

    if (has_scalar) then
       ! check for multiple scalars
       if (size(neko_case%scalars%scalar_fields) .gt. 1) then
          call neko_error('Multiple scalars not supported')
       end if
       ! zero out RHS
       call field_rzero(neko_case%scalars%scalar_fields(1)%scalar%f_Xh)
       ! reset the forward scalar
       call json_get(neko_case%params, &
            'case.scalar.initial_condition.type', string_val)
       call json_get(neko_case%params, &
            'case.scalar.initial_condition', json_subdict)
       if (trim(string_val) .ne. 'user') then
          if (trim(neko_case%scalars%scalar_fields(1)%scalar%name) .eq. &
               'temperature') then
             call set_scalar_ic(neko_case%scalars%scalar_fields(1)%scalar%s, &
                  neko_case%fluid%c_Xh, neko_case%fluid%gs_Xh, string_val, &
                  json_subdict, 0)
          else
             call set_scalar_ic(neko_case%scalars%scalar_fields(1)%scalar%s, &
                  neko_case%fluid%c_Xh, neko_case%fluid%gs_Xh, string_val, &
                  json_subdict, 1)
          end if
       else
          call set_scalar_ic(neko_case%scalars%scalar_fields(1)%scalar%name, &
               neko_case%scalars%scalar_fields(1)%scalar%s, &
               neko_case%scalars%scalar_fields(1)%scalar%c_Xh, &
               neko_case%scalars%scalar_fields(1)%scalar%gs_Xh, &
               neko_case%user%initial_conditions)
       end if
       ! set lags to IC
       call neko_case%scalars%scalar_fields(1)%scalar%slag%set(&
            neko_case%scalars%scalar_fields(1)%scalar%s)
    end if

    ! ------------------------------------------------------------------------ !
    ! Reset the "freeze" parameter of the flow
    ! ------------------------------------------------------------------------ !

    call json_get_or_default(neko_case%params, &
         'case.fluid.freeze_flow', freezeflow, .false.)

    neko_case%fluid%freeze = freezeflow

  end subroutine reset

  !> @brief Reset the adjoint case data structure
  !!
  !! @details This subroutine resets the adjoint case data structure. It is
  !! called at the beginning of each iteration. It is used to reset the time
  !! step counter, the lagged time step parameters, the external BDF
  !! coefficients, the adjoint fluid_adj and adjoint scalar fields.
  !!
  !! @param[inout] adjoint_case Adjoint case data structure.
  !! @param[in] neko_case Primal case.
  subroutine reset_adjoint(adjoint_case, neko_case)
    type(adjoint_case_t), intent(inout) :: adjoint_case
    type(case_t), intent(inout) :: neko_case
    real(kind=rp) :: t
    integer :: i
    character(len=:), allocatable :: string_val
    logical :: has_scalar, freezeflow
    type(field_t), pointer :: u_adj, v_adj, w_adj, p_adj, s_adj
    type(json_file) :: json_subdict

    ! ------------------------------------------------------------------------ !
    ! Setup shorthand notation
    ! ------------------------------------------------------------------------ !

    u_adj => adjoint_case%fluid_adj%u_adj
    v_adj => adjoint_case%fluid_adj%v_adj
    w_adj => adjoint_case%fluid_adj%w_adj
    p_adj => adjoint_case%fluid_adj%p_adj
    if (allocated(adjoint_case%adjoint_scalars)) then
       s_adj => adjoint_case%adjoint_scalars%adjoint_scalar_fields(1)%s_adj
    else
       nullify(s_adj)
    end if

    ! ------------------------------------------------------------------------ !
    ! Reset the timing parameters
    ! ------------------------------------------------------------------------ !

    call adjoint_case%time%reset()
    t = adjoint_case%time%start_time
    do i = 1, size(adjoint_case%time%tlag)
       adjoint_case%time%tlag(i) = t - i*adjoint_case%time%dtlag(i)
    end do

    ! Reset the time step counter
    call adjoint_case%output_controller%set_counter(adjoint_case%time)
    if (adjoint_case%norm_output_enabled) then
       call adjoint_case%norm_output_ctrl%set_counter(adjoint_case%time)
    end if

    ! Reset the external BDF coefficients
    do i = 1, size(adjoint_case%time%dtlag)
       call adjoint_case%fluid_adj%ext_bdf%set_coeffs(adjoint_case%time%dtlag)
    end do

    ! ------------------------------------------------------------------------ !
    ! Reset the adjoint fluid to the initial (final) condition
    ! ------------------------------------------------------------------------ !

    ! don't fallback to the fluid here
    call json_get(neko_case%params, &
         'case.adjoint_fluid.initial_condition.type', string_val)
    call json_get(neko_case%params, 'case.adjoint_fluid.initial_condition', &
         json_subdict)

    if (trim(string_val) .ne. 'user') then
       call set_flow_ic(u_adj, v_adj, w_adj, p_adj, &
            adjoint_case%fluid_adj%c_Xh, adjoint_case%fluid_adj%gs_Xh, &
            string_val, json_subdict)
    else
       call neko_error("adjoint user initial conditions not supported")
    end if

    ! set lags to IC
    call adjoint_case%fluid_adj%ulag%set(u_adj)
    call adjoint_case%fluid_adj%vlag%set(v_adj)
    call adjoint_case%fluid_adj%wlag%set(w_adj)
    ! zero out RHS etc
    select type (f => adjoint_case%fluid_adj)
    type is (adjoint_fluid_pnpn_t)
       call field_rzero(f%abx1)
       call field_rzero(f%aby1)
       call field_rzero(f%abz1)
       call field_rzero(f%abx2)
       call field_rzero(f%aby2)
       call field_rzero(f%abz2)
    end select
    ! zero out all lags etc
    ! (not sure what to do with the abx's_adj)
    call field_rzero(adjoint_case%fluid_adj%f_adj_x)
    call field_rzero(adjoint_case%fluid_adj%f_adj_y)
    call field_rzero(adjoint_case%fluid_adj%f_adj_z)
    ! ------------------------------------------------------------------------ !
    ! Reset the scalar field to the initial condition
    ! ------------------------------------------------------------------------ !

    ! check for a single scalar
    call json_get_or_default(neko_case%params, 'case.scalar.enabled', &
         has_scalar, .false.)

    if (has_scalar) then
       ! check for multiple adjoint_scalars
       if (size(adjoint_case%adjoint_scalars%adjoint_scalar_fields) .gt. 1) then
          call neko_error('Multiple adjoint scalars not supported')
       end if
       ! zero out lag terms
       call field_rzero( &
            adjoint_case%adjoint_scalars%adjoint_scalar_fields(1)%f_Xh)
       ! reset the forward scalar
       call json_get(neko_case%params, &
            'case.adjoint_scalar.initial_condition.type', string_val)
       call json_get(neko_case%params, &
            'case.adjoint_scalar.initial_condition', json_subdict)
       if (trim(string_val) .ne. 'user') then
          if (trim(neko_case%scalars%scalar_fields(1)%scalar%name) .eq. &
               'temperature') then
             call set_scalar_ic( &
                  adjoint_case%adjoint_scalars%adjoint_scalar_fields(1)%s_adj, &
                  adjoint_case%fluid_adj%c_Xh, adjoint_case%fluid_adj%gs_Xh, &
                  string_val, json_subdict, 0)
          else
             call set_scalar_ic( &
                  adjoint_case%adjoint_scalars%adjoint_scalar_fields(1)%s_adj, &
                  adjoint_case%fluid_adj%c_Xh, adjoint_case%fluid_adj%gs_Xh, &
                  string_val, json_subdict, 1)
          end if
       else
          call neko_error("adjoint scalar user IC not supported")
       end if
       ! set lags to IC
       call adjoint_case%adjoint_scalars%adjoint_scalar_fields(1)%s_adj_lag% &
            set(adjoint_case%adjoint_scalars%adjoint_scalar_fields(1)%s_adj)
    end if

    ! ------------------------------------------------------------------------ !
    ! Reset the "freeze" parameter of the flow
    ! ------------------------------------------------------------------------ !

    call json_get_or_default(neko_case%params, &
         'case.adjoint_fluid.freeze_flow', freezeflow, .false.)

    adjoint_case%fluid_adj%freeze = freezeflow

  end subroutine reset_adjoint

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
       if (scalars%scalar_fields(i)%scalar%name .eq. primal_name) then
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
