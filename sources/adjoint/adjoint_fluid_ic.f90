! Copyright (c) 2021, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> Initial flow condition
module adjoint_fluid_ic
  use num_types, only : rp
  use gather_scatter, only : gs_t, GS_OP_ADD
  use neko_config, only : NEKO_BCKND_DEVICE
  use flow_profile, only : blasius_profile, blasius_linear, blasius_cubic, &
       blasius_quadratic, blasius_quartic, blasius_sin
  use device, only: device_memcpy, HOST_TO_DEVICE
  use field, only : field_t
  use utils, only : neko_error
  use coefs, only : coef_t
  use math, only : col2, cfill, cfill_mask
  use device_math, only : device_col2, device_cfill, device_cfill_mask
  use user_intf, only : useric
  use json_module, only : json_file
  use json_utils, only: json_get
  use point_zone, only: point_zone_t
  use point_zone_registry, only: neko_point_zone_registry
  implicit none
  private

  interface set_adjoint_fluid_ic
     module procedure set_adjoint_fluid_ic_int, set_adjoint_fluid_ic_usr
  end interface set_adjoint_fluid_ic

  public :: set_adjoint_fluid_ic

contains

  !> Set initial flow condition (builtin)
  subroutine set_adjoint_fluid_ic_int(uadj, v_adj, w_adj, p_adj, coef, gs, &
     type, params)
    type(field_t), intent(inout) :: uadj
    type(field_t), intent(inout) :: v_adj
    type(field_t), intent(inout) :: w_adj
    type(field_t), intent(inout) :: p_adj
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs
    character(len=*) :: type
    type(json_file), intent(inout) :: params
    real(kind=rp) :: delta
    real(kind=rp), allocatable :: uinf(:)
    real(kind=rp), allocatable :: zone_value(:)
    character(len=:), allocatable :: blasius_approximation
    character(len=:), allocatable :: zone_name

    character(len=:), allocatable :: type_
    call json_get(params, 'type', type_)

    if (trim(type_) .eq. 'uniform') then
       call json_get(params, 'value', uinf)
       call set_adjoint_fluid_ic_uniform(uadj, v_adj, w_adj, uinf)
    else if (trim(type_) .eq. 'blasius') then
       call json_get(params, 'blasius.delta', delta)
       call json_get(params, 'blasius.approximation', &
            blasius_approximation)
       call json_get(params, 'blasius.freestream_velocity', uinf)
       call set_adjoint_fluid_ic_blasius(uadj, v_adj, w_adj, delta, uinf, &
          blasius_approximation)
    else if (trim(type_) .eq. 'point_zone') then
       call json_get(params, 'base_value', uinf)
       call json_get(params, 'zone_name', &
            zone_name)
       call json_get(params, 'zone_value', &
            zone_value)
       call set_adjoint_fluid_ic_point_zone(uadj, v_adj, w_adj, uinf, &
          zone_name, zone_value)
    else
       call neko_error('Invalid initial condition')
    end if

    call set_adjoint_fluid_ic_common(uadj, v_adj, w_adj, p_adj, coef, gs)

  end subroutine set_adjoint_fluid_ic_int

  !> Set intial flow condition (user defined)
  subroutine set_adjoint_fluid_ic_usr(uadj, v_adj, w_adj, p_adj, coef, gs, &
     usr_ic, params)
    type(field_t), intent(inout) :: uadj
    type(field_t), intent(inout) :: v_adj
    type(field_t), intent(inout) :: w_adj
    type(field_t), intent(inout) :: p_adj
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs
    procedure(useric) :: usr_ic
    type(json_file), intent(inout) :: params

    call usr_ic(uadj, v_adj, w_adj, p_adj, params)

    call set_adjoint_fluid_ic_common(uadj, v_adj, w_adj, p_adj, coef, gs)

  end subroutine set_adjoint_fluid_ic_usr

  subroutine set_adjoint_fluid_ic_common(uadj, v_adj, w_adj, p_adj, coef, gs)
    type(field_t), intent(inout) :: uadj
    type(field_t), intent(inout) :: v_adj
    type(field_t), intent(inout) :: w_adj
    type(field_t), intent(inout) :: p_adj
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs
    integer :: n

    n = uadj%dof%size()

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(uadj%x, uadj%x_d, n, &
            HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(v_adj%x, v_adj%x_d, n, &
            HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(w_adj%x, w_adj%x_d, n, &
            HOST_TO_DEVICE, sync = .false.)
    end if

    ! Ensure continuity across elements for initial conditions
    call gs%op(uadj%x, uadj%dof%size(), GS_OP_ADD)
    call gs%op(v_adj%x, v_adj%dof%size(), GS_OP_ADD)
    call gs%op(w_adj%x, w_adj%dof%size(), GS_OP_ADD)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col2(uadj%x_d, coef%mult_d, uadj%dof%size())
       call device_col2(v_adj%x_d, coef%mult_d, v_adj%dof%size())
       call device_col2(w_adj%x_d, coef%mult_d, w_adj%dof%size())
    else
       call col2(uadj%x, coef%mult, uadj%dof%size())
       call col2(v_adj%x, coef%mult, v_adj%dof%size())
       call col2(w_adj%x, coef%mult, w_adj%dof%size())
    end if

  end subroutine set_adjoint_fluid_ic_common

  !> Uniform initial condition
  subroutine set_adjoint_fluid_ic_uniform(uadj, v_adj, w_adj, uinf)
    type(field_t), intent(inout) :: uadj
    type(field_t), intent(inout) :: v_adj
    type(field_t), intent(inout) :: w_adj
    real(kind=rp), intent(in) :: uinf(3)
    integer :: n
    uadj = uinf(1)
    v_adj = uinf(2)
    w_adj = uinf(3)
    n = uadj%dof%size()

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cfill(uadj%x_d, uinf(1), n)
       call device_cfill(v_adj%x_d, uinf(2), n)
       call device_cfill(w_adj%x_d, uinf(3), n)
    else
       call cfill(uadj%x, uinf(1), n)
       call cfill(v_adj%x, uinf(2), n)
       call cfill(w_adj%x, uinf(3), n)
    end if

  end subroutine set_adjoint_fluid_ic_uniform

  !> Set a Blasius profile as initial condition
  !! @note currently limited to axis aligned flow
  subroutine set_adjoint_fluid_ic_blasius(uadj, v_adj, w_adj, delta, uinf, type)
    type(field_t), intent(inout) :: uadj
    type(field_t), intent(inout) :: v_adj
    type(field_t), intent(inout) :: w_adj
    real(kind=rp), intent(in) :: delta
    real(kind=rp), intent(in) :: uinf(3)
    character(len=*), intent(in) :: type
    procedure(blasius_profile), pointer :: bla => null()
    integer :: i

    select case (trim(type))
      case ('linear')
       bla => blasius_linear
      case ('quadratic')
       bla => blasius_quadratic
      case ('cubic')
       bla => blasius_cubic
      case ('quartic')
       bla => blasius_quartic
      case ('sin')
       bla => blasius_sin
      case default
       call neko_error('Invalid Blasius approximation')
    end select

    if ((uinf(1) .gt. 0.0_rp) .and. (uinf(2) .le. 0.0_rp) &
         .and. (uinf(3) .le. 0.0_rp)) then
       do i = 1, uadj%dof%size()
          uadj%x(i,1,1,1) = bla(uadj%dof%z(i,1,1,1), delta, uinf(1))
          v_adj%x(i,1,1,1) = 0.0_rp
          w_adj%x(i,1,1,1) = 0.0_rp
       end do
    else if ((uinf(1) .le. 0.0_rp) .and. (uinf(2) .gt. 0.0_rp) &
         .and. (uinf(3) .le. 0.0_rp)) then
       do i = 1, uadj%dof%size()
          uadj%x(i,1,1,1) = 0.0_rp
          v_adj%x(i,1,1,1) = bla(uadj%dof%x(i,1,1,1), delta, uinf(2))
          w_adj%x(i,1,1,1) = 0.0_rp
       end do
    else if ((uinf(1) .le. 0.0_rp) .and. (uinf(2) .le. 0.0_rp) &
         .and. (uinf(3) .gt. 0.0_rp)) then
       do i = 1, uadj%dof%size()
          uadj%x(i,1,1,1) = 0.0_rp
          v_adj%x(i,1,1,1) = 0.0_rp
          w_adj%x(i,1,1,1) = bla(uadj%dof%y(i,1,1,1), delta, uinf(3))
       end do
    end if

  end subroutine set_adjoint_fluid_ic_blasius

  !> Set the initial condition of the flow based on a point zone.
  !! @details The initial condition is set to the base value and then the
  !! zone is filled with the zone value.
  !! @param uadj The x-component of the adjoint velocity field.
  !! @param v_adj The y-component of the adjoint velocity field.
  !! @param w_adj The z-component of the adjoint velocity field.
  !! @param base_value The base value of the initial condition.
  !! @param zone_name The name of the point zone.
  !! @param zone_value The value of the point zone.
  subroutine set_adjoint_fluid_ic_point_zone(uadj, v_adj, w_adj, base_value, &
     zone_name, zone_value)
    type(field_t), intent(inout) :: uadj
    type(field_t), intent(inout) :: v_adj
    type(field_t), intent(inout) :: w_adj
    real(kind=rp), intent(in), dimension(3) :: base_value
    character(len=*), intent(in) :: zone_name
    real(kind=rp), intent(in) :: zone_value(:)

    ! Internal variables
    class(point_zone_t), pointer :: zone
    integer :: size

    call set_adjoint_fluid_ic_uniform(uadj, v_adj, w_adj, base_value)
    size = uadj%dof%size()

    zone => neko_point_zone_registry%get_point_zone(trim(zone_name))

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cfill_mask(uadj%x_d, zone_value(1), size, &
            zone%mask_d, zone%size)
       call device_cfill_mask(v_adj%x_d, zone_value(2), size, &
            zone%mask_d, zone%size)
       call device_cfill_mask(w_adj%x_d, zone_value(3), size, &
            zone%mask_d, zone%size)
    else
       call cfill_mask(uadj%x, zone_value(1), size, zone%mask, zone%size)
       call cfill_mask(v_adj%x, zone_value(2), size, zone%mask, zone%size)
       call cfill_mask(w_adj%x, zone_value(3), size, zone%mask, zone%size)

    end if
  end subroutine set_adjoint_fluid_ic_point_zone

end module adjoint_fluid_ic
