!> @file objectives_user.f90
!! User-defined inflow for the objective time-window tests.
!!
!! Restores the setup `examples/time_test` used before it became these tests:
!! a paraboloid velocity profile enters the duct, carrying two species split
!! either side of a smooth interface, which leave through the outlet. Every
!! other face is insulated and no-slip.
!!
!! Both halves matter. The velocity profile is shaped to the duct: it vanishes
!! on all four walls, so it agrees with the no-slip condition there instead of
!! being discontinuous at the inlet edges the way a uniform plug is. The
!! scalar influx is what gives the scalar a steady state on the residence time
!! rather than a slow diffusive relaxation, and what makes `scalar_mixing`
!! measure mixing rather than the decay of a trapped blob.
!!
!! The conversion to unit tests replaced the velocity inflow with a plug and
!! the scalar inflow with a zero-flux condition on all six faces plus a
!! `point_zone` initial blob.
module objectives_user
  use num_types, only: rp
  use field, only: field_t
  use field_list, only: field_list_t
  use field_dirichlet, only: field_dirichlet_t
  use time_state, only: time_state_t
  use user_intf, only: user_t
  use neko_config, only: NEKO_BCKND_DEVICE
  use device, only: HOST_TO_DEVICE, DEVICE_TO_HOST, device_memcpy
  implicit none
  private

  public :: objectives_user_setup

  !> Amplitude of the inflow split.
  real(kind=rp), parameter :: split_amplitude = 1.0_rp
  !> Steepness of the split.
  !!
  !! The example used 200 on a 32x8x8 mesh at polynomial order 5. These tests
  !! run 8x4x4 at order 3, which puts thirteen points across the split
  !! direction, so 200 would be a step function in all but name and would ring.
  !! 20 spreads the transition over roughly a fifth of the height, which this
  !! mesh resolves.
  real(kind=rp), parameter :: split_steepness = 20.0_rp
  !> Height at which the two species meet.
  real(kind=rp), parameter :: split_height = 0.5_rp

  !> Peak of the inflow velocity profile.
  !!
  !! 36 is not arbitrary: it makes the mean of `y(y-1)z(z-1)` over the unit
  !! square equal one, so the duct carries unit flow rate.
  real(kind=rp), parameter :: velocity_scale = 36.0_rp

contains

  !> Register the user routines on a case's `user_t`.
  !!
  !! Must be called before the case is initialized: `user_intf_init` only
  !! substitutes its own defaults for pointers that are still null, so
  !! anything assigned here survives.
  !! @param user The user interface to populate.
  subroutine objectives_user_setup(user)
    type(user_t), intent(inout) :: user

    user%dirichlet_conditions => inflow
    user%initial_conditions => scalar_initial_condition
  end subroutine objectives_user_setup

  !> The scalar concentration at a point, as a smooth split in `z`.
  !! @param z Height.
  !! @return The concentration.
  elemental function split_profile(z) result(phi)
    real(kind=rp), intent(in) :: z
    real(kind=rp) :: phi

    phi = split_amplitude / &
         (1.0_rp + exp(-split_steepness * (z - split_height)))
  end function split_profile

  !> The streamwise velocity at a point on the inlet face.
  !!
  !! A paraboloid over the duct cross-section, shaped so that it vanishes on
  !! all four walls and agrees with the no-slip condition applied there.
  !! @param y Spanwise coordinate.
  !! @param z Height.
  !! @return The streamwise velocity.
  pure function velocity_profile(y, z) result(u)
    real(kind=rp), intent(in) :: y
    real(kind=rp), intent(in) :: z
    real(kind=rp) :: u

    u = velocity_scale * y * (y - 1.0_rp) * z * (z - 1.0_rp)
  end function velocity_profile

  !> Impose the inflow profiles on the inlet boundary.
  !!
  !! Called for the velocity fields and for the scalar separately; the first
  !! field's name says which.
  !! @param fields The fields the boundary condition applies to.
  !! @param bc The boundary condition, carrying the mask.
  !! @param time The current time state.
  subroutine inflow(fields, bc, time)
    type(field_list_t), intent(inout) :: fields
    type(field_dirichlet_t), intent(in) :: bc
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: u, v, w, s
    real(kind=rp) :: y, z
    integer :: i, idx

    if (fields%items(1)%ptr%name .eq. 'u') then
       u => fields%get("u")
       v => fields%get("v")
       w => fields%get("w")

       call u%copy_from(DEVICE_TO_HOST, sync = .false.)
       call v%copy_from(DEVICE_TO_HOST, sync = .false.)
       call w%copy_from(DEVICE_TO_HOST, sync = .true.)

       do i = 1, bc%msk(0)
          idx = bc%msk(i)
          y = u%dof%y%x(idx, 1, 1, 1)
          z = u%dof%z%x(idx, 1, 1, 1)
          u%x(idx, 1, 1, 1) = velocity_profile(y, z)
          v%x(idx, 1, 1, 1) = 0.0_rp
          w%x(idx, 1, 1, 1) = 0.0_rp
       end do

       call u%copy_from(HOST_TO_DEVICE, sync = .false.)
       call v%copy_from(HOST_TO_DEVICE, sync = .false.)
       call w%copy_from(HOST_TO_DEVICE, sync = .true.)

       nullify(u, v, w)
    else
       s => fields%get("s")
       call s%copy_from(DEVICE_TO_HOST, sync = .true.)

       do i = 1, bc%msk(0)
          idx = bc%msk(i)
          z = s%dof%z%x(idx, 1, 1, 1)
          s%x(idx, 1, 1, 1) = split_profile(z)
       end do

       call s%copy_from(HOST_TO_DEVICE, sync = .true.)
       nullify(s)
    end if
  end subroutine inflow

  !> Start the scalar from the same split the inflow imposes.
  !!
  !! Starting anywhere else only adds a transient the run then has to sit
  !! through before it can converge.
  !! @param scheme_name The scheme requesting an initial condition.
  !! @param fields The fields to initialize.
  subroutine scalar_initial_condition(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields
    type(field_t), pointer :: s

    if (scheme_name .eq. 'fluid') return

    s => fields%get("s")
    s%x = split_profile(s%dof%z%x)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(s%x, s%x_d, s%size(), HOST_TO_DEVICE, sync = .false.)
    end if
    nullify(s)
  end subroutine scalar_initial_condition

end module objectives_user
