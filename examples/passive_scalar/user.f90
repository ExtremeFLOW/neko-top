! User module for the user defined simulation component
module user
  use user_intf, only: user_t
  use json_module, only: json_file
  use steady_simcomp, only: steady_simcomp_t
  use simcomp_executor, only: neko_simcomps
  use fluid_user_source_term, only: fluid_user_source_term_t
  use num_types, only : rp
  use field, only : field_t
  use field_registry, only : neko_field_registry
  use math, only : rzero, copy, chsign
  use device_math, only: device_copy, device_cmult
  use neko_config, only: NEKO_BCKND_DEVICE
  use operators, only: curl
  use scratch_registry, only : neko_scratch_registry
  implicit none
  !> Case parameters
  ! To define the initial boundary conditions we don't wish to introduce a
  ! discontinuity,
  !
  !    DONT                    DO (but smoother)
  !        _______               ______
  !       |                     /
  !       |                    /
  ! -------             -------
  !
  ! hence,
  !> we use a logistic function defined by
  !! $ \frac{L}{1 + exp{-k(z-z_0)}}$
  real(kind=rp), parameter :: L = 1.0_rp
  real(kind=rp), parameter :: k = 20.0_rp
  real(kind=rp), parameter :: z_0 = 1.0_rp

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%fluid_user_if => user_inflow_eval
    user%scalar_user_bc => scalar_bc
    user%scalar_user_ic => scalar_ic
  end subroutine user_setup

  !> user-defined boundary condition
  subroutine user_inflow_eval(u, v, w, x, y, z, nx, ny, nz, ix, iy, iz, &
       ie, t, tstep)
    real(kind=rp), intent(inout) :: u
    real(kind=rp), intent(inout) :: v
    real(kind=rp), intent(inout) :: w
    real(kind=rp), intent(in) :: x
    real(kind=rp), intent(in) :: y
    real(kind=rp), intent(in) :: z
    real(kind=rp), intent(in) :: nx
    real(kind=rp), intent(in) :: ny
    real(kind=rp), intent(in) :: nz
    integer, intent(in) :: ix
    integer, intent(in) :: iy
    integer, intent(in) :: iz
    integer, intent(in) :: ie
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    ! Inflow velocity profile is a paraboloid
    u = -0.5_rp * (y - 1.0_rp)**2 - 0.5_rp * (z - 1.0_rp)**2 + 1.0_rp
    v = 0._rp
    w = 0._rp
  end subroutine user_inflow_eval

  !> user-defined boundary condition
  subroutine scalar_bc(s, x, y, z, nx, ny, nz, ix, iy, iz, ie, t, tstep)
    real(kind=rp), intent(inout) :: s
    real(kind=rp), intent(in) :: x
    real(kind=rp), intent(in) :: y
    real(kind=rp), intent(in) :: z
    real(kind=rp), intent(in) :: nx
    real(kind=rp), intent(in) :: ny
    real(kind=rp), intent(in) :: nz
    integer, intent(in) :: ix
    integer, intent(in) :: iy
    integer, intent(in) :: iz
    integer, intent(in) :: ie
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    ! Inflow scalar profile is a sigmoid separating the two species
    s = L / (1.0_rp + exp(-k*(z - z_0)))

  end subroutine scalar_bc

  !> user-defined initial condition
  subroutine scalar_ic(s, params)
    type(field_t), intent(inout) :: s
    type(json_file), intent(inout) :: params
    integer :: i

    ! Initial scalar profile is a sigmoid separating the two species
    do i = 1, s%dof%size()
       s%x(i,1,1,1) = L / (1.0_rp + exp(-k*(s%dof%z(i,1,1,1) - z_0)))
    end do

  end subroutine scalar_ic

end module user
