!> @file adjoint_fluid_pnpn.f90
module adjoint_fluid_pnpn
  use perturb_fluid_pnpn, only: perturb_fluid_pnpn_t, &
       perturb_fluid_pnpn_step_adjoint
  use checkpoint, only: chkp_t
  use json_module, only: json_file
  use mesh, only: mesh_t
  use user_intf, only: user_t
  use time_state, only: time_state_t
  use time_step_controller, only: time_step_controller_t
  implicit none
  private

  type, public, extends(perturb_fluid_pnpn_t) :: adjoint_fluid_pnpn_t
   contains
     procedure, pass(this) :: init => adjoint_fluid_pnpn_init
     procedure, pass(this) :: step => adjoint_fluid_pnpn_step
  end type adjoint_fluid_pnpn_t

contains

  subroutine adjoint_fluid_pnpn_init(this, msh, lx, params, user, chkp, &
       if_adjoint)
    class(adjoint_fluid_pnpn_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(in) :: lx
    type(json_file), target, intent(inout) :: params
    type(user_t), target, intent(in) :: user
    type(chkp_t), target, intent(inout) :: chkp
    logical, intent(in) :: if_adjoint

    call this%init_common(msh, lx, params, user, chkp, .true.)

  end subroutine adjoint_fluid_pnpn_init

  subroutine adjoint_fluid_pnpn_step(this, time, dt_controller)
    class(adjoint_fluid_pnpn_t), target, intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(time_step_controller_t), intent(in) :: dt_controller

    call perturb_fluid_pnpn_step_adjoint(this, time, dt_controller)

  end subroutine adjoint_fluid_pnpn_step

end module adjoint_fluid_pnpn
