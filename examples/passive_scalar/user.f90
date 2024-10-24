! User module for the user defined simulation component
module user
  use user_intf, only: user_t, simulation_component_user_settings
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

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%init_user_simcomp => user_simcomp
    user%fluid_user_if => fluid_bc
    user%scalar_user_bc => scalar_bc
    user%scalar_user_ic => scalar_ic
  end subroutine user_setup

  subroutine user_simcomp(params)
    type(json_file), intent(inout) :: params
    type(steady_simcomp_t), allocatable :: steady_comp
    type(json_file) :: simcomp_settings

    ! Allocate a simulation component
    allocate(steady_comp)
    simcomp_settings = simulation_component_user_settings("steady", params)

    call neko_simcomps%add_user_simcomp(steady_comp, simcomp_settings)

  end subroutine user_simcomp

    ! user-defined boundary condition
    subroutine fluid_bc(u, v, w, x, y, z, nx, ny, nz, ix, iy, iz, ie, t, tstep)
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

      ! Casper said he used a parabloid, which is not a solution to NS but
      ! I suppose it will sort itself out with enough distance...
   
      u = -0.5_rp * (y - 1.0_rp)**2 - 0.5_rp * (z - 1.0_rp)**2 + 1.0_rp
      v = 0._rp
      w = 0._rp
   
    end subroutine fluid_bc

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
    real(kind=rp) :: L, k, z_0
    ! TODO
    ! OK here I'm doing something different to Casper.
    ! I feel like since we get a term in the adjoint velocity equation that 
    ! looks like s_adj * grad(s), it's not a good idea to have discontinuity
    ! on this boundary.
    ! ie,
    ! 
    !    DONT                    DO (but smoother)
    !        _______               ______
    !       |                     /
    !       |                    /
    ! -------             -------
    !
    !
    ! A logistic function seems reasonable...
    L = 1.0_rp
    k = 6.0_rp
    z_0 = 1.0_rp

    s = L / (1.0_rp + exp(-k*(z - z_0)))

  end subroutine scalar_bc

  !> User initial condition                                                     
  subroutine scalar_ic(s, params)                                                  
    type(field_t), intent(inout) :: s                                           
    type(json_file), intent(inout) :: params                                    
    integer :: i
    real(kind=rp) :: L, k, z_0

    L = 1.0_rp
    k = 6.0_rp
    z_0 = 1.0_rp
                                                                                
    do i = 1, s%dof%size()                                                      
       s%x(i,1,1,1) = L / (1.0_rp + exp(-k*(s%dof%z(i,1,1,1) - z_0)))
    end do                                                                      
                                                                                
    end subroutine scalar_ic

end module user
