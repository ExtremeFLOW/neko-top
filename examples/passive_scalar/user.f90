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

    if(y.gt.1.0_rp) then
    	s = 1.0_rp
    else
    	s = 0.0_rp
    end if
  end subroutine scalar_bc

  !> User initial condition                                                     
  subroutine scalar_ic(s, params)                                                  
    type(field_t), intent(inout) :: s                                           
    type(json_file), intent(inout) :: params                                    
    integer :: i, e, k, j                                                       
    real(kind=rp) :: rand, z                                                    
                                                                                
    do i = 1, s%dof%size()                                                      
    	if(s%dof%y(i,1,1,1).gt.1.0_rp) then
       s%x(i,1,1,1) = 1.0_rp                                       
      else
       s%x(i,1,1,1) = 0.0_rp                                       
      end if
    end do                                                                      
                                                                                
    end subroutine scalar_ic

end module user
