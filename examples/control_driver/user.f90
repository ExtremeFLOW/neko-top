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
    user%fluid_user_f_vector => adjoint_forcing
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

  !> Forcing
  subroutine adjoint_forcing(f, t)
    class(fluid_user_source_term_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t
    integer :: i
    type(field_t), pointer :: u, v, w
    !type(field_t), pointer :: dudx, dudy, dudz
    !type(field_t), pointer :: dvdx, dvdy, dvdz
    !type(field_t), pointer :: dwdx, dwdy, dwdz
    type(field_t), pointer :: wo1, wo2, wo3, wo4, wo5, wo6
    type(field_t), pointer :: t1 , t2
    integer :: temp_indices(8)
    integer n
    print *, "called forcing"


    n = f%dm%size()

    ! fuck I'm not sure about this... I need a pen and paper
    ! also there should be a way to pre-process this forcing term...
    ! instead of recalculating it every time
    u => neko_field_registry%get_field('u')
    v => neko_field_registry%get_field('v')
    w => neko_field_registry%get_field('w')
    call neko_scratch_registry%request_field(wo1, temp_indices(1))
    call neko_scratch_registry%request_field(wo2, temp_indices(2))
    call neko_scratch_registry%request_field(wo3, temp_indices(3))
    call neko_scratch_registry%request_field(wo4, temp_indices(4))
    call neko_scratch_registry%request_field(wo5, temp_indices(5))
    call neko_scratch_registry%request_field(wo6, temp_indices(6))
    call neko_scratch_registry%request_field(t1 , temp_indices(7))
    call neko_scratch_registry%request_field(t2 , temp_indices(8))

    ! ok we're computing gradients at every timestep... which is stupid...
    ! BUT
    ! if this was unsteady we would have to do this.

    ! this is cheating a little bit...
    ! in strong form, \nabla u . \nabla v =>  v . \nabla^2 u + bdry
    !
    ! we can do this properly later in weak form, ideally using ax_helm or so
    !
    ! for now we'll work in strong form and ignore the bdry
    ! and suffer the double derivative :/
    !
    ! in fact, we'll go even quicker and use
    ! \nabla ^2 u = grad (div (u)) - curl ( curl (u )) and assume divergence free u

    call curl(wo1, wo2, wo3, u, v, w, t1, t2, f%coef)
    call curl(wo4, wo5, wo6, wo1, wo2, wo3, t1, t2, f%coef)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(f%u_d, wo4%x_d, n)
       call device_copy(f%v_d, wo5%x_d, n)
       call device_copy(f%w_d, wo6%x_d, n)
       call device_cmult(f%u_d, -1.0_rp, n)
       call device_cmult(f%v_d, -1.0_rp, n)
       call device_cmult(f%w_d, -1.0_rp, n)
    else
       call copy(f%u, wo4%x, n)
       call copy(f%v, wo5%x, n)
       call copy(f%w, wo6%x, n)
       call chsign(f%u, n)
       call chsign(f%v, n)
       call chsign(f%w, n)
    end if
    ! don't worry... we'll write this MUCH cleaner in the final version

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine adjoint_forcing


end module user
