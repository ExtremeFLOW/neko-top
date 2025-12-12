module user
  use neko
  implicit none

  ! Data streamer
  type(data_streamer_t) :: dstream
  integer :: ipostproc ! frequency of the streaming
  type(coef_t), pointer :: coef

contains

  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%compute => compute
    user%startup => startup
    user%initialize => initialize
    user%finalize => finalize
  end subroutine user_setup

  subroutine startup(params)
    type(json_file), intent(inout) :: params

    call json_get(params, "case.istream", ipostproc)
  end subroutine startup

  ! User-defined initialization called just before time loop starts
  subroutine initialize(time)
    type(time_state_t), intent(in) :: time

    integer :: tstep
    character(len=50) :: mess

    ! read postprocessing interval
    write(mess,*) "streaming steps : ", ipostproc
    call neko_log%message(mess)

    ! Initialize the streamer
    coef => neko_user_access%case%fluid%c_Xh
    call dstream%init(coef)

    ! Stream the mesh
    call dstream%stream(coef%dof%x)
    call dstream%stream(coef%dof%y)
    call dstream%stream(coef%dof%z)

  end subroutine initialize

  ! User-defined routine called at the end of every time step
  subroutine compute(time)
    type(time_state_t), intent(in) :: time

    integer :: ntot, i, n
    real(kind=rp) :: ekin, enst
    type(dofmap_t), pointer :: dof
    type (field_t), pointer :: u, v, w, p

    u => neko_registry%get_field('u')
    v => neko_registry%get_field('v')
    w => neko_registry%get_field('w')

    if (mod(time%tstep, ipostproc) .ne. 0) return

    n = u%dof%size()

    if (NEKO_BCKND_DEVICE .eq .1) then
       ! Ensure host buffers are up to date before streaming from a GPU run.
       call device_memcpy(u%x, u%x_d, n, DEVICE_TO_HOST, sync=.true.)
       call device_memcpy(v%x, v%x_d, n, DEVICE_TO_HOST, sync=.true.)
       call device_memcpy(w%x, w%x_d, n, DEVICE_TO_HOST, sync=.true.)
    end if

    ! Stream the data
    call dstream%stream(u%x)
    call dstream%stream(v%x)
    call dstream%stream(w%x)

  end subroutine compute

  ! User-defined finalization routine called at the end of the simulation
  subroutine finalize(time)
    type(time_state_t), intent(in) :: time

    ! Finalize the stream
    call dstream%free()

  end subroutine finalize

end module user
