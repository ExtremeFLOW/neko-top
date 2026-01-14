module user
  use fld_file_output, only: fld_file_output_t
  use file,  only: file_t
  use csv_file,  only: csv_file_t
  use field_math, only: field_add2s2
  use neko

  implicit none

  ! Data streamer
  type(data_streamer_t) :: dstream
  ! Data streamer (POD modes back)
  type(data_streamer_t) :: dstream_back
  integer :: ipostproc ! frequency of the streaming
  type(coef_t), pointer :: coef
  ! N modes
  integer :: n_modes = 3
  ! modes
  type(field_t), dimension(:), allocatable :: u_list
  type(field_t), dimension(:), allocatable :: v_list
  type(field_t), dimension(:), allocatable :: w_list
  ! writer
  type(fld_file_output_t) :: output
  ! writer to reconstruct
  type(fld_file_output_t) :: output_reconstruct
  ! reader for time coefs
  type(file_t) :: csv_reader
  ! matrix for time coefs
  type(matrix_t)   :: a_time   ! columns: [t, a1, a2, ...]

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

    call json_get(params, "POD.i_stream", ipostproc)
  end subroutine startup

  ! User-defined initialization called just before time loop starts
  subroutine initialize(time)
    type(time_state_t), intent(in) :: time

    integer :: tstep
    character(len=50) :: mess
    integer :: i
    character(len=80) :: str
    type (field_t), pointer :: u, v, w

    ! read postprocessing interval
    write(mess,*) "streaming steps : ", ipostproc
    call neko_log%message(mess)

    ! Initialize the streamer
    coef => neko_user_access%case%fluid%c_Xh
    call dstream%init(coef)


    call output%init(sp, 'yofam', 3 * n_modes)

    allocate(u_list(n_modes))
    allocate(v_list(n_modes))
    allocate(w_list(n_modes))
    do i = 1, n_modes
       write(str, '(A,I0)') "u_fam_", i
       call u_list(i)%init(coef%dof, str)
       call output%fields%assign_to_field(3*(i-1) + 1, u_list(i))

       write(str, '(A,I0)') "v_fam_", i
       call v_list(i)%init(coef%dof, str)
       call output%fields%assign_to_field(3*(i-1) + 2, v_list(i))

       write(str, '(A,I0)') "w_fam_", i
       call w_list(i)%init(coef%dof, str)
       call output%fields%assign_to_field(3*(i-1) + 3, w_list(i))
    end do

    ! Stream the mesh
    call dstream%stream(coef%dof%x)
    call dstream%stream(coef%dof%y)
    call dstream%stream(coef%dof%z)

    ! set up the csv reader
    call csv_reader%init('pod_time_coeffs.csv')

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

    if (NEKO_BCKND_DEVICE .eq. 1) then
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
    integer :: i, j
    integer :: n_lines, nrows, ncols, ierr
    type (field_t), pointer :: u, v, w
    character(len=1024) :: fname
integer :: ctr

    ! Finalize the stream
    call dstream%free()


    ! Give Python time to tear down reader and open the back-writer
    call execute_command_line("sleep 3")

    ! Now start the stream backward
    call dstream_back%init(coef)

    do i = 1, n_modes
      call dstream_back%recieve(u_list(i)%x)
      call dstream_back%recieve(v_list(i)%x)
      call dstream_back%recieve(w_list(i)%x)
    end do

    ! Finalize the stream
    call dstream_back%free()

    ! now we should be certain the csv file is written
    ! Count lines in the CSV
    n_lines = csv_file_count_lines(csv_reader)

    nrows = n_lines - 1

    ! We expect columns = (time + n_modes)
    ncols = 1 + n_modes

    ! init the matrix
    call a_time%init(nrows, ncols)

    call csv_reader%read(a_time)

    ! share with everyone
    call MPI_Allreduce(MPI_IN_PLACE, a_time%x, &
                 a_time%size(), mpi_real_precision, mpi_sum, neko_comm, ierr)

    u => neko_registry%get_field('u')
    v => neko_registry%get_field('v')
    w => neko_registry%get_field('w')

    call output_reconstruct%init(sp, 'reconstruct', 3)
    call output_reconstruct%fields%assign(1, u)
    call output_reconstruct%fields%assign(2, v)
    call output_reconstruct%fields%assign(3, w)
    ! now we can reconstruct
    do i = 1, nrows
    call field_rzero(u)
    call field_rzero(v)
    call field_rzero(w)
      do j = 1, n_modes
        call field_add2s2(u, u_list(j), a_time%x(i, j+1))
        call field_add2s2(v, v_list(j), a_time%x(i, j+1))
        call field_add2s2(w, w_list(j), a_time%x(i, j+1))
      end do
      !if (mod(time%tstep, ipostproc) .eq. 0) then
      call output_reconstruct%sample(a_time%x(i, 1))
      !end if
    end do

    ! and drop the modes
    call output%sample(1.0_rp)


  end subroutine finalize

  function csv_file_count_lines(file_in) result(n)
    class(file_t), intent(in) :: file_in

    integer :: n
    integer :: ierr, file_unit

    open(file = trim(file_in%get_fname()), status = 'old', newunit = file_unit, &
         iostat = ierr)
    if (ierr .ne. 0) then
       call neko_error("Error while opening " // trim(file_in%get_fname()))
    end if
    rewind(file_unit)

    n = 0

    ! Keep reading (ierr = 0) until we reach the end (ierr != 0)
    do
       read (file_unit, *, iostat = ierr)
       if (ierr .ne. 0) exit
       n = n + 1
    end do
    rewind(file_unit)
    close(unit = file_unit)

  end function csv_file_count_lines

end module user
