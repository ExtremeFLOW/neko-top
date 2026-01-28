!> @file POD_state_recover.f90
!! @copyright
!! Copyright (c) 2025, The Neko-TOP Authors
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
!
!> @brief POD-based state recovery using in-situ streaming and reconstruction.
module simulation_POD_state_recover
  use num_types, only: rp, sp, dp
  use case, only: case_t
  use json_file_module, only: json_file
  use json_utils, only: json_get, json_get_or_default
  use field, only: field_t
  use file, only: file_t
  use matrix, only: matrix_t
  use fld_file_output, only: fld_file_output_t
  use coefs, only: coef_t
  use data_streamer, only: data_streamer_t
  use profiler, only: profiler_start_region, profiler_end_region
  use state_recover, only: state_recover_t
  use time_state, only: time_state_t
  use logger, only : neko_log
  use comm, only: neko_comm, mpi_real_precision
  use mpi_f08, only: MPI_Allreduce, MPI_IN_PLACE, MPI_SUM

  use neko_config, only: NEKO_BCKND_DEVICE
  use device, only: DEVICE_TO_HOST, device_memcpy, HOST_TO_DEVICE
  use utils, only: neko_error
  use vector, only: vector_t
  use field_math, only: field_add2s2, field_rzero

  use neko_ctrl_mod, only: ctrl_stream_t, &
       MODE_IDLE, MODE_FORWARD, MODE_ADJOINT, MODE_STOP, &
       PHASE_INIT, PHASE_FWD_RUNNING, PHASE_FWD_DONE, PHASE_ADJ_RUNNING, PHASE_ADJ_DONE
  use, intrinsic :: iso_c_binding, only: c_int, c_double

  implicit none
  private

  !> POD state recovery implementation for forward/adjoint runs.
  type, public, extends(state_recover_t) :: POD_state_recover_t
     private

     logical :: enabled = .true.

     integer :: i_stream
     integer :: n_modes
     integer :: n_flds = 3
     character(len=16) :: dtype = "double"

     type(coef_t), pointer :: coef => null()

     ! SINGLE lifetime streamer (both directions)
     type(data_streamer_t) :: dstream

     ! POD mode storage
     type(field_t), allocatable :: u_modes(:), v_modes(:), w_modes(:)

     ! CSV coeffs
     type(file_t)   :: csv_reader
     type(matrix_t) :: time_coefs
     type(vector_t) :: a_interp

     ! optional output
     type(fld_file_output_t) :: output
     logical :: write_modes = .false.
     logical :: output_reconstruction = .false.
     type(fld_file_output_t) :: recon_output
     character(len=16) :: recon_output_control = "never"
     real(kind=rp) :: recon_output_value = 0.0_rp
     real(kind=rp) :: recon_time_interval = 0.0_rp
     integer :: recon_nsteps = 0

     ! Control state
     type(ctrl_stream_t) :: ctrl

     logical :: have_received_modes = .false.
     logical :: adjoint_started = .false.
     logical :: adj_running_sent = .false.

   contains
     procedure, public, pass(this) :: init => POD_state_recover_init_from_json
     procedure, public, pass(this) :: init_from_json => POD_state_recover_init_from_json
     procedure, public, pass(this) :: init_from_components => POD_state_recover_init_from_components
     procedure, public, pass(this) :: free => POD_state_recover_free
     procedure, public, pass(this) :: reset => POD_state_recover_reset
     procedure, public, pass(this) :: save => POD_state_recover_save
     procedure, public, pass(this) :: restore => POD_state_recover_restore
  end type POD_state_recover_t

contains

  !> Initialize POD state recovery from JSON parameters.
  !! @param[inout] this POD state recovery instance.
  !! @param[inout] neko_case Case data structure.
  !! @param[inout] params JSON parameters.
  subroutine POD_state_recover_init_from_json(this, neko_case, params)
    class(POD_state_recover_t), intent(inout) :: this
    class(case_t), target, intent(inout) :: neko_case
    type(json_file), target, intent(inout) :: params
    integer :: i_stream, n_modes
    logical :: write_modes
    logical :: debug
    logical :: output_reconstruction
    character(len=:), allocatable :: output_precision
    character(len=:), allocatable :: output_control
    real(kind=rp) :: output_value
    integer :: output_prec
    character(len=:), allocatable :: dtype

    call json_get(params, "i_stream", i_stream)
    call json_get(params, "n_modes", n_modes)
    call json_get_or_default(params, "dtype", dtype, "double")
    call json_get_or_default(params, "write_modes", write_modes, .false.)
    call json_get_or_default(params, "debug", debug, .false.)
    call json_get_or_default(params, "output_reconstruction", &
         output_reconstruction, .false.)

    call this%init_from_components(neko_case, i_stream, n_modes, debug)
    this%dtype = adjustl(dtype)
    this%write_modes = write_modes
    this%output_reconstruction = output_reconstruction

    select case (trim(this%dtype))
    case ("single", "SINGLE", "Single")
       if (rp .ne. sp) call neko_error("POD dtype single but code not single precision.")
    case ("double", "DOUBLE", "Double")
       if (rp .ne. dp) call neko_error("POD dtype double but code not double precision.")
    case default
       call neko_error("Unsupported POD dtype: " // trim(this%dtype))
    end select

    if (this%output_reconstruction) then
       call json_get_or_default(neko_case%params, 'case.output_precision', &
            output_precision, 'single')
       if (trim(output_precision) .eq. 'double') then
          output_prec = dp
       else
          output_prec = sp
       end if

       call this%recon_output%init(output_prec, 'pod_reconstruction', &
            this%n_flds, path = trim(neko_case%output_directory))
       call this%recon_output%fields%assign_to_field(1, neko_case%fluid%u)
       call this%recon_output%fields%assign_to_field(2, neko_case%fluid%v)
       call this%recon_output%fields%assign_to_field(3, neko_case%fluid%w)

       call json_get_or_default(neko_case%params, 'case.fluid.output_control', &
            output_control, 'org')
       this%recon_output_control = trim(lower_string(output_control))

       select case (this%recon_output_control)
       case ('org')
          call json_get(neko_case%params, 'case.nsamples', output_value)
          this%recon_output_control = 'nsamples'
          this%recon_output_value = output_value
          if (output_value .gt. 0.0_rp) then
             this%recon_time_interval = (neko_case%time%end_time - &
                  neko_case%time%start_time) / output_value
          end if
       case ('nsamples')
          call json_get(neko_case%params, 'case.fluid.output_value', output_value)
          this%recon_output_value = output_value
          if (output_value .gt. 0.0_rp) then
             this%recon_time_interval = (neko_case%time%end_time - &
                  neko_case%time%start_time) / output_value
          end if
       case ('simulationtime')
          call json_get(neko_case%params, 'case.fluid.output_value', output_value)
          this%recon_output_value = output_value
          this%recon_time_interval = output_value
       case ('tsteps')
          call json_get(neko_case%params, 'case.fluid.output_value', output_value)
          this%recon_output_value = output_value
          this%recon_nsteps = int(output_value)
       case ('never')
          call json_get_or_default(neko_case%params, 'case.fluid.output_value', &
               output_value, 0.0_rp)
          this%recon_output_value = output_value
       case default
          call neko_error('Unsupported output_control for reconstruction: ' // &
               trim(output_control))
       end select
    end if
  end subroutine POD_state_recover_init_from_json


  !> Initialize POD state recovery from explicit components.
  !! @param[inout] this POD state recovery instance.
  !! @param[inout] neko_case Case data structure.
  !! @param[in] i_stream Snapshot stride.
  !! @param[in] n_modes Number of POD modes to keep.
  !! @param[in] debug Optional debug flag.
  subroutine POD_state_recover_init_from_components(this, neko_case, i_stream, n_modes, debug)
    class(POD_state_recover_t), intent(inout), target :: this
    class(case_t), target, intent(inout) :: neko_case
    integer, intent(in) :: i_stream, n_modes
    logical, intent(in), optional :: debug
    integer :: i
    character(len=80) :: str

    this%i_stream = i_stream
    this%n_modes  = n_modes
    this%coef => neko_case%fluid%c_Xh

    allocate(this%u_modes(this%n_modes))
    allocate(this%v_modes(this%n_modes))
    allocate(this%w_modes(this%n_modes))

    call this%output%init(sp, 'POD_modes', this%n_flds * this%n_modes)
    do i = 1, this%n_modes
       write(str, '(A,I0)') "u_mode_", i
       call this%u_modes(i)%init(this%coef%dof, trim(str))
       call this%output%fields%assign(this%n_flds*(i-1) + 1, this%u_modes(i))

       write(str, '(A,I0)') "v_mode_", i
       call this%v_modes(i)%init(this%coef%dof, trim(str))
       call this%output%fields%assign(this%n_flds*(i-1) + 2, this%v_modes(i))

       write(str, '(A,I0)') "w_mode_", i
       call this%w_modes(i)%init(this%coef%dof, trim(str))
       call this%output%fields%assign(this%n_flds*(i-1) + 3, this%w_modes(i))
    end do

    call this%a_interp%init(this%n_modes)
    call this%csv_reader%init('pod_time_coeffs.csv')

    ! Single lifetime streamer: init once, keep open
    call this%dstream%init(this%coef)

    ! Send mesh once
    call this%dstream%stream(this%coef%dof%x)
    call this%dstream%stream(this%coef%dof%y)
    call this%dstream%stream(this%coef%dof%z)

    ! Stream initial condition once (t=0 snapshot)
    block
      type(field_t), pointer :: u, v, w
      integer :: n

      u => neko_case%fluid%u
      v => neko_case%fluid%v
      w => neko_case%fluid%w
      n = u%dof%size()

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_memcpy(u%x, u%x_d, n, DEVICE_TO_HOST, sync=.true.)
         call device_memcpy(v%x, v%x_d, n, DEVICE_TO_HOST, sync=.true.)
         call device_memcpy(w%x, w%x_d, n, DEVICE_TO_HOST, sync=.true.)
      end if

      call this%dstream%stream(u%x)
      call this%dstream%stream(v%x)
      call this%dstream%stream(w%x)
    end block

    ! Control init (use neko_comm%mpi_val – your working pattern)
    if (present(debug)) then
       this%ctrl%debug = debug
    end if
    call this%ctrl%init(int(neko_comm%mpi_val, c_int))

    ! Fire an init tick (python might miss it; harmless)
    call this%ctrl%put(MODE_IDLE, PHASE_INIT, 0_c_int, 0.0_c_double)

  end subroutine POD_state_recover_init_from_components


  !> Release POD state recovery resources.
  !! @param[inout] this POD state recovery instance.
  subroutine POD_state_recover_free(this)
    class(POD_state_recover_t), intent(inout) :: this
    integer :: i

    if (allocated(this%u_modes)) then
       do i = 1, size(this%u_modes); call this%u_modes(i)%free(); end do
       deallocate(this%u_modes)
    end if
    if (allocated(this%v_modes)) then
       do i = 1, size(this%v_modes); call this%v_modes(i)%free(); end do
       deallocate(this%v_modes)
    end if
    if (allocated(this%w_modes)) then
       do i = 1, size(this%w_modes); call this%w_modes(i)%free(); end do
       deallocate(this%w_modes)
    end if

    call this%dstream%free()
    call this%csv_reader%free()
    call this%a_interp%free()
    nullify(this%coef)

    call this%ctrl%free()

    this%enabled = .false.
  end subroutine POD_state_recover_free


  !> Reset POD state recovery to initial control state.
  !! @param[inout] this POD state recovery instance.
  subroutine POD_state_recover_reset(this)
    class(POD_state_recover_t), intent(inout) :: this
    integer :: i

    if (.not. this%enabled) return

    ! If we were in adjoint previously, emit ADJ_DONE once on reset
    if (this%ctrl%inited .and. this%adjoint_started) then
       call this%ctrl%put(MODE_ADJOINT, PHASE_ADJ_DONE, 0_c_int, 0.0_c_double)
    end if

    this%have_received_modes = .false.
    this%adjoint_started = .false.
    this%adj_running_sent = .false.

    call this%set_n_timesteps(0)

    do i = 1, this%n_modes
       call field_rzero(this%u_modes(i))
       call field_rzero(this%v_modes(i))
       call field_rzero(this%w_modes(i))
    end do

  end subroutine POD_state_recover_reset


  !> Stream forward state for POD updates.
  !! @param[inout] this POD state recovery instance.
  !! @param[inout] neko_case Case data structure.
  !! @param[in] time Current time state.
  subroutine POD_state_recover_save(this, neko_case, time)
    class(POD_state_recover_t), intent(inout) :: this
    class(case_t), intent(inout) :: neko_case
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: u, v, w
    integer :: n

    if (.not. this%enabled) return
    if (mod(time%tstep, this%i_stream) .ne. 0) return

    call neko_log%message(" ")
    call neko_log%message("----------------")
    call neko_log%message("Streaming fields")
    call neko_log%message("----------------")
    call neko_log%message(" ")

    u => neko_case%fluid%u
    v => neko_case%fluid%v
    w => neko_case%fluid%w
    n = u%dof%size()

    call profiler_start_region("POD save")

    if (this%ctrl%inited) then
       call this%ctrl%put(MODE_FORWARD, PHASE_FWD_RUNNING, &
            int(time%tstep, c_int), real(time%t, c_double))
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(u%x, u%x_d, n, DEVICE_TO_HOST, sync=.true.)
       call device_memcpy(v%x, v%x_d, n, DEVICE_TO_HOST, sync=.true.)
       call device_memcpy(w%x, w%x_d, n, DEVICE_TO_HOST, sync=.true.)
    end if

    call this%dstream%stream(u%x)
    call this%dstream%stream(v%x)
    call this%dstream%stream(w%x)

    call profiler_end_region("POD save")
  end subroutine POD_state_recover_save


  !> Reconstruct and restore state from POD during adjoint.
  !! @param[inout] this POD state recovery instance.
  !! @param[inout] neko_case Case data structure.
  !! @param[in] time Target time state.
  subroutine POD_state_recover_restore(this, neko_case, time)
    class(POD_state_recover_t), intent(inout) :: this
    class(case_t), target, intent(inout) :: neko_case
    type(time_state_t), intent(in) :: time
    integer :: i, ierr, n_lines, nrows, ncols, n
    integer(c_int) :: mode_cmd, phase_cmd
    type(time_state_t) :: time_out

    if (.not. this%enabled) return

    ! First restore() call is the phase boundary forward->adjoint
    if (.not. this%have_received_modes) then
       call profiler_start_region("POD recieve modes")

       if (this%ctrl%inited) then
          call this%ctrl%put(MODE_FORWARD, PHASE_FWD_DONE, &
               int(time%tstep, c_int), real(time%t, c_double))

          mode_cmd  = MODE_FORWARD
          phase_cmd = PHASE_FWD_DONE

          ! BLOCK until Python says "go adjoint"
          call this%ctrl%wait_cmd(mode_cmd, phase_cmd)

          if (mode_cmd /= MODE_ADJOINT) then
             call neko_error("Expected MODE_ADJOINT from Python at forward->adjoint boundary.")
          end if
       end if

       this%adjoint_started = .true.

       ! Receive modes (Python streams after sending cmd)
       do i = 1, this%n_modes
          call this%dstream%recieve(this%u_modes(i)%x)
          call this%dstream%recieve(this%v_modes(i)%x)
          call this%dstream%recieve(this%w_modes(i)%x)
       end do

       ! Move modes back to GPU
       n = this%u_modes(1)%dof%size()
       if (NEKO_BCKND_DEVICE .eq. 1) then
       do i = 1, this%n_modes
          call device_memcpy(this%u_modes(i)%x, this%u_modes(i)%x_d, n, &
               HOST_TO_DEVICE, sync=.true.)
          call device_memcpy(this%v_modes(i)%x, this%v_modes(i)%x_d, n, &
               HOST_TO_DEVICE, sync=.true.)
          call device_memcpy(this%w_modes(i)%x, this%w_modes(i)%x_d, n, &
               HOST_TO_DEVICE, sync=.true.)
       end do
       end if

       if (this%write_modes) then
          call this%output%sample(0.0_rp)
       end if

       ! Read CSV once
       n_lines = csv_file_count_lines(this%csv_reader)
       nrows = n_lines - 1
       ncols = 1 + this%n_modes
       call this%time_coefs%init(nrows, ncols)
       call this%csv_reader%read(this%time_coefs)

       call MPI_Allreduce(MPI_IN_PLACE, this%time_coefs%x, &
            this%time_coefs%size(), mpi_real_precision, MPI_SUM, neko_comm, ierr)

       this%have_received_modes = .true.
       call profiler_end_region("POD recieve modes")
    end if

    ! Emit ADJ_RUNNING only once (avoid flooding SST)
    if (this%ctrl%inited .and. .not. this%adj_running_sent) then
       call this%ctrl%put(MODE_ADJOINT, PHASE_ADJ_RUNNING, &
            int(time%tstep, c_int), real(time%t, c_double))
       this%adj_running_sent = .true.
    end if

    call profiler_start_region("POD restore")
    call interpolate_time_coeffs_vec(this%a_interp, this%time_coefs, time%t)
    call reconstruct_from_coeffs(this, neko_case, this%a_interp)
    if (this%output_reconstruction) then
       if (recon_should_output(this, time, time_out)) then
          call this%recon_output%sample(time_out%t)
       end if
    end if
    call profiler_end_region("POD restore")
  end subroutine POD_state_recover_restore

  logical function recon_should_output(this, time, time_out)
    class(POD_state_recover_t), intent(in) :: this
    type(time_state_t), intent(in) :: time
    type(time_state_t), intent(out) :: time_out
    real(kind=rp) :: t_rel, tol, interval

    recon_should_output = .false.
    time_out = time
    time_out%t = time%start_time + real(time%tstep, rp) * time%dt

    select case (this%recon_output_control)
    case ('tsteps')
       if (this%recon_nsteps .gt. 0) then
          recon_should_output = mod(time%tstep, this%recon_nsteps) .eq. 0
       end if
    case ('simulationtime', 'nsamples')
       interval = this%recon_time_interval
       if (interval .gt. 0.0_rp) then
          t_rel = time_out%t - time%start_time
          tol = 0.1_rp * abs(time%dt)
          recon_should_output = abs(t_rel - nint(t_rel / interval) * interval) &
               .le. tol
       end if
    case ('never')
       recon_should_output = .false.
    case default
       recon_should_output = .false.
    end select
  end function recon_should_output

  pure function lower_string(str) result(out)
    character(len=*), intent(in) :: str
    character(len=len(str)) :: out
    integer :: i, c

    do i = 1, len(str)
       c = iachar(str(i:i))
       if (c >= iachar('A') .and. c <= iachar('Z')) then
          out(i:i) = achar(c + 32)
       else
          out(i:i) = str(i:i)
       end if
    end do
  end function lower_string


  function csv_file_count_lines(file_in) result(n)
    class(file_t), intent(in) :: file_in
    integer :: n
    integer :: ierr, file_unit

    open(file = trim(file_in%get_fname()), status='old', newunit=file_unit, iostat=ierr)
    if (ierr .ne. 0) call neko_error("Error opening " // trim(file_in%get_fname()))
    rewind(file_unit)

    n = 0
    do
       read(file_unit, *, iostat=ierr)
       if (ierr .ne. 0) exit
       n = n + 1
    end do

    close(unit=file_unit)
  end function csv_file_count_lines


  subroutine find_bracket_time(i0, i1, time_coefs, tq)
    integer, intent(out)       :: i0, i1
    type(matrix_t), intent(in) :: time_coefs
    real(kind=rp), intent(in)  :: tq
    integer :: lo, hi, mid
    real(kind=rp) :: tm

    lo = 1
    hi = time_coefs%get_nrows()

    do while (hi - lo > 1)
       mid = (lo + hi) / 2
       tm = time_coefs%x(mid, 1)
       if (tm <= tq) then
          lo = mid
       else
          hi = mid
       end if
    end do

    i0 = lo
    i1 = lo + 1
  end subroutine find_bracket_time


  subroutine interpolate_time_coeffs_vec(a_out, time_coefs, tq)
    type(vector_t), intent(inout) :: a_out
    type(matrix_t), intent(in)    :: time_coefs
    real(kind=rp), intent(in)     :: tq
    integer :: nrows, ncols, j, i0, i1
    real(kind=rp) :: t0, t1, w

    nrows = time_coefs%get_nrows()
    ncols = time_coefs%get_ncols()

    if (ncols < 2) call neko_error("time_coefs must have (t + coeffs) columns")
    if (size(a_out%x) /= ncols-1) call neko_error("a_out wrong size")

    if (tq <= time_coefs%x(1,1)) then
       do j = 1, ncols-1
          a_out%x(j) = time_coefs%x(1, j+1)
       end do
       return
    end if

    if (tq >= time_coefs%x(nrows,1)) then
       do j = 1, ncols-1
          a_out%x(j) = time_coefs%x(nrows, j+1)
       end do
       return
    end if

    call find_bracket_time(i0, i1, time_coefs, tq)
    t0 = time_coefs%x(i0,1)
    t1 = time_coefs%x(i1,1)

    if (abs(t1 - t0) < tiny(1.0_rp)) then
       w = 0.0_rp
    else
       w = (tq - t0) / (t1 - t0)
    end if

    do j = 1, ncols-1
       a_out%x(j) = (1.0_rp - w) * time_coefs%x(i0, j+1) + w * time_coefs%x(i1, j+1)
    end do
  end subroutine interpolate_time_coeffs_vec


  subroutine reconstruct_from_coeffs(this, neko_case, a)
    class(POD_state_recover_t), intent(inout) :: this
    class(case_t), intent(inout)              :: neko_case
    type(vector_t), intent(in)                :: a
    integer :: j
    type(field_t), pointer :: u, v, w

    u => neko_case%fluid%u
    v => neko_case%fluid%v
    w => neko_case%fluid%w

    call field_rzero(u)
    call field_rzero(v)
    call field_rzero(w)

    do j = 1, this%n_modes
       call field_add2s2(u, this%u_modes(j), a%x(j))
       call field_add2s2(v, this%v_modes(j), a%x(j))
       call field_add2s2(w, this%w_modes(j), a%x(j))
    end do
  end subroutine reconstruct_from_coeffs

end module simulation_POD_state_recover
