!> @file POD_state_recover.f90
!! @copyright
!! Copyright (c) 2025-2026, The Neko-TOP Authors
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
  use field_output, only: field_output_t
  use coefs, only: coef_t
  use data_streamer, only: data_streamer_t
  use profiler, only: profiler_start_region, profiler_end_region
  use state_recover, only: state_recover_t
  use time_state, only: time_state_t
  use simulation, only: simulation_step
  use time_step_controller, only: time_step_controller_t
  use logger, only : neko_log
  use comm, only: neko_comm, mpi_real_precision
  use mpi_f08, only: MPI_Allreduce, MPI_IN_PLACE, MPI_SUM, MPI_WTIME

  use neko_config, only: NEKO_BCKND_DEVICE
  use device, only: DEVICE_TO_HOST, device_memcpy, HOST_TO_DEVICE
  use utils, only: neko_error
  use vector, only: vector_t
  use field_math, only: field_add2s2, field_rzero

  use neko_ctrl_mod, only: ctrl_stream_t, &
       MODE_IDLE, MODE_FORWARD, MODE_ADJOINT, MODE_STOP, &
       PHASE_INIT, PHASE_FWD_RUNNING, PHASE_FWD_DONE, &
       PHASE_ADJ_RUNNING, PHASE_ADJ_DONE
  use, intrinsic :: iso_c_binding, only: c_int, c_double

  implicit none
  private

  ! Control protocol used by save()/restore() and the Python POD driver:
  ! - MODE_IDLE + PHASE_INIT is a startup tick sent after ctrl init.
  ! - MODE_FORWARD + PHASE_FWD_RUNNING is emitted before each streamed forward
  !   snapshot batch.
  ! - MODE_FORWARD + PHASE_FWD_DONE marks the forward-to-adjoint boundary.
  !   Neko then blocks until Python has updated the POD basis and responds with
  !   MODE_ADJOINT + PHASE_ADJ_RUNNING.
  ! - MODE_ADJOINT + PHASE_ADJ_RUNNING means restore() may reconstruct states
  !   from the received POD basis and time coefficients.
  ! - MODE_ADJOINT + PHASE_ADJ_DONE is emitted on reset after an adjoint pass
  !   so that Python can clear its old POD state before the next forward pass.
  ! - MODE_STOP is emitted during free() to terminate the Python control loop.
  !> POD state recovery implementation for forward/adjoint runs.
  type, public, extends(state_recover_t) :: POD_state_recover_t
     private

     logical :: enabled = .true.

     integer :: i_stream
     integer :: n_modes
     integer :: n_flds = 3
     integer :: pod_tstep = 0
     logical :: include_scalar = .false.
     character(len=16) :: dtype = "double"
     real(kind=rp) :: transience_time = 0.0_rp
     real(kind=rp) :: time_shift = 0.0_rp

     type(coef_t), pointer :: coef => null()

     ! SINGLE lifetime streamer (both directions)
     type(data_streamer_t) :: dstream

     ! POD mode storage
     type(field_t), allocatable :: u_modes(:), v_modes(:), w_modes(:)
     type(field_t), allocatable :: s_modes(:)
     type(field_t), pointer :: s => null()

     ! CSV coeffs
     type(file_t)   :: csv_reader
     type(matrix_t) :: time_coefs
     type(vector_t) :: a_interp

     ! optional mode output configuration
     integer :: mode_output_precision = sp
     character(len=32) :: mode_output_format = 'fld'
     character(len=80) :: mode_file_name = 'POD_modes'
     logical :: write_modes = .false.
     logical :: output_reconstruction = .false.
     type(field_output_t) :: recon_output
     integer :: recon_output_precision = sp
     character(len=32) :: recon_output_format = 'fld'
     character(len=80) :: recon_file_name = 'pod_reconstruction'
     character(len=16) :: recon_output_control = "never"
     real(kind=rp) :: recon_output_value = 0.0_rp
     real(kind=rp) :: recon_time_interval = 0.0_rp
     integer :: recon_nsteps = 0

     ! Control state
     type(ctrl_stream_t) :: ctrl

     logical :: have_received_modes = .false.
     logical :: adjoint_started = .false.
     logical :: adj_running_sent = .false.
     logical :: transience_applied = .true.

   contains
     procedure, public, pass(this) :: init => POD_state_recover_init_from_json
     procedure, public, pass(this) :: init_from_json => &
          POD_state_recover_init_from_json
     procedure, public, pass(this) :: init_from_components => &
          POD_state_recover_init_from_components
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
    real(kind=rp) :: transience_time
    character(len=:), allocatable :: recon_output_precision_str
    character(len=:), allocatable :: output_control
    real(kind=rp) :: output_value
    character(len=:), allocatable :: dtype
    character(len=:), allocatable :: mode_output_precision_str
    character(len=:), allocatable :: mode_output_format
    character(len=:), allocatable :: mode_file_name
    character(len=:), allocatable :: recon_output_format
    character(len=:), allocatable :: recon_file_name
    integer :: mode_output_precision
    integer :: recon_output_precision

    call json_get(params, "i_stream", i_stream)
    call json_get(params, "n_modes", n_modes)
    call json_get_or_default(params, "transience_time", transience_time, &
         0.0_rp)
    if (transience_time .lt. 0.0_rp) then
       call neko_error("transience_time must be non-negative.")
    end if
    this%transience_time = transience_time
    call json_get_or_default(params, "dtype", dtype, "double")
    call json_get_or_default(params, "write_modes", write_modes, .false.)
    call json_get_or_default(params, "debug", debug, .false.)
    call json_get_or_default(params, "output_reconstruction", &
         output_reconstruction, .false.)
    call json_get_or_default(params, "output_precision", &
         mode_output_precision_str, "sp")
    call json_get_or_default(params, "output_format", &
         mode_output_format, "fld")
    call json_get_or_default(params, "output_file_name", mode_file_name, &
         "POD_modes")

    mode_output_precision = precision_from_string(mode_output_precision_str, &
         'state_recovery.output_precision')

    recon_output_precision = mode_output_precision
    recon_output_precision_str = trim(mode_output_precision_str)
    recon_output_format = trim(mode_output_format)
    recon_file_name = 'pod_reconstruction'
    output_control = 'never'
    output_value = 0.0_rp

    if (output_reconstruction) then
       if (.not. ('output_precision' .in. params)) then
          call json_get_or_default(neko_case%params, 'case.output_precision', &
               recon_output_precision_str, 'single')
          recon_output_precision = precision_from_string( &
               recon_output_precision_str, &
               'state_recovery.reconstruction_output_precision')
       end if
       if ('reconstruction_output_precision' .in. params) then
          call json_get(params, 'reconstruction_output_precision', &
               recon_output_precision_str)
          recon_output_precision = precision_from_string( &
               recon_output_precision_str, &
               'state_recovery.reconstruction_output_precision')
       end if
       if ('reconstruction_output_format' .in. params) then
          call json_get(params, 'reconstruction_output_format', &
               recon_output_format)
       end if
       call json_get_or_default(params, 'reconstruction_output_file_name', &
            recon_file_name, 'pod_reconstruction')

       call json_get_or_default(neko_case%params, 'case.fluid.output_control', &
            output_control, 'org')
       select case (trim(lower_string(output_control)))
       case ('org')
          call json_get(neko_case%params, 'case.nsamples', output_value)
       case ('nsamples', 'simulationtime', 'tsteps')
          call json_get(neko_case%params, 'case.fluid.output_value', &
               output_value)
       case ('never')
          call json_get_or_default(neko_case%params, &
               'case.fluid.output_value', output_value, 0.0_rp)
       case default
          call json_get_or_default(neko_case%params, &
               'case.fluid.output_value', output_value, 0.0_rp)
       end select
    end if

    call this%init_from_components(neko_case, i_stream, n_modes, dtype, &
         write_modes, output_reconstruction, output_control, output_value, &
         debug, mode_output_precision, &
         mode_output_format, mode_file_name, recon_output_precision, &
         recon_output_format, recon_file_name)
  end subroutine POD_state_recover_init_from_json


  !> Initialize POD state recovery from explicit components.
  !! @param[inout] this POD state recovery instance.
  !! @param[inout] neko_case Case data structure.
  !! @param[in] i_stream Snapshot stride.
  !! @param[in] n_modes Number of POD modes to keep.
  !! @param[in] dtype POD floating-point precision.
  !! @param[in] write_modes Whether Python should persist POD modes.
  !! @param[in] output_reconstruction Whether to output reconstructions.
  !! @param[in] output_control Reconstruction output cadence control.
  !! @param[in] output_value Reconstruction cadence value.
  !! @param[in] debug Optional debug flag.
  !! @param[in] mode_output_precision Optional POD mode output precision.
  !! @param[in] mode_output_format Optional POD mode output format.
  !! @param[in] mode_file_name Optional POD mode output base name.
  !! @param[in] recon_output_precision Optional reconstruction precision.
  !! @param[in] recon_output_format Optional reconstruction output format.
  !! @param[in] recon_file_name Optional reconstruction output base name.
  subroutine POD_state_recover_init_from_components(this, neko_case, &
       i_stream, n_modes, dtype, write_modes, output_reconstruction, &
       output_control, output_value, debug, &
       mode_output_precision, mode_output_format, mode_file_name, &
       recon_output_precision, recon_output_format, recon_file_name)
    class(POD_state_recover_t), intent(inout), target :: this
    class(case_t), target, intent(inout) :: neko_case
    integer, intent(in) :: i_stream, n_modes
    character(len=*), intent(in) :: dtype
    logical, intent(in) :: write_modes
    logical, intent(in) :: output_reconstruction
    character(len=*), intent(in) :: output_control
    real(kind=rp), intent(in) :: output_value
    logical, intent(in), optional :: debug
    integer, intent(in), optional :: mode_output_precision
    character(len=*), intent(in), optional :: mode_output_format
    character(len=*), intent(in), optional :: mode_file_name
    integer, intent(in), optional :: recon_output_precision
    character(len=*), intent(in), optional :: recon_output_format
    character(len=*), intent(in), optional :: recon_file_name
    integer :: i
    character(len=80) :: str

    this%enabled = .true.
    this%i_stream = i_stream
    this%n_modes = n_modes
    if (this%transience_time .lt. 0.0_rp) then
       call neko_error("transience_time must be non-negative.")
    end if
    this%dtype = adjustl(dtype)
    this%write_modes = write_modes
    this%output_reconstruction = output_reconstruction
    this%recon_output_control = 'never'
    this%recon_output_value = 0.0_rp
    this%recon_time_interval = 0.0_rp
    this%recon_nsteps = 0
    this%recon_output_precision = sp
    this%recon_output_format = 'fld'
    this%recon_file_name = 'pod_reconstruction'
    this%mode_output_precision = sp
    this%mode_output_format = 'fld'
    this%mode_file_name = 'POD_modes'
    this%n_flds = 3
    this%pod_tstep = 0
    this%time_shift = 0.0_rp
    this%transience_applied = this%transience_time .le. 0.0_rp
    this%coef => neko_case%fluid%c_Xh

    this%include_scalar = allocated(neko_case%scalars)

    if (present(mode_output_precision)) then
       if (mode_output_precision .ne. sp .and. &
            mode_output_precision .ne. dp) then
          call neko_error('mode_output_precision must be either sp or dp')
       end if
       this%mode_output_precision = mode_output_precision
    end if
    if (present(mode_output_format)) then
       this%mode_output_format = trim(mode_output_format)
    end if
    if (present(mode_file_name)) then
       this%mode_file_name = trim(mode_file_name)
    end if
    if (present(recon_output_precision)) then
      if (recon_output_precision .ne. sp .and. &
           recon_output_precision .ne. dp) then
         call neko_error('recon_output_precision must be either sp or dp')
      end if
      this%recon_output_precision = recon_output_precision
    end if
    if (present(recon_output_format)) then
       this%recon_output_format = trim(recon_output_format)
    end if
    if (present(recon_file_name)) then
       this%recon_file_name = trim(recon_file_name)
    end if

    select case (trim(lower_string(this%dtype)))
    case ('single')
       if (rp .ne. sp) then
          call neko_error("POD dtype single but code not single precision.")
       end if
    case ('double')
       if (rp .ne. dp) then
          call neko_error("POD dtype double but code not double precision.")
       end if
    case default
       call neko_error("Unsupported POD dtype: " // trim(this%dtype))
    end select

    if (this%include_scalar) then
       if (size(neko_case%scalars%scalar_fields) .ne. 1) then
          call neko_error('POD scalar recovery currently supports ' // &
               'exactly one scalar.')
       end if
       this%s => neko_case%scalars%scalar_fields(1)%scalar%s
       this%n_flds = 4
    else
       nullify(this%s)
    end if

    allocate(this%u_modes(this%n_modes))
    allocate(this%v_modes(this%n_modes))
    allocate(this%w_modes(this%n_modes))
    if (this%include_scalar) allocate(this%s_modes(this%n_modes))

    do i = 1, this%n_modes
       write(str, '(A,I0)') "u_mode_", i
       call this%u_modes(i)%init(this%coef%dof, trim(str))

       write(str, '(A,I0)') "v_mode_", i
       call this%v_modes(i)%init(this%coef%dof, trim(str))

       write(str, '(A,I0)') "w_mode_", i
       call this%w_modes(i)%init(this%coef%dof, trim(str))

       if (this%include_scalar) then
          write(str, '(A,I0)') "s_mode_", i
          call this%s_modes(i)%init(this%s%dof, trim(str))
       end if
    end do

    if (this%output_reconstruction) then
       call this%recon_output%init(trim(this%recon_file_name), this%n_flds, &
            precision = this%recon_output_precision, &
            path = trim(neko_case%output_directory), &
            format = trim(this%recon_output_format))
       call this%recon_output%fields%assign_to_field(1, neko_case%fluid%u)
       call this%recon_output%fields%assign_to_field(2, neko_case%fluid%v)
       call this%recon_output%fields%assign_to_field(3, neko_case%fluid%w)
       if (this%include_scalar) then
          call this%recon_output%fields%assign_to_field(4, this%s)
       end if

       this%recon_output_control = trim(lower_string(output_control))
       this%recon_output_value = output_value

       select case (this%recon_output_control)
       case ('org')
          this%recon_output_control = 'nsamples'
          if (output_value .gt. 0.0_rp) then
             this%recon_time_interval = (neko_case%time%end_time - &
                  neko_case%time%start_time) / output_value
          end if
       case ('nsamples')
          if (output_value .gt. 0.0_rp) then
             this%recon_time_interval = (neko_case%time%end_time - &
                  neko_case%time%start_time) / output_value
          end if
       case ('simulationtime')
          this%recon_time_interval = output_value
       case ('tsteps')
          this%recon_nsteps = int(output_value)
       case ('never')
       case default
          call neko_error('Unsupported output_control for reconstruction: ' // &
               trim(output_control))
       end select
    end if

    call this%a_interp%init(this%n_modes)
    call this%csv_reader%init('pod_time_coeffs.csv')

    ! Single lifetime streamer: init once, keep open
    call this%dstream%init(this%coef)

    ! Send mesh once
    call this%dstream%stream(this%coef%dof%x)
    call this%dstream%stream(this%coef%dof%y)
    call this%dstream%stream(this%coef%dof%z)

    ! Stream the initial condition only when sampling starts immediately.
    if (this%transience_applied) then
       call POD_state_recover_stream_fields(this, neko_case)
    end if

    ! Control init for root-to-root MPI coordination with Python.
    if (present(debug)) then
       this%ctrl%debug = debug
    end if
    call this%ctrl%init()

    ! Fire an init tick (python might miss it; harmless)
    call this%ctrl%send(MODE_IDLE, PHASE_INIT, 0_c_int, 0.0_c_double)
    call this%set_n_timesteps(0)

  end subroutine POD_state_recover_init_from_components


  !> Release POD state recovery resources.
  !! @param[inout] this POD state recovery instance.
  subroutine POD_state_recover_free(this)
    class(POD_state_recover_t), intent(inout) :: this
    integer :: i

    if (this%ctrl%inited) then
       call this%ctrl%send(MODE_STOP, PHASE_ADJ_DONE, 0_c_int, &
            0.0_c_double)
    end if

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
    if (allocated(this%s_modes)) then
       do i = 1, size(this%s_modes); call this%s_modes(i)%free(); end do
       deallocate(this%s_modes)
    end if

    call this%recon_output%free()
    call this%dstream%free()
    call this%csv_reader%free()
    call this%a_interp%free()
    nullify(this%coef)

    call this%ctrl%free()

    this%include_scalar = .false.
    this%pod_tstep = 0
    this%time_shift = 0.0_rp
    this%transience_time = 0.0_rp
    this%transience_applied = .true.
    this%mode_output_precision = sp
    this%mode_output_format = 'fld'
    this%mode_file_name = 'POD_modes'
    this%recon_output_precision = sp
    this%recon_output_format = 'fld'
    this%recon_file_name = 'pod_reconstruction'
    nullify(this%s)
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
       call this%ctrl%send(MODE_ADJOINT, PHASE_ADJ_DONE, 0_c_int, &
            0.0_c_double)
    end if

    this%have_received_modes = .false.
    this%adjoint_started = .false.
    this%adj_running_sent = .false.
    this%pod_tstep = 0

    call this%set_n_timesteps(0)

    do i = 1, this%n_modes
       call field_rzero(this%u_modes(i))
       call field_rzero(this%v_modes(i))
       call field_rzero(this%w_modes(i))
       if (this%include_scalar) call field_rzero(this%s_modes(i))
    end do

  end subroutine POD_state_recover_reset

  subroutine POD_state_recover_apply_transience(this, neko_case)
    class(POD_state_recover_t), intent(inout) :: this
    class(case_t), intent(inout) :: neko_case
    type(time_step_controller_t) :: dt_controller
    real(kind=dp) :: loop_start
    real(kind=rp) :: target_time

    if (this%transience_applied) return

    call neko_log%message(" ")
    call neko_log%message("------------------------------")
    call neko_log%message("Advancing POD transience phase")
    call neko_log%message("------------------------------")
    call neko_log%message(" ")

    call profiler_start_region("POD transience")

    call dt_controller%init(neko_case%params)
    target_time = neko_case%time%start_time + this%transience_time
    loop_start = MPI_WTIME()

    do while (neko_case%time%t .lt. target_time)
       call simulation_step(neko_case, dt_controller, loop_start)
    end do

    this%time_shift = neko_case%time%t
    this%pod_tstep = 0
    call this%set_n_timesteps(0)
    this%transience_applied = .true.

    call profiler_end_region("POD transience")
  end subroutine POD_state_recover_apply_transience

  subroutine POD_state_recover_stream_fields(this, neko_case)
    class(POD_state_recover_t), intent(inout) :: this
    class(case_t), intent(inout) :: neko_case
    type(field_t), pointer :: u, v, w, s
    integer :: n

    u => neko_case%fluid%u
    v => neko_case%fluid%v
    w => neko_case%fluid%w
    if (this%include_scalar) s => this%s
    n = u%dof%size()

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(u%x, u%x_d, n, DEVICE_TO_HOST, sync=.true.)
       call device_memcpy(v%x, v%x_d, n, DEVICE_TO_HOST, sync=.true.)
       call device_memcpy(w%x, w%x_d, n, DEVICE_TO_HOST, sync=.true.)
       if (this%include_scalar) then
          call device_memcpy(s%x, s%x_d, n, DEVICE_TO_HOST, sync=.true.)
       end if
    end if

    call this%dstream%stream(u%x)
    call this%dstream%stream(v%x)
    call this%dstream%stream(w%x)
    if (this%include_scalar) call this%dstream%stream(s%x)
  end subroutine POD_state_recover_stream_fields


  !> Stream forward state for POD updates.
  !! @param[inout] this POD state recovery instance.
  !! @param[inout] neko_case Case data structure.
  !! @param[in] time Current time state.
  subroutine POD_state_recover_save(this, neko_case, time)
    class(POD_state_recover_t), intent(inout) :: this
    class(case_t), intent(inout) :: neko_case
    type(time_state_t), intent(in) :: time
    if (.not. this%enabled) return

    if (.not. this%transience_applied) then
       call POD_state_recover_apply_transience(this, neko_case)
       ! Python blocks on the initial field snapshot before it enters the
       ! control loop, so we delay that first snapshot until the transience
       ! has been skipped.
       call POD_state_recover_stream_fields(this, neko_case)
       return
    end if

    if (time%t .le. this%time_shift) return

    this%pod_tstep = this%pod_tstep + 1
    call this%set_n_timesteps(max(this%get_n_timesteps(), this%pod_tstep))
    if (mod(this%pod_tstep, this%i_stream) .ne. 0) return

    call neko_log%message(" ")
    call neko_log%message("----------------")
    call neko_log%message("Streaming fields")
    call neko_log%message("----------------")
    call neko_log%message(" ")

    call profiler_start_region("POD save")

    if (this%ctrl%inited) then
       call this%ctrl%send(MODE_FORWARD, PHASE_FWD_RUNNING, &
            int(this%pod_tstep, c_int), &
            real(neko_case%time%t - this%time_shift, c_double))
    end if

    call POD_state_recover_stream_fields(this, neko_case)

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
    type(time_state_t) :: time_out
    real(kind=rp) :: t_pod

    if (.not. this%enabled) return

    ! First restore() call is the phase boundary forward->adjoint
    if (.not. this%have_received_modes) then
       call POD_state_recover_recieve_modes(this, time)
    end if

    ! Emit ADJ_RUNNING only once (avoid flooding SST)
    if (this%ctrl%inited .and. .not. this%adj_running_sent) then
       call this%ctrl%send(MODE_ADJOINT, PHASE_ADJ_RUNNING, &
            int(time%tstep, c_int), real(time%t, c_double))
       this%adj_running_sent = .true.
    end if

    call profiler_start_region("POD restore")
    t_pod = time%t
    call interpolate_time_coeffs_vec(this%a_interp, this%time_coefs, t_pod)
    call reconstruct_from_coeffs(this, neko_case, this%a_interp)
    if (this%output_reconstruction) then
       if (recon_should_output(this, time, time_out)) then
          call this%recon_output%sample(time_out%t)
       end if
    end if
    call profiler_end_region("POD restore")
  end subroutine POD_state_recover_restore

  !> Receive POD modes at the forward-to-adjoint boundary.
  !! @param[inout] this POD state recovery instance.
  !! @param[in] time Target time state.
  subroutine POD_state_recover_recieve_modes(this, time)
    class(POD_state_recover_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    integer :: i, ierr, n_lines, nrows, ncols, n
    integer(c_int) :: mode_cmd, phase_cmd

    call profiler_start_region("POD recieve modes")

    if (this%ctrl%inited) then
      call this%ctrl%send(MODE_FORWARD, PHASE_FWD_DONE, &
           int(time%tstep, c_int), real(time%t, c_double))

      mode_cmd  = MODE_FORWARD
      phase_cmd = PHASE_FWD_DONE

      ! BLOCK until Python says "go adjoint"
      call this%ctrl%recieve(mode_cmd, phase_cmd)

      if (mode_cmd /= MODE_ADJOINT) then
         call neko_error('Expected MODE_ADJOINT from Python at ' // &
              'forward->adjoint boundary.')
      end if
    end if

    this%adjoint_started = .true.

    ! Receive modes (Python streams after sending cmd)
    do i = 1, this%n_modes
       call this%dstream%recieve(this%u_modes(i)%x)
       call this%dstream%recieve(this%v_modes(i)%x)
       call this%dstream%recieve(this%w_modes(i)%x)
       if (this%include_scalar) then
          call this%dstream%recieve(this%s_modes(i)%x)
       end if
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
          if (this%include_scalar) then
             call device_memcpy(this%s_modes(i)%x, this%s_modes(i)%x_d, &
                  n, HOST_TO_DEVICE, sync=.true.)
          end if
       end do
    end if

    ! Read CSV once
    n_lines = csv_file_count_lines(this%csv_reader)
    nrows = n_lines - 1
    ncols = 1 + this%n_modes
    call this%time_coefs%init(nrows, ncols)
    call this%csv_reader%read(this%time_coefs)

    call MPI_Allreduce(MPI_IN_PLACE, this%time_coefs%x, &
         this%time_coefs%size(), mpi_real_precision, MPI_SUM, &
         neko_comm, ierr)

    this%have_received_modes = .true.
    call profiler_end_region("POD recieve modes")
  end subroutine POD_state_recover_recieve_modes

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

  integer function precision_from_string(str, name) result(precision)
    character(len=*), intent(in) :: str
    character(len=*), intent(in) :: name

    select case (trim(lower_string(str)))
    case ('sp', 'single')
       precision = sp
    case ('dp', 'double')
       precision = dp
    case default
       call neko_error('Invalid ' // trim(name) // '. Expected ''sp'', ' // &
            '''dp'', ''single'', or ''double''.')
    end select
  end function precision_from_string

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

    open(file = trim(file_in%get_fname()), status='old', &
         newunit=file_unit, iostat=ierr)
    if (ierr .ne. 0) then
       call neko_error("Error opening " // trim(file_in%get_fname()))
    end if
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
       a_out%x(j) = (1.0_rp - w) * time_coefs%x(i0, j+1) + &
            w * time_coefs%x(i1, j+1)
    end do
  end subroutine interpolate_time_coeffs_vec


  subroutine reconstruct_from_coeffs(this, neko_case, a)
    class(POD_state_recover_t), intent(inout) :: this
    class(case_t), intent(inout)              :: neko_case
    type(vector_t), intent(in)                :: a
    integer :: j
    type(field_t), pointer :: u, v, w, s

    u => neko_case%fluid%u
    v => neko_case%fluid%v
    w => neko_case%fluid%w
    if (this%include_scalar) s => this%s

    call field_rzero(u)
    call field_rzero(v)
    call field_rzero(w)
    if (this%include_scalar) call field_rzero(s)

    do j = 1, this%n_modes
       call field_add2s2(u, this%u_modes(j), a%x(j))
       call field_add2s2(v, this%v_modes(j), a%x(j))
       call field_add2s2(w, this%w_modes(j), a%x(j))
       if (this%include_scalar) then
          call field_add2s2(s, this%s_modes(j), a%x(j))
       end if
    end do
  end subroutine reconstruct_from_coeffs

end module simulation_POD_state_recover
