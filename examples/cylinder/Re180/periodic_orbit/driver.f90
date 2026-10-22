program usrneko
  use LightKrylov
  use LightKrylov, only: wp => dp
  use LightKrylov_Constants
  use LightKrylov_Logger
  use field_math, only: field_copy
  use field_output, only: field_output_t
  use fld_file, only: fld_file_t
  use json_module, only: json_file
  use json_utils, only: json_get_or_default
  use json_utils_ext, only: json_read_file
  use neko, only: neko_finalize, neko_init
  use case, only: case_t
  use neko_ext, only: reset
  use neko_ext_vector, only: ext_state_vector_t
  use periodic_orbit_system, only: periodic_orbit_jacobian_t, &
       periodic_orbit_system_t
  use neko_stability_propagator, only: copy_state_to_primal, &
       fit_timestep_to_interval
  use neko_top, only: neko_top_register_types
  use neko_vector, only: state_vector_t
  use num_types, only: neko_rp => rp, neko_sp => sp
  use perturb_case, only: perturb_case_t
  use simulation, only: simulation_finalize, simulation_init, simulation_step
  use stability_case_utils, only: stability_case_free, stability_case_init
  use time_step_controller, only: time_step_controller_t
  use utils, only: neko_error
  use mpi_f08, only: MPI_WTIME
  implicit none

  character(len=*), parameter :: this_module = 'Cylinder Re180 periodic orbit'
  type(case_t), target :: base_case
  type(perturb_case_t), target :: lns_case
  type(periodic_orbit_system_t), target :: periodic
  type(periodic_orbit_jacobian_t) :: jacobian
  type(ext_state_vector_t) :: orbit
  type(state_vector_t) :: vector_mold
  type(newton_dp_opts) :: newton_options
  real(wp) :: newton_tol, initial_period, transient_time
  real(wp) :: target_dt
  character(len=:), allocatable :: base_output_base_fname
  integer :: argc, info, movie_frames
  logical :: write_outer_iterations
  character(len=256) :: parameter_file
  character(len=:), allocatable :: movie_output_prefix, outer_output_prefix
  type(json_file) :: parameters

  call neko_init()
  call neko_top_register_types()
  call logger_setup()

  argc = command_argument_count()
  if (argc .lt. 1) call neko_error('Missing parameter file')
  call get_command_argument(1, parameter_file)

  parameters = json_read_file(trim(parameter_file))
  call json_get_or_default(parameters, &
       'case.lightkrylov.periodic_orbit.initial_period', initial_period, &
       5.158_wp)
  call json_get_or_default(parameters, &
       'case.lightkrylov.periodic_orbit.absolute_tolerance', newton_tol, &
       1.0e-6_wp)
  call json_get_or_default(parameters, &
       'case.lightkrylov.periodic_orbit.transient_time', &
       transient_time, 0.0_wp)
  call json_get_or_default(parameters, &
       'case.lightkrylov.periodic_orbit.movie_frames', movie_frames, 32)
  call json_get_or_default(parameters, &
       'case.lightkrylov.periodic_orbit.write_outer_iterations', &
       write_outer_iterations, .false.)
  call json_get_or_default(parameters, &
       'case.lightkrylov.periodic_orbit.outer_iteration_output_prefix', &
       outer_output_prefix, 'outer_periodic_orbit_Re180')
  call json_get_or_default(parameters, &
       'case.lightkrylov.periodic_orbit.movie_output_prefix', &
       movie_output_prefix, 'periodic_orbit_Re180_movie')

  if (initial_period <= 0.0_wp) call neko_error('Initial period must be positive')
  if (transient_time < 0.0_wp) call neko_error('Transient time must be non-negative')
  if (newton_tol <= 0.0_wp) call neko_error('Newton tolerance must be positive')
  if (movie_frames < 2) call neko_error('Periodic orbit movie_frames must be >= 2')

  call stability_case_init(base_case, parameters)
  target_dt = real(base_case%time%dt, kind=wp)
  if (target_dt <= 0.0_wp) call neko_error('Case timestep must be positive')
  base_output_base_fname = trim(base_case%f_out%file_%get_base_fname())
  lns_case%if_adjoint = .false.
  call lns_case%init(base_case)

  if (transient_time > 0.0_wp) then
     call run_forward_transient(transient_time)
  end if

  vector_mold%coef => base_case%fluid%c_Xh
  vector_mold%if_2d = (base_case%fluid%c_Xh%msh%gdim .eq. 2)
  call orbit%q%init_like(vector_mold)
  call field_copy(orbit%q%u, base_case%fluid%u)
  call field_copy(orbit%q%v, base_case%fluid%v)
  call field_copy(orbit%q%w, base_case%fluid%w)
  call field_copy(orbit%q%p, base_case%fluid%p)
  orbit%T = initial_period

  call periodic%init(base_case)
  call periodic%set_fixed_step_count(initial_period)
  if (write_outer_iterations) then
     call periodic%enable_outer_iteration_output(trim(outer_output_prefix))
  end if

  periodic%jacobian = jacobian
  select type (f => periodic%jacobian)
  type is (periodic_orbit_jacobian_t)
     call f%init(periodic, lns_case)
  end select

  newton_options = newton_dp_opts(maxiter=40, ifbisect=.false.)
  call newton(periodic, orbit, gmres_rdp, info, &
       atol=newton_tol, options=newton_options, scheduler=dynamic_tol_dp)
  if (info < 0) call check_info(info, 'newton', module=this_module, &
       procedure='main')

  call write_baseflow_file(orbit%q, orbit%T)

  ! Close the Newton flow-map simulations before starting the independent
  ! integration used to write the movie frames.
  call periodic%free()

  call field_copy(base_case%fluid%u, orbit%q%u)
  call field_copy(base_case%fluid%v, orbit%q%v)
  call field_copy(base_case%fluid%w, orbit%q%w)
  call field_copy(base_case%fluid%p, orbit%q%p)
  call write_periodic_movie(orbit%q, orbit%T, movie_frames, &
       movie_output_prefix)

  call orbit%free()
  call lns_case%free()
  call stability_case_free(base_case)
  call neko_finalize()
contains

  subroutine write_state_fields(name, vec, idx)
    character(len=*), intent(in) :: name
    type(state_vector_t), intent(inout) :: vec
    integer, intent(in) :: idx
    type(field_output_t) :: output

    call execute_command_line('mkdir -p viz')
    call output%init(trim(name), 4, precision=neko_sp, path='viz/', format='nek5000')
    select type (ft => output%file_%file_type)
    type is (fld_file_t)
       ft%skip_pressure = .false.
       ft%skip_velocity = .false.
       ft%skip_temperature = .true.
    end select
    call output%fields%assign_to_field(1, vec%p)
    call output%fields%assign_to_field(2, vec%u)
    call output%fields%assign_to_field(3, vec%v)
    call output%fields%assign_to_field(4, vec%w)
    call output%set_counter(idx - 1)
    call output%sample(real(idx, kind=neko_rp))
    call output%free()
  end subroutine write_state_fields

  subroutine write_baseflow_file(vec, period)
    type(state_vector_t), intent(inout) :: vec
    real(wp), intent(in) :: period
    type(field_output_t) :: output

    call output%init('baseflow', 4, precision=neko_sp, format='nek5000')
    select type (ft => output%file_%file_type)
    type is (fld_file_t)
       ft%skip_pressure = .false.
       ft%skip_velocity = .false.
       ft%skip_temperature = .true.
    end select
    call output%fields%assign_to_field(1, vec%p)
    call output%fields%assign_to_field(2, vec%u)
    call output%fields%assign_to_field(3, vec%v)
    call output%fields%assign_to_field(4, vec%w)
    call output%set_counter(-1)
    call output%sample(real(period, kind=neko_rp))
    call output%free()
  end subroutine write_baseflow_file

  subroutine run_forward_transient(time_interval)
    real(wp), intent(in) :: time_interval
    type(time_step_controller_t) :: dt_controller
    real(wp) :: old_end_time, old_dt, loop_start

    old_end_time = real(base_case%time%end_time, kind=wp)
    old_dt = real(base_case%time%dt, kind=wp)
    base_case%time%end_time = real(time_interval, kind(base_case%time%end_time))
    call fit_timestep_to_interval(base_case%time, time_interval)

    call dt_controller%init(base_case%params)
    call simulation_init(base_case, dt_controller)
    loop_start = MPI_WTIME()
    do while (.not. base_case%time%is_done())
       call simulation_step(base_case, dt_controller, loop_start)
    end do
    call simulation_finalize(base_case)

    base_case%time%end_time = real(old_end_time, kind(base_case%time%end_time))
    base_case%time%dt = real(old_dt, kind(base_case%time%dt))
  end subroutine run_forward_transient

  subroutine write_periodic_movie(initial_state, period, nframes, output_prefix)
    type(state_vector_t), intent(in) :: initial_state
    real(wp), intent(in) :: period
    integer, intent(in) :: nframes
    character(len=*), intent(in) :: output_prefix
    type(time_step_controller_t) :: dt_controller
    type(field_output_t) :: output
    real(wp) :: old_end_time, old_dt, next_sample, sample_interval
    real(wp) :: loop_start
    integer :: next_frame, nsteps

    call execute_command_line('mkdir -p viz')

    old_end_time = real(base_case%time%end_time, kind=wp)
    old_dt = real(base_case%time%dt, kind=wp)
    nsteps = max(1, ceiling(period / old_dt))
    base_case%time%end_time = &
         real(period, kind(base_case%time%end_time))
    base_case%time%dt = &
         real(period / real(nsteps, kind=wp), &
         kind(base_case%time%dt))

    call reset(base_case, 0, base_output_base_fname)
    call copy_state_to_primal(base_case, initial_state)

    call output%init(trim(output_prefix), 4, precision=neko_sp, &
         path='viz/', format='nek5000')
    select type (ft => output%file_%file_type)
    type is (fld_file_t)
       ft%skip_pressure = .false.
       ft%skip_velocity = .false.
       ft%skip_temperature = .true.
    end select
    call output%fields%assign_to_field(1, base_case%fluid%p)
    call output%fields%assign_to_field(2, base_case%fluid%u)
    call output%fields%assign_to_field(3, base_case%fluid%v)
    call output%fields%assign_to_field(4, base_case%fluid%w)

    call dt_controller%init(base_case%params)
    call simulation_init(base_case, dt_controller)

    sample_interval = period / real(nframes - 1, kind=wp)
    next_frame = 0
    call output%set_counter(next_frame - 1)
    call output%sample(real(base_case%time%t, kind=neko_rp))
    next_frame = 1
    next_sample = sample_interval

    loop_start = MPI_WTIME()
    do while (.not. base_case%time%is_done())
       call simulation_step(base_case, dt_controller, loop_start)

       if (next_frame < nframes .and. &
            real(base_case%time%t, kind=wp) >= &
            next_sample - 100.0_wp * epsilon(next_sample)) then
          call output%set_counter(next_frame - 1)
          call output%sample(real(base_case%time%t, kind=neko_rp))
          next_frame = next_frame + 1
          next_sample = sample_interval * real(next_frame, kind=wp)
       end if
    end do

    if (next_frame < nframes) then
      call output%set_counter(next_frame - 1)
      call output%sample(real(base_case%time%t, kind=neko_rp))
    end if

    call simulation_finalize(base_case)
    call output%free()
    base_case%time%end_time = real(old_end_time, kind(base_case%time%end_time))
    base_case%time%dt = real(old_dt, kind(base_case%time%dt))
  end subroutine write_periodic_movie
end program usrneko
