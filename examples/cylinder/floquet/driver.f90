program usrneko
  use LightKrylov
  use LightKrylov, only: wp => dp
  use LightKrylov_AbstractVectors, only: init_like_basis
  use LightKrylov_Constants
  use LightKrylov_Logger
  use field_math, only: field_copy
  use field_output, only: field_output_t
  use file, only: file_t
  use fld_file_data, only: fld_file_data_t
  use fld_file, only: fld_file_t
  use json_module, only: json_file
  use json_utils, only: json_get_or_default
  use json_utils_ext, only: json_read_file
  use neko, only: neko_finalize, neko_init
  use case, only: case_t
  use linop_transient_baseflow, only: linop_transient_baseflow_t
  use mpi_f08, only: MPI_WTIME
  use neko_ext, only: reset, reset_adjoint
  use neko_stability_propagator, only: copy_state_to_primal, &
       copy_state_to_perturb, fit_timestep_to_interval
  use neko_top, only: neko_top_register_types
  use neko_vector, only: state_vector_t
  use num_types, only: neko_rp => rp, neko_sp => sp
  use perturb_case, only: perturb_case_t
  use simulation, only: simulation_finalize, simulation_init, simulation_step
  use simulation_perturb, only: simulation_lns_step, &
       simulation_perturb_finalize, simulation_perturb_init
  use stability_case_utils, only: stability_case_free, stability_case_init
  use time_step_controller, only: time_step_controller_t
  use utils, only: extract_fld_file_index, filename_chsuffix, neko_error
  implicit none

  character(len=*), parameter :: this_module = 'Cylinder Re180 Floquet'
  type(case_t), target :: base_case
  type(perturb_case_t), target :: lns_case
  type(linop_transient_baseflow_t) :: A
  type(state_vector_t), allocatable :: X(:), baseflow
  type(state_vector_t) :: vector_mold
  complex(kind=wp), allocatable :: multiplier(:), exponent(:)
  real(kind=wp), allocatable :: residuals(:)
  real(kind=wp) :: period, eigs_tol, old_dt, field_period
  integer :: argc, nev, kdim, info, i, movie_frames, mode_movies, nsteps
  logical :: write_baseflow_movie, write_mode_movies
  character(len=256) :: parameter_file
  character(len=:), allocatable :: baseflow_file
  type(json_file) :: parameters

  call neko_init()
  call neko_top_register_types()
  call logger_setup()

  argc = command_argument_count()
  if (argc .lt. 1) call neko_error('Missing parameter file')
  call get_command_argument(1, parameter_file)

  parameters = json_read_file(trim(parameter_file))
  call json_get_or_default(parameters, 'case.lightkrylov.floquet.period', &
       period, 0.0_wp)
  call json_get_or_default(parameters, &
       'case.fluid.initial_condition.file_name', baseflow_file, '')
  call json_get_or_default(parameters, &
       'case.lightkrylov.floquet.number_of_eigenvalues', nev, 4)
  call json_get_or_default(parameters, &
       'case.lightkrylov.floquet.krylov_dimension', kdim, 4 * nev)
  call json_get_or_default(parameters, &
       'case.lightkrylov.floquet.tolerance', eigs_tol, rtol_dp)
  call json_get_or_default(parameters, &
       'case.lightkrylov.floquet.movie_frames', movie_frames, 32)
  call json_get_or_default(parameters, &
       'case.lightkrylov.floquet.write_baseflow_movie', &
       write_baseflow_movie, .true.)
  call json_get_or_default(parameters, &
       'case.lightkrylov.floquet.write_mode_movies', write_mode_movies, .true.)
  call json_get_or_default(parameters, &
       'case.lightkrylov.floquet.number_of_mode_movies', mode_movies, nev)

  call stability_case_init(base_case, parameters)
  field_period = read_field_time(baseflow_file)
  if (field_period > 0.0_wp) then
     period = field_period
  else if (period <= 0.0_wp) then
     period = real(base_case%time%end_time, kind=wp)
  end if
  if (period <= 0.0_wp) call neko_error('Floquet period must be positive')
  if (nev < 1) call neko_error('Number of eigenvalues must be positive')
  if (kdim < nev) call neko_error('Krylov dimension must be at least nev')
  if (eigs_tol <= 0.0_wp) call neko_error('Eigensolver tolerance must be positive')
  if (movie_frames < 2) call neko_error('Floquet movie_frames must be >= 2')
  if (mode_movies < 0) call neko_error('number_of_mode_movies must be nonnegative')

  lns_case%if_adjoint = .false.
  call lns_case%init(base_case)

  vector_mold%coef => base_case%fluid%c_Xh
  vector_mold%if_2d = (base_case%fluid%c_Xh%msh%gdim .eq. 2)
  allocate(baseflow)
  call baseflow%init_like(vector_mold)
  call field_copy(baseflow%u, base_case%fluid%u)
  call field_copy(baseflow%v, base_case%fluid%v)
  call field_copy(baseflow%w, base_case%fluid%w)
  call write_state_fields('baseflow_floquet_t0', baseflow, 0)

  old_dt = real(base_case%time%dt, kind=wp)
  nsteps = max(1, ceiling(period / old_dt))
  call A%init(base_case, lns_case)
  call A%set_baseflow(baseflow, period, nsteps)

  allocate(X(nev))
  call init_like_basis(X, vector_mold)
  call initialize_krylov_subspace(X)

  call eigs(A, X, multiplier, residuals, info, kdim=kdim, tolerance=eigs_tol)
  call check_info(info, 'eigs', module=this_module, procedure='main')

  allocate(exponent(size(multiplier)))
  exponent = log(multiplier) / period
  call write_floquet_spectrum(multiplier, exponent, residuals)

  if (write_baseflow_movie) then
     call write_baseflow_period_movie(baseflow, period, movie_frames, nsteps)
  end if

  if (write_mode_movies) then
     do i = 1, min(mode_movies, nev)
        call write_mode_period_movie(i, baseflow, X(i), period, movie_frames, &
             nsteps)
     end do
  end if

  call A%free()
  do i = 1, nev
     call X(i)%free()
  end do
  deallocate(X)
  call baseflow%free()
  deallocate(baseflow)
  if (allocated(multiplier)) deallocate(multiplier)
  if (allocated(exponent)) deallocate(exponent)
  if (allocated(residuals)) deallocate(residuals)
  call lns_case%free()
  call stability_case_free(base_case)
  call neko_finalize()

contains

  real(wp) function read_field_time(field_name) result(field_time)
    character(len=*), intent(in) :: field_name
    type(file_t) :: f
    type(fld_file_data_t) :: fld_data
    character(len=1024) :: fld_name
    integer :: sample_idx

    field_time = 0.0_wp
    if (len_trim(field_name) .eq. 0) return

    sample_idx = extract_fld_file_index(trim(field_name), -1)
    if (sample_idx .eq. -1) return

    call filename_chsuffix(trim(field_name), fld_name, 'fld')
    call f%init(trim(fld_name))
    call f%set_counter(sample_idx)
    call f%read(fld_data)
    field_time = real(fld_data%time, kind=wp)
    call fld_data%free()
    call f%free()
  end function read_field_time

  subroutine write_floquet_spectrum(mu, lambda, res)
    complex(kind=wp), intent(in) :: mu(:), lambda(:)
    real(kind=wp), intent(in) :: res(:)
    integer :: j

    open(unit=1234, file='floquet_spectrum.csv', status='replace', &
         action='write')
    write(1234,'(A)') 'index,multiplier_real,multiplier_imag,' // &
         'multiplier_abs,multiplier_phase,exponent_real,exponent_imag,residual'
    do j = 1, size(mu)
       write(1234,'(I0,7(",",ES24.16))') j, real(mu(j), kind=wp), &
            aimag(mu(j)), abs(mu(j)), atan2(aimag(mu(j)), &
            real(mu(j), kind=wp)), real(lambda(j), kind=wp), &
            aimag(lambda(j)), res(j)
    end do
    close(1234)
  end subroutine write_floquet_spectrum

  subroutine write_state_fields(name, vec, idx)
    character(len=*), intent(in) :: name
    type(state_vector_t), intent(inout) :: vec
    integer, intent(in) :: idx
    type(field_output_t) :: output

    call execute_command_line('mkdir -p viz')
    call output%init(trim(name), 3, precision=neko_sp, path='viz/', &
         format='nek5000')
    select type (ft => output%file_%file_type)
    type is (fld_file_t)
       ft%skip_pressure = .true.
       ft%skip_velocity = .false.
       ft%skip_temperature = .true.
    end select
    call output%fields%assign_to_field(1, vec%u)
    call output%fields%assign_to_field(2, vec%v)
    call output%fields%assign_to_field(3, vec%w)
    call output%set_counter(idx - 1)
    call output%sample(real(idx, kind=neko_rp))
    call output%free()
  end subroutine write_state_fields

  subroutine write_baseflow_period_movie(initial_baseflow, time_interval, &
       nframes, fixed_nsteps)
    type(state_vector_t), intent(in) :: initial_baseflow
    real(wp), intent(in) :: time_interval
    integer, intent(in) :: nframes, fixed_nsteps
    type(time_step_controller_t) :: dt_controller
    type(field_output_t) :: output
    real(wp) :: old_end_time, old_base_dt, next_sample, sample_interval
    real(wp) :: loop_start
    integer :: next_frame

    old_end_time = real(base_case%time%end_time, kind=wp)
    old_base_dt = real(base_case%time%dt, kind=wp)
    base_case%time%end_time = real(time_interval, kind(base_case%time%end_time))
    call fit_timestep_to_interval(base_case%time, time_interval, fixed_nsteps)

    call reset(base_case, 0, trim(base_case%f_out%file_%get_base_fname()))
    call copy_state_to_primal(base_case, initial_baseflow)

    call execute_command_line('mkdir -p viz')
    call output%init('baseflow_floquet_movie', 3, precision=neko_sp, &
         path='viz/', format='nek5000')
    select type (ft => output%file_%file_type)
    type is (fld_file_t)
       ft%skip_pressure = .true.
       ft%skip_velocity = .false.
       ft%skip_temperature = .true.
    end select
    call output%fields%assign_to_field(1, base_case%fluid%u)
    call output%fields%assign_to_field(2, base_case%fluid%v)
    call output%fields%assign_to_field(3, base_case%fluid%w)

    call dt_controller%init(base_case%params)
    call simulation_init(base_case, dt_controller)
    call sample_periodic_output(output, base_case%time%t, 0)

    sample_interval = time_interval / real(nframes - 1, kind=wp)
    next_frame = 1
    next_sample = sample_interval
    loop_start = MPI_WTIME()
    do while (.not. base_case%time%is_done())
       call simulation_step(base_case, dt_controller, loop_start)
       if (next_frame < nframes .and. real(base_case%time%t, kind=wp) >= &
            next_sample - 100.0_wp * epsilon(next_sample)) then
          call sample_periodic_output(output, base_case%time%t, next_frame)
          next_frame = next_frame + 1
          next_sample = sample_interval * real(next_frame, kind=wp)
       end if
    end do

    if (next_frame < nframes) then
       call sample_periodic_output(output, base_case%time%t, next_frame)
    end if

    call simulation_finalize(base_case)
    call output%free()
    base_case%time%end_time = real(old_end_time, kind(base_case%time%end_time))
    base_case%time%dt = real(old_base_dt, kind(base_case%time%dt))
    base_case%time%dtlag = base_case%time%dt
  end subroutine write_baseflow_period_movie

  subroutine write_mode_period_movie(mode_index, initial_baseflow, &
       initial_mode, time_interval, nframes, fixed_nsteps)
    integer, intent(in) :: mode_index
    type(state_vector_t), intent(in) :: initial_baseflow, initial_mode
    real(wp), intent(in) :: time_interval
    integer, intent(in) :: nframes, fixed_nsteps
    type(time_step_controller_t) :: primal_dt_controller
    type(time_step_controller_t) :: perturb_dt_controller
    type(field_output_t) :: output
    real(wp) :: old_base_end, old_lns_end, old_base_dt, old_lns_dt
    real(wp) :: next_sample, sample_interval, loop_start
    real(neko_rp) :: cfl
    integer :: next_frame
    character(len=64) :: output_name

    write(output_name, '("floquet_mode_", I3.3)') mode_index
    old_base_end = real(base_case%time%end_time, kind=wp)
    old_lns_end = real(lns_case%time%end_time, kind=wp)
    old_base_dt = real(base_case%time%dt, kind=wp)
    old_lns_dt = real(lns_case%time%dt, kind=wp)

    base_case%time%end_time = real(time_interval, kind(base_case%time%end_time))
    lns_case%time%end_time = real(time_interval, kind(lns_case%time%end_time))
    call fit_timestep_to_interval(base_case%time, time_interval, fixed_nsteps)
    lns_case%time%dt = base_case%time%dt
    lns_case%time%dtlag = lns_case%time%dt

    call reset(base_case, 0, trim(base_case%f_out%file_%get_base_fname()))
    call reset_adjoint(lns_case, base_case, 0, &
         trim(lns_case%f_out%file_%get_base_fname()))
    call copy_state_to_primal(base_case, initial_baseflow)
    call copy_state_to_perturb(lns_case, initial_mode)

    call execute_command_line('mkdir -p viz')
    call output%init(trim(output_name), 3, precision=neko_sp, path='viz/', &
         format='nek5000')
    select type (ft => output%file_%file_type)
    type is (fld_file_t)
       ft%skip_pressure = .true.
       ft%skip_velocity = .false.
       ft%skip_temperature = .true.
    end select
    call output%fields%assign_to_field(1, lns_case%fluid_pert%u_pert)
    call output%fields%assign_to_field(2, lns_case%fluid_pert%v_pert)
    call output%fields%assign_to_field(3, lns_case%fluid_pert%w_pert)

    call primal_dt_controller%init(base_case%params)
    call perturb_dt_controller%init(base_case%params)
    call simulation_init(base_case, primal_dt_controller)
    call simulation_perturb_init(lns_case, perturb_dt_controller)
    call sample_periodic_output(output, lns_case%time%t, 0)

    cfl = lns_case%fluid_pert%compute_cfl(lns_case%time%dt)
    sample_interval = time_interval / real(nframes - 1, kind=wp)
    next_frame = 1
    next_sample = sample_interval
    loop_start = MPI_WTIME()
    do while (.not. base_case%time%is_done())
       call simulation_step(base_case, primal_dt_controller, loop_start)
       call simulation_lns_step(lns_case, perturb_dt_controller, cfl, &
            loop_start)
       if (next_frame < nframes .and. real(base_case%time%t, kind=wp) >= &
            next_sample - 100.0_wp * epsilon(next_sample)) then
          call sample_periodic_output(output, lns_case%time%t, next_frame)
          next_frame = next_frame + 1
          next_sample = sample_interval * real(next_frame, kind=wp)
       end if
    end do

    if (next_frame < nframes) then
       call sample_periodic_output(output, lns_case%time%t, next_frame)
    end if

    call simulation_perturb_finalize(lns_case)
    call simulation_finalize(base_case)
    call output%free()

    base_case%time%end_time = real(old_base_end, kind(base_case%time%end_time))
    lns_case%time%end_time = real(old_lns_end, kind(lns_case%time%end_time))
    base_case%time%dt = real(old_base_dt, kind(base_case%time%dt))
    lns_case%time%dt = real(old_lns_dt, kind(lns_case%time%dt))
    base_case%time%dtlag = base_case%time%dt
    lns_case%time%dtlag = lns_case%time%dt
  end subroutine write_mode_period_movie

  subroutine sample_periodic_output(output, time_value, idx)
    type(field_output_t), intent(inout) :: output
    real(neko_rp), intent(in) :: time_value
    integer, intent(in) :: idx

    call output%set_counter(idx - 1)
    call output%sample(time_value)
  end subroutine sample_periodic_output

end program usrneko
