program usrneko
  use LightKrylov, only: wp => dp
  use LightKrylov_Logger, only: logger_setup
  use field_math, only: field_copy
  use file, only: file_t
  use fld_file_data, only: fld_file_data_t
  use json_module, only: json_file
  use json_utils, only: json_get_or_default
  use json_utils_ext, only: json_read_file
  use neko, only: neko_finalize, neko_init
  use neko_ext, only: reset
  use case, only: case_t
  use linop_transient_baseflow, only: linop_transient_baseflow_t
  use neko_top, only: neko_top_register_types
  use neko_stability_propagator, only: copy_primal_to_state, &
       copy_state_to_primal
  use neko_vector, only: state_vector_t
  use perturb_case, only: perturb_case_t
  use simulation, only: simulation_step
  use stability_case_utils, only: stability_case_free, stability_case_init
  use time_step_controller, only: time_step_controller_t
  use utils, only: extract_fld_file_index, filename_chsuffix, neko_error
  implicit none

  character(len=*), parameter :: this_module = 'Cylinder adjoint consistency'
  type(case_t), target :: base_case
  type(perturb_case_t), target :: lns_case, adjoint_case
  type(linop_transient_baseflow_t) :: A
  type(state_vector_t) :: vector_mold, baseflow, x, y, Ax, Aty
  type(state_vector_t) :: Ax_repeat, difference
  type(state_vector_t) :: fd_base, fd_perturbed, fd_result, fd_result_pert
  type(state_vector_t) :: fd_difference
  type(time_step_controller_t) :: fd_dt_controller
  real(kind=wp) :: period, old_dt, fitted_dt, test_interval
  real(kind=wp) :: lhs, rhs, abs_error, rel_error, tolerance
  real(kind=wp) :: repeat_error, max_rel_error, max_repeat_error
  real(kind=wp) :: finite_error, max_finite_error, epsilon_fd
  real(kind=wp), parameter :: epsilon_values(3) = [1.0e-3_wp, &
       1.0e-4_wp, 1.0e-5_wp]
  integer :: argc, full_nsteps, test_nsteps, number_random_tests, itest
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
  call json_get_or_default(parameters, &
       'case.fluid.initial_condition.file_name', baseflow_file, '')
  call json_get_or_default(parameters, &
       'case.lightkrylov.adjoint_test.number_of_timesteps', test_nsteps, 4)
  call json_get_or_default(parameters, &
       'case.lightkrylov.adjoint_test.number_of_random_tests', &
       number_random_tests, 3)
  call json_get_or_default(parameters, &
       'case.lightkrylov.adjoint_test.tolerance', tolerance, 1.0e-10_wp)

  call stability_case_init(base_case, parameters)
  if (len_trim(baseflow_file) == 0) &
       call neko_error('Adjoint consistency baseflow file_name is required')

  period = read_field_time(baseflow_file)
  if (period <= 0.0_wp) &
       call neko_error('Adjoint consistency period must be read from baseflow')
  if (test_nsteps < 1) &
       call neko_error('Adjoint consistency number_of_timesteps must be positive')
  if (number_random_tests < 1) &
       call neko_error('Adjoint consistency number_of_random_tests must be positive')
  if (tolerance <= 0.0_wp) &
       call neko_error('Adjoint consistency tolerance must be positive')

  old_dt = real(base_case%time%dt, kind=wp)
  if (old_dt <= 0.0_wp) &
       call neko_error('Adjoint consistency timestep must be positive')
  full_nsteps = max(1, ceiling(period / old_dt))
  if (test_nsteps > full_nsteps) &
       call neko_error('Adjoint consistency test interval exceeds one period')

  fitted_dt = period / real(full_nsteps, kind=wp)
  test_interval = real(test_nsteps, kind=wp) * fitted_dt

  lns_case%if_adjoint = .false.
  call lns_case%init(base_case)
  adjoint_case%if_adjoint = .true.
  call adjoint_case%init(base_case)

  vector_mold%coef => base_case%fluid%c_Xh
  vector_mold%if_2d = (base_case%fluid%c_Xh%msh%gdim .eq. 2)

  call baseflow%init_like(vector_mold)
  call field_copy(baseflow%u, base_case%fluid%u)
  call field_copy(baseflow%v, base_case%fluid%v)
  call field_copy(baseflow%w, base_case%fluid%w)

  call A%init(base_case, lns_case, adjoint_case)

  call x%init_like(vector_mold)
  call y%init_like(vector_mold)
  call Ax_repeat%init_like(vector_mold)
  call difference%init_like(vector_mold)
  call fd_base%init_like(vector_mold)
  call fd_perturbed%init_like(vector_mold)
  call fd_result%init_like(vector_mold)
  call fd_result_pert%init_like(vector_mold)
  call fd_difference%init_like(vector_mold)

  call x%rand(.true.)

  ! Verify one LNS step against the finite difference of one nonlinear NS
  ! step. The NS case is reused sequentially so the registry is unchanged.
  call A%set_baseflow(baseflow, fitted_dt, 1)
  call A%matvec(x, Ax)
  call fd_dt_controller%init(base_case%params)
  max_finite_error = 0.0_wp
  do itest = 1, size(epsilon_values)
     epsilon_fd = epsilon_values(itest)

     call reset(base_case, 0, base_case%f_out%file_%get_base_fname())
     call copy_state_to_primal(base_case, baseflow)
     call simulation_step(base_case, fd_dt_controller)
     call copy_primal_to_state(base_case, fd_base)

     call fd_perturbed%axpby(1.0_wp, baseflow, 0.0_wp)
     call fd_perturbed%axpby(epsilon_fd, x, 1.0_wp)
     call reset(base_case, 0, base_case%f_out%file_%get_base_fname())
     call copy_state_to_primal(base_case, fd_perturbed)
     call simulation_step(base_case, fd_dt_controller)
     call copy_primal_to_state(base_case, fd_result_pert)

     call fd_result%axpby(1.0_wp, fd_result_pert, 0.0_wp)
     call fd_result%axpby(-1.0_wp, fd_base, 1.0_wp)
     call fd_result%scal(1.0_wp / epsilon_fd)
     call fd_difference%axpby(1.0_wp, fd_result, 0.0_wp)
     call fd_difference%axpby(-1.0_wp, Ax, 1.0_wp)
     finite_error = sqrt(max(0.0_wp, fd_difference%dot(fd_difference))) / &
          max(sqrt(max(0.0_wp, Ax%dot(Ax))), epsilon(1.0_wp))
     max_finite_error = max(max_finite_error, finite_error)
     write(*,'(A,ES12.4,A,ES24.16)') 'finite-difference epsilon ', &
          epsilon_fd, ' relative error: ', finite_error
  end do

  ! Rebuild the requested two-step coefficient trajectory for the adjoint
  ! identity and checkpoint/replay tests.
  call A%set_baseflow(baseflow, test_interval, test_nsteps)

  max_rel_error = 0.0_wp
  max_repeat_error = 0.0_wp
  do itest = 1, number_random_tests
     call x%rand(.true.)
     call y%rand(.true.)

     call A%matvec(x, Ax)
     call A%matvec(x, Ax_repeat)
     call A%rmatvec(y, Aty)

     lhs = Ax%dot(y)
     rhs = x%dot(Aty)
     abs_error = abs(lhs - rhs)
     rel_error = abs_error / max(abs(lhs), abs(rhs), epsilon(1.0_wp))
     call difference%axpby(1.0_wp, Ax, 0.0_wp)
     call difference%axpby(-1.0_wp, Ax_repeat, 1.0_wp)
     repeat_error = sqrt(max(0.0_wp, difference%dot(difference))) / &
          max(sqrt(max(0.0_wp, Ax%dot(Ax))), epsilon(1.0_wp))

     max_rel_error = max(max_rel_error, rel_error)
     max_repeat_error = max(max_repeat_error, repeat_error)

     write(*,'(A,I0,A,ES24.16)') 'random test ', itest, &
          ' relative adjoint error: ', rel_error
     write(*,'(A,I0,A,ES24.16)') 'random test ', itest, &
          ' repeatability error: ', repeat_error
  end do

  write(*,'(A)') trim(this_module)
  write(*,'(A,I0)') 'full-period fitted steps: ', full_nsteps
  write(*,'(A,I0)') 'tested steps: ', test_nsteps
  write(*,'(A,ES24.16)') 'fitted dt: ', fitted_dt
  write(*,'(A,ES24.16)') 'test interval: ', test_interval
  write(*,'(A,ES24.16)') 'maximum relative adjoint error: ', max_rel_error
  write(*,'(A,ES24.16)') 'maximum repeatability error: ', max_repeat_error
  write(*,'(A,ES24.16)') 'maximum finite-difference error: ', max_finite_error

  open(unit=1234, file='adjoint_consistency.dat', status='replace', &
       action='write')
  write(1234,'(A)') '# full_nsteps test_nsteps fitted_dt test_interval ' // &
       'max_rel_error max_repeat_error max_finite_error tolerance pass'
  write(1234,'(2(I0,1X),6(ES24.16,1X),L1)') full_nsteps, test_nsteps, &
       fitted_dt, test_interval, max_rel_error, max_repeat_error, &
       max_finite_error, tolerance, &
       max_rel_error <= tolerance .and. &
       max_repeat_error <= tolerance .and. max_finite_error <= tolerance
  close(1234)

  if (max_rel_error > tolerance .or. max_repeat_error > tolerance .or. &
       max_finite_error > tolerance) then
     call neko_error('Adjoint consistency check failed')
  end if

  call A%free()
  call Aty%free()
  call Ax%free()
  call Ax_repeat%free()
  call difference%free()
  call fd_base%free()
  call fd_perturbed%free()
  call fd_result%free()
  call fd_result_pert%free()
  call fd_difference%free()
  call y%free()
  call x%free()
  call baseflow%free()
  call adjoint_case%free()
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

end program usrneko
