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
  use case, only: case_t
  use linop_transient_baseflow, only: linop_transient_baseflow_t
  use neko_top, only: neko_top_register_types
  use neko_vector, only: state_vector_t
  use perturb_case, only: perturb_case_t
  use stability_case_utils, only: stability_case_free, stability_case_init
  use utils, only: extract_fld_file_index, filename_chsuffix, neko_error
  implicit none

  character(len=*), parameter :: this_module = 'Cylinder adjoint consistency'
  type(case_t), target :: base_case
  type(perturb_case_t), target :: lns_case, adjoint_case
  type(linop_transient_baseflow_t) :: A
  type(state_vector_t) :: vector_mold, baseflow, x, y, Ax, Aty
  real(kind=wp) :: period, old_dt, fitted_dt, test_interval
  real(kind=wp) :: lhs, rhs, abs_error, rel_error, tolerance
  integer :: argc, full_nsteps, test_nsteps
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
       'case.lightkrylov.adjoint_test.tolerance', tolerance, 1.0e-10_wp)

  call stability_case_init(base_case, parameters)
  if (len_trim(baseflow_file) == 0) &
       call neko_error('Adjoint consistency baseflow file_name is required')

  period = read_field_time(baseflow_file)
  if (period <= 0.0_wp) &
       call neko_error('Adjoint consistency period must be read from baseflow')
  if (test_nsteps < 1) &
       call neko_error('Adjoint consistency number_of_timesteps must be positive')
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

  call x%init_like(vector_mold)
  call y%init_like(vector_mold)
  call x%rand(.true.)
  call y%rand(.true.)

  call A%init(base_case, lns_case, adjoint_case)
  call A%set_baseflow(baseflow, test_interval, test_nsteps)

  call A%matvec(x, Ax)
  call A%rmatvec(y, Aty)

  lhs = Ax%dot(y)
  rhs = x%dot(Aty)
  abs_error = abs(lhs - rhs)
  rel_error = abs_error / max(abs(lhs), abs(rhs), epsilon(1.0_wp))

  write(*,'(A)') trim(this_module)
  write(*,'(A,I0)') 'full-period fitted steps: ', full_nsteps
  write(*,'(A,I0)') 'tested steps: ', test_nsteps
  write(*,'(A,ES24.16)') 'fitted dt: ', fitted_dt
  write(*,'(A,ES24.16)') 'test interval: ', test_interval
  write(*,'(A,ES24.16)') '<A x, y>: ', lhs
  write(*,'(A,ES24.16)') '<x, A* y>: ', rhs
  write(*,'(A,ES24.16)') 'absolute error: ', abs_error
  write(*,'(A,ES24.16)') 'relative error: ', rel_error

  open(unit=1234, file='adjoint_consistency.dat', status='replace', &
       action='write')
  write(1234,'(A)') '# full_nsteps test_nsteps fitted_dt test_interval ' // &
       'lhs rhs abs_error rel_error tolerance pass'
  write(1234,'(2(I0,1X),6(ES24.16,1X),L1)') full_nsteps, test_nsteps, &
       fitted_dt, test_interval, lhs, rhs, abs_error, rel_error, &
       rel_error <= tolerance
  close(1234)

  if (rel_error > tolerance) then
     call neko_error('Adjoint consistency check failed')
  end if

  call A%free()
  call Aty%free()
  call Ax%free()
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
