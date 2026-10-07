program usrneko
  use LightKrylov
  use LightKrylov, only: wp => dp
  use LightKrylov_Constants
  use LightKrylov_Logger
  use json_module, only: json_file
  use json_utils, only: json_get_or_default
  use json_utils_ext, only: json_read_file
  use neko, only: neko_finalize, neko_init
  use case, only: case_t
  use neko_complex_vector, only: complex_state_vector_t
  use neko_config, only: NEKO_BCKND_DEVICE
  use neko_resolvent, only: resolvent_shifted_operator_t
  use neko_top, only: neko_top_register_types
  use neko_vector, only: state_vector_t
  use field_math, only: field_copy
  use field_output, only: field_output_t
  use fld_file, only: fld_file_t
  use num_types, only: neko_rp => rp, neko_sp => sp
  use perturb_case, only: perturb_case_t
  use stability_case_utils, only: stability_case_free, stability_case_init
  use utils, only: neko_error
  implicit none

  character(len=*), parameter :: this_module = 'Cylinder resolvent'
  type(case_t), target :: base_case
  type(perturb_case_t), target :: lns_case
  type(resolvent_shifted_operator_t) :: shifted
  type(state_vector_t) :: real_mold, baseflow
  type(complex_state_vector_t) :: complex_mold, forcing, response
  real(wp) :: omega, solver_tol, gain
  integer :: argc, info
  character(len=256) :: parameter_file
  type(json_file) :: parameters

  call neko_init()
  call neko_top_register_types()
  call logger_setup()

  argc = command_argument_count()
  if (argc .lt. 1) call neko_error('Missing parameter file')
  call get_command_argument(1, parameter_file)

  parameters = json_read_file(trim(parameter_file))
  call json_get_or_default(parameters, &
       'case.lightkrylov.resolvent.omega', omega, 1.0_wp)
  call json_get_or_default(parameters, &
       'case.lightkrylov.resolvent.tolerance', solver_tol, 1.0e-8_wp)

  if (solver_tol <= 0.0_wp) call neko_error('Resolvent tolerance must be positive')

  call stability_case_init(base_case, parameters)
  lns_case%if_adjoint = .false.
  call lns_case%init(base_case)

  real_mold%coef => base_case%fluid%c_Xh
  real_mold%if_2d = (base_case%fluid%c_Xh%msh%gdim .eq. 2)
  call baseflow%init_like(real_mold)
  call field_copy(baseflow%u, base_case%fluid%u)
  call field_copy(baseflow%v, base_case%fluid%v)
  call field_copy(baseflow%w, base_case%fluid%w)
  call write_state_fields('baseflow_resolvent', baseflow, 0)

  call complex_mold%re%init_like(real_mold)
  call complex_mold%im%init_like(real_mold)

  call forcing%init_like(complex_mold)
  call response%init_like(complex_mold)
  call forcing%zero()
  call response%zero()

  call make_actuator(forcing)
  call write_state_fields('forcing_resolvent_re', forcing%re, 0)
  call write_state_fields('forcing_resolvent_im', forcing%im, 0)

  call shifted%init(base_case, lns_case, omega)
  call gmres(shifted, forcing, response, info, rtol=solver_tol)
  call check_info(info, 'gmres', module=this_module, procedure='main')

  call response%write(0)
  call write_state_fields('response_resolvent_re', response%re, 0)
  call write_state_fields('response_resolvent_im', response%im, 0)
  gain = response%norm()
  open(unit=1234, file='resolvent.txt', status='replace', action='write')
  write(1234, *) omega, gain * gain
  close(1234)

  call shifted%free()
  call forcing%free()
  call response%free()
  call complex_mold%free()
  call baseflow%free()
  call lns_case%free()
  call stability_case_free(base_case)
  call neko_finalize()

contains

  subroutine make_actuator(vec)
    type(complex_state_vector_t), intent(inout) :: vec
    real(neko_rp), parameter :: sigma = 0.1_neko_rp
    real(neko_rp), parameter :: radius = 0.6_neko_rp
    real(neko_rp), parameter :: amplitude = 1.0_neko_rp
    real(neko_rp), parameter :: pi_ = 4.0_neko_rp * atan(1.0_neko_rp)
    real(neko_rp), parameter :: theta = 70.0_neko_rp * pi_ / 180.0_neko_rp
    real(neko_rp) :: sfield1, sfield2, x, y
    integer :: iel, ix, iy, iz

    call vec%zero()
    do iel = 1, vec%re%u%msh%nelv
       do iz = 1, vec%re%u%Xh%lz
          do iy = 1, vec%re%u%Xh%ly
             do ix = 1, vec%re%u%Xh%lx
                x = vec%re%u%dof%x(ix, iy, iz, iel)
                y = vec%re%u%dof%y(ix, iy, iz, iel)
                sfield1 = amplitude / (2.0_neko_rp * pi_ * sigma**2) * &
                     exp(-((x - radius * cos(theta))**2 + &
                     (y - radius * sin(theta))**2) / (2.0_neko_rp * sigma**2))
                sfield2 = amplitude / (2.0_neko_rp * pi_ * sigma**2) * &
                     exp(-((x - radius * cos(-theta))**2 + &
                     (y - radius * sin(-theta))**2) / (2.0_neko_rp * sigma**2))

                vec%re%u%x(ix, iy, iz, iel) = cos(theta) * (sfield1 - sfield2)
                vec%re%v%x(ix, iy, iz, iel) = sin(theta) * (sfield1 + sfield2)
                vec%re%w%x(ix, iy, iz, iel) = 0.0_neko_rp
             end do
          end do
       end do
    end do

    call sync_state(vec%re)
  end subroutine make_actuator

  subroutine sync_state(vec)
    use device, only: device_memcpy, HOST_TO_DEVICE
    type(state_vector_t), intent(inout) :: vec

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(vec%u%x, vec%u%x_d, vec%u%size(), HOST_TO_DEVICE, .false.)
       call device_memcpy(vec%v%x, vec%v%x_d, vec%v%size(), HOST_TO_DEVICE, .false.)
       call device_memcpy(vec%w%x, vec%w%x_d, vec%w%size(), HOST_TO_DEVICE, .true.)
    end if
  end subroutine sync_state

  subroutine write_state_fields(name, vec, idx)
    character(len=*), intent(in) :: name
    type(state_vector_t), intent(inout) :: vec
    integer, intent(in) :: idx
    type(field_output_t) :: output

    call execute_command_line('mkdir -p viz')
    call output%init(trim(name), 3, precision=neko_sp, path='viz/', format='nek5000')
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

end program usrneko
