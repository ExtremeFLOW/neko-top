program checkpointing_test
  use simulation_m, only: simulation_t
  use simulation_checkpoint, only: simulation_checkpoint_t
  use design, only: design_t, design_factory
  use problem, only: problem_t
  use optimizer, only: optimizer_t, optimizer_factory
  use simulation, only: simulation_init, simulation_step, simulation_finalize, &
       simulation_restart

  use num_types, only: dp, rp
  use json_module, only: json_file
  use utils, only: neko_error, neko_warning
  use json_utils_ext, only: json_read_file
  use neko_top, only: neko_top_register_types
  use time_step_controller, only: time_step_controller_t
  use chkp_output, only: chkp_output_t
  use field, only: field_t
  use field_math, only: field_copy, field_glsubnorm, field_cfill, field_glsc2
  use math, only: relcmp, abscmp, NEKO_EPS
  use file, only: file_t
  use comm, only: pe_rank
  use mpi_f08, only: MPI_Init, MPI_Wtime
  use json_utils, only: json_get, json_extract_object
  implicit none

  ! JSON related arguments
  integer :: argc
  type(json_file) :: parameters, design_parameters
  character(len=256) :: parameter_file

  ! MPI parameters
  integer :: ierr

  !> The simulation we are working with
  type(simulation_t) :: sim
  !> The design type
  class(design_t), allocatable :: des
  !> The problem type
  type(problem_t) :: prob
  !> The optimizer (in this case mma)
  class(optimizer_t), allocatable :: opt
  !> The simulation checkpointing object
  type(simulation_checkpoint_t) :: chkp

  !> Log message for errors
  character(len=256) :: log_msg

  integer :: n_timesteps, i
  type(field_t), pointer :: p, u, v, w
  type(field_t), allocatable, dimension(:) :: p_fields, &
       u_fields, v_fields, w_fields
  character(len=256) :: field_name

  ! Objects required for the forward simulation
  type(time_step_controller_t) :: dt_controller
  real(kind=dp) :: loop_start

  ! Objects for the consistency check
  real(kind=rp) :: error_p, error_u, error_v, error_w

  logical :: error

  ! -------------------------------------------------------------------------- !
  ! Initialize the MPI environment

  call MPI_Init(ierr)
  call neko_top_register_types()

  ! -------------------------------------------------------------------------- !
  ! Read the parameters file as the first terminal argument

  argc = command_argument_count()
  if (argc .lt. 1) call neko_error('Missing parameter file')
  call get_command_argument(1, parameter_file)

  ! Read the parameters file
  parameters = json_read_file(trim(parameter_file))
  call json_extract_object(parameters, 'optimization.design', design_parameters)

  ! -------------------------------------------------------------------------- !
  ! Initialization of the components

  call sim%init(parameters)
  call design_factory(des, design_parameters, sim)
  call prob%init(parameters, des, sim)
  call optimizer_factory(opt, parameters, prob, des, sim)

  call chkp%init(parameters, sim%neko_case)

  ! Save some pointers and allocate the fields required for the testing
  p => sim%neko_case%fluid%p
  u => sim%neko_case%fluid%u
  v => sim%neko_case%fluid%v
  w => sim%neko_case%fluid%w

  n_timesteps = int(3.5 * real(chkp%n_saves_memory))
  allocate(p_fields(n_timesteps + 1))
  allocate(u_fields(n_timesteps + 1))
  allocate(v_fields(n_timesteps + 1))
  allocate(w_fields(n_timesteps + 1))

  do i = 1, n_timesteps + 1
     write(field_name, '(A,I0)') 'p_', i
     call p_fields(i)%init(p%dof, field_name)
     write(field_name, '(A,I0)') 'u_', i
     call u_fields(i)%init(u%dof, field_name)
     write(field_name, '(A,I0)') 'v_', i
     call v_fields(i)%init(v%dof, field_name)
     write(field_name, '(A,I0)') 'w_', i
     call w_fields(i)%init(w%dof, field_name)
  end do

  ! -------------------------------------------------------------------------- !
  ! Run the forward simulation and save the resulting u fields in a list

  if (pe_rank .eq. 0) then
     write(*, '(A)') repeat('-', 80)
     write(*, '(A)') 'Running the forward simulation...'
     write(*, '(A,I0)') 'Number of time steps: ', n_timesteps
     write(*, '(A,I0)') 'Number of checkpoints: ', chkp%n_saves_memory
  end if

  call dt_controller%init(sim%neko_case%params)
  call simulation_init(sim%neko_case, dt_controller)

  loop_start = MPI_WTIME()

  chkp%n_saves_disc = 0
  do i = 1, n_timesteps
     call simulation_step(sim%neko_case, dt_controller, loop_start)

     call field_copy(p_fields(i), p)
     call field_copy(u_fields(i), u)
     call field_copy(v_fields(i), v)
     call field_copy(w_fields(i), w)

     call chkp%save(sim%neko_case, sim%neko_case%time)
  end do

  call simulation_finalize(sim%neko_case)

  ! -------------------------------------------------------------------------- !
  ! Restore to each time step and check the consistency of the u fields

  if (pe_rank .eq. 0) then
     write(*, '(A)') repeat('-', 80)
     write(*, '(A)') 'Checking the consistency of the time steps...'
     write(*, '(A,I0)') 'Number of time steps: ', n_timesteps
  end if

  call dt_controller%init(sim%neko_case%params)
  call simulation_init(sim%neko_case, dt_controller)
  chkp%loaded_checkpoint = -1

  do i = n_timesteps, 1, -1
     call chkp%restore(sim%neko_case, i)

     !  Compute the relative error of the fields
     error_p = field_glsubnorm(p_fields(i), p) / &
          max(sqrt(field_glsc2(p_fields(i), p_fields(i))), NEKO_EPS)
     error_u = field_glsubnorm(u_fields(i), u) / &
          max(sqrt(field_glsc2(u_fields(i), u_fields(i))), NEKO_EPS)
     error_v = field_glsubnorm(v_fields(i), v) / &
          max(sqrt(field_glsc2(v_fields(i), v_fields(i))), NEKO_EPS)
     error_w = field_glsubnorm(w_fields(i), w) / &
          max(sqrt(field_glsc2(w_fields(i), w_fields(i))), NEKO_EPS)

     error = .false.

     if (.not. abscmp(error_p, 0.0_rp, NEKO_EPS)) then
        error = .true.
        write(log_msg, '(A,I0,E12.5)') &
             'Inconsistency found in time step pressure: ', i, error_p
        call neko_warning(trim(log_msg))
     end if

     if (.not. abscmp(error_u, 0.0_rp, NEKO_EPS)) then
        error = .true.
        write(log_msg, '(A,I0,E12.5)') &
             'Inconsistency found in time step velocity u: ', i, error_u
        call neko_warning(trim(log_msg))
     end if

     if (.not. abscmp(error_v, 0.0_rp, NEKO_EPS)) then
        error = .true.
        write(log_msg, '(A,I0,E12.5)') &
             'Inconsistency found in time step velocity v: ', i, error_v
        call neko_warning(trim(log_msg))
     end if

     if (.not. abscmp(error_w, 0.0_rp, NEKO_EPS)) then
        error = .true.
        write(log_msg, '(A,I0,E12.5)') &
             'Inconsistency found in time step velocity w: ', i, error_w
        call neko_warning(trim(log_msg))
     end if

     if (error) then
        if (pe_rank .eq. 0) then
           write(*, '(A,E12.5)') '    Error for pressure:   ', error_p
           write(*, '(A,E12.5)') '    Error for velocity u: ', error_u
           write(*, '(A,E12.5)') '    Error for velocity v: ', error_v
           write(*, '(A,E12.5)') '    Error for velocity w: ', error_w
        end if
        call neko_error('Inconsistency found in the time step')
     end if
  end do

  if (pe_rank .eq. 0) then
     write(*, '(A)') 'All time steps are consistent.'
     write(*, '(A)') repeat('-', 80)
  end if

  ! -------------------------------------------------------------------------- !
  ! Clean up the components

  do i = 1, n_timesteps
     call p_fields(i)%free()
     call u_fields(i)%free()
     call v_fields(i)%free()
     call w_fields(i)%free()
  end do
  deallocate(p_fields)
  deallocate(u_fields)
  deallocate(v_fields)
  deallocate(w_fields)

  call opt%free()
  call prob%free()
  call des%free()
  call sim%free()

  if (allocated(des)) deallocate(des)
  if (allocated(opt)) deallocate(opt)

end program checkpointing_test
