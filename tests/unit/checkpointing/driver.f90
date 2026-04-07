module checkpointing_test_mod
  use field, only: field_t
  use field_math, only: field_glsubnorm
  use num_types, only: rp
  implicit none
  private
  public :: rmse

contains

  function rmse(field1, field2) result(result)
    type(field_t), intent(in) :: field1, field2
    real(kind=rp) :: result

    result = sqrt(field_glsubnorm(field1, field2)**2 / real(field1%size(), rp))

  end function rmse

end module checkpointing_test_mod

program checkpointing_test
  use simulation_m, only: simulation_t
  use simulation_checkpoint, only: simulation_checkpoint_t
  use design, only: design_t, design_factory
  use problem, only: problem_t
  use simulation, only: simulation_init, simulation_step, simulation_finalize

  use neko, only: neko_init, neko_finalize
  use num_types, only: dp, rp
  use json_module, only: json_file
  use utils, only: neko_error, neko_warning
  use json_utils, only: json_get
  use json_utils_ext, only: json_read_file
  use neko_top, only: neko_top_register_types
  use time_step_controller, only: time_step_controller_t
  use field, only: field_t
  use field_math, only: field_copy, field_glsubnorm, field_glsc2
  use math, only: abscmp, NEKO_EPS
  use comm, only: pe_rank
  use mpi_f08, only: MPI_Init, MPI_Wtime
  use, intrinsic :: iso_fortran_env, only: stderr => error_unit
  use checkpointing_test_mod, only: rmse
  implicit none

  ! JSON related arguments
  integer :: argc
  type(json_file) :: parameters, design_parameters
  character(len=256) :: parameter_file

  !> The simulation we are working with
  type(simulation_t) :: sim
  !> The design type
  class(design_t), allocatable :: des
  !> The problem type
  type(problem_t) :: prob
  !> The simulation checkpointing object
  type(simulation_checkpoint_t) :: chkp

  ! Parameters for the checkpointing
  integer :: n_saves_memory = 10
  character(len=256) :: algorithm = "linear"
  character(len=256) :: filename = "checkpoint"
  character(len=8) :: fmt = "chkp"
  logical :: keep_checkpoints = .false.

  !> Log message for errors
  ! character(len=256) :: log_msg

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

  ! Tolerances for the consistency check
  real(kind=rp), parameter :: tol_vel = NEKO_EPS
  real(kind=rp), parameter :: tol_p = NEKO_EPS

  ! -------------------------------------------------------------------------- !
  ! Initialize the Neko environment

  call neko_init()
  call neko_top_register_types()

  ! -------------------------------------------------------------------------- !
  ! Read the parameters file as the first terminal argument

  argc = command_argument_count()
  if (argc .lt. 1) call neko_error('Missing parameter file')
  call get_command_argument(1, parameter_file)

  ! Read the parameters file
  parameters = json_read_file(trim(parameter_file))
  call json_get(parameters, 'optimization.design', design_parameters)

  ! -------------------------------------------------------------------------- !
  ! Initialization of the components

  call sim%init(parameters)
  call design_factory(des, design_parameters, sim)
  call prob%init(parameters, des, sim)

  ! Save some pointers and allocate the fields required for the testing
  p => sim%neko_case%fluid%p
  u => sim%neko_case%fluid%u
  v => sim%neko_case%fluid%v
  w => sim%neko_case%fluid%w

  n_timesteps = int(5.5 * real(n_saves_memory))
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
  ! Initialize the checkpointing
  if (pe_rank .eq. 0) then
     write(*, '(A)') repeat('-', 80)
     write(*, '(A, A)') 'Checkpointing algorithm:   ', trim(algorithm)
     write(*, '(A, A)') 'Checkpointing file format: ', trim(fmt)
     write(*, '(A,I0)') 'Number of time steps:      ', n_timesteps
     write(*, '(A,I0)') 'Number of memory saves:    ', n_saves_memory
     write(*, '(A)') repeat('-', 80)
  end if

  call chkp%init(sim%neko_case, algorithm, n_saves_memory = n_saves_memory, &
       filename = filename, fmt = fmt, keep_checkpoints = keep_checkpoints)

  ! -------------------------------------------------------------------------- !
  ! Run the forward simulation and save the resulting u fields in a list

  call dt_controller%init(sim%neko_case%params)
  call simulation_init(sim%neko_case, dt_controller)

  loop_start = MPI_Wtime()

  do i = 1, n_timesteps
     call simulation_step(sim%neko_case, dt_controller, loop_start)

     call field_copy(p_fields(i), p)
     call field_copy(u_fields(i), u)
     call field_copy(v_fields(i), v)
     call field_copy(w_fields(i), w)

     call chkp%save(sim%neko_case)
  end do

  call simulation_finalize(sim%neko_case)

  ! -------------------------------------------------------------------------- !
  ! Restore to each time step and check the consistency of the u fields

  if (pe_rank .eq. 0) then
     write(*, '(A)') 'Checking the consistency of the time steps...'
  end if

  call dt_controller%init(sim%neko_case%params)
  call simulation_init(sim%neko_case, dt_controller)

  error = .false.

  do i = n_timesteps, 1, -1
     call chkp%restore(sim%neko_case, i)

     !  Compute the RMSE between the current fields and the saved ones
     error_p = rmse(p_fields(i), p)
     error_u = rmse(u_fields(i), u)
     error_v = rmse(v_fields(i), v)
     error_w = rmse(w_fields(i), w)

     if (.not. abscmp(error_p, 0.0_rp, tol_p)) error = .true.
     if (.not. abscmp(error_u, 0.0_rp, tol_vel)) error = .true.
     if (.not. abscmp(error_v, 0.0_rp, tol_vel)) error = .true.
     if (.not. abscmp(error_w, 0.0_rp, tol_vel)) error = .true.
     if (error) exit
  end do

  call simulation_finalize(sim%neko_case)

  if (error) then
     if (pe_rank .eq. 0) then
        write(stderr, '(A,I0)') 'Inconsistency found at time step: ', i
        write(stderr, '(A,E12.5)') '    Error for pressure:   ', error_p
        write(stderr, '(A,E12.5)') '    Error for velocity u: ', error_u
        write(stderr, '(A,E12.5)') '    Error for velocity v: ', error_v
        write(stderr, '(A,E12.5)') '    Error for velocity w: ', error_w
     end if
     error stop 'Checkpointing consistency test failed.'
  else if (pe_rank .eq. 0) then
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

  call chkp%free()
  call sim%free()
  call des%free()
  call prob%free()

  if (allocated(des)) deallocate(des)

  ! Finalize the Neko environment
  call neko_finalize()

end program checkpointing_test
