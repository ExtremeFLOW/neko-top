program checkpointing_test
  use simulation_m, only: simulation_t
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
  use math, only: relcmp, abscmp
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

  !> Log message for errors
  character(len=256) :: log_msg

  integer :: n_timesteps, i,j
  type(field_t), pointer :: p, u, v, w
  type(field_t), allocatable, dimension(:) :: p_fields, &
       u_fields, v_fields, w_fields
  character(len=256) :: field_name

  ! Objects required for the forward simulation
  type(time_step_controller_t) :: dt_controller
  real(kind=dp) :: loop_start
  character(len=256) :: chkp_file_name

  ! Objects for the consistency check
  type(file_t) :: chkp_file
  real(kind=rp) :: norm_ref_p, norm_ref_vel, norm_diff_p, norm_diff_vel
  real(kind=rp) :: p_tol = 1.0e-9_rp
  real(kind=rp) :: vel_tol = 1.0e-9_rp

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

  ! Save some pointers and allocate the fields required for the testing
  p => sim%neko_case%fluid%p
  u => sim%neko_case%fluid%u
  v => sim%neko_case%fluid%v
  w => sim%neko_case%fluid%w

  sim%first_valid_timestep = 2

  n_timesteps = int(3.5 * real(sim%n_saves_memory))
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
     write(*, '(A,I0)') 'Number of checkpoints: ', sim%n_saves_memory
  end if

  call dt_controller%init(sim%neko_case%params)
  call simulation_init(sim%neko_case, dt_controller)

  loop_start = MPI_WTIME()

  sim%n_saves_disc = 0
  do i = 1, n_timesteps
     call simulation_step(sim%neko_case, dt_controller, loop_start)

     call field_copy(p_fields(i), p)
     call field_copy(u_fields(i), u)
     call field_copy(v_fields(i), v)
     call field_copy(w_fields(i), w)

     call sim%save_state(sim%neko_case%time)
  end do

  call simulation_finalize(sim%neko_case)

  ! -------------------------------------------------------------------------- !
  ! Load each checkpoint and check the consistency of the u fields

  if (pe_rank .eq. 0) then
     write(*, '(A)') repeat('-', 80)
     write(*, '(A)') 'Checking the consistency of the save states...'
     write(*, '(A,I0)') 'Number of save states: ', sim%n_saves_disc
  end if
  if (sim%n_saves_disc .eq. 0) then
     call neko_error('No save states found.')
  end if

  call dt_controller%init(sim%neko_case%params)
  call simulation_init(sim%neko_case, dt_controller)

  do i = 1, sim%n_saves_disc
     write(chkp_file_name, '(A,I5.5,A)') 'forward_chkp_', i - 1, '.chkp'
     call chkp_file%init(chkp_file_name)
     call chkp_file%read(sim%neko_case%chkp)
     call simulation_restart(sim%neko_case, sim%neko_case%chkp)
     j = sim%neko_case%time%tstep

     ! Compute the l2 norm of the u field and the original one
     norm_ref_p = sqrt(field_glsc2(p_fields(j), p_fields(j)))
     norm_ref_vel = (sqrt(field_glsc2(u_fields(j), u_fields(j))) &
          + sqrt(field_glsc2(v_fields(j), v_fields(j))) &
          + sqrt(field_glsc2(w_fields(j), w_fields(j)))) / 3.0_rp
     norm_diff_p = field_glsubnorm(p_fields(j), p) / norm_ref_p
     norm_diff_vel = (field_glsubnorm(u_fields(j), u) &
          + field_glsubnorm(v_fields(j), v) &
          + field_glsubnorm(w_fields(j), w)) / (3.0_rp * norm_ref_vel)

     if (.not. abscmp(norm_diff_p, 0.0_rp, p_tol)) then
        write(log_msg, '(A,A,E12.5)') &
             'Inconsistency found in save state pressure: ', &
             trim(chkp_file_name), norm_diff_p
        call neko_error(trim(log_msg))
     else if (.not. abscmp(norm_diff_vel, 0.0_rp, vel_tol)) then
        write(log_msg, '(A,A,E12.5)') &
             'Inconsistency found in save state velocity: ', &
             trim(chkp_file_name), norm_diff_vel
        call neko_error(trim(log_msg))
     else
        if (pe_rank .eq. 0) then
           write(*, '(A,A,A)') 'Save state ', trim(chkp_file_name), &
                ' is consistent'
        end if
     end if
  end do

  ! -------------------------------------------------------------------------- !
  ! Restore to each time step and check the consistency of the u fields

  if (pe_rank .eq. 0) then
     write(*, '(A)') repeat('-', 80)
     write(*, '(A)') 'Checking the consistency of the time steps...'
     write(*, '(A,I0)') 'Number of time steps: ', n_timesteps
  end if

  call dt_controller%init(sim%neko_case%params)
  call simulation_init(sim%neko_case, dt_controller)
  sim%loaded_checkpoint = -1

  do i = n_timesteps, sim%first_valid_timestep, -1
     ! do i = sim%first_valid_timestep, n_timesteps
     if (pe_rank .eq. 0) write(*, '(A,I0)') 'Checking time step ', i
     call sim%restore_state(i)

     ! Compute the l2 norm of the u field and the original one
     norm_ref_p = sqrt(field_glsc2(p_fields(i), p_fields(i)))
     norm_ref_vel = (sqrt(field_glsc2(u_fields(i), u_fields(i))) &
          + sqrt(field_glsc2(v_fields(i), v_fields(i))) &
          + sqrt(field_glsc2(w_fields(i), w_fields(i)))) / 3.0_rp
     norm_diff_p = field_glsubnorm(p_fields(i), p)
     norm_diff_vel = (field_glsubnorm(u_fields(i), u) &
          + field_glsubnorm(v_fields(i), v) &
          + field_glsubnorm(w_fields(i), w)) / 3.0_rp


     if (.not. abscmp(norm_diff_p / norm_ref_p, 0.0_rp, p_tol)) then
        write(log_msg, '(A,I0,E12.5)') &
             'Inconsistency found in time step pressure: ', i, norm_diff_p

        if (pe_rank .eq. 0) then
           write(*, '(A,E12.5)') 'Norm difference in pressure: ', norm_diff_p
           write(*, '(A,E12.5)') 'Norm reference pressure: ', norm_ref_p
           write(*, '(A,E12.5)') 'Relative difference: ', &
                norm_diff_p / norm_ref_p
        end if

        call neko_error(trim(log_msg))
     else if (.not. abscmp(norm_diff_vel / norm_ref_vel, 0.0_rp, vel_tol)) then

        write(log_msg, '(A,I0,E12.5)') &
             'Inconsistency found in time step velocity: ', i, norm_diff_vel
        call neko_warning(trim(log_msg))

        if (pe_rank .eq. 0) then
           write(*, '(A,E12.5)') 'Norm difference in velocity: ', norm_diff_vel
           write(*, '(A,E12.5)') 'Norm reference velocity: ', norm_ref_vel
           write(*, '(A,E12.5)') 'Relative difference: ', &
                norm_diff_p / norm_ref_p
        end if

        call neko_error(trim(log_msg))
     else
        if (pe_rank .eq. 0) then
           write(*, '(A,I0,A)') 'Timestep ', i, ' is consistent'
        end if
     end if
  end do

  ! -------------------------------------------------------------------------- !
  ! Clean up the components

  do i = 1, n_timesteps
     call p_fields(i)%free()
     call u_fields(i)%free()
  end do
  deallocate(p_fields)
  deallocate(u_fields)

  call opt%free()
  call prob%free()
  call des%free()
  call sim%free()

  if (allocated(des)) deallocate(des)
  if (allocated(opt)) deallocate(opt)

end program checkpointing_test
