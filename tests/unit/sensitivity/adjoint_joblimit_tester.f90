!> @file adjoint_joblimit_tester.f90
!! Driver asserting that a steady adjoint whose clock stops short of
!! `end_time` still counts as complete.
!!
!! The steady adjoint marches as many steps as the forward run took, with a
!! time step its own CFL sets, so its clock need not reach `end_time`. With
!! `output_at_end` off, the adjoint's finalisation writes a job-limit
!! checkpoint when it judges the run to have stopped short. A steady adjoint
!! always takes all its steps, so any such checkpoint is spurious.
!!
!! The case file must make the adjoint's clock stop short of `end_time`, or
!! the check cannot bite; the driver stops with an error if it does not. It
!! sets a variable timestep whose `timestep` exceeds the step the CFL target
!! allows once the flow has developed. The forward run starts almost at rest,
!! so its first step is the full `timestep`. The adjoint sizes its steps by
!! the developed final flow, so each of them is shorter.
program adjoint_joblimit_tester

  use simulation_m, only: simulation_t
  use brinkman_design, only: brinkman_design_t
  use problem, only: problem_t

  ! Standard modules shared by most of our tests
  use neko, only: neko_init, neko_finalize
  use json_module, only: json_file
  use json_utils, only: json_get, json_get_or_default
  use json_utils_ext, only: json_read_file
  use utils, only: neko_error
  use neko_top, only: neko_top_register_types
  use logger, only: neko_log

  ! Modules specific to this test
  use num_types, only: dp
  use file, only: file_t
  implicit none

  ! JSON related arguments
  integer :: argc
  character(len=256) :: parameter_file
  type(json_file) :: parameters, design_parameters

  !> The simulation we are working with
  type(simulation_t) :: sim
  !> The design type
  type(brinkman_design_t) :: des
  !> The problem type
  type(problem_t) :: prob

  real(kind=dp) :: end_time, forward_time, adjoint_time
  integer :: adjoint_steps
  logical :: output_at_end, joblimit_written
  !> The job-limit checkpoint the adjoint would write, named as it names it.
  type(file_t) :: joblimit_file
  character(len=1024) :: joblimit_fname
  character(len=:), allocatable :: checkpoint_format
  ! Wider than LOG_SIZE: these lines carry full-precision values.
  character(len=256) :: log_buf

  ! -------------------------------------------------------------------------- !
  ! Initialize the Neko environment

  call neko_init()
  call neko_top_register_types()

  ! -------------------------------------------------------------------------- !
  ! Read the parameters file as the first terminal argument

  argc = command_argument_count()
  if (argc .lt. 1) call neko_error('Missing parameter file')
  call get_command_argument(1, parameter_file)

  parameters = json_read_file(trim(parameter_file))
  call json_get(parameters, 'optimization.design', design_parameters)

  call json_get_or_default(parameters, 'case.output_at_end', output_at_end, &
       .true.)
  if (output_at_end) then
     call neko_error('This test needs "output_at_end": false.')
  end if

  ! -------------------------------------------------------------------------- !
  ! Initialization of the components

  call sim%init(parameters)
  call des%init(design_parameters, sim)
  call prob%init(parameters, des, sim)

  if (sim%unsteady) then
     call neko_error('This test needs a steady problem.')
  end if

  ! Build the name the way `simulation_adjoint_joblimit_chkp` does, so that
  ! this looks for the file it would actually write.
  call json_get_or_default(parameters, 'case.checkpoint_format', &
       checkpoint_format, 'chkp')
  if (checkpoint_format .eq. 'hdf5') then
     call joblimit_file%init(sim%neko_case%output_directory // 'joblimit.h5')
  else
     call joblimit_file%init(sim%neko_case%output_directory // &
          'joblimit.chkp')
  end if
  joblimit_fname = joblimit_file%file_type%get_next_output_fname()
  call joblimit_file%free()

  ! Clear any left over from an earlier run of this test.
  call delete_file(joblimit_fname)

  ! -------------------------------------------------------------------------- !
  ! Solve the forward and the adjoint problem once

  call prob%compute(des, sim)
  call prob%compute_sensitivity(des, sim)

  end_time = sim%neko_case%time%end_time
  forward_time = sim%neko_case%time%t
  adjoint_time = sim%adjoint_case%time%t
  adjoint_steps = sim%adjoint_case%time%tstep

  write(log_buf, '(A,I0,A,E22.15,A,I0,A,E22.15,A,E22.15)') &
       'The forward run took ', sim%n_timesteps, ' steps to t = ', &
       forward_time, ', the adjoint ', adjoint_steps, ' steps to t = ', &
       adjoint_time, ', end_time = ', end_time
  call neko_log%message(trim(log_buf))

  ! -------------------------------------------------------------------------- !
  ! Check the adjoint treated its run as complete

  if (adjoint_steps .ne. sim%n_timesteps) then
     call neko_error('The adjoint did not take as many steps as the ' // &
          'forward run')
  end if
  if (adjoint_time .ge. end_time) then
     call neko_error('The adjoint reached end_time, so this case cannot ' // &
          'test a run that stops short of it; retune the case.')
  end if

  inquire(file = trim(joblimit_fname), exist = joblimit_written)
  if (joblimit_written) then
     call neko_error('The adjoint wrote ' // trim(joblimit_fname) // &
          ', so it judged its complete run to have stopped short')
  end if
  call neko_log%message('The adjoint wrote no job-limit checkpoint')

  ! -------------------------------------------------------------------------- !
  ! Clean up

  call prob%free()
  call des%free()
  call sim%free()

  call neko_finalize()

contains

  !> Delete a file, if it exists.
  !! @param fname Name of the file.
  subroutine delete_file(fname)
    character(len=*), intent(in) :: fname
    logical :: exists
    integer :: unit, stat

    inquire(file = trim(fname), exist = exists)
    if (.not. exists) return

    open(newunit = unit, file = trim(fname), iostat = stat, status = 'old')
    if (stat .eq. 0) close(unit, status = 'delete')
  end subroutine delete_file

end program adjoint_joblimit_tester
