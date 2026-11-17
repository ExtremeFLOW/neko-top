!> @file steady_tester.f90
!! Driver running a steady problem once, through `problem_t%compute`.
!!
!! It checks nothing itself. The cases it runs carry a `steady` simulation
!! component, and the checks read what that component wrote to the log; see
!! `check_steady_log.cmake` and the readme in this directory. It differs from
!! the plain Neko-TOP driver only in registering the inflow of
!! `objectives_user` and in building the design and problem, so that the run
!! is the forward solve of a steady optimisation problem.
program steady_tester

  use simulation_m, only: simulation_t
  use brinkman_design, only: brinkman_design_t
  use problem, only: problem_t

  ! Standard modules shared by most of our tests
  use neko, only: neko_init, neko_finalize
  use json_module, only: json_file
  use json_utils, only: json_get
  use json_utils_ext, only: json_read_file
  use utils, only: neko_error
  use neko_top, only: neko_top_register_types

  ! Modules specific to this test
  use objectives_user, only: objectives_user_setup
  implicit none

  integer :: argc
  character(len=256) :: parameter_file
  type(json_file) :: parameters, design_parameters

  !> The simulation we are working with
  type(simulation_t) :: sim
  !> The design type
  type(brinkman_design_t) :: des
  !> The problem type
  type(problem_t) :: prob

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

  ! -------------------------------------------------------------------------- !
  ! Initialization of the components

  ! Register the user-defined inflow before the case is built.
  call objectives_user_setup(sim%neko_case%user)

  call sim%init(parameters)
  call des%init(design_parameters, sim)
  call prob%init(parameters, des, sim)

  if (sim%unsteady) then
     call neko_error('This driver runs steady problems only; ' // &
          'set "unsteady": false in the case file.')
  end if

  ! -------------------------------------------------------------------------- !
  ! Run the forward problem once

  call prob%compute(des, sim)

  ! -------------------------------------------------------------------------- !
  ! Clean up

  call prob%free()
  call des%free()
  call sim%free()

  call neko_finalize()

end program steady_tester
