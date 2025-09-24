program minimum_dissipation_sensitivity

  use simulation_m, only: simulation_t
  use brinkman_design, only: brinkman_design_t
  use problem, only : problem_t

  ! Standard modules shared by most of our tests
  use json_module, only: json_file
  use json_utils, only: json_extract_object
  use json_utils_ext, only: json_read_file
  use utils, only: neko_error
  use neko_top, only: neko_top_register_types
  use mpi_f08, only: MPI_Init
  use mask_ops, only: mask_exterior_const
  use neko_config, only: NEKO_BCKND_DEVICE
  use device, only: device_memcpy, DEVICE_TO_HOST

  ! Modules specific to this test
  use num_types, only: rp
  use vector, only: vector_t
  use matrix, only: matrix_t
  use math, only: abscmp, copy
  use sensitivity, only: compute_sensitivity
  use user, only: user_setup
  implicit none

  ! JSON related arguments
  integer :: argc
  character(len=256) :: parameter_file
  type(json_file) :: parameters, design_parameters

  ! MPI parameters
  integer :: ierr

  !> The simulation we are working with
  type(simulation_t) :: sim
  !> The design type
  type(brinkman_design_t) :: des
  !> The problem type
  type(problem_t) :: prob

  ! Test specific variables
  real(kind=rp) :: tolerance = 1e-5_rp
  real(kind=rp), parameter :: perturbations(8) = [ &
       5e-3_rp, 1e-3_rp, 5e-4_rp, 1e-4_rp, 5e-5_rp, 1e-5_rp, 5e-6_rp, 1e-6_rp]

  type(vector_t) :: sensitivities
  type(matrix_t) :: constraint_sensitivity

  integer :: i_max

  ! True => testing an objective, F => testing a constraint
  logical :: is_objective
  character(len=12) :: nobj_str, ncon_str

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
  call des%init(design_parameters, sim)
  call prob%init(parameters, des, sim)

  ! -------------------------------------------------------------------------- !
  ! Determine if objective or constraint
  if ((prob%get_n_objectives() .gt. 0) .and. &
       (prob%get_n_constraints() .eq. 0)) then
     is_objective = .true.
  else if (prob%get_n_constraints() .eq. 1) then
     ! note we always have a dummy objective
     is_objective = .false.
  else
     write(nobj_str, '(I0)') prob%get_n_objectives()
     write(ncon_str, '(I0)') prob%get_n_constraints()
     call neko_error("Specify a) a single constraint b) multiple " // &
          "objectives. You have" // nobj_str // &
          "objectives and " // ncon_str // " constraints.")
  end if

  ! -------------------------------------------------------------------------- !
  ! Compute the sensitivity with our method

  call prob%compute(des, sim)
  call prob%compute_sensitivity(des, sim)
  if (is_objective) then
     call des%get_sensitivity(sensitivities)
  else
     call prob%get_constraint_sensitivities(constraint_sensitivity)
     call sensitivities%init(constraint_sensitivity%size())
     call copy(sensitivities%x, constraint_sensitivity%x, &
          constraint_sensitivity%size())
  end if

  call des%write(1)
  ! --------------------------------------
  ! Reset the simulation
  call sim%reset()

  if (NEKO_BCKND_DEVICE .eq. 1) then
     call device_memcpy(sensitivities%x, &
          sensitivities%x_d, sensitivities%size(), &
          DEVICE_TO_HOST, .true.)
  end if

  i_max = maxloc(abs(sensitivities%x), dim=1)

  ! -------------------------------------------------------------------------- !
  ! Loop over the perturbations and compare the finite difference estimate with
  ! the sensitivity computed by our method.

  call compute_sensitivity(prob, sim, des, sensitivities, &
       i_max, perturbations, tolerance, trim(parameter_file), is_objective)

  ! -------------------------------------------------------------------------- !
  ! Clean up the components

  call prob%free()
  call des%free()
  call sim%free()

end program minimum_dissipation_sensitivity
