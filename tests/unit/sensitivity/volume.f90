program volume_sensitivity

  use simulation_m, only: simulation_t
  use brinkman_design, only: brinkman_design_t
  use volume_constraint, only: volume_constraint_t

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
  use math, only: abscmp
  use sensitivity, only: compute_sensitivity
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
  !> The volume constraint_object object
  type(volume_constraint_t) :: constraint_object

  ! Test specific variables
  real(kind=rp) :: tolerance = 1e-5_rp
  real(kind=rp), parameter :: perturbations(4) = [ &
       1e-1_rp, 1e-2_rp, 1e-3_rp, 1e-4_rp]

  type(vector_t) :: constraint_sensitivities

  integer :: i_max

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

  ! Initialize our constraint object
  call constraint_object%init_from_components(des, sim, &
       "Volume", "optimization_domain", is_max = .true., limit = 0.2_rp)

  ! -------------------------------------------------------------------------- !
  ! Compute the sensitivity with our method

  call sim%run_forward()
  call constraint_object%update_value(des)
  call sim%run_backward()
  call constraint_object%update_sensitivity(des)
  call sim%reset()
  constraint_sensitivities = constraint_object%get_sensitivity()
  if (NEKO_BCKND_DEVICE .eq. 1) then
     call device_memcpy(constraint_sensitivities%x, &
          constraint_sensitivities%x_d, constraint_sensitivities%size(), &
          DEVICE_TO_HOST, .true.)
  end if

  i_max = maxloc(abs(constraint_sensitivities%x), dim=1)

  ! -------------------------------------------------------------------------- !
  ! Loop over the perturbations and compare the finite difference estimate with
  ! the sensitivity computed by our method.

  call compute_sensitivity(constraint_object, sim, des, &
       constraint_sensitivities, i_max, perturbations, tolerance)

  ! -------------------------------------------------------------------------- !
  ! Clean up the components

  call constraint_object%free()
  call des%free()
  call sim%free()

end program volume_sensitivity
