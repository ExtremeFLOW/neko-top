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

  ! Modules specific to this test
  use num_types, only: rp
  use vector, only: vector_t
  use matrix, only: matrix_t
  use math, only: abscmp
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
  real(kind=rp) :: tolerance = 1e-6_rp
  integer, parameter :: n_perturbations = 6
  real(kind=rp) :: perturbations(n_perturbations) = [ &
       1e-1_rp, 1e-2_rp, 1e-3_rp, 1e-4_rp, 1e-5_rp, 1e-6_rp]

  real(kind=rp) :: fd_estimate, fd_error
  real(kind=rp) :: perturb
  type(vector_t) :: design_vector
  type(vector_t) :: design_perturbed

  real(kind=rp) :: constraint, perturbed_constraint
  type(vector_t) :: constraint_sensitivities

  character(len=*), parameter :: fmt_header = '(4X,A12,4X,A10,6X,A11,5X,A5,10X)'
  character(len=*), parameter :: fmt_data = '(4X,4E15.6E3)'

  integer :: ip, i_max, i, n_elements

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

  constraint = constraint_object%get_value()
  constraint_sensitivities = constraint_object%get_sensitivity()

  call sim%reset()

  ! -------------------------------------------------------------------------- !
  ! Loop over the perturbations and compare the finite difference estimate with
  ! the sensitivity computed by our method.

  n_elements = des%size()

  ! Determine the largest sensitivity to perturb the design variable
  i_max = maxloc(abs(constraint_sensitivities%x), dim=1)

  ! Get the design vector for reference
  ! This is the design vector we will perturb
  design_vector = des%get_values()

  ! do i = 1, n_elements,

  i = i_max

  write(*, '(I0,1X,A,F10.6,1X,A,F10.6,F10.6,F10.6,A)') &
       i, 'Design variable ', design_vector%x(i), &
       'Location [', des%get_x(i), des%get_y(i), des%get_z(i), ']'
  write(*, fmt_header) "Perturbation", "Constraint", "FD Estimate", "Error"
  write(*, fmt_data) 0.0_rp, constraint, constraint_sensitivities%x(i), 0.0_rp

  do ip = 1, n_perturbations
     perturb = perturbations(ip)

     ! Ensure the perturbation stays within the bounds of the design variable
     if (design_vector%x(i) .gt. 0.5_rp) perturb = -perturb

     ! Reset and Perturb the design field by a small amount
     design_perturbed = design_vector
     design_perturbed%x(i) = design_vector%x(i) + perturb
     call des%update_design(design_perturbed)

     ! Compute the objective value of the perturbed design
     call sim%run_forward()
     call constraint_object%update_value(des)
     perturbed_constraint = constraint_object%get_value()
     call sim%reset()

     fd_estimate = perturbed_constraint - constraint
     if (.not. abscmp(fd_estimate, 0.0_rp)) fd_estimate = fd_estimate / perturb

     fd_error = (fd_estimate - constraint_sensitivities%x(i)) / &
          constraint_sensitivities%x(i)

     write(*, fmt_data) perturb, perturbed_constraint, fd_estimate, fd_error

     if (abs(fd_error) .gt. tolerance) then
        call neko_error('Finite difference estimate does not match sensitivity')
     end if
  end do
  ! end do
  ! -------------------------------------------------------------------------- !
  ! Clean up the components

  call constraint_object%free()
  call des%free()
  call sim%free()

end program volume_sensitivity
