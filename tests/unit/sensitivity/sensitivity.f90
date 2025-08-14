module sensitivity
  use simulation_m, only: simulation_t
  use design, only: design_t
  use base_functional, only: base_functional_t
  use utils, only: neko_error
  use num_types, only: rp
  use math, only: abscmp
  use vector, only: vector_t
  use neko_config, only: NEKO_BCKND_DEVICE
  use problem, only : problem_t
  use device, only: device_memcpy, DEVICE_TO_HOST, HOST_TO_DEVICE
  use csv_file, only : csv_file_t
  implicit none

  interface compute_sensitivity
     module procedure compute_sensitivity_i, compute_sensitivity_list, &
        compute_problem_sensitivity_i, compute_problem_sensitivity_list
  end interface compute_sensitivity

contains

  subroutine compute_sensitivity_i(object, sim, des, target_sensitivities, i, &
       perturbations, tolerance)
    class(base_functional_t), intent(inout) :: object
    type(simulation_t), intent(inout) :: sim
    class(design_t), intent(inout) :: des
    type(vector_t), intent(inout) :: target_sensitivities
    integer, intent(in) :: i
    real(kind=rp), intent(in) :: perturbations(:)
    real(kind=rp), intent(in) :: tolerance

    character(len=*), parameter :: fmt_head = '(4X,A12,4X,A10,6X,A11,5X,A5,10X)'
    character(len=*), parameter :: fmt_data = '(4X,4E15.6E3)'

    integer :: n_perturbations, ip
    real(kind=rp) :: perturb, tol
    type(vector_t) :: design_vector, design_perturbed
    real(kind=rp) :: constraint, perturbed_constraint
    real(kind=rp) :: fd_estimate, fd_error

    ! Get the design vector for reference
    ! This is the design vector we will perturb
    design_vector = des%get_values()
    constraint = object%get_value()
    call design_perturbed%init(design_vector%size())

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(design_vector%x, design_vector%x_d, &
            design_vector%size(), DEVICE_TO_HOST, .true.)
       call device_memcpy(target_sensitivities%x, &
            target_sensitivities%x_d, target_sensitivities%size(), &
            DEVICE_TO_HOST, .true.)
    end if

    write(*, '(I0,1X,A,F10.6,1X,A,F10.6,F10.6,F10.6,A)') &
         i, 'Design variable ', design_vector%x(i), &
         'Location [', des%get_x(i), des%get_y(i), des%get_z(i), ']'
    write(*, fmt_head) "Perturbation", "Constraint", "FD Estimate", "Error"
    write(*, fmt_data) 0.0_rp, constraint, target_sensitivities%x(i), 0.0_rp

    n_perturbations = size(perturbations)
    do ip = 1, n_perturbations
       perturb = perturbations(ip)
       tol = perturb / maxval(perturbations) * tolerance

       ! Ensure the perturbation stays within the bounds of the design variable
       if (design_vector%x(i) .gt. 0.5_rp) perturb = -perturb

       ! Reset and Perturb the design field by a small amount
       design_perturbed%x = design_vector%x
       design_perturbed%x(i) = design_vector%x(i) + perturb
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_memcpy(design_perturbed%x, design_perturbed%x_d, &
               design_perturbed%size(), HOST_TO_DEVICE, .true.)
       end if
       call des%update_design(design_perturbed)

       ! Compute the objective value of the perturbed design
       call sim%run_forward()
       call object%update_value(des)
       perturbed_constraint = object%get_value()
       call sim%reset()

       fd_estimate = perturbed_constraint - constraint
       if (.not. abscmp(fd_estimate, 0.0_rp)) then
          fd_estimate = fd_estimate / perturb
       end if

       fd_error = (fd_estimate - target_sensitivities%x(i)) / &
            target_sensitivities%x(i)

       write(*, fmt_data) perturb, perturbed_constraint, fd_estimate, fd_error

       if (abs(fd_error) .gt. tol) then
          !call neko_error('Finite difference estimate does not match ' // &
          !     'sensitivity')
       end if
    end do

  end subroutine compute_sensitivity_i

  subroutine compute_sensitivity_list(object, sim, des, target_sensitivities, &
       list, perturbations, tolerance)
    class(base_functional_t), intent(inout) :: object
    type(simulation_t), intent(inout) :: sim
    class(design_t), intent(inout) :: des
    type(vector_t), intent(inout) :: target_sensitivities
    integer, dimension(:), intent(in) :: list
    real(kind=rp), dimension(:), intent(in) :: perturbations
    real(kind=rp), intent(in) :: tolerance

    integer :: i, n

    n = size(list)
    do i = 1, n
       call compute_sensitivity_i(object, sim, des, target_sensitivities, &
            list(i), perturbations, tolerance)
    end do
  end subroutine compute_sensitivity_list

  subroutine compute_problem_sensitivity_i(problem, sim, des, target_sensitivities, i, &
       perturbations, tolerance, file_name)
    class(problem_t), intent(inout) :: problem
    type(simulation_t), intent(inout) :: sim
    class(design_t), intent(inout) :: des
    type(vector_t), intent(inout) :: target_sensitivities
    integer, intent(in) :: i
    real(kind=rp), intent(in) :: perturbations(:)
    real(kind=rp), intent(in) :: tolerance
    character(len=*), intent(in) :: file_name

    character(len=*), parameter :: fmt_head = '(4X,A12,4X,A10,6X,A11,5X,A5,10X)'
    character(len=*), parameter :: fmt_data = '(4X,4E15.6E3)'

    integer :: n_perturbations, ip
    real(kind=rp) :: perturb, tol
    type(vector_t) :: design_vector, design_perturbed, log_data
    real(kind=rp) :: constraint, perturbed_constraint
    real(kind=rp) :: fd_estimate, fd_error
    type(csv_file_t) :: logger

    ! Get the design vector for reference
    ! This is the design vector we will perturb
    design_vector = des%get_values()
    call problem%get_objective_value(constraint)
    call design_perturbed%init(design_vector%size())

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(design_vector%x, design_vector%x_d, &
            design_vector%size(), DEVICE_TO_HOST, .true.)
       call device_memcpy(target_sensitivities%x, &
            target_sensitivities%x_d, target_sensitivities%size(), &
            DEVICE_TO_HOST, .true.)
    end if

    write(*, '(I0,1X,A,F10.6,1X,A,F10.6,F10.6,F10.6,A)') &
         i, 'Design variable ', design_vector%x(i), &
         'Location [', des%get_x(i), des%get_y(i), des%get_z(i), ']'
    write(*, fmt_head) "Perturbation", "Constraint", "FD Estimate", "Error"
    write(*, fmt_data) 0.0_rp, constraint, target_sensitivities%x(i), 0.0_rp

    ! Init the csv writer
    call logger%init('FD_check_'//trim(file_name)//'.csv')
    call logger%set_header('perturbation,F,dFdx,error')
    call log_data%init(4)

    n_perturbations = size(perturbations)
    do ip = 1, n_perturbations
       perturb = perturbations(ip)
       tol = perturb / maxval(perturbations) * tolerance

       ! Ensure the perturbation stays within the bounds of the design variable
       if (design_vector%x(i) .gt. 0.5_rp) perturb = -perturb

       ! Reset and Perturb the design field by a small amount
       design_perturbed%x = design_vector%x
       design_perturbed%x(i) = design_vector%x(i) + perturb
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_memcpy(design_perturbed%x, design_perturbed%x_d, &
               design_perturbed%size(), HOST_TO_DEVICE, .true.)
       end if
       call des%update_design(design_perturbed)

       ! Compute the objective value of the perturbed design
       call problem%compute(des, sim)
       call problem%get_objective_value(perturbed_constraint)
       call sim%write(ip)
       call sim%reset()

       fd_estimate = perturbed_constraint - constraint
       if (.not. abscmp(fd_estimate, 0.0_rp)) then
          fd_estimate = fd_estimate / perturb
       end if

       fd_error = (fd_estimate - target_sensitivities%x(i)) / &
            target_sensitivities%x(i)

       write(*, fmt_data) perturb, perturbed_constraint, fd_estimate, fd_error

       log_data%x(1) = perturb
       log_data%x(2) = perturbed_constraint
       log_data%x(3) = fd_estimate
       log_data%x(4) = fd_error
       call logger%write(log_data)

       if (abs(fd_error) .gt. tol) then
          !call neko_error('Finite difference estimate does not match ' // &
          !     'sensitivity')
       end if
    end do

  end subroutine compute_problem_sensitivity_i

  subroutine compute_problem_sensitivity_list(problem, sim, des, target_sensitivities, &
       list, perturbations, tolerance, file_name)
    class(problem_t), intent(inout) :: problem
    type(simulation_t), intent(inout) :: sim
    class(design_t), intent(inout) :: des
    type(vector_t), intent(inout) :: target_sensitivities
    integer, dimension(:), intent(in) :: list
    real(kind=rp), dimension(:), intent(in) :: perturbations
    real(kind=rp), intent(in) :: tolerance
    character(len=*), intent(in) :: file_name

    integer :: i, n

    n = size(list)
    do i = 1, n
       call compute_problem_sensitivity_i(problem, sim, des, target_sensitivities, &
            list(i), perturbations, tolerance, file_name)
    end do
  end subroutine compute_problem_sensitivity_list

end module sensitivity
