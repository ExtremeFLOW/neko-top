module sensitivity
  use simulation_m, only: simulation_t
  use design, only: design_t
  use base_functional, only: base_functional_t
  use utils, only: neko_error
  use num_types, only: rp
  use math, only: abscmp
  use vector, only: vector_t
  implicit none

  interface compute_sensitivity
     module procedure compute_sensitivity_i, compute_sensitivity_list
  end interface compute_sensitivity

contains

  subroutine compute_sensitivity_i(object, sim, des, target_sensitivities, i, &
       perturbations, tolerance)
    class(base_functional_t), intent(inout) :: object
    type(simulation_t), intent(inout) :: sim
    class(design_t), intent(inout) :: des
    type(vector_t), intent(in) :: target_sensitivities
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
       design_perturbed = design_vector
       design_perturbed%x(i) = design_vector%x(i) + perturb
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
          call neko_error('Finite difference estimate does not match ' // &
               'sensitivity')
       end if
    end do

  end subroutine compute_sensitivity_i

  subroutine compute_sensitivity_list(object, sim, des, target_sensitivities, &
       list, perturbations, tolerance)
    class(base_functional_t), intent(inout) :: object
    type(simulation_t), intent(inout) :: sim
    class(design_t), intent(inout) :: des
    type(vector_t), intent(in) :: target_sensitivities
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

end module sensitivity
