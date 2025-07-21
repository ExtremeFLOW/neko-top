module sensitivity
  use simulation_m, only: simulation_t
  use design, only: design_t
  use base_functional, only: base_functional_t
  use utils, only: neko_error
  use num_types, only: rp
  use math, only: abscmp
  use vector, only: vector_t
  implicit none

contains

  subroutine compute_sensitivity(object, sim, des, i, perturbations, tolerance)
    class(base_functional_t), intent(inout) :: object
    type(simulation_t), intent(inout) :: sim
    class(design_t), intent(inout) :: des
    integer, intent(in) :: i
    real(kind=rp), intent(in) :: perturbations(:)
    real(kind=rp), intent(in) :: tolerance

    character(len=*), parameter :: fmt_head = '(4X,A12,4X,A10,6X,A11,5X,A5,10X)'
    character(len=*), parameter :: fmt_data = '(4X,4E15.6E3)'

    integer :: n_perturbations, ip
    real(kind=rp) :: perturb
    type(vector_t) :: design_vector, design_perturbed
    real(kind=rp) :: constraint, perturbed_constraint
    type(vector_t) :: constraint_sensitivities
    real(kind=rp) :: fd_estimate, fd_error

    ! Get the design vector for reference
    ! This is the design vector we will perturb
    design_vector = des%get_values()
    constraint = object%get_value()
    constraint_sensitivities = object%get_sensitivity()

    write(*, '(I0,1X,A,F10.6,1X,A,F10.6,F10.6,F10.6,A)') &
         i, 'Design variable ', design_vector%x(i), &
         'Location [', des%get_x(i), des%get_y(i), des%get_z(i), ']'
    write(*, fmt_head) "Perturbation", "Constraint", "FD Estimate", "Error"
    write(*, fmt_data) 0.0_rp, constraint, constraint_sensitivities%x(i), 0.0_rp

    n_perturbations = size(perturbations)
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
       call object%update_value(des)
       perturbed_constraint = object%get_value()
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

  end subroutine compute_sensitivity

end module sensitivity
