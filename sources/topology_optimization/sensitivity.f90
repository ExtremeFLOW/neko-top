! Todo: This module should follow the templates of simulation components. This
! should make it easier to add new components to the code, whenever user
! components gets supported.

!> Sensitivity module.
!! This module contains the sensitivity computation of the topology optimization.
module sensitivity
  use simulation_component, only : simulation_component_t
  use json_file_module, only : json_file
  use case, only : case_t
  use num_types, only : rp
  use field, only : field_t

  implicit none
  private

  type, public :: sensitivity_t

     !> test field
     type(field_t), pointer :: test

   contains

     !> Constructor.
     procedure, pass(this) :: init => sensitivity_init
     !> Destructor.
     procedure, pass(this) :: free => sensitivity_free
     !> Compute the lambda2 field
     procedure, pass(this) :: update => update_sensitivity
  end type sensitivity_t

contains

  !> Actual constructor.
  subroutine sensitivity_init(this, json, neko_case)
    use field_registry, only : neko_field_registry
    use json_utils, only : json_get_or_default
    implicit none

    class(sensitivity_t), intent(inout) :: this
    class(json_file), intent(inout) :: json
    class(case_t), intent(inout) :: neko_case

    ! Allocate new fields if they don't exist yet.
    if (.not. neko_field_registry%field_exists("test")) then
       call neko_field_registry%add_field(neko_case%fluid%u%dof, "test")
    end if

    this%test => neko_field_registry%get_field_by_name("test")

    call neko_case%f_out%fluid%append(this%test)

  end subroutine sensitivity_init

  !> Destructor.
  subroutine sensitivity_free(this)
    class(sensitivity_t), intent(inout) :: this

  end subroutine sensitivity_free

  !> Compute the sensitivities field.
  !! @param t The time value.
  !! @param tstep The current time-step
  subroutine update_sensitivity(this, t, tstep)
    class(sensitivity_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    call compute_sensitivity(this)

  end subroutine update_sensitivity

  ! ========================================================================== !
  ! Each part of the simulation component is implemented as a subroutine.

  !> Compute the sensitivity of our topology optimization.
  subroutine compute_sensitivity(this)
    class(sensitivity_t), intent(inout) :: this

  end subroutine compute_sensitivity

  subroutine compute_adjoint(this)
    class(sensitivity_t), intent(inout) :: this

  end subroutine compute_adjoint

end module sensitivity
