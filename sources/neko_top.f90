module neko_top

  use neko, only: case_t
  use design_module, only: design_t

  implicit none

  private
  public :: neko_top_init, neko_top_solve, neko_top_finalize

  type(design_t) :: design

contains

  !> Register user defined functions (see nekos user_intf.f90)
  subroutine neko_top_init(neko_case)
    use neko, only: neko_init
    use logger, only: neko_log
    use sensitivity, only: sensitivity_t
    use json_module, only: json_file
    use topology_optimization_user_module, only: setup_topology_optimization_user_case

    implicit none

    type(case_t), intent(inout) :: neko_case

    call setup_topology_optimization_user_case(neko_case)

    ! Initialize the neko solver
    call neko_init(neko_case)
    call neko_log%init()

    ! Initialize the design
    call design%init(neko_case)

  end subroutine neko_top_init

  subroutine neko_top_solve(neko_case)
    use neko, only: neko_solve
    use neko_ext, only: setup_iteration
    use json_utils, only: json_get_or_default
    implicit none

    type(case_t), intent(inout) :: neko_case

    integer :: iter, max_iter

    ! Convergence parameter
    logical :: converged = .false.

    call json_get_or_default(neko_case%params, &
                             'case.topology_optimization.max_iter', &
                             max_iter, 4)

    do iter = 1, max_iter
       call setup_iteration(neko_case, iter)

       call neko_solve(neko_case)

       ! Update the permeability
       call design%update(converged)

       if (converged) exit
    end do
  end subroutine neko_top_solve

  subroutine neko_top_finalize(neko_case)
    use neko, only: neko_finalize
    use develop, only: estimate_temperature
    use logger, only: neko_log
    use json_utils, only: json_get_or_default

    type(case_t), intent(inout) :: neko_case

    logical :: temperature_enabled

    ! ---------------------------------------------------------------------- !
    ! Compute the outlet area-weighted average temperature
    ! ---------------------------------------------------------------------- !
    ! Read the case file for options
    call json_get_or_default(neko_case%params, 'case.scalar.enabled', &
                             temperature_enabled, .false.)

    if (temperature_enabled) then
       call estimate_temperature(neko_case)
    end if

    ! ---------------------------------------------------------------------- !
    ! Save the design
    ! ---------------------------------------------------------------------- !

    call neko_finalize(neko_case)
    call design%free()
    call neko_log%end()

  end subroutine neko_top_finalize

end module neko_top

