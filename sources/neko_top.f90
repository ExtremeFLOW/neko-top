module neko_top

  ! External modules
  use case, only: case_t
  use neko, only: neko_init, neko_finalize
  use simcomp_executor, only: simcomp_executor_t
  use json_module, only: json_file, json_value, json_core
  use comm, only: pe_rank
  use utils, only: neko_error, filename_suffix
  use user_intf, only: simulation_component_user_settings
  use num_types, only: rp
  use logger, only: neko_log

  use design_module, only: design_t
  use sensitivity, only: sensitivity_t
  use topology_optimization_user_module, only: neko_user_init
  use, intrinsic :: iso_fortran_env

  implicit none

  private
  public :: neko_top_init, neko_top_solve, neko_top_finalize

  type(simcomp_executor_t) :: topopt_components

contains

  !> Register user defined functions (see nekos user_intf.f90)
  subroutine neko_top_init(neko_case)
    type(case_t), intent(inout) :: neko_case
    type(json_file) :: params, design_params
    type(design_t), allocatable :: design
    type(json_value), pointer :: simcomp_object
    type(json_file) :: comp_subdict
    logical :: found
    type(json_core) :: core

    call neko_log%section('Initialize Topology Optimization')

    ! ------------------------------------------------------------------------ !
    ! Initialize the neko solver
    ! ------------------------------------------------------------------------ !

    call neko_log%section('Neko')

    call neko_user_init(neko_case)
    call neko_init(neko_case)

    call neko_log%end()

    ! ------------------------------------------------------------------------ !
    ! Initialize the topopt components
    ! ------------------------------------------------------------------------ !

    call neko_log%section('Topology Optimization Components')

    call topopt_components%init(neko_case, 'topology_optimization.components')

    call neko_case%params%get_core(core)
    call neko_case%params%get('topology_optimization.components', simcomp_object, found)
    comp_subdict = json_file(simcomp_object)

    ! Allocation of the design
    allocate(design)
    design_params = simulation_component_user_settings('design', comp_subdict)
    call topopt_components%add_user_simcomp(design, design_params)

    call topopt_components%finalize()

    call neko_log%end()

    call neko_log%end()
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

       call topopt_components%preprocess(0.0_rp, 1)

       ! Forward analysis
       call neko_solve(neko_case)



       ! If converged, exit the loop
       if (converged) exit

       ! Call the design update routine


       call topopt_components%compute(0.0_rp, 1)

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
    call neko_log%end()

  end subroutine neko_top_finalize

end module neko_top

