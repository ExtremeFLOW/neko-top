!> Module designed to setup the topology optimization user interface for Neko.
!! This module should initialize, finalize, and run the time dependent
!! operations of our topology optimization problem.
!!
!! This should include, but not be limited to:
!! 1. Setting the material properties (Required by Neko when using the user
!!    interface)
!! 2. Setting the initial conditions for the scalar field
!! 3. Assigning the user defined force for the fluid
!! 4. Define the adjoint and sensitivity analysis
!! 5. Define the objective function
!! 6. Define the constraints
!!
!! Currently the computations should be setup in the "user_check" subroutine.
!! This is a hacky way to run the physics and adjoint physics, but it works.
module topology_optimization_user_module
  use case, only: case_t
  use field, only: field_t
  use json_file_module, only: json_file
  use num_types, only: rp
  use design_module, only: topopt_permeability_force
  use coefs, only: coef_t
  use initial_conditions, only: scalar_z_split_ic

  implicit none

  private
  public :: setup_topology_optimization_user_case

contains

  !> Assign user conditions for the neko case
  !!
  !!
  !! \param[inout] neko_case The neko case to setup the user interface for
  !!
  !! @todo We use a hacky way to run the physics and adjoint physics. This
  !! should be replaced with a more robust way to run the physics and adjoint
  !! physics.
  subroutine setup_topology_optimization_user_case(neko_case)

    implicit none

    type(case_t), intent(inout) :: neko_case

    ! General initialization and finalization
    neko_case%usr%user_init_modules => init_topopt_computation
    neko_case%usr%user_finalize_modules => finalize_topopt_computation

    ! Set the properties for the fluid
    neko_case%usr%material_properties => set_material_properties
    neko_case%usr%fluid_user_f_vector => topopt_permeability_force
    neko_case%usr%scalar_user_ic => scalar_z_split_ic

    ! Hack: Kinda hacky way to run the physics and adjoint physics
    neko_case%usr%user_check => topopt_computation

  end subroutine setup_topology_optimization_user_case

  !> Initialize the topology optimization computation
  subroutine init_topopt_computation(t, u, v, w, p, coef, params)
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params

  end subroutine init_topopt_computation

  !> Finalize the topology optimization computation
  subroutine finalize_topopt_computation(t, params)
    real(kind=rp) :: t
    type(json_file), intent(inout) :: params

  end subroutine finalize_topopt_computation

  !> Run the topology optimization computation
  subroutine topopt_computation(t, tstep, u, v, w, p, coef, params)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params

  end subroutine topopt_computation

  !> Initialize the material properties, unfortunately required from Neko.
  subroutine set_material_properties(t, tstep, rho, mu, cp, lambda, params)
    use json_file_module, only: json_file
    use json_utils, only: json_get
    use num_types, only: rp
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    real(kind=rp), intent(inout) :: rho, mu, cp, lambda
    type(json_file), intent(inout) :: params

    real(kind=rp) :: Re, Pe

    call json_get(params, 'case.fluid.Re', Re)
    call json_get(params, 'case.scalar.Pe', Pe)

    mu = 1.0_rp / Re
    lambda = 1.0_rp / Pe
    rho = 1.0_rp
    cp = 1.0_rp
  end subroutine set_material_properties




end module topology_optimization_user_module
