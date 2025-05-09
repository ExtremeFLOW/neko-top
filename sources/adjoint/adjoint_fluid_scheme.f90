! Copyright (c) 2020-2024, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> Fluid formulations
module adjoint_fluid_scheme
  use gather_scatter, only: gs_t
  use mean_sqr_flow, only: mean_sqr_flow_t
  use neko_config, only: NEKO_BCKND_DEVICE
  use checkpoint, only: chkp_t
  use mean_flow, only: mean_flow_t
  use num_types, only: rp, i8
  use comm, only: NEKO_COMM
  use adjoint_source_term, only: adjoint_source_term_t
  use field, only: field_t
  use space, only: space_t, GLL
  use dofmap, only: dofmap_t
  use zero_dirichlet, only: zero_dirichlet_t
  use krylov, only: ksp_t, krylov_solver_factory, KSP_MAX_ITER
  use coefs, only: coef_t
  use usr_inflow, only: usr_inflow_t, usr_inflow_eval
  use dirichlet, only: dirichlet_t
  use field_dirichlet, only: field_dirichlet_t
  use field_dirichlet_vector, only: field_dirichlet_vector_t
  use jacobi, only: jacobi_t
  use sx_jacobi, only: sx_jacobi_t
  use device_jacobi, only: device_jacobi_t
  use hsmg, only: hsmg_t
  use phmg, only: phmg_t
  use precon, only: pc_t, precon_factory, precon_destroy
  use fluid_stats, only: fluid_stats_t
  use bc, only: bc_t
  use bc_list, only: bc_list_t
  use mesh, only: mesh_t, NEKO_MSH_MAX_ZLBLS, NEKO_MSH_MAX_ZLBL_LEN
  use math, only: cfill, add2s2, glsum
  use device_math, only: device_cfill, device_add2s2
  use time_scheme_controller, only: time_scheme_controller_t
  use operators, only: cfl
  use logger, only: neko_log, LOG_SIZE, NEKO_LOG_VERBOSE
  use field_registry, only: neko_field_registry
  use json_utils, only: json_get, json_get_or_default, json_extract_object, &
       json_extract_item
  use json_module, only: json_file, json_core, json_value
  use scratch_registry, only: scratch_registry_t
  use user_intf, only: user_t, dummy_user_material_properties, &
       user_material_properties
  use utils, only: neko_error, neko_warning
  use field_series, only: field_series_t
  use time_step_controller, only: time_step_controller_t
  use field_math, only: field_cfill, field_add2s2
  use wall_model_bc, only: wall_model_bc_t
  use shear_stress, only: shear_stress_t
  use field_list, only : field_list_t
  use gradient_jump_penalty, only: gradient_jump_penalty_t
  use field_math, only: field_addcol3

  use mpi_f08, only: MPI_INTEGER, MPI_SUM, MPI_Allreduce
  use json_utils_ext, only: json_key_fallback
  use device, only : device_event_sync, glb_cmd_event, DEVICE_TO_HOST, &
       device_memcpy
  implicit none
  private

  !> Base type of all fluid formulations
  type, abstract :: adjoint_fluid_scheme_t
     !> A name that can be used to distinguish this solver in e.g. user routines
     character(len=:), allocatable :: name
     !> x-component of Velocity
     type(field_t), pointer :: u_adj => null()
     !> y-component of Velocity
     type(field_t), pointer :: v_adj => null()
     !> z-component of Velocity
     type(field_t), pointer :: w_adj => null()
     !> Pressure
     type(field_t), pointer :: p_adj => null()
     !> fluid field (lag)
     type(field_series_t) :: ulag, vlag, wlag
     !> Function space \f$ X_h \f$
     type(space_t) :: Xh
     !> Dofmap associated with \f$ X_h \f$
     type(dofmap_t) :: dm_Xh
     !> Gather-scatter associated with \f$ X_h \f$
     type(gs_t) :: gs_Xh
     !> Coefficients associated with \f$ X_h \f$
     type(coef_t) :: c_Xh
     !> The source term for the momentum equation.
     type(adjoint_source_term_t) :: source_term
     !> X-component of the right-hand side.
     type(field_t), pointer :: f_adj_x => null()
     !> Y-component of the right-hand side.
     type(field_t), pointer :: f_adj_y => null()
     !> Z-component of the right-hand side.
     type(field_t), pointer :: f_adj_z => null()

     ! Krylov solvers and settings
     !> Krylov solver for velocity
     class(ksp_t), allocatable :: ksp_vel
     !> Krylov solver for pressure
     class(ksp_t), allocatable :: ksp_prs
     !> Velocity Preconditioner
     class(pc_t), allocatable :: pc_vel
     !> Velocity Preconditioner
     class(pc_t), allocatable :: pc_prs
     !> Size of the projection space for ksp_vel
     integer :: vel_projection_dim
     !> Size of the projection space for ksp_pr
     integer :: pr_projection_dim
     !> Steps to activate projection for ksp_vel
     integer :: vel_projection_activ_step
     !> Steps to activate projection for ksp_pr
     integer :: pr_projection_activ_step
     !> Strict convergence for the velocity solver
     logical :: strict_convergence

     ! List of boundary conditions for pressure
     type(bc_list_t) :: bcs_prs
     ! List of boundary conditions for velocity
     type(bc_list_t) :: bcs_vel
     !> Mesh
     type(mesh_t), pointer :: msh => null()
     !> Checkpoint
     type(chkp_t) :: chkp
     !> Mean flow field
     type(mean_flow_t) :: mean
     !> Fluid statistics
     type(fluid_stats_t) :: stats
     !> Mean squared flow field
     type(mean_sqr_flow_t) :: mean_sqr
     !> Is the flow rate forced?
     logical :: forced_flow_rate = .false.
     !> Freeze velocity at initial condition?
     logical :: freeze = .false.

     !> The turbulent kinematic viscosity field name
     character(len=:), allocatable :: nut_field_name
     !> Is mu varying in time? Currently only due to LES models.
     logical :: variable_material_properties = .false.

     !> Global number of GLL points for the fluid (not unique)
     integer(kind=i8) :: glb_n_points
     !> Global number of GLL points for the fluid (unique)
     integer(kind=i8) :: glb_unique_points
     !> Manager for temporary fields
     type(scratch_registry_t) :: scratch

     !> Density field
     type(field_t) :: rho

     !> The dynamic viscosity
     type(field_t) :: mu

     !> A helper that packs material properties to pass to the user routine.
     type(field_list_t) :: material_properties

     !> User material properties routine
     procedure(user_material_properties), nopass, pointer :: &
          user_material_properties => null()
   contains
     !> Constructor for the base type
     procedure, pass(this) :: init_base => adjoint_fluid_scheme_init_base
     !> Destructor for the base type
     procedure, pass(this) :: scheme_free => adjoint_fluid_scheme_free
     !> Validate that all components are properly allocated
     procedure, pass(this) :: validate => adjoint_fluid_scheme_validate
     !> Apply pressure boundary conditions
     procedure, pass(this) :: bc_apply_vel => adjoint_fluid_scheme_bc_apply_vel
     !> Apply velocity boundary conditions
     procedure, pass(this) :: bc_apply_prs => adjoint_fluid_scheme_bc_apply_prs
     !> Compute the CFL number
     procedure, pass(this) :: compute_cfl => adjoint_compute_cfl
     !> Set rho and mu
     procedure, pass(this) :: set_material_properties => &
          adjoint_fluid_scheme_set_material_properties
     !> Constructor
     procedure(adjoint_fluid_scheme_init_intrf), pass(this), deferred :: init
     !> Destructor
     procedure(adjoint_fluid_scheme_free_intrf), pass(this), deferred :: free
     !> Advance one step in time
     procedure(adjoint_fluid_scheme_step_intrf), pass(this), deferred :: step
     !> Restart from a checkpoint
     procedure(adjoint_fluid_scheme_restart_intrf), &
          pass(this), deferred :: restart
     !> Setup boundary conditions
     procedure(adjoint_fluid_scheme_setup_bcs_intrf), pass(this), deferred :: &
          setup_bcs
     !> Update variable material properties
     procedure, pass(this) :: update_material_properties => &
          adjoint_fluid_scheme_update_material_properties
     !> Linear solver factory, wraps a KSP constructor
     procedure, nopass :: solver_factory => adjoint_fluid_scheme_solver_factory
     !> Preconditioner factory
     procedure, pass(this) :: precon_factory_ => &
          adjoint_fluid_scheme_precon_factory
  end type adjoint_fluid_scheme_t


  !> Abstract interface to initialize an adjoint formulation
  abstract interface
     subroutine adjoint_fluid_scheme_init_intrf(this, msh, lx, params, user, &
          time_scheme)
       import adjoint_fluid_scheme_t
       import json_file
       import mesh_t
       import user_t
       import time_scheme_controller_t
       class(adjoint_fluid_scheme_t), target, intent(inout) :: this
       type(mesh_t), target, intent(inout) :: msh
       integer, intent(in) :: lx
       type(json_file), target, intent(inout) :: params
       type(user_t), target, intent(in) :: user
       type(time_scheme_controller_t), target, intent(in) :: time_scheme
     end subroutine adjoint_fluid_scheme_init_intrf
  end interface

  !> Abstract interface to dealocate an adjoint formulation
  abstract interface
     subroutine adjoint_fluid_scheme_free_intrf(this)
       import adjoint_fluid_scheme_t
       class(adjoint_fluid_scheme_t), intent(inout) :: this
     end subroutine adjoint_fluid_scheme_free_intrf
  end interface

  !> Abstract interface to compute a time-step
  abstract interface
     subroutine adjoint_fluid_scheme_step_intrf(this, t, tstep, dt, ext_bdf, &
          dt_controller)
       import adjoint_fluid_scheme_t
       import time_scheme_controller_t
       import time_step_controller_t
       import rp
       class(adjoint_fluid_scheme_t), target, intent(inout) :: this
       real(kind=rp), intent(in) :: t
       integer, intent(in) :: tstep
       real(kind=rp), intent(in) :: dt
       type(time_scheme_controller_t), intent(in) :: ext_bdf
       type(time_step_controller_t), intent(in) :: dt_controller
     end subroutine adjoint_fluid_scheme_step_intrf
  end interface

  !> Abstract interface to restart an adjoint scheme
  abstract interface
     subroutine adjoint_fluid_scheme_restart_intrf(this, dtlag, tlag)
       import adjoint_fluid_scheme_t
       import rp
       class(adjoint_fluid_scheme_t), target, intent(inout) :: this
       real(kind=rp) :: dtlag(10), tlag(10)

     end subroutine adjoint_fluid_scheme_restart_intrf
  end interface

  !> Abstract interface to setup boundary conditions
  abstract interface
     subroutine adjoint_fluid_scheme_setup_bcs_intrf(this, user, params)
       import adjoint_fluid_scheme_t, user_t, json_file
       class(adjoint_fluid_scheme_t), intent(inout) :: this
       type(user_t), target, intent(in) :: user
       type(json_file), intent(inout) :: params
     end subroutine adjoint_fluid_scheme_setup_bcs_intrf
  end interface

  interface
     !> Initialise an adjoint scheme
     module subroutine adjoint_fluid_scheme_factory(object, type_name)
       class(adjoint_fluid_scheme_t), intent(inout), allocatable :: object
       character(len=*) :: type_name
     end subroutine adjoint_fluid_scheme_factory
  end interface

  public :: adjoint_fluid_scheme_t, adjoint_fluid_scheme_factory

contains

  !> Initialize common data for the current scheme
  subroutine adjoint_fluid_scheme_init_base(this, msh, lx, params, scheme, &
       user, kspv_init)
    class(adjoint_fluid_scheme_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(in) :: lx
    character(len=*), intent(in) :: scheme
    type(json_file), target, intent(inout) :: params
    type(user_t), target, intent(in) :: user
    logical, intent(in) :: kspv_init
    character(len=LOG_SIZE) :: log_buf
    real(kind=rp) :: real_val
    logical :: logical_val
    integer :: integer_val
    character(len=:), allocatable :: string_val1, string_val2
    character(len=:), allocatable :: json_key
    type(json_file) :: json_subdict

    !
    ! SEM simulation fundamentals
    !

    this%msh => msh

    if (msh%gdim .eq. 2) then
       call this%Xh%init(GLL, lx, lx)
    else
       call this%Xh%init(GLL, lx, lx, lx)
    end if

    call this%dm_Xh%init(msh, this%Xh)

    call this%gs_Xh%init(this%dm_Xh)

    call this%c_Xh%init(this%gs_Xh)

    ! Local scratch registry
    this%scratch = scratch_registry_t(this%dm_Xh, 10, 2)

    !
    ! First section of fluid log
    !

    call neko_log%section('Adjoint fluid')
    write(log_buf, '(A, A)') 'Type       : ', trim(scheme)
    call neko_log%message(log_buf)

    !
    ! Material properties
    !
    call this%set_material_properties(params, user)

    ! Projection spaces
    json_key = json_key_fallback(params, &
         'case.adjoint_fluid.velocity_solver.projection_space_size', &
         'case.fluid.velocity_solver.projection_space_size')
    call json_get_or_default(params, json_key, this%vel_projection_dim, 20)
    json_key = json_key_fallback(params, &
         'case.adjoint_fluid.pressure_solver.projection_space_size', &
         'case.fluid.pressure_solver.projection_space_size')
    call json_get_or_default(params, json_key, this%pr_projection_dim, 20)

    json_key = json_key_fallback(params, &
         'case.adjoint_fluid.velocity_solver.projection_hold_steps', &
         'case.fluid.velocity_solver.projection_hold_steps')
    call json_get_or_default(params, json_key, &
         this%vel_projection_activ_step, 5)
    json_key = json_key_fallback(params, &
         'case.adjoint_fluid.pressure_solver.projection_hold_steps', &
         'case.fluid.pressure_solver.projection_hold_steps')
    call json_get_or_default(params, json_key, &
         this%pr_projection_activ_step, 5)

    json_key = json_key_fallback(params, 'case.adjoint_fluid.freeze', &
         'case.fluid.freeze')
    call json_get_or_default(params, json_key, this%freeze, .false.)

    ! TODO
    ! This needs to be discussed... In pricipal, I think if the forward is
    ! forced to a
    ! certain flow rate, then the adjoint should be forced to zero flow rate,
    ! we had a derivation in
    ! https://www.sciencedirect.com/science/article/pii/S0045782522006764
    ! for now I'm commenting this out
    if (params%valid_path("case.fluid.flow_rate_force")) then
       call neko_error('Flow rate forcing not yet implemented')
       this%forced_flow_rate = .true.
    end if


    if (lx .lt. 10) then
       write(log_buf, '(A, I1)') 'Poly order : ', lx-1
    else if (lx .ge. 10) then
       write(log_buf, '(A, I2)') 'Poly order : ', lx-1
    else
       write(log_buf, '(A, I3)') 'Poly order : ', lx-1
    end if
    call neko_log%message(log_buf)
    this%glb_n_points = int(this%msh%glb_nelv, i8)*int(this%Xh%lxyz, i8)
    this%glb_unique_points = int(glsum(this%c_Xh%mult, this%dm_Xh%size()), i8)

    write(log_buf, '(A, I0)') 'GLL points : ', this%glb_n_points
    call neko_log%message(log_buf)
    write(log_buf, '(A, I0)') 'Unique pts.: ', this%glb_unique_points
    call neko_log%message(log_buf)

    call json_get(params, 'case.numerics.dealias', logical_val)
    write(log_buf, '(A, L1)') 'Dealias    : ', logical_val
    call neko_log%message(log_buf)

    write(log_buf, '(A, L1)') 'LES        : ', this%variable_material_properties
    call neko_log%message(log_buf)

    call json_get_or_default(params, 'case.output_boundary', logical_val, &
         .false.)
    write(log_buf, '(A, L1)') 'Save bdry  : ', logical_val
    call neko_log%message(log_buf)

    !
    ! Setup right-hand side fields.
    !
    allocate(this%f_adj_x)
    allocate(this%f_adj_y)
    allocate(this%f_adj_z)
    call this%f_adj_x%init(this%dm_Xh, fld_name = "adjoint_rhs_x")
    call this%f_adj_y%init(this%dm_Xh, fld_name = "adjoint_rhs_y")
    call this%f_adj_z%init(this%dm_Xh, fld_name = "adjoint_rhs_z")

    ! Initialize the source term
    ! TODO
    ! HARRY
    ! Ok, this one I dislike the most... I think the CORRECT solution here is to
    ! modify
    ! the fluid_source_term.f90 so that case.fluid is not hardcoded
    ! whoever coded it originally is probably much smarter than me and has their
    ! reasons
    !
    ! I would envision somehow merging fluid_source_term and scalar_source_term
    ! (maybe they're different because because it's a vector forcing vs a scalar
    ! forcing,
    ! but still...
    !
    ! this is my attempt...
    ! params replaced by adjoint_json
    !call this%source_term%init(params, this%f_adj_x, this%f_adj_y, &
    ! this%f_adj_z, this%c_Xh,&
    !     user)

    call this%source_term%init(this%f_adj_x, this%f_adj_y, this%f_adj_z, &
         this%c_Xh, user)
    call this%source_term%add(params, 'case.adjoint_fluid.source_term')

    ! Todo: HARRY
    ! --------------------------------------------------------------------------
    ! ok here is a chance that we can maybe steal the krylov solvers from the
    ! forward??
    !
    ! something along the lines of
    !
    ! if (.not.params%valid_path('case.adjoint_fluid.velocity_solver')) then
    ! 	this%ksp_vel => case%fluid%ksp_vel
    ! 	this%pc_vel => case%fluid%ksp_vel
    ! else
    !  initialize everything
    ! end if
    !
    ! but I don't know how we can steal the krylov solvers out of the forward,
    ! so for now we'll just initialize two of them...
    !  Initialize velocity solver
    if (kspv_init) then
       call neko_log%section("Adjoint Velocity solver")

       json_key = json_key_fallback(params, &
            'case.adjoint_fluid.velocity_solver.max_iterations', &
            'case.fluid.velocity_solver.max_iterations')
       call json_get_or_default(params, json_key, integer_val, KSP_MAX_ITER)

       json_key = json_key_fallback(params, &
            'case.adjoint_fluid.velocity_solver.type', &
            'case.fluid.velocity_solver.type')
       call json_get(params, json_key, string_val1)

       json_key = json_key_fallback(params, &
            'case.adjoint_fluid.velocity_solver.preconditioner.type', &
            'case.fluid.velocity_solver.preconditioner.type')
       call json_get(params, json_key, string_val2)

       json_key = json_key_fallback(params, &
            'case.adjoint_fluid.velocity_solver.preconditioner', &
            'case.fluid.velocity_solver.preconditioner')
       call json_extract_object(params, json_key, json_subdict)

       json_key = json_key_fallback(params, &
            'case.adjoint_fluid.velocity_solver.absolute_tolerance', &
            'case.fluid.velocity_solver.absolute_tolerance')
       call json_get(params, json_key, real_val)

       json_key = json_key_fallback(params, &
            'case.adjoint_fluid.velocity_solver.monitor', &
            'case.fluid.velocity_solver.monitor')
       call json_get_or_default(params, json_key, logical_val, .false.)

       call neko_log%message('Type       : ('// trim(string_val1) // &
            ', ' // trim(string_val2) // ')')

       write(log_buf, '(A,ES13.6)') 'Abs tol    :', real_val
       call neko_log%message(log_buf)
       call this%solver_factory(this%ksp_vel, this%dm_Xh%size(), &
            string_val1, integer_val, real_val, logical_val)
       call this%precon_factory_(this%pc_vel, this%ksp_vel, &
            this%c_Xh, this%dm_Xh, this%gs_Xh, this%bcs_vel, &
            string_val2, json_subdict)
       call neko_log%end_section()
    end if

    ! Strict convergence for the velocity solver
    call json_get_or_default(params, 'case.fluid.strict_convergence', &
         this%strict_convergence, .false.)

    ! Assign velocity fields
    call neko_field_registry%add_field(this%dm_Xh, 'u_adj')
    call neko_field_registry%add_field(this%dm_Xh, 'v_adj')
    call neko_field_registry%add_field(this%dm_Xh, 'w_adj')
    this%u_adj => neko_field_registry%get_field('u_adj')
    this%v_adj => neko_field_registry%get_field('v_adj')
    this%w_adj => neko_field_registry%get_field('w_adj')

    !! Initialize time-lag fields
    call this%ulag%init(this%u_adj, 2)
    call this%vlag%init(this%v_adj, 2)
    call this%wlag%init(this%w_adj, 2)

    call neko_log%end_section()

  end subroutine adjoint_fluid_scheme_init_base

  ! ========================================================================== !
  ! Todo: This section need to be moved

  ! !> Initialize all components of the current scheme
  ! subroutine adjoint_fluid_scheme_init_all(this, msh, lx, params, kspv_init, &
  !      kspp_init, scheme, user)
  !   implicit none
  !   class(adjoint_fluid_scheme_t), target, intent(inout) :: this
  !   type(mesh_t), target, intent(inout) :: msh
  !   integer, intent(inout) :: lx
  !   type(json_file), target, intent(inout) :: params
  !   type(user_t), target, intent(in) :: user
  !   logical :: kspv_init
  !   logical :: kspp_init
  !   character(len=*), intent(in) :: scheme
  !   real(kind=rp) :: abs_tol
  !   integer :: integer_val, ierr
  !   character(len=:), allocatable :: solver_type, precon_type, json_key
  !   character(len=LOG_SIZE) :: log_buf

  !   call adjoint_fluid_scheme_init_base(this, msh, lx, params, scheme, user, &
  !        kspv_init)

  !   call neko_field_registry%add_field(this%dm_Xh, 'p_adj')
  !   this%p_adj => neko_field_registry%get_field('p_adj')

  !   ! HARRY
  !   ! Holy shit pressure BC's confuse me on a good day...
  !   ! but they're especially confusing for the adjoint.
  !   ! for now... I'm keeping them the same as the primal
  !   ! but I wouldn't be surprised if this is incorrect for the adjoint
  !   !
  !   ! Setup pressure boundary conditions
  !   !
  !   call this%bcs_prs%init()
  !   call this%bc_prs%init_base(this%c_Xh)
  !   call this%bc_prs%mark_zones_from_list(msh%labeled_zones,&
  !        'o', this%bc_labels)
  !   call this%bc_prs%mark_zones_from_list(msh%labeled_zones,&
  !        'on', this%bc_labels)

  !   ! Field dirichlet pressure bc
  !   call this%user_field_bc_prs%init_base(this%c_Xh)
  !   call this%user_field_bc_prs%mark_zones_from_list(msh%labeled_zones,&
  !        'd_pres', this%bc_labels)
  !   call this%user_field_bc_prs%finalize()
  !   call MPI_Allreduce(this%user_field_bc_prs%msk(0), integer_val, 1, &
  !        MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)

  !   if (integer_val .gt. 0) call this%user_field_bc_prs%init_field('d_pres')
  !   call this%bcs_prs%append(this%user_field_bc_prs)
  !   call this%user_field_bc_vel%bc_list%append(this%user_field_bc_prs)

  !   if (msh%outlet%size .gt. 0) then
  !      call this%bc_prs%mark_zone(msh%outlet)
  !   end if
  !   if (msh%outlet_normal%size .gt. 0) then
  !      call this%bc_prs%mark_zone(msh%outlet_normal)
  !   end if

  !   call this%bc_prs%finalize()
  !   call this%bc_prs%set_g(0.0_rp)
  !   call this%bcs_prs%append(this%bc_prs)
  !   call this%bc_dong%init_base(this%c_Xh)
  !   call this%bc_dong%mark_zones_from_list(msh%labeled_zones,&
  !        'o+dong', this%bc_labels)
  !   call this%bc_dong%mark_zones_from_list(msh%labeled_zones,&
  !        'on+dong', this%bc_labels)
  !   call this%bc_dong%finalize()


  !   call this%bc_dong%init(this%c_Xh, params)

  !   call this%bcs_prs%append(this%bc_dong)

  !   ! Pressure solver
  !   ! hopefully we can do the same trick here as we will eventually do for the
  !   ! velocity solver
  !   if (kspp_init) then
  !      call neko_log%section("Pressure solver")

  !      json_key = json_key_fallback(params, &
  !           'case.adjoint_fluid.pressure_solver.max_iterations', &
  !           'case.fluid.pressure_solver.max_iterations')
  !      call json_get_or_default(params, json_key, integer_val, 800)

  !      json_key = json_key_fallback(params, &
  !           'case.adjoint_fluid.pressure_solver.type', &
  !           'case.fluid.pressure_solver.type')
  !      call json_get(params, json_key, solver_type)

  !      json_key = json_key_fallback(params, &
  !           'case.adjoint_fluid.pressure_solver.preconditioner', &
  !           'case.fluid.pressure_solver.preconditioner')
  !      call json_get(params, json_key, precon_type)

  !      json_key = json_key_fallback(params, &
  !           'case.adjoint_fluid.pressure_solver.absolute_tolerance', &
  !           'case.fluid.pressure_solver.absolute_tolerance')
  !      call json_get(params, json_key, abs_tol)

  !      call neko_log%message('Type       : ('// trim(solver_type) // &
  !           ', ' // trim(precon_type) // ')')
  !      write(log_buf, '(A,ES13.6)') 'Abs tol    :', abs_tol
  !      call neko_log%message(log_buf)

  !      call adjoint_fluid_scheme_solver_factory(this%ksp_prs, &
  !           this%dm_Xh%size(), &
  !           solver_type, integer_val, abs_tol)
  !      call adjoint_fluid_scheme_precon_factory(this%pc_prs, this%ksp_prs, &
  !           this%c_Xh, this%dm_Xh, this%gs_Xh, this%bcs_prs, precon_type)

  !      call neko_log%end_section()

  !   end if


  !   call neko_log%end_section()

  ! end subroutine adjoint_fluid_scheme_init_all

  ! End of section to be moved
  ! ========================================================================== !

  !> Deallocate a fluid formulation
  subroutine adjoint_fluid_scheme_free(this)
    class(adjoint_fluid_scheme_t), intent(inout) :: this

    !
    ! Free everything related to field_dirichlet BCs
    !

    call this%Xh%free()

    if (allocated(this%ksp_vel)) then
       call this%ksp_vel%free()
       deallocate(this%ksp_vel)
    end if

    if (allocated(this%ksp_prs)) then
       call this%ksp_prs%free()
       deallocate(this%ksp_prs)
    end if

    if (allocated(this%pc_vel)) then
       call precon_destroy(this%pc_vel)
       deallocate(this%pc_vel)
    end if

    if (allocated(this%pc_prs)) then
       call precon_destroy(this%pc_prs)
       deallocate(this%pc_prs)
    end if

    call this%source_term%free()

    call this%gs_Xh%free()

    call this%c_Xh%free()

    call this%scratch%free()

    nullify(this%u_adj)
    nullify(this%v_adj)
    nullify(this%w_adj)
    nullify(this%p_adj)

    call this%ulag%free()
    call this%vlag%free()
    call this%wlag%free()


    if (associated(this%f_adj_x)) then
       call this%f_adj_x%free()
    end if

    if (associated(this%f_adj_y)) then
       call this%f_adj_y%free()
    end if

    if (associated(this%f_adj_z)) then
       call this%f_adj_z%free()
    end if

    nullify(this%f_adj_x)
    nullify(this%f_adj_y)
    nullify(this%f_adj_z)

    call this%mu%free()

  end subroutine adjoint_fluid_scheme_free

  !> Validate that all fields, solvers etc necessary for
  !! performing time-stepping are defined
  subroutine adjoint_fluid_scheme_validate(this)
    class(adjoint_fluid_scheme_t), target, intent(inout) :: this
    ! Variables for retrieving json parameters

    if ((.not. associated(this%u_adj)) .or. &
         (.not. associated(this%v_adj)) .or. &
         (.not. associated(this%w_adj)) .or. &
         (.not. associated(this%p_adj))) then
       call neko_error('Fields are not registered')
    end if

    if ((.not. allocated(this%u_adj%x)) .or. &
         (.not. allocated(this%v_adj%x)) .or. &
         (.not. allocated(this%w_adj%x)) .or. &
         (.not. allocated(this%p_adj%x))) then
       call neko_error('Fields are not allocated')
    end if

    if (.not. allocated(this%ksp_vel)) then
       call neko_error('No Krylov solver for velocity defined')
    end if

    if (.not. allocated(this%ksp_prs)) then
       call neko_error('No Krylov solver for pressure defined')
    end if

    !
    ! Setup checkpoint structure (if everything is fine)
    !
    call this%chkp%init(this%u_adj, this%v_adj, this%w_adj, this%p_adj)

  end subroutine adjoint_fluid_scheme_validate

  !> Apply all boundary conditions defined for velocity
  !! Here we perform additional gs operations to take care of
  !! shared points between elements that have different BCs, as done in Nek5000.
  !! @todo Why can't we call the interface here?
  subroutine adjoint_fluid_scheme_bc_apply_vel(this, t, tstep, strong)
    class(adjoint_fluid_scheme_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    logical, intent(in) :: strong

    call this%bcs_vel%apply_vector(&
         this%u_adj%x, this%v_adj%x, this%w_adj%x, this%dm_Xh%size(), &
         t, tstep, strong)

  end subroutine adjoint_fluid_scheme_bc_apply_vel

  !> Apply all boundary conditions defined for pressure
  !! @todo Why can't we call the interface here?
  subroutine adjoint_fluid_scheme_bc_apply_prs(this, t, tstep)
    class(adjoint_fluid_scheme_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    call this%bcs_prs%apply(this%p_adj, t, tstep)

  end subroutine adjoint_fluid_scheme_bc_apply_prs

  !> Initialize a linear solver
  !! @note Currently only supporting Krylov solvers
  subroutine adjoint_fluid_scheme_solver_factory(ksp, n, solver, &
       max_iter, abstol, monitor)
    class(ksp_t), allocatable, target, intent(inout) :: ksp
    integer, intent(in), value :: n
    character(len=*), intent(in) :: solver
    integer, intent(in) :: max_iter
    real(kind=rp), intent(in) :: abstol
    logical, intent(in) :: monitor

    call krylov_solver_factory(ksp, n, solver, max_iter, abstol, &
         monitor = monitor)

  end subroutine adjoint_fluid_scheme_solver_factory

  !> Initialize a Krylov preconditioner
  subroutine adjoint_fluid_scheme_precon_factory(this, pc, ksp, coef, dof, gs, &
       bclst, pctype, pcparams)
    class(adjoint_fluid_scheme_t), intent(inout) :: this
    class(pc_t), allocatable, target, intent(inout) :: pc
    class(ksp_t), target, intent(inout) :: ksp
    type(coef_t), target, intent(in) :: coef
    type(dofmap_t), target, intent(in) :: dof
    type(gs_t), target, intent(inout) :: gs
    type(bc_list_t), target, intent(inout) :: bclst
    character(len=*) :: pctype
    type(json_file), intent(inout) :: pcparams

    call precon_factory(pc, pctype)

    select type (pcp => pc)
    type is (jacobi_t)
       call pcp%init(coef, dof, gs)
    type is (sx_jacobi_t)
       call pcp%init(coef, dof, gs)
    type is (device_jacobi_t)
       call pcp%init(coef, dof, gs)
    type is (hsmg_t)
       call pcp%init(coef, bclst, pcparams)
    type is (phmg_t)
       call pcp%init(dof%msh, dof%Xh, coef, dof, gs, bclst)
    end select

    call ksp%set_pc(pc)

  end subroutine adjoint_fluid_scheme_precon_factory

  !> Compute CFL
  ! TODO
  ! HARRY
  ! ok this needs to be discussed.
  ! To me, the CFL is based on the convecting velocity, so I would argue
  ! that the CFL of the adjoint is based on
  ! u, v, w
  ! not
  ! u_adj, v_adj, w_adj
  !
  ! Then again, this function is really only used for varaible timestep,
  ! which we can't use if we're checkpointing...
  !
  ! hmmmmm.... requires more thought. Because people optimal ICs or Arnouldi
  ! may use variable timesteps... or us in a steady case..
  !
  ! for now.... let's ignore it
  function adjoint_compute_cfl(this, dt) result(c)
    class(adjoint_fluid_scheme_t), intent(in) :: this
    real(kind=rp), intent(in) :: dt
    real(kind=rp) :: c

    c = cfl(dt, this%u_adj%x, this%v_adj%x, this%w_adj%x, &
         this%Xh, this%c_Xh, this%msh%nelv, this%msh%gdim)

  end function adjoint_compute_cfl

  ! ========================================================================== !
  ! Todo: This section need to be moved

  ! !> Set boundary types for the diagnostic output.
  ! !! @param params The JSON case file.
  ! subroutine adjoint_fluid_scheme_set_bc_type_output(this, params)
  !   class(adjoint_fluid_scheme_t), intent(inout) :: this
  !   type(json_file), intent(inout) :: params
  !   type(dirichlet_t) :: bdry_mask
  !   logical :: found

  !   !
  !   ! Check if we need to output boundaries
  !   !
  !   call json_get_or_default(params, 'case.output_boundary', found, .false.)

  !   if (found) then
  !      call this%bdry%init(this%dm_Xh, 'bdry')
  !      this%bdry = 0.0_rp

  !      call bdry_mask%init_base(this%c_Xh)
  !      call bdry_mask%mark_zone(this%msh%wall)
  !      call bdry_mask%mark_zones_from_list(this%msh%labeled_zones,&
  !           'w', this%bc_labels)
  !      call bdry_mask%finalize()
  !      call bdry_mask%set_g(1.0_rp)
  !      call bdry_mask%apply_scalar(this%bdry%x, this%dm_Xh%size())
  !      call bdry_mask%free()

  !      call bdry_mask%init_base(this%c_Xh)
  !      call bdry_mask%mark_zone(this%msh%inlet)
  !      call bdry_mask%mark_zones_from_list(this%msh%labeled_zones,&
  !           'v', this%bc_labels)
  !      call bdry_mask%finalize()
  !      call bdry_mask%set_g(2.0_rp)
  !      call bdry_mask%apply_scalar(this%bdry%x, this%dm_Xh%size())
  !      call bdry_mask%free()

  !      call bdry_mask%init_base(this%c_Xh)
  !      call bdry_mask%mark_zone(this%msh%outlet)
  !      call bdry_mask%mark_zones_from_list(this%msh%labeled_zones,&
  !           'o', this%bc_labels)
  !      call bdry_mask%finalize()
  !      call bdry_mask%set_g(3.0_rp)
  !      call bdry_mask%apply_scalar(this%bdry%x, this%dm_Xh%size())
  !      call bdry_mask%free()

  !      call bdry_mask%init_base(this%c_Xh)
  !      call bdry_mask%mark_zone(this%msh%sympln)
  !      call bdry_mask%mark_zones_from_list(this%msh%labeled_zones,&
  !           'sym', this%bc_labels)
  !      call bdry_mask%finalize()
  !      call bdry_mask%set_g(4.0_rp)
  !      call bdry_mask%apply_scalar(this%bdry%x, this%dm_Xh%size())
  !      call bdry_mask%free()

  !      call bdry_mask%init_base(this%c_Xh)
  !      call bdry_mask%mark_zone(this%msh%outlet_normal)
  !      call bdry_mask%mark_zones_from_list(this%msh%labeled_zones,&
  !           'on', this%bc_labels)
  !      call bdry_mask%finalize()
  !      call bdry_mask%set_g(5.0_rp)
  !      call bdry_mask%apply_scalar(this%bdry%x, this%dm_Xh%size())
  !      call bdry_mask%free()

  !      call bdry_mask%init_base(this%c_Xh)
  !      call bdry_mask%mark_zone(this%msh%periodic)
  !      call bdry_mask%finalize()
  !      call bdry_mask%set_g(6.0_rp)
  !      call bdry_mask%apply_scalar(this%bdry%x, this%dm_Xh%size())
  !      call bdry_mask%free()
  !   end if

  ! end subroutine adjoint_fluid_scheme_set_bc_type_output

  ! End of section to be moved
  ! ========================================================================== !



  !> Call user material properties routine and update the values of `mu`
  !! if necessary.
  !! @param this The fluid scheme.
  !! @param t Time value.
  !! @param tstep Current time step.
  subroutine adjoint_fluid_scheme_update_material_properties(this, t, tstep)
    class(adjoint_fluid_scheme_t), intent(inout) :: this
    real(kind=rp),intent(in) :: t
    integer, intent(in) :: tstep
    type(field_t), pointer :: nut

    call this%user_material_properties(t, tstep, this%name, &
         this%material_properties)

    if (len(trim(this%nut_field_name)) > 0) then
       nut => neko_field_registry%get_field(this%nut_field_name)
       call field_addcol3(this%mu, nut, this%rho)
    end if

    ! Since mu, rho is a field_t, and we use the %x(1,1,1,1)
    ! host array data to pass constant density and viscosity
    ! to some routines, we need to make sure that the host
    ! values are also filled
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%rho%x, this%rho%x_d, this%rho%size(), &
            DEVICE_TO_HOST, sync=.false.)
       call device_memcpy(this%mu%x, this%mu%x_d, this%mu%size(), &
            DEVICE_TO_HOST, sync=.false.)
    end if
  end subroutine adjoint_fluid_scheme_update_material_properties

  !> Sets rho and mu
  !! @param this The fluid scheme.
  !! @param params The case paramter file.
  !! @param user The user interface.
  subroutine adjoint_fluid_scheme_set_material_properties(this, params, user)
    class(adjoint_fluid_scheme_t), intent(inout) :: this
    type(json_file), intent(inout) :: params
    type(user_t), target, intent(in) :: user
    character(len=LOG_SIZE) :: log_buf
    ! A local pointer that is needed to make Intel happy
    procedure(user_material_properties), pointer :: dummy_mp_ptr
    real(kind=rp) :: const_mu, const_rho


    dummy_mp_ptr => dummy_user_material_properties

    call this%mu%init(this%dm_Xh, "mu")
    call this%rho%init(this%dm_Xh, "rho")
    call this%material_properties%init(2)
    call this%material_properties%assign_to_field(1, this%rho)
    call this%material_properties%assign_to_field(2, this%mu)

    if (.not. associated(user%material_properties, dummy_mp_ptr)) then

       write(log_buf, '(A)') "Material properties must be set in the user&
       & file!"
       call neko_log%message(log_buf)
       this%user_material_properties => user%material_properties

       call user%material_properties(0.0_rp, 0, this%name, &
            this%material_properties)

    else
       this%user_material_properties => dummy_user_material_properties
       ! Incorrect user input
       if (params%valid_path('case.fluid.Re') .and. &
            (params%valid_path('case.fluid.mu') .or. &
            params%valid_path('case.fluid.rho'))) then
          call neko_error("To set the material properties for the fluid, " // &
               "either provide Re OR mu and rho in the case file.")

       else if (params%valid_path('case.fluid.Re')) then
          ! Non-dimensional case
          write(log_buf, '(A)') 'Non-dimensional fluid material properties &
          & input.'
          call neko_log%message(log_buf, lvl = NEKO_LOG_VERBOSE)
          write(log_buf, '(A)') 'Density will be set to 1, dynamic viscosity to&
          & 1/Re.'
          call neko_log%message(log_buf, lvl = NEKO_LOG_VERBOSE)

          ! Read Re into mu for further manipulation.
          call json_get(params, 'case.fluid.Re', const_mu)
          write(log_buf, '(A)') 'Read non-dimensional material properties'
          call neko_log%message(log_buf)
          write(log_buf, '(A,ES13.6)') 'Re         :', const_mu
          call neko_log%message(log_buf)

          ! Set rho to 1 since the setup is non-dimensional.
          const_rho = 1.0_rp
          ! Invert the Re to get viscosity.
          const_mu = 1.0_rp/const_mu
       else
          ! Dimensional case
          call json_get(params, 'case.fluid.mu', const_mu)
          call json_get(params, 'case.fluid.rho', const_rho)
       end if
    end if

    ! We need to fill the fields based on the parsed const values
    ! if the user routine is not used.
    if (associated(user%material_properties, dummy_mp_ptr)) then
       ! Fill mu and rho field with the physical value
       call field_cfill(this%mu, const_mu)
       call field_cfill(this%rho, const_rho)


       write(log_buf, '(A,ES13.6)') 'rho        :', const_rho
       call neko_log%message(log_buf)
       write(log_buf, '(A,ES13.6)') 'mu         :', const_mu
       call neko_log%message(log_buf)
    end if

    ! Since mu, rho is a field_t, and we use the %x(1,1,1,1)
    ! host array data to pass constant density and viscosity
    ! to some routines, we need to make sure that the host
    ! values are also filled
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%rho%x, this%rho%x_d, this%rho%size(), &
            DEVICE_TO_HOST, sync=.false.)
       call device_memcpy(this%mu%x, this%mu%x_d, this%mu%size(), &
            DEVICE_TO_HOST, sync=.false.)
    end if
  end subroutine adjoint_fluid_scheme_set_material_properties

end module adjoint_fluid_scheme
