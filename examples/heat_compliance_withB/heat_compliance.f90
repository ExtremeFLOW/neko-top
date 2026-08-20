! Copyright (c) 2024, The Neko Authors
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

module heat_compliance
  use num_types, only: rp, sp
  use objective, only: objective_t
  use design, only: design_t
  use constraint, only: constraint_t
  use vector, only: vector_t
  use vector_math, only: vector_rone, vector_cmult, vector_copy
  use field, only: field_t
  use field_math, only: field_rzero, field_col2, field_addcol3, field_rone, &
       field_copy, field_cmult, field_cfill, field_add2, field_sub2
  use registry, only: neko_registry
  use fld_file_output, only: fld_file_output_t
  use mapping_handler, only: mapping_handler_t
  use coefs, only: coef_t
  use scratch_registry, only: neko_scratch_registry
  use point_zone_registry, only: neko_point_zone_registry
  use point_zone, only: point_zone_t
  use json_module, only: json_file
  use json_utils, only: json_get, json_get_or_default
  use neko_ext, only: field_to_vector, vector_to_field
  use optimization_ic, only: set_optimization_ic
  use device_math, only: device_copy, device_cmult, device_glsc2, device_col2, &
       device_cfill, device_subcol3
  use math, only: col2, cmult, copy, glsc2, glsum
  use ax_product, only: ax_t, ax_helm_factory
  use krylov, only: ksp_t, ksp_monitor_t, krylov_solver_factory
  use precon, only: pc_t, precon_allocator, precon_destroy
  use jacobi, only: jacobi_t
  use device_jacobi, only: device_jacobi_t
  use sx_jacobi, only: sx_jacobi_t
  use hsmg, only: hsmg_t
  use bc_list, only: bc_list_t
  use zero_dirichlet, only: zero_dirichlet_t
  use profiler, only: profiler_start_region, profiler_end_region
  use gather_scatter, only: gs_t, GS_OP_ADD
  use dofmap, only: dofmap_t
  use logger, only: neko_log, LOG_SIZE
  use neko_config, only: NEKO_BCKND_DEVICE
  use operators, only: grad
  use mask_ops, only: mask_exterior_const
  use utils, only: neko_error
  use mpi_f08, only: MPI_SUM, MPI_Allreduce, MPI_INTEGER
  use comm, only: pe_rank, NEKO_COMM, pe_size, MPI_REAL_PRECISION
  implicit none
  private

  type, public, extends(objective_t) :: heat_compliance_t
     ! temperature
     type(field_t) :: phi
     class(ax_t), allocatable :: Ax
     type(ksp_monitor_t) :: ksp_results(1)
     class(ksp_t), allocatable :: ksp
     class(pc_t), allocatable :: pc
     type(bc_list_t) :: bclst
     type(zero_dirichlet_t) :: bc_sink
     real(kind=rp) :: abstol = 1.0e-4_rp
     integer :: ksp_max_iter = 10000
     character(len=2) :: ksp_solver = "cg"
     character(len=6) :: precon_type = "jacobi"
     type(coef_t), pointer :: coef
     type(field_t) :: thermal_conductivity
     type(mapping_handler_t) :: mapping
     integer :: ksp_n, n, i
     real(kind=rp) :: writting_counter
     type(fld_file_output_t), private :: output
     real(kind=rp) :: avg_B
     logical :: enable_output
   contains
     procedure, public :: init_from_attributes => heat_compliance_init
     procedure, public, pass(this) :: heat_compliance_init
     procedure, public, pass(this) :: free => heat_compliance_free
     procedure, public, pass(this) :: update_value => &
          heat_compliance_update_value
     procedure, public, pass(this) :: update_sensitivity => &
          heat_compliance_update_sensitivity
  end type heat_compliance_t

  type, extends(design_t), public :: thermal_conductivity_design_t
     private
     type(field_t), pointer :: design_indicator
     type(field_t), pointer :: sensitivity
     class(point_zone_t), pointer :: optimization_domain
     logical :: has_mask
     type(coef_t), pointer :: coef
   contains
     procedure, pass(this) :: get_values      => thermal_conductivity_design_get_design
     procedure, pass(this) :: get_sensitivity => thermal_conductivity_design_get_sensitivity
     procedure, pass(this) :: init            => thermal_conductivity_design_init
     procedure, pass(this) :: update_design   => thermal_conductivity_design_update_design
     procedure, pass(this) :: map_forward     => thermal_conductivity_design_map_forward
     procedure, pass(this) :: map_backward    => thermal_conductivity_design_map_backward
     procedure, pass(this) :: write           => thermal_conductivity_design_write
     procedure, pass(this) :: free            => thermal_conductivity_design_free
  end type thermal_conductivity_design_t

  !> Volume constraint for thermal_conductivity_design_t
  !!   V = (1/|Omega|) ∫ rho dΩ
  !!   g = limit - V       (or -g if is_max = .true.)
  type, public, extends(constraint_t) :: thermal_volume_constraint_t
     private
     logical :: is_max               ! .false. => V > V_min; .true. => V < V_max
     real(kind=rp) :: limit          ! V_min or V_max
     real(kind=rp) :: volume_domain  ! |Omega|
     type(coef_t), pointer :: coef => null()  ! to access B and volume
     type(mapping_handler_t) :: constraint_mapping
     logical :: has_mapping
   contains
     procedure, public, pass(this) :: init_from_components => &
          thermal_volume_constraint_init
     procedure, public, pass(this) :: free => thermal_volume_constraint_free
     procedure, public, pass(this) :: update_value => &
          thermal_volume_constraint_update_value
     procedure, public, pass(this) :: update_sensitivity => &
          thermal_volume_constraint_update_sensitivity
     procedure, private, pass(this) :: compute_volume => &
          thermal_volume_constraint_compute_volume
  end type thermal_volume_constraint_t

contains
  !=========================================================================!
  !  Objective init
  !=========================================================================!
  subroutine heat_compliance_init (this, design, coef, parameters)
    class(heat_compliance_t), intent(inout) :: this
    class(design_t),          intent(in)    :: design
    type(coef_t), target,     intent(in)    :: coef
    type(json_file),          intent(inout), optional :: parameters
    character(len=256), parameter :: name = 'heat_compliance'
    integer :: n, n_global, ierr
    character(len=LOG_SIZE) :: log_buf
    integer, parameter :: DIRICHLET_ZONE_ID = 2
    type(json_file) :: dummy_json
    real(kind=rp) :: glsumB
    logical :: enable_output


    call this%init_base(name, design%size(), 1.0_rp)
    this%coef => coef
    this%writting_counter = 1.0_rp

    call json_get_or_default(parameters, 'optimization.solver.enable_output', &
         enable_output, .true.)
    this%enable_output = enable_output


    ! Design-to-conductivity mapping
    call this%mapping%init_base(coef)
    if (present(parameters)) then
       if ('mapping' .in. parameters) then
          call this%mapping%add(parameters, 'mapping')
       else
          call parameters%print()
          call neko_error("heat_compliance: missing 'mapping' block in parameters")
       end if
    else
       call neko_error("heat_compliance: parameters JSON required (for mapping + BC init)")
    end if

    ! Number of dofs
    n = this%coef%dof%size()

    ! ------------------------------------------------------------------ !
    ! Boundary conditions: Dirichlet φ = 0 on mesh zone with id = 2
    ! ------------------------------------------------------------------ !
    call this%bclst%init()

    ! zero_dirichlet_t%init(coef, json)
    if (present(parameters)) then
       call this%bc_sink%init(this%coef, parameters)
    else
       call dummy_json%initialize()
       call this%bc_sink%init(this%coef, dummy_json)
       call dummy_json%destroy()
    end if

    ! Mark zone index 2 as Dirichlet sink (φ = 0)
    call this%bc_sink%mark_zone(this%coef%msh%labeled_zones(DIRICHLET_ZONE_ID))

    ! Finalize the BC (build masks etc.)
    call this%bc_sink%finalize()

    ! Add this bc to the list used by Ax/KSP
    call this%bclst%append(this%bc_sink)

    ! ------------------------------------------------------------------ !
    ! Helmholtz operator and Krylov solver
    ! ------------------------------------------------------------------ !
    call ax_helm_factory(this%Ax, full_formulation = .false.)

    call krylov_solver_factory(this%ksp, n, this%ksp_solver, &
         this%ksp_max_iter, this%abstol)

    call heat_compliance_precon_factory(this%pc, this%ksp, &
         this%coef, this%coef%dof, this%coef%gs_h, this%bclst, &
         this%precon_type)

    ! Temperature field and conductivity field
    call this%phi%init(this%coef%dof)
    call this%thermal_conductivity%init(this%coef%dof)

    ! Output fields
    select type(design)
    type is (thermal_conductivity_design_t)
       call this%output%init(sp, 'design', 4)
       call this%output%fields%assign_to_field(1, design%design_indicator)
       call this%output%fields%assign_to_field(2, this%thermal_conductivity)
       call this%output%fields%assign_to_field(3, this%phi)
       call this%output%fields%assign_to_field(4, design%sensitivity)
    class default
       call neko_error("heat_compliance_init: design must be thermal_conductivity_design_t")
    end select

    write(log_buf,'(A)') 'heat_compliance: initialized with Dirichlet sink on zone 2.'
    call neko_log%message(log_buf)

    glsumB = glsum(this%coef%B, n)
    n_global = 0
    call MPI_Allreduce(n, n_global, 1, MPI_INTEGER, &
         MPI_SUM, NEKO_COMM, ierr)
    this%avg_B = glsumB/n_global

  end subroutine heat_compliance_init

  !=========================================================================!
  !  Objective free
  !=========================================================================!
  subroutine heat_compliance_free(this)
    class(heat_compliance_t), intent(inout) :: this

    call this%free_base()

    if (allocated(this%Ax)) then
       deallocate(this%Ax)
    end if

    if (allocated(this%ksp)) then
       call this%ksp%free()
       deallocate(this%ksp)
    end if

    if (allocated(this%pc)) then
       call precon_destroy(this%pc)
       deallocate(this%pc)
    end if

    call this%bc_sink%free()
    call this%mapping%free()
    call this%bclst%free()
    call this%phi%free()
    call this%thermal_conductivity%free()

  end subroutine heat_compliance_free

  !=========================================================================!
  !  Forward solve and objective value
  !=========================================================================!
  subroutine heat_compliance_update_value(this, design)
    class(heat_compliance_t), intent(inout) :: this
    class(design_t),         intent(in)    :: design
    integer :: n
    type(field_t), pointer :: RHS, work, delta_phi
    character(len=LOG_SIZE) :: log_buf
    integer :: temp_indices(3)

    select type(design)
    type is (thermal_conductivity_design_t)
       n = this%coef%dof%size()
       call neko_scratch_registry%request_field(RHS, temp_indices(1), .false.)
       call neko_scratch_registry%request_field(work, temp_indices(2), .false.)
       call neko_scratch_registry%request_field(delta_phi, temp_indices(3), .false.)

       ! Optional masking of design outside optimization region
       if (design%has_mask) then
          call mask_exterior_const(design%design_indicator, &
               design%optimization_domain, 0.0_rp)
       end if

       ! Map design indicator -> physical thermal conductivity
       call this%mapping%apply_forward(this%thermal_conductivity, &
            design%design_indicator)

       ! ----------------------------------------------------------------!
       ! RHS: uniform heating throughout the domain
       ! ----------------------------------------------------------------!
       call field_rone(RHS)             ! RHS = 1 everywhere
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_col2(RHS%x_d, this%coef%B_d, n)     ! apply mass matrix B
       else
          call col2(RHS%x, this%coef%B, n)
       end if

       ! Set Helmholtz coefficients: k(x) from mapped design
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_copy(this%coef%h1_d, this%thermal_conductivity%x_d, n)
          call device_cfill(this%coef%h2_d, 0.0_rp, n)
       else
          call copy(this%coef%h1, this%thermal_conductivity%x, n)
          this%coef%h2 = 0.0_rp
       end if
       this%coef%ifh2 = .false.

       ! We solve for delta phi, so use the previous phi field as initial guess
       call this%Ax%compute(work%x, this%phi%x, this%coef, this%coef%msh, &
         this%coef%Xh)

      call field_sub2(RHS, work)

       ! Gather-scatter on RHS
       call this%coef%gs_h%op(RHS, GS_OP_ADD)

       ! Apply strong BCs to RHS (Dirichlet sink on zone 2)
       call this%bclst%apply_scalar(RHS%x, n)

       ! Solve Helmholtz equation for phi
       call profiler_start_region('Forward solve')
       ! Update preconditioner (if needed)
       call this%pc%update()
       this%ksp_results(1) = &
            this%ksp%solve(this%Ax, delta_phi, RHS%x, n, this%coef, &
            this%bclst, this%coef%gs_h)
       call profiler_end_region('Forward solve')

       ! add result
       call field_add2(this%phi, delta_phi)


       call neko_log%message('Forward problem')
       write(log_buf, '(A,A,A)') 'Iterations:   ', &
            'Start residual:     ', 'Final residual:'
       call neko_log%message(log_buf)
       write(log_buf, '(I11,3x, E15.7,5x, E15.7)') this%ksp_results%iter, &
            this%ksp_results%res_start, this%ksp_results%res_final
       call neko_log%message(log_buf)

       call neko_scratch_registry%relinquish_field(temp_indices)

       ! Objective: compliance = ∫ φ dΩ (mass-matrix weighted)
       if (NEKO_BCKND_DEVICE .eq. 1) then
          this%value = device_glsc2(this%phi%x_d, this%coef%B_d, n)
       else
          this%value = glsc2(this%phi%x, this%coef%B, n)
       end if

    class default
       call neko_error("heat_compliance_update_value: requires thermal_conductivity_design_t")
    end select
  end subroutine heat_compliance_update_value

  !=========================================================================!
  !  Sensitivity (gradient wrt design indicator)
  !=========================================================================!
  subroutine heat_compliance_update_sensitivity(this, design)
    class(heat_compliance_t), intent(inout) :: this
    class(design_t),         intent(in)    :: design
    type(field_t), pointer :: grad_phi_x, grad_phi_y, grad_phi_z
    integer :: temp_indices(3)
    integer :: n

    n = this%coef%dof%size()
    call neko_scratch_registry%request_field(grad_phi_x, temp_indices(1), .false.)
    call neko_scratch_registry%request_field(grad_phi_y, temp_indices(2), .false.)
    call neko_scratch_registry%request_field(grad_phi_z, temp_indices(3), .false.)

    ! Self-adjoint problem: dF/dk ~ |grad(phi)|^2 (up to mapping chain rule)
    call grad(grad_phi_x%x, grad_phi_y%x, grad_phi_z%x, this%phi%x, this%coef)
    call field_col2(grad_phi_x, grad_phi_x)
    call field_addcol3(grad_phi_x, grad_phi_y, grad_phi_y)
    call field_addcol3(grad_phi_x, grad_phi_z, grad_phi_z)

    call field_cmult(grad_phi_x, -1.0_rp)
    !call field_cmult(grad_phi_x, this%avg_B)
       if (NEKO_BCKND_DEVICE .eq. 1) then
           call device_col2(grad_phi_x%x_d, this%coef%B_d, n)
       else
           call col2(grad_phi_x%x, this%coef%B, n)
       end if
       !call this%coef%gs_h%op(grad_phi_x, GS_OP_ADD)
       !if (NEKO_BCKND_DEVICE .eq. 1) then
       !   call device_col2(grad_phi_x%x_d, this%coef%Binv_d, n)
       !else
       !   call col2(grad_phi_x%x, this%coef%Binv, n)
       !end if      

    select type(design)
    type is (thermal_conductivity_design_t)
       ! Chain rule through mapping
       call this%mapping%apply_backward(design%sensitivity, grad_phi_x)
       if (design%has_mask) then
          call mask_exterior_const(design%sensitivity, &
               design%optimization_domain, 0.0_rp)
       end if
       call field_to_vector(this%sensitivity, design%sensitivity)
    class default
       call neko_error("heat_compliance_update_sensitivity: requires thermal_conductivity_design_t")
    end select

    call neko_scratch_registry%relinquish_field(temp_indices)
    if (this%enable_output) then
      call this%output%sample(this%writting_counter)
      this%writting_counter = this%writting_counter + 1
    end if
  end subroutine heat_compliance_update_sensitivity

  !=========================================================================!
  !  Preconditioner factory
  !=========================================================================!
  subroutine heat_compliance_precon_factory(pc, ksp, coef, dof, gs, bclst, pctype)
    implicit none
    class(pc_t),   allocatable, target, intent(inout) :: pc
    class(ksp_t),  target,      intent(inout)         :: ksp
    type(coef_t),  target,      intent(in)            :: coef
    type(dofmap_t),target,      intent(in)            :: dof
    type(gs_t),    target,      intent(inout)         :: gs
    type(bc_list_t),target,     intent(inout)         :: bclst
    character(len=*),           intent(in)            :: pctype

    call precon_allocator(pc, pctype)

    select type (pcp => pc)
    type is (jacobi_t)
       call pcp%init(coef, dof, gs)
    type is (sx_jacobi_t)
       call pcp%init(coef, dof, gs)
    type is (device_jacobi_t)
       call pcp%init(coef, dof, gs)
    end select

    call ksp%set_pc(pc)
  end subroutine heat_compliance_precon_factory

  !=========================================================================!
  !  Design: init
  !=========================================================================!
  subroutine thermal_conductivity_design_init(this, parameters, coef, gs)
    class(thermal_conductivity_design_t), intent(inout) :: this
    type(json_file),                      intent(inout) :: parameters
    type(coef_t),              target,    intent(in)    :: coef
    type(gs_t),    target,      intent(inout)         :: gs
    type(json_file) :: json_subdict
    character(len=:), allocatable :: domain_name, domain_type, name

    call json_get_or_default(parameters, 'name', name, 'Thermal Conductivity Design')
    call json_get_or_default(parameters, 'domain.type', domain_type, 'full')

    select case (trim(domain_type))
    case ('full')
       this%has_mask = .false.
    case ('point_zone')
       this%has_mask = .true.
       call json_get(parameters, 'domain.zone_name', domain_name)
       this%optimization_domain => &
            neko_point_zone_registry%get_point_zone(domain_name)
    case default
       call neko_error('Thermal Conductivity design only supports point_zones for optimization domain types')
    end select

    this%coef => coef

    ! Design and sensitivity fields
    call neko_registry%add_field(coef%dof, "design_indicator", .true.)
    call neko_registry%add_field(coef%dof, "sensitivity",     .true.)
    this%design_indicator => neko_registry%get_field("design_indicator")
    this%sensitivity      => neko_registry%get_field("sensitivity")

    call this%init_base(name, this%design_indicator%dof%size())

    if ('initial_distribution' .in. parameters) then
       call json_get(parameters, 'initial_distribution', json_subdict)
       call set_optimization_ic(this%design_indicator, coef, gs, &
             json_subdict)
    else
       call field_cfill(this%design_indicator, 0.5_rp)
    end if
  end subroutine thermal_conductivity_design_init

  !=========================================================================!
  !  Design: free
  !=========================================================================!
  subroutine thermal_conductivity_design_free(this)
    class(thermal_conductivity_design_t), intent(inout) :: this

    call this%free_base()
    call this%design_indicator%free()
    call this%sensitivity%free()
  end subroutine thermal_conductivity_design_free

  !=========================================================================!
  !  Design: update from optimizer vector
  !=========================================================================!
  subroutine thermal_conductivity_design_update_design(this, values)
    class(thermal_conductivity_design_t), intent(inout) :: this
    type(vector_t),                       intent(inout) :: values
    integer :: n

    n = this%size()

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(this%design_indicator%x_d, values%x_d, n)
    else
       call copy(this%design_indicator%x, values%x, n)
    end if
  end subroutine thermal_conductivity_design_update_design

  subroutine thermal_conductivity_design_map_forward(this)
    class(thermal_conductivity_design_t), intent(inout) :: this
    ! No-op: mapping handled by heat_compliance_t%mapping
  end subroutine thermal_conductivity_design_map_forward

  subroutine thermal_conductivity_design_map_backward(this, sensitivity)
    class(thermal_conductivity_design_t), intent(inout) :: this
    type(vector_t),                       intent(in)    :: sensitivity
    ! No-op: mapping back handled in heat_compliance_update_sensitivity
  end subroutine thermal_conductivity_design_map_backward

  subroutine thermal_conductivity_design_write(this, idx)
    class(thermal_conductivity_design_t), intent(inout) :: this
    integer,                             intent(in)     :: idx
    ! Hook for custom writing of design if you want it
  end subroutine thermal_conductivity_design_write

  !=========================================================================!
  !  Design: get current design vector
  !=========================================================================!
  subroutine thermal_conductivity_design_get_design(this, values)
    class(thermal_conductivity_design_t), intent(in)    :: this
    type(vector_t),                       intent(inout) :: values
    integer :: n

    n = this%size()
    call values%init(n)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(values%x_d, this%design_indicator%x_d, n)
    else
       call copy(values%x, this%design_indicator%x, n)
    end if
  end subroutine thermal_conductivity_design_get_design

  !=========================================================================!
  !  Design: get sensitivity vector
  !=========================================================================!
  subroutine thermal_conductivity_design_get_sensitivity(this, values)
    class(thermal_conductivity_design_t), intent(in)    :: this
    type(vector_t),                       intent(inout) :: values
    integer :: n

    n = this%size()
    call values%init(n)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(values%x_d, this%sensitivity%x_d, n)
    else
       call copy(values%x, this%sensitivity%x, n)
    end if
  end subroutine thermal_conductivity_design_get_sensitivity

  !=========================================================================!
  !  Volume constraint: init
  !=========================================================================!
  !> Direct initializer:
  !!   - design: must be thermal_conductivity_design_t
  !!   - coef  : same coef used by the objective
  !!   - name  : constraint name (for logs)
  !!   - is_max: .false. => V > limit; .true. => V < limit
  !!   - limit : target volume fraction
  subroutine thermal_volume_constraint_init(this, design, coef, name, &
       is_max, limit, parameters)
    class(thermal_volume_constraint_t), intent(inout) :: this
    class(design_t),                    intent(in)    :: design
    type(coef_t),           target,     intent(in)    :: coef
    character(len=*),                  intent(in)    :: name
    logical,                           intent(in)    :: is_max
    real(kind=rp),                     intent(in)    :: limit
    type(json_file),          intent(inout), optional :: parameters
    real(kind=rp) :: avg_B
    integer :: n, n_global, ierr
    type(vector_t) :: beforemapping
    ! Base class init (no separate mask_name)
    call this%init_base(name, design%size(), '')

    this%is_max  = is_max
    this%limit   = limit
    this%coef    => coef

    n = design%size()

    ! Domain volume: use the SEM coef volume
    this%volume_domain = this%coef%volume
    if (pe_rank ==0 ) then
       print *, "volume of full domain =", this%volume_domain
    end if
    
    ! Get the PDE filter if exists as mapping
    this%has_mapping = .false.
    call this%constraint_mapping%init_base(coef)
    if (present(parameters)) then
       if ('constraint_mapping' .in. parameters) then
          call this%constraint_mapping%add(parameters, 'constraint_mapping')
          this%has_mapping = .true.
       end if
    end if

    ! Initialize value and sensitivity
    call this%update_value(design)

        ! Sensitivity: d/d rho = -B / |Omega|   (flip sign if is_max)
        if (NEKO_BCKND_DEVICE .eq. 1) then
           call device_copy(this%sensitivity%x_d, this%coef%B_d, n)
           call device_cmult(this%sensitivity%x_d, -1.0_rp / this%volume_domain, n)
           if (this%is_max) then
              call device_cmult(this%sensitivity%x_d, -1.0_rp, n)
           end if
        else
           call copy(this%sensitivity%x, this%coef%B, n)
           call cmult(this%sensitivity%x, -1.0_rp / this%volume_domain, n)
           if (this%is_max) then
              call cmult(this%sensitivity%x, -1.0_rp, n)
           end if
        end if
   ! Sensitivity: d/d rho = -1 / |Omega|   (flip sign if is_max)
   ! scaled by average B
   !call MPI_Allreduce(n, n_global, 1, MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
   !avg_B = glsum(this%coef%B, n)/n_global
   !call vector_rone(this%sensitivity, n)
   !call vector_cmult(this%sensitivity, -avg_B / this%volume_domain, n)
   !if (this%is_max) then
   !   call vector_cmult(this%sensitivity, -1.0_rp, n)
   !end if
   
   if (this%has_mapping) then
      call beforemapping%init(n)
      call vector_copy(beforemapping, this%sensitivity, n)
      call this%constraint_mapping%apply_backward(this%sensitivity, beforemapping)
      call beforemapping%free()
   end if

   

  end subroutine thermal_volume_constraint_init

  subroutine thermal_volume_constraint_free(this)
    class(thermal_volume_constraint_t), intent(inout) :: this
    call this%free_base()
  end subroutine thermal_volume_constraint_free

  !> Recompute g(design) = limit - V/|Omega|  (flip sign if is_max)
  subroutine thermal_volume_constraint_update_value(this, design)
    class(thermal_volume_constraint_t), intent(inout) :: this
    class(design_t),                    intent(in)    :: design
    real(kind=rp) :: volume

    volume = this%compute_volume(design)

    this%value = this%limit - volume / this%volume_domain
    if (this%is_max) this%value = -this%value

  end subroutine thermal_volume_constraint_update_value

  !> Sensitivity is constant in this implementation (set in init).
  subroutine thermal_volume_constraint_update_sensitivity(this, design)
    class(thermal_volume_constraint_t), intent(inout) :: this
    class(design_t),                    intent(in)    :: design
    ! No-op: this%sensitivity is already correct from init.
  end subroutine thermal_volume_constraint_update_sensitivity

  !> Computes V(design) = ∫ rho dΩ over the whole domain
  ! PDE filter is volume conserving so we can ignore that for forward solve.
  ! But if you have any other mapping for volume constaint (Heavyside or anything)
  ! then you should comment this function and use the one commented in the end.
  function thermal_volume_constraint_compute_volume(this, design) result(volume)
    class(thermal_volume_constraint_t), intent(inout) :: this
    class(design_t),                    intent(in)    :: design
    real(kind=rp) :: volume

    type(vector_t) :: values
    integer :: n

    volume = 0.0_rp

    select type(d => design)
    type is (thermal_conductivity_design_t)
       n = d%size()
       call d%get_values(values)

       if (NEKO_BCKND_DEVICE .eq. 1) then
          volume = device_glsc2(values%x_d, this%coef%B_d, n)
       else
          volume = glsc2(values%x, this%coef%B, n)
       end if

       call values%free()

    class default
       call neko_error('thermal_volume_constraint_compute_volume: requires thermal_conductivity_design_t')
    end select

  end function thermal_volume_constraint_compute_volume

  !> Computes V(design) = ∫ rho dΩ over the whole domain
   !   function thermal_volume_constraint_compute_volume(this, design) result(volume)
   !     class(thermal_volume_constraint_t), intent(inout) :: this
   !     class(design_t),                    intent(in)    :: design
   !     real(kind=rp) :: volume

   !     type(vector_t) :: values, values_aftermapping
   !     integer :: n

   !     volume = 0.0_rp

   !     select type(d => design)
   !     type is (thermal_conductivity_design_t)
   !        n = d%size()
   !        call d%get_values(values)
   !        if (this%has_mapping) then
   !           call values_aftermapping%init(n)
   !           call this%constraint_mapping%apply_forward(values_aftermapping, values)
   !           call vector_copy(values, values_aftermapping, n)
   !           call values_aftermapping%free()
   !        end if

   !        if (NEKO_BCKND_DEVICE .eq. 1) then
   !           volume = device_glsc2(values%x_d, this%coef%B_d, n)
   !        else
   !           volume = glsc2(values%x, this%coef%B, n)
   !        end if

   !        call values%free()

   !     class default
   !        call neko_error('thermal_volume_constraint_compute_volume: requires thermal_conductivity_design_t')
   !     end select

   !   end function thermal_volume_constraint_compute_volume

end module heat_compliance
