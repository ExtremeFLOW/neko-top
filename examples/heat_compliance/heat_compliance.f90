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
  use num_types, only: rp

  use objective, only: objective_t

  use design, only: design_t
  use math, only: glsum
  use vector, only: vector_t

  use device, only: device_memcpy, HOST_TO_DEVICE, DEVICE_TO_HOST
  use neko_config, only: NEKO_BCKND_DEVICE

  use mpi_f08, only: mpi_exscan, mpi_sum, MPI_INTEGER, MPI_Allreduce
  use comm, only: pe_rank, pe_size, neko_comm, mpi_real_precision

    use num_types, only: rp, sp
  use field, only: field_t
  use json_module, only: json_file
  use mapping_handler, only: mapping_handler_t
  use coefs, only: coef_t
  use adjoint_fluid_pnpn, only: adjoint_fluid_pnpn_t
  use scratch_registry, only: neko_scratch_registry
  use fld_file_output, only: fld_file_output_t
  use point_zone_registry, only: neko_point_zone_registry
  use point_zone, only: point_zone_t
  use mask_ops, only: mask_exterior_const
  use neko_config, only: NEKO_BCKND_DEVICE
  use device, only: device_memcpy, HOST_TO_DEVICE
  use design, only: design_t
  use math, only: rzero
  use simulation_m, only: simulation_t
  use json_module, only: json_file
  use simple_brinkman_source_term, only: simple_brinkman_source_term_t
  use vector, only: vector_t
  use math, only: copy
  use device_math, only: device_copy
  use field_registry, only: neko_field_registry
  use neko_ext, only: field_to_vector, vector_to_field
  use optimization_ic, only: set_optimization_ic
  use field_math, only: field_rzero
  use json_utils, only: json_get, json_get_or_default, json_get
  use utils, only: neko_error
    use num_types, only: rp
  use json_module, only: json_file
  use field_registry, only: neko_field_registry
  use field, only: field_t
  use coefs, only: coef_t
  use ax_product, only: ax_t, ax_helm_factory
  use krylov, only: ksp_t, ksp_monitor_t, krylov_solver_factory
  use precon, only: pc_t, precon_factory, precon_destroy
  use bc_list, only: bc_list_t
  use neumann, only: neumann_t
  use profiler, only: profiler_start_region, profiler_end_region
  use gather_scatter, only: gs_t, GS_OP_ADD
  use pnpn_residual, only: pnpn_prs_res_t
  use mesh, only: mesh_t, NEKO_MSH_MAX_ZLBLS, NEKO_MSH_MAX_ZLBL_LEN
  use field_registry, only: neko_field_registry
  use mapping, only: mapping_t
  use scratch_registry, only: neko_scratch_registry
  use field_math, only: field_copy, field_add3
  use coefs, only: coef_t
  use logger, only: neko_log, LOG_SIZE
  use neko_config, only: NEKO_BCKND_DEVICE
  use dofmap, only: dofmap_t
  use jacobi, only: jacobi_t
  use device_jacobi, only: device_jacobi_t
  use sx_jacobi, only: sx_jacobi_t
  use utils, only: neko_error
  use device_math, only: device_cfill, device_subcol3, device_cmult
  use json_utils, only: json_get, json_get_or_default
  use math, only: glsc2
  use device_math, only: device_glsc2

  implicit none
  private

  ! ========================================================================== !
  ! Objective: minimum compliance
  type, public, extends(objective_t) :: heat_compliance_t
     ! temperature
     type(field_t) :: phi
     !> Ax
     class(ax_t), allocatable :: Ax
     !> Solver results monitors
     type(ksp_monitor_t) :: ksp_results(1)
     !> Krylov solver
     class(ksp_t), allocatable :: ksp
     !> Preconditioner
     class(pc_t), allocatable :: pc
     !> boundary conditions (they will all be Neumann, so empty)
     type(bc_list_t) :: bclst
     !> tolerance
     real(kind=rp) :: abstol = 0.0000000001_rp
     !> max iterations
     integer :: ksp_max_iter = 200
     !> method for solving PDE
     character(len=5) :: ksp_solver = "gmres"
     ! > preconditioner type
     character(len=4) :: precon_type = "hsmg"
     ! > coef
     type(coef_t), pointer :: coef
     integer :: ksp_n, n, i
   contains
     procedure, public, pass(this) :: heat_compliance_init
     procedure, public, pass(this) :: free => heat_compliance_free
     procedure, public, pass(this) :: update_value => &
          heat_compliance_update_value
     procedure, public, pass(this) :: update_sensitivity => &
          heat_compliance_update_sensitivity
  end type heat_compliance_t

!> A topology optimization design variable
  type, extends(design_t), public :: thermal_conductivity_design_t
     private

     type(field_t), pointer :: design_indicator
     type(field_t), pointer :: thermal_conductivity
     type(field_t), pointer :: sensitivity
     type(mapping_handler_t), public :: mapping
     !> A mask indicating the optimization domain
     class(point_zone_t), pointer :: optimization_domain
     !> A logical if we're restricting the optimization domain
     logical :: has_mask
     type(fld_file_output_t), private :: output
     type(coef_t), pointer :: coef

   contains
     !> Retrieve the design variables
     procedure, pass(this) :: get_values => thermal_conductivity_design_get_design
     !> Retrieve the sensitivity
     procedure, pass(this) :: get_sensitivity => thermal_conductivity_design_get_sensitivity
     !> Initialize the design
     procedure, pass(this) :: init => thermal_conductivity_design_init
     !> Update the design
     procedure, pass(this) :: update_design => thermal_conductivity_design_update_design
     !> map forward
     procedure, pass(this) :: map_forward => thermal_conductivity_design_map_forward
     !> chain rule
     procedure, pass(this) :: map_backward => thermal_conductivity_design_map_backward
     !> chain rule
     procedure, pass(this) :: map_backward_safe => thermal_conductivity_design_map_backward_safe
     ! a writer being called from outside would be nice
     procedure, pass(this) :: write => thermal_conductivity_design_write
     !> Destructor
     procedure, pass(this) :: free => thermal_conductivity_design_free

  end type thermal_conductivity_design_t
contains

  subroutine heat_compliance_init (this, design, coef)
    class(heat_compliance_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(coef_t), target, intent(in) :: coef
    character(len=256), parameter :: name = 'heat_compliance'
    integer :: n

    call this%init_base(name, design%size(), 1.0_rp)
    this%coef => coef

    ! set the number of dofs
    n = this%coef%dof%size()

    ! init the bc list (all Neuman BCs, will remain empty)
    call this%bclst%init()

    ! Setup backend dependent Ax routines
    call ax_helm_factory(this%Ax, full_formulation = .false.)

    ! set up krylov solver
    call krylov_solver_factory(this%ksp, n, this%ksp_solver, &
         this%ksp_max_iter, this%abstol)

    ! set up preconditioner
    call heat_compliance_precon_factory(this%pc, this%ksp, &
         this%coef, this%coef%dof, this%coef%gs_h, this%bclst, &
         this%precon_type)
    
    ! init the temperature field
    call this%phi%init(this%coef%dof)

  end subroutine heat_compliance_init

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

    call this%bclst%free()
    call this%phi%free()

    call this%free_base()
  end subroutine heat_compliance_free

  subroutine heat_compliance_update_value(this, design)
    class(heat_compliance_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    integer :: n, i
    type(field_t), pointer :: RHS
    character(len=LOG_SIZE) :: log_buf
    integer :: temp_indices(1)

select type(design)
  type is (thermal_conductivity_design_t)
    n = this%coef%dof%size()
    call neko_scratch_registry%request_field(RHS, temp_indices(1))

    ! set up Helmholtz operators and RHS
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(this%coef%h1_d, design%thermal_conductivity%x_d, n)
       call device_cmult(this%coef%h1_d, -1.0_rp, n)
       call device_cfill(this%coef%h2_d, 0.0_rp, n)
    else
       ! h1 is now the design variable
       call copy(this%coef%h1, design%thermal_conductivity%x, n)
       ! but negative
       call cmult(this%coef%h1, -1.0_rp, n)
       ! no h2
       this%coef%h2 = 0.0_rp
    end if
    this%coef%ifh2 = .false.
    
    ! RHS is unit forcing
    call field_rone(RHS)

    ! Solve Helmholtz equation for phi
    call profiler_start_region('Forward solve')
    this%ksp_results(1) = &
         this%ksp%solve(this%Ax, this%phi, RHS%x, n, this%coef, &
         this%bclst, this%coef%gs_h)

    call profiler_end_region

    ! update preconditioner (needed?)
    call this%pc%update()

    ! write it all out
    call neko_log%message('Forward problem')

    write(log_buf, '(A,A,A)') 'Iterations:   ',&
         'Start residual:     ', 'Final residual:'
    call neko_log%message(log_buf)
    write(log_buf, '(I11,3x, E15.7,5x, E15.7)') this%ksp_results%iter, &
         this%ksp_results%res_start, this%ksp_results%res_final
    call neko_log%message(log_buf)

    call neko_scratch_registry%relinquish_field(temp_indices)

    ! Now compute the objective
    if (NEKO_BCKND_DEVICE .eq. 1) then
       this%value = device_glsc2(this%phi%x_d, this%coef%B_d, n)
    else
       this%value = glsc2(this%phi%x, this%coef%B, n)
    end if

    class default
     call neko_error("needs a thermal conductivity design type")
  end select
  end subroutine heat_compliance_update_value


  subroutine heat_compliance_update_sensitivity(this, design)
    class(heat_compliance_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(field_t), pointer :: grad_phi_x, grad_phi_y, grad_phi_z
    integer :: temp_indices(3)
    integer :: n


    n = this%coef%dof%size()
    call neko_scratch_registry%request_field(grad_phi_x, temp_indices(1))
    call neko_scratch_registry%request_field(grad_phi_y, temp_indices(2))
    call neko_scratch_registry%request_field(grad_phi_z, temp_indices(3))

    ! dF = grad(phi) . grad(phi_adj)
    ! This problem should be self adjoint?
    call grad(grad_phi_x, grad_phi_y, grad_phi_z, this%phi)
    call field_col2(grad_phi_x, grad_phi_x)
    call field_addcol3(grad_phi_x, grad_phi_y, grad_phi_y)
    call field_addcol3(grad_phi_x, grad_phi_z, grad_phi_z)

    ! now we want to map this backwards through the design
    select type(design)
      type is (thermal_conductivity_design_t)
        call design%map_backward_safe(this%sensitivity, grad_phi_x)
      class default
        call neko_error("Heat compliance requires ad thermal conductivity design type")
    end select

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine heat_compliance_update_sensitivity

  subroutine heat_compliance_precon_factory(pc, ksp, coef, dof, gs, bclst, pctype)

    implicit none
    class(pc_t), allocatable, target, intent(inout) :: pc
    class(ksp_t), target, intent(inout) :: ksp
    type(coef_t), target, intent(in) :: coef
    type(dofmap_t), target, intent(in) :: dof
    type(gs_t), target, intent(inout) :: gs
    type(bc_list_t), target, intent(inout) :: bclst
    character(len=*) :: pctype

    call precon_factory(pc, pctype)

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

  !----------------------------------------------------------------------------!
  !> Initialize the design from a JSON file
  subroutine thermal_conductivity_design_init(this, parameters, coef)
    class(thermal_conductivity_design_t), intent(inout) :: this
    type(json_file), intent(inout) :: parameters
    type(coef_t), target, intent(in) :: coef
    type(json_file) :: json_subdict
    character(len=:), allocatable :: domain_name, domain_type, name
    logical :: dealias

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
       call neko_error('Thermal Conductivity design only supports point_zones for ' // &
            'optimization domain types')

    end select

    this%coef => coef

    ! Initialize the fields
    call neko_field_registry%add_field(coef%dof, "design_indicator", .true.)
      call neko_field_registry%add_field(coef%dof, "thermal_conductivity", .true.)
      call neko_field_registry%add_field(coef%dof, "sensitivity", .true.)
    this%design_indicator => &
         neko_field_registry%get_field("design_indicator")
    this%thermal_conductivity => &
         neko_field_registry%get_field("thermal_conductivity")
    this%sensitivity => &
         neko_field_registry%get_field("sensitivity")

    ! Initialize the output
    call this%output%init(sp, 'design', 3)
    call this%output%fields%assign_to_field(1, this%design_indicator)
    call this%output%fields%assign_to_field(2, this%thermal_conductivity)
    call this%output%fields%assign_to_field(3, this%sensitivity)

    ! Initialize the mapper
      if ('mapping' .in. parameters) then
         call this%mapping%init_base(coef)
         call this%mapping%add(parameters, 'mapping')
      end if

      if ('initial_distribution' .in. parameters) then
         call json_get(parameters, 'initial_distribution', json_subdict)
         call set_optimization_ic(this%design_indicator, coef, coef%gs_h, &
              json_subdict)
      else
         call field_rzero(this%design_indicator)
      end if

    ! Map to the thermal conductivity amplitude
    call this%map_forward()

  end subroutine thermal_conductivity_design_init

  !> Free the design
  subroutine thermal_conductivity_design_free(this)
    class(thermal_conductivity_design_t), intent(inout) :: this

    call this%free_base()
    call this%thermal_conductivity%free()
    call this%design_indicator%free()
    call this%sensitivity%free()

  end subroutine thermal_conductivity_design_free

  subroutine thermal_conductivity_design_update_design(this, values)
    class(thermal_conductivity_design_t), intent(inout) :: this
    type(vector_t), intent(inout) :: values
    integer :: n

    n = this%size()
    call copy(this%design_indicator%x, values%x, n)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(this%design_indicator%x_d, values%x_d, n)
    end if

    call this%map_forward()

    call copy(values%x, this%design_indicator%x, n)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(values%x_d, this%design_indicator%x_d, n)
    end if
  end subroutine thermal_conductivity_design_update_design

  subroutine thermal_conductivity_design_map_forward(this)
    class(thermal_conductivity_design_t), intent(inout) :: this

    if (this%has_mask) then
       call mask_exterior_const(this%design_indicator, &
            this%optimization_domain, 0.0_rp)
    end if

    call this%mapping%apply_forward(this%thermal_conductivity, &
         this%design_indicator)

  end subroutine thermal_conductivity_design_map_forward

  subroutine thermal_conductivity_design_map_backward(this, sensitivity)
    class(thermal_conductivity_design_t), intent(inout) :: this
    type(vector_t), intent(in) :: sensitivity

  end subroutine thermal_conductivity_design_map_backward

  subroutine thermal_conductivity_design_map_backward_safe(this, sens_out, sens_in)
    class(thermal_conductivity_design_t), intent(inout) :: this
    type(vector_t), intent(inout) :: sens_out
    type(field_t), intent(in) :: sens_in
    type(field_t), pointer :: tmp_fld_in, tmp_fld_out
    integer :: temp_indices(2)

    call neko_scratch_registry%request_field(tmp_fld_in, temp_indices(1))
    call neko_scratch_registry%request_field(tmp_fld_out, temp_indices(2))

    call field_copy(tmp_fld_in, sens_in)

    call this%mapping%apply_backward(tmp_fld_out, sens_in)

    if (this%has_mask) then
       call mask_exterior_const(tmp_fld_out, this%optimization_domain, &
            0.0_rp)
    end if

    call field_to_vector(sens_out, tmp_fld_out)

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine thermal_conductivity_design_map_backward_safe

  subroutine thermal_conductivity_design_write(this, idx)
    class(thermal_conductivity_design_t), intent(inout) :: this
    integer, intent(in) :: idx

    call this%output%sample(real(idx, kind=rp))

  end subroutine thermal_conductivity_design_write

  subroutine thermal_conductivity_design_get_design(this, values)
    class(thermal_conductivity_design_t), intent(in) :: this
    type(vector_t), intent(inout) :: values
    integer :: n

    n = this%size()
    call values%init(n)
    call copy(values%x, this%design_indicator%x, n)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(values%x_d, this%design_indicator%x_d, n)
    end if

  end subroutine thermal_conductivity_design_get_design

  subroutine thermal_conductivity_design_get_sensitivity(this, values)
    class(thermal_conductivity_design_t), intent(in) :: this
    type(vector_t), intent(inout) :: values
    integer :: n

    n = this%size()
    call values%init(n)
    call copy(values%x, this%sensitivity%x, n)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(values%x_d, this%sensitivity%x_d, n)
    end if

  end subroutine thermal_conductivity_design_get_sensitivity

end module heat_compliance
