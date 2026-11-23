!> @file adjoint_fluid_scheme.f90
!! @copyright
!! Copyright (c) 2024-2025, The Neko-TOP Authors
!! All rights reserved.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions
!! are met:
!!
!!   * Redistributions of source code must retain the above copyright
!!     notice, this list of conditions and the following disclaimer.
!!
!!   * Redistributions in binary form must reproduce the above
!!     copyright notice, this list of conditions and the following
!!     disclaimer in the documentation and/or other materials provided
!!     with the distribution.
!!
!!   * Neither the name of the authors nor the names of its
!!     contributors may be used to endorse or promote products derived
!!     from this software without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
!! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
!! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
!! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
!! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
!! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
!! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!! POSSIBILITY OF SUCH DAMAGE.
!
module adjoint_fluid_scheme
  use gather_scatter, only: gs_t
  use checkpoint, only: chkp_t
  use num_types, only: rp
  use field, only: field_t
  use space, only: space_t, GL
  use dofmap, only: dofmap_t
  use coefs, only: coef_t
  use dirichlet, only: dirichlet_t
  use precon, only: pc_t, precon_allocator
  use fluid_stats, only: fluid_stats_t
  use bc_list, only: bc_list_t
  use mesh, only: mesh_t, NEKO_MSH_MAX_ZLBL_LEN
  use time_state, only: time_state_t
  use time_scheme_controller, only: time_scheme_controller_t
  use logger, only: LOG_SIZE
  use json_module, only: json_file
  use user_intf, only: user_t, user_material_properties_intf
  use field_series, only: field_series_t
  use time_step_controller, only: time_step_controller_t
  use field_list, only : field_list_t
  use interpolation, only: interpolator_t
  use scratch_registry, only : scratch_registry_t
  use utils, only: neko_error

  implicit none
  private
  public :: adjoint_fluid_scheme_t, adjoint_fluid_scheme_factory

  !> Base type of all fluid formulations.
  type, abstract :: adjoint_fluid_scheme_t
     !> A name that can be used to distinguish this solver in e.g. user routines
     character(len=:), allocatable :: name

     type(space_t) :: Xh !< Function space \f$ X_h \f$
     type(dofmap_t) :: dm_Xh !< Dofmap associated with \f$ X_h \f$
     type(gs_t) :: gs_Xh !< Gather-scatter associated with \f$ X_h \f$
     type(coef_t) :: c_Xh !< Coefficients associated with \f$ X_h \f$

     ! Tim. This will need to be refactored, but since we have so many extra
     ! terms that involve products of variables, we often need to increase our
     ! quadrature. So it's natural to have an over integration coef and a
     ! way of converting between them
     type(space_t) :: Xh_GL !< Function space \f$ X_h_GL \f$
     type(dofmap_t) :: dm_Xh_GL !< Dofmap associated with \f$ X_h_GL \f$
     type(gs_t) :: gs_Xh_GL !< Gather-scatter associated with \f$ X_h_GL \f$
     type(coef_t) :: c_Xh_GL !< Coefficients associated with \f$ X_h_GL \f$
     !> Interpolator between the original and higher-order spaces
     type(interpolator_t) :: GLL_to_GL
     !> Scratch registry on the GL space
     type(scratch_registry_t) :: scratch_GL
     !> Has the Gauss over-integration stack above been built?
     logical :: gauss_initialised = .false.

     type(time_scheme_controller_t), allocatable :: ext_bdf

     !> The velocity field
     type(field_t), pointer :: u_adj => null() !< x-component of Velocity
     type(field_t), pointer :: v_adj => null() !< y-component of Velocity
     type(field_t), pointer :: w_adj => null() !< z-component of Velocity
     type(field_t), pointer :: p_adj => null() !< Pressure
     type(field_series_t) :: ulag, vlag, wlag !< fluid field (lag)

     type(field_t), pointer :: u_b => null() !< x-component of baseflow velocity
     type(field_t), pointer :: v_b => null() !< y-component of baseflow Velocity
     type(field_t), pointer :: w_b => null() !< z-component of baseflow Velocity
     type(field_t), pointer :: p_b => null() !< Baseflow pressure

     !> Checkpoint
     type(chkp_t), pointer :: chkp => null()

     !> X-component of the right-hand side.
     type(field_t), pointer :: f_adj_x => null()
     !> Y-component of the right-hand side.
     type(field_t), pointer :: f_adj_y => null()
     !> Z-component of the right-hand side.
     type(field_t), pointer :: f_adj_z => null()

     !> Boundary conditions
     ! List of boundary conditions for pressure
     type(bc_list_t) :: bcs_prs
     ! List of boundary conditions for velocity
     type(bc_list_t) :: bcs_vel

     type(json_file), pointer :: params !< Parameters
     type(mesh_t), pointer :: msh => null() !< Mesh

     !> Boundary condition labels (if any)
     character(len=NEKO_MSH_MAX_ZLBL_LEN), allocatable :: bc_labels(:)

     !> Density field
     type(field_t) :: rho

     !> The dynamic viscosity
     type(field_t) :: mu

     !> A helper that packs material properties to pass to the user routine.
     type(field_list_t) :: material_properties

     !> Is the fluid frozen at the moment
     logical :: freeze = .false.

     !> User material properties routine
     procedure(user_material_properties_intf), nopass, pointer :: &
          user_material_properties => null()

   contains
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

     !> Set the user inflow
     procedure(validate_intrf), pass(this), deferred :: validate
     !> Compute the CFL number
     procedure(fluid_scheme_base_compute_cfl_intrf), pass(this), deferred :: compute_cfl
     !> Set rho and mu
     procedure(update_material_properties), pass(this), deferred :: update_material_properties
     !> Build the Gauss over-integration stack on first request
     procedure, pass(this), non_overridable :: require_gauss_stack => &
          adjoint_fluid_scheme_require_gauss_stack
  end type adjoint_fluid_scheme_t

  !> Initialize all fields
  abstract interface
     subroutine adjoint_fluid_init_all_intrf(this, msh, lx, params, kspv_init, &
          kspp_init, scheme, user)
       import adjoint_fluid_scheme_t
       import mesh_t
       import json_file
       import user_t
       class(adjoint_fluid_scheme_t), target, intent(inout) :: this
       type(mesh_t), target, intent(inout) :: msh
       integer, intent(inout) :: lx
       type(json_file), target, intent(inout) :: params
       logical, intent(inout) :: kspv_init
       logical, intent(inout) :: kspp_init
       type(user_t), target, intent(in) :: user
       character(len=*), intent(in) :: scheme
     end subroutine adjoint_fluid_init_all_intrf
  end interface

  !> Initialize common data for the current scheme
  abstract interface
     subroutine adjoint_fluid_init_common_intrf(this, msh, lx, params, scheme, &
          user, kspv_init)
       import adjoint_fluid_scheme_t
       import mesh_t
       import json_file
       import user_t
       class(adjoint_fluid_scheme_t), target, intent(inout) :: this
       type(mesh_t), target, intent(inout) :: msh
       integer, intent(inout) :: lx
       character(len=*), intent(in) :: scheme
       type(json_file), target, intent(inout) :: params
       type(user_t), target, intent(in) :: user
       logical, intent(in) :: kspv_init
     end subroutine adjoint_fluid_init_common_intrf
  end interface

  !> Deallocate a fluid formulation
  abstract interface
     subroutine adjoint_fluid_free_intrf(this)
       import adjoint_fluid_scheme_t
       class(adjoint_fluid_scheme_t), intent(inout) :: this
     end subroutine adjoint_fluid_free_intrf
  end interface

  !> Abstract interface to initialize an adjoint formulation
  abstract interface
     subroutine adjoint_fluid_scheme_init_intrf(this, msh, lx, params, user, &
          chkp)
       import adjoint_fluid_scheme_t
       import json_file
       import mesh_t
       import user_t
       import chkp_t
       class(adjoint_fluid_scheme_t), target, intent(inout) :: this
       type(mesh_t), target, intent(inout) :: msh
       integer, intent(in) :: lx
       type(json_file), target, intent(inout) :: params
       type(user_t), target, intent(in) :: user
       type(chkp_t), target, intent(inout) :: chkp
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
     subroutine adjoint_fluid_scheme_step_intrf(this, time, dt_controller)
       import time_state_t
       import adjoint_fluid_scheme_t
       import time_step_controller_t
       class(adjoint_fluid_scheme_t), target, intent(inout) :: this
       type(time_state_t), intent(in) :: time
       type(time_step_controller_t), intent(in) :: dt_controller
     end subroutine adjoint_fluid_scheme_step_intrf
  end interface

  !> Abstract interface to restart an adjoint scheme
  abstract interface
     subroutine adjoint_fluid_scheme_restart_intrf(this, chkp)
       import adjoint_fluid_scheme_t
       import chkp_t
       class(adjoint_fluid_scheme_t), target, intent(inout) :: this
       type(chkp_t), intent(inout) :: chkp
     end subroutine adjoint_fluid_scheme_restart_intrf
  end interface

  !> Abstract interface to setup boundary conditions
  abstract interface
     subroutine adjoint_fluid_scheme_setup_bcs_intrf(this, user, params)
       import adjoint_fluid_scheme_t, user_t, json_file
       class(adjoint_fluid_scheme_t), target, intent(inout) :: this
       type(user_t), target, intent(in) :: user
       type(json_file), intent(inout) :: params
     end subroutine adjoint_fluid_scheme_setup_bcs_intrf
  end interface

  !> Abstract interface to validate the user inflow
  abstract interface
     subroutine validate_intrf(this)
       import adjoint_fluid_scheme_t
       class(adjoint_fluid_scheme_t), target, intent(inout) :: this
     end subroutine validate_intrf
  end interface

  !> Abstract interface to sets rho and mu
  abstract interface
     subroutine update_material_properties(this, time)
       import adjoint_fluid_scheme_t, time_state_t
       class(adjoint_fluid_scheme_t), intent(inout) :: this
       type(time_state_t), intent(in) :: time
     end subroutine update_material_properties
  end interface

  !> Compute the CFL number
  abstract interface
     function fluid_scheme_base_compute_cfl_intrf(this, dt) result(c)
       import adjoint_fluid_scheme_t
       import rp
       class(adjoint_fluid_scheme_t), intent(in) :: this
       real(kind=rp), intent(in) :: dt
       real(kind=rp) :: c
     end function fluid_scheme_base_compute_cfl_intrf
  end interface

  interface
     !> Initialise an adjoint scheme
     module subroutine adjoint_fluid_scheme_factory(object, type_name)
       class(adjoint_fluid_scheme_t), intent(inout), allocatable :: object
       character(len=*) :: type_name
     end subroutine adjoint_fluid_scheme_factory
  end interface

contains

  !> Build the Gauss over-integration stack on first request.
  !! @details Constructs `Xh_GL`, `dm_Xh_GL`, `gs_Xh_GL`, `c_Xh_GL`,
  !! `GLL_to_GL` and `scratch_GL`. The stack is only needed by consumers that
  !! over-integrate, so it is built on demand rather than in `init_base`.
  !! The routine is idempotent: a repeated call is a no-op, leaving any
  !! pointers already captured into the stack valid.
  !! @note `dofmap_t%init`, `gs_t%init` and `coef_t%init` are collective, so
  !! every call site must be reached by all ranks.
  !! @note Must be called after `init_base`, which associates `msh` and
  !! initialises `Xh`; both are read here.
  !! @param this The adjoint fluid scheme owning the stack.
  subroutine adjoint_fluid_scheme_require_gauss_stack(this)
    class(adjoint_fluid_scheme_t), target, intent(inout) :: this
    integer :: lxd

    if (this%gauss_initialised) return

    if (.not. associated(this%msh)) then
       call neko_error('require_gauss_stack called before init_base')
    end if

    ! over intergration order (hard coded now, should be optional)
    lxd = (3 * (this%Xh%lx + 1)) / 2

    call this%Xh_GL%init(GL, lxd, lxd, lxd)

    call this%dm_Xh_GL%init(this%msh, this%Xh_GL)

    call this%gs_Xh_GL%init(this%dm_Xh_GL)

    call this%c_Xh_GL%init(this%gs_Xh_GL)

    call this%GLL_to_GL%init(this%Xh_GL, this%Xh)

    ! Overintegration scratch registry (5 should be sufficient)
    call this%scratch_GL%init(5, 2, this%dm_Xh_GL)

    this%gauss_initialised = .true.

  end subroutine adjoint_fluid_scheme_require_gauss_stack

end module adjoint_fluid_scheme
