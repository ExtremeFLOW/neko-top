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

! Implements the `adjoint_case_t` type.
module adjoint_case
  use num_types, only: rp, dp
  use json_module, only: json_file
  use json_utils, only: json_get, json_get_or_default
  use simulation_component, only: simulation_component_t
  use case, only: case_t
  use field, only: field_t
  use coefs, only: coef_t
  use field_registry, only: neko_field_registry
  use scratch_registry, only: neko_scratch_registry
  use adjoint_pnpn, only: adjoint_pnpn_t
  use adjoint_output, only: adjoint_output_t
  use neko_config, only: NEKO_BCKND_DEVICE
  use field_math, only: field_cfill, field_sub2, field_copy, field_glsc2, &
       field_glsc3
  use field_math, only: field_add2
  use math, only: glsc2, glsc3
  use device_math, only: device_glsc2
  use adv_lin_no_dealias, only: adv_lin_no_dealias_t
  use logger, only: neko_log, LOG_SIZE
  use adjoint_scheme, only: adjoint_scheme_t
  use adjoint_fctry, only: adjoint_scheme_factory
  use time_step_controller, only: time_step_controller_t
  use time_scheme_controller, only: time_scheme_controller_t
  use mpi_f08, only: MPI_WTIME
  use jobctrl, only: jobctrl_time_limit
  use profiler, only: profiler_start, profiler_stop, profiler_start_region, &
       profiler_end_region
  use file, only: file_t
  use num_types, only: rp, sp, dp
  use fluid_scheme, only: fluid_scheme_factory
  use fluid_pnpn, only: fluid_pnpn_t
  use fluid_scheme, only: fluid_scheme_t
  use fluid_output, only: fluid_output_t
  use chkp_output, only: chkp_output_t
  use mean_sqr_flow_output, only: mean_sqr_flow_output_t
  use mean_flow_output, only: mean_flow_output_t
  use fluid_stats_output, only: fluid_stats_output_t
  use mpi_f08, only: MPI_COMM_WORLD
  use mesh_field, only: mesh_fld_t, mesh_field_init, mesh_field_free
  use parmetis, only: parmetis_partmeshkway
  use redist, only: redist_mesh
  use output_controller, only: output_controller_t
  use flow_ic, only: set_flow_ic
  use scalar_ic, only: set_scalar_ic
  use field, only: field_t
  use field_registry, only: neko_field_registry
  use stats, only: stats_t
  use file, only: file_t
  use utils, only: neko_error
  use mesh, only: mesh_t
  use time_scheme_controller, only: time_scheme_controller_t
  use logger, only: neko_log, NEKO_LOG_QUIET, LOG_SIZE
  use jobctrl, only: jobctrl_set_time_limit
  use user_intf, only: user_t
  use scalar_pnpn, only: scalar_pnpn_t
  use json_module, only: json_file, json_core, json_value
  use json_utils, only: json_get, json_get_or_default
  use scratch_registry, only: scratch_registry_t, neko_scratch_registry
  use point_zone_registry, only: neko_point_zone_registry
  use adjoint_ic, only: set_adjoint_ic
  use json_utils, only: json_extract_item
  use json_utils_ext, only: json_key_fallback, json_get_subdict
  use adjoint_scalar_scheme, only: adjoint_scalar_scheme_t
  use adjoint_scalar_pnpn, only: adjoint_scalar_pnpn_t
  use adjoint_scalar_convection_source_term, only: &
  adjoint_scalar_convection_source_term_t
  use usr_scalar, only : usr_scalar_t, usr_scalar_bc_eval
  use, intrinsic :: iso_fortran_env, only: stderr => error_unit
  implicit none
  private
  public :: adjoint_case_t, adjoint_init, adjoint_free

  !> Adjoint case type.
  !! Todo: This should Ideally be a subclass of case_t, however, that is not yet
  !! suppoerted by Neko.
  type :: adjoint_case_t

     ! TODO
     ! I think it would be nicer if this was called `fluid` just like the 
     ! forward.
     !
     ! so you're acessing forward%fluid
     ! and also           adjoint%fluid
     !
     ! and same for the scalar.
     ! I'll update this in a different PR
     class(adjoint_scheme_t), allocatable :: scheme
     ! I don't know why this is the pnpn type, and not the scheme...
     ! but that's how it is in neko
     type(adjoint_scalar_pnpn_t), allocatable :: scalar
     ! this is the extra term for adjoint convection
     type(adjoint_scalar_convection_source_term_t), allocatable :: &
     adjoint_convection_term
     type(case_t), pointer :: case

     ! Fields
     real(kind=rp) :: tol
     type(adjoint_output_t) :: f_out
     type(output_controller_t) :: output_controller

     logical :: have_scalar = .false.

     ! this is to force 'w' on user_bcs
     procedure(usr_scalar_bc_eval), nopass, pointer :: force_user_bc_w => &
     adjoint_force_user_bc_w

  end type adjoint_case_t

  interface adjoint_init
     module procedure adjoint_init_from_json, adjoint_init_from_attributes
  end interface adjoint_init

contains

  ! Constructor from json.
  subroutine adjoint_init_from_json(this, neko_case)
    class(adjoint_case_t), intent(inout) :: this
    type(case_t), target, intent(inout) :: neko_case

    this%case => neko_case

    ! Read the tolerance
    call json_get_or_default(neko_case%params, "tol", this%tol, 1.0e-6_rp)

    ! Check if the scalar field is allocated
    if (allocated(neko_case%scalar)) then
       this%have_scalar = .true.
    end if

    call adjoint_case_init_common(this, neko_case)

  end subroutine adjoint_init_from_json

  ! Constructor from attributes
  subroutine adjoint_init_from_attributes(this, neko_case, tol)
    class(adjoint_case_t), intent(inout) :: this
    class(case_t), intent(inout), target :: neko_case
    real(kind=rp), intent(in) :: tol

    this%case => neko_case
    this%tol = tol

    ! Check if the scalar field is allocated
    ! TODO
    ! Man this is tricky... I agree that in principle, if we have a passive 
    ! scalar in the forward then we should have a passive scalar in the adjoint
    !
    ! But if you look closely in Caspers paper, they have an objective and a
    ! pressure drop constraint.
    !
    ! And the adjoint for the pressure drop constraint doesn't involve the 
    ! passive scalar.
    !
    ! I'm sure with enough thought we can figure out a smart way of handing
    ! this. 
    ! For now...
    ! I think it's ok to assume passive scalar forward => passive scalar adjoint
    !
    ! But just putting this note here so we know to come back and discuss this
    ! point.
    if (allocated(neko_case%scalar)) then
       this%have_scalar = .true.
    end if

    call adjoint_case_init_common(this, neko_case)

  end subroutine adjoint_init_from_attributes

  !> Initialize a neko_case from its (loaded) params object
  subroutine adjoint_case_init_common(this, neko_case)
    class(adjoint_case_t), intent(inout) :: this
    type(case_t), intent(inout) :: neko_case
    character(len=:), allocatable :: output_directory
    integer :: lx = 0
    logical :: scalar = .false.
    type(file_t) :: msh_file, bdry_file, part_file
    logical :: found, logical_val
    integer :: integer_val
    real(kind=rp) :: real_val
    character(len=:), allocatable :: string_val
    real(kind=rp) :: stats_start_time, stats_output_val
    integer :: stats_sampling_interval
    integer :: output_dir_len
    integer :: precision

    ! extra things for json
    type(json_file) :: ic_json
    character(len=:), allocatable :: json_key
    !
    ! Setup fluid scheme
    !
    call json_get(neko_case%params, 'case.fluid.scheme', string_val)
    call adjoint_scheme_factory(this%scheme, trim(string_val))

    call json_get(neko_case%params, 'case.numerics.polynomial_order', lx)
    lx = lx + 1 ! add 1 to get number of gll points
    call this%scheme%init(neko_case%msh, lx, neko_case%params, neko_case%usr, &
         neko_case%ext_bdf)

    !
    ! Setup scalar scheme
    !
    ! TODO
    ! Tims forward vs adjoint JSON
     if (neko_case%params%valid_path('case.scalar')) then
        call json_get_or_default(neko_case%params, 'case.scalar.enabled', &
     scalar,                             .true.)
     end if

     if (scalar) then
        allocate(this%scalar)
        ! TODO
        ! this%scalar%chkp%tlag => this%tlag
        ! this%scalar%chkp%dtlag => this%dtlag
        !
        ! I don't know what this is yet
        call this%scalar%init(neko_case%msh, this%scheme%c_Xh, &
             this%scheme%gs_Xh, neko_case%params, neko_case%usr,&
             this%scheme%ulag, this%scheme%vlag, this%scheme%wlag, &
             neko_case%ext_bdf, this%scheme%rho)
        call this%scheme%chkp%add_scalar(this%scalar%s_adj)

        ! TODO
        ! we don't have checkpoints yet
        ! this%scheme%chkp%abs1 => this%scalar%abx1
        ! this%scheme%chkp%abs2 => this%scalar%abx2
        ! this%scheme%chkp%slag => this%scalar%slag

        ! TODO HUGE HUGE TODO
        ! So if we have a passive scalar we also get a source term entering
        ! the adjoint velocity equation which arises when you linearize the
        ! the convective term in passive scalar equation.
        !
        ! $\phi^\dagger \nabla \phi$
        ! 
        ! I'm SOOOOO worried I have the sign the wrong way around.
        ! We really have to write the adjoint derivation nicely.
        !
        ! for now I'm assuming in our adjoint derivation we ADD all the 
        ! equations together.
        ! - So it starts as being positive on the LHS
        ! - if we treat this term as a source term it goes on the RHS, so now
        !   it's negative on the RHS.
        !
        ! I checked through Casper's adjoint equations and the first term
        ! after the = sign of eq (14) looks like the term I'm talking about.
        ! And his is negative too.
        ! So I THINK this is correct, but we need to double check.

        ! and it should be appended to the adjoint velocity

        ! allocate it
        allocate(this%adjoint_convection_term)
        ! initialize it
        call this%adjoint_convection_term%init_from_components( &
             this%scheme%f_adj_x, this%scheme%f_adj_y, &
             this%scheme%f_adj_z, &
             neko_case%scalar%s, &
             this%scalar%s_adj, this%scheme%c_Xh)
        ! append it to the adjoint velocity equation
        call this%scheme%source_term%add(this%adjoint_convection_term)
     end if

    ! TODO
    ! I don't really know how we should handle this...
    ! I feel like, even if someone puts a strange BC it will be diriclette
    ! so the adjoint will go to 'w'.
    !
    ! What we really need is a robust way to set BCs based on objective 
    ! functions, because that's when we're going to get weird unique BCs for
    ! the adjoint.
    !
    !Setup user defined conditions
    !
    if (neko_case%params%valid_path('case.fluid.inflow_condition')) then
       call json_get(neko_case%params, 'case.fluid.inflow_condition.type', &
            string_val)
       if (trim(string_val) .eq. 'user') then
          call this%scheme%set_usr_inflow(neko_case%usr%fluid_user_if)
       end if
    end if

    ! Setup user boundary conditions for the scalar.
    if (scalar) then
       ! TODO
       ! this is not very clear what's going on, in the standard scalar you
       ! would call:
       ! call this%scalar%set_user_bc(neko_case%usr%scalar_user_bc)
       ! Here, we've defined (the last subroutine) essentially a user_bc but
       ! it always forces to zero.
       call this%scalar%set_user_bc(this%force_user_bc_w)
    end if

    !
    ! Setup initial conditions
    !
    json_key = json_key_fallback(neko_case%params, &
         'case.adjoint.initial_condition', 'case.fluid.initial_condition')

    call json_get(neko_case%params, json_key//'.type', string_val)
    call json_get_subdict(neko_case%params, json_key, ic_json)

    if (trim(string_val) .ne. 'user') then
       call set_adjoint_ic( &
            this%scheme%u_adj, this%scheme%v_adj, this%scheme%w_adj, &
            this%scheme%p_adj, this%scheme%c_Xh, this%scheme%gs_Xh, &
            string_val, ic_json)
    else
       call set_adjoint_ic( &
            this%scheme%u_adj, this%scheme%v_adj, this%scheme%w_adj, &
            this%scheme%p_adj, this%scheme%c_Xh, this%scheme%gs_Xh, &
            neko_case%usr%fluid_user_ic, ic_json)
    end if

    if (scalar) then
       ! TODO
       ! perhaps we should consider: 
       ! `adjoint` -> `adjoint_fluid` in the casefile?
       json_key = json_key_fallback(neko_case%params, &
            'case.adjoint_scalar.initial_condition', &
            'case.scalar.initial_condition')

       call json_get(neko_case%params, json_key//'.type', string_val)
       call json_get_subdict(neko_case%params, json_key, ic_json)
       if (trim(string_val) .ne. 'user') then
          call set_scalar_ic(this%scalar%s_adj, &
            this%scalar%c_Xh, this%scalar%gs_Xh, string_val, ic_json)
       else
       	 ! TODO
       	 ! this is wrong, it's point to the neko case.
       	 ! I think we need to discus the case file before we can do user
       	 ! defined stuff. 
       	 ! But I guess we should ALSO have user functionality for the adjoint
       	 !
       !   call set_scalar_ic(this%scalar%s_adj, &
       !     this%scalar%c_Xh, this%scalar%gs_Xh, &
       !     neko_case%usr%scalar_user_ic, ic_json)
       end if
    end if

    ! Add initial conditions to BDF scheme (if present)
    select type (f => this%scheme)
      type is (adjoint_pnpn_t)
       call f%ulag%set(f%u_adj)
       call f%vlag%set(f%v_adj)
       call f%wlag%set(f%w_adj)
    end select

    !
    ! Validate that the neko_case is properly setup for time-stepping
    !
    call this%scheme%validate

    if (scalar) then
       call this%scalar%slag%set(this%scalar%s_adj)
       call this%scalar%validate
    end if

    !
    ! Setup output precision of the field files
    !
    call json_get_or_default(neko_case%params, 'case.output_precision', &
         string_val, 'single')

    if (trim(string_val) .eq. 'double') then
       precision = dp
    else
       precision = sp
    end if

    !
    ! Setup output_controller
    !
    call this%output_controller%init(neko_case%end_time)
    if (scalar) then
       this%f_out = adjoint_output_t(precision, this%scheme, this%scalar, &
            path = trim(output_directory))
    else
       this%f_out = adjoint_output_t(precision, this%scheme, &
            path = trim(output_directory))
    end if

    call json_get_or_default(neko_case%params, 'case.fluid.output_control',&
         string_val, 'org')

    if (trim(string_val) .eq. 'org') then
       ! yes, it should be real_val below for type compatibility
       call json_get(neko_case%params, 'case.nsamples', real_val)
       call this%output_controller%add(this%f_out, real_val, 'nsamples')
    else if (trim(string_val) .eq. 'never') then
       ! Fix a dummy 0.0 output_value
       call json_get_or_default(neko_case%params, 'case.fluid.output_value', &
            real_val, 0.0_rp)
       call this%output_controller%add(this%f_out, 0.0_rp, string_val)
    else
       call json_get(neko_case%params, 'case.fluid.output_value', real_val)
       call this%output_controller%add(this%f_out, real_val, string_val)
    end if

    ! !
    ! ! Save checkpoints (if nothing specified, default to saving at end of sim)
    ! !
    ! call json_get_or_default(neko_case%params, 'case.output_checkpoints',&
    !      logical_val, .true.)
    ! if (logical_val) then
    !    call json_get_or_default(neko_case%params, 'case.checkpoint_format', &
    !         string_val, "chkp")
    !   neko_case%f_chkp = chkp_output_t(this%scheme%chkp, &
    ! path = output_directory, &
    !         ! fmt = trim(string_val))
    !    call json_get_or_default(neko_case%params, 'case.checkpoint_control', &
    !         string_val, "simulationtime")
    !    call json_get_or_default(neko_case%params, 'case.checkpoint_value', &
    ! real_val,&
    !         1e10_rp)
    !   !  call this%output_controller%add(neko_case%f_chkp, real_val, string_val)
    ! end if

  end subroutine adjoint_case_init_common

  ! Destructor.
  subroutine adjoint_free(this)
    class(adjoint_case_t), intent(inout) :: this

    nullify(this%case)
    call this%scheme%free()
    call this%output_controller%free()

  end subroutine adjoint_free

  subroutine adjoint_force_user_bc_w(s, x, y, z, nx, ny, nz, &                    
                                   ix, iy, iz, ie, t, tstep)                    
       real(kind=rp), intent(inout) :: s                                        
       real(kind=rp), intent(in) :: x                                           
       real(kind=rp), intent(in) :: y                                           
       real(kind=rp), intent(in) :: z                                           
       real(kind=rp), intent(in) :: nx                                          
       real(kind=rp), intent(in) :: ny                                          
       real(kind=rp), intent(in) :: nz                                          
       integer, intent(in) :: ix                                                
       integer, intent(in) :: iy                                                
       integer, intent(in) :: iz                                                
       integer, intent(in) :: ie                                                
       real(kind=rp), intent(in) :: t                                           
       integer, intent(in) :: tstep                                             

       s = 0.0_rp
   end subroutine adjoint_force_user_bc_w
end module adjoint_case

