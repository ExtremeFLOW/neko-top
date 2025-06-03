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
  use num_types, only: rp, dp, sp
  use case, only: case_t
  use adjoint_fluid_scheme, only: adjoint_fluid_scheme_t
  use adjoint_fluid_fctry, only: adjoint_fluid_scheme_factory
  use adjoint_fluid_pnpn, only: adjoint_fluid_pnpn_t
  use adjoint_output, only: adjoint_output_t
  use scalar_ic, only: set_scalar_ic
  use flow_ic, only: set_flow_ic
  use output_controller, only: output_controller_t
  use file, only: file_t
  use json_module, only: json_file
  use json_utils, only: json_get, json_get_or_default, json_extract_object
  use adjoint_scalar_scheme, only: adjoint_scalar_scheme_t
  use adjoint_scalar_pnpn, only : adjoint_scalar_pnpn_t
  use logger, only : neko_log
  use time_state, only : time_state_t
  use utils, only: neko_error
  use adjoint_scalar_convection_source_term, only: &
       adjoint_scalar_convection_source_term_t
  use json_utils_ext, only: json_key_fallback
  implicit none
  private
  public :: adjoint_case_t, adjoint_init, adjoint_free

  !> Adjoint case type.
  !! Todo: This should Ideally be a subclass of case_t, however, that is not yet
  !! supported by Neko.
  type :: adjoint_case_t
     !> Adjoint fluid
     class(adjoint_fluid_scheme_t), allocatable :: fluid_adj
     !> Adjoint scalar
     type(adjoint_scalar_pnpn_t), allocatable :: scalar_adj
     !> Source term coupling the adjoint scalar to the adjoint fluid
     type(adjoint_scalar_convection_source_term_t), allocatable :: &
          adjoint_convection_term
     type(case_t), pointer :: case
     type(time_state_t) :: time

     ! Fields
     real(kind=rp) :: tol
     type(adjoint_output_t) :: f_out
     type(output_controller_t) :: output_controller

     logical :: have_scalar = .false.

  end type adjoint_case_t

  interface adjoint_init
     module procedure adjoint_init_from_json, adjoint_init_from_attributes
  end interface adjoint_init

contains

  ! Constructor from json.
  subroutine adjoint_init_from_json(this, neko_case)
    class(adjoint_case_t), intent(inout) :: this
    type(case_t), target, intent(inout) :: neko_case
    real(kind=rp) :: tol
    logical :: have_scalar

    ! Read the tolerance
    call json_get_or_default(neko_case%params, "tol", tol, 1.0e-6_rp)

    ! I think this is correct.
    ! Maybe there would be a case where we would want a scalar but
    ! no adjoint scalar. So this forces us to prescribe an adjoint scalar.
    call json_get_or_default(neko_case%params, &
         'case.adjoint_scalar.enabled', have_scalar, .false.)

    call adjoint_init_from_attributes(this, neko_case, tol, have_scalar)

  end subroutine adjoint_init_from_json

  ! Constructor from attributes
  subroutine adjoint_init_from_attributes(this, neko_case, tol, have_scalar)
    class(adjoint_case_t), intent(inout) :: this
    class(case_t), intent(inout), target :: neko_case
    real(kind=rp), intent(in) :: tol
    logical :: have_scalar

    this%case => neko_case
    this%tol = tol
    this%have_scalar = have_scalar


    call adjoint_case_init_common(this, neko_case)

  end subroutine adjoint_init_from_attributes

  !> Initialize a neko_case from its (loaded) params object
  subroutine adjoint_case_init_common(this, neko_case)
    class(adjoint_case_t), intent(inout) :: this
    type(case_t), intent(inout) :: neko_case
    integer :: lx = 0
    real(kind=rp) :: real_val = 0.0_rp
    character(len=:), allocatable :: string_val
    integer :: precision

    ! extra things for json
    type(json_file) :: ic_json, numerics_params
    character(len=:), allocatable :: json_key

    !
    ! Setup adjoint fluid
    !
    call json_get(neko_case%params, 'case.fluid.scheme', string_val)
    call adjoint_fluid_scheme_factory(this%fluid_adj, trim(string_val))

    call json_get(neko_case%params, 'case.numerics.polynomial_order', lx)
    lx = lx + 1 ! add 1 to get number of gll points

    select type (f => this%fluid_adj)
    type is (adjoint_fluid_pnpn_t)
       call f%init(neko_case%msh, lx, neko_case%params, &
            neko_case%user, neko_case%chkp)
    end select
    !
    ! Setup adjoint scalar
    !
    ! @todo no scalar_adj factory for now, probably not needed

    ! hmmm should we check for scalar or adjoint scalar?
    ! I'm going to check for adjoint scalar because maybe there would be
    ! a corner case where someone would want the scalar but not the
    ! adjoint scalar?


    if (this%have_scalar) then
       allocate(this%scalar_adj)
       call json_extract_object(neko_case%params, 'case.numerics', &
            numerics_params)

       ! @todo
       ! these tlag and dtlag are new, we likely need to update the standard
       ! fluid in a different PR.
       ! For now I'm commenting them out in the scalar.
       ! this%scalar_adj%chkp%tlag => neko_case%tlag
       ! this%scalar_adj%chkp%dtlag => neko_case%dtlag
       call this%scalar_adj%init(neko_case%msh, neko_case%fluid%c_Xh, &
            neko_case%fluid%gs_Xh, neko_case%params, numerics_params, &
            neko_case%user, neko_case%chkp, &
            neko_case%fluid%ulag, neko_case%fluid%vlag, &
            neko_case%fluid%wlag, neko_case%fluid%ext_bdf, neko_case%fluid%rho)


       ! call neko_case%fluid%chkp%add_scalar(this%scalar_adj%s_adj)

       ! ----------------------------------------------------------------------
       ! @todo
       ! I don't really understand checkpoints or why the fluid would need to
       ! know about the scalar's lag and time integration terms.
       !
       ! Since we won't be using checkpoints, I'm commenting this out, but
       ! leaving a rather large TODO here for when we come back to unsteady.
       ! ----------------------------------------------------------------------
       ! neko_case%fluid%chkp%abs1 => this%scalar_adj%abx1
       ! neko_case%fluid%chkp%abs2 => this%scalar_adj%abx2
       ! neko_case%fluid%chkp%slag => this%scalar_adj%slag

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

       ! allocate the coupling term
       allocate(this%adjoint_convection_term)
       ! initialize the coupling term
       call this%adjoint_convection_term%init_from_components( &
            this%fluid_adj%f_adj_x, this%fluid_adj%f_adj_y, &
            this%fluid_adj%f_adj_z, this%case%scalar%s, &
            this%scalar_adj%s_adj, this%fluid_adj%c_Xh)

       select type (f => this%fluid_adj)
       type is (adjoint_fluid_pnpn_t)
          ! append the coupling term to the adjoint velocity equation
          call f%source_term%add(this%adjoint_convection_term)
       end select
    end if

    !
    ! Time step
    !
    ! call json_get_or_default(neko_case%params, 'case.variable_timestep', &
    !      logical_val, .false.)
    ! if (.not. logical_val) then
    call json_get(neko_case%params, 'case.timestep', this%time%dt)
    ! else
    !    ! randomly set an initial dt to get cfl when dt is variable
    !    this%time%dt = 1.0_rp
    ! end if

    !
    ! End time
    !
    call json_get(neko_case%params, 'case.end_time', this%time%end_time)

    !
    ! Setup user defined conditions
    !
    ! if (neko_case%params%valid_path('case.fluid.inflow_condition')) then
    !    call json_get(neko_case%params, 'case.fluid.inflow_condition.type', &
    !         string_val)
    !    if (trim(string_val) .eq. 'user') then
    !       call neko_case%fluid%set_usr_inflow(neko_case%user%fluid_user_if)
    !    end if
    ! end if

    ! Setup user boundary conditions for the scalar.
    ! if (scalar_adj) then
    !    call neko_case%scalar_adj%set_user_bc(neko_case%user%scalar_user_bc)
    ! end if

    !
    ! Setup initial conditions
    !
    json_key = json_key_fallback(neko_case%params, &
         'case.adjoint_fluid.initial_condition', 'case.fluid.initial_condition')

    call json_extract_object(neko_case%params, json_key, ic_json)
    call json_get(ic_json, 'type', string_val)

    if (trim(string_val) .ne. 'user') then
       call set_flow_ic( &
            this%fluid_adj%u_adj, this%fluid_adj%v_adj, this%fluid_adj%w_adj, &
            this%fluid_adj%p_adj, this%fluid_adj%c_Xh, this%fluid_adj%gs_Xh, &
            string_val, ic_json)
    else
       call set_flow_ic( &
            this%fluid_adj%u_adj, this%fluid_adj%v_adj, this%fluid_adj%w_adj, &
            this%fluid_adj%p_adj, this%fluid_adj%c_Xh, this%fluid_adj%gs_Xh, &
            neko_case%user%fluid_user_ic, neko_case%params)
    end if

    call neko_log%end_section()

    if (this%have_scalar) then

       ! we shouldn't fallback to the primal here.
       call json_get(neko_case%params, &
            'case.adjoint_scalar.initial_condition.type', string_val)
       call json_extract_object(neko_case%params, &
            'case.adjoint_scalar.initial_condition', ic_json)

       !call neko_log%section("Adjoint scalar initial condition ")

       if (trim(string_val) .ne. 'user') then
          call set_scalar_ic(this%scalar_adj%s_adj, this%scalar_adj%c_Xh, &
               this%scalar_adj%gs_Xh, string_val, ic_json)
       else
          call neko_error("user defined ICs not implemented for adjoint scalar")
          ! call set_scalar_ic(this%scalar_adj%s_adj, &
          !      this%scalar_adj%c_Xh, this%scalar_adj%gs_Xh, &
          !      this%usr%scalar_user_ic, neko_case%params)
       end if

       ! call neko_log%end_section()

    end if

    ! Add initial conditions to BDF fluid_adj (if present)
    select type (f => this%fluid_adj)
    type is (adjoint_fluid_pnpn_t)
       call f%ulag%set(f%u_adj)
       call f%vlag%set(f%v_adj)
       call f%wlag%set(f%w_adj)
    end select

    !
    ! Validate that the neko_case is properly setup for time-stepping
    !
    call this%fluid_adj%validate

    if (this%have_scalar) then
       call this%scalar_adj%s_adj_lag%set(this%scalar_adj%s_adj)
       call this%scalar_adj%validate
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
    call this%output_controller%init(neko_case%time%end_time)
    if (this%have_scalar) then
       this%f_out = adjoint_output_t(precision, this%fluid_adj, &
            this%scalar_adj, path = trim(neko_case%output_directory))
    else
       this%f_out = adjoint_output_t(precision, this%fluid_adj, &
            path = trim(neko_case%output_directory))
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
    !   neko_case%f_chkp = chkp_output_t(this%fluid_adj%chkp, &
    ! path = output_directory, &
    !         ! fmt = trim(string_val))
    !    call json_get_or_default(neko_case%params, 'case.checkpoint_control', &
    !         string_val, "simulationtime")
    !    call json_get_or_default(neko_case%params, 'case.checkpoint_value', &
    ! real_val,&
    !         1e10_rp)
    !   call this%output_controller%add(neko_case%f_chkp, real_val, string_val)
    ! end if

    !
    ! Initialize time and step
    !
    this%time%t = 0d0
    this%time%tstep = 0

  end subroutine adjoint_case_init_common

  ! Destructor.
  subroutine adjoint_free(this)
    class(adjoint_case_t), intent(inout) :: this

    nullify(this%case)
    if (allocated(this%scalar_adj)) then
       call this%scalar_adj%free()
       deallocate(this%scalar_adj)
    end if

    if (allocated(this%fluid_adj)) then
       call this%fluid_adj%free()
       deallocate(this%fluid_adj)
    end if
    call this%output_controller%free()

  end subroutine adjoint_free

end module adjoint_case

