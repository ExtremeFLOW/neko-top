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
  use flow_ic, only: set_flow_ic
  use output_controller, only: output_controller_t
  use file, only: file_t
  use json_module, only: json_file
  use json_utils, only: json_get, json_get_or_default, json_extract_object
  use json_utils_ext, only: json_key_fallback
  use comm, only: pe_rank
  use, intrinsic :: iso_fortran_env, only: output_unit
  implicit none
  private
  public :: adjoint_case_t, adjoint_init, adjoint_free

  !> Adjoint case type.
  !! Todo: This should Ideally be a subclass of case_t, however, that is not yet
  !! suppoerted by Neko.
  type :: adjoint_case_t

     class(adjoint_fluid_scheme_t), allocatable :: scheme
     type(case_t), pointer :: case

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

    ! Read the tolerance
    call json_get_or_default(neko_case%params, "tol", tol, 1.0e-6_rp)

    call adjoint_init_from_attributes(this, neko_case, tol)

  end subroutine adjoint_init_from_json

  ! Constructor from attributes
  subroutine adjoint_init_from_attributes(this, neko_case, tol)
    class(adjoint_case_t), intent(inout) :: this
    class(case_t), intent(inout), target :: neko_case
    real(kind=rp), intent(in) :: tol

    this%case => neko_case
    this%tol = tol

    ! Check if the scalar field is allocated
    if (allocated(neko_case%scalar)) then
       this%have_scalar = .true.
    end if

    call adjoint_case_init_common(this)

  end subroutine adjoint_init_from_attributes

  !> Initialize a neko_case from its (loaded) params object
  subroutine adjoint_case_init_common(this)
    class(adjoint_case_t), intent(inout) :: this
    integer :: lx = 0
    logical :: scalar = .false.
    real(kind=rp) :: real_val = 0.0_rp
    character(len=:), allocatable :: string_val
    integer :: precision

    ! extra things for json
    type(json_file) :: ic_json
    character(len=:), allocatable :: json_key



    !
    ! Setup fluid scheme
    !
    call json_get(this%case%params, 'case.fluid.scheme', string_val)
    call adjoint_fluid_scheme_factory(this%scheme, trim(string_val))





    call json_get(this%case%params, 'case.numerics.polynomial_order', lx)
    lx = lx + 1 ! add 1 to get number of gll points
    call this%scheme%init(this%case%msh, lx, this%case%params, this%case%usr, &
         this%case%fluid%ext_bdf)


    !
    ! Setup scalar scheme
    !
    ! @todo Scalar adjoint is not implemented yet
    ! if (this%case%params%valid_path('case.scalar')) then
    !    call json_get_or_default(this%case%params, 'case.scalar.enabled', &
    ! scalar,                             .true.)
    ! end if

    ! if (scalar) then
    !    allocate(this%case%scalar)
    !    call this%case%scalar%init(this%case%msh, this%scheme%c_Xh, &
    !         this%scheme%gs_Xh, this%case%params, this%case%usr,&
    !         this%case%material_properties)
    !    call this%scheme%chkp%add_scalar(this%case%scalar%output_controller)
    !    this%scheme%chkp%abs1 => this%case%scalar%abx1
    !    this%scheme%chkp%abs2 => this%case%scalar%abx2
    !    this%scheme%chkp%slag => this%case%scalar%slag
    ! end if

    !
    ! Setup user defined conditions
    !
    ! if (this%case%params%valid_path('case.fluid.inflow_condition')) then
    !    call json_get(this%case%params, 'case.fluid.inflow_condition.type', &
    !         string_val)
    !    if (trim(string_val) .eq. 'user') then
    !       call this%case%fluid%set_usr_inflow(this%case%usr%fluid_user_if)
    !    end if
    ! end if

    ! Setup user boundary conditions for the scalar.
    ! if (scalar) then
    !    call this%case%scalar%set_user_bc(this%case%usr%scalar_user_bc)
    ! end if

    !
    ! Setup initial conditions
    !
    json_key = 'case.adjoint_fluid.initial_condition'
    ! json_key_fallback(this%case%params, &
    !      , 'case.fluid.initial_condition')

    if (pe_rank .eq. 0) then



       call this%case%params%print()



    end if

    call json_extract_object(this%case%params, json_key, ic_json)
    call json_get(ic_json, 'type', string_val)

    if (trim(string_val) .ne. 'user') then
       call set_flow_ic( &
            this%scheme%u_adj, this%scheme%v_adj, this%scheme%w_adj, &
            this%scheme%p_adj, this%scheme%c_Xh, this%scheme%gs_Xh, &
            string_val, ic_json)
    else
       call set_flow_ic( &
            this%scheme%u_adj, this%scheme%v_adj, this%scheme%w_adj, &
            this%scheme%p_adj, this%scheme%c_Xh, this%scheme%gs_Xh, &
            this%case%usr%fluid_user_ic, this%case%params)
    end if

    ! if (scalar) then
    !    call json_get(this%case%params, 'case.scalar.initial_condition.type', &
    ! string_val)

    !    if (trim(string_val) .ne. 'user') then
    !       call set_scalar_ic(this%case%scalar%output_controller, &
    !         this%case%scalar%c_Xh, this%case%scalar%gs_Xh, string_val, &
    ! this%case%params)
    !    else
    !       call set_scalar_ic(this%case%scalar%output_controller, &
    !         this%case%scalar%c_Xh, this%case%scalar%gs_Xh, &
    ! this%case%usr%scalar_user_ic, this%case%params)
    !    end if
    ! end if

    ! Add initial conditions to BDF scheme (if present)
    select type (f => this%scheme)
    type is (adjoint_fluid_pnpn_t)
       call f%ulag%set(f%u_adj)
       call f%vlag%set(f%v_adj)
       call f%wlag%set(f%w_adj)
    end select

    !
    ! Validate that the this%case is properly setup for time-stepping
    !
    call this%scheme%validate

    ! if (scalar) then
    !    call this%case%scalar%slag%set(this%case%scalar%output_controller)
    !    call this%case%scalar%validate
    ! end if

    !
    ! Setup output precision of the field files
    !
    call json_get_or_default(this%case%params, 'case.output_precision', &
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
    if (scalar) then
       this%f_out = adjoint_output_t(precision, this%scheme, this%case%scalar, &
            path = trim(this%case%output_directory))
    else
       this%f_out = adjoint_output_t(precision, this%scheme, &
            path = trim(this%case%output_directory))
    end if

    call json_get_or_default(this%case%params, 'case.fluid.output_control',&
         string_val, 'org')

    if (trim(string_val) .eq. 'org') then
       ! yes, it should be real_val below for type compatibility
       call json_get(this%case%params, 'case.nsamples', real_val)
       call this%output_controller%add(this%f_out, real_val, 'nsamples')
    else if (trim(string_val) .eq. 'never') then
       ! Fix a dummy 0.0 output_value
       call json_get_or_default(this%case%params, 'case.fluid.output_value', &
            real_val, 0.0_rp)
       call this%output_controller%add(this%f_out, 0.0_rp, string_val)
    else
       call json_get(this%case%params, 'case.fluid.output_value', real_val)
       call this%output_controller%add(this%f_out, real_val, string_val)
    end if

    ! !
    ! ! Save checkpoints (if nothing specified, default to saving at end of sim)
    ! !
    ! call json_get_or_default(this%case%params, 'case.output_checkpoints',&
    !      logical_val, .true.)
    ! if (logical_val) then
    !    call json_get_or_default(this%case%params, 'case.checkpoint_format', &
    !         string_val, "chkp")
    !   this%case%f_chkp = chkp_output_t(this%scheme%chkp, &
    ! path = output_directory, &
    !         ! fmt = trim(string_val))
    !    call json_get_or_default(this%case%params, 'case.checkpoint_control', &
    !         string_val, "simulationtime")
    !    call json_get_or_default(this%case%params, 'case.checkpoint_value', &
    ! real_val,&
    !         1e10_rp)
    !   call this%output_controller%add(this%case%f_chkp, real_val, string_val)
    ! end if

  end subroutine adjoint_case_init_common

  ! Destructor.
  subroutine adjoint_free(this)
    class(adjoint_case_t), intent(inout) :: this

    nullify(this%case)
    call this%scheme%free()
    call this%output_controller%free()

  end subroutine adjoint_free

end module adjoint_case

