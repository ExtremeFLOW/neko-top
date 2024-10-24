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
  use adjoint_scheme, only: adjoint_scheme_t
  use adjoint_fctry, only: adjoint_scheme_factory
  use adjoint_pnpn, only: adjoint_pnpn_t
  use adjoint_output, only: adjoint_output_t
  use adjoint_ic, only: set_adjoint_ic
  use output_controller, only: output_controller_t
  use file, only: file_t
  use json_module, only: json_file
  use json_utils, only: json_get, json_get_or_default
  use json_utils_ext, only: json_key_fallback, json_get_subdict
  implicit none
  private
  public :: adjoint_case_t, adjoint_init, adjoint_free

  !> Adjoint case type.
  !! Todo: This should Ideally be a subclass of case_t, however, that is not yet
  !! suppoerted by Neko.
  type :: adjoint_case_t

     class(adjoint_scheme_t), allocatable :: scheme
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
    if (allocated(neko_case%scalar)) then
       this%have_scalar = .true.
    end if

    call adjoint_case_init_common(this, neko_case)

  end subroutine adjoint_init_from_attributes

  !> Initialize a neko_case from its (loaded) params object
  subroutine adjoint_case_init_common(this, neko_case)
    class(adjoint_case_t), intent(inout) :: this
    type(case_t), intent(inout) :: neko_case
    integer :: lx = 0
    logical :: scalar = .false.
    real(kind=rp) :: real_val
    character(len=:), allocatable :: string_val
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
    ! @todo Scalar adjoint is not implemented yet
    ! if (neko_case%params%valid_path('case.scalar')) then
    !    call json_get_or_default(neko_case%params, 'case.scalar.enabled', &
    ! scalar,                             .true.)
    ! end if

    ! if (scalar) then
    !    allocate(neko_case%scalar)
    !    call neko_case%scalar%init(neko_case%msh, this%scheme%c_Xh, &
    !         this%scheme%gs_Xh, neko_case%params, neko_case%usr,&
    !         neko_case%material_properties)
    !    call this%scheme%chkp%add_scalar(neko_case%scalar%output_controller)
    !    this%scheme%chkp%abs1 => neko_case%scalar%abx1
    !    this%scheme%chkp%abs2 => neko_case%scalar%abx2
    !    this%scheme%chkp%slag => neko_case%scalar%slag
    ! end if

    !
    ! Setup user defined conditions
    !
    ! if (neko_case%params%valid_path('case.fluid.inflow_condition')) then
    !    call json_get(neko_case%params, 'case.fluid.inflow_condition.type', &
    !         string_val)
    !    if (trim(string_val) .eq. 'user') then
    !       call neko_case%fluid%set_usr_inflow(neko_case%usr%fluid_user_if)
    !    end if
    ! end if

    ! Setup user boundary conditions for the scalar.
    ! if (scalar) then
    !    call neko_case%scalar%set_user_bc(neko_case%usr%scalar_user_bc)
    ! end if

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

    ! if (scalar) then
    !    call json_get(neko_case%params, 'case.scalar.initial_condition.type', &
    ! string_val)
    !    if (trim(string_val) .ne. 'user') then
    !       call set_scalar_ic(neko_case%scalar%output_controller, &
    !         neko_case%scalar%c_Xh, neko_case%scalar%gs_Xh, string_val, &
    ! neko_case%params)
    !    else
    !       call set_scalar_ic(neko_case%scalar%output_controller, &
    !         neko_case%scalar%c_Xh, neko_case%scalar%gs_Xh, &
    ! neko_case%usr%scalar_user_ic, neko_case%params)
    !    end if
    ! end if

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

    ! if (scalar) then
    !    call neko_case%scalar%slag%set(neko_case%scalar%output_controller)
    !    call neko_case%scalar%validate
    ! end if

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
       this%f_out = adjoint_output_t(precision, this%scheme, neko_case%scalar, &
            path = trim(neko_case%output_directory))
    else
       this%f_out = adjoint_output_t(precision, this%scheme, &
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

end module adjoint_case

