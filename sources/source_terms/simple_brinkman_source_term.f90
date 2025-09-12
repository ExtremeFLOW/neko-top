! Copyright (c) 2023, The Neko Authors
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
!> Implements the `simple_brinkman_source_term_t` type.
! a term in the form $\chi \mathbf{u}$
module simple_brinkman_source_term
  use num_types, only: rp
  use field_list, only: field_list_t
  use json_module, only: json_file
  use json_utils, only: json_get, json_get_or_default
  use source_term, only: source_term_t
  use coefs, only: coef_t
  use neko_config, only: NEKO_BCKND_DEVICE
  use time_state, only: time_state_t
  use utils, only: neko_error
  use field, only: field_t
  use field_math, only: field_subcol3
  use interpolation, only: interpolator_t
  use space, only: space_t, GL
  use math, only: col2, invcol2, sub2
  use vector_math, only: vector_col3
  use scratch_registry, only: neko_scratch_registry
  use vector, only: vector_t
  implicit none
  private

  !> A simple Brinkman source term.
  ! We have a source term of the form $\chi \mathbf{u}$
  type, public, extends(source_term_t) :: simple_brinkman_source_term_t
     !> the fields corresponding to \f$\chi\f$
     type(field_t), pointer :: chi => null()
     !> the fields corresponding to u
     type(field_t), pointer :: u => null()
     !> the fields corresponding to v
     type(field_t), pointer :: v => null()
     !> the fields corresponding to w
     type(field_t), pointer :: w => null()
     ! --- for over-integration
     !> The original space used in the simulation
     type(space_t), pointer :: Xh_GLL
     !> The additional higher-order space used in dealiasing
     type(space_t), pointer :: Xh_GL
     !> cfs of the higher-order space
     type(coef_t), pointer :: c_Xh_GL
     !> Interpolator between the original and higher-order spaces
     type(interpolator_t), pointer :: GLL_to_GL
     !> if dealiasing should be applied
     logical :: dealias
     !> work arrays
     type(vector_t) :: accumulate, fld_GL, chi_GL

   contains
     !> The common constructor using a JSON object.
     procedure, pass(this) :: init => &
          simple_brinkman_source_term_init_from_json
     !> The constructor from type components.
     procedure, pass(this) :: init_from_components => &
          simple_brinkman_source_term_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => simple_brinkman_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute_ => simple_brinkman_source_term_compute
  end type simple_brinkman_source_term_t

contains
  !> The common constructor using a JSON object.
  !! @param this The source term.
  !! @param json The JSON object for the source.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  !! @param variable_name The name of the variable where the source term acts.
  subroutine simple_brinkman_source_term_init_from_json(this, json, fields, &
       coef, variable_name)
    class(simple_brinkman_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef
    character(len=*), intent(in) :: variable_name


    ! we shouldn't be initializing this from JSON
    ! maybe throw an error?


  end subroutine simple_brinkman_source_term_init_from_json

  !> The constructor from type components.
  !! @param this The source term.
  !! @param f_x, f_y, f_z the RHS of the equation (either primal or adjoint).
  !! @param chi the brinkman amplitude field.
  !! @param u, v, w the velocity field (either primal or adjoint).
  !! @param coef The SEM coeffs.
  !! @param c_Xh_GL The SEM coeffs on the over integration mesh.
  !! @param GLL_to_GL Interpolator between GLL and GL.
  !! @param dealias if dealiasing should be applied.
  subroutine simple_brinkman_source_term_init_from_components(this, &
       f_x, f_y, f_z, chi, u, v, w, coef, c_Xh_GL, GLL_to_GL, dealias)
    class(simple_brinkman_source_term_t), intent(inout) :: this
    type(field_t), pointer, intent(in) :: f_x, f_y, f_z
    type(field_list_t) :: fields
    type(coef_t), intent(in) :: coef
    type(coef_t), intent(in), target :: c_Xh_GL
    type(interpolator_t), intent(in), target :: GLL_to_GL
    logical, intent(in) :: dealias
    real(kind=rp) :: start_time
    real(kind=rp) :: end_time
    type(field_t), intent(in), target :: u, v, w
    type(field_t), intent(in), target :: chi
    integer :: nel, n_GL

    ! I wish you didn't need a start time and end time...
    ! but I'm just going to set a super big number...
    start_time = 0.0_rp
    end_time = 100000000.0_rp

    call this%free()

    ! this is copying the fluid source term init
    ! We package the fields for the source term to operate on in a field list.
    call fields%init(3)
    call fields%assign(1, f_x)
    call fields%assign(2, f_y)
    call fields%assign(3, f_z)

    call this%init_base(fields, coef, start_time, end_time)

    ! point everything in the correct places
    this%u => u
    this%v => v
    this%w => w
    ! and get chi out of the design
    this%chi => chi
    ! for over integration
    this%dealias = dealias
    this%c_Xh_GL => c_Xh_GL
    this%Xh_GL => this%c_Xh_GL%Xh
    this%Xh_GLL => this%coef%Xh
    this%GLL_to_GL => GLL_to_GL

    if (this%dealias) then
       ! allocate work arrays for dealiasing
       nel = this%coef%msh%nelv
       n_GL = nel * this%Xh_GL%lxyz
       call this%accumulate%init(n_GL)
       call this%fld_GL%init(n_GL)
       call this%chi_GL%init(n_GL)
    end if

  end subroutine simple_brinkman_source_term_init_from_components

  !> Destructor.
  subroutine simple_brinkman_source_term_free(this)
    class(simple_brinkman_source_term_t), intent(inout) :: this

    call this%free_base()
    call this%accumulate%free()
    call this%fld_GL%free()
    call this%chi_GL%free()

  end subroutine simple_brinkman_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param this The source term.
  !! @param time The time state.
  subroutine simple_brinkman_source_term_compute(this, time)
    class(simple_brinkman_source_term_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: fu, fv, fw
    type(field_t), pointer :: work
    integer :: temp_indices(1)
    integer :: n_GL, nel

    fu => this%fields%get_by_index(1)
    fv => this%fields%get_by_index(2)
    fw => this%fields%get_by_index(3)

    call neko_scratch_registry%request_field(work, temp_indices(1))

    if (this%dealias) then
       nel = this%coef%msh%nelv
       n_GL = nel * this%Xh_GL%lxyz
       call this%GLL_to_GL%map(this%chi_GL%x, this%chi%x, nel, this%Xh_GL)

       ! u
       call this%GLL_to_GL%map(this%fld_GL%x, this%u%x, nel, this%Xh_GL)
       call vector_col3(this%accumulate, this%chi_GL, this%fld_GL)
       ! Evaluate term on GL and preempt the GLL premultiplication
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_col2(this%accumulate%x_d, this%c_Xh_GL%B_d, n_GL)
          call this%GLL_to_GL%map(work%x, this%accumulate%x, nel, this%Xh_GLL)
          call device_invcol2(work%x_d, this%coef%B_d, work%size())
       else
          call col2(this%accumulate%x, this%c_Xh_GL%B, n_GL)
          call this%GLL_to_GL%map(work%x, this%accumulate%x, nel, this%Xh_GLL)
          call invcol2(work%x, this%coef%B, work%size())
       end if
       call sub2(fu%x, work%x, work%size())

       ! v
       call this%GLL_to_GL%map(this%fld_GL%x, this%v%x, nel, this%Xh_GL)
       call vector_col3(this%accumulate, this%chi_GL, this%fld_GL)
       ! Evaluate term on GL and preempt the GLL premultiplication
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_col2(this%accumulate%x_d, this%c_Xh_GL%B_d, n_GL)
          call this%GLL_to_GL%map(work%x, this%accumulate%x, nel, this%Xh_GLL)
          call device_invcol2(work%x_d, this%coef%B_d, work%size())
       else
          call col2(this%accumulate%x, this%c_Xh_GL%B, n_GL)
          call this%GLL_to_GL%map(work%x, this%accumulate%x, nel, this%Xh_GLL)
          call invcol2(work%x, this%coef%B, work%size())
       end if
       call sub2(fv%x, work%x, work%size())

       ! w
       call this%GLL_to_GL%map(this%fld_GL%x, this%w%x, nel, this%Xh_GL)
       call vector_col3(this%accumulate, this%chi_GL, this%fld_GL)
       ! Evaluate term on GL and preempt the GLL premultiplication
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_col2(this%accumulate%x_d, this%c_Xh_GL%B_d, n_GL)
          call this%GLL_to_GL%map(work%x, this%accumulate%x, nel, this%Xh_GLL)
          call device_invcol2(work%x_d, this%coef%B_d, work%size())
       else
          call col2(this%accumulate%x, this%c_Xh_GL%B, n_GL)
          call this%GLL_to_GL%map(work%x, this%accumulate%x, nel, this%Xh_GLL)
          call invcol2(work%x, this%coef%B, work%size())
       end if
       call sub2(fv%x, work%x, work%size())
    else

       call field_subcol3(fu, this%u, this%chi)
       call field_subcol3(fv, this%v, this%chi)
       call field_subcol3(fw, this%w, this%chi)

    end if

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine simple_brinkman_source_term_compute

end module simple_brinkman_source_term
