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
!> Implements the `adjoint_minimum_dissipation_source_term_t` type.
!
!
! If the objective function \int |\nabla u|^2,
! the corresponding adjoint forcing is \int \nabla v \cdot \nabla u
module adjoint_minimum_dissipation_source_term
  use num_types, only : rp
  use field_list, only : field_list_t
  use json_module, only : json_file
  use json_utils, only: json_get, json_get_or_default
  use source_term, only : source_term_t
  use coefs, only : coef_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use utils, only : neko_error
  use field, only: field_t
  use new_design, only: new_design_t
  use field_math, only: field_subcol3
  use user_intf, only: user_t, simulation_component_user_settings
  use json_module, only: json_file
  use steady_simcomp, only: steady_simcomp_t
  use simcomp_executor, only: neko_simcomps
  use fluid_user_source_term, only: fluid_user_source_term_t
  use num_types, only : rp
  use field, only : field_t
  use field_registry, only : neko_field_registry
  use math, only : rzero, copy, chsign
  use device_math, only: device_copy, device_cmult
  use neko_config, only: NEKO_BCKND_DEVICE
  use operators, only: curl
  use scratch_registry, only : neko_scratch_registry
  implicit none
  private

  !> A constant source term.
  !! The strength is specified with the `values` keyword, which should be an
  !! array, with a value for each component of the source.
  type, public, extends(source_term_t) :: adjoint_minimum_dissipation_source_term_t
   
   ! again, for a mask this is silly... we can fix this later in the week
   type(field_t), pointer :: u,v,w, mask

   contains
     !> The common constructor using a JSON object.
     procedure, pass(this) :: init => adjoint_minimum_dissipation_source_term_init_from_json
     !> The constructor from type components.
     procedure, pass(this) :: init_from_components => &
       adjoint_minimum_dissipation_source_term_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => adjoint_minimum_dissipation_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute_ => adjoint_minimum_dissipation_source_term_compute
  end type adjoint_minimum_dissipation_source_term_t

contains
  !> The common constructor using a JSON object.
  !! @param json The JSON object for the source.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  subroutine adjoint_minimum_dissipation_source_term_init_from_json(this, json, fields, coef)
    class(adjoint_minimum_dissipation_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(inout), target :: fields
    type(coef_t), intent(inout), target :: coef
    real(kind=rp), allocatable :: values(:)
    real(kind=rp) :: start_time, end_time


    ! we shouldn't be initializing this from JSON
    ! maybe throw an error?


  end subroutine adjoint_minimum_dissipation_source_term_init_from_json

  !> The constructor from type components.
  ! NOTE!
  ! u,v,w reffer to the primal, not the adjoint
  subroutine adjoint_minimum_dissipation_source_term_init_from_components(this, f_x, f_y, f_z, &	
                                                    u, v, w, coef)
    class(adjoint_minimum_dissipation_source_term_t), intent(inout) :: this
    type(field_t), pointer, intent(in) :: f_x, f_y, f_z
    type(field_list_t) :: fields
    type(coef_t) :: coef
    real(kind=rp) :: start_time
    real(kind=rp) :: end_time
    type(field_t), intent(in), target :: u, v, w 
    ! TODo
    ! do masks later
    !type(field_t), intent(in), target :: mask

    ! I wish you didn't need a start time and end time...
    ! Tim you're going to hate this... but I'm just going to set a super big number...
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

	 ! Real stuff

    ! point everything in the correct places
    ! NOTE!!!
    ! this is the primal!
    this%u => u
    this%v => v
    this%w => w

    ! TODO
    !this%mask => mask





  end subroutine adjoint_minimum_dissipation_source_term_init_from_components

  !> Destructor.
  subroutine adjoint_minimum_dissipation_source_term_free(this)
    class(adjoint_minimum_dissipation_source_term_t), intent(inout) :: this

    call this%free_base()
  end subroutine adjoint_minimum_dissipation_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine adjoint_minimum_dissipation_source_term_compute(this, t, tstep)
    class(adjoint_minimum_dissipation_source_term_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    integer :: i
    type(field_t), pointer :: u, v, w
    type(field_t), pointer :: fu, fv, fw
    !type(field_t), pointer :: dudx, dudy, dudz
    !type(field_t), pointer :: dvdx, dvdy, dvdz
    !type(field_t), pointer :: dwdx, dwdy, dwdz
    type(field_t), pointer :: wo1, wo2, wo3, wo4, wo5, wo6
    type(field_t), pointer :: t1 , t2
    integer :: temp_indices(8)
    integer n
    print *, "called forcing"


    fu => this%fields%get_by_index(1)
    fv => this%fields%get_by_index(2)
    fw => this%fields%get_by_index(3)

    n = fu%size()

    ! fuck I'm not sure about this... I need a pen and paper
    ! also there should be a way to pre-process this forcing term...
    ! instead of recalculating it every time
    u => this%u
    v => this%v
    w => this%w
    call neko_scratch_registry%request_field(wo1, temp_indices(1))
    call neko_scratch_registry%request_field(wo2, temp_indices(2))
    call neko_scratch_registry%request_field(wo3, temp_indices(3))
    call neko_scratch_registry%request_field(wo4, temp_indices(4))
    call neko_scratch_registry%request_field(wo5, temp_indices(5))
    call neko_scratch_registry%request_field(wo6, temp_indices(6))
    call neko_scratch_registry%request_field(t1 , temp_indices(7))
    call neko_scratch_registry%request_field(t2 , temp_indices(8))

    ! ok we're computing gradients at every timestep... which is stupid...
    ! BUT
    ! if this was unsteady we would have to do this.

    ! this is cheating a little bit...
    ! in strong form, \nabla u . \nabla v =>  v . \nabla^2 u + bdry
    !
    ! we can do this properly later in weak form, ideally using ax_helm or so
    !
    ! for now we'll work in strong form and ignore the bdry
    ! and suffer the double derivative :/
    !
    ! in fact, we'll go even quicker and use
    ! \nabla ^2 u = grad (div (u)) - curl ( curl (u )) and assume divergence free u

    call curl(wo1, wo2, wo3, u, v, w, t1, t2, this%coef)
    call curl(wo4, wo5, wo6, wo1, wo2, wo3, t1, t2, this%coef)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(fu%x_d, wo4%x_d, n)
       call device_copy(fv%x_d, wo5%x_d, n)
       call device_copy(fw%x_d, wo6%x_d, n)
       call device_cmult(fu%x_d, -1.0_rp, n)
       call device_cmult(fv%x_d, -1.0_rp, n)
       call device_cmult(fw%x_d, -1.0_rp, n)
       !TODO
       ! we would also consider a mask here
    else
       call copy(fu%x, wo4%x, n)
       call copy(fv%x, wo5%x, n)
       call copy(fw%x, wo6%x, n)
       call chsign(fu%x, n)
       call chsign(fv%x, n)
       call chsign(fw%x, n)
       !TODO
       ! we would also consider a mask here
    end if
    ! don't worry... we'll write this MUCH cleaner in the final version

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine adjoint_minimum_dissipation_source_term_compute

end module adjoint_minimum_dissipation_source_term
