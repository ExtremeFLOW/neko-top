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
!> Implements the `volume_constraint_t` type.
! $V = \frac{1}{|\Omega_O|}\int_{\Omega_O} \tilde{\rho} d\Omega$
! Either
! $V < V_\text{max}$
! $V > V_\text{min}$
module volume_constraint
  use num_types, only : rp
  use field_list, only : field_list_t
  use json_utils, only: json_get, json_get_or_default
  use source_term, only : source_term_t
  use coefs, only : coef_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use utils, only : neko_error
  use field, only: field_t
  use field_math, only: field_col3, field_addcol3
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
  use operators, only: curl, grad
  use scratch_registry, only : neko_scratch_registry
  use objective_function, only : objective_function_t
  use fluid_scheme, only : fluid_scheme_t
  use adjoint_scheme, only : adjoint_scheme_t
  use fluid_source_term, only: fluid_source_term_t
  use math, only : glsc2
  use field_math, only: field_rone, field_cmult
  use topopt_design, only: topopt_design_t
  implicit none
  private

  !> A constraint on the volume of the design.
  type, public, extends(objective_function_t) :: volume_constraint_t

     !> whether it is minimum or maximum volume
     ! is_max = .false., 	ie V > V_min  		=>		 -V + V_max < 0
     ! is_max = .true. , 	ie V < V_max  		=>		  V - V_max < 0
     ! TODO
     ! this can be done smarter with parameters
     logical :: is_max
     !> Maximum (or minimum) volume
     ! maximum volume prescribed
     real(kind=rp) :: v_max
     !> current volume
     real(kind=rp) :: volume


   contains
     !> The common constructor using a JSON object.
     procedure, pass(this) :: init => volume_constraint_init
     !> Destructor.
     procedure, pass(this) :: free => volume_constraint_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute => volume_constraint_compute
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute_sensitivity => &
          volume_constraint_compute_sensitivity
  end type volume_constraint_t

contains
  !> The common constructor using a JSON object.
  !! @param design the design
  !! @param fluid the fluid scheme
  !! @param adjoint the adjoint scheme
  subroutine volume_constraint_init(this, design, fluid, adjoint)
    class(volume_constraint_t), intent(inout) :: this
    class(fluid_scheme_t), intent(inout) :: fluid
    class(adjoint_scheme_t), intent(inout) :: adjoint
    type(topopt_design_t), intent(inout) :: design

    ! TODO
    ! I don't think there's much to init here
    ! maybe we should include a fmax in here,
    ! ie, if we want f_i(x) < f_i_max
    ! with the MMA notation, this is
    ! f_i(x) - a_i*z - y_i <= f_i_max,

    ! f_i(x) - a_i*z - y_i - f_i_max <= 0,

    ! so we compute
    ! f_i(x) - f_i_max

    ! when we compute the value of the constraint
    !
    ! But I actually think it's better to include the mins and max's
    ! in MMA it's self.
    !
    ! anyway...
    ! here we hard code for now
    this%is_max = .false.
    this%v_max = 0.2

    call this%init_base(fluid%dm_Xh)

  end subroutine volume_constraint_init


  !> Destructor.
  subroutine volume_constraint_free(this)
    class(volume_constraint_t), intent(inout) :: this

    call this%free_base()
  end subroutine volume_constraint_free

  !> The computation of the constraint.
  !! @param design the design
  !! @param fluid the fluid scheme
  subroutine volume_constraint_compute(this, design, fluid)
    class(volume_constraint_t), intent(inout) :: this
    class(fluid_scheme_t), intent(in) :: fluid
    type(topopt_design_t), intent(inout) :: design
    integer :: i
    type(field_t), pointer :: wo1, wo2, wo3
    type(field_t), pointer :: objective_field
    integer :: temp_indices(4)
    integer n

    ! Again, we don't really need to take design and fluid in here...
    n = design%design_indicator%size()
    ! TODO
    ! in the future we should be using the mapped design varaible
    !corresponding to this constraint!!!
    this%volume = glsc2(design%design_indicator%x, fluid%c_xh%B, n)

    ! NOTE
    ! TODO
    ! the definition of the "volume" can get a little tricky once we start
    ! introducing masks for the optimization domain, because then we should
    ! calculate the volume of the optimization domain and take a ratio that way.
    ! point is, a design should have an internal parameter of the volume of
    ! the design domain,
    ! and we should us THAT volume for computing the volume percentage
    this%volume = this%volume/fluid%c_xh%volume

    ! then we need to check min or max
    if (this%is_max) then
       ! max volume
       this%objective_function_value = this%volume - this%v_max
    else
       ! min volume
       this%objective_function_value = -this%volume + this%v_max
    end if


    ! TODo
    ! GPU


  end subroutine volume_constraint_compute

  !> The computation of the sensitivity.
  !! @param design the design
  !! @param fluid the fluid scheme
  !! @param adjoint the adjoint scheme
  subroutine volume_constraint_compute_sensitivity(this, design, fluid, adjoint)
    class(volume_constraint_t), intent(inout) :: this
    type(topopt_design_t), intent(inout) :: design
    class(fluid_scheme_t), intent(in) :: fluid
    class(adjoint_scheme_t), intent(in) :: adjoint

    call field_rone(this%sensitivity_to_coefficient)

    if (this%is_max) then
       ! max volume
       call field_cmult(this%sensitivity_to_coefficient, &
            1.0_rp / fluid%c_xh%volume)
    else
       ! min volume
       call field_cmult(this%sensitivity_to_coefficient, &
            -1.0_rp / fluid%c_xh%volume)
    end if

  end subroutine volume_constraint_compute_sensitivity

end module volume_constraint
