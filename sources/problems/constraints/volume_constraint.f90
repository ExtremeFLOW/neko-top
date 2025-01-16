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
  use constraint, only : constraint_t
  use simulation, only : simulation_t
  use fluid_scheme, only : fluid_scheme_t
  use adjoint_scheme, only : adjoint_scheme_t
  use fluid_source_term, only: fluid_source_term_t
  use math, only : glsc2
  use field_math, only: field_rone, field_cmult
  use design, only: design_t
  use topopt_design, only: topopt_design_t
  use mask_ops, only: mask_exterior_const
  use math_ext, only: glsc2_mask
  implicit none
  private

  !> A constraint on the volume of the design.
  type, public, extends(constraint_t) :: volume_constraint_t
     private

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
     !> volume of the optimization domain
     real(kind=rp) :: volume_domain

     class(fluid_scheme_t), pointer :: fluid

   contains
     !> The common constructor using a JSON object.
     procedure, public, pass(this) :: init => volume_constraint_init
     !> Destructor.
     procedure, public, pass(this) :: free => volume_constraint_free
     !> Computes the source term and adds the result to `fields`.
     procedure, public, pass(this) :: update_value => &
          volume_constraint_update_value
     !> Computes the source term and adds the result to `fields`.
     procedure, public, pass(this) :: update_sensitivity => &
          volume_constraint_update_sensitivity
  end type volume_constraint_t

contains
  !> The common constructor using a JSON object.
  !! @param design the design
  !! @param fluid the fluid scheme
  !! @param adjoint the adjoint scheme
  subroutine volume_constraint_init(this, design, simulation)
    class(volume_constraint_t), intent(inout) :: this
    type(simulation_t), target, intent(inout) :: simulation
    class(design_t), intent(in) :: design
    type(topopt_design_t), pointer :: topopt_design => null()
    type(field_t), pointer :: work
    integer :: temp_indices(1)
    integer :: n

    select type (design)
      type is (topopt_design_t)
       topopt_design => design
      class default
       call neko_error('Volume constraint only works with topopt_design')
    end select

    ! TODO
    ! I don't think there's much to init here
    ! maybe we should include a fmax in here,
    ! ie, if we want f_i(x) < f_i_max
    ! with the MMA notation, this is
    ! f_i(x) - a_i*z - y_i <= f_i_max,

    ! f_i(x) - a_i*z - y_i - f_i_max <= 0,

    ! so we update_value
    ! f_i(x) - f_i_max

    ! when we update_value the value of the constraint
    !
    ! But I actually think it's better to include the mins and max's
    ! in MMA it's self.
    !
    ! anyway...
    ! here we hard code for now
    this%is_max = .false.
    this%v_max = 0.2

    this%fluid => simulation%neko_case%fluid

    ! Now we can extract the mask/if_mask from the design
    n = topopt_design%design_indicator%size()
    if (topopt_design%if_mask) then
       ! init the base
       call this%init_base(topopt_design, topopt_design%if_mask, &
            topopt_design%optimization_domain%name)

       ! calculate the volume of the optimization domain
       call neko_scratch_registry%request_field(work , temp_indices(1))
       call field_rone(work)
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call neko_error('GPU not supported volume constraint')
       else
          this%volume_domain = glsc2_mask(work%x, this%fluid%c_xh%B, &
               n, this%mask%mask, this%mask%size)
       end if
       call neko_scratch_registry%relinquish_field(temp_indices)
    else
       ! init the base
       call this%init_base(topopt_design)

       ! point to the volume of the domain
       this%volume_domain = this%fluid%c_xh%volume
    end if

  end subroutine volume_constraint_init


  !> Destructor.
  subroutine volume_constraint_free(this)
    class(volume_constraint_t), intent(inout) :: this

    call this%free_base()
  end subroutine volume_constraint_free

  !> The computation of the constraint.
  !! @param design the design
  !! @param fluid the fluid scheme
  subroutine volume_constraint_update_value(this, design)
    class(volume_constraint_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(topopt_design_t), pointer :: topopt_design => null()
    integer n

    select type (design)
      type is (topopt_design_t)
       topopt_design => design
      class default
       call neko_error('Volume constraint only works with topopt_design')
    end select

    n = topopt_design%design_indicator%size()
    ! TODO
    ! in the future we should be using the mapped design varaible
    !corresponding to this constraint!!!
    if (topopt_design%if_mask) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call neko_error('GPU not supported volume constraint')
       else
          this%volume = glsc2_mask(topopt_design%design_indicator%x, this%fluid%c_xh%B, &
               n, this%mask%mask, this%mask%size)
       end if
    else
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call neko_error('GPU not supported volume constraint')
       else
          this%volume = glsc2(topopt_design%design_indicator%x, this%fluid%c_xh%B, n)
       end if
    end if

    this%volume = this%volume/this%volume_domain

    ! then we need to check min or max
    if (this%is_max) then
       ! max volume
       this%value = this%volume - this%v_max
    else
       ! min volume
       this%value = -this%volume + this%v_max
    end if


    ! TODo
    ! GPU


  end subroutine volume_constraint_update_value

  !> The computation of the sensitivity.
  !! @param design the design
  !! @param fluid the fluid scheme
  !! @param adjoint the adjoint scheme
  subroutine volume_constraint_update_sensitivity(this, design)
    class(volume_constraint_t), intent(inout) :: this
    class(design_t), intent(in) :: design

    this%sensitivity = 1.0_rp

    if (this%is_max) then
       ! max volume
       this%sensitivity = this%sensitivity * ( 1.0_rp / this%volume_domain)
    else
       ! min volume
       this%sensitivity = this%sensitivity * (-1.0_rp / this%volume_domain)
    end if


    ! TODO
    ! Abbas, don't just mask the sensitivity like I'm doing here, make sure
    ! the only design variables entering MMA are those within the mask.
    ! This way you get the correct N etc.

    ! Look into the `masked_red_copy` function that Martin implemented.
    ! That function will copy from one array to another, but the target
    ! only have the size of the mask, not the full size.
    if (this%if_mask) then
       call mask_exterior_const(&
            this%sensitivity, &
            this%mask, 0.0_rp)
    end if

  end subroutine volume_constraint_update_sensitivity

end module volume_constraint
