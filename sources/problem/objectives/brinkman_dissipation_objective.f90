!> @file brinkman_dissipation_objective.f90
!! @copyright
!! Copyright (c) 2025, The Neko-TOP Authors
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
!> Implements the `brinkman_dissipation_objective_t` type.
!
! F = \int |\nabla u|^2  + K \int \chi \u^2
!
!  |viscous dissipation | |"Brinkman dissipation"|
!
module brinkman_dissipation_objective
  use objective, only: objective_t
  use design, only: design_t
  use brinkman_design, only: brinkman_design_t
  use simulation_m, only: simulation_t
  use adjoint_brinkman_dissipation_source_term, only: &
       adjoint_brinkman_dissipation_source_term_t
  use adjoint_fluid_pnpn, only: adjoint_fluid_pnpn_t
  use num_types, only: rp
  use field, only: field_t
  use scratch_registry, only: neko_scratch_registry, scratch_registry_t
  use neko_config, only: NEKO_BCKND_DEVICE
  use mask_ops, only: mask_exterior_const, compute_masked_volume
  use utils, only: neko_error
  use json_module, only: json_file
  use json_utils, only: json_get_or_default
  use registry, only: neko_registry
  use interpolation, only: interpolator_t
  use space, only: space_t, GL
  use coefs, only: coef_t
  use math, only: glsc2, copy, col2, invcol2
  use device_math, only: device_copy, device_glsc2, device_col2, device_invcol2
  use math_ext, only: glsc2_mask
  use field_math, only: field_col3, field_addcol3, field_cmult, field_col2
  use continuation_scheduler, only: nekotop_continuation
  implicit none
  private

  !> An objective function corresponding to out of plane stresses
  !! \f$ F =  \int_Omega \frac{1}{2} \chi |\mathbf{u}|^2 d \Omega \f$
  type, public, extends(objective_t) :: brinkman_dissipation_objective_t
     private

     !> Pointer to the u field.
     type(field_t), pointer :: u => null()
     !> Pointer to the v field.
     type(field_t), pointer :: v => null()
     !> Pointer to the w field.
     type(field_t), pointer :: w => null()
     !> Pointer to the brinkman amplitude field.
     type(field_t), pointer :: brinkman_amplitude => null()
     !> dealias sensitivity
     logical :: dealias_sensitivity
     !> dealias forcing
     logical :: dealias_forcing
     !> volume of objective domain
     real(kind=rp) :: volume

     ! ---- everything GLL ----
     !> The coefs
     type(coef_t), pointer :: c_Xh_GLL
     !> The original space used in the simulation
     type(space_t), pointer :: Xh_GLL

     ! ---- everything for GL ----
     !> The additional higher-order space used in dealiasing
     type(space_t), pointer :: Xh_GL
     !> coefs of the higher-order space
     type(coef_t), pointer :: c_Xh_GL
     !> Interpolator between the original and higher-order spaces
     type(interpolator_t), pointer :: GLL_to_GL
     !> GL scratch registry
     type(scratch_registry_t), pointer :: scratch_GL
     !> Physical dimension
     integer :: gdim


   contains

     !> The common constructor using a JSON object.
     procedure, public, pass(this) :: init_json_sim => &
          brinkman_dissipation_init_json_sim
     !> The actual constructor.
     procedure, public, pass(this) :: init_from_attributes => &
          brinkman_dissipation_init_attributes
     !> Destructor.
     procedure, public, pass(this) :: free => brinkman_dissipation_free
     !> Computes the value of the objective function.
     procedure, public, pass(this) :: update_value => &
          brinkman_dissipation_update_value
     !> Computes the sensitivity with respect to the coefficient \f$\chi\f$.
     procedure, public, pass(this) :: update_sensitivity => &
          brinkman_dissipation_update_sensitivity

  end type brinkman_dissipation_objective_t

contains

  !> The common constructor using a JSON object.
  !! @param this The objective.
  !! @param json the JSON object.
  !! @param design the design.
  !! @param simulation the simulation.
  subroutine brinkman_dissipation_init_json_sim(this, json, design, simulation)
    class(brinkman_dissipation_objective_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(design_t), intent(in) :: design
    type(simulation_t), target, intent(inout) :: simulation

    character(len=:), allocatable :: mask_name
    character(len=:), allocatable :: name
    real(kind=rp) :: weight
    logical :: dealias_sensitivity, dealias_forcing
    real(kind=rp) :: start_time, end_time

    call nekotop_continuation%json_get_or_register(json, 'weight', &
         this%weight, weight, 1.0_rp)
    call json_get_or_default(json, "mask_name", mask_name, "")
    call json_get_or_default(json, "name", name, "Brinkman dissipation")
    call json_get_or_default(json, "dealias_sensitivity", &
         dealias_sensitivity, .true.)
    call json_get_or_default(json, "dealias_forcing", &
         dealias_forcing, .true.)
    call json_get_or_default(json, "start_time", start_time, 0.0_rp)
    call json_get_or_default(json, "end_time", end_time, huge(0.0_rp))

    call this%init_from_attributes(design, simulation, weight, name, &
         mask_name, dealias_sensitivity, dealias_forcing, start_time, &
         end_time)
  end subroutine brinkman_dissipation_init_json_sim

  !> The actual constructor.
  !! @param this The objective.
  !! @param design the design.
  !! @param simulation the simulation.
  !! @param weight the weight of the objective function.
  !! @param name the name of the objective.
  !! @param mask_name the name of the mask.
  !! @param dealias_sensitivity use dealiasing on the sensitivity.
  !! @param dealias_forcing use dealiasing on the adjoint forcing.
  !! @param start_time start of the integration window.
  !! @param end_time end of the integration window.
  subroutine brinkman_dissipation_init_attributes(this, design, simulation, &
       weight, name, mask_name, dealias_sensitivity, dealias_forcing, &
       start_time, end_time)
    class(brinkman_dissipation_objective_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(simulation_t), target, intent(inout) :: simulation
    real(kind=rp), intent(in) :: weight
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: mask_name
    logical, intent(in) :: dealias_sensitivity
    logical, intent(in) :: dealias_forcing
    real(kind=rp), intent(in) :: start_time
    real(kind=rp), intent(in) :: end_time

    type(adjoint_brinkman_dissipation_source_term_t) :: brinkman_dissipation

    ! Call the base initializer
    call this%init_base(name, design%size(), weight, mask_name, &
         start_time, end_time)

    this%dealias_forcing = dealias_forcing
    this%dealias_sensitivity = dealias_sensitivity

    ! Grab the Brinkman amplitude field.
    select type (design)
    type is (brinkman_design_t)
       this%brinkman_amplitude => &
            neko_registry%get_field("brinkman_amplitude")


    class default
       call neko_error('Brinkman dissipation only works with '// &
            'brinkman_design')
    end select

    this%u => neko_registry%get_field('u')
    this%v => neko_registry%get_field('v')
    this%w => neko_registry%get_field('w')

    ! GLL
    this%c_Xh_GLL => simulation%neko_case%fluid%c_Xh
    this%Xh_GLL => simulation%neko_case%fluid%c_Xh%Xh
    this%gdim = this%c_Xh_GLL%msh%gdim

    ! GL
    this%c_Xh_GL => simulation%adjoint_case%fluid_adj%c_Xh_GL
    this%Xh_GL => this%c_Xh_GL%Xh

    ! GLL to GL
    this%GLL_to_GL => simulation%adjoint_case%fluid_adj%GLL_to_GL

    this%scratch_GL => simulation%adjoint_case%fluid_adj%scratch_GL

    ! compute the volume of the objective domain
    if (this%has_mask) then
       this%volume = compute_masked_volume(this%mask, this%c_Xh_GLL)
    else
       this%volume = this%c_Xh_GLL%volume
    end if

    ! if we have the Brinkman dissipation we initialize and append it too

    associate(f_adj_x => simulation%adjoint_fluid%f_adj_x, &
         f_adj_y => simulation%adjoint_fluid%f_adj_y, &
         f_adj_z => simulation%adjoint_fluid%f_adj_z)

      call brinkman_dissipation%init_from_components(f_adj_x, f_adj_y, &
           f_adj_z, design, this%weight, this%u, this%v, this%w, this%mask, &
           this%has_mask, this%c_Xh_GLL, this%c_Xh_GL, this%GLL_to_GL, &
           this%dealias_forcing, this%volume, this%scratch_GL, this%gdim, &
           this%start_time, this%end_time)

    end associate

    ! append adjoint forcing term based on objective function
    select type (f => simulation%adjoint_fluid)
    type is (adjoint_fluid_pnpn_t)
       call f%source_term%add_source_term(brinkman_dissipation)
    end select

  end subroutine brinkman_dissipation_init_attributes

  !> Destructor.
  subroutine brinkman_dissipation_free(this)
    class(brinkman_dissipation_objective_t), intent(inout) :: this
    call this%free_base()

    this%u => null()
    this%v => null()
    this%w => null()
    this%c_Xh_GLL => null()
    this%brinkman_amplitude => null()
    nullify(this%c_Xh_GL)
    nullify(this%Xh_GL)
    nullify(this%Xh_GLL)
    nullify(this%GLL_to_GL)
    nullify(this%scratch_GL)

  end subroutine brinkman_dissipation_free

  !> Compute the objective function.
  !! @param this The objective.
  !! @param design the design.
  subroutine brinkman_dissipation_update_value(this, design)
    class(brinkman_dissipation_objective_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(field_t), pointer :: work
    integer :: temp_indices(1)

    call neko_scratch_registry%request_field(work, temp_indices(1), .false.)

    call field_col3(work, this%u, this%u)
    call field_addcol3(work, this%v, this%v)
    if (this%gdim .eq. 3) then
       call field_addcol3(work, this%w, this%w)
    end if
    call field_col2(work, this%brinkman_amplitude)

    if (this%has_mask) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          ! note, this could be done more elagantly by writing
          ! device_glsc2_mask
          call mask_exterior_const(work, this%mask, 0.0_rp)
          this%value = device_glsc2(work%x_d, this%c_Xh_GLL%B_d, design%size())
       else
          this%value = glsc2_mask(work%x, this%c_Xh_GLL%B, design%size(), &
               this%mask%mask%get(), this%mask%size)
       end if
    else
       if (neko_bcknd_device .eq. 1) then
          this%value = device_glsc2(work%x_d, this%c_Xh_GLL%B_d, design%size())
       else
          this%value = glsc2(work%x, this%c_Xh_GLL%B, design%size())
       end if
    end if
    this%value = 0.5_rp * this%value / this%volume

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine brinkman_dissipation_update_value

  !> update_value the sensitivity of the objective function with respect to
  !! \f$chi\f$
  !! @param this The objective.
  !! @param design the design.
  subroutine brinkman_dissipation_update_sensitivity(this, design)
    class(brinkman_dissipation_objective_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(field_t), pointer :: work
    integer :: temp_indices(1)
    integer :: n_GL, nel
    type(field_t), pointer :: accumulate, fld_GL
    integer :: temp_indices_GL(2)

    ! The Brinkman dissipation adds an extra term in the sensitivity.

    call neko_scratch_registry%request_field(work, temp_indices(1), .false.)

    if(this%dealias_sensitivity) then
       nel = this%c_Xh_GLL%msh%nelv
       n_GL = nel * this%Xh_GL%lxyz
       call this%scratch_GL%request_field(accumulate, temp_indices_GL(1), &
            .false.)
       call this%scratch_GL%request_field(fld_GL, temp_indices_GL(2), .false.)

       call this%GLL_to_GL%map(fld_GL%x, this%u%x, nel, this%Xh_GL)
       call field_col3(accumulate, fld_GL, fld_GL)
       call this%GLL_to_GL%map(fld_GL%x, this%v%x, nel, this%Xh_GL)
       call field_addcol3(accumulate, fld_GL, fld_GL)
       if (this%gdim .eq. 3) then
          call this%GLL_to_GL%map(fld_GL%x, this%w%x, nel, this%Xh_GL)
          call field_addcol3(accumulate, fld_GL, fld_GL)
       end if
       ! scale
       call field_cmult(accumulate, this%weight * 0.5_rp / this%volume)

       ! Evaluate term on GL and preempt the GLL premultiplication
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_col2(accumulate%x_d, this%c_Xh_GL%B_d, n_GL)
          call this%GLL_to_GL%map(work%x, accumulate%x, nel, this%Xh_GLL)
          call device_invcol2(work%x_d, this%c_Xh_GLL%B_d, work%size())
       else
          call col2(accumulate%x, this%c_Xh_GL%B, n_GL)
          call this%GLL_to_GL%map(work%x, accumulate%x, nel, this%Xh_GLL)
          call invcol2(work%x, this%c_Xh_GLL%B, work%size())
       end if

       call this%scratch_GL%relinquish_field(temp_indices_GL)

    else
       call field_col3(work, this%u, this%u)
       call field_addcol3(work, this%v, this%v)
       if (this%gdim .eq. 3) then
          call field_addcol3(work, this%w, this%w)
       end if
       ! scale
       call field_cmult(work, this%weight * 0.5_rp / this%volume)
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(this%sensitivity%x_d, work%x_d, this%sensitivity%size())
    else
       call copy(this%sensitivity%x, work%x, this%sensitivity%size())
    end if

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine brinkman_dissipation_update_sensitivity

end module brinkman_dissipation_objective
