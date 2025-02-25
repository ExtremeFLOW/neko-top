! Copyright (c) 2023, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!    *  Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!    *  Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!    *  Neither the name of the authors nor the names of its
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
!> This provides a template to be modified in constructing new objectives.
!> For a detailed description see ///
!
!------------------------------------------------------------------------------
! TO BE FILLED: a description of the objective
!------------------------------------------------------------------------------
!
module TEMPLATE_objective
  use num_types, only: rp
  use objective, only: objective_t
  use simulation, only: simulation_t
  use design, only: design_t
  use json_module, only: json_file
  use json_utils, only: json_get_or_default
  !----------------------------------------------------------------------------
  ! TO BE FILLED: Add additional modules
  !----------------------------------------------------------------------------
  implicit none
  private

  !----------------------------------------------------------------------------
  !> TO BE FILLED: a description of the objective
  !----------------------------------------------------------------------------
  type, public, extends(objective_t) :: TEMPLATE_objective_t
     private

     !-------------------------------------------------------------------------
     ! TO BE FILLED: additional private variables used by you objective
     !-------------------------------------------------------------------------
     ! eg,

     ! !> pointers to the primal velocity fields
     ! type(field_t), pointer :: u, v, w

   contains

     !> The common constructor using a JSON object.
     procedure, public, pass(this) :: init_json => TEMPLATE_init_json
     !> The actual constructor.
     procedure, public, pass(this) :: init_from_attributes => &
          TEMPLATE_init_attributes
     !> Destructor.
     procedure, public, pass(this) :: free => TEMPLATE_free
     !> Computes the value of the objective function.
     procedure, public, pass(this) :: update_value => &
          TEMPLATE_update_value
     !> Computes the sensitivity with respect to the coefficient $\chi$.
     procedure, public, pass(this) :: update_sensitivity => &
          TEMPLATE_update_sensitivity

  end type TEMPLATE_objective_t

contains

  !> The common constructor using a JSON object.
  !! @param json The JSON subdictionary corresponding to your objective
  !! @param design The design
  !! @param simulation The simulation
  subroutine TEMPLATE_init_json(this, json, design, simulation)
    class(TEMPLATE_objective_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(design_t), intent(in) :: design
    type(simulation_t), target, intent(inout) :: simulation
    real(kind=rp) :: weight
    character(len=*) :: mask_name

    call json_get_or_default(json, "weight", weight, 1.0_rp)
    call json_get_or_default(json, "mask_name", mask_name, "")
    !--------------------------------------------------------------------------
    ! TO BE FILLED: Read your additional parameters from the JSON
    !--------------------------------------------------------------------------

    ! initialize
    call this%init_from_attributes(design, simulation, weight, mask_name)
  end subroutine TEMPLATE_init_json

  !> The actual constructor.
  !! @param design the design.
  !! @param simulation the simulation.
  !! @param weight the weight of the objective function.
  !! @param mask_name the name of the mask.
  !----------------------------------------------------------------------------
  ! TO BE FILLED: Document your additional parameters
  !----------------------------------------------------------------------------
  subroutine TEMPLATE_init_attributes(this, design, simulation, weight, &
       mask_name)
    class(TEMPLATE_objective_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(simulation_t), target, intent(inout) :: simulation
    real(kind=rp), intent(in) :: weight
    character(len=*), intent(in) :: mask_name

    ! Call the base initializer
    call this%init_base(design%size(), weight, mask_name)

    !--------------------------------------------------------------------------
    ! TO BE FILLED: Associate various fields from the primal in the registry
    !--------------------------------------------------------------------------
    ! eg

    ! this%u => neko_field_registry%get_field('u')
    ! this%v => neko_field_registry%get_field('v')
    ! this%w => neko_field_registry%get_field('w')


    !--------------------------------------------------------------------------
    ! TO BE FILLED: Initialize adjoint source terms
    !--------------------------------------------------------------------------
    ! eg,

    ! associate(f_adj_x => simulation%adjoint_case%scheme%f_adj_x, &
    !      f_adj_y => simulation%adjoint_case%scheme%f_adj_y, &
    !      f_adj_z => simulation%adjoint_case%scheme%f_adj_z, &
    !      c_Xh => simulation%adjoint_case%scheme%c_Xh)

    !   call my_adjoint_source_term%init_from_components(f_adj_x, f_adj_y, &
    !        f_adj_z, design, this%k * this%weight, this%u, this%v, this%w, &
    !        this%mask, this%has_mask, c_Xh)

    ! end associate

    !--------------------------------------------------------------------------
    ! TO BE FILLED: append adjoint forcing term
    !--------------------------------------------------------------------------
    ! eg.

    ! call simulation%adjoint_case%scheme%source_term%add_source_term(&
    !    my_adjoint_source_term)

  end subroutine TEMPLATE_init_attributes

  !> Destructor.
  subroutine TEMPLATE_free(this)
    class(TEMPLATE_objective_t), intent(inout) :: this
    call this%free_base()

    !--------------------------------------------------------------------------
    ! TO BE FILLED: Free everything
    !--------------------------------------------------------------------------
    ! eg,

    ! this%u => null()
    ! this%v => null()
    ! this%w => null()
    ! this%B => null()

  end subroutine TEMPLATE_free

  !> Compute the objective function.
  !! @param design the design.
  subroutine TEMPLATE_update_value(this, design)
    class(TEMPLATE_objective_t), intent(inout) :: this
    class(design_t), intent(in) :: design

    !--------------------------------------------------------------------------
    ! TO BE FILLED: Compute you objective value
    !--------------------------------------------------------------------------
    ! eg,

    ! this%value = my_value

  end subroutine TEMPLATE_update_value

  !> Compute the sensitivity of the objective function with respect to the
  !! design
  !! @param design the design.
  subroutine TEMPLATE_update_sensitivity(this, design)
    class(TEMPLATE_objective_t), intent(inout) :: this
    class(design_t), intent(in) :: design

    !--------------------------------------------------------------------------
    ! TO BE FILLED: Compute your sensitivity with respect to the design
    !--------------------------------------------------------------------------
    ! eg,
    ! call field_rzero(this%sensitivity)

  end subroutine TEMPLATE_update_sensitivity

end module TEMPLATE_objective
