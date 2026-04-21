!> @file mapping.f90
!! @copyright
!! Copyright (c) 2024-2025, The Neko-TOP Authors
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
!> Mappings to be applied to a scalar field
module mapping
  use num_types, only: rp
  use json_module, only: json_file
  use coefs, only: coef_t
  use field, only: field_t
  use field_math, only: field_copy
  use registry, only: neko_registry
  implicit none
  private

  !> Base abstract class for mapping.
  type, abstract, public :: mapping_t
     !> Coefficients for the SEM.
     type(coef_t), pointer :: coef => null()
     !> A copy of the unmapped field (often used for chain rule)
     type(field_t), pointer :: X_in => null()
     !> A name for the field
     character(len=80) :: fld_name = ""

   contains
     !> Constructor for the mapping_t class.
     procedure, pass(this) :: init_base => mapping_init_base
     !> Destructor for the mapping_t (base) class.
     procedure, pass(this) :: free_base => mapping_free_base
     !> Apply the forward mapping.
     procedure, pass(this) :: apply_forward => mapping_apply_forward_wrapper
     !> Apply the backward mapping (ie, chain rule)
     procedure, pass(this) :: apply_backward => mapping_apply_backward_wrapper
     !> The common constructor using a JSON dictionary.
     procedure(mapping_init), pass(this), deferred :: init
     !> Destructor.
     procedure(mapping_free), pass(this), deferred :: free
     !> forward mapping to be computed
     procedure(mapping_forward_mapping), pass(this), deferred :: forward_mapping
     !> Backward mapping to be computed
     procedure(mapping_backward_mapping), pass(this), deferred :: &
          backward_mapping
  end type mapping_t

  !> A helper type that is needed to have an array of polymorphic objects
  type, public :: mapping_wrapper_t
     !> Wrapped polymorphic mapping.
     class(mapping_t), allocatable :: mapping
   contains
     !> Destructor.
     procedure, pass(this) :: free => mapping_wrapper_free
  end type mapping_wrapper_t

  abstract interface
     !> The common constructor using a JSON dictionary.
     !! @param this The mapping object.
     !! @param json The JSON with properties.
     !! @param case The case_t object.
     subroutine mapping_init(this, json, coef)
       import mapping_t, json_file, coef_t
       class(mapping_t), intent(inout) :: this
       type(json_file), intent(inout) :: json
       type(coef_t), intent(inout) :: coef
     end subroutine mapping_init
  end interface

  abstract interface
     !> Destructor.
     subroutine mapping_free(this)
       import mapping_t
       class(mapping_t), intent(inout) :: this
     end subroutine mapping_free
  end interface

  abstract interface
     !> The application of the mapping (\f$\rho \mapsto \tilde{\rho}\f$).
     !! @param this The mapping object
     !! @param X_out The mapped field (\f$\tilde{\rho}\f$)
     !! @param X_in The unmapped field (\f$\rho\f$)
     subroutine mapping_forward_mapping(this, X_out, X_in)
       import mapping_t, field_t
       class(mapping_t), intent(inout) :: this
       type(field_t), intent(in) :: X_in
       type(field_t), intent(inout) :: X_out
     end subroutine mapping_forward_mapping
  end interface

  abstract interface
     !> The application of the mapping backward with chain rule).
     !! \f$\frac{\partial F}{\partial \tilde{\rho}} \mapsto
     !! \frac{\partial F}{\partial \rho}\f$
     !! @param this The mapping object
     !! @param X_in The original input field (\f$\rho\f$)
     !! @param sens_in sensitivity wrt to the mapped field
     !! (\f$\frac{\partial F}{\partial \tilde{\rho}}\f$)
     !! @param sens_out sensitivity wrt to the unmapped field
     !! (\f$\frac{\partial F}{\partial \rho}\f$)
     subroutine mapping_backward_mapping(this, sens_out, sens_in, X_in)
       import mapping_t, field_t
       class(mapping_t), intent(inout) :: this
       type(field_t), intent(in) :: sens_in
       type(field_t), intent(in) :: X_in
       type(field_t), intent(inout) :: sens_out
     end subroutine mapping_backward_mapping
  end interface

  !> mapping factory. Both constructs and initializes the object.
  !! @param object The mapping object.
  !! @param json JSON object initializing the mapping.
  !! @param coef The SEM coefficients.
  interface mapping_factory
     module subroutine mapping_factory(object, json, coef)
       class(mapping_t), allocatable, intent(inout) :: object
       type(json_file), intent(inout) :: json
       type(coef_t), intent(inout) :: coef
     end subroutine mapping_factory
  end interface mapping_factory

  public :: mapping_factory
contains

  !> Constructor for the `mapping_t` (base) class.
  subroutine mapping_init_base(this, json, coef, fld_name)
    class(mapping_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: fld_name
    type(coef_t), intent(inout), target :: coef

    this%coef => coef
    this%fld_name = fld_name
    ! The mapping handler will take care of naming this field
    call neko_registry%add_field(coef%dof, this%fld_name)
    this%X_in => neko_registry%get_field(this%fld_name)

  end subroutine mapping_init_base

  !> Destructor for the `mapping_t` (base) class.
  subroutine mapping_free_base(this)
    class(mapping_t), intent(inout) :: this

    nullify(this%X_in)
    nullify(this%coef)

  end subroutine mapping_free_base

  !> Destructor for the `mapping_wrapper_t` type.
  subroutine mapping_wrapper_free(this)
    class(mapping_wrapper_t), intent(inout) :: this

    if (allocated(this%mapping)) then
       call this%mapping%free()
       deallocate(this%mapping)
    end if
  end subroutine mapping_wrapper_free

  !> Executes `apply_forward`
  !! @param this The mapping object
  !! @param X_out The mapped field (\f$\tilde{\rho}\f$)
  !! @param X_in The unmapped field (\f$\rho\f$)
  subroutine mapping_apply_forward_wrapper(this, X_out, X_in)
    class(mapping_t), intent(inout) :: this
    type(field_t), intent(inout) :: X_out
    type(field_t), intent(in) :: X_in

    call field_copy(this%X_in, X_in)
    call this%forward_mapping(X_out, this%X_in)

  end subroutine mapping_apply_forward_wrapper

  !> Executes `apply_backward`
  !! \f$\frac{\partial F}{\partial \tilde{\rho}} \mapsto
  !! \frac{\partial F}{\partial \rho}\f$
  !! @param this The mapping object
  !! @param sens_in sensitivity wrt to the mapped field
  !! (\f$\frac{\partial F}{\partial \tilde{\rho}}\f$)
  !! @param sens_out sensitivity wrt to the unmapped field
  !! (\f$\frac{\partial F}{\partial \rho}\f$)
  subroutine mapping_apply_backward_wrapper(this, sens_out, sens_in)
    class(mapping_t), intent(inout) :: this
    type(field_t), intent(inout) :: sens_out
    type(field_t), intent(in) :: sens_in
    ! @todo
    ! hmmmm, it would be silly to call mapping backward without mapping forward
    ! but at least this%X_in is certainly initialized.
    ! but it won't contain the correct information unless a map forward has
    ! occured.
    call this%backward_mapping(sens_out, sens_in, this%X_in)

  end subroutine mapping_apply_backward_wrapper

end module mapping
