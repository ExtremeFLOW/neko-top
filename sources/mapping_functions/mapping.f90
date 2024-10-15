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
!
!> Mappings to be applied to a scalar field
module mapping
  use num_types, only : rp
  use json_module, only : json_file
  use coefs, only : coef_t
  use json_utils, only : json_get_or_default, json_get
  use field, only: field_t
  implicit none
  private

  !> Base abstract class for mapping.
  type, abstract, public :: mapping_t
  	  !> Coefficients for the SEM.
     type(coef_t), pointer :: coef => null()
 
   contains
     !> Constructor for the mapping_t class.
     procedure, pass(this) :: init_base => mapping_init_base
     !> Destructor for the mapping_t (base) class.
     procedure, pass(this) :: free_base => mapping_free_base
     !> The common constructor using a JSON dictionary.
     procedure(mapping_init), pass(this), deferred :: init
     !> Destructor.
     procedure(mapping_free), pass(this), deferred :: free
     !> Apply forward
     procedure(mapping_apply), pass(this), deferred :: apply_forward
     !> Apply backwards (with chain rule)
     procedure(mapping_apply_backward), pass(this), deferred :: apply_backward
  end type mapping_t




  abstract interface
     !> The common constructor using a JSON dictionary.
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
     !> The application of the mapping ($\rho \mapsto \tilde{\rho}$).
     !! @param X_out The mapped field ($\tilde{\rho}$)
     !! @param X_in The unmapped field ($\rho$)
     subroutine mapping_apply(this, X_out, X_in)
       import mapping_t, field_t
       class(mapping_t), intent(inout) :: this
       type(field_t), intent(in) ::  X_in
       type(field_t), intent(inout) ::  X_out
     end subroutine mapping_apply
  end interface

  abstract interface
     !> The application of the mapping backward with chain rule).
     !! $\frac{\partial F}{\partial \tilde{\rho}} \mapsto
     !! \frac{\partial F}{\partial \rho}$
     !! @param X_in The original input field ($\rho$)
     !! @param dF_dX_out, sensitivity wrt to the mapped field
     !! ($\frac{\partial F}{\partial \tilde{\rho}}$)
     !! @param dF_dX_in, sensitivity wrt to the unmapped field
     !! ($\frac{\partial F}{\partial \rho}$)
     subroutine mapping_apply_backward(this, dF_dX_in, dF_dX_out, X_in)
       import mapping_t, field_t
       class(mapping_t), intent(inout) :: this
       type(field_t), intent(in) ::  dF_dX_out
       type(field_t), intent(in) ::  X_in
       type(field_t), intent(inout) ::  dF_dX_in
     end subroutine mapping_apply_backward
  end interface
contains

  !> Constructor for the `mapping_t` (base) class.
  subroutine mapping_init_base(this, json, coef)
    class(mapping_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(coef_t), intent(inout), target :: coef
    character(len=:), allocatable :: compute_control, output_control
    real(kind=rp) :: compute_value, output_value
    integer :: order


    this%coef => coef

  end subroutine mapping_init_base

  !> Destructor for the `mapping_t` (base) class.
  subroutine mapping_free_base(this)
    class(mapping_t), intent(inout) :: this

    nullify(this%coef)
  end subroutine mapping_free_base



end module mapping
