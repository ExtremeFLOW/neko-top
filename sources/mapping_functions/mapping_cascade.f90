
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
!> Implements the `brinkman_mapping_t` type.
module brinkman_mapping
  use fluid_user_mapping, only: fluid_user_mapping_t
  use mapping, only: mapping_t
  use mapping_handler, only: mapping_handler_t
  use field, only: field_t
  use field_list, only: field_list_t
  use coefs, only: coef_t
  use user_intf, only: user_t
  implicit none
  private

  !> Wrapper contaning and executing the fluid source terms.
  !! @details
  !! Exists mainly to keep the `fluid_scheme_incompressible_t` type smaller and
  !! also as placeholder for future optimizations.
  type, public, extends(mapping_handler_t) :: brinkman_mapping_t

   contains
     !> Constructor.
     procedure, pass(this) :: init => brinkman_mapping_init
     !> Initialize the user source term.
     procedure, nopass :: init_user_source => fluid_init_user_source

  end type brinkman_mapping_t

contains

  !> Constructor.
  subroutine brinkman_mapping_init(this, coef)
    class(brinkman_mapping_t), intent(inout) :: this
    type(coef_t), target, intent(in) :: coef

    call this%init_base(coef)
  end subroutine brinkman_mapping_init

! This may be stupid.. It was way too heavily inspired by the source terms.
! So Tim, double check if this was the right way of doing things or this has
! just made the mapping
  ! !> Initialize the user source term.
  ! !! @param mapping The allocatable source term to be initialized to a user.
  ! !! @param rhs_fields The field list with the 3 right-hand-side components.
  ! !! @param coef The SEM coefs.
  ! !! @param type The type of the user source term, "user_vector" or
  ! !! "user_poinwise".
  ! !! @param user The user type containing the user source term routines.
  ! subroutine fluid_init_user_source(mapping, rhs_fields, coef, type, user)
  !   class(mapping_t), allocatable, intent(inout) :: mapping
  !   type(field_list_t) :: rhs_fields
  !   type(coef_t), intent(in) :: coef
  !   character(len=*) :: type
  !   type(user_t), intent(in) :: user

  !   allocate(fluid_user_mapping_t::mapping)

  !   select type (mapping)
  !   type is (fluid_user_mapping_t)
  !      call mapping%init_from_components(rhs_fields, coef, type, &
  !           user%fluid_user_f_vector, &
  !           user%fluid_user_f)
  !   end select
  ! end subroutine fluid_init_user_source

end module brinkman_mapping
