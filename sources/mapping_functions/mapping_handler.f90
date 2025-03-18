
! Copyright (c) 2023, The Neko Authors
! All rights reserved.
!
! Redistribution and use in mapping and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of mapping code must retain the above copyright
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
!> Implements the `mapping_handler_t` type.
module mapping_handler
  use neko_config, only: NEKO_BCKND_DEVICE
  use num_types, only: rp
  use mapping, only: mapping_wrapper_t, mapping_t, &
       mapping_factory
  use field, only: field_t
  use field_list, only: field_list_t
  use json_utils, only: json_get, json_extract_item, json_get_or_default
  use json_module, only: json_file
  use coefs, only: coef_t
  use user_intf, only: user_t
  use field_math, only: field_rzero, field_copy
  use math, only: col2
  use device_math, only: device_col2
  use scratch_registry, only: neko_scratch_registry
  use utils, only: neko_warning
  implicit none
  private

  !> Abstract class for handling mapping_cascade.
  !!
  !! @details
  !! This class is responsible for managing the mapping_cascade in a sequential
  !! manor. It is also responsible for using the chain rule to propogate
  !! sensitivity backwards throughout the system.
  type, abstract, public :: mapping_handler_t
     !> Array of mapping_cascade.
     !! @note the order really matter's here since they'll be executed in
     !! sequence.
     class(mapping_wrapper_t), allocatable :: mapping_cascade(:)
     !> The right-hand side.
     type(field_list_t) :: rhs_fields
     !> The coefficients of the (space, mesh) pair.
     type(coef_t), pointer :: coef
     ! @todo, we should handle user mapping_cascade in a different PR
     !!> The user object.
     !!type(user_t), pointer :: user

   contains
     !> Constructor.
     procedure, pass(this) :: init_base => mapping_handler_init_base
     !> Destructor.
     procedure, pass(this) :: free => mapping_handler_free
     !> Cycle through all the mapping_cascade and return the final field
     procedure, pass(this) :: apply_forward => mapping_handler_apply_forward
     !> Cycle backwards through all the mapping_cascade and return the sensitivity
     procedure, pass(this) :: apply_backward => mapping_handler_apply_backward
     !> Generic interface to add a mapping to the list.
     generic :: add => add_mapping, add_json_mappings
     !> Append a new mapping to the mapping_cascade array.
     procedure, pass(this) :: add_mapping => &
          mapping_handler_add_mapping
     !> Read from the json file and initialize the mapping_cascade.
     procedure, pass(this) :: add_json_mappings => &
          mapping_handler_add_json_mappings
     !!> Initialize the user mapping.
     ! procedure(mapping_handler_init_user_mapping), &
     !     nopass, deferred :: init_user_mapping
  end type mapping_handler_t

!   abstract interface
!      subroutine mapping_handler_init_user_mapping(mapping, rhs_fields, &
!           coef, type, user)
!        import :: mapping_t, field_list_t, coef_t, user_t
!        class(mapping_t), allocatable, intent(inout) :: mapping
!        type(field_list_t) :: rhs_fields
!        type(coef_t), intent(in) :: coef
!        character(len=*) :: type
!        type(user_t), intent(in) :: user
!      end subroutine mapping_handler_init_user_mapping
!   end interface

contains

  !> Constructor.
  subroutine mapping_handler_init_base(this, coef)
    class(mapping_handler_t), intent(inout) :: this
    type(coef_t), target, intent(in) :: coef

    call this%free()

    this%coef => coef

  end subroutine mapping_handler_init_base


  !> Destructor.
  subroutine mapping_handler_free(this)
    class(mapping_handler_t), intent(inout) :: this
    integer :: i

    if (allocated(this%mapping_cascade)) then
       do i = 1, size(this%mapping_cascade)
          call this%mapping_cascade(i)%free()
       end do
       deallocate(this%mapping_cascade)
    end if

  end subroutine mapping_handler_free

  !> apply the cascade of mapping_cascade.
  !! @param X_out The mapped field ($\tilde{\rho}$)
  !! @param X_in The unmapped field ($\rho$)
  subroutine mapping_handler_apply_forward(this, X_out, X_in)
    class(mapping_handler_t), intent(inout) :: this
    type(field_t), intent(in) :: X_in
    type(field_t), intent(inout) :: X_out
    integer :: i
    type(field_t), pointer :: tmp_fld_in, tmp_fld_out
    integer :: temp_indices(2)

    call neko_scratch_registry%request_field(tmp_fld_in, temp_indices(1))
    call neko_scratch_registry%request_field(tmp_fld_out, temp_indices(2))

    ! Start by copying the first X_in into the tmp_fld_out to begin the
    ! cascade.
    call field_copy(tmp_fld_out, X_in)

    ! We enter the cascade
    if (allocated(this%mapping_cascade)) then
       do i = 1, size(this%mapping_cascade)
          ! the output from one mapping becomes the input for the next.
          call field_copy(tmp_fld_in, tmp_fld_out)
          ! apply the mapping on temp_fld
          call this%mapping_cascade(i)%mapping%apply_forward(tmp_fld_out, &
               tmp_fld_in)

       end do

    end if

    ! our final mapping should now live in tmp_fld_out
    call field_copy(X_out, tmp_fld_out)

    ! free the scratch.
    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine mapping_handler_apply_forward

  !> apply the cascade of mapping_cascade.
  !! @param X_out The mapped field ($\tilde{\rho}$)
  !! @param X_in The unmapped field ($\rho$)
  subroutine mapping_handler_apply_backward(this, sens_out, sens_in)
    class(mapping_handler_t), intent(inout) :: this
    type(field_t), intent(in) :: sens_in
    type(field_t), intent(inout) :: sens_out
    integer :: i
    type(field_t), pointer :: tmp_fld_in, tmp_fld_out
    integer :: temp_indices(2)

    call neko_scratch_registry%request_field(tmp_fld_in, temp_indices(1))
    call neko_scratch_registry%request_field(tmp_fld_out, temp_indices(2))

    ! Start by copying the first sens_in into the tmp_fld_out to begin the
    ! cascade.
    call field_copy(tmp_fld_out, sens_in)

    ! We enter the cascade
    if (allocated(this%mapping_cascade)) then
       do i = 1, size(this%mapping_cascade)
          ! the output from one mapping becomes the input for the next.
          call field_copy(tmp_fld_in, tmp_fld_out)
          ! apply the mapping on temp_fld
          ! NOTE
          ! all the X_in that is required to map backward should be held
          ! internally by each mapping
          call this%mapping_cascade(i)%mapping%apply_backward(tmp_fld_out, &
               tmp_fld_in)

       end do

    end if

    ! our final mapping should now live in tmp_fld_out
    call field_copy(sens_out, tmp_fld_out)

    ! free the scratch.
    call neko_scratch_registry%relinquish_field(temp_indices)


  end subroutine mapping_handler_apply_backward

  !> Read from the json file and initialize the mapping_cascade.
  subroutine mapping_handler_add_json_mappings(this, json, name)
    class(mapping_handler_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: name

    class(mapping_wrapper_t), dimension(:), allocatable :: temp

    ! A single mapping as its own json_file.
    type(json_file) :: mapping_subdict
    character(len=:), allocatable :: type
    integer :: n_mappings, i, i0

    if (json%valid_path(name)) then
       ! Get the number of mapping_cascade.
       call json%info(name, n_children = n_mappings)

       if (allocated(this%mapping_cascade)) then
          i0 = size(this%mapping_cascade)
          call move_alloc(this%mapping_cascade, temp)
          allocate(this%mapping_cascade(i0 + n_mappings))
          if (allocated(temp)) then
             do i = 1, i0
                call move_alloc(temp(i)%mapping, &
                     this%mapping_cascade(i)%mapping)
             end do
          end if
       else
          i0 = 0
          allocate(this%mapping_cascade(n_mappings))
       end if

       do i = 1, n_mappings
          ! Create a new json containing just the subdict for this mapping.
          call json_extract_item(json, name, i, mapping_subdict)
          call json_get(mapping_subdict, "type", type)
          call mapping_factory(this%mapping_cascade(i + i0)%mapping, &
               mapping_subdict, this%coef)
       end do
    else
       ! I was considering an error, but maybe a warning is more appropriate
       call neko_warning("No mappings selected")
    end if

  end subroutine mapping_handler_add_json_mappings

  !> Add new mapping to the list.
  !! @param mapping The mapping to be added.
  subroutine mapping_handler_add_mapping(this, mapping)
    class(mapping_handler_t), intent(inout) :: this
    class(mapping_t), intent(in) :: mapping
    class(mapping_wrapper_t), dimension(:), allocatable :: temp

    integer :: n_mappings, i

    if (allocated(this%mapping_cascade)) then
       n_mappings = size(this%mapping_cascade)
    else
       n_mappings = 0
    end if

    call move_alloc(this%mapping_cascade, temp)
    allocate(this%mapping_cascade(n_mappings + 1))

    if (allocated(temp)) then
       do i = 1, n_mappings
          call move_alloc(temp(i)%mapping, &
               this%mapping_cascade(i)%mapping)
       end do
    end if

    this%mapping_cascade(n_mappings + 1)%mapping = mapping

  end subroutine mapping_handler_add_mapping
end module mapping_handler
