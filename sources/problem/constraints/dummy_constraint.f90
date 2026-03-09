!> @file dummy_constraint.f90
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
!> Implements the `dummy_constraint_t` type.
! $C = -1$
! $dC/dx = 0$
module dummy_constraint
  use constraint, only: constraint_t
  use json_module, only: json_file

  use design, only: design_t

  use num_types, only: rp
  use vector_math, only: vector_add2, vector_cmult

  implicit none
  private

  !> A dummy constraint
  type, public, extends(constraint_t) :: dummy_constraint_t
     private
   contains
     !> The common constructor using a JSON object.
     procedure, public, pass(this) :: init_json => &
          dummy_constraint_init_json
     !> The direct initializer from attributes.
     procedure, public, pass(this) :: init_from_attributes => &
          dummy_constraint_init_attributes
     !> Destructor.
     procedure, public, pass(this) :: free => dummy_constraint_free
     !> Computes the source term and adds the result to `fields`.
     procedure, public, pass(this) :: update_value => &
          dummy_constraint_update_value
     !> Computes the source term and adds the result to `fields`.
     procedure, public, pass(this) :: update_sensitivity => &
          dummy_constraint_update_sensitivity
     !> Get number of log entries
     procedure, public, pass(this) :: get_log_size => &
          dummy_constraint_get_log_size
     !> Get header labels for log entries
     procedure, public, pass(this) :: get_log_headers => &
          dummy_constraint_get_log_headers
     !> Get values for log entries
     procedure, public, pass(this) :: get_log_values => &
          dummy_constraint_get_log_values
  end type dummy_constraint_t

contains

  !> The common constructor using a JSON object.
  subroutine dummy_constraint_init_json(this, json, design)
    class(dummy_constraint_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(design_t), intent(in) :: design
    call this%init_from_attributes(design)
  end subroutine dummy_constraint_init_json

  !> The direct initializer from attributes.
  subroutine dummy_constraint_init_attributes(this, design)
    class(dummy_constraint_t), intent(inout) :: this
    class(design_t), intent(in) :: design

    ! Initialize the base class
    call this%init_base("dummy_constraint", design%size())

    ! ------------------------------------------------------------------------ !
    ! Initialize the value of constraint

    this%value = -1.0_rp

    ! ------------------------------------------------------------------------ !
    ! Initialize the sensitivity value

    this%sensitivity = 0.0_rp

  end subroutine dummy_constraint_init_attributes

  !> Destructor.
  subroutine dummy_constraint_free(this)
    class(dummy_constraint_t), intent(inout) :: this

    call this%free_base()
  end subroutine dummy_constraint_free

  !> The computation of the constraint.
  subroutine dummy_constraint_update_value(this, design)
    class(dummy_constraint_t), intent(inout) :: this
    class(design_t), intent(in) :: design
  end subroutine dummy_constraint_update_value

  !> The computation of the sensitivity.
  subroutine dummy_constraint_update_sensitivity(this, design)
    class(dummy_constraint_t), intent(inout) :: this
    class(design_t), intent(in) :: design
  end subroutine dummy_constraint_update_sensitivity

  !> Return number of log entries for dummy constraint.
  !! @param[in] this The constraint object.
  !! @return n Number of log entries.
  function dummy_constraint_get_log_size(this) result(n)
    class(dummy_constraint_t), intent(in) :: this
    integer :: n

    n = 0
  end function dummy_constraint_get_log_size

  !> Populate log header labels for dummy constraint.
  !! @param[in] this The constraint object.
  !! @param[out] headers Header labels for each log entry.
  subroutine dummy_constraint_get_log_headers(this, headers)
    class(dummy_constraint_t), intent(in) :: this
    character(len=*), intent(out) :: headers(:)

    if (size(headers) .eq. 0) return
    headers = ""
  end subroutine dummy_constraint_get_log_headers

  !> Populate log values for dummy constraint.
  !! @param[in] this The constraint object.
  !! @param[out] values Values corresponding to the log headers.
  subroutine dummy_constraint_get_log_values(this, values)
    class(dummy_constraint_t), intent(in) :: this
    real(kind=rp), intent(out) :: values(:)

    if (size(values) .eq. 0) return
    values = 0.0_rp
  end subroutine dummy_constraint_get_log_values

end module dummy_constraint
