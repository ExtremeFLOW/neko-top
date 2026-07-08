!> @file base_functional.f90
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

!> Defines the abstract the `base_functional_t` type.
module base_functional
  use design, only: design_t
  use json_module, only: json_file
  use json_utils, only: json_get
  use num_types, only: rp
  use point_zone, only: point_zone_t
  use simulation_m, only: simulation_t
  use time_state, only: time_state_t
  use vector, only: vector_t
  use utils, only: neko_error
  use vector_math, only: vector_copy, vector_add2s1
  implicit none
  private

  !> The base functional type
  !!
  !! This is the base class for objectives and constraints alike.
  !! A base functional should be able to evaluate itself and its sensitivity
  !! with respect to the design variables.
  !!
  !! The base functional is also responsible for managing the adjoint forcing
  !! terms that are required for the adjoint problem, any source terms
  !! simulation components that are required to evaluate the base functional.
  !! All of which should be prepared in the `init` method.
  type, abstract, public :: base_functional_t

     !> Value of the base_functional
     real(kind=rp) :: value = 0.0_rp
     !> Old value for time integration
     real(kind=rp) :: value_old = 0.0_rp
     !> Sensitivity field
     type(vector_t) :: sensitivity
     !> Old sensitivity field for time integration
     type(vector_t) :: sensitivity_old
     !> Name of constraint/objective in the logfile
     character(len=25) :: name = ""
     !> containing a mask
     logical :: has_mask = .false.
     !> A mask for where the objective function is evaluated
     class(point_zone_t), pointer :: mask => null()
     !> Time window start for accumulation
     real(kind=rp) :: start_time = 0.0_rp
     !> Time window end for accumulation
     real(kind=rp) :: end_time = huge(0.0_rp)

   contains

     ! ----------------------------------------------------------------------- !
     ! Derived class interfaces

     !> Constructor interface
     generic :: init => init_json, init_json_sim

     !> Constructor based on json input
     procedure, pass(this) :: init_json => functional_init_json
     !> Constructor based on json input and simulation
     procedure, pass(this) :: init_json_sim => functional_init_json_sim
     !> Destructor
     procedure(functional_free), pass(this), deferred :: free

     !> Update the value of the function
     procedure(functional_update_value), pass(this), deferred :: update_value
     !> Update the sensitivity of the function
     procedure(functional_update_sensitivity), pass(this), deferred :: &
          update_sensitivity

     !> Get the value of the function
     procedure, pass(this) :: get_value => functional_get_value
     !> Get the sensitivity of the function
     procedure, pass(this) :: get_sensitivity => functional_get_sensitivity
     !> Get number of log entries for this functional
     procedure, pass(this) :: get_log_size => functional_get_log_size
     !> Get header labels for this functional's log entries
     procedure, pass(this) :: get_log_headers => functional_get_log_headers
     !> Get values for this functional's log entries
     procedure, pass(this) :: get_log_values => functional_get_log_values
     !> Set the value to zero
     procedure, pass(this) :: reset_value => functional_reset_value
     !> Set the sensitivity to zero
     procedure, pass(this) :: reset_sensitivity => functional_reset_sensitivity
     !> Accumulate the value
     procedure, pass(this) :: accumulate_value => functional_accumulate_value
     !> Accumulate the sensitivity
     procedure, pass(this) :: accumulate_sensitivity => &
          functional_accumulate_sensitivity

  end type base_functional_t

  ! -------------------------------------------------------------------------- !
  ! Interface specifications for the derived types, these are the constructors
  ! for the different types of objective functions.

  abstract interface

     !> Destructor
     subroutine functional_free(this)
       import base_functional_t
       class(base_functional_t), intent(inout) :: this
     end subroutine functional_free

     !> Compute the objective function
     subroutine functional_update_value(this, design)
       import base_functional_t, design_t
       class(base_functional_t), intent(inout) :: this
       class(design_t), intent(in) :: design
     end subroutine functional_update_value

     !> Compute the sensitivity
     subroutine functional_update_sensitivity(this, design)
       import base_functional_t, design_t
       class(base_functional_t), intent(inout) :: this
       class(design_t), intent(in) :: design
     end subroutine functional_update_sensitivity

  end interface

contains

  !> Initialize the objective function
  subroutine functional_init_json(this, json, design)
    class(base_functional_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(design_t), intent(in) :: design
    character(len=:), allocatable :: type

    call json_get(json, 'type', type)
    call neko_error("Functional type: '" // type // &
         "' does not support initialization without simulation")
  end subroutine functional_init_json

  !> Initialize the objective function
  subroutine functional_init_json_sim(this, json, design, simulation)
    class(base_functional_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(design_t), intent(in) :: design
    type(simulation_t), target, intent(inout) :: simulation
    character(len=:), allocatable :: type

    call json_get(json, 'type', type)
    call neko_error("Functional type: '" // type // &
         "' does not support initialization with simulation")
  end subroutine functional_init_json_sim

  !> Get the value of the function
  function functional_get_value(this) result(v)
    class(base_functional_t), intent(in) :: this
    real(kind=rp) :: v

    v = this%value
  end function functional_get_value

  !> Get the sensitivity of the function
  subroutine functional_get_sensitivity(this, sensitivity)
    class(base_functional_t), intent(in) :: this
    type(vector_t), intent(inout) :: sensitivity

    sensitivity = this%sensitivity
  end subroutine functional_get_sensitivity

  !> Return number of log entries for this functional.
  !! @param[in] this The functional object.
  !! @return n Number of log entries.
  function functional_get_log_size(this) result(n)
    class(base_functional_t), intent(in) :: this
    integer :: n

    n = 1
  end function functional_get_log_size

  !> Populate log header labels for this functional.
  !! @param[in] this The functional object.
  !! @param[out] headers Header labels for each log entry.
  subroutine functional_get_log_headers(this, headers)
    class(base_functional_t), intent(in) :: this
    character(len=*), intent(out) :: headers(:)

    if (size(headers) .eq. 0) return
    headers(1) = trim(this%name)
  end subroutine functional_get_log_headers

  !> Populate log values for this functional.
  !! @param[in] this The functional object.
  !! @param[out] values Values corresponding to the log headers.
  subroutine functional_get_log_values(this, values)
    class(base_functional_t), intent(in) :: this
    real(kind=rp), intent(out) :: values(:)

    if (size(values) .eq. 0) return
    values(1) = this%value
  end subroutine functional_get_log_values

  !> Zero value of the function
  subroutine functional_reset_value(this)
    class(base_functional_t), intent(inout) :: this

    this%value = 0.0_rp
  end subroutine functional_reset_value

  !> Zero sensitivity of the function
  subroutine functional_reset_sensitivity(this)
    class(base_functional_t), intent(inout) :: this

    this%sensitivity = 0.0_rp
  end subroutine functional_reset_sensitivity

  !> Accumulate the value of the function
  subroutine functional_accumulate_value(this, design, time)
    class(base_functional_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(time_state_t), intent(in) :: time

    if (time%t .lt. this%start_time .or. time%t .gt. this%end_time) return

    this%value_old = this%value
    call this%update_value(design)

    ! could potentially use higher order trapezoidal/Simpson etc, but this
    ! should suffice
    this%value = this%value_old + this%value * time%dt / &
         (time%end_time - time%start_time)
  end subroutine functional_accumulate_value

  !> Accumulate the sensitivity of the function
  subroutine functional_accumulate_sensitivity(this, design, time)
    class(base_functional_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(time_state_t), intent(in) :: time

    if (time%t .lt. this%start_time .or. time%t .gt. this%end_time) return

    call vector_copy(this%sensitivity_old, this%sensitivity)
    call this%update_sensitivity(design)

    ! could potentially use higher order trapezoidal/Simpson etc, but this
    ! should suffice
    call vector_add2s1(this%sensitivity, this%sensitivity_old, time%dt / &
         (time%end_time - time%start_time))
  end subroutine functional_accumulate_sensitivity
end module base_functional
