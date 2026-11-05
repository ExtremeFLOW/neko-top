!> @file objective.f90
!! @copyright
!! Copyright (c) 2025-2026, The Neko-TOP Authors
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

!> Implements the `objective_t` type.
module objective
  use base_functional, only: base_functional_t
  use simulation_m, only: simulation_t
  use design, only: design_t
  use num_types, only: rp, xp
  use point_zone_registry, only: neko_point_zone_registry
  use json_module, only: json_file
  use time_state, only: time_state_t
  implicit none
  private

  public :: objective_t, objective_factory, objective_wrapper_t

  !> The abstract objective type.
  !!
  !! This is the base class for objectives, which is a type of base functional.
  !! Each objective contain a weight that is used to scale the objective value.
  type, abstract, extends(base_functional_t) :: objective_t
     !> Weight of the objective in the overall cost function
     real(kind=rp) :: weight = 1.0_rp
     !> Scale of the objective's adjoint forcing. The forcing source term
     !! is a copy held by the adjoint, and the objective is moved after
     !! building it, so the two share this scale through an allocated
     !! pointer. It is one outside an unsteady adjoint run, during which
     !! `update_forcing_scale` sets it before every adjoint step.
     real(kind=rp), pointer :: forcing_scale => null()

   contains

     !> Initializer for the base class
     procedure, pass(this) :: init_base => objective_init_base
     !> Destructor of the base class
     procedure, pass(this) :: free_base => objective_free_base
     !> Get the weight of the objective
     procedure, pass(this) :: get_weight => objective_get_weight
     !> Get the time window to give the adjoint forcing source term
     procedure, pass(this) :: get_forcing_window => &
          objective_get_forcing_window
     !> Set the adjoint forcing scale for one unsteady adjoint step
     procedure, pass(this) :: update_forcing_scale => &
          objective_update_forcing_scale
     !> Reset the adjoint forcing scale to one
     procedure, pass(this) :: reset_forcing_scale => &
          objective_reset_forcing_scale
     !> Get number of log entries for this objective
     procedure, pass(this) :: get_log_size => objective_get_log_size
     !> Get header labels for this objective's log entries
     procedure, pass(this) :: get_log_headers => objective_get_log_headers
     !> Get values for this objective's log entries
     procedure, pass(this) :: get_log_values => objective_get_log_values

  end type objective_t

  !> Wrapper for objectives for use in lists.
  type :: objective_wrapper_t
     class(objective_t), allocatable :: objective
   contains
     procedure, pass(this) :: free => objective_wrapper_free
  end type objective_wrapper_t

  ! -------------------------------------------------------------------------- !
  ! Explicit interfaces

  !> Factory function
  !! Allocates and initializes an objective function object
  !! @param object The objective function object to be created
  !! @param type The type of the objective function
  !! @param design The design object
  !! @param simulation The simulation object
  interface objective_factory
     module subroutine objective_factory(object, json, design, simulation)
       class(objective_t), allocatable, intent(inout) :: object
       type(json_file), intent(inout) :: json
       class(design_t), intent(in) :: design
       type(simulation_t), target, optional, intent(inout) :: simulation
     end subroutine objective_factory
  end interface objective_factory

contains

  ! -------------------------------------------------------------------------- !
  ! Implementations for the base class

  !> Initialize the objective base class.
  !! @param this The objective.
  !! @param name The name of the objective.
  !! @param design_size The number of design variables.
  !! @param weight The weight of the objective function.
  !! @param mask_name The name design the mask. [optional]
  !! @param start_time start of the integration window. [optional]
  !! @param end_time end of the integration window. [optional]
  subroutine objective_init_base(this, name, design_size, weight, mask_name, &
       start_time, end_time)
    class(objective_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer, intent(in) :: design_size
    real(kind=rp), intent(in) :: weight
    character(len=*), intent(in), optional :: mask_name
    real(kind=rp), intent(in), optional :: start_time
    real(kind=rp), intent(in), optional :: end_time

    call this%free_base()

    this%name = trim(name)
    call this%sensitivity%init(design_size)

    this%weight = weight

    allocate(this%forcing_scale)
    this%forcing_scale = 1.0_rp

    if (present(mask_name)) then
       if (mask_name .ne. "") then
          this%has_mask = .true.
          this%mask => neko_point_zone_registry%get_point_zone(mask_name)
       end if
    end if

    if (present(start_time)) then
       this%start_time = start_time
    else
       this%start_time = 0.0_rp
    end if

    if (present(end_time)) then
       this%end_time = end_time
    else
       this%end_time = huge(0.0_rp)
    end if

  end subroutine objective_init_base

  !> Free the base class
  subroutine objective_free_base(this)
    class(objective_t), target, intent(inout) :: this

    if (allocated(this%name)) deallocate(this%name)
    this%weight = 1.0_rp
    if (associated(this%forcing_scale)) deallocate(this%forcing_scale)

    this%value = 0.0_rp
    this%value_weight = 0.0_rp
    this%sensitivity_weight = 0.0_rp
    this%start_time = 0.0_rp
    this%end_time = huge(0.0_rp)
    call this%sensitivity%free()

    this%has_mask = .false.
    if (associated(this%mask)) nullify(this%mask)

  end subroutine objective_free_base

  !> Return number of log entries for this objective.
  !! @param[in] this The objective object.
  !! @return n Number of log entries.
  function objective_get_log_size(this) result(n)
    class(objective_t), intent(in) :: this
    integer :: n

    n = 2
  end function objective_get_log_size

  !> Populate log header labels for this objective.
  !! @param[in] this The objective object.
  !! @param[out] headers Header labels for each log entry.
  subroutine objective_get_log_headers(this, headers)
    class(objective_t), intent(in) :: this
    character(len=*), intent(out) :: headers(:)
    character(len=:), allocatable :: prefix

    if (size(headers) .lt. 1) return
    prefix = this%get_name()
    headers(1) = prefix
    if (size(headers) .lt. 2) return
    headers(2) = prefix // '.weight'
  end subroutine objective_get_log_headers

  !> Populate log values for this objective.
  !! @param[in] this The objective object.
  !! @param[out] values Values corresponding to the log headers.
  subroutine objective_get_log_values(this, values)
    class(objective_t), intent(in) :: this
    real(kind=rp), intent(out) :: values(:)

    if (size(values) .lt. 1) return
    values(1) = this%value
    if (size(values) .lt. 2) return
    values(2) = this%weight
  end subroutine objective_get_log_values

  ! -------------------------------------------------------------------------- !
  ! Implementations for the wrapper

  !> Free the objective wrapper.
  subroutine objective_wrapper_free(this)
    class(objective_wrapper_t), intent(inout) :: this
    if (allocated(this%objective)) then
       call this%objective%free()
       deallocate(this%objective)
    end if
  end subroutine objective_wrapper_free

  !> Get the weight of the objective
  pure function objective_get_weight(this) result(w)
    class(objective_t), intent(in) :: this
    real(kind=rp) :: w
    w = this%weight
  end function objective_get_weight

  !> Get the time window to give the adjoint forcing source term.
  !!
  !! A source term compares its window with the adjoint clock, which on the
  !! unsteady path runs one step behind the primal state the forcing
  !! consumes. There `update_forcing_scale` gates the forcing instead, and
  !! the source term gets an open window.
  !! @param this The objective.
  !! @param simulation The simulation.
  !! @param start_time Start of the source term's window.
  !! @param end_time End of the source term's window.
  subroutine objective_get_forcing_window(this, simulation, start_time, &
       end_time)
    class(objective_t), intent(in) :: this
    type(simulation_t), intent(in) :: simulation
    real(kind=rp), intent(out) :: start_time
    real(kind=rp), intent(out) :: end_time

    if (simulation%unsteady) then
       start_time = 0.0_rp
       end_time = huge(0.0_rp)
    else
       start_time = real(this%start_time, kind=rp)
       end_time = real(this%end_time, kind=rp)
    end if
  end subroutine objective_get_forcing_window

  !> Set the adjoint forcing scale for one unsteady adjoint step.
  !!
  !! The objective is the mean over its window, of length `value_weight`,
  !! and its forcing is the plain derivative of the integrand. The design
  !! coupling term is the mean over the whole run instead. The forcing is
  !! therefore scaled by the run length over the window length, and is off
  !! when the primal sample the step consumes lies outside the window. For
  !! an objective without a window the scale is exactly one.
  !! @param this The objective.
  !! @param time The time state of the primal sample the step consumes.
  !! @param run_time Length of the forward run, summed as `value_weight`.
  subroutine objective_update_forcing_scale(this, time, run_time)
    class(objective_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    real(kind=xp), intent(in) :: run_time

    if (this%in_window(time) .and. this%value_weight .gt. 0.0_xp) then
       this%forcing_scale = real(run_time / this%value_weight, kind=rp)
    else
       this%forcing_scale = 0.0_rp
    end if
  end subroutine objective_update_forcing_scale

  !> Reset the adjoint forcing scale to one.
  !!
  !! The unsteady path calls this after its last adjoint step, so that the
  !! scale set for that step does not outlive the run. Everywhere else, as
  !! on the steady path, the forcing is applied unscaled.
  !! @param this The objective.
  subroutine objective_reset_forcing_scale(this)
    class(objective_t), intent(inout) :: this

    this%forcing_scale = 1.0_rp
  end subroutine objective_reset_forcing_scale

end module objective

