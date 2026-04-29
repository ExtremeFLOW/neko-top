!> @file continuation_scheduler.f90
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

!> Continuation scheduler for the optimization loop
module continuation_scheduler
  use num_types, only: rp
  use json_module, only: json_file
  use json_utils, only: json_get_or_default, json_get_or_lookup, json_get
  use utils, only: neko_error
  implicit none
  private

  public :: continuation_scheduler_t

  !------------------------------------------------------------------
  ! Continuation parameter
  !------------------------------------------------------------------
  type continuation_parameter_t
     character(len=:), allocatable :: name ! Parameter name
     real(rp), pointer :: target => null() ! Parameter to update
     real(rp), allocatable :: values(:) ! Values to step through
     integer :: iterations_per_value = 0 ! How many iterations per value
   contains
     procedure :: update => continuation_parameter_update
     procedure :: free => continuation_parameter_free
  end type continuation_parameter_t

  !------------------------------------------------------------------
  ! Continuation scheduler
  !------------------------------------------------------------------
  type continuation_scheduler_t
     type(continuation_parameter_t), allocatable :: params(:)
     integer :: default_iterations = 1
   contains
     procedure :: init => init_scheduler
     procedure :: free => free_scheduler
     procedure :: json_get_or_register
     procedure :: register_parameter
     procedure :: update
     procedure :: get_param_name
     procedure :: get_n_params
  end type continuation_scheduler_t

  !------------------------------------------------------------------
  ! Global continuation object
  !------------------------------------------------------------------
  type(continuation_scheduler_t), public, target :: nekotop_continuation

contains

  !> Initialize the continuation scheduler from JSON defaults
  subroutine init_scheduler(this, json)
    class(continuation_scheduler_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    integer :: n_default

    call json_get_or_default(json, &
         'optimization.solver.continuation_iterations', n_default, 1)

    this%default_iterations = n_default

    ! Initialize empty parameter array
    if (allocated(this%params)) deallocate(this%params)
  end subroutine init_scheduler

  !> free continuation scheduler
  subroutine free_scheduler(this)
    class(continuation_scheduler_t), intent(inout) :: this
    integer :: i

    if (.not. allocated(this%params)) return

    do i = 1, size(this%params)
       call this%params(i)%free()
    end do

    deallocate(this%params)
  end subroutine free_scheduler

  !> Read a parameter from JSON and optionally register it for continuation.
  !! If the parameter is given as a scalar, its value is returned.
  !! If it is given as an array, the first entry is used as the initial value
  !! and the full array is registered for continuation.
  !! If the parameter is not found and no default_value is provided,
  !! the routine terminates with an error.
  !! @param this The continuation scheduler instance.
  !! @param json The json_file object.
  !! @param name The parameter name.
  !! @param target Variable that may be registered for continuation.
  !! @param value Returned value of the parameter.
  !! @param default_value Optional default value if the parameter is not found.
  subroutine json_get_or_register(this, json, name, target, value, &
       default_value)
    class(continuation_scheduler_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: name
    real(kind=rp), target, intent(inout) :: target
    real(kind=rp), intent(out) :: value
    real(kind=rp), intent(in), optional :: default_value
    real(kind=rp), allocatable :: values(:)
    real(kind=rp) :: scalar_parameter
    logical :: found
    integer :: var_type, iter

    ! Inspect JSON for the parameter
    call json%info(name, found=found, var_type=var_type)
    if (.not. found) then
       if (present(default_value)) then
          value = default_value
       else
          call neko_error(trim(name)//" not found and no default provided")
       end if
    else if (var_type == 6) then ! scalar
       if (present(default_value)) then
          call json_get_or_default(json, name, scalar_parameter, default_value)
       else
          call json_get(json, name, scalar_parameter)
       end if
       value = scalar_parameter
    else if (var_type == 3) then ! array
       call json_get_or_lookup(json, name, values)
       value = values(1)
       ! Read optional iterations per parameter step
       call json_get_or_default(json, trim(name)//'_iterations', iter, &
            this%default_iterations)
       ! Register the parameter for continuation
       call this%register_parameter(name, target, values, iter)
    else
       call neko_error(trim(name)//" can only be real variable or real array!")
    end if

  end subroutine json_get_or_register

  !> Register a continuation parameter
  subroutine register_parameter(this, name, target, values, iterations)
    class(continuation_scheduler_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    real(rp), target, intent(inout) :: target
    real(rp), intent(in) :: values(:)
    integer, optional, intent(in) :: iterations

    integer :: n, old_size
    type(continuation_parameter_t), allocatable :: tmp(:)
    integer :: iter_val

    ! Determine iterations per value
    if (present(iterations)) then
       iter_val = iterations
    else
       iter_val = this%default_iterations
    end if

    ! Append new parameter
    if (.not. allocated(this%params)) then
       allocate(this%params(1))
       n = 1
    else
       old_size = size(this%params)
       allocate(tmp(old_size))
       tmp = this%params
       deallocate(this%params)
       allocate(this%params(old_size + 1))
       this%params(1:old_size) = tmp
       deallocate(tmp)
       n = old_size + 1
    end if

    this%params(n)%name = name
    this%params(n)%target => target
    allocate(this%params(n)%values(size(values)))
    this%params(n)%values = values
    this%params(n)%iterations_per_value = iter_val
  end subroutine register_parameter

  !> Update all registered parameters for the current iteration
  subroutine update(this, iter)
    class(continuation_scheduler_t), intent(inout) :: this
    integer, intent(in) :: iter
    integer :: i

    if (.not. allocated(this%params)) return

    do i = 1, size(this%params)
       call this%params(i)%update(iter)
    end do
  end subroutine update

  !> Get the name ith param in the scheduler
  character(len=32) function get_param_name(this, i)
    class(continuation_scheduler_t), intent(in) :: this
    integer, intent(in) :: i

    get_param_name = ''
    if (.not. allocated(this%params)) return
    if (i < 1 .or. i > size(this%params)) return
    get_param_name = this%params(i)%name
  end function get_param_name

  !> Get the number of params in the scheduler
  function get_n_params(this) result(n)
    class(continuation_scheduler_t), intent(in) :: this
    integer :: n
    if (allocated(this%params)) then
       n = size(this%params)
    else
       n = 0
    end if
  end function get_n_params

  !> Update a single parameter based on the iteration
  subroutine continuation_parameter_update(this, iter)
    class(continuation_parameter_t), intent(inout) :: this
    integer, intent(in) :: iter
    integer :: idx

    if (.not. allocated(this%values)) return
    if (this%iterations_per_value <= 0) return

    idx = (iter - 1) / this%iterations_per_value + 1
    idx = min(idx, size(this%values))

    this%target = this%values(idx)
  end subroutine continuation_parameter_update

  !> Free a single continuation parameter
  subroutine continuation_parameter_free(this)
    class(continuation_parameter_t), intent(inout) :: this
    if (allocated(this%name)) deallocate(this%name)
    if (allocated(this%values)) deallocate(this%values)
    nullify(this%target)
  end subroutine continuation_parameter_free

end module continuation_scheduler
