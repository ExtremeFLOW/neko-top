!> @file optimizer_factory.f90
!! @copyright (c) 2025, The Neko-TOP Authors
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

submodule(optimizer) optimizer_factory_mod
  use json_utils, only: json_get
  use utils, only: neko_type_error
  use mma_optimizer, only: mma_optimizer_t


  implicit none

  !> Known function types
  character(len=25), parameter :: KNOWN_TYPES(1) = [ character(len=25) :: &
       "mma"]

contains

  ! -------------------------------------------------------------------------- !
  ! Interfaces for the factory functions

  !> @brief Factory function for the optimizer.
  !! @details This function creates an optimizer object based on the type
  !! specified in the JSON file under "optimizer.type".
  !!
  !! @param object The optimizer object to be created.
  !! @param parameters The JSON file containing the optimizer parameters.
  !! @param problem The problem object.
  !! @param design The design object.
  !! @param simulation The simulation object.
  module subroutine optimizer_factory(object, parameters, problem, design, &
       simulation)
    class(optimizer_t), allocatable, intent(inout) :: object
    type(json_file), intent(inout) :: parameters
    class(problem_t), intent(inout) :: problem
    class(design_t), intent(in) :: design
    type(simulation_t), optional, intent(in) :: simulation

    character(len=:), allocatable :: type

    if (allocated(object)) then
       call object%free()
       deallocate(object)
    end if

    ! Get the type of the optimizer
    call json_get(parameters, "optimization.solver.type", type)

    ! Select the optimizer type
    select case (trim(type))
    case ("mma")
       allocate(mma_optimizer_t::object)
    case default
       call neko_type_error("Optimizer", type, KNOWN_TYPES)
    end select

    call object%init_from_json(parameters, problem, design, simulation)

  end subroutine optimizer_factory


end submodule optimizer_factory_mod
