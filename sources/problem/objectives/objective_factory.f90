!> @file objective_factory.f90
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

!> Submodule for the objective function factory
submodule (objective) objective_factory_mod
  use json_utils, only: json_get
  use utils, only: neko_type_error

  ! Import the objective function types
  use viscous_dissipation_objective, only: viscous_dissipation_objective_t
  use brinkman_dissipation_objective, only: brinkman_dissipation_objective_t
  use scalar_mixing_objective, only: scalar_mixing_objective_t

  implicit none

  !> Known function types
  character(len=25), parameter :: KNOWN_TYPES(3) = [ character(len=25) :: &
       "viscous_dissipation", &
       "scalar_mixing", &
       "brinkman_dissipation"]

contains

  ! -------------------------------------------------------------------------- !
  ! Factory function

  !> Factory function
  module subroutine objective_factory(object, json, design, simulation)
    class(objective_t), allocatable, intent(inout) :: object
    type(json_file), intent(inout) :: json
    class(design_t), intent(in) :: design
    type(simulation_t), target, optional, intent(inout) :: simulation
    character(len=:), allocatable :: type

    if (allocated(object)) then
       call object%free()
       deallocate(object)
    end if

    call json_get(json, "type", type)
    select case (trim(type))
    case ("viscous_dissipation")
       allocate(viscous_dissipation_objective_t::object)
    case ("scalar_mixing")
       allocate(scalar_mixing_objective_t::object)
    case ("brinkman_dissipation")
       allocate(brinkman_dissipation_objective_t::object)

    case default
       call neko_type_error("Objective", type, KNOWN_TYPES)
    end select

    call object%init(json, design, simulation)
  end subroutine objective_factory

end submodule objective_factory_mod
