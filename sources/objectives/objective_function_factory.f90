! Copyright (c) 2025, The Neko Authors
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

!> Submodule for the objective function factory
submodule (objective_function) objective_function_factory_mod

  use utils, only: neko_type_error
  use volume_constraint, only: volume_constraint_t
  use minimum_dissipation_objective_function, only: minimum_dissipation_objective_function_t
  use topopt_design, only: topopt_design_t
  use simulation, only: simulation_t

  implicit none

  character(len=80), parameter :: KNOWN_TYPES(2) = [ character(len=80) :: &
       'minimum_dissipation', 'volume_constraint']


contains

  ! -------------------------------------------------------------------------- !
  ! Factory function

  !> Factory function
  !! Allocates and initializes an objective function object
  !! @param object The objective function object to be created
  !! @param type The type of the objective function
  !! @param design The design object
  !! @param simulation The simulation object
  module subroutine objective_function_factory(object, type, design, simulation)
    class(objective_function_t), allocatable, intent(inout) :: object
    character(len=*), intent(in) :: type
    type(topopt_design_t), intent(inout) :: design
    type(simulation_t), target, intent(inout) :: simulation

    select case(type)
      case('minimum_dissipation')
       allocate(minimum_dissipation_objective_function_t::object)
      case('volume_constraint')
       allocate(volume_constraint_t::object)
      case default
       call neko_type_error('objective_function', type, KNOWN_TYPES)

    end select

    call object%init(design, simulation)
  end subroutine objective_function_factory

end submodule objective_function_factory_mod
