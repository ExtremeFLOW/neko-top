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

!> Submodule for the design function factory
submodule (design) design_factory_mod
  use json_utils, only: json_get
  use utils, only: neko_type_error

  ! Import the design function types
  use brinkman_design, only: brinkman_design_t
  use simple_design, only: simple_design_t

  implicit none

  !> Known function types
  character(len=25), parameter :: KNOWN_TYPES(2) = [ character(len=25) :: &
       "brinkman", &
       "simple"]

contains

  ! -------------------------------------------------------------------------- !
  ! Factory function

  !> Factory function
  !! Allocates and initializes an design function object
  !! @param object The design function object to be created
  !! @param type The type of the design function
  !! @param design The design object
  !! @param simulation The simulation object
  module subroutine design_factory_w_simulation(object, parameters, simulation)
    class(design_t), allocatable, intent(inout) :: object
    type(json_file), intent(inout) :: parameters
    type(simulation_t), intent(inout) :: simulation
    character(len=:), allocatable :: type

    if (allocated(object)) then
       call object%free()
       deallocate(object)
    end if

    call json_get(parameters, "optimization.design.type", type)
    select case (trim(type))
    case ("brinkman")
       allocate(brinkman_design_t::object)
    case ("simple")
       allocate(simple_design_t::object)

    case default
       call neko_type_error("design", type, KNOWN_TYPES)
    end select

    call object%init_from_json(parameters, simulation)
  end subroutine design_factory_w_simulation

  !> Factory function
  !! Allocates and initializes an design function object
  !! @param object The design function object to be created
  !! @param type The type of the design function
  !! @param design The design object
  !! @param simulation The simulation object
  module subroutine design_factory_wo_simulation(object, parameters)
    class(design_t), allocatable, intent(inout) :: object
    type(json_file), intent(inout) :: parameters
    character(len=:), allocatable :: type

    if (allocated(object)) then
       call object%free()
       deallocate(object)
    end if

    call json_get(parameters, "optimization.design.type", type)
    select case (trim(type))
    case ("brinkman")
       allocate(brinkman_design_t::object)
    case ("simple")
       allocate(simple_design_t::object)

    case default
       call neko_type_error("design", type, KNOWN_TYPES)
    end select

    call object%init_from_json(parameters)
  end subroutine design_factory_wo_simulation

end submodule design_factory_mod
