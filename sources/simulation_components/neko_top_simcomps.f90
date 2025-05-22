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

!> Neko-TOP Simulation component register.
submodule(neko_top) neko_top_simcomps
  use simulation_component, only: simulation_component_t, &
       simulation_component_allocate, register_simulation_component

  ! Our user-defined simulation components
  use steady_simcomp, only: steady_simcomp_t

contains

  !> @brief Register the known simulation components from Neko-TOP.
  module subroutine register_simcomps
    procedure(simulation_component_allocate), pointer :: steady

    ! Assign the pointers
    steady => steady_simcomp_allocate

    ! Register the simulation components
    call register_simulation_component('steady', steady)
  end subroutine register_simcomps

  ! ========================================================================== !
  ! Definitions of the simulation component allocators

  !> Allocator for the steady simulation component.
  subroutine steady_simcomp_allocate(obj)
    class(simulation_component_t), allocatable, intent(inout) :: obj
    allocate(steady_simcomp_t::obj)
  end subroutine steady_simcomp_allocate

end submodule neko_top_simcomps
