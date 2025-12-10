!> @file user_module.f90
!! @copyright
!! Copyright (c) 2024-2025, The Neko-TOP Authors
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

!> Module designed to setup the topology optimization user interface for Neko.
!! This module should initialize, finalize, and run the time dependent
!! operations of our topology optimization problem.
!!
!! This should include, but not be limited to:
!! 1. Setting the material properties (Required by Neko when using the user
!!    interface)
!! 2. Setting the initial conditions for the scalar field
!! 3. Assigning the user defined force for the fluid
!! 4. Define the adjoint and sensitivity analysis
!! 5. Define the objective function
!! 6. Define the constraints
!!
!! Currently the computations should be setup in the "user_check" subroutine.
!! This is a hacky way to run the physics and adjoint physics, but it works.
module topology_optimization_user_module
  use case, only: case_t
  use json_file_module, only: json_file
  use json_utils, only: json_get
  use num_types, only: rp
  use user_initial_conditions, only: scalar_z_split_ic

  implicit none

  private
  public :: neko_user_init

contains

  !> Assign user conditions for the neko case
  !!
  !!
  !! \param[inout] neko_case The neko case to setup the user interface for
  !!
  !! @todo
  !! We use a hacky way to run the physics and adjoint physics. This
  !! should be replaced with a more robust way to run the physics and adjoint
  !! physics.
  subroutine neko_user_init(neko_case)
    type(case_t), intent(inout) :: neko_case

    ! Set the properties for the fluid
    ! neko_case%user%initial_conditions => scalar_z_split_ic

  end subroutine neko_user_init



end module topology_optimization_user_module
