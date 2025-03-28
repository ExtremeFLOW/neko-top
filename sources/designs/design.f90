! Copyright (c) 2024-2025, The Neko Authors
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

!> Implements the `design_t`.
module design
  use json_module, only: json_file
  use simulation_m, only: simulation_t
  use vector, only: vector_t
  use utils, only: neko_error
  implicit none
  private

  !> An abstract design type.
  !!
  !! This type is used to define the interface for the design variables. These
  !! design are used to define the optimization problem. A given design can be
  !! initialized from the factory and is responsible for providing any
  !! alterations to the simulation based on the design variables.
  type, abstract :: design_t
     private

     !> The number of design variables
     integer :: n = 0

   contains

     ! ----------------------------------------------------------------------- !
     ! Interfaces

     !> Initialize the design.
     !! @details
     !! This method is used to initialize the design from a JSON file. The
     !! design object will be initialized based on the parameters provided in
     !! the JSON file. The simulation object is also provided to allow for any
     !! additional setup.
     !!
     !! @param this The design object.
     !! @param parameters The JSON parameters.
     !! @param simulation The simulation object.
     procedure, pass(this) :: init_from_json_sim => design_init_from_json_sim

     !> Initialize the design.
     !! @details
     !! This method is used to initialize the design from a JSON file. The
     !! design object will be initialized based on the parameters provided in
     !! the JSON file.
     !!
     !! @param this The design object.
     !! @param parameters The JSON parameters.
     procedure, pass(this) :: init_from_json => design_init_from_json

     !> Free the design.
     procedure(design_free), public, pass(this), deferred :: free

     !> Retrieve the design variables.
     procedure(design_get_design), public, pass(this), deferred :: get_design
     !> Update the design variables.
     procedure(design_update_design), public, pass(this), deferred :: &
          update_design

     !> Run the forward mapping of the design
     procedure(design_map_forward), public, pass(this), deferred :: &
          map_forward
     !> Run the backward mapping of the design
     procedure(design_map_backward), public, pass(this), deferred :: &
          map_backward
     !> Write the design
     procedure(design_write), public, pass(this), deferred :: write

     ! ----------------------------------------------------------------------- !
     ! Methods

     !> Initialize the base design
     procedure, pass(this) :: init_base => design_init_base
     !> Free the base design
     procedure, pass(this) :: free_base => design_free_base
     !> Return the number of design variables
     procedure, public, pass(this) :: size => design_size

  end type design_t

  ! ========================================================================== !
  ! Interface for the factory function

  !> Factory function for the design object.
  !! @details
  !! This function is used to create a new design object based on the provided
  !! JSON file. The simulation object is also provided to allow for any
  !! additional setup.
  !!
  !! @param object The design object.
  !! @param parameters The JSON file.
  !! @param simulation The simulation object [Optional].
  interface design_factory
     module subroutine design_factory(object, parameters, simulation)
       class(design_t), allocatable, intent(inout) :: object
       type(json_file), intent(inout) :: parameters
       type(simulation_t), intent(inout), optional :: simulation
     end subroutine design_factory
  end interface design_factory


  ! ========================================================================== !
  ! Public interface for the deferred methods

  abstract interface
     subroutine design_free(this)
       import design_t
       class(design_t), intent(inout) :: this
     end subroutine design_free

     function design_get_design(this) result(x)
       import design_t, vector_t
       class(design_t), intent(in) :: this
       type(vector_t) :: x
     end function design_get_design

     subroutine design_update_design(this, x)
       import design_t, vector_t
       class(design_t), intent(inout) :: this
       type(vector_t), intent(inout) :: x
     end subroutine design_update_design

     subroutine design_map_forward(this)
       import design_t
       class(design_t), intent(inout) :: this
     end subroutine design_map_forward

     subroutine design_map_backward(this, sensitivity)
       import design_t, vector_t
       class(design_t), intent(inout) :: this
       type(vector_t), intent(in) :: sensitivity
     end subroutine design_map_backward

     subroutine design_write(this, idx)
       import design_t
       class(design_t), intent(inout) :: this
       integer, intent(in) :: idx
     end subroutine design_write
  end interface

  public :: design_t, design_factory
contains

  !> Dummy initialization from JSON
  !! @param this The design object.
  !! @param parameters The JSON parameters.
  subroutine design_init_from_json(this, parameters)
    class(design_t), intent(inout) :: this
    type(json_file), intent(inout) :: parameters

    call neko_error("Design type does not support initialization " // &
         "without simulation")
  end subroutine design_init_from_json

  !> Dummy initialization from JSON
  !! @param this The design object.
  !! @param parameters The JSON parameters.
  !! @param simulation The simulation object.
  subroutine design_init_from_json_sim(this, parameters, simulation)
    class(design_t), intent(inout) :: this
    type(json_file), intent(inout) :: parameters
    type(simulation_t), intent(inout) :: simulation

    call neko_error("Design type does not support initialization " // &
         "with simulation")
  end subroutine design_init_from_json_sim

  !> Initialize the base design
  !! @param this The design object.
  !! @param n The number of design variables.
  subroutine design_init_base(this, n)
    class(design_t), intent(inout) :: this
    integer, intent(in) :: n
    this%n = n
  end subroutine design_init_base

  !> Free the base design
  !! @param this The design object.
  subroutine design_free_base(this)
    class(design_t), intent(inout) :: this
    this%n = 0
  end subroutine design_free_base

  !> Return the number of design variables
  !! @param this The design object.
  !! @return The number of design variables.
  pure function design_size(this) result(n)
    class(design_t), intent(in) :: this
    integer :: n
    n = this%n
  end function design_size

end module design
