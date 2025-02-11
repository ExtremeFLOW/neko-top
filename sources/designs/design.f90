! Copyright (c) 2024, The Neko Authors
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

 ! Implements the `topopt_design_t` type.
module design
  use json_module, only: json_file
  use simulation, only: simulation_t
  implicit none
  private

  !> A topology optimization design variable
  type, abstract, public :: design_t
     private

     !> The number of design variables
     integer :: n = 0

   contains

     !> Initialize the design
     procedure(design_init), pass(this), public, deferred :: init

     !> Initialize the base design
     procedure, pass(this), public :: init_base => design_init_base

     !> Return the number of design variables
     procedure, pass(this), public :: size => design_size

  end type design_t

  abstract interface
     subroutine design_init(this, parameters, simulation)
       import design_t, simulation_t, json_file
       class(design_t), intent(inout) :: this
       type(json_file), intent(inout) :: parameters
       type(simulation_t), intent(inout) :: simulation
     end subroutine design_init
  end interface

contains

  !> Initialize the base design
  subroutine design_init_base(this, n)
    class(design_t), intent(inout) :: this
    integer, intent(in) :: n
    this%n = n
  end subroutine design_init_base

  !> Return the number of design variables
  pure function design_size(this) result(n)
    class(design_t), intent(in) :: this
    integer :: n
    n = this%n
  end function design_size

end module design
