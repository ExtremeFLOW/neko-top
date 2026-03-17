!> @file design.f90
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

!> Implements the `design_t`.
module design
  use json_module, only: json_file
  use simulation_m, only: simulation_t
  use num_types, only: rp
  use vector, only: vector_t
  use utils, only: neko_error, filename_suffix
  use point_zone, only: point_zone_t
  use comm, only: neko_comm
  use mpi_f08, only: MPI_Allreduce, MPI_INTEGER, MPI_SUM
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
     !> The name of the design
     character(len=:), allocatable :: name

     !> The number of design variables
     integer :: n = 0
     !> The global number of design variables
     integer :: n_global = 0

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
     procedure(design_get_values), public, pass(this), deferred :: get_values
     !> Retrieve the x location of the design variables.
     generic :: get_x => design_get_x
     generic :: x => design_get_x_i
     !> Retrieve the y location of the design variables.
     generic :: get_y => design_get_y
     generic :: y => design_get_y_i
     !> Retrieve the z location of the design variables.
     generic :: get_z => design_get_z
     generic :: z => design_get_z_i

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

     !> Save the design to a checkpoint file
     procedure, public, pass(this) :: save_checkpoint => design_save_checkpoint
     !> Load the design from a checkpoint file
     procedure, public, pass(this) :: load_checkpoint => design_load_checkpoint
     !> Set the design output counter
     procedure, public, pass(this) :: set_output_counter => &
          design_set_output_counter
     ! ----------------------------------------------------------------------- !
     ! Methods

     !> Initialize the base design
     procedure, pass(this) :: init_base => design_init_base
     !> Free the base design
     procedure, pass(this) :: free_base => design_free_base

     !> Get the name of the design.
     procedure, public, pass(this) :: get_name => design_get_name
     !> Return the number of design variables
     procedure, public, pass(this) :: size => design_size
     !> Return the number of global design variables
     procedure, public, pass(this) :: size_global => design_size_global

     !> Getter of x vector
     procedure, pass(this) :: design_get_x
     !> Getter of i'th element of x
     procedure, pass(this) :: design_get_x_i
     !> Getter of y vector
     procedure, pass(this) :: design_get_y
     !> Getter of i'th element of y
     procedure, pass(this) :: design_get_y_i
     !> Getter of z vector
     procedure, pass(this) :: design_get_z
     !> Getter of i'th element of z
     procedure, pass(this) :: design_get_z_i

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

     subroutine design_get_values(this, values)
       import design_t, vector_t
       class(design_t), intent(in) :: this
       type(vector_t), intent(inout) :: values
     end subroutine design_get_values

     subroutine design_update_design(this, values)
       import design_t, vector_t
       class(design_t), intent(inout) :: this
       type(vector_t), intent(inout) :: values
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

     subroutine set_output_counter(this, idx)
       import design_t
       class(design_t), intent(inout) :: this
       integer, intent(in) :: idx
     end subroutine set_output_counter
  end interface

  ! ========================================================================== !
  ! Module subroutine implementations

  interface
     module subroutine design_save_checkpoint_hdf5(this, filename, overwrite)
       class(design_t), intent(in) :: this
       character(len=*), intent(in) :: filename
       logical, intent(in), optional :: overwrite
     end subroutine design_save_checkpoint_hdf5

     module subroutine design_load_checkpoint_hdf5(this, filename)
       class(design_t), intent(inout) :: this
       character(len=*), intent(in) :: filename
     end subroutine design_load_checkpoint_hdf5
  end interface

  public :: design_t, design_factory

contains

  ! ========================================================================== !
  ! Initializers and destructors

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
  !! @param name The name of the design.
  !! @param n The number of design variables.
  subroutine design_init_base(this, name, n)
    class(design_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    integer :: ierr

    this%name = name
    this%n = n
    call MPI_Allreduce(n, this%n_global, 1, MPI_INTEGER, MPI_SUM, &
         neko_comm, ierr)
  end subroutine design_init_base

  !> Free the base design
  !! @param this The design object.
  subroutine design_free_base(this)
    class(design_t), intent(inout) :: this
    this%name = ""
    this%n = 0
    this%n_global = 0
  end subroutine design_free_base

  ! ========================================================================== !
  ! IO methods

  !> Save the design to a checkpoint file
  !! @param this The design object.
  !! @param filename The filename to save the checkpoint to.
  !! @param overwrite Whether to overwrite the file if it exists.
  subroutine design_save_checkpoint(this, filename, overwrite)
    class(design_t), intent(in) :: this
    character(len=*), intent(in) :: filename
    logical, intent(in), optional :: overwrite
    character(len=12) :: file_ext

    ! Determine the file extension
    call filename_suffix(filename, file_ext)

    select case (trim(file_ext))
    case ('h5', 'hdf5', 'hf5')
       call design_save_checkpoint_hdf5(this, filename, overwrite)
    case default
       call neko_error('design_save_checkpoint: Unsupported file format: ' // &
            trim(file_ext))
    end select

  end subroutine design_save_checkpoint

  !> Load the design from a checkpoint file
  !! @param this The design object.
  !! @param filename The filename to load the checkpoint from.
  subroutine design_load_checkpoint(this, filename)
    class(design_t), intent(inout) :: this
    character(len=*), intent(in) :: filename
    character(len=12) :: file_ext

    ! Determine the file extension
    call filename_suffix(filename, file_ext)

    select case (trim(file_ext))
    case ('h5', 'hdf5', 'hf5')
       call design_load_checkpoint_hdf5(this, filename)
    case default
       call neko_error('design_load_checkpoint: Unsupported file format: ' // &
            trim(file_ext))
    end select

  end subroutine design_load_checkpoint

  !> Set design output counters.
  !! Default implementation is a no-op.
  subroutine design_set_output_counter(this, idx)
    class(design_t), intent(inout) :: this
    integer, intent(in) :: idx

  end subroutine design_set_output_counter

  ! ========================================================================== !
  ! Getter methods

  !> Get the name of the design.
  !! @param this The design object.
  !! @return The name of the design.
  function design_get_name(this) result(name)
    class(design_t), intent(in) :: this
    character(len=:), allocatable :: name
    name = this%name
  end function design_get_name

  !> Return the number of design variables
  !! @param this The design object.
  !! @return The number of design variables.
  pure function design_size(this) result(n)
    class(design_t), intent(in) :: this
    integer :: n
    n = this%n
  end function design_size

  !> Return the number of global design variables
  pure function design_size_global(this) result(n)
    class(design_t), intent(in) :: this
    integer :: n
    n = this%n_global
  end function design_size_global

  subroutine design_get_x(this, x)
    class(design_t), intent(in) :: this
    type(vector_t), intent(inout) :: x
    call neko_error("Design type does not support x retrieval")
  end subroutine design_get_x

  function design_get_x_i(this, i) result(x_i)
    class(design_t), intent(in) :: this
    integer, intent(in) :: i
    real(kind=rp) :: x_i
    x_i = -huge(x_i)
    call neko_error("Design type does not support x retrieval")
  end function design_get_x_i

  subroutine design_get_y(this, y)
    class(design_t), intent(in) :: this
    type(vector_t), intent(inout) :: y
    call neko_error("Design type does not support y retrieval")
  end subroutine design_get_y

  function design_get_y_i(this, i) result(y_i)
    class(design_t), intent(in) :: this
    integer, intent(in) :: i
    real(kind=rp) :: y_i
    y_i = -huge(y_i)
    call neko_error("Design type does not support y retrieval")
  end function design_get_y_i

  subroutine design_get_z(this, z)
    class(design_t), intent(in) :: this
    type(vector_t), intent(inout) :: z
    call neko_error("Design type does not support z retrieval")
  end subroutine design_get_z

  function design_get_z_i(this, i) result(z_i)
    class(design_t), intent(in) :: this
    integer, intent(in) :: i
    real(kind=rp) :: z_i
    z_i = -huge(z_i)
    call neko_error("Design type does not support z retrieval")
  end function design_get_z_i

  ! ========================================================================= !
  ! Dummy implementations for module procedures

#if !HAVE_HDF5
  module subroutine design_save_checkpoint_hdf5(this, filename, overwrite)
    class(design_t), intent(in) :: this
    character(len=*), intent(in) :: filename
    logical, intent(in), optional :: overwrite
    call neko_error('design: HDF5 support not enabled rebuild with HAVE_HDF5')
  end subroutine design_save_checkpoint_hdf5

  module subroutine design_load_checkpoint_hdf5(this, filename)
    class(design_t), intent(inout) :: this
    character(len=*), intent(in) :: filename
    call neko_error('design: HDF5 support not enabled rebuild with HAVE_HDF5')
  end subroutine design_load_checkpoint_hdf5
#endif

end module design
