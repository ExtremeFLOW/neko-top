!> @file design_simplefield.f90
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

! Implements the `simplefield_design_t` type.
module simplefield_design
  use num_types, only: rp, sp
  use field, only: field_t
  use json_module, only: json_file
  use mapping, only: mapping_t
  use PDE_filter, only: PDE_filter_t
  use RAMP_mapping, only: RAMP_mapping_t
  use coefs, only: coef_t
  use scratch_registry, only: neko_scratch_registry
  use point_zone_registry, only: neko_point_zone_registry
  use point_zone, only: point_zone_t
  use mask_ops, only: mask_exterior_const
  use neko_config, only: NEKO_BCKND_DEVICE
  use device, only: device_memcpy, HOST_TO_DEVICE
  use device_math, only: device_copy
  use design, only: design_t
  use math, only: rzero
  use simulation_m, only: simulation_t
  use json_module, only: json_file
  use json_utils, only: json_get
  use simple_brinkman_source_term, only: simple_brinkman_source_term_t
  use vector, only: vector_t
  use math, only: copy
  use utils, only: neko_error

  use fld_file_output, only: fld_file_output_t

  implicit none
  private

  !> A topology optimization design variable
  type, extends(design_t), public :: simplefield_design_t
     private

     type(field_t) :: designfield
     type(vector_t) :: x_coord
     type(vector_t) :: y_coord
     type(vector_t) :: z_coord

     ! needed to write the design
     type(fld_file_output_t) :: output


   contains

     ! ----------------------------------------------------------------------- !
     ! Initializations

     !> Initialize the design
     generic, public :: init => init_from_components
     !> Initialize the design from components
     procedure, pass(this) :: init_from_components => &
          design_simple_init_from_components

     !> Add mappings to the design
     procedure, pass(this) :: add_mapping => design_simple_add_mapping

     !> Retrieve the design variables
     procedure, pass(this) :: get_values => design_simple_get_values

     ! Overrides of base class deferred procedures
     procedure, pass(this) :: design_get_x => design_simple_get_x
     procedure, pass(this) :: design_get_y => design_simple_get_y
     procedure, pass(this) :: design_get_z => design_simple_get_z

     !> Update the design
     procedure, pass(this) :: update_design => design_simple_update_design

     !> map (this will include everything from mapping
     procedure, pass(this) :: map_forward => design_simple_map_forward

     !> this will contain chain rule for going backwards
     procedure, pass(this) :: map_backward => design_simple_map_backward

     !> Write the design
     procedure, pass(this) :: write => design_simple_write

     !> Destructor
     procedure, pass(this) :: free => design_simple_free

  end type simplefield_design_t

contains

  subroutine design_simple_init_from_components(this, n, x, y, z, neko_field)
    class(simplefield_design_t), intent(inout) :: this
    integer, intent(in) :: n
    type(vector_t), intent(in) :: x, y, z
    type(field_t) :: neko_field

    call this%free()
    call this%init_base('simplefield_design', n)

    this%x_coord = x
    this%y_coord = y
    this%z_coord = z
    this%designfield = neko_field
    call this%output%init(sp, 'design', 1)
    call this%output%fields%assign_to_field(1, this%designfield)
  end subroutine design_simple_init_from_components

  !> Free the design
  subroutine design_simple_free(this)
    class(simplefield_design_t), intent(inout) :: this

    call this%free_base()
    call this%x_coord%free()
    call this%y_coord%free()
    call this%z_coord%free()

    call this%designfield%free()
  end subroutine design_simple_free

  !> Add mappings to the design
  subroutine design_simple_add_mapping(this, parameters, simulation)
    class(simplefield_design_t), intent(inout) :: this
    type(json_file), intent(inout) :: parameters
    type(simulation_t), intent(inout) :: simulation

  end subroutine design_simple_add_mapping


  subroutine design_simple_map_forward(this)
    class(simplefield_design_t), intent(inout) :: this


  end subroutine design_simple_map_forward

  subroutine design_simple_get_values(this, values)
    class(simplefield_design_t), intent(in) :: this
    type(vector_t), intent(inout) :: values
    integer :: n

    n = this%size()
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(values%x_d, this%designfield%x_d, n)
    else
       call copy(values%x, this%designfield%x, n)
    end if
  end subroutine design_simple_get_values

  subroutine design_simple_get_x(this, x)
    class(simplefield_design_t), intent(in) :: this
    type(vector_t), intent(inout) :: x

    if (this%size() .ne. x%size()) then
       call neko_error('Get x: size mismatch')
    end if
    x = this%x_coord
  end subroutine design_simple_get_x

  subroutine design_simple_get_y(this, y)
    class(simplefield_design_t), intent(in) :: this
    type(vector_t), intent(inout) :: y

    if (this%size() .ne. y%size()) then
       call neko_error('Get y: size mismatch')
    end if
    y = this%y_coord
  end subroutine design_simple_get_y

  subroutine design_simple_get_z(this, z)
    class(simplefield_design_t), intent(in) :: this
    type(vector_t), intent(inout) :: z

    if (this%size() .ne. z%size()) then
       call neko_error('Get z: size mismatch')
    end if
    z = this%z_coord
  end subroutine design_simple_get_z



  subroutine design_simple_update_design(this, values)
    class(simplefield_design_t), intent(inout) :: this
    type(vector_t), intent(inout) :: values
    integer :: n

    n = this%size()
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(this%designfield%x_d, values%x_d, n)
    else
       call copy(this%designfield%x, values%x, n)
    end if
  end subroutine design_simple_update_design

  subroutine design_simple_map_backward(this, sensitivity)
    class(simplefield_design_t), intent(inout) :: this
    type(vector_t), intent(in) :: sensitivity

  end subroutine design_simple_map_backward

  subroutine design_simple_write(this, idx)
    class(simplefield_design_t), intent(inout) :: this
    integer, intent(in) :: idx
    call this%output%sample(real(idx, kind=rp))
  end subroutine design_simple_write

end module simplefield_design
