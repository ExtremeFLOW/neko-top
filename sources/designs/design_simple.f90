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

! Implements the `simple_design_t` type.
module simple_design
  use num_types, only: rp, sp
  use field, only: field_t
  use json_module, only: json_file
  use mapping, only: mapping_t
  use PDE_filter, only: PDE_filter_t
  use RAMP_mapping, only: RAMP_mapping_t
  use coefs, only: coef_t
  use scratch_registry, only: neko_scratch_registry
  use fld_file_output, only: fld_file_output_t
  use point_zone_registry, only: neko_point_zone_registry
  use point_zone, only: point_zone_t
  use mask_ops, only: mask_exterior_const
  use neko_config, only: NEKO_BCKND_DEVICE
  use device, only: device_memcpy, HOST_TO_DEVICE
  use design, only: design_t
  use math, only: rzero
  use simulation, only: simulation_t
  use json_module, only: json_file
  use json_utils, only: json_get
  use simple_brinkman_source_term, only: simple_brinkman_source_term_t
  use vector, only: vector_t
  use math, only: copy
  use field_registry, only: neko_field_registry
  implicit none
  private

  !> A topology optimization design variable
  type, extends(design_t), public :: simple_design_t
     private

     type(vector_t) :: x

   contains

     ! ----------------------------------------------------------------------- !
     ! Initializations

     !> Initialize the design
     generic, public :: init => init_from_json, init_from_components
     !> Initialize the design from a JSON file
     procedure, pass(this), public :: init_from_json => &
          design_simple_init_from_json
     !> Initialize the design from components
     procedure, pass(this), public :: init_from_components => &
          design_simple_init_from_components

     !> Add mappings to the design
     procedure, pass(this) :: add_mapping => design_simple_add_mapping

     !> Retrieve the design variables
     procedure, pass(this) :: get_design => design_simple_get_design

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

  end type simple_design_t

contains

  !> Initialize the design from a JSON file
  subroutine design_simple_init_from_json(this, parameters, simulation)
    class(simple_design_t), intent(inout) :: this
    type(json_file), intent(inout) :: parameters
    type(simulation_t), intent(inout), optional :: simulation
    integer :: n

    call json_get(parameters, 'design.n', n)

    call this%init_from_components(n, simulation)

  end subroutine design_simple_init_from_json

  !> Free the design
  subroutine design_simple_free(this)
    class(simple_design_t), intent(inout) :: this

    call this%free_base()
    call this%x%free()

  end subroutine design_simple_free

  subroutine design_simple_init_from_components(this, n, simulation)
    class(simple_design_t), intent(inout) :: this
    integer, intent(in) :: n
    type(simulation_t), intent(inout) :: simulation

    call this%init_base(n)
    call this%x%init(n)
    this%x = 0.0_rp

  end subroutine design_simple_init_from_components

  !> Add mappings to the design
  subroutine design_simple_add_mapping(this, parameters, simulation)
    class(simple_design_t), intent(inout) :: this
    type(json_file), intent(inout) :: parameters
    type(simulation_t), intent(inout) :: simulation

  end subroutine design_simple_add_mapping


  subroutine design_simple_map_forward(this)
    class(simple_design_t), intent(inout) :: this


  end subroutine design_simple_map_forward

  function design_simple_get_design(this) result(x)
    class(simple_design_t), intent(in) :: this
    type(vector_t) :: x

    x = this%x
  end function design_simple_get_design

  subroutine design_simple_update_design(this, x)
    class(simple_design_t), intent(inout) :: this
    type(vector_t), intent(inout) :: x

    this%x = x

  end subroutine design_simple_update_design

  subroutine design_simple_map_backward(this, sensitivity)
    class(simple_design_t), intent(inout) :: this
    type(vector_t), intent(in) :: sensitivity

  end subroutine design_simple_map_backward

  subroutine design_simple_write(this, idx)
    class(simple_design_t), intent(inout) :: this
    integer, intent(in) :: idx

  end subroutine design_simple_write

end module simple_design
