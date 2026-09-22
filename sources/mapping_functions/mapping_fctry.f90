!> @file mapping_fctry.f90
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
!
!> Defines a factory subroutine for mapping functions.
submodule (mapping) mapping_fctry
  use utils, only: neko_type_error
  use linear_mapping, only : linear_mapping_t
  use PDE_filter_mapping, only: PDE_filter_t
  use RAMP_mapping, only: RAMP_mapping_t
  use Borrvall_Petersson_mapping, only: Borrvall_Petersson_mapping_t
  use heaviside_mapping, only: heaviside_mapping_t
  use SIMP_mapping, only: SIMP_mapping_t
  use json_utils, only : json_get
  use utils, only : concat_string_array, neko_error
  implicit none

  ! List of all possible types created by the factory routine
  character(len=20) :: MAPPING_KNOWN_TYPES(6) = [character(len=20) :: &
       "linear", &
       "PDE_filter", &
       "RAMP", &
       "Borrvall_Petersson", &
       "SIMP", &
       "heaviside_mapping"]

contains

  !> Mapping factory. Both constructs and initializes the object.
  !! @param json JSON object initializing the mapping.
  !! @param coef The SEM coefficients.
  module subroutine mapping_factory(object, json, coef)
    class(mapping_t), allocatable, intent(inout) :: object
    type(json_file), intent(inout) :: json
    type(coef_t), intent(inout) :: coef
    character(len=:), allocatable :: type_name

    call json_get(json, "type", type_name)

    select case (trim(type_name))
    case ("linear")
       allocate(linear_mapping_t::object)
    case ("PDE_filter")
       allocate(PDE_filter_t::object)
    case ("RAMP")
       allocate(RAMP_mapping_t::object)
    case ("Borrvall_Petersson")
       allocate(Borrvall_Petersson_mapping_t::object)
    case ("SIMP")
       allocate(SIMP_mapping_t::object)
    case ("heaviside_mapping", "heaviside_projection")
       allocate(heaviside_mapping_t::object)
    case default
       call neko_type_error("Mapping function", type_name, MAPPING_KNOWN_TYPES)
    end select

    ! Initialize
    call object%init(json, coef)

  end subroutine mapping_factory

end submodule mapping_fctry
