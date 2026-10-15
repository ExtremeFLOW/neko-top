!> @file initial_conditions.f90
!! @copyright
!! Copyright (c) 2023-2025, The Neko-TOP Authors
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

module user_initial_conditions
  use field, only: field_t
  use json_file_module, only: json_file
  use json_utils, only: json_get_or_default
  use device, only: HOST_TO_DEVICE
  use num_types, only: rp
  implicit none

contains

  !> Set the initial condition for the scalar field
  !! @details This function will initialize the scalar field with a two part
  !! uniform value. Above z=z0 the scalar field will be 0.0 and below z=z0 the
  !! scalar field will be 1.0.
  !! z0 is read from the JSON file under the key
  !! 'case.scalar.initial_condition.value' or set to 0.0 if not found.
  !!
  !! @param[inout] s Scalar field
  !! @param[in] split_value
  !! @param[in] value_below Value assigned below the split value
  !! @param[in] value_above Value assigned above the split value
  subroutine scalar_z_split_ic(s, split_value, value_below, value_above)
    type(field_t), intent(inout) :: s
    real(kind=rp), intent(in) :: split_value
    real(kind=rp), intent(in) :: value_below, value_above

    real(kind=rp) :: z_value
    integer :: i

    do i = 1, s%dof%size()
       z_value = s%dof%z%x(i, 1, 1, 1)

       if (z_value .gt. split_value) then
          s%x(i, 1, 1, 1) = value_above
       else
          s%x(i, 1, 1, 1) = value_below
       end if

    end do

    call s%copy_from(HOST_TO_DEVICE, .true.)

  end subroutine scalar_z_split_ic

end module user_initial_conditions
