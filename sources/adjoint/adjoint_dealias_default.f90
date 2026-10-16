!> @file adjoint_dealias_default.f90
!! @copyright
!! Copyright (c) 2026, The Neko-TOP Authors
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
!> Shared default for the adjoint over-integration (dealiasing) flags.
!!
!! The adjoint machinery exposes several independent `dealias` switches -- on
!! the Brinkman design, on the adjoint scalar coupling term, and on the
!! objectives' forcing and sensitivity contributions. Historically each of
!! them defaulted to `.true.` in isolation, so a case that disabled
!! over-integration through `case.numerics.dealias` still over-integrated
!! every adjoint contribution. This module provides the single fallback those
!! switches now share: the case's own `case.numerics.dealias`.
module adjoint_dealias_default
  use json_module, only: json_file
  use json_utils, only: json_get_or_default
  implicit none
  private

  public :: dealias_default

contains

  !> Default value for the adjoint over-integration flags.
  !!
  !! Returns the case's `case.numerics.dealias`, so that an adjoint switch
  !! which the case file does not set explicitly follows whatever the forward
  !! problem does.
  !!
  !! The `.true.` fallback applies only when the key is absent altogether.
  !! `adjoint_fluid_scheme_incompressible.f90` reads `case.numerics.dealias`
  !! unconditionally, and its `init` runs before the design or the
  !! objectives are constructed, so on every supported Neko-TOP
  !! configuration this branch is not expected to be exercised; it merely
  !! preserves the previous behaviour on any path which skips that check.
  !!
  !! @param case_params The root case JSON, i.e. the object holding `case`
  !! and `optimization` as its children.
  !! @return The default to use for an unset adjoint dealiasing flag.
  logical function dealias_default(case_params) result(dealias)
    type(json_file), intent(inout) :: case_params

    call json_get_or_default(case_params, 'case.numerics.dealias', &
         dealias, .true.)
  end function dealias_default

end module adjoint_dealias_default
