!> @file neko_top_source_terms.f90
!! @copyright
!! Copyright (c) 2025-2026, The Neko-TOP Authors
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

!> Neko-TOP Source term register.
submodule(neko_top) neko_top_source_terms
  use source_term, only: source_term_allocate, register_source_term

  ! Our user-defined source terms
  use adjoint_brinkman_dissipation_source_term, only: &
       adjoint_brinkman_dissipation_source_term_allocate
  use adjoint_viscous_dissipation_source_term, only: &
       adjoint_viscous_dissipation_source_term_allocate
  use adjoint_mixing_scalar_source_term, only: &
       adjoint_mixing_scalar_source_term_allocate
  use adjoint_scalar_convection_source_term, only: &
       adjoint_scalar_convection_source_term_allocate
  use simple_brinkman_source_term, only: simple_brinkman_source_term_allocate

contains

  !> @brief Register the known source terms from Neko-TOP in the Neko system.
  module subroutine register_source_terms()
    procedure(source_term_allocate), pointer :: adjoint_brinkman_dissipation
    procedure(source_term_allocate), pointer :: adjoint_viscous_dissipation
    procedure(source_term_allocate), pointer :: adjoint_mixing_scalar
    procedure(source_term_allocate), pointer :: adjoint_scalar_convection
    procedure(source_term_allocate), pointer :: simple_brinkman

    ! Assign the pointers
    adjoint_brinkman_dissipation => &
         adjoint_brinkman_dissipation_source_term_allocate
    adjoint_viscous_dissipation => &
         adjoint_viscous_dissipation_source_term_allocate
    adjoint_mixing_scalar => adjoint_mixing_scalar_source_term_allocate
    adjoint_scalar_convection => adjoint_scalar_convection_source_term_allocate
    simple_brinkman => simple_brinkman_source_term_allocate

    ! Register the source terms
    call register_source_term('adjoint_brinkman_dissipation', &
         adjoint_brinkman_dissipation)
    call register_source_term('adjoint_viscous_dissipation', &
         adjoint_viscous_dissipation)
    call register_source_term('adjoint_mixing_scalar', adjoint_mixing_scalar)
    call register_source_term('adjoint_scalar_convection', &
         adjoint_scalar_convection)
    call register_source_term('simple_brinkman', simple_brinkman)
  end subroutine register_source_terms

end submodule neko_top_source_terms
