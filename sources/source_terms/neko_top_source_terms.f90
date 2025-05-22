! Copyright (c) 2025, The Neko Authors
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

!> Neko-TOP Source term register.
submodule(neko_top) neko_top_source_terms
  use source_term, only: source_term_t, source_term_allocate, &
       register_source_term

  ! Our user-defined source terms
  use adjoint_lube_source_term, only: adjoint_lube_source_term_t
  use adjoint_minimum_dissipation_source_term, only: &
       adjoint_minimum_dissipation_source_term_t
  use adjoint_mixing_scalar_source_term, only: &
       adjoint_mixing_scalar_source_term_t
  use adjoint_scalar_convection_source_term, only: &
       adjoint_scalar_convection_source_term_t
  use simple_brinkman_source_term, only: simple_brinkman_source_term_t

contains

  !> @brief Register the known source terms from Neko-TOP.
  module subroutine register_source_terms()
    procedure(source_term_allocate), pointer :: adjoint_lube
    procedure(source_term_allocate), pointer :: adjoint_minimum_dissipation
    procedure(source_term_allocate), pointer :: adjoint_mixing_scalar
    procedure(source_term_allocate), pointer :: adjoint_scalar_convection
    procedure(source_term_allocate), pointer :: simple_brinkman

    ! Assign the pointers
    adjoint_lube => adjoint_lube_source_term_allocate
    adjoint_minimum_dissipation => adjoint_minimum_dissipation_source_term_allocate
    adjoint_mixing_scalar => adjoint_mixing_scalar_source_term_allocate
    adjoint_scalar_convection => adjoint_scalar_convection_source_term_allocate
    simple_brinkman => simple_brinkman_source_term_allocate

    ! Register the source terms
    call register_source_term('adjoint_lube', adjoint_lube)
    call register_source_term('adjoint_minimum_dissipation', adjoint_minimum_dissipation)
    call register_source_term('adjoint_mixing_scalar', adjoint_mixing_scalar)
    call register_source_term('adjoint_scalar_convection', adjoint_scalar_convection)
    call register_source_term('simple_brinkman', simple_brinkman)
  end subroutine register_source_terms

  ! ========================================================================== !
  ! Definitions of the source term allocators

  !> Allocator for the adjoint lube source term.
  subroutine adjoint_lube_source_term_allocate(obj)
    class(source_term_t), allocatable, intent(inout) :: obj
    allocate(adjoint_lube_source_term_t::obj)
  end subroutine adjoint_lube_source_term_allocate

  !> Allocator for the adjoint minimum dissipation source term.
  subroutine adjoint_minimum_dissipation_source_term_allocate(obj)
    class(source_term_t), allocatable, intent(inout) :: obj
    allocate(adjoint_minimum_dissipation_source_term_t::obj)
  end subroutine adjoint_minimum_dissipation_source_term_allocate

  !> Allocator for the adjoint mixing scalar source term.
  subroutine adjoint_mixing_scalar_source_term_allocate(obj)
    class(source_term_t), allocatable, intent(inout) :: obj
    allocate(adjoint_mixing_scalar_source_term_t::obj)
  end subroutine adjoint_mixing_scalar_source_term_allocate

  !> Allocator for the adjoint scalar convection source term.
  subroutine adjoint_scalar_convection_source_term_allocate(obj)
    class(source_term_t), allocatable, intent(inout) :: obj
    allocate(adjoint_scalar_convection_source_term_t::obj)
  end subroutine adjoint_scalar_convection_source_term_allocate

  !> Allocator for the simple brinkman source term.
  subroutine simple_brinkman_source_term_allocate(obj)
    class(source_term_t), allocatable, intent(inout) :: obj
    allocate(simple_brinkman_source_term_t::obj)
  end subroutine simple_brinkman_source_term_allocate

end submodule neko_top_source_terms
