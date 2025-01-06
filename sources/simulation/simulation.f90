! Copyright (c) 2023, The Neko Authors
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
!
!> Implements the `steady_problem_t` type.
! Here, we simply march forward to steady state solutions
module simulation
  use num_types, only: rp, sp
  use problem, only: problem_t
  use adjoint_case, only: adjoint_case_t, adjoint_init, adjoint_free
  use simulation_adjoint, only: solve_adjoint
  use neko, only: neko_init, neko_finalize, neko_solve
  use case, only: case_t
  use design, only: design_t
  use topopt_design, only: topopt_design_t
  use json_file_module, only: json_file
  use fld_file_output, only: fld_file_output_t
  use steady_simcomp, only: steady_simcomp_t
  use simple_brinkman_source_term, only: simple_brinkman_source_term_t
  use json_utils_ext, only: json_key_fallback, json_get_subdict
  use json_value_module, only: json_value
  use utils, only: neko_error
  use user_intf, only: simulation_component_user_settings
  use simcomp_executor, only: neko_simcomps
  use neko_ext, only: reset
  use field_math, only: field_rzero
  use json_utils, only: json_extract_item
  implicit none
  private

  type, public :: simulation_t
     !> and primal case
     type(case_t), public :: neko_case
     !> and adjoint case
     type(adjoint_case_t), public :: adjoint_case

   contains
     !> Initialize the simulation
     procedure, pass(this) :: init => simulation_init
     !> Free the simulation
     procedure, pass(this) :: free => simulation_free
     !> Run the simulation
     procedure, pass(this) :: run => simulation_run
     !> Reset the simulation
     procedure, pass(this) :: reset => simulation_reset
  end type simulation_t

contains

  !> Initialize the simulation
  subroutine simulation_init(this)
    class(simulation_t), intent(inout) :: this
    type(steady_simcomp_t), allocatable :: steady_comp
    type(json_file) :: simcomp_settings

    ! initialize the primal
    call neko_init(this%neko_case)
    ! initialize the adjoint
    call adjoint_init(this%adjoint_case, this%neko_case)

    !> Initialize the steady state simulation component
    allocate(steady_comp)
    call json_extract_item(this%neko_case%params, &
         "case.simulation_components", 1, simcomp_settings)

    call steady_comp%init(simcomp_settings, this%neko_case)

    call neko_simcomps%add_user_simcomp(steady_comp)

  end subroutine simulation_init

  !> Free the simulation
  subroutine simulation_free(this)
    class(simulation_t), intent(inout) :: this

    call adjoint_free(this%adjoint_case)
    call neko_finalize(this%neko_case)

  end subroutine simulation_free

  !> Run the simulation
  subroutine simulation_run(this)
    class(simulation_t), intent(inout) :: this

    ! run the primal
    call neko_solve(this%neko_case)
    ! run the adjoint
    call solve_adjoint(this%adjoint_case)

  end subroutine simulation_run

  !> Reset the simulation
  subroutine simulation_reset(this)
    class(simulation_t), intent(inout) :: this

    call reset(this%neko_case)

    ! TODO
    ! reset for the adjoint
    ! call reset(this%adjoint_case)
    call field_rzero(this%adjoint_case%scheme%u_adj)
    call field_rzero(this%adjoint_case%scheme%v_adj)
    call field_rzero(this%adjoint_case%scheme%w_adj)

  end subroutine simulation_reset

end module simulation
