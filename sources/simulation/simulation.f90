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
  use case, only: case_t
  use neko, only: neko_init, neko_finalize, neko_solve
  use adjoint_case, only: adjoint_case_t, adjoint_init, adjoint_free
  use fluid_scheme_incompressible, only: fluid_scheme_incompressible_t
  use fluid_pnpn, only: fluid_pnpn_t
  use simulation_adjoint, only: solve_adjoint
  use fld_file_output, only: fld_file_output_t
  use steady_simcomp, only: steady_simcomp_t
  use simcomp_executor, only: neko_simcomps
  use neko_ext, only: reset
  use field_math, only: field_rzero
  use json_file_module, only: json_file
  use json_utils, only: json_extract_item
  use num_types, only: rp, sp
  implicit none
  private

  type :: simulation_t
     !> and primal case
     type(case_t), public :: neko_case
     !> and adjoint case
     type(adjoint_case_t), public :: adjoint_case

     !> The fluid
     class(fluid_scheme_incompressible_t), public, pointer :: &
          fluid_scheme => null()

     !> An output sampler for the problem. This should probably be an output
     !! controller at some point instead.
     type(fld_file_output_t), public :: output

   contains
     !> Initialize the simulation
     procedure, pass(this) :: init => simulation_init
     !> Free the simulation
     procedure, pass(this) :: free => simulation_free
     !> Run the simulation
     procedure, pass(this) :: run_forward => simulation_run_forward
     !> Run the simulation
     procedure, pass(this) :: run_backward => simulation_run_backward
     !> Reset the simulation
     procedure, pass(this) :: reset => simulation_reset
     !> Write current state of the simulation to disk
     procedure, pass(this) :: write => simulation_write
  end type simulation_t

  public :: simulation_t
contains

  !> Initialize the simulation
  subroutine simulation_init(this, parameters)
    class(simulation_t), intent(inout) :: this
    type(json_file), intent(inout) :: parameters
    type(steady_simcomp_t), allocatable :: steady_comp
    type(json_file) :: simcomp_settings

    ! initialize the primal
    call neko_init(this%neko_case)
    ! initialize the adjoint
    call adjoint_init(this%adjoint_case, this%neko_case)

    select type (fluid => this%neko_case%fluid)
    type is (fluid_pnpn_t)
       this%fluid_scheme => fluid

    end select

    !> Initialize the steady state simulation component
    allocate(steady_comp)
    call json_extract_item(parameters, &
         "case.simulation_components", 1, simcomp_settings)

    call steady_comp%init(simcomp_settings, this%neko_case)

    call neko_simcomps%add_user_simcomp(steady_comp)

    ! init the sampler
    !---------------------------------------------------------
    ! TODO
    ! obviously when we do the mappings properly, to many coefficients, we'll
    ! also have to modify this
    ! for now:
    ! - forward (p,u,v,w)                      1,2,3,4           p,vx,vy,vz
    ! - adjoint (p,u,v,w)                      5,6,7,8           t,s1,s2,s3

    ! Allocate the output type
    call this%output%init(sp, 'optimization', 8)
    call this%output%fields%assign(1, this%fluid_scheme%p)
    call this%output%fields%assign(2, this%fluid_scheme%u)
    call this%output%fields%assign(3, this%fluid_scheme%v)
    call this%output%fields%assign(4, this%fluid_scheme%w)
    call this%output%fields%assign(5, this%adjoint_case%scheme%p_adj)
    call this%output%fields%assign(6, this%adjoint_case%scheme%u_adj)
    call this%output%fields%assign(7, this%adjoint_case%scheme%v_adj)
    call this%output%fields%assign(8, this%adjoint_case%scheme%w_adj)

  end subroutine simulation_init

  !> Free the simulation
  subroutine simulation_free(this)
    class(simulation_t), intent(inout) :: this

    call adjoint_free(this%adjoint_case)
    call neko_finalize(this%neko_case)

  end subroutine simulation_free

  !> Run the simulation
  subroutine simulation_run_forward(this)
    class(simulation_t), intent(inout) :: this

    ! run the primal
    call neko_solve(this%neko_case)

  end subroutine simulation_run_forward

  !> Run the simulation
  subroutine simulation_run_backward(this)
    class(simulation_t), intent(inout) :: this

    ! run the adjoint
    call solve_adjoint(this%adjoint_case)

  end subroutine simulation_run_backward

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

  !> Write current state of the simulation to disk
  subroutine simulation_write(this, idx)
    class(simulation_t), intent(inout) :: this
    integer, intent(in) :: idx

    call this%output%sample(real(idx, kind=rp))

  end subroutine simulation_write

end module simulation
