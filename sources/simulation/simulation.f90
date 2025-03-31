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
module simulation_m
  use case, only: case_t
  use neko, only: neko_init, neko_finalize, neko_solve
  use adjoint_case, only: adjoint_case_t, adjoint_init, adjoint_free
  use fluid_scheme_incompressible, only: fluid_scheme_incompressible_t
  use adjoint_fluid_scheme, only: adjoint_fluid_scheme_t
  use adjoint_fluid_pnpn, only: adjoint_fluid_pnpn_t
  use scalar_pnpn, only: scalar_pnpn_t
  use adjoint_scalar_pnpn, only: adjoint_scalar_pnpn_t
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
  use user_intf, only: user_t, simulation_component_user_settings
  implicit none
  private

  type :: simulation_t
     !> and primal case
     type(case_t), public :: neko_case
     !> and adjoint case
     type(adjoint_case_t), public :: adjoint_case
     !> The fluid
     class(fluid_scheme_incompressible_t), public, pointer :: fluid => null()
     !> The scalar
     type(scalar_pnpn_t), public, pointer :: scalar => null()
     !> The adjoint fluid
     class(adjoint_fluid_scheme_t), public, pointer :: adjoint_fluid => null()
     !> The adjoint scalar
     type(adjoint_scalar_pnpn_t), public, pointer :: adjoint_scalar => null()
     !> An output sampler for the forward problem.
     !! This should probably be an output controller at some point instead.
     type(fld_file_output_t), public :: output_forward
     !> An output sampler for the adjoint problem.
     !! This should probably be an output controller at some point instead.
     type(fld_file_output_t), public :: output_adjoint

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

  subroutine user_simcomp(params)
    type(json_file), intent(inout) :: params
    type(steady_simcomp_t), allocatable :: steady_comp
    type(json_file) :: simcomp_settings

    ! Allocate a simulation component
    allocate(steady_comp)
    simcomp_settings = simulation_component_user_settings("steady", params)
    call neko_simcomps%add_user_simcomp(steady_comp, simcomp_settings)

  end subroutine user_simcomp

  !> Initialize the simulation
  subroutine simulation_init(this, parameters)
    class(simulation_t), intent(inout), target :: this
    type(json_file), intent(inout) :: parameters

    this%neko_case%usr%init_user_simcomp => user_simcomp

    ! initialize the primal
    call neko_init(this%neko_case)
    ! initialize the adjoint
    call adjoint_init(this%adjoint_case, this%neko_case)

    select type (fluid => this%neko_case%fluid)
    type is (fluid_pnpn_t)
       this%fluid => fluid

    end select

    select type (adjoint_fluid => this%adjoint_case%fluid_adj)
    type is (adjoint_fluid_pnpn_t)
       this%adjoint_fluid => adjoint_fluid

    end select

    if (allocated(this%neko_case%scalar)) then
       this%scalar => this%neko_case%scalar
    end if

    if (allocated(this%adjoint_case%scalar_adj)) then
       this%adjoint_scalar => this%adjoint_case%scalar_adj
    end if

    ! init the sampler
    !---------------------------------------------------------
    ! Allocate the output type
    if (allocated(this%neko_case%scalar)) then
       call this%output_forward%init(sp, 'forward_fields', 5)
       call this%output_forward%fields%assign(5, this%scalar%s)
    else
       call this%output_forward%init(sp, 'forward_fields', 4)
    end if

    call this%output_forward%fields%assign(1, this%fluid%p)
    call this%output_forward%fields%assign(2, this%fluid%u)
    call this%output_forward%fields%assign(3, this%fluid%v)
    call this%output_forward%fields%assign(4, this%fluid%w)

    if (allocated(this%adjoint_case%scalar_adj)) then
       call this%output_adjoint%init(sp, 'adjoint_fields', 5)
       call this%output_adjoint%fields%assign(5, this%adjoint_scalar%s_adj)
    else
       call this%output_adjoint%init(sp, 'adjoint_fields', 4)
    end if
    call this%output_adjoint%fields%assign(1, this%adjoint_fluid%p_adj)
    call this%output_adjoint%fields%assign(2, this%adjoint_fluid%u_adj)
    call this%output_adjoint%fields%assign(3, this%adjoint_fluid%v_adj)
    call this%output_adjoint%fields%assign(4, this%adjoint_fluid%w_adj)

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
    call field_rzero(this%adjoint_case%fluid_adj%u_adj)
    call field_rzero(this%adjoint_case%fluid_adj%v_adj)
    call field_rzero(this%adjoint_case%fluid_adj%w_adj)
    if (allocated(this%neko_case%scalar)) then
       call field_rzero(this%adjoint_case%scalar_adj%s_adj)
    end if

  end subroutine simulation_reset

  !> Write current state of the simulation to disk
  subroutine simulation_write(this, idx)
    class(simulation_t), intent(inout) :: this
    integer, intent(in) :: idx

    call this%output_forward%sample(real(idx, kind=rp))
    call this%output_adjoint%sample(real(idx, kind=rp))

  end subroutine simulation_write

end module simulation_m
