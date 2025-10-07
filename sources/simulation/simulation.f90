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
  use adjoint_scalars, only: adjoint_scalars_t
  use scalars, only: scalars_t
  use scalar_scheme, only: scalar_scheme_t
  use fluid_pnpn, only: fluid_pnpn_t
  use time_step_controller, only: time_step_controller_t
  use time_state, only: time_state_t
  use fld_file_output, only: fld_file_output_t
  use chkp_output, only: chkp_output_t
  use simcomp_executor, only: neko_simcomps
  use neko_ext, only: reset
  use field, only: field_t
  use field_registry, only: neko_field_registry
  use field_math, only: field_rzero, field_copy
  use checkpoint, only: chkp_t
  use file, only: file_t
  use utils, only: neko_warning, neko_error
  use comm, only: pe_rank
  use json_file_module, only: json_file
  use json_utils, only: json_get, json_get_or_default
  use num_types, only: rp, sp, dp
  use logger, only: LOG_SIZE, neko_log
  use mpi_f08, only: MPI_WTIME
  use jobctrl, only: jobctrl_time_limit
  use profiler, only: profiler_start, profiler_stop, &
       profiler_start_region, profiler_end_region
  use simulation_adjoint, only: simulation_adjoint_init, &
       simulation_adjoint_step, simulation_adjoint_finalize
  use simulation, only: simulation_init, simulation_step, simulation_finalize, &
       simulation_restart
  use simulation_checkpoint, only: simulation_checkpoint_t
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
     type(scalars_t), public, pointer :: scalars => null()
     !> The adjoint fluid
     class(adjoint_fluid_scheme_t), public, pointer :: adjoint_fluid => null()
     !> The adjoint scalar
     type(adjoint_scalars_t), public, pointer :: adjoint_scalars => null()
     !> An output sampler for the forward problem.
     !! This should probably be an output controller at some point instead.
     type(fld_file_output_t), public :: output_forward
     !> An output sampler for the adjoint problem.
     !! This should probably be an output controller at some point instead.
     type(fld_file_output_t), public :: output_adjoint

     logical :: have_scalar = .false.
     integer :: n_timesteps = 0

     ! ----------------------------------------------------------------------- !
     ! Checkpoint system

     !> The checkpoint system data
     type(simulation_checkpoint_t) :: checkpoint

   contains
     !> Initialize the simulation
     procedure, pass(this) :: init => simulation_initialize
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
  subroutine simulation_initialize(this, parameters)
    class(simulation_t), intent(inout), target :: this
    type(json_file), intent(inout) :: parameters
    type(json_file) :: checkpoint_params
    integer :: i, n_scalars

    ! initialize the primal
    call neko_init(this%neko_case)
    ! initialize the adjoint
    call adjoint_init(this%adjoint_case, this%neko_case)

    ! Start the profiler
    call profiler_start

    select type (fluid => this%neko_case%fluid)
    type is (fluid_pnpn_t)
       this%fluid => fluid
    end select

    select type (adjoint_fluid => this%adjoint_case%fluid_adj)
    type is (adjoint_fluid_pnpn_t)
       this%adjoint_fluid => adjoint_fluid
    end select

    if (allocated(this%neko_case%scalars)) then
       this%scalars => this%neko_case%scalars
    end if

    if (allocated(this%adjoint_case%adjoint_scalars)) then
       this%adjoint_scalars => this%adjoint_case%adjoint_scalars
    end if

    ! init the sampler
    !---------------------------------------------------------
    ! Allocate the output type
    n_scalars = 0
    if (allocated(this%neko_case%scalars)) then
       n_scalars = size(this%neko_case%scalars%scalar_fields)
    end if
    call this%output_forward%init(sp, 'forward_fields', 4 + n_scalars)

    call this%output_forward%fields%assign(1, this%fluid%p)
    call this%output_forward%fields%assign(2, this%fluid%u)
    call this%output_forward%fields%assign(3, this%fluid%v)
    call this%output_forward%fields%assign(4, this%fluid%w)

    ! Assign all scalar fields
    if (allocated(this%neko_case%scalars)) then
       do i = 1, n_scalars
          call this%output_forward%fields%assign(4 + i, &
               this%scalars%scalar_fields(i)%s)
       end do
    end if

    n_scalars = 0
    if (allocated(this%adjoint_case%adjoint_scalars)) then
       n_scalars = size(this%adjoint_case%adjoint_scalars%adjoint_scalar_fields)
    end if
    call this%output_adjoint%init(sp, 'adjoint_fields', 4 + n_scalars)
    call this%output_adjoint%fields%assign(1, this%adjoint_fluid%p_adj)
    call this%output_adjoint%fields%assign(2, this%adjoint_fluid%u_adj)
    call this%output_adjoint%fields%assign(3, this%adjoint_fluid%v_adj)
    call this%output_adjoint%fields%assign(4, this%adjoint_fluid%w_adj)

    ! Assign all scalar fields
    if (allocated(this%adjoint_case%adjoint_scalars)) then
       do i = 1, n_scalars
          call this%output_adjoint%fields%assign(4 + i, &
               this%adjoint_scalars%adjoint_scalar_fields(i)%s_adj)
       end do
    end if

    if ("checkpoints" .in. parameters) then
       call json_get(parameters, 'checkpoints', checkpoint_params)
       call this%checkpoint%init(this%neko_case, checkpoint_params)
    end if

  end subroutine simulation_initialize

  !> Free the simulation
  subroutine simulation_free(this)
    class(simulation_t), intent(inout) :: this

    ! Stop the profiler
    call profiler_stop

    call this%checkpoint%free()
    call adjoint_free(this%adjoint_case)
    call neko_finalize(this%neko_case)

  end subroutine simulation_free

  !> Run the simulation
  subroutine simulation_run_forward(this)
    class(simulation_t), intent(inout) :: this
    type(time_step_controller_t) :: dt_controller
    real(kind=dp) :: loop_start

    call dt_controller%init(this%neko_case%params)

    call simulation_init(this%neko_case, dt_controller)

    loop_start = MPI_WTIME()
    call profiler_start_region("Forward simulation")
    do while (this%neko_case%time%t .lt. this%neko_case%time%end_time)
       this%n_timesteps = this%n_timesteps + 1

       call simulation_step(this%neko_case, dt_controller, loop_start)

       call this%checkpoint%save(this%neko_case)
    end do
    call profiler_end_region("Forward simulation")

    call simulation_finalize(this%neko_case)

  end subroutine simulation_run_forward

  !> Run the simulation
  subroutine simulation_run_backward(this)
    class(simulation_t), intent(inout) :: this
    type(time_step_controller_t) :: dt_controller
    real(kind=dp) :: loop_start
    real(kind=rp) :: cfl
    integer :: i

    call dt_controller%init(this%neko_case%params)

    call simulation_adjoint_init(this%adjoint_case, dt_controller)

    cfl = this%adjoint_case%fluid_adj%compute_cfl(this%adjoint_case%time%dt)
    loop_start = MPI_WTIME()

    call profiler_start_region("Adjoint simulation")
    do i = this%n_timesteps, 1, -1
       call this%checkpoint%restore(this%neko_case, i)

       call simulation_adjoint_step(this%adjoint_case, dt_controller, cfl, &
            loop_start)
    end do
    call profiler_end_region("Adjoint simulation")

    call simulation_adjoint_finalize(this%adjoint_case)

  end subroutine simulation_run_backward

  !> Reset the simulation
  subroutine simulation_reset(this)
    class(simulation_t), intent(inout) :: this
    integer :: i, n_scalars

    call reset(this%neko_case)

    ! TODO
    ! reset for the adjoint
    ! call reset(this%adjoint_case)
    this%adjoint_case%time%t = 0.0_rp
    this%adjoint_case%time%tstep = 0
    this%n_timesteps = 0

    call field_rzero(this%adjoint_case%fluid_adj%u_adj)
    call field_rzero(this%adjoint_case%fluid_adj%v_adj)
    call field_rzero(this%adjoint_case%fluid_adj%w_adj)
    n_scalars = 0
    if (allocated(this%adjoint_case%adjoint_scalars)) then
       n_scalars = size(this%adjoint_case%adjoint_scalars%adjoint_scalar_fields)
       do i = 1, n_scalars
          call field_rzero(&
               this%adjoint_case%adjoint_scalars%adjoint_scalar_fields(i)%s_adj)
       end do
    end if

    call this%checkpoint%reset()
  end subroutine simulation_reset

  !> Write current state of the simulation to disk
  subroutine simulation_write(this, idx)
    class(simulation_t), intent(inout) :: this
    integer, intent(in) :: idx

    call this%output_forward%sample(real(idx, kind=rp))
    call this%output_adjoint%sample(real(idx, kind=rp))

  end subroutine simulation_write

end module simulation_m
