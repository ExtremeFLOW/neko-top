!> @file simulation.f90
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
!
!> Implements the `steady_problem_t` type.
! Here, we simply march forward to steady state solutions
module simulation_m
  use case, only: case_t
  use user_access_singleton, only: neko_user_access
  use adjoint_case, only: adjoint_case_t
  use fluid_scheme_incompressible, only: fluid_scheme_incompressible_t
  use adjoint_fluid_scheme, only: adjoint_fluid_scheme_t
  use adjoint_fluid_pnpn, only: adjoint_fluid_pnpn_t
  use adjoint_scalars, only: adjoint_scalars_t
  use scalars, only: scalars_t
  use fluid_pnpn, only: fluid_pnpn_t
  use time_step_controller, only: time_step_controller_t
  use field_output, only: field_output_t
  use simcomp_executor, only: neko_simcomps
  use neko_ext, only: reset, reset_adjoint
  use utils, only: neko_error
  use json_file_module, only: json_file
  use json_utils, only: json_get, json_get_or_default
  use num_types, only: rp, sp, dp
  use mpi_f08, only: MPI_WTIME
  use profiler, only: profiler_start, profiler_stop, &
       profiler_start_region, profiler_end_region
  use simulation_adjoint, only: simulation_adjoint_init, &
       simulation_adjoint_step, simulation_adjoint_finalize
  use simulation, only: simulation_init, simulation_step, simulation_finalize, &
       simulation_restart
  use state_recover, only: state_recover_t
  use state_recover_fctry, only: state_recover_factory
  use runtime_stats, only: neko_rt_stats
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
     type(field_output_t), public :: output_forward
     !> An output sampler for the adjoint problem.
     !! This should probably be an output controller at some point instead.
     type(field_output_t), public :: output_adjoint
     !> Base output filename for Neko's own forward field output
     !! (`neko_case%f_out`), captured before any design-iteration tag is
     !! spliced into it, so that tag can be rebuilt fresh on every reset.
     character(len=:), allocatable :: forward_field_base_fname
     !> Base output filename for Neko's own adjoint field output
     !! (`adjoint_case%f_out`), captured before any design-iteration tag is
     !! spliced into it, so that tag can be rebuilt fresh on every reset.
     character(len=:), allocatable :: adjoint_field_base_fname
     !> The current design iteration. Used so that `reset` can keep Neko's
     !! own field output (`neko_case%f_out` / `adjoint_case%f_out`) from
     !! overwriting the previous design iteration's output.
     integer :: current_design_iteration = 0
     !> Whether the simulation is steady or unsteady
     logical :: unsteady = .false.

     logical :: have_scalar = .false.
     integer :: n_timesteps = 0

     ! ----------------------------------------------------------------------- !
     ! State recovery system

     !> The state recovery system data
     class(state_recover_t), allocatable :: state_recover

   contains
     !> Initialize the simulation
     procedure, pass(this) :: init => simulation_initialize
     !> Free the simulation
     procedure, pass(this) :: free => simulation_free
     !> Run the simulation
     procedure, pass(this) :: run_forward => simulation_run_forward
     !> Run the adjoint simulation
     procedure, pass(this) :: run_backward => simulation_run_backward
     !> Reset the simulation
     procedure, pass(this) :: reset => simulation_reset
     !> Set simulation output counters.
     procedure, pass(this) :: set_output_counter => &
          simulation_set_output_counter
     !> Set the current design iteration.
     procedure, pass(this) :: set_design_iteration => &
          simulation_set_design_iteration
     !> Write current state of the simulation to disk
     procedure, pass(this) :: write => simulation_write
     !> Write current state of the forward simulation to disk
     procedure, pass(this) :: write_forward => simulation_write_forward
     !> Write current state of the adjoint simulation to disk
     procedure, pass(this) :: write_adjoint => simulation_write_adjoint

  end type simulation_t
  public :: simulation_t
contains

  !> Initialize the simulation
  subroutine simulation_initialize(this, parameters)
    class(simulation_t), intent(inout), target :: this
    type(json_file), intent(inout) :: parameters
    type(json_file) :: state_recovery_params
    integer :: i, n_scalars
    character(len=:), allocatable :: output_directory, precision_s, file_format
    integer :: precision
    logical :: unsteady, subdivide

    ! initialize the primal Neko objects
    call this%neko_case%init(parameters)
    call neko_user_access%init(this%neko_case)

    call neko_rt_stats%init(parameters)
    call neko_simcomps%init(this%neko_case)

    ! initialize the adjoint
    call this%adjoint_case%init(this%neko_case)

    ! Capture Neko's own field output filenames before any design-iteration
    ! tag ever gets spliced into them (see `reset`/`neko_ext::reset`).
    this%forward_field_base_fname = &
         trim(this%neko_case%f_out%file_%get_base_fname())
    this%adjoint_field_base_fname = &
         trim(this%adjoint_case%f_out%file_%get_base_fname())

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

    !---------------------------------------------------------
    ! Initialize the output types

    ! Read settings for the output types, with some reasonable defaults
    call json_get_or_default(parameters, 'case.output_directory', &
         output_directory, '')
    call json_get_or_default(parameters, 'case.output_precision', precision_s, &
         'single')
    call json_get_or_default(parameters, 'case.fluid.output_format', &
         file_format, 'fld')
    call json_get_or_default(parameters, 'case.fluid.output_subdivide', &
         subdivide, .false.)

    if (trim(precision_s) .eq. 'double') then
       precision = dp
    else
       precision = sp
    end if

    ! Allocate the output type
    n_scalars = 0
    if (allocated(this%neko_case%scalars)) then
       n_scalars = size(this%neko_case%scalars%scalar_fields)
    end if
    call this%output_forward%init('forward_fields', 4 + n_scalars, &
         precision = precision, &
         path = trim(output_directory), &
         format = trim(file_format))
    call this%output_forward%file_%set_subdivide(subdivide)

    call this%output_forward%fields%assign(1, this%fluid%p)
    call this%output_forward%fields%assign(2, this%fluid%u)
    call this%output_forward%fields%assign(3, this%fluid%v)
    call this%output_forward%fields%assign(4, this%fluid%w)

    ! Assign all scalar fields
    if (allocated(this%neko_case%scalars)) then
       do i = 1, n_scalars
          call this%output_forward%fields%assign(4 + i, &
               this%scalars%scalar_fields(i)%scalar%s)
       end do
    end if

    ! Read settings for the output types, with some reasonable defaults
    call json_get_or_default(parameters, &
         'case.adjoint_fluid.output_format', file_format, 'fld')
    call json_get_or_default(parameters, &
         'case.adjoint_fluid.output_subdivide', subdivide, .false.)

    n_scalars = 0
    if (allocated(this%adjoint_case%adjoint_scalars)) then
       n_scalars = size(this%adjoint_case%adjoint_scalars%adjoint_scalar_fields)
    end if
    call this%output_adjoint%init('adjoint_fields', 4 + n_scalars, &
         precision = precision, &
         path = trim(output_directory), &
         format = trim(file_format))
    call this%output_adjoint%file_%set_subdivide(subdivide)

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

    ! Check if the simulation is steady or unsteady
    call json_get_or_default(parameters, "unsteady", unsteady, .false.)
    this%unsteady = unsteady

    ! State recovery is only needed for unsteady runs.
    if (this%unsteady) then
       call json_get(parameters, 'state_recovery', state_recovery_params)
       call state_recover_factory(this%state_recover, this%neko_case, &
            state_recovery_params)
    end if


  end subroutine simulation_initialize

  !> Free the simulation
  subroutine simulation_free(this)
    class(simulation_t), intent(inout) :: this

    ! Stop the profiler
    call profiler_stop

    if (allocated(this%state_recover)) then
       call this%state_recover%free()
       deallocate(this%state_recover)
    end if

    ! Free the objects
    call this%neko_case%free()
    call this%adjoint_case%free()
    call this%output_forward%free()
    call this%output_adjoint%free()

    ! Nullify pointers
    nullify(this%fluid)
    nullify(this%scalars)
    nullify(this%adjoint_fluid)
    nullify(this%adjoint_scalars)

    ! Reset flags and counters
    this%unsteady = .false.
    this%have_scalar = .false.
    this%n_timesteps = 0
    this%current_design_iteration = 0
    if (allocated(this%forward_field_base_fname)) then
       deallocate(this%forward_field_base_fname)
    end if
    if (allocated(this%adjoint_field_base_fname)) then
       deallocate(this%adjoint_field_base_fname)
    end if

    ! Close global objects
    call neko_simcomps%free()

  end subroutine simulation_free

  !> Run the simulation
  subroutine simulation_run_forward(this)
    class(simulation_t), intent(inout) :: this
    type(time_step_controller_t) :: dt_controller
    real(kind=dp) :: loop_start

    call dt_controller%init(this%neko_case%params)

    call this%neko_case%time%reset()
    call simulation_init(this%neko_case, dt_controller)

    if (this%unsteady) then
       if (.not. allocated(this%state_recover)) then
          call neko_error("State recovery not initialized.")
       end if
    end if

    call profiler_start_region("Forward simulation")
    loop_start = MPI_WTIME()
    this%n_timesteps = 0
    do while (.not. this%neko_case%time%is_done())
       this%n_timesteps = this%n_timesteps + 1

       call simulation_step(this%neko_case, dt_controller, loop_start)

       if (this%unsteady) then
          call this%state_recover%save()
       end if
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

    if (this%unsteady) then
       if (.not. allocated(this%state_recover)) then
          call neko_error("State recovery not initialized.")
       end if
    end if

    call profiler_start_region("Adjoint simulation")
    cfl = this%adjoint_case%fluid_adj%compute_cfl(this%adjoint_case%time%dt)
    loop_start = MPI_WTIME()
    do i = this%n_timesteps, 1, -1
       if (this%unsteady) then
          call this%state_recover%restore(i)
       end if

       call simulation_adjoint_step(this%adjoint_case, dt_controller, cfl, &
            loop_start)
    end do
    call profiler_end_region("Adjoint simulation")

    call simulation_adjoint_finalize(this%adjoint_case)

  end subroutine simulation_run_backward

  !> Reset the simulation
  subroutine simulation_reset(this)
    class(simulation_t), intent(inout) :: this

    call reset(this%neko_case, this%current_design_iteration, &
         this%forward_field_base_fname)
    call reset_adjoint(this%adjoint_case, this%neko_case, &
         this%current_design_iteration, this%adjoint_field_base_fname)
    if (this%unsteady) then
       if (.not. allocated(this%state_recover)) then
          call neko_error("State recovery not initialized.")
       end if
       call this%state_recover%reset()
    end if

  end subroutine simulation_reset

  !> Set the current design iteration, so that the next `reset` tags Neko's
  !! own field output (forward and adjoint) with it instead of overwriting
  !! the previous design iteration's output.
  !! @param this The simulation object.
  !! @param iteration The current design iteration.
  subroutine simulation_set_design_iteration(this, iteration)
    class(simulation_t), intent(inout) :: this
    integer, intent(in) :: iteration

    this%current_design_iteration = iteration

  end subroutine simulation_set_design_iteration

  subroutine simulation_set_output_counter(this, idx)
    class(simulation_t), intent(inout) :: this
    integer, intent(in) :: idx

    call this%output_forward%set_counter(idx)
    call this%output_adjoint%set_counter(idx)

  end subroutine simulation_set_output_counter

  !> Write current state of the simulation to disk
  subroutine simulation_write(this, idx)
    class(simulation_t), intent(inout) :: this
    integer, intent(in) :: idx

    call this%output_forward%sample(real(idx, kind=rp))
    call this%output_adjoint%sample(real(idx, kind=rp))

  end subroutine simulation_write

  !> Write current state of the forward simulation to disk
  subroutine simulation_write_forward(this, idx)
    class(simulation_t), intent(inout) :: this
    integer, intent(in) :: idx

    call this%output_forward%sample(real(idx, kind=rp))

  end subroutine simulation_write_forward

  !> Write current state of the adjoint simulation to disk
  subroutine simulation_write_adjoint(this, idx)
    class(simulation_t), intent(inout) :: this
    integer, intent(in) :: idx

    call this%output_adjoint%sample(real(idx, kind=rp))

  end subroutine simulation_write_adjoint

end module simulation_m
