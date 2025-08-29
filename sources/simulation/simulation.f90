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
  use json_utils, only: json_extract_item, json_get_or_default
  use num_types, only: rp, sp, dp
  use logger, only: LOG_SIZE, neko_log
  use mpi_f08, only: MPI_WTIME
  use jobctrl, only: jobctrl_time_limit
  use profiler, only: profiler_start, profiler_stop
  use simulation_adjoint, only: simulation_adjoint_init, &
       simulation_adjoint_step, simulation_adjoint_finalize
  use simulation, only: simulation_init, simulation_step, simulation_finalize, &
       simulation_restart
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

     ! ----------------------------------------------------------------------- !
     ! Unsteady simulation parameters

     ! This is used to save the state of the simulation at certain time steps
     ! to be able to restart the simulation from there. This is used for the
     ! adjoint simulation.
     logical :: unsteady = .false.
     integer :: n_saves_memory = 10
     integer :: n_saves_disc = 0
     integer :: n_timesteps = 0
     integer :: first_valid_timestep = 2
     integer :: loaded_checkpoint = -1
     type(chkp_output_t) :: chkp_output
     type(field_t), dimension(:), allocatable :: p_list
     type(field_t), dimension(:), allocatable :: u_list
     type(field_t), dimension(:), allocatable :: v_list
     type(field_t), dimension(:), allocatable :: w_list
     type(field_t), dimension(:), allocatable :: s_list

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

     !> Save the current state of the simulation to disk
     procedure, pass(this) :: save_state => simulation_save_state
     !> Restore the forward simulation state
     procedure, pass(this) :: restore_state => simulation_restore_state
  end type simulation_t

  public :: simulation_t
contains

  !> Initialize the simulation
  subroutine simulation_initialize(this, parameters)
    class(simulation_t), intent(inout), target :: this
    type(json_file), intent(inout) :: parameters
    integer :: i, j, n_scalars
    class(scalar_scheme_t), pointer :: scalar_i
    character(len=10) :: str
    character(len=:), allocatable :: chkp_file_name

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


    call json_get_or_default(parameters, "unsteady.enable", &
         this%unsteady, .false.)

    if (this%unsteady) then
       ! --------------------------------------------------------------------- !
       ! Handle unsteady simulation
       ! this needs some refactoring, but I envision the json going:
       ! unsteady.strategy.type
       ! 
       ! type has options
       ! - RAM (implying we just save everything to RAM, what we have rn)
       ! - DISK (implying we save uniformly to disk)
       ! - REVOLVE (self explanitory)
       ! - POD (sneaky POD idea)
       !
       ! Based on this we can read parameters specific to the checkpointing
       ! strategy.
       !
       ! then the functions 
       ! - simulation_save_state 
       ! - simulation_restore_state
       !
       ! point to unique ways in which each of these checkpointing strategies
       ! save and restore u,v,w. Maybe we can even have a dedicated checkpoint
       ! object with proceedure save and restore, and we use a factory to
       ! build the various kinds.

       ! Read options related to check pointing
       call json_get_or_default(parameters, "unsteady.n_memory", &
            this%n_saves_memory, 10)
       call json_get_or_default(parameters, "unsteady.filename", &
            chkp_file_name, "forward_chkp_")

       call this%chkp_output%init(this%neko_case%chkp, chkp_file_name)

       ! Allocate the RAM Checkpoints
       allocate(this%p_list(this%n_saves_memory))
       allocate(this%u_list(this%n_saves_memory))
       allocate(this%v_list(this%n_saves_memory))
       allocate(this%w_list(this%n_saves_memory))
       if (this%have_scalar) then
          allocate(this%s_list(this%n_saves_memory * n_scalars))
       end if

       do i = 1, this%n_saves_memory
          write(str, '(A,I0)') "p_chkp_", i
          call this%p_list(i)%init(this%neko_case%fluid%p%dof, str)
          write(str, '(A,I0)') "u_chkp_", i
          call this%u_list(i)%init(this%neko_case%fluid%u%dof, str)
          write(str, '(A,I0)') "v_chkp_", i
          call this%v_list(i)%init(this%neko_case%fluid%v%dof, str)
          write(str, '(A,I0)') "w_chkp_", i
          call this%w_list(i)%init(this%neko_case%fluid%w%dof, str)
          if (this%have_scalar) then
             do j = 1, n_scalars
                write(str, '(A,I0,A,I0)') "s_chkp_", i, "_", j
                scalar_i => this%neko_case%scalars%scalar_fields(j)
                call this%s_list((i - 1) * n_scalars + j)%init(scalar_i%s%dof, str)
             end do
          end if
       end do
    else
       ! --------------------------------------------------------------------- !
       ! Handle steady simulation
       ! this would be a good place to append the steady simcomp
    end if

  end subroutine simulation_initialize

  !> Free the simulation
  subroutine simulation_free(this)
    class(simulation_t), intent(inout) :: this
    integer :: i, n_scalars

    ! Free the RAM Checkpoints
    if (this%unsteady) then
       do i = 1, this%n_saves_memory
          call this%p_list(i)%free()
          call this%u_list(i)%free()
          call this%v_list(i)%free()
          call this%w_list(i)%free()
       end do

       if (this%have_scalar) then
          n_scalars = size(this%s_list)
          do i = 1, n_scalars * this%n_saves_memory
             call this%s_list(i)%free()
          end do
       end if

       if (allocated(this%p_list)) deallocate(this%p_list)
       if (allocated(this%u_list)) deallocate(this%u_list)
       if (allocated(this%v_list)) deallocate(this%v_list)
       if (allocated(this%w_list)) deallocate(this%w_list)
       if (allocated(this%s_list)) deallocate(this%s_list)

       this%unsteady = .false.
       this%n_saves_memory = 0
       this%n_saves_disc = 0
       this%n_timesteps = 0
       this%loaded_checkpoint = -1
    end if

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

    call profiler_start
    loop_start = MPI_WTIME()
    do while (this%neko_case%time%t .lt. this%neko_case%time%end_time)
       call simulation_step(this%neko_case, dt_controller, loop_start)

       if (this%unsteady) call this%save_state(this%neko_case%time)
    end do
    call profiler_stop

    call simulation_finalize(this%neko_case)

  end subroutine simulation_run_forward

  !> Run the simulation
  subroutine simulation_run_backward(this)
    class(simulation_t), intent(inout) :: this
    type(time_step_controller_t) :: dt_controller
    real(kind=dp) :: loop_start
    real(kind=rp) :: cfl
    integer :: i, total_timesteps

    call dt_controller%init(this%neko_case%params)

    call simulation_adjoint_init(this%adjoint_case, dt_controller)

    call profiler_start
    cfl = this%adjoint_case%fluid_adj%compute_cfl(this%adjoint_case%time%dt)
    loop_start = MPI_WTIME()

    total_timesteps = int(this%neko_case%time%end_time / this%adjoint_case%time%dt)
    do i = total_timesteps, 1, -1
       if (this%unsteady) call this%restore_state(i)

       call simulation_adjoint_step(this%adjoint_case, dt_controller, cfl, &
            loop_start)
    end do
    call profiler_stop

    call simulation_adjoint_finalize(this%adjoint_case)

  end subroutine simulation_run_backward

  !> Reset the simulation
  subroutine simulation_reset(this)
    class(simulation_t), intent(inout) :: this
    integer :: i, n_scalars, j

    call reset(this%neko_case)

    ! TODO
    ! reset for the adjoint
    ! call reset(this%adjoint_case)
    this%adjoint_case%time%t = 0.0_rp
    this%adjoint_case%time%tstep = 0

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

    ! Reset our unsteady
    if (this%unsteady) then
       this%n_saves_disc = 0
       this%loaded_checkpoint = -1
       this%n_timesteps = 0

       do i = 1, this%n_saves_memory
          call field_rzero(this%p_list(i))
          call field_rzero(this%u_list(i))
          call field_rzero(this%v_list(i))
          call field_rzero(this%w_list(i))
          if (this%have_scalar) then
             n_scalars = size(this%s_list)
             do j = 1, n_scalars
                call field_rzero(this%s_list((i - 1) * n_scalars + j))
             end do
          end if
       end do
    end if

  end subroutine simulation_reset

  !> Save the current state of the simulation to disk
  subroutine simulation_save_state(this, time)
    class(simulation_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    ! Update the number of recorded timesteps
    this%n_timesteps = this%n_timesteps + 1

    ! Sample the checkpoint if needed
    if (modulo(time%tstep, this%n_saves_memory) .eq. 0 .or. &
         time%tstep .le. this%first_valid_timestep) then

       call this%chkp_output%set_counter(time%tstep)
       call this%chkp_output%sample(time%t)
       this%n_saves_disc = this%n_saves_disc + 1
    end if

  end subroutine simulation_save_state

  !> Restore the forward simulation state
  subroutine simulation_restore_state(this, tstep)
    class(simulation_t), intent(inout) :: this
    integer, intent(in) :: tstep
    type(time_step_controller_t) :: dt_controller
    type(file_t) :: chkp_file
    character(len=256) :: chkp_file_name
    real(kind=dp) :: loop_start
    integer :: j, k, previous_save_state, next_save_state
    integer :: i_scalars, n_scalars
    type(field_t), pointer :: u, v, w, p, s

    loop_start = MPI_WTIME()

    u => this%neko_case%fluid%u
    v => this%neko_case%fluid%v
    w => this%neko_case%fluid%w
    p => this%neko_case%fluid%p
    s => null()

    ! Determine the nearest save states
    previous_save_state = tstep - modulo(tstep, this%n_saves_memory)
    next_save_state = previous_save_state + this%n_saves_memory

    if (previous_save_state .le. this%first_valid_timestep) then
       previous_save_state = min(tstep, this%first_valid_timestep)
    end if

    if (previous_save_state .lt. this%first_valid_timestep) then
       next_save_state = previous_save_state + 1
    end if

    if (this%loaded_checkpoint .ne. previous_save_state) then
       write(chkp_file_name, '(A,I5.5,A)') 'forward_chkp_', &
            previous_save_state, '.chkp'

       call chkp_file%init(chkp_file_name)
       call chkp_file%read(this%neko_case%chkp)

       call dt_controller%init(this%neko_case%params)
       call simulation_restart(this%neko_case, this%neko_case%chkp)
       this%neko_case%time%tstep = previous_save_state
       this%loaded_checkpoint = this%neko_case%time%tstep

       do k = previous_save_state, min(next_save_state - 1, this%n_timesteps)

          if (k .ne. previous_save_state) then
             if (this%neko_case%time%t .ge. this%neko_case%time%end_time) exit
             call simulation_step(this%neko_case, dt_controller, loop_start)
          end if

          j = modulo(k, this%n_saves_memory) + 1
          call field_copy(this%p_list(j), p)
          call field_copy(this%u_list(j), u)
          call field_copy(this%v_list(j), v)
          call field_copy(this%w_list(j), w)
          if (this%have_scalar) then
             n_scalars = size(this%neko_case%scalars%scalar_fields)
             do i_scalars = 1, n_scalars
                s => this%neko_case%scalars%scalar_fields(i_scalars)%s
                call field_copy(this%s_list((j - 1) * n_scalars + i_scalars), s)
             end do
          end if
       end do

    end if

    j = modulo(tstep, this%n_saves_memory) + 1

    call field_copy(p, this%p_list(j))
    call field_copy(u, this%u_list(j))
    call field_copy(v, this%v_list(j))
    call field_copy(w, this%w_list(j))
    if (this%have_scalar) then
       n_scalars = size(this%neko_case%scalars%scalar_fields)
       do i_scalars = 1, n_scalars
          s => this%neko_case%scalars%scalar_fields(i_scalars)%s
          call field_copy(s, this%s_list((j - 1) * n_scalars + i_scalars))
       end do
    end if


  end subroutine simulation_restore_state


  !> Write current state of the simulation to disk
  subroutine simulation_write(this, idx)
    class(simulation_t), intent(inout) :: this
    integer, intent(in) :: idx

    call this%output_forward%sample(real(idx, kind=rp))
    call this%output_adjoint%sample(real(idx, kind=rp))

  end subroutine simulation_write

end module simulation_m
