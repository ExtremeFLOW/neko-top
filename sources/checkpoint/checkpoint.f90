!> @file checkpoint.f90
!! @copyright
!! Copyright (c) 2025, The Neko-TOP Authors
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
module simulation_checkpoint
  use num_types, only: rp, sp, dp
  use case, only: case_t
  use json_file_module, only: json_file
  use json_utils, only: json_get_or_default
  use scalar_scheme, only: scalar_scheme_t
  use time_state, only: time_state_t
  use chkp_output, only: chkp_output_t
  use field, only: field_t
  use mpi_f08, only: MPI_WTIME
  use utils, only: neko_error
  use field_math, only: field_copy, field_rzero
  use profiler, only: profiler_start_region, profiler_end_region
  implicit none
  private

  type, public :: simulation_checkpoint_t
     private

     ! ----------------------------------------------------------------------- !
     ! User parameters

     !> Whether checkpointing is enabled
     logical :: enabled = .false.
     !> The checkpointing algorithm to use
     character(len=256) :: algorithm = "linear"
     !> The name of the checkpoint file
     character(len=256) :: filename = "checkpoint"
     !> The format of the checkpoint file
     character(len=8) :: fmt = "chkp"
     !> Number of checkpoints to keep in memory
     integer :: n_saves_memory = 10
     !> Whether to keep the checkpoint files on disk after the simulation ends
     logical :: keep_checkpoints = .true.

     ! Internal parameters
     integer :: n_saves_disc = 0
     integer :: n_timesteps = 0
     integer :: first_valid_timestep = 2
     integer :: loaded_checkpoint = -1

     ! Structures to hold the checkpoint data
     type(chkp_output_t) :: chkp_output
     type(field_t), dimension(:), allocatable :: p_list
     type(field_t), dimension(:), allocatable :: u_list
     type(field_t), dimension(:), allocatable :: v_list
     type(field_t), dimension(:), allocatable :: w_list

     integer :: n_scalars = 0
     type(field_t), dimension(:), allocatable :: s_list

   contains
     !> Initialization
     generic, public :: init => init_from_json, init_from_components
     !> Initialization from a JSON file
     procedure, public, pass(this) :: init_from_json => &
          checkpoint_init_from_json
     !> Initialization from components
     procedure, public, pass(this) :: init_from_components => &
          checkpoint_init_from_components
     !> Free
     procedure, public, pass(this) :: free => checkpoint_free
     !> Reset the checkpoint data
     procedure, public, pass(this) :: reset => checkpoint_reset
     !> Save the current state of the simulation to disk
     procedure, public, pass(this) :: save => checkpoint_save
     !> Restore the forward simulation state
     procedure, public, pass(this) :: restore => checkpoint_restore

  end type simulation_checkpoint_t

  ! ========================================================================== !
  ! Module procedures for our algorithm implementations.


  interface
     !> Save the current state of the simulation in a linear fashion
     module subroutine checkpoint_save_linear(this, neko_case)
       class(simulation_checkpoint_t), intent(inout) :: this
       class(case_t), intent(inout) :: neko_case
     end subroutine checkpoint_save_linear

     !> Restore the forward simulation state in a linear fashion
     module subroutine checkpoint_restore_linear(this, neko_case, tstep)
       class(simulation_checkpoint_t), intent(inout) :: this
       class(case_t), target, intent(inout) :: neko_case
       integer, intent(in) :: tstep
     end subroutine checkpoint_restore_linear
  end interface

contains

  ! ========================================================================== !
  ! Initialization and deallocation

  !> Initialization
  subroutine checkpoint_init_from_json(this, neko_case, params)
    class(simulation_checkpoint_t), intent(inout) :: this
    class(case_t), target, intent(inout) :: neko_case
    type(json_file), target, intent(inout) :: params
    integer :: n_saves_memory
    character(len=:), allocatable :: filename, algorithm, fmt
    logical :: enabled, keep_checkpoints

    call json_get_or_default(params, "enabled", enabled, .false.)
    if (.not. enabled) return

    call json_get_or_default(params, "algorithm", algorithm, "linear")
    call json_get_or_default(params, "n_memory", n_saves_memory, 10)
    call json_get_or_default(params, "filename", filename, "checkpoint")
    call json_get_or_default(params, "format", fmt, "chkp")
    call json_get_or_default(params, "keep_checkpoints", keep_checkpoints, &
         .true.)

    call this%init_from_components(neko_case, algorithm, n_saves_memory, &
         filename, fmt, keep_checkpoints)
  end subroutine checkpoint_init_from_json

  !> Initialization from components
  subroutine checkpoint_init_from_components(this, neko_case, algorithm, &
       n_saves_memory, filename, fmt, keep_checkpoints)
    class(simulation_checkpoint_t), intent(inout), target :: this
    class(case_t), target, intent(inout) :: neko_case
    character(len=*), optional, intent(in) :: algorithm
    integer, optional, intent(in) :: n_saves_memory
    character(len=*), optional, intent(in) :: filename
    character(len=*), optional, intent(in) :: fmt
    logical, optional, intent(in) :: keep_checkpoints

    class(scalar_scheme_t), pointer :: scalar_i
    integer :: i, j
    character(len=80) :: str

    call this%free()

    ! Set internal parameters
    this%enabled = .true.
    if (present(algorithm)) this%algorithm = algorithm
    if (present(filename)) this%filename = filename
    if (present(n_saves_memory)) this%n_saves_memory = n_saves_memory
    if (present(fmt)) this%fmt = fmt
    if (present(keep_checkpoints)) this%keep_checkpoints = keep_checkpoints

    if (allocated(neko_case%scalars)) then
       this%n_scalars = size(neko_case%scalars%scalar_fields)
    end if


    ! Initialize the Neko checkpoint output
    call this%chkp_output%init(neko_case%chkp, this%filename, fmt = this%fmt, &
         overwrite = .true.)

    ! Allocate the RAM Checkpoints
    allocate(this%p_list(this%n_saves_memory))
    allocate(this%u_list(this%n_saves_memory))
    allocate(this%v_list(this%n_saves_memory))
    allocate(this%w_list(this%n_saves_memory))
    if (this%n_scalars .gt. 0) then
       this%n_scalars = size(neko_case%scalars%scalar_fields)
       allocate(this%s_list(this%n_saves_memory * this%n_scalars))
    end if

    do i = 1, this%n_saves_memory
       write(str, '(A,I0)') "p_chkp_", i
       call this%p_list(i)%init(neko_case%fluid%p%dof, str)
       write(str, '(A,I0)') "u_chkp_", i
       call this%u_list(i)%init(neko_case%fluid%u%dof, str)
       write(str, '(A,I0)') "v_chkp_", i
       call this%v_list(i)%init(neko_case%fluid%v%dof, str)
       write(str, '(A,I0)') "w_chkp_", i
       call this%w_list(i)%init(neko_case%fluid%w%dof, str)
       if (this%n_scalars .gt. 0) then
          do j = 1, this%n_scalars
             write(str, '(A,I0,A,I0)') "s_chkp_", i, "_", j
             scalar_i => neko_case%scalars%scalar_fields(j)%scalar
             call this%s_list((i - 1) * this%n_scalars + j)%init(scalar_i%s%dof, str)
          end do
       end if
    end do

  end subroutine checkpoint_init_from_components

  !> Free
  subroutine checkpoint_free(this)
    class(simulation_checkpoint_t), intent(inout) :: this
    integer :: i
    character(len=1024) :: file_name
    logical :: exists
    integer :: stat, unit

    ! Free the RAM Checkpoints
    do i = 1, this%n_saves_memory
       if (allocated(this%p_list)) call this%p_list(i)%free()
       if (allocated(this%u_list)) call this%u_list(i)%free()
       if (allocated(this%v_list)) call this%v_list(i)%free()
       if (allocated(this%w_list)) call this%w_list(i)%free()
    end do

    if (allocated(this%s_list)) then
       do i = 1, size(this%s_list)
          call this%s_list(i)%free()
       end do
       this%n_scalars = 0
    end if

    if (allocated(this%p_list)) deallocate(this%p_list)
    if (allocated(this%u_list)) deallocate(this%u_list)
    if (allocated(this%v_list)) deallocate(this%v_list)
    if (allocated(this%w_list)) deallocate(this%w_list)
    if (allocated(this%s_list)) deallocate(this%s_list)

    ! Delete the checkpoint file list
    if (.not. this%keep_checkpoints) then
       do i = this%n_timesteps, 1, -1
          call this%chkp_output%set_counter(i)
          file_name = this%chkp_output%file_%get_fname()
          inquire(file = trim(file_name), exist = exists)
          if (exists) then
             open(newunit = unit, file = trim(file_name), iostat = stat, &
                  status='old')
             if (stat .eq. 0) close(unit, status = 'delete')
          end if
       end do
    end if

    ! Reset to default values
    this%enabled = .false.
    this%filename = "checkpoint"
    this%fmt = "chkp"
    this%algorithm = "linear"
    this%n_saves_memory = 10
    this%keep_checkpoints = .true.

    this%n_saves_disc = 0
    this%n_timesteps = 0
    this%first_valid_timestep = 2
    this%loaded_checkpoint = -1

  end subroutine checkpoint_free

  ! ========================================================================== !
  ! Saving and Restoring

  !> Save the current state of the simulation to disk
  subroutine checkpoint_save(this, neko_case)
    class(simulation_checkpoint_t), intent(inout) :: this
    class(case_t), intent(inout) :: neko_case

    if (.not. this%enabled) return

    call profiler_start_region("Checkpoint save")

    ! Update the number of recorded timesteps
    this%n_timesteps = this%n_timesteps + 1

    select case (this%algorithm)
    case ("linear")
       call checkpoint_save_linear(this, neko_case)
    case default
       call neko_error("Unknown checkpoint algorithm: " // this%algorithm)
    end select

    call profiler_end_region("Checkpoint save")
  end subroutine checkpoint_save

  !> Restore the forward simulation state
  subroutine checkpoint_restore(this, neko_case, tstep)
    class(simulation_checkpoint_t), intent(inout) :: this
    class(case_t), target, intent(inout) :: neko_case
    integer, intent(in) :: tstep
    character(len=256) :: msg

    if (.not. this%enabled) return

    call profiler_start_region("Checkpoint restore")

    if (tstep .lt. 1 .or. tstep .gt. this%n_timesteps) then
       write(msg, '(A,I0,A,I0,A)') "Requested timestep ", tstep, &
            " is out of range [1, ", this%n_timesteps, "]"
       call neko_error(trim(msg))
    end if

    select case (this%algorithm)
    case ("linear")
       call checkpoint_restore_linear(this, neko_case, tstep)
    case default
       call neko_error("Unknown checkpoint algorithm: " // this%algorithm)
    end select

    call profiler_end_region("Checkpoint restore")
  end subroutine checkpoint_restore

  ! ========================================================================== !
  ! Meta handling

  !> Reset the checkpoint data
  subroutine checkpoint_reset(this)
    class(simulation_checkpoint_t), intent(inout) :: this
    integer :: i

    if (.not. this%enabled) return

    ! Reset our checkpoints
    this%loaded_checkpoint = -1
    this%n_saves_disc = 0
    this%n_timesteps = 0

    do i = 1, this%n_saves_memory
       call field_rzero(this%p_list(i))
       call field_rzero(this%u_list(i))
       call field_rzero(this%v_list(i))
       call field_rzero(this%w_list(i))
    end do

    if (allocated(this%s_list)) then
       do i = 1, size(this%s_list)
          call field_rzero(this%s_list(i))
       end do
    end if
  end subroutine checkpoint_reset

end module simulation_checkpoint
