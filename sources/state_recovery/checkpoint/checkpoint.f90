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
!> @brief Checkpoint-based state recovery for adjoint runs.
module simulation_checkpoint
  use num_types, only: rp
  use case, only: case_t
  use json_file_module, only: json_file
  use json_utils, only: json_get, json_get_or_default
  use time_state, only: time_state_t
  use chkp_output, only: chkp_output_t
  use field, only: field_t
  use field_list, only: field_list_t
  use logger, only: neko_log, LOG_SIZE, NEKO_LOG_DEBUG
  use utils, only: neko_error
  use math, only: copy, rzero
  use profiler, only: profiler_start_region, profiler_end_region
  use state_recover, only: state_recover_t
  use comm, only: pe_rank, NEKO_COMM
  use neko_config, only: NEKO_BCKND_DEVICE
  use device, only: device_memcpy, DEVICE_TO_HOST, HOST_TO_DEVICE
  use registry, only: neko_registry
  use mpi_f08, only: MPI_WTIME, MPI_Barrier
  implicit none
  private

  type, public, extends(state_recover_t) :: simulation_checkpoint_t
     private

     ! ----------------------------------------------------------------------- !
     ! User parameters

     !> Whether checkpointing is enabled
     logical :: enabled = .false.
     !> The checkpointing algorithm to use
     character(len=256) :: algorithm = "linear"
     !> The name of the checkpoint file
     character(len=256) :: filename = "forward_checkpoint"
     !> The path to the checkpoint file (directory)
     character(len=256) :: path = "checkpoints/"
     !> The format of the checkpoint file
     character(len=8) :: fmt = "chkp"
     !> Number of checkpoints to keep in memory
     integer :: n_saves_memory = 10
     !> Whether to keep the checkpoint files on disk after the simulation ends
     logical :: keep_checkpoints = .false.

     ! Internal parameters
     integer :: n_saves_disc = 0
     integer :: first_valid_timestep = 2
     integer :: loaded_checkpoint = -1

     ! Field pointers
     type(field_list_t) :: state_list
     type(host_array), dimension(:,:), allocatable :: state_storage

     ! Structures to hold the checkpoint data
     type(chkp_output_t) :: chkp_output

   contains
     !> Initialization from a JSON file
     procedure, public, pass(this) :: init => checkpoint_init_from_json
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

     !> Save current data to the ram checkpoint at index
     procedure, pass(this) :: save_data => checkpoint_save_data
     !> Restore data from the ram checkpoint at index to the current state
     procedure, pass(this) :: load_data => checkpoint_load_data
  end type simulation_checkpoint_t

  type :: host_array
     real(kind=rp), allocatable :: data(:)
     integer :: size = 0
   contains
     procedure, pass(this) :: init => host_array_init
     procedure, pass(this) :: free => host_array_free
     procedure, pass(this) :: is_allocated => host_array_is_allocated
  end type host_array

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
  !! @param[inout] this Checkpointing implementation.
  !! @param[inout] neko_case Case data structure.
  !! @param[inout] params JSON parameters.
  subroutine checkpoint_init_from_json(this, neko_case, params)
    class(simulation_checkpoint_t), intent(inout) :: this
    class(case_t), target, intent(inout) :: neko_case
    type(json_file), target, intent(inout) :: params
    integer :: n_saves_memory
    character(len=:), allocatable :: path, filename, algorithm, fmt
    character(len=256), dimension(:), allocatable :: extra_field_names
    type(field_list_t) :: extra_fields
    type(field_t), pointer :: fi
    integer :: i
    logical :: enabled, keep_checkpoints

    call json_get_or_default(params, "enabled", enabled, .false.)
    if (.not. enabled) return

    call json_get_or_default(params, "algorithm", algorithm, "linear")
    call json_get_or_default(params, "n_memory", n_saves_memory, 10)
    call json_get_or_default(params, "path", path, "checkpoints/")
    call json_get_or_default(params, "filename", filename, "checkpoint")
    call json_get_or_default(params, "format", fmt, "chkp")
    call json_get_or_default(params, "keep_checkpoints", keep_checkpoints, &
         .false.)

    select case (trim(algorithm))
    case ("linear", "LINEAR", "Linear")
    case default
       call neko_error("Only the linear checkpoint strategy is supported.")
    end select

    if ("extra_fields" .in. params) then
       allocate(extra_field_names(0))
       call json_get(params, "extra_fields", extra_field_names)
       call extra_fields%init(size(extra_field_names))
       do i = 1, size(extra_field_names)
          fi => neko_registry%get_field(extra_field_names(i))
          call extra_fields%assign(i, fi)
       end do
       ! Create a field list for the extra fields
       call this%init_from_components(neko_case, algorithm, n_saves_memory, &
            path, filename, fmt, keep_checkpoints, extra_fields)
    else
       ! Create a field list without the extra fields
       call this%init_from_components(neko_case, algorithm, n_saves_memory, &
            path, filename, fmt, keep_checkpoints)
    end if

  end subroutine checkpoint_init_from_json

  !> Initialization from components
  !! @param[inout] this Checkpointing implementation.
  !! @param[inout] neko_case Case data structure.
  !! @param[in] algorithm Checkpointing algorithm identifier.
  !! @param[in] n_saves_memory Number of checkpoints in memory.
  !! @param[in] path Output directory for checkpoint files.
  !! @param[in] filename Checkpoint base filename.
  !! @param[in] fmt Checkpoint file format.
  !! @param[in] keep_checkpoints Whether to keep checkpoint files on disk.
  !! @param[inout] extra_fields Additional fields to include in checkpoints.
  subroutine checkpoint_init_from_components(this, neko_case, algorithm, &
       n_saves_memory, path, filename, fmt, keep_checkpoints, extra_fields)
    class(simulation_checkpoint_t), intent(inout), target :: this
    class(case_t), target, intent(inout) :: neko_case
    character(len=*), optional, intent(in) :: algorithm
    integer, optional, intent(in) :: n_saves_memory
    character(len=*), optional, intent(in) :: path
    character(len=*), optional, intent(in) :: filename
    character(len=*), optional, intent(in) :: fmt
    logical, optional, intent(in) :: keep_checkpoints
    type(field_list_t), optional, intent(inout) :: extra_fields
    type(field_t), pointer :: si
    character(len=LOG_SIZE) :: msg
    integer :: i, n_states
    logical :: exists

    call this%free()

    ! Assign parameters from arguments or defaults
    this%enabled = .true.
    if (present(algorithm)) this%algorithm = algorithm
    if (present(n_saves_memory)) this%n_saves_memory = n_saves_memory
    if (present(path)) this%path = trim(path)
    if (present(filename)) this%filename = trim(filename)
    if (present(fmt)) this%fmt = trim(fmt)
    if (present(keep_checkpoints)) this%keep_checkpoints = keep_checkpoints

    inquire(file = trim(this%path), exist = exists)
    if (.not. exists) then
       call MPI_Barrier(NEKO_COMM)
       if (pe_rank .eq. 0) then
          call execute_command_line("mkdir -p '" // trim(this%path) // "'")
       end if
       call MPI_Barrier(NEKO_COMM)
    end if

    ! Initialize the Neko checkpoint output
    call this%chkp_output%init(neko_case%chkp, this%filename, &
         fmt = this%fmt, path = this%path, overwrite = .true.)

    n_states = 4
    if (allocated(neko_case%scalars)) then
       n_states = n_states + size(neko_case%scalars%scalar_fields)
    end if
    if (present(extra_fields)) then
       n_states = n_states + extra_fields%size()
    end if

    call this%state_list%init(n_states)

    ! Assign fluid pointers
    call this%state_list%assign(1, neko_case%fluid%p)
    call this%state_list%assign(2, neko_case%fluid%u)
    call this%state_list%assign(3, neko_case%fluid%v)
    call this%state_list%assign(4, neko_case%fluid%w)
    n_states = 4

    ! Assign scalar pointers
    if (allocated(neko_case%scalars)) then
       do i = 1, size(neko_case%scalars%scalar_fields)
          si => neko_case%scalars%scalar_fields(i)%scalar%s
          call this%state_list%assign(n_states + i, si)
       end do
       n_states = n_states + size(neko_case%scalars%scalar_fields)
    end if

    ! Assign any extra fields specified by the user
    if (present(extra_fields)) then
       do i = 1, extra_fields%size()
          si => extra_fields%get_by_index(i)
          call this%state_list%assign(n_states + i, si)
       end do
       n_states = n_states + extra_fields%size()
    end if

    ! Allocate the storage for the RAM checkpoints
    allocate(this%state_storage(this%n_saves_memory, this%state_list%size()))

    ! Write a status message with the parameters set
    call neko_log%section("Checkpointing")

    write(msg, '(A, A)') "Algorithm:                    ", trim(this%algorithm)
    call neko_log%message(trim(msg))
    write(msg, '(A,I0)') "Number of checkpoints in RAM: ", this%n_saves_memory
    call neko_log%message(trim(msg))
    write(msg, '(A, A)') "Checkpoint file path:         ", trim(this%path)
    call neko_log%message(trim(msg))
    write(msg, '(A, A)') "Checkpoint file name:         ", trim(this%filename)
    call neko_log%message(trim(msg))
    write(msg, '(A, A)') "Checkpoint file format:       ", trim(this%fmt)
    call neko_log%message(trim(msg))

    if (.not. this%keep_checkpoints) then
       call neko_log%message("Checkpoint files will be deleted.")
    else
       call neko_log%message("Checkpoint files will be kept.")
    end if

    call neko_log%message("Fields in checkpoint:", NEKO_LOG_DEBUG)
    do i = 1, this%state_list%size()
       si => this%state_list%get(i)
       call neko_log%message("  - " // trim(si%name), NEKO_LOG_DEBUG)
    end do

    call neko_log%end_section()

  end subroutine checkpoint_init_from_components

  !> Free
  !> Free checkpointing resources.
  !! @param[inout] this Checkpointing implementation.
  subroutine checkpoint_free(this)
    class(simulation_checkpoint_t), intent(inout) :: this
    integer :: i, j
    character(len=1024) :: file_name
    logical :: exists
    integer :: stat, unit

    ! Free the RAM Checkpoints
    if (allocated(this%state_storage)) then
       do i = 1, this%n_saves_memory
          do j = 1, this%state_list%size()
             call this%state_storage(i, j)%free()
          end do
       end do
    end if

    call this%state_list%free()
    if (allocated(this%state_storage)) deallocate(this%state_storage)

    ! Delete the checkpoint file list
    if (.not. this%keep_checkpoints .and. pe_rank .eq. 0) then
       do i = this%get_n_timesteps(), 1, -1
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
    call MPI_Barrier(NEKO_COMM)

    ! Reset to default values
    this%enabled = .false.
    this%filename = "checkpoint"
    this%fmt = "chkp"
    this%algorithm = "linear"
    this%n_saves_memory = 10
    this%keep_checkpoints = .false.

    this%n_saves_disc = 0
    call this%set_n_timesteps(0)
    this%first_valid_timestep = 2
    this%loaded_checkpoint = -1

  end subroutine checkpoint_free

  ! ========================================================================== !
  ! Saving and Restoring

  !> Save the current state of the simulation to disk
  !> Save forward state.
  !! @param[inout] this Checkpointing implementation.
  !! @param[inout] neko_case Case data structure.
  !! @param[in] time Current time state.
  subroutine checkpoint_save(this, neko_case, time)
    class(simulation_checkpoint_t), intent(inout) :: this
    class(case_t), intent(inout) :: neko_case
    type(time_state_t), intent(in) :: time

    if (.not. this%enabled) return

    call profiler_start_region("Checkpoint save")

    ! Update the number of recorded timesteps
    call this%set_n_timesteps(max(this%get_n_timesteps(), time%tstep))

    select case (this%algorithm)
    case ("linear")
       call checkpoint_save_linear(this, neko_case)
    case default
       call neko_error("Unknown checkpoint algorithm: " // this%algorithm)
    end select

    call profiler_end_region("Checkpoint save")
  end subroutine checkpoint_save

  !> Restore the forward simulation state
  !> Restore forward state for adjoint.
  !! @param[inout] this Checkpointing implementation.
  !! @param[inout] neko_case Case data structure.
  !! @param[in] time Target time state.
  subroutine checkpoint_restore(this, neko_case, time)
    class(simulation_checkpoint_t), intent(inout) :: this
    class(case_t), target, intent(inout) :: neko_case
    type(time_state_t), intent(in) :: time
    character(len=256) :: msg
    integer :: tstep

    if (.not. this%enabled) return

    call profiler_start_region("Checkpoint restore")

    tstep = time%tstep
    if (tstep .lt. 1 .or. tstep .gt. this%get_n_timesteps()) then
       write(msg, '(A,I0,A,I0,A)') "Requested timestep ", tstep, &
            " is out of range [1, ", this%get_n_timesteps(), "]"
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

  !> Save current data to the RAM checkpoint at the specified index.
  !! @param this The checkpoint object.
  !! @param index The index in the RAM checkpoint to save to.
  subroutine checkpoint_save_data(this, index)
    class(simulation_checkpoint_t), intent(inout) :: this
    integer, intent(in) :: index
    type(field_t), pointer :: si !< Pointer to the i'th state field
    integer :: i
    character(len=1024) :: msg

    if (index .lt. 1 .or. index .gt. this%n_saves_memory) then
       write(msg, '(A,I0,A,I0,A)') "Checkpoint save index ", index, &
            " is out of range [1, ", this%n_saves_memory, "]"
       call neko_error(trim(msg))
    end if

    ! Allocate the RAM checkpoint if not already allocated
    do i = 1, this%state_list%size()
       if (.not. this%state_storage(index, i)%is_allocated()) then
          si => this%state_list%get(i)
          call this%state_storage(index, i)%init(si%size())
       end if
    end do

    ! Save the current iterates to memory
    if (NEKO_BCKND_DEVICE .eq. 0) then
       do i = 1, this%state_list%size()
          si => this%state_list%get(i)
          call copy(this%state_storage(index, i)%data, si%x, si%size())
       end do
    else
       do i = 1, this%state_list%size()
          si => this%state_list%get(i)
          call device_memcpy(this%state_storage(index, i)%data, si%x_d, &
               si%size(), DEVICE_TO_HOST, this%state_list%size() .eq. i)
       end do
    end if

    nullify(si)
  end subroutine checkpoint_save_data

  !> Restore data from the RAM checkpoint at the specified index to the current
  !! state.
  !! @param this The checkpoint object.
  !! @param index The index in the RAM checkpoint to restore from.
  subroutine checkpoint_load_data(this, index)
    class(simulation_checkpoint_t), intent(inout) :: this
    integer, intent(in) :: index
    type(field_t), pointer :: si
    character(len=1024) :: msg
    integer :: i

    if (index .lt. 1 .or. index .gt. this%n_saves_memory) then
       write(msg, '(A,I0,A,I0,A)') "Checkpoint save index ", index, &
            " is out of range [1, ", this%n_saves_memory, "]"
       call neko_error(trim(msg))
    end if

    ! Save the current iterates to memory
    if (NEKO_BCKND_DEVICE .eq. 0) then
       do i = 1, this%state_list%size()
          si => this%state_list%get(i)
          call copy(si%x, this%state_storage(index, i)%data, si%size())
       end do
    else
       do i = 1, this%state_list%size()
          si => this%state_list%get(i)
          call device_memcpy(this%state_storage(index, i)%data, si%x_d, &
               si%size(), HOST_TO_DEVICE, this%state_list%size() .eq. i)
       end do
    end if

    nullify(si)
  end subroutine checkpoint_load_data

  ! ========================================================================== !
  ! Meta handling

  !> Reset the checkpoint data
  !> Reset checkpointing state.
  !! @param[inout] this Checkpointing implementation.
  subroutine checkpoint_reset(this)
    class(simulation_checkpoint_t), intent(inout) :: this
    integer :: i, j

    if (.not. this%enabled) return

    ! Reset our checkpoints
    this%loaded_checkpoint = -1
    this%n_saves_disc = 0
    call this%set_n_timesteps(0)

    do i = 1, size(this%state_storage, 1)
       do j = 1, size(this%state_storage, 2)
          if (this%state_storage(i, j)%is_allocated()) then
             call rzero(this%state_storage(i, j)%data, &
                  this%state_storage(i, j)%size)
          end if
       end do
    end do

  end subroutine checkpoint_reset

  ! -------------------------------------------------------------------------- !
  ! Host array routines

  subroutine host_array_init(this, size)
    class(host_array), intent(inout) :: this
    integer, intent(in) :: size

    call this%free()
    this%size = size
    allocate(this%data(size))
    call rzero(this%data, this%size)

  end subroutine host_array_init

  subroutine host_array_free(this)
    class(host_array), intent(inout) :: this

    this%size = 0
    if (allocated(this%data)) deallocate(this%data)

  end subroutine host_array_free

  pure function host_array_is_allocated(this) result(is_alloc)
    class(host_array), intent(in) :: this
    logical :: is_alloc

    is_alloc = allocated(this%data)

  end function host_array_is_allocated

end module simulation_checkpoint
