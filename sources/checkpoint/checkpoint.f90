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
!
module simulation_checkpoint
  use num_types, only: rp, sp, dp
  use case, only: case_t
  use json_file_module, only: json_file
  use json_utils, only: json_get_or_default
  use scalar_scheme, only: scalar_scheme_t
  use time_step_controller, only: time_step_controller_t
  use time_state, only: time_state_t
  use chkp_output, only: chkp_output_t
  use field, only: field_t
  use file, only: file_t
  use mpi_f08, only: MPI_WTIME
  use simulation, only: simulation_step, simulation_restart
  use field_math, only: field_copy, field_rzero
  use utils, only: neko_error
  implicit none
  private

  type, public :: simulation_checkpoint_t

     character(len=256) :: algorithm = "linear"
     character(len=256) :: filename = "checkpoint_"
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

     logical :: have_scalar = .false.
     type(field_t), dimension(:), allocatable :: s_list

   contains
     !> Initialization
     generic :: init => init_from_json, init_from_components
     procedure, pass(this) :: init_from_json => checkpoint_init_from_json
     procedure, pass(this) :: init_from_components => &
          checkpoint_init_from_components
     !> Free
     procedure, pass(this) :: free => checkpoint_free
     !> Reset the checkpoint data
     procedure, pass(this) :: reset => checkpoint_reset
     !> Save the current state of the simulation to disk
     procedure, pass(this) :: save => checkpoint_save
     !> Restore the forward simulation state
     procedure, pass(this) :: restore => checkpoint_restore

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
  subroutine checkpoint_init_from_json(this, neko_case, parameters)
    class(simulation_checkpoint_t), intent(inout) :: this
    class(case_t), target, intent(inout) :: neko_case
    type(json_file), target, intent(inout) :: parameters
    integer :: n_saves_memory
    character(len=:), allocatable :: filename, algorithm

    call json_get_or_default(parameters, "algorithm", algorithm, "linear")
    call json_get_or_default(parameters, "n_memory", n_saves_memory, 10)
    call json_get_or_default(parameters, "filename", filename, "checkpoint_")

    call this%init_from_components(neko_case, algorithm, n_saves_memory, &
         filename)
  end subroutine checkpoint_init_from_json

  !> Initialization from components
  subroutine checkpoint_init_from_components(this, neko_case, algorithm, &
       n_saves_memory, chkp_file_name)
    class(simulation_checkpoint_t), intent(inout), target :: this
    class(case_t), target, intent(inout) :: neko_case
    character(len=*), intent(in) :: algorithm
    integer, intent(in) :: n_saves_memory
    character(len=*), intent(in) :: chkp_file_name
    class(scalar_scheme_t), pointer :: scalar_i
    integer :: i, j, n_scalars
    character(len=10) :: str

    call this%free()

    ! Set internal parameters
    this%algorithm = algorithm
    this%filename = chkp_file_name
    this%n_saves_memory = n_saves_memory
    this%have_scalar = allocated(neko_case%scalars)

    ! Initialize the Neko checkpoint output
    call this%chkp_output%init(neko_case%chkp, this%filename)

    ! Allocate the RAM Checkpoints
    allocate(this%p_list(this%n_saves_memory))
    allocate(this%u_list(this%n_saves_memory))
    allocate(this%v_list(this%n_saves_memory))
    allocate(this%w_list(this%n_saves_memory))
    if (this%have_scalar) then
       n_scalars = size(neko_case%scalars%scalar_fields)
       allocate(this%s_list(this%n_saves_memory * n_scalars))
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
       if (this%have_scalar) then
          do j = 1, n_scalars
             write(str, '(A,I0,A,I0)') "s_chkp_", i, "_", j
             scalar_i => neko_case%scalars%scalar_fields(j)
             call this%s_list((i - 1) * n_scalars + j)%init(scalar_i%s%dof, str)
          end do
       end if
    end do

  end subroutine checkpoint_init_from_components

  !> Free
  subroutine checkpoint_free(this)
    class(simulation_checkpoint_t), intent(inout) :: this
    integer :: i, n_scalars

    ! Free the RAM Checkpoints
    do i = 1, this%n_saves_memory
       if (allocated(this%p_list)) call this%p_list(i)%free()
       if (allocated(this%u_list)) call this%u_list(i)%free()
       if (allocated(this%v_list)) call this%v_list(i)%free()
       if (allocated(this%w_list)) call this%w_list(i)%free()
    end do

    if (allocated(this%s_list)) then
       n_scalars = size(this%s_list)
       do i = 1, n_scalars
          call this%s_list(i)%free()
       end do
    end if

    if (allocated(this%p_list)) deallocate(this%p_list)
    if (allocated(this%u_list)) deallocate(this%u_list)
    if (allocated(this%v_list)) deallocate(this%v_list)
    if (allocated(this%w_list)) deallocate(this%w_list)
    if (allocated(this%s_list)) deallocate(this%s_list)

  end subroutine checkpoint_free

  ! ========================================================================== !
  ! Meta handling

  !> Reset the checkpoint data
  subroutine checkpoint_reset(this)
    class(simulation_checkpoint_t), intent(inout) :: this
    integer :: i

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

    do i = 1, size(this%s_list)
       call field_rzero(this%s_list(i))
    end do
  end subroutine checkpoint_reset

  !> Save the current state of the simulation to disk
  subroutine checkpoint_save(this, neko_case)
    class(simulation_checkpoint_t), intent(inout) :: this
    class(case_t), intent(inout) :: neko_case

    ! Update the number of recorded timesteps
    this%n_timesteps = this%n_timesteps + 1

    select case (this%algorithm)
    case ("linear")
       call checkpoint_save_linear(this, neko_case)
    case default
       call neko_error("Unknown checkpoint algorithm: " // this%algorithm)
    end select
  end subroutine checkpoint_save

  !> Save the current state of the simulation in a linear fashion
  module subroutine checkpoint_save_linear(this, neko_case)
    class(simulation_checkpoint_t), intent(inout) :: this
    class(case_t), intent(inout) :: neko_case

    ! Sample the checkpoint if needed
    if (modulo(neko_case%time%tstep, this%n_saves_memory) .eq. 0 .or. &
         neko_case%time%tstep .le. this%first_valid_timestep) then

       call this%chkp_output%set_counter(neko_case%time%tstep)
       call this%chkp_output%sample(neko_case%time%t)
       this%n_saves_disc = this%n_saves_disc + 1
    end if
  end subroutine checkpoint_save_linear

  !> Restore the forward simulation state
  subroutine checkpoint_restore(this, neko_case, tstep)
    class(simulation_checkpoint_t), intent(inout) :: this
    class(case_t), target, intent(inout) :: neko_case
    integer, intent(in) :: tstep

    select case (this%algorithm)
    case ("linear")
       call checkpoint_restore_linear(this, neko_case, tstep)
    case default
       call neko_error("Unknown checkpoint algorithm: " // this%algorithm)
    end select
  end subroutine checkpoint_restore

  !> Restore the forward simulation state in a linear fashion
  module subroutine checkpoint_restore_linear(this, neko_case, tstep)
    class(simulation_checkpoint_t), intent(inout) :: this
    class(case_t), target, intent(inout) :: neko_case
    integer, intent(in) :: tstep
    type(time_step_controller_t) :: dt_controller
    type(file_t) :: chkp_file
    character(len=256) :: chkp_file_name
    real(kind=dp) :: loop_start
    integer :: j, k, previous_save, next_save
    integer :: i_scalars, n_scalars
    type(field_t), pointer :: u, v, w, p, s

    loop_start = MPI_WTIME()

    u => neko_case%fluid%u
    v => neko_case%fluid%v
    w => neko_case%fluid%w
    p => neko_case%fluid%p
    s => null()

    ! Determine the nearest save states
    previous_save = tstep - modulo(tstep, this%n_saves_memory)
    next_save = previous_save + this%n_saves_memory

    if (previous_save .le. this%first_valid_timestep) then
       previous_save = min(tstep, this%first_valid_timestep)
    end if

    if (previous_save .lt. this%first_valid_timestep) then
       next_save = previous_save + 1
    end if

    if (this%loaded_checkpoint .ne. previous_save) then
       write(chkp_file_name, '(A,I5.5,A)') trim(this%filename), &
            previous_save, '.chkp'

       call chkp_file%init(chkp_file_name)
       call chkp_file%read(neko_case%chkp)

       call dt_controller%init(neko_case%params)
       call simulation_restart(neko_case, neko_case%chkp)
       neko_case%time%tstep = previous_save
       this%loaded_checkpoint = neko_case%time%tstep

       do k = previous_save, min(next_save - 1, this%n_timesteps)

          if (k .ne. previous_save) then
             if (neko_case%time%t .ge. neko_case%time%end_time) exit
             call simulation_step(neko_case, dt_controller, loop_start)
          end if

          j = modulo(k, this%n_saves_memory) + 1
          call field_copy(this%p_list(j), p)
          call field_copy(this%u_list(j), u)
          call field_copy(this%v_list(j), v)
          call field_copy(this%w_list(j), w)
          if (this%have_scalar) then
             n_scalars = size(neko_case%scalars%scalar_fields)
             do i_scalars = 1, n_scalars
                s => neko_case%scalars%scalar_fields(i_scalars)%s
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
       n_scalars = size(neko_case%scalars%scalar_fields)
       do i_scalars = 1, n_scalars
          s => neko_case%scalars%scalar_fields(i_scalars)%s
          call field_copy(s, this%s_list((j - 1) * n_scalars + i_scalars))
       end do
    end if

  end subroutine checkpoint_restore_linear

end module simulation_checkpoint
