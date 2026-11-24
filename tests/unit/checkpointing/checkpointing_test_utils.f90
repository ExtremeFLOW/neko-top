module checkpointing_test_utils
  use simulation_m, only: simulation_t
  use state_recover_checkpoint, only: state_recover_checkpoint_t
  use design, only: design_t, design_factory
  use problem, only: problem_t
  use neko, only: neko_init
  use num_types, only: rp
  use json_module, only: json_file
  use json_utils, only: json_get
  use json_utils_ext, only: json_read_file
  use neko_top, only: neko_top_register_types
  use field, only: field_t
  use field_array, only: field_array_t
  use field_math, only: field_glsubnorm
  use registry, only: neko_registry
  use scratch_registry, only: neko_scratch_registry
  use vector, only: vector_t
  use device, only: HOST_TO_DEVICE, DEVICE_TO_HOST
  use comm, only: NEKO_COMM
  use mpi_f08, only: MPI_Barrier
  implicit none
  private

  logical, save :: runtime_initialized = .false.

  !> Design that only stores its values, for the design checkpoint tests.
  type, extends(design_t), public :: stored_design_t
     !> The values of the design variables on this rank.
     type(vector_t) :: values
   contains
     procedure, pass(this) :: init => stored_design_init
     procedure, pass(this) :: free => stored_design_free
     procedure, pass(this) :: get_values => stored_design_get_values
     procedure, pass(this) :: update_design => stored_design_update_design
     procedure, pass(this) :: map_forward => stored_design_map_forward
     procedure, pass(this) :: map_backward => stored_design_map_backward
     procedure, pass(this) :: write => stored_design_write
  end type stored_design_t

  public :: stop_on_every_rank
  public :: initialize_runtime_once
  public :: initialize_case_from_file
  public :: finalize_case_objects
  public :: cleanup_runtime_state
  public :: allocate_scalar_snapshot_fields
  public :: free_scalar_snapshot_fields
  public :: allocate_snapshot_fields
  public :: free_snapshot_fields
  public :: rmse

contains

  subroutine initialize_runtime_once()
    if (runtime_initialized) return

    call neko_init()
    call neko_top_register_types()

    runtime_initialized = .true.
  end subroutine initialize_runtime_once

  ! Shared case bootstrapping used by multiple checkpoint tests.
  subroutine initialize_case_from_file(parameter_file, sim, des, prob)
    character(len=*), intent(in) :: parameter_file
    type(simulation_t), intent(inout), target :: sim
    class(design_t), allocatable, intent(inout) :: des
    type(problem_t), intent(inout) :: prob
    type(json_file) :: parameters, design_parameters

    parameters = json_read_file(trim(parameter_file))
    call json_get(parameters, 'optimization.design', design_parameters)

    call sim%init(parameters)
    call design_factory(des, design_parameters, sim)
    call prob%init(parameters, des, sim)
  end subroutine initialize_case_from_file

  subroutine cleanup_runtime_state()
    call neko_registry%free()
    call neko_scratch_registry%free()
    call neko_registry%init()
    call neko_scratch_registry%init()
  end subroutine cleanup_runtime_state

  function rmse(field1, field2) result(result)
    type(field_t), intent(in) :: field1
    type(field_t), intent(in) :: field2
    real(kind=rp) :: result

    result = sqrt(field_glsubnorm(field1, field2)**2 / real(field1%size(), rp))
  end function rmse

  subroutine allocate_snapshot_series(field, n_timesteps, snapshots)
    type(field_t), intent(in) :: field
    integer, intent(in) :: n_timesteps
    type(field_array_t), intent(inout) :: snapshots
    integer :: i

    call snapshots%init(n_timesteps)
    do i = 1, n_timesteps
       call snapshots%assign(i, field)
    end do
  end subroutine allocate_snapshot_series

  subroutine allocate_snapshot_fields(p, u, v, w, n_timesteps, &
       p_fields, u_fields, v_fields, w_fields)
    type(field_t), intent(in) :: p, u, v, w
    integer, intent(in) :: n_timesteps
    type(field_array_t), intent(inout) :: p_fields, u_fields, v_fields, &
         w_fields

    call allocate_snapshot_series(p, n_timesteps, p_fields)
    call allocate_snapshot_series(u, n_timesteps, u_fields)
    call allocate_snapshot_series(v, n_timesteps, v_fields)
    call allocate_snapshot_series(w, n_timesteps, w_fields)
  end subroutine allocate_snapshot_fields

  subroutine allocate_scalar_snapshot_fields(sim, n_timesteps, scalar_fields)
    type(simulation_t), intent(inout) :: sim
    integer, intent(in) :: n_timesteps
    type(field_array_t), allocatable, intent(inout) :: scalar_fields(:)

    type(field_t), pointer :: scalar
    integer :: i, n_scalars

    call free_scalar_snapshot_fields(scalar_fields)

    if (.not. allocated(sim%neko_case%scalars)) return

    n_scalars = size(sim%neko_case%scalars%scalar_fields)
    if (n_scalars .le. 0) return

    allocate(scalar_fields(n_scalars))

    do i = 1, n_scalars
       scalar => sim%neko_case%scalars%scalar_fields(i)%scalar%s
       call allocate_snapshot_series(scalar, n_timesteps, scalar_fields(i))
    end do
  end subroutine allocate_scalar_snapshot_fields

  subroutine free_scalar_snapshot_fields(scalar_fields)
    type(field_array_t), allocatable, intent(inout) :: scalar_fields(:)

    integer :: i

    if (.not. allocated(scalar_fields)) return

    do i = 1, size(scalar_fields)
       call scalar_fields(i)%free()
    end do

    deallocate(scalar_fields)
  end subroutine free_scalar_snapshot_fields

  subroutine free_snapshot_fields(p_fields, u_fields, v_fields, w_fields)
    type(field_array_t), intent(inout) :: p_fields, u_fields, v_fields, w_fields

    call p_fields%free()
    call u_fields%free()
    call v_fields%free()
    call w_fields%free()
  end subroutine free_snapshot_fields

  subroutine finalize_case_objects(sim, des, prob, chkp)
    type(simulation_t), intent(inout) :: sim
    class(design_t), allocatable, intent(inout) :: des
    type(problem_t), intent(inout) :: prob
    type(state_recover_checkpoint_t), optional, intent(inout) :: chkp

    if (present(chkp)) call chkp%free()
    call sim%free()
    if (allocated(des)) then
       call des%free()
       deallocate(des)
    end if
    call prob%free()
  end subroutine finalize_case_objects

  ! ========================================================================== !
  ! Error handling for tests that expect a run to stop

  !> Error handler that stops only once every rank has raised an error.
  !! @details neko_error prints its message before calling this, so every
  !! message reaches the output before any rank stops. An error raised on
  !! some ranks only never passes the barrier, which shows as a hang.
  !! @param filename The file that raised the error.
  !! @param line_number The line that raised the error.
  !! @param message The error message, if any.
  subroutine stop_on_every_rank(filename, line_number, message)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: line_number
    character(len=*), optional, intent(in) :: message

    call MPI_Barrier(NEKO_COMM)
    error stop 1

  end subroutine stop_on_every_rank

  ! ========================================================================== !
  ! stored_design_t

  !> Initialise the design with the values given on this rank.
  !! @param x0 The values.
  subroutine stored_design_init(this, x0)
    class(stored_design_t), intent(inout) :: this
    real(kind=rp), intent(in) :: x0(:)

    call this%free()
    call this%init_base('stored_design', size(x0))
    this%values = x0

  end subroutine stored_design_init

  !> Free the design.
  subroutine stored_design_free(this)
    class(stored_design_t), intent(inout) :: this

    call this%values%free()
    call this%free_base()

  end subroutine stored_design_free

  !> Copy the stored values into `values`, on the host and the device.
  !! @param values The vector to fill.
  subroutine stored_design_get_values(this, values)
    class(stored_design_t), intent(in) :: this
    type(vector_t), intent(inout) :: values

    if (this%size() .gt. 0) values%x(1:this%size()) = this%values%x(:)
    call values%copy_from(HOST_TO_DEVICE, .true.)

  end subroutine stored_design_get_values

  !> Store `values`, whose device side is current.
  !! @param values The new values.
  subroutine stored_design_update_design(this, values)
    class(stored_design_t), intent(inout) :: this
    type(vector_t), intent(inout) :: values

    if (this%size() .eq. 0) return
    call values%copy_from(DEVICE_TO_HOST, .true.)
    this%values%x(:) = values%x(1:this%size())
    call this%values%copy_from(HOST_TO_DEVICE, .true.)

  end subroutine stored_design_update_design

  !> No mapping: the design has no simulation to map to.
  subroutine stored_design_map_forward(this)
    class(stored_design_t), intent(inout) :: this
  end subroutine stored_design_map_forward

  !> No mapping: the design has no simulation to map from.
  !! @param sensitivity The sensitivity, unused.
  subroutine stored_design_map_backward(this, sensitivity)
    class(stored_design_t), intent(inout) :: this
    type(vector_t), intent(in) :: sensitivity
  end subroutine stored_design_map_backward

  !> No output.
  !! @param idx The output index, unused.
  subroutine stored_design_write(this, idx)
    class(stored_design_t), intent(inout) :: this
    integer, intent(in) :: idx
  end subroutine stored_design_write

end module checkpointing_test_utils
