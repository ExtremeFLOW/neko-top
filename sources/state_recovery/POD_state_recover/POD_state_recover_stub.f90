!> @file POD_state_recover_stub.f90
!! @brief Stub POD state recovery used when ADIOS2 is unavailable.
module simulation_POD_state_recover
  use case, only: case_t
  use json_file_module, only: json_file
  use state_recover, only: state_recover_t
  use time_state, only: time_state_t
  use utils, only: neko_error
  implicit none
  private

  type, public, extends(state_recover_t) :: POD_state_recover_t
   contains
     procedure, public, pass(this) :: init => POD_state_recover_init_from_json
     procedure, public, pass(this) :: init_from_json => &
          POD_state_recover_init_from_json
     procedure, public, pass(this) :: init_from_components => &
          POD_state_recover_init_from_components
     procedure, public, pass(this) :: free => POD_state_recover_free
     procedure, public, pass(this) :: reset => POD_state_recover_reset
     procedure, public, pass(this) :: save => POD_state_recover_save
     procedure, public, pass(this) :: restore => POD_state_recover_restore
  end type POD_state_recover_t

contains

  subroutine POD_state_recover_init_from_json(this, neko_case, params)
    class(POD_state_recover_t), intent(inout) :: this
    class(case_t), target, intent(inout) :: neko_case
    type(json_file), target, intent(inout) :: params

    call POD_state_recover_unavailable()
  end subroutine POD_state_recover_init_from_json

  subroutine POD_state_recover_init_from_components(this, neko_case, &
       i_stream, n_modes, dtype, write_modes, output_reconstruction, &
       output_control, output_value, debug, &
       mode_output_precision, mode_output_format, mode_file_name, &
       recon_output_precision, recon_output_format, recon_file_name)
    class(POD_state_recover_t), intent(inout) :: this
    class(case_t), target, intent(inout) :: neko_case
    integer, intent(in) :: i_stream, n_modes
    character(len=*), intent(in) :: dtype
    logical, intent(in) :: write_modes
    logical, intent(in) :: output_reconstruction
    character(len=*), intent(in) :: output_control
    real, intent(in) :: output_value
    logical, intent(in), optional :: debug
    integer, intent(in), optional :: mode_output_precision
    character(len=*), intent(in), optional :: mode_output_format
    character(len=*), intent(in), optional :: mode_file_name
    integer, intent(in), optional :: recon_output_precision
    character(len=*), intent(in), optional :: recon_output_format
    character(len=*), intent(in), optional :: recon_file_name

    call POD_state_recover_unavailable()
  end subroutine POD_state_recover_init_from_components

  subroutine POD_state_recover_free(this)
    class(POD_state_recover_t), intent(inout) :: this
  end subroutine POD_state_recover_free

  subroutine POD_state_recover_reset(this)
    class(POD_state_recover_t), intent(inout) :: this
  end subroutine POD_state_recover_reset

  subroutine POD_state_recover_save(this, neko_case)
    class(POD_state_recover_t), intent(inout) :: this
    class(case_t), intent(inout) :: neko_case

    call POD_state_recover_unavailable()
  end subroutine POD_state_recover_save

  subroutine POD_state_recover_restore(this, neko_case, time)
    class(POD_state_recover_t), intent(inout) :: this
    class(case_t), target, intent(inout) :: neko_case
    type(time_state_t), intent(in) :: time

    call POD_state_recover_unavailable()
  end subroutine POD_state_recover_restore

  subroutine POD_state_recover_unavailable()
    call neko_error('POD state recovery requires an ADIOS2-enabled build.')
  end subroutine POD_state_recover_unavailable

end module simulation_POD_state_recover
