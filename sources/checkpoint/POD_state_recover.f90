!> @file POD_state_recover.f90
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
module simulation_POD_state_recover
  use num_types, only: rp, sp
  use case, only: case_t
  use json_file_module, only: json_file
  use json_utils, only: json_get
  use field, only: field_t
  use file, only: file_t
  use matrix, only: matrix_t
  use fld_file_output, only: fld_file_output_t
  use coefs, only: coef_t
  use data_streamer, only: data_streamer_t
  use profiler, only: profiler_start_region, profiler_end_region
  use state_recover, only: state_recover_t
  use time_state, only: time_state_t
  use comm, only: neko_comm, mpi_real_precision
  use mpi_f08, only: MPI_Allreduce, MPI_IN_PLACE, MPI_SUM
  use neko_config, only: NEKO_BCKND_DEVICE
  use device, only: DEVICE_TO_HOST, device_memcpy
  use utils, only: neko_error
  use vector, only: vector_t
  use field_math, only: field_add2s2, field_rzero
  implicit none
  private

  type, public, extends(state_recover_t) :: POD_state_recover_t
     private

     ! ----------------------------------------------------------------------- !
     ! User parameters

     !> enabled
     logical :: enabled = .true.
     !> Has the csv been read and modes sent back?
     logical :: data_sent_back = .false.
     !> Data streamer Neko -> pySEMTools
     type(data_streamer_t) :: dstream_forward
     ! Data streamer pySEMTools -> Neko
     type(data_streamer_t) :: dstream_back
     !> Stream every ith timestep
     integer :: i_stream
     !> Number of modes to retain
     integer :: n_modes
     !> Reader for time coefs
     type(file_t) :: csv_reader
     !> Matrix for time coefs, columns: [t, a1, a2, ...]
     type(matrix_t) :: time_coefs
     !> Interpolated time coefficients at a single time instance
     type(vector_t) :: a_interp

     !> Structures to hold the POD mode
     type(field_t), dimension(:), allocatable :: u_modes
     type(field_t), dimension(:), allocatable :: v_modes
     type(field_t), dimension(:), allocatable :: w_modes
     !> in case we have scalars etc
     integer :: n_flds = 3

     !> writer to write modes
     type(fld_file_output_t) :: output

     !> Coefficients
     type(coef_t), pointer :: coef => null()

   contains
     !> Initialization from a JSON file
     procedure, public, pass(this) :: init => POD_state_recover_init_from_json
     procedure, public, pass(this) :: init_from_json => &
          POD_state_recover_init_from_json
     !> Initialization from components
     procedure, public, pass(this) :: init_from_components => &
          POD_state_recover_init_from_components
     !> Free
     procedure, public, pass(this) :: free => POD_state_recover_free
     !> Reset the POD_state_recover data
     procedure, public, pass(this) :: reset => POD_state_recover_reset
     !> Save the current state of the simulation to disk
     procedure, public, pass(this) :: save => POD_state_recover_save
     !> Restore the forward simulation state
     procedure, public, pass(this) :: restore => POD_state_recover_restore

  end type POD_state_recover_t

contains

  ! ========================================================================== !
  ! Initialization and deallocation

  !> Initialization
  subroutine POD_state_recover_init_from_json(this, neko_case, params)
    class(POD_state_recover_t), intent(inout) :: this
    class(case_t), target, intent(inout) :: neko_case
    type(json_file), target, intent(inout) :: params
    integer :: i_stream, n_modes

    call json_get(params, "i_stream", i_stream)
    call json_get(params, "n_modes", n_modes)

    call this%init_from_components(neko_case, i_stream, n_modes)

  end subroutine POD_state_recover_init_from_json

  !> Initialization from components
  subroutine POD_state_recover_init_from_components(this, neko_case, i_stream, &
       n_modes)
    class(POD_state_recover_t), intent(inout), target :: this
    class(case_t), target, intent(inout) :: neko_case
    integer, intent(in) :: i_stream
    integer, intent(in) :: n_modes
    integer :: i
    character(len=80) :: str
    

    this%i_stream = i_stream
    this%n_modes = n_modes
    this%coef => neko_case%fluid%c_Xh

    ! Allocate the modes
    allocate(this%u_modes(this%n_modes))
    allocate(this%v_modes(this%n_modes))
    allocate(this%w_modes(this%n_modes))

    ! Initialize writer
    call this%output%init(sp, 'POD_modes', this%n_flds * this%n_modes)
    do i = 1, this%n_modes
       write(str, '(A,I0)') "u_mode_", i
       call this%u_modes(i)%init(this%coef%dof, trim(str))
       call this%output%fields%assign(this%n_flds*(i-1) + 1, this%u_modes(i))

       write(str, '(A,I0)') "v_mode_", i
       call this%v_modes(i)%init(this%coef%dof, trim(str))
       call this%output%fields%assign(this%n_flds*(i-1) + 2, this%v_modes(i))

       write(str, '(A,I0)') "w_mode_", i
       call this%w_modes(i)%init(this%coef%dof, trim(str))
       call this%output%fields%assign(this%n_flds*(i-1) + 3, this%w_modes(i))
    end do

    ! Initialize the data streamer
    call this%dstream_forward%init(this%coef)

    ! Stream the mesh
    call this%dstream_forward%stream(this%coef%dof%x)
    call this%dstream_forward%stream(this%coef%dof%y)
    call this%dstream_forward%stream(this%coef%dof%z)

    ! set up the csv reader
    call this%csv_reader%init('pod_time_coeffs.csv')

    ! init vector
    call this%a_interp%init(this%n_modes)

  end subroutine POD_state_recover_init_from_components

  !> Free
  subroutine POD_state_recover_free(this)
    class(POD_state_recover_t), intent(inout) :: this
    integer :: i

    if (allocated(this%u_modes)) then
       do i = 1, size(this%u_modes)
          call this%u_modes(i)%free()
       end do
       deallocate(this%u_modes)
    end if
    if (allocated(this%v_modes)) then
       do i = 1, size(this%v_modes)
          call this%v_modes(i)%free()
       end do
       deallocate(this%v_modes)
    end if
    if (allocated(this%w_modes)) then
       do i = 1, size(this%w_modes)
          call this%w_modes(i)%free()
       end do
       deallocate(this%w_modes)
    end if

    call this%dstream_forward%free()
    call this%dstream_back%free()
    call this%csv_reader%free()
    nullify(this%coef)
    call this%a_interp%free()

    this%enabled = .false.
    this%i_stream = 0
    this%n_modes = 0

  end subroutine POD_state_recover_free

  ! ========================================================================== !
  ! Saving and Restoring

  !> Save the current state of the simulation to disk
  subroutine POD_state_recover_save(this, neko_case, time)
    class(POD_state_recover_t), intent(inout) :: this
    class(case_t), intent(inout) :: neko_case
    type(time_state_t), intent(in) :: time
    integer :: n
    type (field_t), pointer :: u, v, w

    if (.not. this%enabled) return
    if (mod(time%tstep, this%i_stream) .ne. 0) return

    u => neko_case%fluid%u
    v => neko_case%fluid%v
    w => neko_case%fluid%w
    n = u%dof%size()

    call profiler_start_region("POD save")

    if (NEKO_BCKND_DEVICE .eq. 1) then
       ! Ensure host buffers are up to date before streaming from a GPU run.
       call device_memcpy(u%x, u%x_d, n, DEVICE_TO_HOST, sync=.true.)
       call device_memcpy(v%x, v%x_d, n, DEVICE_TO_HOST, sync=.true.)
       call device_memcpy(w%x, w%x_d, n, DEVICE_TO_HOST, sync=.true.)
    end if

    ! Stream the data
    call this%dstream_forward%stream(u%x)
    call this%dstream_forward%stream(v%x)
    call this%dstream_forward%stream(w%x)


    call profiler_end_region("POD save")
  end subroutine POD_state_recover_save

  !> Restore the forward simulation state
  subroutine POD_state_recover_restore(this, neko_case, time)
    class(POD_state_recover_t), intent(inout) :: this
    class(case_t), target, intent(inout) :: neko_case
    type(time_state_t), intent(in) :: time
    character(len=256) :: msg
    integer :: i, ierr, n_lines, nrows, ncols

    if (.not. this%enabled) return

    if (.not. this%data_sent_back) then
        ! this is the first time a field is asked for, so we must request the
        ! modes and read the csv file

        ! Finalize the forward stream
        call this%dstream_forward%free()

        ! Give Python time to tear down reader and open the back-writer
        call execute_command_line("sleep 3")

        ! Now start the stream backward
        call this%dstream_back%init(this%coef)

        do i = 1, this%n_modes
            call this%dstream_back%recieve(this%u_modes(i)%x)
            call this%dstream_back%recieve(this%v_modes(i)%x)
            call this%dstream_back%recieve(this%w_modes(i)%x)
        end do

        ! Finalize the backward stream
        call this%dstream_back%free()

        ! Count lines in the CSV
        n_lines = csv_file_count_lines(this%csv_reader)

        nrows = n_lines - 1

        ! We expect columns = (time + n_modes)
        ncols = 1 + this%n_modes

        ! init the matrix
        call this%time_coefs%init(nrows, ncols)

        call this%csv_reader%read(this%time_coefs)

        ! share with everyone
        call MPI_Allreduce(MPI_IN_PLACE, this%time_coefs%x, &
             this%time_coefs%size(), mpi_real_precision, MPI_SUM, &
             neko_comm, ierr)

        this%data_sent_back = .true.
    end if

    call profiler_start_region("POD restore")

    ! interpolate the current time coeficients
    call interpolate_time_coeffs_vec(this%a_interp, this%time_coefs, time%t)
    ! reconstruct spatial fields
    call reconstruct_from_coeffs(this, neko_case, this%a_interp)

    call profiler_end_region("POD restore")
  end subroutine POD_state_recover_restore

  ! ========================================================================== !
  ! Meta handling

  !> Reset the POD_state_recover data
  subroutine POD_state_recover_reset(this)
    class(POD_state_recover_t), intent(inout) :: this
    integer :: i

    if (.not. this%enabled) return

    call this%set_n_timesteps(0)
    do i = 1, this%n_modes
       call field_rzero(this%u_modes(i))
       call field_rzero(this%v_modes(i))
       call field_rzero(this%w_modes(i))
    end do
    this%data_sent_back = .false.


  end subroutine POD_state_recover_reset

    function csv_file_count_lines(file_in) result(n)
    class(file_t), intent(in) :: file_in

    integer :: n
    integer :: ierr, file_unit

    open(file = trim(file_in%get_fname()), status = 'old', newunit = file_unit, &
         iostat = ierr)
    if (ierr .ne. 0) then
       call neko_error("Error while opening " // trim(file_in%get_fname()))
    end if
    rewind(file_unit)

    n = 0

    ! Keep reading (ierr = 0) until we reach the end (ierr != 0)
    do
       read (file_unit, *, iostat = ierr)
       if (ierr .ne. 0) exit
       n = n + 1
    end do
    rewind(file_unit)
    close(unit = file_unit)

  end function csv_file_count_lines

   !> Find indices i0,i1 such that t(i0) <= tq <= t(i1), with i1=i0+1.
  !! Assumes time column is sorted ascending.
  subroutine find_bracket_time(i0, i1, time_coefs, tq)
    integer, intent(out)       :: i0, i1
    type(matrix_t), intent(in) :: time_coefs
    real(kind=rp), intent(in)  :: tq

    integer :: lo, hi, mid
    real(kind=rp) :: tm

    lo = 1
    hi = time_coefs%get_nrows()

    do while (hi - lo > 1)
       mid = (lo + hi) / 2
       tm = time_coefs%x(mid, 1)
       if (tm <= tq) then
          lo = mid
       else
          hi = mid
       end if
    end do

    i0 = lo
    i1 = lo + 1
  end subroutine find_bracket_time


  !> Interpolate POD coefficients at time tq into a_out.
  !! time_coefs columns: [t, a1, a2, ...]
  !! Clamps to endpoints if tq outside range.
  subroutine interpolate_time_coeffs_vec(a_out, time_coefs, tq)
    type(vector_t), intent(inout) :: a_out
    type(matrix_t), intent(in)    :: time_coefs
    real(kind=rp), intent(in)     :: tq

    integer :: nrows, ncols, j, i0, i1
    real(kind=rp) :: t0, t1, w

    nrows = time_coefs%get_nrows()
    ncols = time_coefs%get_ncols()

    if (ncols < 2) call neko_error("time_coefs must have (t + coeffs) columns")
    if (size(a_out%x) /= ncols-1) call neko_error("a_out has wrong size")

    ! Clamp outside range
    if (tq <= time_coefs%x(1,1)) then
       do j = 1, ncols-1
          a_out%x(j) = time_coefs%x(1, j+1)
       end do
       return
    end if

    if (tq >= time_coefs%x(nrows,1)) then
       do j = 1, ncols-1
          a_out%x(j) = time_coefs%x(nrows, j+1)
       end do
       return
    end if

    call find_bracket_time(i0, i1, time_coefs, tq)

    t0 = time_coefs%x(i0,1)
    t1 = time_coefs%x(i1,1)

    if (abs(t1 - t0) < tiny(1.0_rp)) then
       w = 0.0_rp
    else
       w = (tq - t0) / (t1 - t0)
    end if

    do j = 1, ncols-1
       a_out%x(j) = (1.0_rp - w) * time_coefs%x(i0, j+1) + &
                     w          * time_coefs%x(i1, j+1)
    end do
  end subroutine interpolate_time_coeffs_vec


  !> Reconstruct u,v,w from POD modes using coefficients a.
  subroutine reconstruct_from_coeffs(this, neko_case, a)
    class(POD_state_recover_t), intent(inout) :: this
    class(case_t), intent(inout)              :: neko_case
    type(vector_t), intent(in)                :: a

    integer :: j
    type(field_t), pointer :: u, v, w

    u => neko_case%fluid%u
    v => neko_case%fluid%v
    w => neko_case%fluid%w

    call field_rzero(u)
    call field_rzero(v)
    call field_rzero(w)

    do j = 1, this%n_modes
       call field_add2s2(u, this%u_modes(j), a%x(j))
       call field_add2s2(v, this%v_modes(j), a%x(j))
       call field_add2s2(w, this%w_modes(j), a%x(j))
    end do
  end subroutine reconstruct_from_coeffs


end module simulation_POD_state_recover
