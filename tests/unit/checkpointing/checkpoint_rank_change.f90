!> Driver for check_checkpoint_rank_change.cmake: write or load the MMA and
!! the design HDF5 checkpoints on the ranks of the current run.
!! @details
!! Usage: `checkpoint_rank_change <mode> <file> [uneven]`, where `mode` is
!! `save` (both checkpoints into one file), `load_mma` or `load_design`.
!! The 6 design variables are spread as evenly as possible over the ranks,
!! or, with `uneven`, with one moved from the second rank to the first.
!! A load that does not match the save must stop with an error on every
!! rank; see stop_on_every_rank. A completed mode prints one line per rank.
program checkpoint_rank_change
  use, intrinsic :: iso_fortran_env, only: output_unit
  use num_types, only: rp
  use neko, only: neko_init, neko_finalize
  use comm, only: pe_rank, pe_size
  use utils, only: neko_error, throw_error
  use device, only: HOST_TO_DEVICE
  use vector, only: vector_t
  use mma, only: mma_t
  use checkpointing_test_utils, only: stored_design_t, stop_on_every_rank
  implicit none

  !> Global number of design variables.
  integer, parameter :: n_global = 6
  !> Number of constraints.
  integer, parameter :: m = 1
  character(len=16) :: mode, layout
  character(len=256) :: filename
  type(mma_t) :: opt
  type(stored_design_t) :: des
  type(vector_t) :: x
  real(kind=rp), allocatable :: xmin(:), xmax(:)
  real(kind=rp) :: a(m), c(m), d(m)
  integer :: n

  call get_command_argument(1, mode)
  call get_command_argument(2, filename)
  call get_command_argument(3, layout)

  call neko_init()
  throw_error => stop_on_every_rank

  n = n_global / pe_size
  if (pe_rank .lt. mod(n_global, pe_size)) n = n + 1
  if (trim(layout) .eq. 'uneven' .and. pe_size .gt. 1) then
     if (pe_rank .eq. 0) n = n + 1
     if (pe_rank .eq. 1) n = n - 1
  end if

  call x%init(n)
  if (n .gt. 0) x%x = 0.5_rp
  call x%copy_from(HOST_TO_DEVICE, sync = .true.)

  allocate(xmin(n), xmax(n))
  xmin = 0.0_rp
  xmax = 1.0_rp
  a = 0.0_rp
  c = 1000.0_rp
  d = 1.0_rp
  call opt%init_from_components(x, n, m, 1.0_rp, a, c, d, xmin, xmax, &
       bcknd = 'cpu', subsolver = 'dip')
  ! Only the layout of the design values matters here
  call des%init(xmax)

  select case (trim(mode))
  case ('save')
     call opt%save_checkpoint(trim(filename))
     call des%save_checkpoint(trim(filename))
  case ('load_mma')
     call opt%load_checkpoint(trim(filename))
  case ('load_design')
     call des%load_checkpoint(trim(filename))
  case default
     call neko_error('checkpoint_rank_change: unknown mode "' // &
          trim(mode) // '"')
  end select

  write(output_unit, '(A,I0)') 'checkpoint_rank_change: ' // trim(mode) // &
       ' completed on rank ', pe_rank
  flush(output_unit)

  call des%free()
  call opt%free()
  call x%free()
  deallocate(xmin, xmax)
  call neko_finalize()

end program checkpoint_rank_change
