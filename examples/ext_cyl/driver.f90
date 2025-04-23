program usrneko

!  use stdlib_io_npy, only : save_npy
  use LightKrylov
  use LightKrylov, only: wp => dp
  use LightKrylov_Logger
  use cylinder
  use global_coef, only: global_coef_t, global_coef_getter
  implicit none

  character(len=128), parameter :: this_module = 'Example cylinder'

  !------------------------------------------------
  !-----     LINEAR OPERATOR INVESTIGATED     -----
  !------------------------------------------------

  !> Exponential propagator.
  type(neko_propagator), allocatable :: A
  !> Sampling time.
  real(kind=wp) :: tau
  !> A way to access coef globally
  type(global_coef_t), target :: my_global_coef_getter

  !---------------------------------------------------
  !-----     KRYLOV-BASED EIGENDECOMPOSITION     -----
  !---------------------------------------------------

  !> Number of eigenvalues we wish to converge.
  integer, parameter :: nev = 8
  !> Krylov subspace.
  type(state_vector), allocatable :: X(:)
  !> Eigenvalues.
  complex(kind=wp), allocatable :: lambda(:)
  !> Residual.
  real(kind=wp), allocatable    :: residuals(:)
  !> Information flag.
  integer          :: info

  !> Miscellaneous.
  integer :: i, j

  !=============================================================================

  !----------------------------------
  !-----     INITIALIZATION     -----
  !----------------------------------

  !> Set up logging
  call logger_setup()

  !> Initialize exponential propagator.
  allocate(A)
  call A%init()

  !> Get the integration time
  tau = real(A%neko_case%time%end_time,kind=wp)

  !> Extract the global coef from neko
  my_global_coef_getter%global_coef = A%neko_case%fluid%c_Xh
  global_coef_getter => my_global_coef_getter

  !> Initialize Krylov subspace.
  allocate(X(nev))
  ! (this should also initialize them)
  call zero_basis(X)

  !------------------------------------------
  !-----     EIGENVALUE COMPUTATION     -----
  !------------------------------------------

  !> Call to LightKrylov.
  call eigs(A, X, lambda, residuals, info)
  call check_info(info, 'eigs', module=this_module, procedure='main')

  !> Transform eigenspectrum from unit-disk to standard complex plane.
  lambda = log(lambda) / tau

  !--------------------------------
  !-----     SAVE TO DISK     -----
  !--------------------------------

  !> just write them out...
  print *, "SPECTRA"
  do i = 1, nev
    print *, i, lambda(i)
  enddo

  !> Save eigenvectors to disk.
  do i = 1, nev
    ! this is also so stupid, we need a better way to write these out
    do j = 1, nev - i + 1
      call X( nev - i + 1)%write(j)
    end do
  enddo



  !> Clean up
  call A%free()
  do i = 1, nev
    call X(i)%free()
  enddo



  !> Save the eigenspectrum.
  ! call save_eigenspectrum(lambda, residuals, "example/ginzburg_landau/eigenspectrum.npy")


end program usrneko
