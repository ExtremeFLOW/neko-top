program usrneko

!  use stdlib_io_npy, only : save_npy
  use LightKrylov
  use LightKrylov, only: wp => dp
  use LightKrylov_Logger
  use cylinder
  implicit none

  character(len=128), parameter :: this_module = 'Example cylinder'

  !------------------------------------------------
  !-----     LINEAR OPERATOR INVESTIGATED     -----
  !------------------------------------------------

  !> Exponential propagator.
  type(neko_propagator), allocatable :: A
  !> Sampling time.
  real(kind=wp) :: tau

  !---------------------------------------------------
  !-----     KRYLOV-BASED EIGENDECOMPOSITION     -----
  !---------------------------------------------------

  !> Number of eigenvalues we wish to converge.
  integer, parameter :: nev = 32
  !> Krylov subspace.
  type(state_vector), allocatable :: X(:)
  !> Eigenvalues.
  complex(kind=wp), allocatable :: lambda(:)
  !> Residual.
  real(kind=wp), allocatable    :: residuals(:)
  !> Information flag.
  integer          :: info

  !> Miscellaneous.
  integer       :: i
  !complex(wp) :: eigenvectors(nx, nev)

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

  !> Initialize Krylov subspace.
  allocate(X(nev))
  !> Use coef to initialize all the fields
  do i = 1, nev
     call X(i)%init(A%neko_case%fluid%C_Xh)
  end do
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

  !> Save the eigenspectrum.
  call save_eigenspectrum(lambda, residuals, "example/ginzburg_landau/eigenspectrum.npy")

  !> Reconstruct the leading eigenvectors from the Krylov basis.
  ! do i = 1, nev
  !   eigenvectors(:, i) = X(i)%state
  ! enddo

  !> Save eigenvectors to disk.
!  call save_npy("example/ginzburg_landau/eigenvectors.npy", eigenvectors)


end program usrneko
