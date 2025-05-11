program usrneko

!  use stdlib_io_npy, only : save_npy
  use LightKrylov
  use LightKrylov, only: wp => dp, rtol => rtol_dp
  use LightKrylov_Logger
  use LightKrylov_Constants
  use cylinder, only: neko_propagator, state_vector, my_eigs
  use global_coef, only: global_coef_t, global_coef_getter
  use LightKrylov_IterativeSolvers, only: write_results_cdp
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
  integer, parameter :: nev = 15
  !> Size of Krylov subspace.
  integer, parameter :: kdim = 126
  !> Krylov subspace.
  type(state_vector), allocatable :: X(:)
  !> Eigenvalues.
  complex(kind=wp), allocatable :: lambda(:)
  !> Residual.
  real(kind=wp), allocatable    :: residuals(:)
  !> Information flag.
  integer          :: info

  !> writer
  type(state_vector), allocatable :: X_writer

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
  allocate(X(nev)); call initialize_krylov_subspace(X)

  !> initialize writer
  allocate(X_writer); call X_writer%init()


  !---------------------------------------
  !-----     CHECK LINEAR SOLVER     -----
  !---------------------------------------
  ! check just a massive forward run:
  ! call X(1)%rand()
  ! do i = 1, 100
  !    call A%matvec(X(1), X(2))
  !    call X(1)%axpby(1.0_wp, X(2), 0.0_wp)
  !    call X(1)%scal(1.0_wp/x(1)%norm())
  ! end do

  !------------------------------------------
  !-----     EIGENVALUE COMPUTATION     -----
  !------------------------------------------

  !> Call to LightKrylov.
  ! call eigs(A, X, lambda, residuals, info)
  call eigs(A, X, lambda, residuals, info, kdim = kdim, tolerance = rtol)
  
  ! call my_eigs(A, X, lambda, residuals, info, kdim = kdim, tolerance = rtol)

  call check_info(info, 'eigs', module=this_module, procedure='main')
  !> Transform eigenspectrum from unit-disk to standard complex plane.
  lambda = log(lambda) / tau

  !--------------------------------
  !-----     SAVE TO DISK     -----
  !--------------------------------

  !> Save eigenvectors to disk.
  do i = 1, nev
      call X_writer%axpby(1.0_wp, X(i), 0.0_wp)
      call X_writer%write(i)
  enddo

  !> Clean up
  call A%free()
  do i = 1, nev
    call X(i)%free()
  enddo

  !> Save the eigenspectrum.
  ! call save_eigenspectrum(lambda, residuals, "example/ginzburg_landau/eigenspectrum.npy")


end program usrneko
