! Re-implementation of tests/regression/mma (1D cantilever beam: weight
! objective, tip-deflection + 10 stress constraints) driving the UNMODIFIED
! neko-top mma_t through the same calls mma_optimizer.f90 makes.
program beam_driver
  use mpi_f08
  use num_types, only: rp
  use comm, only: NEKO_COMM, MPI_REAL_PRECISION
  use vector, only: vector_t
  use matrix, only: matrix_t
  use mma, only: mma_t
  implicit none
  integer, parameter :: n = 221184, m = 11, ncon = 10
  real(rp), parameter :: L_total = 2.0_rp, b = 0.02_rp, rho = 7800.0_rp
  real(rp), parameter :: E = 210.0e9_rp, P = 1000.0_rp, h_min = 0.005_rp
  real(rp), parameter :: h_max = 0.05_rp, u_tip_max = 0.25_rp
  real(rp), parameter :: sigma_max = 250e6_rp
  type(mma_t) :: opt
  type(vector_t) :: x, xold, df0dx, fval
  type(matrix_t) :: dfdx
  real(rp) :: f0, a(m), c(m), d(m), xmin(n), xmax(n), Delta(n), Le, maxdx
  real(rp) :: asyinit, asyincr, asydecr, cval, move_limit
  integer :: idx(ncon), iter, ierr, k, nit, u
  character(len=64) :: subsolver, outname, arg

  call MPI_Init(ierr)
  NEKO_COMM = MPI_COMM_WORLD
  MPI_REAL_PRECISION = MPI_DOUBLE_PRECISION

  ! args: subsolver outname nit asyinit asyincr asydecr c move_limit
  call get_command_argument(1, subsolver)
  call get_command_argument(2, outname)
  call get_command_argument(3, arg); read(arg, *) nit
  call get_command_argument(4, arg); read(arg, *) asyinit
  call get_command_argument(5, arg); read(arg, *) asyincr
  call get_command_argument(6, arg); read(arg, *) asydecr
  call get_command_argument(7, arg); read(arg, *) cval
  call get_command_argument(8, arg); read(arg, *) move_limit

  call x%init(n); call xold%init(n); call df0dx%init(n); call fval%init(m)
  call dfdx%init(m, n)
  x = 0.5_rp

  call fill_constraint_indices(idx, ncon, ncon, n)
  Le = L_total / real(n, rp)
  do k = 1, n
     Delta(k) = ((L_total - Le*real(k-1, rp))**3 - &
          (L_total - Le*real(k, rp))**3) / 3.0_rp
  end do

  a = 0.0_rp; c = cval; d = 1.0_rp; xmin = 0.0_rp; xmax = 1.0_rp
  call opt%init(x, n, m, 1.0_rp, a, c, d, xmin, xmax, max_iter = 100, &
       asyinit = asyinit, asyincr = asyincr, asydecr = asydecr, &
       bcknd = 'cpu', subsolver = trim(subsolver), move_limit = move_limit)

  call evaluate(x%x, f0, df0dx%x, fval%x, dfdx%x)
  call opt%kkt(x, df0dx, fval, dfdx)   ! as mma_optimizer_initialize does

  open(newunit=u, file=trim(outname)//'.csv', status='replace')
  write(u, '(A)') 'iter,f0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,kktmax,kktnorm,maxdx'
  write(u, '(I0,15(",",ES24.16))') 0, f0, fval%x, 0.0_rp, 0.0_rp, 0.0_rp
  open(newunit=k, file=trim(outname)//'.x', access='stream', &
       form='unformatted', status='replace')
  do iter = 1, nit
     xold = x
     call opt%update(iter, x, df0dx, fval, dfdx)
     call evaluate(x%x, f0, df0dx%x, fval%x, dfdx%x)
     call opt%kkt(x, df0dx, fval, dfdx)
     maxdx = maxval(abs(x%x - xold%x))
     write(u, '(I0,15(",",ES24.16))') iter, f0, fval%x, &
          opt%get_residumax(), opt%get_residunorm(), maxdx
     write(k) x%x
  end do
  close(u); close(k)
  call MPI_Finalize(ierr)

contains

  subroutine evaluate(xv, f0, df0, g, dg)
    real(rp), intent(in) :: xv(n)
    real(rp), intent(out) :: f0, df0(n), g(m), dg(m, n)
    real(rp) :: h(n), Ie, ce, Me, xe
    integer :: i, j
    h = h_min + (h_max - h_min) * xv
    f0 = rho * b * Le * sum(h)
    df0 = rho * b * Le * (h_max - h_min)
    g(1) = sum(Delta / (b * h**3 / 12.0_rp) * (P / E)) / u_tip_max - 1.0_rp
    dg = 0.0_rp
    dg(1, :) = Delta * (P * (-36.0_rp) * (h_max - h_min) / (E * b)) / h**4 &
         / u_tip_max
    do i = 1, ncon
       j = idx(i)
       xe = Le * real(j - 1, rp)
       Ie = b * h(j)**3 / 12.0_rp
       ce = h(j) / 2.0_rp
       Me = P * (L_total - xe)
       g(1+i) = (Me * ce / Ie) / sigma_max - 1.0_rp
       dg(1+i, j) = Me * ((1.0_rp / (2.0_rp * Ie)) - &
            (ce * 3.0_rp * b * h(j)**2 / 12.0_rp) / (Ie**2)) * &
            (h_max - h_min) / sigma_max
    end do
  end subroutine evaluate

  ! Copied from tests/regression/mma/driver.f90
  subroutine fill_constraint_indices(stress_global_indices, num_constraints, &
       num_constraint_partitions, design_size)
    integer, intent(in) :: num_constraints, num_constraint_partitions, &
         design_size
    integer, intent(out) :: stress_global_indices(num_constraints)
    integer :: base_size, remainder_size, constraints_per_partition, &
         remainder_constraints, i, j, idx2, partition_start, partition_end, &
         partition_constraints, partition_size, cum_start
    constraints_per_partition = num_constraints / num_constraint_partitions
    remainder_constraints = mod(num_constraints, num_constraint_partitions)
    base_size = design_size / num_constraint_partitions
    remainder_size = mod(design_size, num_constraint_partitions)
    idx2 = 1
    cum_start = 1
    do i = 1, num_constraint_partitions
       if (i <= remainder_constraints) then
          partition_constraints = constraints_per_partition + 1
       else
          partition_constraints = constraints_per_partition
       endif
       if (i <= remainder_size) then
          partition_size = base_size + 1
       else
          partition_size = base_size
       endif
       partition_start = cum_start
       partition_end = partition_start + partition_size - 1
       do j = 1, partition_constraints
          if (idx2 > num_constraints) exit
          stress_global_indices(idx2) = min(partition_start + (j-1), &
               partition_end)
          idx2 = idx2 + 1
       end do
       cum_start = partition_end + 1
    end do
  end subroutine fill_constraint_indices
end program beam_driver
