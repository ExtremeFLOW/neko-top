submodule (mma) mma_cpu

contains
  module subroutine mma_gensub_cpu(this, iter, x, df0dx, fval, dfdx)
    ! ----------------------------------------------------- !
    ! Generate the approximation sub problem by computing   !
    ! the lower and upper asymtotes and the other necessary !
    ! parameters (alpha, beta, p0j, q0j, pij, qij, ...).    !
    ! ----------------------------------------------------- !
    class(mma_t), intent(inout) :: this
    real(kind=rp), dimension(this%n), intent(in) :: x
    real(kind=rp), dimension(this%n), intent(in) :: df0dx(:)
    real(kind=rp), dimension(this%m), intent(in) :: fval(:)
    real(kind=rp), dimension(this%m, this%n), intent(in) :: dfdx(:,:)
    integer, intent(in) :: iter
    integer :: i, j, ierr
    real(kind=rp), dimension(this%m) :: globaltmp_m

    if (iter .lt. 3) then
       do j = 1, this%n
          this%low(j) = x(j) - this%asyinit * (this%xmax(j) - &
               this%xmin(j))
          this%upp(j) = x(j) + this%asyinit * (this%xmax(j) - &
               this%xmin(j))
       end do
    else
       !Move asymptotes low and upp
       do j = 1, this%n
          if ((x(j) - this%xold1(j))*(this%xold1(j) - this%xold2(j)) &
               .lt. 0) then
             this%low(j) = x(j) - &
                  this%asydecr * (this%xold1(j) - this%low(j))
             this%upp(j) = x(j) + &
                  this%asydecr * (this%upp(j) - this%xold1(j))

          else if ((x(j) - this%xold1(j))* &
               (this%xold1(j) - this%xold2(j)) .gt. 0) then
             this%low(j) = x(j) - &
                  this%asyincr * (this%xold1(j) - this%low(j))
             this%upp(j) = x(j) + &
                  this%asyincr * (this%upp(j) - this%xold1(j))
          else
             this%low(j) = x(j) - (this%xold1(j) - this%low(j))
             this%upp(j) = x(j) + (this%upp(j) - this%xold1(j))
          end if

          ! setting a minimum and maximum for the low and upp
          ! asymptotes (eq3.9)
          this%low(j) = max(this%low(j), &
               x(j) - 10*(this%xmax(j) - this%xmin(j)))
          this%low(j) = min(this%low(j), &
               x(j) - 0.01*(this%xmax(j) - this%xmin(j)))

          this%upp(j) = min(this%upp(j), &
               x(j) + 10*(this%xmax(j) - this%xmin(j)))
          this%upp(j) = max(this%upp(j), &
               x(j) + 0.01*(this%xmax(j) - this%xmin(j)))
       end do
    end if
    ! we can move alpha and beta out of the following loop if needed as:
    ! this%alpha = max(this%xmin, this%low + &
    !     0.1*(this%x- this%low), this%x - 0.5*(this%xmax - this%xmin))
    ! this%beta = min(this%xmax, this%upp -  &
    !     0.1*(this%upp - this%x), this%x + 0.5*(this%xmax - this%xmin))
    do j = 1, this%n
       ! set the the bounds and coefficients for the approximation
       ! the move bounds (alpha and beta )are slightly more restrictive
       ! than low and upp. This is done based on eq(3.6)--eq(3.10).
       ! also check
       ! https://comsolyar.com/wp-content/uploads/2020/03/gcmma.pdf
       ! eq (2.8) and (2.9)
       this%alpha(j) = max(this%xmin(j), this%low(j) + &
            0.1*(x(j)- this%low(j)), &
            x(j) - 0.5*(this%xmax(j) - this%xmin(j)))
       this%beta(j) = min(this%xmax(j), this%upp(j) - &
            0.1*(this%upp(j) - x(j)), &
            x(j) + 0.5*(this%xmax(j) - this%xmin(j)))

       !Calculate p0j, q0j, pij, qij
       !where j = 1,2,...,n and i = 1,2,...,m  (eq(2.3)-eq(2.5))
       this%p0j(j) = (this%upp(j) - x(j))**2 * &
            (1.001*max(df0dx(j),0.0) + &
            0.001*max(-df0dx(j),0.0) + &
            (0.00001/(max(0.00001, &
            (this%xmax(j) - this%xmin(j))))))

       this%q0j(j) = (x(j) - this%low(j))**2 * &
            (0.001*max(df0dx(j),0.0) + &
            1.001*max(-df0dx(j),0.0) + &
            (0.00001/(max(0.00001, &
            (this%xmax(j) - this%xmin(j))))))

       do i = 1, this%m
          this%pij(i,j) = (this%upp(j) - x(j))**2 * &
               (1.001*max(dfdx(i,j),0.0) + &
               0.001*max(-dfdx(i,j),0.0) + &
               (0.00001/(max(0.00001, &
               (this%xmax(j) - this%xmin(j))))))
          this%qij(i,j) = (x(j) - this%low(j))**2 * &
               (0.001*max(dfdx(i, j), 0.0) + &
               1.001*max(-dfdx(i, j), 0.0) + &
               (0.00001/(max(0.00001, &
               (this%xmax(j) - this%xmin(j))))))
       end do
    end do

    !computing bi as defined in page 5
    this%bi = 0_rp
    do i = 1, this%m
       !MPI: here this%n is the global n
       do j = 1, this%n
          this%bi(i) = this%bi(i) + &
               this%pij(i,j) / (this%upp(j) - x(j)) + &
               this%qij(i,j) / (x(j) - this%low(j))
       end do
    end do


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!Showing that for double precision, bi will be different when!!!!!!!!
    !!!!!!!!!!!computed in parallel compare to sequential!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! this%bi = 0_rp
    ! longbi = 0.0
    ! do i = 1, this%m
    !     !MPI: here this%n is the global n
    !     do j = 1, this%n
    !         this%bi(i) = this%bi(i) + &
    !                     this%pij(i,j)/ (this%upp(j) - x(j)) + &
    !                     this%qij(i,j)/(x(j) - this%low(j))
    !         longbi(i) = longbi(i) + &
    !                     this%pij(i,j)/ (this%upp(j) - x(j)) + &
    !                     this%qij(i,j)/(x(j) - this%low(j))
    !     end do
    ! end do
    ! print *, "bi =  ", this%bi, "this%n = ", this%n
    ! print *, "longbi =  ", longbi
    ! ierr = 2160
    ! longbi = 0.0
    ! this%bi = 0_rp
    ! do i = 1, this%m
    !     do j = 1, ierr
    !         this%bi(i) = this%bi(i) + &
    !                     this%pij(i,j)/ (this%upp(j) - x(j)) + &
    !                     this%qij(i,j)/(x(j) - this%low(j))
    !         longbi(i) = longbi(i) + &
    !                     this%pij(i,j)/ (this%upp(j) - x(j)) + &
    !                     this%qij(i,j)/(x(j) - this%low(j))
    !     end do
    ! end do
    ! print *, "bi =  ", this%bi, "first batch(1-ierr)"
    ! print *, "longbi =  ", longbi, "first batch(1-ierr)"
    ! longbiglobal = longbi
    ! longbi = 0.0
    ! globaltmp_m = this%bi
    ! this%bi = 0_rp
    ! do i = 1, this%m
    !     do j = ierr+1, this%n
    !         this%bi(i) = this%bi(i) + &
    !                     this%pij(i,j)/ (this%upp(j) - x(j)) + &
    !                     this%qij(i,j)/(x(j) - this%low(j))
    !         longbi(i) = longbi(i) + &
    !                     this%pij(i,j)/ (this%upp(j) - x(j)) + &
    !                     this%qij(i,j)/(x(j) - this%low(j))
    !     end do
    ! end do
    ! print *, "bi =  ", this%bi, "second batch(ierr+1:end)"
    ! print *, "longbi =  ", longbi, "second batch(ierr+1:end)"
    ! print *, "bi =  ", this%bi+globaltmp_m, "first + second"
    ! print *, "longbi =  ", longbi+longbiglobal, "first + second"
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    globaltmp_m = 0.0_rp
    call MPI_Allreduce(this%bi, globaltmp_m, this%m, &
         mpi_real_precision, mpi_sum, neko_comm, ierr)
    this%bi = globaltmp_m - fval

  end subroutine mma_gensub_cpu

end submodule mma_cpu
