module mma_optimizer
  use optimizer, only: optimizer_t
  use steady_state_problem, only : steady_state_problem_t
  use mma, only: mma_t
  use problem, only: problem_t
  use num_types, only : rp
  use utils, only : neko_error

  !only to print nglobal when running in parallel
  use comm, only: neko_comm
  use mpi_f08, only: MPI_INTEGER, mpi_sum, MPI_Allreduce

  
  use neko_config, only: NEKO_BCKND_DEVICE
  ! Inclusions from external dependencies and standard libraries
  use, intrinsic :: iso_fortran_env, only: stderr => error_unit

  use math, only: copy, cmult
  use field_math, only: field_rzero
  use neko_ext, only: reset
  use mask_ops, only: mask_exterior_const


  use csv_file, only: csv_file_t
  use matrix, only: matrix_t

  implicit none
  private
  public :: mma_optimizer_t

  ! Concrete type for MMA optimizer
  type, extends(optimizer_t) :: mma_optimizer_t

      type(mma_t) :: mma

      !Pointer to the specific problem_type(steady_state_problem_t in this case)
      type(steady_state_problem_t), pointer :: steady_state_prob => null()
      ! type(unsteady_problem_t), pointer :: unsteady_prob => null()

  contains
      ! Override the deferred methods
      procedure :: init => mma_optimizer_init
      procedure :: run => mma_optimizer_run
      procedure :: free => mma_optimizer_free
  end type mma_optimizer_t

contains

  !> Initialize the MMA optimizer, associate it with a specific problem
  subroutine mma_optimizer_init(this, prob)
    class(mma_optimizer_t), intent(inout) :: this
    class(problem_t), target, intent(inout) :: prob
    real(kind=rp), allocatable :: xmax(:), xmin(:)


    ! Associate the problem with the optimizer
    this%prob => prob

    !setting xmax, xmin using neko_scratch_registry
    allocate(xmax(prob%design%design_indicator%size()))
    allocate(xmin(prob%design%design_indicator%size()))
    xmin=0_rp
    xmax=1_rp

    !>scaling fval and dfdx
    this%scale = 100
    this%auto_scale = .true.

    ! Initialize MMA solver
    ! Check the type of the problem using select type
    select type(prob)
    type is (steady_state_problem_t)
      ! Now we know prob is of type steady_state_problem_t, assign the pointer
      this%steady_state_prob => prob
      print *, "Initializing mma_optimizer with steady_state_problem_t."
      ! mma_init_json(this, x, n, m, a0, a, c, d, xmin, xmax, json)
      call this%mma%init_json( prob%design%design_indicator%x, &
        prob%design%design_indicator%size(), 1, 1.0_rp, [0.0_rp], [100.0_rp], &
        [0.0_rp], xmin, xmax, prob%C%params)


    class default
      !Unknown problem
      call neko_error('Unknown problem type in the mma_optimizer_init')
    end select
  end subroutine mma_optimizer_init

  ! Define the optimization loop for MMA
  subroutine mma_optimizer_run(this, tolerance, max_iter)
    class(mma_optimizer_t), intent(inout) :: this
    integer, intent(in) :: max_iter
    real(kind=rp), intent(in) :: tolerance
    integer :: iter, rank, ierr, nglobal
    real(kind=rp) :: scalingfactor

    character(len=1024) :: header_line
    type(csv_file_t) :: opt_outputf
    type(matrix_t) :: opt_data


    ! call MPI_Comm_rank(neko_comm, rank, ierr)
    call MPI_Allreduce(this%mma%get_n(), nglobal, 1, &
        MPI_INTEGER, mpi_sum, neko_comm, ierr)

    !>initializing the scaling factor
    scalingfactor = 1_rp

    call opt_data%init(max_iter+1,6)

    !check if there is a drived type for prob
    !if not, then prob is of abstract type problem_t and therefore we get errors
    !later on we will add other types of problem here as well:
    ! if (.not. (associated(this%steady_state_prob) .and. &
    !           (associated(this%unsteady_prob) .and. ... &
    !     ......  ))
    if (.not. associated(this%steady_state_prob)) then
      call neko_error('Undefined problem type initialized in mma_optimizer_run')
    endif

    call this%prob%compute()
    print *, "initial objective function value=" , &
      this%steady_state_prob%volume_constraint%objective_function_value
    print *, "size(this%prob%design%design_indicator%x)=", &
      size(this%steady_state_prob%design%design_indicator%x)
    print *, "size(&
      &this%prob%volume_constraint%sensitivity_to_coefficient%x)=",&
      size(&
     this%steady_state_prob%volume_constraint%sensitivity_to_coefficient%x)
          

    !Writing the optimization data in a separate file
    open(1368, file="optimization_data.txt", status="replace")

    ! Loop to write labeled integer and real values

    associate(x => this%steady_state_prob%design%design_indicator%x, &
      f0val => &
        this%steady_state_prob%objective_function%objective_function_value, &
      fval => &
        this%steady_state_prob%volume_constraint%objective_function_value, &
      df0dx => &
        this%steady_state_prob%design%sensitivity%x, &
      dfdx => &
        this%steady_state_prob%volume_constraint%sensitivity_to_coefficient%x)


 
    write(1368, '("n=", I10, ", m=", I10)') nglobal, this%mma%get_m()
    write(1368, '("iter=", I3, ", f0val=", ES25.17, &
      & ", fval(1)=", ES25.17, ", KKTmax=", ES25.17, ", tolerance=", ES25.17,&
      & ", KKTnorm2=", ES25.17)') 0, f0val, fval, &
      & this%mma%get_residumax(), tolerance, this%mma%get_residunorm()

    opt_data%x(1,1) = 0
    opt_data%x(1,2) = f0val
    opt_data%x(1,3) = fval
    opt_data%x(1,4) = this%mma%get_residumax()
    opt_data%x(1,5) = this%mma%get_residunorm()


    do iter = 1, max_iter
      if (this%mma%get_residumax() .lt. tolerance) exit
      !Scaling
      if (this%auto_scale .eqv. .true.) then
        scalingfactor = abs(this%scale/fval)
      else
        scalingfactor = abs(this%scale)
      end if
      ! call cmult(dfdx, 100.0_rp, n)
      ! fval(1) = fval(1)*100.0_rp

      ! mma_update_cpu(this, iter, x, df0dx, fval, dfdx)
      ! call this%mma%mma_update_cpu( iter, x, df0dx, &
      !   reshape([fval*100.0_rp],[this%mma%get_m()]) , dfdx*100.0_rp)


      if (NEKO_BCKND_DEVICE .eq. 0) then
        call this%mma%mma_update_cpu( iter, x, df0dx, &
          reshape([fval*scalingfactor],[this%mma%get_m()]) , dfdx*scalingfactor)
      else
        write(stderr, *) "Device not supported in mma_optimizer.f90."
        error stop
      end if
      ! call this%mma%update( iter, x, df0dx, &
      !   reshape([fval*scalingfactor],[this%mma%get_m()]) , dfdx*scalingfactor)


      call this%steady_state_prob%compute()
      call this%steady_state_prob%compute_sensitivity()
      if (this%steady_state_prob%design%if_mask) then
        call mask_exterior_const(&
          this%steady_state_prob%volume_constraint%sensitivity_to_coefficient, &
          this%steady_state_prob%design%optimization_domain, 0.0_rp)
      end if

      call this%mma%KKT(x,df0dx,reshape([fval],[this%mma%get_m()]),dfdx)

      print *, 'iter=', iter,&
        '-------,f0val= ', f0val, ',   fval= ', fval, &
        ',  KKTmax=', this%mma%get_residumax(), ', KKTnorm2=',&
        this%mma%get_residunorm()

      write(1368, '("iter=", I3, ", f0val=", ES25.17, &
        & ", fval(1)=", ES25.17, ", KKTmax=", ES25.17, &
        & ", KKTnorm2=", ES25.17, ", scalingfactor=", ES25.17)') iter, &
        f0val, fval, this%mma%get_residumax(), this%mma%get_residunorm(), &
        scalingfactor
      ! Flush the buffer to write the data during the run
      flush(1368)


      opt_data%x(iter+1,1) = iter
      opt_data%x(iter+1,2) = f0val
      opt_data%x(iter+1,3) = fval
      opt_data%x(iter+1,4) = this%mma%get_residumax()
      opt_data%x(iter+1,5) = this%mma%get_residunorm()
      opt_data%x(iter+1,6) = scalingfactor

      call this%steady_state_prob%sample(real(iter, rp))


      call this%steady_state_prob%design%map_forward()
      call reset(this%steady_state_prob%C)
      ! TODO
      ! reset for the adjoint
      call field_rzero(this%steady_state_prob%adj%scheme%u_adj)
      call field_rzero(this%steady_state_prob%adj%scheme%v_adj)
      call field_rzero(this%steady_state_prob%adj%scheme%w_adj)
      this%steady_state_prob%C%fluid%freeze = .false.
    end do
    end associate


    close(1368)

    write(header_line, '(A,A,I10,A,I10)') "iter,f0val,fval1,KKTmax,KKTnorm2,", &
      "Scalingfactor,n=",&
      nglobal, ", m=", this%mma%get_m()
    call opt_outputf%init("optimization_data.csv")
    call opt_outputf%set_header(header_line)
    call opt_outputf%write(opt_data)


    ! Final state after optimization
    print*, "MMA Optimization completed after", iter-1, "iterations."

  end subroutine mma_optimizer_run

  ! Free resources associated with the MMA optimizer
  subroutine mma_optimizer_free(this)
    class(mma_optimizer_t), intent(inout) :: this

    ! Free MMA-specific data
    call this%mma%free()
  end subroutine mma_optimizer_free

end module mma_optimizer

