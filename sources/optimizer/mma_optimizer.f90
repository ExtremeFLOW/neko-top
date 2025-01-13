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

  contains
      ! Override the deferred methods
      procedure :: init => mma_optimizer_init
      procedure :: run => mma_optimizer_run
      procedure :: free => mma_optimizer_free

      procedure, pass(this) :: run_ss => mma_optimizer_run_steady_state_prob
  end type mma_optimizer_t

contains

  !> Initialize the MMA optimizer, associate it with a specific problem
  subroutine mma_optimizer_init(this, prob)
    class(mma_optimizer_t), intent(inout) :: this
    class(problem_t), intent(inout) :: prob

    ! Initialize MMA solver
    ! Check the type of the problem using select type
    select type(prob)
    type is (steady_state_problem_t)
      ! Now we know prob is of type steady_state_problem_t, assign the pointer
      ! this%steady_state_prob => prob
      print *, "Initializing mma_optimizer with steady_state_problem_t."
      ! mma_init_json( x, n, json, auto_scale, scale)
      call this%mma%init_json( prob%design%design_indicator%x, &
        prob%design%design_indicator%size(), prob%C%params, this%scale, &
        this%auto_scale)
      print *, "scale = ", this%scale
      print *, "auto_scale = ", this%auto_scale   
    class default
      !Unknown problem
      call neko_error('Unknown problem type in the mma_optimizer_init')
    end select
  end subroutine mma_optimizer_init

  ! Define the optimization loop for MMA
  subroutine mma_optimizer_run(this, prob, tolerance)
    class(mma_optimizer_t), intent(inout) :: this
    class(problem_t), intent(inout) :: prob
    real(kind=rp), intent(in) :: tolerance
    
    ! Check the type of the problem using select type
    select type(prob)
    type is (steady_state_problem_t)
      ! Now we know prob is of type steady_state_problem_t, call the run_ss
      call this%run_ss(prob, tolerance)
      ! steady_state_prob => prob
    class default
      !Unknown problem
      call neko_error('Unknown problem type in the mma_optimizer_run')
    end select
  end subroutine mma_optimizer_run

  subroutine mma_optimizer_run_steady_state_prob(this, prob, tolerance)
    class(mma_optimizer_t), intent(inout) :: this
    class(steady_state_problem_t), intent(inout) :: prob
    real(kind=rp), intent(in) :: tolerance
    integer :: max_iter
    integer :: iter, rank, ierr, nglobal
    real(kind=rp) :: scalingfactor

    character(len=1024) :: header_line
    type(csv_file_t) :: opt_outputf
    type(matrix_t) :: opt_data


    max_iter = this%mma%get_max_iter()
    ! call MPI_Comm_rank(neko_comm, rank, ierr)
    call MPI_Allreduce(this%mma%get_n(), nglobal, 1, &
        MPI_INTEGER, mpi_sum, neko_comm, ierr)

    !>initializing the scaling factor
    scalingfactor = 1.0_rp
    print *, "max_iter for the optimization loop = ", max_iter
    call opt_data%init(max_iter+1,6)

    call prob%compute()
    print *, "initial objective function value = " , &
      prob%volume_constraint%objective_function_value
    print *, "size(prob%design%design_indicator%x) = ", &
      size(prob%design%design_indicator%x)
    print *, "size(&
      &prob%volume_constraint%sensitivity_to_coefficient%x) = ",&
      size(&
     prob%volume_constraint%sensitivity_to_coefficient%x)
          

    !Writing the optimization data in a separate file
    open(1368, file = "optimization_data.txt", status = "replace")

    ! Loop to write labeled integer and real values

    associate(x => prob%design%design_indicator%x, &
      f0val => &
        prob%objective_function%objective_function_value, &
      fval => &
        prob%volume_constraint%objective_function_value, &
      df0dx => &
        prob%design%sensitivity%x, &
      dfdx => &
        prob%volume_constraint%sensitivity_to_coefficient%x)


 
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

      if (NEKO_BCKND_DEVICE .eq. 0) then
        call this%mma%mma_update_cpu( iter, x, df0dx, &
          reshape([fval*scalingfactor],[this%mma%get_m()]) , dfdx*scalingfactor)
      else
        write(stderr, *) "Device not supported in mma_optimizer.f90."
        error stop
      end if

      call prob%compute()
      call prob%compute_sensitivity()
      if (prob%design%if_mask) then
        call mask_exterior_const(&
          prob%volume_constraint%sensitivity_to_coefficient, &
          prob%design%optimization_domain, 0.0_rp)
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

      call prob%sample(real(iter, rp))


      call prob%design%map_forward()
      call reset(prob%C)
      ! TODO
      ! reset for the adjoint
      call field_rzero(prob%adj%scheme%u_adj)
      call field_rzero(prob%adj%scheme%v_adj)
      call field_rzero(prob%adj%scheme%w_adj)
      prob%C%fluid%freeze = .false.
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
  end subroutine mma_optimizer_run_steady_state_prob

  ! Free resources associated with the MMA optimizer
  subroutine mma_optimizer_free(this)
    class(mma_optimizer_t), intent(inout) :: this

    ! Free MMA-specific data
    call this%mma%free()
  end subroutine mma_optimizer_free

end module mma_optimizer

