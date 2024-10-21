program usrneko
  use neko, only: neko_init, neko_solve, neko_finalize
  use case, only: case_t
  use adjoint_case, only: adjoint_case_t, adjoint_init, adjoint_free
  use simulation_adjoint, only: solve_adjoint


  ! here are all the extra things i'll add to "problem"
  use source_term, only: source_term_t
  use source_term_handler, only: source_term_handler_t
  use topopt_design, only: topopt_design_t
  use simple_brinkman_source_term, only: simple_brinkman_source_term_t
  use objective_function, only: objective_function_t
  use field, only:field_t
  use scratch_registry, only : neko_scratch_registry
  use num_types, only : rp, sp, dp, qp
  use field_math, only: field_rzero, field_rone, field_cmult
  use volume_constraint, only: volume_constraint_t
  use fld_file_output, only : fld_file_output_t
  use steady_state_problem, only : steady_state_problem_t
  use mma, only: mma_t
  use neko_ext, only: reset
  use math, only: copy, cmult


  !> a problem type
  type(steady_state_problem_t) :: problem
  !> a design type
  type(topopt_design_t) :: design
  !> an optimizer (in this case mma)
  type(mma_t) :: optimizer

  ! these are some things needed for MMA/ work arrays (all these will become
  ! redundant when we do this properly)
  type(field_t), pointer :: wo1, wo2, wo3
  integer :: temp_indices(3)
  integer :: n, optimization_iteration
  real(kind=rp), dimension(1) :: fval
  real(kind=rp), allocatable :: x_switch(:)




  ! init the problem (base)
  call problem%init_base()

  ! init the design
  call design%init(problem%C%params, problem%C%fluid%c_Xh)

  ! init the problem, with the design
  call problem%init_design(design)

  ! TODO
  ! call `optimizer%init(problem, design)`
!------------------------------------------------------------------------------
  ! this all has to be wrapped up again!
  ! so that we have an `optimizer_t` class, where it would pick out MMA.
  ! so this would be:


  ! - `problem` knows how many constraints etc there are
  ! - `design` knows that we're doing topology optimization
  ! here we've already read through the JSON to find the number of constraints

  ! here we hard code MMA
  ! Then we need to tell MMA how many constraints we need etc
  ! work arrays for xmin and xmax
  call neko_scratch_registry%request_field(wo1, temp_indices(1))
  call neko_scratch_registry%request_field(wo2, temp_indices(2))
  call field_rzero(wo1)
  call field_rone (wo2)
  ! obviously do this properly in the future...
  n = design%design_indicator%size()
  call optimizer%init_json(design%design_indicator%x, n, &
  !    m, a0         a_i          c_i           d_i
       1, 0.0_rp, [0.0_rp], [100.0_rp], [0.0_rp], wo1%x, wo2%x, &
       problem%C%params)
  ! -------------------------------------------------------------------!
  !      Internal parameters for MMA                                   !
  !      Minimize  f_0(x) + a_0*z + sum( c_i*y_i + 0.5*d_i*(y_i)^2 )   !
  !    subject to  f_i(x) - a_i*z - y_i <= 0,  i = 1,...,m             !
  !                xmin_j <= x_j <= xmax_j,    j = 1,...,n             !
  !                z >= 0,   y_i >= 0,         i = 1,...,m             !
  ! -------------------------------------------------------------------!

  call neko_scratch_registry%relinquish_field(temp_indices)
!------------------------------------------------------------------------------



  ! TODO
  ! call `optimizer%compute(problem, design)`
!------------------------------------------------------------------------------
  optimization_iteration = 1
  do while (optimization_iteration .lt. 100)

     ! compute objective function
     call problem%compute()

     ! in this case it's MMA so we need gradient information
     call problem%compute_sensitivity()

     ! now we have the optimizer act on the design field.

     fval(1) = problem%volume_constraint%objective_function_value
     associate(x => design%design_indicator%x, &
          df0dx => design%sensitivity%x, &
          dfdx => problem%volume_constraint%sensitivity_to_coefficient%x)
       ! TODO
       ! reshape everything and use the propper "update"

       ! TODO
       ! I'm CERTAIN the mass matrix comes in here, but I need to sit down
       ! with a pen and paper

       ! TODO
       ! for heuristic reasons, it's important to rescale the dfdx etc so
       ! they're within a reasonable range unfortunately this takes a bit of
       ! trial and error, but we should include a subroutine that
       ! allows us to rescale
       ! call cmult(df0dx,0.01_rp,n)
       call cmult(dfdx, 100.0_rp, n)
       fval(1) = fval(1)*100.0_rp
       ! (and also prints out some norms to make the trial and error and
       ! bit easier)

       call problem%sample(real(optimization_iteration, rp))

       !call optimizer%mma_update_cpu( &
       !     optimization_iteration, x, df0dx, fval, dfdx)

       ! TODO
       ! this is a really dumb way of handling the reshaping..
       if ( .not. allocated(x_switch) ) then
       allocate(x_switch(optimizer%get_n()))
       end if

       x_switch = reshape(x,[optimizer%get_n()])

       call optimizer%mma_update_cpu( &
            optimization_iteration, &
            x_switch, &
            reshape(df0dx,[optimizer%get_n()]), &
            fval, &
            reshape(dfdx,[optimizer%get_m(), optimizer%get_n()]))

       call copy(x, x_switch, optimizer%get_n())

       ! TODO
       ! do a KKT check and do a propper convergence check..
       call optimizer%mma_KKT_cpu(x, df0dx, fval, dfdx)

       ! TODO
       ! on the MMA side, we need residunorm public
       print *, 'KKT', optimizer%get_residunorm()

     end associate

     ! TODO
     ! instead of writing the objective function etc to the log file,
     ! we should come up with our own "optimization log" file, that writes out
     ! iterations, KKT, objective function values etc

     optimization_iteration = optimization_iteration + 1
     call design%map_forward()

     call reset(problem%C)
     ! TODO
     ! reset for the adjoint
     call field_rzero(problem%adj%scheme%u_adj)
     call field_rzero(problem%adj%scheme%v_adj)
     call field_rzero(problem%adj%scheme%w_adj)

     ! don't forget to unfreeze the fluid!
     problem%C%fluid%freeze = .false.



  end do
!------------------------------------------------------------------------------

  deallocate(x_switch)
  call problem%free()
  call design%free()
  ! TODO
  call optimizer%free()

end program usrneko
