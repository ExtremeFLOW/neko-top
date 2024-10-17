program usrneko
  use neko, only: neko_init, neko_solve, neko_finalize
  use case, only: case_t
  use user, only: user_setup
  use adjoint_case, only: adjoint_case_t, adjoint_init, adjoint_free
  use simulation_adjoint, only: solve_adjoint


  ! here are all the extra things i'll add to "problem"
  use source_term, only: source_term_t
  use source_term_handler, only: source_term_handler_t
  use topopt_design, only: topopt_design_t
  use simple_brinkman_source_term, only: simple_brinkman_source_term_t
  use adjoint_minimum_dissipation_source_term, &
       only: adjoint_minimum_dissipation_source_term_t
  use objective_function, only: objective_function_t
  use minimum_dissipation_objective_function, &
       only: minimum_dissipation_objective_function_t
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
  use adjoint_mixing_scalar_source_term, only: &
  adjoint_mixing_scalar_source_term_t
  use scalar_mixing_objective_function, only: &
  scalar_mixing_objective_function_t


  !> a problem type
  type(steady_state_problem_t) :: problem
  !> a design type
  type(topopt_design_t) :: design
  !> an optimizer (in this case mma)
  type(mma_t) :: optimizer

  ! these are some things needed for MMA/ work arrays (all these will become
  ! redundant when we do this properly)
  type(field_t), pointer :: wo1, wo2, wo3
  integer :: temp_indices(2)
  integer :: n, optimization_iteration
  real(kind=rp), dimension(1) :: fval
  real(kind=rp), allocatable :: x_switch(:)

  ! these would live in "problem"
  !> and primal case                                                         
  type(case_t) :: C                                                  
  !> and adjoint case                                                        
  type(adjoint_case_t) :: adj
  ! the brinkman source terms
  type(simple_brinkman_source_term_t) :: forward_brinkman, adjoint_brinkman
  ! the objective function
  type(scalar_mixing_objective_function_t) :: objective


  type(fld_file_output_t) :: output
  ! init the problem (base)
 !call problem%init_base()
!!------------------------------------------------------------------------------
    call user_setup(C%usr)
    ! initialize the primal
    call neko_init(C)
    ! initialize the adjoint
    call adjoint_init(adj, C)


 ! init the design
 ! call design%init(problem%C%params, problem%C%fluid%c_Xh)
!!------------------------------------------------------------------------------
  call design%init(C%params, C%fluid%c_Xh)


  ! init the problem, with the design
!  call problem%init_design(design)
!!------------------------------------------------------------------------------
   ! init the simple brinkman term for the forward problem
    call forward_brinkman%init_from_components( &
         C%fluid%f_x, C%fluid%f_y, C%fluid%f_z, &
         design, &
         C%fluid%u, C%fluid%v, C%fluid%w, &
         C%fluid%c_Xh)
    ! append brinkman source term to the forward problem
    call C%fluid%source_term%add(forward_brinkman)

    ! init the simple brinkman term for the adjoint
    call adjoint_brinkman%init_from_components( &
         adj%scheme%f_adj_x, adj%scheme%f_adj_y, &
         adj%scheme%f_adj_z, &
         design, &
         adj%scheme%u_adj, adj%scheme%v_adj, adj%scheme%w_adj, &
         adj%scheme%c_Xh)
    ! append brinkman source term based on design
    call adj%scheme%source_term%add(adjoint_brinkman)





    ! NOW here is where we would initialize our new objective functions!
    print *, "about to append the objective"
    call objective%init(design, C, adj)
    print *, "done"

    print *, "check apended one", allocated(adj%scalar%source_term%source_terms), &   
    size(adj%scalar%source_term%source_terms)

    ! init the sampler
    !---------------------------------------------------------
    ! TODO
    ! obviously when we do the mappings properly, to many coeficients, we'll
    ! also have to modify this
    ! for now:
    ! - forward (p,u,v,w)                      1,2,3,4           p,vx,vy,vz
    ! - design (\rho)                          5                 temperature
    ! - adjoint (u,v,w,p)                      6,7,8,9           s1,s2,s3,s4
    ! - mapped (\chi)                          10                s5
    ! - sensitivity (dF/d\chi and dC/d\chi)    11, 12            s6,s7
    ! - sensitivity (dF/d\rho and dC/d\rho)    13, 14            s8,s9

    call output%init(sp, 'optimization', 14)
    call output%fields%assign(1, C%fluid%p)
    call output%fields%assign(2, C%fluid%u)
    call output%fields%assign(3, C%fluid%v)
    call output%fields%assign(4, C%fluid%w)
    ! I don't know why these ones need assign_to_field?
    call output%fields%assign_to_field(5, design%design_indicator)
    call output%fields%assign(6, adj%scheme%u_adj)
    call output%fields%assign(7, adj%scheme%v_adj)
    call output%fields%assign(8, adj%scheme%w_adj)
    call output%fields%assign(9, adj%scheme%p_adj)
    call output%fields%assign_to_field(10, design%brinkman_amplitude)
    call output%fields%assign_to_field(11, &
         objective%sensitivity_to_coefficient)
    call output%fields%assign_to_field(12, design%sensitivity)
    call output%fields%assign(13, C%scalar%s)
    call output%fields%assign(14, adj%scalar%s_adj)



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
       0, 0.0_rp, (/0.0_rp/), (/100.0_rp/), (/0.0_rp/), wo1%x, wo2%x, C%params)
  !m, a0         a_i          c_i           d_i
  ! -------------------------------------------------------------------!
  !      Internal parameters for MMA                                   !
  !      Minimize  f_0(x) + a_0*z + sum( c_i*y_i + 0.5*d_i*(y_i)^2 )   !
  !    subject to  f_i(x) - a_i*z - y_i <= 0,  i = 1,...,m             !
  !                xmin_j <= x_j <= xmax_j,    j = 1,...,n             !
  !                z >= 0,   y_i >= 0,         i = 1,...,m             !
  ! -------------------------------------------------------------------!

  call neko_scratch_registry%relinquish_field(temp_indices)



  optimization_iteration = 1
  do while (optimization_iteration .lt. 100)

!     call problem%compute()
!!------------------------------------------------------------------------------
    print *, "check apended two", allocated(adj%scalar%source_term%source_terms), &   
    size(adj%scalar%source_term%source_terms)
      call neko_solve(C)
      call objective%compute(design,C)
      print *,'OBJECTIVE', objective%objective_function_value
    print *, "check apended three", allocated(adj%scalar%source_term%source_terms), &   
    size(adj%scalar%source_term%source_terms)

!     call problem%compute_sensitivity()
!!------------------------------------------------------------------------------
		call solve_adjoint(adj)

    call objective%compute_sensitivity(&
         design, C, adj)

    ! do the adjoint mapping
    call design%map_backward(&
         objective%sensitivity_to_coefficient)
     fval(1) = problem%volume_constraint%objective_function_value

     !call design%sample(real(optimization_iteration,kind=rp))
    call output%sample(real(optimization_iteration,kind=rp))

     associate(x => design%design_indicator%x, &
          df0dx=>design%sensitivity%x, &
          dfdx => design%sensitivity%x)
       ! TODO
       ! this is a really dumb way of handling the reshaping..
       if( .not. allocated(x_switch)) then
       allocate(x_switch(optimizer%get_n()))
       end if

       x_switch = reshape(x,[optimizer%get_n()])

       call cmult(dfdx,100.0_rp,n)
       fval(1) = fval(1)*100.0_rp
       ! (and also prints out some norms to make the trial and error and
       ! bit easier)

       call optimizer%mma_update_cpu( &
            optimization_iteration, &
            x_switch, &
            reshape(df0dx,[optimizer%get_n()]), &
            fval, &
            reshape(dfdx,[optimizer%get_m(), optimizer%get_n()]))

		 call copy(x,x_switch,optimizer%get_n())

       ! TODO
       ! do a KKT check and do a propper convergence check..
       call optimizer%kkt(x, df0dx, fval, dfdx)

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

     call reset(C)
     ! TODO
     ! reset for the adjoint
     call field_rzero(adj%scheme%u_adj)
     call field_rzero(adj%scheme%v_adj)
     call field_rzero(adj%scheme%w_adj)

     ! don't forget to unfreeze the fluid!
     C%fluid%freeze = .false.



  end do


  deallocate(x_switch)
  call design%free()
  call optimizer%free()

end program usrneko
