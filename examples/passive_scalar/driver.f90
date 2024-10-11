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

  ! these would live in "problem"
  !> and primal case                                                         
  type(case_t) :: C                                                  
  !> and adjoint case                                                        
  type(adjoint_case_t) :: adj
  ! the brinkman source terms
  type(simple_brinkman_source_term_t) :: forward_brinkman, adjoint_brinkman
  ! A term for the objective function
  type(adjoint_mixing_scalar_source_term_t) :: obj_source


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
    call obj_source%init_from_components( & 
         adj%scalar%f_Xh, & 
         C%scalar%s, &
         adj%scheme%c_Xh)
    ! append brinkman source term based on design
    call adj%scalar%source_term%add(obj_source)



!     call problem%compute()
!!------------------------------------------------------------------------------
      call neko_solve(C)

!     call problem%compute_sensitivity()
!!------------------------------------------------------------------------------
		call solve_adjoint(adj)

end program usrneko
