program usrneko
  use neko, only: neko_init, neko_solve, neko_finalize
  use case, only: case_t
  use user, only: user_setup
  use adjoint_case, only: adjoint_case_t, adjoint_init, adjoint_free
  use simulation_adjoint, only: solve_adjoint


  ! here are all the extra things i'll add to "problem"
  use source_term, only:  source_term_t
  use source_term_handler, only: source_term_handler_t
  use new_design, only: new_design_t
  use simple_brinkman_source_term, only: simple_brinkman_source_term_t
  use adjoint_minimum_dissipation_source_term, only: adjoint_minimum_dissipation_source_term_t
  use objective_function, only: objective_function_t
  use minimum_dissipation_objective_function, only: minimum_dissipation_objective_function_t


  ! use topopt, only: topopt_init, topopt_finalize, topopt_t

  ! I suppose in the future, all of this init stuff would be held inside of "problem"
  ! So that we make derived problem types, ie, include heat transfer etc.
  ! In which case the "problem" should own all the fluid and adjoint solvers etc,
  type(case_t) :: C
  type(adjoint_case_t) :: adj
  ! and this driver!
  ! For now we hard code everything.
  !
  ! So a problem would own a design
  type(new_design_t) :: design
  ! I would also argue that either the design or the "problem" should contain a way to dump all of 
  ! the mapped coeficients into both the forward and adjoint.
  ! Right now, we only have the brinkman term, so it's simple.
  ! But if we have conjugate heat transfer etc, there are more equations and more coefficients to map.
  ! So for now... again... hard coding just these two. 

  ! use derived type after!
  type(simple_brinkman_source_term_t) :: forward_brinkman
  type(simple_brinkman_source_term_t) :: adjoint_brinkman


  ! 
  ! And an objective type will also be necessary that contains
  ! TODO
  ! in the future we need a factory or so that generates these.
  ! for now, I'm picking a specific one
  type(minimum_dissipation_objective_function_t) :: objective_function

  ! - how the objective function is computed
  ! - the adjoint forcing
  ! - possibly also boundary condition modifications
  ! - a way to compute sensistivity to coeficcients, ie dF/d\chi in our case
  !
  ! (again I will hard code it)
  ! it would require 
  ! - the design (to compute chain rules for filter's etc)
  ! - the fluid/adjoint and potentially other equations like passive scalars
  ! - the objective 

  call user_setup(C%usr)
! In the future I guess all of this will be contained in the "problem" type.
! And this would be a 


! call problem%init()
!--------------------------------------------------------------------------------------------------------------------
  call neko_init(C)
  ! init the design field
  ! - init the initial design
  ! - filter it
  ! - map to brinkman
  call design%init(C%fluid%dm_Xh)

  ! init the objective function 
  ! - somehow append a user_check 
  ! TODO:
  ! Tim, I loved what you did with with the source term handler. I'm hoping when you get a chance you can do something
  ! similar with simulation components?
  ! as in, this kind of post processing isn't just one function, but a list of post processing modules that can be appended
  ! to a simulation (and then of course, we could append others to our adjoint!)
  !
  ! The thing is, because right now we're doing steady calculations, so the computations of:
  ! - The objective function:	performed at the end of the steady run
  ! - The sensitivity:			performed at the end of the adjoint run
  !
  ! but when we move to unsteady calculations we'll have: 
  ! - The objective function:	accumulated DURING the forward run
  ! - The sensitivity:			accumulated DURING the adjoint run
  !
  ! So they'll have to be simcomps that get appended to C and adj.
  ! I trust you can whip that up lickity split!





  ! manually init the forward brinkman term
  	! init the simple brinkman term
  	call forward_brinkman%init_from_components(C%fluid%f_x, C%fluid%f_y, C%fluid%f_z, design, &
  																C%fluid%u, C%fluid%v, C%fluid%w, &
  																C%fluid%c_Xh)
  	! append brinkman source term based on design
  	call C%fluid%source_term%add(forward_brinkman)

  ! append "user_check" for computing objective function

  

  ! initialize the adjoint
  call adjoint_init(adj, C)

  ! manually add to adjoint
  ! init the simple brinkman term for the adjoint
  	call adjoint_brinkman%init_from_components(adj%scheme%f_adj_x, adj%scheme%f_adj_y, adj%scheme%f_adj_z, design, &
  																adj%scheme%u_adj, adj%scheme%v_adj, adj%scheme%w_adj, &
  																adj%scheme%c_Xh)
  ! append brinkman source term based on design
  call adj%scheme%source_term%add(adjoint_brinkman)

  !! append a source term based on objective function
  call objective_function%init(C%fluid,adj%scheme)

  !! init the simple brinkman term for the adjoint
  !	call adjoint_forcing%init_from_components(adj%scheme%f_adj_x, adj%scheme%f_adj_y, adj%scheme%f_adj_z, &
  !																C%fluid%u, C%fluid%v, C%fluid%w, &
  !																adj%scheme%c_Xh)
  !! append brinkman source term based on design
  !call adj%scheme%source_term%add(adjoint_forcing)

!--------------------------------------------------------------------------------------------------------------------


  call neko_solve(C)
  ! TODO
  ! In the future, the objective_function_t will potentially include simulation components so that we can 
  ! accumulate the objective function during the run...
  ! here, we just compute it on the last step
  call objective_function%compute(C%fluid)
  print *, 'OBJECTIVE FUNCTION',  objective_function%objective_function_value


  ! TODO
  ! the steady simcomp only works for the forward, we either hardcode another one, or we have a forward adjoint flag
  ! or probably the smartest thing to do would be to accept a list of fields in the registry that will be checked...
  ! anyhow, I don't like the termination condition based on simulation time anyway...
  call solve_adjoint(adj)

  ! again, in the future, the objective_function_t will potentially include simulation components so that we can 
  ! accumulate the sensitivity during the run...
  ! here, we just compute it on the last step
  call objective_function%compute_sensitivity(C%fluid, adj%scheme)
  ! it would be nice to visualize this

  ! do the adjoint mapping
  call design%map_backward(objective_function%sensitivity_to_coefficient)



  call adjoint_free(adj)
  call neko_finalize(C)
end program usrneko
