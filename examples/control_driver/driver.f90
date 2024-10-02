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
  use field, only:field_t
  use scratch_registry, only : neko_scratch_registry
  use num_types, only : rp, sp, dp, qp
  use field_math, only: field_rzero, field_rone
  use volume_constraint, only: volume_constraint_t
  use fld_file_output, only : fld_file_output_t


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
  type(volume_constraint_t) :: volume_constraint

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


  ! these are some things needed for MMA/ work arrays (all these will become redundant when we do this properly)
  type(field_t), pointer :: wo1, wo2, wo3
  integer :: temp_indices(3)
  integer :: n, optimization_iteration
  real(kind=rp), dimension(1) :: fval

  ! Fuck the internal samplers, they don't make sense in this context, we should have our own.
  ! - design (\rho)
  ! - mapped (\chi)
  ! - forward (u,v,w,p)
  ! - adjoint (u,v,w,p)
  ! - sensitivity to coefficients (dF/d\chi and dC/d\chi)   (maybe this is redundant... but I want it for debugging)
  ! - sensitivity (dF/d\rho and dC/d\rho)
   type(fld_file_output_t) :: output


  call user_setup(C%usr)
! In the future I guess all of this will be contained in the "problem" type.
! And this would be a 


! call problem%init()
!###################################################################################################################
! init_from_json
!--------------------------------------------------------------------------------------------------------------------
! - here we would read what type of optimization we're doing, ie, passive scalar y/n?
! - steady/ unsteady? etc

! - also the type of objective function and number/type of constraints
!--------------------------------------------------------------------------------------------------------------------






!--------------------------------------------------------------------------------------------------------------------
! Here we init the solvers needed
!---------------------------------------------------------
  call neko_init(C)
  ! initialize the adjoint
  call adjoint_init(adj, C)

  ! init the design field
!---------------------------------------------------------
  ! - init the initial design
  ! - filter it
  ! - map to brinkman

  call design%init(C%fluid%dm_Xh)
  ! manually init the forward brinkman term
  	! init the simple brinkman term
  	call forward_brinkman%init_from_components(C%fluid%f_x, C%fluid%f_y, C%fluid%f_z, design, &
  																C%fluid%u, C%fluid%v, C%fluid%w, &
  																C%fluid%c_Xh)
  	! append brinkman source term based on design
  	call C%fluid%source_term%add(forward_brinkman)

  ! manually add to adjoint
  ! init the simple brinkman term for the adjoint
  	call adjoint_brinkman%init_from_components(adj%scheme%f_adj_x, adj%scheme%f_adj_y, adj%scheme%f_adj_z, design, &
  																adj%scheme%u_adj, adj%scheme%v_adj, adj%scheme%w_adj, &
  																adj%scheme%c_Xh)
  ! append brinkman source term based on design
  call adj%scheme%source_term%add(adjoint_brinkman)

  ! here we've already read through the JSON to find the number of constraints etc
!---------------------------------------------------------
  ! Then we need to tell MMA how many constraints we need etc
  ! work arrays for xmin and xmax
  call neko_scratch_registry%request_field(wo1, temp_indices(1))
  call neko_scratch_registry%request_field(wo2, temp_indices(2))
  call field_rzero(wo1)
  call field_rone (wo2)
  ! obviously do this properly in the future...
  n = design%design_indicator%size()
  call design%optimizer%init_json(design%design_indicator%x, n, &
   1, 0.0_rp, (/0.0_rp/), (/10.0_rp/), (/1.0_rp/), wo1%x, wo2%x, C%params)
  !m, a0         a_i          c_i           d_i
  ! -------------------------------------------------------------------!
  !      Internal parameters for MMA                                   !
  !      Minimize  f_0(x) + a_0*z + sum( c_i*y_i + 0.5*d_i*(y_i)^2 )   !
  !    subject to  f_i(x) - a_i*z - y_i <= 0,  i = 1,...,m             !
  !                xmin_j <= x_j <= xmax_j,    j = 1,...,n             !
  !                z >= 0,   y_i >= 0,         i = 1,...,m             !
  ! -------------------------------------------------------------------!

	call neko_scratch_registry%relinquish_field(temp_indices)


  ! TODO
  ! Note, Tim, while you're reading this I'm sure you can already see we need to unmangle a lot of this.
  ! for instance,
  ! FIRST we have to read what our objective is and all of our constraints
  ! now if an objective involve the fluid (which it will) THIS will tell us we need to init a fluid and an adjoint
  ! SECOND now we know how many coefficients we need to map in the design
  ! THIRD we can start adding adjoint forcing etc...
  ! so the order is a little fucked up here...
  ! technically 


  ! init the objective function 
!---------------------------------------------------------
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
  !
  ! in the future, the "problem" init will have already read out all the types of objective functions & constraints
  ! so in this place, we would already know from the JSON what objectives/constraints we need to init
  !
  ! So we need something akin to the "source_term_t" where we can have a list of them.
  !
  ! for this test we'll have 2
  ! minimum dissipation objective function
  call objective_function%init(design, C%fluid, adj%scheme)
  ! volume constraint
  call volume_constraint%init(design, C%fluid, adj%scheme)

  ! init the sampler
!---------------------------------------------------------
	! TODO
   ! obviously when we do the mappings properly, to many coeficients, we'll also have to modify this
   ! for now:
  ! - forward (p,u,v,w)								1,2,3,4				p,vx,vy,vz
  ! - design (\rho)									5						temperature
  ! - adjoint (u,v,w,p)								6,7,8,9				s1,s2,s3,s4
  ! - mapped (\chi)									10						s5
  ! - sensitivity (dF/d\chi and dC/d\chi)		11, 12				s6,s7
  ! - sensitivity (dF/d\rho and dC/d\rho)		13, 14				s8,s9

   call output%init(sp,'optimization',13)
   call output%fields%assign(1,  C%fluid%p)
   call output%fields%assign(2,  C%fluid%u)
   call output%fields%assign(3,  C%fluid%v)
   call output%fields%assign(4,  C%fluid%w)
   ! I don't know why these ones need assign_to_field?
   call output%fields%assign_to_field(5,  design%design_indicator)
   call output%fields%assign(6,  adj%scheme%u_adj)
   call output%fields%assign(7,  adj%scheme%v_adj)
   call output%fields%assign(8,  adj%scheme%w_adj)
   call output%fields%assign(9,  adj%scheme%p_adj)
   call output%fields%assign_to_field(10, design%brinkman_amplitude)
   call output%fields%assign_to_field(11, objective_function%sensitivity_to_coefficient)
   call output%fields%assign_to_field(12, volume_constraint%sensitivity_to_coefficient)
   call output%fields%assign_to_field(13, design%sensitivity)
   ! TODO
   ! I still haven't done the design%sensitivity as a field list!
   ! so it will eventually be 
   ! call this%output%fields%assign(13, design%sensitivity(1))
   ! call this%output%fields%assign(14, design%sensitivity(2))
   ! or something to this effect

!--------------------------------------------------------------------------------------------------------------------
!###################################################################################################################












! call problem%solve()
!###################################################################################################################
	optimization_iteration = 1
	do while (optimization_iteration.lt.100)


  call neko_solve(C)
  ! TODO
  ! In the future, the objective_function_t will potentially include simulation components so that we can 
  ! accumulate the objective function during the run...
  ! here, we just compute it on the last step
  ! TODO
  ! We would presumable have a list that holds all of objective functions and constraints, such that this would be a
  ! objectives%compute()
  call objective_function%compute(design, C%fluid)
  call volume_constraint%compute(design, C%fluid)
  print *, 'OBJECTIVE FUNCTION',  objective_function%objective_function_value
  print *, 'VOLUME CONSTRAINT', volume_constraint%objective_function_value, volume_constraint%volume


  ! TODO
  ! the steady simcomp only works for the forward, we either hardcode another one, or we have a forward adjoint flag
  ! or probably the smartest thing to do would be to accept a list of fields in the registry that will be checked...
  ! anyhow, I don't like the termination condition based on simulation time anyway...
  !
  call solve_adjoint(adj)

  ! again, in the future, the objective_function_t will potentially include simulation components so that we can 
  ! accumulate the sensitivity during the run...
  ! here, we just compute it on the last step

  ! TODO
  ! now that I look at this, we could have easily had design, fluid and adjoint passed in on init, and then have
  ! pointers pointing to what we're reffering to. 
  ! So then compute would take in t and t-step.
  ! this would probably be smarter than passing the whole thing down!!!
  ! TODO
  ! We would presumable have a list that holds all of objective functions and constraints, such that this would be a
  ! objectives%compute_sensitivity()
  ! and it would cycled through the list.
  call objective_function%compute_sensitivity(design, C%fluid, adj%scheme)
  call volume_constraint%compute_sensitivity(design, C%fluid, adj%scheme)
  ! it would be nice to visualize this

  ! do the adjoint mapping
  call design%map_backward(objective_function%sensitivity_to_coefficient)
  ! ok now you've fucked up the whole "list of sensitivity fields" aspect...
  ! we somehow need to populate the list 


	! this an appropriate place to sample, it's the current iteration + the sensitivity
	call output%sample(real(optimization_iteration,rp))
	! TODO
	! you've done this very incorrectly for the future,
	! to do this properly you need:
	! - a field list in the design holding sensitivity
	! - a list of mappings from rho -> whatever enters the PDE
	! - a list of adjoint mappings 
	! (the only reason this works is because we're considering the volume constraint on the unmapped field)
	! - a subroutine in design%update that compiles all this stuff into a way readable for MMA

	fval(1) = volume_constraint%objective_function_value
	associate(x => design%design_indicator%x, df0dx=>design%sensitivity%x, &
				dfdx => volume_constraint%sensitivity_to_coefficient%x)
	! TODO
	! reshape everything and use the propper "update"

	! TODO
	! I'm CERTAIN the mass matrix comes in here, but I need to sit down with a pen and paper

	! TODO
	! for heuristic reasons, it's important to rescale the dfdx etc so they're within a reasonable range
	! unfortunately this takes a bit of trial and error, but we should include a subroutine that 
	! allows us to rescale (and also prints out some norms to make the trial and error and bit easier)

	call design%optimizer%mma_update_cpu(optimization_iteration, x, df0dx, fval, dfdx) 

	! TODO
	! do a KKT check and do a propper convergence check..
	call design%optimizer%kkt(x, df0dx, fval, dfdx)

	! TODO
	! on the MMA side, we need residunorm public
	print *, 'KKT', design%optimizer%residunorm

	end associate

	! TODO
	! instead of writing the objective function etc to the log file, we should come up with our own
	! "optimization log" file, that writes out iterations, KKT, objective function values etc

	optimization_iteration = optimization_iteration + 1
	call design%map_forward()

	! TODO
	! use the propper restart functionality here!
	! but I can't find it on this branch :/
	call field_rzero(C%fluid%u)
	call field_rzero(C%fluid%v)
	call field_rzero(C%fluid%w)
	! don't forget to unfreeze the fluid!
	C%fluid%freeze = .false.

	call field_rzero(adj%scheme%u_adj)
	call field_rzero(adj%scheme%v_adj)
	call field_rzero(adj%scheme%w_adj)


	enddo
!###################################################################################################################





! call problem free
!###################################################################################################################
  call adjoint_free(adj)
  call neko_finalize(C)
  !
!###################################################################################################################
end program usrneko
