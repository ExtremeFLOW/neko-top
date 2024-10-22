# Problem  {#problem}

Problem

aim to solve one step of : min F
				  					s.t. C <= 0

	proceedure sensitivity (both F and C)


	Contains:
	current design (design x)
	Objective - (object) general
		- determine objective function value (function) F
		- Compute dFdx
	
	instantiated from objective:
			easy_one - (derived from objective)
				- determine objective function value (function) F
				- Compute dFdx (from design variable)
				eg, volume constraints

			fluid/adjoint - (derived from objective)
			contains:
			fluid
				- init: compute objective
			adjoint
				- init: adjoint forcing

			proceedures
				- determine objective function value (function) F
				- Compute dFdx
					* determine adjoint forcing (source term)
					* compute sensitivity (proceedure) 
	
	We also need something similar for constraints

	objective F
	constraints C[:]


	Should be able to:
	compute F
	compute C

	compute dF/dx
	compute dC[i]/dx


	"find min F s.t. ..."





  type(case_t) :: C
  type(adjoint_obj) :: adj

  !
  call init_problem(P)
  !
  			call user_setup(C%usr)
  			call neko_init(C)


  			call init_problem_from_json(P)
  				! init design
  				! init F
  				 		call adj%init(C)
  				! init C
  				! loop all constraints
  					know C[i]

! compute F and dF/Dx 
! compute C[i] and dC[i]/dx
!------------------------
call solve()
			"steady"
  			call neko_solve(C)
  			call solve_adjoint(adj)
  			call sensitivity
  			or 
  			something else....
!----------------------

  call neko_finalize(C)

  call mma(P)




	
