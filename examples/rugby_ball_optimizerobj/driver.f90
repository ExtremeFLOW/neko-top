program usrneko
  ! use neko, only: neko_init, neko_solve, neko_finalize
  ! use case, only: case_t
  ! use adjoint_case, only: adjoint_case_t, adjoint_init, adjoint_free
  ! use simulation_adjoint, only: solve_adjoint

  ! use optimizer, only : optimizer_t
  use mma_optimizer, only : mma_optimizer_t 
  use steady_state_problem, only : steady_state_problem_t
  use topopt_design, only: topopt_design_t
  use num_types, only : rp
  
  use mask_ops, only: mask_exterior_const

  real(kind=rp) :: tolerance
  integer :: max_iter

  !> a problem type
  type(steady_state_problem_t) :: problem
  !> a design type
  type(topopt_design_t) :: design
  !> an optimizer (in this case mma)
  type(mma_optimizer_t) :: optimizer

  ! init the problem (base)
  call problem%init_base()

  ! init the design
  call design%init(problem%C%params, problem%C%fluid%c_Xh)

  ! init the problem, with the design
  call problem%init_design(design)

  ! init the optimizer
  call optimizer%init(problem)


  tolerance = 1.0e-3_rp
  max_iter = 100
  call optimizer%run(tolerance,max_iter)
!------------------------------------------------------------------------------

  ! deallocate(x_switch)
  call problem%free()
  call design%free()
  ! TODO
  call optimizer%free()

end program usrneko
