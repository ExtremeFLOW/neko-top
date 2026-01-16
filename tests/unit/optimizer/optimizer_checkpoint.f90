module optimizer_checkpointing_test_mod
  use vector, only: vector_t
  use vector_math, only: vector_glsubnorm
  use num_types, only: rp
  implicit none
  private
  public :: rmse

contains

  function rmse(v1, v2) result(result)
    type(vector_t), intent(in) :: v1, v2
    real(kind=rp) :: result

    result = sqrt(vector_glsubnorm(v1, v2)**2 / real(v1%size(), rp))

  end function rmse

end module optimizer_checkpointing_test_mod

program optimizer_checkpointing
  use mpi_f08, only: MPI_Init
  use, intrinsic :: iso_fortran_env, only: error_unit
  use simulation_m, only: simulation_t
  use design, only: design_t, design_factory
  use problem, only: problem_t
  use optimizer, only: optimizer_t, optimizer_factory

  use neko, only: neko_init, neko_finalize
  use json_module, only: json_file
  use json_utils, only: json_get
  use json_utils_ext, only: json_read_file
  use neko_top, only: neko_top_register_types
  use num_types, only: rp
  use profiler, only: profiler_start_region, profiler_end_region
  use utils, only: neko_error
  use vector, only: vector_t
  use registry, only: neko_registry
  use scratch_registry, only: neko_scratch_registry

  use optimizer_checkpointing_test_mod, only: rmse

  implicit none

  ! JSON related arguments
  integer :: argc
  character(len=1024) :: parameter_file
  type(json_file) :: parameters, design_parameters

  !> The simulation we are working with
  type(simulation_t) :: sim
  !> The design type
  class(design_t), allocatable :: des
  !> The problem type
  type(problem_t) :: prob
  !> The optimizer (in this case mma)
  class(optimizer_t), allocatable :: opt

  real(kind=rp) :: f_ref, f_half, f_restored
  type(vector_t) :: x_ref, x_half, x_restored

  ! Internal variables for the optimizer loop
  integer :: iter, iter0
  logical :: converged

  ! Set test parameters
  integer, parameter :: max_iter = 16
  real(kind=rp), parameter :: tol = 1.0e-8_rp

  ! -------------------------------------------------------------------------- !
  ! Initialize the Neko environment

  call neko_init()
  call neko_top_register_types()

  ! -------------------------------------------------------------------------- !
  ! Read the parameters file as the first terminal argument

  argc = command_argument_count()
  if (argc .lt. 1) call neko_error('Missing parameter file')
  call get_command_argument(1, parameter_file)

  ! Read the parameters file
  parameters = json_read_file(trim(parameter_file))
  call json_get(parameters, 'optimization.design', design_parameters)

  ! -------------------------------------------------------------------------- !
  ! Initialization of the components

  call sim%init(parameters)
  call design_factory(des, design_parameters, sim)
  call prob%init(parameters, des, sim)
  call optimizer_factory(opt, parameters, prob, des, sim)

  call x_ref%init(des%size())
  call x_half%init(des%size())
  call x_restored%init(des%size())

  ! -------------------------------------------------------------------------- !
  ! Run the optimization and save a checkpoint halfway

  iter0 = 0
  call opt%initialize(prob, des, sim)
  call opt%write(iter0, prob)

  do iter = iter0 + 1, max_iter
     converged = opt%step(iter, prob, des, sim)
     call opt%write(iter, prob)

     if (iter .eq. max_iter / 2 ) then
        call opt%save_checkpoint('optimizer_checkpoint.h5', iter, des, .true.)
        call prob%get_objective_value(f_half)
        call des%get_values(x_half)
     end if
  end do

  call prob%get_objective_value(f_ref)
  call des%get_values(x_ref)

  call opt%free()
  call prob%free()
  call des%free()
  call sim%free()
  deallocate(des, opt)

  ! Free global registries
  call neko_registry%free()
  call neko_scratch_registry%free()

  ! -------------------------------------------------------------------------- !
  ! Restore from checkpoint and verify

  call neko_registry%init()
  call neko_scratch_registry%init()

  call sim%init(parameters)
  call design_factory(des, design_parameters, sim)
  call prob%init(parameters, des, sim)
  call optimizer_factory(opt, parameters, prob, des, sim)

  call opt%restore_checkpoint('optimizer_checkpoint.h5', iter0, des)
  call opt%initialize(prob, des, sim)
  call opt%write(iter0, prob)

  call prob%get_objective_value(f_restored)
  call des%get_values(x_restored)

  if (iter0 .ne. max_iter / 2) then
     write(error_unit, *) 'Restored iteration: ', iter0, &
          ' Expected iteration: ', max_iter / 2
     call neko_error('Restored iteration does not match expected iteration')
  else if (abs(f_half - f_restored) .gt. tol) then
     write(error_unit, *) 'Objective at checkpoint: ', f_half, &
          ' Restored objective: ', f_restored, ' Difference: ', abs(f_half - f_restored)
     call neko_error('Restored objective does not match objective at checkpoint')
  else if (rmse(x_half, x_restored) .gt. tol) then
     write(error_unit, *) 'RMSE between original and restored design: ', &
          rmse(x_ref, x_restored)
     call neko_error('Restored design does not match original design')
  end if

  ! -------------------------------------------------------------------------- !
  ! Continue optimization from restored state

  do iter = iter0 + 1, max_iter
     converged = opt%step(iter, prob, des, sim)
     call opt%write(iter, prob)
  end do

  call prob%get_objective_value(f_restored)
  call des%get_values(x_restored)

  ! -------------------------------------------------------------------------- !
  ! Verify that the restored design matches the original one

  if (abs(f_ref - f_restored) .gt. tol) then
     write(error_unit, *) 'Final objective: ', f_restored, &
          ' Original objective: ', f_ref, ' Difference: ', abs(f_ref - f_restored)
     call neko_error('Final objective does not match original objective')
  else if (rmse(x_ref, x_restored) .gt. tol) then
     write(error_unit, *) 'RMSE between original and restored design: ', &
          rmse(x_ref, x_restored)
     call neko_error('Final design does not match original design')
  end if

  ! -------------------------------------------------------------------------- !
  ! Clean up the components

  call x_ref%free()
  call x_half%free()
  call x_restored%free()

  call opt%free()
  call prob%free()
  call des%free()
  call sim%free()

  if (allocated(des)) deallocate(des)
  if (allocated(opt)) deallocate(opt)

  ! Finalize the Neko environment
  call neko_finalize()

end program optimizer_checkpointing
