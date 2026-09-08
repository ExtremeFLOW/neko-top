!> @file topopt.f90
!! @copyright
!! Copyright (c) 2024-2025, The Neko-TOP Authors
!! All rights reserved.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions
!! are met:
!!
!!   * Redistributions of source code must retain the above copyright
!!     notice, this list of conditions and the following disclaimer.
!!
!!   * Redistributions in binary form must reproduce the above
!!     copyright notice, this list of conditions and the following
!!     disclaimer in the documentation and/or other materials provided
!!     with the distribution.
!!
!!   * Neither the name of the authors nor the names of its
!!     contributors may be used to endorse or promote products derived
!!     from this software without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
!! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
!! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
!! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
!! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
!! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
!! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!! POSSIBILITY OF SUCH DAMAGE.

program topopt
  use neko, only: neko_init, neko_finalize, neko_job_info
  use simulation_m, only: simulation_t
  use design, only: design_t, design_factory
  use problem, only: problem_t
  use optimizer, only: optimizer_t, optimizer_factory

  use json_module, only: json_file
  use json_utils, only: json_get
  use utils, only: neko_error
  use json_utils_ext, only: json_read_file
  use neko_top, only: neko_top_register_types
  use continuation_scheduler, only: nekotop_continuation

  implicit none

  ! JSON related arguments
  integer :: argc
  character(len=256) :: parameter_file
  type(json_file) :: parameters, design_parameters
  character(10) :: time
  character(8) :: date

  !> The simulation we are working with
  type(simulation_t), target :: sim
  !> The design type
  class(design_t), allocatable :: des
  !> The problem type
  type(problem_t) :: prob
  !> The optimizer (in this case mma)
  class(optimizer_t), allocatable :: opt

  ! -------------------------------------------------------------------------- !
  ! Initialize the Neko environment

  call date_and_time(time = time, date = date)
  call neko_init()
  call neko_job_info(date, time)
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

  ! initialize the global continuation_scheduler object (nekotop_continuation)
  call nekotop_continuation%init(parameters)

  ! initialize the simulation
  call sim%init(parameters)

  ! initialize the design
  call design_factory(des, design_parameters, sim)

  ! initialize the problem
  call prob%init(parameters, des, sim)

  ! initialize the optimizer
  call optimizer_factory(opt, parameters, prob, des, sim)

  ! -------------------------------------------------------------------------- !
  ! Execute the optimization

  call opt%run(prob, des, sim)

  ! -------------------------------------------------------------------------- !
  ! Clean up the components

  call opt%free()
  call prob%free()
  call des%free()
  call sim%free()
  call nekotop_continuation%free()

  if (allocated(des)) deallocate(des)
  if (allocated(opt)) deallocate(opt)

  ! Finalize the Neko environment
  call neko_finalize()

end program topopt
