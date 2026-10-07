module fixed_point_system
   use LightKrylov, only: abstract_jacobian_linop_rdp, abstract_vector_rdp
   use LightKrylov, only: wp => dp
   use field_output, only: field_output_t
   use fld_file, only: fld_file_t
   use linop_fixed_baseflow, only: linop_fixed_baseflow_t
   use mpi_f08, only: MPI_WTIME
   use neko_ext, only: reset
   use neko_stability_propagator, only: copy_primal_to_state, &
        copy_state_to_primal, fit_timestep_to_interval
   use neko_vector, only: state_vector_t
   use num_types, only: neko_rp => rp, neko_sp => sp, rp
   use perturb_case, only: perturb_case_t
   use stability_solver_tolerances, only: solver_tolerance_guard_t
   use simulation, only: simulation_finalize, simulation_init, simulation_step
   use system_base, only: stability_system_free, stability_system_init, &
        stability_system_t
   use time_step_controller, only: time_step_controller_t

   implicit none

   type, extends(stability_system_t), public :: fixed_point_system_t
      logical :: write_outer_iterations = .false.
      integer :: outer_write_counter = 0
      character(len=:), allocatable :: outer_write_name
   contains
      private
      procedure, pass(self), public :: response => fixed_point_response
      procedure, pass(self), public :: init => fixed_point_system_init
      procedure, pass(self), public :: free => fixed_point_system_free
      procedure, pass(self), public :: enable_outer_iteration_output => &
           fixed_point_enable_outer_output
      procedure, pass(self), public :: disable_outer_iteration_output => &
           fixed_point_disable_outer_output
      procedure, pass(self) :: write_outer_state
   end type fixed_point_system_t

   type, extends(abstract_jacobian_linop_rdp), public :: fixed_point_jacobian_t
      type(linop_fixed_baseflow_t) :: propagator
   contains
      private
      procedure, pass(self), public :: matvec => fixed_point_jacobian_matvec
      procedure, pass(self), public :: rmatvec => fixed_point_jacobian_matvec
      procedure, pass(self), public :: init => fixed_point_jacobian_init
      procedure, pass(self), public :: free => fixed_point_jacobian_free
   end type fixed_point_jacobian_t

contains

   subroutine fixed_point_response(self, vec_in, vec_out, atol)
      class(fixed_point_system_t), intent(inout) :: self
      class(abstract_vector_rdp), intent(in) :: vec_in
      class(abstract_vector_rdp), intent(inout) :: vec_out
      real(wp), intent(in) :: atol

      select type (vec_in)
      type is (state_vector_t)
         select type (vec_out)
         type is (state_vector_t)
            call copy_state_to_primal(self%base_case, vec_in)
            if (self%write_outer_iterations) call self%write_outer_state()

            call propagate_nonlinear(self, vec_in, &
                 real(self%base_case%time%end_time, kind=wp), vec_out, &
                 solver_tolerance=atol)
            call vec_out%sub(vec_in)
         end select
      end select
   end subroutine fixed_point_response

   subroutine propagate_nonlinear(self, q0, time_interval, qT, fixed_nsteps, &
        solver_tolerance)
      class(fixed_point_system_t), intent(inout) :: self
      type(state_vector_t), intent(in) :: q0
      real(wp), intent(in) :: time_interval
      type(state_vector_t), intent(inout) :: qT
      integer, optional, intent(in) :: fixed_nsteps
      real(wp), optional, intent(in) :: solver_tolerance
      type(time_step_controller_t) :: dt_controller
      type(solver_tolerance_guard_t) :: tolerance_guard
      real(kind=rp) :: old_end_time, old_dt
      real(kind=wp) :: loop_start

      old_end_time = self%base_case%time%end_time
      old_dt = self%base_case%time%dt
      self%base_case%time%end_time = real(time_interval, kind=rp)
      call fit_timestep_to_interval(self%base_case%time, time_interval, &
           fixed_nsteps)

      call reset(self%base_case, 0, self%base_output_base_fname)
      call copy_state_to_primal(self%base_case, q0)
      if (present(solver_tolerance)) then
         call tolerance_guard%set(self%base_case, solver_tolerance)
      end if

      call dt_controller%init(self%base_case%params)
      call simulation_init(self%base_case, dt_controller)

      loop_start = MPI_WTIME()
      do while (.not. self%base_case%time%is_done())
         call simulation_step(self%base_case, dt_controller, loop_start)
      end do

      call qT%init_like(q0)
      call copy_primal_to_state(self%base_case, qT)

      call simulation_finalize(self%base_case)
      call tolerance_guard%restore()
      self%base_case%time%end_time = old_end_time
      self%base_case%time%dt = old_dt
   end subroutine propagate_nonlinear

   subroutine fixed_point_system_init(self, base_case)
      use case, only: case_t
      class(fixed_point_system_t), intent(inout) :: self
      type(case_t), target, intent(inout) :: base_case

      call stability_system_init(self, base_case)
      self%write_outer_iterations = .false.
      self%outer_write_counter = 0
   end subroutine fixed_point_system_init

   subroutine fixed_point_enable_outer_output(self, name)
      class(fixed_point_system_t), intent(inout) :: self
      character(len=*), intent(in) :: name

      if (allocated(self%outer_write_name)) deallocate(self%outer_write_name)
      self%outer_write_name = trim(name)
      self%outer_write_counter = 0
      self%write_outer_iterations = .true.
   end subroutine fixed_point_enable_outer_output

   subroutine fixed_point_disable_outer_output(self)
      class(fixed_point_system_t), intent(inout) :: self

      self%write_outer_iterations = .false.
      self%outer_write_counter = 0
      if (allocated(self%outer_write_name)) deallocate(self%outer_write_name)
   end subroutine fixed_point_disable_outer_output

   subroutine write_outer_state(self)
      class(fixed_point_system_t), intent(inout) :: self
      type(field_output_t) :: output
      integer :: idx

      if (.not. allocated(self%outer_write_name)) return

      idx = self%outer_write_counter
      call execute_command_line('mkdir -p viz')
      call output%init(trim(self%outer_write_name), 3, precision=neko_sp, &
           path='viz/', format='nek5000')
      select type (ft => output%file_%file_type)
      type is (fld_file_t)
         ft%skip_pressure = .true.
         ft%skip_velocity = .false.
         ft%skip_temperature = .true.
      end select
      call output%fields%assign_to_field(1, self%base_case%fluid%u)
      call output%fields%assign_to_field(2, self%base_case%fluid%v)
      call output%fields%assign_to_field(3, self%base_case%fluid%w)
      call output%set_counter(idx - 1)
      call output%sample(real(idx, kind=neko_rp))
      call output%free()
      self%outer_write_counter = self%outer_write_counter + 1
   end subroutine write_outer_state

   subroutine fixed_point_system_free(self)
      class(fixed_point_system_t), intent(inout) :: self

      call stability_system_free(self)
      call self%disable_outer_iteration_output()
   end subroutine fixed_point_system_free

   subroutine fixed_point_jacobian_matvec(self, vec_in, vec_out)
      class(fixed_point_jacobian_t), intent(inout) :: self
      class(abstract_vector_rdp), intent(in) :: vec_in
      class(abstract_vector_rdp), intent(inout) :: vec_out
      type(state_vector_t) :: tangent

      select type (base => self%X)
      type is (state_vector_t)
         select type (vec_in)
         type is (state_vector_t)
            select type (vec_out)
            type is (state_vector_t)
               call tangent%init_like(base)
               call self%propagator%set_baseflow(base)
               call self%propagator%matvec(vec_in, tangent)
               call vec_out%axpby(1.0_wp, tangent, 0.0_wp)
               call vec_out%sub(vec_in)
               call tangent%free()
            end select
         end select
      end select
   end subroutine fixed_point_jacobian_matvec

   subroutine fixed_point_jacobian_init(self, base_case, lns_case, adjoint_case)
      use case, only: case_t
      class(fixed_point_jacobian_t), intent(inout) :: self
      type(case_t), intent(inout), target :: base_case
      type(perturb_case_t), intent(inout), target :: lns_case
      type(perturb_case_t), optional, intent(inout), target :: adjoint_case

      if (present(adjoint_case)) then
         call self%propagator%init(base_case, lns_case, adjoint_case)
      else
         call self%propagator%init(base_case, lns_case)
      end if
   end subroutine fixed_point_jacobian_init

   subroutine fixed_point_jacobian_free(self)
      class(fixed_point_jacobian_t), intent(inout) :: self

      call self%propagator%free()
   end subroutine fixed_point_jacobian_free

end module fixed_point_system
