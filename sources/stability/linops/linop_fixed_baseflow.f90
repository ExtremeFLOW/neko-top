module linop_fixed_baseflow
   use LightKrylov, only: abstract_vector_rdp
   use LightKrylov, only: wp => dp
   use linop_base, only: linop_base_t, linop_base_require_adjoint, &
        linop_base_require_baseflow, linop_base_require_lns
   use mpi_f08, only: MPI_WTIME
   use neko_ext, only: reset_adjoint
   use neko_stability_propagator, only: copy_perturb_to_state, &
        copy_state_to_perturb, copy_state_to_primal, fit_timestep_to_interval
   use neko_vector, only: state_vector_t
   use num_types, only: rp
   use simulation_perturb, only: simulation_adjoint_step, &
        simulation_lns_step, simulation_perturb_finalize, &
        simulation_perturb_init
   use time_step_controller, only: time_step_controller_t

   implicit none
   private

   type, extends(linop_base_t), public :: linop_fixed_baseflow_t
   contains
      private
      procedure, pass(self), public :: matvec => linop_fixed_baseflow_matvec
      procedure, pass(self), public :: rmatvec => linop_fixed_baseflow_rmatvec
   end type linop_fixed_baseflow_t

contains

   subroutine linop_fixed_baseflow_matvec(self, vec_in, vec_out)
      class(linop_fixed_baseflow_t), intent(inout) :: self
      class(abstract_vector_rdp), intent(in) :: vec_in
      class(abstract_vector_rdp), intent(inout) :: vec_out

      call linop_base_require_lns(self)
      select type (x => vec_in)
      type is (state_vector_t)
         select type (y => vec_out)
         type is (state_vector_t)
            call apply_fixed_lns(self, x, y)
         end select
      end select
   end subroutine linop_fixed_baseflow_matvec

   subroutine linop_fixed_baseflow_rmatvec(self, vec_in, vec_out)
      class(linop_fixed_baseflow_t), intent(inout) :: self
      class(abstract_vector_rdp), intent(in) :: vec_in
      class(abstract_vector_rdp), intent(inout) :: vec_out

      call linop_base_require_adjoint(self)
      select type (x => vec_in)
      type is (state_vector_t)
         select type (y => vec_out)
         type is (state_vector_t)
            call apply_fixed_adjoint(self, x, y)
         end select
      end select
   end subroutine linop_fixed_baseflow_rmatvec

   subroutine apply_fixed_lns(self, vec_in, vec_out)
      class(linop_fixed_baseflow_t), intent(inout) :: self
      type(state_vector_t), intent(in) :: vec_in
      type(state_vector_t), intent(inout) :: vec_out
      type(time_step_controller_t) :: dt_controller
      real(kind=wp) :: loop_start
      real(kind=rp) :: cfl, old_dt

      call linop_base_require_baseflow(self)
      call copy_state_to_primal(self%base_case, self%baseflow)

      old_dt = self%lns_case%time%dt
      call fit_timestep_to_interval(self%lns_case%time, &
           real(self%lns_case%time%end_time, kind=wp))

      call reset_adjoint(self%lns_case, self%base_case, 0, &
           self%lns_output_base_fname)
      call copy_state_to_perturb(self%lns_case, vec_in)

      call dt_controller%init(self%base_case%params)
      call simulation_perturb_init(self%lns_case, dt_controller)

      cfl = self%lns_case%fluid_pert%compute_cfl(self%lns_case%time%dt)
      loop_start = MPI_WTIME()
      do while (self%lns_case%time%t .lt. self%lns_case%time%end_time)
         call simulation_lns_step(self%lns_case, dt_controller, cfl, &
              loop_start)
      end do

      call vec_out%init_like(vec_in)
      call copy_perturb_to_state(self%base_case, self%lns_case, vec_out)

      call simulation_perturb_finalize(self%lns_case)
      self%lns_case%time%dt = old_dt
   end subroutine apply_fixed_lns

   subroutine apply_fixed_adjoint(self, vec_in, vec_out)
      class(linop_fixed_baseflow_t), intent(inout) :: self
      type(state_vector_t), intent(in) :: vec_in
      type(state_vector_t), intent(inout) :: vec_out
      type(time_step_controller_t) :: dt_controller
      real(kind=wp) :: loop_start
      real(kind=rp) :: cfl, final_time, old_dt

      call linop_base_require_baseflow(self)
      call copy_state_to_primal(self%base_case, self%baseflow)

      old_dt = self%adjoint_case%time%dt
      call fit_timestep_to_interval(self%adjoint_case%time, &
           real(self%adjoint_case%time%end_time, kind=wp))

      call reset_adjoint(self%adjoint_case, self%base_case, 0, &
           self%adjoint_output_base_fname)
      call copy_state_to_perturb(self%adjoint_case, vec_in)

      call dt_controller%init(self%base_case%params)
      call simulation_perturb_init(self%adjoint_case, dt_controller)

      cfl = self%adjoint_case%fluid_pert%compute_cfl(self%adjoint_case%time%dt)
      final_time = self%adjoint_case%time%end_time
      loop_start = MPI_WTIME()
      do while (self%adjoint_case%time%t .lt. &
           self%adjoint_case%time%end_time)
         call simulation_adjoint_step(self%adjoint_case, dt_controller, cfl, &
              loop_start, final_time)
      end do

      call vec_out%init_like(vec_in)
      call copy_perturb_to_state(self%base_case, self%adjoint_case, vec_out)

      call simulation_perturb_finalize(self%adjoint_case)
      self%adjoint_case%time%dt = old_dt
   end subroutine apply_fixed_adjoint

end module linop_fixed_baseflow
