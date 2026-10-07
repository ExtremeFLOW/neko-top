module linop_transient_baseflow
   use LightKrylov, only: abstract_vector_rdp
   use LightKrylov, only: wp => dp
   use case, only: case_t
   use json_file_module, only: json_file
   use json_utils, only: json_get
   use linop_base, only: linop_base_t, linop_base_init, &
        linop_base_require_adjoint, linop_base_require_baseflow, &
        linop_base_require_lns, linop_base_store_baseflow, linop_base_free
   use mpi_f08, only: MPI_WTIME
   use neko_ext, only: reset, reset_adjoint
   use neko_stability_propagator, only: copy_perturb_to_state, &
        copy_primal_to_state, copy_state_to_perturb, copy_state_to_primal, &
        fit_timestep_to_interval
   use neko_vector, only: state_vector_t
   use num_types, only: rp
   use perturb_case, only: perturb_case_t
   use simulation, only: simulation_finalize, simulation_init, simulation_step
   use simulation_checkpoint, only: simulation_checkpoint_t
   use simulation_perturb, only: simulation_adjoint_step, &
        simulation_lns_step, simulation_perturb_finalize, &
        simulation_perturb_init
   use time_step_controller, only: time_step_controller_t
   use utils, only: neko_error, neko_warning

   implicit none
   private

   type, extends(linop_base_t), public :: linop_transient_baseflow_t
      type(simulation_checkpoint_t) :: checkpoint
      real(wp) :: time_interval = 0.0_wp
      integer :: fixed_step_count = 0
      integer :: n_timesteps = 0
      logical :: checkpoint_initialized = .false.
   contains
      private
      procedure, pass(self), public :: init => linop_transient_baseflow_init
      procedure, pass(self), public :: matvec => linop_transient_baseflow_matvec
      procedure, pass(self), public :: rmatvec => linop_transient_baseflow_rmatvec
      procedure, pass(self), public :: set_baseflow => &
           linop_transient_baseflow_set_baseflow
      procedure, pass(self), public :: free => linop_transient_baseflow_free
   end type linop_transient_baseflow_t

contains

   subroutine linop_transient_baseflow_init(self, base_case, lns_case, &
        adjoint_case)
      class(linop_transient_baseflow_t), intent(inout) :: self
      type(case_t), target, intent(inout) :: base_case
      type(perturb_case_t), optional, target, intent(inout) :: lns_case
      type(perturb_case_t), optional, target, intent(inout) :: adjoint_case
      type(json_file) :: checkpoint_params

      call linop_base_init(self, base_case, lns_case, adjoint_case)
      self%time_interval = 0.0_wp
      self%fixed_step_count = 0
      self%n_timesteps = 0
      self%checkpoint_initialized = .false.

      if ("checkpoints" .in. base_case%params) then
         call json_get(base_case%params, 'checkpoints', checkpoint_params)
         call self%checkpoint%init(base_case, checkpoint_params)
         self%checkpoint_initialized = .true.
      else
         call neko_warning('Transient-baseflow linop initialized without ' // &
              'checkpoints; rmatvec will be unavailable')
      end if
   end subroutine linop_transient_baseflow_init

   subroutine linop_transient_baseflow_set_baseflow(self, baseflow, &
        time_interval, fixed_step_count)
      class(linop_transient_baseflow_t), intent(inout) :: self
      type(state_vector_t), intent(in) :: baseflow
      real(wp), optional, intent(in) :: time_interval
      integer, optional, intent(in) :: fixed_step_count

      call linop_base_store_baseflow(self, baseflow)
      self%time_interval = 0.0_wp
      self%fixed_step_count = 0
      if (present(time_interval)) self%time_interval = time_interval
      if (present(fixed_step_count)) self%fixed_step_count = fixed_step_count
   end subroutine linop_transient_baseflow_set_baseflow

   subroutine linop_transient_baseflow_matvec(self, vec_in, vec_out)
      class(linop_transient_baseflow_t), intent(inout) :: self
      class(abstract_vector_rdp), intent(in) :: vec_in
      class(abstract_vector_rdp), intent(inout) :: vec_out

      call linop_base_require_lns(self)
      select type (x => vec_in)
      type is (state_vector_t)
         select type (y => vec_out)
         type is (state_vector_t)
            call apply_lns(self, x, y)
         end select
      end select
   end subroutine linop_transient_baseflow_matvec

   subroutine linop_transient_baseflow_rmatvec(self, vec_in, vec_out)
      class(linop_transient_baseflow_t), intent(inout) :: self
      class(abstract_vector_rdp), intent(in) :: vec_in
      class(abstract_vector_rdp), intent(inout) :: vec_out

      call linop_base_require_adjoint(self)
      if (.not. self%checkpoint_initialized) then
         call neko_error('Transient-baseflow adjoint requires checkpoints')
      end if

      select type (x => vec_in)
      type is (state_vector_t)
         select type (y => vec_out)
         type is (state_vector_t)
            call apply_adjoint(self, x, y)
         end select
      end select
   end subroutine linop_transient_baseflow_rmatvec

   subroutine apply_lns(self, vec_in, vec_out)
      class(linop_transient_baseflow_t), intent(inout) :: self
      type(state_vector_t), intent(in) :: vec_in
      type(state_vector_t), intent(inout) :: vec_out
      type(state_vector_t) :: mapped
      type(time_step_controller_t) :: primal_dt_controller
      type(time_step_controller_t) :: perturb_dt_controller
      real(kind=rp) :: old_base_end, old_lns_end, old_base_dt, old_lns_dt
      real(kind=rp) :: cfl
      real(kind=wp) :: interval, loop_start

      call linop_base_require_baseflow(self)
      interval = propagation_interval(self)

      old_base_end = self%base_case%time%end_time
      old_lns_end = self%lns_case%time%end_time
      old_base_dt = self%base_case%time%dt
      old_lns_dt = self%lns_case%time%dt

      self%base_case%time%end_time = real(interval, kind=rp)
      self%lns_case%time%end_time = real(interval, kind=rp)
      if (self%fixed_step_count > 0) then
         call fit_timestep_to_interval(self%base_case%time, interval, &
              self%fixed_step_count)
      else
         call fit_timestep_to_interval(self%base_case%time, interval)
      end if
      self%lns_case%time%dt = self%base_case%time%dt
      self%lns_case%time%dtlag = self%lns_case%time%dt

      call reset(self%base_case, 0, self%base_output_base_fname)
      call reset_adjoint(self%lns_case, self%base_case, 0, &
           self%lns_output_base_fname)
      call copy_state_to_primal(self%base_case, self%baseflow)
      call copy_state_to_perturb(self%lns_case, vec_in)

      call primal_dt_controller%init(self%base_case%params)
      call perturb_dt_controller%init(self%base_case%params)
      call simulation_init(self%base_case, primal_dt_controller)
      call simulation_perturb_init(self%lns_case, perturb_dt_controller)

      cfl = self%lns_case%fluid_pert%compute_cfl(self%lns_case%time%dt)
      loop_start = MPI_WTIME()
      self%n_timesteps = 0
      do while (.not. self%base_case%time%is_done())
         self%n_timesteps = self%n_timesteps + 1
         call simulation_step(self%base_case, primal_dt_controller, loop_start)
         call simulation_lns_step(self%lns_case, perturb_dt_controller, &
              cfl, loop_start)
      end do

      call mapped%init_like(self%baseflow)
      call vec_out%init_like(vec_in)
      call copy_primal_to_state(self%base_case, mapped)
      call copy_perturb_to_state(self%base_case, self%lns_case, vec_out)

      call simulation_perturb_finalize(self%lns_case)
      call simulation_finalize(self%base_case)

      self%base_case%time%end_time = old_base_end
      self%lns_case%time%end_time = old_lns_end
      self%base_case%time%dt = old_base_dt
      self%lns_case%time%dt = old_lns_dt
      call mapped%free()
   end subroutine apply_lns

   subroutine apply_adjoint(self, vec_in, vec_out)
      class(linop_transient_baseflow_t), intent(inout) :: self
      type(state_vector_t), intent(in) :: vec_in
      type(state_vector_t), intent(inout) :: vec_out
      real(kind=rp) :: old_base_end, old_adj_end, old_base_dt, old_adj_dt
      real(kind=wp) :: interval

      call linop_base_require_baseflow(self)
      interval = propagation_interval(self)

      old_base_end = self%base_case%time%end_time
      old_adj_end = self%adjoint_case%time%end_time
      old_base_dt = self%base_case%time%dt
      old_adj_dt = self%adjoint_case%time%dt

      call run_base_with_checkpoints(self, interval)
      call run_checkpointed_adjoint(self, vec_in, vec_out, interval)

      self%base_case%time%end_time = old_base_end
      self%adjoint_case%time%end_time = old_adj_end
      self%base_case%time%dt = old_base_dt
      self%adjoint_case%time%dt = old_adj_dt
   end subroutine apply_adjoint

   subroutine run_base_with_checkpoints(self, interval)
      class(linop_transient_baseflow_t), intent(inout) :: self
      real(wp), intent(in) :: interval
      type(time_step_controller_t) :: dt_controller
      real(kind=wp) :: loop_start

      self%base_case%time%end_time = real(interval, kind=rp)
      if (self%fixed_step_count > 0) then
         call fit_timestep_to_interval(self%base_case%time, interval, &
              self%fixed_step_count)
      else
         call fit_timestep_to_interval(self%base_case%time, interval)
      end if

      call reset(self%base_case, 0, self%base_output_base_fname)
      call self%checkpoint%reset()
      call copy_state_to_primal(self%base_case, self%baseflow)

      call dt_controller%init(self%base_case%params)
      call simulation_init(self%base_case, dt_controller)

      loop_start = MPI_WTIME()
      self%n_timesteps = 0
      do while (.not. self%base_case%time%is_done())
         self%n_timesteps = self%n_timesteps + 1
         call simulation_step(self%base_case, dt_controller, loop_start)
         call self%checkpoint%save(self%base_case)
      end do

      call simulation_finalize(self%base_case)
   end subroutine run_base_with_checkpoints

   subroutine run_checkpointed_adjoint(self, vec_in, vec_out, interval)
      class(linop_transient_baseflow_t), intent(inout) :: self
      type(state_vector_t), intent(in) :: vec_in
      type(state_vector_t), intent(inout) :: vec_out
      real(wp), intent(in) :: interval
      type(time_step_controller_t) :: dt_controller
      real(kind=wp) :: loop_start
      real(kind=rp) :: cfl, total_time
      integer :: i

      self%adjoint_case%time%end_time = real(interval, kind=rp)
      self%adjoint_case%time%dt = self%base_case%time%dt
      self%adjoint_case%time%dtlag = self%adjoint_case%time%dt

      call reset_adjoint(self%adjoint_case, self%base_case, 0, &
           self%adjoint_output_base_fname)
      call copy_state_to_perturb(self%adjoint_case, vec_in)

      call dt_controller%init(self%base_case%params)
      call simulation_perturb_init(self%adjoint_case, dt_controller)

      cfl = self%adjoint_case%fluid_pert%compute_cfl(self%adjoint_case%time%dt)
      total_time = self%n_timesteps * self%adjoint_case%time%dt
      loop_start = MPI_WTIME()
      do i = self%n_timesteps, 1, -1
         call self%checkpoint%restore(self%base_case, i)
         call simulation_adjoint_step(self%adjoint_case, dt_controller, cfl, &
              loop_start, total_time)
      end do

      call vec_out%init_like(vec_in)
      call copy_perturb_to_state(self%base_case, self%adjoint_case, vec_out)

      call simulation_perturb_finalize(self%adjoint_case)
   end subroutine run_checkpointed_adjoint

   real(wp) function propagation_interval(self) result(interval)
      class(linop_transient_baseflow_t), intent(in) :: self

      interval = self%time_interval
      if (interval <= 0.0_wp) then
         interval = real(self%base_case%time%end_time, kind=wp)
      end if
   end function propagation_interval

   subroutine linop_transient_baseflow_free(self)
      class(linop_transient_baseflow_t), intent(inout) :: self

      self%time_interval = 0.0_wp
      self%fixed_step_count = 0
      self%n_timesteps = 0
      self%checkpoint_initialized = .false.
      call self%checkpoint%free()
      call linop_base_free(self)
   end subroutine linop_transient_baseflow_free

end module linop_transient_baseflow
