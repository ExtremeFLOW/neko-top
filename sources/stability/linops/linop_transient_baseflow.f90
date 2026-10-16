module linop_transient_baseflow
   use LightKrylov, only: abstract_vector_rdp
   use LightKrylov, only: wp => dp
   use case, only: case_t
   use json_file_module, only: json_file
   use json_utils, only: json_get
   use linop_base, only: linop_base_t, linop_base_init, &
        linop_base_require_adjoint, linop_base_require_lns, &
        linop_base_store_baseflow, linop_base_free
   use mpi_f08, only: MPI_WTIME
   use neko_ext, only: reset, reset_adjoint
   use neko_stability_propagator, only: copy_perturb_to_state, &
        copy_state_to_perturb, copy_state_to_primal, fit_timestep_to_interval
   use neko_vector, only: state_vector_t
   use num_types, only: rp
   use perturb_case, only: perturb_case_t
   use simulation, only: simulation_init, simulation_step
   use simulation_checkpoint, only: simulation_checkpoint_t
   use simulation_perturb, only: simulation_adjoint_step, &
        simulation_lns_step, simulation_perturb_init
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
        adjoint_case, baseflow, time_interval, fixed_step_count)
      class(linop_transient_baseflow_t), intent(inout) :: self
      type(case_t), target, intent(inout) :: base_case
      type(perturb_case_t), optional, target, intent(inout) :: lns_case
      type(perturb_case_t), optional, target, intent(inout) :: adjoint_case
      type(state_vector_t), optional, intent(in) :: baseflow
      real(wp), optional, intent(in) :: time_interval
      integer, optional, intent(in) :: fixed_step_count
      type(json_file) :: checkpoint_params
      type(time_step_controller_t) :: base_dt_controller
      type(time_step_controller_t) :: perturb_dt_controller

      call linop_base_init(self, base_case, lns_case, adjoint_case, baseflow)
      self%time_interval = 0.0_wp
      self%fixed_step_count = 0
      self%n_timesteps = 0
      self%checkpoint_initialized = .false.

      if ("checkpoints" .in. base_case%params) then
         call json_get(base_case%params, 'checkpoints', checkpoint_params)
         call self%checkpoint%init(base_case, checkpoint_params)
         self%checkpoint_initialized = .true.
      else if (self%adjoint_initialized) then
         call neko_error('Transient-baseflow adjoint requires checkpoints')
      else
         call neko_warning('Transient-baseflow linop initialized without ' // &
              'checkpoints; rmatvec will be unavailable')
      end if

      call base_dt_controller%init(self%base_case%params)
      call simulation_init(self%base_case, base_dt_controller)

      if (self%lns_initialized) then
         call perturb_dt_controller%init(self%base_case%params)
         call simulation_perturb_init(self%lns_case, perturb_dt_controller)
      end if

      if (self%adjoint_initialized) then
         call perturb_dt_controller%init(self%base_case%params)
         call simulation_perturb_init(self%adjoint_case, perturb_dt_controller)
      end if

      if (present(baseflow) .or. present(time_interval) .or. &
           present(fixed_step_count)) then
         if (.not. present(baseflow)) then
            call neko_error('Transient-baseflow set_baseflow requires baseflow')
         end if
         call self%set_baseflow(baseflow, time_interval, fixed_step_count)
      end if
   end subroutine linop_transient_baseflow_init

   subroutine linop_transient_baseflow_set_baseflow(self, baseflow, &
        time_interval, fixed_step_count)
      class(linop_transient_baseflow_t), intent(inout) :: self
      type(state_vector_t), intent(in) :: baseflow
      real(wp), optional, intent(in) :: time_interval
      integer, optional, intent(in) :: fixed_step_count
      type(time_step_controller_t) :: dt_controller
      real(wp) :: dt, loop_start
      integer :: istep

      call linop_base_store_baseflow(self, baseflow)

      if (present(time_interval)) then
         self%time_interval = time_interval
      else
         call neko_error('Transient-baseflow set_baseflow requires time_interval')
      end if
      if (self%time_interval <= 0.0_wp) then
         call neko_error('Transient-baseflow time interval must be positive')
      end if

      if (present(fixed_step_count)) then
         if (fixed_step_count <= 0) then
            call neko_error('Transient-baseflow step count must be positive')
         end if
         self%fixed_step_count = fixed_step_count
         self%n_timesteps = fixed_step_count
      else
         dt = real(self%base_case%time%dt, kind=wp)
         if (dt <= 0.0_wp) then
            call neko_error('Transient-baseflow timestep must be positive')
         end if
         self%n_timesteps = max(1, ceiling(self%time_interval / dt))
         self%fixed_step_count = self%n_timesteps
      end if

      self%base_case%time%end_time = real(self%time_interval, kind=rp)
      call fit_timestep_to_interval(self%base_case%time, self%time_interval, &
           self%n_timesteps)

      if (self%lns_initialized) then
         self%lns_case%time%end_time = real(self%time_interval, kind=rp)
         self%lns_case%time%dt = self%base_case%time%dt
         self%lns_case%time%dtlag = self%base_case%time%dtlag
      end if
      if (self%adjoint_initialized) then
         self%adjoint_case%time%end_time = real(self%time_interval, kind=rp)
         self%adjoint_case%time%dt = self%base_case%time%dt
         self%adjoint_case%time%dtlag = self%base_case%time%dtlag

         call reset(self%base_case, 0, self%base_output_base_fname)
         call self%checkpoint%reset()
         call copy_state_to_primal(self%base_case, self%baseflow)

         call dt_controller%init(self%base_case%params)
         loop_start = MPI_WTIME()

         do istep = 1, self%n_timesteps
            call simulation_step(self%base_case, dt_controller, loop_start)
            call self%checkpoint%save(self%base_case)
         end do

         call reset(self%base_case, 0, self%base_output_base_fname)
         call copy_state_to_primal(self%base_case, self%baseflow)
      end if
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
      type(time_step_controller_t) :: primal_dt_controller
      type(time_step_controller_t) :: perturb_dt_controller
      real(kind=rp) :: cfl
      real(kind=wp) :: loop_start
      integer :: istep

      call reset(self%base_case, 0, self%base_output_base_fname)
      call reset_adjoint(self%lns_case, self%base_case, 0, &
           self%lns_output_base_fname)
      call copy_state_to_primal(self%base_case, self%baseflow)
      call copy_state_to_perturb(self%lns_case, vec_in)

      call primal_dt_controller%init(self%base_case%params)
      call perturb_dt_controller%init(self%base_case%params)

      cfl = self%lns_case%fluid_pert%compute_cfl(self%lns_case%time%dt)
      loop_start = MPI_WTIME()
      do istep = 1, self%n_timesteps
         call simulation_step(self%base_case, primal_dt_controller, loop_start)
         call simulation_lns_step(self%lns_case, perturb_dt_controller, &
              cfl, loop_start)
      end do

      call vec_out%init_like(vec_in)
      call copy_perturb_to_state(self%base_case, self%lns_case, vec_out)
   end subroutine apply_lns

   subroutine apply_adjoint(self, vec_in, vec_out)
      class(linop_transient_baseflow_t), intent(inout) :: self
      type(state_vector_t), intent(in) :: vec_in
      type(state_vector_t), intent(inout) :: vec_out
      type(time_step_controller_t) :: dt_controller
      real(kind=wp) :: loop_start
      real(kind=rp) :: cfl, total_time
      integer :: i

      call reset_adjoint(self%adjoint_case, self%base_case, 0, &
           self%adjoint_output_base_fname)
      call copy_state_to_perturb(self%adjoint_case, vec_in)

      call dt_controller%init(self%base_case%params)

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
   end subroutine apply_adjoint

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
