module periodic_orbit_system
   use LightKrylov, only: abstract_vector_rdp, abstract_jacobian_linop_rdp
   use LightKrylov, only: wp => dp
   use case, only: case_t
   use field, only: field_t
   use field_output, only: field_output_t
   use fld_file, only: fld_file_t
   use linop_transient_baseflow, only: linop_transient_baseflow_t
   use mpi_f08, only: MPI_WTIME
   use neko_ext, only: reset
   use neko_ext_vector, only: ext_state_vector_t
   use neko_stability_propagator, only: copy_primal_to_state, &
        copy_state_to_primal, fit_timestep_to_interval
   use neko_vector, only: state_vector_t
   use num_types, only: neko_rp => rp, neko_sp => sp, rp
   use perturb_case, only: perturb_case_t
   use stability_solver_tolerances, only: solver_tolerance_guard_t
   use simulation, only: simulation_finalize, simulation_init, simulation_step
   use system_base, only: stability_system_t
   use time_step_controller, only: time_step_controller_t
   use utils, only: neko_error

   implicit none

   type, extends(stability_system_t), public :: periodic_orbit_system_t
      type(state_vector_t) :: reference
      type(state_vector_t) :: phase_direction
      logical :: has_reference = .false.
      logical :: has_phase_direction = .false.
      logical :: write_outer_iterations = .false.
      logical :: suppress_outer_output = .false.
      integer :: outer_write_counter = 0
      integer :: fixed_step_count = 0
      character(len=:), allocatable :: outer_write_name
   contains
      private
      procedure, pass(self), public :: response => periodic_orbit_map
      procedure, pass(self), public :: init => periodic_orbit_system_init
      procedure, pass(self), public :: set_reference
      procedure, pass(self), public :: set_phase_direction
      procedure, pass(self), public :: enable_outer_iteration_output => &
           periodic_orbit_enable_outer_output
      procedure, pass(self), public :: disable_outer_iteration_output => &
           periodic_orbit_disable_outer_output
      procedure, pass(self), public :: set_fixed_step_count => &
           periodic_orbit_set_fixed_step_count
      procedure, pass(self), public :: free => periodic_orbit_system_free
      procedure, pass(self) :: flow_map
      procedure, pass(self) :: compute_fdot
      procedure, pass(self) :: timestep_for_period
      procedure, pass(self) :: write_outer_state
   end type periodic_orbit_system_t

   type, extends(abstract_jacobian_linop_rdp), public :: periodic_orbit_jacobian_t
      type(periodic_orbit_system_t), pointer :: system => null()
      type(linop_transient_baseflow_t) :: propagator
   contains
      private
      procedure, pass(self), public :: matvec => periodic_orbit_jacobian_matvec
      procedure, pass(self), public :: rmatvec => periodic_orbit_jacobian_rmatvec
      procedure, pass(self), public :: init => periodic_orbit_jacobian_init
      procedure, pass(self), public :: free => periodic_orbit_jacobian_free
   end type periodic_orbit_jacobian_t

contains

   subroutine periodic_orbit_system_init(self, base_case)
      class(periodic_orbit_system_t), intent(inout) :: self
      type(case_t), target, intent(inout) :: base_case

      self%base_case => base_case
      self%base_output_base_fname = trim(base_case%f_out%file_%get_base_fname())
      self%if_2d = (self%base_case%fluid%c_Xh%msh%gdim .eq. 2)
      self%write_outer_iterations = .false.
      self%suppress_outer_output = .false.
      self%outer_write_counter = 0
      self%fixed_step_count = 0
   end subroutine periodic_orbit_system_init

   subroutine set_reference(self, reference)
      class(periodic_orbit_system_t), intent(inout) :: self
      type(state_vector_t), intent(in) :: reference

      call self%reference%init_like(reference)
      call self%reference%axpby(1.0_wp, reference, 0.0_wp)
      self%has_reference = .true.

      if (.not. self%has_phase_direction) then
         call self%phase_direction%init_like(reference)
         call self%phase_direction%axpby(1.0_wp, reference, 0.0_wp)
         self%has_phase_direction = .true.
      end if
   end subroutine set_reference

   subroutine set_phase_direction(self, phase_direction)
      class(periodic_orbit_system_t), intent(inout) :: self
      type(state_vector_t), intent(in) :: phase_direction

      call self%phase_direction%init_like(phase_direction)
      call self%phase_direction%axpby(1.0_wp, phase_direction, 0.0_wp)
      self%has_phase_direction = .true.
   end subroutine set_phase_direction

   subroutine periodic_orbit_map(self, vec_in, vec_out, atol)
      class(periodic_orbit_system_t), intent(inout) :: self
      class(abstract_vector_rdp), intent(in) :: vec_in
      class(abstract_vector_rdp), intent(inout) :: vec_out
      real(wp), intent(in) :: atol
      type(state_vector_t) :: mapped

      select type (vec_in)
      type is (ext_state_vector_t)
         select type (vec_out)
         type is (ext_state_vector_t)
            call vec_out%init_like(vec_in)
            call mapped%init_like(vec_in%q)

            if (self%write_outer_iterations .and. &
                 .not. self%suppress_outer_output) then
               call self%write_outer_state(vec_in%q, vec_in%T)
            end if

            call self%flow_map(vec_in%q, vec_in%T, mapped, &
                 solver_tolerance=atol)

            call vec_out%q%axpby(1.0_wp, mapped, 0.0_wp)
            call vec_out%q%axpby(-1.0_wp, vec_in%q, 1.0_wp)
            vec_out%T = 0.0_wp

            call mapped%free()
         end select
      end select
   end subroutine periodic_orbit_map

   subroutine flow_map(self, q, period, mapped, solver_tolerance)
      class(periodic_orbit_system_t), intent(inout) :: self
      type(state_vector_t), intent(in) :: q
      real(wp), intent(in) :: period
      type(state_vector_t), intent(inout) :: mapped
      real(wp), optional, intent(in) :: solver_tolerance
      type(time_step_controller_t) :: dt_controller
      type(solver_tolerance_guard_t) :: tolerance_guard
      real(kind=rp) :: old_end_time, old_dt
      real(kind=wp) :: loop_start

      old_end_time = self%base_case%time%end_time
      old_dt = self%base_case%time%dt
      self%base_case%time%end_time = real(period, kind=rp)
      if (self%fixed_step_count > 0) then
         call fit_timestep_to_interval(self%base_case%time, period, &
              self%fixed_step_count)
      else
         call fit_timestep_to_interval(self%base_case%time, period)
      end if

      call reset(self%base_case, 0, self%base_output_base_fname)
      call copy_state_to_primal(self%base_case, q)
      if (present(solver_tolerance)) then
         call tolerance_guard%set(self%base_case, solver_tolerance)
      end if

      call dt_controller%init(self%base_case%params)
      call simulation_init(self%base_case, dt_controller)

      loop_start = MPI_WTIME()
      do while (.not. self%base_case%time%is_done())
         call simulation_step(self%base_case, dt_controller, loop_start)
      end do

      call mapped%init_like(q)
      call copy_primal_to_state(self%base_case, mapped)

      call simulation_finalize(self%base_case)
      call tolerance_guard%restore()
      self%base_case%time%end_time = old_end_time
      self%base_case%time%dt = old_dt
   end subroutine flow_map

   subroutine compute_fdot(self, q, fdot, timestep)
      class(periodic_orbit_system_t), intent(inout) :: self
      type(state_vector_t), intent(in) :: q
      type(state_vector_t), intent(inout) :: fdot
      real(wp), optional, intent(in) :: timestep
      type(state_vector_t) :: q_next
      real(wp) :: dt
      integer :: saved_fixed_step_count

      if (present(timestep)) then
         dt = timestep
      else
         dt = real(self%base_case%time%dt, kind=wp)
      end if
      if (dt <= 0.0_wp) call neko_error( &
           'Cannot compute fdot with non-positive timestep')

      call fdot%init_like(q)
      call q_next%init_like(q)
      saved_fixed_step_count = self%fixed_step_count
      self%fixed_step_count = 1
      call self%flow_map(q, dt, q_next)
      self%fixed_step_count = saved_fixed_step_count
      call fdot%axpby(1.0_wp, q_next, 0.0_wp)
      call fdot%axpby(-1.0_wp, q, 1.0_wp)
      call fdot%scal(1.0_wp / dt)
      call q_next%free()
   end subroutine compute_fdot

   subroutine periodic_orbit_set_fixed_step_count(self, period)
      class(periodic_orbit_system_t), intent(inout) :: self
      real(wp), intent(in) :: period
      real(wp) :: dt

      dt = real(self%base_case%time%dt, kind=wp)
      if (period <= 0.0_wp) call neko_error( &
           'Cannot set periodic-orbit step count from non-positive period')
      if (dt <= 0.0_wp) call neko_error( &
           'Cannot set periodic-orbit step count from non-positive timestep')

      self%fixed_step_count = max(1, ceiling(period / dt))
   end subroutine periodic_orbit_set_fixed_step_count

   real(wp) function timestep_for_period(self, period) result(dt)
      class(periodic_orbit_system_t), intent(in) :: self
      real(wp), intent(in) :: period

      if (self%fixed_step_count > 0) then
         dt = period / real(self%fixed_step_count, kind=wp)
      else
         dt = real(self%base_case%time%dt, kind=wp)
      end if
   end function timestep_for_period

   subroutine periodic_orbit_enable_outer_output(self, name)
      class(periodic_orbit_system_t), intent(inout) :: self
      character(len=*), intent(in) :: name

      if (allocated(self%outer_write_name)) deallocate(self%outer_write_name)
      self%outer_write_name = trim(name)
      self%outer_write_counter = 0
      self%write_outer_iterations = .true.
   end subroutine periodic_orbit_enable_outer_output

   subroutine periodic_orbit_disable_outer_output(self)
      class(periodic_orbit_system_t), intent(inout) :: self

      self%write_outer_iterations = .false.
      self%suppress_outer_output = .false.
      self%outer_write_counter = 0
      if (allocated(self%outer_write_name)) deallocate(self%outer_write_name)
   end subroutine periodic_orbit_disable_outer_output

   subroutine write_outer_state(self, q, period)
      class(periodic_orbit_system_t), intent(inout) :: self
      type(state_vector_t), intent(in) :: q
      real(wp), intent(in) :: period
      type(field_output_t) :: output
      integer :: idx
      integer :: unit
      character(len=:), allocatable :: period_file

      if (.not. allocated(self%outer_write_name)) return

      call copy_state_to_primal(self%base_case, q)

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

      period_file = 'viz/' // trim(self%outer_write_name) // '_period.csv'
      if (idx .eq. 0) then
         open(newunit=unit, file=period_file, status='replace', action='write')
         write(unit, '(A)') 'iteration,T'
      else
         open(newunit=unit, file=period_file, status='old', position='append', &
              action='write')
      end if
      write(unit, '(I0,A,ES24.16E3)') idx, ',', period
      close(unit)

      self%outer_write_counter = self%outer_write_counter + 1
   end subroutine write_outer_state

   subroutine periodic_orbit_system_free(self)
      class(periodic_orbit_system_t), intent(inout) :: self

      if (self%has_reference) call self%reference%free()
      if (self%has_phase_direction) call self%phase_direction%free()
      if (allocated(self%jacobian)) then
         select type (jac => self%jacobian)
         type is (periodic_orbit_jacobian_t)
            call jac%free()
         end select
         deallocate(self%jacobian)
      end if
      nullify(self%base_case)
      self%has_reference = .false.
      self%has_phase_direction = .false.
      self%if_2d = .false.
      self%fixed_step_count = 0
      if (allocated(self%base_output_base_fname)) then
         deallocate(self%base_output_base_fname)
      end if
      call self%disable_outer_iteration_output()
   end subroutine periodic_orbit_system_free

   subroutine periodic_orbit_jacobian_init(self, system, lns_case, adjoint_case)
      class(periodic_orbit_jacobian_t), intent(inout) :: self
      type(periodic_orbit_system_t), target, intent(inout) :: system
      type(perturb_case_t), target, intent(inout) :: lns_case
      type(perturb_case_t), target, intent(inout) :: adjoint_case

      self%system => system
      call self%propagator%init(system%base_case, lns_case, adjoint_case)
   end subroutine periodic_orbit_jacobian_init

   subroutine periodic_orbit_jacobian_matvec(self, vec_in, vec_out)
      class(periodic_orbit_jacobian_t), intent(inout) :: self
      class(abstract_vector_rdp), intent(in) :: vec_in
      class(abstract_vector_rdp), intent(inout) :: vec_out
      type(state_vector_t) :: mapped, tangent, fdot_endpoint, fdot_initial

      select type (x => self%X)
      type is (ext_state_vector_t)
         select type (v => vec_in)
         type is (ext_state_vector_t)
            select type (y => vec_out)
            type is (ext_state_vector_t)
               call y%init_like(x)

               call mapped%init_like(x%q)
               call tangent%init_like(x%q)
               call fdot_endpoint%init_like(x%q)
               call fdot_initial%init_like(x%q)

               call self%system%flow_map(x%q, x%T, mapped)
               if (self%system%fixed_step_count > 0) then
                  call self%propagator%set_baseflow(x%q, &
                       time_interval=x%T, &
                       fixed_step_count=self%system%fixed_step_count)
               else
                  call self%propagator%set_baseflow(x%q, time_interval=x%T)
               end if
               call self%propagator%matvec(v%q, tangent)

               call y%q%axpby(1.0_wp, tangent, 0.0_wp)
               call y%q%axpby(-1.0_wp, v%q, 1.0_wp)

               call self%system%compute_fdot(mapped, fdot_endpoint, &
                    self%system%timestep_for_period(x%T))
               call y%q%axpby(v%T, fdot_endpoint, 1.0_wp)

               call self%system%compute_fdot(x%q, fdot_initial, &
                    self%system%timestep_for_period(x%T))
               y%T = v%q%dot(fdot_initial)

               call mapped%free()
               call tangent%free()
               call fdot_endpoint%free()
               call fdot_initial%free()
            end select
         end select
      end select
   end subroutine periodic_orbit_jacobian_matvec

   subroutine periodic_orbit_jacobian_rmatvec(self, vec_in, vec_out)
      class(periodic_orbit_jacobian_t), intent(inout) :: self
      class(abstract_vector_rdp), intent(in) :: vec_in
      class(abstract_vector_rdp), intent(inout) :: vec_out
      type(state_vector_t) :: mapped, adjoint_image, fdot_endpoint, fdot_initial

      select type (x => self%X)
      type is (ext_state_vector_t)
         select type (v => vec_in)
         type is (ext_state_vector_t)
            select type (y => vec_out)
            type is (ext_state_vector_t)
               call y%init_like(x)
               call mapped%init_like(x%q)
               call adjoint_image%init_like(x%q)
               call fdot_endpoint%init_like(x%q)
               call fdot_initial%init_like(x%q)

               call self%system%flow_map(x%q, x%T, mapped)
               if (self%system%fixed_step_count > 0) then
                  call self%propagator%set_baseflow(x%q, &
                       time_interval=x%T, &
                       fixed_step_count=self%system%fixed_step_count)
               else
                  call self%propagator%set_baseflow(x%q, time_interval=x%T)
               end if
               call self%propagator%rmatvec(v%q, adjoint_image)

               call y%q%axpby(1.0_wp, adjoint_image, 0.0_wp)
               call y%q%axpby(-1.0_wp, v%q, 1.0_wp)

               call self%system%compute_fdot(x%q, fdot_initial, &
                    self%system%timestep_for_period(x%T))
               call y%q%axpby(v%T, fdot_initial, 1.0_wp)

               call self%system%compute_fdot(mapped, fdot_endpoint, &
                    self%system%timestep_for_period(x%T))
               y%T = v%q%dot(fdot_endpoint)

               call mapped%free()
               call adjoint_image%free()
               call fdot_endpoint%free()
               call fdot_initial%free()
            end select
         end select
      end select
   end subroutine periodic_orbit_jacobian_rmatvec

   subroutine periodic_orbit_jacobian_free(self)
      class(periodic_orbit_jacobian_t), intent(inout) :: self

      nullify(self%system)
      call self%propagator%free()
   end subroutine periodic_orbit_jacobian_free

   subroutine z_plane_fix(fld)
      type(field_t), intent(inout) :: fld
      integer :: iel, iz, iy, ix

      do iel = 1, fld%msh%nelv
         do iz = 2, fld%xh%lz
            do iy = 1, fld%xh%ly
               do ix = 1, fld%xh%lx
                  fld%x(ix, iy, iz, iel) = fld%x(ix, iy, 1, iel)
               end do
            end do
         end do
      end do
   end subroutine z_plane_fix

end module periodic_orbit_system
