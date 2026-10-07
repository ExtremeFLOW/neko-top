module stability_solver_tolerances
   use case, only: case_t
   use fluid_pnpn, only: fluid_pnpn_t
   use LightKrylov, only: wp => dp
   use logger, only: neko_log
   use num_types, only: rp
   use utils, only: neko_error

   implicit none
   private

   type, public :: solver_tolerance_guard_t
      class(case_t), pointer :: base_case => null()
      real(kind=rp) :: old_velocity_tolerance = 0.0_rp
      real(kind=rp) :: old_pressure_tolerance = 0.0_rp
      logical :: active = .false.
   contains
      procedure, pass(self), public :: set => solver_tolerance_guard_set
      procedure, pass(self), public :: restore => solver_tolerance_guard_restore
   end type solver_tolerance_guard_t

contains

   subroutine solver_tolerance_guard_set(self, base_case, tolerance)
      class(solver_tolerance_guard_t), intent(inout) :: self
      type(case_t), target, intent(inout) :: base_case
      real(kind=wp), intent(in) :: tolerance
      real(kind=rp) :: tolerance_rp
      character(len=128) :: log_buf

      if (tolerance <= 0.0_wp) call neko_error( &
           'Cannot set dynamic solver tolerance from non-positive value')

      if (self%active) call self%restore()

      tolerance_rp = real(tolerance, kind=rp)
      select type (fluid => base_case%fluid)
      type is (fluid_pnpn_t)
         if (.not. allocated(fluid%ksp_vel)) call neko_error( &
              'Velocity solver is not initialized for dynamic tolerance')
         if (.not. allocated(fluid%ksp_prs)) call neko_error( &
              'Pressure solver is not initialized for dynamic tolerance')

         self%base_case => base_case
         self%old_velocity_tolerance = fluid%ksp_vel%abs_tol
         self%old_pressure_tolerance = fluid%ksp_prs%abs_tol
         fluid%ksp_vel%abs_tol = tolerance_rp
         fluid%ksp_prs%abs_tol = tolerance_rp
         self%active = .true.

         write(log_buf, '(A,ES13.6)') &
              'Dynamic velocity/pressure abs tol:', tolerance_rp
         call neko_log%message(log_buf)
      class default
         call neko_error('Dynamic solver tolerance is only implemented for pnpn')
      end select
   end subroutine solver_tolerance_guard_set

   subroutine solver_tolerance_guard_restore(self)
      class(solver_tolerance_guard_t), intent(inout) :: self

      if (.not. self%active) return
      if (.not. associated(self%base_case)) then
         self%active = .false.
         return
      end if

      select type (fluid => self%base_case%fluid)
      type is (fluid_pnpn_t)
         if (allocated(fluid%ksp_vel)) then
            fluid%ksp_vel%abs_tol = self%old_velocity_tolerance
         end if
         if (allocated(fluid%ksp_prs)) then
            fluid%ksp_prs%abs_tol = self%old_pressure_tolerance
         end if
      end select

      nullify(self%base_case)
      self%old_velocity_tolerance = 0.0_rp
      self%old_pressure_tolerance = 0.0_rp
      self%active = .false.
   end subroutine solver_tolerance_guard_restore

end module stability_solver_tolerances
