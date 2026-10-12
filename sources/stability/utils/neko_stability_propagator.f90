module neko_stability_propagator
   use LightKrylov, only: wp => dp
   use case, only: case_t
   use field, only: field_t
   use field_math, only: field_copy, field_rzero
   use fluid_pnpn, only: fluid_pnpn_t
   use neko_vector, only: state_vector_t
   use num_types, only: rp
   use perturb_case, only: perturb_case_t
   use perturb_fluid_pnpn, only: perturb_fluid_pnpn_t
   use time_state, only: time_state_t
   use utils, only: neko_error

   implicit none
   private

   public :: copy_state_to_primal, copy_primal_to_state
   public :: copy_state_to_perturb, copy_perturb_to_state
   public :: fit_timestep_to_interval

contains

   subroutine copy_state_to_primal(neko_case, q)
      type(case_t), intent(inout) :: neko_case
      type(state_vector_t), intent(in) :: q

      call field_copy(neko_case%fluid%u, q%u)
      call field_copy(neko_case%fluid%v, q%v)
      call field_copy(neko_case%fluid%w, q%w)
      call field_copy(neko_case%fluid%p, q%p)
      call neko_case%fluid%ulag%set(neko_case%fluid%u)
      call neko_case%fluid%vlag%set(neko_case%fluid%v)
      call neko_case%fluid%wlag%set(neko_case%fluid%w)

      select type (fluid => neko_case%fluid)
      type is (fluid_pnpn_t)
         call field_copy(fluid%u_e, fluid%u)
         call field_copy(fluid%v_e, fluid%v)
         call field_copy(fluid%w_e, fluid%w)
      end select
   end subroutine copy_state_to_primal

   subroutine copy_primal_to_state(neko_case, q)
      type(case_t), intent(inout) :: neko_case
      type(state_vector_t), intent(inout) :: q

      call field_copy(q%u, neko_case%fluid%u)
      call field_copy(q%v, neko_case%fluid%v)
      call field_copy(q%w, neko_case%fluid%w)
      call field_copy(q%p, neko_case%fluid%p)
      call enforce_dimension(q, neko_case%fluid%c_Xh%msh%gdim)
   end subroutine copy_primal_to_state

   subroutine copy_state_to_perturb(perturb_case, q)
      type(perturb_case_t), intent(inout) :: perturb_case
      type(state_vector_t), intent(in) :: q

      call field_copy(perturb_case%fluid_pert%u_pert, q%u)
      call field_copy(perturb_case%fluid_pert%v_pert, q%v)
      call field_copy(perturb_case%fluid_pert%w_pert, q%w)
      call field_copy(perturb_case%fluid_pert%p_pert, q%p)
      call perturb_case%fluid_pert%ulag%set(perturb_case%fluid_pert%u_pert)
      call perturb_case%fluid_pert%vlag%set(perturb_case%fluid_pert%v_pert)
      call perturb_case%fluid_pert%wlag%set(perturb_case%fluid_pert%w_pert)

      select type (fluid => perturb_case%fluid_pert)
      class is (perturb_fluid_pnpn_t)
         call field_copy(fluid%u_pert_e, fluid%u_pert)
         call field_copy(fluid%v_pert_e, fluid%v_pert)
         call field_copy(fluid%w_pert_e, fluid%w_pert)
      end select
   end subroutine copy_state_to_perturb

   subroutine copy_perturb_to_state(neko_case, perturb_case, q)
      type(case_t), intent(inout) :: neko_case
      type(perturb_case_t), intent(inout) :: perturb_case
      type(state_vector_t), intent(inout) :: q

      call field_copy(q%u, perturb_case%fluid_pert%u_pert)
      call field_copy(q%v, perturb_case%fluid_pert%v_pert)
      call field_copy(q%w, perturb_case%fluid_pert%w_pert)
      call field_copy(q%p, perturb_case%fluid_pert%p_pert)
      call enforce_dimension(q, neko_case%fluid%c_Xh%msh%gdim)
   end subroutine copy_perturb_to_state

   subroutine fit_timestep_to_interval(time, time_interval, fixed_nsteps)
      type(time_state_t), intent(inout) :: time
      real(wp), intent(in) :: time_interval
      integer, optional, intent(in) :: fixed_nsteps
      real(wp) :: dt, interval
      integer :: nsteps

      interval = time_interval
      dt = real(time%dt, kind=wp)

      if (interval <= 0.0_wp) call neko_error( &
           'Cannot propagate over a non-positive time interval')
      if (dt <= 0.0_wp) call neko_error( &
           'Cannot fit a non-positive timestep to a time interval')

      if (present(fixed_nsteps)) then
         if (fixed_nsteps <= 0) call neko_error( &
              'Cannot fit time interval with non-positive step count')
         nsteps = fixed_nsteps
      else
         nsteps = max(1, ceiling(interval / dt))
      end if
      time%dt = real(interval / real(nsteps, kind=wp), kind=rp)
      time%dtlag = time%dt
   end subroutine fit_timestep_to_interval

   subroutine enforce_dimension(q, gdim)
      type(state_vector_t), intent(inout) :: q
      integer, intent(in) :: gdim

      if (gdim .eq. 2) then
         call z_plane_fix(q%u)
         call z_plane_fix(q%v)
         call field_rzero(q%w)
      end if
   end subroutine enforce_dimension

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

end module neko_stability_propagator
