module neko_resolvent
   use LightKrylov, only: abstract_linop_cdp, abstract_vector_cdp
   use LightKrylov, only: wp => dp
   use case, only: case_t
   use neko_complex_vector, only: complex_state_vector_t
   use linop_fixed_baseflow, only: linop_fixed_baseflow_t
   use perturb_case, only: perturb_case_t

   implicit none

   type, extends(abstract_linop_cdp), public :: resolvent_shifted_operator_t
      type(linop_fixed_baseflow_t) :: propagator
      real(wp) :: omega = 0.0_wp
      real(wp) :: tau = 1.0_wp
   contains
      private
      procedure, pass(self), public :: matvec => resolvent_shifted_matvec
      procedure, pass(self), public :: rmatvec => resolvent_shifted_matvec
      procedure, pass(self), public :: init => resolvent_shifted_operator_init
      procedure, pass(self), public :: free => resolvent_shifted_operator_free
   end type resolvent_shifted_operator_t

contains

   subroutine resolvent_shifted_operator_init(self, base_case, lns_case, omega)
      class(resolvent_shifted_operator_t), intent(inout) :: self
      type(case_t), target, intent(inout) :: base_case
      type(perturb_case_t), target, intent(inout) :: lns_case
      real(wp), intent(in) :: omega

      call self%propagator%init(base_case, lns_case)
      self%omega = omega
      self%tau = real(base_case%time%end_time, kind=wp)
      if (self%tau <= 0.0_wp) self%tau = 1.0_wp
   end subroutine resolvent_shifted_operator_init

   subroutine resolvent_shifted_matvec(self, vec_in, vec_out)
      class(resolvent_shifted_operator_t), intent(inout) :: self
      class(abstract_vector_cdp), intent(in) :: vec_in
      class(abstract_vector_cdp), intent(inout) :: vec_out
      type(complex_state_vector_t) :: phi

      select type (x => vec_in)
      type is (complex_state_vector_t)
         select type (y => vec_out)
         type is (complex_state_vector_t)
            call y%init_like(x)
            call phi%init_like(x)

            call self%propagator%matvec(x%re, phi%re)
            call self%propagator%matvec(x%im, phi%im)

            ! y = i omega x - (Phi_tau x - x) / tau
            call y%re%axpby(1.0_wp / self%tau, x%re, 0.0_wp)
            call y%re%axpby(-1.0_wp / self%tau, phi%re, 1.0_wp)
            call y%re%axpby(-self%omega, x%im, 1.0_wp)

            call y%im%axpby(1.0_wp / self%tau, x%im, 0.0_wp)
            call y%im%axpby(-1.0_wp / self%tau, phi%im, 1.0_wp)
            call y%im%axpby(self%omega, x%re, 1.0_wp)

            call phi%free()
         end select
      end select
   end subroutine resolvent_shifted_matvec

   subroutine resolvent_shifted_operator_free(self)
      class(resolvent_shifted_operator_t), intent(inout) :: self

      call self%propagator%free()
      self%omega = 0.0_wp
      self%tau = 1.0_wp
   end subroutine resolvent_shifted_operator_free

end module neko_resolvent
