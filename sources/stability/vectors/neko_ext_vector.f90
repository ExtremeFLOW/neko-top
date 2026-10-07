module neko_ext_vector
   use LightKrylov, only: abstract_vector_rdp
   use LightKrylov, only: wp => dp
   use neko_vector, only: state_vector_t

   implicit none

   type, extends(abstract_vector_rdp), public :: ext_state_vector_t
      type(state_vector_t) :: q
      real(wp) :: T = 0.0_wp
   contains
      private
      procedure, pass(self), public :: zero
      procedure, pass(self), public :: dot
      procedure, pass(self), public :: scal
      procedure, pass(self), public :: axpby
      procedure, pass(self), public :: rand
      procedure, pass(self), public :: get_size
      procedure, pass(self), public :: init_like => ext_state_vector_init_like
      procedure, pass(self), public :: free => ext_state_vector_free
      procedure, pass(self), public :: write => ext_state_vector_write
   end type ext_state_vector_t

contains

   subroutine ext_state_vector_init_like(self, mold)
      class(ext_state_vector_t), intent(inout) :: self
      class(abstract_vector_rdp), intent(in) :: mold

      select type (mold)
      type is (ext_state_vector_t)
         call self%q%init_like(mold%q)
      type is (state_vector_t)
         call self%q%init_like(mold)
      end select

      self%is_initialized = self%q%is_initialized
      self%owns_data = .true.
   end subroutine ext_state_vector_init_like

   subroutine zero(self)
      class(ext_state_vector_t), intent(inout) :: self

      call self%q%zero()
      self%T = 0.0_wp
   end subroutine zero

   real(wp) function dot(self, vec) result(alpha)
      class(ext_state_vector_t), intent(in) :: self
      class(abstract_vector_rdp), intent(in) :: vec

      alpha = 0.0_wp
      select type (vec)
      type is (ext_state_vector_t)
         alpha = self%q%dot(vec%q) + self%T * vec%T
      end select
   end function dot

   subroutine scal(self, alpha)
      class(ext_state_vector_t), intent(inout) :: self
      real(wp), intent(in) :: alpha

      call self%q%scal(alpha)
      self%T = alpha * self%T
   end subroutine scal

   subroutine axpby(alpha, vec, beta, self)
      class(ext_state_vector_t), intent(inout) :: self
      class(abstract_vector_rdp), intent(in) :: vec
      real(wp), intent(in) :: alpha, beta

      select type (vec)
      type is (ext_state_vector_t)
         if (.not. self%is_initialized) call self%init_like(vec)
         call self%q%axpby(alpha, vec%q, beta)
         self%T = alpha * vec%T + beta * self%T
      end select
   end subroutine axpby

   subroutine rand(self, ifnorm)
      class(ext_state_vector_t), intent(inout) :: self
      logical, optional, intent(in) :: ifnorm
      real(wp) :: alpha

      call self%q%rand(.false.)
      self%T = 0.0_wp

      if (present(ifnorm)) then
         if (ifnorm) then
            alpha = self%norm()
            if (alpha > 0.0_wp) call self%scal(1.0_wp / alpha)
         end if
      else
         alpha = self%norm()
         if (alpha > 0.0_wp) call self%scal(1.0_wp / alpha)
      end if
   end subroutine rand

   integer function get_size(self) result(N)
      class(ext_state_vector_t), intent(in) :: self

      N = self%q%get_size() + 1
   end function get_size

   subroutine ext_state_vector_free(self)
      class(ext_state_vector_t), intent(inout) :: self

      call self%q%free()
      self%T = 0.0_wp
      self%is_initialized = .false.
      self%owns_data = .true.
   end subroutine ext_state_vector_free

   subroutine ext_state_vector_write(self, idx)
      class(ext_state_vector_t), intent(inout) :: self
      integer, intent(in) :: idx

      call self%q%write(idx)
   end subroutine ext_state_vector_write

end module neko_ext_vector
