module neko_complex_vector
   use LightKrylov, only: abstract_vector, abstract_vector_cdp
   use LightKrylov, only: wp => dp
   use neko_vector, only: state_vector_t

   implicit none

   type, extends(abstract_vector_cdp), public :: complex_state_vector_t
      type(state_vector_t) :: re
      type(state_vector_t) :: im
   contains
      private
      procedure, pass(self), public :: zero
      procedure, pass(self), public :: dot
      procedure, pass(self), public :: scal
      procedure, pass(self), public :: axpby
      procedure, pass(self), public :: rand
      procedure, pass(self), public :: get_size
      procedure, pass(self), public :: init_like => complex_state_vector_init_like
      procedure, pass(self), public :: assign => complex_state_vector_assign
      procedure, pass(self), public :: swap => complex_state_vector_swap
      procedure, pass(self), public :: free => complex_state_vector_free
      procedure, pass(self), public :: write => complex_state_vector_write
   end type complex_state_vector_t

contains

   subroutine complex_state_vector_init_like(self, mold)
      class(complex_state_vector_t), intent(inout) :: self
      class(abstract_vector_cdp), intent(in) :: mold

      select type (mold)
      type is (complex_state_vector_t)
         call self%re%init_like(mold%re)
         call self%im%init_like(mold%im)
      end select

      self%is_initialized = self%re%is_initialized .and. self%im%is_initialized
      self%owns_data = .true.
   end subroutine complex_state_vector_init_like

   subroutine complex_state_vector_assign(self, from)
      class(complex_state_vector_t), intent(inout) :: self
      class(abstract_vector), intent(in) :: from

      select type (from)
      type is (complex_state_vector_t)
         call self%re%assign(from%re)
         call self%im%assign(from%im)
         self%is_initialized = self%re%is_initialized .and. &
              self%im%is_initialized
         self%owns_data = .true.
      class default
         error stop 'complex_state_vector_assign: from must be a complex_state_vector_t'
      end select
   end subroutine complex_state_vector_assign

   subroutine complex_state_vector_swap(self, other)
      class(complex_state_vector_t), intent(inout) :: self
      class(abstract_vector), intent(inout) :: other
      logical :: logical_tmp

      select type (other)
      type is (complex_state_vector_t)
         call self%re%swap(other%re)
         call self%im%swap(other%im)

         logical_tmp = self%is_initialized
         self%is_initialized = other%is_initialized
         other%is_initialized = logical_tmp

         logical_tmp = self%owns_data
         self%owns_data = other%owns_data
         other%owns_data = logical_tmp
      class default
         error stop 'complex_state_vector_swap: other must be a complex_state_vector_t'
      end select
   end subroutine complex_state_vector_swap

   subroutine zero(self)
      class(complex_state_vector_t), intent(inout) :: self

      call self%re%zero()
      call self%im%zero()
   end subroutine zero

   complex(wp) function dot(self, vec) result(alpha)
      class(complex_state_vector_t), intent(in) :: self
      class(abstract_vector_cdp), intent(in) :: vec
      real(wp) :: rr, ri, ir, ii

      alpha = cmplx(0.0_wp, 0.0_wp, kind=wp)
      select type (vec)
      type is (complex_state_vector_t)
         rr = self%re%dot(vec%re)
         ri = self%re%dot(vec%im)
         ir = self%im%dot(vec%re)
         ii = self%im%dot(vec%im)
         alpha = cmplx(rr + ii, ri - ir, kind=wp)
      end select
   end function dot

   subroutine scal(self, alpha)
      class(complex_state_vector_t), intent(inout) :: self
      complex(wp), intent(in) :: alpha
      type(state_vector_t) :: re_old, im_old

      call re_old%init_like(self%re)
      call im_old%init_like(self%im)
      call re_old%axpby(1.0_wp, self%re, 0.0_wp)
      call im_old%axpby(1.0_wp, self%im, 0.0_wp)

      call self%re%zero()
      call self%re%axpby(real(alpha, kind=wp), re_old, 0.0_wp)
      call self%re%axpby(-aimag(alpha), im_old, 1.0_wp)

      call self%im%zero()
      call self%im%axpby(aimag(alpha), re_old, 0.0_wp)
      call self%im%axpby(real(alpha, kind=wp), im_old, 1.0_wp)

      call re_old%free()
      call im_old%free()
   end subroutine scal

   subroutine axpby(alpha, vec, beta, self)
      class(complex_state_vector_t), intent(inout) :: self
      class(abstract_vector_cdp), intent(in) :: vec
      complex(wp), intent(in) :: alpha, beta
      type(state_vector_t) :: y_re, y_im

      select type (vec)
      type is (complex_state_vector_t)
         call self%init_like(vec)
         call y_re%init_like(self%re)
         call y_im%init_like(self%im)
         call y_re%axpby(1.0_wp, self%re, 0.0_wp)
         call y_im%axpby(1.0_wp, self%im, 0.0_wp)

         call self%re%zero()
         call self%re%axpby(real(alpha, kind=wp), vec%re, 0.0_wp)
         call self%re%axpby(-aimag(alpha), vec%im, 1.0_wp)
         call self%re%axpby(real(beta, kind=wp), y_re, 1.0_wp)
         call self%re%axpby(-aimag(beta), y_im, 1.0_wp)

         call self%im%zero()
         call self%im%axpby(aimag(alpha), vec%re, 0.0_wp)
         call self%im%axpby(real(alpha, kind=wp), vec%im, 1.0_wp)
         call self%im%axpby(aimag(beta), y_re, 1.0_wp)
         call self%im%axpby(real(beta, kind=wp), y_im, 1.0_wp)

         call y_re%free()
         call y_im%free()
      end select
   end subroutine axpby

   subroutine rand(self, ifnorm)
      class(complex_state_vector_t), intent(inout) :: self
      logical, optional, intent(in) :: ifnorm
      real(wp) :: alpha

      call self%re%rand(.false.)
      call self%im%rand(.false.)
      if (present(ifnorm)) then
         if (ifnorm) then
            alpha = self%norm()
            if (alpha > 0.0_wp) call self%scal(cmplx(1.0_wp / alpha, 0.0_wp, kind=wp))
         end if
      end if
   end subroutine rand

   integer function get_size(self) result(N)
      class(complex_state_vector_t), intent(in) :: self

      N = self%re%get_size()
   end function get_size

   subroutine complex_state_vector_free(self)
      class(complex_state_vector_t), intent(inout) :: self

      call self%re%free()
      call self%im%free()
      self%is_initialized = .false.
      self%owns_data = .true.
   end subroutine complex_state_vector_free

   subroutine complex_state_vector_write(self, idx)
      class(complex_state_vector_t), intent(inout) :: self
      integer, intent(in) :: idx

      call self%re%write(2 * idx)
      call self%im%write(2 * idx + 1)
   end subroutine complex_state_vector_write

end module neko_complex_vector
