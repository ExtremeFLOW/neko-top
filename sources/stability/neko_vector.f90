module neko_vector
   use LightKrylov, only: abstract_vector_rdp
   use LightKrylov, only: wp => dp
   use stdlib_optval, only: optval
   use num_types, only: rp, sp
   use field, only: field_t
   use coefs, only: coef_t
   use fld_file_output, only: fld_file_output_t
   use neko_config, only: NEKO_BCKND_DEVICE
   use field_math, only: field_rzero, field_cmult, field_add3s2
   use math, only: glsc3
   use device_math, only: device_glsc3
   use device, only: device_associated, device_memcpy, HOST_TO_DEVICE
   use gather_scatter, only: GS_OP_ADD
   use user_access_singleton, only: neko_user_access
   use, intrinsic :: iso_c_binding, only: c_null_ptr

   implicit none

   type, extends(abstract_vector_rdp), public :: state_vector_t
      ! velocity and pressure fields
      type(field_t) :: u, v, w, p
      ! we need the mass matrix to integrate fields..
      type(coef_t), pointer :: coef => null()
      logical :: if_2d = .false.
    contains
      private
      procedure, pass(self), public :: zero
      procedure, pass(self), public :: dot
      procedure, pass(self), public :: scal
      procedure, pass(self), public :: axpby
      procedure, pass(self), public :: rand
      procedure, pass(self), public :: get_size
      ! we also want some other things
      procedure, pass(self), public :: init => state_vector_init
      procedure, pass(self), public :: free => state_vector_free
      procedure, pass(self), public :: write => state_vector_write
   end type state_vector_t
 
 contains
 
   subroutine state_vector_init(self)
      class(state_vector_t), intent(inout) :: self

      call state_vector_attach_coef(self)
      if (state_vector_is_initialized(self)) return

      call state_vector_sanitize_field(self%u)
      call state_vector_sanitize_field(self%v)
      call state_vector_sanitize_field(self%w)
      call state_vector_sanitize_field(self%p)

      call self%u%init(self%coef%dof, fld_name='state_u')
      call self%v%init(self%coef%dof, fld_name='state_v')
      call self%w%init(self%coef%dof, fld_name='state_w')
      call self%p%init(self%coef%dof, fld_name='state_p')
   end subroutine state_vector_init
 
   subroutine zero(self)
      class(state_vector_t), intent(inout) :: self

      call field_rzero(self%u)
      call field_rzero(self%v)
      call field_rzero(self%w)
      call field_rzero(self%p)
   end subroutine zero
 
   real(kind=wp) function dot(self, vec) result(alpha)
      class(state_vector_t), intent(in) :: self
      class(abstract_vector_rdp), intent(in) :: vec
      integer :: n
      real(kind=rp) :: alpha_rp

      select type (vec)
      type is (state_vector_t)
         n = self%u%size()
         alpha_rp = 0.0_rp

         if (NEKO_BCKND_DEVICE .eq. 1) then
            alpha_rp = device_glsc3(self%u%x_d, vec%u%x_d, self%coef%B_d, n)
            alpha_rp = alpha_rp + &
                 device_glsc3(self%v%x_d, vec%v%x_d, self%coef%B_d, n)
            if (.not. self%if_2d) then
               alpha_rp = alpha_rp + &
                    device_glsc3(self%w%x_d, vec%w%x_d, self%coef%B_d, n)
            end if
         else
            alpha_rp = glsc3(self%u%x, vec%u%x, self%coef%B, n)
            alpha_rp = alpha_rp + glsc3(self%v%x, vec%v%x, self%coef%B, n)
            if (.not. self%if_2d) then
               alpha_rp = alpha_rp + glsc3(self%w%x, vec%w%x, self%coef%B, n)
            end if
         end if

         alpha = real(alpha_rp, wp)
      end select
   end function dot
 
   subroutine scal(self, alpha)
      class(state_vector_t), intent(inout) :: self
      real(kind=wp), intent(in) :: alpha
      real(kind=rp) :: alpha_rp

      alpha_rp = real(alpha, rp)
      call field_cmult(self%u, alpha_rp)
      call field_cmult(self%v, alpha_rp)
      call field_cmult(self%w, alpha_rp)
      if (self%if_2d) call field_rzero(self%w)
      call field_cmult(self%p, alpha_rp)
   end subroutine scal
 
   subroutine axpby(alpha, vec, beta, self)
      class(state_vector_t), intent(inout) :: self
      class(abstract_vector_rdp), intent(in) :: vec
      real(kind=wp), intent(in) :: alpha, beta
      real(kind=rp) :: alpha_rp, beta_rp

      select type (vec)
      type is (state_vector_t)
         if (.not. associated(self%coef)) self%coef => vec%coef
         call self%init()

         alpha_rp = real(alpha, rp)
         beta_rp = real(beta, rp)

         call field_add3s2(self%u, self%u, vec%u, beta_rp, alpha_rp)
         call field_add3s2(self%v, self%v, vec%v, beta_rp, alpha_rp)
         call field_add3s2(self%w, self%w, vec%w, beta_rp, alpha_rp)
         call field_add3s2(self%p, self%p, vec%p, beta_rp, alpha_rp)
      end select
   end subroutine axpby
 
   integer function get_size(self) result(N)
      class(state_vector_t), intent(in) :: self

      N = self%u%size() + self%v%size() + self%w%size()
   end function get_size
 
   subroutine rand(self, ifnorm)
      class(state_vector_t), intent(inout) :: self
      logical, optional, intent(in) :: ifnorm
      logical :: normalize
      real(kind=wp) :: alpha

      normalize = optval(ifnorm, .true.)
      call rand_ic(self%u, self%v, self%w, self%if_2d)

      call self%coef%gs_h%op(self%u, GS_OP_ADD)
      call self%coef%gs_h%op(self%v, GS_OP_ADD)
      call self%coef%gs_h%op(self%w, GS_OP_ADD)

      if (normalize) then
         alpha = self%norm()
         call self%scal(1.0_wp / alpha)
      end if
   end subroutine rand

   subroutine state_vector_free(self)
      class(state_vector_t), intent(inout) :: self

      call state_vector_sanitize_field(self%u)
      call state_vector_sanitize_field(self%v)
      call state_vector_sanitize_field(self%w)
      call state_vector_sanitize_field(self%p)

      call self%u%free()
      call self%v%free()
      call self%w%free()
      call self%p%free()

      nullify(self%coef)
      self%if_2d = .false.
   end subroutine state_vector_free

  ! User defined initial condition
  subroutine rand_ic(u, v, w, if_2d)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    logical, intent(in) :: if_2d
    integer :: iel, ix, iy, iz
    real(kind=rp) :: fcoeff(3), xl(2)



    do iel = 1, u%msh%nelv
       do iz = 1, u%Xh%lz
          do iy = 1, u%Xh%ly
             do ix = 1, u%Xh%lx
                xl(1) = u%dof%x(ix, iy, iz, iel)
                xl(2) = u%dof%y(ix, iy, iz, iel)
                fcoeff(1) = 3.0e4_rp
                fcoeff(2) = -1.5e3_rp
                fcoeff(3) = 0.5e5_rp
                u%x(ix, iy, iz, iel) = math_ran_dst(ix, iy, 1, iel, xl, &
                     fcoeff) * 1.0e-08_rp
                fcoeff(1) = 2.3e4_rp
                fcoeff(2) = 2.3e3_rp
                fcoeff(3) = -2.0e5_rp
                v%x(ix, iy, iz, iel) = math_ran_dst(ix, iy, 1, iel, xl, &
                     fcoeff) * 1.0e-08_rp
                if (if_2d) then
                   w%x(ix, iy, iz, iel) = 0.0_rp
                else
                fcoeff(1) = 0.5e4_rp
                fcoeff(2) = -2.3e3_rp
                fcoeff(3) = 2.0e5_rp
                w%x(ix, iy, iz, iel) = math_ran_dst(ix, iy, 1, iel, xl, &
                     fcoeff) * 1.0e-08_rp
                end if
             end do
          end do
       end do
    end do

    if (if_2d) then
       call z_plane_fix(u)
       call z_plane_fix(v)
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
      call device_memcpy(u%x, u%x_d, u%size(), host_to_device, .true.)
      call device_memcpy(v%x, v%x_d, v%size(), host_to_device, .true.)
      call device_memcpy(w%x, w%x_d, w%size(), host_to_device, .true.)
    end if

  end subroutine rand_ic


  ! The original Nek5000 random number generator is implemented
  ! in @ref ran1. This totally ad-hoc random number generator below
  ! could be preferable to the original one for the simple reason that it
  ! gives the same initial condition independent of the number of
  real(kind=rp) function math_ran_dst(ix, iy, iz, ieg, xl, fcoeff)
    implicit none
    integer ix, iy, iz, ieg
    real(kind=rp) :: fcoeff(3), xl(2)

    math_ran_dst = fcoeff(1)*(ieg+xl(1)*sin(xl(2))) + &
         fcoeff(2)*ix*iy + fcoeff(3)*ix
    math_ran_dst = 1.0e3_rp * sin(math_ran_dst)
    math_ran_dst = 1.0e3_rp * sin(math_ran_dst)
    math_ran_dst = cos(math_ran_dst)

    return
  end function math_ran_dst

   subroutine state_vector_write(self, idx)
      class(state_vector_t), intent(inout) :: self
      integer, intent(in) :: idx
      type(fld_file_output_t) :: output

      call output%init(sp, 'state', 4)
      call output%fields%assign_to_field(1, self%p)
      call output%fields%assign_to_field(2, self%u)
      call output%fields%assign_to_field(3, self%v)
      call output%fields%assign_to_field(4, self%w)
      call output%set_counter(idx - 1)
      call output%sample(real(idx, kind=rp))
      call output%free()
   end subroutine state_vector_write

   subroutine state_vector_attach_coef(self)
      class(state_vector_t), intent(inout) :: self

      if (.not. associated(neko_user_access%case)) then
         if (.not. associated(self%coef)) then
            error stop 'Neko user access is not initialized!'
         end if
      else
         self%coef => neko_user_access%case%fluid%c_Xh
      end if

      self%if_2d = (self%coef%msh%gdim .eq. 2)
   end subroutine state_vector_attach_coef

   logical function state_vector_is_initialized(self)
      class(state_vector_t), intent(inout) :: self

      state_vector_is_initialized = .false.

      if (.not. associated(self%coef)) return
      if (.not. allocated(self%u%x)) return
      if (.not. allocated(self%v%x)) return
      if (.not. allocated(self%w%x)) return
      if (.not. allocated(self%p%x)) return
      if (.not. associated(self%u%dof, self%coef%dof)) return
      if (.not. associated(self%v%dof, self%coef%dof)) return
      if (.not. associated(self%w%dof, self%coef%dof)) return
      if (.not. associated(self%p%dof, self%coef%dof)) return

      if (NEKO_BCKND_DEVICE .eq. 1) then
         if (.not. device_associated(self%u%x)) return
         if (.not. device_associated(self%v%x)) return
         if (.not. device_associated(self%w%x)) return
         if (.not. device_associated(self%p%x)) return
      end if

      state_vector_is_initialized = .true.
   end function state_vector_is_initialized

   subroutine state_vector_sanitize_field(fld)
      type(field_t), intent(inout) :: fld
      logical :: x_device_associated

      if (NEKO_BCKND_DEVICE .ne. 1) return

      if (.not. allocated(fld%x)) then
         fld%x_d = c_null_ptr
         fld%internal_dofmap = .false.
         fld%name = ''
         nullify(fld%dof)
         nullify(fld%xh)
         nullify(fld%msh)
         return
      end if

      x_device_associated = device_associated(fld%x)
      if (.not. x_device_associated) fld%x_d = c_null_ptr
   end subroutine state_vector_sanitize_field
   
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
 
 end module neko_vector
