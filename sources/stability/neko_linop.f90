module neko_linop
   use LightKrylov, only: abstract_linop_rdp, abstract_vector_rdp
   use num_types, only: rp
   use simulation_m, only: simulation_t
   use field, only: field_t
   use field_math, only: field_copy, field_rzero
   use neko_vector, only: state_vector_t

   implicit none
   type, extends(abstract_linop_rdp), public :: linear_propagator_t
      type(simulation_t), pointer, public :: simulation => null()
      logical :: if_2d = .false.
    contains
      private
      procedure, pass(self), public :: matvec => direct_solver
      procedure, pass(self), public :: rmatvec => direct_solver
      procedure, pass(self), public :: init => neko_propagator_init
      procedure, pass(self), public :: free => neko_propagator_free
   end type linear_propagator_t
 
 contains

   subroutine direct_solver(self, vec_in, vec_out)
      class(linear_propagator_t), intent(inout) :: self
      class(abstract_vector_rdp), intent(in) :: vec_in
      class(abstract_vector_rdp), intent(inout) :: vec_out

      select type (vec_in)
      type is (state_vector_t)
         select type (vec_out)
         type is (state_vector_t)
            call self%simulation%reset_adjoint()

            call field_copy(self%simulation%adjoint_case%fluid_adj%u_adj, &
                 vec_in%u)
            call field_copy(self%simulation%adjoint_case%fluid_adj%v_adj, &
                 vec_in%v)
            call field_copy(self%simulation%adjoint_case%fluid_adj%w_adj, &
                 vec_in%w)
            call field_copy(self%simulation%adjoint_case%fluid_adj%p_adj, &
                 vec_in%p)

            call self%simulation%run_backward()

            call vec_out%init()
            call field_copy(vec_out%u, &
                 self%simulation%adjoint_case%fluid_adj%u_adj)
            call field_copy(vec_out%v, &
                 self%simulation%adjoint_case%fluid_adj%v_adj)
            call field_copy(vec_out%w, &
                 self%simulation%adjoint_case%fluid_adj%w_adj)
            call field_copy(vec_out%p, &
                 self%simulation%adjoint_case%fluid_adj%p_adj)

            if (self%if_2d) then
               call z_plane_fix(vec_out%u)
               call z_plane_fix(vec_out%v)
               call field_rzero(vec_out%w)
               call z_plane_fix(vec_out%p)
            end if
         end select
      end select
   end subroutine direct_solver

   subroutine neko_propagator_init(self, simulation)
     class(linear_propagator_t), intent(inout) :: self
     type(simulation_t), intent(in), target :: simulation
     real(kind=rp) :: T_fin, dt
     integer :: n_steps

     self%simulation => simulation
     T_fin = self%simulation%neko_case%time%end_time
     dt = self%simulation%neko_case%time%dt
     n_steps = int(T_fin/dt)
     self%simulation%n_timesteps = n_steps
     if (self%simulation%neko_case%fluid%c_Xh%msh%gdim .eq. 2) then
        self%if_2d = .true.
     end if
   end subroutine neko_propagator_init

   subroutine neko_propagator_free(self)
     class(linear_propagator_t), intent(inout) :: self

     nullify(self%simulation)
     self%if_2d = .false.
   end subroutine neko_propagator_free

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
 
 end module neko_linop
