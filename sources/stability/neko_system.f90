module neko_system
   use LightKrylov, only: abstract_system_rdp, abstract_vector_rdp
   use LightKrylov, only: wp => dp
   use simulation_m, only: simulation_t
   use field, only: field_t
   use field_math, only: field_copy, field_rzero
   use neko_vector, only: state_vector_t

   implicit none
   type, extends(abstract_system_rdp), public :: non_linear_propagator_t
      type(simulation_t), pointer, public :: simulation => null()
      logical :: if_2d = .false.
   contains
      private
      procedure, pass(self), public :: response => nonlinear_map
      procedure, pass(self), public :: init => non_linear_propagator_init
      procedure, pass(self), public :: free => non_linear_propagator_free
   end type non_linear_propagator_t
 
 contains
 
   subroutine nonlinear_map(self, vec_in, vec_out, atol)
     class(non_linear_propagator_t), intent(inout) :: self
     class(abstract_vector_rdp), intent(in) :: vec_in
     class(abstract_vector_rdp) , intent(inout) :: vec_out
     real(wp), intent(in) :: atol

     select type (vec_in)
     type is (state_vector_t)
        select type (vec_out)
        type is (state_vector_t)
           call self%simulation%reset_forward()

           call field_copy(self%simulation%neko_case%fluid%u, vec_in%u)
           call field_copy(self%simulation%neko_case%fluid%v, vec_in%v)
           call field_copy(self%simulation%neko_case%fluid%w, vec_in%w)
           call field_copy(self%simulation%neko_case%fluid%p, vec_in%p)

           call self%simulation%run_forward()

           call vec_out%init()
           call field_copy(vec_out%u, self%simulation%neko_case%fluid%u)
           call field_copy(vec_out%v, self%simulation%neko_case%fluid%v)
           call field_copy(vec_out%w, self%simulation%neko_case%fluid%w)
           call field_copy(vec_out%p, self%simulation%neko_case%fluid%p)

           if (self%if_2d) then
              call z_plane_fix(vec_out%u)
              call z_plane_fix(vec_out%v)
              call z_plane_fix(vec_out%p)
              call field_rzero(vec_out%w)
           end if

           call vec_out%sub(vec_in)
        end select
     end select
   end subroutine nonlinear_map

   subroutine non_linear_propagator_init(self, simulation)
     class(non_linear_propagator_t), intent(inout) :: self
     type(simulation_t), target, intent(in) :: simulation

     self%simulation => simulation
     if (self%simulation%neko_case%fluid%c_Xh%msh%gdim .eq. 2) then
        self%if_2d = .true.
     end if
   end subroutine non_linear_propagator_init

   subroutine non_linear_propagator_free(self)
     class(non_linear_propagator_t), intent(inout) :: self

     nullify(self%simulation)
     self%if_2d = .false.
   end subroutine non_linear_propagator_free

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
 
 end module neko_system
