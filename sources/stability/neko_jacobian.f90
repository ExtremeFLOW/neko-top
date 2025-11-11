module neko_jacobian
   use LightKrylov, only: abstract_vector_rdp, abstract_jacobian_linop_rdp
   use field, only: field_t
   use field_math, only: field_copy, field_rzero
   use neko_vector, only: state_vector_t
   use simulation_m, only: simulation_t

   implicit none
   type, extends(abstract_jacobian_linop_rdp), public :: jacobian_t
      type(simulation_t), pointer, public :: simulation
      logical :: if_2d = .false.
    contains
      private
      procedure, pass(self), public :: matvec => linear_map
      procedure, pass(self), public :: rmatvec => linear_map
      procedure, pass(self), public :: init => linear_propagator_init
      procedure, pass(self), public :: free => linear_propagator_free
   end type jacobian_t

 
 contains


   subroutine linear_map(self, vec_in, vec_out)
     ! Linear Operator.
     class(jacobian_t), intent(inout)  :: self
     ! Input vector.
     class(abstract_vector_rdp) , intent(in)  :: vec_in
     ! Output vector.
     class(abstract_vector_rdp) , intent(out) :: vec_out
 
     select type(vec_in)
     type is(state_vector_t)
        select type(vec_out)
        type is(state_vector_t)
           ! Reset propagator.
           call self%simulation%reset_adjoint()

           ! Get state vector.
           call field_copy(self%simulation%adjoint_case%fluid_adj%u_adj, vec_in%u)
           call field_copy(self%simulation%adjoint_case%fluid_adj%v_adj, vec_in%v)
           call field_copy(self%simulation%adjoint_case%fluid_adj%w_adj, vec_in%w)
           call field_copy(self%simulation%adjoint_case%fluid_adj%p_adj, vec_in%p)

           ! Integrate forward in time.
           call self%simulation%run_backward()

           ! Pass-back the state vector.
           call vec_out%init()
           call field_copy(vec_out%u, self%simulation%adjoint_case%fluid_adj%u_adj)
           call field_copy(vec_out%v, self%simulation%adjoint_case%fluid_adj%v_adj)
           call field_copy(vec_out%w, self%simulation%adjoint_case%fluid_adj%w_adj)
           call field_copy(vec_out%p, self%simulation%adjoint_case%fluid_adj%p_adj)

           ! quasi 2D hack...
           if (self%if_2d) then
           call z_plane_fix(vec_out%u)
           call z_plane_fix(vec_out%v)
           call field_rzero(vec_out%w)
           call z_plane_fix(vec_out%p)
           end if

           ! Evaluate residual F(X) - X.
           call vec_out%sub(vec_in)
        end select
     end select
     return
   end subroutine linear_map

   subroutine linear_propagator_init(self, simulation)
     ! Linear Operator.
     class(jacobian_t), intent(inout)  :: self
     type(simulation_t), intent(in), target :: simulation
     self%simulation => simulation
     if(self%simulation%neko_case%fluid%c_Xh%msh%gdim .eq. 2) then
             self%if_2d = .true.
     end if
     return
   end subroutine linear_propagator_init

   subroutine linear_propagator_free(self)
     ! Linear Operator.
     class(jacobian_t), intent(inout)  :: self
     nullify(self%simulation)
     return
   end subroutine linear_propagator_free

   subroutine z_plane_fix(fld)
  type(field_t), intent(inout) :: fld
  integer :: iel, iz, iy, ix, nel
  ! note this wont work on GPUs

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
 
 end module neko_jacobian
