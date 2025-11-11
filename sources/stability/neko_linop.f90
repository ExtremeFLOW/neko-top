module neko_linop
   use LightKrylov, only: abstract_linop_rdp, abstract_vector_rdp
   use num_types, only : rp, sp
   use simulation_m, only : simulation_t
   use field, only : field_t
   use field_math, only : field_copy
   use fld_file_output, only : fld_file_output_t
   use neko_vector, only : state_vector_t

   implicit none
   type, extends(abstract_linop_rdp), public :: linear_propagator_t
      type(simulation_t), pointer, public :: simulation
      type(fld_file_output_t), public :: output_primal
      type(fld_file_output_t), public :: output_linear
      type(fld_file_output_t), public :: output_adjoint
    contains
      private
      procedure, pass(self), public :: matvec => direct_solver
      procedure, pass(self), public :: rmatvec => direct_solver
      procedure, pass(self), public :: init => neko_propagator_init
      procedure, pass(self), public :: free => neko_propagator_free
      procedure, pass(self), public :: write_linear => write_linear_wrapper
      procedure, pass(self), public :: write_adjoint => write_adjoint_wrapper
   end type linear_propagator_t
 
 contains

   subroutine direct_solver(self, vec_in, vec_out)
     class(linear_propagator_t), intent(inout)  :: self
     class(abstract_vector_rdp) , intent(in)  :: vec_in
     class(abstract_vector_rdp) , intent(out) :: vec_out
 
     select type(vec_in)
     type is(state_vector_t)
        select type(vec_out)
        type is(state_vector_t)
           ! Reset propagator.
           call self%simulation%reset_adjoint()

           ! Set state vector.
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
           call z_plane_fix(vec_out%u)
           call z_plane_fix(vec_out%v)
           call z_plane_fix(vec_out%w)
           call z_plane_fix(vec_out%p)
        end select
     end select
     return
   end subroutine direct_solver

   subroutine neko_propagator_init(self, simulation)
     ! Linear Operator.
     class(linear_propagator_t), intent(inout)  :: self
     type(simulation_t), intent(in), target :: simulation
     real(kind=rp) :: T_fin, dt 
     integer :: n_steps

     self%simulation => simulation
     T_fin = self%simulation%neko_case%time%end_time
     dt = self%simulation%neko_case%time%dt
     n_steps = int(T_fin/dt)
     self%simulation%n_timesteps = n_steps
     
     ! NOTE baseflow should be loaded via IC in .case file, but let's double
     ! check
          call self%output_primal%init(sp, 'checking_base', 4)
          call self%output_primal%fields%assign(2, self%simulation%neko_case%fluid%u)
          call self%output_primal%fields%assign(3, self%simulation%neko_case%fluid%v)
          call self%output_primal%fields%assign(4, self%simulation%neko_case%fluid%w)
          call self%output_primal%fields%assign(1, self%simulation%neko_case%fluid%p)
          call self%output_primal%sample(0.0_rp)
     ! Assign samplers for the forward and adjoint in case we want to look at
     ! them (debugging)
          call self%output_linear%init(sp, 'checking_linear', 4)
          call self%output_linear%fields%assign(2, self%simulation%adjoint_case%fluid_adj%u_adj)
          call self%output_linear%fields%assign(3, self%simulation%adjoint_case%fluid_adj%v_adj)
          call self%output_linear%fields%assign(4, self%simulation%adjoint_case%fluid_adj%w_adj)
          call self%output_linear%fields%assign(1, self%simulation%adjoint_case%fluid_adj%p_adj)
          call self%output_adjoint%init(sp, 'checking_adjoint', 4)
          call self%output_adjoint%fields%assign(2, self%simulation%adjoint_case%fluid_adj%u_adj)
          call self%output_adjoint%fields%assign(3, self%simulation%adjoint_case%fluid_adj%v_adj)
          call self%output_adjoint%fields%assign(4, self%simulation%adjoint_case%fluid_adj%w_adj)
          call self%output_adjoint%fields%assign(1, self%simulation%adjoint_case%fluid_adj%p_adj)

     return
   end subroutine neko_propagator_init

   subroutine neko_propagator_free(self)
     ! Linear Operator.
     class(linear_propagator_t), intent(inout)  :: self

     call self%simulation%free()
     return
   end subroutine neko_propagator_free

   ! silly little wrapper to ignore intent.
   subroutine write_linear_wrapper(self, idx)
     ! Linear Operator.
     class(linear_propagator_t) :: self
     integer :: idx

     call self%output_linear%sample(real(idx, kind=rp))

     return
   end subroutine write_linear_wrapper

   ! silly little wrapper to ignore intent.
   subroutine write_adjoint_wrapper(self, idx)
     ! Linear Operator.
     class(linear_propagator_t) :: self
     integer :: idx

     call self%output_adjoint%sample(real(idx, kind=rp))

     return
   end subroutine write_adjoint_wrapper

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
 
 end module neko_linop
