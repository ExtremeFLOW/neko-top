module cylinder
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov, only: wp => dp
   ! Standard Library.
   use stdlib_math, only : linspace
   use stdlib_optval, only : optval
   ! Additional neko-top libraries (to solve linearized and adjoint)
   use num_types, only : rp
   use neko, only: neko_init, neko_finalize, neko_solve
   use case, only : case_t
   use adjoint_case, only: adjoint_case_t, adjoint_init, adjoint_free
   use simulation_adjoint, only: solve_adjoint
   use json_module, only: json_file
   ! Specific to the state vector
   use field, only: field_t
   use coefs, only: coef_t
   ! various imports for some specific neko operations
   use neko_config, only : NEKO_BCKND_DEVICE
   use scratch_registry, only: neko_scratch_registry
   use field_math, only: field_col3, field_rzero, field_cmult, field_add3s2, &
      field_copy
   use math, only: glsc2
   use device_math, only: device_glsc2
   use device, only : device_memcpy, HOST_TO_DEVICE
   implicit none
 
   character*128, parameter, private :: this_module = 'cylinder'
 
   !------------------------------
   !-----     PARAMETERS     -----
   !------------------------------

   ! generally parameters enter via the .case file

   !-------------------------------------------
   !-----     LIGHTKRYLOV VECTOR TYPE     -----
   !-------------------------------------------
 
   type, extends(abstract_vector_rdp), public :: state_vector
      ! adjoint velocity fields
      type(field_t), allocatable :: u
      type(field_t), allocatable :: v
      type(field_t), allocatable :: w
      ! we need the mass matrix to integrate fields..
     type(coef_t), pointer :: coef => null()
    contains
      private
      procedure, pass(self), public :: zero
      procedure, pass(self), public :: dot
      procedure, pass(self), public :: scal
      procedure, pass(self), public :: axpby
      procedure, pass(self), public :: rand
      procedure, pass(self), public :: get_size
      ! we also want a clean way to initialize and free
      procedure, pass(self), public :: init => state_vector_init
      procedure, pass(self), public :: free => state_vector_free
   end type state_vector
 
   !-----------------------------------
   !-----     NEKO PROPAGATOR     -----
   !-----------------------------------
 
   type, extends(abstract_linop_rdp), public :: neko_propagator
      ! The primal case
      type(case_t), public :: neko_case
      ! The linear case 
      ! (I know the name is stupid here... I've had to hack around the fact we
      ! only have an adjoint type, not "perturb" as in Nek5000)
      type(adjoint_case_t), public :: linear_case
      ! The adjoint case
      type(adjoint_case_t), public :: adjoint_case
    contains
      private
      procedure, pass(self), public :: matvec => direct_solver
      procedure, pass(self), public :: rmatvec => adjoint_solver
      ! we also want a clean way to initialize, free and reset
      procedure, pass(self), public :: init => neko_propagator_init
      procedure, pass(self), public :: free => neko_propagator_free
   end type neko_propagator
 
 contains
 
   !======================================================================
   !======================================================================
   !=====                                                            =====
   !=====     PHYSICAL MODEL : Linearize NAVIER-STOKES EQUATIONS     =====
   !=====                                                            =====
   !======================================================================
   !======================================================================
 
    subroutine state_vector_init(self, coef)
     class(state_vector), intent(inout) :: self
     type(coef_t), intent(in), target :: coef

     ! pass coef so we can integrate
     self%coef => coef

     ! allocate fields
     allocate(self%u)
     allocate(self%v)
     allocate(self%w)
     call self%u%init(self%coef%dof, fld_name = "state_u")
     call self%v%init(self%coef%dof, fld_name = "state_v")
     call self%w%init(self%coef%dof, fld_name = "state_w")

     return
   end subroutine state_vector_init
 
   !=========================================================
   !=========================================================
   !=====                                               =====
   !=====     LIGHTKRYLOV MANDATORY IMPLEMENTATIONS     =====
   !=====                                               =====
   !=========================================================
   !=========================================================
 
   !----------------------------------------------------
   !-----     TYPE-BOUND PROCEDURE FOR VECTORS     -----
   !----------------------------------------------------
 
   subroutine zero(self)
     class(state_vector), intent(inout) :: self
     call field_rzero(self%u)
     call field_rzero(self%v)
     call field_rzero(self%w)
     return
   end subroutine zero
 
   real(kind=wp) function dot(self, vec) result(alpha)
     class(state_vector)   , intent(in) :: self
     class(abstract_vector_rdp), intent(in) :: vec
     ! a working array
     type(field_t), pointer :: work
     integer :: temp_indices(1)
     integer :: n
     real(kind=rp) :: tmp_real
     select type(vec)
     type is(state_vector)
        ! here we're going to take an energy norm I guess...
        call neko_scratch_registry%request_field(work, temp_indices(1))
        call field_rzero(work)
        call field_col3(work, self%u, vec%u)
        call field_col3(work, self%v, vec%v)
        call field_col3(work, self%w, vec%w)
        ! integrate 
        ! (I guess the GPU backend doesn't matter so much here...)
        n = work%size()
        if (NEKO_BCKND_DEVICE .eq. 1) then
           tmp_real = device_glsc2(work%x_d, self%coef%B_d, n)
       else
           tmp_real = glsc2(work%x, self%coef%B, n)
       end if
       tmp_real = tmp_real * 0.5_rp
       alpha = real(tmp_real, wp)

       call neko_scratch_registry%relinquish_field(temp_indices)
     end select
     return
   end function dot
 
   subroutine scal(self, alpha)
     class(state_vector), intent(inout) :: self
     real(kind=wp)      , intent(in)    :: alpha
     real(kind=rp) :: tmp_real
     tmp_real = real(alpha, rp)
     call field_cmult(self%u, tmp_real)
     call field_cmult(self%v, tmp_real)
     call field_cmult(self%w, tmp_real)
     return
   end subroutine scal
 
   subroutine axpby(self, alpha, vec, beta)
     class(state_vector)   , intent(inout) :: self
     class(abstract_vector_rdp), intent(in)    :: vec
     real(kind=wp)         , intent(in)    :: alpha, beta
     select type(vec)
     type is(state_vector)
        call field_add3s2(self%u, self%u, vec%u, alpha, beta)
        call field_add3s2(self%v, self%v, vec%v, alpha, beta)
        call field_add3s2(self%w, self%w, vec%w, alpha, beta)
     end select
     return
   end subroutine axpby
 
   integer function get_size(self) result(N)
     class(state_vector), intent(in) :: self
     ! hmmm self is a bit confusing. I assume you mean the TOTAL size, ie,
     ! number of GLL pts for all 3 components...
     N = self%u%size() + self%v%size() + self%w%size()
     return
   end function get_size
 
   subroutine rand(self, ifnorm)
     class(state_vector), intent(inout) :: self
     logical, optional,   intent(in)    :: ifnorm
     ! internals
     logical :: normalize
     real(kind=wp) :: alpha
     normalize = optval(ifnorm,.true.)
     call rand_ic(self%u, self%v, self%w)
     if (normalize) then
       alpha = self%norm()
       call self%scal(1.0_wp/alpha)
     endif
     return
   end subroutine rand

   !---------------------------------------------------
   !-----     EXTRA PROCEDURES FOR THE VECTOR     -----
   !---------------------------------------------------

   subroutine state_vector_free(self)
     class(state_vector), intent(inout) :: self

     call self%u%free()
     call self%v%free()
     call self%w%free()
     call self%coef%free()

     return
   end subroutine state_vector_free

  ! User defined initial condition
  subroutine rand_ic(u, v, w)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
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
                u%x(ix, iy, iz, iel) = math_ran_dst(ix, iy, iz, iel, xl, &
                     fcoeff) * 1.0e-08_rp
                fcoeff(1) = 2.3e4_rp
                fcoeff(2) = 2.3e3_rp
                fcoeff(3) = -2.0e5_rp
                v%x(ix, iy, iz, iel) = math_ran_dst(ix, iy, iz, iel, xl, &
                     fcoeff) * 1.0e-08_rp
                w%x(ix, iy, iz, iel) = 0.0_rp
             end do
          end do
       end do
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
      call device_memcpy(u%x, u%x_d, u%size(), host_to_device, .true.)
      call device_memcpy(v%x, v%x_d, v%size(), host_to_device, .true.)
      call device_memcpy(w%x, w%x_d, w%size(), host_to_device, .true.)
    end if

  end subroutine rand_ic


  ! The original Nek5000 random number generator is implementted
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
 
   !------------------------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURES FOR THE EXPONENTIAL PROPAGATOR     -----
   !------------------------------------------------------------------------
 
   subroutine direct_solver(self, vec_in, vec_out)
     ! Linear Operator.
     class(neko_propagator), intent(in)  :: self
     ! Input vector.
     class(abstract_vector_rdp) , intent(in)  :: vec_in
     ! Output vector.
     class(abstract_vector_rdp) , intent(out) :: vec_out
 
     select type(vec_in)
     type is(state_vector)
        select type(vec_out)
        type is(state_vector)
           ! Reset propagator.
           ! @todo
           ! Get state vector.
           ! (again... the naming with "adjoint" isn't smart here)
           call field_copy(self%linear_case%fluid_adj%u_adj, vec_in%u)
           call field_copy(self%linear_case%fluid_adj%v_adj, vec_in%v)
           call field_copy(self%linear_case%fluid_adj%w_adj, vec_in%w)
           ! Integrate forward in time.
           call solve_wrapper(self%linear_case)
           ! Pass-back the state vector.
           call field_copy(vec_out%u, self%linear_case%fluid_adj%u_adj)
           call field_copy(vec_out%v, self%linear_case%fluid_adj%v_adj)
           call field_copy(vec_out%w, self%linear_case%fluid_adj%w_adj)
        end select
     end select
     return
   end subroutine direct_solver
 
   subroutine adjoint_solver(self, vec_in, vec_out)
     ! Linear Operator.
     class(neko_propagator), intent(in)  :: self
     ! Input vector.
     class(abstract_vector_rdp) , intent(in)  :: vec_in
     ! Output vector.
     class(abstract_vector_rdp) , intent(out) :: vec_out
 
     select type(vec_in)
     type is(state_vector)
        select type(vec_out)
        type is(state_vector)
           ! Reset propagator.
           ! @todo
           ! Get state vector.
           ! (again... the naming with "adjoint" isn't smart here)
           call field_copy(self%adjoint_case%fluid_adj%u_adj, vec_in%u)
           call field_copy(self%adjoint_case%fluid_adj%v_adj, vec_in%v)
           call field_copy(self%adjoint_case%fluid_adj%w_adj, vec_in%w)
           ! Integrate forward in time.
           call solve_wrapper(self%linear_case)
           ! Pass-back the state vector.
           call field_copy(vec_out%u, self%adjoint_case%fluid_adj%u_adj)
           call field_copy(vec_out%v, self%adjoint_case%fluid_adj%v_adj)
           call field_copy(vec_out%w, self%adjoint_case%fluid_adj%w_adj)
        end select
     end select
     return
   end subroutine adjoint_solver

   !-------------------------------------------------------------------
   !-----     EXTRA PROCEDURES FOR THE EXPONENTIAL PROPAGATOR     -----
   !-------------------------------------------------------------------

   subroutine neko_propagator_init(self)
     ! Linear Operator.
     class(neko_propagator), intent(inout)  :: self
     ! type(json_file), intent(inout) :: parameters

     ! initialize the "baseflow"
     call neko_init(self%neko_case)
     ! initialize the linear
     call adjoint_init(self%linear_case, self%neko_case)
     ! initialize the adjoint
     call adjoint_init(self%adjoint_case, self%neko_case)

     ! TODO
     ! Load in the baseflow

     return
   end subroutine neko_propagator_init

   subroutine neko_propagator_free(self)
     ! Linear Operator.
     class(neko_propagator), intent(inout)  :: self

     call adjoint_free(self%adjoint_case)
     call adjoint_free(self%linear_case)
     call neko_finalize(self%neko_case)
     return
   end subroutine neko_propagator_free

   ! silly little wrapper to ignore intent.
   subroutine solve_wrapper(self)
     ! Linear Operator.
     class(adjoint_case_t)  :: self

     call solve_adjoint(self)

     return
   end subroutine solve_wrapper
 
 end module cylinder
