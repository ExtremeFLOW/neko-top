module Ginzburg_Landau
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov, only: wp => dp
   ! Standard Library.
   use stdlib_math, only : linspace
   use stdlib_optval, only : optval
   ! Additional neko-top libraries
   use simulation_m, only: simulation_t
   use json_module, only: json_file
   use field, only: field_t
   ! specific to some neko operations
   use coef, only: coef_t
   use scratch_registry, only: neko_scratch_registry
   use field_math, only: field_col3, field_rzero, field_cmult, field_add3s2, &
      field_copy
   use math, only: glsc2
   use device_math, only: device_glsc2
   implicit none
 
   character*128, parameter, private :: this_module = 'cylinder'
 
   public :: initialize_parameters
 
   !------------------------------
   !-----     PARAMETERS     -----
   !------------------------------

   ! generally parameters enter via the .case file

   !-------------------------------------------
   !-----     LIGHTKRYLOV VECTOR TYPE     -----
   !-------------------------------------------
 
   type, extends(abstract_vector_rdp), public :: state_vector
      ! adjoint velocity fields
      type(field_t) :: u
      type(field_t) :: v
      type(field_t) :: w
      ! we need the mass matrix to integrate fields..
     type(coef_t), pointer :: coef => null()
    contains
      private
      procedure, pass(this), public :: zero
      procedure, pass(this), public :: dot
      procedure, pass(this), public :: scal
      procedure, pass(this), public :: axpby
      procedure, pass(this), public :: rand
      procedure, pass(this), public :: get_size
      ! we also want a clean way to initialize and free
      procedure, pass(this), public :: init => state_vector_init
      procedure, pass(this), public :: free => state_vector_free
   end type state_vector
 
   !-----------------------------------
   !-----     NEKO PROPAGATOR     -----
   !-----------------------------------
 
   type, extends(abstract_linop_rdp), public :: neko_propagator
      type(simulation_t) :: simulation
    contains
      private
      procedure, pass(this), public :: matvec => direct_solver
      procedure, pass(this), public :: rmatvec => adjoint_solver
      ! we also want a clean way to initialize, free and reset
      procedure, pass(this), public :: init => neko_propagator_init
      procedure, pass(this), public :: free => neko_propagator_free
      procedure, pass(this), public :: reset => neko_propagator_reset
   end type neko_propagator
 
 contains
 
   !======================================================================
   !======================================================================
   !=====                                                            =====
   !=====     PHYSICAL MODEL : Linearize NAVIER-STOKES EQUATIONS     =====
   !=====                                                            =====
   !======================================================================
   !======================================================================
 
   subroutine initialize_parameters(parameters)
     implicit none
     ! json containing the .case file
     type(json_file) :: parameters

    ! initialize the simulation
    call this%state%init(parameters)

    ! load the baseflow
    ! @todo
 
     return
   end subroutine initialize_parameters
 
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
 
   subroutine zero(this)
     class(state_vector), intent(inout) :: this
     call field_rzero(this%u)
     call field_rzero(this%v)
     call field_rzero(this%w)
     return
   end subroutine zero
 
   real(kind=wp) function dot(this, vec) result(alpha)
     class(state_vector)   , intent(in) :: this
     class(abstract_vector_cdp), intent(in) :: vec
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
        call field_col3(work, this%u, vec%u)
        call field_col3(work, this%v, vec%v)
        call field_col3(work, this%w, vec%w)
        ! integrate 
        ! (I guess the GPU backend doesn't matter so much here...)
        n = work%size()
        if (NEKO_BCKND_DEVICE .eq. 1) then
           tmp_real = device_glsc2(work%x_d, this%coef%C_Xh%B_d, n)
       else
           tmp_real = glsc2(work%x, this%coef%C_Xh%B, n)
       end if
       tmp_real = tmp_real * 0.5_rp
       alpha = real(tmp_real, wp)

       call neko_scratch_registry%relinquish_field(temp_indices)
     end select
     return
   end function dot
 
   subroutine scal(this, alpha)
     class(state_vector), intent(inout) :: this
     complex(kind=wp)      , intent(in)    :: alpha
     real(kind=rp) :: tmp_real
     tmp_real = real(alpha, rp)
     call field_cmult(this%u, tmp_real)
     call field_cmult(this%v, tmp_real)
     call field_cmult(this%w, tmp_real)
     return
   end subroutine scal
 
   subroutine axpby(this, alpha, vec, beta)
     class(state_vector)   , intent(inout) :: this
     class(abstract_vector_cdp), intent(in)    :: vec
     complex(kind=wp)         , intent(in)    :: alpha, beta
     select type(vec)
     type is(state_vector)
        call field_add3s2(this%u, this%u, vec%u, alpha, beta)
        call field_add3s2(this%v, this%v, vec%v, alpha, beta)
        call field_add3s2(this%w, this%w, vec%w, alpha, beta)
     end select
     return
   end subroutine axpby
 
   integer function get_size(this) result(N)
     class(state_vector), intent(in) :: this
     ! hmmm this is a bit confusing. I assume you mean the TOTAL size, ie,
     ! number of GLL pts for all 3 components...
     N = this%u%size() + this%v%size() + this%w%size()
     return
   end function get_size
 
   subroutine rand(this, ifnorm)
     class(state_vector), intent(inout) :: this
     logical, optional,   intent(in)    :: ifnorm
     real(kind=wp) :: tmp(nx, 2)
     ! internals
     logical :: normalize
     real(kind=wp) :: alpha
     normalize = optval(ifnorm,.true.)
     call rand_ic(this%u, this%v, this%w)
     if (normalize) then
       alpha = this%norm()
       call this%scal(1.0_wp/alpha)
     endif
     return
   end subroutine rand

   !---------------------------------------------------
   !-----     EXTRA PROCEDURES FOR THE VECTOR     -----
   !---------------------------------------------------

   subroutine state_vector_init(this, coef)
     class(state_vector), intent(inout) :: this
     type(coef_t), intent(in), target :: coef

     ! pass coef so we can integrate
     this%coef => coef

     ! allocate fields
     allocate(this%u)
     allocate(this%v)
     allocate(this%w)
     call this%u%init(this%coef%dm_Xh, fld_name = "state_u")
     call this%v%init(this%coef%dm_Xh, fld_name = "state_v")
     call this%w%init(this%coef%dm_Xh, fld_name = "state_w")

     return
   end subroutine state_vector_init

   subroutine state_vector_free(this)
     class(state_vector), intent(inout) :: this

     call this%u%free()
     call this%v%free()
     call this%w%free()
     call this%coef%free()

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
 
   subroutine direct_solver(this, vec_in, vec_out)
     ! Linear Operator.
     class(neko_propagator), intent(in)  :: this
     ! Input vector.
     class(abstract_vector_cdp) , intent(in)  :: vec_in
     ! Output vector.
     class(abstract_vector_cdp) , intent(out) :: vec_out
 
     ! Time-integrator.
     type(rks54_class) :: prop
     real(kind=wp)     :: dt = 1.0_wp
     real(kind=wp)     :: state_ic(2*nx), state_fc(2*nx)
 
     select type(vec_in)
     type is(state_vector)
        select type(vec_out)
        type is(state_vector)
           ! Reset propagator.
           call this%simulation%reset()
           ! Get state vector.
           call field_copy(this%simulation%adjoint_fluid%u_adj, vec_in%u)
           call field_copy(this%simulation%adjoint_fluid%v_adj, vec_in%v)
           call field_copy(this%simulation%adjoint_fluid%w_adj, vec_in%w)
           ! Integrate forward in time.
           call thi


           ! Get state vector.
           state_ic(:nx) = vec_in%state%re
           state_ic(nx+1:) = vec_in%state%im
           ! Initialize propagator.
           call prop%initialize(n=2*nx, f=rhs)
           ! Integrate forward in time.
           call prop%integrate(0.0_wp, state_ic, dt, this%tau, state_fc)
           ! Pass-back the state vector.
           vec_out%state%re = state_fc(:nx)
           vec_out%state%im = state_fc(nx+1:)
        end select
     end select
     return
   end subroutine direct_solver
 
   subroutine adjoint_solver(this, vec_in, vec_out)
     ! Linear Operator.
     class(neko_propagator), intent(in)  :: this
     ! Input vector.
     class(abstract_vector_cdp) , intent(in)  :: vec_in
     ! Output vector.
     class(abstract_vector_cdp) , intent(out) :: vec_out
 
     ! Time-integrator.
     type(rks54_class) :: prop
     real(kind=wp)     :: dt = 1.0_wp
     real(kind=wp)     :: state_ic(2*nx), state_fc(2*nx)
 
     select type(vec_in)
     type is(state_vector)
        select type(vec_out)
        type is(state_vector)
           ! Get the state.
           state_ic(:nx) = vec_in%state%re
           state_fc(nx+1:) = vec_in%state%im
           ! Initialize propagator.
           call prop%initialize(n=2*nx, f=adjoint_rhs)
           ! Integrate forward in time.
           call prop%integrate(0.0_wp, state_ic, dt, this%tau, state_fc)
           ! Pass-back the state.
           vec_out%state%re = state_fc(:nx)
           vec_out%state%im = state_fc(nx+1:)
        end select
     end select
     return
   end subroutine adjoint_solver

   ! some extra stuff


 
 end module Ginzburg_Landau
