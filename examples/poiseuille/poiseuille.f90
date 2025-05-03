module poiseuille
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov, only: wp => dp
   ! Standard Library.
   use stdlib_math, only : linspace
   use stdlib_optval, only : optval
   ! Additional neko-top libraries (to solve linearized and adjoint)
   use num_types, only : rp, dp
   use neko, only: neko_init, neko_finalize, neko_solve
   use case, only : case_t
   use adjoint_case, only: adjoint_case_t, adjoint_init, adjoint_free
   use simulation_adjoint, only: solve_adjoint, adjoint_reset
   use json_module, only: json_file
   ! Specific to the state vector
   use field, only: field_t
   use coefs, only: coef_t
   use dofmap, only : dofmap_t
   use fld_file_output, only: fld_file_output_t
   ! various imports for some specific neko operations
   use neko_config, only : NEKO_BCKND_DEVICE
   use scratch_registry, only: neko_scratch_registry
   use field_math, only: field_addcol3, field_rzero, field_cmult, &
      field_add2s2, field_copy
   use math, only: glsc3
   use device_math, only: device_glsc3
   use device, only : device_memcpy, HOST_TO_DEVICE
   use gather_scatter, only: gs_t, GS_OP_ADD
   use comm, only: pe_rank
   ! This one is silly... But we need coef to initialize fields
   use global_coef, only: global_coef_t, global_coef_getter
   ! sponges
   use sponge_source_term, only: sponge_source_term_t
   ! debugging, my own eigs
   use LightKrylov_Constants
   use LightKrylov_Utils
       use LightKrylov_Logger, only: log_warning, log_error, log_message, log_information, &
    &                             log_debug, stop_error, check_info, type_error
   implicit none
 
   character*128, parameter, private :: this_module = 'poiseuille'
 
   !------------------------------
   !-----     PARAMETERS     -----
   !------------------------------

   ! generally parameters enter via the .case file

   !-------------------------------------------
   !-----     LIGHTKRYLOV VECTOR TYPE     -----
   !-------------------------------------------
 
   type, extends(abstract_vector_rdp), public :: state_vector
      ! adjoint velocity fields
      type(field_t), pointer :: u => null()
      type(field_t), pointer :: v => null()
      type(field_t), pointer :: w => null()
      type(field_t), pointer :: p => null()
      ! we need the mass matrix to integrate fields..
      type(coef_t), pointer :: coef => null()
      logical :: initialized = .false.
      ! a way to spy on the vector (mostly for debugging)
      type(fld_file_output_t), public :: output
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
      procedure, pass(self), public :: copy => state_vector_copy
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
      ! a way to spy on the vector (mostly for debugging)
      type(fld_file_output_t), public :: output_primal
      type(fld_file_output_t), public :: output_linear
      type(fld_file_output_t), public :: output_adjoint
    contains
      private
      procedure, pass(self), public :: matvec => direct_solver
      procedure, pass(self), public :: rmatvec => adjoint_solver
      ! we also want a clean way to initialize, free and reset
      procedure, pass(self), public :: init => neko_propagator_init
      procedure, pass(self), public :: free => neko_propagator_free
      procedure, pass(self), public :: write_linear => write_linear_wrapper
      procedure, pass(self), public :: write_adjoint => write_adjoint_wrapper
   end type neko_propagator
 
 contains
 
   !======================================================================
   !======================================================================
   !=====                                                            =====
   !=====     PHYSICAL MODEL : Linearize NAVIER-STOKES EQUATIONS     =====
   !=====                                                            =====
   !======================================================================
   !======================================================================
 
    subroutine state_vector_init(self)
     class(state_vector), intent(inout) :: self

     ! Check for global coef
     if (.not. self%initialized) then
       if (.not. associated(global_coef_getter)) then
          error stop "No global coef set!"
       end if
          ! Take the global coef
          self%coef => global_coef_getter%global_coef

          ! initialize fields
          allocate(self%u)
          allocate(self%v)
          allocate(self%w)
          allocate(self%p)
          call self%u%init(self%coef%dof, fld_name = "state_u")
          call self%v%init(self%coef%dof, fld_name = "state_v")
          call self%w%init(self%coef%dof, fld_name = "state_w")
          call self%p%init(self%coef%dof, fld_name = "state_p")

          ! initialize the sampler so we can spy in if needed
          call self%output%init(sp, 'state', 4)
          call self%output%fields%assign(2, self%u)
          call self%output%fields%assign(3, self%v)
          call self%output%fields%assign(4, self%w)
          call self%output%fields%assign(1, self%p)

          ! done
          self%initialized = .true.
     end if

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
     
     ! always try initializing
     call self%init()

     call field_rzero(self%u)
     call field_rzero(self%v)
     call field_rzero(self%w)
     call field_rzero(self%p)
     return
   end subroutine zero
 
   real(kind=wp) function dot(self, vec) result(alpha)
     class(state_vector)   , intent(in) :: self
     class(abstract_vector_rdp), intent(in) :: vec
     integer :: n
     real(kind=rp) :: alpha_rp
     select type(vec)
     type is(state_vector)
        ! always try initializing
        call state_vector_init_wrapper(self)
        call state_vector_init_wrapper(vec)

        ! here we're going to take an energy norm I guess...
        n = self%u%size()
        if (NEKO_BCKND_DEVICE .eq. 1) then
           alpha_rp = device_glsc3(self%u%x_d, vec%u%x_d, self%coef%B_d, n)
           alpha_rp = alpha_rp + device_glsc3(self%v%x_d, vec%v%x_d, self%coef%B_d, n)
           alpha_rp = alpha_rp + device_glsc3(self%w%x_d, vec%w%x_d, self%coef%B_d, n)
       else
           alpha_rp = glsc3(self%u%x, vec%u%x, self%coef%B, n)
           alpha_rp = alpha_rp + glsc3(self%v%x, vec%v%x, self%coef%B, n)
           alpha_rp = alpha_rp + glsc3(self%w%x, vec%w%x, self%coef%B, n)
       end if
       alpha_rp = alpha_rp * 0.5_rp
       alpha = real(alpha_rp, wp)

     end select
     return
   end function dot
 
   subroutine scal(self, alpha)
     class(state_vector), intent(inout) :: self
     real(kind=wp)      , intent(in)    :: alpha
     real(kind=rp) :: alpha_rp
     ! always try initializing
     call state_vector_init_wrapper(self)

     alpha_rp = real(alpha, rp)
     call field_cmult(self%u, alpha_rp)
     call field_cmult(self%v, alpha_rp)
     call field_cmult(self%w, alpha_rp)
     call field_cmult(self%p, alpha_rp)
     return
   end subroutine scal
 
   subroutine axpby(alpha, vec, beta, self)
     class(state_vector)   , intent(inout) :: self
     class(abstract_vector_rdp), intent(in)    :: vec
     real(kind=wp)         , intent(in)    :: alpha, beta
     real(kind=rp) :: alpha_rp, beta_rp
     select type(vec)
     type is(state_vector)
        ! always try initializing
        call state_vector_init_wrapper(self)
        call state_vector_init_wrapper(vec)

        alpha_rp = real(alpha, kind=rp)
        beta_rp = real(beta, kind=rp)

        if (pe_rank.eq.0) then
        print *, "axbpy, alpha = ", alpha_rp, " beta = ", beta_rp
        end if

        ! be careful with the order here !
        ! notice the axpby(alpha, vec, beta, self)
        ! in Ginzberg_landau we have:
        ! self%state = beta*self%state + alpha*vec%state
        call self%scal(beta_rp)
        call field_add2s2(self%u, vec%u, alpha_rp)
        call field_add2s2(self%v, vec%v, alpha_rp)
        call field_add2s2(self%w, vec%w, alpha_rp)
        call field_add2s2(self%p, vec%p, alpha_rp)
     end select
     return
   end subroutine axpby
 
   integer function get_size(self) result(N)
     class(state_vector), intent(in) :: self
     ! always try initializing
     call state_vector_init_wrapper(self)
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
     
     ! always try initializing
     call self%init()

     normalize = optval(ifnorm,.true.)
     call rand_ic(self%u, self%v, self%w)

    ! enforce continuity across the field
    call self%coef%gs_h%op(self%u, GS_OP_ADD)
    call self%coef%gs_h%op(self%v, GS_OP_ADD)
    call self%coef%gs_h%op(self%w, GS_OP_ADD)

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
     call self%p%free()
     call self%coef%free()

     return
   end subroutine state_vector_free

   subroutine state_vector_copy(self, vec)
     class(state_vector), intent(inout) :: self
     class(state_vector), intent(in) :: vec

     call field_copy(self%u, vec%u)
     call field_copy(self%v, vec%v)
     call field_copy(self%w, vec%w)
     call field_copy(self%p, vec%p)

     return
   end subroutine state_vector_copy

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
                fcoeff(1) = 2.3e4_rp
                fcoeff(2) = 2.3e3_rp
                fcoeff(3) = -2.0e5_rp
                w%x(ix, iy, iz, iel) = math_ran_dst(ix, iy, iz, iel, xl, &
                     fcoeff) * 1.0e-08_rp

                if (xl(1) .lt. -12.0_rp) then
                   ! this is not needed I was just double checking my BCs
                   u%x(ix, iy, iz, iel) = 0.0_rp
                   v%x(ix, iy, iz, iel) = 0.0_rp
                end if
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

  subroutine state_vector_write(self, idx)
     class(state_vector), intent(inout) :: self
     integer :: idx

     call self%output%sample(real(idx, kind=rp))

     return
   end subroutine state_vector_write
 
  ! silly wrapper to ignore intent
  subroutine state_vector_init_wrapper(self)
     class(state_vector) :: self

     call self%init()

     return
   end subroutine state_vector_init_wrapper
   !------------------------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURES FOR THE EXPONENTIAL PROPAGATOR     -----
   !------------------------------------------------------------------------
 
   subroutine direct_solver(self, vec_in, vec_out)
     ! Linear Operator.
     class(neko_propagator), intent(inout)  :: self
     ! Input vector.
     class(abstract_vector_rdp) , intent(in)  :: vec_in
     ! Output vector.
     class(abstract_vector_rdp) , intent(out) :: vec_out
 
     select type(vec_in)
     type is(state_vector)
        select type(vec_out)
        type is(state_vector)
           ! Reset propagator.
           call reset_wrapper(self%linear_case)
           ! Get state vector.
           ! (again... the naming with "adjoint" isn't smart here)
           call field_copy(self%linear_case%fluid_adj%u_adj, vec_in%u)
           call field_copy(self%linear_case%fluid_adj%v_adj, vec_in%v)
           call field_copy(self%linear_case%fluid_adj%w_adj, vec_in%w)
           call field_copy(self%linear_case%fluid_adj%p_adj, vec_in%p)
           ! Integrate forward in time.
           call self%write_linear(0)
           call solve_wrapper(self%linear_case)
           ! call self%write_linear(1)
           ! There is a chance that vec_out isn't initialized!
           call init_wrapper(vec_out)
           ! Pass-back the state vector.
           call field_copy(vec_out%u, self%linear_case%fluid_adj%u_adj)
           call field_copy(vec_out%v, self%linear_case%fluid_adj%v_adj)
           call field_copy(vec_out%w, self%linear_case%fluid_adj%w_adj)
           call field_copy(vec_out%p, self%linear_case%fluid_adj%p_adj)
        end select
     end select
     return
   end subroutine direct_solver
 
   subroutine adjoint_solver(self, vec_in, vec_out)
     ! Linear Operator.
     class(neko_propagator), intent(inout)  :: self
     ! Input vector.
     class(abstract_vector_rdp) , intent(in)  :: vec_in
     ! Output vector.
     class(abstract_vector_rdp) , intent(out) :: vec_out
 
     select type(vec_in)
     type is(state_vector)
        select type(vec_out)
        type is(state_vector)
           ! Reset propagator.
           call reset_wrapper(self%adjoint_case)
           ! Get state vector.
           ! (again... the naming with "adjoint" isn't smart here)
           call field_copy(self%adjoint_case%fluid_adj%u_adj, vec_in%u)
           call field_copy(self%adjoint_case%fluid_adj%v_adj, vec_in%v)
           call field_copy(self%adjoint_case%fluid_adj%w_adj, vec_in%w)
           ! Integrate forward in time.
           !call self%write_adjoint(0)
           call solve_wrapper(self%linear_case)
           !call self%write_adjoint(0)
           ! There is a chance that vec_out isn't initialized!
           call init_wrapper(vec_out)
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
     type(sponge_source_term_t) :: linear_sponge, adjoint_sponge
     integer :: i, n
     real :: y_c

     ! initialize the "baseflow"
     call neko_init(self%neko_case)
     ! set the baseflow
     n = self%neko_case%fluid%u%size()
     do i = 1, n
        y_c = self%neko_case%fluid%dm_xh%y(i,1,1,1)
        self%neko_case%fluid%u%x(i,1,1,1) = -(y_c -1)*(y_c + 1)
        self%neko_case%fluid%v%x(i,1,1,1) = 0.0_rp
        self%neko_case%fluid%w%x(i,1,1,1) = 0.0_rp
     end do
     ! initialize the linear
     ! This is hacky, we'll work on a cleaner way soon :-)
     self%linear_case%if_adjoint = .false.
     call adjoint_init(self%linear_case, self%neko_case)
     ! initialize the adjoint
     call adjoint_init(self%adjoint_case, self%neko_case)

     ! Also hacky... but throw on a sponge
     call linear_sponge%init_from_components(self%linear_case%fluid_adj%f_adj_x, &
        self%linear_case%fluid_adj%f_adj_y, self%linear_case%fluid_adj%f_adj_z, &
        self%linear_case%fluid_adj%u_adj, self%linear_case%fluid_adj%v_adj, &
        self%linear_case%fluid_adj%w_adj, self%linear_case%fluid_adj%c_Xh)
     call self%linear_case%fluid_adj%source_term%add(linear_sponge)

     ! NOTE baseflow should be loaded via IC in .case file, but let's double
     ! check
          call self%output_primal%init(sp, 'checking_base', 4)
          call self%output_primal%fields%assign(2, self%neko_case%fluid%u)
          call self%output_primal%fields%assign(3, self%neko_case%fluid%v)
          call self%output_primal%fields%assign(4, self%neko_case%fluid%w)
          call self%output_primal%fields%assign(1, self%neko_case%fluid%p)
          call self%output_primal%sample(0.0_rp)
     ! Assign samplers for the forward and adjoint in case we want to look at
     ! them (debugging)
          call self%output_linear%init(sp, 'checking_linear', 4)
          call self%output_linear%fields%assign(2, self%linear_case%fluid_adj%u_adj)
          call self%output_linear%fields%assign(3, self%linear_case%fluid_adj%v_adj)
          call self%output_linear%fields%assign(4, self%linear_case%fluid_adj%w_adj)
          call self%output_linear%fields%assign(1, self%linear_case%fluid_adj%p_adj)
          call self%output_adjoint%init(sp, 'checking_adjoint', 4)
          call self%output_adjoint%fields%assign(2, self%adjoint_case%fluid_adj%u_adj)
          call self%output_adjoint%fields%assign(3, self%adjoint_case%fluid_adj%v_adj)
          call self%output_adjoint%fields%assign(4, self%adjoint_case%fluid_adj%w_adj)
          call self%output_adjoint%fields%assign(1, self%adjoint_case%fluid_adj%p_adj)

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
     ! Case
     class(adjoint_case_t)  :: self

     call solve_adjoint(self)

     return
   end subroutine solve_wrapper

   ! silly little wrapper to ignore intent.
   subroutine reset_wrapper(self)
     ! Linear Operator.
     class(adjoint_case_t)  :: self

     call adjoint_reset(self)

     return
   end subroutine reset_wrapper

   ! silly little wrapper to ignore intent.
   subroutine write_linear_wrapper(self, idx)
     ! Linear Operator.
     class(neko_propagator) :: self
     integer :: idx

     call self%output_linear%sample(real(idx, kind=rp))

     return
   end subroutine write_linear_wrapper

   ! silly little wrapper to ignore intent.
   subroutine write_adjoint_wrapper(self, idx)
     ! Linear Operator.
     class(neko_propagator) :: self
     integer :: idx

     call self%output_adjoint%sample(real(idx, kind=rp))

     return
   end subroutine write_adjoint_wrapper

    ! silly little wrapper to ignore intent.
   subroutine init_wrapper(self)
     class(state_vector) :: self

     call self%init()

     return
   end subroutine init_wrapper

 end module poiseuille
