module user
  use neko
  use field_math, only: field_col3
  implicit none

  ! Global user variables
  type(field_t) :: divergence

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%user_check => user_calc_quantities
    user%fluid_user_ic => user_ic
    user%user_init_modules => initialize
  end subroutine user_setup

  ! Initialize user variables or external objects
  subroutine initialize(t, u, v, w, p, coef, params)
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params
  
    call neko_field_registry%add_field(u%dof, "dudx", .true.)
    call neko_field_registry%add_field(u%dof, "dudy", .true.)
    call neko_field_registry%add_field(u%dof, "dvdx", .true.)
    call neko_field_registry%add_field(u%dof, "dvdy", .true.)
  
  end subroutine initialize
  



  ! User-defined routine called at the end of every time step
  subroutine user_calc_quantities(t, tstep, u, v, w, p, coef, params)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(field_t), pointer :: brinkman
    type(field_t), pointer :: dudx, dudy, dudz
    type(field_t), pointer :: dvdx, dvdy, dvdz
    type(field_t), pointer :: work
    integer :: temp_indices(1)

    integer ntot
    real(kind=rp) :: leakage, lift, drag
    real(kind=rp) :: div_tot

    ! This may be really inefficent... but we're going to need to probe the 
    ! gradient if we want the viscous contribution to the drag/lift.
    dudx => neko_field_registry%get_field("dudx")
    dudy => neko_field_registry%get_field("dvdx")
    dvdx => neko_field_registry%get_field("dudy")
    dvdy => neko_field_registry%get_field("dvdy")

    call dudxyz(dudx%x, u%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call dudxyz(dudy%x, u%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call dudxyz(dvdx%x, v%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call dudxyz(dvdy%x, v%x, coef%drdy, coef%dsdy, coef%dtdy, coef)

    ! We need the Brinkman term in the registry 
    call neko_field_registry%add_field(u%dof, "brinkman_indicator", .true.)
    brinkman => neko_field_registry%get_field("brinkman_indicator")

    ntot = u%dof%size()

    ! Another good metric inspired by
    ! A. Ghasemi & A. Elham (2019)
    ! "FLOW TOPOLOGY OPTIMIZATION IN PERIODIC DOMAINS WITH APPLICATION TO 
    ! MICRO HEAT EXCHANGER OPTIMIZATION"

    ! I think they are aluding to that: we only have one source term (brinkman)
    ! in the equation, so by integrating this source term we effectively are
    ! computing the lift and drag. (without the need to know the interface
    ! location)
    ! I think that's quite clever!

    call neko_scratch_registry%request_field(work, temp_indices(1))
    call field_col3(work, brinkman, u)
    drag = glsc2(work%x, coef%B, ntot)
    call field_col3(work, brinkman, v)
    lift = glsc2(work%x, coef%B, ntot)





    ! calculate the leakage
    leakage = leak(brinkman%x,u%x,v%x,w%x,coef%b,ntot)
    if (pe_rank .eq. 0) then
       print *, 'Leakage = ', leakage, ',  ', t
       print *, 'Lift = ', lift, ',  ', t
       print *, 'Drag = ', drag, ',  ', t
    endif








  end subroutine user_calc_quantities

  function leak(brink, u, v, w, B, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: brink
    real(kind=rp), dimension(n), intent(in) :: u
    real(kind=rp), dimension(n), intent(in) :: v
    real(kind=rp), dimension(n), intent(in) :: w
    real(kind=rp), dimension(n), intent(in) :: B
    real(kind=rp) :: leak, tmp
    integer :: i, ierr

    tmp = 0.0_rp
    do i = 1, n
       tmp = tmp + brink(i) *sqrt( u(i)**2 + v(i)**2 + w(i)**2) * B(i)
    end do

    call mpi_allreduce(tmp, leak, 1, &
         mpi_real_precision, mpi_sum, neko_comm, ierr)

  end function leak


  ! User-defined initial condition
  subroutine user_ic(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: params
    integer :: i, ntot

    ntot = u%dof%size()
    do i = 1, ntot
       u%x(i,1,1,1) = 1.0_rp
       w%x(i,1,1,1) = 1.0_rp

       ! just to break the symmetry and induce shedding quicker
       if(abs(u%dof%y(i,1,1,1)).lt.4.0_rp) then
          v%x(i,1,1,1) = 0.1_rp
       else
          v%x(i,1,1,1) = 0.0_rp
       endif

    end do
    p = 0._rp
  end subroutine user_ic


end module user
