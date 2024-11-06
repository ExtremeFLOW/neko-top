module user
  use neko
  use field_math, only: field_col3
  use point_zone, only: point_zone_t
  use point_zone_registry, only: neko_point_zone_registry
  use math, only: cfill_mask, glsc2
  implicit none

  ! Global user variables
  type(field_t) :: divergence

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%user_check => user_calc_quantities
    user%fluid_user_ic => user_ic
  end subroutine user_setup

  ! User-defined routine called at the end of every time step
  subroutine user_calc_quantities(t, tstep, u, v, w, p, coef, params)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params
    type(field_t), intent(inout) :: u, v, w, p
    type(field_t), pointer :: brinkman, mapped_brinkman
    type(field_t), pointer :: work
    integer :: temp_indices(2)

    integer :: ntot
    real(kind=rp) :: leakage, lift, drag

    if (neko_field_registry%field_exists("brinkman_indicator")) then
       brinkman => neko_field_registry%get_field("brinkman_indicator")
       mapped_brinkman => neko_field_registry%get_field("brinkman")

       ! Another good metric inspired by
       ! A. Ghasemi & A. Elham (2019)
       ! "FLOW TOPOLOGY OPTIMIZATION IN PERIODIC DOMAINS WITH APPLICATION TO
       ! MICRO HEAT EXCHANGER OPTIMIZATION"
       ! I think they are aluding to that: we only have one source term (brinkman)
       ! in the equation, so by integrating this source term we effectively are
       ! computing the lift and drag. (without the need to know the interface
       ! location)
       ! I think that's quite clever!

       ntot = u%dof%size()
       call neko_scratch_registry%request_field(work, temp_indices(1))
       call field_col3(work, mapped_brinkman, u)
       drag = glsc2(work%x, coef%B, ntot)
       call field_col3(work, mapped_brinkman, v)
       lift = glsc2(work%x, coef%B, ntot)
       call neko_scratch_registry%relinquish_field(temp_indices)

       ! calculate the leakage
       leakage = leak(brinkman%x, u%x, v%x, w%x, coef%b, ntot)
       if (pe_rank .eq. 0) then
          print *, 'Leakage = ', leakage, ',  ', t
          print *, 'Lift = ', lift, ',  ', t
          print *, 'Drag = ', drag, ',  ', t
       end if
    end if

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
    class(point_zone_t), pointer :: zone
    integer :: i, ntot

    ntot = u%dof%size()
    do i = 1, ntot
       u%x(i,1,1,1) = 1.0_rp
       w%x(i,1,1,1) = 0.0_rp

       ! just to break the symmetry and induce shedding quicker
       if (abs(u%dof%y(i,1,1,1)) .lt. 4.0_rp) then
          v%x(i,1,1,1) = 0.1_rp
       else
          v%x(i,1,1,1) = 0.0_rp
       end if

    end do
    p = 0.0_rp

    ! Get the cylinder mask
    if (neko_point_zone_registry%point_zone_exists("cylinder")) then
       zone => neko_point_zone_registry%get_point_zone("cylinder")

       ! Set the velocity to zero inside the cylinder
       call cfill_mask(u%x, 0.0_rp, u%size(), zone%mask, zone%size)
       call cfill_mask(v%x, 0.0_rp, v%size(), zone%mask, zone%size)
       call cfill_mask(w%x, 0.0_rp, w%size(), zone%mask, zone%size)
    end if
  end subroutine user_ic


end module user
