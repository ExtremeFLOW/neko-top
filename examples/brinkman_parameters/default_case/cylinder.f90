module user
  use user_intf, only: user_t
  use field, only: field_t
  use field_list, only: field_list_t
  use registry, only: neko_registry
  use user_access_singleton, only: neko_user_access
  use time_state, only: time_state_t
  use comm, only: neko_comm, pe_rank, mpi_real_precision
  use mpi_f08, only: mpi_allreduce, mpi_sum
  use num_types, only: rp
  use coefs, only: coef_t
  use device, only: HOST_TO_DEVICE, DEVICE_TO_HOST
  use point_zone, only: point_zone_t
  use point_zone_registry, only: neko_point_zone_registry
  use math, only: cfill_mask, glsc3
  implicit none

  ! Global user variables
  type(field_t) :: divergence

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%compute => user_calc_quantities
    user%initial_conditions => user_ic
  end subroutine user_setup

  ! User-defined routine called at the end of every time step
  subroutine user_calc_quantities(time)
    type(time_state_t), intent(in) :: time
    type(coef_t), pointer :: coef
    type(field_t), pointer :: u, v, w
    type(field_t), pointer :: brinkman, mapped_brinkman

    integer :: ntot
    real(kind=rp) :: leakage, lift, drag

    if (neko_registry%field_exists("brinkman_indicator")) then
       coef => neko_user_access%case%fluid%c_Xh
       u => neko_registry%get_field("u")
       v => neko_registry%get_field("v")
       w => neko_registry%get_field("w")
       brinkman => neko_registry%get_field("brinkman_indicator")
       mapped_brinkman => neko_registry%get_field("brinkman_permeability")

       ! The quantities below are computed on the host
       call u%copy_from(DEVICE_TO_HOST, sync = .false.)
       call v%copy_from(DEVICE_TO_HOST, sync = .false.)
       call w%copy_from(DEVICE_TO_HOST, sync = .false.)
       call brinkman%copy_from(DEVICE_TO_HOST, sync = .false.)
       call mapped_brinkman%copy_from(DEVICE_TO_HOST, sync = .true.)

       ! Another good metric inspired by
       ! A. Ghasemi & A. Elham (2019)
       ! "FLOW TOPOLOGY OPTIMIZATION IN PERIODIC DOMAINS WITH APPLICATION TO
       ! MICRO HEAT EXCHANGER OPTIMIZATION"
       ! I think they are aluding to that: we only have one source term
       ! (brinkman)
       ! in the equation, so by integrating this source term we effectively are
       ! computing the lift and drag. (without the need to know the interface
       ! location)
       ! I think that's quite clever!

       ntot = u%dof%size()
       drag = glsc3(mapped_brinkman%x, u%x, coef%B, ntot)
       lift = glsc3(mapped_brinkman%x, v%x, coef%B, ntot)

       ! calculate the leakage
       leakage = leak(brinkman%x, u%x, v%x, w%x, coef%B, ntot)
       if (pe_rank .eq. 0) then
          print *, 'Leakage = ', leakage, ',  ', time%t
          print *, 'Lift = ', lift, ',  ', time%t
          print *, 'Drag = ', drag, ',  ', time%t
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
       tmp = tmp + brink(i) * sqrt(u(i)**2 + v(i)**2 + w(i)**2) * B(i)
    end do

    call mpi_allreduce(tmp, leak, 1, &
         mpi_real_precision, mpi_sum, neko_comm, ierr)

  end function leak

  ! User-defined initial condition
  subroutine user_ic(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields
    type(field_t), pointer :: u, v, w, p
    class(point_zone_t), pointer :: zone
    integer :: i, ntot

    if (trim(scheme_name) .ne. 'fluid') return

    u => fields%get("u")
    v => fields%get("v")
    w => fields%get("w")
    p => fields%get("p")

    ntot = u%dof%size()
    do i = 1, ntot
       u%x(i,1,1,1) = sqrt(1.0_rp - 0.1_rp**2)
       w%x(i,1,1,1) = 0.0_rp

       ! just to break the symmetry and induce shedding quicker
       if (abs(u%dof%y%x(i,1,1,1)) .lt. 4.0_rp) then
          v%x(i,1,1,1) = 0.1_rp
       else
          v%x(i,1,1,1) = 0.0_rp
       end if

    end do
    p%x = 0.0_rp

    ! Get the cylinder mask
    if (neko_point_zone_registry%point_zone_exists("cylinder")) then
       zone => neko_point_zone_registry%get_point_zone("cylinder")

       ! Set the velocity to zero inside the cylinder
       call cfill_mask(u%x, 0.0_rp, u%size(), zone%mask%get(), zone%size)
       call cfill_mask(v%x, 0.0_rp, v%size(), zone%mask%get(), zone%size)
       call cfill_mask(w%x, 0.0_rp, w%size(), zone%mask%get(), zone%size)
    end if

    ! Neko does not copy user initial conditions to the device
    call u%copy_from(HOST_TO_DEVICE, sync = .false.)
    call v%copy_from(HOST_TO_DEVICE, sync = .false.)
    call w%copy_from(HOST_TO_DEVICE, sync = .false.)
    call p%copy_from(HOST_TO_DEVICE, sync = .true.)
  end subroutine user_ic


end module user
