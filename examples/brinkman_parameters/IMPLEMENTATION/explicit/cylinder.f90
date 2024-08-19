module user
  use neko
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
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(field_t), pointer :: brinkman
    type(field_t), pointer :: divergence

    integer ntot
    real(kind=rp) :: leakage
    real(kind=rp) :: div_tot

    brinkman => neko_field_registry%get_field("brinkman_indicator")
    divergence => neko_field_registry%get_field_by_name('div_fam')

    ntot = u%dof%size()


   ! calculate the leakage
    leakage = leak(brinkman%x,u%x,v%x,w%x,coef%b,ntot)
   ! calculate the divergence
   call div(divergence%x, u%x, v%x, w%x, coef)
   div_tot = glsc3(divergence%x,divergence%x,coef%b,ntot)
        if (pe_rank .eq. 0) then
        	 print *, 'Leakage = ', leakage, ',  ', t
        	 print *, 'Divergence = ', sqrt(div_tot), ', ', t
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
        tmp = tmp + brink(i) *sqrt( u(i)**2 + v(i)**2 +  w(i)**2) * B(i)
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
