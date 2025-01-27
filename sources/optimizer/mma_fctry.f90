!> Defines MMA factory for the mma_t
submodule (mma) mma_fctry
  use mma, only : mma_t
  use mma_cpu, only : mma_cpu_t
  use mma_device, only : mma_device_t
  use neko_config, only : NEKO_BCKND_DEVICE, NEKO_BCKND_SX
  ! use pnpn_res_device, only : pnpn_prs_res_device_t, pnpn_vel_res_device_t
  ! use pnpn_res_cpu, only : pnpn_prs_res_cpu_t, pnpn_vel_res_cpu_t
  ! use pnpn_res_sx, only : pnpn_prs_res_sx_t, pnpn_vel_res_sx_t
  implicit none

contains

  !> Factory for the mma_t
  !! @details Only selects the compute backend.
  !! @param object The object to be allocated by the factory.
  module subroutine mma_factory(object,backnd)
    class(mma_t), allocatable, intent(inout) :: object
    character(len=:), allocatable :: backnd

    if (allocated(object)) then
       deallocate(object)
    end if

    if (backnd .eq. "cpu") then
      allocate(mma_cpu_t::object)
      print *,"mma_cpu_t is set as the mma_t drived type in mma_factory!"
    else if (backnd .eq. "cuda") then
      allocate(mma_device_t::object)
      print *,"mma_device_t is set as the mma_t drived type in mma_factory!"
    else
      allocate(mma_cpu_t::object)
      print *,"mma_cpu_t is set as the mma_t drived type in mma_factory!"
    end if

    ! if (NEKO_BCKND_SX .eq. 1) then
    ! else if (NEKO_BCKND_DEVICE .eq. 1) then
    !   print *,"mma_device_t is set as the mma_t drived type in mma_factory!"
    ! else
    !   allocate(mma_cpu_t::object)
    !   print *,"mma_cpu_t is set as the mma_t drived type in mma_factory!"
    ! end if

    ! if (NEKO_BCKND_SX .eq. 1) then
    !    allocate(pnpn_prs_res_sx_t::object)
    ! else if (NEKO_BCKND_DEVICE .eq. 1) then
    !    allocate(pnpn_prs_res_device_t::object)
    ! else
    !    allocate(pnpn_prs_res_cpu_t::object)
    ! end if

  end subroutine mma_factory

end submodule mma_fctry
