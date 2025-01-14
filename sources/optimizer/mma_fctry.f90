!> Defines MMA factory for the mma_t
submodule (mma) mma_fctry
  use mma, only : mma_t
  use mma_cpu, only : mma_cpu_t
  use neko_config, only : NEKO_BCKND_DEVICE, NEKO_BCKND_SX
  ! use pnpn_res_device, only : pnpn_prs_res_device_t, pnpn_vel_res_device_t
  ! use pnpn_res_cpu, only : pnpn_prs_res_cpu_t, pnpn_vel_res_cpu_t
  ! use pnpn_res_sx, only : pnpn_prs_res_sx_t, pnpn_vel_res_sx_t
  implicit none

contains

  !> Factory for the mma_t
  !! @details Only selects the compute backend.
  !! @param object The object to be allocated by the factory.
  module subroutine mma_factory(object)
    class(mma_t), allocatable, intent(inout) :: object

    if (allocated(object)) then
       deallocate(object)
    end if

    if (NEKO_BCKND_SX .eq. 1) then
    else if (NEKO_BCKND_DEVICE .eq. 1) then
    else
      allocate(mma_cpu_t::object)
    end if

    ! if (NEKO_BCKND_SX .eq. 1) then
    !    allocate(pnpn_prs_res_sx_t::object)
    ! else if (NEKO_BCKND_DEVICE .eq. 1) then
    !    allocate(pnpn_prs_res_device_t::object)
    ! else
    !    allocate(pnpn_prs_res_cpu_t::object)
    ! end if

  end subroutine mma_factory

end submodule mma_fctry
