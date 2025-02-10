!> Defines MMA factory for the mma_t
submodule (mma) mma_fctry
  use mma, only : mma_t
  use mma_cpu, only : mma_cpu_t
  use mma_device, only : mma_device_t
  use neko_config, only : NEKO_BCKND_DEVICE, NEKO_BCKND_SX
  use comm, only : pe_rank
  
  implicit none

contains

  !> Factory for the mma_t
  !! @details Only selects the compute backend.
  !! @param object The object to be allocated by the factory.
  module subroutine mma_factory(object,backnd)
    class(mma_t), allocatable, intent(inout) :: object
    character(len=:), allocatable, intent(in) :: backnd

    if (allocated(object)) then
       deallocate(object)
    end if

    if (backnd .eq. "cpu") then
      allocate(mma_cpu_t::object)
      if (pe_rank .eq. 0) then
         print *,"mma_cpu_t is set as the mma_t drived type in mma_factory!"
      end if
    else if (backnd .eq. "cuda") then
      allocate(mma_device_t::object)
      if (pe_rank .eq. 0) then
         print *,"mma_device_t is set as the mma_t drived type in mma_factory!"
      end if
    else
      allocate(mma_cpu_t::object)
      if (pe_rank .eq. 0) then
         print *,"mma_cpu_t is set as the mma_t drived type in mma_factory!"
      end if
    end if

  end subroutine mma_factory

end submodule mma_fctry
