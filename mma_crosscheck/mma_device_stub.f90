! The device submodule is not built in this harness; provide the two module
! procedures so the type-bound bindings in mma.f90 link.
submodule (mma) mma_device_stub
contains
  module subroutine mma_update_device(this, iter, x, df0dx, fval, dfdx)
    class(mma_t), intent(inout) :: this
    integer, intent(in) :: iter
    type(c_ptr), intent(inout) :: x
    type(c_ptr), intent(in) :: df0dx, fval, dfdx
    call neko_error('device backend not built in harness')
  end subroutine mma_update_device
  module subroutine mma_KKT_device(this, x, df0dx, fval, dfdx)
    class(mma_t), intent(inout) :: this
    type(c_ptr), intent(in) :: x, df0dx, fval, dfdx
    call neko_error('device backend not built in harness')
  end subroutine mma_KKT_device
end submodule mma_device_stub
