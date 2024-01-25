module math_ext
  use num_types, only: rp
  implicit none

contains

  !> @brief Add a constant to a masked vector.
  !! \f$ a_i = a_i + c, for i in mask \f$
  subroutine cadd_mask(a, c, size, mask, mask_size)
    real(kind=rp), dimension(size), intent(inout) :: a
    real(kind=rp), intent(in) :: c
    integer, intent(in) :: size
    integer, dimension(mask_size), intent(in) :: mask
    integer, intent(in) :: mask_size
    integer :: i

    do i = 1, mask_size
       a(mask(i)) = a(mask(i)) + c
    end do

  end subroutine cadd_mask

  !> @brief Invert a masked vector.
  !! \f$ a_i = 1/a_i, for i in mask \f$
  subroutine invcol1_mask(a, size, mask, mask_size)
    real(kind=rp), dimension(size), intent(inout) :: a
    integer, intent(in) :: size
    integer, dimension(mask_size), intent(in) :: mask
    integer, intent(in) :: mask_size
    integer :: i

    do i = 1, mask_size
       a(mask(i)) = 1.0_rp / a(mask(i))
    end do

  end subroutine invcol1_mask

  !> @brief Multiply a masked vector by a constant.
  !! \f$ a_i = c * a_i, for i in mask \f$
  subroutine cmult_mask(a, c, size, mask, mask_size)
    real(kind=rp), dimension(size), intent(inout) :: a
    real(kind=rp), intent(in) :: c
    integer, intent(in) :: size
    integer, dimension(mask_size), intent(in) :: mask
    integer, intent(in) :: mask_size
    integer :: i

    do i = 1, mask_size
       a(mask(i)) = c * a(mask(i))
    end do

  end subroutine cmult_mask

  !> @brief Multiply 2 masked vectors. Save the result in a new vector.
  !! \f$ a_i = b_i * c_i, for i in mask \f$
  subroutine col2_mask(a, b, size, mask, mask_size)
    real(kind=rp), dimension(size), intent(inout) :: a
    real(kind=rp), dimension(size), intent(in) :: b
    integer, intent(in) :: size
    integer, dimension(mask_size), intent(in) :: mask
    integer, intent(in) :: mask_size
    integer :: i

    do i = 1, mask_size
       a(mask(i)) = a(mask(i)) * b(mask(i))
    end do

  end subroutine col2_mask

  !> @brief Multiply 2 masked vectors. Save the result in a new vector.
  !! \f$ a_i = b_i * c_i, for i in mask \f$
  subroutine col3_mask(a, b, c, size, mask, mask_size)
    real(kind=rp), dimension(size), intent(inout) :: a
    real(kind=rp), dimension(size), intent(in) :: b, c
    integer, intent(in) :: size
    integer, dimension(mask_size), intent(in) :: mask
    integer, intent(in) :: mask_size
    integer :: i

    do i = 1, mask_size
       a(mask(i)) = b(mask(i)) * c(mask(i))
    end do

  end subroutine col3_mask

end module math_ext
