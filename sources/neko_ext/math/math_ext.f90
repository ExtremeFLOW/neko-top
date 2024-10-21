module math_ext
  use num_types, only: rp, xp
  use comm
  implicit none

contains

  !> @brief Add a constant to a masked vector.
  !! \f$ a_i = a_i + c, for i in mask \f$
  subroutine cadd_mask(a, c, size, mask, mask_size)
    integer, intent(in) :: size, mask_size
    real(kind=rp), dimension(size), intent(inout) :: a
    real(kind=rp), intent(in) :: c
    integer, dimension(mask_size), intent(in) :: mask
    integer :: i

    do i = 1, mask_size
       a(mask(i)) = a(mask(i)) + c
    end do

  end subroutine cadd_mask

  !> @brief Invert a masked vector.
  !! \f$ a_i = 1/a_i, for i in mask \f$
  subroutine invcol1_mask(a, size, mask, mask_size)
    integer, intent(in) :: size, mask_size
    real(kind=rp), dimension(size), intent(inout) :: a
    integer, dimension(mask_size), intent(in) :: mask
    integer :: i

    do i = 1, mask_size
       a(mask(i)) = 1.0_rp / a(mask(i))
    end do

  end subroutine invcol1_mask

  !> @brief Multiply a masked vector by a constant.
  !! \f$ a_i = c * a_i, for i in mask \f$
  subroutine cmult_mask(a, c, size, mask, mask_size)
    integer, intent(in) :: size, mask_size
    real(kind=rp), dimension(size), intent(inout) :: a
    real(kind=rp), intent(in) :: c
    integer, dimension(mask_size), intent(in) :: mask
    integer :: i

    do i = 1, mask_size
       a(mask(i)) = c * a(mask(i))
    end do

  end subroutine cmult_mask

  !> @brief Multiply 2 masked vectors. Save the result in a new vector.
  !! \f$ a_i = b_i * c_i, for i in mask \f$
  subroutine col2_mask(a, b, size, mask, mask_size)
    integer, intent(in) :: size, mask_size
    real(kind=rp), dimension(size), intent(inout) :: a
    real(kind=rp), dimension(size), intent(in) :: b
    integer, dimension(mask_size), intent(in) :: mask
    integer :: i

    do i = 1, mask_size
       a(mask(i)) = a(mask(i)) * b(mask(i))
    end do

  end subroutine col2_mask

  !> @brief Multiply 2 masked vectors. Save the result in a new vector.
  !! \f$ a_i = b_i * c_i, for i in mask \f$
  subroutine col3_mask(a, b, c, size, mask, mask_size)
    integer, intent(in) :: size, mask_size
    real(kind=rp), dimension(size), intent(inout) :: a
    real(kind=rp), dimension(size), intent(in) :: b, c
    integer, dimension(mask_size), intent(in) :: mask
    integer :: i

    do i = 1, mask_size
       a(mask(i)) = b(mask(i)) * c(mask(i))
    end do

  end subroutine col3_mask

  !> @brief Subtract 2 masked vectors. Save the result in a new vector.
  !! \f$ a_i = b_i - c_i, for i in mask \f$
  subroutine sub3_mask(a, b, c, size, mask, mask_size)
    integer, intent(in) :: size, mask_size
    real(kind=rp), dimension(size), intent(inout) :: a
    real(kind=rp), dimension(size), intent(in) :: b, c
    integer, dimension(mask_size), intent(in) :: mask
    integer :: i

    do i = 1, mask_size
       a(mask(i)) = b(mask(i)) - c(mask(i))
    end do

  end subroutine sub3_mask

  !> @brief Weighted inner product 
  !! \f$ a^T b \f$ for indices in the mask
  function glsc2_mask(a, b, size, mask, mask_size)
    integer, intent(in) :: size, mask_size
    real(kind=rp), dimension(size), intent(in) :: a
    real(kind=rp), dimension(size), intent(in) :: b
    integer, dimension(mask_size), intent(in) :: mask
    real(kind=rp) :: glsc2_mask
    real(kind=xp) :: tmp
    integer :: i, ierr

    tmp = 0.0_xp
    do i = 1, mask_size
       tmp = tmp + a(mask(i)) * b(mask(i))
    end do

    call MPI_Allreduce(MPI_IN_PLACE, tmp, 1, &
         MPI_EXTRA_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    glsc2_mask = tmp
  end function glsc2_mask
end module math_ext
