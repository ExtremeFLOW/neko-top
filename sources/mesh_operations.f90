! Copyright (c) 2024, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> Mesh operations and transformations
module mesh_operators
  use point, only: point_t
  use mesh, only: mesh_t
  use tri_mesh, only: tri_mesh_t
  use tet_mesh, only: tet_mesh_t
  use num_types, only: dp
  use utils, only: neko_error

  implicit none

  !> @brief Translate a given mesh.
  interface translate
     module procedure translate_array, translate_point, &
       translate_mesh, translate_tri_mesh, translate_tet_mesh
  end interface translate

  !> @brief Scale a given mesh.
  interface scale
     module procedure scale_array, scale_point, &
       scale_mesh, scale_tri_mesh, scale_tet_mesh, &
       scale_mesh_uniform, scale_tri_mesh_uniform, scale_tet_mesh_uniform
  end interface scale

  !> @brief Apply an affine transformation to a given mesh.
!   interface transform
!      module procedure array_transform, point_transform, &
!        mesh_transform, tri_mesh_transform, tet_mesh_transform
!   end interface transform

contains

  ! ========================================================================== !
  ! Actual Implementations

  !> @brief Translate a list of points by a given vector
  !> @param[inout] points The list of points to be translated
  !> @param[in] vec The translation vector
  subroutine translate_array(points, vec)
    real(kind=dp), dimension(:,:), intent(inout) :: points
    real(kind=dp), dimension(3), intent(in) :: vec
    integer :: i

    if (size(points, 1) /= 3) then
       call neko_error('translate_array: points must be a 3xN array')
    end if

    do i = 1, size(points, 2)
       points(:,i) = points(:,i) + vec
    end do

  end subroutine translate_array

  !> @brief Scale a list of points by a uniform factor.
  !> @param[inout] points The list of points to be scaled.
  !> @param[in] factor The uniform scaling factor for all dimensions.
  subroutine scale_array_uniform(points, factor)
    real(kind=dp), dimension(:,:), intent(inout) :: points
    real(kind=dp), intent(in) :: factor
    integer :: i

    if (size(points, 1) /= 3) then
       call neko_error('scale_array_uniform: points must be a 3xN array')
    end if

    do i = 1, size(points, 2)
       points(:,i) = points(:,i) * factor
    end do

  end subroutine scale_array_uniform

  !> @brief Scale a list of points by a list of factors.
  !> @param[inout] points The list of points to be scaled.
  !> @param[in] factor The scaling non-uniform factor.
  subroutine scale_array(points, factor)
    real(kind=dp), dimension(:,:), intent(inout) :: points
    real(kind=dp), dimension(3), intent(in) :: factor
    integer :: i

    if (size(points, 1) /= 3) then
       call neko_error('scale_array: points must be a 3xN array')
    end if

    do i = 1, size(points, 2)
       points(:,i) = points(:,i) * factor
    end do

  end subroutine scale_array

  !> @brief Apply an affine transformation to a list of points.
  !> @param[inout] points The list of points to be transformed.
  !> @param[in] matrix The transformation matrix.
  subroutine transform_array(points, matrix)
    real(kind=dp), dimension(:,:), intent(inout) :: points
    real(kind=dp), dimension(4,4), intent(in) :: matrix
    real(kind=dp), dimension(:,:), allocatable :: x

    if (size(points, 1) /= 3) then
       call neko_error('transform_array: points must be a 3xN array')
    end if

    allocate(x(4,size(points,2)))

    x(1:3,:) = points
    x(4,:) = 1.0_dp

    x = matmul(matrix, x)

    points = x(1:3, :)
  end subroutine transform_array

  !> @brief Build an affine transformation matrix for a list of transforms.
  !> @param[in, optional] translation The translation vector.
  !> @param[in, optional] scaling The scaling vector.
  !> @param[in, optional] rotation The rotation matrix.
  !> @param[in, optional] shear The shear matrix.
  !> @param[in, optional] reflect The boolean list of axis to reflect along.
  !> @returns matrix The resulting affine transformation matrix.
  pure function build_transform_matrix( &
    & translation, scaling, rotation, shear, reflect ) &
    & result(matrix)
    real(kind=dp), dimension(3), intent(in), optional :: translation, scaling
    real(kind=dp), dimension(3,3), intent(in), optional :: rotation, shear
    logical, dimension(3), intent(in), optional :: reflect
    real(kind=dp), dimension(4,4) :: matrix
    integer :: i, j

    matrix = reshape((/ 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                        0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
                        0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
                        0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp /), shape(matrix))

    if (present(translation)) then
       matrix(1:3,4) = translation
    end if

    if (present(scaling)) then
       matrix(1, 1) = scaling(1)
       matrix(2, 2) = scaling(2)
       matrix(3, 3) = scaling(3)
    end if

    if (present(rotation)) then
       matrix(1:3,1:3) = matrix(1:3,1:3) * rotation
    end if

    if (present(shear)) then
       matrix(1:3,1:3) = matrix(1:3,1:3) * shear
    end if

    if (present(reflect)) then
       do i = 1, 3
          if (reflect(i)) then
             do j = 1, 3
                matrix(j,i) = -matrix(j,i)
             end do
          end if
       end do
    end if

  end function build_transform_matrix

  ! ========================================================================== !
  ! Translation operations

  !> @brief Translate a list of points by a given vector
  !> @param[inout] points The list of points to be translated
  !> @param[in] vec The translation vector
  subroutine translate_point(points, vec)
    type(point_t), dimension(:), intent(inout) :: points
    real(kind=dp), dimension(3), intent(in) :: vec

    real(kind=dp), dimension(:,:), allocatable :: x
    integer :: i

    allocate(x(3, size(points)))

    do i = 1, size(points)
       x(:,i) = points(i)%x
    end do

    call translate_array(x, vec)

    do i = 1, size(points)
       points(i)%x = x(:,i)
    end do

  end subroutine translate_point

  !> @brief Translate a mesh by a given vector
  !> @param[in] object The mesh to be translated
  !> @param[in] vec The translation vector
  subroutine translate_mesh(object, vec)
    class(mesh_t), intent(inout) :: object
    real(kind=dp), dimension(3), intent(in) :: vec

    call translate_point(object%points, vec)

  end subroutine translate_mesh

  !> @brief Translate a triangular mesh by a given vector
  !> @param[in] object The mesh to be translated
  !> @param[in] vec The translation vector
  subroutine translate_tri_mesh(object, vec)
    class(tri_mesh_t), intent(inout) :: object
    real(kind=dp), dimension(3), intent(in) :: vec

    call translate_point(object%points, vec)

  end subroutine translate_tri_mesh

  !> @brief Translate a tetrahedral mesh by a given vector
  !> @param[in] object The mesh to be translated
  !> @param[in] vec The translation vector
  subroutine translate_tet_mesh(object, vec)
    class(tet_mesh_t), intent(inout) :: object
    real(kind=dp), dimension(3), intent(in) :: vec

    call translate_point(object%msh%points, vec)

  end subroutine translate_tet_mesh

  ! ========================================================================== !
  ! Scaling operations

  !> @brief Scale a point by a given factor.
  !> @param[inout] object The point to be scaled.
  !> @param[in] factor The scaling factor.
  subroutine scale_point(object, factor)
    type(point_t), intent(inout) :: object
    real(kind=dp), dimension(3), intent(in) :: factor

    object%x = object%x * factor

  end subroutine scale_point

  !> @brief Scale a mesh by a given factor.
  !> @param[in] object The mesh to be scaled.
  !> @param[in] factor The uniform scaling factor for all dimensions.
  subroutine scale_mesh_uniform(object, factor)
    class(mesh_t), intent(inout) :: object
    real(kind=dp), intent(in) :: factor
    integer :: i

    do i = 1, object%mpts
       object%points(i)%x = object%points(i)%x * factor
    end do

  end subroutine scale_mesh_uniform

  !> @brief Scale a triangular mesh by a given factor.
  !> @param[in] object The mesh to be scaled.
  !> @param[in] factor The uniform scaling factor for all dimensions.
  subroutine scale_tri_mesh_uniform(object, factor)
    class(tri_mesh_t), intent(inout) :: object
    real(kind=dp), intent(in) :: factor
    integer :: i

    do i = 1, object%mpts
       object%points(i)%x = object%points(i)%x * factor
    end do

  end subroutine scale_tri_mesh_uniform

  !> @brief Scale a tetrahedral mesh by a given factor.
  !> @param[in] object The mesh to be scaled.
  !> @param[in] factor The uniform scaling factor for all dimensions.
  subroutine scale_tet_mesh_uniform(object, factor)
    class(tet_mesh_t), intent(inout) :: object
    real(kind=dp), intent(in) :: factor
    integer :: i

    do i = 1, object%msh%mpts
       object%msh%points(i)%x = object%msh%points(i)%x * factor
    end do

  end subroutine scale_tet_mesh_uniform


  !> @brief Scale a mesh by a given factor.
  !> @param[in] object The mesh to be scaled.
  !> @param[in] factor The scaling factor.
  subroutine scale_mesh(object, factor)
    class(mesh_t), intent(inout) :: object
    real(kind=dp), dimension(3), intent(in) :: factor
    integer :: i

    do i = 1, object%mpts
       object%points(i)%x = object%points(i)%x * factor
    end do

  end subroutine scale_mesh

  !> @brief Scale a triangular mesh by a given factor.
  !> @param[in] object The mesh to be scaled.
  !> @param[in] factor The scaling factor.
  subroutine scale_tri_mesh(object, factor)
    class(tri_mesh_t), intent(inout) :: object
    real(kind=dp), dimension(3), intent(in) :: factor
    integer :: i

    do i = 1, object%mpts
       object%points(i)%x = object%points(i)%x * [factor]
    end do

  end subroutine scale_tri_mesh

  !> @brief Scale a tetrahedral mesh by a given factor.
  !> @param[in] object The mesh to be scaled.
  !> @param[in] factor The scaling factor.
  subroutine scale_tet_mesh(object, factor)
    class(tet_mesh_t), intent(inout) :: object
    real(kind=dp), dimension(3), intent(in) :: factor
    integer :: i

    do i = 1, object%msh%mpts
       object%msh%points(i)%x = object%msh%points(i)%x * factor
    end do

  end subroutine scale_tet_mesh

  ! ========================================================================== !
  ! Affine transformation operations

  !> @brief Apply an affine transformation to a list of points.
  !> @param[in] points The list of points to be transformed.
  !> @param[in] matrix The transformation matrix.
  subroutine transform_points(points, matrix)
    real(kind=dp), dimension(:,:), intent(inout) :: points
    real(kind=dp), dimension(4,4), intent(in) :: matrix
    real(kind=dp), dimension(:,:), allocatable :: x

    allocate(x(4,size(points,2)))

    x(1:3,:) = points
    x(4,:) = 1.0_dp

    x = matmul(matrix, x)

    points = x(1:3, :)
  end subroutine transform_points

  !> @brief Apply an affine transformation to a mesh.
  !> @param[in] object The mesh to be transformed.
  !> @param[in] matrix The transformation matrix.
  subroutine mesh_transform(object, matrix)
    class(mesh_t), intent(inout) :: object
    real(kind=dp), dimension(4,4), intent(in) :: matrix
    real(kind=dp), dimension(:,:), allocatable :: x
    integer :: i

    allocate(x(4,object%mpts))
    do i = 1, object%mpts
       x(1:3,i) = object%points(i)%x
       x(4,i) = 1.0_dp
    end do

    x = matmul(matrix, x)

    do i = 1, object%mpts
       object%points(i)%x = x(1:3,i)
    end do

    deallocate(x)
  end subroutine mesh_transform

end module mesh_operators
