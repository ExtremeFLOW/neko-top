! ============================================================================ !
! Supporting functions
! ============================================================================ !
module develop

  integer, public, parameter :: facet_type_interior = 0
  integer, public, parameter :: facet_type_outlet = 1
  integer, public, parameter :: facet_type_inlet = 2
  integer, public, parameter :: facet_type_wall = 3
  integer, public, parameter :: facet_type_sympln = 4

contains

  !> Compute cross product of two vectors
  !> @param[in] a First vector
  !> @param[in] b Second vector
  !> @return Cross product \f$ a \times b \f$
  pure function cross(a, b) result(c)
    use num_types
    implicit none

    real(kind=rp), dimension(3), intent(in) :: a
    real(kind=rp), dimension(3), intent(in) :: b
    real(kind=rp), dimension(3) :: c

    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)

  end function cross

  !> Compute dot product of two vectors.
  !> @param[in] a First vector.
  !> @param[in] b Second vector.
  !> @return dot product \f$ a \cdot b \f$.
  pure function dot(a, b) result(d)
    use num_types
    implicit none

    real(kind=rp), dimension(3), intent(in) :: a(3)
    real(kind=rp), dimension(3), intent(in) :: b(3)
    real(kind=rp) :: d

    d = a(1) * b(1) + a(2) * b(2) + a(3) * b(3)

  end function dot

  !> Compute area of a polygon.
  !> @param[in] nv Number of vertices.
  !> @param[in] vertices Vertices.
  !> @return Area of polygon.
  !> @note The vertices must be ordered in order either clockwise or
  !> counter-clockwise around the polygon.
  pure function area(nv, vertices) result(a)
    use num_types
    implicit none
    integer, intent(in) :: nv
    real(kind=rp), dimension(3, nv), intent(in) :: vertices
    real(kind=rp) :: a

    ! Internal variables
    real(kind=rp), dimension(3) :: v1
    real(kind=rp) :: cp(3)
    integer :: v

    ! Initialize the origin to the first vertex, which reduce the
    ! computational cost by allowing us to skip the first and last cross
    ! products.
    v1 = vertices(:, 1)
    cp = 0.0
    do v = 2, nv - 1
       cp = cp + 0.5 * cross((vertices(:, v) - v1), (vertices(:, v + 1) - v1))
    end do

    a = sqrt(dot(cp, cp))

  end function area

  !> Construct a list of facets.
  !> @param[in] C Case data structure.
  !> @param[out] out List of facets, size (2, N_facets), first row is element.
  !> ID, second row is facet ID.
  !> @param[in] facet_type Type of facet, optional,
  !> (0 = interior, 1 = outlet [default]).
  subroutine get_facets(C, out, facet_type)
    use case
    implicit none

    type(case_t), intent(in) :: C
    integer, dimension(:, :), allocatable, intent(out) :: out
    integer, intent(in), optional :: facet_type

    ! Local variables
    integer :: e, f, facet_id
    integer :: N_elem, N_facets
    integer :: facet_type_

    if (present(facet_type)) then
       facet_type_ = facet_type
    else
       facet_type_ = 1
    end if

    ! Total number of elements
    N_elem = C%msh%nelv

    ! Count number of output facets (Facet_type is 1)
    N_facets = 0
    do e = 1, N_elem
       do f = 1, 6
          if (C%msh%facet_type(f, e) .eq. facet_type_) then
             N_facets = N_facets + 1
          end if
       end do
    end do

    ! Allocate the outlet_elements array
    allocate (out(2, N_facets))

    ! Loop across the elements and store the outlet elements
    facet_id = 1
    do e = 1, N_elem
       do f = 1, 6
          if (C%msh%facet_type(f, e) .eq. facet_type_) then
             out(1, facet_id) = e
             out(2, facet_id) = f
             facet_id = facet_id + 1
          end if
       end do
    end do

  end subroutine get_facets

  !> Construct a list of facets.
  !> @param[in] C Case data structure.
  !> @param[out] out List of outlet facets, size (2, N_facets),
  !> first row is element ID, second row is facet ID.
  subroutine get_facets_outlet(C, out)
    use case, only: case_t
    implicit none

    type(case_t), intent(in) :: C
    integer, dimension(:, :), allocatable, intent(out) :: out

    call get_facets(C, out, facet_type_outlet)
  end subroutine get_facets_outlet

  !> Compute the weighted sum of an array.
  !> @param[in] a Array.
  !> @param[in] w Weights, optional.
  !> @return Weighted sum of array.
  !> @details This function computes the weighted sum of an array. The
  !> weights must be the same size as the array. If no weights are supplied, the
  !> unweighted sum is computed.
  pure function sum_weighted(a, w) result(s)
    use num_types
    implicit none

    real(kind=rp), dimension(:), intent(in) :: a
    real(kind=rp), dimension(:), intent(in), optional :: w
    real(kind=rp) :: s

    if (present(w)) then
       s = sum(a * w)
    else
       s = sum(a)
    end if

  end function sum_weighted

  !> Compute the weighted average of an array.
  !> @param[in] a Array.
  !> @param[in] w Weights, optional.
  !> @return Weighted average of array.
  !> @details This function computes the weighted average of an array. The
  !> weights must be the same size as the array.
  pure function average_weighted(a, w) result(m)
    use num_types
    implicit none

    real(kind=rp), dimension(:), intent(in) :: a
    real(kind=rp), dimension(:), intent(in), optional :: w
    real(kind=rp) :: m

    if (present(w)) then
       m = sum_weighted(a, w) / sum(w)
    else
       m = sum(a) / size(a)
    end if

  end function average_weighted

  !> Construct a list of nodes for a facet.
  !> @param[in] C Case data structure.
  !> @param[in] element_id Element ID.
  !> @param[in] facet_id Facet ID.
  !> @param[out] nodes Node list.
  subroutine get_facet_nodes(C, element_id, facet_id, nodes)
    use case
    use num_types
    use tuple
    use hex
    implicit none

    type(case_t), intent(in) :: C
    integer, intent(in) :: element_id
    integer, intent(in) :: facet_id
    real(kind=rp), dimension(:, :), allocatable, intent(out) :: nodes

    ! Local variables
    integer :: n, v
    integer :: N_nodes
    type(tuple4_i4_t) :: t_hex

    ! Determine the number of nodes in the facet
    N_nodes = 0
    select type (ele => C%msh%elements(element_id)%e)
      type is (hex_t)
       N_nodes = 4
    end select

    if (N_nodes .eq. 0) then
       call neko_error("Error: the facet shape is not supported.")
    end if

    ! Allocate the nodes array
    if (allocated(nodes) .and. size(nodes, 2) .eq. N_nodes) then
       nodes = 0.0
    else if (.not. allocated(nodes)) then
       allocate (nodes(3, N_nodes))
    else
       call neko_error("Error: the nodes array is not the correct size.")
    end if

    ! Get the nodes
    select type (ele => C%msh%elements(element_id)%e)
      type is (hex_t)
       call ele%facet_order(t_hex, facet_id)
       do n = 1, N_nodes
          v = t_hex%x(n)
          nodes(:, n) = C%msh%points(v)%x
       end do
    end select

  end subroutine get_facet_nodes

  subroutine estimate_temperature(neko_case)
    use case, only: case_t
    use global_interpolation, only: global_interpolation_t
    use json_utils, only: json_get_or_default
    use logger, only: neko_log, LOG_SIZE
    use num_types, only: rp
    use tuple, only: tuple4_i4_t
    implicit none

    type(case_t), intent(inout) :: neko_case
    type(global_interpolation_t) :: interpolator

    real(kind=rp) :: target_temperature
    real(kind=rp) :: temperature_mean
    real(kind=rp), dimension(:), allocatable :: temperature_local

    integer, allocatable :: facet_list(:, :)

    real(kind=rp), dimension(:, :), allocatable :: facet_nodes
    real(kind=rp), dimension(:, :), allocatable :: facet_centers
    real(kind=rp), dimension(:), allocatable :: facet_area

    integer :: element_id, facet_id
    integer :: N_facets
    integer :: f, n

    character(len=LOG_SIZE) :: log_buf

    ! Read the case file for options
    call json_get_or_default(neko_case%params, 'topopt.target_temperature', &
                             target_temperature, 0.5_rp)

    ! Initialize the global interpolation
    call interpolator%init(neko_case%scalar%dm_xh)

    ! Get the list of outlet facets
    call get_facets(neko_case, facet_list)

    ! Allocate the facet_list array
    N_facets = size(facet_list, 2)
    allocate (facet_centers(3, N_facets))
    allocate (facet_area(N_facets))
    allocate (temperature_local(N_facets))

    ! Loop across the outlet elements and print the scalar field
    do f = 1, N_facets
       element_id = facet_list(1, f)
       facet_id = facet_list(2, f)

       ! Get the nodes for the facet
       call get_facet_nodes(neko_case, element_id, facet_id, facet_nodes)

       ! Compute the center of the facet
       facet_centers(:, f) = 0.0
       do n = 1, size(facet_nodes, 2)
          facet_centers(:, f) = facet_centers(:, f) + facet_nodes(:, n)
       end do
       facet_centers(:, f) = facet_centers(:, f) / size(facet_nodes, 2)

       ! Compute the area of the facet
       facet_area(f) = area(size(facet_nodes, 2), facet_nodes)
    end do

    ! Find the outlet temperature at the supplied list of points
    call interpolator%find_points_xyz(facet_centers, N_facets)
    call interpolator%evaluate(temperature_local, neko_case%scalar%s%x)

    temperature_mean = average_weighted( &
                                         temperature_local - target_temperature, &
                                         facet_area)

    write (log_buf, '(a,f15.7)') &
         "Outlet area-weighted average temperature deviation: ", &
         temperature_mean
    call neko_log%message(log_buf)

    call interpolator%free()

  end subroutine estimate_temperature

end module develop
