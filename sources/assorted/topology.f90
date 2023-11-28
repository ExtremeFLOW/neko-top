module topology_module
   use neko

   implicit none
   private
   public :: topopt_permeability_force

   real(kind=rp), dimension(:), allocatable :: resistance
   type(c_ptr) :: resistance_d = c_null_ptr

   type, public :: topology_t
      !> @brief The lower left front corner
      real(kind=rp), dimension(3) :: llf = 0.0_rp
      !> @brief The upper right corner
      real(kind=rp), dimension(3) :: urb = 0.0_rp
      !> @brief The number of cells in each direction
      integer, dimension(3) :: n = 0
      !> @brief array describing the topology
      real(kind=rp), dimension(:, :, :), allocatable :: design
      !> @brief Pointer to the design domain
      class(point_zone_t), pointer :: design_domain => null()

      !Todo: We might need to include something to suport GPU here, such as a
      !      pointer to the GPU array of the design

      ! Limits of the permeability
      ! Todo: These should be either set by the user or based on raynolds number
      real(kind=rp) :: min_permeability, max_permeability

   contains

      !> @brief Initialize the topology
      !><
      !> @param[in] llf The lower left front corner
      !> @param[in] urc The upper right corner
      !> @param[in] n The number of cells in each direction
      procedure, pass(this) :: init => init_topology

      !> @brief Free the topology
      procedure, pass(this) :: free => free_topology

      !> @brief Update the topology
      procedure, pass(this) :: update => update_topology
   end type topology_t

contains

   ! ========================================================================= !
   ! Public routines
   ! ========================================================================= !

   !> @brief Initialize the topology
   !>
   !> @param[in] lower_left_front The lower left front corner
   !> @param[in] upper_right_back The upper right corner
   !> @param[in] resolution The number of cells in each direction
   subroutine init_topology(this, neko_case, resolution, reynolds)
      use json_utils, only: json_get
      implicit none

      class(topology_t), intent(inout) :: this
      type(case_t), intent(in) :: neko_case
      integer, dimension(3), intent(in) :: resolution
      real(kind=rp), optional :: reynolds

      integer :: i
      real(kind=rp) :: x, y, z
      integer :: n

      ! Allocate the design array
      this%n = resolution
      allocate (this%design(this%n(1), this%n(2), this%n(3)))
      this%design = 0.0_rp

      this%max_permeability = 0.0_rp
      if (present(reynolds)) then
         this%min_permeability = -reynolds
      else
         this%min_permeability = -1000.0_rp
      end if

      ! ---------------------------------------------------------------------- !
      ! Set the bounding box for the topology
      ! ---------------------------------------------------------------------- !
      ! Assign infinity to the corners
      this%llf = huge(this%llf)
      this%urb = -huge(this%urb)

      ! Loop through and detect the bounding box
      this%design_domain => neko_point_zone_registry%get_point_zone("design_domain")
      do i = 1, this%design_domain%size
         x = neko_case%fluid%u%dof%x(this%design_domain%mask(i), 1, 1, 1)
         y = neko_case%fluid%u%dof%y(this%design_domain%mask(i), 1, 1, 1)
         z = neko_case%fluid%u%dof%z(this%design_domain%mask(i), 1, 1, 1)

         this%llf(1) = min(this%llf(1), x)
         this%llf(2) = min(this%llf(2), y)
         this%llf(3) = min(this%llf(3), z)

         this%urb(1) = max(this%urb(1), x)
         this%urb(2) = max(this%urb(2), y)
         this%urb(3) = max(this%urb(3), z)
      end do

      ! ---------------------------------------------------------------------- !
      ! Initialize the resistance array
      ! ---------------------------------------------------------------------- !

      n = neko_case%fluid%dm_Xh%size()
      allocate (resistance(n))
      call rzero(resistance, n)

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_map(resistance, resistance_d, n)
      end if

   end subroutine init_topology

   !> @brief Free the topology
   subroutine free_topology(this)
      class(topology_t), intent(inout) :: this

      deallocate (this%design)
      deallocate (resistance)


   end subroutine free_topology

   !> @brief Update the topology
   subroutine update_topology(this, neko_case)
      class(topology_t), intent(inout) :: this
      type(case_t), intent(in) :: neko_case

      this%design = this%design + 0.3_rp

      this%design = min(this%design, 1.0_rp)
      this%design = max(this%design, 0.0_rp)

      call update_permeability(this, neko_case)
   end subroutine update_topology

   !> @brief Compute the permeability force term
   !! @details This function computes the permeability force term. This is
   !!          done by computing the permeability at each point in the domain
   !!          and then multiplying it by the velocity at that point. This
   !!          is then added to the force term.
   !!
   !! @param[inout] f The force term
   !! @param[in] t The current time
   subroutine topopt_permeability_force(f, t)
      class(fluid_user_source_term_t), intent(inout) :: f
      real(kind=rp), intent(in) :: t
      type(field_t), pointer :: u, v, w

      u => neko_field_registry%get_field('u')
      v => neko_field_registry%get_field('v')
      w => neko_field_registry%get_field('w')

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_col3(f%u_d, u%x_d, resistance_d, f%dm%size())
         call device_col3(f%v_d, v%x_d, resistance_d, f%dm%size())
         call device_col3(f%w_d, w%x_d, resistance_d, f%dm%size())
      else
         call col3(f%u, u%x, resistance, f%dm%size())
         call col3(f%v, v%x, resistance, f%dm%size())
         call col3(f%w, w%x, resistance, f%dm%size())
      end if
   end subroutine topopt_permeability_force

   ! ========================================================================= !
   ! Private routines
   ! ========================================================================= !

   !> @brief Update the permeability
   !! @details This function updates the permeability based on the current
   !!          design. This is done by interpolating the design to the
   !!          permeability array.
   !!
   !! @param[in] neko_case The neko case
   subroutine update_permeability(this, neko_case)
      type(topology_t), intent(inout) :: this
      type(case_t), intent(in) :: neko_case

      integer :: i, current_id, total_size
      real(kind=rp) :: x, y, z
      real(kind=rp), dimension(:,:,:,:), pointer :: x_ptr, y_ptr, z_ptr

      x_ptr => neko_case%fluid%u%dof%x
      y_ptr => neko_case%fluid%u%dof%y
      z_ptr => neko_case%fluid%u%dof%z

      do i = 1, this%design_domain%size
         current_id = this%design_domain%mask(i)
         x = x_ptr(current_id, 1, 1, 1)
         y = y_ptr(current_id, 1, 1, 1)
         z = z_ptr(current_id, 1, 1, 1)

         resistance(current_id) = interpolate_permeability(this, x, y, z)
      end do

      ! Copy to device
      if (NEKO_BCKND_DEVICE .eq. 1) then
         total_size = neko_case%fluid%dm_Xh%size()
         call device_memcpy(resistance, resistance_d, total_size, host_to_device)
      end if

   end subroutine update_permeability

   !> @brief Interpolate the permeability
   !! @details This function interpolates the permeability from the topology
   !!          array. The permeability is a function of the topology array
   !!          and the permeability bounds. Currently we just look up the
   !!          design value in the topology array and return the permeability
   !!          based on that.
   !!
   !! @param[in] x The x coordinate
   !! @param[in] y The y coordinate
   !! @param[in] z The z coordinate
   !! @return The permeability
   pure function interpolate_permeability(this, x, y, z) result(per)
      class(topology_t), intent(in) :: this
      real(kind=rp), intent(in) :: x
      real(kind=rp), intent(in) :: y
      real(kind=rp), intent(in) :: z
      real(kind=rp) :: per
      real(kind=rp) :: dx, dy, dz
      integer :: i, j, k
      real(kind=rp) :: local_design

      !Todo: Implement filter operations here.
      !      We might want to use some dropoff filtering outside the domain to
      !      ensure that the permeability is not too discontinuous.

      if (x < this%llf(1) .or. x > this%urb(1) .or. &
          y < this%llf(2) .or. y > this%urb(2) .or. &
          z < this%llf(3) .or. z > this%urb(3)) then
         per = this%max_permeability
      else
         dx = (x - this%llf(1))/(this%urb(1) - this%llf(1))
         dy = (y - this%llf(2))/(this%urb(2) - this%llf(2))
         dz = (z - this%llf(3))/(this%urb(3) - this%llf(3))

         i = min(max(floor(dx*this%n(1)), 1), this%n(1))
         j = min(max(floor(dy*this%n(2)), 1), this%n(2))
         k = min(max(floor(dz*this%n(3)), 1), this%n(3))

         local_design = this%design(i, j, k)

         per = local_design*this%min_permeability + (1.0_rp - local_design)*this%max_permeability
      end if
   end function interpolate_permeability

end module topology_module
