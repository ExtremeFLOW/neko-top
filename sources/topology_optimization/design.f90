!> @file topology_module.f90
!! @brief Module containing the topology type
!! @details This module contains the topology type, which is used to describe
!!          the designs in topology optimization.

module design_module

  use num_types, only: rp
  use case, only: case_t
  use fluid_user_source_term, only: fluid_user_source_term_t
  use field, only: field_t
  use field_registry, only: neko_field_registry
  use point_zone, only: point_zone_t

  use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr

  implicit none

  ! ========================================================================= !
  ! Module interface
  ! ========================================================================= !
  private

  public :: design_t, topopt_permeability_force

  ! ========================================================================= !
  ! Global variables
  ! ========================================================================= !

  !> @brief The resistance array
  real(kind=rp), dimension(:), allocatable :: resistance
  !> @brief Pointer to the resistance array on the device
  type(c_ptr) :: resistance_d = c_null_ptr
  !> @brief Pointer to the design domain mask
  integer, dimension(:), pointer :: design_domain_mask
  !> @brief Pointer to the design domain mask on the device
  type(c_ptr) :: design_domain_mask_d = c_null_ptr

  ! ========================================================================= !
  ! Topology type
  ! ========================================================================= !

  !> @brief The topology type, which is used to describe the designs in
  !! topology optimization
  type :: design_t

     !> @brief array describing the topology.
     class(field_t), private, pointer :: design_field
     !> @brief Pointer to the design domain.
     class(point_zone_t), private, pointer :: design_domain => null()
     !> @brief Pointer to the neko case.
     class(case_t), private, pointer :: neko_case => null()

     ! Limits of the permeability.
     real(kind=rp), private :: perm_0, perm_1, perm_penalty

     !> Size parameters
     integer :: total_size, design_size

   contains

     !> @brief Initialize the topology
     procedure, public, pass(this) :: init => init_design

     !> @brief Free the topology
     procedure, public, pass(this) :: free => free_design

     !> @brief Update the topology
     procedure, public, pass(this) :: update => update_design
  end type design_t

contains

  ! ========================================================================= !
  ! Public routines
  ! ========================================================================= !

  !> @brief Initialize the topology
  !>
  !> @param[inout] this The topology
  !> @param[in] neko_case The neko case
  !> @param[in] resolution The resolution of the topology
  subroutine init_design(this, neko_case)
    use neko_config, only: NEKO_BCKND_DEVICE
    use device, only: device_map, device_memcpy, HOST_TO_DEVICE
    use math, only: rzero
    use logger, only: neko_log
    use device_math, only: device_rzero
    use point_zone_registry, only: neko_point_zone_registry
    use json_utils, only: json_get, json_get_or_default
    implicit none

    ! Subroutine arguments
    class(design_t), intent(inout) :: this
    type(case_t), intent(inout), target :: neko_case

    ! Local variables
    real(kind=rp) :: reynolds

    ! Save the neko case object.
    this%neko_case => neko_case

    ! ---------------------------------------------------------------------- !
    ! Assign local variables
    ! ---------------------------------------------------------------------- !

    call json_get_or_default(this%neko_case%params, "topopt.perm_penalty", &
                             this%perm_penalty, 1.0_rp)
    call json_get_or_default(this%neko_case%params, "topopt.perm_0", &
                             this%perm_0, 1000.0_rp)
    call json_get_or_default(this%neko_case%params, "topopt.perm_1", &
                             this%perm_1, 0.0_rp)
    call json_get(this%neko_case%params, 'case.fluid.Re', reynolds)

    ! ---------------------------------------------------------------------- !
    ! Set the design domain
    ! ---------------------------------------------------------------------- !

    if (neko_point_zone_registry%point_zone_exists("design_domain") ) then
       this%design_domain => &
         neko_point_zone_registry%get_point_zone("design_domain")
    else
       call neko_log%error("design_domain point zone does not exist")
    end if

    design_domain_mask => this%design_domain%mask
    design_domain_mask_d = this%design_domain%mask_d
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_log%warning("We assume that NEKO is using the old design " // &
                             "domain masks. Please check the status of PR " // &
                             "on masks if results seem wonky.")

       design_domain_mask = design_domain_mask -1
       call device_memcpy(design_domain_mask, design_domain_mask_d,&
                          this%design_domain%size, HOST_TO_DEVICE, &
                          .true. &
                          )
       design_domain_mask = design_domain_mask + 1

    end if
    ! ---------------------------------------------------------------------- !
    ! Initialize the design field
    ! ---------------------------------------------------------------------- !

    if (.not. neko_field_registry%field_exists("design")) then
       call neko_field_registry%add_field(this%neko_case%fluid%dm_Xh, "design")
    end if

    this%design_field => neko_field_registry%get_field_by_name("design")

    call this%neko_case%f_out%fluid%append(this%design_field)

    ! Initialize the design
    call rzero(this%design_field%x, this%design_field%dof%size())
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_rzero(this%design_field%x_d, this%design_field%dof%size())
    end if

    ! ---------------------------------------------------------------------- !
    ! Initialize the resistance array
    ! ---------------------------------------------------------------------- !

    this%total_size = this%design_field%dof%size()
    this%design_size = this%design_domain%size

    allocate (resistance(this%total_size))
    call rzero(resistance, this%total_size)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(resistance, resistance_d, this%total_size)
    end if

  end subroutine init_design

  !> @brief Free the topology
  subroutine free_design(this)
    use device, only: device_free
    implicit none

    class(design_t), intent(inout) :: this

    deallocate (resistance)
    call device_free(resistance_d)

  end subroutine free_design

  !> @brief Update the topology
  !!
  !! @todo This is currently just a dummy function. We need to implement
  !!       the actual topology optimization algorithm here.
  subroutine update_design(this, converged)
    use neko_config, only: NEKO_BCKND_DEVICE
    use math, only: cadd
    use math_ext, only: cadd_mask
    use device_math_ext, only: device_cadd_mask
    use device_math, only: device_cadd
    use device, only: device_memcpy, host_to_device
    implicit none

    class(design_t), intent(inout) :: this
    logical, intent(out) :: converged
    real(kind=rp) :: max_norm

    ! ---------------------------------------------------------------------- !
    ! Update the design
    ! ---------------------------------------------------------------------- !

    associate(x => this%design_field%x, &
              x_d => this%design_field%x_d, &
              mask => this%design_domain%mask, &
              mask_d => this%design_domain%mask_d &
              )

      if (.not. NEKO_BCKND_DEVICE .eq. 1) then
         call cadd_mask(x, 0.3_rp, this%total_size, mask, this%design_size)
      else
         call device_cadd_mask(x_d, 0.3_rp, this%total_size, &
                               mask_d, this%design_size)
      end if
    end associate

    ! ---------------------------------------------------------------------- !
    ! Compute the maximum norm of the change in the design and check if we
    ! have converged.
    ! ---------------------------------------------------------------------- !

    max_norm = 1.0_rp
    converged = max_norm < 0.01_rp

    ! ---------------------------------------------------------------------- !
    ! Finish the update by updating the permeability field
    ! ---------------------------------------------------------------------- !

    call update_permeability(this)

  end subroutine update_design

  !> @brief Compute the permeability force term
  !! @details This function computes the permeability force term. This is
  !!          done by computing the permeability at each point in the domain
  !!          and then multiplying it by the velocity at that point. This
  !!          is then added to the force term.
  !!
  !! @param[inout] f The force term
  !! @param[in] t The current time
  subroutine topopt_permeability_force(f, t)
    use neko_config, only: NEKO_BCKND_DEVICE
    use math, only: col3
    use math_ext, only: col3_mask
    use device_math, only: device_col3
    use device_math_ext, only: device_col3_mask
    use field, only: field_t
    use field_registry, only: neko_field_registry

    implicit none

    class(fluid_user_source_term_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t
    type(field_t), pointer :: u, v, w
    integer :: total_size, design_size

    u => neko_field_registry%get_field('u')
    v => neko_field_registry%get_field('v')
    w => neko_field_registry%get_field('w')

    design_size = size(design_domain_mask)
    total_size = f%dm%size()

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col3_mask(f%u_d, u%x_d, resistance_d, total_size, &
                             design_domain_mask_d, design_size)
       call device_col3_mask(f%v_d, v%x_d, resistance_d, total_size, &
                             design_domain_mask_d, design_size)
       call device_col3_mask(f%w_d, w%x_d, resistance_d, total_size, &
                             design_domain_mask_d, design_size)
    else
       call col3_mask(f%u, u%x, resistance, total_size, &
                      design_domain_mask, design_size)
       call col3_mask(f%v, v%x, resistance, total_size, &
                      design_domain_mask, design_size)
       call col3_mask(f%w, w%x, resistance, total_size, &
                      design_domain_mask, design_size)
    end if
  end subroutine topopt_permeability_force

  ! ========================================================================= !
  ! Private routines
  ! ========================================================================= !

  !> @brief Update the permeability
  !! @details This function updates the permeability based on the current
  !! design. This is done by interpolating the design and updating the
  !! pointwise permeability force.
  !!
  !! From Andreasen et al. (2008) we have that the permeability is given by
  !! \f[
  !!    \kappa = \kappa_{0} + (\kappa_{1} - \kappa_{0})\phi(\mathbf{x})
  !!    \frac{\beta + 1}{\beta + \phi(\mathbf{x})},
  !! \f]
  !! where \f$\kappa_{0}\f$ and \f$\kappa_{1}\f$ are the permeability at solid
  !! and fluid, respectively, \f$\phi(\mathbf{x})\f$ is the
  !! design variable at point \f$\mathbf{x}\f$, and \f$\beta\f$ is a penalty
  !! parameter.
  !!
  !! @param[in] neko_case The neko case
  subroutine update_permeability(this)
    use neko_config, only: NEKO_BCKND_DEVICE
    use math, only: col2, cmult, invcol1, cadd, add2, cfill
    use device_math, only: device_col2, device_cmult, device_invcol1, &
      device_cadd, device_add2, device_cfill

    implicit none

    ! Arguments
    type(design_t), intent(inout) :: this

    ! Local variables
    real(kind=rp) :: constant

    constant = (this%perm_0 - this%perm_1) &
      * (this%perm_penalty + 1.0_rp)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cfill(resistance_d, this%perm_penalty, this%total_size)
       call device_add2(resistance_d, this%design_field%x_d, this%total_size)
       call device_invcol1(resistance_d, this%total_size)
       call device_col2(resistance_d, this%design_field%x_d, this%total_size)
       call device_cmult(resistance_d, constant, this%total_size)
       call device_cadd(resistance_d, this%perm_1, this%total_size)

       ! Multiply by -1, since we wish to affect the velocity negatively.
       call device_cmult(resistance_d, -1.0_rp, this%total_size)
    else
       call cfill(resistance, this%perm_penalty, this%total_size)
       call add2(resistance, this%design_field%x, this%total_size)
       call invcol1(resistance, this%total_size)
       call col2(resistance, this%design_field%x, this%total_size)
       call cmult(resistance, constant, this%total_size)
       call cadd(resistance, this%perm_1, this%total_size)

       ! Multiply by -1, since we wish to affect the velocity negatively.
       call cmult(resistance, -1.0_rp, this%total_size)
    end if

  end subroutine update_permeability

end module design_module
