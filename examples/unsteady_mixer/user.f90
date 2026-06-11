! User module for the user defined simulation component
module user
  use user_intf, only: user_t
  use json_module, only: json_file
  use steady_simcomp, only: steady_simcomp_t
  use simcomp_executor, only: neko_simcomps
  use num_types, only : rp
  use field, only : field_t
  use vector, only : vector_t
  use field_list, only : field_list_t
  use field_dirichlet, only : field_dirichlet_t
  use time_state, only : time_state_t
  use registry, only : neko_registry
  use math, only : rzero, copy, chsign, masked_scatter_copy_0
  use device_math, only: device_copy, device_cmult, device_masked_scatter_copy_0
  use neko_config, only: NEKO_BCKND_DEVICE
  use operators, only: curl
  use scratch_registry, only : neko_scratch_registry
  use device, only: HOST_TO_DEVICE, device_memcpy
  implicit none
  !> Case parameters
  ! To define the initial boundary conditions we don't wish to introduce a
  ! discontinuity,
  !
  !    DONT                    DO (but smoother)
  !        _______               ______
  !       |                     /
  !       |                    /
  ! -------             -------
  !
  ! hence,
  !> we use a logistic function defined by
  !! $ \frac{L}{1 + exp{-k(z-z_0)}}$
  real(kind=rp), parameter :: L = 1.0_rp
  real(kind=rp), parameter :: k = 20.0_rp
  real(kind=rp), parameter :: z_0 = 0.5_rp

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%dirichlet_conditions => user_bc
    user%initial_conditions => scalar_ic
  end subroutine user_setup

  !> user-defined boundary condition
  subroutine user_bc(fields, bc, time)
    type(field_list_t), intent(inout) :: fields
    type(field_dirichlet_t), intent(in) :: bc
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: u, v, w, s
    type(vector_t), pointer :: bc_u, bc_v, bc_w, bc_s
    real(kind=rp) :: x, y, z
    integer :: i, idx
    logical :: is_fluid
    logical, save :: fluid_initialized = .false.
    logical, save :: scalar_initialized = .false.

    if (bc%msk(0) .eq. 0) return ! No boundary points to apply BCs to

    is_fluid = (fields%items(1)%ptr%name .eq. 'u')

    if (is_fluid) then

       u => fields%get("u")
       v => fields%get("v")
       w => fields%get("w")

       if (.not. fluid_initialized) then
          call neko_registry%add_vector(bc%msk(0), "user_bc_u", .true.)
          call neko_registry%add_vector(bc%msk(0), "user_bc_v", .true.)
          call neko_registry%add_vector(bc%msk(0), "user_bc_w", .true.)

          bc_u => neko_registry%get_vector("user_bc_u")
          bc_v => neko_registry%get_vector("user_bc_v")
          bc_w => neko_registry%get_vector("user_bc_w")

          do i = 1, bc%msk(0)
             idx = bc%msk(i)
             x = u%dof%x(idx, 1, 1, 1)
             y = u%dof%y(idx, 1, 1, 1)
             z = u%dof%z(idx, 1, 1, 1)

             ! Inflow velocity profile is a paraboloid
             bc_u%x(i) = - (y - 0.5_rp)**2 - (z - 0.5_rp)**2 + 1.0_rp
             bc_v%x(i) = 0.0_rp
             bc_w%x(i) = 0.0_rp
          end do

          if (NEKO_BCKND_DEVICE .eq. 1) then
             call bc_u%copy_from(HOST_TO_DEVICE, sync=.false.)
             call bc_v%copy_from(HOST_TO_DEVICE, sync=.false.)
             call bc_w%copy_from(HOST_TO_DEVICE, sync=.true.)
          end if

          nullify(bc_u)
          nullify(bc_v)
          nullify(bc_w)
          fluid_initialized = .true.
       end if

       bc_u => neko_registry%get_vector("user_bc_u")
       bc_v => neko_registry%get_vector("user_bc_v")
       bc_w => neko_registry%get_vector("user_bc_w")

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_masked_scatter_copy_0(u%x_d, bc_u%x_d, bc%msk_d, &
               u%size(), bc%msk(0))
          call device_masked_scatter_copy_0(v%x_d, bc_v%x_d, bc%msk_d, &
               v%size(), bc%msk(0))
          call device_masked_scatter_copy_0(w%x_d, bc_w%x_d, bc%msk_d, &
               w%size(), bc%msk(0))
       else
          call masked_scatter_copy_0(u%x, bc_u%x, bc%msk, u%size(), bc%msk(0))
          call masked_scatter_copy_0(v%x, bc_v%x, bc%msk, v%size(), bc%msk(0))
          call masked_scatter_copy_0(w%x, bc_w%x, bc%msk, w%size(), bc%msk(0))
       end if

    else
       s => fields%get("s")

       if (.not. scalar_initialized) then
          call neko_registry%add_vector(bc%msk(0), "user_bc_s", .true.)
          bc_s => neko_registry%get_vector("user_bc_s")

          do i = 1, bc%msk(0)
             idx = bc%msk(i)
             z = s%dof%z(idx, 1, 1, 1)
             ! Inflow scalar profile is a sigmoid separating the two species
             bc_s%x(i) = L / (1.0_rp + exp(-k*(z - z_0)))
          end do
          if (NEKO_BCKND_DEVICE .eq. 1) then
             call bc_s%copy_from(HOST_TO_DEVICE, sync=.true.)
          end if

          nullify(bc_s)
          scalar_initialized = .true.
       end if

       bc_s => neko_registry%get_vector("user_bc_s")

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_masked_scatter_copy_0(s%x_d, bc_s%x_d, bc%msk_d, &
               s%size(), bc%msk(0))
       else
          call masked_scatter_copy_0(s%x, bc_s%x, bc%msk, s%size(), bc%msk(0))
       end if

    end if
  end subroutine user_bc

  !> user-defined initial condition
  subroutine scalar_ic(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields
    type(field_t), pointer :: s
    integer :: i

    if (scheme_name .eq. 'fluid') return

    ! Initial scalar profile is a sigmoid separating the two species
    s => fields%get("s")
    do i = 1, s%dof%size()
       s%x(i,1,1,1) = L / (1.0_rp + exp(-k*(s%dof%z(i,1,1,1) - z_0)))
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(s%x, s%x_d, s%size(), HOST_TO_DEVICE, sync=.false.)
    end if


  end subroutine scalar_ic

end module user
