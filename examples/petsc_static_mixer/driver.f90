program usrneko
  use user, only: user_setup
  use neko, only: neko_init, neko_solve, neko_finalize
  use field, only: field_t
  use field_registry, only: neko_field_registry
  use case, only: case_t
  use file, only: file_t
  use initial_conditions, only: scalar_z_split_ic
  use tri_mesh, only: tri_mesh_t
  use json_utils, only: json_get, json_get_or_default
  use signed_distance, only: signed_distance_field
  use device, only: device_memcpy, HOST_TO_DEVICE
  use utils, only: neko_error
  use neko_config, only: NEKO_BCKND_DEVICE

  use filters, only: smooth_step_field, permeability_field, step_function_field
  use num_types, only: rp

  implicit none

  type(case_t) :: C

  ! Options
  character(len=:), allocatable :: base
  character(len=:), allocatable :: mesh_file_name
  character(len=:), allocatable :: distance_transform
  character(len=:), allocatable :: filter_type
  real(kind=rp), dimension(:), allocatable :: brinkman_limits
  real(kind=rp) :: brinkman_penalty

  real(kind=rp) :: json_read_scalar

  type(file_t) :: mesh_file
  type(tri_mesh_t) :: boundary_mesh

  type(field_t), pointer :: brinkman

  call user_setup(C%usr)
  call neko_init(C)

  ! ------------------------------------------------------------------------ !
  ! Read options from the input file

  base = 'case.test.'

  ! Read the mesh file name
  call json_get(C%params, base // &
    & 'mesh_file', mesh_file_name)

  ! Read the options for the permeability field
  call json_get(C%params, base // &
    & 'brinkman.limits', brinkman_limits)
  call json_get_or_default(C%params, base // &
    & 'brinkman.penalty', brinkman_penalty, 1.0_rp)

  ! Settings on how to transform the distance field to a design field
  call json_get(C%params, base // &
    & 'distance_transform.type', distance_transform)

  ! Settings on how to filter the design field
  call json_get_or_default(C%params, base // &
    & 'filter.type', filter_type, 'none')

  ! Ensure the options are valid
  if (size(brinkman_limits) .ne. 2) then
     call neko_error('brinkman_limits must be a 2 element array')
  end if

  ! ------------------------------------------------------------------------ !
  ! Load the immersed boundary mesh

  mesh_file = file_t(mesh_file_name)

  call mesh_file%read(boundary_mesh)

  if (boundary_mesh%nelv .eq. 0) then
     call neko_error('No elements in the boundary mesh')
  end if

  ! ------------------------------------------------------------------------ !
  ! Allocate the permeability field

  if (.not. neko_field_registry%field_exists('brinkman')) then
     call neko_field_registry%add_field(C%fluid%u%dof, 'brinkman')
  end if

  brinkman => neko_field_registry%get_field_by_name('brinkman')
  call C%f_out%fluid%append(brinkman)

  ! ------------------------------------------------------------------------ !
  ! Compute the permeability field

  ! Assign the signed distance field to all GLL points in the permeability
  ! field. Initally we just run a brute force loop over all GLL points and
  ! compute the signed distance function. This should be replaced with a
  ! more efficient method, such as a tree search.

  ! Select how to transform the distance field to a design field
  select case (distance_transform)
    case ('smooth_step')
     call json_get(C%params, base // &
       & 'distance_transform.value', json_read_scalar)

     call signed_distance_field(brinkman, boundary_mesh, json_read_scalar)
     call smooth_step_field(brinkman, 0.0_rp, json_read_scalar)

    case ('step')

     call json_get(C%params, base // &
       & 'distance_transform.value', json_read_scalar)

     call signed_distance_field(brinkman, boundary_mesh, json_read_scalar)
     call step_function_field(brinkman, json_read_scalar, 1.0_rp, 0.0_rp)

    case default
     call neko_error('Unknown distance transform')
  end select

  ! Run filter on the permeability field to smooth it out.
  ! This filter should initially be the classic heaviside filter, but we wish
  ! to alter it to be a PDE based filter to avoid MPI communication.
  ! The "helmholtz" solver in Neko should be able to solve the PDE filter.

  select case (filter_type)
    case ('none')
     ! Do nothing
    case default
     call neko_error('Unknown filter type')
  end select

  ! ------------------------------------------------------------------------ !
  ! Compute the permeability field and copy to the device

  call permeability_field(brinkman, &
    & brinkman_limits(1), brinkman_limits(2), brinkman_penalty)

  if (NEKO_BCKND_DEVICE .eq. 1) then
     call device_memcpy(brinkman%x, brinkman%x_d, &
                        brinkman%dof%size(), HOST_TO_DEVICE, .true.)
  end if

  call neko_solve(C)
  call neko_finalize(C)


end program usrneko
