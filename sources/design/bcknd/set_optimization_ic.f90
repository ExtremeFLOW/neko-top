!> @file set_optimization_ic.f90
!! @copyright
!! Copyright (c) 2025, The Neko-TOP Authors
!! All rights reserved.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions
!! are met:
!!
!!   * Redistributions of source code must retain the above copyright
!!     notice, this list of conditions and the following disclaimer.
!!
!!   * Redistributions in binary form must reproduce the above
!!     copyright notice, this list of conditions and the following
!!     disclaimer in the documentation and/or other materials provided
!!     with the distribution.
!!
!!   * Neither the name of the authors nor the names of its
!!     contributors may be used to endorse or promote products derived
!!     from this software without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
!! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
!! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
!! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
!! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
!! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
!! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!! POSSIBILITY OF SUCH DAMAGE.
!
!> Optimization initial condition
module optimization_ic
  use gather_scatter, only : gs_t, GS_OP_ADD
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use device_math, only : device_col2
  use device, only : device_memcpy, HOST_TO_DEVICE
  use field, only : field_t
  use utils, only : neko_error, filename_chsuffix, filename_suffix, &
       neko_warning, NEKO_FNAME_LEN, extract_fld_file_index
  use coefs, only : coef_t
  use math, only : col2, cfill, cfill_mask
  use json_module, only : json_file
  use json_utils, only: json_get, json_get_or_default
  use point_zone, only: point_zone_t
  use point_zone_registry, only: neko_point_zone_registry
  use registry, only: neko_registry
  use nekotop_logger, only: nekotop_log, LOG_SIZE
  use fld_file_data, only: fld_file_data_t
  use fld_file, only: fld_file_t
  use checkpoint, only: chkp_t
  use file, only: file_t
  use global_interpolation, only: global_interpolation_t
  use interpolation, only: interpolator_t
  use space, only: space_t, GLL
  implicit none
  private

  interface set_optimization_ic
     module procedure set_optimization_ic_int
  end interface set_optimization_ic

  public :: set_optimization_ic

contains

  !> Set optimization initial design field (builtin)
  !! @details Set optimization initial condition using one of the builtin types
  !! currently supported:
  !! - uniform
  !! - point zone
  !! - field
  !! It should be noted that point zones currently contain functionality to
  !! perform boolean operations and will soon be extended to allow for point
  !! zones defined by boundary meshes, allowing for complex initial material
  !! distributions to be assigned.
  !! @param fld design field.
  !! @param coef Coefficient.
  !! @param gs Gather-Scatter object.
  !! @param params JSON parameters.
  subroutine set_optimization_ic_int(fld, coef, gs, params)
    type(field_t), intent(inout) :: fld
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs
    type(json_file), intent(inout) :: params

    ! Variables for retrieving JSON parameters
    character(len=:), allocatable :: type
    real(kind=rp) :: ic_value
    character(len=:), allocatable :: read_str
    character(len=:), allocatable :: fname, mesh_fname
    real(kind=rp) :: zone_value, tol
    logical :: interpolate

    call json_get(params, "type", type)

    select case (trim(type))
    case ('uniform')
       call json_get(params, 'value', ic_value)
       call set_optimization_ic_uniform(fld, ic_value)

    case ('point_zone')
       call json_get(params, 'base_value', ic_value)
       call json_get(params, 'zone_name', read_str)
       call json_get(params, 'zone_value', zone_value)

       call set_optimization_ic_point_zone(fld, ic_value, read_str, zone_value)

    case ('field')
       call json_get(params, 'file_name', fname)
       call json_get_or_default(params, 'interpolate', interpolate, .false.)
       call json_get_or_default(params, 'tolerance', tol, 0.000001_rp)
       call json_get_or_default(params, 'mesh_file_name', mesh_fname, "none")

       call set_optimization_ic_fld(fld, fname, interpolate, tol, mesh_fname)

    case default
       call neko_error('Invalid initial condition')

    end select

    call set_optimization_ic_common(fld, coef, gs)

  end subroutine set_optimization_ic_int

  !> Set design initial condition (common)
  !! @details Finalize design initial condition by distributing the initial
  !! condition across elements and multiplying by the coefficient (if any).
  !! @param fld Scalar field.
  !! @param coef Coefficient.
  !! @param gs Gather-Scatter object.
  subroutine set_optimization_ic_common(fld, coef, gs)
    type(field_t), intent(inout) :: fld
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs
    integer :: n

    n = fld%dof%size()
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(fld%x, fld%x_d, n, &
            HOST_TO_DEVICE, sync = .false.)
    end if

    ! Ensure continuity across elements for initial conditions
    call gs%op(fld%x, n, GS_OP_ADD)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col2(fld%x_d, coef%mult_d, n)
    else
       call col2(fld%x, coef%mult, n)
    end if

  end subroutine set_optimization_ic_common

  !> Uniform initial condition
  !! @details Set design initial condition to a uniform value across the domain.
  !! @param fld Scalar field.
  !! @param ic_value Desired value of the design field.
  subroutine set_optimization_ic_uniform(fld, ic_value)
    type(field_t), intent(inout) :: fld
    real(kind=rp), intent(in) :: ic_value
    integer :: n
    character(len=LOG_SIZE) :: log_buf

    call nekotop_log%message("Type : uniform")
    write (log_buf, '(A,ES12.6)') "Value: ", ic_value
    call nekotop_log%message(log_buf)

    fld = ic_value
    n = fld%dof%size()
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call cfill(fld%x, ic_value, n)
    end if

  end subroutine set_optimization_ic_uniform

  !> Point zone initial condition
  !! @details Set design initial condition to a uniform value across a point
  !! zone.
  !! @param fld Scalar field.
  !! @param base_value Base value of the design field.
  !! @param zone_name Name of the point zone.
  !! @param zone_value Desired value of the design field in the point zone.
  subroutine set_optimization_ic_point_zone(fld, base_value, zone_name, &
       zone_value)
    type(field_t), intent(inout) :: fld
    real(kind=rp), intent(in) :: base_value
    character(len=*), intent(in) :: zone_name
    real(kind=rp), intent(in) :: zone_value

    ! Internal variables
    character(len=LOG_SIZE) :: log_buf
    class(point_zone_t), pointer :: zone
    integer :: size

    call nekotop_log%message("Type      : point_zone")
    write (log_buf, '(A,ES12.6)') "Base value: ", base_value
    call nekotop_log%message(log_buf)
    call nekotop_log%message("Zone name : " // trim(zone_name))
    write (log_buf, '(A,ES12.6)') "Zone value: ", zone_value
    call nekotop_log%message(log_buf)

    size = fld%dof%size()
    zone => neko_point_zone_registry%get_point_zone(trim(zone_name))

    call set_optimization_ic_uniform(fld, base_value)
    call cfill_mask(fld%x, zone_value, size, zone%mask%get(), zone%size)

  end subroutine set_optimization_ic_point_zone

  !> Set the initial condition of the design based on a field.
  !! @detail The field is read from an `fld` file. If enabled, interpolation
  !! is also possible. In that case, the mesh coordinates can be read from
  !! another file in the `fld` field series.
  !! @param fld The design field.
  !! @param file_name The name of the "fld" file series.
  !! @param sample_idx index of the field file .f000* to read, default is
  !! -1.
  !! @param interpolate Flag to indicate wether or not to interpolate the
  !! values onto the current mesh.
  !! @param tolerance If interpolation is enabled, tolerance for finding the
  !! points in the mesh.
  !! @param sample_mesh_idx If interpolation is enabled, index of the field
  !! file where the mesh coordinates are located.
  subroutine set_optimization_ic_fld(fld, file_name, &
       interpolate, tolerance, mesh_file_name)
    type(field_t), intent(inout) :: fld
    character(len=*), intent(in) :: file_name
    logical, intent(in) :: interpolate
    real(kind=rp), intent(in) :: tolerance
    character(len=*), intent(inout) :: mesh_file_name

    character(len=LOG_SIZE) :: log_buf
    integer :: sample_idx, sample_mesh_idx
    type(fld_file_data_t) :: fld_data
    type(file_t) :: f
    logical :: mesh_mismatch

    ! ---- For the mesh to mesh interpolation
    type(global_interpolation_t) :: global_interp
    ! -----

    ! ---- For space to space interpolation
    type(space_t) :: prev_Xh
    type(interpolator_t) :: space_interp
    ! ----

    call nekotop_log%message("Type         : field")
    call nekotop_log%message("File name    : " // trim(file_name))
    write (log_buf, '(A,L1)') "Interpolation: ", interpolate
    call nekotop_log%message(log_buf)

    ! Extract sample index from the file name
    sample_idx = extract_fld_file_index(file_name, -1)

    if (sample_idx .eq. -1) then
       call neko_error("Invalid file name for the initial condition. The " // &
            "file format must be e.g. 'mean0.f00001'")
    end if

    ! Change from "field0.f000*" to "field0.fld" for the fld reader
    call filename_chsuffix(file_name, trim(file_name), 'fld')

    call fld_data%init()
    call f%init(trim(file_name))

    if (interpolate) then

       ! If no mesh file is specified, use the default file name
       if (mesh_file_name .eq. "none") then
          mesh_file_name = trim(file_name)
          sample_mesh_idx = sample_idx
       else

          ! Extract sample index from the mesh file name
          sample_mesh_idx = extract_fld_file_index(mesh_file_name, -1)

          if (sample_mesh_idx .eq. -1) then
             call neko_error("Invalid file name for the initial condition." // &
                  " The file format must be e.g. 'mean0.f00001'")
          end if

          write (log_buf, '(A,ES12.6)') "Tolerance     : ", tolerance
          call nekotop_log%message(log_buf)
          write (log_buf, '(A,A)') "Mesh file     : ", &
               trim(mesh_file_name)
          call nekotop_log%message(log_buf)

       end if ! if mesh_file_name .eq. none

       ! Read the mesh coordinates if they are not in our fld file
       if (sample_mesh_idx .ne. sample_idx) then
          call f%set_counter(sample_mesh_idx)
          call f%read(fld_data)
       end if

    end if

    ! Read the field file containing (u,v,w,p)
    call f%set_counter(sample_idx)
    call f%read(fld_data)

    !
    ! Check that the data in the fld file matches the current case.
    ! Note that this is a safeguard and there are corner cases where
    ! two different meshes have the same dimension and same # of elements
    ! but this should be enough to cover obvious cases.
    !
    mesh_mismatch = (fld_data%glb_nelv .ne. fld%msh%glb_nelv .or. &
         fld_data%gdim .ne. fld%msh%gdim)

    if (mesh_mismatch .and. .not. interpolate) then
       call neko_error("The fld file must match the current mesh! " // &
            "Use 'interpolate': 'true' to enable interpolation.")
    else if (.not. mesh_mismatch .and. interpolate) then
       call nekotop_log%warning("You have activated interpolation but you " // &
            "might still be using the same mesh.")
    end if


    ! Mesh interpolation if specified
    if (interpolate) then
       ! Issue a warning if the mesh is in single precision
       select type (ft => f%file_type)
       type is (fld_file_t)
          if (.not. ft%dp_precision) then
             call neko_warning("The coordinates read from the field file " // &
                  "are in single precision.")
             call nekotop_log%message("It is recommended to use a mesh in " // &
                  "double precision for better interpolation results.")
             call nekotop_log%message("If the interpolation does not " // &
                  "work, you can try to increase the tolerance.")
          end if
       class default
       end select

       ! Generates an interpolator object and performs the point search
       call fld_data%generate_interpolator(global_interp, fld%dof, fld%msh, &
            tolerance)

       ! Evaluate design
       ! NOTE. Is is currently the convention in the design0.f* files that
       ! the unfiltered design should be stored in the u component.
       call global_interp%evaluate(fld%x, fld_data%u%x, on_host=.false.)
       call global_interp%free()

    else ! No interpolation, just potentially from different spaces

       ! Build a space_t object from the data in the fld file
       call prev_Xh%init(GLL, fld_data%lx, fld_data%ly, fld_data%lz)
       call space_interp%init(fld%Xh, prev_Xh)

       ! Do the space-to-space interpolation
       call space_interp%map_host(fld%x, fld_data%t%x, fld_data%nelv, fld%Xh)

       call space_interp%free()

    end if

    call fld_data%free()

  end subroutine set_optimization_ic_fld

end module optimization_ic
