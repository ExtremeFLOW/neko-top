!> @file perturb_output.f90
!! @copyright
!! Copyright (c) 2024-2025, The Neko-TOP Authors
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
!> Defines output for a perturbation field.
module perturb_output
  use num_types, only: rp
  use perturb_fluid_scheme, only: perturb_fluid_scheme_t
  use scalar_scheme, only: scalar_scheme_t
  use perturb_scalar_scheme, only: perturb_scalar_scheme_t
  use field_list, only: field_list_t
  use neko_config, only: NEKO_BCKND_DEVICE
  use device, only: device_memcpy, DEVICE_TO_HOST
  use output, only: output_t
  use perturb_scalars, only : perturb_scalars_t
  use fld_file, only: fld_file_t
  implicit none
  private

  !> Perturbation output.
  type, public, extends(output_t) :: perturb_output_t
     type(field_list_t) :: perturb
     logical :: always_write_mesh = .false.
   contains
     !> Constructor
     procedure, pass(this) :: init => perturb_output_init
     !> Destructor
     procedure, pass(this) :: free => perturb_output_free
     !> Sample, i.e. extract the values of the fields and write.
     procedure, pass(this) :: sample => perturb_output_sample
  end type perturb_output_t

contains

  !> Constructor.
  !! @details initialize the output.
  !! @param[inout] this The output sampler.
  !! @param[in] precision The precision of the output fields.
  !! @param[in] perturb The perturbation fluid scheme.
  !! @param[in] perturb_scalars The perturbation scalar schemes.
  !! @param[in] name The name of the .fld file.
  !! @param[in] path The path to save the .fld files.
  !! @param[in] fmt The format of the .fld files.
  !! @param[in] layout The layout of the .fld files.
  !! @param[in] always_write_mesh Whether to always write the mesh.
  subroutine perturb_output_init(this, precision, perturb, perturb_scalars, &
       name, path, fmt, layout, always_write_mesh)
    class(perturb_output_t), intent(inout) :: this
    integer, intent(in) :: precision
    class(perturb_fluid_scheme_t), intent(in), target :: perturb
    class(perturb_scalars_t), intent(in), optional, target :: perturb_scalars
    character(len=*), intent(in), optional :: name
    character(len=*), intent(in), optional :: path
    character(len=*), intent(in), optional :: fmt
    logical, intent(in), optional :: always_write_mesh
    integer, intent(in), optional :: layout
    character(len=1024) :: fname, suffix
    integer :: i, n_scalars

    call this%free()

    suffix = '.fld'
    if (present(fmt)) then
       if (fmt .eq. 'adios2') then
          suffix = '.bp'
       else if (fmt .eq. 'vtkhdf') then
          suffix = '.vtkhdf'
       end if
    end if

    if (present(always_write_mesh)) then
       this%always_write_mesh = always_write_mesh
    end if

    if (present(name) .and. present(path)) then
       fname = trim(path) // trim(name) // trim(suffix)
    else if (present(name)) then
       fname = trim(name) // trim(suffix)
    else if (present(path)) then
       fname = trim(path) // 'field' // trim(suffix)
    else
       fname = 'field' // trim(suffix)
    end if

    if (present(layout)) then
       call this%init_base(fname, precision, layout)
    else
       call this%init_base(fname, precision)
    end if

    ! Calculate total number of fields
    n_scalars = 0
    if (present(perturb_scalars)) then
       n_scalars = size(perturb_scalars%perturb_scalar_fields)
    end if

    ! Initialize field list with appropriate size
    call this%perturb%init(4 + n_scalars)

    call this%perturb%assign(1, perturb%p_pert)
    call this%perturb%assign(2, perturb%u_pert)
    call this%perturb%assign(3, perturb%v_pert)
    call this%perturb%assign(4, perturb%w_pert)

    ! Assign all scalar fields
    if (present(perturb_scalars)) then
       do i = 1, n_scalars
          call this%perturb%assign(4 + i, &
               perturb_scalars%perturb_scalar_fields(i)%s_pert)
       end do
    end if

  end subroutine perturb_output_init

  !> Destructor
  subroutine perturb_output_free(this)
    class(perturb_output_t), intent(inout) :: this

    call this%free_base()
    call this%perturb%free()

  end subroutine perturb_output_free

  !> Sample a perturbation solution at time @a t
  !! @param[inout] this The output sampler.
  !! @param[in] t The time.
  subroutine perturb_output_sample(this, t)
    class(perturb_output_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer :: i

    if (NEKO_BCKND_DEVICE .eq. 1) then

       associate(fields => this%perturb%items)
         do i = 1, size(fields)
            call device_memcpy(fields(i)%ptr%x, fields(i)%ptr%x_d, &
                 fields(i)%ptr%dof%size(), DEVICE_TO_HOST, &
                 sync = (i .eq. size(fields))) ! Sync on the last field
         end do
       end associate

    end if

    select type (ft => this%file_%file_type)
       ! Only fld files have the option to write the mesh at command
    type is (fld_file_t)
       ft%write_mesh = this%always_write_mesh
       call ft%write(this%perturb, t)
    class default
       call ft%write(this%perturb, t)
    end select

  end subroutine perturb_output_sample

end module perturb_output
