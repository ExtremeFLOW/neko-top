!> @file mapping_handler.f90
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
!> Implements the `mapping_handler_t` type.
module mapping_handler
  use neko_config, only: NEKO_BCKND_DEVICE
  use num_types, only: rp, sp, dp
  use mapping, only: mapping_wrapper_t, mapping_t, &
       mapping_factory
  use field, only: field_t
  use field_output, only: field_output_t
  use json_utils, only: json_extract_item
  use json_module, only: json_file
  use coefs, only: coef_t
  use field_math, only: field_copy
  use math, only: col2
  use device_math, only: device_col2
  use scratch_registry, only: neko_scratch_registry
  use utils, only: neko_warning, neko_error
  use vector, only:vector_t
  use neko_ext, only: field_to_vector, vector_to_field
  use gather_scatter, only : gs_op_add
  implicit none
  private

  !> Abstract class for handling mapping_cascade.
  !!
  !! @details
  !! This class is responsible for managing the mapping_cascade in a sequential
  !! manor. It is also responsible for using the chain rule to propagate
  !! sensitivity backwards throughout the system.
  type, public :: mapping_handler_t
     !> Array of mapping_cascade.
     !! @note the order really matter's here since they'll be executed in
     !! sequence.
     class(mapping_wrapper_t), allocatable :: mapping_cascade(:)
     !> The coefficients of the (space, mesh) pair.
     type(coef_t), pointer :: coef
     !> Output for design and intermediate forward-mapping stages.
     type(field_output_t) :: design_output
     !> Output for sensitivity and intermediate backward-mapping stages.
     type(field_output_t) :: sensitivity_output
     !> Field pointers used to initialize design/sensitivity outputs.
     type(field_t), pointer :: design_in => null()
     type(field_t), pointer :: design_out => null()
     type(field_t), pointer :: sensitivity_out => null()
     !> Internal storage of sensitivity propagation stages.
     type(field_t), allocatable :: sensitivity_stages(:)
     !> Flags controlling output initialization.
     logical :: output_fields_set = .false.
     logical :: outputs_initialized = .false.
     !> Flag controlling output verbose design output (default = false).
     logical :: verbose_design = .false.
     !> Flag controlling output verbose sensitivity output (default = false).
     logical :: verbose_sensitivity = .false.
     !> Precision used for design/sensitivity field outputs.
     integer :: output_precision = sp
     !> Format used for design/sensitivity field outputs.
     character(len=32) :: output_format = 'fld'
     !> Base file name for forward mapping output.
     character(len=80) :: forward_file_name = ''
     !> Base file name for sensitivity output.
     character(len=80) :: sensitivity_file_name = ''

   contains
     !> Constructor.
     procedure, pass(this) :: init_base => mapping_handler_init_base
     !> Destructor.
     procedure, pass(this) :: free => mapping_handler_free
     !> Cycle through all the mapping_cascade and return the final field
     generic :: apply_forward => mapping_handler_apply_forward_field, &
          mapping_handler_apply_forward_vector
     procedure, pass(this) :: mapping_handler_apply_forward_field
     procedure, pass(this) :: mapping_handler_apply_forward_vector
     !> Cycle backwards through all the mapping_cascade and return the
     !! sensitivity
     generic :: apply_backward => mapping_handler_apply_backward_field, &
          mapping_handler_apply_backward_vector
     procedure, pass(this) :: mapping_handler_apply_backward_field
     procedure, pass(this) :: mapping_handler_apply_backward_vector
     !> Generic interface to add a mapping to the list.
     generic :: add => add_mapping, add_json_mappings
     !> Append a new mapping to the mapping_cascade array.
     procedure, pass(this) :: add_mapping => &
          mapping_handler_add_mapping
     !> Read from the json file and initialize the mapping_cascade.
     procedure, pass(this) :: add_json_mappings => &
          mapping_handler_add_json_mappings
     !> Set names of the stored mapping-stage fields.
     procedure, pass(this) :: set_stage_names => &
          mapping_handler_set_stage_names
     !> Force a field to be continuous.
     procedure, pass(this) :: make_cts => mapping_handler_make_cts
     !> Configure output fields and initialize mapping-stage writers.
     procedure, pass(this) :: init_output_fields => &
          mapping_handler_init_output_fields
     !> Write design and forward-mapping stages.
     procedure, pass(this) :: write_design => mapping_handler_write_design
     !> Write sensitivity and backward-mapping stages.
     procedure, pass(this) :: write_sensitivity => &
          mapping_handler_write_sensitivity
     !> Reset design and sensitivity output counters.
     procedure, pass(this) :: set_output_counter => &
          mapping_handler_set_output_counter
  end type mapping_handler_t

contains

  !> Setup output objects after mappings are known.
  !! @param this The handler object.
  subroutine mapping_handler_setup_outputs(this)
    class(mapping_handler_t), intent(inout) :: this
    integer :: i, n_mappings, n_fields

    if (.not. this%output_fields_set) then
       call neko_error('Mapping output fields have not been configured')
    end if
    if (.not. associated(this%coef)) then
       call neko_error('Mapping handler coefficients are not initialized')
    end if
    if (.not. associated(this%design_out)) then
       call neko_error('Mapped design output field is not associated')
    end if
    if (.not. associated(this%design_in)) then
       call neko_error('Unmapped design input field is not associated')
    end if
    if (.not. associated(this%sensitivity_out)) then
       call neko_error('Sensitivity output field is not associated')
    end if

    n_mappings = 0
    if (allocated(this%mapping_cascade)) n_mappings = size(this%mapping_cascade)

    if (this%verbose_design) then
       ! all mapping inputs + mapped output
       n_fields = n_mappings + 1
       if (n_fields .le. 0) n_fields = 1
    else
       ! only first unmapped stage and mapped output
       if (n_mappings .gt. 0) then
          n_fields = 2
       else
          n_fields = 1
       end if
    end if

    call this%set_stage_names()

    call this%design_output%init(trim(this%forward_file_name), n_fields, &
         precision = this%output_precision, &
         format = trim(this%output_format))

    if (this%verbose_design) then
       do i = 1, n_mappings
          call this%design_output%fields%assign_to_field(i, &
               this%mapping_cascade(i)%mapping%X_in)
       end do
       call this%design_output%fields%assign_to_field(n_fields, this%design_out)
    else if (n_mappings .eq. 0) then
       call this%design_output%fields%assign_to_field(1, this%design_out)
    else
       call this%design_output%fields%assign_to_field(1, &
            this%mapping_cascade(1)%mapping%X_in)
       call this%design_output%fields%assign_to_field(2, this%design_out)
    end if

    if (this%verbose_sensitivity) then
       ! sensitivity_in + chain-rule stages + final sensitivity
       n_fields = n_mappings + 2
       allocate(this%sensitivity_stages(n_mappings + 1))
       do i = 1, n_mappings + 1
          call this%sensitivity_stages(i)%init(this%coef%dof)
       end do
       call this%sensitivity_output%init( &
            trim(this%sensitivity_file_name), n_fields, &
            precision = this%output_precision, &
            format = trim(this%output_format))
       do i = 1, n_mappings + 1
          call this%sensitivity_output%fields%assign_to_field(i, &
               this%sensitivity_stages(i))
       end do
       call this%sensitivity_output%fields%assign_to_field(n_fields, &
            this%sensitivity_out)
    else
       call this%sensitivity_output%init( &
            trim(this%sensitivity_file_name), 1, &
            precision = this%output_precision, &
            format = trim(this%output_format))
       call this%sensitivity_output%fields%assign_to_field(1, &
            this%sensitivity_out)
    end if

    this%outputs_initialized = .true.

  end subroutine mapping_handler_setup_outputs

  !> Constructor.
  subroutine mapping_handler_init_base(this, coef)
    class(mapping_handler_t), intent(inout) :: this
    type(coef_t), target, intent(in) :: coef

    call this%free()

    this%coef => coef

  end subroutine mapping_handler_init_base


  !> Destructor.
  subroutine mapping_handler_free(this)
    class(mapping_handler_t), intent(inout) :: this
    integer :: i

    if (allocated(this%mapping_cascade)) then
       do i = 1, size(this%mapping_cascade)
          call this%mapping_cascade(i)%free()
       end do
       deallocate(this%mapping_cascade)
    end if
    if (allocated(this%sensitivity_stages)) then
       do i = 1, size(this%sensitivity_stages)
          call this%sensitivity_stages(i)%free()
       end do
       deallocate(this%sensitivity_stages)
    end if

    call this%design_output%free()
    call this%sensitivity_output%free()

    this%outputs_initialized = .false.
    this%output_fields_set = .false.
    this%output_precision = sp
    this%output_format = 'fld'
    this%forward_file_name = ''
    this%sensitivity_file_name = ''
    if (associated(this%design_in)) nullify(this%design_in)
    if (associated(this%design_out)) nullify(this%design_out)
    if (associated(this%sensitivity_out)) nullify(this%sensitivity_out)
    if (associated(this%coef)) nullify(this%coef)

  end subroutine mapping_handler_free

  !> apply the cascade of mapping_cascade.
  !! @param this The handler object
  !! @param X_out The mapped field (\f$\tilde{\rho}\f$)
  !! @param X_in The unmapped field (\f$\rho\f$)
  subroutine mapping_handler_apply_forward_field(this, X_out, X_in)
    class(mapping_handler_t), intent(inout) :: this
    type(field_t), intent(in) :: X_in
    type(field_t), intent(inout) :: X_out
    integer :: i
    type(field_t), pointer :: tmp_fld_in, tmp_fld_out
    integer :: temp_indices(2)

    call neko_scratch_registry%request_field(tmp_fld_in, temp_indices(1), &
         .false.)
    call neko_scratch_registry%request_field(tmp_fld_out, temp_indices(2), &
         .false.)

    ! Start by copying the first X_in into the tmp_fld_out to begin the
    ! cascade.
    call field_copy(tmp_fld_out, X_in)

    ! enforce continuity in the field
    call this%make_cts(tmp_fld_out)

    ! We enter the cascade
    if (allocated(this%mapping_cascade)) then
       do i = 1, size(this%mapping_cascade)
          ! the output from one mapping becomes the input for the next.
          call field_copy(tmp_fld_in, tmp_fld_out)
          ! apply the mapping on temp_fld
          call this%mapping_cascade(i)%mapping%apply_forward(tmp_fld_out, &
               tmp_fld_in)

       end do

    end if

    ! our final mapping should now live in tmp_fld_out
    call field_copy(X_out, tmp_fld_out)

    ! free the scratch.
    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine mapping_handler_apply_forward_field

  !> apply the cascade of mapping_cascade.
  !! @param this The handler object
  !! @param X_out The mapped vector (\f$\tilde{\rho}\f$)
  !! @param X_in The unmapped vector (\f$\rho\f$)
  subroutine mapping_handler_apply_forward_vector(this, X_out, X_in)
    class(mapping_handler_t), intent(inout) :: this
    type(vector_t), intent(in) :: X_in
    type(vector_t), intent(inout) :: X_out
    type(field_t), pointer :: tmp_fld_in, tmp_fld_out
    integer :: temp_indices(2)

    call neko_scratch_registry%request_field(tmp_fld_in, temp_indices(1), &
         .false.)
    call neko_scratch_registry%request_field(tmp_fld_out, temp_indices(2), &
         .false.)

    call vector_to_field(tmp_fld_in, X_in)
    call mapping_handler_apply_forward_field(this, tmp_fld_out, tmp_fld_in)
    call field_to_vector(X_out, tmp_fld_out)

    ! free the scratch.
    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine mapping_handler_apply_forward_vector

  !> Apply the cascade of mapping_cascade.
  !! @param this The handler object
  !! @param sens_out The sensitivity after applying the chain rule
  !! (\f$\frac{\partial F}{\partial \rho}\f$)
  !! @param sens_in The sensitivity before applying the chain rule
  !! (\f$\frac{\partial F}{\partial \tilde{\rho}}\f$)
  subroutine mapping_handler_apply_backward_field(this, sens_out, sens_in)
    class(mapping_handler_t), intent(inout) :: this
    type(field_t), intent(inout) :: sens_out
    type(field_t), intent(in) :: sens_in
    integer :: i
    type(field_t), pointer :: tmp_fld_in, tmp_fld_out
    integer :: temp_indices(2)

    call neko_scratch_registry%request_field(tmp_fld_in, temp_indices(1), &
         .false.)
    call neko_scratch_registry%request_field(tmp_fld_out, temp_indices(2), &
         .false.)

    ! Start by copying the first sens_in into the tmp_fld_out to begin the
    ! cascade.
    call field_copy(tmp_fld_out, sens_in)
    if (allocated(this%sensitivity_stages)) then
       call field_copy(this%sensitivity_stages(1), sens_in)
    end if

    ! enforce continuity in the field
    call this%make_cts(tmp_fld_out)

    ! We enter the cascade
    if (allocated(this%mapping_cascade)) then
       do i = size(this%mapping_cascade), 1, -1
          ! the output from one mapping becomes the input for the next.
          call field_copy(tmp_fld_in, tmp_fld_out)
          ! apply the mapping on temp_fld
          ! NOTE
          ! all the X_in that is required to map backward should be held
          ! internally by each mapping
          call this%mapping_cascade(i)%mapping%apply_backward(tmp_fld_out, &
               tmp_fld_in)
          if (allocated(this%sensitivity_stages)) then
             call field_copy(this%sensitivity_stages( &
                  size(this%mapping_cascade) - i + 2), tmp_fld_out)
          end if

       end do

    end if

    ! our final mapping should now live in tmp_fld_out
    call field_copy(sens_out, tmp_fld_out)

    ! free the scratch.
    call neko_scratch_registry%relinquish_field(temp_indices)


  end subroutine mapping_handler_apply_backward_field

  !> apply the cascade of mapping_cascade.
  !! @param this The handler object
  !! @param X_out The sensitivity after applying the chain rule
  !! (\f$\frac{\partial F}{\partial \rho}\f$)
  !! @param X_in The sensitivity before applying the chain rule
  !! (\f$\frac{\partial F}{\partial \tilde{\rho}}\f$)
  subroutine mapping_handler_apply_backward_vector(this, X_out, X_in)
    class(mapping_handler_t), intent(inout) :: this
    type(vector_t), intent(in) :: X_in
    type(vector_t), intent(inout) :: X_out
    type(field_t), pointer :: tmp_fld_in, tmp_fld_out
    integer :: temp_indices(2)

    call neko_scratch_registry%request_field(tmp_fld_in, temp_indices(1), &
         .false.)
    call neko_scratch_registry%request_field(tmp_fld_out, temp_indices(2), &
         .false.)

    call vector_to_field(tmp_fld_in, X_in)
    call mapping_handler_apply_backward_field(this, tmp_fld_out, tmp_fld_in)
    call field_to_vector(X_out, tmp_fld_out)

    ! free the scratch.
    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine mapping_handler_apply_backward_vector

  !> Read from the json file and initialize the mapping_cascade.
  subroutine mapping_handler_add_json_mappings(this, json, name)
    class(mapping_handler_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: name

    class(mapping_wrapper_t), dimension(:), allocatable :: temp

    ! A single mapping as its own json_file.
    type(json_file) :: mapping_subdict
    integer :: n_mappings, i, i0

    if (json%valid_path(name)) then
       ! Get the number of mapping_cascade.
       call json%info(name, n_children = n_mappings)

       if (allocated(this%mapping_cascade)) then
          i0 = size(this%mapping_cascade)
          call move_alloc(this%mapping_cascade, temp)
          allocate(this%mapping_cascade(i0 + n_mappings))
          if (allocated(temp)) then
             do i = 1, i0
                call move_alloc(temp(i)%mapping, &
                     this%mapping_cascade(i)%mapping)
             end do
          end if
       else
          i0 = 0
          allocate(this%mapping_cascade(n_mappings))
       end if

       do i = 1, n_mappings
          ! Create a new json containing just the subdict for this mapping.
          call json_extract_item(json, name, i, mapping_subdict)
          call mapping_factory(this%mapping_cascade(i + i0)%mapping, &
               mapping_subdict, this%coef)
       end do
       if (this%output_fields_set) call mapping_handler_setup_outputs(this)
    else
       ! I was considering an error, but maybe a warning is more appropriate
       call neko_warning("No mappings selected")
    end if

  end subroutine mapping_handler_add_json_mappings

  !> Add new mapping to the list.
  !! @param this The handler object
  !! @param mapping The mapping to be added.
  subroutine mapping_handler_add_mapping(this, mapping)
    class(mapping_handler_t), intent(inout) :: this
    class(mapping_t), intent(in) :: mapping
    class(mapping_wrapper_t), dimension(:), allocatable :: temp

    integer :: n_mappings, i

    if (allocated(this%mapping_cascade)) then
       n_mappings = size(this%mapping_cascade)
    else
       n_mappings = 0
    end if

    call move_alloc(this%mapping_cascade, temp)
    allocate(this%mapping_cascade(n_mappings + 1))

    if (allocated(temp)) then
       do i = 1, n_mappings
          call move_alloc(temp(i)%mapping, this%mapping_cascade(i)%mapping)
       end do
    end if

    this%mapping_cascade(n_mappings + 1)%mapping = mapping
    if (this%output_fields_set) call mapping_handler_setup_outputs(this)

  end subroutine mapping_handler_add_mapping

  !> Prefix a field name to indicate a derivative with respect to it.
  pure function mapping_handler_derivative_name(field_name) result(name)
    character(len=*), intent(in) :: field_name
    character(len=80) :: name

    name = 'dFd_' // trim(field_name)

  end function mapping_handler_derivative_name

  !> Set the names of the stored mapping and sensitivity stages.
  !! @param this The handler object.
  subroutine mapping_handler_set_stage_names(this)
    class(mapping_handler_t), intent(inout) :: this
    integer :: i, j
    character(len=80) :: previous_name, current_name

    ! Since each mapping stores their unmapped field (for chain rule purposes)
    ! and stores the name for the mapped field, we need to loop through and
    ! name them correctly
    previous_name = trim(this%design_in%name)
    if (allocated(this%mapping_cascade)) then
       do i = 1, size(this%mapping_cascade)
          current_name = trim(this%mapping_cascade(i)%mapping%fld_name)
          this%mapping_cascade(i)%mapping%X_in%name = trim(previous_name)
          previous_name = trim(current_name)
       end do
    end if

    ! same thing for sensitivities
    this%sensitivity_out%name = mapping_handler_derivative_name( &
         this%design_in%name)
    if (allocated(this%sensitivity_stages)) then
       this%sensitivity_stages(1)%name = mapping_handler_derivative_name( &
            this%design_out%name)
       if (allocated(this%mapping_cascade)) then
          do i = size(this%mapping_cascade), 1, -1
             j = size(this%mapping_cascade) - i + 2
             this%sensitivity_stages(j)%name = &
                  mapping_handler_derivative_name( &
                  this%mapping_cascade(i)%mapping%X_in%name)
          end do
       end if
    end if

  end subroutine mapping_handler_set_stage_names

  !> Configure fields used by design/sensitivity output writers.
  !! @param this The handler object.
  !! @param design_in Initial unmapped design field.
  !! @param design_out Final mapped design field.
  !! @param sensitivity_out Final backward-mapped sensitivity field.
  !! @param[in] verbose_design If true, output all forward cascade stages.
  !! @param[in] verbose_sensitivity If true, output all backward stages.
  !! @param[in] output_precision Output precision (sp or dp).
  !! @param[in] output_format Output format for design/sensitivity fields.
  !! @param[in] forward_file_name Base file name for forward mapping output.
  !! @param[in] sensitivity_file_name Base file name for sensitivity output.
  subroutine mapping_handler_init_output_fields(this, design_in, design_out, &
       sensitivity_out, verbose_design, verbose_sensitivity, &
       output_precision, output_format, forward_file_name, &
       sensitivity_file_name)
    class(mapping_handler_t), intent(inout) :: this
    type(field_t), target, intent(inout) :: design_in
    type(field_t), target, intent(inout) :: design_out
    type(field_t), target, intent(inout) :: sensitivity_out
    logical, intent(in) :: verbose_design
    logical, intent(in) :: verbose_sensitivity
    integer, intent(in) :: output_precision
    character(len=*), intent(in) :: output_format
    character(len=*), intent(in) :: forward_file_name
    character(len=*), intent(in) :: sensitivity_file_name

    this%design_in => design_in
    this%design_out => design_out
    this%sensitivity_out => sensitivity_out
    this%verbose_design = verbose_design
    this%verbose_sensitivity = verbose_sensitivity
    if (output_precision .ne. sp .and. output_precision .ne. dp) then
       call neko_error('output_precision must be either sp or dp')
    end if
    this%output_precision = output_precision
    this%output_format = trim(output_format)
    this%forward_file_name = trim(forward_file_name)
    this%sensitivity_file_name = trim(sensitivity_file_name)
    this%output_fields_set = .true.

    call mapping_handler_setup_outputs(this)

  end subroutine mapping_handler_init_output_fields

  !> Write design-related output fields.
  !! @param this The handler object.
  !! @param[in] idx Output sample index.
  subroutine mapping_handler_write_design(this, idx)
    class(mapping_handler_t), intent(inout) :: this
    integer, intent(in) :: idx

    if (.not. this%outputs_initialized) then
       call neko_error('Mapping design output is not initialized')
    end if
    call this%design_output%sample(real(idx, kind=rp))

  end subroutine mapping_handler_write_design

  !> Write sensitivity-related output fields.
  !! @param this The handler object.
  !! @param[in] idx Output sample index.
  subroutine mapping_handler_write_sensitivity(this, idx)
    class(mapping_handler_t), intent(inout) :: this
    integer, intent(in) :: idx

    if (.not. this%outputs_initialized) then
       call neko_error('Mapping sensitivity output is not initialized')
    end if
    call this%sensitivity_output%sample(real(idx, kind=rp))

  end subroutine mapping_handler_write_sensitivity

  !> Reset design-related output counters.
  !! @param this The handler object.
  !! @param[in] idx Output sample index.
  subroutine mapping_handler_set_output_counter(this, idx)
    class(mapping_handler_t), intent(inout) :: this
    integer, intent(in) :: idx

    if (.not. this%outputs_initialized) return

    call this%design_output%set_counter(idx)
    call this%sensitivity_output%set_counter(idx)

  end subroutine mapping_handler_set_output_counter

  !> Force a field to be continuous.
  !! This can be done in many ways, here it is a simple average.
  !! @param this The handler object
  !! @param fld The field to be made continuous.
  subroutine mapping_handler_make_cts(this, fld)
    class(mapping_handler_t), intent(inout) :: this
    type(field_t), intent(inout) :: fld

    call this%coef%gs_h%op(fld, gs_op_add)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col2(fld%x_d, this%coef%mult_d, fld%size())
    else
       call col2(fld%x, this%coef%mult, fld%size())
    end if

  end subroutine mapping_handler_make_cts
end module mapping_handler
