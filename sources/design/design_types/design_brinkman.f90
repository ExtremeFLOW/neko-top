!> @file design_brinkman.f90
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

! Implements the `brinkman_design_t` type.
module brinkman_design
  use num_types, only: rp, sp, dp
  use field, only: field_t
  use json_module, only: json_file
  use mapping_handler, only: mapping_handler_t
  use coefs, only: coef_t
  use adjoint_fluid_pnpn, only: adjoint_fluid_pnpn_t
  use scratch_registry, only: neko_scratch_registry
  use point_zone_registry, only: neko_point_zone_registry
  use point_zone, only: point_zone_t
  use mask_ops, only: mask_exterior_const
  use neko_config, only: NEKO_BCKND_DEVICE
  use design, only: design_t
  use simulation_m, only: simulation_t
  use simple_brinkman_source_term, only: simple_brinkman_source_term_t
  use vector, only: vector_t
  use matrix, only: matrix_t
  use vector_math, only: vector_cmult
  use matrix_math, only: matrix_cmult
  use math, only: copy, col2
  use device_math, only: device_copy, device_col2
  use registry, only: neko_registry
  use neko_ext, only: field_to_vector, vector_to_field
  use optimization_ic, only: set_optimization_ic
  use field_math, only: field_rzero
  use json_utils, only: json_get, json_get_or_default
  use utils, only: neko_error
  implicit none
  private

  !> A topology optimization design variable
  type, extends(design_t), public :: brinkman_design_t
     private

     ! TODO
     ! in the future make this a derived type of a `design_variable`
     ! type, public, extends(design_variable_t) :: brinkman_design_t
     !> the unfilitered design
     type(field_t), pointer :: design_indicator

     !> the mapped coefficient (Brinkman term)
     ! TODO
     ! NOTE: Tim, right now we only have the brinkman term
     ! in the future we may map to other coeeficients for other equations...
     ! in which case this should be a field list
     !
     ! or as I describe below, we also have multiple constraints,
     ! so a list-of-lists may be the correct way forward
     type(field_t), pointer :: brinkman_amplitude

     ! NOTE:
     ! again, we have to be so clear with nomenclature.
     ! If we have an objective function F.
     ! I also like to use \rho to describe the design_indicator
     !
     ! and let's say, for completness, we're doing conjugate heat transfer,
     ! so we have to map
     ! to 3 coefficients, \chi, C and \kappa.
     !
     ! then we will perform 3 mapping,
     ! \rho -> \chi
     ! \rho -> C
     ! \rho -> \kappa
     !
     ! Then during the forward/adjoint looping there will be an additional
     ! object, the "objective_function" object that will be responsible for
     ! computing the sensitivity of the objective function with respect to the
     ! coefficients
     ! ie,
     ! dF/d\chi, dF/dC and dF/d\kappa
     !
     ! What I'm calling "sensitivity" here, is the sensitivity with respect to
     ! the design indicator
     ! so dF/d\rho
     !
     ! so the proceedure "map_backwards" will take in the field list
     ! dF/d\chi, dF/dC and dF/d\kappa
     !and chain rule it's way back to
     ! dF/d\rho
     ! and store it here          v
     type(field_t), pointer :: sensitivity
     ! have a field list here
     ! type(filed_list_t), public :: constraint_sensitivity
     ! HOWEVER !
     ! What is a bit confusing to me... is how we'll deal with constraints.
     !
     ! Because in principle we could have constraints C1, C2, C3 etc
     ! Implying we also want dC1/d\rho, dC2/d\rho etc
     ! So this "sensitivity" should also be a list...
     !
     ! So I bet you'll have a nice abstract way of constructing this in the
     ! future but I think a sort of list-of-lists would be nice.
     !
     ! For F:
     ! \rho -> \tild{\rho} -> \chi
     ! \rho -> \tild{\rho} -> C
     ! \rho -> \tild{\rho} -> \kappa
     !
     ! For C1:
     ! \rho -> \tild{\rho}  (eg for a volume constraint or so)
     !
     ! For C2:
     ! \rho -> \tild{\rho}  (for something else)
     !
     ! etc..
     !
     ! So now we have multiple instances of the "objective" type,
     ! each for F, C1, C2 etc
     ! each of them can pass their dF\d_coefficents to the design
     ! and we continue from there.
     !
     ! perhaps even the "objective" type is defined as a list of objectives.

     ! Let's say way have a chain of two mappings

     ! This is still up for discussion where the mapping will eventually live.
     ! Maybe it makes more sense to keep the design as clean as possible,
     ! and then the resposibility of converting the design into coefficients
     ! used in the PDE being solved would be the responsibility of the module
     ! used for solving the PDE, ie, the simulation.
     ! Then the simulation would responsiblity for asking "do we have a scalar?"
     ! ok, let me figure out how to map this design into a coeffient in the
     ! scalar equation.
     !
     ! I guess then the mapping backwards gets confusing to me.
     ! anyway, for now, let's keep things how they are and I suppose the
     ! intention is to have different types of designs for different types of
     ! problems.

     !> A mapper to map the design into coefficients in the PDE.
     !! @todo It is currently assumed to be the Brinkman amplitude in the fluid
     !! equations. Mapping to multiple coeficients is currently not supported,
     !! ie, CHT.
     type(mapping_handler_t) :: mapping

     !> A mask indicating the optimization domain
     class(point_zone_t), pointer :: optimization_domain
     !> A logical if we're restricting the optimization domain
     logical :: has_mask
     !> SEM coefficients
     class(coef_t), pointer :: coef
     !> Average mass matrix
     real(kind=rp) :: avg_B

     ! TODO
     ! you also had logicals for convergence etc,
     ! I feel they should live in "problem"

   contains

     ! ----------------------------------------------------------------------- !
     ! Initializations

     !> Initialize the design
     generic, public :: init => init_from_json_sim, init_from_components
     !> Initialize the design from a JSON file
     procedure, pass(this), public :: init_from_json_sim => &
          brinkman_design_init_from_json_sim
     !> Initialize the design from components
     procedure, pass(this), public :: init_from_components => &
          brinkman_design_init_from_components

     !> Retrieve the design variables
     procedure, pass(this) :: get_values => brinkman_design_get_design
     !> Retrieve the sensitivity
     procedure, pass(this) :: get_sensitivity => brinkman_design_get_sensitivity

     !> Retrieve the x location of the design variables
     procedure, pass(this) :: design_get_x => brinkman_design_get_x
     !> Retrieve the x location of the i'th design variable
     procedure, pass(this) :: design_get_x_i => brinkman_design_get_x_i
     !> Retrieve the y location of the design variables
     procedure, pass(this) :: design_get_y => brinkman_design_get_y
     !> Retrieve the y location of the i'th design variable
     procedure, pass(this) :: design_get_y_i => brinkman_design_get_y_i
     !> Retrieve the z location of the design variables
     procedure, pass(this) :: design_get_z => brinkman_design_get_z
     !> Retrieve the z location of the i'th design variable
     procedure, pass(this) :: design_get_z_i => brinkman_design_get_z_i

     !> Update the design
     procedure, pass(this) :: update_design => brinkman_design_update_design

     !> map (this will include everything from mapping
     ! design_indicator -> filtering -> chi
     ! and ultimately handle mapping different coeficients!
     procedure, pass(this) :: map_forward => brinkman_design_map_forward
     !> this will contain chain rule for going backwards
     ! d_design_indicator <- d_filtering <- d_chi
     ! and ultimately handle mapping different coeficients!
     procedure, pass(this) :: map_backward => brinkman_design_map_backward
     !> Project a gradient field onto all GLL points, recovering a vector of
     !! directional derivatives
     procedure, pass(this) :: convert_to_directional_derivative => &
          brinkman_design_convert_to_directional_derivative
     !> Project a vector sensitivity before optimization.
     procedure, pass(this) :: project_sensitivity_vector => &
          brinkman_design_project_sensitivity_vector
     !> Project a matrix sensitivity before optimization.
     procedure, pass(this) :: project_sensitivity_matrix => &
          brinkman_design_project_sensitivity_matrix
     generic, public :: project_sensitivity => &
          project_sensitivity_vector, project_sensitivity_matrix
     ! TODO
     ! maybe it would have been smarter to have a "coeficient" type,
     ! which is just a scalar field and set of mappings going from
     ! design_indicator -> coeficient and their corresponding chain rules
     ! maybe also some information about what equation they live in...

     ! a writer being called from outside would be nice
     procedure, pass(this) :: write => brinkman_design_write
     !> Reset output counters after a restart.
     procedure, pass(this) :: set_output_counter => &
          brinkman_design_set_output_counter

     !> Destructor
     procedure, pass(this) :: free => brinkman_design_free
     ! TODO
     ! I'm not sure who owns the optimizer...
     ! but it would make sense to have it in here so you provide it
     ! with dF/d_design_indicator and it updates itself.
     ! procedure, pass(this) :: update => brinkman_design_update_design
  end type brinkman_design_t

contains

  !> Initialize the design from a JSON file
  subroutine brinkman_design_init_from_json_sim(this, parameters, simulation)
    class(brinkman_design_t), intent(inout) :: this
    type(json_file), intent(inout) :: parameters
    type(simulation_t), intent(inout) :: simulation
    type(json_file) :: json_subdict
    character(len=:), allocatable :: domain_name, domain_type, name
    character(len=:), allocatable :: output_format_str, output_precision_str
    logical :: dealias, verbose_design, verbose_sensitivity
    integer :: output_precision

    call json_get_or_default(parameters, 'name', name, 'Brinkman Design')
    call json_get_or_default(parameters, 'domain.type', domain_type, 'full')
    call json_get_or_default(parameters, 'dealias', dealias, .true.)
    call json_get_or_default(parameters, 'verbose_design', verbose_design, &
         .false.)
    call json_get_or_default(parameters, 'verbose_sensitivity', &
         verbose_sensitivity, .false.)
    call json_get_or_default(parameters, 'output_precision', &
         output_precision_str, 'sp')
    call json_get_or_default(parameters, 'output_format', &
         output_format_str, 'fld')

    select case (trim(output_precision_str))
    case ('sp', 'SP', 'single')
       output_precision = sp
    case ('dp', 'DP', 'double')
       output_precision = dp
    case default
       call neko_error('Invalid output_precision in design.brinkman. ' // &
            'Expected "sp" or "dp".')
    end select

    select case (trim(domain_type))
    case ('full')
       this%has_mask = .false.
    case ('point_zone')
       this%has_mask = .true.
       call json_get(parameters, 'domain.zone_name', domain_name)
       this%optimization_domain => &
            neko_point_zone_registry%get_point_zone(domain_name)

    case default
       call neko_error('brinkman design only supports point_zones for ' // &
            'optimization domain types')

    end select

    ! Initialize and inject into the simulation
    call this%init_from_components(name, simulation, dealias)

    ! Initialize the mapper
    associate(coef => simulation%neko_case%fluid%c_Xh, &
         gs => simulation%neko_case%fluid%gs_Xh)

      if ('mapping' .in. parameters) then
         call this%mapping%init_base(coef)
         call this%mapping%add(parameters, 'mapping')
         call this%mapping%init_output_fields(this%design_indicator, &
              this%brinkman_amplitude, this%sensitivity, verbose_design, &
              verbose_sensitivity, output_precision, output_format_str, &
              'design', 'sensitivity')
      end if

      if ('initial_distribution' .in. parameters) then
         call json_get(parameters, 'initial_distribution', json_subdict)
         call set_optimization_ic(this%design_indicator, coef, gs, &
              json_subdict)
      else
         call field_rzero(this%design_indicator)
      end if
    end associate

    ! Map to the Brinkman amplitude
    call this%map_forward()

  end subroutine brinkman_design_init_from_json_sim

  !> Free the design
  subroutine brinkman_design_free(this)
    class(brinkman_design_t), intent(inout) :: this

    call this%free_base()
    call this%mapping%free()
    nullify(this%brinkman_amplitude)
    nullify(this%design_indicator)
    nullify(this%sensitivity)

  end subroutine brinkman_design_free

  subroutine brinkman_design_init_from_components(this, name, simulation, &
       dealias)
    class(brinkman_design_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    type(simulation_t), intent(inout) :: simulation
    logical, intent(in) :: dealias
    integer :: n
    type(simple_brinkman_source_term_t) :: forward_brinkman, adjoint_brinkman

    associate(dof => simulation%neko_case%fluid%dm_Xh)

      call neko_registry%add_field(dof, "design_indicator", .true.)
      call neko_registry%add_field(dof, "brinkman_amplitude", .true.)
      call neko_registry%add_field(dof, "sensitivity", .true.)

    end associate

    this%design_indicator => neko_registry%get_field("design_indicator")
    this%brinkman_amplitude => neko_registry%get_field("brinkman_amplitude")
    this%sensitivity => neko_registry%get_field("sensitivity")

    this%coef => simulation%fluid%c_Xh

    ! TODO
    ! this is where we steal basically everything in
    ! brinkman_source_term regarding loading initial fields
    ! for now, make it a cylinder by hand
    this%design_indicator = 0.0_rp
    this%brinkman_amplitude = 0.0_rp
    this%design_indicator = 0.0_rp

    ! TODO
    ! Regarding masks and filters,
    ! I suppose there are two ways of thinking about it:
    ! 1) Mask first, then filter
    ! 2) filter first, then mask
    !
    ! Each one can be used to achieve different results, and when we do complex
    ! non-linear filter cascades the choice here has implications for exactly
    ! how we define "minimum size control".
    !
    ! On top of this, I've only ever used spatial convolution filters before,
    ! not PDE based filters.
    ! With these filters you need to think of how your domain is "padded" when
    ! you move the kernel over the boundary. Again, you can achieve different
    ! effects depending on the choice of padding.
    ! I don't really know if there's a similar implication for PDE based
    ! based filters, I bet there is. (We should ask Niels and Casper) I bet it
    ! means we need to consider the boundary of the mask as the boundary we
    ! enforce to be Nuemann, not the boundary of the computational domain.
    ! That sounds REALLY hard to implement...I hope it doesn't come down to
    ! that.
    !
    ! We can always change this decision later, but I'm going with (1)
    ! mask first, then filter.
    ! The reason being, one of the purposes of the filtering is to avoid sharp
    ! interfaces, if we filter first then mask there's a chance we have a sharp
    ! interface on the boundary of the optimization domain.
    ! if we mask first then filter, at least all the boundaries will be smooth.
    if (this%has_mask) then
       call mask_exterior_const(this%design_indicator, &
            this%optimization_domain, 0.0_rp)
    end if

    n = this%design_indicator%dof%size()
    call this%init_base(name, n)

    ! compute the average mass matrix
    this%avg_B = this%coef%volume / real(simulation%fluid%glb_unique_points)

    ! init the simple brinkman term for the forward problem
    call forward_brinkman%init_from_components( &
         simulation%fluid%f_x, &
         simulation%fluid%f_y, &
         simulation%fluid%f_z, &
         this%brinkman_amplitude, &
         simulation%fluid%u, &
         simulation%fluid%v, &
         simulation%fluid%w, &
         simulation%fluid%c_Xh, &
         simulation%adjoint_fluid%c_Xh_GL, &
         simulation%adjoint_fluid%GLL_to_GL, &
         dealias, simulation%adjoint_fluid%scratch_GL)
    ! append brinkman source term to the forward problem
    call simulation%fluid%source_term%add(forward_brinkman)

    ! init the simple brinkman term for the adjoint
    call adjoint_brinkman%init_from_components( &
         simulation%adjoint_fluid%f_adj_x, &
         simulation%adjoint_fluid%f_adj_y, &
         simulation%adjoint_fluid%f_adj_z, &
         this%brinkman_amplitude, &
         simulation%adjoint_fluid%u_adj, &
         simulation%adjoint_fluid%v_adj, &
         simulation%adjoint_fluid%w_adj, &
         simulation%adjoint_fluid%c_Xh, &
         simulation%adjoint_fluid%c_Xh_GL, &
         simulation%adjoint_fluid%GLL_to_GL, &
         dealias, simulation%adjoint_fluid%scratch_GL)
    ! append brinkman source term based on design

    select type (f => simulation%adjoint_fluid)
    type is (adjoint_fluid_pnpn_t)
       call f%source_term%add(adjoint_brinkman)
    class default
    end select

  end subroutine brinkman_design_init_from_components


  subroutine brinkman_design_map_forward(this)
    class(brinkman_design_t), intent(inout) :: this

    ! TODO, see previous todo about mask first, then mapping
    if (this%has_mask) then
       call mask_exterior_const(this%design_indicator, &
            this%optimization_domain, 0.0_rp)
    end if

    call this%mapping%apply_forward(this%brinkman_amplitude, &
         this%design_indicator)

  end subroutine brinkman_design_map_forward

  subroutine brinkman_design_get_design(this, values)
    class(brinkman_design_t), intent(in) :: this
    type(vector_t), intent(inout) :: values
    integer :: n

    n = this%size()
    if (n .ne. values%size()) then
       call neko_error('Get design: size mismatch')
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(values%x_d, this%design_indicator%x_d, n)
    else
       call copy(values%x, this%design_indicator%x, n)
    end if

  end subroutine brinkman_design_get_design

  subroutine brinkman_design_get_sensitivity(this, values)
    class(brinkman_design_t), intent(in) :: this
    type(vector_t), intent(inout) :: values
    integer :: n

    n = this%size()
    if (n .ne. values%size()) then
       call neko_error('Get sensitivity: size mismatch')
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(values%x_d, this%sensitivity%x_d, n)
    else
       call copy(values%x, this%sensitivity%x, n)
    end if

  end subroutine brinkman_design_get_sensitivity

  ! The interpretation of the sensitivity DF is a function that satisfies
  ! <DF,phi> = delta F(x;phi)
  ! To recover directional derivatives we should project the sensitivity onto
  ! the basis functions, equivelent to multiplying by the mass matrix.
  subroutine brinkman_design_convert_to_directional_derivative(this, vec_in)
    class(brinkman_design_t), intent(in) :: this
    type(vector_t), intent(inout) :: vec_in
    integer :: n

    n = this%size()
    if (n .ne. vec_in%size()) then
       call neko_error('convert_to_directional_derivative: size mismatch')
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col2(vec_in%x_d, this%coef%B_d, n)
    else
       call col2(vec_in%x, this%coef%B, n)
    end if

  end subroutine brinkman_design_convert_to_directional_derivative

  ! to aleviate lumpy designs due to non uniform spacing, we instead project
  ! onto a uniform mass matrix with average weight, as opposed to using the
  ! directional derivative
  subroutine brinkman_design_project_sensitivity_vector(this, sensitivity)
    class(brinkman_design_t), intent(inout) :: this
    type(vector_t), intent(inout) :: sensitivity

    call vector_cmult(sensitivity, this%avg_B)
  end subroutine brinkman_design_project_sensitivity_vector

  subroutine brinkman_design_project_sensitivity_matrix(this, sensitivity)
    class(brinkman_design_t), intent(inout) :: this
    type(matrix_t), intent(inout) :: sensitivity

    call matrix_cmult(sensitivity, this%avg_B)
  end subroutine brinkman_design_project_sensitivity_matrix

  subroutine brinkman_design_get_x(this, x)
    class(brinkman_design_t), intent(in) :: this
    type(vector_t), intent(inout) :: x
    integer :: n

    n = this%size()
    if (n .ne. x%size()) then
       call neko_error('Get x: size mismatch')
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(x%x_d, this%design_indicator%dof%x_d, n)
    else
       call copy(x%x, this%design_indicator%dof%x, n)
    end if

  end subroutine brinkman_design_get_x

  function brinkman_design_get_x_i(this, i) result(x_i)
    class(brinkman_design_t), intent(in) :: this
    integer, intent(in) :: i
    real(kind=rp) :: x_i
    integer :: n

    n = this%size()
    if (i .lt. 1 .or. i .gt. n) then
       call neko_error('brinkman_design_get_x_i: index out of bounds')
    end if

    x_i = this%design_indicator%dof%x(i,1,1,1)

  end function brinkman_design_get_x_i

  subroutine brinkman_design_get_y(this, y)
    class(brinkman_design_t), intent(in) :: this
    type(vector_t), intent(inout) :: y
    integer :: n

    n = this%size()
    if (n .ne. y%size()) then
       call neko_error('Get y: size mismatch')
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(y%x_d, this%design_indicator%dof%y_d, n)
    else
       call copy(y%x, this%design_indicator%dof%y, n)
    end if

  end subroutine brinkman_design_get_y

  function brinkman_design_get_y_i(this, i) result(y_i)
    class(brinkman_design_t), intent(in) :: this
    integer, intent(in) :: i
    real(kind=rp) :: y_i
    integer :: n

    n = this%size()
    if (i .lt. 1 .or. i .gt. n) then
       call neko_error('brinkman_design_get_y_i: index out of bounds')
    end if

    y_i = this%design_indicator%dof%y(i,1,1,1)

  end function brinkman_design_get_y_i

  subroutine brinkman_design_get_z(this, z)
    class(brinkman_design_t), intent(in) :: this
    type(vector_t), intent(inout) :: z
    integer :: n

    n = this%size()
    if (n .ne. z%size()) then
       call neko_error('Get z: size mismatch')
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(z%x_d, this%design_indicator%dof%z_d, n)
    else
       call copy(z%x, this%design_indicator%dof%z, n)
    end if

  end subroutine brinkman_design_get_z

  function brinkman_design_get_z_i(this, i) result(z_i)
    class(brinkman_design_t), intent(in) :: this
    integer, intent(in) :: i
    real(kind=rp) :: z_i
    integer :: n

    n = this%size()
    if (i .lt. 1 .or. i .gt. n) then
       call neko_error('brinkman_design_get_z_i: index out of bounds')
    end if

    z_i = this%design_indicator%dof%z(i,1,1,1)

  end function brinkman_design_get_z_i

  subroutine brinkman_design_update_design(this, values)
    class(brinkman_design_t), intent(inout) :: this
    type(vector_t), intent(inout) :: values
    integer :: n

    n = this%size()
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(this%design_indicator%x_d, values%x_d, n)
    else
       call copy(this%design_indicator%x, values%x, n)
    end if

    call this%map_forward()

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(values%x_d, this%design_indicator%x_d, n)
    else
       call copy(values%x, this%design_indicator%x, n)
    end if

  end subroutine brinkman_design_update_design

  subroutine brinkman_design_map_backward(this, sensitivity)
    class(brinkman_design_t), intent(inout) :: this
    type(vector_t), intent(in) :: sensitivity
    type(field_t), pointer :: tmp_fld
    integer :: temp_index

    call neko_scratch_registry%request(tmp_fld, temp_index, .false.)

    call vector_to_field(tmp_fld, sensitivity)

    call this%mapping%apply_backward(this%sensitivity, tmp_fld)

    if (this%has_mask) then
       call mask_exterior_const(this%sensitivity, this%optimization_domain, &
            0.0_rp)
    end if

    call neko_scratch_registry%relinquish_field(temp_index)

  end subroutine brinkman_design_map_backward

  subroutine brinkman_design_write(this, idx)
    class(brinkman_design_t), intent(inout) :: this
    integer, intent(in) :: idx

    call this%mapping%write_design(idx)
    call this%mapping%write_sensitivity(idx)

  end subroutine brinkman_design_write

  subroutine brinkman_design_set_output_counter(this, idx)
    class(brinkman_design_t), intent(inout) :: this
    integer, intent(in) :: idx

    call this%mapping%set_output_counter(idx)

  end subroutine brinkman_design_set_output_counter

end module brinkman_design
