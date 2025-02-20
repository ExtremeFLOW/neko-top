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

! Implements the `topopt_design_t` type.
module topopt_design
  use num_types, only: rp, sp
  use field, only: field_t
  use json_module, only: json_file
  use mapping, only: mapping_t
  use PDE_filter, only: PDE_filter_t
  use RAMP_mapping, only: RAMP_mapping_t
  use coefs, only: coef_t
  use scratch_registry, only: neko_scratch_registry
  use fld_file_output, only: fld_file_output_t
  use point_zone_registry, only: neko_point_zone_registry
  use point_zone, only: point_zone_t
  use mask_ops, only: mask_exterior_const
  use neko_config, only: NEKO_BCKND_DEVICE
  use device, only: device_memcpy, HOST_TO_DEVICE
  use design, only: design_t
  use math, only: rzero
  use simulation, only: simulation_t
  use json_module, only: json_file
  use simple_brinkman_source_term, only: simple_brinkman_source_term_t
  use vector, only: vector_t
  use math, only: copy
  use field_registry, only: neko_field_registry
  implicit none
  private

  !> A topology optimization design variable
  type, extends(design_t), public :: topopt_design_t
     private

     ! TODO
     ! in the future make this a derived type of a `design_variable`
     ! type, public, extends(design_variable_t) :: topopt_design_t
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
     type(PDE_filter_t) :: filter
     type(RAMP_mapping_t) :: mapping

     ! and we need to hold onto a field for the chain of mappings
     type(field_t), pointer :: filtered_design


     !> A mask indicating the optimization domain
     class(point_zone_t), pointer :: optimization_domain
     !> A logical if we're restricting the optimization domain
     logical :: if_mask

     ! TODO
     ! you also had logicals for convergence etc,
     ! I feel they should live in "problem"

     ! afield writer would be nice too
     type(fld_file_output_t), private :: output

   contains

     ! ----------------------------------------------------------------------- !
     ! Initializations

     !> Initialize the design
     generic, public :: init => init_from_json, &
          init_from_components
     !> Initialize the design from a JSON file
     procedure, pass(this), public :: init_from_json => &
          topopt_design_init_from_json
     !> Initialize the design from components
     procedure, pass(this), public :: init_from_components => &
          topopt_design_init_from_components

     !> Retrieve the design variables
     procedure, pass(this) :: get_design => topopt_design_get_design

     !> Update the design
     procedure, pass(this) :: update_design => topopt_design_update_design

     !> map (this will include everything from mapping
     ! design_indicator -> filtering -> chi
     ! and ultimately handle mapping different coeficients!
     procedure, pass(this) :: map_forward => topopt_design_map_forward
     !> this will contain chain rule for going backwards
     ! d_design_indicator <- d_filtering <- d_chi
     ! and ultimately handle mapping different coeficients!
     procedure, pass(this) :: map_backward => topopt_design_map_backward
     ! TODO
     ! maybe it would have been smarter to have a "coeficient" type,
     ! which is just a scalar field and set of mappings going from
     ! design_indicator -> coeficient and their corresponding chain rules
     ! maybe also some information about what equation they live in...

     ! a writer being called from outside would be nice
     procedure, pass(this) :: write => topopt_design_write

     !> Destructor
     procedure, pass(this) :: free => topopt_design_free
     ! TODO
     ! I'm not sure who owns the optimizer...
     ! but it would make sense to have it in here so you provide it
     ! with dF/d_design_indicator and it updates itself.
     ! procedure, pass(this) :: update => topopt_design_update_design
  end type topopt_design_t


contains

  !> Initialize the design from a JSON file
  subroutine topopt_design_init_from_json(this, parameters, simulation)
    class(topopt_design_t), intent(inout) :: this
    type(json_file), intent(inout) :: parameters
    type(simulation_t), intent(inout) :: simulation

    call this%init_from_components(simulation)

    ! Todo: This need to be read from the parameters in the JSON
    associate(coef => simulation%neko_case%fluid%c_Xh)
      call this%filter%init(parameters, coef)
      call this%mapping%init(parameters, coef)
    end associate

    ! and then we would map for the first one
    call this%map_forward()

  end subroutine topopt_design_init_from_json

  !> Free the design
  subroutine topopt_design_free(this)
    class(topopt_design_t), intent(inout) :: this

    call this%free_base()
    call this%brinkman_amplitude%free()
    call this%design_indicator%free()
    call this%filtered_design%free()
    call this%sensitivity%free()

  end subroutine topopt_design_free

  subroutine topopt_design_init_from_components(this, simulation)
    class(topopt_design_t), intent(inout) :: this
    type(simulation_t), intent(inout) :: simulation
    character(len=:), allocatable :: optimization_domain_zone_name
    integer :: n, i
    type(simple_brinkman_source_term_t) :: forward_brinkman, adjoint_brinkman

    associate(dof => simulation%neko_case%fluid%dm_Xh)

      call neko_field_registry%add_field(dof, "design_indicator", .true.)
      call neko_field_registry%add_field(dof, "brinkman_amplitude", .true.)
      call neko_field_registry%add_field(dof, "sensitivity", .true.)
      call neko_field_registry%add_field(dof, "filtered_design", .true.)

    end associate

    this%design_indicator => &
         neko_field_registry%get_field("design_indicator")
    this%brinkman_amplitude => &
         neko_field_registry%get_field("brinkman_amplitude")
    this%sensitivity => &
         neko_field_registry%get_field("sensitivity")
    this%filtered_design => &
         neko_field_registry%get_field("filtered_design")

    ! TODO
    ! this is where we steal basically everything in
    ! brinkman_source_term regarding loading initial fields
    ! for now, make it a cylinder by hand
    this%design_indicator = 0.0_rp
    this%brinkman_amplitude = 0.0_rp
    this%design_indicator%x = 0.0_rp

    n = this%design_indicator%dof%size()
    do i = 1, n
       if (sqrt((this%design_indicator%dof%x(i,1,1,1) - 0.5_rp)**2 + &
            (this%design_indicator%dof%y(i,1,1,1) &
            - 0.5_rp)**2) .lt. 0.25_rp) then
          this%design_indicator%x(i,1,1,1) = 1.0_rp
       end if
    end do

    ! again this will be handled better in the future...
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%design_indicator%x, &
            this%design_indicator%x_d, n, &
            HOST_TO_DEVICE, sync = .false.)
    end if

    ! TODO, of course when we move all of Tim's stuff for initialization of
    ! the initial design field we'll be reading things properly from the JSON.
    ! call json_get(parameters, 'name', optimization_domain_zone_name)
    ! Right now, I'm hardcoding the name of the point zone.
    this%if_mask = .true.
    optimization_domain_zone_name = "optimization_domain"

    ! Initialize the mask
    if (this%if_mask) then
       this%optimization_domain => &
            neko_point_zone_registry%get_point_zone(&
            optimization_domain_zone_name)
    end if

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
    if (this%if_mask) then
       call mask_exterior_const(this%design_indicator, &
            this%optimization_domain, 0.0_rp)
    end if

    ! a field writer would be nice to output
    ! - design indicator (\rho)
    ! - mapped design (\chi)
    ! - sensitivity (dF/d\chi)
    ! TODO
    ! do this properly with JSON
    ! TODO
    ! obviously when we do the mappings properly, to many coeficients,
    ! we'll also have to modify this
    call this%output%init(sp, 'design', 3)
    call this%output%fields%assign_to_field(1, this%design_indicator)
    call this%output%fields%assign_to_field(2, this%brinkman_amplitude)
    call this%output%fields%assign_to_field(3, this%sensitivity)

    call this%init_base(n)

    ! init the simple brinkman term for the forward problem
    call forward_brinkman%init_from_components( &
         simulation%fluid_scheme%f_x, &
         simulation%fluid_scheme%f_y, &
         simulation%fluid_scheme%f_z, &
         this%brinkman_amplitude, &
         simulation%fluid_scheme%u, &
         simulation%fluid_scheme%v, &
         simulation%fluid_scheme%w, &
         simulation%fluid_scheme%c_Xh)
    ! append brinkman source term to the forward problem
    call simulation%fluid_scheme%source_term%add(forward_brinkman)

    ! init the simple brinkman term for the adjoint
    call adjoint_brinkman%init_from_components( &
         simulation%adjoint_case%scheme%f_adj_x, &
         simulation%adjoint_case%scheme%f_adj_y, &
         simulation%adjoint_case%scheme%f_adj_z, &
         this%brinkman_amplitude, &
         simulation%adjoint_case%scheme%u_adj, &
         simulation%adjoint_case%scheme%v_adj, &
         simulation%adjoint_case%scheme%w_adj, &
         simulation%adjoint_case%scheme%c_Xh)
    ! append brinkman source term based on design
    call simulation%adjoint_case%scheme%source_term%add(adjoint_brinkman)

  end subroutine topopt_design_init_from_components


  subroutine topopt_design_map_forward(this)
    class(topopt_design_t), intent(inout) :: this

    ! TODO, see previous todo about mask first, then mapping
    if (this%if_mask) then
       call mask_exterior_const(this%design_indicator, &
            this%optimization_domain, 0.0_rp)
    end if

    ! TODO
    ! this should be somehow deffered so we can pick different mappings!!!
    ! so this would be:
    ! call mapper%forward(fld_out, fld_in)

    call this%filter%apply_forward(this%filtered_design, &
         this%design_indicator)

    call this%mapping%apply_forward(this%brinkman_amplitude, &
         this%filtered_design)

  end subroutine topopt_design_map_forward

  function topopt_design_get_design(this) result(x)
    class(topopt_design_t), intent(in) :: this
    type(vector_t) :: x
    integer :: n

    n = this%size()
    call x%init(n)
    call copy(x%x, this%design_indicator%x, n)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(x%x_d, this%design_indicator%x_d, n)
    end if

  end function topopt_design_get_design

  subroutine topopt_design_update_design(this, x)
    class(topopt_design_t), intent(inout) :: this
    type(vector_t), intent(inout) :: x
    integer :: n

    n = this%size()
    call copy(this%design_indicator%x, x%x, n)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(this%design_indicator%x_d, x%x_d, n)
    end if

    call this%map_forward()

    call copy(x%x, this%design_indicator%x, n)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(x%x_d, this%design_indicator%x_d, n)
    end if

  end subroutine topopt_design_update_design

  subroutine topopt_design_map_backward(this, sensitivity)
    class(topopt_design_t), intent(inout) :: this
    type(vector_t), intent(in) :: sensitivity
    type(field_t), pointer :: df_dchi
    type(field_t), pointer :: dF_dfiltered_design
    integer :: temp_indices(2)

    ! it would be nice to visualize this

    call neko_scratch_registry%request_field(df_dchi, temp_indices(1))
    call copy(df_dchi%x, sensitivity%x, this%size())

    ! TODO
    ! again..
    ! so this would be:
    ! call mapper%backward(fld_out, fld_in)
    call neko_scratch_registry%request_field(dF_dfiltered_design, &
         temp_indices(2))

    call this%mapping%apply_backward(dF_dfiltered_design, df_dchi, &
         this%filtered_design)

    call this%filter%apply_backward(this%sensitivity, dF_dfiltered_design, &
         this%filtered_design)

    ! TODO
    ! DELETE THIS LATER
    !
    ! When Abbas writes the interface for the optimization
    ! module this may be a moot point, because we would only really collect
    ! the sensitivity of the design variables inside the mask.
    !
    ! Note for Abbas,
    ! I'm NOT doing this because I'm too lazy and I just need masks so I can
    ! test something in the passive scalar.
    if (this%if_mask) then
       call mask_exterior_const(this%sensitivity, this%optimization_domain, &
            0.0_rp)
    end if

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine topopt_design_map_backward

  subroutine topopt_design_write(this, idx)
    class(topopt_design_t), intent(inout) :: this
    integer, intent(in) :: idx

    call this%output%sample(real(idx, kind=rp))

  end subroutine topopt_design_write

end module topopt_design
