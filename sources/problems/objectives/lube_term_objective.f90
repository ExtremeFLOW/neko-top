! Copyright (c) 2023, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!    *  Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!    *  Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!    *  Neither the name of the authors nor the names of its
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
!
!> Implements the `lube_term_objective_t` type.
!
! I promise I'll write this document properly in the future...
!
! But the Borval Peterson (I think) paper had an objective function
! that had 2 terms, dissipation and this term they claimed represented
! out of plane stresses.
! I never really understood that extra term, I also don't think it
! applies to 3D cases, but everyone includes it anyway.
!
! It appears to me to be basically a heuristic penality that targets
! non-binary designs
!
! so let's call
!
! F = \int |\nabla u|^2  + K \int \chi \u^2
!
!      | dissipation |     |"lube term"|
!
! I say "lube term" because they said it came from lubrication theory...
! Anyway, we can change all this later (especially the names!)

! If the objective function \int |\nabla u|^2,
! the corresponding adjoint forcing is \int \nabla v \cdot \nabla u
!
! for the lube term, the adjoint forcing is \chi u
!
! This has always annoyed me...
! because now I see one objective and one constraint
!
module lube_term_objective
  use num_types, only: rp
  use field, only: field_t
  use field_math, only: field_col3, field_addcol3, field_cmult, field_add2s2
  use operators, only: grad
  use scratch_registry, only: neko_scratch_registry
  use objective, only: objective_t
  use simulation, only: simulation_t
  use fluid_scheme_incompressible, only: fluid_scheme_incompressible_t
  use adjoint_scheme, only: adjoint_scheme_t
  use neko_config, only: NEKO_BCKND_DEVICE
  use math, only: glsc2, copy
  use device_math, only: device_copy
  use design, only: design_t
  use topopt_design, only: topopt_design_t
  use adjoint_lube_source_term, only: adjoint_lube_source_term_t
  use point_zone, only: point_zone_t
  use mask_ops, only: mask_exterior_const
  use math_ext, only: glsc2_mask
  use utils, only: neko_error
  use json_module, only: json_file
  use json_utils, only: json_get_or_default
  use field_registry, only: neko_field_registry
  use logger, only: neko_log
  implicit none
  private

  !> An objective function corresponding to minimum dissipation
  !! $ F =  \int_\Omega |\nabla u|^2 d \Omega + K \int_Omega \frac{1}{2} \chi
  !! |\mathbf{u}|^2 d \Omega $
  type, public, extends(objective_t) :: lube_term_objective_t
     private

     !> The coefficient for the lube term.
     real(kind=rp) :: K

     ! Internal references to the simulation components.

     type(field_t), pointer :: u, v, w
     real(kind=rp), pointer :: B(:,:,:,:)
     type(field_t), pointer :: brinkman_amplitude

   contains

     !> The common constructor using a JSON object.
     procedure, public, pass(this) :: init_json => lube_term_init_json
     !> The actual constructor.
     procedure, public, pass(this) :: init_from_attributes => &
          lube_term_init_attributes
     !> Destructor.
     procedure, public, pass(this) :: free => lube_term_free
     !> Computes the value of the objective function.
     procedure, public, pass(this) :: update_value => &
          lube_term_update_value
     !> Computes the sensitivity with respect to the coefficient $\chi$.
     procedure, public, pass(this) :: update_sensitivity => &
          lube_term_update_sensitivity

  end type lube_term_objective_t

contains

  !> The common constructor using a JSON object.
  !! @param design the design.
  !! @param fluid the fluid scheme.
  !! @param adjoint the fluid adjoint.
  subroutine lube_term_init_json(this, json, design, simulation)
    class(lube_term_objective_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(design_t), intent(in) :: design
    type(simulation_t), target, intent(inout) :: simulation
    character(len=:), allocatable :: mask_name
    character(len=:), allocatable :: name
    real(kind=rp) :: weight, K

    call json_get_or_default(json, "weight", weight, 1.0_rp)
    call json_get_or_default(json, "mask_name", mask_name, "")
    call json_get_or_default(json, "name", name, "Out of plane stresses")
    call json_get_or_default(json, "K", K, 1.0_rp)

    call this%init_from_attributes(design, simulation, weight, name, &
         mask_name, K)
  end subroutine lube_term_init_json

  !> The actual constructor.
  !! @param design the design.
  !! @param simulation the simulation.
  !! @param weight the weight of the objective function.
  !! @param mask_name the name of the mask.
  !! @param K the coefficient for the lube term.
  subroutine lube_term_init_attributes(this, design, simulation, weight, &
       name, mask_name, K)
    class(lube_term_objective_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(simulation_t), target, intent(inout) :: simulation
    real(kind=rp), intent(in) :: weight
    character(len=*), intent(in) :: mask_name
    character(len=*), intent(in) :: name
    real(kind=rp), intent(in) :: K
    type(adjoint_lube_source_term_t) :: lube_term

    ! Call the base initializer
    call this%init_base(name, design%size(), weight, mask_name)

    ! Set the coefficient for the lube term
    this%K = K

    ! Grab the brinkman amplitude for the lube term
    select type (design)
    type is (topopt_design_t)
       this%brinkman_amplitude => design%brinkman_amplitude

    class default
       call neko_error('Minimum dissipation only works with topopt_design')
    end select

    this%u => neko_field_registry%get_field('u')
    this%v => neko_field_registry%get_field('v')
    this%w => neko_field_registry%get_field('w')
    this%B => simulation%neko_case%fluid%c_Xh%B

    ! if we have the lube term we need to initialize and append that too

    associate(f_adj_x => simulation%adjoint_case%scheme%f_adj_x, &
         f_adj_y => simulation%adjoint_case%scheme%f_adj_y, &
         f_adj_z => simulation%adjoint_case%scheme%f_adj_z, &
         c_Xh => simulation%adjoint_case%scheme%c_Xh)

      call lube_term%init_from_components(f_adj_x, f_adj_y, f_adj_z, design, &
           this%k * this%weight, &
           this%u, this%v, this%w, &
           this%mask, this%has_mask, c_Xh)

    end associate

    ! append adjoint forcing term based on objective function
    call simulation%adjoint_case%scheme%source_term%add_source_term(lube_term)

  end subroutine lube_term_init_attributes

  !> Destructor.
  subroutine lube_term_free(this)
    class(lube_term_objective_t), intent(inout) :: this
    call this%free_base()

    this%u => null()
    this%v => null()
    this%w => null()
    this%B => null()
    this%brinkman_amplitude => null()

  end subroutine lube_term_free

  !> Compute the objective function.
  !! @param design the design.
  !! @param fluid the fluid scheme.
  !! @param adjoint the fluid adjoint.
  subroutine lube_term_update_value(this, design)
    class(lube_term_objective_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(field_t), pointer :: work
    integer :: temp_indices(1)

    call neko_scratch_registry%request_field(work, temp_indices(1))

    ! it's becoming so stupid to pass the whole fluid and adjoint and
    ! design through
    ! I feel like every objective function should have internal pointers to
    ! u,v,w and u_adj, v_adj, w_adj and perhaps the design
    ! (the whole design, so we get all the coeffients)
    call field_col3(work, this%u, this%brinkman_amplitude)
    call field_addcol3(work, this%v, this%brinkman_amplitude)
    call field_addcol3(work, this%w, this%brinkman_amplitude)

    if (this%has_mask) then
       this%value = glsc2_mask(work%x, this%B, design%size(), &
            this%mask%mask, this%mask%size)
    else
       this%value = glsc2(work%x, this%B, design%size())
    end if
    this%value = 0.5 * this%K * this%value

    !TODO
    ! GPUS

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine lube_term_update_value

  !> update_value the sensitivity of the objective function with respect to $\chi$
  !! @param design the design.
  subroutine lube_term_update_sensitivity(this, design)
    class(lube_term_objective_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(field_t), pointer :: work
    integer :: temp_indices(1)

    ! if we have the lube term we also get an extra term in the sensitivity
    ! K * u^2
    ! TODO
    ! omfg be so careful with non-dimensionalization etc
    ! I bet this is scaled a smidge wrong (ie, track if it's 1/2 or not etc)
    ! do this later

    call neko_scratch_registry%request_field(work, temp_indices(1))

    call field_col3(work, this%u, this%u)
    call field_addcol3(work, this%v, this%v)
    call field_addcol3(work, this%w, this%w)
    call field_cmult(work, this%K)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(this%sensitivity%x_d, work%x_d, this%sensitivity%size())
    else
       call copy(this%sensitivity%x, work%x, this%sensitivity%size())
    end if

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine lube_term_update_sensitivity

end module lube_term_objective
