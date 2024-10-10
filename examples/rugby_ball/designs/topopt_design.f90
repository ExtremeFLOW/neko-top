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
  use num_types, only: rp, dp
  use json_module, only: json_file
  use json_utils, only: json_get, json_get_or_default
  use simulation_component, only: simulation_component_t
  use case, only: case_t
  use field, only: field_t
  use coefs, only: coef_t
  use field_registry, only: neko_field_registry
  use scratch_registry, only: neko_scratch_registry
  use adjoint_pnpn, only: adjoint_pnpn_t
  use adjoint_output, only: adjoint_output_t
  use neko_config, only: NEKO_BCKND_DEVICE
  use field_math, only: field_cfill, field_sub2, field_copy, field_glsc2, &
       field_glsc3
  use field_math, only: field_add2
  use math, only: glsc2, glsc3
  use device_math, only: device_glsc2
  use adv_lin_no_dealias, only: adv_lin_no_dealias_t
  use logger, only: neko_log, LOG_SIZE
  use adjoint_scheme, only: adjoint_scheme_t
  use adjoint_fctry, only: adjoint_scheme_factory
  use time_step_controller, only: time_step_controller_t
  use time_scheme_controller, only: time_scheme_controller_t
  use mpi_f08, only: MPI_WTIME
  use jobctrl, only: jobctrl_time_limit
  use profiler, only: profiler_start, profiler_stop, profiler_start_region, &
       profiler_end_region
  use file, only: file_t
  use num_types, only : rp, sp, dp
  use fluid_scheme, only : fluid_scheme_factory
  use fluid_pnpn, only : fluid_pnpn_t
  use fluid_scheme, only : fluid_scheme_t
  use fluid_output, only : fluid_output_t
  use chkp_output, only : chkp_output_t
  use mean_sqr_flow_output, only : mean_sqr_flow_output_t
  use mean_flow_output, only : mean_flow_output_t
  use fluid_stats_output, only : fluid_stats_output_t
  use mpi_f08
  use mesh_field, only : mesh_fld_t, mesh_field_init, mesh_field_free
  use parmetis, only : parmetis_partmeshkway
  use redist, only : redist_mesh
  use sampler, only : sampler_t
  use flow_ic, only : set_flow_ic
  use scalar_ic, only : set_scalar_ic
  use field, only : field_t
  use field_registry, only : neko_field_registry
  use stats, only : stats_t
  use file, only : file_t
  use utils, only : neko_error
  use mesh, only : mesh_t
  use comm
  use time_scheme_controller, only : time_scheme_controller_t
  use logger, only : neko_log, NEKO_LOG_QUIET, LOG_SIZE
  use jobctrl, only : jobctrl_set_time_limit
  use user_intf, only : user_t
  use scalar_pnpn, only : scalar_pnpn_t
  use json_module, only : json_file, json_core, json_value
  use json_utils, only : json_get, json_get_or_default
  use scratch_registry, only : scratch_registry_t, neko_scratch_registry
  use point_zone_registry, only: neko_point_zone_registry
  use adjoint_ic, only : set_adjoint_ic
  use json_utils, only : json_extract_item
  use json_utils_ext, only: json_key_fallback, json_get_subdict
  use dofmap, only : dofmap_t
  use filters, only: permeability_field
  use mma, only: mma_t
  use fld_file_output, only : fld_file_output_t
  use linear_mapping, only: linear_mapping_t
  use RAMP_mapping, only: RAMP_mapping_t
  use Casper_mapping, only: Casper_mapping_t
  use PDE_filter, only: PDE_filter_t
  !use design_variable, only: design_variable_t
  implicit none
  private

  !> A topology optimization design variable
  type, public :: topopt_design_t

     ! TODO
     ! in the future make this a derived type of a `design_variable`
     ! type, public, extends(design_variable_t) :: topopt_design_t
     !> the unfilitered design
     type(field_t), public :: design_indicator

     !> the mapped coefficient (Brinkman term)
     ! TODO
     ! NOTE: Tim, right now we only have the brinkman term
     ! in the future we may map to other coeeficients for other equations...
     ! in which case this should be a field list
     !
     ! or as I describe below, we also have multiple constraints,
     ! so a list-of-lists may be the correct way forward
     type(field_t), public :: brinkman_amplitude

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
     type(field_t), public :: sensitivity
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
     type(Casper_mapping_t) :: mapping
     ! and we need to hold onto a field for the chain of mappings
     type(field_t) :: filtered_design


     !
     ! TODO
     ! you also had lots of masks etc that was a nice idea,
     ! but we'll cross that bridge later

     ! TODO
     ! you also had logicals for convergence etc,
     ! I feel they should live in "problem"

     ! afield writer would be nice too
     type(fld_file_output_t), private :: output


   contains
     !> init (will make this legit at some point)
     procedure, pass(this) :: init => topopt_design_init
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

     ! a sampler being called from outside would be nice
     procedure, pass(this) :: sample => topopt_design_sample

     !> Destructor
     procedure, pass(this) :: free => topopt_design_free
     ! TODO
     ! I'm not sure who owns the optimizer...
     ! but it would make sense to have it in here so you provide it
     ! with dF/d_design_indicator and it updates itself.
     ! procedure, pass(this) :: update => topopt_design_update
  end type topopt_design_t


contains

  subroutine topopt_design_init(this, json, coef)
    class(topopt_design_t), target, intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(coef_t), intent(inout) :: coef
    integer :: n, i
    ! init the fields
    call this%design_indicator%init(coef%dof, "design_indicator")
    call this%brinkman_amplitude%init(coef%dof, "brinkman_amplitude")
    call this%sensitivity%init(coef%dof, "sensitivity")
    call this%filtered_design%init(coef%dof, "filtered_design")

    ! TODO
    ! this is where we steal basically everything in
    ! brinkman_source_term regarding loading initial fields
    ! for now, make it a cylinder by hand
    this%design_indicator = 0.0_rp
    this%brinkman_amplitude = 0.0_rp

    n = this%design_indicator%dof%size()
    do i = 1, n
       if(sqrt((this%design_indicator%dof%x(i,1,1,1) - 0.5_rp)**2 + &
            (this%design_indicator%dof%y(i,1,1,1) &
            - 0.5_rp)**2).lt.0.25_rp) then
          this%design_indicator%x(i,1,1,1) = 1.0_rp
       endif
    enddo

    ! TODO
    ! we would also need to make a mapping type that reads in
    ! parameters etc about filtering and mapping
    ! ie,
    ! call mapper%init(this woud be from JSON)
    call this%filter%init(json, coef)
    call this%mapping%init(json, coef)

    ! and then we would map for the first one
    call this%map_forward()


    ! a field writer would be nice to output
    ! - design indicator (\rho)
    ! - mapped design (\chi)
    ! - sensitivity (dF/d\chi)
    ! TODO
    ! do this properly with JSON
    ! TODO
    ! obviously when we do the mappings properly, to many coeficients,
    ! we'll also have to modify this
    call this%output%init(sp,'design',3)
    call this%output%fields%assign(1, this%design_indicator)
    call this%output%fields%assign(2, this%brinkman_amplitude)
    call this%output%fields%assign(3, this%sensitivity)


  end subroutine topopt_design_init


  subroutine topopt_design_map_forward(this)
    class(topopt_design_t), target, intent(inout) :: this


    ! TODO
    ! this should be somehow deffered so we can pick different mappings!!!
    ! so this would be:
    ! call mapper%forward(fld_out, fld_in)

    call this%filter%apply_forward(this%filtered_design, &
    this%design_indicator)

    call this%mapping%apply_forward(this%brinkman_amplitude, &
    this%filtered_design)


  end subroutine topopt_design_map_forward

  subroutine topopt_design_map_backward(this, df_dchi)
    class(topopt_design_t), target, intent(inout) :: this
    type(field_t), intent(in) :: df_dchi
    type(field_t), pointer :: dF_dfiltered_design
    integer :: temp_indices(1)

    ! TODO
    ! again..
    ! so this would be:
    ! call mapper%backward(fld_out, fld_in)
    call neko_scratch_registry%request_field(dF_dfiltered_design, &
    temp_indices(1))

    call this%mapping%apply_backward(dF_dfiltered_design, df_dchi, &
    this%filtered_design)

    call this%filter%apply_backward(this%sensitivity, dF_dfiltered_design, &
    this%filtered_design)

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine topopt_design_map_backward

  subroutine topopt_design_free(this)
    class(topopt_design_t), target, intent(inout) :: this
    call this%brinkman_amplitude%free()
    call this%design_indicator%free()

  end subroutine topopt_design_free

  subroutine topopt_design_sample(this,t)
    class(topopt_design_t), target, intent(inout) :: this
    real(kind=rp), intent(in) :: t

    call this%output%sample(t)

  end subroutine topopt_design_sample
end module topopt_design
