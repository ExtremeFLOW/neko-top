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

 ! Implements the `adjoint_t` type.
module simcomp_example
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
  use field_math, only: field_cfill, field_sub2, field_copy, field_glsc2, field_glsc3
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
  use fluid_fctry, only : fluid_scheme_factory
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
  use material_properties, only : material_properties_t
  implicit none
  private

  ! An empty user defined simulation component.
  ! This is a simple example of a user-defined simulation component.
  type, public, extends(simulation_component_t) :: adjoint_t

     class(adjoint_scheme_t), allocatable :: scheme

     ! Fields
     type(field_t) :: u_old, v_old, w_old, p_old, s_old
     type(field_t), pointer :: u_adj, v_adj, w_adj, p_adj, s_adj
     real(kind=rp) :: tol
     type(adjoint_output_t) :: f_out
     type(sampler_t) :: s

     logical :: have_scalar = .false.
     logical :: converged = .false.
     logical :: computed = .false.

   contains
     ! Constructor from json, wrapping the actual constructor.
     procedure, pass(this) :: init => simcomp_test_init_from_json
     ! Actual constructor.
     procedure, pass(this) :: init_from_attributes => &
          simcomp_test_init_from_attributes
     ! Destructor.
     procedure, pass(this) :: free => simcomp_test_free
     ! Compute the simcomp_test field.
     procedure, pass(this) :: compute_ => simcomp_test_compute
  end type adjoint_t

contains

  ! Constructor from json.
  subroutine simcomp_test_init_from_json(this, json, case)
    class(adjoint_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case

    call this%init_from_attributes()
    call this%init_base(json, case)

    call adjoint_case_init_common(this, case)

    ! Read the tolerance
    call json_get_or_default(json, "tol", this%tol, 1.0e-6_rp)

    ! Point local fields to the scratch fields
    this%u_old = case%fluid%u
    this%v_old = case%fluid%v
    this%w_old = case%fluid%w
    this%p_old = case%fluid%p

    call field_cfill(this%u_old, 0.0_rp)
    call field_cfill(this%v_old, 0.0_rp)
    call field_cfill(this%w_old, 0.0_rp)
    call field_cfill(this%p_old, 0.0_rp)

    ! Check if the scalar field is allocated
    if (allocated(case%scalar)) then
       this%have_scalar = .true.
       this%s_old = case%scalar%s
       call field_cfill(this%s_old, 0.0_rp)
    end if

    ! Point to the adjoint fields
    this%u_adj => this%scheme%u_adj
    this%v_adj => this%scheme%v_adj
    this%w_adj => this%scheme%w_adj
    this%p_adj => this%scheme%p_adj

  end subroutine simcomp_test_init_from_json

  !> Initialize a case from its (loaded) params object
  subroutine adjoint_case_init_common(this, C)
    class(adjoint_t), intent(inout) :: this
    type(case_t), target, intent(inout) :: C
    character(len=:), allocatable :: output_directory
    integer :: lx = 0
    logical :: scalar = .false.
    type(file_t) :: msh_file, bdry_file, part_file
    logical :: found, logical_val
    integer :: integer_val
    real(kind=rp) :: real_val
    character(len=:), allocatable :: string_val
    real(kind=rp) :: stats_start_time, stats_output_val
    integer :: stats_sampling_interval
    integer :: output_dir_len
    integer :: precision
    type(field_t), pointer :: u_b, v_b, w_b

    !
    ! Setup fluid scheme
    !
    call json_get(C%params, 'case.fluid.scheme', string_val)
    call adjoint_scheme_factory(this%scheme, trim(string_val))

    call json_get(C%params, 'case.numerics.polynomial_order', lx)
    lx = lx + 1 ! add 1 to get number of gll points
    call this%scheme%init(C%msh, lx, C%params, C%usr, C%material_properties)
    ! this%scheme%chkp%tlag => C%tlag
    ! this%scheme%chkp%dtlag => C%dtlag
    select type (f => this%scheme)
      type is (adjoint_pnpn_t)
       !  f%chkp%abx1 => f%abx1
       !  f%chkp%abx2 => f%abx2
       !  f%chkp%aby1 => f%aby1
       !  f%chkp%aby2 => f%aby2
       !  f%chkp%abz1 => f%abz1
       !  f%chkp%abz2 => f%abz2
    end select

    ! !
    ! ! Setup scalar scheme
    ! !
    ! ! @todo no scalar factory for now, probably not needed
    ! if (C%params%valid_path('case.scalar')) then
    !    call json_get_or_default(C%params, 'case.scalar.enabled', scalar,&
    !                             .true.)
    ! end if

    ! if (scalar) then
    !    allocate(C%scalar)
    !    call C%scalar%init(C%msh, this%scheme%c_Xh, this%scheme%gs_Xh, C%params, C%usr,&
    !                       C%material_properties)
    !    call this%scheme%chkp%add_scalar(C%scalar%s)
    !    this%scheme%chkp%abs1 => C%scalar%abx1
    !    this%scheme%chkp%abs2 => C%scalar%abx2
    !    this%scheme%chkp%slag => C%scalar%slag
    ! end if

    !
    ! Setup user defined conditions
    !
    if (C%params%valid_path('case.fluid.inflow_condition')) then
       call json_get(C%params, 'case.fluid.inflow_condition.type',&
            string_val)
       if (trim(string_val) .eq. 'user') then
          call this%scheme%set_usr_inflow(C%usr%fluid_user_if)
       end if
    end if

    ! ! Setup user boundary conditions for the scalar.
    ! if (scalar) then
    !    call C%scalar%set_user_bc(C%usr%scalar_user_bc)
    ! end if

    !
    ! Setup initial conditions
    !
    call json_get(C%params, 'case.fluid.initial_condition.type',&
         string_val)
    if (trim(string_val) .ne. 'user') then
       call set_flow_ic(this%scheme%u_adj, this%scheme%v_adj, this%scheme%w_adj, this%scheme%p_adj, &
            this%scheme%c_Xh, this%scheme%gs_Xh, string_val, C%params)
    else
       call set_flow_ic(this%scheme%u_adj, this%scheme%v_adj, this%scheme%w_adj, this%scheme%p_adj, &
            this%scheme%c_Xh, this%scheme%gs_Xh, C%usr%fluid_user_ic, C%params)
    end if

    ! if (scalar) then
    !    call json_get(C%params, 'case.scalar.initial_condition.type', string_val)
    !    if (trim(string_val) .ne. 'user') then
    !       call set_scalar_ic(C%scalar%s, &
    !         C%scalar%c_Xh, C%scalar%gs_Xh, string_val, C%params)
    !    else
    !       call set_scalar_ic(C%scalar%s, &
    !         C%scalar%c_Xh, C%scalar%gs_Xh, C%usr%scalar_user_ic, C%params)
    !    end if
    ! end if

    ! Add initial conditions to BDF scheme (if present)
    select type (f => this%scheme)
      type is (adjoint_pnpn_t)
       call f%ulag%set(f%u_adj)
       call f%vlag%set(f%v_adj)
       call f%wlag%set(f%w_adj)

		! baseflow is solution to forward problem
       u_b => neko_field_registry%get_field('u')
       v_b => neko_field_registry%get_field('v')
       w_b => neko_field_registry%get_field('w')

       !!
       !! Setup initial baseflow
       !!
       !call json_get(C%params, 'case.fluid.baseflow.type', string_val)

       !if (trim(string_val) .ne. 'user') then
       !   call set_baseflow(u_b, v_b, w_b, this%scheme%c_Xh, this%scheme%gs_Xh, &
       !        string_val, C%params)
       !else
       !   call set_baseflow(u_b, v_b, w_b, this%scheme%c_Xh, this%scheme%gs_Xh, &
       !        C%usr%baseflow_user, C%params)
       !end if



		! Tim what is this for?
       call field_cfill(f%u_b, 0.0_rp)
       call field_cfill(f%v_b, 0.0_rp)
       call field_cfill(f%w_b, 0.0_rp)


    end select

    !
    ! Validate that the case is properly setup for time-stepping
    !
    call this%scheme%validate

    ! if (scalar) then
    !    call C%scalar%slag%set(C%scalar%s)
    !    call C%scalar%validate
    ! end if

    !
    ! Setup output precision of the field files
    !
    call json_get_or_default(C%params, 'case.output_precision', string_val,&
         'single')

    if (trim(string_val) .eq. 'double') then
       precision = dp
    else
       precision = sp
    end if

    !
    ! Setup sampler
    !
    call this%s%init(C%end_time)
    if (scalar) then
       this%f_out = adjoint_output_t(precision, this%scheme, C%scalar, &
            path = trim(output_directory))
    else
       this%f_out = adjoint_output_t(precision, this%scheme, &
            path = trim(output_directory))
    end if

    call json_get_or_default(C%params, 'case.fluid.output_control',&
         string_val, 'org')

    if (trim(string_val) .eq. 'org') then
       ! yes, it should be real_val below for type compatibility
       call json_get(C%params, 'case.nsamples', real_val)
       call this%s%add(this%f_out, real_val, 'nsamples')
    else if (trim(string_val) .eq. 'never') then
       ! Fix a dummy 0.0 output_value
       call json_get_or_default(C%params, 'case.fluid.output_value', real_val, &
            0.0_rp)
       call this%s%add(this%f_out, 0.0_rp, string_val)
    else
       call json_get(C%params, 'case.fluid.output_value', real_val)
       call this%s%add(this%f_out, real_val, string_val)
    end if

    ! !
    ! ! Save checkpoints (if nothing specified, default to saving at end of sim)
    ! !
    ! call json_get_or_default(C%params, 'case.output_checkpoints',&
    !      logical_val, .true.)
    ! if (logical_val) then
    !    call json_get_or_default(C%params, 'case.checkpoint_format', &
    !         string_val, "chkp")
    !   !   C%f_chkp = chkp_output_t(this%scheme%chkp, path = output_directory, &
    !         ! fmt = trim(string_val))
    !    call json_get_or_default(C%params, 'case.checkpoint_control', &
    !         string_val, "simulationtime")
    !    call json_get_or_default(C%params, 'case.checkpoint_value', real_val,&
    !         1e10_rp)
    !   !  call this%s%add(C%f_chkp, real_val, string_val)
    ! end if

  end subroutine adjoint_case_init_common

  ! Actual constructor.
  subroutine simcomp_test_init_from_attributes(this)
    class(adjoint_t), intent(inout) :: this
  end subroutine simcomp_test_init_from_attributes

  ! Destructor.
  subroutine simcomp_test_free(this)
    class(adjoint_t), intent(inout) :: this

    call this%u_old%free()
    call this%v_old%free()
    call this%w_old%free()
    call this%p_old%free()
    if (this%have_scalar) then
       call this%s_old%free()
    end if

    call this%scheme%free()
    call this%free_base()

  end subroutine simcomp_test_free

  ! Compute the simcomp_test field.
  subroutine simcomp_test_compute(this, t, tstep)
    class(adjoint_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    real(kind=rp), dimension(5) :: normed_diff
    type(field_t), pointer :: u, v, w, p, s
    character(len=256) :: msg

    real(kind=rp) :: t_adj, cfl
    real(kind=dp) :: start_time_org, start_time, end_time
    character(len=LOG_SIZE) :: log_buf
    integer :: tstep_adj
    character(len=:), allocatable :: restart_file
    logical :: output_at_end, found
    ! for variable_tsteping
    real(kind=rp) :: cfl_avrg = 0.0_rp
    type(time_step_controller_t) :: dt_controller
    integer :: idx

    ! ------------------------------------------------------------------------ !
    ! Computation of the maximal normed difference.
    !
    ! Our goal is to freeze the simulation when the normed difference between
    ! the old and new fields is below a certain tolerance.
    ! @todo: This should be refactored into a separate function.
    !
    ! Tim I'm making the following changes:
    ! - no u_b

    u => neko_field_registry%get_field("u")
    v => neko_field_registry%get_field("v")
    w => neko_field_registry%get_field("w")
    p => neko_field_registry%get_field("p")

    if (this%have_scalar) then
       s => neko_field_registry%get_field("s")
    else
       s => null()
    end if


    if (.not. this%converged) then

       ! Compute the difference between the old and new fields
       call field_sub2(this%u_old, u)
       call field_sub2(this%v_old, v)
       call field_sub2(this%w_old, w)
       call field_sub2(this%p_old, p)
       if (this%have_scalar) then
          call field_sub2(this%s_old, s)
       end if

       ! Compute the normed difference
       normed_diff(1) = field_energy_norm(this%u_old, this%v_old, this%w_old, &
        this%case%fluid%c_Xh)
       ! divide by timestep
       normed_diff(1) = normed_diff(1)/this%case%dt
       print*, 'NORM', normed_diff(1)

       ! If the normed difference is below the tolerance, we consider the
       ! simulation to have converged. Otherwise, we copy the new fields to the
       ! old fields and continue the simulation.
       if (maxval(normed_diff) .gt. this%tol) then
          call field_copy(this%u_old, u)
          call field_copy(this%v_old, v)
          call field_copy(this%w_old, w)
          call field_copy(this%p_old, p)
          if (this%have_scalar) then
             call field_copy(this%s_old, s)
          end if

       else
          this%case%fluid%freeze = .true.
          this%converged = .true.
       end if
    end if

    ! Return if the simulation has not converged
    if (.not. this%converged .or. this%computed) then
       return
    end if

    ! ------------------------------------------------------------------------ !
    ! Computation of the adjoint field.
    !
    ! This is where the actual computation of the adjoint field should be
    ! implemented. This will so far just be a steady state field based on the
    ! current state of the fluid fields.

    ! ! Hack: Lets modify settings of the case.
    ! this%case%f_out%file_%file_type%fname = "adjoint.fld"
    ! call this%case%f_out%set_counter(0)
    ! do idx = 1, size(this%case%f_out%fluid%items)
    !    if (this%case%f_out%fluid%items(idx)%ptr%name == "u") then
    !       this%case%f_out%fluid%items(idx)%ptr => this%u_adj
    !    end if
    !    if (this%case%f_out%fluid%items(idx)%ptr%name == "v") then
    !       this%case%f_out%fluid%items(idx)%ptr => this%v_adj
    !    end if
    !    if (this%case%f_out%fluid%items(idx)%ptr%name == "w") then
    !       this%case%f_out%fluid%items(idx)%ptr => this%w_adj
    !    end if
    !    if (this%case%f_out%fluid%items(idx)%ptr%name == "p") then
    !       this%case%f_out%fluid%items(idx)%ptr => this%p_adj
    !    end if
    ! end do


    ! ------------------------------------------------------------------------ !
    ! Full copy of the `simulation.f90` file from the Neko source code.
    t_adj = 0d0
    tstep_adj = 0
    call neko_log%section('Starting adjoint')
    write(log_buf,'(A, E15.7,A,E15.7,A)') 'T  : [', 0d0,',',this%case%end_time,')'
    call neko_log%message(log_buf)
    call dt_controller%init(this%case%params)
    if (.not. dt_controller%if_variable_dt) then
       write(log_buf,'(A, E15.7)') 'dt :  ', this%case%dt
       call neko_log%message(log_buf)
    else
       write(log_buf,'(A, E15.7)') 'CFL :  ', dt_controller%set_cfl
       call neko_log%message(log_buf)
    end if

    !> Call stats, samplers and user-init before time loop
    call neko_log%section('Postprocessing')
    call this%case%q%eval(t_adj, this%case%dt, tstep_adj)
    call this%s%sample(t_adj, tstep_adj)

    call this%case%usr%user_init_modules(t_adj, this%scheme%u_adj, this%scheme%v_adj, this%scheme%w_adj,&
         this%scheme%p_adj, this%scheme%c_Xh, this%case%params)
    call neko_log%end_section()
    call neko_log%newline()

    call profiler_start
    cfl = this%scheme%compute_cfl(this%case%dt)
    start_time_org = MPI_WTIME()

    do while (t_adj .lt. this%case%end_time .and. (.not. jobctrl_time_limit()))
       call profiler_start_region('Time-Step')
       tstep_adj = tstep_adj + 1
       start_time = MPI_WTIME()
       if (dt_controller%dt_last_change .eq. 0) then
          cfl_avrg = cfl
       end if
       call dt_controller%set_dt(this%case%dt, cfl, cfl_avrg, tstep_adj)
       !calculate the cfl after the possibly varied dt
       cfl = this%scheme%compute_cfl(this%case%dt)

       call neko_log%status(t_adj, this%case%end_time)
       write(log_buf, '(A,I6)') 'Time-step: ', tstep_adj
       call neko_log%message(log_buf)
       call neko_log%begin()

       write(log_buf, '(A,E15.7,1x,A,E15.7)') 'CFL:', cfl, 'dt:', this%case%dt
       call neko_log%message(log_buf)
       call simulation_settime(t_adj, this%case%dt, this%case%ext_bdf, this%case%tlag, this%case%dtlag, tstep_adj)

       call neko_log%section('Fluid')
       call this%scheme%step(t_adj, tstep_adj, this%case%dt, this%case%ext_bdf, dt_controller)
       end_time = MPI_WTIME()
       write(log_buf, '(A,E15.7,A,E15.7)') &
            'Elapsed time (s):', end_time-start_time_org, ' Step time:', &
            end_time-start_time
       call neko_log%end_section(log_buf)

       ! Scalar step
       if (allocated(this%case%scalar)) then
          start_time = MPI_WTIME()
          call neko_log%section('Scalar')
          call this%case%scalar%step(t_adj, tstep_adj, this%case%dt, this%case%ext_bdf, dt_controller)
          end_time = MPI_WTIME()
          write(log_buf, '(A,E15.7,A,E15.7)') &
               'Elapsed time (s):', end_time-start_time_org, ' Step time:', &
               end_time-start_time
          call neko_log%end_section(log_buf)
       end if

       call neko_log%section('Postprocessing')

       call this%case%q%eval(t_adj, this%case%dt, tstep_adj)
       call this%s%sample(t_adj, tstep_adj)

       ! Update material properties
       call this%case%usr%material_properties(t_adj, tstep_adj, this%case%material_properties%rho,&
            this%case%material_properties%mu, &
            this%case%material_properties%cp, &
            this%case%material_properties%lambda, &
            this%case%params)

       call neko_log%end_section()

       call neko_log%end()
       call profiler_end_region
    end do
    call profiler_stop

    call json_get_or_default(this%case%params, 'case.output_at_end',&
         output_at_end, .true.)
    call this%s%sample(t_adj, tstep_adj, output_at_end)

    if (.not. (output_at_end) .and. t_adj .lt. this%case%end_time) then
       call simulation_joblimit_chkp(this%case, t_adj)
    end if

    call this%case%usr%user_finalize_modules(t_adj, this%case%params)

    call neko_log%end_section('Normal end.')
    this%computed = .true.

  end subroutine simcomp_test_compute

  subroutine simulation_settime(t, dt, ext_bdf, tlag, dtlag, step)
    real(kind=rp), intent(inout) :: t
    real(kind=rp), intent(in) :: dt
    type(time_scheme_controller_t), intent(inout) :: ext_bdf
    real(kind=rp), dimension(10) :: tlag
    real(kind=rp), dimension(10) :: dtlag
    integer, intent(in) :: step
    integer :: i


    do i = 10, 2, -1
       tlag(i) = tlag(i-1)
       dtlag(i) = dtlag(i-1)
    end do

    dtlag(1) = dt
    tlag(1) = t
    if (ext_bdf%ndiff .eq. 0) then
       dtlag(2) = dt
       tlag(2) = t
    end if

    t = t + dt

    call ext_bdf%set_coeffs(dtlag)

  end subroutine simulation_settime

  !> Restart a case @a C from a given checkpoint
  subroutine simulation_restart(C, t)
    implicit none
    type(case_t), intent(inout) :: C
    real(kind=rp), intent(inout) :: t
    integer :: i
    type(file_t) :: chkpf, previous_meshf
    character(len=LOG_SIZE) :: log_buf
    character(len=:), allocatable :: restart_file
    character(len=:), allocatable :: restart_mesh_file
    real(kind=rp) :: tol
    logical :: found

    call C%params%get('case.restart_file', restart_file, found)
    call C%params%get('case.restart_mesh_file', restart_mesh_file,&
         found)

    if (found) then
       previous_meshf = file_t(trim(restart_mesh_file))
       call previous_meshf%read(C%fluid%chkp%previous_mesh)
    end if

    call C%params%get('case.mesh2mesh_tolerance', tol,&
         found)

    if (found) C%fluid%chkp%mesh2mesh_tol = tol

    chkpf = file_t(trim(restart_file))
    call chkpf%read(C%fluid%chkp)
    C%dtlag = C%fluid%chkp%dtlag
    C%tlag = C%fluid%chkp%tlag

    !Free the previous mesh, dont need it anymore
    call C%fluid%chkp%previous_mesh%free()
    do i = 1, size(C%dtlag)
       call C%ext_bdf%set_coeffs(C%dtlag)
    end do

    call C%fluid%restart(C%dtlag, C%tlag)
    if (allocated(C%scalar)) call C%scalar%restart( C%dtlag, C%tlag)

    t = C%fluid%chkp%restart_time()
    call neko_log%section('Restarting from checkpoint')
    write(log_buf,'(A,A)') 'File :   ', &
         trim(restart_file)
    call neko_log%message(log_buf)
    write(log_buf,'(A,E15.7)') 'Time : ', t
    call neko_log%message(log_buf)
    call neko_log%end_section()

    call C%s%set_counter(t)
  end subroutine simulation_restart

  !> Write a checkpoint at joblimit
  subroutine simulation_joblimit_chkp(C, t)
    type(case_t), intent(inout) :: C
    real(kind=rp), intent(inout) :: t
    type(file_t) :: chkpf
    character(len=:), allocatable :: chkp_format
    character(len=LOG_SIZE) :: log_buf
    character(len=10) :: format_str
    logical :: found

    call C%params%get('case.checkpoint_format', chkp_format, found)
    call C%fluid%chkp%sync_host()
    format_str = '.chkp'
    if (found) then
       if (chkp_format .eq. 'hdf5') then
          format_str = '.h5'
       end if
    end if
    chkpf = file_t('joblimit'//trim(format_str))
    call chkpf%write(C%fluid%chkp, t)
    write(log_buf, '(A)') '! saving checkpoint >>>'
    call neko_log%message(log_buf)

  end subroutine simulation_joblimit_chkp

    function field_energy_norm(u, v, w, coef, n) result(norm)
    ! I'm new to field maths... but the problem is u,v,w are fields
    ! and coef%B is an array.
    ! But I would assume the norm we're interested in is some measure of 
    ! |du/dt|
    ! which is approximated by
    ! |(u_{n} - u_{n-1})/(t_{n} - t_{n_1})|
    !
    ! So I would imagine we would want to take:
    !
    ! sqrt(int (du**2 + dv**2 + dw**2))
    ! 
    ! and then divide by the timestep outside
    integer, intent(in), optional :: n
    type(coef_t), intent(in) :: coef
    type(field_t), intent(in) :: u,v,w
    real(kind=rp) :: norm
    integer :: size

    if (present(n)) then
       size = n
    else
       size = u%size()
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       ! fill this in later
    else
       norm = energy_norm(u%x, v%x, w%x, coef%B, size)
    end if

  	 end function field_energy_norm

  function energy_norm(a, b, c, d, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), dimension(n), intent(in) :: c
    real(kind=rp), dimension(n), intent(in) :: d
    real(kind=rp) :: energy_norm, tmp
    integer :: i, ierr

    tmp = 0.0_rp
    do i = 1, n
       tmp = tmp + (a(i)**2 +  b(i)**2 + c(i)**2) * d(i)
    end do

    call MPI_Allreduce(tmp, energy_norm, 1, &
                       MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    energy_norm = sqrt(energy_norm)
    

  end function energy_norm



end module simcomp_example

