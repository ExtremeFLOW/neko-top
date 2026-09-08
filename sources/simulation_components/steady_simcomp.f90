!> @file steady_simcomp.f90
!! @copyright
!! Copyright (c) 2024-2026, The Neko-TOP Authors
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

 !> Implements the `steady_simcomp_t` type.
module steady_simcomp
  use simulation_component, only: simulation_component_t
  use num_types, only: rp, dp
  use field, only: field_t
  use json_module, only: json_file
  use json_utils, only: json_get_or_default
  use case, only: case_t
  use field_math, only: field_sub2, field_copy
  use math, only: glsc3
  use device_math, only: device_glsc3
  use coefs, only: coef_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use csv_file, only : csv_file_t
  use vector, only: vector_t
  use time_state, only: time_state_t
  use utils, only: neko_error
  implicit none
  private
  public :: steady_simcomp_allocate

  !> The `steady_simcomp_t` type is a simulation component that terminates a
  !! simulation when the normed difference between the old and new fields is
  !! below a certain tolerance. This allow us to compute steady-state solutions
  !! without modifying the simulation loop.
  type, public, extends(simulation_component_t) :: steady_simcomp_t
     private

     type(field_t), pointer :: u, v, w, p, s

     ! Old fields
     type(field_t) :: u_old, v_old, w_old, p_old, s_old
     real(kind=dp) :: tol

     !> If true, the fluid and scalar fields are considered to be coupled, and
     !! the simulation will only be considered converged if both the fluid and
     !! scalar fields have converged.
     logical :: scalar_coupled = .false.

     !> A file writer to document the convergence of the steady state
     type(csv_file_t) :: logger
     !> Log every ___ iterations
     integer :: log_frequency
     ! vector to store the data being logged
     type(vector_t) :: log_data

   contains
     ! Constructor from json, wrapping the actual constructor.
     procedure, public, pass(this) :: init => steady_simcomp_init_from_json
     ! Actual constructor.
     procedure, public, pass(this) :: init_from_attributes => &
          steady_simcomp_init_from_attributes
     ! Destructor.
     procedure, public, pass(this) :: free => steady_simcomp_free
     ! Compute the steady_simcomp field.
     procedure, public, pass(this) :: compute_ => steady_simcomp_compute
  end type steady_simcomp_t

contains

  !> Allocator for the steady simulation component.
  subroutine steady_simcomp_allocate(obj)
    class(simulation_component_t), allocatable, intent(inout) :: obj
    allocate(steady_simcomp_t::obj)
  end subroutine steady_simcomp_allocate

  ! Constructor from json.
  subroutine steady_simcomp_init_from_json(this, json, case)
    class(steady_simcomp_t), intent(inout), target :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    real(kind=dp) :: tol
    integer :: log_frequency

    call this%init_base(json, case)

    ! Read the tolerance
    call json_get_or_default(json, "tol", tol, 1.0e-6_dp)
    ! Read the log frequency
    call json_get_or_default(json, "log_frequency", log_frequency, 50)
    call json_get_or_default(json, "scalar_coupled", this%scalar_coupled, .false.)

    call this%init_from_attributes(tol, log_frequency)

  end subroutine steady_simcomp_init_from_json

  ! Actual constructor.
  subroutine steady_simcomp_init_from_attributes(this, tol, log_frequency)
    class(steady_simcomp_t), intent(inout) :: this
    real(kind=dp), intent(in) :: tol
    integer, intent(in) :: log_frequency

    this%tol = tol
    this%log_frequency = log_frequency

    this%u => this%case%fluid%u
    this%v => this%case%fluid%v
    this%w => this%case%fluid%w
    this%p => this%case%fluid%p
    if (allocated(this%case%scalars)) then
       if (size(this%case%scalars%scalar_fields) .gt. 1) then
          call neko_error('steady simcomp only works for a single scalar')
       end if
       this%s => this%case%scalars%scalar_fields(1)%scalar%s
    else
       nullify(this%s)
    end if

    ! initialize the logger
    call this%logger%init('steady_state_data.csv')
    call this%logger%set_header('iter,time,u,v,w,p,scalar')
    call this%log_data%init(7)

    ! Point local fields to the scratch fields
    call this%u_old%init(this%u%dof)
    call this%v_old%init(this%v%dof)
    call this%w_old%init(this%w%dof)
    call this%p_old%init(this%p%dof)

    ! Check if the scalar field is allocated
    if (associated(this%s)) then
       call this%s_old%init(this%s%dof)
    end if

  end subroutine steady_simcomp_init_from_attributes

  ! Destructor.
  subroutine steady_simcomp_free(this)
    class(steady_simcomp_t), intent(inout) :: this

    if (associated(this%u)) nullify(this%u)
    if (associated(this%v)) nullify(this%v)
    if (associated(this%w)) nullify(this%w)
    if (associated(this%p)) nullify(this%p)
    if (associated(this%s)) nullify(this%s)

    call this%u_old%free()
    call this%v_old%free()
    call this%w_old%free()
    call this%p_old%free()
    call this%s_old%free()
    call this%log_data%free()

    call this%free_base()

  end subroutine steady_simcomp_free

  ! Compute the steady_simcomp field.
  subroutine steady_simcomp_compute(this, time)
    class(steady_simcomp_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    integer :: tstep
    real(kind=rp) :: t
    real(kind=rp) :: dt

    real(kind=rp), dimension(4) :: normed_diff_fluid
    real(kind=rp), dimension(1) :: normed_diff_scalar
    logical :: converged_fluid, converged_scalar, converged

    ! A frozen field is not interesting to compute differences for.
    if (this%case%fluid%freeze) return

    ! Get the current time and time step
    t = time%t
    tstep = time%tstep
    dt = this%case%time%dt

    ! ------------------------------------------------------------------------ !
    ! Computation of the maximal normed difference.
    !
    ! Our goal is to freeze the simulation when the normed difference between
    ! the old and new fields is below a certain tolerance.

    ! Compute the difference between the old and new fields
    call field_sub2(this%u_old, this%u)
    call field_sub2(this%v_old, this%v)
    call field_sub2(this%w_old, this%w)
    call field_sub2(this%p_old, this%p)
    if (associated(this%s)) then
       call field_sub2(this%s_old, this%s)
    end if

    ! Here we compute the squared difference between the old and new fields
    ! and store the result in the `normed_diff` array.
    normed_diff_fluid(1) = energy_norm(this%u_old, this%case%fluid%C_Xh, dt)
    normed_diff_fluid(2) = energy_norm(this%v_old, this%case%fluid%C_Xh, dt)
    normed_diff_fluid(3) = energy_norm(this%w_old, this%case%fluid%C_Xh, dt)
    normed_diff_fluid(4) = energy_norm(this%p_old, this%case%fluid%C_Xh, dt)
    if (associated(this%s)) then
       normed_diff_scalar(1) = energy_norm(this%s_old, this%case%fluid%C_Xh, dt)
    else
       normed_diff_scalar(1) = 0.0_rp
    end if

    converged_fluid = maxval(normed_diff_fluid) .le. this%tol
    if (associated(this%s)) then
       converged_scalar = maxval(normed_diff_scalar) .le. this%tol
    else
       converged_scalar = .true.
    end if
    converged = converged_fluid .and. converged_scalar

    if (this%scalar_coupled) then
       converged_fluid = converged
       converged_scalar = converged
    end if

    ! If the normed difference is below the tolerance, we consider the
    ! simulation to have converged. Otherwise, we copy the new fields to the
    ! old fields and continue the simulation.
    if (.not. converged) then
       call field_copy(this%u_old, this%u)
       call field_copy(this%v_old, this%v)
       call field_copy(this%w_old, this%w)
       call field_copy(this%p_old, this%p)
       if (associated(this%s)) then
          call field_copy(this%s_old, this%s)
       end if

       ! this could be a csv, but I think it's nicer in the logfile, it can
       ! always be greped out later if you want to plot.
       if (mod(tstep, this%log_frequency) .eq. 0) then
          this%log_data%x(1) = real(tstep, kind=rp)
          this%log_data%x(2) = t
          this%log_data%x(3) = normed_diff_fluid(1)
          this%log_data%x(4) = normed_diff_fluid(2)
          this%log_data%x(5) = normed_diff_fluid(3)
          this%log_data%x(6) = normed_diff_fluid(4)
          this%log_data%x(7) = normed_diff_scalar(1)
          call this%logger%write(this%log_data)
       end if
    end if

    if (converged_fluid) then
       this%case%fluid%freeze = .true.
       ! @todo
       ! we should implement a freeze for the scalar too.

       ! The scalar freeze should only happen after the fluid has converged.
       ! if (converged_scalar) this%case%scalar%freeze = .true.
    end if
  end subroutine steady_simcomp_compute

  !> A norm for a scalar field based on the energy.
  !! Assuming the time derivative $\frac{\partial u}{\partial t}
  !! \approx \frac{\Delta u}{\Delta t}$, this norm computes the spatial
  !! integral normalized by the domain volume.
  function energy_norm(delta_fld, coef, dt)
    type(field_t), intent(in) :: delta_fld
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(in) :: dt
    real(kind=rp) :: energy_norm, tmp
    integer :: n

    n = delta_fld%size()
    if (NEKO_BCKND_DEVICE .eq. 1) then
       tmp = device_glsc3(delta_fld%x_d, delta_fld%x_d, coef%B_d, n)
    else
       tmp = glsc3(delta_fld%x, delta_fld%x, coef%B, n)
    end if
    energy_norm = sqrt(tmp) / dt / coef%volume

  end function energy_norm


end module steady_simcomp
