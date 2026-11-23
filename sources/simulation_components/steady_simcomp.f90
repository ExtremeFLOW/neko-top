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
  use time_based_controller, only: time_based_controller_t
  use utils, only: neko_error
  use logger, only: neko_log
  implicit none
  private
  public :: steady_simcomp_allocate

  !> The `steady_simcomp_t` type is a simulation component that terminates a
  !! simulation when the normed difference between the old and new fields is
  !! below a certain tolerance. This allow us to compute steady-state solutions
  !! without modifying the simulation loop.
  !!
  !! Only the fluid is frozen. The scalar keeps evolving, so its residual is
  !! still followed, and a warning is issued when the fluid is frozen while
  !! the scalar has not converged, as well as at the end of such a run if it
  !! still has not.
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

     !> Whether this component froze the fluid. Only read while the fluid is
     !! frozen, so a reset that unfreezes the fluid needs no reset here.
     logical :: froze_fluid = .false.
     !> Whether the scalar was not converged when the fluid was frozen, so
     !! that it is to be checked again at the end of the run.
     logical :: check_scalar_at_end = .false.

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

    ! A frozen field is not interesting to compute differences for. The
    ! scalar is not frozen with the fluid, though, so if it was this component
    ! that froze the fluid, the scalar is still followed.
    if (this%case%fluid%freeze) then
       if (this%froze_fluid .and. associated(this%s)) then
          call steady_simcomp_follow_scalar(this, time)
       end if
       return
    end if

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
          call steady_simcomp_write_row(this, tstep, t, normed_diff_fluid, &
               normed_diff_scalar(1))
       end if
    end if

    if (converged_fluid) then
       this%case%fluid%freeze = .true.
       this%froze_fluid = .true.
       this%check_scalar_at_end = .false.
       ! @todo
       ! we should implement a freeze for the scalar too.

       ! The scalar freeze should only happen after the fluid has converged.
       ! if (converged_scalar) this%case%scalar%freeze = .true.

       if (associated(this%s)) then
          ! Nothing was copied above if everything converged, which leaves the
          ! difference in the old scalar. Following the scalar from here on
          ! needs the scalar itself.
          if (converged) call field_copy(this%s_old, this%s)

          ! Without coupling the scalar took no part in the decision.
          if (.not. this%scalar_coupled .and. .not. converged_scalar) then
             call neko_log%warning('steady: fluid frozen at t = ' // &
                  real_str(time%t) // ' with the scalar not converged ' // &
                  '(residual ' // real_str(real(normed_diff_scalar(1), dp)) &
                  // ' > tol ' // real_str(this%tol) // '). The scalar ' // &
                  'is not frozen and keeps evolving, so any scalar ' // &
                  'objective is evaluated on a scalar that is not steady. ' &
                  // 'Set "scalar_coupled": true to freeze the fluid only ' &
                  // 'once the scalar has converged too.')
             this%check_scalar_at_end = .true.
          end if
       end if
    end if
  end subroutine steady_simcomp_compute

  !> Follow the scalar once this component has frozen the fluid.
  !! @param time The current time state.
  !! @details The scalar is not frozen with the fluid and keeps evolving, so
  !! its residual is still computed and logged. The fluid columns of the log
  !! hold zero, the exact increment of a frozen field. If the scalar was not
  !! converged when the fluid was frozen, it is checked again on the last
  !! computation of the run.
  subroutine steady_simcomp_follow_scalar(this, time)
    class(steady_simcomp_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    real(kind=rp), dimension(4), parameter :: frozen_fluid = 0.0_rp
    real(kind=rp) :: t, dt, residual

    t = time%t
    dt = this%case%time%dt

    call field_sub2(this%s_old, this%s)
    residual = energy_norm(this%s_old, this%case%fluid%C_Xh, dt)
    call field_copy(this%s_old, this%s)

    if (mod(time%tstep, this%log_frequency) .eq. 0) then
       call steady_simcomp_write_row(this, time%tstep, t, frozen_fluid, &
            residual)
    end if

    if (.not. this%check_scalar_at_end) return
    if (.not. steady_simcomp_is_last_compute(this, time)) return
    this%check_scalar_at_end = .false.

    ! The last check need not be on the final step, so the time given is
    ! that of the check.
    if (residual .le. this%tol) then
       call neko_log%message('steady: run ending, scalar converged since ' &
            // 'the fluid was frozen (t = ' // real_str(time%t) // &
            ', residual ' // real_str(real(residual, dp)) // ' <= tol ' // &
            real_str(this%tol) // ').')
    else
       call neko_log%warning('steady: run ending, scalar still not ' // &
            'converged at the last check (t = ' // real_str(time%t) // &
            ', residual ' // real_str(real(residual, dp)) // ' > tol ' // &
            real_str(this%tol) // '). Any scalar objective is evaluated ' // &
            'on a scalar that is not steady. Set "scalar_coupled": true, ' // &
            'with an end_time long enough for the scalar to converge.')
    end if
  end subroutine steady_simcomp_follow_scalar

  !> Whether the run ends before this component is computed again.
  !! @param time The current time state.
  !! @return True if this is the last computation of the run.
  !! @details The component's own compute controller decides, so that the
  !! answer follows its schedule whatever the control mode. It is asked on
  !! copies of itself and of the time state, which leaves the originals
  !! untouched: the copied controller registers the current computation, as
  !! the real one does after it, and the copied time is advanced as the time
  !! loop advances it, up to the step at which the run is done.
  function steady_simcomp_is_last_compute(this, time) result(is_last)
    class(steady_simcomp_t), intent(in) :: this
    type(time_state_t), intent(in) :: time
    logical :: is_last
    type(time_based_controller_t) :: controller
    type(time_state_t) :: ahead

    controller = this%compute_controller
    call controller%register_execution(time)

    ahead = time
    is_last = .true.
    do while (is_last .and. .not. ahead%is_done())
       ahead%tstep = ahead%tstep + 1
       ahead%t = ahead%t + ahead%dt
       is_last = .not. controller%check(ahead)
    end do

    call controller%free()

  end function steady_simcomp_is_last_compute

  !> Write one row of `steady_state_data.csv`.
  !! @param tstep The time step.
  !! @param t The time.
  !! @param fluid The normed increments of `u`, `v`, `w` and `p`.
  !! @param scalar The normed increment of the scalar.
  subroutine steady_simcomp_write_row(this, tstep, t, fluid, scalar)
    class(steady_simcomp_t), intent(inout) :: this
    integer, intent(in) :: tstep
    real(kind=rp), intent(in) :: t
    real(kind=rp), dimension(4), intent(in) :: fluid
    real(kind=rp), intent(in) :: scalar

    this%log_data%x(1) = real(tstep, kind=rp)
    this%log_data%x(2) = t
    this%log_data%x(3:6) = fluid
    this%log_data%x(7) = scalar
    call this%logger%write(this%log_data)

  end subroutine steady_simcomp_write_row

  !> Format a real number for a log message.
  !! @param x The number.
  !! @return The number in scientific notation, without padding.
  function real_str(x) result(str)
    real(kind=dp), intent(in) :: x
    character(len=:), allocatable :: str
    character(len=16) :: buf

    write(buf, '(ES12.5)') x
    str = trim(adjustl(buf))

  end function real_str

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
