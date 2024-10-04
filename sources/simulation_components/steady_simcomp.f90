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

 !> Implements the `steady_simcomp_t` type.
module steady_simcomp
  use simulation_component, only: simulation_component_t
  use num_types, only: rp, dp
  use field, only: field_t
  use json_module, only: json_file
  use json_utils, only: json_get_or_default
  use case, only: case_t
  use coefs, only: coef_t
  use field_math, only: field_sub2, field_glsc2, field_copy
  implicit none
  private

  !> The `steady_simcomp_t` type is a simulation component that terminates a
  !! simulation when the normed difference between the old and new fields is
  !! below a certain tolerance. This allow us to compute steady-state solutions
  !! without modifying the simulation loop.
  type, public, extends(simulation_component_t) :: steady_simcomp_t
     private

     ! Old fields
     type(field_t) :: u_old, v_old, w_old, p_old, s_old
     real(kind=dp) :: tol

     logical :: have_scalar = .false.

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

  ! Constructor from json.
  subroutine steady_simcomp_init_from_json(this, json, case)
    class(steady_simcomp_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    real(kind=dp) :: tol

    call this%init_base(json, case)

    ! Read the tolerance
    call json_get_or_default(json, "tol", tol, 1.0e-6_dp)

    call this%init_from_attributes(tol)

  end subroutine steady_simcomp_init_from_json

  ! Actual constructor.
  subroutine steady_simcomp_init_from_attributes(this, tol)
    class(steady_simcomp_t), intent(inout) :: this
    real(kind=dp), intent(in) :: tol

    this%tol = tol

    ! Point local fields to the scratch fields
    call this%u_old%init(this%case%fluid%u%dof)
    call this%v_old%init(this%case%fluid%v%dof)
    call this%w_old%init(this%case%fluid%w%dof)
    call this%p_old%init(this%case%fluid%p%dof)

    ! Check if the scalar field is allocated
    if (allocated(this%case%scalar)) then
       this%have_scalar = .true.
       call this%s_old%init(this%case%scalar%s%dof)
    end if

  end subroutine steady_simcomp_init_from_attributes

  ! Destructor.
  subroutine steady_simcomp_free(this)
    class(steady_simcomp_t), intent(inout) :: this

    call this%u_old%free()
    call this%v_old%free()
    call this%w_old%free()
    call this%p_old%free()
    if (this%have_scalar) then
       call this%s_old%free()
    end if

    call this%free_base()

  end subroutine steady_simcomp_free

  ! Compute the steady_simcomp field.
  subroutine steady_simcomp_compute(this, t, tstep)
    class(steady_simcomp_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    real(kind=rp), dimension(5) :: normed_diff
    type(field_t), pointer :: u, v, w, p, s

    ! A frozen field is not interesting to compute differences for.
    if (this%case%fluid%freeze) return

    ! ------------------------------------------------------------------------ !
    ! Computation of the maximal normed difference.
    !
    ! Our goal is to freeze the simulation when the normed difference between
    ! the old and new fields is below a certain tolerance.

    u => this%case%fluid%u
    v => this%case%fluid%v
    w => this%case%fluid%w
    p => this%case%fluid%p

    if (this%have_scalar) then
       s => this%case%scalar%s
    else
       s => null()
    end if

    ! Compute the difference between the old and new fields
    call field_sub2(this%u_old, u)
    call field_sub2(this%v_old, v)
    call field_sub2(this%w_old, w)
    call field_sub2(this%p_old, p)
    if (this%have_scalar) then
       call field_sub2(this%s_old, s)
    end if

    ! Here we compute the squared difference between the old and new fields
    ! and store the result in the `normed_diff` array.
    normed_diff(1) = field_glsc2(this%u_old, this%u_old)
    normed_diff(2) = field_glsc2(this%v_old, this%v_old)
    normed_diff(3) = field_glsc2(this%w_old, this%w_old)
    normed_diff(4) = field_glsc2(this%p_old, this%p_old)
    if (this%have_scalar) then
       normed_diff(5) = field_glsc2(this%s_old, this%s_old)
    else
       normed_diff(5) = 0.0_rp
    end if

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
    end if


  end subroutine steady_simcomp_compute

end module steady_simcomp

