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
  use num_types, only: rp
  use json_module, only: json_file
  use json_utils, only: json_get_or_default
  use simulation_component, only: simulation_component_t
  use case, only: case_t
  use field, only: field_t
  use field_registry, only: neko_field_registry
  use scratch_registry, only: neko_scratch_registry
  use neko_config, only: NEKO_BCKND_DEVICE
  use field_math, only: field_cfill, field_sub2, field_copy, field_glsc2
  use field_math, only: field_add2
  use math, only: glsc2
  use device_math, only: device_glsc2
  use adv_lin_no_dealias, only: adv_lin_no_dealias_t
  use logger, only: neko_log
  use comm
  implicit none
  private

  ! An empty user defined simulation component.
  ! This is a simple example of a user-defined simulation component.
  type, public, extends(simulation_component_t) :: adjoint_t

     type(field_t) :: u_old, v_old, w_old, p_old, s_old
     real(kind=rp) :: tol

     logical :: have_scalar = .false.

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

  end subroutine simcomp_test_init_from_json

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
    call this%free_base()

  end subroutine simcomp_test_free

  ! Compute the simcomp_test field.
  subroutine simcomp_test_compute(this, t, tstep)
    class(adjoint_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    real(kind=rp), dimension(5) :: normed_diff
    type(field_t), pointer :: u, v, w, p, s
    type(field_t), pointer :: u_base, v_base, w_base, p_base, s_base
    character(len=256) :: msg

    logical :: converged = .false.

    ! ------------------------------------------------------------------------ !
    ! Computation of the maximal normed difference.
    !
    ! Our goal is to freeze the simulation when the normed difference between
    ! the old and new fields is below a certain tolerance.
    ! @todo: This should be refactored into a separate function.

    if (.not. this%case%fluid%freeze) then
       u => neko_field_registry%get_field("u")
       v => neko_field_registry%get_field("v")
       w => neko_field_registry%get_field("w")
       p => neko_field_registry%get_field("p")
       if (this%have_scalar) then
          s => neko_field_registry%get_field("s")
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

       ! Compute the normed difference
       normed_diff(1) = field_glsc2(this%u_old, this%u_old)
       normed_diff(2) = field_glsc2(this%v_old, this%v_old)
       normed_diff(3) = field_glsc2(this%w_old, this%w_old)
       normed_diff(4) = field_glsc2(this%p_old, this%p_old)
       if (this%have_scalar) then
          normed_diff(5) = field_glsc2(this%s_old, this%s_old)
       end if

       write(msg, '(A,ES8.2)' ) "normed_diff = ", maxval(normed_diff)
       call neko_log%message(msg)

       if (maxval(normed_diff) .gt. this%tol) then
          ! Copy the new fields to the old fields
          call field_copy(this%u_old, u)
          call field_copy(this%v_old, v)
          call field_copy(this%w_old, w)
          call field_copy(this%p_old, p)
          if (this%have_scalar) then
             call field_copy(this%s_old, s)
          end if

          return
          !  else if (.not. converged) then
          !     u_base => neko_field_registry%get_field("u_b")
          !     v_base => neko_field_registry%get_field("v_b")
          !     w_base => neko_field_registry%get_field("w_b")
          !     p_base => neko_field_registry%get_field("p_b")
          !     if (this%have_scalar) then
          !        s_base => neko_field_registry%get_field("s_b")
          !     else
          !        s_base => null()
          !     end if

          !     call field_add2(u_base, this%u_old)
          !     call field_add2(v_base, this%v_old)
          !     call field_add2(w_base, this%w_old)
          !     call field_add2(p_base, this%p_old)
          !     if (this%have_scalar) then
          !        call field_add2(s_base, this%s_old)
          !     end if

          !     this%case%fluid%toggle_adjoint = .true.
          !     converged = .true.

          !  else
          !     this%case%fluid%freeze = .true.
       end if
    end if

    ! ------------------------------------------------------------------------ !
    ! Computation of the adjoint field.
    !
    ! This is where the actual computation of the adjoint field should be
    ! implemented. This will so far just be a steady state field based on the
    ! current state of the fluid fields.





  end subroutine simcomp_test_compute

end module simcomp_example

