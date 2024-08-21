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

 ! Implements the `mma_comp_t` type.
module mma_simcomp
  use num_types, only: rp
  use case, only: case_t
  use json_module, only: json_file
  use json_utils, only: json_get_or_default
  use simulation_component, only: simulation_component_t
  use field, only: field_t
  use logger, only: neko_log
  use mma, only: mma_t
  implicit none
  private

  ! An empty user defined simulation component.
  ! This is a simple example of a user-defined simulation component.
  type, public, extends(simulation_component_t) :: mma_comp_t

     real(kind=rp) :: tol !< Just some dummy variable to show it working.
     type(field_t) :: tmp !< Just some dummy field to show it working.


     type(field_t) :: designx !< Just some dummy field to show it working.
     type(field_t) :: xmax !< Just some dummy field to show it working.
     type(field_t) :: xmin !< Just some dummy field to show it working.

     integer :: m !< Just some dummy variable to show it working.
     real(kind=rp) :: a0_const !< Just some dummy variable to show it working.
     real(kind=rp) :: a_const !< Just some dummy variable to show it working.
     real(kind=rp) :: c_const !< Just some dummy variable to show it working.
     real(kind=rp) :: d_const !< Just some dummy variable to show it working.

     type(mma_t) :: mma !< The actual MMA simulation component.
    
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
  end type mma_comp_t

contains

  ! Constructor from json.
  subroutine simcomp_test_init_from_json(this, json, case)
    class(mma_comp_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case

    call this%tmp%init(case%msh, case%fluid%Xh, "tmp")
    call this%designx%init(case%msh, case%fluid%Xh, "designx")
    call this%xmax%init(case%msh, case%fluid%Xh, "xmax")
    call this%xmin%init(case%msh, case%fluid%Xh, "xmin")

    ! Read the tolerance
    call json_get_or_default(json, "tol", this%tol, 1.0e-6_rp)
    call json_get_or_default(json, "m", this%m, 1)
    call json_get_or_default(json, "a0_const", this%a0_const, 1.0_rp)
    call json_get_or_default(json, "a_const", this%a_const, 0.0_rp)
    call json_get_or_default(json, "c_const", this%c_const, 1000.0_rp)
    call json_get_or_default(json, "d_const", this%d_const, 1.0_rp)

    call this%init_from_attributes()
    call this%init_base(json, case)
  end subroutine simcomp_test_init_from_json

  ! Actual constructor.
  subroutine simcomp_test_init_from_attributes(this)
    class(mma_comp_t), intent(inout) :: this

    ! integer :: n
    ! real(kind=rp), dimension(n) ::  xmax, xmin
    real(kind=rp), allocatable ::a(:), c(:), d(:)
    real(kind=rp) :: a0

    allocate(a(this%m))
    allocate(c(this%m))
    allocate(d(this%m))
    a0= this%a0_const
    a= this%a_const
    c= this%c_const
    d= this%d_const
    ! a0= 1.0
    ! a= 0
    ! c= 1000.0
    ! d= 1.0

    ! n=this%designx%dof%size()


    call this%mma%init(this%designx%x, this%designx%dof%size(), this%m, a0, a, c, d, this%xmin%x, this%xmax%x)
    print *, "yooooooooooooo"
    print *, this%designx%dof%size()
    ! print *, size(this%designx%x)
    ! print *, size(this%designx%x,1)
    ! print *, size(this%designx%x,2)
    ! print *, size(this%designx%x,3)
    ! print *, size(this%designx%x,4)
    ! print *, size(this%designx%msh%elements)
    ! print *, this%designx%msh%elements(160)%e%pts(1)%p%x

    ! print *, this%designx%x(10)

  end subroutine simcomp_test_init_from_attributes

  ! Destructor.
  subroutine simcomp_test_free(this)
    class(mma_comp_t), intent(inout) :: this

    call this%tmp%free()
    call this%mma%free()

    call this%free_base()
  end subroutine simcomp_test_free

  ! Computations.
  subroutine simcomp_test_compute(this, t, tstep)
    class(mma_comp_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    !     ! Declarations
    ! integer :: n, nrhs, lda, ldb, info
    ! integer, allocatable :: ipiv(:)
    ! double precision, allocatable :: A(:,:), B(:,:)

    ! ! Matrix and vector dimensions
    ! n = 3     ! Size of the matrix A (n x n)
    ! nrhs = 1  ! Number of right-hand sides, i.e., number of columns of B
    ! lda = n   ! Leading dimension of A
    ! ldb = n   ! Leading dimension of B

    ! ! Allocate arrays
    ! allocate(A(n, n), B(n, nrhs), ipiv(n))

    ! ! Define the matrix A (example: A = [3, 1, 2; 1, 4, 0; 2, 0, 5])
    ! A = reshape([3.0d0, 1.0d0, 2.0d0, &
    !              1.0d0, 4.0d0, 0.0d0, &
    !              2.0d0, 0.0d0, 5.0d0], &
    !              shape(A))

    ! ! Define the right-hand side vector B (example: B = [10; 10; 10])
    ! B = reshape([10.0d0, 10.0d0, 10.0d0], shape(B))

    ! ! Call DGESV to solve the system A * X = B
    ! call dgesv(n, nrhs, A, lda, ipiv, B, ldb, info)

    ! ! Check if the solution was successful
    ! if (info == 0) then
    !     print *, "Solution found:"
    !     print *, "X = ", B
    ! else
    !     print *, "DGESV failed with error code ", info
    ! end if

    real(kind=rp), allocatable ::df0dx(:), fval(:), dfdx(:,:)
   
    allocate(df0dx(this%designx%dof%size()))
    allocate(fval(this%m))
    allocate(dfdx(this%m,this%designx%dof%size()))

    df0dx=0
    fval=0
    dfdx=0
    ! update(iter, x, df0dx, fval, dfdx)
    call this%mma%update(this%m, this%designx%x, df0dx, fval, dfdx)

    print *, "The stuff is computed and this%designx%x is updated"

  end subroutine simcomp_test_compute

end module mma_simcomp

