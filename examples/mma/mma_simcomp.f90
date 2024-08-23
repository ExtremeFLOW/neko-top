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
  use neko_config
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!Functions to read/write a filed_t object from/to a vector!!!!!!!!!!!!!!
  subroutine write_field_to_vector(field, vector, size)
    class(field_t), intent(inout) :: field
    integer, intent(in) :: size
    real(kind=rp), intent(inout), dimension(size) :: vector

    integer :: i, j, k, l, idx
    idx=1
    if (neko_bcknd_device .eq. 1) then
      !call device_cfill(this%x_d, a, this%size())
      !stuf for gpus
    else
      do i = 1, field%msh%nelv
          do l = 1, field%Xh%lz
            do k = 1, field%Xh%ly
                do j = 1, field%Xh%lx
                  vector(idx)=field%x(j, k, l, i)
                  idx = idx + 1
                end do
            end do
          end do
      end do
    end if
  end subroutine write_field_to_vector

  subroutine write_vector_to_field(field, vector, size)
    class(field_t), intent(inout) :: field
    integer, intent(in) :: size
    real(kind=rp), intent(inout), dimension(size) :: vector

    integer :: i, j, k, l, idx
    idx=1
    if (neko_bcknd_device .eq. 1) then
      !call device_cfill(this%x_d, a, this%size())
      !stuf for gpus
    else
      do i = 1, field%msh%nelv
          do l = 1, field%Xh%lz
            do k = 1, field%Xh%ly
                do j = 1, field%Xh%lx
                  field%x(j, k, l, i) = vector(idx)
                  idx = idx + 1
                end do
            end do
          end do
      end do
    end if
  end subroutine write_vector_to_field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    call json_get_or_default(json, "m", this%m, 2)
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

    integer :: n
    real(kind=rp), allocatable ::a(:), c(:), d(:), x(:)
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
    this%designx%x=0
    n=this%designx%dof%size()
    allocate(x(n))
    call write_field_to_vector(this%designx, x, n)
    call this%mma%init(x, n, this%m, a0, a, c, d, this%xmin%x, this%xmax%x)

    print *, "yooooooooooooo"
    ! print *, this%designx%dof%size()

    print *, size(this%designx%x)
    print *, size(this%designx%x,1)
    print *, size(this%designx%x,2)
    print *, size(this%designx%x,3)
    print *, size(this%designx%x,4)
    ! print *, size(this%designx%msh%elements)
    ! ! print *, this%designx%msh%elements(160)%e%pts(1)%p%x

    ! ! print *, this%designx%x(10)

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
    !!!!!!! TEST_CASE_3
    real(kind=rp), parameter :: L = 4 !sqrt(1.25)+ sqrt(9.25)
    real(kind=rp), parameter :: x0 = 1.0_rp 
    real(kind=rp), parameter :: xL= 3.0_rp 
    integer :: iter, i
    real(kind=rp) :: start_time, end_time
    real(kind=rp), dimension(this%mma%n) :: x
    
    ! integer :: n, m
    ! ! n=this%designx%dof%size()
    ! n=this%mma%n
    ! m = this%mma%m
    real(kind=rp), dimension(this%mma%m) :: fval
    real(kind=rp), dimension(this%mma%m,this%mma%n) :: dfdx
    real(kind=rp) :: f0val, tol, f0valeps
    real(kind=rp), dimension(this%mma%n) :: df0dx
    character(len=50) :: filename
    ! print *, "this%designx%dof%size()=", this%designx%dof%size()
    ! print *, "size(this%mma%x)=", size(this%mma%x)
    
    this%mma%xmax = 3.0_rp
    this%mma%xmin = 0.0_rp
    ! print *, "size(this%xmin%x)=", size(this%xmin%x)
    ! print *, "size(this%mma%xmin)=", size(this%mma%xmin)
    
!     x=0.00001
    ! call cpu_time(start_time)
    call write_field_to_vector(this%designx,x,this%mma%n) 

    call func1 (this%mma%n, this%mma%m, L, x, x0, xL, f0val, df0dx, fval , dfdx)
    print *, f0val


    filename = 'Catenary.csv'
    ! Open the file in append mode to keep adding data in each iteration
    open(unit=10, file=filename, status='replace', action='write')
    ! Write the header row with column names
    write(10, '(A)', advance='no') 'Iteration,'
    do i = 0, this%mma%n+3
        if (i < this%mma%n+2) then
            write(10, '(A)', advance='no') 'x(' // trim(adjustl(to_string(i))) // '),'
        else if (i .eq. this%mma%n+2) then
            write(10, '(A)', advance='no') 'f0val,'
        else
            write(10, '(A)') 'fval'
        end if
    end do
!write the initial design
    write(10, '(I0, A)', advance='no') 0, ','
    do i = 0, this%mma%n+3
        if (i .eq. 0) then
            write(10, '(F12.8, A)', advance='no') x0, ','
        else if (i < this%mma%n+1) then
            write(10, '(F12.8, A)', advance='no') x(i), ','
        else if (i .eq. this%mma%n+1) then
            write(10, '(F12.8, A)', advance='no') xL, ','
        else if (i .eq. this%mma%n+2) then
            write(10, '(F12.8, A)', advance='no') f0val, ','
        else
            write(10, '(F12.8)') fval(1)
        end if
    end do
  ! The optimization loop
    do iter = 1, 1000
      call this%mma%update(iter, x, df0dx, fval, dfdx)
      ! this%designx%x = this%mma%x
      call func1 (this%mma%n, this%mma%m, L, x, x0, xL, f0val, df0dx, fval , dfdx)

      call this%mma%KKT(x,df0dx,fval,dfdx)
            print *, 'iter=', iter,&
            '-------,f0val= ', f0val, ',   fval= ', fval(1), &
            ',  KKTmax=', this%mma%residumax, ', KKTnorm2=', this%mma%residunorm


            ! Write the iteration number and all x values in the same row
      write(10, '(I0, A)', advance='no') iter, ','
      do i = 0, this%mma%n+3
          if (i .eq. 0) then
              write(10, '(F12.8, A)', advance='no') x0, ','
          else if (i < this%mma%n+1) then
              write(10, '(F12.8, A)', advance='no') x(i), ','
          else if (i .eq. this%mma%n+1) then
              write(10, '(F12.8, A)', advance='no') xL, ','
          else if (i .eq. this%mma%n+2) then
              write(10, '(F12.8, A)', advance='no') f0val, ','
          else
              write(10, '(F12.8)') fval(1)
          end if
      end do

      if (this%mma%residunorm .lt. 1.0e-8_rp) exit
    end do

    ! print *, "f0val=", f0val, "fval=", fval
    call cpu_time(end_time)


    print *, 'Elapsed Time: ', end_time - start_time, ' seconds'

    close(10)
    print *, 'Data written to ', trim(filename)

    call write_vector_to_field(this%designx, x, this%mma%n)
    print *, "The stuff is computed and this%designx%x is updated"

  end subroutine simcomp_test_compute

pure function to_string(i) result(str)
    integer, intent(in) :: i
    character(len=11) :: str
    write(str, '(I0)') i
end function to_string

subroutine func1 (n, m, L, x, x0, xL, f0val, df0dx, fval , dfdx)
    implicit none

    integer, intent(in) :: n, m
    real(kind=rp), intent(in) :: x0, xL, L
    real(kind=rp), dimension(n), intent(in) :: x
    real(kind=rp), intent(inout) :: f0val
    real(kind=rp), dimension(n), intent(inout) :: df0dx
    real(kind=rp), dimension(m), intent(inout) :: fval
    real(kind=rp), dimension(m,n), intent(inout) :: dfdx
    real(kind=rp) :: h, u_j, v_j, u_jm1, v_jm1
    real(kind=rp), dimension(n+2) :: xx
    integer :: j
    ! ----------------------------------------------------------- !
    !  This file calculates function values and gradients         !
    !  for "toy problem 3":                                       !
    !                                                             !
    !    minimize sum_(j=1,..,n) {(x[j]+x[j-1])/2 *               !
    !                    sqrt(1+((x[j]-x[j-1])/h))^2}             !
    !  subject to sum_(j=1,..,n) h*sqrt(1+((x[j]-x[j-1])/h))^2-L=0!
    !              -L =< x(j) =< L, for j=1,...,n                 !
    ! ----------------------------------------------------------- !
    xx(2:n+1)=x
    xx(1)= x0
    xx(n+2)= xL
    ! print *, "xx=", xx
    h=1/(1.0*n+1)
    f0val=0
    do j=2, n+2
        f0val = f0val + h*(xx(j)+xx(j-1))/2*sqrt(1+((xx(j)-xx(j-1))/h)**2)
    end do

    !!!!! differentiating the discrete f=h*u_j*v_j wrt xj
    !!!!! note that the term for xj repeats two times in the sum
    do j=1, n
        ! Compute u and v for j
        u_j = (xx(j+1) + xx(j)) / 2.0
        v_j = sqrt(1.0 + ((xx(j+1) - xx(j)) / h)**2)


        ! for first time xj showup in the equations
        df0dx(j) = h * (0.5 * v_j + &
                u_j * (xx(j+1) - xx(j)) / (h**2 * v_j))

        ! if (j .lt. n) then
        ! Compute u and v for the second time xj shows up in the sum
        u_jm1 = (xx(j+2) + xx(j+1)) / 2.0
        v_jm1 = sqrt(1.0 + ((xx(j+2) - xx(j+1)) / h)**2)
        df0dx(j) = df0dx(j) + h * (0.5 * v_jm1 - &
            u_jm1  * (xx(j+2) - xx(j+1)) / (h**2 * v_jm1))
        ! end if

    end do

    fval(1)=0
    do j=1, n+1
        fval(1) = fval(1) + h*sqrt(1+((xx(j+1)-xx(j))/h)**2)
    end do
    fval(1) = fval(1) -L

    !!!!! differentiating the discrete f=h*v_j wrt xj
    !!!!! note that the term for xj repeats two times in the sum
    do j=1, n
        ! Compute u and v for j
        v_j = sqrt(1.0 + ((xx(j+1) - xx(j)) / h)**2)

        ! for first time xj showup in the equations
        dfdx(1,j)= h *((xx(j+1) - xx(j)) / (h**2 * v_j))

        v_jm1 = sqrt(1.0 + ((xx(j+2) - xx(j+1)) / h)**2)
        dfdx(1,j)= dfdx(1,j) - h * ((xx(j+2) - xx(j+1)) / (h**2 * v_jm1))
    end do

    fval(2) = -fval(1)
    dfdx(2,:) = - dfdx(1,:)
end subroutine func1

end module mma_simcomp

