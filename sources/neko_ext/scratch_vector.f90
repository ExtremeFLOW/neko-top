! Copyright (c) 2018-2023, The Neko Authors
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
!
!> Defines a registry for storing and requesting temporary vectors
!! This can be used when you have a function that will be called
!! often and you don't want to create temporary vectors (work arrays) inside
!! it on each call.
module scratch_vector
  use vector, only : vector_t, vector_ptr_t
  implicit none
  private


  type, public :: scratch_vector_t
     !> list of scratch vectors
     type(vector_ptr_t), private, allocatable :: vectors(:)
     !> Tracks which vectors are used
     logical, private, allocatable :: inuse(:)
     !> number of registered vectors
     integer, private :: nvectors
     !> number of vectors in use
     integer, private :: nvectors_inuse
     !> the registry_size the vectors array is increased by upon reallocation
     integer, private :: expansion_size
     !> Dofmap
     integer :: vector_size
   contains
     procedure, private, pass(this) :: expand
     !> destructor
     procedure, pass(this) :: free => scratch_registry_free
     !> getter for nvectors
     procedure, pass(this) :: get_nvectors
     !> getter for nvectors_inuse
     procedure, pass(this) :: get_nvectors_inuse
     !> getter for expansion_size
     procedure, pass(this) :: get_expansion_size
     !> return registry_size of allocated vectors
     procedure, pass(this) :: get_size
     !> get value of inuse for a given index
     procedure, pass(this) :: get_inuse
     !> get a new scratch vector
     procedure, pass(this) :: request_vector
     procedure, pass(this) :: relinquish_vector_single
     procedure, pass(this) :: relinquish_vector_multiple
     !> free a vector for later reuse
     generic :: relinquish_vector => relinquish_vector_single, &
          relinquish_vector_multiple
  end type scratch_vector_t

  interface scratch_vector_t
     procedure :: init
  end interface scratch_vector_t

  !> Global scratch registry
  type(scratch_vector_t), public, target :: neko_scratch_vector

contains

  !> Constructor, optionally taking initial registry and expansion
  !! registry_size as argument
  type(scratch_vector_t) function init(vector_size, registry_size, expansion_size) result(this)
    integer, intent(in) :: vector_size
    integer, optional, intent(in) :: registry_size
    integer, optional, intent(in) :: expansion_size
    integer :: i

    this%vector_size = vector_size

    if (present(registry_size)) then
       allocate (this%vectors(registry_size))
       do i= 1, registry_size
          allocate(this%vectors(i)%ptr)
       end do
       allocate (this%inuse(registry_size))
    else
       allocate (this%vectors(10))
       allocate (this%inuse(10))
    end if

    this%inuse(:) = .false.
    if (present(expansion_size)) then
       this%expansion_size = expansion_size
    else
       this%expansion_size = 10
    end if

    this%nvectors = 0
    this%nvectors_inuse = 0
  end function init

  !> Destructor
  subroutine scratch_registry_free(this)
    class(scratch_vector_t), intent(inout):: this
    integer :: i

    if (allocated(this%vectors)) then
       do i=1, this%nvectors
          call this%vectors(i)%ptr%free()
          deallocate(this%vectors(i)%ptr)
       end do

       deallocate(this%vectors)
       deallocate(this%inuse)
    end if

  end subroutine scratch_registry_free


  !> Get the number of vectors stored in the registry
  pure function get_nvectors(this) result(n)
    class(scratch_vector_t), intent(in) :: this
    integer :: n

    n = this%nvectors
  end function get_nvectors

  pure function get_nvectors_inuse(this) result(n)
    class(scratch_vector_t), intent(in) :: this
    integer :: n, i

    n = 0
    do i=1,this%get_size()
       if (this%inuse(i)) n = n + 1
    end do
  end function get_nvectors_inuse

  !> Get the registry_size of the vectors array
  pure function get_size(this) result(n)
    class(scratch_vector_t), intent(in) :: this
    integer :: n

    n = size(this%vectors)
  end function get_size

  !> Get the expansion registry_size
  pure function get_expansion_size(this) result(n)
    class(scratch_vector_t), intent(in) :: this
    integer :: n

    n = this%expansion_size
  end function get_expansion_size

  subroutine expand(this)
    class(scratch_vector_t), intent(inout) :: this
    type(vector_ptr_t), allocatable :: temp(:)
    logical, allocatable :: temp2(:)
    integer :: i

    allocate(temp(this%get_size() + this%expansion_size))
    temp(1:this%nvectors) = this%vectors(1:this%nvectors)

    do i=this%nvectors +1, size(temp)
       allocate(temp(i)%ptr)
    enddo

    call move_alloc(temp, this%vectors)

    allocate(temp2(this%get_size() + this%expansion_size))
    temp2(1:this%nvectors) = this%inuse(1:this%nvectors)
    temp2(this%nvectors+1:) = .false.
    this%inuse = temp2
  end subroutine expand


  !> Get a vector from the registry by assigning it to a pointer
  subroutine request_vector(this, f, index)
    class(scratch_vector_t), target, intent(inout) :: this
    type(vector_t), pointer, intent(inout) :: f
    integer, intent(inout) :: index !< The index of the vector in the inuse array


    associate(nvectors => this%nvectors, nvectors_inuse => this%nvectors_inuse)

      do index=1,this%get_size()
         if (.not. this%inuse(index)) then

            if (.not. allocated(this%vectors(index)%ptr%x)) then
               call this%vectors(index)%ptr%init(this%vector_size)
               nvectors = nvectors + 1
            end if
            f => this%vectors(index)%ptr
            this%inuse(index) = .true.
            this%nvectors_inuse = this%nvectors_inuse + 1
            return
         end if
      end do
      ! all existing vectors in use, we need to expand to add a new one
      index = nvectors +1
      call this%expand()
      nvectors = nvectors + 1
      nvectors_inuse = nvectors_inuse + 1
      this%inuse(nvectors) = .true.
      call this%vectors(nvectors)%ptr%init(this%vector_size)
      f => this%vectors(nvectors)%ptr

    end associate
  end subroutine request_vector

  !> Relinquish the use of a vector in the registry
  subroutine relinquish_vector_single(this, index)
    class(scratch_vector_t), target, intent(inout) :: this
    integer, intent(inout) :: index !< The index of the vector to free

    this%inuse(index) = .false.
    this%nvectors_inuse = this%nvectors_inuse - 1
  end subroutine relinquish_vector_single

  subroutine relinquish_vector_multiple(this, indices)
    class(scratch_vector_t), target, intent(inout) :: this
    integer, intent(inout) :: indices(:) !< The indices of the vector to free
    integer :: i

    do i=1, size(indices)
       this%inuse(indices(i)) = .false.
    end do
    this%nvectors_inuse = this%nvectors_inuse - size(indices)
  end subroutine relinquish_vector_multiple

  logical function get_inuse(this, index)
    class(scratch_vector_t), target, intent(inout) :: this
    integer, intent(inout) :: index !< The index of the vector to check

    get_inuse = this%inuse(index)
  end function get_inuse

end module scratch_vector
