module test
contains

  subroutine split_string(string, word_list)
    use iso_fortran_env, only: iostat_end
    implicit none

    character(len=*), intent(in) :: string
    character(len=255), dimension(:), allocatable, intent(inout):: word_list

    character(len=255), dimension(:), allocatable:: tmp


    integer :: nr_words
    integer :: status, i

    nr_words = 1
    do
       if (allocated(word_list)) then
          deallocate(word_list)
       end if
       allocate(character(len=255) :: word_list(nr_words))

       read(string, *, iostat=status) word_list

       if (status == iostat_end) exit

       nr_words = nr_words + 1
    end do

    ! Delete any empty strings
    tmp = word_list
    do i = 1, size(tmp)
       tmp(i) = trim(tmp(i))
    end do

    deallocate(word_list)
    allocate(character(len=255) :: word_list(count(tmp /= '')))
    do i = 1, size(tmp)
       if (tmp(i) /= '') then
          word_list(i) = tmp(i)
       end if
    end do

  end subroutine split_string



! This is a simple test program used to develop our vtk reader system.
  subroutine read_binary_vtk_file(filename)
    use iso_fortran_env, only: iostat_end

    implicit none
    character(len=*), intent(in) :: filename
    integer :: unit, status
    integer :: nr_words

    real, dimension(:,:,:), allocatable :: data

    ! File reading temporary variables
    character(len=255) :: line
    character(len=255), dimension(:), allocatable :: word_list
    character(len=255) :: word

    ! VTK File header variables
    real :: vtk_version
    logical :: is_binary
    character(len=255) :: data_name, data_type

    integer:: i

    ! Open the file
    open(newunit=unit, file=filename, status='old', action='read')

    ! ======================================================================== !
    ! Read Header, this will also determine the VTK version

    read(unit, '(A23,F3.2)') word, vtk_version
    if (word /= '# vtk DataFile Version') then
       print *, 'Error: ' // &
         'VTK File error in header, expected "# vtk DataFile Version"'
       stop
    end if

    ! Read the title
    read(unit, '(A255)') data_name

    ! Read the file type
    read(unit, '(A)') line
    select case (word_list(1))
      case ('BINARY')
       is_binary = .true.
      case ('ASCII')
       is_binary = .false.
      case default
       print *, 'Error: ' // &
         'VTK File error in header line 3, expected "BINARY" or "ASCII"'
       stop
    end select

    print *, "Reading VTK file: ", filename
    print *, "VTK version: ", vtk_version
    print *, "Data name: ", data_name
    if (is_binary) then
       print *, "Data type: binary"
    else
       print *, "Data type: ascii"
    end if


    !  ! ======================================================================== !
    !  ! Read the actual data

    !  read_lines: do
    !     read(unit, '(A)', iostat=status) line
    !     if (status == iostat_end) exit read_lines

    !     ! Split the line into words
    !     call split_string(line, word_list)
    !     nr_words = size(word_list)

    !     select case (word_list(1))


    !        ! Read the data types in the file
    !       case ('DATASET')

    !        do i = 1, nr_words
    !           print *, trim(word_list(i))
    !        end do
    !       case ('DIMENSIONS')

    !        do i = 1, nr_words
    !           print *, trim(word_list(i))
    !        end do
    !       case ('ORIGIN')

    !        do i = 1, nr_words
    !           print *, trim(word_list(i))
    !        end do
    !       case ('SPACING')

    !        do i = 1, nr_words
    !           print *, trim(word_list(i))
    !        end do
    !       case ('POINT_DATA')

    !        do i = 1, nr_words
    !           print *, trim(word_list(i))
    !        end do
    !       case ('SCALARS')

    !        do i = 1, nr_words
    !           print *, trim(word_list(i))
    !        end do
    !       case ('LOOKUP_TABLE')

    !        do i = 1, nr_words
    !           print *, trim(word_list(i))
    !        end do
    !        !  case default
    !        !   print *, trim(word(1))

    !     end select
    !  end do read_lines


!   ! Read the grid dimensions
!   read(unit, '(A,I5,I5,I5)') line, nx, ny, nz

!   ! Allocate the data array
!   allocate(data(nx, ny, nz))

!   ! Skip the data description line
!   read(unit, '(A)') line

!   ! Read the data
!   do k = 1, nz
!      do j = 1, ny
!         do i = 1, nx
!            read(unit) data(i, j, k)
!         end do
!      end do
!   end do

    ! Close the file
    close(unit)
  end subroutine read_binary_vtk_file

end module test

program reader_test
  use test
  implicit none

  call read_binary_vtk_file('data_local/petsc_static_mixer_lowest.vtk')

end program reader_test
