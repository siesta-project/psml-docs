module sets_m
!
! Utilities to handle set information in arrays
!
implicit none
private
integer, parameter :: MAX_NSETS = 10

type alloc_t
   integer, allocatable :: ind(:)
end type alloc_t

type, public :: set_info_t
   private
   integer :: nsets
   integer :: set(MAX_NSETS)
   integer :: nels(MAX_NSETS)
   type(alloc_t) :: indexes(MAX_NSETS)
end type set_info_t

public :: sort_sets, display
public :: nsets, set_id, set_indexes
CONTAINS

subroutine clear(set_info)
type(set_info_t),  intent(inout) :: set_info

integer :: j

set_info%set(:) = 0
set_info%nels(:) = 0
do j = 1, MAX_NSETS
   if (allocated(set_info%indexes(j)%ind)) then
      deallocate(set_info%indexes(j)%ind)
   endif
enddo
set_info%nsets = 0
end subroutine clear

subroutine display(set_info)
type(set_info_t),  intent(in) :: set_info

integer :: j, i

do j = 1, set_info%nsets
      print *, "Set: ", set_info%set(j)
      do i = 1, size(set_info%indexes(j)%ind)
         print *, "index: ", set_info%indexes(j)%ind(i)
      enddo
enddo
end subroutine display
!
!

subroutine sort_sets(n,a,set_info)
integer, intent(in) :: n
integer, intent(in) :: a(:)
type(set_info_t),  intent(inout) :: set_info

integer :: i, j, dim, nsets, current_set, set

call clear(set_info)

! First pass. Assumes bunched sets
current_set = 0
nsets = 0
do i = 1, n
   set = a(i)
   if (set /= current_set) then
      nsets = nsets + 1
      set_info%set(nsets) = set
      set_info%nels(nsets) = 1
      current_set = set
   else
      set_info%nels(nsets) = set_info%nels(nsets) + 1
   endif
enddo
set_info%nsets = nsets
! Allocate
do j = 1, nsets
   dim = set_info%nels(j)
   allocate(set_info%indexes(j)%ind(dim))
enddo
! Second pass.
current_set = 0
nsets = 0
do i = 1, n
   set = a(i)
   if (set /= current_set) then
      nsets = nsets + 1
      if (nsets > MAX_NSETS) stop "set overflow"
      set_info%nels(nsets) = 1
      set_info%indexes(nsets)%ind(1) = i
      current_set = set
   else
      set_info%nels(nsets) = set_info%nels(nsets) + 1
      set_info%indexes(nsets)%ind(set_info%nels(nsets)) = i
   endif
enddo
end subroutine sort_sets

function nsets(set_info) result(n)
  type(set_info_t), intent(in) :: set_info
  integer :: n

  n = set_info%nsets

end function nsets

subroutine set_indexes(set_info,j,idx)
  type(set_info_t), intent(in) :: set_info
  integer, intent(in) :: j
  integer, allocatable, intent(inout) :: idx(:)

  integer :: n

  if (j > set_info%nsets) then
     n = 0
  else
     n = set_info%nels(j)
  endif

  if (allocated(idx)) deallocate(idx)
  allocate(idx(n))
  idx(:) = set_info%indexes(j)%ind(:)

end subroutine set_indexes

function set_id(set_info,j) 
  type(set_info_t), intent(in) :: set_info
  integer, intent(in) :: j
  integer :: set_id

  if (j > set_info%nsets) then
     set_id = 0
  else
     set_id = set_info%set(j)
  endif
end function set_id

end module sets_m

#ifdef __TEST__
program ts
  use sets_m
  type(set_info_t) :: set_info
  integer :: a(20)
  integer :: n

  print *, "Enter n"
  read *, n
  print *, "enter a"
  read *, a(1:n)

  call sort_sets(n,a,set_info)
  call display(set_info)
end program ts
#endif
