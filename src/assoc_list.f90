#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module m_psml_assoc_list
! First version, with fixed (initial) length,
! and fixed-length fields.
! Alberto Garcia, Sept 2014
!
!-----------------------------------------------------------
type, public :: assoc_list_t
  private
  integer                               :: nslots = 0
  integer                               :: nitems = 0
  character(len=200), allocatable       :: key(:)
  character(len=1000), allocatable      :: value(:)
end type assoc_list_t

public :: assoc_list_init
public :: assoc_list_reset
public :: assoc_list_insert
public :: assoc_list_nitems
public :: assoc_list_get_key
public :: assoc_list_get_value
public :: assoc_list_print

! Now, define and export an EMPTY_ASSOC_LIST
! This feature might be deprecated in future versions. Client programs can
! use 'assoc_list_reset' to generate an empty association list.


! NOTE: A module variable of a type with default initialization must
! have the SAVE attribute by the Fortran 2003 or earlier standards
! The 'protected' attribute forbids changing the value of the variable
! from outside the module.

! NOTE: A 'parameter' could be used instead of a variable, with explicit initialization
! of the form:
! type(assoc_list_t), parameter :: PAR = assoc_list_type_t(0,0,null(),null())
! ... but compiler support for this might be spotty.

type(assoc_list_t), protected, save   :: EMPTY_ASSOC_LIST 


interface assoc_list_get_value
  module procedure assoc_list_get_value_by_index
  module procedure assoc_list_get_value_of_key
end interface

CONTAINS

subroutine assoc_list_init(a,n,stat)
type(assoc_list_t), intent(inout) :: a
integer, intent(in)               :: n
integer, intent(out)              :: stat

  if (allocated(a%key)) then
     deallocate(a%key)
  endif
  if (allocated(a%value)) then
     deallocate(a%value)
  endif
  a%nslots = n
  a%nitems = 0
  allocate(a%key(n),a%value(n),stat=stat)

end subroutine assoc_list_init

subroutine assoc_list_reset(a)
type(assoc_list_t), intent(inout) :: a

  if (allocated(a%key)) then
     deallocate(a%key)
  endif
  if (allocated(a%value)) then
     deallocate(a%value)
  endif
  a%nslots = 0
  a%nitems = 0

end subroutine assoc_list_reset

subroutine assoc_list_insert(a,key,value,stat)
type(assoc_list_t), intent(inout) :: a
character(len=*), intent(in)      :: key, value
integer, intent(out)              :: stat

integer :: i
character(len=1000), allocatable  :: b(:)

if (.not. allocated(a%key)) then
   call assoc_list_init(a,4,stat)
   if (stat /= 0) return
endif

! Replace if key exists already
do i = 1, a%nitems
   if (a%key(i) == key) then
      a%value(i) = value
      stat = 0
      return
   endif
enddo
!
! Add at the end
a%nitems = a%nitems + 1
if (a%nitems > a%nslots)  then
   ! Enlarge a%data by 4 slots
   !
   allocate(b(a%nslots))
   !
   b = a%key
   deallocate(a%key)
   allocate(a%key(a%nslots + 4))
   a%key(1:a%nslots) = b
   !
   b = a%value
   deallocate(a%value)
   allocate(a%value(a%nslots + 4))
   a%value(1:a%nslots) = b
   deallocate(b)
   !
   a%nslots = size(a%key)
endif
i = a%nitems
a%key(i) = key
a%value(i) = value
stat = 0

end subroutine assoc_list_insert

function assoc_list_nitems(a) result(n)
type(assoc_list_t), intent(in)    :: a
integer                           :: n

n = a%nitems
end function assoc_list_nitems

subroutine assoc_list_get_key(a,i,key,stat)
type(assoc_list_t), intent(in)    :: a
integer, intent(in)               :: i
character(len=*), intent(out)     :: key
integer, intent(out)              :: stat

if (i > a%nitems) then
   stat = -1
   return
endif
key = a%key(i)
stat = 0
end subroutine assoc_list_get_key

subroutine assoc_list_get_value_of_key(a,key,value,stat)
type(assoc_list_t), intent(in) :: a
character(len=*), intent(in)      :: key
character(len=*), intent(out)     :: value
integer, intent(out)              :: stat

integer :: i

do i = 1, a%nitems
   if (a%key(i) == key) then
      value = a%value(i) 
      stat = 0
      return
   endif
enddo
stat = -1
end subroutine assoc_list_get_value_of_key

subroutine assoc_list_get_value_by_index(a,i,value,stat)
type(assoc_list_t), intent(in) :: a
integer, intent(in)            :: i
character(len=*), intent(out)     :: value
integer, intent(out)              :: stat

if (i <= a%nitems) then
   value = a%value(i) 
   stat = 0
else
   stat = -1
   value =""
endif

end subroutine assoc_list_get_value_by_index

subroutine assoc_list_print(a)
type(assoc_list_t), intent(in) :: a

integer :: i

if (a%nitems > 0) then
   print "(a)", "---------------------"
   do i = 1, a%nitems
      print "(3x,a,a,a)", trim(a%key(i)), " : ", trim(a%value(i))
   enddo
   print "(a)", "---------------------"
endif
end subroutine assoc_list_print

end module m_psml_assoc_list
