#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module class_Grid

  use assoc_list, only: ps_annotation_t => assoc_list_t
  use assoc_list, only: ps_clean_annotation => assoc_list_reset
  
  implicit none

  character(len=*), parameter :: mod_name= "Grid"

  integer, parameter :: dp = selected_real_kind(10,100)

  !
  type Grid_
    integer :: refCount = 0
    character(len=36)   :: id = "null_id"
    !----------------------
    character(len=256)   :: name = "null Grid"
    integer                        :: npts = 0
    real(dp), pointer              :: grid_data(:) => null()
    type(ps_annotation_t)          :: annotation 
  end type Grid_

  type Grid
    type(Grid_), pointer :: data => null()
  end type Grid

  public  :: newGrid, print_type, valGrid, annotationGrid, sizeGrid

  interface print_type
    module procedure printGrid
  end interface

!========================
#define TYPE_NAME Grid
#include "basic_type.inc"
!========================

     subroutine delete_Data(gd)
      type(Grid_) :: gd
      if (associated(gd%grid_data)) then
         deallocate(gd%grid_data)
         gd%grid_data => null()
      endif
      call ps_clean_annotation(gd%annotation)
     end subroutine delete_Data


  subroutine newGrid(this,n,name)

   type(Grid), intent(inout)  :: this
   integer, intent(in)        :: n
   character(len=*), intent(in), optional  :: name

   integer :: stat

   ! We release the previous incarnation
   ! This means that we relinquish access to the previous
   ! memory location. It will be deallocated when nobody
   ! else is using it.

   call init(this)

   if (present(name)) then
      this%data%name = trim(name)
   else
      this%data%name = "(Grid from n)"
   endif

   allocate(this%data%grid_data(n))

   call tag_new_object(this)

  end subroutine newGrid

  function valGrid(this) result(p)
   type(Grid), intent(in)  :: this
   real(dp), pointer       :: p(:)

   nullify(p)
   p => this%data%grid_data
 end function valGrid

  function annotationGrid(this) result(p)
   type(Grid), intent(in)         :: this
   type(ps_annotation_t) , pointer   :: p

   nullify(p)
   p => this%data%annotation
 end function annotationGrid

 function sizeGrid(this) result(n)
   type(Grid), intent(in)  :: this
   integer                 :: n

   if (.not. initialized(this)) then
      n = 0
   else
      n = size(this%data%grid_data)
   endif
 end function sizeGrid

 subroutine printGrid(this)
   type(Grid), intent(in)  :: this

    integer :: n, m

   if (.not. associated(this%data)) then
      print "(a)", "Grid Not Associated"
      RETURN
   endif

    n = size(this%data%grid_data)

   print "(a,i0,a,i0,a,i0,a)", "  <grid:" // trim(this%data%name) // &
                               " n=",  n,   &
                               ", refcount: ", refcount(this),">"
 end subroutine printGrid

end module class_Grid
