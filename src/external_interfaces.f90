#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module external_interfaces
!
! The user must provide external subroutines
! with these interfaces
!
public :: psml_die  
interface
   ! Called to terminate the program, printing a message
      subroutine psml_die(str)
      character(len=*), intent(in)   :: str
      end subroutine psml_die
end interface

end module external_interfaces
