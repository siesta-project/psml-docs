#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module external_interfaces
  !
  interface
     ! Called to terminate the program, printing a message
     subroutine psml_die_interf(str)
       character(len=*), intent(in)   :: str
     end subroutine psml_die_interf
  end interface

  ! This initialization to non-null() is a F2008 feature
  procedure(psml_die_interf),pointer ::   &
      psml_die => simple_die_routine

  public :: set_die_routine
  public :: psml_die

CONTAINS

  subroutine set_die_routine(func)

    procedure(psml_die_interf) :: func

    psml_die => func

  end subroutine set_die_routine

  subroutine simple_die_routine(str)
    character(len=*), intent(in)   :: str
    write(0,'(a,a)') "[libpsml default error handler]: " // trim(str)
    write(6,'(a,a)') "[libpsml default error handler]: " // trim(str)
    write(6,'(a)') "Use ps_set_error_handler(custom_handler) to configure"
    stop
  end subroutine simple_die_routine

end module external_interfaces
