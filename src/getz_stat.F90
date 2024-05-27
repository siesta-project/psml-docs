program getz

!
  ! Example of extraction of information from PSML file
  ! Just the atomic number, for scripts
!
use m_psml
use m_getopts

external :: custom_psml_die

! Example of use of the provided kind specification
real(kind=ps_real_kind) :: atnum

character(len=1), dimension(0:4) :: sym = (/ "s", "p", "d", "f", "g" /)

type(ps_t)   :: ps

      character(len=200) :: filename
      logical            :: debug
      character(len=200) :: opt_arg, mflnm, ref_line
      character(len=10)  :: opt_name 
      integer :: nargs, iostat, n_opts, nlabels
      
      integer :: i, j, l, n, num, nfun, set, seq
!
!     Process options
!
      n_opts = 0
      debug = .false.

      do
         call getopts('d',opt_name,opt_arg,n_opts,iostat)
         if (iostat /= 0) exit
         select case(opt_name)
           case ('d')
              debug = .true.
           case ('?',':')
             write(0,*) "Invalid option: ", opt_arg(1:1)
             write(0,*) "Usage: getz [-d]  FILE"
             write(0,*) " -d   : debug"
             STOP
          end select
       enddo

       nargs = command_argument_count()
       nlabels = nargs - n_opts + 1
       if (nlabels /= 1)  then
          write(0,*) "Usage: getz FILE"
          write(0,*) " -d for debug output"
          STOP
       endif

       call get_command_argument(n_opts,value=filename,status=iostat)
       if (iostat /= 0) then
          STOP "Cannot get filename"
       endif
!
       call ps_set_error_handler(custom_psml_die)
       
if (debug) print "(a)", "Processing: " // trim(filename)
call psml_reader(filename,ps,debug=debug,stat=iostat)
if (iostat /=0) then
   print *, "PSML reader error. Iostat: ", iostat
   stop 
endif

call ps_PseudoAtomSpec_Get(ps,atomic_number=atnum)
print *,  int(atnum)

end program getz

