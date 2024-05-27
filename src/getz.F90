program getz

!
! Example of PSML data extraction
! 
! + Tests of possible use of custom error handler, depending on the PSML library version
!
use m_psml
use m_getopts

#ifdef PSML_HAS_PP
external :: custom_psml_die
#endif

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
      integer :: version_number

      version_number = ps_GetLibPSMLVersion()
      print "(a,i5)", "libpsml Version: ", version_number
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
#ifdef PSML_HAS_PP       
       call ps_set_error_handler(custom_psml_die)
#endif

       if (debug) print "(a)", "Processing: " // trim(filename)

#ifdef PSML_READER_HAS_STAT
! Released versions with procedure pointers also have the 'stat' optional argument
call psml_reader(filename,ps,debug=debug,stat=iostat)
if (iostat /=0) then
   print *, "PSML reader error. Iostat: ", iostat
   stop 
endif
#else
call psml_reader(filename,ps,debug=debug)
#endif


call ps_PseudoAtomSpec_Get(ps,atomic_number=atnum)
print *,  int(atnum)

end program getz

