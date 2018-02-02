program normalize

  !
  ! Parses a PSML file and dumps the contents of the resulting
  ! ps_t object in a PSML 1.1 file.
  !
  ! If it detects a 1.0 PSML file, it inserts
  ! a new provenance element, but keeps the original uuid.
  !
  ! It will also add 'record-number' attributes to provenance elements.
  !
use m_psml
use m_getopts

integer, parameter :: dp = selected_real_kind(10,100)

type(ps_t)   :: ps

      character(len=200) :: filename, output_filename
      logical            :: debug
      character(len=200) :: opt_arg, mflnm, ref_line
      character(len=10)  :: opt_name 
      integer :: nargs, iostat, n_opts, nlabels
      
      integer :: i, j, l, n, num, nfun, set, seq
      character(len=20) :: date, version
!
!     Process options
!
      n_opts = 0
      debug = .false.
      output_filename = "PSML_DUMP"
      do
         call getopts('do:',opt_name,opt_arg,n_opts,iostat)
         if (iostat /= 0) exit
         select case(opt_name)
           case ('d')
              debug = .true.
           case ('o')
              read(opt_arg,*) output_filename
           case ('?',':')
             write(0,*) "Invalid option: ", opt_arg(1:1)
             write(0,*) "Usage: test_dump [-o output_file] PSML_FILE"
             STOP
          end select
       enddo

       nargs = command_argument_count()
       nlabels = nargs - n_opts + 1
       if (nlabels /= 1)  then
             write(0,*) "Invalid option: ", opt_arg(1:1)
             write(0,*) "Usage: test_dump [-o output_file] PSML_FILE"
             STOP
       endif

       call get_command_argument(n_opts,value=filename,status=iostat)
       if (iostat /= 0) then
          STOP "Cannot get filename"
       endif
!
call ps_destroy(ps)
if (debug) print "(a)", "Processing: " // trim(filename)
call date_and_time(date)
call psml_reader(filename,ps,debug=debug)
call ps_RootAttributes_Get(ps,version=version)

call ps_RootAttributes_Set(ps,version="1.1",&
     namespace="http://esl.cecam.org/PSML/ns/1.1")

if (trim(version) == "1.0") then
   if (nitems_annotation(ps%annotation)>0) then
      call ps_Provenance_Add(ps,creator="1.0-to-1.1-conversion",&
           date=trim(date),  annotation=ps%annotation)
   else
      call ps_Provenance_Add(ps,creator="1.0-to-1.1-conversion",&
           date=trim(date), annotation=EMPTY_ANNOTATION)
   endif
endif

call ps_DumpToPSMLFile(ps,trim(output_filename))
write(0,"(a)") "Written PSML 1.1 file: " // trim(output_filename)
end program normalize

