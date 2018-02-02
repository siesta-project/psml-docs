program show_psml

!
! Example of extraction of information from PSML file
! This is quite a "low-level" example, in which the
! program queries all the contents and prints the information.
!
! Examples closer to the use cases in electronic-structure
! codes are in preparation
!
! No examples of extraction of "j" or "s" values yet.
!
use m_psml
use m_getopts

implicit none

integer, parameter :: dp = selected_real_kind(10,100)
character(len=1), dimension(0:4) :: sym = (/ "s", "p", "d", "f", "g" /)

type(ps_t)   :: ps
type(ps_annotation_t)  :: annot

      character(len=200) :: filename
      logical            :: debug, plot, plot_raw
      character(len=200) :: opt_arg, mflnm, ref_line
      character(len=10)  :: opt_name 
      character(len=10)  :: interp_type
      integer            :: interp_quality
      integer :: nargs, iostat, n_opts, nlabels, iorb, ikb
      
      integer :: i, j, l, n, num, nfun, set, seq, ncont
      real(dp) :: ekb, rc, jval

      real(dp), allocatable, dimension(:) :: r
      real(dp) :: Rmax, delta, val
      integer  :: npts, ir
      character(len=50) :: outname,  target_set, set_str
      integer, allocatable   :: setidx(:), idxl(:)

      integer :: npots, nprojs, nwfs
      real(dp), allocatable   :: raw_r(:), raw_data(:)
      type(ps_radfunc_t) :: rfunc

      character(len=10) :: relativity
      logical           :: core_corrections

      character(len=100) :: creator, date
      integer            :: prov_depth, ninput
      character(len=100) :: uuid, version, namespace
      character(len=100) :: name, type
      real(dp)           :: weight
      integer            :: code
      logical            :: has_chlocal
      logical            :: use_effective_range = .true.
      
      ! There is a copy of dpnint1 inside the library
      ! This external declaration is to support calls
      ! to ps_SetInterpolator
      external :: interpolate_drh

      ! Uncomment the following when a new interpolator
      ! is added to the source tree. The one used previously
      ! was removed due to license issues.
!!      external :: interpolate_other  

!
!     Process options
!
      n_opts = 0
      debug = .false.
      plot = .false.
      plot_raw = .false.

      Rmax = 4.0_dp   ! default maximum radius for plots
      npts = 200      ! default number of points

      interp_type = "drh"
      interp_quality = 7  ! this is npoint=2 in OTHER interpolator
      do
         call getopts('dR:n:pi:q:tr',opt_name,opt_arg,n_opts,iostat)
         if (iostat /= 0) exit
         select case(opt_name)
           case ('d')
              debug = .true.
           case ('t')
              use_effective_range = .false.
           case ('p')
              plot = .true.
           case ('r')
              plot_raw = .true.
           case ('R')
              read(opt_arg,*) Rmax
           case ('n')
              read(opt_arg,*) npts
           case ('i')
              read(opt_arg,*) interp_type
           case ('q')
              read(opt_arg,*) interp_quality
           case ('?',':')
             write(0,*) "Invalid option: ", opt_arg(1:1)
             write(0,*) "Usage: show_psml [ -d -p -r -R Rmax -n npts ] FILE"
             write(0,*) " -d for debug output"
             write(0,*) " -p for output of evaluated data"
             write(0,*) " -r for output of raw data"
             write(0,*) " -R Rmax : maximum range of output grid (def: 4 bohr)"
             write(0,*) " -n npts : number of points of output grid (def: 200)"
             write(0,*) " -i type : interpolator type (drh) ('other' not implemented)"
             write(0,*) " -q nq   : interpolator quality (npoly for drh)"
             STOP
          end select
       enddo

       nargs = command_argument_count()
       nlabels = nargs - n_opts + 1
       if (nlabels /= 1)  then
             write(0,*) "Usage: show_psml [ -d -p -r -R Rmax -n npts ] FILE"
             write(0,*) " -d for debug output"
             write(0,*) " -p for output of evaluated data"
             write(0,*) " -r for output of raw data"
             write(0,*) " -R Rmax : maximum range of output grid (def: 4 bohr)"
             write(0,*) " -n npts : number of points of output grid (def: 200)"
             write(0,*) " -i type : interpolator type (drh) ('other' not implemented)"
             write(0,*) " -q nq   : interpolator quality (npoly for drh)"
          STOP
       endif

       call get_command_argument(n_opts,value=filename,status=iostat)
       if (iostat /= 0) then
          STOP "Cannot get filename"
       endif
!
call ps_destroy(ps)
print "(a)", "Processing: " // trim(filename)
call psml_reader(filename,ps,debug=debug)

call ps_SetEvaluatorOptions(use_effective_range=use_effective_range)
#ifdef __NO_PROC_POINTERS__
if (trim(interp_type)=="other") then
   print "(a)", "No support for other interpolators if no proc pointers"
   STOP 
endif
call ps_SetEvaluatorOptions(quality_index=interp_quality,debug=debug)
print "(a,i3)", "Using DRH (default) interpolator with npoly:", interp_quality

#else

if (trim(interp_type)=="other") then
!   call ps_SetInterpolator(interpolate_other,interp_quality)
!   print "(a,i3)", "Using OTHER interpolator with quality index:",interp_quality
   print "(a)", "Please implement 'other' interpolator first"
   STOP 
else if (trim(interp_type)=="drh") then
   call ps_SetEvaluatorOptions(custom_interpolator=interpolate_drh,&
                               quality_level=interp_quality,debug=debug)
   print "(a,i3)", "Using DRH interpolator with npoly:",interp_quality
else
   STOP "unknown interpolator"
endif

#endif

call ps_RootAttributes_Get(ps,uuid=uuid,version=version,&
                 namespace=namespace)
print "(a,1x,a,1x,a)", "PSML version, uuid:", version, uuid
print "(a,1x,a,1x,a)", "PSML namespace:", trim(namespace)

prov_depth = ps_Provenance_Depth(ps)
do i = prov_depth, 1, -1
   call ps_Provenance_Get(ps,i,creator=creator,date=date,annotation=annot,&
        number_of_input_files=ninput)
   print "(a,1x,i1,1x,a,1x,a,1x,i1)", "Prov: (depth)", i, trim(creator), trim(date), &
        ninput
   call print_annotation(annot)
end do
!
! Set up our grid
allocate(r(npts))
delta = Rmax / (npts - 1)
do ir = 1, npts
   r(ir) = (ir-1)*delta
enddo

call ps_PseudoAtomSpec_Get(ps,relativity=relativity,&
     core_corrections=core_corrections,annotation=annot)
call print_annotation(annot)

print *, "Relativity: ", trim(relativity)

call ps_ExchangeCorrelation_Get(ps,annotation=annot,n_libxc_functionals=nfun)
print *, "Exchange-correlation info:"

print *, "Number of functionals: ", nfun
do i = 1, nfun
   call ps_LibxcFunctional_Get(ps,i, &
           name=name,code=code,weight=weight,type=type)
  print "(a,1x,a,1x,a,1x,i4,f10.6)", "name, type, id, weight:", &
    trim(name), trim(type), code, weight
enddo
call print_annotation(annot)
!
! ================================================
!
call ps_Potential_Filter(ps,set=SET_ALL,indexes=setidx,number=npots)
print *, "There are ", npots, " semi-local pots"

do i = 1, npots
   call ps_Potential_Get(ps,setidx(i), &
                         l=l,n=n,j=jval,rc=rc,set=set,annotation=annot,func=rfunc)
   set_str = str_of_set(set)
   if (set == SET_LJ) then
      print "(a,1x,i1,a1,1x,f3.1,f10.3)", trim(set_str), n, sym(l), jval, rc
   else
      print "(a,1x,i1,a1,f10.3)", trim(set_str), n, sym(l), rc
   endif

   if (plot) then
   write(outname,"(a,i0,a1)") "Vsl." //trim(set_str) // ".", n, sym(l)
   open(4,file=outname,form="formatted",status="unknown",&
        action="write",position="rewind")
   do ir = 1, npts
      val = ps_Potential_Value(ps,setidx(i),r(ir))
      write(4,*) r(ir), val
   enddo
   close(4)
   endif
    if (plot_raw) then
      call ps_GetRawData(rfunc,raw_r,raw_data)
      write(outname,"(a)") trim(outname) // ".raw"
      open(4,file=outname,form="formatted",status="unknown",&
           action="write",position="rewind")
      do ir = 1, size(raw_data)
         write(4,*) raw_r(ir), raw_data(ir)
      enddo
      close(4)
    endif

enddo

if (.not. ps_HasLocalPotential(ps)) then
   print *, "This file does not have a local potential"
else
   call ps_LocalPotential_Get(ps,annotation=annot,type=type,&
                              has_local_charge=has_chlocal,func=rfunc)
   call print_annotation(annot)
   print *, "Vlocal type: " // trim(type)
   print *, "Has Local Charge: ", has_chlocal
   
   if (plot) then
   write(outname,"(a)") "Vlocal"
   open(4,file=outname,form="formatted",status="unknown",&
        action="write",position="rewind")
   do ir = 1, npts
      val = ps_LocalPotential_Value(ps,r(ir))
      write(4,*) r(ir), val
   enddo
   close(4)
   endif
    if (plot_raw) then
      call ps_GetRawData(rfunc,raw_r,raw_data)
      write(outname,"(a)") trim(outname) // ".raw"
      open(4,file=outname,form="formatted",status="unknown",&
           action="write",position="rewind")
      do ir = 1, size(raw_data)
         write(4,*) raw_r(ir), raw_data(ir)
      enddo
      close(4)
    endif
endif

if (.not. ps_HasProjectors(ps)) then
   print *, "This file does not have non-local projectors"
else

   call ps_Projector_Filter(ps,set=SET_ALL,indexes=setidx,number=nprojs)
   print *, "There are ", nprojs, " projectors"

   print "(a20,a8,a4,a12)", "set", "l (j)", "seq", "Ekb"
   do l = 0, 4
      call ps_Projector_Filter(ps,indexes_in=setidx,&
                   l=l,indexes=idxl,number=num)

   do j = 1, num
      call ps_Projector_Get(ps,idxl(j),ekb=ekb,set=set,seq=seq,j=jval,func=rfunc)
      set_str = str_of_set(set)
      if (set == SET_LJ) then
         print "(a20,i4,1x,f3.1,i4,f12.6)", set_str, l, jval, seq, ekb
         write(outname,"(a,i0,a1,f3.1,a1,i0)") &
                   "Proj." //trim(set_str) // ".", &
                    l, ".", jval,".", seq
      else
         print "(a20,i4,i4,f12.6)", set_str, l, seq, ekb
         write(outname,"(a,i0,a1,i0)") &
              "Proj." //trim(set_str) // ".", l, ".", seq
      endif

      if (plot) then

      open(4,file=outname,form="formatted",status="unknown",&
           action="write",position="rewind")
      do ir = 1, npts
         val = ps_Projector_Value(ps,idxl(j),r(ir))
         write(4,*) r(ir), val
      enddo
      close(4)
    if (plot_raw) then
      call ps_GetRawData(rfunc,raw_r,raw_data)
      write(outname,"(a)") trim(outname) // ".raw"
      open(4,file=outname,form="formatted",status="unknown",&
           action="write",position="rewind")
      do ir = 1, size(raw_data)
         write(4,*) raw_r(ir), raw_data(ir)
      enddo
      close(4)
    endif
   endif

   enddo
enddo
endif

! Valence charge density
call ps_ValenceCharge_Get(ps,annotation=annot)
call print_annotation(annot)

   if (plot) then
      write(outname,"(a)") "Valence.charge"
      open(4,file=outname,form="formatted",status="unknown",&
           action="write",position="rewind")
      do ir = 1, npts
         val = ps_ValenceCharge_Value(ps,r(ir))
         write(4,*) r(ir), val
      enddo
      close(4)
   endif

! Pseudo-Core charge

   if (core_corrections) then
      call ps_CoreCharge_Get(ps,annotation=annot,&
                             nderivs=ncont)
      call print_annotation(annot)
      print *, "No of continuous derivs: ", ncont
      if (plot) then
         write(outname,"(a)") "Core.charge"
         open(4,file=outname,form="formatted",status="unknown",&
              action="write",position="rewind")
         do ir = 1, npts
            val = ps_CoreCharge_Value(ps,r(ir))
            write(4,"(1p,2e21.13)") r(ir), val
         enddo
         close(4)
      endif
   endif
   
call ps_PseudoWaveFunctions_Filter(ps,set=SET_ALL,indexes=setidx,number=nwfs)
print *, "There are ", nwfs, " pseudo wave-functions"
!call display_annotation(ps,"pseudo-wavefunctions")
!
do i = 1, nwfs
   call ps_PseudoWf_Get(ps,setidx(i),l=l,n=n,set=set,j=jval)
   set_str = str_of_set(set)
   if (set == SET_LJ) then
      print "(a,1x,i1,a1,1x,f3.1)", trim(set_str), n, sym(l), jval
   else
      print "(a,1x,i1,a1)", trim(set_str), n, sym(l)
   endif

   if (plot) then
   write(outname,"(a,i0,a1)") "Pswf." //trim(set_str) // ".", n, sym(l)
   open(4,file=outname,form="formatted",status="unknown",&
        action="write",position="rewind")
   do ir = 1, npts
      val = ps_PseudoWf_Value(ps,setidx(i),r(ir))
      write(4,*) r(ir), val
   enddo
   close(4)
   endif
enddo
   
call ps_destroy(ps)

end program show_psml

