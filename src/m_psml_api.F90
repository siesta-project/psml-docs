!+ libPSML API implementation
module m_psml_api
!+ graph: false
!+ Procedures to handle the PSML pseudopotential format.

use m_psml_core   ! For basic structures

use assoc_list, only: ps_annotation_t => assoc_list_t
use assoc_list, only: EMPTY_ANNOTATION => EMPTY_ASSOC_LIST

use external_interfaces, only: die => psml_die
use class_Grid

implicit none

integer, parameter    :: dp = selected_real_kind(14)
logical               :: global_debug = .false.
character(len=1), dimension(0:4) :: sym = (/ "s", "p", "d", "f", "g" /)

public :: ps_GetLibPSMLVersion
public :: ps_SetDebug
#ifndef  __NO_PROC_POINTERS__
public :: ps_SetInterpolator
#endif
public :: ps_SetInterpolatorQuality
!
! Accessor list
!
public :: ps_AtomicSymbol
public :: ps_AtomicLabel
public :: ps_AtomicNumber
public :: ps_ZPseudo
public :: ps_GenerationZval
public :: ps_PseudoFlavor

public :: ps_GetUUID
public :: ps_GetPSMLVersion
public :: ps_Creator
public :: ps_Date

public :: ps_GetAnnotation

public :: ps_NLibxcFunctionals
public :: ps_LibxcName
public :: ps_LibxcId
public :: ps_LibxcWeight

public :: ps_ValidLibxc
public :: ps_XCAnnotation

public :: ps_Relativity
public :: ps_IsSpinPolarized
public :: ps_HasCoreCorrections

public :: ps_NValenceShells
public :: ps_ValenceShellL
public :: ps_ValenceShellN
public :: ps_ValenceShellOccupation
!
! Semilocal potentials
!
public :: ps_HasSemilocalPotentials
public :: ps_Get_Potential_Indexes
public :: ps_Number_Of_Potentials
public :: ps_Potential_L
public :: ps_Potential_J
public :: ps_Potential_N
public :: ps_Potential_Rc
public :: ps_Potential_Set
public :: ps_Potential_Value
public :: ps_Potential_GetRawData
!
public :: ps_HasPsOperator
!
! Vlocal 
!
public :: ps_HasLocalPotential
public :: ps_LocalPotential_Value
public :: ps_LocalPotential_Type
public :: ps_LocalPotential_GetRawData
!
! Local Charge
!
public :: ps_HasLocalCharge
public :: ps_LocalCharge_Value
!
! Projectors
!
public :: ps_HasProjectors
public :: ps_Number_Of_Projectors
public :: ps_Get_Projector_Indexes
public :: ps_Get_Projector_Indexes_byL
public :: ps_Projector_L
public :: ps_Projector_J
public :: ps_Projector_Seq
public :: ps_Projector_Ekb
public :: ps_Projector_Type
public :: ps_Projector_Set
public :: ps_Projector_Value
public :: ps_Projector_GetRawData
!
! Pseudo-wave-functions
!
public :: ps_Number_Of_PseudoWfs
public :: ps_Get_PseudoWf_Indexes
public :: ps_PseudoWf_L
public :: ps_PseudoWf_J
public :: ps_PseudoWf_N
public :: ps_PseudoWf_Set
public :: ps_PseudoWf_Value
public :: ps_PseudoWf_GetRawData
!
public :: ps_ValenceCharge_Value
public :: ps_ValenceCharge_GetRawData
!
public :: ps_CoreCharge_Value
public :: ps_CoreCharge_GetRawData
public :: ps_CoreCharge_MatchingRadius
public :: ps_CoreCharge_NumberOfKeptDerivatives
!

private

CONTAINS !===============================================

!>  Returns the library version in integer format
function ps_GetLibPSMLVersion() result(v)
integer :: v
  v = PSML_LIBRARY_VERSION
end function ps_GetLibPSMLVersion

!>  Sets the global debug flag
!> @param debug:logical
subroutine ps_SetDebug(debug)
logical, intent(in) :: debug
  global_debug = debug
end subroutine ps_SetDebug

#ifndef __NO_PROC_POINTERS__

!> Sets the default interpolator and
!> its quality parameter
subroutine ps_SetInterpolator(func,nquality)
use m_interp, only: interpolator, nq

! Parameter for interpolator's quality
! It might mean different things for different
! interpolators
integer, intent(in) :: nquality

interface
   subroutine func(nquality,x,y,npts,r,val,debug)

     integer, parameter :: dp = selected_real_kind(10,100)

     integer, intent(in)  :: nquality  ! Quality parameter
     real(dp), intent(in) :: x(*), y(*)
     integer, intent(in)  :: npts    ! Size of x, y arrays
     real(dp), intent(in) :: r
     real(dp), intent(out):: val
     logical, intent(in) :: debug
   end subroutine func
end interface

  interpolator => func
  nq = nquality

end subroutine ps_SetInterpolator

#endif

!> Sets the quality parameter of the current
!> default interpolator. Useful when we do
!> not care about the type of evaluator, but
!> want to compare different qualities
subroutine ps_SetInterpolatorQuality(nquality)
use m_interp, only: nq

! Parameter for interpolator's quality
! It might mean different things for different
! interpolators
integer, intent(in) :: nquality

  nq = nquality
end subroutine ps_SetInterpolatorQuality

!> Returns the number of non-empty valence shells
!> in the ps generation configuration
!> @param ps is a handle to the psml information
function ps_NValenceShells(ps) result(nshells)
type(ps_t), intent(in) :: ps
integer                :: nshells

nshells = ps%config_val%nshells

end function ps_NValenceShells

!> Returns the angular momentum of the i'th valence shell
!> in the ps generation configuration
!> @param ps is a handle to the psml information
!> @param i is the index of the shell
!> @note i should be within range
function ps_ValenceShellL(ps,i) result(l)
type(ps_t), intent(in) :: ps
integer, intent(in)    :: i
integer                :: l

character(len=1) :: str

call check_index(i,ps%config_val%nshells,"valence shell")
l = l_of_sym(ps%config_val%l(i),"valence shell")

end function ps_ValenceShellL

!>  Returns the principal quantum number of the i'th valence shell
!> in the ps generation configuration
!> @author Alberto Garcia
!> @date 2014
!> @param ps is a handle to the psml information
!> @param i is the index of the shell
!> @note i should be within range
function ps_ValenceShellN(ps,i) result(n)
type(ps_t), intent(in) :: ps
integer, intent(in)    :: i
integer                :: n

call check_index(i,ps%config_val%nshells,"valence shell")
n = ps%config_val%n(i)

end function ps_ValenceShellN

!>  Returns the occupation of the i'th valence shell
!> in the ps generation configuration
!> @author Alberto Garcia
!> @date 2014
!> @param ps is a handle to the psml information
!> @param i is the index of the shell
!> @param channel is an optional parameter for spin-polarized 
!> calculations ("u" or "d"). 
!> @note i should be within range
!> @note If "channel" is present, the occupation returned
!> corresponds to the given channel
function ps_ValenceShellOccupation(ps,i,channel) result(occ)
type(ps_t), intent(in) :: ps
integer, intent(in)    :: i
character(len=1), intent(in), optional :: channel
real(dp)                :: occ

call check_index(i,ps%config_val%nshells,"valence shell")

if (present(channel)) then
   if (ps_IsSpinPolarized(ps)) then
      if (channel == "u") then
         occ = ps%config_val%occ_up(i)
      else if (channel == "d") then
         occ = ps%config_val%occ_down(i)
      else
         call die("Wrong channel in ValShellOccupation")
      endif
   else
      call die("Cannot speficy channel in ValShellOccupation")
   endif
else
   occ = ps%config_val%occ(i)
endif

end function ps_ValenceShellOccupation
!
!-------------------------------------------------------
!>  Computes the value of the valence charge at r
!> @param ps is a handle to the psml information
!> @param r is the radius
!> It returns the valence charge density integrated over
!> solid angle, so that Q_val = int{ val*r*r }
!> 
function ps_ValenceCharge_Value(ps,r) result(val)
type(ps_t), intent(in) :: ps
real(dp), intent(in)       :: r
real(dp)                   :: val

if (r > max_range(ps%valence_charge%rho_val)) then
   val = 0.0_dp
else
   val = eval_radfunc(ps%valence_charge%rho_val,r,debug=global_debug)
endif
end function ps_ValenceCharge_Value

subroutine ps_ValenceCharge_GetRawData(ps,raw_r,raw_data)
type(ps_t), intent(in) :: ps
real(dp), allocatable, intent(out)  :: raw_r(:), raw_data(:)

call get_raw_radfunc(ps%valence_charge%rho_val,raw_r,raw_data)

end subroutine ps_ValenceCharge_GetRawData

!>  Computes the value of the pseudo-core charge at r
!> @param ps is a handle to the psml information
!> @param r is the radius
!> It returns the pseudo-core charge density integrated over
!> solid angle, so that Q_core = int{ val*r*r }
!> 
function ps_CoreCharge_Value(ps,r) result(val)
type(ps_t), intent(in) :: ps
real(dp), intent(in)       :: r
real(dp)                   :: val

if (r > max_range(ps%core_charge%rho_core)) then
   val = 0.0_dp
else
   val = eval_radfunc(ps%core_charge%rho_core,r,debug=global_debug)
endif
end function ps_CoreCharge_Value

subroutine ps_CoreCharge_GetRawData(ps,raw_r,raw_data)
type(ps_t), intent(in) :: ps
real(dp), allocatable, intent(out)  :: raw_r(:), raw_data(:)

call get_raw_radfunc(ps%core_charge%rho_core,raw_r,raw_data)

end subroutine ps_CoreCharge_GetRawData


!>  Returns the matching radius used in the fit of the
!>  pseudo-core charge
!> @param ps is a handle to the psml information
!> It returns -1.0 if the information is not provided in the file
function ps_CoreCharge_MatchingRadius(ps) result(rmatch)
type(ps_t), intent(in) :: ps
real(dp)               :: rmatch

rmatch = ps%core_charge%rcore

end function ps_CoreCharge_MatchingRadius

!>  Returns the number of derivatives kept in the fit of the
!>  pseudo-core charge
!> @param ps is a handle to the psml information
!> It returns -1 if the information is not provided in the file
function ps_CoreCharge_NumberOfKeptDerivatives(ps) result(n)
type(ps_t), intent(in) :: ps
integer                :: n

n = ps%core_charge%n_cont_derivs

end function ps_CoreCharge_NumberOfKeptDerivatives

!>  Returns the atomic symbol
!> @param ps is a handle to the psml information
!
function ps_AtomicSymbol(ps) result(name)
type(ps_t), intent(in) :: ps
character(len=2) :: name
name = ps%header%atomic_label(1:2)
end function ps_AtomicSymbol
!
!>  Returns the atomic label
!> @param ps is a handle to the psml information
!> The label is more general than the atomic symbol
!>
function ps_AtomicLabel(ps) result(name)
type(ps_t), intent(in) :: ps
character(len=len_trim(ps%header%atomic_label)) :: name
name = trim(ps%header%atomic_label)
end function ps_AtomicLabel
!
!>  Returns the atomic number
!> @param ps is a handle to the psml information
!>
function ps_AtomicNumber(ps) result(z)
type(ps_t), intent(in) :: ps
real(dp) :: z
 z = ps%header%z
end function ps_AtomicNumber
!
!>  Returns the PSML version used in the file
!
function ps_GetPSMLVersion(ps) result(version)
type(ps_t), intent(in) :: ps
character(len=10) :: version
version = ps%version
end function ps_GetPSMLVersion
!
!>  Returns the unique uuid of the PSML file
!
function ps_GetUUID(ps) result(name)
type(ps_t), intent(in) :: ps
character(len=36) :: name
name = ps%uuid
end function ps_GetUUID
!
!>  Returns the creator of the PSML file, as
!> written in the provenance element.
!> @param ps is a handle to the psml information
!> The format is arbitrary
!
function ps_Creator(ps) result(name)
type(ps_t), intent(in) :: ps
character(len=len_trim(ps%provenance%creator)) :: name
name = trim(ps%provenance%creator)
end function ps_Creator
!
!>  Returns the date appearing in the provenance element.
!> @param ps is a handle to the psml information
!> The format is arbitrary
!
function ps_Date(ps) result(str)
type(ps_t), intent(in) :: ps
character(len=len_trim(ps%provenance%date)) :: str
str = trim(ps%provenance%date)
end function ps_Date
!
! **AG**
!>
!
function ps_PseudoFlavor(ps) result(str)
!+ category: needs_work
!+  Returns the pseudization flavor appearing in the header element.
  !*### NOTE
  ! Generalize to accept an index argument for potential
  ! In principle, different flavors for different channels are
  ! possible in the PSML file, but there is not support in the API for it.
type(ps_t), intent(in) :: ps
!+ Handle for PSML information
character(len=len_trim(ps%header%flavor)) :: str
str = trim(ps%header%flavor)
end function ps_PseudoFlavor

function ps_ZPseudo(ps) result(zpseudo)
!*  Returns the effective valence of the pseudo-atom,
! i.e., the atomic number minus the number of "core" electrons.
  type(ps_t), intent(in) :: ps
!+ Handle for PSML information
real(dp)                   :: zpseudo
zpseudo = ps%header%zpseudo
end function ps_ZPseudo

!>  Returns the total valence charge density in the
!> atomic configuration used to generate the pseudopotential.
!> !> @param ps is a handle to the psml information
!>
function ps_GenerationZval(ps) result(zval)
type(ps_t), intent(in) :: ps
real(dp)                   :: zval
zval = ps%config_val%total_charge
end function ps_GenerationZval
!
!>  Returns the number of libxc functionals that
!> would correspond to the exchange-correlation scheme
!> used in the generation code.
!> @param ps is a handle to the psml information
!>
function ps_NLibxcFunctionals(ps) result(xc_n)
type(ps_t), intent(in) :: ps
integer                :: xc_n

xc_n =  ps%xc_info%n_functs_libxc 
end function ps_NLibxcFunctionals
!
function ps_LibxcName(ps,i) result(xc_name)
type(ps_t), intent(in) :: ps
integer, intent(in)    :: i
character(len=50)      :: xc_name

call check_index(i,ps%xc_info%n_functs_libxc,"libxc functional")
xc_name = ps%xc_info%libxc_name(i)
end function ps_LibxcName
!
function ps_LibxcId(ps,i) result(xc_id)
type(ps_t), intent(in) :: ps
integer, intent(in)    :: i
integer                :: xc_id

call check_index(i,ps%xc_info%n_functs_libxc,"libxc functional")
xc_id = ps%xc_info%libxc_id(i)
end function ps_LibxcId
!
function ps_LibxcWeight(ps,i) result(xc_weight)
type(ps_t), intent(in) :: ps
integer, intent(in)    :: i
real(dp)               :: xc_weight

call check_index(i,ps%xc_info%n_functs_libxc,"libxc functional")
xc_weight = ps%xc_info%libxc_weight(i)
end function ps_LibxcWeight
!
function ps_LibxcIdArray(ps) result(xc_id_array)
type(ps_t), intent(in) :: ps
integer                :: xc_id_array(2)

xc_id_array(:) = ps%xc_info%libxc_id(:)
end function ps_LibxcIdArray
!
function ps_ValidLibxc(ps) result(libxc_ok)
type(ps_t), intent(in), target :: ps
logical                :: libxc_ok

integer, pointer  :: xc_id_array(:)

xc_id_array => ps%xc_info%libxc_id(:)
libxc_ok = .true.
if (any (xc_id_array(:) <= 0)) then
   libxc_ok = .false.
endif
end function ps_ValidLibxc
!
!=================================================
!>  Returns the annotation associated to a
!> given element. For grids, if a radial function
!> handle is given, the annotation for that 
!> radial function's grid is returned. Otherwise,
!> the return value is the annotation for the global grid.
!> If there is no appropriate annotation, an empty
!> structure is returned.
!> @note Support for radfunc handles is not yet implemented
!> @param ps is a handle to the psml information
!> @param name is the element name or a common alias
!> @param radfunc is a handle to a radial function structure
!>
function ps_GetAnnotation(ps,name) result(annotation)
                             !! ,radfunc)  To be implemented
 type(ps_t), intent(in)                 :: ps
 character(len=*), intent(in)           :: name
!!  type(radfunc_t), intent(in), optional  :: radfunc
 
 type(ps_annotation_t)  :: annotation
 type(ps_annotation_t), pointer  :: annotation_p

select case (name)

   case ("psml","PSML","top-level","global")
      annotation = ps%annotation
   case ("provenance")
         annotation = ps%provenance%annotation
   case ("exchange-correlation","xc","XC")
         annotation = ps%xc_info%annotation
   case ("valence-configuration")
         annotation = ps%config_val%annotation

   case ("grid")

!!$      if (present(radfunc)) then
!!$         ! We are told to get the grid annotation
!!$         ! for a specific radial function
!!$
!!$         if (.not. associated(radfunc%grid)) then
!!$            call die("get_annotation: Invalid radial function")
!!$         endif
!!$         annotation = radfunc%grid%annotation
!!$
!!$      else

         ! This is the global grid annotation
         if (.not. initialized(ps%global_grid)) then
            annotation = EMPTY_ANNOTATION
         else
            annotation_p => annotationGrid(ps%global_grid)
            annotation = annotation_p
         endif

!!$      endif

   case ("semilocal-potentials")
         annotation = ps%semilocal%annotation
   case ("nonlocal-projectors")
         annotation = ps%nonlocal%annotation
   case ("local-potential")
         annotation = ps%local%annotation
   case ("pseudo-wavefunctions")
         annotation = ps%pswfs%annotation
   case ("valence-charge")
         annotation = ps%valence_charge%annotation
   case ("core-charge")
         annotation = ps%core_charge%annotation
   case default
      call die("Unrecognized annotation name: "//trim(name))

   end select

 end function ps_GetAnnotation
!
! The following two functions are retained for
! backwards compatibility
!
function ps_XCAnnotation(ps) result(xc_annotation)
type(ps_t), intent(in) :: ps
type(ps_annotation_t)  :: xc_annotation
xc_annotation = ps%xc_info%annotation
end function ps_XCAnnotation
!
!=================================================
function ps_Relativity(ps) result(rel)
type(ps_t), intent(in) :: ps
character(len=6)       :: rel
rel = ps%header%relativity
end function ps_Relativity
!
function ps_IsSpinPolarized(ps) result(pol)
type(ps_t), intent(in) :: ps
logical                    :: pol
pol = ps%header%polarized
end function ps_IsSpinPolarized
!
function ps_HasCoreCorrections(ps) result(cc)
type(ps_t), intent(in) :: ps
logical                    :: cc
cc = (ps%header%core_corrections == "yes")
end function ps_HasCoreCorrections
!=================================================
!
!----- Convenience set handling routines
!
subroutine ps_Get_Potential_Indexes(ps,set,indexes)
type(ps_t), intent(in)                 :: ps
integer, intent(in)                    :: set
integer, allocatable, intent(inout)    :: indexes(:)

integer :: i, n

n = 0
do i = 1, size(ps%sl_table)
   if (iand(ps%sl_table(i)%p%set,set) /= 0) n = n+1
enddo

if (allocated(indexes))  deallocate(indexes)
allocate(indexes(n))

n = 0
do i = 1, size(ps%sl_table)
   if (iand(ps%sl_table(i)%p%set,set) /= 0) then
      n = n+1
      indexes(n) = i
   end if
enddo

end subroutine ps_Get_Potential_Indexes

subroutine ps_Get_Projector_Indexes(ps,set,indexes)
type(ps_t), intent(in)                 :: ps
integer, intent(in)                    :: set
integer, allocatable, intent(inout)    :: indexes(:)

integer :: i, n

n = 0
do i = 1, size(ps%nl_table)
   if (iand(ps%nl_table(i)%p%set,set) /= 0) n = n+1
enddo

if (allocated(indexes))  deallocate(indexes)
allocate(indexes(n))

n = 0
do i = 1, size(ps%nl_table)
   if (iand(ps%nl_table(i)%p%set,set) /= 0) then
      n = n+1
      indexes(n) = i
   end if
enddo

end subroutine ps_Get_Projector_Indexes

subroutine ps_Get_PseudoWf_Indexes(ps,set,indexes)
type(ps_t), intent(in)                 :: ps
integer, intent(in)                    :: set
integer, allocatable, intent(inout)    :: indexes(:)

call GetSetIndexes(ps%pswfs%npswfs,ps%pswfs%set,set,indexes)
end subroutine ps_Get_PseudoWf_Indexes

function ps_Number_Of_Potentials(ps,set) result(n)
type(ps_t), intent(in)                 :: ps
integer, intent(in)                    :: set
integer                         :: n

integer, allocatable  :: idx(:)
call ps_Get_Potential_Indexes(ps,set,idx)
n = size(idx)
deallocate(idx)
end function ps_Number_Of_Potentials

function ps_Number_Of_Projectors(ps,set) result(n)
type(ps_t), intent(in)                 :: ps
integer, intent(in)                    :: set
integer                         :: n

integer, allocatable  :: idx(:)
call ps_Get_Projector_Indexes(ps,set,idx)
n = size(idx)
deallocate(idx)
end function ps_Number_Of_Projectors

function ps_Number_Of_PseudoWfs(ps,set) result(n)
type(ps_t), intent(in)                 :: ps
integer, intent(in)                    :: set
integer                         :: n

integer, allocatable  :: idx(:)
call ps_Get_PseudoWf_Indexes(ps,set,idx)
n = size(idx)
deallocate(idx)
end function ps_Number_Of_PseudoWfs
!
!====================================================
! Semilocal potentials
!
function ps_Potential_Set(ps,i) result(set)
type(ps_t), intent(in) :: ps
integer,   intent(in)      :: i
integer                    :: set

call check_index(i,size(ps%sl_table),"SL pot")
set = ps%sl_table(i)%p%set
!
end function ps_Potential_Set
!
function ps_Potential_L(ps,i) result(l)
type(ps_t), intent(in) :: ps
integer,   intent(in)      :: i
integer                    :: l

call check_index(i,size(ps%sl_table),"SL pot")
l = l_of_sym(ps%sl_table(i)%p%l,"SL pot")

end function ps_Potential_L
!
function ps_Potential_J(ps,i) result(j)
type(ps_t), intent(in) :: ps
integer,   intent(in)      :: i
real(dp)                   :: j

call check_index(i,size(ps%sl_table),"SL pot")
if (  iand(ps%sl_table(i)%p%set, &
           SET_LJ+SET_USER1+SET_USER2)  &
       == 0) then
   call die("Attempt to get j from wrong set")
endif
j = ps%sl_table(i)%p%j

end function ps_Potential_J
!
function ps_Potential_Rc(ps,i) result(rc)
type(ps_t), intent(in) :: ps
integer,   intent(in)      :: i
real(dp)                   :: rc

call check_index(i,size(ps%sl_table),"SL pot")
rc = ps%sl_table(i)%p%rc

end function ps_Potential_Rc
!
function ps_Potential_N(ps,i) result(n)
type(ps_t), intent(in) :: ps
integer,   intent(in)      :: i
integer                    :: n

call check_index(i,size(ps%sl_table),"SL pot")
n = ps%sl_table(i)%p%n

end function ps_Potential_N
!
!> Evaluator by storage index
function ps_Potential_Value(ps,i,r) result(val)
type(ps_t), intent(in) :: ps
integer,   intent(in)      :: i
real(dp),  intent(in)      :: r
real(dp)                   :: val

logical           :: coulomb_tail

call check_index(i,size(ps%sl_table),"SL pot")

select case (ps%sl_table(i)%p%set)
   case (SET_SO, SET_SPINDIFF)
      coulomb_tail = .false.
   case default
      coulomb_tail = .true.
end select
   
if (r> max_range(ps%sl_table(i)%p%V)) then
   if (coulomb_tail) then
      val = - ps_ZPseudo(ps)/r
   else
      val = 0.0_dp
   endif
else
   val = eval_radfunc(ps%sl_table(i)%p%V,r, debug=global_debug)
endif

end function ps_Potential_Value

subroutine ps_Potential_GetRawData(ps,i,raw_r,raw_data)
type(ps_t), intent(in) :: ps
integer, intent(in)    :: i
real(dp), allocatable, intent(out)  :: raw_r(:), raw_data(:)

call get_raw_radfunc(ps%sl_table(i)%p%V,raw_r,raw_data)

end subroutine ps_Potential_GetRawData

!
!====================================================
! PseudoWavefunctions
!
! Basic accessors
!
function ps_PseudoWf_L(ps,i) result(l)
type(ps_t), intent(in) :: ps
integer,   intent(in)      :: i
integer                    :: l

call check_index(i,ps%pswfs%npswfs,"pswf")
l = l_of_sym(ps%pswfs%l(i),"pswf")

end function ps_PseudoWf_L
!
function ps_PseudoWf_J(ps,i) result(j)
type(ps_t), intent(in) :: ps
integer,   intent(in)      :: i
real(dp)                   :: j

call check_index(i,ps%pswfs%npswfs,"pswf")
if (  iand(ps%pswfs%set(i), &
           SET_LJ+SET_USER1+SET_USER2)  &
       == 0) then
   call die("Attempt to get j from wrong set")
endif
j = ps%pswfs%j(i)

end function ps_PseudoWf_J
!
function ps_PseudoWf_N(ps,i) result(n)
type(ps_t), intent(in) :: ps
integer,   intent(in)      :: i
integer                    :: n

call check_index(i,ps%pswfs%npswfs,"pswf")
n = ps%pswfs%n(i)

end function ps_PseudoWf_N
!
function ps_PseudoWf_Set(ps,i) result(set)
type(ps_t), intent(in) :: ps
integer,   intent(in)      :: i
integer              :: set

call check_index(i,ps%pswfs%npswfs,"pswf")
set = ps%pswfs%set(i)
!
end function ps_PseudoWf_Set
!
!
function ps_PseudoWf_Value(ps,i,r) result(val)
type(ps_t), intent(in) :: ps
integer,   intent(in)  :: i
real(dp),  intent(in)  :: r
real(dp)               :: val

call check_index(i,ps%pswfs%npswfs,"pswf")
if (r> max_range(ps%pswfs%Phi(i))) then
   val = 0.0_dp
else
   val = eval_radfunc(ps%pswfs%Phi(i),r, debug=global_debug)
endif

end function ps_PseudoWf_Value
!
!
subroutine ps_PseudoWf_GetRawData(ps,i,raw_r,raw_data)
type(ps_t), intent(in) :: ps
integer, intent(in)    :: i
real(dp), allocatable, intent(out)  :: raw_r(:), raw_data(:)

call get_raw_radfunc(ps%pswfs%Phi(i),raw_r,raw_data)

end subroutine ps_PseudoWf_GetRawData
!
!====================================================
! Pseudopotential operator (Vlocal, LocalCharge, projectors)
!
! Basic accessors
!
function ps_HasProjectors(ps) result(p)
type(ps_t), intent(in) :: ps
logical                    :: p
!
p = (size(ps%nl_table) > 0)

end function ps_HasProjectors

function ps_HasLocalPotential(ps) result(p)
type(ps_t), intent(in) :: ps
logical                    :: p
!
p = (initialized(ps%local%Vlocal%grid))

end function ps_HasLocalPotential

function ps_HasLocalCharge(ps) result(p)
type(ps_t), intent(in) :: ps
logical                    :: p
!
p = (initialized(ps%local%chlocal%grid))

end function ps_HasLocalCharge

function ps_HasPSOperator(ps) result(psop)
type(ps_t), intent(in) :: ps
logical                    :: psop
!
psop = (ps_HasProjectors(ps) .and. ps_HasLocalPotential(ps))

end function ps_HasPSOperator
!
function ps_HasSemilocalPotentials(ps) result(p)
type(ps_t), intent(in) :: ps
logical                    :: p
!
p = (associated(ps%semilocal))

end function ps_HasSemilocalPotentials

!
function ps_LocalPotential_Value(ps,r) result(val)
type(ps_t), intent(in) :: ps
real(dp), intent(in)       :: r
real(dp)                   :: val

if (r > max_range(ps%local%vlocal)) then

   ! There should be a sanity check upon parsing
   ! to guarantee that the coulomb behavior has
   ! been reached at the end of the grid range

   val = - ps_ZPseudo(ps)/r
else
   val = eval_radfunc(ps%local%vlocal,r,debug=global_debug)
endif
end function ps_LocalPotential_Value

function ps_LocalCharge_Value(ps,r) result(val)
type(ps_t), intent(in) :: ps
real(dp), intent(in)       :: r
real(dp)                   :: val

if (r > max_range(ps%local%chlocal)) then

   ! There should be a sanity check upon parsing
   ! to guarantee that the coulomb behavior has
   ! been reached at the end of the grid range

   val = 0.0_dp
else
   val = eval_radfunc(ps%local%chlocal,r,debug=global_debug)
endif
end function ps_LocalCharge_Value
!
function ps_LocalPotential_Type(ps) result(type)
type(ps_t), intent(in) :: ps
character(len=40)          :: type

type = ps%local%vlocal_type
!
end function ps_LocalPotential_Type
!
subroutine ps_LocalPotential_GetRawData(ps,raw_r,raw_data)
type(ps_t), intent(in) :: ps
real(dp), allocatable, intent(out)  :: raw_r(:), raw_data(:)

call get_raw_radfunc(ps%local%vlocal,raw_r,raw_data)

end subroutine ps_LocalPotential_GetRawData

function ps_Projector_L(ps,i) result(l)
type(ps_t), intent(in) :: ps
integer,   intent(in)      :: i
integer                    :: l

call check_index(i,size(ps%nl_table),"proj")
l = l_of_sym(ps%nl_table(i)%p%l,"proj")

end function ps_Projector_L
!
function ps_Projector_J(ps,i) result(j)
type(ps_t), intent(in) :: ps
integer,   intent(in)      :: i
real(dp)                   :: j

call check_index(i,size(ps%nl_table),"proj")
if (  iand(ps%nl_table(i)%p%set, &
           SET_LJ+SET_USER1+SET_USER2)  &
       == 0) then
   call die("Attempt to get j from wrong set")
endif
j = ps%nl_table(i)%p%j

end function ps_Projector_J
!
function ps_Projector_Ekb(ps,i) result(ekb)
type(ps_t), intent(in) :: ps
integer,   intent(in)  :: i
real(dp)                   :: ekb

call check_index(i,size(ps%nl_table),"proj")
ekb = ps%nl_table(i)%p%ekb

end function ps_Projector_Ekb
!
function ps_Projector_Seq(ps,i) result(seq)
type(ps_t), intent(in) :: ps
integer,   intent(in)      :: i
integer                    :: seq

call check_index(i,size(ps%nl_table),"proj")
seq = ps%nl_table(i)%p%seq

end function ps_Projector_Seq
!
function ps_Projector_Set(ps,i) result(set)
type(ps_t), intent(in) :: ps
integer,   intent(in)      :: i
integer                 :: set

call check_index(i,size(ps%nl_table),"proj")
set = ps%nl_table(i)%p%set
!
end function ps_Projector_Set
!
function ps_Projector_Type(ps,i) result(type)
type(ps_t), intent(in) :: ps
integer,   intent(in)      :: i
character(len=40)          :: type

call check_index(i,size(ps%nl_table),"proj")
type = ps%nl_table(i)%p%type
!
end function ps_Projector_Type

function ps_Projector_Value(ps,i,r) result(val)
!+ display: private
type(ps_t), intent(in) :: ps
integer,   intent(in)  :: i
real(dp),  intent(in)  :: r
real(dp)               :: val

call check_index(i,size(ps%nl_table),"proj")
if (r> max_range(ps%nl_table(i)%p%proj)) then
   val = 0.0_dp
else
   val = eval_radfunc(ps%nl_table(i)%p%proj,r, debug=global_debug)
endif

end function ps_Projector_Value

subroutine ps_Projector_GetRawData(ps,i,raw_r,raw_data)
!+ deprecated: true
type(ps_t), intent(in) :: ps
integer, intent(in)    :: i
real(dp), allocatable, intent(out)  :: raw_r(:), raw_data(:)

call get_raw_radfunc(ps%nl_table(i)%p%proj,raw_r,raw_data)

end subroutine ps_Projector_GetRawData
!
subroutine ps_Get_Projector_Indexes_byL(ps,l,idxset,idxl)
!* Subset of projectors with given l
! Note that this function takes an array of indexes
! and returns another array of indexes
! There is currently no way to check that the idxset
! really corresponds to projectors...
type(ps_t), intent(in) :: ps
integer, intent(in)    :: l
integer, intent(in)    :: idxset(:)
integer, allocatable, intent(inout)   :: idxl(:)

integer                    :: n

integer :: n_in_set, i, l_i

n_in_set = size(idxset)
n = 0
do i = 1, n_in_set
  l_i = ps_Projector_L(ps,idxset(i))
  if (l_i == l) n = n + 1
enddo
if (allocated(idxl)) deallocate(idxl)
allocate(idxl(n))
n = 0
do i = 1, n_in_set
  l_i = ps_Projector_L(ps,idxset(i))
  if (l_i == l) then
     n = n + 1
     idxl(n) = idxset(i)
  endif
enddo

end subroutine ps_Get_Projector_Indexes_byL
!
!====================================================
! Low-level routines
!
!
subroutine GetSetIndexes(nitems,setarray,set,indexes)
integer, intent(in)                    :: nitems
integer, intent(in)                    :: setarray(:)
integer, intent(in)                    :: set
integer, allocatable, intent(inout)    :: indexes(:)

integer ::  n, i

n = 0
do i = 1, nitems
   if (iand(setarray(i),set) /= 0) n = n+1
enddo
if (allocated(indexes))  deallocate(indexes)
allocate(indexes(n))
n = 0
do i = 1, nitems
   if (iand(setarray(i),set) /= 0) then
      n = n+1
      indexes(n) = i
   end if
enddo
end subroutine GetSetIndexes
!
subroutine check_index(i,n,str)
integer, intent(in) :: i, n
character(len=*), intent(in) :: str

call assert( (i <=  n), "Index overflow in "//trim(str))
call assert( (i >  0), "Non-positive index in "//trim(str))
end subroutine check_index
!
function l_of_sym(str,name) result(l)
character(len=*), intent(in) :: str, name
integer                :: l
!
! This routine will disappear once we store
! l as integer in the data structure
!
do l = 0,4
   if (str == sym(l)) RETURN
enddo
call die("Wrong l symbol in "//trim(name))
end function l_of_sym

!>  Returns the maximum radius in a radfunc's grid
function max_range(f) result(range)
type(radfunc_t), intent(in) :: f
real(dp)                  :: range

real(dp),  pointer :: a(:)

integer :: npts

npts = sizeGrid(f%grid)
a => valGrid(f%grid)
range = a(npts)

end function max_range
!----------
function eval_radfunc(f,r,debug) result(val)
use m_interp, only: interpolator, nq

type(radfunc_t), intent(in) :: f
real(dp), intent(in)      :: r
real(dp)                  :: val
logical, intent(in)       :: debug

real(dp), pointer :: x(:) => null(), y(:) => null()

x => valGrid(f%grid)
y => f%data(:)

call interpolator(nq,x,y,size(x),r,val,debug)

end function eval_radfunc

subroutine get_raw_radfunc(f,raw_r,raw_data)
!
type(radfunc_t), intent(in) :: f
real(dp), intent(out), allocatable  :: raw_r(:)
real(dp), intent(out), allocatable  :: raw_data(:)

integer npts
real(dp), pointer :: a(:) 

npts = sizeGrid(f%grid)
allocate(raw_r(npts), raw_data(npts))

a => valGrid(f%grid)
raw_r(:) = a(:)
raw_data(:) = f%data(:)

end subroutine get_raw_radfunc
!
      FUNCTION atomic_number(SYMBOL) result(z)

! Given the atomic symbol, it returns the atomic number
! Based on code by J. Soler

      character(len=2), intent(in)    :: SYMBOL  ! Atomic symbol
      integer                         :: Z       ! Atomic number

      integer, parameter  :: NZ=103
      character(len=2), parameter :: NAME(NZ) =  &
               (/'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne', &
                 'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca', &
                 'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn', &
                 'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr', &
                 'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn', &
                 'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd', &
                 'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', &
                 'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg', &
                 'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', &
                 'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm', &
                 'Md','No','Lr'/)

     do z = 1, NZ
        if (SYMBOL == NAME(Z)) then
           RETURN
        endif
     enddo
     call die("Cannot find atomic number for " // symbol)
        
   end FUNCTION atomic_number

   subroutine assert(cond,message)
     logical, intent(in) :: cond
     character(len=*) message

     if (.not. cond) call die(message)
   end subroutine assert

 end module m_psml_api






