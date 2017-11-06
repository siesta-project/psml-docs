!+ libPSML API implementation
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module m_psml_api
!+ graph: false
!+ Procedures to handle the PSML pseudopotential format.

use m_psml_core
use m_aux_aliases, only: ps_radfunc_t

use assoc_list, only: ps_annotation_t => assoc_list_t
use assoc_list, only: EMPTY_ANNOTATION => EMPTY_ASSOC_LIST

use external_interfaces, only: die => psml_die
use class_Grid

implicit none

integer, parameter    :: dp = selected_real_kind(14)
logical               :: global_debug = .false.
logical               :: global_use_effective_range = .true.
character(len=1), dimension(0:4) :: sym = (/ "s", "p", "d", "f", "g" /)

! Library operation routines
public :: ps_GetLibPSMLVersion
public :: ps_SetEvaluatorOptions
!
! Accessor list
!
public :: ps_RootAttributes_Get
public :: ps_Provenance_Depth
public :: ps_Provenance_Get

public :: ps_PseudoAtomSpec_Get
public :: ps_ExchangeCorrelation_Get
public :: ps_LibxcFunctional_Get

public :: ps_ValenceConfiguration_Get
public :: ps_ValenceShell_Get
!
public :: ps_ValenceCharge_Value
public :: ps_ValenceCharge_Get
!
public :: ps_CoreCharge_Value
public :: ps_CoreCharge_Get
!
! Semilocal potentials
!
public :: ps_HasSemilocalPotentials
public :: ps_SemilocalPotentials_Filter
public :: ps_Potential_Get
public :: ps_Potential_Value

! Pseudo Operator
!
public :: ps_HasPsOperator
!
! Vlocal 
!
public :: ps_HasLocalPotential
public :: ps_LocalPotential_Get
public :: ps_LocalPotential_Value
!
! Projectors
!
public :: ps_HasProjectors
public :: ps_NonlocalProjectors_Filter
public :: ps_Projector_Get
public :: ps_Projector_Value
!
! Pseudo-wave-functions
!
public :: ps_PseudoWavefunctions_Filter
public :: ps_PseudoWf_Get
public :: ps_PseudoWf_Value

! Grid Annotations
public :: ps_GridAnnotation 
!
! Exported low-level routines
!
public :: ps_GetValue
public :: ps_GetRawData
!
! Aliases
!
interface ps_Potential_Filter
   module procedure ps_SemiLocalPotentials_Filter
end interface ps_Potential_Filter
public :: ps_Potential_Filter

interface ps_Projector_Filter
   module procedure ps_NonLocalProjectors_Filter
end interface ps_Projector_Filter
public :: ps_Projector_Filter

interface ps_PseudoWf_Filter
   module procedure ps_PseudoWaveFunctions_Filter
end interface ps_PseudoWf_Filter
public :: ps_PseudoWf_Filter

private

CONTAINS 

!
! ==============================================================
!
!>  Returns the library version in integer format
function ps_GetLibPSMLVersion() result(v)
integer :: v
  v = PSML_LIBRARY_VERSION
end function ps_GetLibPSMLVersion

!> Sets various parameters for the operation of
!> the evaluator
subroutine ps_SetEvaluatorOptions(quality_level,debug,&
     use_effective_range,&
     custom_interpolator)

use m_psml_interp, only: nq, interpolator

! Parameter for interpolator's quality
! It might mean different things for different
! interpolators
integer, intent(in), optional :: quality_level

logical, intent(in), optional :: debug
logical, intent(in), optional :: use_effective_range

interface
   subroutine interp(nquality,x,y,npts,r,val,debug)

     integer, parameter :: dp = selected_real_kind(10,100)

     integer, intent(in)  :: nquality  ! Quality parameter
     real(dp), intent(in) :: x(*), y(*)
     integer, intent(in)  :: npts    ! Size of x, y arrays
     real(dp), intent(in) :: r
     real(dp), intent(out):: val
     logical, intent(in) :: debug
   end subroutine interp
end interface

procedure(interp), optional :: custom_interpolator

if (present(quality_level)) then
   nq = quality_level
endif
if (present(debug)) then
   global_debug = debug
endif
if (present(use_effective_range)) then
   global_use_effective_range = use_effective_range
endif
if (present(custom_interpolator)) then
   interpolator => custom_interpolator
endif

end subroutine ps_SetEvaluatorOptions
!
! ==============================================================
!
subroutine ps_RootAttributes_Get(ps,uuid,version,namespace,annotation)
  type(ps_t), intent(in) :: ps
  
  character(len=*), intent(out), optional :: uuid
  character(len=*), intent(out), optional :: version
  character(len=*), intent(out), optional :: namespace
  type(ps_annotation_t), intent(out), optional :: annotation

  if (present(uuid)) then
     uuid = trim(ps%uuid)
  endif
  if (present(version)) then
     version = trim(ps%version)
  endif
  if (present(annotation)) then
     annotation = ps%annotation
  endif
  if (present(namespace)) then
     namespace = ps%namespace
  endif

end subroutine ps_RootAttributes_Get
!
! ==============================================================
!
function ps_Provenance_Depth(ps) result(n)
  type(ps_t), intent(in) :: ps
  integer :: n
  
  type(provenance_t), pointer  :: p

  n = 0
  p => ps%provenance
  do while (associated(p))
     n = n + 1
     p => p%next
  enddo
end function ps_Provenance_Depth
  
subroutine ps_Provenance_Get(ps,level,creator,date,&
     annotation,number_of_input_files)
  type(ps_t), intent(in) :: ps
  integer   , intent(in) :: level

  character(len=*), intent(out), optional :: creator
  character(len=*), intent(out), optional :: date
  type(ps_annotation_t), intent(out), optional :: annotation
  integer, intent(out), optional :: number_of_input_files

  type(provenance_t), pointer  :: p
  logical :: found_level

  ! Here "level" means "record_number", with
  ! the oldest having a value of 1.
  
  found_level = .false.
  p => ps%provenance
  do while (associated(p))
     if (p%record_number == level) then
        found_level = .true.
        exit
     endif
     p => p%next
  enddo
  
  if (.not. found_level) call die("Cannot reach provenance level")

  if (present(creator)) then
     creator = p%creator
  endif
  if (present(date)) then
     date = p%date
  endif
  if (present(number_of_input_files)) then
     number_of_input_files = p%n_input_files
  endif
  if (present(annotation)) then
     annotation = p%annotation
  endif
end subroutine ps_Provenance_Get

! To be implemented
!subroutine ps_Provenance_InputFile(ps,level,file_index,&
!                             filename,file_content)
!
! ===================================================================
!
subroutine ps_PseudoAtomSpec_Get(ps,atomic_symbol, atomic_label, &
     atomic_number, z_pseudo, pseudo_flavor,&
     relativity, spin_dft, core_corrections, annotation)
  
  type(ps_t), intent(in) :: ps
  
  character(len=*), intent(out), optional :: atomic_symbol
  character(len=*), intent(out), optional :: atomic_label
  real(dp), intent(out), optional         :: atomic_number
  real(dp), intent(out), optional         :: z_pseudo
  character(len=*), intent(out), optional :: pseudo_flavor
  character(len=*), intent(out), optional :: relativity
  logical, intent(out), optional          :: spin_dft
  logical, intent(out), optional          :: core_corrections
  type(ps_annotation_t), intent(out), optional :: annotation

  if (present(atomic_symbol)) then
     atomic_symbol = ps%header%atomic_label(1:2)
  endif
  
  if (present(atomic_label)) then
     atomic_label = trim(ps%header%atomic_label)
  endif

  if (present(atomic_number)) then
     atomic_number = ps%header%z
  endif
  
  if (present(z_pseudo)) then
     z_pseudo = ps%header%zpseudo
  endif
  
  if (present(pseudo_flavor)) then
     pseudo_flavor = trim(ps%header%flavor)
  endif

  if (present(relativity)) then
     relativity = trim(ps%header%relativity)
  endif

  if (present(spin_dft)) then
     spin_dft = ps%header%polarized
  endif

  if (present(core_corrections)) then
     core_corrections = (ps%header%core_corrections == "yes")
  endif

  if (present(annotation)) then
     annotation = ps%header%annotation
  endif

end subroutine ps_PseudoAtomSpec_Get
!
! ===================================================================
!
subroutine ps_ValenceConfiguration_Get(ps,nshells,charge,annotation)
  type(ps_t), intent(in) :: ps
  integer, intent(out), optional  :: nshells
  real(dp), intent(out), optional :: charge
  type(ps_annotation_t), intent(out), optional :: annotation

  if (present(nshells)) then
     nshells = ps%config_val%nshells
  endif
  if (present(charge)) then
     charge = ps%config_val%total_charge
  endif
  if (present(annotation)) then
     annotation = ps%config_val%annotation
  endif
end subroutine ps_ValenceConfiguration_Get

subroutine ps_ValenceShell_Get(ps,i,n,l,occupation,occ_up,occ_down)
  type(ps_t), intent(in) :: ps
  integer, intent(in)    :: i
  
  integer, intent(out), optional  :: n
  integer, intent(out), optional  :: l
  real(dp), intent(out), optional :: occupation 
  real(dp), intent(out), optional :: occ_up
  real(dp), intent(out), optional :: occ_down

  call check_index(i,ps%config_val%nshells,"valence shell")
  if (present(n)) then
     n =  ps%config_val%n(i)
  endif
  if (present(l)) then
     l = l_of_sym(ps%config_val%l(i),"valence shell")
  endif
  if (present(occupation)) then
     occupation = ps%config_val%occ(i)
  endif
  if (present(occ_up)) then
     if (ps%header%polarized) then
        occ_up = ps%config_val%occ_up(i)
     else
        call die("Cannot get per spin occupation")
     endif
  endif
  if (present(occ_down)) then
     if (ps%header%polarized) then
        occ_down = ps%config_val%occ_down(i)
     else
        call die("Cannot get per spin occupation")
     endif
  endif
end subroutine ps_ValenceShell_Get
!
! ===================================================================
!
subroutine ps_ValenceCharge_Get(ps,total_charge,&
     is_unscreening_charge, rescaled_to_z_pseudo,&
     annotation,func)

type(ps_t), intent(in) :: ps

real(dp), intent(out), optional              :: total_charge
character(len=*), intent(out), optional      :: is_unscreening_charge
character(len=*), intent(out), optional      :: rescaled_to_z_pseudo
type(ps_annotation_t), intent(out), optional :: annotation
type(ps_radfunc_t), intent(out), optional    :: func

if (present(total_charge)) then
   total_charge = ps%valence_charge%total_charge
endif

if (present(is_unscreening_charge)) then
   is_unscreening_charge = ps%valence_charge%is_unscreening_charge
endif

if (present(rescaled_to_z_pseudo)) then
   rescaled_to_z_pseudo = ps%valence_charge%rescaled_to_z_pseudo
endif

if (present(annotation)) then
   annotation = ps%valence_charge%annotation
endif

if (present(func)) then
   func = ps%valence_charge%rho_val
endif

end subroutine ps_ValenceCharge_Get
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

val = ps_GetValue(ps%valence_charge%rho_val,r)

end function ps_ValenceCharge_Value
!
! ===================================================================
!
subroutine ps_CoreCharge_Get(ps,rc,nderivs,annotation,func)

type(ps_t), intent(in) :: ps

real(dp), intent(out), optional              :: rc
integer , intent(out), optional              :: nderivs
type(ps_annotation_t), intent(out), optional :: annotation
type(ps_radfunc_t), intent(out), optional    :: func

if (present(rc)) then
   rc = ps%core_charge%rcore
endif
if (present(nderivs)) then
   nderivs = ps%core_charge%n_cont_derivs
endif
if (present(annotation)) then
   annotation = ps%core_charge%annotation
endif
if (present(func)) then
   func = ps%core_charge%rho_core
endif

end subroutine ps_CoreCharge_Get

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

val = ps_GetValue(ps%core_charge%rho_core,r)

end function ps_CoreCharge_Value
!
! ===================================================================
!
subroutine ps_ExchangeCorrelation_Get(ps,annotation,n_libxc_functionals)
type(ps_t), intent(in) :: ps

type(ps_annotation_t), intent(out), optional :: annotation
integer, intent(out), optional               :: n_libxc_functionals

if (present(annotation)) then
   annotation = ps%xc_info%annotation
endif
if (present(n_libxc_functionals)) then
   n_libxc_functionals =  ps%xc_info%n_functs_libxc
endif
end subroutine ps_ExchangeCorrelation_Get

subroutine ps_LibxcFunctional_Get(ps,i,name,code,type,weight)
  type(ps_t), intent(in) :: ps

  integer, intent(in)      :: i
  character(len=*), intent(out), optional :: name
  integer, intent(out), optional          :: code
  character(len=*), intent(out), optional :: type
  real(dp), intent(out), optional         :: weight

  call check_index(i,ps%xc_info%n_functs_libxc,"libxc functional")

  if (present(name)) then
     name = ps%xc_info%libxc_name(i)
  endif
  if (present(code)) then
     code = ps%xc_info%libxc_id(i)
  endif
  if (present(type)) then
     type = ps%xc_info%libxc_type(i)
  endif
  if (present(weight)) then
     weight = ps%xc_info%libxc_weight(i)
  endif
end subroutine ps_LibxcFunctional_Get
!
!=================================================
!> Returns the annotation associated to a grid.
!> If a radial function
!> handle is given, the annotation for that 
!> radial function's grid is returned. Otherwise,
!> the return value is the annotation for the global grid.
!> If there is no appropriate annotation, an empty
!> structure is returned.
!> @param ps is a handle to the psml information
!> @param radfunc is a handle to a radial function structure
!>
function ps_GridAnnotation(ps,radfunc) result(annotation)

 type(ps_t), intent(in)                 :: ps
 type(ps_radfunc_t), intent(in), optional  :: radfunc
 
 type(ps_annotation_t)  :: annotation
 type(ps_annotation_t), pointer  :: annotation_p

 if (present(radfunc)) then
    ! We are told to get the grid annotation
    ! for a specific radial function
    
    if (.not. initialized(radfunc%grid)) then
       call die("get_annotation: Invalid radial function")
    endif
    annotation_p => annotationGrid(radfunc%grid)
    ! If npts_data /= npts_grid, add a record to reflect this
    annotation = annotation_p

 else

    ! This is the global grid annotation
    if (.not. initialized(ps%global_grid)) then
       annotation = EMPTY_ANNOTATION
    else
       annotation_p => annotationGrid(ps%global_grid)
       annotation = annotation_p
    endif

 endif
end function ps_GridAnnotation

!
!====================================================
! Semilocal potentials
!
subroutine ps_Potential_Get(ps,i,&
            l,j,n,rc,set,flavor,eref,annotation,func)
type(ps_t), intent(in), target            :: ps
integer,   intent(in)             :: i
integer, intent(out), optional    :: set
integer, intent(out), optional    :: l
real(dp), intent(out), optional    :: j
integer, intent(out), optional    :: n
real(dp), intent(out), optional    :: rc
real(dp), intent(out), optional    :: eref
character(len=*), intent(out), optional    :: flavor
type(ps_annotation_t), intent(out), optional  :: annotation
type(ps_radfunc_t), intent(out), optional    :: func

type(sl_table_t), pointer :: q(:)
q => ps%sl_table

call check_index(i,size(q),"SL pot")
if (present(set)) then
   set = q(i)%p%set
endif
if (present(l)) then
   l = l_of_sym(q(i)%p%l,"SL pot")
endif
if (present(n)) then
   n = q(i)%p%n
endif
if (present(j)) then
   if (q(i)%p%j < 0.0) then
      j = -1.0_dp
      ! Maybe optional status flag raised?
   else
      j = q(i)%p%j
   endif
endif
if (present(rc)) then
   rc = q(i)%p%rc
endif

if (present(eref)) then
   ! Will return a very large positive value if the attribute eref is not
   ! present in the file
   eref = q(i)%p%eref
endif

if (present(flavor)) then
   flavor = q(i)%p%flavor
endif
if (present(annotation)) then
   annotation = q(i)%p%parent_group%annotation
endif
if (present(func)) then
   func = q(i)%p%V
endif
end subroutine ps_Potential_Get
!
subroutine ps_SemilocalPotentials_Filter(ps,&
                    indexes_in, &
                    l,j,n,set, &
                    indexes,number)
type(ps_t), intent(in), target   :: ps
integer, intent(in), optional    :: indexes_in(:)
integer, intent(in), optional    :: set
integer, intent(in), optional    :: l
real(dp), intent(in), optional   :: j
integer, intent(in), optional    :: n
integer,   intent(out), allocatable, optional  :: indexes(:)
integer,   intent(out), optional  :: number

integer :: i, ii, num
integer, allocatable :: idx(:), range(:)

type(sl_table_t), pointer :: q(:)
q => ps%sl_table

!
!  If we specify an initial domain, we
!  only loop over it for the rest of
!  the operations
!

if (present(indexes_in)) then
   allocate(range(size(indexes_in)))
   range(:) = indexes_in(:)
else
   allocate(range(size(q)))
   do ii = 1, size(q)
      range(ii) = ii
   enddo
endif

num = 0
do ii = 1, size(range)
   i = range(ii)
   if (present(set)) then
      if (iand(q(i)%p%set,set) == 0) cycle
   endif
   if (present(l)) then
      if (l_of_sym(q(i)%p%l,"pot") /= l) cycle
   endif
   if (present(j)) then
      ! A negative p%j signals the absence of j...
      if (q(i)%p%j < 0.0) cycle
      ! Perform integer comparison
      if (nint(10*q(i)%p%j) /= nint(10*j)) cycle
   endif
   if (present(n)) then
      if (q(i)%p%n /= n) cycle
   endif
   num = num + 1
enddo

allocate(idx(num))

num = 0
do ii = 1, size(range)
   i = range(ii)
   if (present(set)) then
      if (iand(q(i)%p%set,set) == 0) cycle
   endif
   if (present(j)) then
      ! A negative p%j signals the absence of j...
      if (q(i)%p%j < 0.0) cycle
      ! Perform integer comparison
      if (nint(10*q(i)%p%j) /= nint(10*j)) cycle
   endif
   if (present(l)) then
      if (l_of_sym(q(i)%p%l,"pot") /= l) cycle
   endif
   if (present(n)) then
      if (q(i)%p%n /= n) cycle
   endif
   num = num + 1
   idx(num) = i
enddo

if (present(indexes)) then
   if (allocated(indexes))  deallocate(indexes)
   allocate(indexes(num))
   indexes(:) = idx(:)
endif
if (present(number)) then
   number = num
endif

deallocate(idx,range)

end subroutine ps_SemilocalPotentials_Filter
!
function ps_GetValue(f,r) result(val)
type(ps_radfunc_t), intent(in)     ::  f
real(dp),  intent(in)      :: r
real(dp)                   :: val

if (r> max_range(f)) then
   if (f%has_coulomb_tail) then
      val = f%tail_factor / r
   else
      val = 0.0_dp
   endif
else
   val = eval_radfunc(f,r, debug=global_debug)
endif
end function ps_GetValue

!> Evaluator by storage index
function ps_Potential_Value(ps,i,r) result(val)
type(ps_t), intent(in)     ::  ps
integer,   intent(in)      :: i
real(dp),  intent(in)      :: r
real(dp)                   :: val

call check_index(i,size(ps%sl_table),"SL pot")
val = ps_GetValue(ps%sl_table(i)%p%V,r)

end function ps_Potential_Value

subroutine ps_GetRawData(f,raw_r,raw_data)
type(ps_radfunc_t), intent(in) :: f
real(dp), allocatable, intent(out)  :: raw_r(:), raw_data(:)

integer npts_data
real(dp), pointer :: a(:) 

! Cover the case in which the data set uses only
! a first section of the grid

npts_data = size(f%data)
allocate(raw_r(npts_data), raw_data(npts_data))

a => valGrid(f%grid)
raw_r(1:npts_data) = a(1:npts_data)
raw_data(:) = f%data(:)

end subroutine ps_GetRawData

!
!====================================================
! PseudoWavefunctions
!
! Basic accessors
!
subroutine ps_PseudoWaveFunctions_Filter(ps,&
                    indexes_in,&
                    l,j,n,set, &
                    indexes,number)
type(ps_t), intent(in), target   :: ps
integer,  intent(in), optional   :: indexes_in(:)
integer, intent(in), optional    :: set
integer, intent(in), optional    :: l
real(dp), intent(in), optional   :: j
integer, intent(in), optional    :: n
integer,   intent(out), allocatable, optional  :: indexes(:)
integer,   intent(out), optional  :: number

integer :: i, ii, num
integer, allocatable :: idx(:), range(:)

type(wf_table_t), pointer :: q(:)
q => ps%wf_table

!
!  If we specify an initial domain, we
!  only loop over it for the rest of
!  the operations
!

if (present(indexes_in)) then
   allocate(range(size(indexes_in)))
   range(:) = indexes_in(:)
else
   allocate(range(size(q)))
   do ii = 1, size(q)
      range(ii) = ii
   enddo
endif

num = 0
do ii = 1, size(range)
   i = range(ii)
   if (present(set)) then
      if (iand(q(i)%p%set,set) == 0) cycle
   endif
   if (present(l)) then
      if (l_of_sym(q(i)%p%l,"wf") /= l) cycle
   endif
   if (present(j)) then
      ! A negative p%j signals the absence of j...
      if (q(i)%p%j < 0.0) cycle
      ! Perform integer comparison
      if (nint(10*q(i)%p%j) /= nint(10*j)) cycle
   endif
   if (present(n)) then
      if (q(i)%p%n /= n) cycle
   endif
   num = num + 1
enddo

allocate(idx(num))

num = 0
do ii = 1, size(range)
   i = range(ii)
   if (present(set)) then
      if (iand(q(i)%p%set,set) == 0) cycle
   endif
   if (present(j)) then
      ! A negative p%j signals the absence of j...
      if (q(i)%p%j < 0.0) cycle
      ! Perform integer comparison
      if (nint(10*q(i)%p%j) /= nint(10*j)) cycle
   endif
   if (present(l)) then
      if (l_of_sym(q(i)%p%l,"wf") /= l) cycle
   endif
   if (present(n)) then
      if (q(i)%p%n /= n) cycle
   endif
   num = num + 1
   idx(num) = i
enddo

if (present(indexes)) then
   if (allocated(indexes))  deallocate(indexes)
   allocate(indexes(num))
   indexes(:) = idx(:)
endif
if (present(number)) then
   number = num
endif

deallocate(idx,range)

end subroutine ps_PseudoWaveFunctions_Filter

subroutine ps_PseudoWf_Get(ps,i,&
            l,j,n,set,energy_level,annotation,func)
type(ps_t), intent(in), target    :: ps
integer,   intent(in)             :: i
integer, intent(out), optional    :: set
integer, intent(out), optional    :: l
real(dp), intent(out), optional   :: j
integer, intent(out), optional    :: n
real(dp), intent(out), optional   :: energy_level
type(ps_annotation_t), intent(out), optional :: annotation
type(ps_radfunc_t), intent(out), optional    :: func

type(wf_table_t), pointer :: q(:)
q => ps%wf_table

call check_index(i,size(q),"wf")
if (present(set)) then
   set = q(i)%p%set
endif
if (present(l)) then
   l = l_of_sym(q(i)%p%l,"wf")
endif
if (present(n)) then
   n = q(i)%p%n
endif
if (present(j)) then
   if (q(i)%p%j < 0.0) then
      j = -1.0_dp
      ! Maybe optional status flag raised?
   else
      j = q(i)%p%j
   endif
endif

if (present(energy_level)) then
   ! Will return a very large positive value if the attribute energy_level is not
   ! present in the file
   energy_level = q(i)%p%energy_level
endif


if (present(annotation)) then
   annotation = q(i)%p%parent_group%annotation
endif

if (present(func)) then
   func = q(i)%p%Phi
endif
end subroutine ps_PseudoWf_Get

!
function ps_PseudoWf_Value(ps,i,r) result(val)
type(ps_t), intent(in) :: ps
integer,   intent(in)  :: i
real(dp),  intent(in)  :: r
real(dp)               :: val

call check_index(i,size(ps%wf_table),"Wf")
val = ps_GetValue(ps%wf_table(i)%p%Phi,r)

end function ps_PseudoWf_Value
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

subroutine ps_LocalPotential_Get(ps,type,annotation,func,&
                                 has_local_charge,func_local_charge)
type(ps_t), intent(in) :: ps

character(len=*), intent(out), optional      :: type
type(ps_annotation_t), intent(out), optional :: annotation
type(ps_radfunc_t), intent(out), optional    :: func
logical, intent(out), optional               :: has_local_charge
type(ps_radfunc_t), intent(out), optional    :: func_local_charge

if (present(type)) then
   type = ps%local%vlocal_type
endif
if (present(annotation)) then
   annotation = ps%local%annotation
endif
if (present(func)) then
   func = ps%local%vlocal
endif
if (present(has_local_charge)) then
   has_local_charge = (initialized(ps%local%chlocal%grid))
endif
if (present(func_local_charge)) then
   ! No checks
   func_local_charge = ps%local%chlocal
endif
end subroutine ps_LocalPotential_Get
!
function ps_LocalPotential_Value(ps,r) result(val)
type(ps_t), intent(in) :: ps
real(dp), intent(in)       :: r
real(dp)                   :: val

val = ps_GetValue(ps%local%vlocal,r)

end function ps_LocalPotential_Value
!
! ==========================================================
!
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

!================================ Projectors
subroutine ps_NonlocalProjectors_Filter(ps,&
                    indexes_in,&
                    l,j,seq,set, &
                    indexes,number)
type(ps_t), intent(in), target   :: ps
integer, intent(in), optional    :: indexes_in(:)
integer, intent(in), optional    :: set
integer, intent(in), optional    :: l
real(dp), intent(in), optional   :: j
integer, intent(in), optional    :: seq
integer,   intent(out), allocatable, optional  :: indexes(:)
integer,   intent(out), optional  :: number

integer :: i, ii, num
integer, allocatable :: idx(:), range(:)

type(nl_table_t), pointer :: q(:)
q => ps%nl_table

!
!  If we specify an initial domain, we
!  only loop over it for the rest of
!  the operations
!

if (present(indexes_in)) then
   allocate(range(size(indexes_in)))
   range(:) = indexes_in(:)
else
   allocate(range(size(q)))
   do ii = 1, size(q)
      range(ii) = ii
   enddo
endif

num = 0
do ii = 1, size(range)
   i = range(ii)
   if (present(set)) then
      if (iand(q(i)%p%set,set) == 0) cycle
   endif
   if (present(l)) then
      if (l_of_sym(q(i)%p%l,"pot") /= l) cycle
   endif
   if (present(j)) then
      ! A negative p%j signals the absence of j...
      if (q(i)%p%j < 0.0) cycle
      ! Perform integer comparison
      if (nint(10*q(i)%p%j) /= nint(10*j)) cycle
   endif
   if (present(seq)) then
      if (q(i)%p%seq /= seq) cycle
   endif
   num = num + 1
enddo

allocate(idx(num))

num = 0
do ii = 1, size(range)
   i = range(ii)
   if (present(set)) then
      if (iand(q(i)%p%set,set) == 0) cycle
   endif
   if (present(j)) then
      ! A negative p%j signals the absence of j...
      if (q(i)%p%j < 0.0) cycle
      ! Perform integer comparison
      if (nint(10*q(i)%p%j) /= nint(10*j)) cycle
   endif
   if (present(l)) then
      if (l_of_sym(q(i)%p%l,"pot") /= l) cycle
   endif
   if (present(seq)) then
      if (q(i)%p%seq /= seq) cycle
   endif
   num = num + 1
   idx(num) = i
enddo

if (present(indexes)) then
   if (allocated(indexes))  deallocate(indexes)
   allocate(indexes(num))
   indexes(:) = idx(:)
endif
if (present(number)) then
   number = num
endif

deallocate(idx,range)

end subroutine ps_NonlocalProjectors_Filter

subroutine ps_Projector_Get(ps,i,&
            l,j,seq,set,ekb,eref,type,annotation,func)
type(ps_t), intent(in), target    :: ps
integer,   intent(in)             :: i
integer, intent(out), optional    :: set
integer, intent(out), optional    :: l
real(dp), intent(out), optional    :: j
integer, intent(out), optional    :: seq
real(dp), intent(out), optional    :: ekb
real(dp), intent(out), optional    :: eref
character(len=*), intent(out), optional    :: type
type(ps_annotation_t), intent(out), optional :: annotation
type(ps_radfunc_t), intent(out), optional    :: func

type(nl_table_t), pointer :: q(:)
q => ps%nl_table

call check_index(i,size(q),"NL pot")
if (present(set)) then
   set = q(i)%p%set
endif
if (present(l)) then
   l = l_of_sym(q(i)%p%l,"NL pot")
endif
if (present(seq)) then
   seq = q(i)%p%seq
endif
if (present(j)) then
   if (q(i)%p%j < 0.0) then
      j = -1.0_dp
      ! Maybe optional status flag raised?
   else
      j = q(i)%p%j
   endif
endif
if (present(ekb)) then
   ekb = q(i)%p%ekb
endif
if (present(eref)) then
   ! Will return a very large positive value if the attribute eref is not
   ! present in the file
   eref = q(i)%p%eref
endif
if (present(type)) then
   type = q(i)%p%type
endif
if (present(annotation)) then
   annotation = q(i)%p%parent_group%annotation
endif
if (present(func)) then
   func = q(i)%p%proj
endif
end subroutine ps_Projector_Get


function ps_Projector_Value(ps,i,r) result(val)
!+ display: private
type(ps_t), intent(in) :: ps
integer,   intent(in)  :: i
real(dp),  intent(in)  :: r
real(dp)               :: val

call check_index(i,size(ps%nl_table),"proj")
val = ps_GetValue(ps%nl_table(i)%p%proj,r)

end function ps_Projector_Value

!====================================================
! Low-level routines
!
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

!>  Returns the maximum radius in a radfunc's data
function max_range(f) result(range)
type(radfunc_t), intent(in) :: f
real(dp)                  :: range

real(dp),  pointer :: a(:)
integer :: npts_data

a => valGrid(f%grid)

!
if (global_use_effective_range) then
   ! We use the effective end_of_range
   range = f%rcut_eff
else
   ! Use the nominal range
   ! This covers the case in which the data set uses only
   ! a first section of the grid
   npts_data = size(f%data)
   range = a(npts_data)
endif

end function max_range

!----------
function eval_radfunc(f,r,debug) result(val)
use m_psml_interp, only: interpolator, nq

type(radfunc_t), intent(in) :: f
real(dp), intent(in)      :: r
logical, intent(in)       :: debug

real(dp)                  :: val

real(dp), pointer :: x(:) => null(), y(:) => null()
integer :: npts

x => valGrid(f%grid)
y => f%data(:)

!Note size(y) to cover the case in which the data set uses only
! a first section of the grid
npts = size(y)
call interpolator(nq,x,y,npts,r,val,debug)

end function eval_radfunc
!
!------
!
   subroutine assert(cond,message)
     logical, intent(in) :: cond
     character(len=*) message

     if (.not. cond) call die(message)
   end subroutine assert

end module m_psml_api






