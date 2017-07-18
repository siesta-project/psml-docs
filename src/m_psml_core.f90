  !+ graph: false
  !+ author: Alberto Garcia

  !+ Data structures to handle the PSML pseudopotential format.

module m_psml_core

use iso_varying_string, only: varying_string, var_str

use assoc_list, only: ps_annotation_t => assoc_list_t
use assoc_list, only: ps_clean_annotation => assoc_list_reset
use class_Grid

use external_interfaces, only: die => psml_die

implicit none

private

!
! Update this. Up to 99...
integer, parameter    :: PATCH_LEVEL = 4
!
! Only update 1000 when changing major/minor version
integer, parameter, public  :: PSML_LIBRARY_VERSION = 1100 + PATCH_LEVEL
!
!  Simple sanity checks while the format evolves
!  This version is able to read v1.0 and v1.1 PSML files
!
!  Note that the version is really given by the generators.
!
! The "hi" value is intended to auto-revoke the library, but
! it is neither completely foolproof nor flexible enough.
!
real, parameter, public  :: PSML_TARGET_VERSION_LO = 1.00
real, parameter, public  :: PSML_TARGET_VERSION_HI = 1.10

!----------------------------------------------------------------
! Hardwired parameters (to be made dynamical in a later version)

! Maximum number of valence shells (including semicore) or core shells:
integer, parameter, private    :: MAXN_SHELLS = 20
!----------------------------------------------------------------

integer, parameter, private    :: dp = selected_real_kind(14)
!
!-----------------------------------------------------------

type, public :: input_file_t
        character(len=100)      :: name = "-----"
        type(varying_string)    :: buffer
end type input_file_t

!------
type, public :: provenance_t
   type(provenance_t), pointer  :: prev  => null()
        integer                 :: record_number
        character(len=100)      :: creator = "-----"
        character(len=60)       :: date    = "-----"
        integer                 :: n_input_files = 0  ! Max 1 for now !!
        type(input_file_t)      :: input_file
        type(ps_annotation_t)   :: annotation
   type(provenance_t), pointer  :: next => null()
end type provenance_t
!------
type, public :: header_t
   ! This is the 'pseudo-atom-spec' section
   
        character(len=100)      :: atomic_label    !! generalized symbol
        real(kind=dp)           :: z  !! atomic number (might be non-integer)
        real(kind=dp)           :: zpseudo !! Z - ncore-electrons
        character(len=100)      :: flavor  !! pseudization method
        character(len=6)        :: relativity !! "no|scalar|dirac"
        logical                 :: polarized !! is spin-DFT?
        !
        character(len=3)        :: core_corrections !! are there NLCC's?
      !
      type(ps_annotation_t)   :: annotation
end type header_t
!------
type, public :: config_val_t
      integer                          :: nshells
      real(kind=dp)                    :: total_charge
      integer, dimension(MAXN_SHELLS)  :: n
      character(len=1), dimension(MAXN_SHELLS) :: l
      real(dp), dimension(MAXN_SHELLS) :: occ
      real(dp), dimension(MAXN_SHELLS) :: occ_up
      real(dp), dimension(MAXN_SHELLS) :: occ_down
      !
      type(ps_annotation_t)   :: annotation
end type config_val_t
!------
type, public :: xc_t
        integer                         :: n_functs_libxc = 0
        character(len=200), allocatable :: libxc_name(:)
        character(len=100), allocatable :: libxc_type(:)
        integer, allocatable            :: libxc_id(:)
        real(dp), allocatable           :: libxc_weight(:)
        type(ps_annotation_t)           :: annotation
end type xc_t
!------
type, public :: radfunc_t
      type(Grid)                              :: grid
      real(kind=dp), dimension(:), pointer    :: data => null()
      logical                                 :: has_coulomb_tail
      real(dp)                                :: tail_factor = 0.0_dp
      integer                                 :: nnz ! # of 'non-zero' values
      real(dp)                                :: rcut_eff ! effective end of range
 
end type radfunc_t      
!
!===============================================
type, public :: slps_t
      integer           :: n
      character(len=1)  :: l
      real(dp)          :: j = -1.0_dp
      integer           :: set
      character(len=100):: flavor
      real(dp)          :: rc
      real(dp)          :: eref  ! Reference energy
      type(radfunc_t)   :: V
      type(semilocal_t), pointer :: parent_group => null()
   type(slps_t), pointer :: next => null()

end type slps_t

type, public :: sl_table_t
   type(slps_t), pointer :: p => null()
end type sl_table_t

type, public :: semilocal_t
   type(slps_t), pointer :: pot => null()
   integer           :: set
   !
   ! Optional private grid
   !
   type(Grid)             :: grid
   type(ps_annotation_t)  :: annotation
   !
   type(semilocal_t), pointer     :: next => null()

end type semilocal_t
!===============================================

type, public :: local_t
   !
      type(ps_annotation_t)                     :: annotation

   ! Optional private grid
   !
      type(Grid)                               :: grid

      type(radfunc_t)                          :: Vlocal
      character(len=100)                       :: vlocal_type

      type(radfunc_t)                          :: Chlocal
end type local_t

!===============================================
type, public :: nlpj_t
      integer           :: seq
      character(len=1)  :: l
      real(dp)          :: j = -1.0_dp
      integer           :: set
      character(len=100):: type
      real(dp)          :: ekb
      real(dp)          :: eref  ! Reference energy
      type(radfunc_t)   :: proj

      type(nonlocal_t), pointer :: parent_group => null()
   type(nlpj_t), pointer :: next => null()

end type nlpj_t

type, public :: nl_table_t
   type(nlpj_t), pointer :: p => null()
end type nl_table_t

type, public :: nonlocal_t
   type(nlpj_t), pointer :: proj => null()
   integer           :: set
   !
   ! Optional private grid
   !
   type(Grid)                               :: grid
   type(ps_annotation_t)  :: annotation
   !
   type(nonlocal_t), pointer     :: next => null()

end type nonlocal_t
!===============================================
! Wavefunctions
!
type, public :: wf_t
      integer           :: n
      character(len=1)  :: l
      integer           :: set
      real(dp)          :: j = -1.0_dp
      real(dp)          :: energy_level
      type(radfunc_t)   :: Phi
      type(wfns_t), pointer :: parent_group => null()
   type(wf_t), pointer :: next => null()

end type wf_t

type, public :: wf_table_t
   type(wf_t), pointer :: p => null()
end type wf_table_t

type, public :: wfns_t
   type(wf_t), pointer :: wf => null()
   integer             :: set
   character(len=100)  :: type = ""
   !
   ! Optional private grid
   !
   type(Grid)             :: grid
   type(ps_annotation_t)  :: annotation
   !
   type(wfns_t), pointer     :: next => null()

end type wfns_t

type, public :: valence_charge_t
      real(dp)        :: total_charge
      character(len=3):: is_unscreening_charge = ""
      character(len=3):: rescaled_to_z_pseudo  = ""

      type(radfunc_t) :: rho_val
      type(ps_annotation_t)   :: annotation
end type valence_charge_t

type, public :: core_charge_t
      integer         :: n_cont_derivs
      real(dp)        :: rcore
      type(radfunc_t) :: rho_core

      type(ps_annotation_t)           :: annotation
end type core_charge_t


type, public :: ps_t
!! Main derived type to hold the PSML information
      character(len=10)                  :: version     = ""
      character(len=40)                  :: energy_unit = ""
      character(len=40)                  :: length_unit = ""
      character(len=36)                  :: uuid = ""
      character(len=200)                 :: namespace = ""
      type(ps_annotation_t)              :: annotation   ! V1.0 only
      type(provenance_t), pointer        :: provenance => null()
      type(header_t)                     :: header     ! pseudo-atom-spec
      type(config_val_t)                 :: config_val
      type(config_val_t)                 :: config_core  ! extension
      type(xc_t)                         :: xc_info
      type(Grid)                         :: global_grid
      type(local_t)                      :: local
      type(semilocal_t), pointer         :: semilocal => null()
      type(nonlocal_t), pointer          :: nonlocal => null()
      type(wfns_t), pointer              :: wavefunctions => null()
      !
      type(valence_charge_t)             :: valence_charge
      type(core_charge_t)                :: core_charge
      !
      ! index tables
      !
      type(sl_table_t), allocatable      :: sl_table(:)
      type(nl_table_t), allocatable      :: nl_table(:)
      type(wf_table_t), allocatable      :: wf_table(:)

   end type ps_t

   integer,  parameter, public   &
                          :: SET_NULL     =   0, &
                             SET_SREL     =   1, &
                             SET_NONREL   =   2, &
                             SET_SO       =   4, &
                             SET_LJ       =   8, &
                             SET_UP       =  16, &
                             SET_DOWN     =  32, &
                             SET_SPINAVE  =  64, &
                             SET_SPINDIFF = 128, &         ! 2^7
                             SET_USER1    = 256, &         ! 2^8
                             SET_USER2    = 512            ! 2^9

   integer, parameter, public    :: SET_ALL =  2**10 -1

 public  :: ps_destroy
 public  :: str_of_set

 public  :: setcode_of_string    ! utility function, not for client normal use

 public  :: destroy_local
 public  :: destroy_nonlocal
 public  :: destroy_wavefunctions
 
 CONTAINS

subroutine ps_destroy(ps)
!! Cleans the ps object
type(ps_t), intent(inout)     :: ps

integer :: i

call ps_clean_annotation(ps%annotation)

call destroy_provenance(ps%provenance)

call ps_clean_annotation(ps%header%annotation)
call ps_clean_annotation(ps%config_val%annotation)
call ps_clean_annotation(ps%config_core%annotation)
call destroy_xc(ps%xc_info)
!
! Note that freshly declared objects must have
! npots = 0 and npswfs = 0 !
!
call destroy_semilocal(ps%semilocal)
call destroy_nonlocal(ps%nonlocal)
!
call destroy_local(ps%local)
!
call destroy_wavefunctions(ps%wavefunctions)
!
call destroy_radfunc(ps%valence_charge%rho_val)
call ps_clean_annotation(ps%valence_charge%annotation)
!
call destroy_radfunc(ps%core_charge%rho_core)
call ps_clean_annotation(ps%core_charge%annotation)

!
call delete(ps%global_grid)

end subroutine ps_destroy

subroutine destroy_provenance(p)
type(provenance_t), pointer :: p

type(provenance_t), pointer :: q

do while (associated(p))
   call ps_clean_annotation(p%annotation)
   ! clean buffers for input files?
   q => p%next
   deallocate(p)
   p => q
enddo

end subroutine destroy_provenance

!==================================================
subroutine destroy_semilocal(p)
type(semilocal_t), pointer :: p

type(semilocal_t), pointer :: q

do while (associated(p))
   call ps_clean_annotation(p%annotation)
   call destroy_slps(p%pot)
   call delete(p%grid)
   q => p%next
   deallocate(p)
   p => q
enddo

end subroutine destroy_semilocal
!
subroutine destroy_slps(p)
type(slps_t), pointer :: p

type(slps_t), pointer :: q

do while (associated(p))
   call destroy_radfunc(p%V)
   q => p%next
   deallocate(p)
   p => q
enddo

end subroutine destroy_slps
!==================================================
subroutine destroy_local(p)
  type(local_t)  :: p
  
call ps_clean_annotation(p%annotation)
call delete(p%grid)
call destroy_radfunc(p%vlocal)
call destroy_radfunc(p%chlocal)
end subroutine destroy_local
!==================================================
subroutine destroy_nonlocal(p)
type(nonlocal_t), pointer :: p

type(nonlocal_t), pointer :: q

do while (associated(p))
   call ps_clean_annotation(p%annotation)
   call destroy_nlpj(p%proj)
   call delete(p%grid)
   q => p%next
   deallocate(p)
   p => q
enddo

end subroutine destroy_nonlocal
!
subroutine destroy_nlpj(p)
type(nlpj_t), pointer :: p

type(nlpj_t), pointer :: q

do while (associated(p))
   call destroy_radfunc(p%proj)
   q => p%next
   deallocate(p)
   p => q
enddo

end subroutine destroy_nlpj
!
!==================================================
subroutine destroy_wavefunctions(p)
type(wfns_t), pointer :: p

type(wfns_t), pointer :: q

do while (associated(p))
   call ps_clean_annotation(p%annotation)
   call destroy_pswf(p%wf)
   call delete(p%grid)
   q => p%next
   deallocate(p)
   p => q
enddo

end subroutine destroy_wavefunctions
subroutine destroy_pswf(p)
type(wf_t), pointer :: p

type(wf_t), pointer :: q

do while (associated(p))
   call destroy_radfunc(p%Phi)
   q => p%next
   deallocate(p)
   p => q
enddo

end subroutine destroy_pswf

!==================================================
subroutine destroy_radfunc(rp)
type(radfunc_t) :: rp

call delete(rp%grid)
if (associated(rp%data)) then
   deallocate(rp%data)
   rp%data => null()
endif
end subroutine destroy_radfunc

!
subroutine destroy_xc(xp)
type(xc_t), intent(inout) :: xp

if (allocated(xp%libxc_name)) deallocate(xp%libxc_name)
if (allocated(xp%libxc_type)) deallocate(xp%libxc_type)
if (allocated(xp%libxc_id)) deallocate(xp%libxc_id)
if (allocated(xp%libxc_weight)) deallocate(xp%libxc_weight)
call ps_clean_annotation(xp%annotation)

end subroutine destroy_xc

function setcode_of_string(str) result(code)
       character(len=*), intent(in) :: str
       integer                      :: code

       select case (trim(str))
       case ("non_relativistic")
          code = SET_NONREL
       case ("scalar_relativistic")
          code = SET_SREL
       case ("spin_orbit")
          code = SET_SO
       case ("lj")
          code = SET_LJ
       case ("spin_up")
          code = SET_UP
       case ("spin_down")
          code = SET_DOWN
       case ("spin_average")
          code = SET_SPINAVE
       case ("spin_difference")
          code = SET_SPINDIFF
       case ("user_extension1")
          code = SET_USER1
       case ("user_extension2")
          code = SET_USER2
       case ("all","any")
          code = SET_ALL
       case ("invalid","INVALID")
          code = SET_NULL
       case default
          call die("Wrong set string: "//trim(str))
       end select

end function setcode_of_string

function str_of_set(code) result(str)
       integer, intent(in)          :: code
       character(len=40)            :: str

       character(len=100) :: msg

       select case (code)
       case (SET_NONREL)
          str ="non_relativistic"
       case (SET_SREL)
          str ="scalar_relativistic"
       case (SET_SO)
          str ="spin_orbit"
       case (SET_LJ)
          str ="lj"
       case (SET_UP)
          str ="spin_up"
       case (SET_DOWN)
          str ="spin_down"
       case (SET_SPINAVE)
          str ="spin_average"
       case (SET_SPINDIFF)
          str ="spin_difference"
       case (SET_USER1)
          str ="user_extension1"
       case (SET_USER2)
          str ="user_extension2"
       case (SET_ALL)
          str ="all"
       case (SET_NULL)
          str ="invalid"
       case default
          write(msg,"(a,i4)") "Wrong set code: ", code
          call die(msg)
       end select

end function str_of_set

end module m_psml_core






