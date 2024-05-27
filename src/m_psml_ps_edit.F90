!> Functions to edit the PSML ps_t structure
!!
!> @author Alberto Garcia
!
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module m_psml_ps_edit

use m_psml_core   ! For basic structures
use m_psml_external_interfaces, only: die => psml_die

use m_psml_assoc_list, ps_annotation_t => assoc_list_t

implicit none

public :: ps_RootAttributes_Set
public :: ps_SetUUID
public :: ps_SetPSMLVersion
public :: ps_AddProvenanceRecord
public :: ps_Delete_LocalPotential
public :: ps_Delete_NonLocalProjectors

! Aliases
interface ps_Provenance_Add
   module procedure ps_AddProvenanceRecord
end interface ps_Provenance_Add
public :: ps_Provenance_Add

interface ps_NonLocalProjectors_Delete
   module procedure ps_Delete_NonLocalProjectors
end interface ps_NonLocalProjectors_Delete
public :: ps_NonLocalProjectors_Delete

interface ps_LocalPotential_Delete
   module procedure ps_Delete_LocalPotential
end interface ps_LocalPotential_Delete
public :: ps_LocalPotential_Delete

!public :: ps_AddSLBlock
private

CONTAINS !===============================================
  subroutine ps_RootAttributes_Set(ps,version,uuid,namespace) 
    type(ps_t), intent(inout) :: ps
    character(len=*), intent(in), optional :: version
    character(len=*), intent(in), optional :: uuid
    character(len=*), intent(in), optional :: namespace

    if (present(version)) then
       ps%version = version
    endif

    if (present(uuid)) then
       ps%uuid = uuid
    endif
    
    if (present(namespace)) then
       ps%namespace = namespace
    endif
    
  end  subroutine ps_RootAttributes_Set
  
  subroutine ps_SetPSMLVersion(ps,version) 
    type(ps_t), intent(inout) :: ps
    character(len=*), intent(in) :: version
    ps%version = version
  end subroutine ps_SetPSMLVersion
  
  subroutine ps_SetUUID(ps,id) 
    type(ps_t), intent(inout) :: ps
    character(len=36), intent(in) :: id
    ps%uuid = id
  end subroutine ps_SetUUID
  
  subroutine ps_AddProvenanceRecord(ps,creator,date,annotation)

    type(ps_t), intent(inout) :: ps
    character(len=*), intent(in) :: creator
    character(len=*), intent(in) :: date
    type(ps_annotation_t), intent(in), target :: annotation

type(provenance_t), pointer :: p
type(provenance_t), pointer :: q
integer :: depth

allocate(p)
!
! Find the depth of the provenance stack
!
depth = 0
q => ps%provenance
do while (associated(q))
   depth = depth + 1
   q => q%next
enddo

q => ps%provenance
if (associated(q)) then
   p%next => q
   q%prev => p
endif
ps%provenance => p

p%record_number = depth + 1
p%creator = trim(creator)
p%date = trim(date)
p%annotation= annotation

end subroutine ps_AddProvenanceRecord
!
subroutine ps_Delete_NonlocalProjectors(ps)
type(ps_t), intent(inout) :: ps

call destroy_nonlocal(ps%nonlocal)
end subroutine ps_Delete_NonlocalProjectors

subroutine ps_Delete_LocalPotential(ps)
type(ps_t), intent(inout) :: ps

call destroy_local(ps%local)
end subroutine ps_Delete_LocalPotential

end module m_psml_ps_edit
   



