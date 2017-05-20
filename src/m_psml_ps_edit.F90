!> Functions to edit the PSML ps_t structure
!!
!> @author Alberto Garcia
!
module m_psml_ps_edit

use m_psml_core   ! For basic structures
use external_interfaces, only: die => psml_die

use assoc_list, ps_annotation_t => assoc_list_t

implicit none

public :: ps_SetUUID
public :: ps_SetPSMLVersion
public :: ps_AddProvenanceRecord
public :: ps_Delete_LocalPotential
public :: ps_Delete_NonLocalProjectors
!public :: ps_AddSLBlock
private

CONTAINS !===============================================
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

allocate(p)

q => ps%provenance
if (associated(q)) then
   p%next => q
   q%prev => p
endif
ps%provenance => p

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
   



