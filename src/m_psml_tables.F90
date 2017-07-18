!> Functions to handle PSML pseudopotential format structures
!!
!! (Table generators)
!> @author Alberto Garcia
!
module m_psml_tables

use m_psml_core   ! For basic structures

implicit none

public :: ps_GenerateTables

private

CONTAINS !===============================================


!-------------------------------------
subroutine generate_table_sl(ps)

  type(ps_t), intent(inout), target :: ps

  type(semilocal_t), pointer :: slp
  type(slps_t), pointer :: slvp

  integer :: npots

  if (allocated(ps%sl_table)) then
     deallocate(ps%sl_table)
  endif

  npots = 0
  slp => ps%semilocal
  do while (associated(slp))
     slvp => slp%pot
     do while (associated(slvp))
        npots = npots + 1
        slvp => slvp%next
     enddo
     slp => slp%next
  enddo

  allocate(ps%sl_table(npots))
  
  npots = 0
  slp => ps%semilocal
  do while (associated(slp))
     slvp => slp%pot
     do while (associated(slvp))
        npots = npots + 1
        ps%sl_table(npots)%p => slvp
        slvp => slvp%next
     enddo
     slp => slp%next
  enddo

end subroutine generate_table_sl

!-------------------------------------
subroutine generate_table_nl(ps)

  type(ps_t), intent(inout), target :: ps

  type(nonlocal_t), pointer :: nlp
  type(nlpj_t), pointer :: nlpp

  integer :: nprojs

  if (allocated(ps%nl_table)) then
     deallocate(ps%nl_table)
  endif

  nprojs = 0
  nlp => ps%nonlocal
  do while (associated(nlp))
     nlpp => nlp%proj
     do while (associated(nlpp))
        nprojs = nprojs + 1
        nlpp => nlpp%next
     enddo
     nlp => nlp%next
  enddo

  allocate(ps%nl_table(nprojs))
  
  nprojs = 0
  nlp => ps%nonlocal
  do while (associated(nlp))
     nlpp => nlp%proj
     do while (associated(nlpp))
        nprojs = nprojs + 1
        ps%nl_table(nprojs)%p => nlpp
        nlpp => nlpp%next
     enddo
     nlp => nlp%next
  enddo

end subroutine generate_table_nl

!------------------------------------

subroutine generate_table_wf(ps)

  type(ps_t), intent(inout), target :: ps

  type(wfns_t), pointer :: wfp
  type(wf_t), pointer   :: wfpp

  integer :: nwfns

  if (allocated(ps%wf_table)) then
     deallocate(ps%wf_table)
  endif

  nwfns = 0
  wfp => ps%wavefunctions
  do while (associated(wfp))
     wfpp => wfp%wf
     do while (associated(wfpp))
        nwfns = nwfns + 1
        wfpp => wfpp%next
     enddo
     wfp => wfp%next
  enddo

  allocate(ps%wf_table(nwfns))
  
  nwfns = 0
  wfp => ps%wavefunctions
  do while (associated(wfp))
     wfpp => wfp%wf
     do while (associated(wfpp))
        nwfns = nwfns + 1
        ps%wf_table(nwfns)%p => wfpp
        wfpp => wfpp%next
     enddo
     wfp => wfp%next
  enddo

end subroutine generate_table_wf

subroutine ps_GenerateTables(ps)

  type(ps_t), intent(inout), target :: ps

  call generate_table_sl(ps)
  call generate_table_nl(ps)
  call generate_table_wf(ps)

end subroutine ps_GenerateTables

end module m_psml_tables

