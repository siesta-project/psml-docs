!> Functions to handle PSML pseudopotential format structures
!!
!! (Dumpers)
!> @author Alberto Garcia
!
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module m_psml_dump

use m_psml_core   ! For basic structures

use m_psml_assoc_list, only: ps_annotation_t => assoc_list_t
!use m_psml_assoc_list, only: EMPTY_ANNOTATION => EMPTY_ASSOC_LIST

use m_psml_class_grid, Grid_t => Grid
use m_psml_external_interfaces, only: die => psml_die

implicit none

integer, parameter    :: dp = selected_real_kind(14)
logical               :: global_debug = .false.
character(len=1), dimension(0:4) :: sym = (/ "s", "p", "d", "f", "g" /)

public :: ps_DumpToPSMLFile

private

CONTAINS !===============================================

subroutine ps_DumpToPSMLFile(ps,fname,indent)

  use xmlf90_wxml

  type(ps_t), intent(in) :: ps
  character(len=*), intent(in) :: fname
  logical, intent(in), optional :: indent

  type(xmlf_t)  :: xf

  call xml_OpenFile(trim(fname),xf, indent)

  call xml_AddXMLDeclaration(xf,"UTF-8")

  call xml_NewElement(xf,"psml")
  call my_add_attribute(xf,"version",trim(ps%version))
  call my_add_attribute(xf,"energy_unit",trim(ps%energy_unit))
  call my_add_attribute(xf,"length_unit",trim(ps%length_unit))
  call my_add_attribute(xf,"uuid",ps%uuid)
      call my_add_attribute(xf,"xmlns",ps%namespace)

! No top-level annotations in V1.1      
!      call dump_annotation(xf,ps%annotation)

  call dump_provenance(xf,ps%provenance)

      call dump_pseudo_atom_spec(xf,ps)

  if (initialized(ps%global_grid)) then
     call dump_grid(xf,ps%global_grid)
  endif

  call dump_valence_charge(xf,ps%valence_charge,ps%global_grid)
  if (trim(ps%header%core_corrections) == "yes") then
     call dump_core_charge(xf,ps%core_charge,ps%global_grid)
  endif

  if (trim(ps%header%meta_gga) == "yes") then
     call dump_valence_kinetic_energy_density(xf,ps%valence_kinetic_energy_density,ps%global_grid)
     if (trim(ps%header%core_corrections) == "yes") then
        call dump_core_kinetic_energy_density(xf,ps%core_kinetic_energy_density,ps%global_grid)
     endif
  endif

  call dump_semilocal_potentials(xf,ps)
  call dump_local_potential(xf,ps)
  call dump_nonlocal_projectors(xf,ps)
  call dump_pseudo_wavefunctions(xf,ps)

      
  call xml_EndElement(xf,"psml")
  call xml_Close(xf)

end subroutine ps_DumpToPSMLFile

subroutine dump_provenance(xf,p)

  use xmlf90_wxml
  
  type(xmlf_t), intent(inout) :: xf
  type(provenance_t), pointer :: p

  integer :: depth
  type(provenance_t), pointer :: q

  depth = 0
  q => p
  do while (associated(q))
     depth = depth + 1
     q => q%next
  enddo
  
  do while (associated(p))
     call xml_NewElement(xf,"provenance")
     call my_add_attribute(xf,"record-number",str(p%record_number))
     call my_add_attribute(xf,"creator",trim(p%creator))
     call my_add_attribute(xf,"date",trim(p%date))
     call dump_annotation(xf,p%annotation)
     if (len(p%input_file%buffer) > 0) then
        call xml_NewElement(xf,"input-file")
        call my_add_attribute(xf,"name",trim(p%input_file%name))
        call xml_AddCDataSection(xf,p%input_file%buffer, &
                                 line_feed=.true.)
        call xml_EndElement(xf,"input-file")
     endif
     call xml_EndElement(xf,"provenance")
     p => p%next
  end do

end subroutine dump_provenance


subroutine dump_xc_info(xf,p)

  use xmlf90_wxml

  type(xmlf_t), intent(inout) :: xf
  type(xc_t), intent(in) :: p

  integer :: i
  
  call xml_NewElement(xf,"exchange-correlation")
  call dump_annotation(xf,p%annotation)
  call xml_NewElement(xf,"libxc-info")
  call my_add_attribute(xf,"number-of-functionals",str(p%n_functs_libxc))
  do i = 1, p%n_functs_libxc
     call xml_NewElement(xf,"functional")
     call my_add_attribute(xf,"name",trim(p%libxc_name(i)))
     if (trim(p%libxc_type(i)) /= "UNKNOWN") then
        call my_add_attribute(xf,"type",trim(p%libxc_type(i)))
     endif
     call my_add_attribute(xf,"id",str(p%libxc_id(i)))
     if (p%libxc_weight(i) /= 1.0_dp) then
        call my_add_attribute(xf,"weight",str(p%libxc_weight(i)))
     endif
     call xml_EndElement(xf,"functional")
  enddo
  call xml_EndElement(xf,"libxc-info")
  call xml_EndElement(xf,"exchange-correlation")
     
end subroutine dump_xc_info

subroutine dump_config_val(xf,p)

  use xmlf90_wxml

  type(xmlf_t), intent(inout) :: xf
  type(config_val_t), intent(in) :: p

  integer :: i

  call xml_NewElement(xf,"valence-configuration")
  call my_add_attribute(xf,"total-valence-charge",str(p%total_charge))
  call dump_annotation(xf,p%annotation)
  do i = 1, p%nshells
     call xml_NewElement(xf,"shell")
     call my_add_attribute(xf,"n",str(p%n(i)))
     call my_add_attribute(xf,"l",p%l(i))
     call my_add_attribute(xf,"occupation",str(p%occ(i)))
     if ((p%occ_up(i) + p%occ_down(i)) /= 0.0_dp) then
        call my_add_attribute(xf,"occupation-down",str(p%occ_down(i)))
        call my_add_attribute(xf,"occupation-up",str(p%occ_up(i)))
     endif
     call xml_EndElement(xf,"shell")
  enddo
  call xml_EndElement(xf,"valence-configuration")
     
end subroutine dump_config_val


subroutine dump_pseudo_atom_spec(xf,ps)

  use xmlf90_wxml

  type(xmlf_t), intent(inout) :: xf
  type(ps_t), intent(in), target :: ps

  type(header_t), pointer :: h

  h => ps%header
  call xml_NewElement(xf,"pseudo-atom-spec")
  call my_add_attribute(xf,"atomic-label",trim(h%atomic_label))
  call my_add_attribute(xf,"atomic-number",str(h%z))
  call my_add_attribute(xf,"z-pseudo",str(h%zpseudo))
  call my_add_attribute(xf,"flavor",trim(h%flavor))
  call my_add_attribute(xf,"relativity",trim(h%relativity))
  if (h%polarized) then
     call my_add_attribute(xf,"spin-dft","yes")
  else
     call my_add_attribute(xf,"spin-dft","no")
  endif
  call my_add_attribute(xf,"core-corrections",trim(h%core_corrections))
  call my_add_attribute(xf,"meta-gga",trim(h%meta_gga))

  call dump_annotation(xf,h%annotation)

  call dump_xc_info(xf,ps%xc_info)
  call dump_config_val(xf,ps%config_val)

  call xml_EndElement(xf,"pseudo-atom-spec")

end subroutine dump_pseudo_atom_spec

subroutine dump_radfunc(xf,rf,parent_grid)

  use xmlf90_wxml

  type(xmlf_t), intent(inout) :: xf
  type(radfunc_t), intent(in) :: rf
  type(Grid_t)                  :: parent_grid ! Only one level for now

  if (.not. initialized(rf%grid)) return

  call xml_NewElement(xf,"radfunc")
  if (same(rf%grid,parent_grid)) then
     ! do nothing
  else
     call dump_grid(xf,rf%grid)
  endif
  call xml_NewElement(xf,"data")
  ! Cover the case in which the data uses only an
  ! initial section of the grid
  if (size(rf%data) < sizeGrid(rf%grid)) then
     call my_add_attribute(xf,"npts",str(size(rf%data)))
  endif
  call xml_AddArray(xf,rf%data(:))
  call xml_EndElement(xf,"data")
  call xml_EndElement(xf,"radfunc")
end subroutine dump_radfunc

subroutine dump_valence_charge(xf,val,parent_grid)

  use xmlf90_wxml

  type(xmlf_t), intent(inout) :: xf
  type(valence_charge_t), intent(in) :: val
  type(Grid_t)           :: parent_grid

  call xml_NewElement(xf,"valence-charge")
  call my_add_attribute(xf,"total-charge",str(val%total_charge))
  call dump_annotation(xf,val%annotation)
  call dump_radfunc(xf,val%rho_val,parent_grid)
  call xml_EndElement(xf,"valence-charge")
end subroutine dump_valence_charge

subroutine dump_core_charge(xf,core,parent_grid)

  use xmlf90_wxml

  type(xmlf_t), intent(inout) :: xf
  type(core_charge_t), intent(in) :: core
  type(Grid_t)                 :: parent_grid

  call xml_NewElement(xf,"pseudocore-charge")
  if (core%rcore >= 0.0_dp) then
     call my_add_attribute(xf,"matching-radius",str(core%rcore))
  endif
  if (core%n_cont_derivs >= 0 ) then
     call my_add_attribute(xf,"number-of-continuous-derivatives",str(core%n_cont_derivs))
  endif
  call dump_annotation(xf,core%annotation)
  call dump_radfunc(xf,core%rho_core,parent_grid)
  call xml_EndElement(xf,"pseudocore-charge")
end subroutine dump_core_charge
!

subroutine dump_valence_kinetic_energy_density(xf,val,parent_grid)

  use xmlf90_wxml

  type(xmlf_t), intent(inout) :: xf
  type(valence_kinetic_energy_density_t), intent(in) :: val
  type(Grid_t)           :: parent_grid

  call xml_NewElement(xf,"valence-kinetic-energy-density")
  call my_add_attribute(xf,"is-unscreening-tau",trim(val%is_unscreening_tau))
  call dump_annotation(xf,val%annotation)
  call dump_radfunc(xf,val%kin_edens_val,parent_grid)
  call xml_EndElement(xf,"valence-kinetic-energy-density")
end subroutine dump_valence_kinetic_energy_density

subroutine dump_core_kinetic_energy_density(xf,core,parent_grid)

  use xmlf90_wxml

  type(xmlf_t), intent(inout) :: xf
  type(core_kinetic_energy_density_t), intent(in) :: core
  type(Grid_t)                 :: parent_grid

  call xml_NewElement(xf,"pseudocore-kinetic-energy-density")
  if (core%rcore >= 0.0_dp) then
     call my_add_attribute(xf,"matching-radius",str(core%rcore))
  endif
  if (core%n_cont_derivs >= 0 ) then
     call my_add_attribute(xf,"number-of-continuous-derivatives",str(core%n_cont_derivs))
  endif
  call dump_annotation(xf,core%annotation)
  call dump_radfunc(xf,core%kin_edens_core,parent_grid)
  call xml_EndElement(xf,"pseudocore-charge")
end subroutine dump_core_kinetic_energy_density
!

subroutine dump_semilocal_potentials(xf,ps)

  use xmlf90_wxml

  type(xmlf_t), intent(inout) :: xf
  type(ps_t), intent(in), target :: ps

  type(semilocal_t), pointer :: slp
  type(slps_t), pointer :: slvp
  type(Grid_t)            :: parent_grid
  integer :: i, j, set

  slp => ps%semilocal
  do while (associated(slp))
     set = slp%set
     call xml_NewElement(xf,"semilocal-potentials")
     if (set /= SET_NULL) then
        ! Group set was specified
        call my_add_attribute(xf,"set",str_of_set(set))
     else
        call die("Set not specified in SemiLocalPotentials block")
     endif
     call dump_annotation(xf,slp%annotation)
     if (initialized(slp%grid)) then
        parent_grid = slp%grid
        call dump_grid(xf,slp%grid)
     else
        parent_grid = ps%global_grid
     endif
     slvp => slp%pot
     do while (associated(slvp))
        call xml_NewElement(xf,"slps")
        call my_add_attribute(xf,"n",str(slvp%n))
        call my_add_attribute(xf,"l",slvp%l)
        call my_add_attribute(xf,"rc",str(slvp%rc))
        ! If eref has a physical value, output it
        if (slvp%eref < 0.1*huge(1.0_dp)) then
           call my_add_attribute(xf,"eref",str(slvp%eref))
        endif
        call my_add_attribute(xf,"flavor",slvp%flavor)
        if (set == SET_LJ) then
           call my_add_attribute(xf,"j",str(slvp%j,format="(f3.1)"))
        endif
        call dump_radfunc(xf,slvp%V,parent_grid)
        call xml_EndElement(xf,"slps")
        slvp => slvp%next
     enddo
     call xml_EndElement(xf,"semilocal-potentials")
     slp => slp%next
  enddo
  call delete(parent_grid)
end subroutine dump_semilocal_potentials

subroutine dump_local_potential(xf,ps)

  use xmlf90_wxml

  type(xmlf_t), intent(inout) :: xf
  type(ps_t), intent(in), target :: ps

  type(local_t), pointer :: lop
  type(Grid_t)           :: parent_grid

  logical :: has_vlocal, has_local_charge

  lop => ps%local
  
  has_vlocal = associated(lop%Vlocal%data)
  has_local_charge = associated(lop%Chlocal%data)

  if (has_vlocal) then
     call xml_NewElement(xf,"local-potential")
     call my_add_attribute(xf,"type",lop%vlocal_type)
     call dump_annotation(xf,lop%annotation)
     ! No processing of grids here
     if (initialized(lop%grid)) then
        parent_grid = lop%grid
        call dump_grid(xf,lop%grid)
     else
        parent_grid = ps%global_grid
     endif
     call dump_radfunc(xf,lop%Vlocal,parent_grid)
     if (has_local_charge) then
        call xml_NewElement(xf,"local-charge")
        call dump_radfunc(xf,lop%chlocal,parent_grid)
        call xml_EndElement(xf,"local-charge")
     endif
     call xml_EndElement(xf,"local-potential")
     call delete(parent_grid)
  endif
end subroutine dump_local_potential

!----------------------------------------------------------
subroutine dump_nonlocal_projectors(xf,ps)

  use xmlf90_wxml

  type(xmlf_t), intent(inout) :: xf
  type(ps_t), intent(in), target :: ps

  type(nonlocal_t), pointer :: nlp
  type(nlpj_t), pointer :: nlpp
  type(Grid_t)            :: parent_grid
  integer :: set

  nlp => ps%nonlocal
  do while (associated(nlp))
     set = nlp%set
     call xml_NewElement(xf,"nonlocal-projectors")
     if (set /= SET_NULL) then
        ! Group set was specified
        call my_add_attribute(xf,"set",str_of_set(set))
     else
        call die("Set not specified in NonLocalProjectors block")
     endif
     call dump_annotation(xf,nlp%annotation)
     if (initialized(nlp%grid)) then
        parent_grid = nlp%grid
        call dump_grid(xf,nlp%grid)
     else
        parent_grid = ps%global_grid
     endif
     nlpp => nlp%proj
     do while (associated(nlpp))
        call xml_NewElement(xf,"proj")
        call my_add_attribute(xf,"l",nlpp%l)
        if (set == SET_LJ) then
           call my_add_attribute(xf,"j",str(nlpp%j,format="(f3.1)"))
        endif
        call my_add_attribute(xf,"seq",str(nlpp%seq))
        call my_add_attribute(xf,"ekb",str(nlpp%ekb))
        ! If eref has a physical value, output it
        if (nlpp%eref < 0.1*huge(1.0_dp)) then
           call my_add_attribute(xf,"eref",str(nlpp%eref))
        endif
        call my_add_attribute(xf,"type",nlpp%type)
        call dump_radfunc(xf,nlpp%proj,parent_grid)
        call xml_EndElement(xf,"proj")
        nlpp => nlpp%next
     enddo
     call xml_EndElement(xf,"nonlocal-projectors")
     nlp => nlp%next
  enddo
  call delete(parent_grid)
end subroutine dump_nonlocal_projectors
!----------------------------------------------------------
subroutine dump_pseudo_wavefunctions(xf,ps)

  use xmlf90_wxml

  type(xmlf_t), intent(inout) :: xf
  type(ps_t), intent(in), target :: ps

  type(wfns_t), pointer :: wfp
  type(wf_t), pointer   :: wfpp
  type(Grid_t)             :: parent_grid

  integer :: set
  
  wfp => ps%wavefunctions
  do while (associated(wfp))
     set = wfp%set
     call xml_NewElement(xf,"pseudo-wave-functions")
     if (set /= SET_NULL) then
        ! Group set was specified
        call my_add_attribute(xf,"set",str_of_set(set))
     else
        call die("Set not specified in Wavefunctions block")
     endif
!     if (wfp%type /= "") then
!        call my_add_attribute(xf,"type",wfp%type)
!     endif
     call dump_annotation(xf,wfp%annotation)
     if (initialized(wfp%grid)) then
        parent_grid = wfp%grid
        call dump_grid(xf,wfp%grid)
     else
        parent_grid = ps%global_grid
     endif
     wfpp => wfp%wf
     do while (associated(wfpp))
        call xml_NewElement(xf,"pswf")
        call my_add_attribute(xf,"n",str(wfpp%n))
        call my_add_attribute(xf,"l",wfpp%l)
        if (set == SET_LJ) then
           call my_add_attribute(xf,"j",str(wfpp%j,format="(f3.1)"))
        endif
        ! If energy_level has a physical value, output it
        if (wfpp%energy_level < 0.1*huge(1.0_dp)) then
           call my_add_attribute(xf,"energy_level",str(wfpp%energy_level))
        endif

        call dump_radfunc(xf,wfpp%Phi,parent_grid)
        call xml_EndElement(xf,"pswf")
        wfpp => wfpp%next
     enddo
     call xml_EndElement(xf,"pseudo-wave-functions")
     wfp => wfp%next
  enddo
  call delete(parent_grid)

end subroutine dump_pseudo_wavefunctions

subroutine dump_grid(xf,agrid)

  use xmlf90_wxml

  type(xmlf_t), intent(inout) :: xf
  type(Grid_t), intent(in) :: agrid

  if (.not. initialized(agrid)) return
  
  call xml_NewElement(xf,"grid")
  call my_add_attribute(xf,"npts",str(sizeGrid(agrid)))
  call dump_annotation(xf,annotationGrid(agrid))
  call xml_NewElement(xf,"grid-data")
  call xml_AddArray(xf,valGrid(agrid))
  call xml_EndElement(xf,"grid-data")
  call xml_EndElement(xf,"grid")
end subroutine dump_grid

subroutine dump_annotation(xf,annotation)

  use xmlf90_wxml

  use m_psml_assoc_list, only: ps_annotation_t => assoc_list_t
  use m_psml_assoc_list, only: nitems_annotation => assoc_list_nitems
  use m_psml_assoc_list, only: get_annotation_key => assoc_list_get_key
  use m_psml_assoc_list, only: get_annotation_value => assoc_list_get_value

  type(xmlf_t), intent(inout) :: xf
  type(ps_annotation_t), intent(in)  :: annotation

integer :: n_items, i, stat
character(len=256) :: key, val

n_items = nitems_annotation(annotation)
if (n_items > 0) then

   call xml_NewElement(xf,"annotation")
   do i = 1, n_items
      call get_annotation_key(annotation,i,key,stat)
      call get_annotation_value(annotation,i,val,stat)
      call my_add_attribute(xf,trim(key),trim(val))
   enddo
   call xml_EndElement(xf,"annotation")
endif
end subroutine dump_annotation

subroutine my_add_attribute(xf,name,value)

  use xmlf90_wxml

  type(xmlf_t), intent(inout)   :: xf
  character(len=*), intent(in)  :: name
  character(len=*), intent(in)  :: value

  call xml_AddAttribute(xf,name,trim(value))
end subroutine my_add_attribute

end module m_psml_dump

