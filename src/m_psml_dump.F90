!> Functions to handle PSML pseudopotential format structures
!!
!! (Dumpers)
!> @author Alberto Garcia
!
module m_psml_dump

use m_psml_core   ! For basic structures

use assoc_list, only: ps_annotation_t => assoc_list_t
!use assoc_list, only: EMPTY_ANNOTATION => EMPTY_ASSOC_LIST

use class_grid, Grid_t => Grid
use external_interfaces, only: die => psml_die

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

      call dump_annotation(xf,ps%annotation)

      call dump_provenance(xf,ps%provenance)

      call dump_header(xf,ps)

      if (initialized(ps%global_grid)) then
         call dump_grid(xf,ps%global_grid)
      endif

      call dump_semilocal_potentials(xf,ps)
      call dump_valence_charge(xf,ps%valence_charge,ps%global_grid)
      call dump_local_potential(xf,ps)
      call dump_nonlocal_projectors(xf,ps)
      call dump_pseudo_wavefunctions(xf,ps)

      if (trim(ps%header%core_corrections) == "yes") then
         call dump_core_charge(xf,ps%core_charge,ps%global_grid)
      endif
      
      call xml_EndElement(xf,"psml")
      call xml_Close(xf)

end subroutine ps_DumpToPSMLFile

subroutine dump_provenance(xf,p)
  use xmlf90_wxml
  use iso_varying_string, only: put, len, char
  
  type(xmlf_t), intent(inout) :: xf
  type(provenance_t), pointer :: p

  do while (associated(p))
     call xml_NewElement(xf,"provenance")
     call my_add_attribute(xf,"creator",trim(p%creator))
     call my_add_attribute(xf,"date",trim(p%date))
     if (len(p%input_file%buffer) > 0) then
        call xml_NewElement(xf,"input-file")
        call my_add_attribute(xf,"name",trim(p%input_file%name))
        call xml_AddCDataSection(xf,char(p%input_file%buffer), &
                                 line_feed=.true.)
        call xml_EndElement(xf,"input-file")
     endif
     call dump_annotation(xf,p%annotation)
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
  call dump_annotation(xf,p%annotation)
  call xml_EndElement(xf,"valence-configuration")
     
end subroutine dump_config_val


subroutine dump_header(xf,ps)
  use xmlf90_wxml
  type(xmlf_t), intent(inout) :: xf
  type(ps_t), intent(in), target :: ps

  type(header_t), pointer :: h

  h => ps%header
  call xml_NewElement(xf,"header")
  call my_add_attribute(xf,"atomic-label",trim(h%atomic_label))
  call my_add_attribute(xf,"z-pseudo",str(h%zpseudo))
  call my_add_attribute(xf,"atomic-number",str(h%z))
  call my_add_attribute(xf,"flavor",trim(h%flavor))
  call my_add_attribute(xf,"relativity",trim(h%relativity))
  if (h%polarized) then
     call my_add_attribute(xf,"polarized","yes")
  else
     call my_add_attribute(xf,"polarized","no")
  endif
  call my_add_attribute(xf,"core-corrections",trim(h%core_corrections))

  call dump_xc_info(xf,ps%xc_info)
  call dump_config_val(xf,ps%config_val)

  call xml_EndElement(xf,"header")

end subroutine dump_header

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
  call my_add_attribute(xf,"matching-radius",str(core%rcore))
  call my_add_attribute(xf,"number-of-continuous-derivatives",str(core%n_cont_derivs))
  call dump_annotation(xf,core%annotation)
  call dump_radfunc(xf,core%rho_core,parent_grid)
  call xml_EndElement(xf,"pseudocore-charge")
end subroutine dump_core_charge
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
        call my_add_attribute(xf,"flavor",slvp%flavor)
        if (set == SET_LJ) then
           call my_add_attribute(xf,"j",str(slvp%j))
        endif
        ! Group set was not specified
        if (set == SET_NULL) then
           call my_add_attribute(xf,"set",str_of_set(slvp%set))
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

  logical :: has_vlocal

  lop => ps%local
  
  has_vlocal = associated(lop%Vlocal%data)

  if (has_vlocal) then
     call xml_NewElement(xf,"local-potential")
     call my_add_attribute(xf,"type",lop%vlocal_type)
     ! No processing of grids here
     call dump_radfunc(xf,lop%Vlocal,ps%global_grid)
     call dump_radfunc(xf,lop%chlocal,ps%global_grid)
     call xml_EndElement(xf,"local-potential")
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
           call my_add_attribute(xf,"j",str(nlpp%j))
        endif
        call my_add_attribute(xf,"seq",str(nlpp%seq))
        call my_add_attribute(xf,"ekb",str(nlpp%ekb))
        call my_add_attribute(xf,"type",nlpp%type)
        ! Group set was not specified
        if (set == SET_NULL) then
           call my_add_attribute(xf,"set",str_of_set(nlpp%set))
        endif
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
  use sets_m

  type(xmlf_t), intent(inout) :: xf
  type(ps_t), intent(in), target :: ps

  type(pswfs_t), pointer :: wfp
  type(Grid_t)             :: parent_grid
  type(set_info_t) :: set_info

  integer :: i, j, set, n_grid
  integer, allocatable :: idx(:)

  wfp => ps%pswfs
  call sort_sets(wfp%npswfs,wfp%set,set_info)
  do j = 1, nsets(set_info)
     call set_indexes(set_info,j,idx)
     set = set_id(set_info,j)
     call xml_NewElement(xf,"pseudo-wave-functions")
     call my_add_attribute(xf,"set",str_of_set(set))
     !
     ! Check grids and decide whether to include a <grid> element
     ! First, check whether any wfn is not using the global grid.
     ! If all the wfns are in that case, assume that a mid-level grid was
     ! specified, dump a grid in a new <grid> element, and pass it 
     ! to the radfunc dumper. In the worst case scenario, we will
     ! have chosen a truly "radfunc-private" grid and there would be
     ! some replication of data. To avoid this, one would need to
     ! classify the grids. Maybe in a new version.
     !
     parent_grid = ps%global_grid
     n_grid = 0
     do i = 1, size(idx)
        if (.not. same(wfp%Phi(idx(i))%grid,ps%global_grid)) then
           n_grid = n_grid + 1
        endif
     enddo
     if (n_grid == size(idx)) then
        call dump_grid(xf,wfp%Phi(idx(1))%grid)
        parent_grid = wfp%Phi(idx(1))%grid
     endif
     do i = 1, size(idx)
        call xml_NewElement(xf,"pswf")
        call my_add_attribute(xf,"n",str(wfp%n(idx(i))))
        call my_add_attribute(xf,"l",wfp%l(idx(i)))
        if (set == SET_LJ) then
           call my_add_attribute(xf,"j",str(wfp%j(idx(i))))
        endif
        call dump_radfunc(xf,wfp%Phi(idx(i)),parent_grid)
        call xml_EndElement(xf,"pswf")
     enddo
     call xml_EndElement(xf,"pseudo-wave-functions")
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
  use assoc_list, only: ps_annotation_t => assoc_list_t
  use assoc_list, only: nitems_annotation => assoc_list_nitems
  use assoc_list, only: get_annotation_key => assoc_list_get_key
  use assoc_list, only: get_annotation_value => assoc_list_get_value

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
