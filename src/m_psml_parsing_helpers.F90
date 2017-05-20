module m_psml_parsing_helpers
!
!  This module reads a pseudopotential file written in the PSML format
!  A full example of the building up of a data structure using
!  the SAX paradigm.
!
 use m_psml_core              ! For data types and basic utilities
 use external_interfaces, only: die => psml_die
 use class_Grid
 use assoc_list, only: ps_annotation_t => assoc_list_t

 
implicit none

private

!
! It defines the routines that are called from xml_parser in response
! to particular events.
!
public  :: begin_element, end_element, pcdata_chunk
public  :: cdata_section_chunk

!
! The data will be stored in this public variable
! There are some design issues to decide:
! -- Should this be a pointer, associated by the client
!    program to its own variable? In that case, the
!    client should make sure that the variable is "clean"
!    before calling this routine, as some fields will be
!    allocated here.
! -- Perhaps it should be a pointer allocated here (and
!    then destroyed when done by the client). It should be
!    allocated at the beginning of processing, maybe detected
!    with a (default) "begin_Document" handler, or by 
!    "begin_Element" after checking for association.
!    This is the cleanest option, as the caller might want
!    to keep several instances alive at the same time...
!    (... but this should be handled by the user)
! -- If "pseudo" here is a normal variable, it should also
!    be "cleaned" before the next use. The current usage
!    in Abinit falls in this category: psxml is a pointer
!    associated to "pseudo", and cleaned after use.
!
!    We implement the first option now

type(ps_t), pointer, public, save :: pseudo => null() 
logical, public, save             :: debug_parsing = .false.

logical, private, save  :: in_psml = .false.
logical, private, save  :: in_slps = .false. , in_radfunc = .false.
logical, private, save  :: in_semilocal = .false. , in_header = .false.
logical, private, save  :: in_coreCharge = .false. , in_data = .false.
logical, private, save  :: in_grid_data = .false. , in_grid = .false.
logical, private, save  :: in_valenceCharge = .false.
logical, private, save  :: in_provenance = .false., in_input_file = .false.
logical, private, save  :: in_valence_config = .false.
logical, private, save  :: in_xc = .false., in_libxc_info = .false.
logical, private, save  :: in_pseudowavefun = .false. , in_pswf = .false.
logical, private, save  :: in_chlocal = .false., in_nonlocal = .false.
logical, private, save  :: in_proj = .false. , in_local_potential = .false.
logical, private, save  :: got_explicit_grid_data

integer, private, save  :: ndata, ndata_grid
integer, private, save  :: n_funct

character(len=20), private, save :: current_wf_set
character(len=20), private, save :: current_sl_set
character(len=20), private, save :: current_proj_set
character(len=40), private, save :: top_flavor

integer, parameter, private    :: dp = selected_real_kind(14)
real(dp), private, save        :: zval_generation

type(Grid), private, save                  :: tmp_grid
real(dp), private, save, pointer           :: gdata(:) => null()
type(ps_annotation_t), private, save, pointer :: gannot => null()
!
! Pointers to make it easier to manage the data
!
type(provenance_t), private, pointer      :: pp => null()
type(provenance_t), private, pointer      :: qp => null()
type(input_file_t), private, pointer      :: ifp => null()
type(header_t), private, pointer          :: hp => null()
type(config_val_t), private, pointer      :: cp => null()
type(xc_t), private, pointer              :: xp => null()
type(pswfs_t), private, pointer           :: wfp => null()
type(semilocal_t), private, pointer       :: slp => null()
type(semilocal_t), private, pointer       :: qslp => null()
type(slps_t), private, pointer            :: slvp => null()
type(slps_t), private, pointer            :: qslvp => null()
type(local_t), private, pointer           :: lop => null()
type(nonlocal_t), private, pointer        :: nlp => null()
type(nonlocal_t), private, pointer        :: qnlp => null()
type(nlpj_t), private, pointer            :: nlpp => null()
type(nlpj_t), private, pointer            :: qnlpp => null()
type(valence_charge_t), private, pointer  :: valp => null()
type(core_charge_t), private, pointer     :: corep => null()
type(radfunc_t), private, pointer         :: rp => null()

character(len=100), private, save :: parent_element=""

CONTAINS  !===========================================================

!----------------------------------------------------------------------
#ifndef PSML_USE_FOX
subroutine begin_element(name,attributes)
use xmlf90_sax, only: dictionary_t, get_value
#else
subroutine begin_element(namespaceURI,localName,name,attributes)
use Fox_sax, only: dictionary_t
use fox_extra, only: get_value=>get_value_by_key
character(len=*), intent(in)    :: namespaceURI
character(len=*), intent(in)    :: localName
#endif
character(len=*), intent(in)    :: name
type(dictionary_t), intent(in)  :: attributes

character(len=100)  :: value, msg
real                :: version_number
integer             :: status

integer             :: i, npts

if (debug_parsing) print *, "Element: ", trim(name)

select case(name)

      case ("psml")

         in_psml = .true.
         ! Make sure that pseudo is pointing to something

         if (.not. associated(pseudo)) then
            call die("ps_t object not initialized by client")
         endif

         call get_value(attributes,"version",value,status)
         if (status /= 0 ) call die("No psml version")
         read(value,fmt=*) version_number
         if ( (version_number < PSML_TARGET_VERSION_LO) .or. &
              (version_number > PSML_TARGET_VERSION_HI)) then
            write(msg,"('[',f4.2,',',f4.2,']')")  &
                       PSML_TARGET_VERSION_LO, &
                       PSML_TARGET_VERSION_HI
            call die("This version of the library can "// &
                     "process PSML files with version in "//trim(msg))
         endif
         pseudo%version = value
         call get_value(attributes,"energy_unit",pseudo%energy_unit,status)
         if (status /= 0 ) call die("No energy unit")
         call get_value(attributes,"length_unit",pseudo%length_unit,status)
         if (status /= 0 ) call die("No length unit")

         call get_value(attributes,"uuid",pseudo%uuid,status)
         if (status /= 0 ) pseudo%uuid = "no-uuid-specified"

         ! Initialize counters

         pseudo%pswfs%npswfs = 0

      case ("provenance")
         in_provenance = .true.
         ! This will gather all element provenances at the top level
            allocate(pp)
            if (associated(pseudo%provenance)) then
               qp => pseudo%provenance
               do while (associated(qp%next))
                  qp => qp%next
               enddo
               qp%next => pp
            else
               pseudo%provenance => pp
            endif

            call get_value(attributes,"creator",pp%creator,status)
            if (status /= 0 ) pp%creator="unknown"
 
            call get_value(attributes,"date",pp%date,status)
            if (status /= 0 ) pp%date="unknown"

      case ("input-file")
         if (.not. in_provenance) call die("<input-file> outside <provenance>")
         in_input_file = .true.
         ifp => pp%input_file
         call get_value(attributes,"name",ifp%name,status)
         if (status /= 0 ) ifp%name="unknown"
         
      case ("header")
         in_header = .true.
         hp => pseudo%header
         
         call get_value(attributes,"atomic-label",hp%atomic_label,status)
         if (status /= 0 ) call die("Cannot determine atomic-label")

         call get_value(attributes,"z-pseudo",value,status)
         if (status /= 0 ) call die("Cannot determine z-pseudo")
         read(unit=value,fmt=*) hp%zpseudo

         call get_value(attributes,"atomic-number",value,status)
         if (status /= 0 ) call die("Cannot determine atomic number")
         read(unit=value,fmt=*) hp%z

         call get_value(attributes,"flavor",hp%flavor,status)
         if (status /= 0 ) hp%flavor="not-unique"

         call get_value(attributes,"relativity",hp%relativity,status)
         if (status /= 0 ) call die("Cannot determine relativity scheme")

         call get_value(attributes,"polarized",value,status)
         if (status /= 0 ) value = "no"
         hp%polarized = (value == "yes")
         if (hp%polarized .and. trim(hp%relativity)=="dirac") then
            call die("Cannot be polarized and fully relativistic at the same time")
         endif

         call get_value(attributes,"core-corrections", &
                                    hp%core_corrections,status)
         if (status /= 0 ) hp%core_corrections = "no"

      case ("exchange-correlation")
         in_xc = .true.
         xp => pseudo%xc_info

      case ("libxc-info")
         if (.not. in_xc) call die("Orphan <libxc-info>")
         in_libxc_info = .true.
         call get_value(attributes,"number-of-functionals", &
                                    value,status)
         if (status /= 0 ) call die("Error reading number of libxc functs")
         read(unit=value,fmt=*)  xp%n_functs_libxc 
         n_funct = xp%n_functs_libxc
         allocate(xp%libxc_name(n_funct), xp%libxc_id(n_funct))
         allocate(xp%libxc_weight(n_funct),xp%libxc_type(n_funct))
         xp%libxc_weight(1:n_funct) = 1.0_dp
         xp%libxc_type(1:n_funct) = "UNKNOWN"

         n_funct = 0   ! for checking the counting on the fly

      case ("functional")
         if (.not. in_libxc_info) call die("Orphan <functional>")
         n_funct = n_funct + 1
         if (n_funct > xp%n_functs_libxc) &
           call die("Too many <functional> elements in <libxc-info>")

         call get_value(attributes,"name", &
                                    xp%libxc_name(n_funct),status)
         if (status /= 0 ) call die("Error reading libxc name")

         call get_value(attributes,"id", value, status)
         if (status /= 0 ) call die("Error reading libxc id")
         read(unit=value,fmt=*)  xp%libxc_id(n_funct) 

         ! optional attribute(s)
         call get_value(attributes,"weight", value, status)
         if (status == 0 ) then
            read(unit=value,fmt=*)  xp%libxc_weight(n_funct) 
         endif
         call get_value(attributes,"type", value, status)
         if (status == 0 ) then
            xp%libxc_type(n_funct)  = trim(value)
         endif

      case ("valence-configuration")
         in_valence_config = .true.

         pseudo%config_val%nshells = 0
         pseudo%config_val%occ_up(:) = 0.0_dp
         pseudo%config_val%occ_down(:) = 0.0_dp
         pseudo%config_val%occ(:) = 0.0_dp
         call get_value(attributes,"total-valence-charge",value,status)
         if (status /= 0 ) call die("Cannot determine total-valence-charge")
         read(unit=value,fmt=*) pseudo%config_val%total_charge

      case ("shell")

         if (in_valence_config) then
            cp => pseudo%config_val
!         else if (in_core_config) then
!            cp => pseudo%config_core
         else
            call die("Orphan <shell> element")
         endif

         cp%nshells = cp%nshells + 1
         call get_value(attributes,"l",cp%l(cp%nshells),status)
         if (status /= 0 ) call die("Cannot determine l for shell")

         call get_value(attributes,"n",value,status)
         if (status /= 0 ) call die("Cannot determine n for shell")
         read(unit=value,fmt=*) cp%n(cp%nshells)

         call get_value(attributes,"occupation",value,status)
         if (status /= 0 ) call die("Cannot determine occupation for shell")
         read(unit=value,fmt=*) cp%occ(cp%nshells)

         call get_value(attributes,"occupation-up",value,status)
         if (status == 0 ) then
            read(unit=value,fmt=*) cp%occ_up(cp%nshells)
         endif
         call get_value(attributes,"occupation-down",value,status)
         if (status == 0 ) then
            read(unit=value,fmt=*) cp%occ_down(cp%nshells)
         endif

      case ("slps")
         in_slps = .true.
         if (.not. in_semilocal) call die("Orphan <slps> element")
         allocate(slvp)

         ! Append to end of list  !! call append(slp%pot,slvp)
         if (associated(slp%pot)) then
            qslvp => slp%pot
            do while (associated(qslvp%next))
               qslvp => qslvp%next
            enddo
            qslvp%next => slvp
         else
            !First link
            slp%pot => slvp
         endif

         rp => slvp%V

         slvp%parent_group => slp   ! current semilocal-potentials element
         
         call get_value(attributes,"set",value,status)
         if (status /= 0 ) then
            value = current_sl_set
         endif
         slvp%set = setcode_of_string(value)

         call get_value(attributes,"l",slvp%l,status)
         if (status /= 0 ) call die("Cannot determine l for SL potential")

         call get_value(attributes,"n",value,status)
         if (status /= 0 ) call die("Cannot determine n for SL potential")
         read(unit=value,fmt=*) slvp%n

         call get_value(attributes,"rc",value,status)
         if (status /= 0 ) call die("Cannot determine rc for SL potential")
         read(unit=value,fmt=*) slvp%rc

         call get_value(attributes,"j",value,status)
         if (status /= 0 ) then
            if (slvp%set == SET_LJ) &
                 call die("Cannot determine j for SLPS in set " // str_of_set(SET_LJ))
         else
            read(unit=value,fmt=*) slvp%j
         endif

         call get_value(attributes,"flavor",slvp%flavor,status)
         if (status /= 0 ) then
            slvp%flavor = top_flavor
         endif

      case ("proj")
         in_proj = .true.
         if (.not. in_nonlocal) call die("Orphan <proj> element")
         allocate(nlpp)
         
         ! Append to end of list  !! call append(nlp%proj,nlpp)
         if (associated(nlp%proj)) then
            qnlpp => nlp%proj
            do while (associated(qnlpp%next))
               qnlpp => qnlpp%next
            enddo
            qnlpp%next => nlpp
         else
            !First link
            nlp%proj => nlpp
         endif

         rp => nlpp%proj
         nlpp%parent_group => nlp   ! current nonlocal-projectors element

         call get_value(attributes,"set",value,status)
         if (status /= 0 ) then
            value = current_proj_set
         endif
         nlpp%set = setcode_of_string(value)

         call get_value(attributes,"l",nlpp%l,status)
         if (status /= 0 ) call die("Cannot determine l for proj")

         call get_value(attributes,"seq",value,status)
         if (status /= 0 ) call die("Cannot determine seq number for proj")
         read(unit=value,fmt=*) nlpp%seq

         call get_value(attributes,"ekb",value,status)
         if (status /= 0 ) call die("Cannot determine Ekb for proj")
         read(unit=value,fmt=*) nlpp%ekb

         call get_value(attributes,"j",value,status)
         if (status /= 0 ) then
            if (nlpp%set == SET_LJ) &
                 call die("Cannot determine j for Proj in set " // str_of_set(SET_LJ))
         else
            read(unit=value,fmt=*) nlpp%j
         endif

         call get_value(attributes,"type",nlpp%type,status)
         if (status /= 0 ) call die("Cannot determine type of proj")

      case ("pswf")

         if (.not. in_pseudowavefun) call die("Orphan <pswf> element")
         in_pswf = .true.

         wfp => pseudo%pswfs
         wfp%npswfs = wfp%npswfs + 1
         i = wfp%npswfs
         rp => wfp%Phi(i)

         call get_value(attributes,"set",value,status)
         if (status /= 0 ) then
            value = current_wf_set
         endif
         wfp%set(i) = setcode_of_string(value)

         call get_value(attributes,"l",wfp%l(i),status)
         if (status /= 0 ) call die("Cannot determine l for PSwf")
                                                                              
         call get_value(attributes,"n",value,status)
         if (status /= 0 ) call die("Cannot determine n for PSwf")
         read(unit=value,fmt=*) wfp%n(i)

         call get_value(attributes,"j",value,status)
         if (status /= 0 ) then
            if (wfp%set(i) == SET_LJ) &
                 call die("Cannot determine j for PSwf in set " // str_of_set(SET_LJ))
         else
            read(unit=value,fmt=*) wfp%j(i)
         endif

      case ("grid")
         in_grid = .true.

         got_explicit_grid_data = .false.

         ! This attribute is mandatory
         call get_value(attributes,"npts",value,status)
         if (status /= 0 ) call die("Cannot determine grid npts")
         read(unit=value,fmt=*) npts
         if (npts == 0) call die("Grid size not specified correctly")

         ! Create working object and associate inner sections
         ! while the parsing is active
         call newGrid(tmp_grid,npts)
         gdata => valGrid(tmp_grid)
         gannot => annotationGrid(tmp_grid)
         !
         ! In this way we allow for a private grid for each radfunc,
         ! or for a global grid specification
         !
         if (in_radfunc) then
            if (debug_parsing) print *, "Found grid in radfunc"
            if (initialized(rp%grid)) then
               call die("psml: Two grids specified for a radfunc")
            endif
            rp%grid = tmp_grid

         ! We check whether we are at the top level,
         ! or at an intermediate grouping level that allows a grid

         else if (in_nonlocal) then

            if (initialized(nlp%grid)) then
                call die("psml: Two grids in same nonlocal block")
            endif
            if (debug_parsing) print *, "Found nonlocal grid"
            nlp%grid = tmp_grid

         else if (in_local_potential) then

            if (initialized(lop%grid)) then
                call die("psml: Two grids in same local block")
            endif
            if (debug_parsing) print *, "Found local grid"
            lop%grid = tmp_grid

         else if (in_semilocal) then

            if (initialized(slp%grid)) then
                call die("psml: Two grids in same semilocal block")
            endif
            if (debug_parsing) print *, "Found semilocal grid"
            slp%grid = tmp_grid

         else if (in_pseudowavefun) then

            if (initialized(wfp%grid)) then
               !call die("psml: Two pseudo-wavefunction grids specified")
            endif
            if (debug_parsing) print *, "Found pseudo-wavefunction grid"
            wfp%grid = tmp_grid

         else  ! We are at the top level

            if (debug_parsing) print *, "Found grid at the top level"
            if (initialized(pseudo%global_grid)) then
               ! Maybe allow this in the future
               call die("psml: Two global grids specified")
            endif
            pseudo%global_grid = tmp_grid
         endif

      case ("data")
         if (.not. in_radfunc) then
            call die("<data> element outside <rad_func> element")
         endif
         in_data = .true.

	 ! The following blocks are a bit more verbose than needed since
         ! the Intel compiler seems to be trying to evaluate all the
         ! clauses joined by an .and. operator, instead of stopping if
         ! the first clause is .false.
       
         if (.not. initialized(rp%grid)) then

            ! Try regional grids first

            if (in_nonlocal) then
               if (initialized(nlp%grid)) then
                  rp%grid = nlp%grid
                  if (debug_parsing) print *,"Associated proj grid with nl parent grid"
               endif
            else if (in_local_potential) then
               if (initialized(lop%grid)) then
                  rp%grid = lop%grid
                  if (debug_parsing) print *,"Associated grid with vlocal parent grid"
               endif
            else if (in_semilocal) then
               if (initialized(slp%grid)) then
                  rp%grid = slp%grid
                  if (debug_parsing) print *,"Associated slps grid with sl parent grid"
               endif
            else if (in_pseudowavefun) then
               if (initialized(wfp%grid)) then
                  rp%grid = wfp%grid
               endif
            endif

         endif

         ! If the parent block does not include a grid, try the global grid

         if (.not. initialized(rp%grid)) then
               if (initialized(pseudo%global_grid)) then
                  rp%grid = pseudo%global_grid
                  if (debug_parsing) print *, "Associated grid with global grid"
               endif
         endif

         ! Now give up
         if (.not. initialized(rp%grid)) call die("Cannot find grid data for radfunc")

         allocate(rp%data(sizeGrid(rp%grid)))
         ndata = 0             ! To start the build up

      case ("grid-data")
         if (.not. in_grid) call die("Grid_data element outside grid element")
         in_grid_data = .true.
         got_explicit_grid_data = .true.
         if (size(gdata) == 0) call die("Grid npts attribute faulty")
         ndata_grid = 0             ! To start the build up

      case ("radfunc")

         ! We need to make sure that a radfunc is allowed at this level
         ! 
         ! For example, if an old-style file with <vps> is used, the finding
         ! of a <vps> element will not increase npots, and there will be a
         ! segmentation fault when trying to store the data in the rp pointer,
         ! which would be non-associated if the vps section is the first
         ! case of radfuncs in the file
         !
         ! Actually, it gets worse: if there is a <core-charge> element
         ! before the old-style <vps> section, the rp pointer used for 
         ! core-charge will be reused and assigned the data for the first
         ! and subsequent vps elements.
         !
         if (     in_slps .or. in_coreCharge .or. in_valenceCharge &
              .or. in_pswf .or. in_proj .or. in_local_potential  &
              .or. in_chlocal) then           
            in_radfunc = .true.
         else
            call die("<radfunc> element found under unallowed <"// &
                     trim(parent_element) // ">")
         endif

      case ("pseudocore-charge")
         in_coreCharge = .true.
         corep => pseudo%core_charge
         rp => corep%rho_core

         call get_value(attributes,"matching-radius",value,status)
         if (status == 0 )  then
            read(unit=value,fmt=*) corep%rcore
         else
            corep%rcore = -1.0_dp
         endif
                                                                              
         call get_value(attributes,"number-of-continuous-derivatives", &
                                   value,status)
         if (status == 0 )  then
            read(unit=value,fmt=*) corep%n_cont_derivs
         else
            corep%n_cont_derivs = -1
         endif

      case ("valence-charge")
         in_valenceCharge = .true.
         valp => pseudo%valence_charge
         rp => valp%rho_val

         call get_value(attributes,"total-charge",value,status)
         if (status /= 0 ) call die("Cannot determine total valence charge")
         read(unit=value,fmt=*) valp%total_charge
                                                                              
      case ("semilocal-potentials")
         in_semilocal = .true.
         allocate(slp)
         
         if (associated(pseudo%semilocal)) then
            qslp => pseudo%semilocal
            do while (associated(qslp%next))
               qslp => qslp%next
            enddo
            qslp%next => slp
         else
            pseudo%semilocal => slp
         endif
            
         current_sl_set = "invalid"
         slp%set = SET_NULL
         call get_value(attributes,"set",value,status)
         if (status == 0 ) then
            current_sl_set = value
            slp%set = setcode_of_string(value)
         endif
         if (debug_parsing) print *, "Found semilocal-potentials set: ", trim(current_sl_set)

         top_flavor = pseudo%header%flavor
         call get_value(attributes,"flavor",value,status)
         if (status == 0 ) then
            top_flavor = value
         endif


      case ("nonlocal-projectors")
         in_nonlocal = .true.

	 ! Allocate new node and add to the end of the linked list
         allocate(nlp)

         if (associated(pseudo%nonlocal)) then
            qnlp => pseudo%nonlocal
            do while (associated(qnlp%next))
               qnlp => qnlp%next
            enddo
            qnlp%next => nlp
         else
            pseudo%nonlocal => nlp
         endif

         current_proj_set = "invalid"
         nlp%set = SET_NULL
         call get_value(attributes,"set",value,status)
         if (status == 0 ) then
            current_proj_set = value
         nlp%set = setcode_of_string(value)
         endif
         if (debug_parsing) print *, "Found nonlocal-projectors set: ", trim(current_proj_set)
         nlp%set = setcode_of_string(value)

      case ("local-potential")
         in_local_potential = .true.
         lop => pseudo%local
         rp => lop%vlocal

         call get_value(attributes,"type",lop%vlocal_type,status)
         if (status /= 0 ) call die("Cannot determine type of local potential")
                                                                              
      case ("local-charge")
         if (.not. in_local_potential) call die("<local-charge> outside <local-potential>")
         in_chlocal = .true.
         lop => pseudo%local
         rp => lop%chlocal

         ! Future expansion: chlocal attributes
                                                                              
      case ("pseudo-wave-functions")
         in_pseudowavefun = .true. 

         wfp => pseudo%pswfs

         current_wf_set = "invalid"
         wfp%set = SET_NULL
         call get_value(attributes,"set",value,status)
         if (status == 0 ) then
            current_wf_set = value
            wfp%set = setcode_of_string(value)
         endif
         if (debug_parsing) print *, "Found pseudo-wavefunction set: ", trim(current_wf_set)

      case ("annotation")
         if (in_provenance) then
            call save_annotation(attributes,pp%annotation)
         else if (in_grid) then
            call save_annotation(attributes,gannot)
         else if (in_xc) then
            call save_annotation(attributes,xp%annotation)
         else if (in_valence_config) then
            call save_annotation(attributes,cp%annotation)
         else if (in_semilocal) then
            call save_annotation(attributes,slp%annotation)
         else if (in_nonlocal) then
            call save_annotation(attributes,nlp%annotation)
         else if (in_local_potential) then
            call save_annotation(attributes,lop%annotation)
         else if (in_pseudowavefun) then
            call save_annotation(attributes,wfp%annotation)
         else if (in_valenceCharge) then
            call save_annotation(attributes,valp%annotation)
         else if (in_coreCharge) then
            call save_annotation(attributes,corep%annotation)
         else if (in_psml) then  ! It must be at the top level
            call save_annotation(attributes,pseudo%annotation)
         else
            ! Do nothing instead of dying
            ! call die("Misplaced <annotation> element")
         endif
                  
end select

parent_element = name

end subroutine begin_element
!----------------------------------------------------------------------

#ifndef PSML_USE_FOX
subroutine end_element(name)
#else
subroutine end_element(namespaceURI,localName,name)
character(len=*), intent(in)    :: namespaceURI
character(len=*), intent(in)    :: localName
#endif

character(len=*), intent(in)     :: name

integer :: i

if (debug_parsing) print *, "-- end Element: ", trim(name)

select case(name)

      case ("radfunc")
         in_radfunc = .false.
         if (.not. associated(rp%data)) then
            call die("No data for radfunc!")
         endif

      case ("grid")
         in_grid = .false.
         !
         if (.not. got_explicit_grid_data) then
            call die("Need explicit grid data!")
         endif
         call delete(tmp_grid)

      case ("data")
      !
      ! We are done filling up the radfunc data
      ! Check that we got the advertised number of items
      !
         in_data = .false.
         if (ndata /= size(rp%data)) then
            call die("npts mismatch in radfunc data")
         endif

      case ("grid-data")
      !
      ! We are done filling up the grid data
      ! Check that we got the advertised number of items
      !
         in_grid_data = .false.
         if (ndata_grid /= size(gdata)) then
            call die("npts mismatch in grid")
         endif
         if (debug_parsing) print *, "Got grid data: ", got_explicit_grid_data

      case ("pseudocore-charge")
         in_coreCharge = .false.

      case ("valence-charge")
         in_valenceCharge = .false.

      case ("semilocal-potentials")
         in_semilocal = .false.
         slp => null()

      case ("nonlocal-projectors")
         in_nonlocal = .false.
         nlp => null()

      case ("slps")
         in_slps = .false.

      case ("proj")
         in_proj = .false.

      case ("local-potential")
         in_local_potential = .false.
         lop => null()

      case ("local-charge")
         in_chlocal = .false.

      case ("pseudo-wave-functions")
         in_pseudowavefun = .false. 

      case ("pswf")
         in_pswf = .false.

      case ("valence-configuration")
         in_valence_config = .false.

      case ("exchange-correlation")
         in_xc = .false.

      case ("libxc-info")
         in_libxc_info = .false.
         if (n_funct /= xp%n_functs_libxc) &
            call die("Too few <functional> elements in <libxc-info>")

      case ("provenance")
         in_provenance = .false.

      case ("input-file")
         in_input_file = .false.
         
      case ("header")
         in_header = .false.

      case ("psml")
         in_psml = .false.
!         call dump_pseudo(pseudo)

end select

end subroutine end_element
!----------------------------------------------------------------------

subroutine pcdata_chunk(chunk)
#ifdef PSML_USE_FOX
  use fox_extra, only: build_data_array
#else
  use xmlf90_sax, only: build_data_array
#endif
  use iso_varying_string, only: operator(//)
  
character(len=*), intent(in) :: chunk

if (len_trim(chunk) == 0) RETURN     ! skip empty chunk

if (in_data) then
!
! Note that we know where we need to put it through the pointer rp...
!
      call build_data_array(chunk,rp%data,ndata)

else if (in_grid_data) then
!
!     Fill the explicit grid data pointer
      call build_data_array(chunk,gdata,ndata_grid)

else if (in_input_file) then

      ifp%buffer = ifp%buffer // chunk 
   
else if (in_header) then
      !
      ! There should not be any pcdata in header in this version...

!      print *, "Header data:"
!      print *, trim(chunk)

endif

end subroutine pcdata_chunk
!
subroutine cdata_section_chunk(chunk)
  use iso_varying_string, only: operator(//)
  
character(len=*), intent(in) :: chunk

if (len_trim(chunk) == 0) RETURN     ! skip empty chunk

if (in_input_file) then
   ifp%buffer = ifp%buffer // chunk 
endif

end subroutine cdata_section_chunk
!----------------------------------------------------------------------

     ! Annotations are encoded as an association list
     ! in a couple of arrays
     ! ( (key "value") (key "value") ...)
     ! 
subroutine save_annotation(atts,annotation)
  use assoc_list, ps_annotation_t => assoc_list_t
#ifdef PSML_USE_FOX
  use Fox_common, only: dictionary_t, len=>getLength
  use fox_extra, only: get_value=>get_value_by_index, get_key
#else
  use xmlf90_sax, only: dictionary_t, get_value, get_key, len
#endif
       type(dictionary_t), intent(in) :: atts
       type(ps_annotation_t), intent(out) :: annotation
       
       integer :: n, i, status
       character(len=300) :: key, value

       n = len(atts)
       call assoc_list_init(annotation,n,status)
       if (status /= 0) call die("Failed to init annotation object")
       do i = 1, n
          call get_key(atts,i,key,status)
          if (status /= 0) call die("cannot get key in atts dict")
          call get_value(atts,i,value,status)
          if (status /= 0) call die("cannot get value in atts dict")
          call assoc_list_insert(annotation,key,value,status)
          if (status /= 0) call die("Failed to insert annotation pair")
       enddo
     end subroutine save_annotation

end module m_psml_parsing_helpers












