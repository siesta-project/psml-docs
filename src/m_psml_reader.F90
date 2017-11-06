#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module m_psml_reader

  public :: psml_reader

  CONTAINS

  subroutine psml_reader(fname,ps,debug)

  use m_psml_core,            only: ps_t, ps_destroy
  use m_psml_tables,          only: ps_GenerateTables

  use m_psml_parsing_helpers, only: begin_element, end_element, pcdata_chunk
  use m_psml_parsing_helpers, only: cdata_section_chunk
  use m_psml_parsing_helpers, only: pseudo, debug_parsing

  use external_interfaces,    only: die => psml_die
  use m_psml_interp,               only: set_default_interpolator

  use xmlf90_sax,        only: xml_t, open_xmlfile, xml_parse, close_xmlfile

  implicit none 

  character(len=*), intent(in)        :: fname
  type(ps_t), intent(inout), target   :: ps 
  logical, intent(in), optional       :: debug

  type(xml_t)                     :: fxml
  integer :: iostat

  ! Clean the object's internal data
  ! Note that the inout intent allow us
  ! to do this, and avoid having ps being
  ! reset by the compiler

  call ps_destroy(ps)

  ! Associate module pointer, so that the parsed data
  ! is written to ps
  pseudo => ps
  if (present(debug)) then
     debug_parsing = debug
  else
     debug_parsing = .false.
  endif

 ! Allocate internal structures here...
 
 call open_xmlfile(fname,fxml,iostat)
 if (iostat /=0) call die("Cannot open PSML file: " // trim(fname))
 call xml_parse(fxml, begin_element,end_element,pcdata_chunk, &
                cdata_section_handler=cdata_section_chunk,verbose=.false.)
 call close_xmlfile(fxml)

 ! Clean up association of module pointer
 pseudo => null()

 call ps_GenerateTables(ps)

 !
 ! Set default interpolator
 !
 call set_default_interpolator()

end subroutine psml_reader

end module m_psml_reader
