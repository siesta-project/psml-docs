module m_psml

!
  use m_psml_reader, only: psml_reader 

  use m_psml_core, only: ps_t, ps_destroy

  use m_psml_core, only: str_of_set
  use m_psml_core, only:     SET_SREL     , &
                             SET_NONREL   , &
                             SET_SO       , &
                             SET_LJ       , &
                             SET_UP       , &
                             SET_DOWN     , &
                             SET_SPINAVE  , &
                             SET_SPINDIFF , &
                             SET_USER1    , &
                             SET_USER2    , &
                             SET_ALL

  ! Temporarily needed for psop
  use m_psml_core, only: nonlocal_t, nlpj_t

  use m_psml_api
#ifndef PSML_NO_OLD_API  
  use m_psml_old_api
#endif
  
  use m_psml_dump
  use m_psml_ps_edit

  ! Exported utility types
  use m_psml_core, only: ps_radfunc_t => radfunc_t
  use assoc_list, only: ps_annotation_t => assoc_list_t

  use assoc_list, only: nitems_annotation => assoc_list_nitems
  use assoc_list, only: get_annotation_key => assoc_list_get_key
  use assoc_list, only: get_annotation_value => assoc_list_get_value
  use assoc_list, only: insert_annotation_pair => assoc_list_insert
  use assoc_list, only: init_annotation => assoc_list_init
  use assoc_list, only: reset_annotation => assoc_list_reset
  use assoc_list, only: print_annotation => assoc_list_print
  use assoc_list, only: EMPTY_ANNOTATION => EMPTY_ASSOC_LIST 

  ! Remove common names from namespace
  use class_Grid, dummy_id => id, dummy_name => name

  public

end module m_psml
