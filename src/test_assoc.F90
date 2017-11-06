#if defined HAVE_CONFIG_H
#include "config.h"
#endif

program test_assoc
  use assoc_list
  
  type(assoc_list_t) :: a, b
  character(len=100) :: k, v, val
  integer :: stat

  do i = 1, 20
     write(k,"(a,i0)") "key_", i
     write(v,"(a,i0)") "val_", i
     call assoc_list_insert(a,k,v,stat)
  end do
  
  call assoc_list_get_value(a,10,val,stat)
  b = a
  print *, trim(val)
  call assoc_list_get_value(b,20,val,stat)
  print *, trim(val)
end program test_assoc
