subroutine psml_die(str)
  character(len=*), intent(in) :: str

  write(0,"(a)") str
  stop

end subroutine psml_die
