! The PSML library calls a "psml_die" procedure pointer when it encounters
! an error. The procedure pointer can be made to point to a given error handler
! by calling the 'ps_set_error_handler' routine.

! The handler should take care of carrying out any needed cleaning and
! terminating the program.  As the details would vary with the client
! program, each program should provide its own and call the
! 'ps_set_error_handler' routine.
! 
! This is an example implementation that could be used for serial programs.

subroutine custom_psml_die(str)
  character(len=*), intent(in) :: str

  write(0,"(a,a)") "Custom psml_die:: ", str
  write(6,"(a,a)") "Custom psml_die:: ", str
  STOP

end subroutine custom_psml_die
