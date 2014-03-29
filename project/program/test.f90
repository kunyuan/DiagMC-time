include "mylib.f90"
program main
      use string_basic
      use logging_module
      implicit none
      integer :: i
      character(len=20) :: a
      type(logging) :: log_loop
      i=4
      a=str(i)
      write(*,*) a
      i=8888888
      write(*,*) str(i),a
      !log_loop=logging("test.log","hello")
      call log_loop%SetFileName("test.log")
      print *,log_loop%FileName
      print *,log_loop%Content
end program
