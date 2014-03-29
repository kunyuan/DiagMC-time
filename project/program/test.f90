include "mylib.f90"
program main
      use string_basic
      implicit none
      integer :: i
      i=4
      write(*,*) str(i)
end program
