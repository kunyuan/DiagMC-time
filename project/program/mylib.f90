module string_basic
  implicit none
contains
  function str(num) result(outstr)
      implicit none
      integer, intent(in) :: num
      character(len=100) :: outstr
      write(outstr,*) num
  end function
end module

