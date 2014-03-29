module string_basic
  implicit none
contains
  function str(num) result(outstr)
      implicit none
      integer, intent(in) :: num
      character(len=20) :: outstr
      write(outstr,*) num
  end function
end module

module logging_module
  implicit none
  type, public :: logging
      character(len=256) :: FileName
      character(len=1024) :: Content
      integer :: level
    contains
      procedure :: SetFileName => set_filename
  end type logging

contains
  subroutine set_filename(this,str)
      class(logging) :: this
      character(len=*),intent(in) :: str
      this%filename=str
  end subroutine
end module
