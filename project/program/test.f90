include "mylib.f90"
program main
      use string_basic
      use logging_module
      implicit none
      integer :: i
      double precision :: f
      character(len=20) :: a
      type(logging) :: LogLoop
      call LogLoop%init("test.log")
      call LogLoop%SetLevel('e')
      call LogLoop%AddLine("================================================")
      call LogLoop%AddLine("================================================")
      call LogLoop%ClearContent()
      !call LogLoop%AddLine("================================================")
      call LogLoop%AddLine('')
      call LogLoop%AddLine("MC steps: "//str(234,'(i5)'))
      call LogLoop%AddLine("Efficiency: "//str(13.6d0))
      call LogLoop%Write()
end program
