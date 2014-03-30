include "mylib/mylib.f90"
program main
      use string_basic
      use logging_module
      implicit none
      integer :: i
      double precision :: f
      character(len=20) :: a
      type(logging) :: LogLoop
      !call LogLoop%init("test.log")
      call LogLoop%init("*")
      call LogLoop%SetLevel('e')


      call LogLoop%AddLine("================================================")
      call LogLoop%AddLine("================================================")
      call LogLoop%ClearContent()
      
      call LogLoop%AddLine("================================================")
      call LogLoop%AddLine("MC steps: "//str(234,'(i5)'))
      call LogLoop%AddLine("Efficiency: "//str(13.6d0))
      call LogLoop%Add("    This will also work.")
      call LogLoop%Add("I am a double number: "//trim(str(13.638734737d19)))
      call LogLoop%Write()

      call LogLoop%QuickLog("I am QuickLog.")
      call LogLoop%QuickLog("I am QuickLog,too",'e')
end program
