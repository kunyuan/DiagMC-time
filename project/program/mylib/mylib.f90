module string_basic
  implicit none
  interface str
    module procedure str_int,str_double,str_logic,str_double_complex
  end interface

  interface operator (+)
    module procedure add_str
  end interface

contains

  function add_str(str1,str2) result (str3)
    implicit none
      character(len=*),intent(in) :: str1
      character(len=*),intent(in) :: str2
      character(len=1024) :: str3
      str3=trim(str1)//trim(str2)
  end function add_str

  function str_int(num,format) result(outstr)
      implicit none
      integer, intent(in) :: num
      character(len=15) :: outstr
      character(len=*),optional,intent(in) :: format
      if(present(format)) then
        write(outstr,trim(adjustl(format))) num
      else
        write(outstr,'(i10)') num
      endif
      outstr=" "//adjustl(outstr)
  end function

  function str_logic(logic) result(outstr)
      implicit none
      logical, intent(in) :: logic 
      character(len=10) :: outstr
      if(logic==.true.) then
        outstr=" .True." 
      else
        outstr=" .False." 
      endif
  end function

  function str_double(num,format) result(outstr)
      implicit none
      double precision, intent(in) :: num
      character(len=32) :: outstr
      character(len=*),optional,intent(in) :: format
      if(present(format)) then
        write(outstr,trim(adjustl(format))) num
      else
        !write(outstr,'(d16.8)') num
        write(outstr,*) num
      endif
      outstr=" "//adjustl(outstr)
  end function

  function str_double_complex(num,format) result(outstr)
      implicit none
      complex*16, intent(in) :: num
      character(len=64) :: outstr
      character(len=*),optional,intent(in) :: format
      if(present(format)) then
        write(outstr,trim(adjustl(format))) num
        outstr=" "//adjustl(outstr)
      else
        outstr=" ("+str(real(num))+","+str(dimag(num))+")"
      endif
  end function
end module

module logging_module
  use string_basic
  implicit none
  private
  character(len=10),parameter :: LevelName(5) &
    & =(/'DEBUG','INFO','WARNING','ERROR','CRITICAL'/)

  public :: logging

  type :: logging
    private
    character(len=256) :: FileName="*"
    character(len=1024) :: Content=''
    integer :: Level=2
    logical :: status=.true.
  contains
    private
    procedure :: set_level_num
    procedure :: set_level_char
    procedure :: set_level_null
    procedure :: quicklog_num
    procedure :: quicklog_char
    procedure :: quicklog_null
    procedure :: write_log
    procedure :: write_log_num
    procedure :: write_log_char
    procedure :: write_stamp
    procedure :: write_stamp_num
    procedure :: write_stamp_char
    generic,public :: SetLevel => set_level_num,set_level_char,set_level_null
    generic,public :: QuickLog => quicklog_num,quicklog_char,quicklog_null
    generic,public :: Write => write_log_num,write_log_char, write_log
    procedure,public :: Init
    procedure,public :: SetFileName
    procedure,public :: ClearContent
    procedure,public :: Add
    procedure,public :: AddLine

    generic,public :: WriteStamp => write_stamp,write_stamp_num,write_stamp_char
    procedure,public :: WriteLine
    procedure,public :: WriteClose
  end type logging

contains
  subroutine Init(this,filename,content,level)
      implicit none
      class(logging) :: this 
      character(len=*),optional,intent(in) :: filename,content
      integer,optional :: level

      if(present(filename)) then
        call this%SetFileName(filename)
      else
        call this%SetFileName('*')
      endif

      if(present(content)) then
        this%Content=content
      else
        this%Content=''
      endif

      if(present(level)) then
        call this%SetLevel(level)
      else
        call this%SetLevel(2)
      endif
      this%status=.false.
  end subroutine

  subroutine SetFilename(this,str)
      implicit none
      class(logging) :: this
      character(len=*),intent(in) :: str
      this%filename=trim(adjustl(str))
  end subroutine

  subroutine set_level_num(this,num)
      implicit none
      class(logging) :: this
      integer,intent(in) :: num
      if(num<=1) then
        this%Level=1
      elseif(num>=5) then
        this%Level=5
      else
        this%Level=num
      endif
  end subroutine

  subroutine set_level_char(this,str)
      implicit none
      class(logging) :: this
      character,intent(in) :: str
      if(str=='d' .or. str=='D') then
        this.Level=1
      elseif(str=='i' .or. str=='I') then
        this.Level=2
      elseif(str=='w' .or. str=='W') then
        this.Level=3
      elseif(str=='e' .or. str=='E') then
        this.Level=4
      else
        this.Level=5
      endif
  end subroutine

  subroutine set_level_null(this)
      implicit none
      class(logging) :: this
      this%Level=2
  end subroutine

  subroutine ClearContent(this)
      implicit none
      class(logging) :: this 
      this%Content=''
  end subroutine

  subroutine Add(this, str)
      implicit none
      class(logging) :: this
      character(len=*),intent(in) :: str
      this%Content=trim(this%Content)//trim(str)
  end subroutine

  subroutine AddLine(this, str)
      implicit none
      class(logging) :: this
      character(len=*),intent(in) :: str
      this%Content=trim(this%Content)//"    "//trim(str)//char(10)
  end subroutine

  subroutine write_log_num(this,num)
      implicit none
      class(logging) :: this
      integer,intent(in) :: num
      call this%SetLevel(num)
      call this%Write()
  end subroutine

  subroutine write_log_char(this,c)
      implicit none
      class(logging) :: this
      character,intent(in) :: c
      call this%SetLevel(c)
      call this%Write()
  end subroutine

  subroutine write_log(this)
      implicit none
      class(logging) :: this
      integer,dimension(8) :: values
      character(len=100) :: timestr
      integer :: tail
      call date_and_time(VALUES=values)
      timestr=trim(adjustl(str(values(1)-2000,'(i2.2)')))//"/"  &
          &  //trim(adjustl(str(values(2),'(i2.2)')))//'/'//trim(adjustl(str(values(3),'(i2.2)'))) &
          &  //' '//trim(adjustl(str(values(5),'(i2.2)')))//':' &
          &  //trim(adjustl(str(values(6),'(i2.2)')))//':'//trim(adjustl(str(values(7),'(i2.2)')))
      tail=len(this%Content)
      if(this%FileName=='*') then
        write(*,101,advance='no') &
            & trim(timestr),trim(LevelName(this%Level)),char(10)//trim(this%Content)
        if(this%Content(tail:tail)/=char(10)) write(*,*)
        !write(*,"(A)",advance='no') "[job]["//trim(timestr)//"]["//trim(LevelName(this%Level)) &
            !& //"]"//trim(this%Content)
        !print *,"hello",this%Content
      else
        if(this%status) then
          close(36)
        endif
        open(36,file=this%FileName,access='append')
        this%status=.true.
        write(36,101,advance='no') &
            & trim(timestr),trim(LevelName(this%Level)),char(10)//trim(this%Content)
        if(this%Content(tail:tail)/=char(10)) write(36,*)
        close(36)
        this%status=.false.
      endif

      call this%ClearContent()
      call this%SetLevel('i')
101 format('[job]','[',A,']','[',A,']:',A)
  end subroutine

  subroutine write_stamp_num(this,num)
      implicit none
      class(logging) :: this
      integer,intent(in) :: num
      call this%SetLevel(num)
      call this%WriteStamp()
  end subroutine

  subroutine write_stamp_char(this,c)
      implicit none
      class(logging) :: this
      character,intent(in) :: c
      call this%SetLevel(c)
      call this%WriteStamp()
  end subroutine

  subroutine write_stamp(this)
      implicit none
      class(logging) :: this
      integer,dimension(8) :: values
      character(len=100) :: timestr
      integer :: tail
      call date_and_time(VALUES=values)
      timestr=trim(adjustl(str(values(1)-2000,'(i2.2)')))//"/"  &
          &  //trim(adjustl(str(values(2),'(i2.2)')))//'/'//trim(adjustl(str(values(3),'(i2.2)'))) &
          &  //' '//trim(adjustl(str(values(5),'(i2.2)')))//':' &
          &  //trim(adjustl(str(values(6),'(i2.2)')))//':'//trim(adjustl(str(values(7),'(i2.2)')))
      if(this%FileName=='*') then
        write(*,101,advance='no') &
            & trim(timestr),trim(LevelName(this%Level))
        write(*,*)
      else
        open(36,file=this%FileName,access='append')
        this%status=.true.
        write(36,101,advance='no') &
            & trim(timestr),trim(LevelName(this%Level))
        write(36,*)
      endif
      call this%ClearContent()
      call this%SetLevel('i')
101 format('[job]','[',A,']','[',A,']:')
  end subroutine

  subroutine WriteLine(this,str)
      implicit none
      class(logging) :: this
      character(len=*),intent(in) :: str
      if(this%FileName=='*') then
        write(*,"(A)",advance='no') "    "//trim(str)//char(10) 
      else
        if(this%status) then
          write(36,"(A)",advance='no') "    "//trim(str)//char(10)  
        else
          call this%ClearContent
          call this%QuickLog("Please open"+' '//this.FileName+" before logging!")
        endif
      endif
  end subroutine

  subroutine WriteClose(this)
      implicit none
      class(logging) :: this
      if(this%status) then
        close(36)
      else
        call this%ClearContent
        call this%QuickLog("Please open"+' '//this.FileName+" before close!")
      endif
      call this%ClearContent
      this%status=.false.
  end subroutine

  subroutine quicklog_num(this,str,num)
      implicit none
      class(logging) :: this
      character(len=*),intent(in) :: str
      integer,intent(in) :: num
      this%Content=''
      call this%SetLevel(num)
      call this%AddLine(str)
      call this%Write
  end subroutine

  subroutine quicklog_char(this,str,c)
      implicit none
      class(logging) :: this
      character(len=*),intent(in) :: str
      character,intent(in) :: c
      this%Content=''
      call this%SetLevel(c)
      call this%AddLine(str)
      call this%Write
  end subroutine

  subroutine quicklog_null(this,str)
      implicit none
      class(logging) :: this
      character(len=*),intent(in) :: str
      this%Content='' 
      this%Level=2
      call this%AddLine(str)
      call this%Write
  end subroutine

end module

