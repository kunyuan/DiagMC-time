module string_basic
  implicit none
  interface str
    module procedure str_int,str_double
  end interface

contains
  function str_int(num,format) result(outstr)
      implicit none
      integer, intent(in) :: num
      character(len=10) :: outstr
      character(len=*),optional,intent(in) :: format
      if(present(format)) then
        write(outstr,trim(adjustl(format))) num
      else
        write(outstr,'(i10)') num
      endif
      outstr=adjustl(outstr)
  end function

  function str_double(num,format) result(outstr)
      implicit none
      double precision, intent(in) :: num
      character(len=25) :: outstr
      character(len=*),optional,intent(in) :: format
      if(present(format)) then
        write(outstr,trim(adjustl(format))) num
      else
        !write(outstr,'(d16.8)') num
        write(outstr,*) num
      endif
      outstr=adjustl(outstr)
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
  contains
    private
    procedure :: set_level_num
    procedure :: set_level_char
    procedure :: set_level_null
    procedure :: quicklog_num
    procedure :: quicklog_char
    procedure :: quicklog_null
    generic,public :: SetLevel => set_level_num,set_level_char,set_level_null
    generic,public :: QuickLog => quicklog_num,quicklog_char,quicklog_null
    procedure,public :: Init
    procedure,public :: SetFileName
    procedure,public :: ClearContent
    procedure,public :: Add
    procedure,public :: AddLine
    procedure,public :: Write => write_log
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
        open(36,file=this%FileName,access='append')
        write(36,101,advance='no') &
            & trim(timestr),trim(LevelName(this%Level)),char(10)//trim(this%Content)
        if(this%Content(tail:tail)/=char(10)) write(36,*)
        close(36)
      endif
101 format('[job]','[',A,']','[',A,']:',A)
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

