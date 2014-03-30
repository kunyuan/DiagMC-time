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
    character(len=256) :: FileName
    character(len=1024) :: Content
    integer :: Level
  contains
    private
    procedure :: set_level_num
    procedure :: set_level_char
    procedure,public :: init
    procedure,public :: SetFileName => set_filename
    generic,public :: SetLevel => set_level_num,set_level_char
    procedure,public :: ClearContent => clear_content
    procedure,public :: Add => add_string
    procedure,public :: AddLine => add_line
    procedure,public :: Write => write_log
  end type logging

contains
  subroutine init(this,filename,content,level)
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
        call this%SetLevel(1)
      endif
  end subroutine

  subroutine ClearContent(this)
      class(logging) :: this 
      this%Content=''
  end subroutine

  subroutine set_filename(this,str)
      class(logging) :: this
      character(len=*),intent(in) :: str
      this%filename=trim(adjustl(str))
  end subroutine

  subroutine set_level_num(this,num)
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

  subroutine clear_content(this)
      class(logging) :: this 
      this%Content=''
  end subroutine

  subroutine add_string(this, str)
      class(logging) :: this
      character(len=*),intent(in) :: str
      this%Content=trim(this%Content)//trim(adjustl(str))
  end subroutine

  subroutine add_line(this, str)
      class(logging) :: this
      character(len=*),intent(in) :: str
      this%Content=trim(this%Content)//"    "//trim(adjustl(str))//char(10)
  end subroutine

  subroutine write_log(this)
      class(logging) :: this
      integer,dimension(8) :: values
      character(len=100) :: timestr
      call date_and_time(VALUES=values)
      timestr=trim(adjustl(str(values(1)-2000,'(i2.2)')))//"/"  &
          &  //trim(adjustl(str(values(2),'(i2.2)')))//'/'//trim(adjustl(str(values(3),'(i2.2)'))) &
          &  //' '//trim(adjustl(str(values(5),'(i2.2)')))//':' &
          &  //trim(adjustl(str(values(6),'(i2.2)')))//':'//trim(adjustl(str(values(7),'(i2.2)')))
      if(this%FileName=='*') then
        write(*,101,advance='no') &
            & trim(timestr),trim(LevelName(this%Level)),trim(adjustl(this%Content))
        !write(*,"(A)",advance='no') "[job]["//trim(timestr)//"]["//trim(LevelName(this%Level)) &
            !& //"]"//trim(this%Content)
      else
        open(36,file=this%FileName,access='append')
        write(36,101,advance='no') &
            & trim(timestr),trim(LevelName(this%Level)),trim(adjustl(this%Content))
        write(36,*)
        close(36)
      endif
101 format('[job]','[',A,']','[',A,']:',A)
  end subroutine
end module

