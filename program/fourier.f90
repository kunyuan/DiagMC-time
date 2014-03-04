PROGRAM MAIN
  implicit none
  integer  :: i, k, dir
  double precision  :: j
  character*100 :: title
  double precision, parameter :: pi=3.141592653
  double precision, parameter :: beta=0.5d0
  integer, parameter :: lent = 128 
  integer, parameter :: lenomega = 128 

  double precision, dimension(0:lent-1) :: tau
  complex*16, dimension(0:lent-1):: MT
  complex*16, dimension(-lenomega:lenomega):: MOmega
  double precision :: tmp(1:3)
  integer :: col

  do i = 0, lent-1
    tau(i) = beta*real(i)/real(lent)
    !tau(i) = beta*real(i)/real(lent)
  enddo

  print *, "title"
  read (*,'(a100)') title
  print *, "direction"
  read (*, *) dir
  print *, "colum"
  read (*, *) col


  if(dir==1) then
    open(11, file=trim(title)//".dat", status="old")
    open(19, file=trim(title)//"_t.dat", access="append")

    do i = -lenomega, lenomega
      read(11, *) j, k, tmp(1:3)
      MOmega(i) = dcmplx(tmp(col), 0.d0)
    enddo

    call fourier(dir)

    do k =0, lent-1
      write(19, *) tau(k), real(MT(k)), dimag(MT(k))
    enddo
    
    close(11)
    close(19)

  else if(dir==-1) then
    open(21, file=trim(title)//".dat", status="old")
    open(29, file=trim(title)//"_omega.dat", status='replace')

    do i = 0, lent-1
      read(21, *) j, MT(i)
    enddo

    call fourier(dir)

    do k = 0, lenomega-1
      write(29, *) k, real(MOmega(k)), dimag(MOmega(k))
    enddo

    close(21)
    close(29)
  endif


CONTAINS

SUBROUTINE fourier(flag)
  implicit none
  integer, intent(in) :: flag
  integer :: t, omega
  integer :: i
  double precision :: omegac 

  if(flag==1) then
    do t = 0, lent-1
      MT(t) = (0.d0, 0.d0)
      do omega = -lenomega, lenomega
        !omegac = (2.d0*omega+1.d0)*pi/beta
        omegac = (2.d0*omega)*pi/beta
        MT(t) = MT(t)+ cdexp((0.d0, -1.d0)*omegac*tau(t))*MOmega(omega)/beta
      enddo
    enddo
  else if(flag==-1) then
    do omega = 0, lenomega-1
      omegac = (2.d0*omega+1.d0)*pi/beta
      MOmega(omega) = (0.d0, 0.d0)
      do t = 0, lent-1
        MOmega(omega) = MOmega(omega)+ cdexp((0.d0, 1.d0)*omegac*tau(t))*MT(t)
      enddo
    enddo
  endif
  
  return
END SUBROUTINE fourier

END PROGRAM MAIN

