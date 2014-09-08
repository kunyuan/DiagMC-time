INCLUDE "vrbls_mc.f90"
INCLUDE "mylib/mylib.f90"
program main
    USE string_basic
    USE vrbls_mc
    use logging_module
    integer :: ix,iy

    L(1:3) = 2
    Beta=1.d0

    Vol = 1.d0
    do i = 1, D
      dVol(i) = Vol

      Vol = Vol *L(i)
      dL(i) = Floor(L(i)/2.d0)
    enddo

    do i = D+1, 3
      dVol(i) = Vol
    enddo

    allocate(W(NTypeW, 0:Vol-1, 0:MxT-1))
    allocate(Gam(NTypeGam, 0:Vol-1, 0:MxT-1, 0:MxT-1))

    print *, "initialize G, W, Gamma"
    print *, "initialize W"
    call ini_W
    print *, "initialize G"
    call ini_G
    print *, "initialize Gamma"
    call ini_Gamma

    call print_G_table
    call print_G_int

    call print_W_table
    call print_W_int

    call print_Gam_table
    call print_Gam_int


CONTAINS
INCLUDE "model_dependent.f90"

SUBROUTINE ini_W
  implicit none
  integer :: x, y, z, tau, site

  do tau = 0, MxT-1
    W(1:NtypeW, :, tau) = dcmplx(((real(tau)+0.5d0)*Beta/MxT)**2.d0+1.d0, 0.d0)
  enddo
  return
END SUBROUTINE ini_W

SUBROUTINE ini_G
  implicit none
  integer :: tau

  do tau = 0, MxT-1
    G(1:NtypeG, tau) =W(1, 0, tau)*cdexp(dcmplx(0.d0, Pi*(real(tau)+0.5d0)/real(MxT)))
  enddo
  return
END SUBROUTINE ini_G

SUBROUTINE ini_Gamma
  implicit none
  integer :: tau1, tau2

  do tau1 = 0, MxT-1
    do tau2 = 0, MxT-1
      Gam(1:NtypeGam, :, tau1, tau2) = dcmplx((real(tau1)*Beta/MxT-0.5d0*Beta)**2.d0 &
        & +(real(tau2)*Beta/MxT-0.5d0*Beta)**2.d0, 0.d0)*cdexp(dcmplx(0.d0,  &
        & Pi*(real(tau1+tau2))/real(MxT)))
    enddo
  enddo
  return
END SUBROUTINE ini_Gamma


SUBROUTINE print_W_table
  implicit none
  integer :: tau

  open(10, status="replace", file="W_table.dat")
  do tau = -MxT, MxT-1
    if(tau<0)  then
      write(10, *) (tau+0.5d0)*Beta/MxT, W(1, 0, tau+MxT)
    else
      write(10, *) (tau+0.5d0)*Beta/MxT, W(1, 0, tau)
    endif
  enddo
  write(10, *)
  close(10)
  return
END SUBROUTINE print_W_table

SUBROUTINE print_G_table
  implicit none
  integer :: tau

  open(10, status="replace", file="G_table.dat")
  do tau = -MxT, MxT-1
    if(tau<0)  then
      write(10, *) (tau+0.5d0)*Beta/MxT, -G(1, tau+MxT)
    else
      write(10, *) (tau+0.5d0)*Beta/MxT, G(1, tau)
    endif
  enddo
  write(10, *)
  close(10)
  return
END SUBROUTINE print_G_table

SUBROUTINE print_Gam_table
  implicit none
  integer :: tau1, tau2

  open(10, status="replace", file="Gam_table.dat")
  do tau1 = -MxT, MxT-1
    !do tau2 = -MxT, MxT-1
    tau2 = 0
      if(tau1<0 .and. tau2<0)  then
        write(10, *) tau1*Beta/MxT, Gam(1, 0, tau1+MxT, tau2+MxT)
      else if(tau1<0 .and. tau2>=0) then
        write(10, *) tau1*Beta/MxT, -Gam(1, 0, tau1+MxT, tau2)
      else if(tau2<0 .and. tau1>=0) then
        write(10, *) tau1*Beta/MxT, -Gam(1, 0, tau1, tau2+MxT)
      else
        write(10, *) tau1*Beta/MxT, Gam(1, 0, tau1, tau2)
      endif
    !enddo
    !write(10, *)
  enddo
  write(10, *)
  close(10)
  return
END SUBROUTINE print_Gam_table

SUBROUTINE print_W_int
  implicit none
  integer :: t
  double precision :: tau

  open(10, status="replace", file="W_interpolate.dat")
  do t = -2*MxT, 2*MxT-1
    tau = t*Beta/2.d0/MxT
    write(10, *) tau, weight_W(1, 0, tau)
  enddo
  write(10, *)
  close(10)
  return
END SUBROUTINE print_W_int

SUBROUTINE print_G_int
  implicit none
  integer :: t
  double precision :: tau

  open(10, status="replace", file="G_interpolate.dat")
  do t = -2*MxT, 2*MxT-1
    tau = t*Beta/2.d0/MxT
    write(10, *) tau, weight_G(1, tau)
  enddo
  write(10, *)
  close(10)
  return
END SUBROUTINE print_G_int

SUBROUTINE print_Gam_int
  implicit none
  integer :: t1, t2
  double precision :: tau1, tau2

  open(10, status="replace", file="Gam_interpolate.dat")
  do t1 = -2*MxT, 2*MxT-1
    tau1 = t1*Beta/2.d0/MxT
    !do t2 = -2*MxT, 2*MxT-1
      !tau2 = t2*Beta/2.d0/MxT
      write(10, *) tau1, weight_Gam(1, 0, tau1, 0.01d0)
    !enddo
    !write(10, *)
  enddo
  write(10, *)
  close(10)
  return
END SUBROUTINE print_Gam_int
END PROGRAM
