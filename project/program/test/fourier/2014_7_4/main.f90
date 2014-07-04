
MODULE test
    implicit none
    integer,parameter :: Omega=2
end module

INCLUDE "vrbls_mc.f90"
INCLUDE "mylib/mylib.f90"
program main
    use test
    USE string_basic
    USE vrbls_mc
    use logging_module
    integer :: ix,iy

    L(1:3) = 2
    Vol = 1.d0
    Beta=1.d0

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


    !print *, "W: FFT r -> k"
    !print *, "Step1 ==> print W(r)"
    !call print_W_space

    !print *, "Step2 ==> FFT"
    !call transfer_W_r(1)

    !print *, "Step3 ==> print W(k)"
    !call print_W_space

    !print *, "W: FFT k -> r"
    !print *, "Step4 ==> anti-FFT"
    !call transfer_W_r(-1)

    !print *, "Step5 ==> print W(r)"
    !call print_W_space


    print *, "W: FFT tau -> omega"
    print *, "Step1 ==> print W(tau)"
    call print_W_time

    print *, "Step2 ==> FFT"
    call transfer_W_t(1)

    print *, "Step3 ==> print W(omega)"
    call print_W_time

    print *, "W: FFT omega -> tau"

    print *, "Step4 ==> anti-FFT"
    call transfer_W_t(-1)

    print *, "Step5 ==> print W(tau)"
    call print_W_time


    print *, " "

    print *, "G: FFT tau -> omega"
    print *, "Step1 ==> print G(tau)"
    call print_G_time

    print *, "Step2 ==> FFT"
    call transfer_G_t(1)

    print *, "Step3 ==> print G(omega)"
    call print_G_time

    print *, "G: FFT omega -> tau"

    print *, "Step4 ==> anti-FFT"
    call transfer_G_t(-1)

    print *, "Step5 ==> print G(tau)"
    call print_G_time

    !print *, "Gamma: FFT tau -> omega"
    !print *, "Gamma(tau1, tau2)"
    !print *, "Gamma(omega1, omega2)"
    !print *, "Gamma: FFT omega -> tau"
    !print *, "Gamma(omega1, omega2)"
    !print *, "Gamma(tau1, tau2)"

CONTAINS
INCLUDE "sc_functions.f90"

SUBROUTINE ini_W
  implicit none
  integer :: x, y, z, tau, site

  do x = 0, L(1)-1
    do y = 0, L(2)-1
      do z = 0, L(3)-1
        do tau = 0, MxT-1
          site = get_site_from_cord(D, (/x,y,z/))
          W(1:NtypeW, site, tau) = dcmplx((((real(tau)+0.5d0)*Beta/MxT)**2.d0+1.d0) &
            & /(dsqrt(x**2.d0+y**2.d0+z**2.d0)+1.d0), 0.d0)
        enddo
      enddo
    enddo
  enddo
  return
END SUBROUTINE ini_W

SUBROUTINE ini_G
  implicit none
  integer :: tau

  do tau = 0, MxT-1
    G(1:NtypeG, tau) =W(1,0, tau)*cdexp(dcmplx(0.d0, Pi*(real(tau)+0.5d0)/real(MxT)))
  enddo
  return
END SUBROUTINE ini_G

SUBROUTINE print_W_space
  implicit none
  integer :: x, y, z, site

  open(10, access="append", file="W_space.dat")
  do x = 0, L(1)-1
    do y = 0, x
      do z = 0, y
        site = get_site_from_cord(D, (/x,y,z/))
        write(10, *) x, y, z, W(1, site, 0)
      enddo
    enddo
  enddo
  write(10, *)
  close(10)
  return
END SUBROUTINE print_W_space

SUBROUTINE print_W_time
  implicit none
  integer :: tau

  open(10, access="append", file="W_time.dat")
  do tau = 0, MxT-1
    write(10, *) tau, W(1, 0, tau)
  enddo
  write(10, *)
  close(10)
  return
END SUBROUTINE print_W_time

SUBROUTINE print_G_time
  implicit none
  integer :: tau

  open(10, access="append", file="G_time.dat")
  do tau = 0, MxT-1
    write(10, *) tau, G(1,  tau)
  enddo
  write(10, *)
  close(10)
  return
END SUBROUTINE print_G_time
END PROGRAM
