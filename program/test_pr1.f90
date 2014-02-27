PROGRAM MAIN
  implicit none
  integer :: it, iomega
  integer, parameter :: Nt=16
  double precision, parameter :: Pi = 3.141592653
  complex*16 :: Gt(0:Nt-1), Gomega(0:Nt-1)

  open(10, file="G_t.dat")
  open(11, file="G_omega_direct.dat")

  do it = 0, Nt-1
    read(10, *) Gt(it)
  enddo

  Gomega(:) = (0.d0, 0.d0)
  do iomega = 0, Nt-1
    do it = 0, Nt-1
      Gomega(iomega) = Gomega(iomega) + Gt(it)*cdexp((0.d0,-1.d0)*(real(it) &
        & *(2.d0*real(iomega)+1.d0)*Pi/real(Nt)))
    enddo
    write(11, *) Gomega(iomega)
  enddo

  close(10)
  close(11)
  return
END PROGRAM
