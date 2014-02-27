PROGRAM MAIN
  implicit none
  integer :: it, it1, it2
  integer :: iloop
  integer, parameter :: Nt = 16
  complex*16 :: G0(0:Nt-1), W0(0:Nt-1), Gam0(0:Nt-1, 0:Nt-1)
  complex*16 :: G(0:Nt-1)

  !-------- initialize ---------
  open(10, file="G0_t.dat")
  do it = 0, Nt-1
    read(10, *) G0(it)
  enddo
  close(10)

  W0(:) = (0.d0, 0.d0)
  W0(0) = (0.25d0, 0.d0)

  Gam0(:,:) = (0.d0, 0.d0)
  Gam0(0,0) = (1.d0, 0.d0)

  !-------- integral ------------
  do iloop = 1, 10
    do it = 0, Nt-1
      G(it) = G0(it)
      do it1 = 0, it
        G(it) = G(it) - 3.d0*G0(it1)*G(0)*G(it-it1)*W0(0)*Gam0(0,0)
      enddo
    enddo
  enddo

  open(11, file="G_t_test.dat")
  do it = 0, Nt-1
    write(11, *) G(it)
  enddo
  close(11)

  return
END PROGRAM
      



