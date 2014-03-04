PROGRAM MAIN
  implicit none
  integer :: it, it1, it2
  integer :: iloop
  integer, parameter :: Nt = 32
  double precision :: ratio
  complex*16 :: G0(0:Nt-1), W0(0:Nt-1), Gam0(0:Nt-1, 0:Nt-1)
  complex*16 :: W(0:Nt-1), Polar(0:Nt-1)

  !-------- initialize ---------
  open(10, file="G0_t.dat")
  do it = 0, Nt-1
    read(10, *) G0(it)
  enddo
  close(10)

  W0(:) = (0.25d0, 0.d0)
  Gam0(:,:) = (0.d0, 0.d0)
  Gam0(0,0) = (1.d0, 0.d0)

  !-------- integral ------------
  ratio = -2.d0
  do it = 0, Nt-1
    Polar(it) = Polar(it) - ratio*G0(it)*G0(-it+Nt)
  enddo

  do it = 0, Nt-1
    do it1 =  0, Nt-1
      do it2 = 0, Nt-1
        if(it2-it1>=0 .and. it-it2>=0) then
          W(it) = W(it) + W0(it1)*Polar(it2-it1)*W0(it-it2)
        else if(it2-it1<0 .and. it-it2>=0) then
          W(it) = W(it) + W0(it1)*Polar(Nt+it2-it1)*W0(it-it2)
        else if(it2-it1>=0 .and. it-it2<0) then
          W(it) = W(it) + W0(it1)*Polar(it2-it1)*W0(Nt+it-it2)
        else 
          W(it) = W(it) + W0(it1)*Polar(Nt+it2-it1)*W0(Nt+it-it2)
        endif

      enddo
    enddo
  enddo


  open(11, file="Polar_t_test.dat")
  do it = 0, Nt-1
    write(11, *) Polar(it)
  enddo
  close(11)

  open(12, file="W_t_test.dat")
  do it = 0, Nt-1
    write(12, *) W(it)
  enddo
  close(12)

  return
END PROGRAM
      



