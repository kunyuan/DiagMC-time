PROGRAM MAIN
  implicit none
  integer :: it, it1, it2
  integer :: iloop
  integer, parameter :: Nt = 16
  double precision :: ratio
  complex*16 :: G0(0:Nt-1), W0(0:Nt-1), Gam0(0:Nt-1, 0:Nt-1)
  complex*16 :: Sigma(0:Nt-1)

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
  ratio = 3.d0*(0.5d0)**2.d0/real(Nt)**2.d0
  do it = 0, Nt-1
    do it1 = 0, Nt-1
      do it2 = 0, Nt-1
        if(it1-it2>=0 .and. it-it2>=0) then
          Sigma(it) = Sigma(it) + ratio*G0(it1)*W0(it2)*Gam0(it1-it2,it-it2)
        else if(it1-it2>=0 .and. it-it2<0) then
          Sigma(it) = Sigma(it) + ratio*G0(it1)*W0(it2)*Gam0(it1-it2,it-it2+Nt)
        else if(it1-it2<0 .and. it-it2>=0) then
          Sigma(it) = Sigma(it) + ratio*G0(it1)*W0(it2)*Gam0(it1-it2+Nt,it-it2)
        else 
          Sigma(it) = Sigma(it) + ratio*G0(it1)*W0(it2)*Gam0(it1-it2+Nt,it-it2+Nt)
        endif
      enddo
    enddo
  enddo


  open(11, file="Sigma_t_test.dat")
  do it = 0, Nt-1
    write(11, *) Sigma(it)
  enddo
  close(11)

  return
END PROGRAM
      



