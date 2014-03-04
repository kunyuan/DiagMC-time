PROGRAM MAIN
  implicit none
  integer :: it, it1, it2
  integer :: iloop, i, j 
  integer :: ix, iy, ix1, iy1, dx, dy, dit
  integer, parameter :: Nt = 128
  integer, parameter :: L = 4
  double precision :: Beta = 0.5d0
  double precision :: ratio
  double precision :: rtmp, itmp
  complex*16 :: Polar(0:Nt-1), W0(0:L-1, 0:L-1) 
  complex*16 :: W(0:L-1, 0:L-1, 0:Nt-1)
  complex*16 :: Wnew(0:L-1, 0:L-1, 0:Nt-1)

  !-------- initialize ---------
  open(10, file="Polar_t.dat")
  do it = 0, Nt-1
    read(10, *) i, rtmp, itmp
    Polar(it) = dcmplx(rtmp, itmp)
  enddo
  close(10)

  open(11, file="W_t.dat")
  do iy = 0, L-1
    do ix = 0, L-1
      read(11, *) i, j, rtmp, itmp
      W(ix, iy, :) = dcmplx(rtmp, itmp)
    enddo
  enddo
  close(11)
  Wnew(:,:,:) = (0.d0, 0.d0)

  W0(:, :) = (0.d0, 0.d0)
  W0(1, 0) = (0.25d0, 0.d0)
  W0(3, 0) = (0.25d0, 0.d0)
  W0(0, 1) = (0.25d0, 0.d0)
  W0(0, 3) = (0.25d0, 0.d0)

  !-------- integral ------------

  do iloop = 1, 20
    W(:, :, :) = Wnew(:, :, :)
    Wnew(:,:,:) = (0.d0, 0.d0)
    do it = 0, Nt-1
      do ix = 0, L-1
        do iy = 0, L-1

          do it1 = 0, Nt-1
            do ix1 = 0, L-1
              do iy1 = 0, L-1
                dx = ix-ix1
                if(dx<0) dx = dx+L
                dy = iy-iy1
                if(dy<0) dy = dy+L
                dit = it-it1
                if(dit<0) dit = dit+Nt
                Wnew(ix,iy,it) = Wnew(ix,iy,it)+(Nt/Beta)*W0(ix1,iy1)*Polar(it1)*W(dx,dy,dit)
              enddo
            enddo
          enddo

          do ix1 = 0, L-1
            do iy1 = 0, L-1
              dx = ix-ix1
              if(dx<0) dx = dx+L
              dy = iy-iy1
              if(dy<0) dy = dy+L
              Wnew(ix, iy, it) = Wnew(ix, iy, it)+(Nt/Beta)**2.d0*W0(ix1, iy1)*Polar(it)*W0(dx, dy)
            enddo
          enddo
        enddo
      enddo
    enddo

  enddo

  open(13, file="W_t_test.dat")
  write(13, *) "r=(0, 0)"
  do it = 0, Nt-1
    write(13, *) (it+0.5d0)*Beta/Nt, W(0, 0, it)
  enddo
  write(13, *) "r=(1, 0)"
  do it = 0, Nt-1
    write(13, *) (it+0.5d0)*Beta/Nt, W(1, 0, it)
  enddo
  write(13, *) "r=(1, 1)"
  do it = 0, Nt-1
    write(13, *) (it+0.5d0)*Beta/Nt, W(1, 1, it)
  enddo
  write(13, *) "new"
  write(13, *) "r=(0, 0)"
  do it = 0, Nt-1
    write(13, *) (it+0.5d0)*Beta/Nt, Wnew(0, 0, it)
  enddo
  write(13, *) "r=(1, 0)"
  do it = 0, Nt-1
    write(13, *) (it+0.5d0)*Beta/Nt, Wnew(1, 0, it)
  enddo
  write(13, *) "r=(1, 1)"
  do it = 0, Nt-1
    write(13, *) (it+0.5d0)*Beta/Nt, Wnew(1, 1, it)
  enddo
  close(13)

  return
END PROGRAM
      



