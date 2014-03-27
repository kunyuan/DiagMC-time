  program main
    implicit none
    integer :: ix, iy, ityp, ibasis, ibasis2, iomega1, iomega2
    integer :: i, itot, iorder
    integer :: MCOrder
    double precision :: ime(200), nme
    integer :: ifile
    integer, parameter :: MxOmega = 16
    double precision :: Gam1MR, GamNorm(200)
    double precision :: GamMC(200, 0:10, 3, 0:3, 0:3, -16:16, -16:16)
    double precision :: GamSqMC(200, 0:10, 3, 0:3, 0:3, -16:16, -16:16)
    double precision :: Gam1, Gam2, Gam3
    character*100 title
    character*1 :: order

    print *, 'title'
    read  *, title
    order = title(lnblnk(title):lnblnk(title))
    read(order, *) MCOrder
    

    open(101, status="old", file="data_collapse/"//trim(title)//"_monte_carlo_data.dat")

    open(300, file=trim(title)//"_coll_monte_carlo_data.dat")
    open(301, file=trim(title)//"_coll_Gamma.dat")

    GamNorm(:) = 0.d0
    ime(:) = 0.d0
    GamMC(:,:,:,:,:,:,:) = 0.d0
    GamSqMC(:,:,:,:,:,:,:) = 0.d0

    itot = 1
    do 
      read(101, *, End=100) ime(itot), GamNorm(itot), Gam1MR
      do iorder = 0, MCOrder
        do ityp = 1, 3 
          do ix = 0, 3
            do iy = 0, 3
              do iomega1 = -MxOmega, MxOmega
                do iomega2 = -MxOmega, MxOmega
                  read(101, *)  GamMC(itot, iorder, ityp, ix, iy, iomega1, iomega2)
                  read(101, *)  GamSqMC(itot, iorder, ityp, ix, iy, iomega1, iomega2)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      write(*, *) itot, ime(itot), GamNorm(itot)
      itot = itot + 1
    enddo
100 continue
    itot = itot - 1
    write(*, *) itot, SUM(ime(:))

    write(300, *) SUM(ime(:)), SUM(GamNorm(:)), Gam1MR
    do iorder = 0, MCOrder
      do ityp = 1, 3
        do ix = 0, 3
          do iy = 0, 3
            do iomega1 = -MxOmega, MxOmega
              do iomega2 = -MxOmega, MxOmega
                write(300, *)  SUM(GamMC(:, iorder, ityp, ix, iy, iomega1, iomega2))
                write(300, *)  SUM(GamSqMC(:, iorder, ityp, ix, iy, iomega1, iomega2))
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo

    write(301, *) "================================================="
    write(301, *) SUM(ime(:)), SUM(GamNorm(:)), Gam1MR
    do iorder = 1, MCOrder
      write(301, *) "Order", iorder
      write(301, *) "(0,0)              (1,0)                 (1,1)"
      do iomega1 = -MxOmega, MxOmega
        Gam1 = SUM(GamMC(:,iorder,1,0,0,iomega1, 0))*Gam1MR/SUM(GamNorm(:))
        Gam2 = SUM(GamMC(:,iorder,1,0,0,0, iomega1))*Gam1MR/SUM(GamNorm(:))
        write(301, *) iomega1, Gam1, Gam2
      enddo
      write(301, *)
    enddo

    close(101)
    close(300)
    close(301)
  END program main

