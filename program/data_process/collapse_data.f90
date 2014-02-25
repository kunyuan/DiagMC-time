  program main
    implicit none
    integer :: itopo, ix, iy, ityp, ibasis, ibasis2, iomega1, iomega2
    integer :: i, itot, iorder
    integer :: MCOrder
    double precision :: ime(100), nme
    integer :: ifile
    integer, parameter :: MxOmega = 16
    double precision :: Gam1MR, GamNorm(100)
    double precision :: GamMC1(100, 0:10, 3, 0:3, 0:3, -16:16, -16:16)
    double precision :: GamSqMC1(100, 0:10, 3, 0:3, 0:3, -16:16, -16:16)
    double precision :: GamMC2(100, 0:10, 3, 0:3, 0:3, -16:16, -16:16)
    double precision :: GamSqMC2(100, 0:10, 3, 0:3, 0:3, -16:16, -16:16)
    double precision :: Gam2Topo(100, 2, -16:16)
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
    GamMC1(:,:,:,:,:,:,:) = 0.d0
    GamSqMC1(:,:,:,:,:,:,:) = 0.d0
    GamMC2(:,:,:,:,:,:,:) = 0.d0
    GamSqMC2(:,:,:,:,:,:,:) = 0.d0

    itot = 1
    do 
      read(101, *, End=100) ime(itot), GamNorm(itot), Gam1MR
      do iorder = 0, MCOrder
        do ityp = 1, 3 
          do ix = 0, 3
            do iy = 0, 3
              do iomega1 = -MxOmega, MxOmega
                do iomega2 = -MxOmega, MxOmega
                  read(101, *)  GamMC1(itot, iorder, ityp, ix, iy, iomega1, iomega2)
                  read(101, *)  GamSqMC1(itot, iorder, ityp, ix, iy, iomega1, iomega2)
                enddo
              enddo
            enddo
          enddo
        enddo

        do ityp = 1, 3 
          do ix = 0, 3
            do iy = 0, 3
              do iomega1 = -MxOmega, MxOmega
                do iomega2 = -MxOmega, MxOmega
                  read(101, *)  GamMC2(itot, iorder, ityp, ix, iy, iomega1, iomega2)
                  read(101, *)  GamSqMC2(itot, iorder, ityp, ix, iy, iomega1, iomega2)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

      do iorder = 1, 2
        do iomega1 = -MxOmega, MxOmega
          read(101, *)  Gam2Topo(itot, iorder, iomega1)
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
                write(300, *)  SUM(GamMC1(:, iorder, ityp, ix, iy, iomega1, iomega2))
                write(300, *)  SUM(GamSqMC1(:, iorder, ityp, ix, iy, iomega1, iomega2))
              enddo
            enddo
          enddo
        enddo
      enddo

      do ityp = 1, 3
        do ix = 0, 3
          do iy = 0, 3
            do iomega1 = -MxOmega, MxOmega
              do iomega2 = -MxOmega, MxOmega
                write(300, *)  SUM(GamMC2(:, iorder, ityp, ix, iy, iomega1, iomega2))
                write(300, *)  SUM(GamSqMC2(:, iorder, ityp, ix, iy, iomega1, iomega2))
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    do iorder = 1, 2
      do iomega1 = -MxOmega, MxOmega
        write(300, *)  SUM(Gam2Topo(:, iorder, iomega1))
      enddo
    enddo

    write(301, *) "================================================="
    write(301, *) SUM(ime(:)), SUM(GamNorm(:)), Gam1MR
    do iorder = 1, MCOrder
      write(301, *) "Order", iorder
      do iomega1 = -MxOmega, MxOmega
        Gam1 = SUM(GamMC1(:,iorder,1,0,0,iomega1,iomega1))*Gam1MR/SUM(GamNorm(:))
        Gam2 = SUM(GamMC2(:,iorder,1,1,0,iomega1,iomega1))*Gam1MR/SUM(GamNorm(:))
        write(301, *) iomega1, Gam1, Gam2
      enddo
      write(301, *)
    enddo

    close(101)
    close(300)
    close(301)
  END program main

