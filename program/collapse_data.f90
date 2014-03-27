INCLUDE "vrbls_mc.f90"
PROGRAM MAIN
  USE vrbls_mc
  program main
    implicit none
    integer :: ix, iy, ityp, it1, it2
    integer :: i, itot, iorder
    integer :: ifile
    integer, parameter :: Mxjobs = 200
    double precision :: imctmp
    complex*16 :: iGam
    complex*16 :: iGamNorm, iGamNormWeight
    double precision :: iReGamSq
    double precision :: iImGamSq

    character*100 title
    character*1 :: order

    print *, 'title'
    read  *, title
    order = title(lnblnk(title):lnblnk(title))
    read(order, *) MCOrder

    title_mc = trim(adjustl(title))//"_coll"
    
    open(101, status="old", file="data_collapse/"//trim(title)//"_monte_carlo_data.dat", form="binary")

    itot = 1
    read(101, *) Lx, Ly

    allocate(GamMC(0:MCOrder,1:NTypeGam/2, 0:Lx-1, 0:Ly-1, 0:MxT-1, 0:MxT-1))
    allocate(ReGamSqMC(0:MCOrder,1:NTypeGam/2, 0:Lx-1, 0:Ly-1, 0:MxT-1, 0:MxT-1))
    allocate(ImGamSqMC(0:MCOrder,1:NTypeGam/2, 0:Lx-1, 0:Ly-1, 0:MxT-1, 0:MxT-1))

    GamNorm = 0.d0
    imc = 0.d0

    GamMC(:,:,:,:,:,:) = 0.d0
    ReGamSqMC(:,:,:,:,:,:) = 0.d0
    ImGamSqMC(:,:,:,:,:,:) = 0.d0

    do 
      read(101, *) imctmp, iGamNorm, iGamNormWeight
      imc = imc + imctmp
      GamNorm = GamNorm + iGamNorm
      GamNormWeight = iGamNormWeight
      do iorder = 0, MCOrder
        do ityp = 1, 3 
          do ix = 0, Lx-1
            do iy = 0, Ly-1
              do it1 = 0, MxT-1
                do it2 = 0, MxT-1
                  read(101, *)  iGam
                  GamMC(iorder, ityp, ix, iy, it1, it2) = GamMC(iorder, ityp, ix, iy, it1, it2) +  iGam
                  read(101, *)  iReGamSq
                  ReGamSqMC(iorder,ityp,ix,iy,it1,it2)=ReGamSqMC(iorder,ityp,ix,iy,it1,it2)+iReGamSq
                  read(101, *)  iImGamSq
                  ImGamSqMC(iorder,ityp,ix,iy,it1,it2)=ImGamSqMC(iorder,ityp,ix,iy,it1,it2)+iImGamSq
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      read(101, *, End=100)
      itot = itot + 1
    enddo
100 continue
    itot = itot - 1
    write(*, *) itot, imc


    call write_monte_carlo_data

    close(101)

  CONTAINS
  INCLUDE "read_write_data.f90"
  END program main

