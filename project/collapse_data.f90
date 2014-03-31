INCLUDE "program/vrbls_mc.f90"
INCLUDE "program/mylib/mylib.f90"
PROGRAM MAIN
  USE vrbls_mc
  USE string_basic
  USE logging_module
    implicit none
    integer :: ix, iy, ityp, it1, it2
    integer :: i, itot, iorder
    integer :: ifile,ios
    integer, parameter :: Mxjobs = 200
    double precision :: imctmp
    complex*16 :: iGam
    complex*16 :: iGamNorm, iGamNormWeight
    double precision :: iReGamSq
    double precision :: iImGamSq
    character(len=100),dimension(Mxjobs) :: title_file
    character*1 :: order_coll
    character*6 :: title_coll
    logical :: flag,alive

    call LogFile%Initial("project.log")
    call LogFile%QuickLog("Reading data files list form read_list.dat...")
    inquire(file="read_list.dat",exist=alive)
    if(.not. alive) then
      call LogFile%QuickLog("read_list.dat is not exist!")
      stop
    endif
    open(10, file="read_list.dat")
    i = 1
    do
      read(10, *, End=100) title_file(i)
      i = i+1
    enddo
100 continue
    itot = i-1

    if(itot==0) stop

    order_coll = title_file(1)(8:8)
    read(order_coll, *) MCOrder
    title_coll = title_file(1)(3:8)
    title_mc = trim(adjustl(title_coll))//"_coll"


    flag=.true.
    do i=1,itot
      inquire(file=title_file(i),exist=alive)
      if(alive) then
        open(101, status="old", file=title_file(i), form="binary", iostat=ios)
        read(101, iostat=ios) Lx, Ly
        if(ios==0) then
          flag=.false.
          close(101)
          exit
        endif
        close(101)
      endif
    enddo 
    if(flag) then
        call LogFile%QuickLog("I can't find any data file contains Lx,Ly information!")
        stop
    endif


    allocate(GamMC(0:MCOrder,1:NTypeGam/2, 0:Lx-1, 0:Ly-1, 0:MxT-1, 0:MxT-1))
    allocate(ReGamSqMC(0:MCOrder,1:NTypeGam/2, 0:Lx-1, 0:Ly-1, 0:MxT-1, 0:MxT-1))
    allocate(ImGamSqMC(0:MCOrder,1:NTypeGam/2, 0:Lx-1, 0:Ly-1, 0:MxT-1, 0:MxT-1))

    GamNorm = 0.d0
    imc = 0.d0

    GamMC(:,:,:,:,:,:) = 0.d0
    ReGamSqMC(:,:,:,:,:,:) = 0.d0
    ImGamSqMC(:,:,:,:,:,:) = 0.d0

    do i = 1, itot
      open(101, status="old", file=title_file(i), form="binary")

      read(101) Lx, Ly
      read(101) imctmp, iGamNorm, iGamNormWeight

      imc = imc + imctmp
      GamNorm = GamNorm + iGamNorm
      GamNormWeight = iGamNormWeight

      do iorder = 0, MCOrder
        do ityp = 1, 3 
          do ix = 0, Lx-1
            do iy = 0, Ly-1
              do it1 = 0, MxT-1
                do it2 = 0, MxT-1
                  read(101)  iGam
                  GamMC(iorder, ityp, ix, iy, it1, it2) = GamMC(iorder, ityp, ix, iy, it1, it2) +  iGam
                  read(101)  iReGamSq
                  ReGamSqMC(iorder,ityp,ix,iy,it1,it2)=ReGamSqMC(iorder,ityp,ix,iy,it1,it2)+iReGamSq
                  read(101)  iImGamSq
                  ImGamSqMC(iorder,ityp,ix,iy,it1,it2)=ImGamSqMC(iorder,ityp,ix,iy,it1,it2)+iImGamSq
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      close(101)
    enddo

    write(*, *) itot, imc
    call write_monte_carlo_data


  CONTAINS

  SUBROUTINE write_monte_carlo_data
    implicit none
    integer :: iorder, itopo, ix, iy, ityp, it1, it2

    !=========== write into files =========================================
    open(104, status="replace", &
      & file=trim(title_mc)//"_monte_carlo_data.bin.dat",form="binary")

    write(104) Lx, Ly
    write(104) imc, GamNorm, GamNormWeight
    do it2 = 0, MxT-1
      do it1 = 0, MxT-1
        do iy = 0, Ly-1
          do ix = 0, Lx-1
            do ityp = 1, NtypeGam/2
              do iorder = 0, MCOrder
                write(104)  GamMC(iorder,  ityp, ix, iy, it1, it2)
                write(104)  ReGamSqMC(iorder,  ityp, ix, iy, it1, it2)
                write(104)  ImGamSqMC(iorder,  ityp, ix, iy, it1, it2)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    !=========  write on the screen ========================================

    
    close(104)
  END SUBROUTINE write_monte_carlo_data

  END program main

