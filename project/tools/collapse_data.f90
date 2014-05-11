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
    integer :: ibin, ibasis
    double precision :: imctmp, Ztmp
    complex*16 :: iGam
    complex*16 :: iGamNorm, iGamNormWeight
    double precision :: iReGamSq
    double precision :: iImGamSq
    character(len=100),dimension(Mxjobs) :: title_file
    logical :: flag,alive
    complex*16, allocatable :: GamTmp(:,:,:,:,:,:)
    complex*16, allocatable :: GamBasisTmp(:,:,:,:,:,:)
    double precision, allocatable :: ReGamSqTmp(:,:,:,:,:,:)
    double precision, allocatable :: ImGamSqTmp(:,:,:,:,:,:)
    integer :: EffectiveSamp

    call LogFile%Initial("project.log","loop.collapse")
    call LogFile%QuickLog("Reading data files list form read_list.dat...")
    inquire(file="read_list.dat",exist=alive)
    if(.not. alive) then
      call LogFile%QuickLog("read_list.dat is not exist!",'e')
      stop -1
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

    flag=.true.
    do i=1,itot
      inquire(file=title_file(i),exist=alive)
      if(alive) then
        open(101, status="old", file=title_file(i), form="binary")
        read(101, iostat=ios) Beta, J2, MCOrder, L(1), L(2)
        if(ios==0) then
          flag=.false.
          close(101)
          exit
        endif
        close(101)
      endif
    enddo 
    if(flag) then
        call LogFile%QuickLog("I can't find any data file contains L(1),L(2) information!",'e')
        stop -1
    endif

    title_mc = str(Beta,'(f4.2)')+'_'+str(J2,'(f4.2)')+'_'+str(MCOrder,'(i1)')+'_coll'

    allocate(GamMC(0:MCOrder,1:NTypeGam/2, 0:L(1)-1, 0:L(2)-1, 0:MxT-1, 0:MxT-1))
    allocate(GamMCBasis(0:MCOrder,1:NTypeGam/2, 0:L(1)-1, 0:L(2)-1, 1:NBinGam, 1:NBasisGam))
    allocate(ReGamSqMC(0:MCOrder,1:NTypeGam/2, 0:L(1)-1, 0:L(2)-1, 0:MxT-1, 0:MxT-1))
    allocate(ImGamSqMC(0:MCOrder,1:NTypeGam/2, 0:L(1)-1, 0:L(2)-1, 0:MxT-1, 0:MxT-1))

    allocate(GamTmp(0:MCOrder,1:NTypeGam/2, 0:L(1)-1, 0:L(2)-1, 0:MxT-1, 0:MxT-1))
    allocate(ReGamSqTmp(0:MCOrder,1:NTypeGam/2, 0:L(1)-1, 0:L(2)-1, 0:MxT-1, 0:MxT-1))
    allocate(ImGamSqTmp(0:MCOrder,1:NTypeGam/2, 0:L(1)-1, 0:L(2)-1, 0:MxT-1, 0:MxT-1))
    allocate(GamBasisTmp(0:MCOrder,1:NTypeGam/2, 0:L(1)-1, 0:L(2)-1, 1:NBinGam, 1:NBasisGam))

    GamNorm = 0.d0
    imc = 0.d0
    EffectiveSamp=0
    Z_normal = 0.d0

    GamMC(:,:,:,:,:,:) = 0.d0
    GamMCBasis(:,:,:,:,:,:) = 0.d0
    ReGamSqMC(:,:,:,:,:,:) = 0.d0
    ImGamSqMC(:,:,:,:,:,:) = 0.d0

    do i = 1, itot
      ios=0
      inquire(file=title_file(i),exist=alive)
      if(alive) then
        open(101, status="old", file=title_file(i), form="binary")

        read(101) Beta, J2, MCOrder, L(1), L(2)
        read(101,iostat=ios) imctmp, iGamNorm, iGamNormWeight
        read(101) Ztmp, ratioerr
        do it2 = 0, MxT-1
          do it1 = 0, MxT-1
            do iy = 0, L(2)-1
              do ix = 0, L(1)-1
                do ityp = 1, NtypeGam/2
                  do iorder = 0, MCOrder
                    read(101,iostat=ios)  GamTmp(iorder, ityp, ix, iy, it1, it2)
                    read(101,iostat=ios)  ReGamSqTmp(iorder,ityp,ix,iy,it1,it2)
                    read(101,iostat=ios)  ImGamSqTmp(iorder,ityp,ix,iy,it1,it2)
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo

        do ibasis = 1, NBasisGam
          do ibin = 1, NBinGam
            do iy = 0, L(2)-1
              do ix = 0, L(1)-1
                do ityp = 1, NtypeGam/2
                  do iorder = 0, MCOrder
                    read(101,iostat=ios)  GamBasisTmp(iorder, ityp, ix, iy, ibin, ibasis)
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
        close(101)
      endif
      if(ios==0) then
        EffectiveSamp=EffectiveSamp+1
        imc = imc + imctmp
        Z_normal = Z_normal+Ztmp
        GamNorm = GamNorm + iGamNorm
        GamNormWeight = iGamNormWeight

        GamMC=GamMC+GamTmp
        ReGamSqMC=ReGamSqMC+ReGamSqTmp
        ImGamSqMC=ImGamSqMC+ImGamSqTmp
        GamMCBasis = GamMCBasis+GamBasisTmp
      endif
    enddo

    call LogFile%QuickLog("File Num:"+str(itot)+" ,Effective Samples Num:"+str(EffectiveSamp))
    if(EffectiveSamp>itot/2) then
      call LogFile%QuickLog("Writing collpased data...")
      call write_monte_carlo_data
      call LogFile%QuickLog("Writing GamMC data...")
      call write_GamMC
      call LogFile%QuickLog("Collpasing data is done!")
      stop 0
    else
      call LogFile%QuickLog("Too few Effective Samples, self consistent loop aborted!",'w')
      stop -1
    endif
  CONTAINS

  SUBROUTINE write_monte_carlo_data
    implicit none
    integer :: iorder, itopo, ix, iy, ityp, it1, it2

    !=========== write into files =========================================
    open(104, status="replace", &
      & file=trim(title_mc)//"_monte_carlo_data.bin.dat",form="binary")

    write(104) Beta, J2, MCOrder, L(1), L(2)
    write(104) imc, GamNorm, GamNormWeight
    write(104) Z_normal, ratioerr
    do it2 = 0, MxT-1
      do it1 = 0, MxT-1
        do iy = 0, L(2)-1
          do ix = 0, L(1)-1
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
    do ibasis = 1, NBasisGam
      do ibin = 1, NBinGam
        do iy = 0, L(2)-1
          do ix = 0, L(1)-1
            do ityp = 1, NtypeGam/2
              do iorder = 0, MCOrder
                write(104)  GamMCBasis(iorder,  ityp, ix, iy, ibin, ibasis)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    !=========  write on the screen ========================================

    
    close(104)
  END SUBROUTINE write_monte_carlo_data

  SUBROUTINE write_GamMC
    implicit none
    integer :: iorder, itopo, ix, iy, ityp, it1, it2
    complex*16 :: gam

    !=========== write into files =========================================
    open(104, status="replace", file=trim(title_mc)//"_GamMC.dat")
    write(104, *) Beta, J2, MCOrder, L(1), L(2)
    write(104, *) imc, GamNorm, GamNormWeight
    write(104, *) Z_normal, ratioerr
    
    do it1 = 0, MxT-1
      it2 =  it1
      gam = GamMC(1, 1, 0, 0, it1, it2)*GamNormWeight/GamNorm
      write(104, '(i3,E20.10E3,"    +i",E20.10E3)') it1, real(gam), dimag(gam)
    enddo
    write(104, *)
    close(104)
  END SUBROUTINE write_GamMC

  END program main

