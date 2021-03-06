INCLUDE "program/vrbls_mc.f90"
INCLUDE "program/mylib/mylib.f90"
PROGRAM MAIN
  USE vrbls_mc
  USE string_basic
  USE logging_module
    implicit none
    integer :: ix, iy, ityp, it1, it2, ir
    integer :: i, itot, iorder
    integer :: ifile,ios
    integer, parameter :: Mxjobs = 200
    integer :: ibin, ibasis
    double precision :: imctmp, Ztmp
    double precision :: targetBeta
    complex*16 :: iGam
    complex*16 :: iGamNorm, iGamNormWeight
    double precision :: iReGamSq
    double precision :: iImGamSq
    character(len=100),dimension(Mxjobs) :: title_file
    logical :: flag,alive

    complex*16, allocatable :: GamMCTmp(:,:,:)
    double precision, allocatable :: ReGamSqMCTmp(:,:,:)
    double precision, allocatable :: ImGamSqMCTmp(:,:,:)

    complex*16, allocatable :: GamMCBasisTmp(:,:,:,:,:)
    double precision, allocatable :: ReGamSqTmp(:,:,:,:,:)
    double precision, allocatable :: ImGamSqTmp(:,:,:,:,:)
    integer :: EffectiveSamp

    open(10, status='old', file="input.inp")
    read(10, *) targetBeta
    close(10)

    call LogFile%Initial("project.log","loop.collapse")
    call LogFile%QuickLog("Reading data files list from read_list.dat...")
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
        read(101, iostat=ios) finalBeta, Beta, MCOrder, L(1:D)
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

    Vol = 1
    do i = 1, D
      Vol =  Vol *L(i)
      dL(i) = L(i)/2
    enddo

    call LogFile%QuickLog("Volume: "+str(Vol))

    title_mc = str(finalBeta,'(f4.2)')+'_'+str(MCOrder,'(i1)')+'_coll'

    allocate(GamMC(0:MCOrder, 0:Vol-1, 0:MxT-1))
    allocate(ReGamSqMC(0:MCOrder, 0:Vol-1, 0:MxT-1))
    allocate(ImGamSqMC(0:MCOrder, 0:Vol-1, 0:MxT-1))

    allocate(GamMCTmp(0:MCOrder, 0:Vol-1, 0:MxT-1))
    allocate(ReGamSqMCTmp(0:MCOrder, 0:Vol-1, 0:MxT-1))
    allocate(ImGamSqMCTmp(0:MCOrder, 0:Vol-1, 0:MxT-1))

    allocate(GamMCBasis(0:MCOrder,1:NTypeGam/2, 0:Vol-1, 1:NBinGam, 1:NBasisGam))
    allocate(ReGamSqBasis(0:MCOrder,1:NTypeGam/2, 0:Vol-1, 1:NBinGam, 1:NBasisGam))
    allocate(ImGamSqBasis(0:MCOrder,1:NTypeGam/2, 0:Vol-1, 1:NBinGam, 1:NBasisGam))

    allocate(GamMCBasisTmp(0:MCOrder,1:NTypeGam/2, 0:Vol-1, 1:NBinGam, 1:NBasisGam))
    allocate(ReGamSqTmp(0:MCOrder,1:NTypeGam/2, 0:Vol-1, 1:NBinGam, 1:NBasisGam))
    allocate(ImGamSqTmp(0:MCOrder,1:NTypeGam/2, 0:Vol-1, 1:NBinGam, 1:NBasisGam))

    GamNorm = 0.d0
    imc = 0.d0
    EffectiveSamp=0
    Z_normal = 0.d0

    GamMC(:,:,:) = 0.d0
    ReGamSqMC(:,:,:) = 0.d0
    ImGamSqMC(:,:,:) = 0.d0
    GamMCBasis(:,:,:,:,:) = 0.d0
    ReGamSqBasis(:,:,:,:,:) = 0.d0
    ImGamSqBasis(:,:,:,:,:) = 0.d0

    do i = 1, itot
      ios=0
      inquire(file=title_file(i),exist=alive)
      if(alive) then
        open(101, status="old", file=title_file(i), form="binary")

        read(101) finalBeta, Beta, MCOrder, L(1:D)
        read(101,iostat=ios) imctmp, iGamNorm, iGamNormWeight
        read(101) Ztmp, ratioerr
        do it1 = 0, MxT-1
          do ir = 0, Vol-1
            do iorder = 0, MCOrder
              read(101,iostat=ios)  GamMCTmp(iorder, ir, it1)
              read(101,iostat=ios)  ReGamSqMCTmp(iorder,ir,it1)
              read(101,iostat=ios)  ImGamSqMCTmp(iorder,ir,it1)
            enddo
          enddo
        enddo

        do ibasis = 1, NBasisGam
          do ibin = 1, NBinGam
            do ir = 0, Vol-1
              do ityp = 1, NtypeGam/2
                do iorder = 0, MCOrder
                  read(101,iostat=ios)  GamMCBasisTmp(iorder, ityp, ir, ibin, ibasis)
                  read(101,iostat=ios)  ReGamSqTmp(iorder, ityp, ir, ibin, ibasis)
                  read(101,iostat=ios)  ImGamSqTmp(iorder, ityp, ir, ibin, ibasis)
                enddo
              enddo
            enddo
          enddo
        enddo
        close(101)
      endif
      if(ios==0 .and. Beta==targetBeta) then
        EffectiveSamp=EffectiveSamp+1
        imc = imc + imctmp
        Z_normal = Z_normal+Ztmp
        GamNorm = GamNorm + iGamNorm
        GamNormWeight = iGamNormWeight

        GamMC=GamMC+GamMCTmp
        ReGamSqMC=ReGamSqMC+ReGamSqMCTmp
        ImGamSqMC=ImGamSqMC+ImGamSqMCTmp
        GamMCBasis = GamMCBasis+GamMCBasisTmp
        ReGamSqBasis = ReGamSqBasis+ReGamSqTmp
        ImGamSqBasis = ImGamSqBasis+ImGamSqTmp
      endif
    enddo

    call LogFile%QuickLog("File Num:"+str(itot)+" ,Effective Samples Num:"+str(EffectiveSamp))
    if(EffectiveSamp>itot/2) then
      call logfile%quicklog("writing collpased data...")
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
    integer :: iorder, itopo, ix, iy, ityp, it1, it2, ir

    !=========== write into files =========================================
    open(104, status="replace", &
      & file=trim(title_mc)//"_monte_carlo_data.bin.dat",form="binary")

    write(104) finalBeta, Beta, MCOrder, L(1:D)
    write(104) imc, GamNorm, GamNormWeight
    write(104) Z_normal, ratioerr
    do it1 = 0, MxT-1
      do ir = 0, Vol-1
        do iorder = 0, MCOrder
          write(104)  GamMC(iorder,  ir, it1)
          write(104)  ReGamSqMC(iorder, ir, it1)
          write(104)  ImGamSqMC(iorder, ir, it1)
        enddo
      enddo
    enddo
    do ibasis = 1, NBasisGam
      do ibin = 1, NBinGam
        do ir = 0, Vol-1
          do ityp = 1, NtypeGam/2
            do iorder = 0, MCOrder
              write(104)  GamMCBasis(iorder,  ityp, ir, ibin, ibasis)
              write(104)  ReGamSqBasis(iorder,  ityp, ir, ibin, ibasis)
              write(104)  ImGamSqBasis(iorder,  ityp, ir, ibin, ibasis)
            enddo
          enddo
        enddo
      enddo
    enddo
    close(104)

    !=========  write on the screen ========================================
    !open(105, status="replace", &
      !& file=trim(title_mc)//"_monte_carlo_data.dat")

    !write(105, *) finalBeta, Beta, MCOrder, L(1:D)
    !write(105, *) imc, GamNorm, GamNormWeight
    !write(105, *) Z_normal, ratioerr
    !do it1 = 0, MxT-1
      !do ir = 0, Vol-1
        !do iorder = 0, MCOrder
          !write(105, *)  GamMC(iorder,  ir, it1)
          !write(105, *)  ReGamSqMC(iorder, ir, it1)
          !write(105, *)  ImGamSqMC(iorder, ir, it1)
        !enddo
      !enddo
    !enddo
    !do ibasis = 1, NBasisGam
      !do ibin = 1, NBinGam
        !do ir = 0, Vol-1
          !do ityp = 1, NtypeGam/2
            !do iorder = 0, MCOrder
              !write(105, *)  GamMCBasis(iorder,  ityp, ir, ibin, ibasis)
              !write(105, *)  ReGamSqBasis(iorder,  ityp, ir, ibin, ibasis)
              !write(105, *)  ImGamSqBasis(iorder,  ityp, ir, ibin, ibasis)
            !enddo
          !enddo
        !enddo
      !enddo
    !enddo
    !close(105)
  END SUBROUTINE write_monte_carlo_data

  SUBROUTINE write_GamMC
    implicit none
    integer :: iorder, itopo, ix, iy, ityp, it1, it2
    complex*16 :: gam

    !=========== write into files =========================================
    open(104, status="replace", file=trim(title_mc)//"_GamMC.dat")
    write(104, *) finalBeta, MCOrder, L(1:D)
    write(104, *) imc, GamNorm, GamNormWeight
    write(104, *) Z_normal, ratioerr
    
    do it1 = 0, MxT-1
      it2 =  it1
      gam = GamMC(1, 0, it1)*GamNormWeight/GamNorm
      write(104, '(i3,E20.10E3,"    +i",E20.10E3)') it1, real(gam), dimag(gam)
    enddo
    write(104, *)
    close(104)
  END SUBROUTINE write_GamMC

  END program main

