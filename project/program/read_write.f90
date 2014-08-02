
!===============================================================
!================= read and write infile and input =============
!===============================================================

SUBROUTINE read_infile(id, isub)
  implicit none
  integer, intent(out) :: id, isub
  open(unit=11, file=trim(infile), action='read')
    read(11,*) ID
    read(11,*) L(1:D)
    read(11,*) Jcp
    read(11,*) iniBeta
    read(11,*) dBeta
    read(11,*) finalBeta
    read(11,*) iniMCOrder
    read(11,*) IsLoad
    read(11,*) isub

    if(ISub==2) then
      !for MC we need the statistics numbers
      read(11,*) Ntoss
      read(11,*) Nsamp
      read(11,*) Nstep
      read(11,*) Seed

      !if we want to use some old statistics, we need the title of the statistics file
      read(11,*) title

      !for MC we need the initial reweighting ratios
      read(11,*) CoefOfWorm
      CoefOfWeight(0) = 1.d0
      read(11,*) CoefOfWeight(1:iniMCOrder)

    elseif(ISub==1 .or. ISub==4) then
      !for GW and output, only the title of Gamma file are needed
      read(11,*) title
    endif
    
    print *, isub, "Read from infile done!"

  close(11)
  return
END SUBROUTINE read_infile



SUBROUTINE read_input(ifmust)
  implicit none
  integer :: ios
  logical :: ifmust

  open(11, status="old", iostat=ios, file="input.inp")
  read(11, *) Beta
  read(11, *) MCOrder
  read(11, *) file_version
  close(11)

  if(ios/=0) then
    if(ifmust) then
      call LogFile%QuickLog(str(ios)+"Fail to read the input file! Stop the code!")
      stop -1
    else
      call LogFile%QuickLog(str(ios)+"Fail to read the input file! continue...")
      return
    endif
  endif

  return
END SUBROUTINE read_input


SUBROUTINE write_input(changeBeta, iffirst)
  implicit none
  integer :: ios 
  logical, intent(in) :: changeBeta, iffirst


  if(.not. iffirst) then
    call read_input(.true.)

    open(10, status='replace', iostat= ios, file='input.inp')
    if(changeBeta) then
      call LogFile%QuickLog("Changing Beta:")
      call LogFile%QuickLog("old Beta:"+str(Beta))
      if(Beta+dBeta<=finalBeta) then
        write(10, *) Beta+dBeta
        call LogFile%QuickLog("change to new Beta:"+str(Beta+dBeta))
      else if(Beta==finalBeta) then
        write(10, *) Beta
        call LogFile%QuickLog("already the target Beta!"+str(Beta))
      else
        call LogFile%QuickLog("error! Beta:"+str(Beta)+", dBeta: "+str(dBeta), 'e')
        stop -1
      endif
    else 
      write(10, *) Beta
    endif

    write(10, *) MCOrder
    write(10, *) file_version+1

  else 
    open(10, status='replace', iostat= ios, file='input.inp')
    write(10, *) Beta
    write(10, *) MCOrder
    write(10, *) 0
  endif

  close(10)

  if(ios/=0) then
    call LogFile%QuickLog(str(ios)+"Fail to write the input into input.inp!")
    stop -1
  endif
  return
END SUBROUTINE write_input


SUBROUTINE write_readlist
  implicit none
  logical :: FEXIST

  INQUIRE(DIRECTORY="readfile",EXIST=FEXIST)
  IF(.not.FEXIST) THEN
    call system("mkdir readfile")
  ENDIF

  !write the pid into a readlist
  open(10,access="append", iostat=ios, file="readfile/"+trim(adjustl(title2))+".dat")
  write(10, *) trim(adjustl(title_mc))//"_monte_carlo_data.bin.dat"
  close(10)
  if(ios/=0) call LogFile%QuickLog("readfile iostat nonzero!",'e')

  return
END SUBROUTINE write_readlist

!!================================================================
!!================= READ/WRITE CONFIGURATIONS ====================
!!================================================================

SUBROUTINE read_GWGamma
  implicit none

  call LogFile%QuickLog("Reading G, W, Gamma...")
  if(read_GW()) call LogFile%QuickLog("Read G, W done!")
  if(read_Gamma()) call LogFile%QuickLog("Read Gamma done!")
  return
END SUBROUTINE read_GWGamma

Logical Function read_GW
  implicit none
  integer :: isite, ityp, it1, it2, ios
  logical :: alive
  complex*16, allocatable :: Gtmp(:,:)
  complex*16, allocatable :: Wtmp(:,:,:)

  inquire(file=trim(title_loop)//"_G_file.dat",exist=alive)
  if(.not. alive) then
    call LogFile%QuickLog("There is no G file yet!",'e')
    stop -1
  endif
  inquire(file=trim(title_loop)//"_W_file.dat",exist=alive)
  if(.not. alive) then
    call LogFile%QuickLog("There is no W file yet!",'e')
    stop -1
  endif

  allocate(Gtmp(1:NtypeG, 0:MxT-1))
  allocate(Wtmp(1:NtypeW, 0:Vol-1, 0:MxT-1))

  open(100, status="old", file=trim(title_loop)//"_G_file.dat")
  open(101, status="old", file=trim(title_loop)//"_W_file.dat")

  do it1 = 0, MxT-1
    do ityp = 1, NTypeG
      read(100, *,iostat=ios) Gtmp(ityp, it1)
    enddo
  enddo

  do it1 = 0, MxT-1
    do isite = 0, Vol-1
      do ityp = 1, NTypeW
        read(101, *,iostat=ios) Wtmp(ityp, isite, it1)
      enddo
    enddo
  enddo

  if(ISub==2) then
    if(ios/=0) then
      call LogFile%QuickLog("Failed to read G,W information!",'e')
      read_GW = .false.
    else 
      G = Gtmp;         deallocate(Gtmp)
      W = Wtmp;         deallocate(Wtmp)
      read_GW = .true.
    endif
  else 
    if(ios/=0) then
      call LogFile%QuickLog("Failed to read G,W information!",'e')
      read_GW = .false.
      close(100)
      close(101)
      stop -1
    else
      G = Gtmp;         deallocate(Gtmp)
      W = Wtmp;         deallocate(Wtmp)
      read_GW = .true.
    endif
  endif

  close(100)
  close(101)
  return
END Function read_GW

LOGICAL FUNCTION read_Gamma
  implicit none
  integer :: isite, ityp, it1, it2, ios, ibin, ibasis
  logical :: alive
  complex*16, allocatable :: GamBasistmp(:,:,:,:)

  inquire(file=trim(title_loop)//"_Gam_file.dat",exist=alive)
  if(.not. alive) then
    call LogFile%QuickLog("There is no Gam file yet!",'e')
    stop -1
  endif

  allocate(GamBasistmp(1:NtypeGam, 0:Vol-1, 1:NbinGam, 1:NbasisGam))

  open(102, status="old", file=trim(title_loop)//"_Gam_file.dat")

  do ibasis = 1, NBasisGam
    do ibin = 1, NbinGam
      do isite = 0, Vol-1
        do ityp = 1, NTypeGam
          read(102, *,iostat=ios) GamBasistmp(ityp, isite, ibin, ibasis)
        enddo
      enddo
    enddo
  enddo

  if(ISub==2) then
    if(ios/=0) then
      call LogFile%QuickLog("Failed to read Gam information!",'e')
      read_Gamma = .false.
      deallocate(GamBasistmp)
    else 
      GamBasis = GamBasistmp;         deallocate(GamBasistmp)
      read_Gamma = .true.
    endif
  else 
    if(ios/=0) then
      call LogFile%QuickLog("Failed to read Gam information!",'e')
      read_Gamma = .false.
      close(102)
      stop -1
    else
      GamBasis = GamBasistmp;         deallocate(GamBasistmp)
      read_Gamma = .true.
    endif
  endif
  close(102)

  if(read_Gamma) then
    call Gam_basis2matrix
  endif

  return
END FUNCTION read_Gamma



SUBROUTINE write_GWGamma
  implicit none
  integer :: ix, iy, isite, ityp
  integer :: it1, it2
  integer :: ibasis, ibin
  character*26 ich


  open(100, status="replace", file=trim(title_loop)//"_G_file.dat")
  open(101, status="replace", file=trim(title_loop)//"_W_file.dat")
  open(102, status="replace", file=trim(title_loop)//"_Gam_file.dat")

  do it1 = 0, MxT-1
    do ityp = 1, NTypeG
      write(100, *) G(ityp, it1)
    enddo
  enddo

  do it1 = 0, MxT-1
    do isite = 0, Vol-1
      do ityp = 1, NTypeW
        write(101, *) W(ityp, isite, it1)
      enddo
    enddo
  enddo

  do ibasis = 1, NBasisGam
    do ibin = 1, NbinGam
      do isite = 0, Vol-1
        do ityp = 1, NTypeGam
          write(102, *) GamBasis(ityp, isite, ibin, ibasis)
        enddo
      enddo
    enddo
  enddo


  close(100)
  close(101)
  close(102)
  return
END SUBROUTINE write_GWGamma




SUBROUTINE read_Gamma_MC(changeBeta, mcBeta)
  IMPLICIT none
  double precision, intent(out) :: mcBeta
  logical, intent(out) :: changeBeta

  call LogFile%QuickLog("read monte carlo data from collapse")
  call read_monte_carlo_data(mcBeta)
  call LogFile%QuickLog("get the new Gamma function...")
  call Gam_mc2matrix_mc(changeBeta)
  call LogFile%QuickLog("read_Gamma_MC done!")
  return
END SUBROUTINE


!!================================================================
!!================= READ/WRITE CONFIGURATIONS ====================
!!================================================================
SUBROUTINE write_monte_carlo_data
  implicit none
  integer :: iorder, itopo, ir, irr, ityp, it1, it2
  double precision :: rgam2, rerr
  integer :: ibin, ibasis
  complex*16 :: gam1, normal

  ratioerr = 1.d0

  !=========== write into files =========================================
  open(104, status="replace", &
    & file=trim(title_mc)//"_monte_carlo_data.bin.dat",form="binary")

  write(104) finalBeta, Beta, MCOrder, L(1:D)
  write(104) imc, GamNorm, GamNormWeight
  write(104) Z_normal, ratioerr
  do it1 = 0, MxT-1
    do ir = 0, Vol-1
      do iorder = 0, MCOrder
        write(104)  GamMC(iorder, ir, it1)
        write(104)  ReGamSqMC(iorder, ir, it1)
        write(104)  ImGamSqMC(iorder, ir, it1)
      enddo
    enddo
  enddo

  do ibasis = 1, NBasisGam
    do ibin = 1, NbinGam
      do ir = 0, Vol-1
        do ityp = 1, NtypeGam/2
          do iorder = 0, MCOrder
            write(104) GamMCBasis(iorder, ityp, ir, ibin, ibasis)
            write(104) ReGamSqBasis(iorder, ityp, ir, ibin, ibasis)
            write(104) ImGamSqBasis(iorder, ityp, ir, ibin, ibasis)
          enddo
        enddo
      enddo
    enddo
  enddo
  
  close(104)

END SUBROUTINE write_monte_carlo_data


SUBROUTINE read_monte_carlo_data(mcBeta)
  implicit none
  double precision, intent(out) :: mcBeta

  integer :: iorder, ir, ityp, it1, it2, itopo,ios
  integer :: ibin, ibasis
  integer :: iL(1:D)
  double precision :: fBeta, mcsteps, oldZ_normal
  complex*16 :: oldGamNormWeight, oldGamNorm
  logical :: alive
  

  inquire(file=trim(title)//"_monte_carlo_data.bin.dat",exist=alive)
  if(.not. alive) then
    call LogFile%QuickLog("There is no monte carlo binary data yet!",'e')
    stop -1
  endif

  open(105, status="old", file=trim(title)//"_monte_carlo_data.bin.dat",form="binary")
  read(105,iostat=ios) fBeta, mcBeta, iorder, iL(1:D)
  read(105,iostat=ios) imc, GamNorm, GamNormWeight
  read(105,iostat=ios) Z_normal, ratioerr

  do it1 = 0, MxT-1
    do ir = 0, Vol-1
      do iorder = 0, MCOrder
        read(105,iostat=ios)  GamMC(iorder, ir, it1)
        read(105,iostat=ios)  ReGamSqMC(iorder,  ir, it1)
        read(105,iostat=ios)  ImGamSqMC(iorder,  ir, it1)
      enddo
    enddo
  enddo

  do ibasis = 1, NBasisGam
    do ibin = 1, NbinGam
      do ir = 0, Vol-1
        do ityp = 1, NtypeGam/2
          do iorder = 0, MCOrder
            read(105) GamMCBasis(iorder, ityp, ir, ibin, ibasis)
            read(105) ReGamSqBasis(iorder, ityp, ir, ibin, ibasis)
            read(105) ImGamSqBasis(iorder, ityp, ir, ibin, ibasis)
          enddo
        enddo
      enddo
    enddo
  enddo

  close(105)
  if(ios/=0) then
    call LogFile%QuickLog("The monte carlo binary data is broken?",'e')
    stop -1
  endif
END SUBROUTINE read_monte_carlo_data



  
!!================================================================
!!========== PRINT OUT THE DATA FILES ============================
!!================================================================

SUBROUTINE output_Quantities
  implicit none
  integer :: ityp, it1, it2, iorder, itt1, itt2
  integer :: dr(D), dx, dy, it, isite, site, ip, iomega
  complex*16 :: gam1
  double precision :: normal, ratio
  integer :: ibin, ibasis
  double precision :: tau1, tau2

  !open(104, access='append', file=trim(title_loop)//"_quantities.dat") 
  open(104, status='replace', file=trim(title_loop)//"_quantities.dat") 

  write(104, *) "##################################Gamma"
  write(104, *) "#tau1:", MxT, ",tau2:", MxT
  write(104, *) "#Beta", Beta, "L", L(1), "Order", MCOrder
  do it2 = 0, MxT-1
    do it1 = 0, MxT-1
      write(104, *)  real(Gam(1, 0, it1, it2)), dimag(Gam(1, 0, it1, it2))
    enddo
  enddo
  write(104, *)

  write(104, *) "##################################GammaR"
  write(104, *) "#r:", Vol
  write(104, *) "#Beta", Beta, "L", L(1), "Order", MCOrder
  do isite = 0, Vol-1
    write(104, *)  real(sum(Gam(1, isite, :, :))/(MxT)**2.d0),  &
      &  dimag(sum(Gam(1, isite, :, :))/(MxT)**2.d0)
  enddo
  write(104, *)

  normal = GamNormWeight/GamNorm
  do iorder = 1, MCOrder
    write(104, *) "##################################GammaDiag",trim(adjustl(str(iorder)))
    write(104, *) "#tau1:", MxT
    write(104, *) "#Beta", Beta, "L", L(1), "Order", MCOrder
    do it1 = 0, MxT-1
      write(104, *) real(GamMC(iorder,0,it1)*normal), dimag(GamMC(iorder,0,it1)*normal)
    enddo
    write(104, *)
  enddo

  do iorder = 1, MCOrder
    write(104, *) "##################################Gamma",trim(adjustl(str(iorder)))
    write(104, *) "#tau1:", MxT, ",tau2:", MxT
    write(104, *) "#Beta", Beta, "L", L(1), "Order", MCOrder
    do it2 = 0, MxT-1
      do it1 = 0, MxT-1
        gam1 = Gam_basis(it1, it2, GamMCBasis(iorder, 1, 0, :, :))
        write(104, *) real(gam1*normal), dimag(gam1*normal)
      enddo
    enddo
    write(104, *)
  enddo

  do iorder = 1, MCOrder
    write(104, *) "##################################GammaR",trim(adjustl(str(iorder)))
    write(104, *) "#r:", Vol
    write(104, *) "#Beta", Beta, "L", L(1), "Order", MCOrder
    do isite = 0, Vol-1
      it1 = MxT/2
      it2 = MxT/2 
      gam1 = Gam_basis(it1, it2, GamMCBasis(iorder, 1, isite, :, :))
      write(104, *) real(gam1*normal), dimag(gam1*normal)
    enddo
    write(104, *)
  enddo

  write(104, *) "##################################G"
  write(104, *) "#tau:", MxT
  write(104, *) "#Beta", Beta, "L", L(1), "Order", MCOrder
  do it1 = 0, MxT-1
    write(104, *)  real(G(1, it1)), dimag(G(1,it1))
  enddo
  write(104, *)

  write(104, *) "##################################G0"
  write(104, *) "#tau:", MxT
  write(104, *) "#Beta", Beta, "L", L(1), "Order", MCOrder
  do it1 = 0, MxT-1
    write(104, *)  real(G0F(it1)), dimag(G0F(it1))
  enddo
  write(104, *)

  write(104, *) "##################################W"
  write(104, *) "#r:", Vol, ",tau:", MxT
  write(104, *) "#Beta", Beta, "L", L(1), "Order", MCOrder
  do it1 = 0, MxT-1
    do isite = 0, Vol-1
      write(104, *)  real(W(1, isite, it1)), dimag(W(1, isite, it1))
    enddo
  enddo
  write(104, *)

  write(104, *) "##################################Chi"
  write(104, *) "#r:", Vol
  write(104, *) "#Beta", Beta, "L", L(1), "Order", MCOrder

  ratio = 1.d0/dble(MxT)
  do isite = 0, Vol-1
    write(104, *) ratio*real(SUM(Chi(isite, :))),ratio*dimag(SUM(Chi(isite, :)))
  enddo
  
  write(104, *) "##################################Polar"
  write(104, *) "#tau:", MxT
  write(104, *) "#Beta", Beta, "L", L(1), "Order", MCOrder
  do it = 0, MxT-1
    write(104, *) real(Polar(0, it)), dimag(Polar(0, it)) 
  enddo

  write(104, *) "##################################Sigma"
  write(104, *) "#tau:", MxT
  write(104, *) "#Beta", Beta, "L", L(1), "Order", MCOrder
  do it = 0, MxT-1
    write(104, *) Vol*(MxT/Beta)**2.d0*real(Sigma(it)),  &
      & Vol*(MxT/Beta)**2.d0*dimag(Sigma(it))
  enddo

  write(104, *) "##################################SUMChi"
  write(104, *) "#tau:", MxT
  write(104, *) "#Beta", Beta, "L", L(1), "Order", MCOrder
  do it = 0, MxT-1
    write(104, *) real(SUM(Chi(:, it))),dimag(SUM(Chi(:, it)))
  enddo

  write(104, *) "##################################Denom"
  write(104, *) "#p:", Vol, ",omega:", MxT
  write(104, *) "#Beta", Beta, "L", L(1), "Order", MCOrder
  do iomega = 0, MxT-1
    do ip = 0, Vol-1
      write(104, *) real(Denom(ip, iomega)),dimag(Denom(ip, iomega))
    enddo
  enddo

  call transfer_Chi_r(1)
  write(104, *) "##################################ChiK"
  write(104, *) "#k:", Vol
  write(104, *) "#Beta", Beta, "L", L(1), "Order", MCOrder
  ratio = 1.d0/dble(MxT)
  do isite = 0, Vol-1
    write(104, *) ratio*real(SUM(Chi(isite, :))),ratio*dimag(SUM(Chi(isite, :)))
  enddo

  write(104, *) "##################################ChiKt0"
  write(104, *) "#k:", Vol
  write(104, *) "#Beta", Beta, "L", L(1), "Order", MCOrder
  do isite = 0, Vol-1
    write(104, *) real(Chi(isite, 0)), dimag(Chi(isite, 0))
  enddo

  call transfer_Chi_r(-1)

  close(104)

  call output_denominator
END SUBROUTINE output_Quantities

SUBROUTINE output_denominator
  implicit none
  integer :: ip

  ip = get_site_from_cord(D, L(1:D)/2)
  open(104, access='append', file=trim(title_loop)//"_denom.dat") 
  write(104, *) real(Denom(ip, 0)),dimag(Denom(ip, 0))
  close(104)
END SUBROUTINE output_denominator

SUBROUTINE output_Gam1
  implicit none
  integer :: ityp, it1, it2, ibin, ibasis
  integer :: dx, dy, it
  complex*16 :: gam1
  double precision :: tau1, tau2

  open(104, status='replace', file=trim(title_loop)//"_Gam1.dat")

  write(104, *) "##################################Gamma"
  write(104, *) "#tau1:", MxT, ",tau2:", MxT
  write(104, *) "#Beta", Beta, "L", L(1), L(2), "Order", MCOrder
  do it2 = 0, MxT-1
    do it1 = 0, MxT-1
      write(104, *)  real(GamOrder1(1, it1, it2)), dimag(GamOrder1(1,it1,it2))
    enddo
  enddo
  write(104, *)

  !write(104, *) "##################################GammaBasis"
  !write(104, *) "#tau1:", MxT, ",tau2:", MxT
  !write(104, *) "#Beta", Beta, "L", L(1), L(2), "Order", MCOrder
  !do it2 = 0, MxT-1
    !do it1 = 0, MxT-1

      !ibin = get_bin_Gam(it1, it2)

      !tau1 = dble(it1)*Beta/dble(MxT)
      !tau2 = dble(it2)*Beta/dble(MxT)

      !if(IsBasis2D(ibin)) then
        !gam1 = (0.d0, 0.d0)
        !do  ibasis = 1, NBasisGam
          !gam1 = gam1 + GamMCBasis(1, 1, 0, ibin, ibasis)* weight_basis_Gam( &
            !& CoefGam(0:BasisOrderGam,0:BasisOrderGam,ibasis,ibin), tau1, tau2)
        !enddo
      !else
        !gam1 = (0.d0, 0.d0)
        !do  ibasis = 1, NBasis
          !gam1 = gam1 + GamMCBasis(1, 1, 0, ibin, ibasis)* weight_basis( &
            !& CoefGam(0:BasisOrder,0,ibasis,ibin), tau1)
        !enddo
      !endif

      !write(104, *) real(gam1), dimag(gam1)
    !enddo
  !enddo
  !write(104, *)
  close(104)
  return
END SUBROUTINE output_Gam1

  
!!================================================================
!!========== PRINT OUT THE DATA FILES ============================
!!================================================================

SUBROUTINE output_test
  implicit none
  integer :: iorder
  !open(104, status='replace', file=trim(title_mc)//"_Gam_Order_test.dat")
  open(104, access='append', file=trim(title_mc)//"_Gam_Order_test.dat")
  do iorder = 1, MCOrder
    write(104, *) iorder, Quan(iorder+2)/Quan(2), Error(iorder+2)/Quan(2)
  enddo
  write(104, *)

  write(104, *) "Type 1,2 <==> 3,4"
  do iorder = 0, MCOrder
    BalenceCheck(iorder,1,3)=(BalenceCheck(iorder,1,1)-BalenceCheck(iorder,1,2)) &
      & /sqrt(BalenceCheck(iorder,1,1))
    write(104, '(i3,3f17.5)') iorder, BalenceCheck(iorder,1,:)
  enddo
  write(104, *)
  close(104)
END SUBROUTINE

!!================================================================
!!================================================================
!!================================================================



!!================================================================
!!================================================================
!!================================================================


SUBROUTINE write_monte_carlo_conf
  implicit none
  integer :: i

  open(103, status="replace", &
   & file=trim(title_mc)//"_monte_carlo_conf.bin.dat",form="binary")

  write(103)  mc_version, Z_normal, Z_worm, StatNum
  write(103)  ProbCall(:,:), ProbProp(:,:), ProbAcc(:,:)
  write(103)  GamOrder(:), GamWormOrder(:)
  write(103)  CoefOfWorm, CoefOfWeight(:)
  write(103)  TestData(:)
  write(103)  TailLn, TailVertex
  write(103)  NextLn(:)
  write(103)  NextVertex(:)
  write(103)  Order, MeasureGam, SignFermiLoop, IsWormPresent 
  write(103)  Ira, Masha, SpinMasha, kMasha
  write(103)  StatusVertex(:)
  write(103)  TypeVertex(:)
  write(103)  GRVertex(:)
  write(103)  WRVertex(:)
  write(103)  TVertex(1, :), TVertex(2, :), TVertex(3, :)
  write(103)  DirecVertex(:)             
  write(103)  IsDeltaVertex(:)
  write(103)  IsDeltaLn(:)
  write(103)  kLn(:)
  write(103)  KindLn(:)

  do i = 1, 2
    write(103)  NeighLn(i,:)
  enddo

  do i = 1, 3
    write(103)  NeighVertex(i,:) 
  enddo

  close(103)
END SUBROUTINE write_monte_carlo_conf



SUBROUTINE read_monte_carlo_conf
  implicit none
  integer :: i,ios
  integer :: ikey, ikeyG, ikeyW
  integer :: iGin, iGout, iW
  integer :: iGam, jGam
  complex*16 :: ComCurrent
  double precision :: tau
  logical :: alive

  inquire(file=trim(title)//"_monte_carlo_conf.bin.dat",exist=alive)
  if(.not. alive) then
    call LogFile%QuickLog("There is no monte carlo binary configuration yet!",'e')
    stop -1
  endif
  open(106, status="old", file=trim(title)//"_monte_carlo_conf.bin.dat",form="binary")

  read(106,iostat=ios)  mc_version, Z_normal, Z_worm, StatNum
  read(106,iostat=ios)  ProbCall(:,:), ProbProp(:,:), ProbAcc(:,:)
  read(106,iostat=ios)  GamOrder(:), GamWormOrder(:)
  read(106,iostat=ios)  CoefOfWorm, CoefOfWeight(:)
  read(106,iostat=ios)  TestData(:)
  read(106,iostat=ios)  TailLn, TailVertex
  read(106,iostat=ios)  NextLn(:)
  read(106,iostat=ios)  NextVertex(:)
  read(106,iostat=ios)  Order, MeasureGam, SignFermiLoop, IsWormPresent 
  read(106,iostat=ios)  Ira, Masha, SpinMasha, kMasha
  read(106,iostat=ios)  StatusVertex(:)
  read(106,iostat=ios)  TypeVertex(:)
  read(106,iostat=ios)  GRVertex(:)
  read(106,iostat=ios)  WRVertex(:)
  read(106,iostat=ios)  TVertex(1, :), TVertex(2, :), TVertex(3, :)
  read(106,iostat=ios)  DirecVertex(:)             
  read(106,iostat=ios)  IsDeltaVertex(:)
  read(106,iostat=ios)  IsDeltaLn(:)
  read(106,iostat=ios)  kLn(:)
  read(106,iostat=ios)  KindLn(:)

  do i = 1, 2
    read(106,iostat=ios)  NeighLn(i,:)
  enddo

  do i = 1, 3
    read(106,iostat=ios)  NeighVertex(i,:) 
  enddo
  close(106)

  if(ios/=0) then
    call LogFile%QuickLog("The monte carlo binary configuration is broken?",'e')
    stop -1
  endif

  NVertex = 2*(Order+1)
  NGLn = NVertex
  NWLn = Order+1
  ComCurrent = (1.d0, 0.d0)

  ikey = 1
  do i = 1, MxNVertex
    if(StatusVertex(i)==-1) cycle
    VertexValue2Key(i) = ikey
    VertexKey2Value(ikey) = i
    ikey = ikey + 1
  enddo

  if(ikey/=NVertex+1) then
    call LogFile%QuickLog("read_monte_carlo_conf: Number of Vertex Error!",'e')
    stop
  endif

  StatusLn(:)=-1
  do ikey = 1, NVertex
    i = VertexKey2Value(ikey)
    if(StatusVertex(i)==-1) then
      call LogFile%QuickLog("read_monte_carlo_conf: Status of Vertex Error!",'e')
      stop
    endif

    iGin = NeighVertex(1, i)
    iGout = NeighVertex(2, i)
    iW = NeighVertex(3, i)

    if(StatusLn(iGin)==-1) then
      StatusLn(iGin) = gline_stat(StatusVertex(NeighLn(1, iGin)), StatusVertex(NeighLn(2, iGin)))
      KindLn(iGin) = 1
    endif
    if(StatusLn(iGout)==-1) then
      StatusLn(iGout) = gline_stat(StatusVertex(NeighLn(1, iGout)), StatusVertex(NeighLn(2, iGout)))
      KindLn(iGout) = 1
    endif
    if(StatusLn(iW)==-1) then
      StatusLn(iW) = wline_stat(StatusVertex(NeighLn(1, iW)), StatusVertex(NeighLn(2, iW)))
      KindLn(iW) = 2
    endif

    if(TypeVertex(i)==1 .or. TypeVertex(i)==4) then
      SpInVertex(:, i) = 1
    else if(TypeVertex(i)==2 .or. TypeVertex(i)==3) then
      SpInVertex(:, i) = 2
    else if(TypeVertex(i)==5 .or. TypeVertex(i)==6) then
      SpInVertex(2, i) = Mod(TypeVertex(i), 2)+1
      SpInVertex(1, i) = 3 - SpInVertex(2, i)
    endif

    WeightVertex(i) = weight_vertex(StatusVertex(i), IsDeltaVertex(i), &
      & diff_r(D, GRVertex(i),WRVertex(i)), &
      & TVertex(3,i)-TVertex(2,i), TVertex(1,i)-TVertex(3,i), TypeVertex(i))

    ComCurrent = ComCurrent* WeightVertex(i)
  enddo

  ikeyG = 1
  ikeyW = 1
  do i = 1, MxNLn
    if(StatusLn(i)==-1) cycle
    if(KindLn(i)==1) then
      LnValue2Key(i) = ikeyG
      GLnKey2Value(ikeyG) = i
      ikeyG = ikeyG + 1
    else if(KindLn(i)==2) then
      LnValue2Key(i) = ikeyW
      WLnKey2Value(ikeyW) = i
      ikeyW = ikeyW + 1
    endif
  enddo

  if(ikeyG/=NGLn+1) then
    call LogFile%QuickLog("read_monte_carlo_conf: Number of Glines Error!",'e')
    stop
  endif

  if(ikeyW/=NWLn+1) then
    call LogFile%QuickLog("read_monte_carlo_conf: Number of Glines Error!",'e')
    stop
  endif

  do ikeyG = 1, NGLn
    i = GLnKey2Value(ikeyG)
    if(Is_k_valid_for_G(kLn(i))) then
      call add_Hash4G(kLn(i),i)
    else
      call LogFile%QuickLog("read_monte_carlo_conf: k of G Error!",'e')
      stop
    endif
    TypeLn(i) = mod(TypeVertex(NeighLn(2,i)), 2)
    if(TypeLn(i)==0)  TypeLn(i) = 2
    tau = TVertex(2, NeighLn(2, i))-TVertex(1,NeighLn(1,i))
    WeightLn(i) = weight_gline(StatusLn(i), tau, TypeLn(i))
    ComCurrent = ComCurrent* WeightLn(i)
  enddo

  do ikeyW = 1, NWLn
    i = WLnKey2Value(ikeyW)
    if(Is_k_valid_for_W(kLn(i))) then
      call add_Hash4W(kLn(i),i)
    else
      call LogFile%QuickLog("read_monte_carlo_conf: k of W Error!",'e')
      stop
    endif
    iGam = NeighLn(1, i)
    jGam = NeighLn(2, i)
    if(StatusLn(i)<=1) then
      TypeLn(i) = TypeGam2W(TypeVertex(iGam), TypeVertex(jGam))
    else
      TypeLn(i) = 1
    endif

    tau = TVertex(3, iGam)-TVertex(3, jGam)
    WeightLn(i) = weight_wline(StatusLn(i), IsDeltaLn(i), diff_r(D, WRVertex(jGam),WRVertex(iGam)), &
      & tau, TypeLn(i))
    ComCurrent = ComCurrent* WeightLn(i)
  enddo

  WeightCurrent = abs(ComCurrent) *(1.d0/Beta)**Order *SignFermiLoop
  Phase = ComCurrent/abs(ComCurrent) *SignFermiLoop

  return
END SUBROUTINE read_monte_carlo_conf
!!================================================================
!!================================================================
!!================================================================


!================================================================
!================= print/draw information    ====================
!================================================================

!==================================================================
!===================== PRINT STATUS ================================
!==================================================================
!If you want to log staff when runing markov, log here!
SUBROUTINE print_status
    implicit none
    integer :: iorder,i
    character*30 :: updatename(Nupdate)

    call LogFile%WriteStamp()
    call LogFile%WriteLine("MC steps:"+str(imc))
    call time_elapse
    call LogFile%WriteLine("Printing interval:"+str(t_elap,'(f12.3)')+'s')
    call LogFile%WriteLine("Efficiency: "+str(imc/t_elap,'(f12.0)')+"steps per second.")
    call LogFile%WriteLine('Statistics Number ='+str(StatNum))

    do i=0,MCOrder
      if(Norm(i)>1e-6) then
        write(logstr,"(i2, A,f15.6,'+/-',f15.6)") i,QuanName(i),Quan(i)/Norm(i),Error(i) 
        call LogFile%WriteLine(logstr)
      endif
    enddo

    call LogFile%WriteLine("------------------------------------------------")

    updatename(1)= " 1: create worm along wline"
    updatename(2)= " 2: delete worm along wline"
    updatename(3)= " 3: create worm along gline"
    updatename(4)= " 4: delete worm along gline"
    updatename(5)= " 5: move worm along wline"
    updatename(6)= " 6: move worm along gline"
    updatename(7)= " 7: add interaction"
    updatename(8)= " 8: remove interaction"
    updatename(9)= " 9: add interaction cross"
    updatename(10)= "10: remove interaction cross"
    updatename(11)= "11: reconnect"
    updatename(12)= "12: shift gline in space"
    updatename(13)= "13: shift wline in space"
    updatename(14)= "14: change Gamma type"
    updatename(15)= "15: move measuring index"
    updatename(16)= "16: change Gamma time"
    updatename(17)= "17: change wline isdelta"
    updatename(18)= "18: change Gamma isdelta"

    do iorder = 0, MCOrder
      call LogFile%WriteLine("Order"+str(iorder))
      do i = 1, Nupdate
        if(ProbProp(iorder, i)/=0.d0) then
          write(logstr, '(A,2f13.1,f13.5,f13.1,f13.5)') updatename(i), ProbCall(iorder, i), ProbProp(iorder, i), &
            & ProbProp(iorder, i)/ProbCall(iorder, i), ProbAcc(iorder, i), &
            & ProbAcc(iorder, i)/ProbProp(iorder, i)
          call LogFile%WriteLine(logstr)
        endif
      enddo
    enddo

    call LogFile%WriteLine("------------------------------------------------")
    do iorder = 0, MCOrder
      call LogFile%WriteLine("Order"+str(iorder))
      BalenceCheck(iorder,1,3)=(BalenceCheck(iorder,1,1)-BalenceCheck(iorder,1,2)) &
        & /sqrt(BalenceCheck(iorder,1,1))
      write(logstr, '(A,3f17.5)') "Gamma Type 1,2 <==> 3,4:",BalenceCheck(iorder,1,:)
      call LogFile%WriteLine(logstr)
    enddo
    call LogFile%WriteLine("------------------------------------------------")
    call LogFile%WriteLine("Reducibility ratio "+str(TestData(2)/TestData(1))+" : "+str(TestData(3)/TestData(1)))

END SUBROUTINE print_status

!====================================================================
!===================== PRINT CONFIGURATION ==========================
!====================================================================


SUBROUTINE print_config
  implicit none
  integer :: i, iln, iv
  
  open(108, access='append', file=trim(title_mc)//"_mc.conf")
  
  write(108, *) "============================================================"
  write(108, *) imc, IsWormPresent, iupdate

  if(IsWormPresent .eqv. .true.) then
    write(108, *) "Ira", Ira, "Masha", Masha, "SpinMasha", SpinMasha
    write(108, *) "kMasha", kMasha
  endif

  write(108, *) "Order", Order
  write(108, *) "SignFermiLoop", SignFermiLoop

  write(108, *) "Measuring Gamma", MeasureGam
  write(108, *) "Phase", Phase
  write(108, *) "Weight", WeightCurrent

  do i = 1, NGLn
    iln = GLnKey2Value(i)
    if(StatusLn(iln) <0) cycle
    write(108, 10) iln, KindLn(iln), IsDeltaLn(iln), TypeLn(iln), kLn(iln), StatusLn(iln), NeighLn(1:2,iln)
  enddo

  do i = 1, NWLn
    iln = WLnKey2Value(i)
    if(StatusLn(iln) <0) cycle
    write(108, 10) iln, KindLn(iln), IsDeltaLn(iln), TypeLn(iln), kLn(iln), StatusLn(iln), NeighLn(1:2,iln)
  enddo

  do i = 1, NVertex
    iv = VertexKey2Value(i)
    if(StatusVertex(iv) <0) cycle
    write(108, 12) iv,IsDeltaVertex(iv), TypeVertex(iv), SpInVertex(1, iv),SpInVertex(2, iv), &
      & GRVertex(iv),WRVertex(iv), TVertex(1, iv), TVertex(2, iv),  &
      & TVertex(3, iv), DirecVertex(iv), StatusVertex(iv), NeighVertex(:,iv)
  enddo
  write(108, *) "============================================================"

  10 format(' Line:',i2,2x,'kind:',i2,2x,'isdelta:',i2,2x,'type:',i2,2x,'k:',i8,2x,'stat:',i2, 2x,&
    & 'neigh:',i6,i6)
  12 format('Gamma:',i2,2x,'isdelta:',i2,2x,'type:',i2,2x,'typein:',i2,2x,'typeout:',i2,2x,&
    & 'gr:(',i4,'), wr:(',i4,')', 't:(', f7.4, f7.4, f7.4, ')',2x, &
    & 'direction:', i2,2x, 'stat:',i2, 2x,'neigh:', i6,i6,i6)

  close(108)
  call DRAW
END SUBROUTINE print_config

!=================== VISUALIZATION  ==================================
SUBROUTINE DRAW
    IMPLICIT NONE
    DOUBLE PRECISION :: x1,y1, x2,y2, y3, scx, scy, sgn   
    DOUBLE PRECISION :: scydash, ca1,ca2,ra,a1,a2, radian
    DOUBLE PRECISION :: phi1, phi2, pi2, pi4,ini,seg,theta2
    INTEGER :: i
    INTEGER :: FEXIST, RES
    integer :: iWLn,iGLn,Vertex1,Vertex2,Vertex3
    character(10) :: imcstr
    character(100) :: tempstr
    
    pi2=dasin(1.d0)
    pi4=pi2/2.d0
    theta2=dasin(0.2d0)
    !theta2=pi4
    radian=90.d0/pi2
    scx=500/beta
    scy=400./Vol
    x1=scx*beta
    y1=scy*Vol
    scydash=scy/40.

    INQUIRE(DIRECTORY="graph",EXIST=FEXIST)
    IF(.not.FEXIST) THEN
      call system("mkdir graph")
    ENDIF

    write(imcstr,'(i10)') int(imc)
    open(11, file='graph/graph_'//trim(adjustl(title2))//'_'//trim(adjustl(imcstr))//'.eps')
    write(11,*) '%!'
    write(11,*) '%BoundingBox: 0 0 ', x1, y1
    write(11,*) '%%EndComments'
    write(11,*) '%%BeginProlog'
    write(11,*) '/L { lineto stroke} def'
    write(11,*) '/M { moveto } def'
    write(11,*) '/N {newpath } def'
    write(11,*) '/Ndashed {[5 5] 0 setdash newpath } def'
    write(11,*) '/Nsolid {[] 0 setdash newpath } def'
    write(11,*) '/Y { 0 360 arc closepath gsave fill'
    write(11,*) '     grestore stroke } def'
    write(11,*) '/YL { 0 360 arc stroke} def' 
    write(11,*) '/Y45 { arc stroke} def'
    write(11,*) '/C {/Times-Roman findfont 15 scalefont setfont show} def'
    write(11,*) '% Put an arrowhead at point x2 y2,'
    write(11,*) '% pointing away from x1 y1'
    write(11,*) '% Replace x2 y2 with coordinates of arrowbase:'
    write(11,*) '% the point to connect lines to'
    write(11,*) '% ArrowHeadSize gives the size of the arrow'
    write(11,*) '/ArrowHeadSize 10 def'
    write(11,*) '/ahead {'
    write(11,*) '    1 index 4 index sub'
    write(11,*) '    1 index 4 index sub'
    write(11,*) '    exch atan'
    write(11,*) '    ArrowHeadSize -.8 mul'
    write(11,*) '    dup'
    write(11,*) '    2 index cos mul 4 index add'
    write(11,*) '    exch'
    write(11,*) '    2 index sin mul 3 index add'
    write(11,*) '    5 2 roll'
    write(11,*) '    gsave'
    write(11,*) '        3 1 roll'
    write(11,*) '        translate'
    write(11,*) '        rotate'
    write(11,*) '        newpath'
    write(11,*) '        0 0 moveto'
    write(11,*) '        ArrowHeadSize dup neg exch .25 mul'
    write(11,*) '        2 copy lineto'
    write(11,*) '        ArrowHeadSize -.8 mul 0'
    write(11,*) '        2 copy'
    write(11,*) '        6 4 roll'
    write(11,*) '        neg curveto'
    write(11,*) '        closepath fill'
    write(11,*) '    grestore'
    write(11,*) '} bind def'
    write(11,*) ''
    write(11,*) '%%EndProlog'
    write(11,*) '%%BeginSetup'
    write(11,*) '1 setlinewidth'
    write(11,*) '5 140 translate'
    write(11,*) '1 1 scale'
    write(11,*) '%%EndSetup'

    ini=640.0
    seg=15.0
    write(11,*) '0 0 0 setrgbcolor'
    write(tempstr, *) "Beta: ",beta,"    L: ",L(1:D)
    write(11,803) 0.,ini,tempstr
    write(tempstr,*)  "Jcp: ",Jcp," Seed: ",Seed
    ini=ini-seg
    write(11,803) 0.,ini,tempstr
    write(tempstr, *) "imc: ",imc," Is Worm Here:", IsWormPresent
    ini=ini-seg
    write(11,803) 0.,ini,tempstr
    if(IsWormPresent .eqv. .true.) then
      write(tempstr, *) "Ira: ", Ira, "Masha: ", Masha
      ini=ini-seg
      write(11,803) 0.,ini,tempstr
      write(tempstr, *) "Spin of Masha: ", SpinMasha, "k of Masha: ", kMasha
      ini=ini-seg
      write(11,803) 0.,ini,tempstr
    endif

    write(tempstr, *) "Order: ", Order
    ini=ini-seg
    write(11,803) 0.,ini,tempstr
    write(tempstr, *) "Fermi Loop Sign: ", SignFermiLoop
    ini=ini-seg
    write(11,803) 0.,ini,tempstr

    write(tempstr, *) "Measuring Gamma: ", MeasureGam
    ini=ini-seg
    write(11,803) 0.,ini,tempstr
    write(tempstr, *) "Phase: ", Phase
    ini=ini-seg
    write(11,803) 0.,ini,tempstr
    write(tempstr, *) "Weight: ", WeightCurrent
    ini=ini-seg
    write(11,803) 0.,ini,tempstr

    ini=455.0
    write(11,*) '0 1 0 setrgbcolor'
    write(11,777) 500.0, ini, scy/20.
    write(11,803) 520.0, ini-5.0, "MeasureGam"
    ini=ini-seg
    write(11,*) '1 0 0 setrgbcolor'
    write(11,777) 500.0, ini, scy/30.
    write(11,803) 520.0, ini-5.0, "Ira"
    ini=ini-seg
    write(11,*) '0 0 1 setrgbcolor'
    write(11,777) 500.0, ini, scy/30.
    write(11,803) 520.0, ini-5.0, "Masha"

    write(11,*) '0 0 0 setrgbcolor'
    do i=1,NWLn;
      iWLn=WLnKey2Value(i)
      Vertex1=NeighLn(1,iWLn)
      x1=scx*TVertex(3, Vertex1)
      y1=scy*GRVertex(Vertex1)
      Vertex2=NeighLn(2,iWLn)
      x2=scx*TVertex(3, Vertex2)
      y2=scy*GRVertex(Vertex2)

      ra=dsqrt((x2-x1)**2+(y2-y1)**2)/2.d0
      if(dabs(ra)<1e-6) then
        write(11,792) x1, y1+scy/5. , scy/5. 
        cycle
      endif
      phi1=(y2-y1)/(2.d0*ra)
      phi2=(x2-x1)/(2.d0*ra)
      ca1=(x1+x2)/2.d0
      ca2=(y1+y2)/2.d0;
      ra=ra/sin(theta2)*cos(theta2)
      ca1=ca1+phi1*ra
      ca2=ca2-phi2*ra; 
      ra=ra/cos(theta2)
      a1=pi2-theta2+dasin(phi1)
      a2=a1+2.0*theta2
      IF(phi2.lt.0) then
        a2=4*pi2-pi2+theta2-dasin(phi1)
        a1=a2-2.0*theta2
      endif
      a1=a1*radian
      a2=a2*radian 

      write(11,791)  ca1, ca2, ra, a1, a2
    enddo
  
    do i=1,NGLn;
      iGLn=GLnKey2Value(i)
      Vertex1=NeighLn(1,iGLn)
      x1=scx*TVertex(3, Vertex1)
      y1=scy*GRVertex(Vertex1)
      Vertex2=NeighLn(2,iGLn)
      x2=scx*TVertex(3, Vertex2)
      y2=scy*GRVertex(Vertex2)

      if(TypeLn(iGLn)==1) then
        write(11,*) '1 0 0 setrgbcolor'
      else
        write(11,*) '0 0 1 setrgbcolor'
      endif


      if(Vertex1/=Vertex2) then

        ra=dsqrt((x2-x1)**2+(y2-y1)**2)/2.d0
        if(dabs(ra)<1e-6) then
          iWLn=NeighVertex(3,Vertex1)
          Vertex3=NeighLn(3-DirecVertex(Vertex1),iWLn)
          y3=scy*GRVertex(Vertex3)
          sgn=1.d0
          IF(y3>y1) sgn=-1.d0   
          write(11,780) x1, y1+sgn*scy/5. , scy/5. 
          cycle
        endif
        phi1=(y2-y1)/(2.d0*ra)
        phi2=(x2-x1)/(2.d0*ra)
        ca1=(x1+x2)/2.d0
        ca2=(y1+y2)/2.d0;
        ca1=ca1+phi1*ra
        ca2=ca2-phi2*ra; 
        ra=ra*dsqrt(2.d0)
        a1=pi4+dasin(phi1)
        a2=a1+pi2
        IF(phi2.lt.0) then
          a2=4*pi2-pi4-dasin(phi1)
          a1=a2-pi2
        endif
        a1=a1*radian
        a2=a2*radian 
        write(11,781)  ca1, ca2, ra, a1, a2  ! propagator lines - arcs 

        ca1=x1+(phi2-phi1)*scy/2.;     
        ca2=y1+(phi2+phi1)*scy/2. 

        write(11,790) ca1, ca2, x1, y1       ! arrows
      else
        iWLn=NeighVertex(3,Vertex1)
        Vertex3=NeighLn(3-DirecVertex(Vertex1),iWLn)
        y3=scy*GRVertex(Vertex3)
        sgn=1.d0
        IF(y3>y1) sgn=-1.d0   
        write(11,780) x1, y1+sgn*scy/5. , scy/5. 
      endif
    enddo

    write(11,*) '0 0 0 setrgbcolor'
    write(11,"(f6.1,x,f6.1,x,' M (Gamma info) C')") 520.0,350.0
    ini=350.0-seg
    ! write Gamma info
    do i=1,NVertex
      Vertex1=VertexKey2Value(i)
      x1=scx*TVertex(3, Vertex1)
      y1=scy*GRVertex(Vertex1)
      if(Vertex1==MeasureGam) then
        write(11,*) '0 1 0 setrgbcolor'
        write(11,777) x1, y1, scy/20.
      endif
      if(Vertex1==Ira .and. IsWormPresent) then
        write(11,*) '1 0 0 setrgbcolor'
      elseif(Vertex1==Masha .and. IsWormPresent) then
        write(11,*) '0 0 1 setrgbcolor'
      else
        write(11,*) '0 0 0 setrgbcolor'
      endif
      write(11,777) x1, y1, scy/30.
      write(11,801) x1-5., y1+7., Vertex1
      if(Vertex1==MeasureGam) then
        write(11,*) '0 1 0 setrgbcolor'
      endif
      write(11,802) 500.0,ini,Vertex1,TVertex(3, Vertex1), &
         & GRVertex(Vertex1)
      ini=ini-seg
    enddo

    write(11,*) '0 0 0 setrgbcolor'
    write(11,"(f6.1,x,f6.1,x,' M (G info) C')") 530.0,ini
    ini=ini-seg
    do i=1,NGLn
      iGLn=GLnKey2Value(i)
      Vertex1=NeighLn(1,iGLn)
      Vertex2=NeighLn(2,iGLn)
      if(TypeLn(iGLn)==1) then
        write(11,*) '1 0 0 setrgbcolor'
      else
        write(11,*) '0 0 1 setrgbcolor'
      endif
      write(11,804) 500.0, ini, iGLn, Vertex2,Vertex1
      ini=ini-seg
    enddo
    write(11,*) '0 0 0 setrgbcolor'
    write(11,"(f6.1,x,f6.1,x,' M (W info) C')") 530.0,ini
    ini=ini-seg
    do i=1,NWLn;
      iWLn=WLnKey2Value(i)
      Vertex1=NeighLn(1,iWLn)
      Vertex2=NeighLn(2,iWLn)
      write(11,805) 500.0, ini, iWLn, Vertex1,Vertex2
      ini=ini-seg
    enddo


    write(11,*) ''
    write(11,*) 'stroke showpage'
    write(11,*) '%%Trailer'

    close (11)
  return
!******************************************************************
 777  format ('Nsolid ',f6.1,x,f6.1,x,f9.3,' Y')
 778  format ('N ',f6.1,x,f6.1,x,' M')
 779  format (     f6.1,x,f6.1,x,' L')
 780  format ('Nsolid ',f6.1,x,f6.1,x,f9.3,' YL')
 781  format ('Nsolid ',f6.1,x,f6.1,x,f6.1,x,f6.1,x,f6.1,' Y45')
 790  format ('Nsolid ',f6.1,x,f6.1,x,f6.1,x,f6.1,x,' ahead')
 !791  format ('Ndashed ',f6.1,x,f6.1,' M',f6.1,x,f6.1,' L')
 791  format ('Ndashed ',f6.1,x,f6.1,x,f6.1,x,f6.1,x,f6.1,' Y45')
 792  format ('Ndashed ',f6.1,x,f6.1,x,f9.3,' YL')
 801  format (f6.1,x,f6.1,x,' M (',i2,') C')
 802  format (f6.1,x,f6.1,x,' M (',i2,':',f6.3,'; (',i3,')) C')
 803  format (f6.1,x,f6.1,x,' M (',A,') C')
 804  format (f6.1,x,f6.1,x,' M (G ',i2,' : ',i2,' --->',i2') C')
 805  format (f6.1,x,f6.1,x,' M (W ',i2,' : ',i2,' <=>',i2') C')
 
END SUBROUTINE DRAW

!__________________________________________ 
!   	'0 0 0 setrgbcolor'    black
!     '0 0 1 setrgbcolor'    blue
!     '1 0 0 setrgbcolor'    red
!     '0 1 0 setrgbcolor'    green

!====================================================================
!====================================================================
