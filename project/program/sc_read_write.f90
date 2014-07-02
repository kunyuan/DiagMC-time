!!================================================================
!!================= READ/WRITE CONFIGURATIONS ====================
!!================================================================


SUBROUTINE read_GW
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
    else 
      G = Gtmp;         deallocate(Gtmp)
      W = Wtmp;         deallocate(Wtmp)
      call read_Gamma
      call update_WeightCurrent
      mc_version = file_version
    endif
  else 
    if(ios/=0) then
      call LogFile%QuickLog("Failed to read G,W information!",'e')
      close(100)
      close(101)
      stop -1
    endif
  endif

  close(100)
  close(101)
  return
END SUBROUTINE read_GW

SUBROUTINE write_GW
  implicit none
  integer :: ix, iy, isite, ityp
  integer :: it1, it2
  character*26 ich

  open(100, status="replace", file=trim(title_loop)//"_G_file.dat")
  open(101, status="replace", file=trim(title_loop)//"_W_file.dat")

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

  close(100)
  close(101)
  return
END SUBROUTINE write_GW



  
!!================================================================
!!========== PRINT OUT THE DATA FILES ============================
!!================================================================

SUBROUTINE output_Quantities
  implicit none
  integer :: ityp, it1, it2, iorder, itt1, itt2
  integer :: dr(D), dx, dy, it, isite, ip, iomega
  complex*16 :: gam1
  double precision :: normal, ratio
  integer :: ibin, ibasis
  double precision :: tau1, tau2

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
      write(104, *) real(GamMC(iorder,0,it1))*normal, dimag(GamMC(iorder,0,it1))*normal
    enddo
    write(104, *)
  enddo

  do iorder = 1, MCOrder
    write(104, *) "##################################Gamma",trim(adjustl(str(iorder)))
    write(104, *) "#tau1:", MxT, ",tau2:", MxT
    write(104, *) "#Beta", Beta, "L", L(1), "Order", MCOrder
    do it2 = 0, MxT-1
      do it1 = 0, MxT-1
        ibin = get_bin_Gam(it1, it2)

        tau1 = dble(it1)*Beta/dble(MxT)
        tau2 = dble(it2)*Beta/dble(MxT)

        gam1 = (0.d0, 0.d0)
        if(IsBasis2D(ibin)) then
          do ibasis = 1, NBasisGam
            gam1 = gam1 + GamBasis(iorder, 1, 0, ibin, ibasis)*  &
              & weight_basis_Gam(CoefGam(0:BasisOrderGam,0:BasisOrderGam,ibasis,ibin), tau1, tau2)
          enddo
        else
          do ibasis = 1, NBasis
            gam1 = gam1 + GamBasis(iorder, 1, 0, ibin, ibasis)* weight_basis( &
              & CoefGam(0:BasisOrder,0,ibasis,ibin), tau1)
          enddo
        endif

        write(104, *) real(gam1)*normal, dimag(gam1)*normal
      enddo
    enddo
    write(104, *)
  enddo

  do iorder = 1, MCOrder
    write(104, *) "##################################GammaR",trim(adjustl(str(iorder)))
    write(104, *) "#r:", VolFold
    write(104, *) "#Beta", Beta, "L", L(1), "Order", MCOrder
    do isite = 0, VolFold-1
      it1 = MxT/2
      it2 = MxT/2 
      ibin = get_bin_Gam(it1, it2)

      tau1 = dble(it1)*Beta/dble(MxT)
      tau2 = dble(it2)*Beta/dble(MxT)

      gam1 = (0.d0, 0.d0)
      if(IsBasis2D(ibin)) then
        do ibasis = 1, NBasisGam
          gam1 = gam1 + GamBasis(iorder, 1, isite, ibin, ibasis)*  &
            & weight_basis_Gam(CoefGam(0:BasisOrderGam,0:BasisOrderGam,ibasis,ibin), tau1, tau2)
        enddo
      else
        do ibasis = 1, NBasis
          gam1 = gam1 + GamBasis(iorder, 1, isite, ibin, ibasis)* weight_basis( &
            & CoefGam(0:BasisOrder,0,ibasis,ibin), tau1)
        enddo
      endif

      write(104, *) real(gam1)*normal, dimag(gam1)*normal
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
  call transfer_Chi_r(-1)

  close(104)
END SUBROUTINE output_Quantities



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
!!================================================================
!!================================================================

