!!================================================================
!!================= READ/WRITE CONFIGURATIONS ====================
!!================================================================

!SUBROUTINE read_first_Gamma
  !implicit none
  !integer :: ix, iy, ityp, iomega1, iomega2

  !open(104, status='old', file=trim(title1)//"_first_order_Gamma.dat")
  !read(104, *)  GamNormWeight
  !do ix = 0, dLx
  !do iy = 0, dLy
    !do ityp = 1, ntypGa
      !do iomega1 = -MxOmegaDiag, MxOmegaDiag
        !do iomega2 = 0, MxOmegaDiag
          !read(104, *) Gam1MR(ix, iy, ityp, iomega1, iomega2) 
        !enddo
      !enddo
    !enddo
  !enddo
  !enddo
  !close(104)
!END SUBROUTINE read_first_Gamma



!SUBROUTINE write_first_Gamma
  !implicit none
  !integer :: ix, iy, ityp, iomega1, iomega2

  !open(104, status='replace', file=trim(title1)//"_first_order_Gamma.dat")
  !write(104, *)  GamNormWeight
  !do ix = 0, dLx
  !do iy = 0, dLy
    !do ityp = 1, ntypGa
      !do iomega1 = -MxOmegaDiag, MxOmegaDiag
        !do iomega2 = 0, MxOmegaDiag
          !write(104, *) Gam1MR(ix, iy, ityp, iomega1, iomega2) 
        !enddo
      !enddo
    !enddo
  !enddo
  !enddo
  !close(104)
!END SUBROUTINE write_first_Gamma


SUBROUTINE read_GWGamma
  implicit none
  integer :: ix, iy, ityp, it1, it2, ios

  open(100, status="old", file=trim(title1)//"_G_file.dat",iostat=ios)
  open(101, status="old", file=trim(title1)//"_W_file.dat",iostat=ios)
  open(102, status="old", file=trim(title1)//"_Gamma_file.dat",iostat=ios)

  if(.not. (ios==0)) then
    write(*,*) "You have to run self consistent loop first!"
    stop
  endif

  do it1 = 0, MxT-1
    do ityp = 1, NTypeG
      read(100, *) G(ityp, it1)
    enddo
  enddo

  do it1 = 0, MxT-1
    do iy = 0, Ly-1
      do ix = 0, Lx-1
        do ityp = 1, NTypeW
          read(101, *) W(ityp, ix, iy, it1)
        enddo
      enddo
    enddo
  enddo

  do it2 = 0, MxT-1
    do it1 = 0, MxT-1
      do iy = 0, Ly-1
        do ix = 0, Lx-1
          do ityp = 1, NTypeG
            read(102, *) Gam(ityp, ix, iy, it1, it2)
          enddo
        enddo
      enddo
    enddo
  enddo



  close(100)
  close(101)
  close(102)
  return
END SUBROUTINE read_GWGamma

SUBROUTINE write_GWGamma
  implicit none
  integer :: ix, iy, ityp
  integer :: it1, it2
  character*26 ich

  open(100, status="replace", file=trim(title1)//"_G_file.dat")
  open(101, status="replace", file=trim(title1)//"_W_file.dat")
  open(102, status="replace", file=trim(title1)//"_Gamma_file.dat")

  do it1 = 0, MxT-1
    do ityp = 1, NTypeG
      write(100, *) G(ityp, it1)
    enddo
  enddo

  do it1 = 0, MxT-1
    do iy = 0, Ly-1
      do ix = 0, Lx-1
        do ityp = 1, NTypeW
          write(101, *) W(ityp, ix, iy, it1)
        enddo
      enddo
    enddo
  enddo

  do it2 = 0, MxT-1
    do it1 = 0, MxT-1
      do iy = 0, Ly-1
        do ix = 0, Lx-1
          do ityp = 1, NTypeG
            write(102, *) Gam(ityp, ix, iy, it1, it2)
          enddo
        enddo
      enddo
    enddo
  enddo

  close(100)
  close(101)
  close(102)
  return
END SUBROUTINE write_GWGamma




SUBROUTINE write_monte_carlo_data
  implicit none
  integer :: iorder, itopo, ix, iy, ityp, it1, it2

  !=========== write into files =========================================
  open(104, status="replace", &
    & file=trim(title3)//"_monte_carlo_data.bin.dat",form="binary")

  write(104) imc, GamNorm, GamNormWeight
  do it2 = 0, MxT-1
    do it1 = 0, MxT-1
      do iy = 0, Ly-1
        do ix = 0, Lx-1
          do ityp = 1, NtypeGam/2
            do itopo = 0, 1
              do iorder = 0, MCOrder
                write(104)  GamMC(iorder, itopo, ityp, ix, iy, it1, it2)
                write(104)  GamSqMC(iorder, itopo, ityp, ix, iy, it1, it2)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
  !=========  write on the screen ========================================

  
  close(104)
END SUBROUTINE write_monte_carlo_data

SUBROUTINE write_monte_carlo_conf
  implicit none
  integer :: i

  open(103, status="replace", &
   & file=trim(title3)//"_monte_carlo_conf.bin.dat",form="binary")

  write(103)  Order, NGLn, NWLn, NVertex, MeasureGam, SignFermiLoop, IsWormPresent 
  write(103)  Ira, Masha, SpinMasha, kMasha
  write(103)  TailLn, TailVertex
  write(103)  NextLn(:)
  write(103)  NextVertex(:)
  write(103)  WeightCurrent, Phase
  write(103)  StatusLn(:)
  write(103)  GLnKey2Value(1:NGLn), WLnKey2Value(1:NWLn)
  write(103)  LnValue2Key(:)
  write(103)  StatusVertex(:)
  write(103)  VertexKey2Value(:)
  write(103)  VertexValue2Key(:)
  write(103)  KindLn(:)
  write(103)  TypeLn(:)
  write(103)  SpInVertex(:, :)
  write(103)  TypeVertex(:)
  write(103)  kLn(:)
  write(103)  Hash4G(:) 
  write(103)  Hash4W(:)
  write(103)  GXVertex(:), GYVertex(:)
  write(103)  WXVertex(:), WYVertex(:)
  write(103)  T1Vertex(:), T2Vertex(:), T3Vertex(:)
  write(103)  DirecVertex(:)             
  write(103)  WeightLn(:) 
  write(103)  WeightVertex(:) 

  do i = 1, 2
    write(103)  NeighLn(i,:)
  enddo

  do i = 1, 3
    write(103)  NeighVertex(i,:) 
  enddo

  close(103)
END SUBROUTINE write_monte_carlo_conf


SUBROUTINE read_monte_carlo_data
  implicit none
  integer :: iorder, ix, iy, ityp, it1, it2, itopo

  open(105, status="old", file=trim(title)//"_monte_carlo_data.bin.dat",form="binary")

  read(105) imc, GamNorm, GamNormWeight
  do it2 = 0, MxT-1
    do it1 = 0, MxT-1
      do iy = 0, Ly-1
        do ix = 0, Lx-1
          do ityp = 1, NtypeGam/2
            do itopo = 0, 1
              do iorder = 0, MCOrder
                read(105)  GamMC(iorder, itopo, ityp, ix, iy, it1, it2)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
  

  close(105)
END SUBROUTINE read_monte_carlo_data



SUBROUTINE read_monte_carlo_conf
  implicit none
  integer :: i

  open(106, status="old", file=trim(title)//"_monte_carlo_conf.bin.dat",form="binary")

  read(106)  Order, NGLn, NWLn, NVertex, MeasureGam, SignFermiLoop, IsWormPresent 
  read(106)  Ira, Masha, SpinMasha, kMasha
  read(106)  TailLn, TailVertex
  read(106)  NextLn(:)
  read(106)  NextVertex(:)
  read(106)  WeightCurrent, Phase
  read(106)  StatusLn(:)
  read(106)  GLnKey2Value(1:NGLn), WLnKey2Value(1:NWLn)
  read(106)  LnValue2Key(:)
  read(106)  StatusVertex(:)
  read(106)  VertexKey2Value(:)
  read(106)  VertexValue2Key(:)
  read(106)  KindLn(:)
  read(106)  TypeLn(:)
  read(106)  SpInVertex(:, :)
  read(106)  TypeVertex(:)
  read(106)  kLn(:)
  read(106)  Hash4G(:) 
  read(106)  Hash4W(:)
  read(106)  GXVertex(:), GYVertex(:)
  read(106)  WXVertex(:), WYVertex(:)
  read(106)  T1Vertex(:), T2Vertex(:), T3Vertex(:)
  read(106)  DirecVertex(:)             
  read(106)  WeightLn(:) 
  read(106)  WeightVertex(:) 

  do i = 1, 2
    read(106)  NeighLn(i,:)
  enddo

  do i = 1, 3
    read(106)  NeighVertex(i,:) 
  enddo

  close(106)
END SUBROUTINE read_monte_carlo_conf
!!================================================================
!!================================================================
!!================================================================





  
!!================================================================
!!========== PRINT OUT THE DATA FILES ============================
!!================================================================

SUBROUTINE output_Quantities
  implicit none
  integer :: ifile, iorder
  integer :: omega, omega1, i, ityp, dx, dy
  integer :: omega2
  integer :: it1, it2, it
  double precision :: GGI
  double precision :: GaR1, GaR2, GaR3
  double precision :: WWR1, WWR2, WWR3, WWR4

  open(15, file=trim(title1)//"_Chi.dat")
  !open(16, file=trim(title1)//"_Sigma.dat")

  !open(11, access="append", file=trim(title1)//"_G.dat")
  !open(12, access="append", file=trim(title1)//"_G_omega_0.dat")
  !open(13, access="append", file=trim(title1)//"_W.dat")
  !open(14, access="append", file=trim(title1)//"_W_omega_0.dat")
  !open(15, access="append", file=trim(title1)//"_Chi.dat")
  !open(16, access="append", file=trim(title1)//"_Chi_omega_0.dat")
  !open(17, access="append", file=trim(title1)//"_Gamma.dat")
  !open(18, access="append", file=trim(title1)//"_Gamma_omega_0.dat")
  !open(19, access="append", file=trim(title1)//"_Sigma.dat")
  !open(20, access="append", file=trim(title1)//"_Sigma_omega_0.dat")
  !open(21, access="append", file=trim(title1)//"_Pi.dat")
  !open(22, access="append", file=trim(title1)//"_Pi_omega_0.dat")
  !open(23, access="append", file=trim(title1)//"_Gamma_matrix.dat")

  do it = 0, MxT-1
    do dy = 0, dLy
      do dx = 0, dLx
        write(15, *) dx, dy, it, Chi(dx, dy, it)
      enddo
    enddo
  enddo

  close(15)

  !do it = 0, MxT-1
    !write(16, *) it, real(Sigma(it)), dimag(Sigma(it))
  !enddo

  !close(16)

  !do ifile = 11, 22
    !write(ifile, *) "=============================================================="
    !write(ifile, *) "Beta", Beta, "Lx, Ly", Lx, Ly, "Order", MCOrder, "Seed",Seed
    !write(ifile, *) ime, GamNorm, GamNormWeight
  !enddo

  !do omega = -MxOmegaG2, MxOmegaG2
    !GGI = weight_G(omega, 2)
    !write(11, *) omega, GGI
    !if(omega==0)      write(12, *) GGI
  !enddo

  !close(11)
  !close(12)

  !do ityp = 1, 2, 5
    !do omega = -MxOmegaW2, MxOmegaW2
      !WWR1 = weight_W( 0, 0, omega, 1)
      !WWR2 = weight_W( 0, 1, omega, 1)
      !WWR3 = weight_W( 1, 1, omega, 1)
      !write(13, *) ityp, omega, WWR1, WWR2, WWR3

      !if(omega==0)      write(14, *) ityp, WWR1, WWR2, WWR3
    !enddo
  !enddo

  !close(13)
  !close(14)

  !do dx = 0, dLx
    !do dy = 0, dLy
      !do omega = -MxOmegaChi, MxOmegaChi
        !write(15, *) Beta, dx, dy, omega, -3.d0*trChiR(dx, dy, omega)/Beta
        !if(omega==0)      write(16, *) Beta, dx, dy, -3.d0*trChiR(dx, dy, 0)/Beta
      !enddo
    !enddo
  !enddo

  !close(15)
  !close(16)

  !do omega1 = -MxOmegaDiag, MxOmegaDiag
    !do omega2 = -MxOmegaDiag, MxOmegaDiag
      !GaR1 = weight_Gamma(1, 0, omega1, omega2, 1)
      !write(23, '(f14.8)', advance='no') GaR1 
    !enddo

    !write(23, *)
    !GaR1 = weight_Gamma(0, 0, omega1, omega1, 1)
    !GaR2 = weight_Gamma(1, 0, omega1, omega1, 1)
    !GaR3 = weight_Gamma(1, 1, omega1, omega1, 1)
    !write(17, *) GaR1-1.d0, GaR2, GaR3
    !if(omega1==0)      write(18, *) GaR1-1.d0, GaR2, GaR3
  !enddo

  !close(17)
  !close(18)
  !close(23)

  !do iorder = 0, MCOrder
    !write(19, *) "Order:", iorder
    !do omega = -MxOmegaSigma, MxOmegaSigma
      !write(19, *) Beta, omega, SigmaI(iorder, 1, omega)
      !if(omega==0)      write(20, *) Beta, iorder, SigmaI(iorder, 1, 0)
    !enddo
  !enddo

  !close(19)
  !close(20)

  !do ityp = 1, ntypPi
    !do dx = 0, dLx
      !do dy = 0, dLy
        !do omega = -MxOmegaChi, MxOmegaChi
          !write(21, *) Beta, dx, dy, omega, PiR(ityp, dx, dy, omega)
          !if(omega==0)      write(22, *) Beta, dx, dy, PiR(ityp, dx, dy, 0)
        !enddo
      !enddo
    !enddo
  !enddo

  !close(21)
  !close(22)
END SUBROUTINE output_Quantities



!SUBROUTINE output_prob_MC
  !implicit none
  !integer :: omega, omega1, omega2, ityp, iorder, iloop
  !integer :: iconf, i, j, iobs, iomega
  !double precision :: Ga0R,GaMR1, GaCR1, GaMR2, GaMR3
  !double precision :: Norm

  !open(33, access="append", file=trim(title3)//"_prob_MC.dat")
  !!open(34, access="append", file=trim(title3)//"_over_complete_1_2.dat")
  !!open(35, access="append", file=trim(title3)//"_over_complete_3_4.dat")
  !!open(36, access="append", file=trim(title3)//"_over_complete_7_8.dat")
  !!open(37, access="append", file=trim(title3)//"_over_complete_9_10.dat")

  !write(33, *) "============================================"
  !write(33, *) "Beta", Beta, "Lx, Ly", Lx, Ly, "Order", MCOrder, "Seed",Seed
  !write(33, *) " 1: create worm along wline"
  !write(33, *) " 2: delete worm along wline"
  !write(33, *) " 3: create worm along gline"
  !write(33, *) " 4: delete worm along gline"
  !write(33, *) " 5: move worm along wline"
  !write(33, *) " 6: move worm along gline"
  !write(33, *) " 7: add interaction"
  !write(33, *) " 8: remove interaction"
  !write(33, *) " 9: add interaction cross"
  !write(33, *) "10: remove interaction cross"
  !write(33, *) "11: reconnect"
  !write(33, *) "12: shift gline in space"
  !write(33, *) "13: shift wline in space"
  !write(33, *) "14: change Gamma type"
  !write(33, *) "15: move measuring index"

  !write(33, *) "MC", imc, "Measure", ime
  !do iorder = 0, MCOrder
    !write(33, *) "Order", iorder
    !do i = 1, Nupdate
      !if(ProbProp(iorder, i)/=0.d0) then
        !write(33, *) i, ProbProp(iorder, i), ProbAcc(iorder, i), ProbAcc(iorder, i)/ProbProp(iorder, i)
      !endif
    !enddo
    !write(33, *)
  !enddo
  !write(33, *)

  !write(33, *)  "Physical configurations for different orders"
  !do iorder = 0, MCOrder
    !!if(iorder ==0) then
      !!write(33, *) iorder, GamOrder(iorder), GamWormOrder(iorder), Gam0Bubble
    !!else
      !!write(33, *) iorder, GamOrder(iorder), GamWormOrder(iorder)
    !!endif
    !write(33, *) iorder, GamOrder(iorder), GamWormOrder(iorder)
  !enddo
  !write(33, *) 

  !write(33, *)  "Average Weight Ratios"
  !do iorder = 0, MCOrder
    !write(33, *) iorder, 2, AveWeightRatio(iorder, 1)/ProbProp(iorder,2)
    !write(33, *) iorder, 4, AveWeightRatio(iorder, 2)/ProbProp(iorder,4)
  !enddo
  !write(33, *) 

  !!write(33, *)  "Histogram for omega on W"
  !!do omega = -MxOmega, MxOmega
    !!if(HistoOmegaW(omega)/=0.d0) then
      !!write(33, *) omega, HistoOmegaW(omega)
    !!endif
  !!enddo

  !!if(ProbAcc(1)/=0.d0 .and. ProbAcc(2)/=0.d0) then
    !!write(34, *) ime, imc, (ProbAcc(1)-ProbAcc(2))/sqrt(ProbAcc(1)+ProbAcc(2))
  !!endif

  !!if(ProbAcc(3)/=0.d0 .and. ProbAcc(4)/=0.d0) then
    !!write(35, *) ime, imc, (ProbAcc(3)-ProbAcc(4))/sqrt(ProbAcc(3)+ProbAcc(4))
  !!endif

  !!if(ProbAcc(7)/=0.d0 .and. ProbAcc(8)/=0.d0) then
    !!write(36, *) ime, imc, (ProbAcc(7)-ProbAcc(8))/sqrt(ProbAcc(7)+ProbAcc(8))
  !!endif

  !!if(ProbAcc(9)/=0.d0 .and. ProbAcc(10)/=0.d0) then
    !!write(37, *) ime, imc, (ProbAcc(9)-ProbAcc(10))/sqrt(ProbAcc(10)+ProbAcc(9))
  !!endif

  !close(33)
  !!close(34)
  !!close(35)
  !!close(36)
  !!close(37)
!END SUBROUTINE output_prob_MC


!SUBROUTINE output_GamMC
  !implicit none
  !integer :: iconf, i, j, iorder, iomega, iomega2
  !double precision :: err, percenterr, norm
  !double precision :: gam1, gam2, gamn

  !open(34, access="append", file=trim(title3)//"_Gamma_MC_matrix.dat")
  !open(35, access="append", file=trim(title3)//"_Gamma_MC.dat")
  !open(36, access="append", file=trim(title3)//"_Gamma_order2_MC.dat")

  !write(34, *) "============================================"
  !write(34, *) "Beta", Beta, "Lx, Ly", Lx, Ly, "Order", MCOrder, "Seed",Seed
  !write(34, *) imc, ime, GamNormWeight, GamNorm

  !write(35, *) "============================================"
  !write(35, *) "Beta", Beta, "Lx, Ly", Lx, Ly, "Order", MCOrder, "Seed",Seed
  !write(35, *) imc, ime, GamNormWeight, GamNorm

  !write(36, *) "============================================"
  !write(36, *) "Beta", Beta, "Lx, Ly", Lx, Ly, "Order", MCOrder, "Seed",Seed
  !write(36, *) imc, ime, GamNormWeight, GamNorm

  !norm = GamNormWeight*ime/GamNorm

  !write(34, *) "Order 1, dx=0, dy=0"
  !do iomega = -MxOmegaDiag, MxOmegaDiag
    !do iomega2 = -MxOmegaDiag, MxOmegaDiag
      !gam1 = GamMC(1,0,1, 0, 0, iomega, iomega2)/ime
      !gamn = gam1*norm
      !write(34, '(f14.8)', advance='no') gamn
    !enddo
    !write(34, *) 
  !enddo
  !write(34, *) 

  !write(34, *) "Order 2, dx=0, dy=0"
  !do iomega = -MxOmegaDiag, MxOmegaDiag
    !do iomega2 = -MxOmegaDiag, MxOmegaDiag
      !gam1 = GamMC(2, 0, 1, 0, 0, iomega, iomega2)/ime
      !gamn = gam1*norm
      !write(34, '(f14.8)', advance='no') gamn
    !enddo
    !write(34, *) 
  !enddo
  !write(34, *) 

  !write(34, *) "Order 2, dx=1, dy=0"
  !do iomega = -MxOmegaDiag, MxOmegaDiag
    !do iomega2 = -MxOmegaDiag, MxOmegaDiag
      !gam1 = GamMC(2, 1, 1, 1, 0, iomega, iomega2)/ime
      !gamn = gam1*norm
      !write(34, '(f14.8)', advance='no') gamn
    !enddo
    !write(34, *) 
  !enddo
  !write(34, *) 

  !write(34, *) "Order 3, dx=0, dy=0"
  !do iomega = -MxOmegaDiag, MxOmegaDiag
    !do iomega2 = -MxOmegaDiag, MxOmegaDiag
      !gam1 = GamMC(3, 0, 1, 0, 0, iomega, iomega2)/ime
      !gamn = gam1*norm
      !write(34, '(f14.8)', advance='no') gamn
    !enddo
    !write(34, *) 
  !enddo
  !write(34, *) 

  !write(34, *) "Order 3, dx=1, dy=0"
  !do iomega = -MxOmegaDiag, MxOmegaDiag
    !do iomega2 = -MxOmegaDiag, MxOmegaDiag
      !gam1 = GamMC(3, 1, 1, 1, 0, iomega, iomega2)/ime
      !gamn = gam1*norm
      !write(34, '(f14.8)', advance='no') gamn
    !enddo
    !write(34, *) 
  !enddo
  !write(34, *) 

  !write(34, *) "total, dx=1, dy=0"
  !do iomega = -MxOmegaDiag, MxOmegaDiag
    !do iomega2 = -MxOmegaDiag, MxOmegaDiag
      !gam1 = SUM(GamMC(0:MCOrder, 1, 1, 1, 0, iomega, iomega2))/ime
      !gamn = gam1*norm
      !write(34, '(f14.8)', advance='no') gamn
    !enddo
    !write(34, *) 
  !enddo
  !write(34, *) 

  !do iorder = 1, MCOrder
    !write(35, *) "Order", iorder
    !write(35, *) "dx = 0, dy = 0"
    !do iomega = -MxOmegaDiag, MxOmegaDiag
      !gam1 = GamMC(iorder, 0, 1, 0, 0, iomega, iomega)/ime
      !gam2 = GamSqMC(iorder,0, 1, 0, 0, iomega, iomega)/ime
      !err = sqrt(gam2-gam1**2.d0)/sqrt(ime-1)
      !if(abs(gam1)<1.d-30) then
        !percenterr = 0.d0
      !else
        !percenterr = err/abs(gam1)
      !endif
      !gamn = gam1*norm
      !write(35, *) iomega, gamn, percenterr*gamn, percenterr
    !enddo

    !write(35, *) "dx = 1, dy = 0, omega1=omega2"
    !do iomega = -MxOmegaDiag, MxOmegaDiag
      !gam1 = GamMC(iorder, 1, 1, 1, 0, iomega, iomega)/ime
      !gam2 = GamSqMC(iorder,1, 1, 1, 0, iomega, iomega)/ime
      !err = sqrt(gam2-gam1**2.d0)/sqrt(ime-1)
      !if(abs(gam1)<1.d-30) then
        !percenterr = 0.d0
      !else
        !percenterr = err/abs(gam1)
      !endif
      !gamn = gam1*norm
      !write(35, *) iomega, gamn, percenterr*gamn, percenterr
    !enddo

    !write(35, *) "dx = 1, dy = 0, omega2=-1"
    !do iomega = -MxOmegaDiag, MxOmegaDiag
      !gam1 = GamMC(iorder, 1, 1, 1, 0, iomega, -1)/ime
      !gam2 = GamSqMC(iorder,1, 1, 1, 0, iomega, -1)/ime
      !err = sqrt(gam2-gam1**2.d0)/sqrt(ime-1)
      !if(abs(gam1)<1.d-30) then
        !percenterr = 0.d0
      !else
        !percenterr = err/abs(gam1)
      !endif
      !gamn = gam1*norm
      !write(35, *) iomega, gamn, percenterr*gamn, percenterr
    !enddo
    !write(35, *) 
    !write(35, *) 
  !enddo

  !write(36, *) "topo 1:"
  !do iomega = -MxOmegaDiag, MxOmegaDiag
    !gam1 = Gam2Topo(1,iomega)/ime
    !gamn = gam1*norm
    !write(36, *) iomega, gamn
  !enddo

  !write(36, *) "topo 2:"
  !do iomega = -MxOmegaDiag, MxOmegaDiag
    !gam1 = Gam2Topo(2,iomega)/ime
    !gamn = gam1*norm
    !write(36, *) iomega, gamn
  !enddo


  !!write(36, *) "Order 2, typ=1, iloop = 2"
  !!do iomega = -MxOmegaDiag, MxOmegaDiag
    !!do iomega2 = -MxOmegaDiag, MxOmegaDiag
      !!gamn = Gam2MR(1, iomega, iomega2)
      !!write(36, '(f14.8)', advance='no') gamn
    !!enddo
    !!write(36, *) 
  !!enddo


  !close(34)
  !close(35)
  !close(36)
!END SUBROUTINE output_GamMC

SUBROUTINE write_monte_carlo_test
  implicit none
  !write(*, *) "conf(IsDelta(3)=1):", TestData(1)/TestData(0)
  !write(*, *) "conf(Status(3)=measure):", TestData(2)/TestData(0)
  return
END SUBROUTINE

!!================================================================
!!================================================================
!!================================================================

