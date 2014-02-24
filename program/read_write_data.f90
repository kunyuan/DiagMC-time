
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
  integer :: ix, iy, ityp, ibasis, ibasis2, iomega1, iomega2

  open(100, status="old", file=trim(title1)//"_G_file.dat")
  open(101, status="old", file=trim(title1)//"_W_file.dat")
  open(102, status="old", file=trim(title1)//"_Gamma_file.dat")


  do ityp = 1, NtypG
    do ibasis = 1, Nbasis
      read(100, *) GITailN(ityp, ibasis)
    enddo
    do iomega1 = -MxOmegaG1, MxOmegaG1
      read(100, *) GI(ityp, iomega1)
    enddo
    do ibasis = 1, Nbasis
      read(100, *) GITailP(ityp, ibasis)
    enddo
  enddo

  do ityp = 1, NtypW
    do ix = 0, Lx-1
      do iy = 0, Ly-1
        read(101, *) WRTailC(ityp, ix, iy)
        do iomega1 = 0, MxOmegaW1
          read(101, *) WR(ityp, ix, iy, iomega1)
        enddo
        do ibasis = 0, Nbasis
          read(101, *) WRTailP(ityp, ix, iy, ibasis)
        enddo
      enddo
    enddo
  enddo


  do ityp = 1, ntypGa
    do ix = 0, Lx-1
      do iy = 0, Ly-1
        read(102, *) GamRTailC(ityp, ix, iy)
        do iomega1 = -MxOmegaGamG1, MxOmegaGamG1
          do iomega2 = -MxOmegaGamG1, MxOmegaGamG1
            read(102, *) GamR(ityp, ix, iy, iomega1, iomega2)
          enddo
        enddo

        do ibasis = 1, nbasis
          do iomega1 = -MxOmegaGamG1, MxOmegaGamG1
            read(102, *) GamRTailMP(ityp, ix, iy, iomega1, ibasis)
            read(102, *) GamRTailMN(ityp, ix, iy, iomega1, ibasis)
          enddo
        enddo

        do ibasis = 1, nbasis
          do iomega2 = -MxOmegaGamG1, MxOmegaGamG1
            read(102, *) GamRTailPM(ityp, ix, iy, ibasis, iomega2)
            read(102, *) GamRTailNM(ityp, ix, iy, ibasis, iomega2)
          enddo
        enddo

        do ibasis = 1, nbasis
          read(102, *) GamRTailDiagP(ityp, ix, iy, ibasis)
          read(102, *) GamRTailDiagN(ityp, ix, iy, ibasis)
        enddo

        do ibasis = 1, nbasisGamma
          read(102, *) GamRTailPN(ityp, ix, iy, ibasis)
          read(102, *) GamRTailNP(ityp, ix, iy, ibasis)
          read(102, *) GamRTailPPR(ityp, ix, iy, ibasis)
          read(102, *) GamRTailNNR(ityp, ix, iy, ibasis)
          read(102, *) GamRTailPPL(ityp, ix, iy, ibasis)
          read(102, *) GamRTailNNL(ityp, ix, iy, ibasis)
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
  integer :: ix, iy, ityp, ibasis, ibasis2, iomega1, iomega2
  character*26 ich

  open(100, status="replace", file=trim(title1)//"_G_file.dat")
  open(101, status="replace", file=trim(title1)//"_W_file.dat")
  open(102, status="replace", file=trim(title1)//"_Gamma_file.dat")

  do ityp = 1, NtypG
    do ibasis = 1, Nbasis
      write(100, *) GITailN(ityp, ibasis)
    enddo
    do iomega1 = -MxOmegaG1, MxOmegaG1
      write(100, *) GI(ityp, iomega1)
    enddo
    do ibasis = 1, Nbasis
      write(100, *) GITailP(ityp, ibasis)
    enddo
  enddo

  do ityp = 1, NtypW
    do ix = 0, Lx-1
      do iy = 0, Ly-1
        write(101, *) WRTailC(ityp, ix, iy)
        do iomega1 = 0, MxOmegaW1
          write(101, *) WR(ityp, ix, iy, iomega1)
        enddo
        do ibasis = 0, Nbasis
          write(101, *) WRTailP(ityp, ix, iy, ibasis)
        enddo
      enddo
    enddo
  enddo

  do ityp = 1, ntypGa
    do ix = 0, Lx-1
      do iy = 0, Ly-1
        write(102, *) GamRTailC(ityp, ix, iy)
        do iomega1 = -MxOmegaGamG1, MxOmegaGamG1
          do iomega2 = -MxOmegaGamG1, MxOmegaGamG1
            write(102, *) GamR(ityp, ix, iy, iomega1, iomega2)
          enddo
        enddo

        do ibasis = 1, nbasis
          do iomega1 = -MxOmegaGamG1, MxOmegaGamG1
            write(102, *) GamRTailMP(ityp, ix, iy, iomega1, ibasis)
            write(102, *) GamRTailMN(ityp, ix, iy, iomega1, ibasis)
          enddo
        enddo

        do ibasis = 1, nbasis
          do iomega2 = -MxOmegaGamG1, MxOmegaGamG1
            write(102, *) GamRTailPM(ityp, ix, iy, ibasis, iomega2)
            write(102, *) GamRTailNM(ityp, ix, iy, ibasis, iomega2)
          enddo
        enddo

        do ibasis = 1, nbasis
          write(102, *) GamRTailDiagP(ityp, ix, iy, ibasis)
          write(102, *) GamRTailDiagN(ityp, ix, iy, ibasis)
        enddo

        do ibasis = 1, nbasisGamma
          write(102, *) GamRTailPN(ityp, ix, iy, ibasis)
          write(102, *) GamRTailNP(ityp, ix, iy, ibasis)
          write(102, *) GamRTailPPR(ityp, ix, iy, ibasis)
          write(102, *) GamRTailNNR(ityp, ix, iy, ibasis)
          write(102, *) GamRTailPPL(ityp, ix, iy, ibasis)
          write(102, *) GamRTailNNL(ityp, ix, iy, ibasis)
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
  integer :: iorder, itopo, ix, iy, ityp, iomega1, iomega2

  open(104, status="replace", file=trim(title3)//"_monte_carlo_data.dat")

  write(104, *) ime, GamNorm, GamNormWeight
  do iorder = 0, MCOrder
    do itopo = 0, 1
      do ityp = 1, ntypGa/2
        do ix = 0, Lx-1
          do iy = 0, Ly-1
            do iomega1 = -MxOmegaDiag, MxOmegaDiag
              do iomega2 = -MxOmegaDiag, MxOmegaDiag
                write(104, *)  GamMC(iorder, itopo, ityp, ix, iy, iomega1, iomega2)
                write(104, *)  GamSqMC(iorder, itopo, ityp, ix, iy, iomega1, iomega2)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
  
  do itopo = 1, 2
    do iomega1 = -MxOmegaDiag, MxOmegaDiag
      write(104, *) Gam2Topo(itopo, iomega1)
    enddo
  enddo

  close(104)
END SUBROUTINE write_monte_carlo_data

SUBROUTINE write_monte_carlo_conf
  implicit none
  integer :: ix, iy, ityp, ibasis, ibasis2, iomega1, iomega2
  integer :: i
  character*34 ich

  open(103, status="replace", file=trim(title3)//"_monte_carlo_conf.dat")

  write(103, *)  Order, NGLn, NWLn, NGam, MeasGamma, NFermiLoop, IsWormPresent 
  write(103, *)  Ira, Masha, SpinMasha, OmegaMasha, kMasha
  write(103, *)  TailLn, TailGam
  write(103, *)  NextLn(:)
  write(103, *)  NextVertex(:)
  write(103, *)  WeightCurrent, Phase
  write(103, *)  StatusLn(:)
  write(103, *)  Ln4GList(1:NGLn), Ln4WList(1:NWLn)
  write(103, *)  List4Ln(:)
  write(103, *)  StatusVertex(:)
  write(103, *)  Vertex4GamList(:)
  write(103, *)  List4Vertex(:)
  write(103, *)  KindLn(:)
  write(103, *)  TypeLn(:)
  write(103, *)  TypeVertexIn(:)
  write(103, *)  TypeVertexOut(:)
  write(103, *)  TypeVertex(:)
  write(103, *)  OmegaLn(:)
  write(103, *)  kLn(:)
  write(103, *)  Hash4G(:) 
  write(103, *)  Hash4W(:)
  write(103, *)  GXVertex(:), GYVertex(:)
  write(103, *)  WXVertex(:), WYVertex(:)
  write(103, *)  DirecVertex(:)             
  write(103, *)  WeightLn(:) 
  write(103, *)  WeightVertex(:) 

  do i = 1, 2
    write(103, *)  NeighLn(i,:)
  enddo

  do i = 1, 3
    write(103, *)  NeighVertex(i,:) 
  enddo

  close(103)
END SUBROUTINE write_monte_carlo_conf


SUBROUTINE read_monte_carlo_data
  implicit none
  integer :: iorder, ix, iy, ityp, iomega1, iomega2, itopo

  open(105, status="old", file=trim(title)//"_monte_carlo_data.dat")

  read(105, *) ime, GamNorm, GamNormWeight
  do iorder = 0, MCOrder
    do itopo = 0, 1
      do ityp = 1, ntypGa/2
        do ix = 0, Lx-1
          do iy = 0, Ly-1
            do iomega1 = -MxOmegaDiag, MxOmegaDiag
              do iomega2 = -MxOmegaDiag, MxOmegaDiag
                read(105, *)  GamMC(iorder, itopo, ityp, ix, iy, iomega1, iomega2)
                read(105, *)  GamSqMC(iorder, itopo, ityp, ix, iy, iomega1, iomega2)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  do itopo = 1, 2
    do iomega1 = -MxOmegaDiag, MxOmegaDiag
      read(105, *) Gam2Topo(itopo, iomega1)
    enddo
  enddo

  close(105)
END SUBROUTINE read_monte_carlo_data



SUBROUTINE read_monte_carlo_conf
  implicit none
  integer :: ix, iy, ityp, ibasis, ibasis2, iomega1, iomega2
  integer :: i

  open(106, status="old", file=trim(title)//"_monte_carlo_conf.dat")

  read(106, *)  Order, NGLn, NWLn, NGam, MeasGamma, NFermiLoop, IsWormPresent 
  read(106, *)  Ira, Masha, SpinMasha, OmegaMasha, kMasha
  read(106, *)  TailLn, TailGam
  read(106, *)  NextLn(:)
  read(106, *)  NextVertex(:)
  read(106, *)  WeightCurrent, Phase
  read(106, *)  StatusLn(:)
  read(106, *)  Ln4GList(1:NGLn), Ln4WList(1:NWLn)
  read(106, *)  List4Ln(:)
  read(106, *)  StatusVertex(:)
  read(106, *)  Vertex4GamList(:)
  read(106, *)  List4Vertex(:)
  read(106, *)  KindLn(:)
  read(106, *)  TypeLn(:)
  read(106, *)  TypeVertexIn(:)
  read(106, *)  TypeVertexOut(:)
  read(106, *)  TypeVertex(:)
  read(106, *)  OmegaLn(:)
  read(106, *)  kLn(:)
  read(106, *)  Hash4G(:) 
  read(106, *)  Hash4W(:)
  read(106, *)  GXVertex(:), GYVertex(:)
  read(106, *)  WXVertex(:), WYVertex(:)
  read(106, *)  DirecVertex(:)             
  read(106, *)  WeightLn(:) 
  read(106, *)  WeightVertex(:) 

  do i = 1, 2
    read(106, *)  NeighLn(i,:)
  enddo

  do i = 1, 3
    read(106, *)  NeighVertex(i,:) 
  enddo

  close(106)
END SUBROUTINE read_monte_carlo_conf

