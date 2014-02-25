  PROGRAM MAIN
    implicit none
    integer :: iblck, omega, ityp, iobs, Nblck
    integer :: Lx, Ly, MCOrder, Seed
    integer, parameter :: MxOmega = 16
    double precision :: Beta, temp
    integer, parameter :: Nobs_b = 67
    double precision :: GamMC(1:Nobs_b-1, 12000), GamNorm(12000), GamNormWeight
    double precision :: Obs(1:Nobs_b, 12000), Ave(1:Nobs_b), Dev(1:Nobs_b), Cor(1:Nobs_b)
    logical :: prt
    character*100 :: str, title

    print *, 'title' 
    read *, title 

    open(11, file="data_collapse/"//trim(title)//"_Gamma_blck_MC.dat")

    GamMC(:,:) = 0.d0
    GamNorm(:) = 0.d0
    Obs(:,:) = 0.d0
    Ave(:) = 0.d0
    Dev(:) = 0.d0
    Cor(:) = 0.d0

    iblck = 1
    do 
      read(11,*,End=100)
      read(11,*) str, Beta, str, str, Lx, Ly, str, MCOrder, str, Seed
      read(11,*) temp, temp, GamNormWeight, temp
      read(11,*)
      read(11,*)
      iobs = 0
      do omega = -MxOmega, MxOmega
        iobs = iobs + 1 
        read(11, *) temp, GamMC(iobs, iblck)
        Obs(iobs, iblck) = GamMC(iobs, iblck)
      enddo
      read(11,*)
      do omega = -MxOmega, MxOmega
        iobs = iobs + 1 
        read(11, *) temp, GamMC(iobs, iblck)
        Obs(iobs, iblck) = GamMC(iobs, iblck)
      enddo
      iobs = iobs + 1
      read(11, *) GamNorm(iblck)
      Obs(iobs, iblck) = GamNorm(iblck)
      read(11, *)
      iblck = iblck + 1
    enddo
100 close(11)
    Nblck = iblck - 1
    write(*, *) Nblck

    call stat_analy

    call norm_data
    call output_GamMC

  CONTAINS
    
  !==============Statistics ==========================================
  !! THIS IS PROJECT-INDEPENDENT 
  SUBROUTINE stat_analy
    implicit none
    integer          :: j, k, k0
    double precision :: devn, devp, nor

    ! -- calculate average -------------------------------------------
    nor  = 1.d0/(NBlck*1.d0)
    do j = 1, NObs_b
      Ave(j) = nor*Sum(Obs(j,1:NBlck))
    enddo

    Coarsen: do
      ! -- calculate error and t=1 correlation for basics obs.--------
      prt = .true.
      DO j = 1, NObs_b
        devp = 0.d0;  Cor(j) = 0.d0;  Dev(j) = 0.d0
        do k = 1,  NBlck
          devn   = Obs(j,k)-Ave(j)
          Dev(j) = Dev(j)+devn*devn
          Cor(j) = Cor(j)+devn*devp
          devp   = devn
        enddo 
        Dev(j)   = Dev(j)*nor;        Cor(j) = Cor(j)*nor
        if(Dev(j)>1.d-30)                Cor(j) = Cor(j)/Dev(j)
        Dev(j)   = dsqrt(Dev(j)/(NBlck-1.d0))
        !if(dabs(Cor(j))>tol) prt = .false.
      ENDDO 

      IF(prt)                         EXIT Coarsen 
      IF(NBlck<=64)    THEN
        prt = .false.;                EXIT Coarsen 
      ENDIF

      ! -- coarsen blocking ------------------------------------------
      NBlck = NBlck/2;      nor = nor*2.d0
      DO j = 1, NObs_b
        k0 = 1
        do k   = 1, NBlck
          Obs(j,k) = (Obs(j,k0)+Obs(j,k0+1))*0.5d0
          k0 = k0 +2
        enddo 
      ENDDO 
    enddo Coarsen 

    return
  END SUBROUTINE stat_analy
  !===================================================================

  SUBROUTINE norm_data
    implicit none
    integer          :: j, k, k0
    double precision :: nor

    nor = GamNormWeight/Ave(NObs_b)

    ! -- calculate deviation --------------------
    do j = 1, NObs_b-1
      Ave(j) = Ave(j)*nor
      Dev(j) = abs(Dev(j)*nor)+abs(Ave(j)*Dev(NObs_b)/Ave(NObs_b))
    enddo

    return
  END SUBROUTINE norm_data

  SUBROUTINE output_GamMC
    implicit none
    integer :: iconf, i, j, iobs, iomega
    double precision :: nor
    character*34 ich

    open(34, access="append", file=trim(title)//"_Gamma_MC.dat")

    write(34, *) "============================================"
    write(34, *) "Beta", Beta, "Lx, Ly", Lx, Ly, "Order", MCOrder, "Seed",Seed
    write(34, *) GamNormWeight, SUM(GamNorm(1:Nblck))
    write(34, *) "Nblck", Nblck
    if(prt) then
      iobs = 0
      write(34, *) "dx = 0, dy = 0"
      do iomega = -MxOmega, MxOmega
        iobs = iobs + 1
        write(34, *) iomega, Ave(iobs), Dev(iobs), Cor(iobs)
      enddo
      write(34, *) "dx = 1, dy = 0"
      do iomega = -MxOmega, MxOmega
        iobs = iobs + 1
        write(34, *) iomega, Ave(iobs), Dev(iobs), Cor(iobs)
      enddo
      write(34, *)
    endif

    close(34)
  END SUBROUTINE output_GamMC

END PROGRAM

