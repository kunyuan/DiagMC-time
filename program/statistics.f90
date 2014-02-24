
  !--- PROJECT-DEPENDENT -------------------------------------------
  !==============Collect data ========================================
  !! THIS IS PROJECT-INDEPENDENT 
  SUBROUTINE coll_data(iblck)
    implicit none
    integer, intent(in) :: iblck
    integer             :: j, iobs
    double precision :: nor
    do iobs = 1, NObs
      Obs(iobs, iblck) = Quan(iobs)
    enddo
  END SUBROUTINE coll_data 

  !==============Statistics ==========================================
  !! THIS IS PROJECT-INDEPENDENT 
  SUBROUTINE stat_analy
    implicit none
    integer          :: j, k, k0
    double precision :: devn, devp, nor

    ! -- calculate average -------------------------------------------
    nor  = 1.d0/(NBlck*1.d0)
    do j = 1, NObs
      Ave(j) = nor*Sum(Obs(j,1:NBlck))
    enddo

    Coarsen: do
      ! -- calculate error and t=1 correlation for basics obs.--------
      prt = .true.
      DO j = 1, NObs
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
      DO j = 1, NObs
        k0 = 1
        do k   = 1, NBlck
          Obs(j,k) = (Obs(j,k0)+Obs(j,k0+1))*0.5d0
          k0 = k0 +2
        enddo 
      ENDDO 
    enddo Coarsen 

    ! -- define auxillary variables and average of composite obs.-----
    !call cal_Obs_comp

    ! -- calculate error and t=1 correlation for composite obs.-----
    !do j = 1+NObs, NObs
      !devp = 0.d0;  Cor(j) = 0.d0;  Dev(j) = 0.d0
      !DO k = 1,  NBlck
        !devn   = Obs(j,k)-Ave(j)
        !Dev(j) = Dev(j)+devn*devn
        !Cor(j) = Cor(j)+devn*devp
        !devp   = devn
      !ENDDO
      !Dev(j)   = Dev(j)*nor;        Cor(j) = Cor(j)*nor
      !IF(Dev(j)>eps)                Cor(j) = Cor(j)/Dev(j)
      !Dev(j)   = dsqrt(Dev(j)/(NBlck-1.d0))
    !enddo
    return
  END SUBROUTINE stat_analy
  !===================================================================


  SUBROUTINE norm_data
    implicit none
    integer          :: j, k, k0
    double precision :: nor

    nor = GamNormWeight/Ave(67)

    ! -- calculate deviation --------------------
    do j = 1, 66
      Ave(j) = Ave(j)*nor
      Dev(j) = abs(Dev(j)*nor)+abs(Ave(j)*Dev(67)/Ave(67))
    enddo

    return
  END SUBROUTINE norm_data
