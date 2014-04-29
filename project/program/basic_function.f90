!============== BASIC PROCEDURES IN UPDATE ==========================

!============== translate G,W,Gamma state ===========================

logical FUNCTION is_mea_worm_near_G(stat)
  implicit none
  integer :: stat
  if(stat==0)then
    is_mea_worm_near_G=.false.
  elseif(stat==1) then
    is_mea_worm_near_G=.true.
  else
    call LogFile%WriteStamp('e')
    call LogFile%WriteLine("The number of update: "+str(imc)+",imc: "+str(imc))
    call LogFile%WriteLine("line status error!"+str(stat))
    call print_config
    stop
  endif
end FUNCTION

logical FUNCTION is_mea_near_W(stat)
  implicit none
  integer :: stat
  if(stat==1 .or. stat==3)then
    is_mea_near_W=.true.
  else
    is_mea_near_W=.false.
  endif
end FUNCTION

logical FUNCTION is_worm_near_W(stat)
  implicit none
  integer :: stat
  if(stat==2)then
    is_worm_near_W=.true.
  else
    is_worm_near_W=.false.
  endif
end FUNCTION

logical FUNCTION is_mea_near_Gamma(stat)
  implicit none
  integer :: stat
  if(stat==1 .or. stat==3)then
    is_mea_near_Gamma=.true.
  else
    is_mea_near_Gamma=.false.
  endif
end FUNCTION

logical FUNCTION is_worm_near_Gamma(stat)
  implicit none
  integer :: stat
  if(stat==2)then
    is_worm_near_Gamma=.true.
  else
    is_worm_near_Gamma=.false.
  endif
end FUNCTION

!====================================================================
!============================== WEIGHTS =============================
!====================================================================

!-------- the weight of a gline -------------------------
! tau = tau2-tau1
COMPLEX*16 FUNCTION weight_gline(stat, tau, typ)
  implicit none
  integer :: stat, typ
  double precision :: tau
  integer :: t, flag

  t = Floor(tau*MxT/Beta)

  if(stat == 0) then
    weight_gline = weight_G(typ, t)
  else if(stat == 1) then
    !  Have measuring vertex around
    weight_gline = weight_meas_G(t)
  else
    call LogFile%WriteStamp('e')
    call LogFile%WriteLine("The number of update: "+str(iupdate))
    call LogFile%WriteLine("line status error!"+str(stat))
    call print_config
    stop
  endif

  !---------------------- test1: fake function ------------------------------
  !t = Floor(tau*MxT/Beta)

  !flag = 0
  !if(t<0) then
    !tau = tau +Beta
    !flag = 1
  !endif

  !if(stat == 0) then

    !!weight_gline = cdexp(-(0.d0, 1.d0)*pi*tau/(2.d0*Beta))
    !weight_gline = (1.d0, 0.d0)
    !!weight_gline = dcmplx(Beta-tau, 0.d0)
    !if(flag==1) then
      !weight_gline = -1.d0*weight_gline
    !endif

  !else if(stat==1) then
    !weight_gline = weight_meas_G(t)
  !else
    !call LogFile%WriteStamp('e')
    !call LogFile%WriteLine("The number of update: "+str(iupdate))
    !call LogFile%WriteLine("line status error!"+str(stat))
    !stop
  !endif

  !---------------------- test2: uniform function -----------------------
  !if(stat==0 .or. stat==1) then
    !weight_gline = weight_meas_G(t)
  !else
    !call LogFile%WriteStamp('e')
    !call LogFile%WriteLine("The number of update: "+str(iupdate))
    !call LogFile%WriteLine("line status error!"+str(stat))
    !stop
  !endif
  !------------------------ end -----------------------------------------

  return
END FUNCTION weight_gline



!-------- the weight of a wline -------------------------

! dx = x2-x1;  dy = y2- y1; tau = tau2-tau1
COMPLEX*16 FUNCTION weight_wline(stat, isdelta, dr0, tau, typ)
  implicit none
  integer :: stat, isdelta, dx, dy, typ
  integer :: dr0(2)
  integer :: dr(2)
  double precision :: tau
  integer :: t

  t = Floor(tau*MxT/Beta)
  call diff_r(dr0, dr)

  if(stat == 0) then
    if(isdelta==0) weight_wline = weight_W(typ, dr, t)
    if(isdelta==1) weight_wline = weight_W0(typ, dr)
  else if(stat == 2) then
    ! Have Ira or Masha around 
    if(isdelta==0) weight_wline = weight_W(1, dr, t)
    if(isdelta==1) weight_wline = weight_W0(1, dr)
  else if(stat == 1 .or. stat==3) then
    ! Have measuring vertex around
    if(isdelta==0) weight_wline = weight_meas_W(dr, t)
    if(isdelta==1) weight_wline = (0.d0, 0.d0)
  else
    call LogFile%WriteStamp('e')
    call LogFile%WriteLine("The number of update: "+str(iupdate))
    call LogFile%WriteLine("line status error!"+str(stat))
    call print_config
    stop
  endif

  !---------------------- test1: fake function ----------------------------
  !if(stat >= 0 .and. stat<=3) then
    !if(isdelta==0 .and. dr(1)==0 .and. dr(2)==0) then
      !weight_wline = weight_meas_W(dr, t)
    !else 
      !weight_wline = (0.d0, 0.d0)
    !endif

    !!if(isdelta==0) weight_wline = weight_meas_W(dr, t)
    !!if(isdelta==1) weight_wline = weight_meas_W(dr, 0)
  !else
    !call LogFile%WriteStamp('e')
    !call LogFile%WriteLine("The number of update: "+str(iupdate))
    !call LogFile%WriteLine("line status error!"+str(stat))
    !stop
  !endif

  !---------------------- test2: uniform function ----------------------------
  !if(stat >= 0 .and. stat<=3) then
    !if(isdelta==0 .and. dr(1)==0 .and. dr(2)==0) then
      !weight_wline = weight_meas_W(dr, t)
    !else 
      !weight_wline = (0.d0, 0.d0)
    !endif

  !else
    !call LogFile%WriteStamp('e')
    !call LogFile%WriteLine("The number of update: "+str(iupdate))
    !call LogFile%WriteLine("line status error!"+str(stat))
    !stop
  !endif
  !------------------------ end -----------------------------------------

  return
END FUNCTION weight_wline

!-------- the weight of a vertex ------------------------
 !dx = xg-xw;  dy = yg-yw; dtau1 = tau3-tau2; dtau2 = tau1-tau3

COMPLEX*16 FUNCTION weight_vertex(stat, isdelta, dr0, dtau1, dtau2, typ)
  implicit none
  integer :: stat, dr(2), t1, t2, typ, isdelta
  integer :: dr0(2), flag
  double precision :: weight
  double precision :: dtau1, dtau2

  t1 = Floor(dtau1*MxT/Beta)
  t2 = Floor(dtau2*MxT/Beta)

  call diff_r(dr0, dr)

  if(stat==0) then
    if(isbold) then
      !----------------- for bold Gamma ------------------------------
      if(isdelta==0) weight_vertex = weight_Gam(typ, dr, t1, t2)
      if(isdelta==1) weight_vertex = weight_Gam0(typ, dr)
    else
      !----------------- for bare Gamma ------------------------------
      if(isdelta==0) weight_vertex = (0.d0, 0.d0)
      if(isdelta==1) weight_vertex = weight_Gam0(typ, dr)
    endif

  else if(stat==2) then
    if(isbold) then
      !----------------- for bold Gamma ------------------------------
      if(isdelta==0) weight_vertex = weight_Gam(typ, dr, t1, t2)
      if(isdelta==1) weight_vertex = weight_Gam0(typ, dr)
    else 
      !----------------- for bare Gamma ------------------------------
      if(isdelta==0) weight_vertex = (0.d0, 0.d0)
      if(isdelta==1) weight_vertex = weight_Gam0(typ, dr)
    endif

  else if(stat==1 .or. stat==3) then
    if(isdelta==0) weight_vertex = (0.d0, 0.d0)
    if(isdelta==1) weight_vertex = weight_meas_Gam0(typ, dr)
  else
    call LogFile%WriteStamp('e')
    call LogFile%WriteLine("The number of update: "+str(iupdate))
    call LogFile%WriteLine("vertex status error!"+str(stat))
    stop
  endif

  !---------------------- test1: fake function -----------------------------------
  !flag = 0
  !if(t1<0) then
    !dtau1 = dtau1 + Beta
    !flag = flag + 1
  !endif
  !if(t2<0) then
    !dtau2 = dtau2 + Beta
    !flag = flag + 1
  !endif

  !if(stat==0 .or. stat==2) then
    !if(isdelta==1) weight_vertex = weight_meas_Gam0(typ, dr)
    !if(isdelta==0) then
      !weight_vertex = (0.d0, 0.d0)
      !if(dr(1)==0 .and. dr(2)==0) then
        !if(typ==1 .or. typ==2 .or. typ==5 .or. typ==6) then
          !if(mod(flag, 2)==0) the
            !weight_vertex = dcmplx(dtau1**2.d0+dtau2**2.d0+1.d0, 0.d0)
            !!weight_vertex = (1.d0, 0.d0)
          !else 
            !weight_vertex = dcmplx(-1.d0*(dtau1**2.d0+dtau2**2.d0+1.d0), 0.d0)
            !!weight_vertex = (-1.d0, 0.d0)
          !endif
        !endif
      !endif
    !endif

  !else if(stat==1 .or. stat==3) then
    !if(isdelta==1) weight_vertex = weight_meas_Gam0(typ, dr)
    !if(isdelta==0) weight_vertex = (0.d0, 0.d0)

  !else
    !call LogFile%WriteStamp('e')
    !call LogFile%WriteLine("The number of update: "+str(iupdate))
    !call LogFile%WriteLine("vertex status error!"+str(stat))
    !stop
  !endif

  !---------------------- test2: uniform function ---------------------
  !if(stat>=0 .and. stat<=3) then
    !if(isdelta==1) weight_vertex = weight_meas_Gam0(typ, dr)
    !if(isdelta==0) weight_vertex = weight_meas_Gam0(1,   dr)

  !else
    !call LogFile%WriteStamp('e')
    !call LogFile%WriteLine("The number of update: "+str(iupdate))
    !call LogFile%WriteLine("vertex status error!"+str(stat))
    !stop
  !endif
  !------------------------ end -----------------------------------------
  return
END FUNCTION weight_vertex


!!======================== WEIGHT EXTRACTING =========================
!! most basic interface to the matrix element

!--------- weight for bare propagator ----------------
Complex*16 FUNCTION weight_G0(typ, t)
  implicit none
  integer, intent(in)  :: typ, t  
  double precision     :: tau
  complex(kind=8)      :: muc  

  muc = dcmplx(0.d0, Mu(1)*pi/(2.d0*Beta))
  tau = real(t)*Beta/MxT
  if(tau>=0) then
    weight_G0 = cdexp(muc*tau)/(1.d0, 1.d0) 
  else
    weight_G0 = -cdexp(muc*(tau+Beta))/(1.d0, 1.d0) 
  endif
  return
END FUNCTION weight_G0


!!--------- calculate weight for bare interaction ----
Complex*16 FUNCTION weight_W0(typ, dr)
  implicit none
  integer, intent(in) :: dr(2), typ
  integer :: dx1, dy1
  double precision :: ratio

  ratio = Jcp

  dx1 = dr(1);       dy1 = dr(2)
  if(dx1>=0  .and. dx1<L(1) .and. dy1>=0 .and. dy1<L(2)) then
    if(dx1>dL(1))     dx1 = L(1)-dx1
    if(dy1>dL(2))     dy1 = L(2)-dy1

    weight_W0 = (0.d0, 0.d0)

    if((dx1==1.and.dy1==0).or.(dx1==0.and.dy1==1)) then
      if(typ ==1 .or. typ == 2) then
        weight_W0 = dcmplx(0.25d0*ratio, 0.d0)
      else if(typ == 3 .or. typ == 4) then
        weight_W0 = dcmplx(-0.25d0*ratio, 0.d0)
      else if(typ == 5 .or. typ == 6) then
        weight_W0 = dcmplx(0.5d0*ratio, 0.d0)
      endif
    endif
  else
    call LogFile%QuickLog("Weight_W"+str(dx1)+str(dy1)+"dx, dy bigger than system size!")
    stop
  endif
END FUNCTION weight_W0


!!--------- calculate weight for bare Gamma ---------
COMPLEX*16 FUNCTION weight_Gam0(typ, dr)
  implicit none
  integer, intent(in)  :: dr(2), typ
  double precision :: ratio

  weight_Gam0 = (0.d0, 0.d0)

  if(dr(1)>=0 .and. dr(1)<L(1) .and. dr(2)>=0 .and. dr(2)<L(2)) then
    if(dr(1)==0.and.dr(2)==0) then
      if(typ==1 .or. typ==2 .or. typ==5 .or. typ==6) then
        weight_Gam0 = (1.d0, 0.d0)
      endif
    endif
  else
    call logFile%QuickLog("Weight_Gam"//str(dr(1))//str(dr(2))//"dx, dy bigger than system size!")
    stop
  endif
END FUNCTION weight_Gam0

!!--------- extract weight for G ---------
COMPLEX*16 FUNCTION weight_G(typ1, t1)
  implicit none
  integer, intent(in)  :: t1, typ1

  if(t1>=0) then
    weight_G = G(typ1, t1)
  else
    weight_G = -G(typ1, t1+MxT)
  endif
END FUNCTION weight_G

!!--------- extract weight for W ---------
COMPLEX*16 FUNCTION weight_W(typ1, dr, t1)
  implicit none
  integer, intent(in)  :: dr(2), t1, typ1

  if(t1>=0) then
    weight_W = W(typ1, dr(1), dr(2), t1)
  else
    weight_W = W(typ1, dr(1), dr(2), t1+MxT)
  endif
END FUNCTION weight_W

!!--------- extract weight for Gamma ---------
COMPLEX*16 FUNCTION weight_Gam(typ1, dr, t1, t2)
  implicit none
  integer, intent(in)  :: dr(2), t1, t2, typ1
  double precision :: GaR

  if(t1>=0 .and. t2>=0) then
    weight_Gam = Gam(typ1, dr(1), dr(2), t1, t2)
  else if(t1<0 .and. t2>=0) then
    weight_Gam = -Gam(typ1, dr(1), dr(2), t1+MxT, t2)
  else if(t1>=0 .and. t2<0) then
    weight_Gam = -Gam(typ1, dr(1), dr(2), t1, t2+MxT)
  else
    weight_Gam = Gam(typ1, dr(1), dr(2), t1+MxT, t2+MxT)
  endif
END FUNCTION weight_Gam

!!====================================================================
!=====================================================================
!==================FAKE WEIGHT FOR MEASURING AND WORM=================
!=====================================================================


!--------- worm weight function  ---------
DOUBLE PRECISION FUNCTION weight_worm(drg, drw, dtau)
  implicit none
  integer, intent(in)  :: drg(2), drw(2)
  double precision, intent(in) :: dtau

  weight_worm = 1.d0
  !weight_worm = weight_worm*dexp(-dtau)
  !weight_worm = weight_worm*dexp(-0.5d0*(abs(dxg)+abs(dyg)))
  !weight_worm = weight_worm*dexp(-0.5d0*(abs(dxw)+abs(dyw)))

  return
END FUNCTION weight_worm


!--------- measure weight function for Gamma ---------
complex*16 FUNCTION weight_meas_G(t1)
  implicit none
  integer, intent(in)  :: t1
  weight_meas_G = (1.d0, 0.d0)
  return
END FUNCTION weight_meas_G

complex*16 FUNCTION weight_meas_W(dr, t1)
  implicit none
  integer, intent(in)  :: dr(2)
  integer, intent(in) :: t1

  weight_meas_W = (1.d0, 0.d0)
  return
END FUNCTION weight_meas_W

!complex*16 FUNCTION weight_meas_Gam(ityp, dr, t1, t2)
  !implicit none
  !integer, intent(in)  :: dr(2), ityp, t1, t2

  !weight_meas_Gam = (0.d0, 0.d0)
  !if(dr(1)==0 .and. dr(2)==0) then
    !weight_meas_Gam = (1.d0, 0.d0)
  !endif
  !return
!END FUNCTION weight_meas_Gam

complex*16 FUNCTION weight_meas_Gam0(ityp, dr)
  implicit none
  integer, intent(in)  :: dr(2), ityp

  !----------------------------------------------------------------------------
  !                       (out)
  !                 (out) / b
  !                   d /||
  !                ====  ||
  !                   c \||
  !                  (in) \ a
  !                       (in)
  !----- type = 1:  a = up;   b = up;    c = up;   d = up  -------------------
  !----- type = 2:  a = down; b = down;  c = down; d = down  -----------------
  !----- type = 3:  a = up;   b = up;    c = down; d = down  -----------------
  !----- type = 4:  a = down; b = down;  c = up;   d = up  -------------------
  !----- type = 5:  a = up;   b = down;  c = up;   d = down  -------------------
  !----- type = 6:  a = down; b = up;    c = down; d = up  -----------------
  !----------------------------------------------------------------------------

  weight_meas_Gam0 = (0.d0, 0.d0)
  if(dr(1)==0 .and. dr(2)==0) then
    if(ityp==1 .or. ityp==2 .or. ityp==5 .or. ityp==6) then
      weight_meas_Gam0 = (1.d0, 0.d0)
    endif
  endif
  return
END FUNCTION weight_meas_Gam0



!!------------- definition of the system symmetry ----------------
SUBROUTINE def_symmetry
  implicit none
  integer :: i, j, omega

  CoefOfSymmetry(:,:) = 2.d0   !typ 1,2; 3,4; 5,6

  !do i = 1, L(1)-1
    !CoefOfSymmetry(i, :) = 2.d0* CoefOfSymmetry(i, :)
  !enddo

  !do j = 1, L(2)-1
    !CoefOfSymmetry(:, j) = 2.d0* CoefOfSymmetry(:, j)
  !enddo
  return
END SUBROUTINE def_symmetry

SUBROUTINE update_WeightCurrent
  implicit none
  integer :: i, ikey, Gam1, Gam2 
  double precision :: tau, tau1, tau2
  Complex*16 :: weight
  Complex*16 :: wln, wgam

  weight = (1.d0, 0.d0)
  do ikey = 1, NGLn
    i = GLnKey2Value(ikey)
    Gam1 = NeighLn(1,i);       Gam2 = NeighLn(2,i)
    tau = TVertex(2, Gam2)-TVertex(1,Gam1)
    wln = weight_gline(StatusLn(i),tau,TypeLn(i))
    WeightLn(i) = wln
    weight = weight *wln
  enddo

  do ikey = 1, NWLn
    i = WLnKey2Value(ikey)
    Gam1 = NeighLn(1,i);       Gam2 = NeighLn(2,i)
    wln = weight_wline(StatusLn(i),IsDeltaLn(i), WRVertex(:, Gam1)-WRVertex(:, Gam2), &
      &  TVertex(3, Gam2)-TVertex(3,Gam1), TypeLn(i))
    WeightLn(i) = wln
    weight = weight *wln
  enddo


  do ikey = 1, NVertex
    i = VertexKey2Value(ikey)
    tau1 = TVertex(3, i)-TVertex(2, i)
    tau2 = TVertex(1, i)-TVertex(3, i)
    wgam = weight_vertex(StatusVertex(i), IsDeltaVertex(i), GRVertex(:, i)-WRVertex(:, i), &
      & tau1, tau2, TypeVertex(i))
    WeightVertex(i) = wgam
    weight = weight *wgam
  enddo

  weight = weight*(-1.d0)**Order *SignFermiLoop

  WeightCurrent = abs(weight)
  Phase = weight/WeightCurrent

  return
END SUBROUTINE update_WeightCurrent




!!====================================================================
!!============== ELEMENTS GENERATORS =================================
!!====================================================================

!---------- double precision tau -------------------------
DOUBLE PRECISION FUNCTION generate_tau()
  implicit none
  double precision :: tau
  generate_tau = rn()*Beta
  return
END FUNCTION generate_tau

DOUBLE PRECISION FUNCTION prob_tau(tau)
  implicit none 
  double precision, intent(in) :: tau
  prob_tau = 1.d0/Beta
  return
END FUNCTION prob_tau

!---------- int x y -------------------------
SUBROUTINE generate_xy(CurrentR,NewR,dR,Weight,Flag)
!Please make sure Weight is already initialized before calling!
!Flag=.true.: generate new X,Y and new dX,dY
!Flag=.false.: generate new X,Y according to input dX,dY
  implicit none
  logical :: Flag
  integer :: NewR(2),CurrentR(2),dR(2),i
  double precision :: Weight,rand
  !dR: CurrentR --> NewR;   dR': CurrentR <-- NewR
  !CurrentR==NewR, dR=0, dR'=0
  !CurrentR>NewR, dR=NewR-CurrentR, dR'=CurrentR-NewR+L
  !CurrentR<NewR, dR=NewR-CurrentR+L, dR'=CurrentR-NewR

  do i=1,2
    if(Flag) then
      rand=rn()
      dR(i)=0.5d0*dexp(rand*logL(i))
      IF(rn()>0.5d0) dR(i)=L(i)-1-dR(i)
    endif
    if(dR(i)/=0) then
      Weight=Weight*SpatialWeight(i,dR(i))
      Weight=Weight/SpatialWeight(i,L(i)-dR(i))
    endif

    NewR(i) = CurrentR(i) + dR(i)
    if(NewR(i)>=L(i)) then
      NewR(i) = NewR(i) - L(i)
    endif
  enddo
  return
END SUBROUTINE generate_xy

SUBROUTINE diff_r(dr0, dr)
  implicit none
  integer, intent(in) :: dr0(2)
  integer, intent(out) :: dr(2)
  integer :: i
  do i = 1, 2
    dr(i) = dr0(i)
    if(dr(i)<0)     dr(i) = L(i)+ dr(i)
  enddo
END SUBROUTINE diff_r

!!---------- int k -------------------------
INTEGER FUNCTION generate_k()
  implicit none
  generate_k = Floor(rn()*(2*MxK))-MxK
  return
END FUNCTION generate_k

INTEGER FUNCTION add_k(k1, k2)
  implicit none
  integer :: k1, k2
  add_k = k1 + k2
  if(add_k>=MxK) then
    add_k = add_k -2*MxK
  else if(add_k<-MxK) then
    add_k = add_k +2*MxK
  endif
  return
END FUNCTION add_k

LOGICAL FUNCTION Is_k_valid(k)
  implicit none
  integer,intent(in) :: k
  if(k<MxK .and. k>=-MxK)  then
    Is_k_valid = .true.
  else
    Is_k_valid = .false.
  endif
END FUNCTION Is_k_valid


!----------- randomly pick a gline -------------
INTEGER FUNCTION generate_gline()
  implicit none
  integer :: rand
  rand = Floor(rn()*NGLn)+1
  generate_gline = GLnKey2Value(rand)
  return
END FUNCTION generate_gline

!----------- randomly pick a wline -------------
INTEGER FUNCTION generate_wline()
  implicit none
  integer :: rand
  rand = Floor(rn()*NWLn)+1
  generate_wline = WLnKey2Value(rand)
  return
END FUNCTION generate_wline

!----------- randomly pick a gamma -------------
INTEGER FUNCTION generate_vertex()
  implicit none
  integer :: rand
  rand = Floor(rn()*NVertex)+1
  generate_vertex = VertexKey2Value(rand)
  return
END FUNCTION generate_vertex
!!!=======================================================================
!!!=======================================================================
!!!=======================================================================



!!=======================================================================
!!================= Hash Table operations     ===========================
!!=======================================================================

!-------- update G Hash table -----------------------------------
SUBROUTINE update_Hash4G(newk, GLn)
  !the newk and oldk is alway about the same G,
  !PLEASE MAKE SURE THAT!!!
  !ALSO MAKE SURE Hash4G(newk)==0 BEFORE CALL THIS SUBROUTINE
  implicit none
  integer, intent(in) :: newk, GLn
  integer :: oldk
  oldk=kLn(GLn)

  if(CheckG .and. (Hash4G(newk)/=0 .or. Hash4G(oldk)/=GLn)) then
      call LogFile%WriteStamp('e')
      call LogFile%WriteLine("Oops, update_Hash4G found a bug!")
      call LogFile%WriteLine("IsWormPresent:"+str(IsWormPresent)+", update number:"+str(imc))
      call LogFile%WriteLine("G Hash table for old k"+str(newk)+" is not 1!!"+str(Hash4G(newk)))
      call print_config
  else
    Hash4G(newk)=Hash4G(oldk)
    Hash4G(oldk)=0
  endif
  return
END SUBROUTINE update_Hash4G

SUBROUTINE add_Hash4G(newk, GLn)
  implicit none
  integer, intent(in) :: newk, GLn
  !MAKE SURE Hash4G(newk)==0 BEFORE CALL THIS SUBROUTINE
  if(CheckG .and. Hash4G(newk)/=0) then
      call LogFile%WriteStamp('e')
      call LogFile%WriteLine("Oops, add_Hash4G found a bug!")
      call LogFile%WriteLine("IsWormPresent:"+str(IsWormPresent)+", update number:"+str(imc))
      call LogFile%WriteLine("G Hash table for old k"+str(newk)+" is not 1!!"+str(Hash4G(newk)))
      call print_config
  else
    Hash4G(newk)=GLn
  endif

END SUBROUTINE add_Hash4G

SUBROUTINE delete_Hash4G(oldk, GLn)
  implicit none
  integer, intent(in) :: oldk, GLn
  if(Hash4G(oldk)==GLn) then
    Hash4G(oldk)=0
  else
    if(CheckG) then
      call LogFile%WriteStamp('e')
      call LogFile%WriteLine("Oops, delete_Hash4G found a bug!")
      call LogFile%WriteLine("IsWormPresent:"+str(IsWormPresent)+", update number:"+str(imc))
      call LogFile%WriteLine("G Hash table for old k"+str(oldk)+" is not 1!!"+str(Hash4G(oldk)))
      call print_config
      stop
    endif
  endif
  return
END SUBROUTINE delete_Hash4G


!-------- update W Hash table -----------------------------------
SUBROUTINE update_Hash4W(newk, WLn)
  implicit none
  integer, intent(in) :: newk, WLn
  integer :: aoldk, anewk
  !COME ON, MAKE SURE Hash4W(aoldk)==0 BEFORE YOU CALL IT
  aoldk = abs(kLn(WLn))
  anewk = abs(newk)
  Hash4W(anewk)=Hash4W(aoldk)
  Hash4W(aoldk)=0
  return
END SUBROUTINE update_Hash4W

SUBROUTINE add_Hash4W(newk, WLn)
  implicit none
  integer, intent(in) :: newk, WLn
  integer :: anewk
  anewk = abs(newk)
  Hash4W(anewk)=WLn
END SUBROUTINE add_Hash4W

SUBROUTINE delete_Hash4W(oldk, WLn)
  implicit none
  integer, intent(in) :: oldk, WLn
  integer :: aoldk

  aoldk = abs(oldk)

  Hash4W(aoldk)=0
  return
END SUBROUTINE delete_Hash4W


!!=======================================================================
!!================= IRREDUCIBILITY CHECK ================================
!!=======================================================================

!--------- check the irreducibility for G -----------------------
LOGICAL FUNCTION Is_k_valid_for_G(k)
  implicit none
  integer,intent(in) :: k
  if(CheckG .and. Hash4G(k)/=0) then
    Is_k_valid_for_G=.false.
  else
    Is_k_valid_for_G=.true.
  endif
END FUNCTION

LOGICAL FUNCTION Is_reducible_G_Gam(GLn)
  implicit none
  integer, intent(in) :: GLn
  integer :: NeighGLn1, NeighGLn2
  integer :: nG, Gam1, Gam2, W1, W2, nW
  integer :: newk, kG 
  integer :: i, nnk1, nnk2

  !Is_reducible_G_Gam=Is_reducible_Gam()
  !return
  Is_reducible_G_Gam = .false.
  if(CheckGam==.false.) return

  Gam1 = NeighLn(1, GLn)
  W1 = NeighVertex(3, Gam1)
  NeighGLn1 = NeighVertex(1, Gam1)
  Gam2 = NeighLn(2, GLn)
  W2 = NeighVertex(3, Gam2)
  NeighGLn2 = NeighVertex(2, Gam2)

  newk = kLn(GLn)

  do i = 1, NGLn
    nG = GLnKey2Value(i)
    if(nG==GLn) cycle  !don't test GLn itself
    kG = kLn(nG)
    nW=Hash4W(abs(add_k(newk,-kG)))
    if(nW/=0) then
      if(nW==W1 .and. nG==NeighGLn1) cycle
      if(nW==W2 .and. nG==NeighGLn2) cycle
      Is_reducible_G_Gam = .true.
      return
    endif
  enddo
  return
END FUNCTION Is_reducible_G_Gam

!--------- check the irreducibility for W -----------------------

LOGICAL FUNCTION Is_k_valid_for_W(k)
  implicit none
  integer,intent(in) :: k
  integer :: ak
  ak=abs(k)
  if(CheckW .and. (ak==0 .or. Hash4W(ak)/=0)) then
    Is_k_valid_for_W=.false.
  else
    Is_k_valid_for_W=.true.
  endif
END FUNCTION

!--------- check the irreducibility for W -----------------------
LOGICAL FUNCTION Is_reducible_W_Gam(WLn)
  implicit none
  integer, intent(in) :: WLn
  integer :: absk, i, Gam1, Gam2, G1, G2, G3, G4
  integer :: nG, kG, pkGG, nkGG, pk, nk
  logical :: test

  Is_reducible_W_Gam = .false.
  if(CheckGam==.false.) return

  absk = abs(kLn(WLn))
  Gam1 = NeighLn(1, WLn)
  Gam2 = NeighLn(2, WLn)
  G1 = NeighVertex(1, Gam1)
  G2 = NeighVertex(2, Gam1)
  G3 = NeighVertex(1, Gam2)
  G4 = NeighVertex(2, Gam2)

  do i = 1, NGLn
    nG = GLnKey2Value(i)
    kG = kLn(nG)
    pkGG=Hash4G(add_k(kG, absk))
    if(pkGG/=0) then
      if((nG/=G1 .or. pkGG/=G2) .and. (nG/=G3 .or. pkGG/=G4) .and. (nG/=G2 .or. pkGG/=G1) .and. (nG/=G4 .or. pkGG/=G3)) then
        Is_reducible_W_Gam=.true.
        return
      endif
    endif

    nkGG=Hash4G(add_k(kG, -absk))
    if(nkGG/=0) then
      if((nG/=G1 .or. nkGG/=G2) .and. (nG/=G3 .or. nkGG/=G4) .and. (nG/=G2 .or. nkGG/=G1) .and. (nG/=G4 .or. nkGG/=G3)) then
        Is_reducible_W_Gam=.true.
        return
      endif
    endif
  enddo
END FUNCTION Is_reducible_W_Gam

!------------- check the irreducibility for add interaction operation --------------------

LOGICAL FUNCTION Is_reducible_add_interaction(kW, kIA, kMB, kIC, kMD)
  implicit none
  integer,intent(in) :: kW, kIA, kMB, kIC, kMD
  integer :: kG
  integer :: ktemp

  Is_reducible_add_interaction = .false.
  if(CheckGam==.false.) return

  do i = 1, NGLn
    kG = kLn(GLnKey2Value(i))
    !check the W which is the new interaction line added
    ktemp = add_k(kG, kW)
    !ktemp is G here
    if((kG/=kIA .or. ktemp/=kIC) .and. (kG/=kIC .or. ktemp/=kIA) &
      & .and. (kG/=kMB .or. ktemp/=kMD) .and. (kG/=kMD .or. ktemp/=kMB)) then
      if(Hash4G(ktemp)/=0) then
          Is_reducible_add_interaction=.true.
          return
      endif
    endif

    ktemp = add_k(kG, -kW)
    !ktemp is G here
    if((kG/=kIA .or. ktemp/=kIC) .and. (kG/=kIC .or. ktemp/=kIA) &
      & .and. (kG/=kMB .or. ktemp/=kMD) .and. (kG/=kMD .or. ktemp/=kMB)) then
      if(Hash4G(ktemp)/=0) then
          Is_reducible_add_interaction=.true.
          return
      endif
    endif

    if(kG/=kIA) then
    !check GIA
      !ktemp is W here
      ktemp = abs(add_k(kG, -kIA))
      if(Hash4W(ktemp)/=0) then
        if(ktemp/=abs(kW)) then
          Is_reducible_add_interaction=.true.
          return
        endif
      endif
    endif

    if(kG/=kMB) then
    !check GMB
      !ktemp is W here
      ktemp = abs(add_k(kG, -kMB))
      if(Hash4W(ktemp)/=0) then
        if(ktemp/=abs(kW)) then
          Is_reducible_add_interaction=.true.
          return
        endif
      endif
    endif
  enddo
END FUNCTION Is_reducible_add_interaction

SUBROUTINE test_reduciblility(is_reducible, mcname)
  implicit none
  logical :: is_reducible
  character(len=*) :: mcname
  if(is_reducible/=Is_reducible_Gam()) then
    call LogFile%QuickLog(mcname+", Reducibility check is wrong!"+" imc: "+str(imc), "e")
    call print_config
    stop
  endif
end SUBROUTINE

logical FUNCTION Is_reducible_Gam()
  implicit none
  integer :: i,j, k, Gi, Gj, Wk, Gam1, Gam2
  Is_reducible_Gam = .false.
  if(CheckGam==.false.) return

  do i = 1, NGLn
    Gi = GLnKey2Value(i)
    do j = i+1, NGLn
      Gj = GLnKey2Value(j)
      !if(NeighLn(1,Gi)==NeighLn(2,Gj) .or. NeighLn(2,Gi)==NeighLn(1,Gj)) cycle
      do k = 1, NWLn
        Wk = WLnKey2Value(k)
        if(abs(add_k(kLn(Gi), -kLn(Gj)))==abs(kLn(Wk))) then
          Gam1=NeighLn(1,Wk)
          Gam2=NeighLn(2,Wk)
          if(NeighVertex(1,Gam1)==Gi .and. NeighVertex(2,Gam1)==Gj) cycle
          if(NeighVertex(2,Gam1)==Gi .and. NeighVertex(1,Gam1)==Gj) cycle
          if(NeighVertex(1,Gam2)==Gi .and. NeighVertex(2,Gam2)==Gj) cycle
          if(NeighVertex(2,Gam2)==Gi .and. NeighVertex(1,Gam2)==Gj) cycle
          Is_reducible_Gam=.true.
          !if(imc==935670) print *,imc, Gi,kLn(Gi), Gj,kLn(Gj), Wk,kLn(Wk)
          return
        endif
      enddo
    enddo
  enddo
end FUNCTION
!!=======================================================================
!!=======================================================================
!!=======================================================================


!!=======================================================================
!!================= SUBROUTINES IN UPDATES ==============================
!!=======================================================================

!-- change the status from normal or measuring to i/m or measuring+i/m --
INTEGER FUNCTION add_ira_stat(stat)
  implicit none
  integer, intent(in) :: stat
  if(stat == -1) then
    call LogFile%QuickLog("add_ira_stat, stat = -1",'e')
    call print_config
    stop
  endif
  if(stat<=1) then
    add_ira_stat = stat + 2
  else
    add_ira_stat = stat
  endif
  if(add_ira_stat>=4) then
    call LogFile%QuickLog(str(imc)+" add_ira_stat(stat>3): "+str(add_ira_stat),'e')
    call print_config
    stop
  endif
END FUNCTION add_ira_stat

!-- change the status from i/m or measuring+i/m to normal or measuring--
INTEGER FUNCTION delete_ira_stat(stat)
  implicit none
  integer, intent(in) :: stat
  if(stat == -1) then
    call LogFile%QuickLog("del_ira_stat, stat = -1",'e')
    call print_config
    stop
  endif
  delete_ira_stat = stat - 2
  if(delete_ira_stat == -2) delete_ira_stat = 0
  if(delete_ira_stat == -1) delete_ira_stat = 1
  if(delete_ira_stat<0) then
    call LogFile%QuickLog(str(imc)+" del_ira_stat(stat<0): "+str(delete_ira_stat),'e')
    call print_config
    stop
  endif
END FUNCTION delete_ira_stat

!-- change the status from normal or i/m to measuring or measuring+i/m --
INTEGER FUNCTION add_mea_stat(stat)
  implicit none
  integer, intent(in) :: stat
  if(stat == -1) then
    call LogFile%QuickLog("add_mea, stat = -1",'e')
    call print_config
    stop
  endif

  add_mea_stat = stat + 1
  if(add_mea_stat == 2)   add_mea_stat = 1

  if(add_mea_stat>=4) then
    call LogFile%QuickLog(str(imc)+" add_mea(stat>3): "+str(add_mea_stat),'e')
    call print_config
    stop
  endif
END FUNCTION add_mea_stat

!-- change the status from measuring or measuring+i/m to normal or i/m--
INTEGER FUNCTION delete_mea_stat(stat)
  implicit none
  integer,intent(in) :: stat
  if(stat == -1) then
    call LogFile%QuickLog("del_mea_stat, stat = -1",'e')
    call print_config
    stop
  endif
  delete_mea_stat = stat - 1
  if(delete_mea_stat<0) then
    call LogFile%QuickLog(str(imc)+"del_mea_stat(stat<0): "+str(delete_mea_stat),'e')
    call print_config
    stop
  endif
END FUNCTION delete_mea_stat

!--- calculate the status of a wline according to the neighbor vertexes --
INTEGER FUNCTION wline_stat(stat1, stat2)
  implicit none
  integer, intent(in) :: stat1, stat2
  if(stat1 == -1) stop
  if(stat2 == -1) stop
  wline_stat = stat1 + stat2
  if(stat1 == stat2) then
    wline_stat = stat1
  else if(wline_stat == 5) then
    wline_stat = 3
  endif

  if(wline_stat>3) then
    call LogFile%WriteStamp('e')
    call LogFile%WriteLine("wline_stat error!")
    call LogFile%WriteLine(str(wline_stat)+'  1: '+str(stat1)+'  2: '+str(stat2))
    call print_config
    stop
  endif
  return
END FUNCTION wline_stat

!--- calculate the status of a line according to the neighbor vertexes --
INTEGER FUNCTION gline_stat(stat1, stat2)
  implicit none
  integer, intent(in) :: stat1, stat2
  integer :: mstat1, mstat2
  if(stat1 == -1) stop
  if(stat2 == -1) stop
  mstat1 = Mod(stat1, 2)
  mstat2 = Mod(stat2, 2)
  gline_stat = mstat1 + mstat2
  if(mstat1 == mstat2) then
    gline_stat = mstat1
  endif

  if(gline_stat>1) then
    call LogFile%WriteStamp('e')
    call LogFile%WriteLine("gline_stat error!")
    call LogFile%WriteLine(str(gline_stat)+'  1: '+str(stat1)+'  2: '+str(stat2))
    call print_config
    stop
  endif
  return
END FUNCTION gline_stat

!------------- update a line  -------------------------------
!mainly to add the new k to the hash table
SUBROUTINE update_line(CurLine, k, knd)
  implicit none
  integer, intent(in) :: CurLine, k, knd

  if(knd==1) then
    call update_Hash4G(k, CurLine)
  else
    call update_Hash4W(k, CurLine)
  endif
  kLn(CurLine)=k
end SUBROUTINE 

SUBROUTINE undo_update_line(CurLine, k, knd)
  implicit none
  integer, intent(in) :: CurLine, k, knd
  if(knd==1) then
    call update_Hash4G(k, CurLine)
  else
    call update_Hash4W(k, CurLine)
  endif
  kLn(CurLine)=k
end SUBROUTINE


!------------- insert a line to the link -------------------
SUBROUTINE insert_line(newline, isdelta, k, knd, typ, stat, weigh)
  implicit none
  integer, intent(out) :: newline
  integer, intent(in) :: isdelta, k, knd, typ, stat
  complex*16, intent(in) :: weigh

  newline = TailLn
  TailLn  = NextLn(TailLn)
  if(StatusLn(TailLn)>=0) then
    call LogFile%QuickLog("update: "+str(imc)+", insert_line error!!!", 'e')
    call print_config
    stop
  endif

  if(TailLn == -1) then
    call LogFile%QuickLog("insert_line error! tail=-1!", 'e')
    call print_config
    stop
  endif

  if(knd==1) then
    NGLn = NGLn + 1
    GLnKey2Value(NGLn) = newline
    LnValue2Key(newline) = NGLn
    call add_Hash4G(k,newline)
  else
    NWLn = NWLn + 1
    WLnKey2Value(NWLn) = newline
    LnValue2Key(newline) = NWLn
    call add_Hash4W(k,newline)
  endif

  kLn(newline) = k
  IsDeltaLn(newline) = isdelta
  KindLn(newline) = knd
  TypeLn(newline) = typ
  StatusLn(newline) = stat
  WeightLn(newline) = weigh
  return
END SUBROUTINE insert_line

!------------- undo_insert a line from the link -------------------
SUBROUTINE undo_insert_line(occline, knd)
  implicit none
  integer, intent(in) :: occline, knd
  integer :: tmp


  if(StatusLn(occline)==-1) then
    call LogFile%QuickLog("update: "+str(imc)+" , undo_insert_line error!!!", 'e')
    call print_config
    stop
  endif

  if(knd==1) then

    tmp = GLnKey2Value(NGLn)
    GLnKey2Value(NGLn) = 0
    GLnKey2Value(LnValue2Key(occline)) = tmp
    LnValue2Key(tmp) = LnValue2Key(occline)
    NGLn = NGLn -1
    call delete_Hash4G(kLn(occline),occline)
  else
    tmp = WLnKey2Value(NWLn)
    WLnKey2Value(NWLn) = 0
    WLnKey2Value(LnValue2Key(occline)) = tmp
    LnValue2Key(tmp) = LnValue2Key(occline)
    NWLn = NWLn -1
    call delete_Hash4W(kLn(occline),occline)
  endif

  NextLn(occline) = TailLn
  StatusLn(occline) = -1
  TailLn = occline


  if(TailLn == -1) then
    call LogFile%QuickLog("undo_insert_line error! tail=-1!", 'e')
    stop
  endif


  return
END SUBROUTINE undo_insert_line


!------------- insert a gamma to the link -------------------
SUBROUTINE insert_gamma(newgamma, isdelta, gx, gy, wx, wy, t1, t2, t3, dir, typ, stat, weigh)
  implicit none
  integer, intent(out) :: newgamma
  integer, intent(in) :: gx, gy, wx, wy, isdelta, dir, typ, stat
  complex*16, intent(in) :: weigh
  double precision, intent(in) :: t1, t2, t3

  newgamma = TailVertex
  TailVertex = NextVertex(TailVertex)
  if(StatusVertex(TailVertex)>=0) then
    call LogFile%QuickLog("update: "+str(imc)+" , insert_gamma error!!!", 'e')
    call print_config
    stop
  endif

  if(TailVertex == -1) then
    call LogFile%QuickLog("insert_gamma error! tail=-1!", 'e')
    call print_config
    stop
  endif
   
  NVertex = NVertex + 1
  VertexKey2Value(NVertex) = newgamma
  VertexValue2Key(newgamma) = NVertex

  IsDeltaVertex(newgamma) = isdelta
  GRVertex(1, newgamma) = gx
  GRVertex(2, newgamma) = gy
  WRVertex(1, newgamma) = wx
  WRVertex(2, newgamma) = wy
  TVertex(1, newgamma) = t1
  TVertex(2, newgamma) = t2
  TVertex(3, newgamma) = t3
  DirecVertex(newgamma) = dir
  TypeVertex(newgamma) = typ
  SpInVertex(1, newgamma) = typ
  SpInVertex(2, newgamma) = typ
  StatusVertex(newgamma) = stat
  WeightVertex(newgamma) = weigh
  return
END SUBROUTINE insert_gamma

!------------- delete a line from the link -------------------
SUBROUTINE delete_line(occline, knd)
  implicit none
  integer, intent(in) :: occline, knd
  integer :: tmp

  
  if(knd==1) then
    call delete_Hash4G(kLn(occline), occline)
  else
    call delete_Hash4W(kLn(occline), occline)
  endif

  if(StatusLn(occline)==-1) then
    call LogFile%QuickLog("update: "+str(imc)+" , delete_gamma error!!!", 'e')
    call print_config
    stop
  endif

  NextLn(occline) = TailLn
  StatusLn(occline) = -1
  TailLn = occline

  if(knd==1) then

    tmp = GLnKey2Value(NGLn)
    GLnKey2Value(NGLn) = 0
    GLnKey2Value(LnValue2Key(occline)) = tmp
    LnValue2Key(tmp) = LnValue2Key(occline)
    NGLn = NGLn -1
  else

    tmp = WLnKey2Value(NWLn)
    WLnKey2Value(NWLn) = 0
    WLnKey2Value(LnValue2Key(occline)) = tmp
    LnValue2Key(tmp) = LnValue2Key(occline)
    NWLn = NWLn -1
  endif

  if(TailLn == -1) then
    call LogFile%QuickLog("delete_gamma error! tail=-1!", 'e')
    stop
  endif
  return
END SUBROUTINE delete_line

!------------- insert a gamma to the link -------------------
SUBROUTINE undo_delete_line(newline, knd, stat, k)
  implicit none
  integer, intent(in) :: newline, knd, stat, k

  if(TailLn/=newline)    then
    call LogFile%QuickLog("update: "+str(imc)+"  ,undo_delete_line error!!!", 'e')
    call print_config
    stop
  endif
  StatusLn(newline) = stat
  TailLn = NextLn(newline)
  if(TailLn == -1) then
    call LogFile%QuickLog("undo_delete_line error! tail=-1!", 'e')
    call print_config
    stop
  endif

  if(knd==1) then
    NGLn = NGLn + 1
    GLnKey2Value(NGLn) = newline
    LnValue2Key(newline) = NGLn
    call add_Hash4G(k, newline)
  else
    NWLn = NWLn + 1
    WLnKey2Value(NWLn) = newline
    LnValue2Key(newline) = NWLn
    call add_Hash4W(k, newline)
  endif
  return
END SUBROUTINE undo_delete_line


!------------- delete a gamma from the link -------------------
SUBROUTINE delete_gamma(occgamma)
  implicit none
  integer, intent(in) :: occgamma
  integer :: tmp

  if(StatusVertex(occgamma)==-1) then
    call LogFile%QuickLog("update: "+str(imc)+" , delete_gamma error!!!", 'e')
    call print_config
    stop
  endif
  NextVertex(occgamma) = TailVertex
  StatusVertex(occgamma) = -1
  TailVertex = occgamma

  tmp = VertexKey2Value(NVertex)
  VertexKey2Value(NVertex) = 0
  VertexKey2Value(VertexValue2Key(occgamma)) = tmp
  VertexValue2Key(tmp) = VertexValue2Key(occgamma)
  NVertex = NVertex -1

  if(TailVertex == -1) then
    call LogFile%QuickLog("delete_gamma error! tail=-1!", 'e')
    stop
  endif
  return
END SUBROUTINE delete_gamma


SUBROUTINE undo_delete_gamma(newgamma)
  implicit none
  integer, intent(in) :: newgamma

  if(TailVertex/=newgamma)    then
    call LogFile%QuickLog("update: "+str(imc)+" , undo_delete_gamma error!!!", 'e')
    stop
  endif
  StatusVertex(newgamma) = 0

  TailVertex = NextVertex(newgamma)
  if(TailVertex == -1) then
    call LogFile%QuickLog("undo_delete_gamma error! tail=-1!", 'e')
    stop
  endif
   
  NVertex = NVertex + 1
  VertexKey2Value(NVertex) = newgamma
  VertexValue2Key(newgamma) = NVertex

  return
END SUBROUTINE undo_delete_gamma
!!=======================================================================
!!=======================================================================
!!=======================================================================



!======================= Complex Operations =============================
COMPLEX*16 FUNCTION d_times_cd(nd, ncd)
  implicit none 
  double precision, intent(in) :: nd
  complex*16, intent(in) :: ncd
  d_times_cd = dcmplx(nd*real(ncd), nd*dimag(ncd))
  return
END FUNCTION d_times_cd

!!=======================================================================
!!======================= Fourier Transformation ========================
!!=======================================================================
!MODULE test
    !implicit none
    !integer,parameter :: L(1)=4,L(2)=4
    !integer,parameter :: Omega=2
    !integer,parameter :: NtypeW=7
    !double precision :: WR(NtypeW,0:L(1)-1,0:L(2)-1,-Omega:Omega,-Omega:Omega)
!end module

!program main
    !use test
    !integer :: ix,iy
    !print *, "Start"
    !WR(:,:,:,:,:)=1
    !WR(1,0,0,:,:)=3
    !write(*,*) "Before"
    !do iy=0,L(2)-1
      !write(*,*) WR(1,:,iy,2,2)
    !enddo
    !call FFT(WR,NtypW,(2*Omega+1)**2,1)
    !print *, "First FFT"
    !do iy=0,L(2)-1
      !write(*,*) WR(1,:,iy,2,2)
    !enddo
    !call FFT(WR,NtypW,(2*Omega+1)**2,-1)
    !print *, "Second FFT"
    !do iy=0,L(2)-1
      !write(*,*) WR(1,:,iy,2,2)
    !enddo

!CONTAINS

SUBROUTINE transfer_W0_r(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    call FFT_r(W0PF,1,MxT,BackForth)
END SUBROUTINE

SUBROUTINE transfer_W0_t(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    call FFT_tau_single(W0PF,1,L(1)*L(2),BackForth)
END SUBROUTINE

SUBROUTINE transfer_Gam0_r(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer :: it1, it2

    call FFT_r(Gam0PF,1,MxT**2,BackForth)
end SUBROUTINE

SUBROUTINE transfer_Gam0_t(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer :: it1, it2

    if(BackForth/=-1) then
      do it2 = 0, MxT-1
        do it1 = 0, MxT-1
          Gam0PF(:,:,it1,it2) = Gam0PF(:,:,it1,it2)* cdexp(dcmplx(0.d0,-Pi*real(it1)/real(MxT)))
          Gam0PF(:,:,it1,it2) = Gam0PF(:,:,it1,it2)* cdexp(dcmplx(0.d0,-Pi*real(it2)/real(MxT)))
        enddo
      enddo

      call FFT_tau_double(Gam0PF,1,L(1)*L(2),BackForth)
    else if(BackForth ==-1) then
      call FFT_tau_double(Gam0PF,1,L(1)*L(2),BackForth)

      do it2 = 0, MxT-1
        do it1 = 0, MxT-1
          Gam0PF(:,:,it1,it2) = Gam0PF(:,:,it1,it2)* cdexp(dcmplx(0.d0,Pi*real(it1)/real(MxT)))
          Gam0PF(:,:,it1,it2) = Gam0PF(:,:,it1,it2)* cdexp(dcmplx(0.d0,Pi*real(it2)/real(MxT)))
        enddo
      enddo
    endif
END SUBROUTINE

SUBROUTINE transfer_G0(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer :: it
    if(BackForth/=-1) then
      do it = 0, MxT-1
        G0F(it) = G0F(it)* cdexp(dcmplx(0.d0,-Pi/real(MxT))*real(it))
      enddo

      call FFT_tau_single(G0F,1,1,BackForth)
    else if(BackForth ==-1) then
      call FFT_tau_single(G0F,1,1,BackForth)

      do it = 0, MxT-1
        G0F(it) = G0F(it)* cdexp(dcmplx(0.d0,Pi/real(MxT))*real(it))
      enddo
    endif
END SUBROUTINE

SUBROUTINE transfer_GamOrder1_t(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer ::it1, it2
    if(BackForth/=-1) then
      do it2 = 0, MxT-1
        do it1 = 0, MxT-1
          GamOrder1(:,it1,it2) = GamOrder1(:,it1,it2)* cdexp(dcmplx(0.d0,-Pi*real(it1)/real(MxT)))
          GamOrder1(:,it1,it2) = GamOrder1(:,it1,it2)* cdexp(dcmplx(0.d0,-Pi*real(it2)/real(MxT)))
        enddo
      enddo

      call FFT_tau_double(GamOrder1,NtypeGam,1,BackForth)
    else if(BackForth ==-1) then
      call FFT_tau_double(GamOrder1,NtypeGam,1,BackForth)

      do it2 = 0, MxT-1
        do it1 = 0, MxT-1
          GamOrder1(:,it1,it2) = GamOrder1(:,it1,it2)* cdexp(dcmplx(0.d0,Pi*real(it1)/real(MxT)))
          GamOrder1(:,it1,it2) = GamOrder1(:,it1,it2)* cdexp(dcmplx(0.d0,Pi*real(it2)/real(MxT)))
        enddo
      enddo
    endif
END SUBROUTINE


SUBROUTINE transfer_G_t(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer            :: it
    if(BackForth/=-1) then
      do it = 0, MxT-1
        G(:,it) = G(:,it)* cdexp(dcmplx(0.d0,-Pi*real(it)/real(MxT)))
      enddo

      call FFT_tau_single(G,NtypeG,1,BackForth)
    else if(BackForth ==-1) then
      call FFT_tau_single(G,NtypeG,1,BackForth)

      do it = 0, MxT-1
        G(:,it) = G(:,it)* cdexp(dcmplx(0.d0, Pi*real(it)/real(MxT)))
      enddo
    endif
END SUBROUTINE

SUBROUTINE transfer_W_r(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer :: ix, iy

    call FFT_r(W,NtypeW,MxT,BackForth)

END SUBROUTINE

SUBROUTINE transfer_W_t(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    call FFT_tau_single(W,NtypeW,L(1)*L(2),BackForth)
END SUBROUTINE

SUBROUTINE transfer_Gam_r(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    call FFT_r(Gam,NtypeGam,MxT**2,BackForth)
END SUBROUTINE

SUBROUTINE transfer_Gam_t(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer :: it1, it2
    if(BackForth/=-1) then
      do it2 = 0, MxT-1
        do it1 = 0, MxT-1
          Gam(:,:,:,it1,it2) = Gam(:,:,:,it1,it2)* cdexp(dcmplx(0.d0,-Pi*real(it1)/real(MxT)))
          Gam(:,:,:,it1,it2) = Gam(:,:,:,it1,it2)* cdexp(dcmplx(0.d0,-Pi*real(it2)/real(MxT)))
        enddo
      enddo

      call FFT_tau_double(Gam,NtypeGam,L(1)*L(2),BackForth)
    else if(BackForth ==-1) then
      call FFT_tau_double(Gam,NtypeGam,L(1)*L(2),BackForth)

      do it2 = 0, MxT-1
        do it1 = 0, MxT-1
          Gam(:,:,:,it1,it2) = Gam(:,:,:,it1,it2)* cdexp(dcmplx(0.d0,Pi*real(it1)/real(MxT)))
          Gam(:,:,:,it1,it2) = Gam(:,:,:,it1,it2)* cdexp(dcmplx(0.d0,Pi*real(it2)/real(MxT)))
        enddo
      enddo
    endif
END SUBROUTINE

SUBROUTINE transfer_Chi_r(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    call FFT_r(Chi,1,MxT,BackForth)
END SUBROUTINE

SUBROUTINE transfer_Chi_t(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    call FFT_tau_single(Chi,1,L(1)*L(2),BackForth)
END SUBROUTINE

SUBROUTINE transfer_Sigma_t(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer       :: it

    if(BackForth/=-1) then
      do it = 0, MxT-1
        Sigma(it) = Sigma(it)* cdexp(dcmplx(0.d0,-Pi*real(it)/real(MxT)))
      enddo
      call FFT_tau_single(Sigma,1,1,BackForth)
    else if(BackForth==-1) then
      call FFT_tau_single(Sigma,1,1,BackForth)
      do it = 0, MxT-1
        Sigma(it) = Sigma(it)* cdexp(dcmplx(0.d0, Pi*real(it)/real(MxT)))
      enddo
    endif
    return
END SUBROUTINE

!SUBROUTINE transfer_Polar_r(BackForth)
    !implicit none
    !integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    !call FFT_r(Polar,1,MxT,BackForth)
!END SUBROUTINE

!SUBROUTINE transfer_Polar_t(BackForth)
    !implicit none
    !integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    !call FFT_tau_single(Polar,1,L(1)*L(2),BackForth)
!END SUBROUTINE

SUBROUTINE transfer_r(Backforth)
  implicit none
  integer,intent(in) :: Backforth

  call transfer_W_r(Backforth)
  call transfer_Gam_r(Backforth)
  return
END SUBROUTINE

SUBROUTINE transfer_t(Backforth)
  implicit none
  integer,intent(in) :: Backforth

  call transfer_G_t(Backforth)
  call transfer_W_t(Backforth)
  call transfer_Gam_t(Backforth)
  return
END SUBROUTINE


SUBROUTINE FFT_r(XR,Ntype,Nz,BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer,intent(in) :: Nz
    integer,intent(in) :: Ntype
    complex*16    :: XR(Ntype,0:L(1)-1,0:L(2)-1,0:Nz-1)
    integer :: Power,Noma
    integer :: ix,iy,iz,it
    double precision,allocatable ::   Real1(:), Im1(:)

    do iy=0,L(2)-1
      do ix=0,L(1)-1
        if(ix>L(1)/2 .and. iy<=L(2)/2) then
          XR(:,ix,iy,:)=XR(:,L(1)-ix,iy,:)
        else if(ix<=L(1)/2 .and. iy>L(2)/2) then
          XR(:,ix,iy,:)=XR(:,ix,L(2)-iy,:)
        else if(ix>L(1)/2 .and. iy>L(2)/2)  then
          XR(:,ix,iy,:)=XR(:,L(1)-ix,L(2)-iy,:)
        endif
      enddo
    enddo

    !do FFT in x direction
    Power=log(L(1)*1.d0)/log(2.d0)+1.d-14  
    Noma=2**power
    allocate(Real1(0:Noma-1),Im1(0:Noma-1))
    do iz=0,Nz-1
      do iy=0,L(2)-1
        do it=1,Ntype
           Real1(:)=real(XR(it,:,iy,iz))
           Im1(:)=dimag(XR(it,:,iy,iz))
           call sffteu(Real1, Im1, Noma, Power, BackForth)
           XR(it,:,iy,iz)=dcmplx(Real1(:), Im1(:))
         enddo
      enddo
    enddo
    deallocate(Real1,Im1)

    !do FFT in y direction
    power=log(L(2)*1.d0)/log(2.d0)+1.d-14  
    Noma=2**Power
    allocate(Real1(0:Noma-1),Im1(0:Noma-1))
    do iz=0,Nz-1
      do ix=0,L(1)-1
        do it=1,Ntype
           Real1(:)=real(XR(it,ix,:,iz))
           Im1(:)=dimag(XR(it,ix,:,iz))
           call sffteu(Real1, Im1, Noma, Power, BackForth)
           XR(it,ix,:,iz)=dcmplx(Real1(:), Im1(:))
         enddo
      enddo
    enddo
    deallocate(Real1,Im1)

end SUBROUTINE FFT_r

SUBROUTINE FFT_tau_single(XR,Ntype,Nz,BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer,intent(in) :: Nz
    integer,intent(in) :: Ntype
    complex(kind=8)    :: XR(Ntype,0:Nz-1,0:MxT-1)
    integer :: Power,Noma
    integer :: itau,iomega,iz,it
    double precision,allocatable ::   Real1(:), Im1(:)

    !do FFT in tau 
    Power=log(MxT*1.d0)/log(2.d0)+1.d-14  
    Noma=2**power
    allocate(Real1(0:Noma-1),Im1(0:Noma-1))
    do iz=0,Nz-1
      do it=1,Ntype
        Real1(:)=real(XR(it,iz,:))
        Im1(:)=dimag(XR(it,iz,:))
        call sffteu(Real1, Im1, Noma, Power, BackForth)
        XR(it,iz,:)=dcmplx(Real1(:), Im1(:))
       enddo
    enddo
    deallocate(Real1,Im1)

end SUBROUTINE FFT_tau_single
    
SUBROUTINE FFT_tau_double(XR,Ntype,Nz,BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer,intent(in) :: Nz
    integer,intent(in) :: Ntype
    complex(kind=8)    :: XR(Ntype,0:Nz-1,0:MxT-1,0:MxT-1)
    integer :: Power,Noma
    integer :: itau1,itau2,iz,it
    double precision,allocatable ::   Real1(:), Im1(:)

    !do FFT in tau1 direction
    Power=log(MxT*1.d0)/log(2.d0)+1.d-14  
    Noma=2**power
    allocate(Real1(0:Noma-1),Im1(0:Noma-1))
    do itau2=0,MxT-1
      do iz=0,Nz-1
        do it=1,Ntype
          Real1(:)=real(XR(it,iz,:,itau2))
          Im1(:)=dimag(XR(it,iz,:,itau2))
          call sffteu(Real1, Im1, Noma, Power, BackForth)
          XR(it,iz,:,itau2)=dcmplx(Real1(:),Im1(:))
        enddo
      enddo
    enddo
    deallocate(Real1,Im1)

    !do FFT in tau2 direction
    Power=log(MxT*1.d0)/log(2.d0)+1.d-14  
    Noma=2**power
    allocate(Real1(0:Noma-1),Im1(0:Noma-1))
    do itau1=0,MxT-1
      do iz=0,Nz-1
        do it=1,Ntype
            Real1(:)=real(XR(it,iz,itau1,:))
            Im1(:)=dimag(XR(it,iz,itau1,:))
            call sffteu(Real1, Im1, Noma, Power, BackForth)
            XR(it,iz,itau1,:)=dcmplx(Real1(:),Im1(:))
        enddo
      enddo
    enddo
    deallocate(Real1,Im1)

end SUBROUTINE FFT_tau_double
    
!-------------------------------------------------------------c
!                                                             c
!  Subroutine sffteu( x, y, n, m, itype )                     c
!                                                             c
!  This routine is a slight modification of a complex split   c
!  radix FFT routine presented by C.S. Burrus.  The original  c
!  program header is shown below.                             c
!                                                             c
!  Arguments:                                                 c
!     x - real array containing real parts of transform       c
!              sequence (in/out)                              c
!     y - real array containing imag parts of transform       c
!              sequence (in/out)                              c
!     n - integer length of transform (in)                    c
!     m - integer such that n = 2**m  (in)                    c
!     itype - integer job specifier (in)                      c
!              itype .ne. -1 --> foward transform             c
!              itype .eq. -1 --> backward transform           c
!                                                             c
!  The forward transform computes                             c
!     X(k) = sum_{j=0}^{N-1} x(j)*exp(-2ijk*pi/N)             c
!                                                             c
!  The backward transform computes                            c
!     x(j) = (1/N) * sum_{k=0}^{N-1} X(k)*exp(2ijk*pi/N)      c
!                                                             c
!  Requires standard FORTRAN functions - sin, cos             c
!                                                             c
!  Steve Kifowit, 9 July 1997                                 c
!                                                             c
!-------------------------------------------------------------C
!  A Duhamel-Hollman Split-Radix DIF FFT                      C
!  Reference:  Electronics Letters, January 5, 1984           C
!  Complex input and output in data arrays X and Y            C
!  Length is N = 2**M                                         C
!                                                             C
!  C.S. Burrus          Rice University         Dec 1984      C
!-------------------------------------------------------------C

subroutine SFFTEU( X, Y, N, M, ITYPE )
integer  N, M, ITYPE
real(kind=8) ::  X(*), Y(*)
integer  I, J, K, N1, N2, N4, IS, ID, I0, I1, I2, I3
real(kind=8) ::  E, A, A3, CC1, SS1, CC3, SS3
real(kind=8) ::  R1, R2, S1, S2, S3, XT
!      INTRINSIC  SIN, COS
real(kind=8), parameter :: TWOPI = 6.2831853071795864769 

if ( N .eq. 1 ) return

if ( ITYPE .eq. -1 ) then
   do I = 1, N
      Y(I) = - Y(I)
   end do
end if

N2 = 2 * N
do K = 1, M-1
   N2 = N2 / 2
   N4 = N2 / 4
   E = TWOPI / N2
   A = 0.0
   do J = 1, N4
      A3 = 3 * A
      CC1 = DCOS( A )
      SS1 = DSIN( A )
      CC3 = DCOS( A3 )
      SS3 = DSIN( A3 )
      A = J * E
      IS = J
      ID = 2 * N2
40        do I0 = IS, N-1, ID
         I1 = I0 + N4
         I2 = I1 + N4
         I3 = I2 + N4
         R1 = X(I0) - X(I2)
         X(I0) = X(I0) + X(I2)
         R2 = X(I1) - X(I3)
         X(I1) = X(I1) + X(I3)
         S1 = Y(I0) - Y(I2)
         Y(I0) = Y(I0) + Y(I2)
         S2 = Y(I1) - Y(I3)
         Y(I1) = Y(I1) + Y(I3)
         S3 = R1 - S2
         R1 = R1 + S2
         S2 = R2 - S1
         R2 = R2 + S1
         X(I2) = R1 * CC1 - S2 * SS1
         Y(I2) = - S2 * CC1 - R1 * SS1
         X(I3) = S3 * CC3 + R2 * SS3
         Y(I3) = R2 * CC3 - S3 * SS3
      end do
      IS = 2 * ID - N2 + J
      ID = 4 * ID
      if ( IS .lt. N ) goto 40
   end do
end do

!--------LAST STAGE, LENGTH-2 BUTTERFLY ----------------------C

IS = 1
ID = 4
50  do I0 = IS, N, ID
   I1 = I0 + 1
   R1 = X(I0)
   X(I0) = R1 + X(I1)
   X(I1) = R1 - X(I1)
   R1 = Y(I0)
   Y(I0) = R1 + Y(I1)
   Y(I1) = R1 - Y(I1)
end do
IS = 2 * ID - 1
ID = 4 * ID
if ( IS .lt. N ) goto 50

!-------BIT REVERSE COUNTER-----------------------------------C

100 J = 1
N1 = N - 1
do I = 1, N1
   if ( I .ge. J ) goto 101
   XT = X(J)
   X(J) = X(I)
   X(I) = XT
   XT = Y(J)
   Y(J) = Y(I)
   Y(I) = XT
101    K = N / 2
102    if ( K .ge. J ) goto 103
   J = J - K
   K = K / 2
   goto 102
103    J = J + K
end do

if ( ITYPE .eq. -1 ) then
   do I = 1, N
      X(I) = X(I) / N
      Y(I) = - Y(I) / N
   end do
end if

return

end  subroutine SFFTEU
!end program
!!=======================================================================
!!=======================================================================
!!=======================================================================



!!=======================================================================
!!================ RANDOM NUMBER GENERATORS =============================
!!=======================================================================

!---------- shift register random number generator ------
SUBROUTINE set_RNG
  implicit none
  integer :: i_r, k_r, k1_r, iseed

  nrannr = mxrn
  iseed  = iabs(Seed)+1
  k_r    = 3**18+2*iseed
  k1_r   = 1313131*iseed
  k1_r   = k1_r-(k1_r/mod2)*mod2

  do i_r = 1, len1
    k_r  = k_r *mult
    k1_r = k1_r*mul2
    k1_r = k1_r-(k1_r/mod2)*mod2
    ir1(i_r) = k_r+k1_r*8193
  enddo

  do i_r = 1, len2
    k_r  = k_r *mult
    k1_r = k1_r*mul2
    k1_r = k1_r-(k1_r/mod2)*mod2
    ir2(i_r) = k_r+k1_r*4099
  enddo

  do i_r = 1, len1
    inxt1(i_r) = i_r+1
  enddo
  inxt1(len1) = 1
  ipnt1 = 1
  ipnf1 = ifd1+1

  do i_r = 1, len2
    inxt2(i_r) = i_r+1
  enddo
  inxt2(len2) = 1
  ipnt2 = 1
  ipnf2 = ifd2 + 1
END SUBROUTINE set_RNG 

!---------- calculate next random number --------------
DOUBLE PRECISION FUNCTION rn()
  implicit none
  integer   :: i_r, l_r, k_r
  nrannr = nrannr +1
  if(nrannr>=mxrn) then
    nrannr = 1
    do i_r= 1, mxrn
      l_r = ieor(ir1(ipnt1),ir1(ipnf1))
      k_r = ieor(ir2(ipnt2),ir2(ipnf2))
      irn(i_r) = ieor(k_r,l_r)
      ir1(ipnt1)=l_r
      ipnt1 = inxt1(ipnt1)
      ipnf1 = inxt1(ipnf1)
      ir2(ipnt2) = k_r
      ipnt2 = inxt2(ipnt2)
      ipnf2 = inxt2(ipnf2)
    enddo
  endif 
  rn = irn(nrannr)*tm32+0.5d0
END FUNCTION rn

!===================================================================



!==============Trace elapsed time ==================================
!! THIS IS PROJECT-INDEPENDENT 
SUBROUTINE set_time_elapse
  implicit none
  !-- read and calculate time (in seconds) -------------------------
  call date_and_time(date, time, zone, tval)
  t_curr = tval(5)*3600.d0+tval(6)*60.d0+tval(7)+tval(8)*0.001d0 
  h_curr = tval(5)
  t_prev = t_curr
  h_prev = h_curr
  return
END SUBROUTINE set_time_elapse
  

!==============Trace elapsed time ==================================
!! THIS IS PROJECT-INDEPENDENT 
SUBROUTINE time_elapse
  implicit none
  
  !-- read and calculate time (in seconds) -------------------------
  call date_and_time(date, time, zone, tval)
  t_curr = tval(5)*3600.d0+tval(6)*60.d0+tval(7)+tval(8)*0.001d0 
  h_curr = tval(5)

  t_elap = t_curr-t_prev
  if(h_curr<h_prev) t_elap = t_elap+24*3600.d0
  t_prev = t_curr
  h_prev = h_curr 
  return
END SUBROUTINE time_elapse



!============   initialize diagrams =================================
SUBROUTINE first_order_diagram
  implicit none
  integer :: i, dr0(2)
  double precision :: ratio
  complex*16 :: Anew

  !-------------- 1-order diagram ------------------------
  Order = 1
  ! the index of measuring gamma
  NGLn = 4;  NWLn = 2;  NVertex = 4
  ! the number of glines, wlines, gamma
  MeasureGam = 1
  ! number of fermi loops
  SignFermiLoop = 1.d0
  ! the phase of the diagram
  Phase = (1.d0, 0.d0)
  ! if Ira and Masha are present
  IsWormPresent = .false.

  !status for lines: 0: normal; 1: with measuring Gamma
  StatusLn(1) = 1
  StatusLn(2) = 1
  StatusLn(3) = 1
  StatusLn(4) = 0
  StatusLn(5) = 0
  StatusLn(6) = 0
  TailLn = 7

  GLnKey2Value(1)= 1;     GLnKey2Value(2)= 2;     GLnKey2Value(3)= 4
  GLnKey2Value(4)= 5
  WLnKey2Value(1)= 3;     WLnKey2Value(2)= 6

  LnValue2Key(1)  = 1;     LnValue2Key(2)  = 2;     LnValue2Key(3)  = 1
  LnValue2Key(4)  = 3;     LnValue2Key(5)  = 4;     LnValue2Key(6)  = 2


  !kind of lines: 1: gline;  2: wline
  KindLn(1) = 1;     KindLn(2) = 1;     KindLn(3) = 2
  KindLn(4) = 1;     KindLn(5) = 1;     KindLn(6) = 2

  !type of glines: 1: spin up; 2: spin down
  TypeLn(1) = 1;     TypeLn(2) = 1
  TypeLn(4) = 1;     TypeLn(5) = 1

  !IsDelta: 1: delta function(W); 0: normal function(W); -1: G
  IsDeltaLn(1:6) = -1
  IsDeltaLn(3) = 0
  IsDeltaLn(6) = 0

  !type of gamma inside spins: 1: spin up; 2: spin down
  SpInVertex(:, :) = 1

  !type of Gamma: 1: gin, 2: gout, 3: win, 4: wout
  TypeVertex(1) = 1
  TypeVertex(2) = 1
  TypeVertex(3) = 1
  TypeVertex(4) = 1

  !type of wlines
  TypeLn(3) = TypeGam2W(TypeVertex(1), TypeVertex(2))
  TypeLn(6) = TypeGam2W(TypeVertex(3), TypeVertex(4))


  ! the momentum on line 1, 2, and 3
  kLn(1) = 100
  call add_Hash4G(kLn(1),1)
  kLn(2) = 32
  call add_Hash4G(kLn(2),2)

  kLn(3) = add_k(kLn(2), -kLn(1))
  call add_Hash4W(kLn(3),3)

  kLn(4) = 23
  call add_Hash4G(kLn(4),4)


  kLn(5) = add_k(kLn(3), kLn(4))
  call add_Hash4G(kLn(5),5)

  kLn(6) = add_k(kLn(1), -kLn(4))
  call add_Hash4W(kLn(6),6)


  ! NeighLn(i, j): the ith neighbor gamma of line j
  NeighLn(1,1) = 1;                NeighLn(2,1) = 3
  NeighLn(1,2) = 4;                NeighLn(2,2) = 1
  NeighLn(1,3) = 1;                NeighLn(2,3) = 2
  NeighLn(1,4) = 3;                NeighLn(2,4) = 2
  NeighLn(1,5) = 2;                NeighLn(2,5) = 4
  NeighLn(1,6) = 3;                NeighLn(2,6) = 4

  StatusVertex(1) = 1
  StatusVertex(2) = 0
  StatusVertex(3) = 0
  StatusVertex(4) = 0

  do i = 1, 4
    VertexKey2Value(i) = i
    VertexValue2Key(i)  = i
  enddo
  TailVertex = 5

  !IsDeltaVertex: 1: delta function(Gam); 0: normal function(Gam)
  IsDeltaVertex(1:4) = 1
  !IsDeltaVertex(1:4) = 0
  !IsDeltaVertex(1) = 1

  ! the site variables of Gamma 1, 2, 3 and 4
  GRVertex(1, 1:4) = 0;            GRVertex(2, 1:4) = 0
  WRVertex(1, 1:4) = 0;            WRVertex(2, 1:4) = 0

  ! the time variables of Gamma 1, 2, 3 and 4
  TVertex(1, 1)=0.d0; TVertex(2, 1)=0.d0; TVertex(3, 1)=0.d0
  TVertex(1, 2)=0.d0; TVertex(2, 2)=0.d0; TVertex(3, 2)=0.d0
  TVertex(1, 3)=0.2d0; TVertex(2, 3)=0.2d0; TVertex(3, 3)=0.2d0
  TVertex(1, 4)=0.3d0; TVertex(2, 4)=0.3d0; TVertex(3, 4)=0.3d0

  ! Direction of Gamma: 1: left of W;  2: right of W
  DirecVertex(1) = 1;                DirecVertex(2) = 2
  DirecVertex(3) = 1;                DirecVertex(4) = 2

  ! NeighVertex(i, j): the ith neighbor line of gamma j
  NeighVertex(1,1) = 2;        NeighVertex(2,1) = 1;        NeighVertex(3,1) = 3
  NeighVertex(1,2) = 4;        NeighVertex(2,2) = 5;        NeighVertex(3,2) = 3
  NeighVertex(1,3) = 1;        NeighVertex(2,3) = 4;        NeighVertex(3,3) = 6
  NeighVertex(1,4) = 5;        NeighVertex(2,4) = 2;        NeighVertex(3,4) = 6

  ! weights for lines and vertexes
  WeightLn(1) = weight_gline(StatusLn(1),TVertex(2,NeighLn(2,1))-TVertex(1,NeighLn(1,1)),TypeLn(1))
  WeightLn(2) = weight_gline(StatusLn(2),TVertex(2,NeighLn(2,2))-TVertex(1,NeighLn(1,2)),TypeLn(2))
  WeightLn(4) = weight_gline(StatusLn(4),TVertex(2,NeighLn(2,4))-TVertex(1,NeighLn(1,4)),TypeLn(4))
  WeightLn(5) = weight_gline(StatusLn(5),TVertex(2,NeighLn(2,5))-TVertex(1,NeighLn(1,5)),TypeLn(5))

  dr0(:) = 0
  WeightLn(3) = weight_wline(StatusLn(3),IsDeltaLn(3),dr0,TVertex(3,NeighLn(2,3))-TVertex(3,NeighLn(1,3)),TypeLn(3))
  WeightLn(6) = weight_wline(StatusLn(6),IsDeltaLn(6),dr0,TVertex(3,NeighLn(2,6))-TVertex(3,NeighLn(1,6)),TypeLn(6))

  WeightVertex(1) = weight_vertex(StatusVertex(1), IsDeltaVertex(1), dr0, 0.d0, 0.d0, TypeVertex(1))
  WeightVertex(2) = weight_vertex(StatusVertex(2), IsDeltaVertex(2), dr0, 0.d0, 0.d0, TypeVertex(2))
  WeightVertex(3) = weight_vertex(StatusVertex(3), IsDeltaVertex(3), dr0, 0.d0, 0.d0, TypeVertex(3))
  WeightVertex(4) = weight_vertex(StatusVertex(4), IsDeltaVertex(4), dr0, 0.d0, 0.d0, TypeVertex(4))


  !ratio = (1.d0/Beta)**Order *SignFermiLoop
  ratio = (-1.d0)**Order *SignFermiLoop
  Anew = d_times_cd(ratio, WeightLn(1)*WeightLn(2)*WeightLn(3)*WeightLn(4)* &
    & WeightLn(5)*WeightLn(6)*WeightVertex(1)*WeightVertex(2)*WeightVertex(3)* &
    & WeightVertex(4))

  WeightCurrent = abs(Anew)
  Phase = Anew/WeightCurrent
  !-------------------------------------------------------
end SUBROUTINE first_order_diagram




SUBROUTINE first_order_diagram_with_bubble
  implicit none
  integer :: i, dr0(2)
  double precision :: ratio
  complex*16 :: Anew

  !-------------- 1-order diagram with bubble ------------
  Order = 1
  ! the index of measuring gamma
  NGLn = 4;  NWLn = 2;  NVertex = 4
  ! the number of glines, wlines, gamma
  MeasureGam = 1
  ! number of fermi loops
  SignFermiLoop = -1.d0
  ! the phase of the diagram
  Phase = (1.d0, 0.d0)
  ! if Ira and Masha are present
  IsWormPresent = .false.

  !status for lines: 0: normal; 1: with measuring Gamma
  StatusLn(1) = 1
  StatusLn(2) = 1
  StatusLn(3) = 1
  StatusLn(4) = 0
  StatusLn(5) = 0
  StatusLn(6) = 0
  TailLn = 7

  GLnKey2Value(1)= 1;     GLnKey2Value(2)= 2;     GLnKey2Value(3)= 4
  GLnKey2Value(4)= 5
  WLnKey2Value(1)= 3;     WLnKey2Value(2)= 6

  LnValue2Key(1)  = 1;     LnValue2Key(2)  = 2;     LnValue2Key(3)  = 1
  LnValue2Key(4)  = 3;     LnValue2Key(5)  = 4;     LnValue2Key(6)  = 2


  !kind of lines: 1: gline;  2: wline
  KindLn(1) = 1;     KindLn(2) = 1;     KindLn(3) = 2
  KindLn(4) = 1;     KindLn(5) = 1;     KindLn(6) = 2

  !type of glines: 1: spin up; 2: spin down
  TypeLn(1) = 1;     TypeLn(2) = 1
  TypeLn(4) = 1;     TypeLn(5) = 1

  !IsDelta: 1: delta function(W); 0: normal function(W); -1: G
  IsDeltaLn(1:6) = -1
  IsDeltaLn(3) = 0
  IsDeltaLn(6) = 0

  !type of gamma inside spins: 1: spin up; 2: spin down
  SpInVertex(:, :) = 1

  !type of Gamma: 1: gin, 2: gout, 3: win, 4: wout
  TypeVertex(1) = 1
  TypeVertex(2) = 1
  TypeVertex(3) = 1
  TypeVertex(4) = 1

  !type of wlines
  TypeLn(3) = TypeGam2W(TypeVertex(1), TypeVertex(2))
  TypeLn(6) = TypeGam2W(TypeVertex(3), TypeVertex(4))


  ! the momentum on line 1, 2, and 3
  !TODO: UPDATE THE HASH TABLE
  !kLn(1) = 100
  !Hash4G(kLn(1)) = 1
  !kLn(2) = 32
  !Hash4G(kLn(2)) = 1
  !kLn(3) = add_k(kLn(1), -kLn(1))
  !Hash4W(abs(kLn(3))) = 1
  !kLn(4) = 23
  !Hash4G(kLn(4)) = 1
  !kLn(5) = add_k(kLn(4), kLn(3))
  !Hash4G(kLn(5)) = 1
  !kLn(6) = add_k(kLn(4), -kLn(2))
  !Hash4W(abs(kLn(6))) = 1


  ! NeighLn(i, j): the ith neighbor gamma of line j
  NeighLn(1,1) = 1;                NeighLn(2,1) = 1
  NeighLn(1,2) = 3;                NeighLn(2,2) = 4
  NeighLn(1,3) = 1;                NeighLn(2,3) = 2
  NeighLn(1,4) = 2;                NeighLn(2,4) = 3
  NeighLn(1,5) = 4;                NeighLn(2,5) = 2
  NeighLn(1,6) = 3;                NeighLn(2,6) = 4

  StatusVertex(1) = 1
  StatusVertex(2) = 0
  StatusVertex(3) = 0
  StatusVertex(4) = 0

  do i = 1, 4
    VertexKey2Value(i) = i
    VertexValue2Key(i)  = i
  enddo
  TailVertex = 5

  !IsDeltaVertex: 1: delta function(Gam); 0: normal function(Gam)
  IsDeltaVertex(1:4) = 1

  ! the site variables of Gamma 1, 2, 3 and 4
  GRVertex(1, 1:4) = 0;            GRVertex(2, 1:4) = 0
  WRVertex(1, 1:4) = 0;            WRVertex(2, 1:4) = 0

  ! the time variables of Gamma 1, 2, 3 and 4
  TVertex(1, 1)=0.d0; TVertex(2, 1)=0.d0; TVertex(3, 1)=0.d0
  TVertex(1, 2)=0.1; TVertex(2, 2)=0.1; TVertex(3, 2)=0.1
  TVertex(1, 3)=0.2; TVertex(2, 3)=0.2; TVertex(3, 3)=0.2
  TVertex(1, 4)=0.3; TVertex(2, 4)=0.3; TVertex(3, 4)=0.3

  ! Direction of Gamma: 1: left of W;  2: right of W
  DirecVertex(1) = 1;                DirecVertex(2) = 2
  DirecVertex(3) = 1;                DirecVertex(4) = 2

  ! NeighVertex(i, j): the ith neighbor line of gamma j
  NeighVertex(1,1) = 1;        NeighVertex(2,1) = 1;        NeighVertex(3,1) = 3
  NeighVertex(1,2) = 5;        NeighVertex(2,2) = 4;        NeighVertex(3,2) = 3
  NeighVertex(1,3) = 4;        NeighVertex(2,3) = 2;        NeighVertex(3,3) = 6
  NeighVertex(1,4) = 2;        NeighVertex(2,4) = 5;        NeighVertex(3,4) = 6

  ! weights for lines and vertexes
  WeightLn(1) = weight_gline(StatusLn(1),TVertex(2,NeighLn(2,1))-TVertex(1,NeighLn(1,1)),TypeLn(1))
  WeightLn(2) = weight_gline(StatusLn(2),TVertex(2,NeighLn(2,2))-TVertex(1,NeighLn(1,2)),TypeLn(2))
  WeightLn(4) = weight_gline(StatusLn(4),TVertex(2,NeighLn(2,4))-TVertex(1,NeighLn(1,4)),TypeLn(4))
  WeightLn(5) = weight_gline(StatusLn(5),TVertex(2,NeighLn(2,5))-TVertex(1,NeighLn(1,5)),TypeLn(5))

  dr0(:) = 0
  WeightLn(3) = weight_wline(StatusLn(3),IsDeltaLn(3),dr0,TVertex(3,NeighLn(2,3))-TVertex(3,NeighLn(1,3)),TypeLn(3))
  WeightLn(6) = weight_wline(StatusLn(6),IsDeltaLn(6),dr0,TVertex(3,NeighLn(2,6))-TVertex(3,NeighLn(1,6)),TypeLn(6))

  WeightVertex(1) = weight_vertex(StatusVertex(1), 1, dr0, 0.d0, 0.d0, TypeVertex(1))
  WeightVertex(2) = weight_vertex(StatusVertex(2), 1, dr0, 0.d0, 0.d0, TypeVertex(2))
  WeightVertex(3) = weight_vertex(StatusVertex(3), 1, dr0, 0.d0, 0.d0, TypeVertex(3))
  WeightVertex(4) = weight_vertex(StatusVertex(4), 1, dr0, 0.d0, 0.d0, TypeVertex(4))


  !ratio = (1.d0/Beta)**Order *SignFermiLoop
  ratio = (-1.d0)**Order *SignFermiLoop
  Anew = d_times_cd(ratio, WeightLn(1)*WeightLn(2)*WeightLn(3)*WeightLn(4)* &
    & WeightLn(5)*WeightLn(6)*WeightVertex(1)*WeightVertex(2)*WeightVertex(3)* &
    & WeightVertex(4))

  WeightCurrent = abs(Anew)
  Phase = Anew/WeightCurrent
  !-------------------------------------------------------
end SUBROUTINE first_order_diagram_with_bubble




SUBROUTINE second_order_diagram
  implicit none
  integer :: i, dr0(2)
  double precision :: ratio
  complex*16 :: Anew

  !-------------- 2-order diagram ------------------------
  Order = 2
  ! the index of measuring gamma
  NGLn = 6;  NWLn = 3;  NVertex = 6
  ! the number of glines, wlines, gamma
  MeasureGam = 1
  ! sign of fermi loops
  SignFermiLoop = -1.d0
  ! the phase of the diagram
  Phase = (1.d0, 0.d0)
  ! if Ira and Masha are present
  IsWormPresent = .false.

  !status for lines: 0: normal; 1: with measuring Gamma
  StatusLn(1) = 1
  StatusLn(2) = 1
  StatusLn(3) = 1
  StatusLn(4) = 0
  StatusLn(5) = 0
  StatusLn(6) = 0
  StatusLn(7) = 0
  StatusLn(8) = 0
  StatusLn(9) = 0
  TailLn = 10

  GLnKey2Value(1)= 1;     GLnKey2Value(2)= 2;     GLnKey2Value(3)= 4
  GLnKey2Value(4)= 5;     GLnKey2Value(5)= 7;     GLnKey2Value(6)= 8
  WLnKey2Value(1)= 3;     WLnKey2Value(2)= 6;     WLnKey2Value(3)= 9

  LnValue2Key(1)  = 1;     LnValue2Key(2)  = 2;     LnValue2Key(3)  = 1
  LnValue2Key(4)  = 3;     LnValue2Key(5)  = 4;     LnValue2Key(6)  = 2
  LnValue2Key(7)  = 5;     LnValue2Key(8)  = 6;     LnValue2Key(9)  = 3


  !kind of lines: 1: gline;  2: wline
  KindLn(1) = 1;     KindLn(2) = 1;     KindLn(3) = 2
  KindLn(4) = 1;     KindLn(5) = 1;     KindLn(6) = 2
  KindLn(7) = 1;     KindLn(8) = 1;     KindLn(9) = 2

  !type of glines: 1: spin up; 2: spin down
  TypeLn(1) = 1;     TypeLn(2) = 1
  TypeLn(4) = 1;     TypeLn(5) = 1
  TypeLn(7) = 1;     TypeLn(8) = 1

  !IsDelta: 1: delta function(W); 0: normal function(W); -1: G
  IsDeltaLn(1:9) = -1
  IsDeltaLn(3) = 0
  IsDeltaLn(6) = 0
  IsDeltaLn(9) = 0

  !type of gamma inside spins: 1: spin up; 2: spin down
  SpInVertex(:, :) = 1

  !type of Gamma: 1: gin, 2: gout, 3: win, 4: wout
  TypeVertex(:) = 1

  !type of wlines
  TypeLn(3) = TypeGam2W(TypeVertex(1), TypeVertex(2))
  TypeLn(6) = TypeGam2W(TypeVertex(3), TypeVertex(4))
  TypeLn(9) = TypeGam2W(TypeVertex(5), TypeVertex(6))


  ! the momentum on line 1, 2, and 3
  !TODO: UPDATE THE HASH TABLE
  !kLn(1) = 10
  !Hash4G(kLn(1)) = 1
  !kLn(2) = 21
  !Hash4G(kLn(2)) = 1
  !kLn(3) = add_k(kLn(1), -kLn(2))
  !Hash4W(abs(kLn(3))) = 1
  !kLn(4) = 33
  !Hash4G(kLn(4)) = 1
  !kLn(5) = add_k(kLn(4), -kLn(3))
  !Hash4G(kLn(5)) = 1
  !kLn(6) = 46
  !Hash4W(abs(kLn(6))) = 1
  !kLn(7) = add_k(kLn(4), -kLn(6))
  !Hash4G(kLn(7)) = 1
  !kLn(8) = add_k(kLn(2), kLn(6))
  !Hash4G(kLn(8)) = 1
  !kLn(9) = add_k(kLn(1), -kLn(8))
  !Hash4W(abs(kLn(9))) = 1


  ! NeighLn(i, j): the ith neighbor gamma of line j
  NeighLn(1,1) = 6;                NeighLn(2,1) = 1
  NeighLn(1,2) = 1;                NeighLn(2,2) = 4
  NeighLn(1,4) = 2;                NeighLn(2,4) = 3
  NeighLn(1,5) = 5;                NeighLn(2,5) = 2
  NeighLn(1,7) = 3;                NeighLn(2,7) = 5
  NeighLn(1,8) = 4;                NeighLn(2,8) = 6

  NeighLn(1,3) = 1;                NeighLn(2,3) = 2
  NeighLn(1,6) = 3;                NeighLn(2,6) = 4
  NeighLn(1,9) = 5;                NeighLn(2,9) = 6

  StatusVertex(1) = 1
  StatusVertex(2:6) = 0

  do i = 1, 6
    VertexKey2Value(i) = i
    VertexValue2Key(i)  = i
  enddo
  TailVertex = 7

  !IsDeltaVertex: 1: delta function(Gam); 0: normal function(Gam)
  IsDeltaVertex(1:6) = 1

  ! the site variables of Gamma 1, 2, 3 and 4
  GRVertex(1, 1:6) = 0;            GRVertex(2, 1:6) = 0
  WRVertex(1, 1:6) = 0;            WRVertex(2, 1:6) = 0

  ! the time variables of Gamma 1, 2, 3 and 4
  TVertex(1, 1)=0.d0;  TVertex(2, 1)=0.d0;  TVertex(3, 1)=0.d0
  TVertex(1, 2)=0.1d0; TVertex(2, 2)=0.1d0; TVertex(3, 2)=0.1d0
  TVertex(1, 3)=0.2d0; TVertex(2, 3)=0.2d0; TVertex(3, 3)=0.2d0
  TVertex(1, 4)=0.3d0; TVertex(2, 4)=0.3d0; TVertex(3, 4)=0.3d0
  TVertex(1, 5)=0.4d0; TVertex(2, 5)=0.4d0; TVertex(3, 5)=0.4d0
  TVertex(1, 6)=0.5d0; TVertex(2, 6)=0.5d0; TVertex(3, 6)=0.5d0

  ! Direction of Gamma: 1: left of W;  2: right of W
  DirecVertex(1) = 1;                DirecVertex(2) = 2
  DirecVertex(3) = 1;                DirecVertex(4) = 2
  DirecVertex(5) = 1;                DirecVertex(6) = 2

  ! NeighVertex(i, j): the ith neighbor line of gamma j
  NeighVertex(1,1) = 1;        NeighVertex(2,1) = 2;        NeighVertex(3,1) = 3
  NeighVertex(1,2) = 5;        NeighVertex(2,2) = 4;        NeighVertex(3,2) = 3
  NeighVertex(1,3) = 4;        NeighVertex(2,3) = 7;        NeighVertex(3,3) = 6
  NeighVertex(1,4) = 2;        NeighVertex(2,4) = 8;        NeighVertex(3,4) = 6
  NeighVertex(1,5) = 7;        NeighVertex(2,5) = 5;        NeighVertex(3,5) = 9
  NeighVertex(1,6) = 8;        NeighVertex(2,6) = 1;        NeighVertex(3,6) = 9

  ! weights for lines and vertexes
  WeightLn(1) = weight_gline(StatusLn(1),TVertex(2,NeighLn(2,1))-TVertex(1,NeighLn(1,1)),TypeLn(1))
  WeightLn(2) = weight_gline(StatusLn(2),TVertex(2,NeighLn(2,2))-TVertex(1,NeighLn(1,2)),TypeLn(2))
  WeightLn(4) = weight_gline(StatusLn(4),TVertex(2,NeighLn(2,4))-TVertex(1,NeighLn(1,4)),TypeLn(4))
  WeightLn(5) = weight_gline(StatusLn(5),TVertex(2,NeighLn(2,5))-TVertex(1,NeighLn(1,5)),TypeLn(5))
  WeightLn(7) = weight_gline(StatusLn(7),TVertex(2,NeighLn(2,7))-TVertex(1,NeighLn(1,7)),TypeLn(7))
  WeightLn(8) = weight_gline(StatusLn(8),TVertex(2,NeighLn(2,8))-TVertex(1,NeighLn(1,8)),TypeLn(8))

  dr0(:) = 0
  WeightLn(3) = weight_wline(StatusLn(3),IsDeltaLn(3),dr0,TVertex(3,NeighLn(2,3))-TVertex(3,NeighLn(1,3)),TypeLn(3))
  WeightLn(6) = weight_wline(StatusLn(6),IsDeltaLn(6),dr0,TVertex(3,NeighLn(2,6))-TVertex(3,NeighLn(1,6)),TypeLn(6))
  WeightLn(9) = weight_wline(StatusLn(9),IsDeltaLn(9),dr0,TVertex(3,NeighLn(2,9))-TVertex(3,NeighLn(1,9)),TypeLn(9))

  WeightVertex(1) = weight_vertex(StatusVertex(1), 1, dr0, 0.d0, 0.d0, TypeVertex(1))
  WeightVertex(2) = weight_vertex(StatusVertex(2), 1, dr0, 0.d0, 0.d0, TypeVertex(2))
  WeightVertex(3) = weight_vertex(StatusVertex(3), 1, dr0, 0.d0, 0.d0, TypeVertex(3))
  WeightVertex(4) = weight_vertex(StatusVertex(4), 1, dr0, 0.d0, 0.d0, TypeVertex(4))
  WeightVertex(5) = weight_vertex(StatusVertex(5), 1, dr0, 0.d0, 0.d0, TypeVertex(5))
  WeightVertex(6) = weight_vertex(StatusVertex(6), 1, dr0, 0.d0, 0.d0, TypeVertex(6))


  !ratio = (1.d0/Beta)**Order *SignFermiLoop
  ratio = (-1.d0)**Order *SignFermiLoop
  Anew = d_times_cd(ratio, WeightLn(1)*WeightLn(2)*WeightLn(3)*WeightLn(4)* &
    & WeightLn(5)*WeightLn(6)*WeightLn(7) *WeightLn(8) *WeightLn(9) *WeightVertex(1) &
    & *WeightVertex(2)*WeightVertex(3)*WeightVertex(4)*WeightVertex(5)*WeightVertex(6))

  WeightCurrent = abs(Anew)
  Phase = Anew/WeightCurrent
  !!-------------------------------------------------------
end SUBROUTINE second_order_diagram
!====================================================================
!====================================================================
!====================================================================
