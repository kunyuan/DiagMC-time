!========= MODEL OR DIMENSION DEPENDENT ==============================

!!------------- definition of W0 function ----------------
!!! return 0, 1, 2 (0: W0=0; 1: nearest neighbor; 2: next nearest neighbor)
Integer Function is_W0_nonzero(dims, site)
  implicit none
  integer, intent(in) :: dims, site
  integer, allocatable :: cord(:)
  integer :: dx, dy, dz 
  
  allocate(cord(1:dims))
  cord = get_cord_from_site(D, L, site)
  
  if(D==2) then
    dx = cord(1)
    dy = cord(2)
    is_W0_nonzero = 0
    if(dx>=0  .and. dx<L(1) .and. dy>=0 .and. dy<L(2)) then
      if(dx>dL(1))     dx = L(1)-dx
      if(dy>dL(2))     dy = L(2)-dy
      cord = (/dx, dy/)
    else
      call logFile%QuickLog("Weight_W dx, dy bigger than system size!")
      stop
    endif

    if(cord(1)==1 .and. cord(2)==0) is_W0_nonzero = 1
    if(cord(1)==0 .and. cord(2)==1) is_W0_nonzero = 1

    if(Is_J1J2) then
      if(cord(1)==1 .and. cord(2)==1) is_W0_nonzero = 2
    endif
  else if(D==3) then

    dx = cord(1)
    dy = cord(2)
    dz = cord(3)
    is_W0_nonzero = 0

    if(dx>=0  .and. dx<L(1) .and. dy>=0 .and. dy<L(2) .and. dz>=0 .and. dz<L(3)) then

      if(dx>dL(1))     dx = L(1)-dx
      if(dy>dL(2))     dy = L(2)-dy
      if(dz>dL(3))     dz = L(3)-dz
      cord = (/dx, dy, dz/)
    else
      call logFile%QuickLog("Weight_W dx, dy, dz bigger than system size!")
      stop
    endif

    if(cord(1)==1 .and. cord(2)==0 .and. cord(3)==0) is_W0_nonzero = 1
    if(cord(1)==0 .and. cord(2)==1 .and. cord(3)==0) is_W0_nonzero = 1
    if(cord(1)==0 .and. cord(2)==0 .and. cord(3)==1) is_W0_nonzero = 1

    if(Is_J1J2) then
      !========3-d J1-J2 ============
      if(cord(1)==1 .and. cord(2)==1 .and. cord(3)==0) is_W0_nonzero = 2
      if(cord(1)==0 .and. cord(2)==1 .and. cord(3)==1) is_W0_nonzero = 2
      if(cord(1)==1 .and. cord(2)==0 .and. cord(3)==1) is_W0_nonzero = 2
    endif
  endif
  return
END FUNCTION is_W0_nonzero

!!------------- definition of Gam0 function ----------------
Logical Function is_Gam0_nonzero(dims, site)
  implicit none
  integer, intent(in) :: dims, site
  integer, allocatable :: cord(:)
  integer :: dx, dy, dz
  
  allocate(cord(1:dims))
  cord = get_cord_from_site(D, L, site)

  if(D==2) then
    !========2-d ============
    dx = cord(1)
    dy = cord(2)
    is_Gam0_nonzero = .false.

    if(dx>=0  .and. dx<L(1) .and. dy>=0 .and. dy<L(2)) then

      if(dx>dL(1))     dx = L(1)-dx
      if(dy>dL(2))     dy = L(2)-dy
      cord = (/dx, dy/)
    else
      call logFile%QuickLog("Weight_Gam dx, dy bigger than system size!")
      stop
    endif
    if(cord(1)==0 .and. cord(2)==0) then
      is_Gam0_nonzero = .true.
    endif

  else

    !========3-d ============
    dx = cord(1)
    dy = cord(2)
    dz = cord(3)
    is_Gam0_nonzero = .false.

    if(dx>=0  .and. dx<L(1) .and. dy>=0 .and. dy<L(2) .and. dz>=0 .and. dz<L(3)) then

      if(dx>dL(1))     dx = L(1)-dx
      if(dy>dL(2))     dy = L(2)-dy
      if(dz>dL(3))     dz = L(3)-dz
      cord = (/dx, dy, dz/)
    else
      call logFile%QuickLog("Weight_Gam dx, dy, dz bigger than system size!")
      stop
    endif

    if(cord(1)==0 .and. cord(2)==0 .and. cord(3)==0) then
      is_Gam0_nonzero = .true.
    endif
  endif

  return
END FUNCTION is_Gam0_nonzero

SUBROUTINE generate_r(dims,CurrentSite,NewSite,dsite,Weight,Flag)
!Please make sure Weight is already initialized before calling!
!Flag=.true.: generate new R and new dR
!Flag=.false.: generate new R according to input dR
  implicit none
  logical, intent(in) :: Flag
  integer, intent(in) :: dims, CurrentSite
  integer, intent(out) :: NewSite
  integer :: dSite
  integer :: i, newr(dims), currentr(dims), dr(dims)
  double precision :: Weight,rand

  !dR: CurrentR --> NewR;   dR': CurrentR <-- NewR
  !CurrentR==NewR, dR=0, dR'=0
  !CurrentR>NewR, dR=NewR-CurrentR, dR'=CurrentR-NewR+L
  !CurrentR<NewR, dR=NewR-CurrentR+L, dR'=CurrentR-NewR

  currentr = get_cord_from_site(dims, L, CurrentSite)
  if(Flag .eqv. .false.)  dr = get_cord_from_site(dims, L, dsite)

  newr(:) = 0

  do i=1,dims
    if(Flag) then
      rand=rn()
      dr(i)=0.5d0*dexp(rand*logL(i))
      IF(rn()>0.5d0) dr(i)=L(i)-1-dr(i)
    endif

    if(dr(i)/=0) then
      Weight=Weight*SpatialWeight(i,dr(i))
      Weight=Weight/SpatialWeight(i,L(i)-dr(i))
    endif

    newr(i) = currentr(i) + dr(i)
    if(newr(i)>=L(i)) then
      newr(i) = newr(i) - L(i)
    endif
  enddo

  dSite = get_site_from_cord(dims, L, dr)
  NewSite = get_site_from_cord(dims, L, newr)

  return
END SUBROUTINE generate_r



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

  if(stat == 0) then
    weight_gline = weight_G(typ, tau)
  else if(stat == 1) then
    !  Have measuring vertex around
    weight_gline = weight_meas_G(tau)
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
COMPLEX*16 FUNCTION weight_wline(stat, isdelta, dr, tau, typ)
  implicit none
  integer :: stat, isdelta, dx, dy, typ
  integer :: dr
  double precision :: tau
  integer :: t

  if(stat == 0) then
    if(isdelta==0) weight_wline = weight_W(typ, dr, tau)
    if(isdelta==1) weight_wline = weight_W0(typ, dr)
  else if(stat == 2) then
    ! Have Ira or Masha around 
    if(isdelta==0) weight_wline = weight_W(1, dr, tau)
    if(isdelta==1) weight_wline = weight_W0(1, dr)
  else if(stat == 1 .or. stat==3) then
    ! Have measuring vertex around
    if(isdelta==0) weight_wline = weight_meas_W(dr, tau)
    if(isdelta==1) weight_wline = (0.d0, 0.d0)
  else
    call LogFile%WriteStamp('e')
    call LogFile%WriteLine("The number of update: "+str(iupdate))
    call LogFile%WriteLine("line status error!"+str(stat))
    call print_config
    stop
  endif

  !---------------------- test1: fake function ----------------------------

  !---------------------- test2: uniform function ----------------------------
  !if(stat >= 0 .and. stat<=3) then
    !if(isdelta==0 .and. dr==0 ) then
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

COMPLEX*16 FUNCTION weight_vertex(stat, isdelta, dr, dtau1, dtau2, typ)
  implicit none
  integer :: stat, dr, t1, t2, typ, isdelta
  integer :: flag
  double precision :: weight
  double precision :: dtau1, dtau2

  if(stat==0) then
    if(IS_BOLD) then
      !----------------- for bold Gamma ------------------------------
      if(isdelta==0) weight_vertex = weight_Gam(typ, dr, dtau1, dtau2)
      if(isdelta==1) weight_vertex = weight_Gam0(typ, dr)
    else
      !----------------- for bare Gamma ------------------------------
      if(isdelta==0) weight_vertex = (0.d0, 0.d0)
      if(isdelta==1) weight_vertex = weight_Gam0(typ, dr)
    endif

  else if(stat==2) then
    if(IS_BOLD) then
      !----------------- for bold Gamma ------------------------------
      if(isdelta==0) weight_vertex = weight_Gam(typ, dr, dtau1, dtau2)
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
      !if(dr==0 ) then
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



!!====================================================================
!=====================================================================
!==================FAKE WEIGHT FOR MEASURING AND WORM=================
!=====================================================================


!--------- worm weight function  ---------
DOUBLE PRECISION FUNCTION weight_worm(drg, drw, dtau)
  implicit none
  integer, intent(in)  :: drg, drw
  double precision, intent(in) :: dtau

  weight_worm = 1.d0
  !weight_worm = weight_worm*dexp(-dtau)
  !weight_worm = weight_worm*dexp(-0.5d0*(abs(dxg)+abs(dyg)))
  !weight_worm = weight_worm*dexp(-0.5d0*(abs(dxw)+abs(dyw)))

  return
END FUNCTION weight_worm


!--------- measure weight function for Gamma ---------
complex*16 FUNCTION weight_meas_G(tau1)
  implicit none
  double precision, intent(in)  :: tau1
  weight_meas_G = (1.d0, 0.d0)
  return
END FUNCTION weight_meas_G

complex*16 FUNCTION weight_meas_W(dr, tau1)
  implicit none
  integer, intent(in)  :: dr
  double precision, intent(in) :: tau1

  weight_meas_W = (1.d0, 0.d0)
  return
END FUNCTION weight_meas_W


complex*16 FUNCTION weight_meas_Gam0(ityp, dr)
  implicit none
  integer, intent(in)  :: dr, ityp

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
  if(dr==0) then
    if(ityp==1 .or. ityp==2 .or. ityp==5 .or. ityp==6) then
      weight_meas_Gam0 = (1.d0, 0.d0)
    endif
  endif
  return
END FUNCTION weight_meas_Gam0


!!--------- extract weight for G ---------
COMPLEX*16 FUNCTION weight_G(typ1, tau1)
  implicit none
  double precision, intent(in)  :: tau1
  integer, intent(in) :: typ1
  call interpolation_GW(tau1, G(typ1, :), weight_G)
  if(tau1<=1.d-14)  weight_G = -1.d0* weight_G
END FUNCTION weight_G


!!--------- extract weight for W ---------
COMPLEX*16 FUNCTION weight_W(typ1, site, tau1)
  implicit none
  integer, intent(in)  :: site, typ1
  double precision, intent(in) :: tau1
  call interpolation_GW(tau1, W(typ1, site, :), weight_W)
END FUNCTION weight_W


!!--------- extract weight for Gamma ---------
COMPLEX*16 FUNCTION weight_Gam(typ1, site, tau1, tau2)
  implicit none
  integer, intent(in)  :: site, typ1
  double precision, intent(in) :: tau1, tau2
  integer :: ityp, isite
  double precision :: GaR
  ityp = (typ1+1)/2
  isite = diff_r(D, site, 0)
  if(if_inside_Vol(D, dGamL(1:D), isite)) then
    call interpolation_Gam(tau1, tau2, GamMCInput(ityp,isite,:,:), weight_Gam)
  else
    weight_Gam = (0.d0, 0.d0)
  endif

  if(tau1<=1.d-14 .and. tau2>1.d-14) then
    weight_Gam = -1.d0 *weight_Gam
  else if(tau1>1.d-14 .and. tau2<=1.d-14) then
    weight_Gam = -1.d0 *weight_Gam
  endif
END FUNCTION weight_Gam



SUBROUTINE interpolation_GW(tau, Table, Val)
  implicit none
  double precision, intent(in) :: tau
  complex*16, dimension(0:MxT-1), intent(in) :: Table
  complex*16, intent(out) :: Val
  integer :: ts, tb
  double precision :: ptau, taus, taub

  ptau = shift_tau(tau)

  ts = Floor((ptau-0.5d0*Beta/MxT)*MxT/Beta)
  tb = ts + 1

  if(ts<0) then
    ts=0
    tb=1
  else if(tb>=MxT) then
    ts=MxT-2
    tb=MxT-1
  endif

  taus = (real(ts)+0.5d0)*Beta/MxT
  taub = (real(tb)+0.5d0)*Beta/MxT

  Val = interpolate(ptau, taus, Table(ts), taub, Table(tb))
  return
END SUBROUTINE interpolation_GW


SUBROUTINE interpolation_Gam(tau1, tau2, Table, Val)
  implicit none
  double precision, intent(in) :: tau1, tau2
  complex*16, dimension(0:MxT-1, 0:MxT-1), intent(in) :: Table
  complex*16, intent(out) :: Val
  integer :: i, ts(1:2), tb(1:2)
  double precision :: ptau(1:2), taus(1:2), taub(1:2)
  complex*16 :: Val1, Val2

  ptau(1) = shift_tau(tau1)
  ptau(2) = shift_tau(tau2)

  do i = 1, 2
    ts(i) = Floor((ptau(i)-0.5d0*Beta/MxT)*MxT/Beta)
    if(ts(i)<0) then
      ts(i) = 0
    else if(ts(i)>=MxT-1) then
      ts(i) = MxT-2
    endif

    tb(i) = ts(i) + 1

    taus(i) = (real(ts(i))+0.5d0)*Beta/MxT
    taub(i) = (real(tb(i))+0.5d0)*Beta/MxT
  enddo

  Val1 = interpolate(ptau(1), taus(1), Table(ts(1), ts(2)), taub(1), Table(tb(1), ts(2)))
  Val2 = interpolate(ptau(1), taus(1), Table(ts(1), tb(2)), taub(1), Table(tb(1), tb(2)))

  Val  = interpolate(ptau(2), taus(2), Val1, taub(2), Val2)
  return
END SUBROUTINE interpolation_Gam

complex*16 FUNCTION interpolate(tau, taus, Fs, taub, Fb)
  double precision, intent(in) :: tau
  double precision, intent(in) :: taus, taub
  complex*16, intent(in) :: Fs, Fb
  interpolate = ((tau-taus)*Fb+(taub-tau)*Fs)*MxT/Beta
  return
END FUNCTION interpolate

DOUBLE PRECISION FUNCTION shift_tau(tau)
  implicit none
  double precision, intent(in) :: tau
  shift_tau = tau
  if(tau<=1.d-14) then
    if(abs(tau)<=1.d-14) then
      shift_tau = Beta
    else
      shift_tau = tau + Beta
    endif
  endif
  return
END FUNCTION shift_tau
