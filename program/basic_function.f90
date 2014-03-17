!============== BASIC PROCEDURES IN UPDATE ==========================


!=====================================================================
!==================FAKE WEIGHT FOR MEASURING AND WORM=================
!=====================================================================


SUBROUTINE calculate_GamNormWeight
  implicit none
  integer :: omega1, omegaW, omega2, ityp
  double precision :: Gam0

  GamNormWeight = 1.d0
  !do omega1 = -MxOmegaDiag, MxOmegaDiag
    !do omega2 = -MxOmegaDiag, MxOmegaDiag
      !omegaW = omega1-omega2
      !do ityp = 1, 6
        !if(ityp==3 .or. ityp==4) cycle
        !Gam0 = weight_meas_W(0,0,omegaW,ityp)
        !Gam0 = Gam0 *weight_meas_Gamma(0,0,0,omega1,omega2,1)
        !Gam0 = Gam0 *weight_Gamma0(0,0,omega1,omega2,1)
        !GamNormWeight = GamNormWeight + Gam0
      !enddo
    !enddo
  !enddo
  return
END SUBROUTINE calculate_GamNormWeight

!--------- worm weight function  ---------
DOUBLE PRECISION FUNCTION weight_worm(dxg, dyg, dxw, dyw, dtau)
  implicit none
  integer, intent(in)  :: dxg, dyg, dxw, dyw
  double precision, intent(in) :: dtau

  weight_worm = 1.d0
  !weight_worm = 1.d0/(1.d0+abs(domega)**2.d0)
  !weight_worm = weight_worm*exp(-0.5d0*(abs(dxg)+abs(dyg)))
  !weight_worm = weight_worm*exp(-0.5d0*(abs(dxw)+abs(dyw)))

  return
END FUNCTION weight_worm


!--------- worm weight function for W ---------
complex*16 FUNCTION weight_worm_W(ityp, dx, dy, t1)
  implicit none
  integer, intent(in)  :: dx, dy, t1, ityp

  !weight_worm_W = weight_W(1, dx, dy, t1)
  weight_worm_W = weight_meas_W(1, dx, dy, t1)
  return
END FUNCTION weight_worm_W

!--------- worm weight function for Gamma ---------
complex*16 FUNCTION weight_worm_Gam(ityp, dx, dy, t1, t2)
  implicit none
  integer, intent(in)  :: dx, dy, t1, t2, ityp

  !weight_worm_Gam = weight_Gam0(1, dx, dy)
  weight_worm_Gam = weight_meas_Gam(1, dx, dy, t1, t2)
  return
END FUNCTION weight_worm_Gam

!--------- measure weight function for Gamma ---------
complex*16 FUNCTION weight_meas_G(ityp, t1)
  implicit none
  integer, intent(in)  :: t1, ityp
  weight_meas_G = 1.d0
  return
END FUNCTION weight_meas_G

complex*16 FUNCTION weight_meas_W(ityp, dx, dy, t1)
  implicit none
  integer, intent(in)  :: dx, dy, t1, ityp
  weight_meas_W = 1.d0
  return
END FUNCTION weight_meas_W

complex*16 FUNCTION weight_meas_Gam(ityp, dx, dy, t1, t2)
  implicit none
  integer, intent(in)  :: dx, dy, t1, t2, ityp

  weight_meas_Gam = 0.d0
  if(dx==0 .and. dy==0) then
    if(ityp ==1 .or. ityp == 2) then
      weight_meas_Gam = 1.d0
    else if(ityp == 3 .or. ityp == 4) then
      weight_meas_Gam = 0.d0
    else if(ityp == 5 .or. ityp == 6) then
      weight_meas_Gam = 1.d0
    endif
  endif

  !================   test ================
  !=======  r independent ================
  !weight_meas_Gam=0.d0
  !if(ityp ==1 .or. ityp == 2) then
    !weight_meas_Gam = 1.d0
  !else if(ityp == 3 .or. ityp == 4) then
    !weight_meas_Gam = 0.d0
  !else if(ityp == 5 .or. ityp == 6) then
    !weight_meas_Gam = 1.d0
  !endif
  !======================================

  return
END FUNCTION weight_meas_Gam



!!------------- definition of the system symmetry ----------------
SUBROUTINE def_symmetry
  implicit none
  integer :: i, j, omega

  CoefOfSymmetry(:,:) = 2.d0

  !do i = 1, Lx-1
    !CoefOfSymmetry(i, :) = 2.d0* CoefOfSymmetry(i, :)
  !enddo

  !do j = 1, Ly-1
    !CoefOfSymmetry(:, j) = 2.d0* CoefOfSymmetry(:, j)
  !enddo
  
  return
END SUBROUTINE def_symmetry

!SUBROUTINE update_WeightCurrent
  !implicit none
  !integer :: i, Gam1, Gam2, G, W
  !double precision :: weight
  !double precision :: wln(MxNLn), wgam(MxNVertex)

  !weight = 1.d0
  !do i = 1, MxNLn
    !if(StatusLn(i)<0)  cycle
    !if(KindLn(i)==1) then
      !wln(i) = weight_line(StatusLn(i),1,0,0,OmegaLn(i),TypeLn(i))
    !else
      !Gam1 = NeighLn(1,i);       Gam2 = NeighLn(2,i)
      !wln(i) = weight_line(StatusLn(i),2,WXVertex(Gam1)-WXVertex(Gam2), &
        !& WYVertex(Gam1)-WYVertex(Gam2), OmegaLn(i), TypeLn(i))
    !endif
    !WeightLn(i) = wln(i)
    !weight = weight *wln(i)
  !enddo

  !do i = 1, MxNVertex
    !if(StatusVertex(i)<0)  cycle
    !G = NeighVertex(1, i);           W = NeighVertex(3, i)
    !wgam(i) = weight_vertex(Order, StatusVertex(i),GXVertex(i)-WXVertex(i), &
      !& GYVertex(i)-WYVertex(i), OmegaLn(G), OmegaLn(W), TypeVertex(i))
    !WeightVertex(i) = wgam(i)
    !weight = weight *wgam(i)
  !enddo

  !weight = weight*CoefOfWeight(Order)*(1.d0/Beta)**Order *(-1.d0)**NFermiLoop

  !WeightCurrent = weight

  !return
!END SUBROUTINE update_WeightCurrent




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
SUBROUTINE generate_x(CurrentX,NewX,Weight)
  implicit none
  integer :: NewX,CurrentX,dX
  double precision :: Weight,rand
  rand=rn()
  dX=0.5d0*dexp(rand*logLx)
  Weight=SpatialWeight(1,dX)
  IF(rn()>0.5d0) dX=Lx-1-dX

  NewX = CurrentX + dX
  if(NewX>=Lx) then
    NewX = NewX - Lx
  endif
  return
END SUBROUTINE generate_x

SUBROUTINE generate_y(CurrentY,NewY,Weight)
  implicit none
  integer :: NewY,CurrentY,dY
  double precision :: Weight,rand
  rand=rn()
  dY=0.5d0*dexp(rand*logLy)
  Weight=SpatialWeight(2,dY)
  IF(rn()>0.5d0) dY=Ly-1-dY

  NewY = CurrentY + dY
  if(NewY>=Ly) then
    NewY = NewY - Ly
  endif
  return
END SUBROUTINE generate_y

DOUBLE PRECISION FUNCTION prob_dx(x)
  implicit none 
  integer :: dx,x
  if(dx<0)then 
    dx=x+Lx
  else 
    dx=x
  endif
  prob_dx = SpatialWeight(1,dx)
  return
END FUNCTION prob_dx

DOUBLE PRECISION FUNCTION prob_dy(y)
  implicit none 
  integer :: dy,y
  if(dy<0)then
    dy=y+Ly;
  else
    dy=y
  endif
  prob_dy = SpatialWeight(2,dy)
  return
END FUNCTION prob_dy

INTEGER FUNCTION diff_x(dx)
  implicit none
  integer, intent(in) :: dx
  diff_x = dx
  if(diff_x<0)     diff_x = Lx+diff_x
END FUNCTION diff_x

INTEGER FUNCTION diff_y(dy)
  implicit none
  integer, intent(in) :: dy
  diff_y = dy
  if(diff_y<0)     diff_y = Ly+diff_y
END FUNCTION diff_y

LOGICAL FUNCTION Is_x_valid(x1, x2)
  implicit none
  integer, intent(in) :: x1, x2
  integer :: dx
  dx = x1 - x2
  if(dx<0)     dx = Lx+dx
  if(dx>dLx)   dx = Lx-dx
  if(abs(dx)<=1) then
    Is_x_valid = .true.
  else
    Is_x_valid = .false.
  endif
END FUNCTION Is_x_valid

LOGICAL FUNCTION Is_y_valid(y1, y2)
  implicit none
  integer, intent(in) :: y1, y2
  integer :: dy
  dy = y1 - y2
  if(dy<0)     dy = Ly+dy
  if(dy>dLy)   dy = Ly-dy
  if(abs(dy)<=1) then
    Is_y_valid = .true.
  else
    Is_y_valid = .false.
  endif
END FUNCTION Is_y_valid


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
INTEGER FUNCTION generate_gamma()
  implicit none
  integer :: rand
  rand = Floor(rn()*NVertex)+1
  generate_gamma = VertexKey2Value(rand)
  return
END FUNCTION generate_gamma
!!!=======================================================================
!!!=======================================================================
!!!=======================================================================



!!=======================================================================
!!================= IRREDUCIBILITY CHECK ================================
!!=======================================================================

!-------- update G Hash table -----------------------------------
SUBROUTINE update_Hash4G(oldk, newk)
  implicit none
  integer, intent(in) :: oldk, newk

  if(Hash4G(oldk)==1) then
    Hash4G(oldk)=0
  else
    if(CheckG) then
      write(*, *) "===================================="
      write(*, *) "Oops, update_Hash4G found a bug!"
      write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      write(*, *) "G Hash table for old k", oldk, " is not 1!!", Hash4G(oldk)
      write(*, *) "===================================="
      call print_config
      stop
    endif
  endif

  if(Hash4G(newk)==0) then
    Hash4G(newk)=1
  else
    if(CheckG) then
      write(*, *) "===================================="
      write(*, *) "Oops, update_Hash4G found a bug!"
      write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      write(*, *) "G Hash table for new k", newk, " is not 0!!", Hash4G(newk)
      write(*, *) "===================================="
      call print_config
      stop
    endif
  endif

  return
END SUBROUTINE update_Hash4G

SUBROUTINE add_Hash4G(newk)
  implicit none
  integer, intent(in) :: newk
  if(Hash4G(newk)==0) then
    Hash4G(newk)=1
  else
    if(CheckG) then
      write(*, *) "===================================="
      write(*, *) "Oops, add_Hash4G found a bug!"
      write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      write(*, *) "G Hash table for new k", newk, " is not 0!!", Hash4G(newk)
      write(*, *) "===================================="
      call print_config
      stop
    endif
  endif
  return
END SUBROUTINE add_Hash4G

SUBROUTINE delete_Hash4G(oldk)
  implicit none
  integer, intent(in) :: oldk
  if(Hash4G(oldk)==1) then
    Hash4G(oldk)=0
  else
    if(CheckG) then
      write(*, *) "===================================="
      write(*, *) "Oops, delete_Hash4G found a bug!"
      write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      write(*, *) "G Hash table for old k", oldk, " is not 1!!", Hash4G(oldk)
      write(*, *) "===================================="
      call print_config
      stop
    endif
  endif
  return
END SUBROUTINE delete_Hash4G


!-------- update W Hash table -----------------------------------
SUBROUTINE update_Hash4W(oldk, newk)
  implicit none
  integer, intent(in) :: oldk, newk
  integer :: aoldk, anewk
  aoldk = abs(oldk)
  anewk = abs(newk)

  if(Hash4W(aoldk)==1) then
    Hash4W(aoldk)=0
  else
    if(CheckW) then
      write(*, *) "===================================="
      write(*, *) "Oops, update_Hash4W found a bug!"
      write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      write(*, *) "W Hash table for old k", aoldk, " is not 1!!", Hash4W(aoldk)
      write(*, *) "===================================="
      call print_config
      stop
    endif
  endif

  if(Hash4W(anewk)==0) then
    Hash4W(anewk)=1
  else
    if(CheckW) then
      write(*, *) "===================================="
      write(*, *) "Oops, update_Hash4W found a bug!"
      write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      write(*, *) "W Hash table for new k", anewk, " is not 0!!", Hash4W(anewk)
      write(*, *) "===================================="
      call print_config
      stop
    endif
  endif
  return
END SUBROUTINE update_Hash4W

SUBROUTINE add_Hash4W(newk)
  implicit none
  integer, intent(in) :: newk
  integer :: anewk
  anewk = abs(newk)

  if(Hash4W(anewk)==0) then
    Hash4W(anewk)=1
  else
    if(CheckW) then
      write(*, *) "===================================="
      write(*, *) "Oops, add_Hash4W found a bug!"
      write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      write(*, *) "W Hash table for new k", anewk, " is not 0!!", Hash4W(anewk)
      write(*, *) "===================================="
      call print_config
      stop
    endif
  endif
  return
END SUBROUTINE add_Hash4W

SUBROUTINE delete_Hash4W(oldk)
  implicit none
  integer, intent(in) :: oldk
  integer :: aoldk
  aoldk = abs(oldk)

  if(Hash4W(aoldk)==1) then
    Hash4W(aoldk)=0
  else
    if(CheckW) then
      write(*, *) "===================================="
      write(*, *) "Oops, delete_Hash4W found a bug!"
      write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      write(*, *) "W Hash table for old k", aoldk, " is not 1!!", Hash4W(aoldk)
      write(*, *) "===================================="
      call print_config
      stop
    endif
  endif
  return
END SUBROUTINE delete_Hash4W



!--------- check the irreducibility for G -----------------------
LOGICAL FUNCTION Is_reducible_G(GLn)
  implicit none
  integer, intent(in) :: GLn
  integer :: newk

  Is_reducible_G = .false.
  newk = kLn(GLn)

  if(Is_k_valid(newk) .eqv. .false.) then
    Is_reducible_G = .true.
    return
  endif

  if(CheckG) then
    if(Hash4G(newk)==1) then
      Is_reducible_G = .true.
      return
    else if(Hash4G(newk)/=0) then
      write(*, *) "===================================="
      write(*, *) "Oops, Is_reducible_G found a bug!"
      write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      write(*, *) "G Hash table for new k", newk, " is not 0 or 1!!", Hash4G(newk)
      write(*, *) "===================================="
      call print_config
      stop
    endif
  endif
    
  return
END FUNCTION Is_reducible_G


LOGICAL FUNCTION Is_reducible_G_Gam(GLn)
  implicit none
  integer, intent(in) :: GLn
  integer :: nG, Gam1, Gam2, W1, W2, nW
  integer :: newk, kG 
  integer :: i, nnk1, nnk2

  Is_reducible_G_Gam = .false.
  newk = kLn(GLn)

  Gam1 = NeighLn(1, GLn)
  Gam2 = NeighLn(2, GLn)
  W1 = NeighVertex(3, Gam1)
  W2 = NeighVertex(3, Gam2)

  if(CheckGam) then
    do i = 1, NGLn
      nG = GLnKey2Value(i)
      kG = kLn(nG)
      if(nG==GLn)   cycle             !! rule out the line itself
      if(nG==NeighVertex(1, Gam1) .or. nG==NeighVertex(2, Gam2)) cycle
      if(abs(add_k(newk,-kG))==abs(kLn(W1)) .or. &
        & abs(add_k(newk,-kG))==abs(kLn(W2))) cycle !! rule out the neighbor wlines of newk
      nnk1 = kLn(NeighVertex(3, NeighLn(1, nG)))
      nnk2 = kLn(NeighVertex(3, NeighLn(2, nG)))
      if(abs(add_k(newk, -kG))==abs(nnk1) .or. &
        & abs(add_k(newk, -kG))==abs(nnk2)) cycle !! rule out the neighbor wlines of neighbor of newk

      if(Hash4W(abs(add_k(newk,-kG)))==1) then
        Is_reducible_G_Gam = .true.
      endif
    enddo
  endif
    
  return
END FUNCTION Is_reducible_G_Gam




!--------- check the irreducibility for W -----------------------
LOGICAL FUNCTION Is_reducible_W(WLn)
  implicit none
  integer, intent(in) :: WLn
  integer :: absk 

  absk = abs(kLn(WLn))
  Is_reducible_W = .false.

  if(Is_k_valid(absk) .eqv. .false. ) then
    Is_reducible_W = .true.
    return
  endif

  if(CheckW) then
    if(absk==0) then
      Is_reducible_W = .true.
      return
    endif

    if(Hash4W(absk)==1) then
      Is_reducible_W = .true.
      return
    else if(Hash4W(absk)/=0) then
      write(*, *) "===================================="
      write(*, *) "Oops, Is_reducible_W found a bug!"
      write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      write(*, *) "W Hash table for new k", absk, " is not 0 or 1!!", Hash4W(absk)
      write(*, *) "===================================="
      call print_config
      stop
    endif
  endif

  return
END FUNCTION Is_reducible_W



!--------- check the irreducibility for W -----------------------
LOGICAL FUNCTION Is_reducible_W_Gam(WLn)
  implicit none
  integer, intent(in) :: WLn
  integer :: absk, i, Gam1, Gam2, G1, G2, G3, G4, kG5, kG6
  integer :: nG, kG, pkGG, nkGG

  absk = abs(kLn(WLn))
  Is_reducible_W_Gam = .false.
  Gam1 = NeighLn(1, WLn)
  Gam2 = NeighLn(2, WLn)
  G1 = NeighVertex(1, Gam1)
  G2 = NeighVertex(2, Gam1)
  G3 = NeighVertex(1, Gam2)
  G4 = NeighVertex(2, Gam2)

  if(CheckGam) then
    do i = 1, NGLn
      nG = GLnKey2Value(i)
      kG5 = kLn(NeighVertex(1, NeighLn(1,nG)))
      kG6 = kLn(NeighVertex(2, NeighLn(2,nG)))
      kG = kLn(nG)
      if(nG==G1 .or. nG==G2 .or. nG==G3 .or. nG==G4) cycle
      pkGG = add_k(kG, absk)
      nkGG = add_k(kG, -absk)
      if(pkGG==kLn(G1) .or. pkGG==kLn(G2) .or. pkGG==kLn(G3) .or. pkGG==kLn(G4)) cycle
      if(nkGG==kLn(G1) .or. nkGG==kLn(G2) .or. nkGG==kLn(G3) .or. nkGG==kLn(G4)) cycle
      if(pkGG==kG5 .or. pkGG==kG6 .or. nkGG==kG5 .or. nkGG==kG6) cycle

      if(Hash4G(pkGG)==1 .or. Hash4G(nkGG)==1) then
        Is_reducible_W_Gam = .true.
      endif
    enddo
  endif
  return
END FUNCTION Is_reducible_W_Gam

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
    write(*, *)  "add_ira_stat, stat = -1"
    call print_config
    stop
  endif
  if(stat<=1) then
    add_ira_stat = stat + 2
  else
    add_ira_stat = stat
  endif
  if(add_ira_stat>=4) then
    write(*, *)  IsWormPresent, iupdate, "add_ira_stat, stat > 3", add_ira_stat
    call print_config
    stop
  endif
END FUNCTION add_ira_stat

!-- change the status from i/m or measuring+i/m to normal or measuring--
INTEGER FUNCTION delete_ira_stat(stat)
  implicit none
  integer, intent(in) :: stat
  if(stat == -1) then
    write(*, *)  IsWormPresent, iupdate, "del_ira_stat, stat = -1"
    call print_config
    stop
  endif
  delete_ira_stat = stat - 2
  !if(delete_ira_stat == -2) delete_ira_stat = 0
  !if(delete_ira_stat == -1) delete_ira_stat = 1
  if(delete_ira_stat<0) then
    write(*, *)  IsWormPresent, iupdate, "del_ira_stat, stat < 0", delete_ira_stat
    call print_config
    stop
  endif
END FUNCTION delete_ira_stat

!-- change the status from normal or i/m to measuring or measuring+i/m --
INTEGER FUNCTION add_mea_stat(stat)
  implicit none
  integer, intent(in) :: stat
  if(stat == -1) then
    write(*, *)  "add_mea_stat, stat = -1"
    call print_config
    stop
  endif

  add_mea_stat = stat + 1
  if(add_mea_stat == 2)   add_mea_stat = 1

  if(add_mea_stat>=4) then
    write(*, *)  IsWormPresent, iupdate, "add_mea_stat, stat > 3", add_mea_stat
    call print_config
    stop
  endif
END FUNCTION add_mea_stat

!-- change the status from measuring or measuring+i/m to normal or i/m--
INTEGER FUNCTION delete_mea_stat(stat)
  implicit none
  integer,intent(in) :: stat
  if(stat == -1) then
    write(*, *)  "del_mea_stat, stat = -1"
    call print_config
    stop
  endif
  delete_mea_stat = stat - 1
  if(delete_mea_stat<0) then
    write(*, *)  IsWormPresent, iupdate, "del_mea_stat, stat < 0", delete_mea_stat
    call print_config
    stop
  endif
END FUNCTION delete_mea_stat

!--- calculate the status of a line according to the neighbor vertexes --
INTEGER FUNCTION line_stat(stat1, stat2)
  implicit none
  integer :: stat1, stat2
  if(stat1 == -1) stop
  if(stat2 == -1) stop
  line_stat = stat1 + stat2
  if(stat1 == stat2) then
    line_stat = stat1
  else if(line_stat == 5) then
    line_stat = 3
  endif

  if(line_stat>3) then
    write(*, *) "line_stat error!", line_stat, stat1, stat2
    call print_config
    stop
  endif
  return
END FUNCTION line_stat


!------------- insert a line to the link -------------------
SUBROUTINE insert_line(newline, isdelta, k, knd, typ, stat, weigh)
  implicit none
  integer, intent(out) :: newline
  integer, intent(in) :: isdelta, k, knd, typ, stat
  double precision, intent(in) :: weigh

  newline = TailLn
  TailLn  = NextLn(TailLn)
  if(StatusLn(TailLn)>=0) then
    write(*, *) IsWormPresent, iupdate, "insert_line error!!!"
    call print_config
    stop
  endif

  if(TailLn == -1) then
    write(*, *) "Tail=-1! Too many lines!"
    call print_config
    stop
  endif

  if(knd==1) then
    NGLn = NGLn + 1
    GLnKey2Value(NGLn) = newline
    LnValue2Key(newline) = NGLn
  else
    NWLn = NWLn + 1
    WLnKey2Value(NWLn) = newline
    LnValue2Key(newline) = NWLn
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
    write(*, *) IsWormPresent, iupdate, "delete_line error!!!"
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
    write(*, *) "Tail=-1! Too many lines!"
    stop
  endif
  return
END SUBROUTINE undo_insert_line


!------------- insert a gamma to the link -------------------
SUBROUTINE insert_gamma(newgamma, isdelta, gx, gy, wx, wy, t1, t2, t3, dir, typ, stat, weigh)
  implicit none
  integer, intent(out) :: newgamma
  integer, intent(in) :: gx, gy, wx, wy, isdelta, dir, typ, stat
  double precision, intent(in) :: weigh, t1, t2, t3

  newgamma = TailVertex
  TailVertex = NextVertex(TailVertex)
  if(StatusVertex(TailVertex)>=0) then
    write(*, *) IsWormPresent, iupdate, "insert_gamma error!!!"
    call print_config
    stop
  endif

  if(TailVertex == -1) then
    write(*, *) "Tail=-1! Too many gammas!"
    stop
  endif
   
  NVertex = NVertex + 1
  VertexKey2Value(NVertex) = newgamma
  VertexValue2Key(newgamma) = NVertex

  IsDeltaVertex(newgamma) = isdelta
  GXVertex(newgamma) = gx
  GYVertex(newgamma) = gy
  WXVertex(newgamma) = wx
  WYVertex(newgamma) = wy
  T1Vertex(newgamma) = t1
  T2Vertex(newgamma) = t2
  T3Vertex(newgamma) = t3
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
    !call delete_Hash4G(kLn(occline))
  else
    call delete_Hash4W(kLn(occline))
  endif

  if(StatusLn(occline)==-1) then
    write(*, *) IsWormPresent, iupdate, "delete_line error!!!"
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
    write(*, *) "Tail=-1! Too many lines!"
    stop
  endif
  return
END SUBROUTINE delete_line

!------------- insert a gamma to the link -------------------
SUBROUTINE undo_delete_line(newline, knd, stat)
  implicit none
  integer, intent(in) :: newline, knd, stat

  if(TailLn/=newline)    then
    write(*, *) "undo_delete_line error!"
    call print_config
    stop
  endif
  StatusLn(newline) = stat
  TailLn = NextLn(newline)
  if(TailLn == -1) then
    write(*, *) "Tail=-1! Too many glines!"
    call print_config
    stop
  endif

  if(knd==1) then
    NGLn = NGLn + 1
    GLnKey2Value(NGLn) = newline
    LnValue2Key(newline) = NGLn
    !call add_Hash4G(kLn(newline))
  else
    NWLn = NWLn + 1
    WLnKey2Value(NWLn) = newline
    LnValue2Key(newline) = NWLn
    call add_Hash4W(kLn(newline))
  endif
  return
END SUBROUTINE undo_delete_line


!------------- delete a gamma from the link -------------------
SUBROUTINE delete_gamma(occgamma)
  implicit none
  integer, intent(in) :: occgamma
  integer :: tmp

  if(StatusVertex(occgamma)==-1) then
    write(*, *) IsWormPresent, iupdate, occgamma, StatusVertex(occgamma), "delete_gamma error!!!"
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
    write(*, *) "Tail=-1! Too many vertexes!"
    stop
  endif
  return
END SUBROUTINE delete_gamma


SUBROUTINE undo_delete_gamma(newgamma)
  implicit none
  integer, intent(in) :: newgamma

  if(TailVertex/=newgamma)    then
    write(*, *) "undo_delete_gamma error!"
    stop
  endif
  StatusVertex(newgamma) = 0

  TailVertex = NextVertex(newgamma)
  if(TailVertex == -1) then
    write(*, *) "Tail=-1! Too many gammas!"
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





!!============== GRAM-SCHMIDT BASIS ==================================

!!-------------- generate the Gram-Schmidt basis --------------
!SUBROUTINE calculate_basis(nbasis, OmegaMin, OmegaMax, EP, EN)
  !implicit none
  !integer, intent(in):: nbasis, OmegaMin, OmegaMax
  !double precision :: EP(nbasis, nbasis), EN(nbasis, nbasis)
  !double precision :: VP(nbasis, nbasis), VN(nbasis, nbasis)
  !double precision :: norm
  !integer :: i, j, k

  !EP(:,:) = 0.d0
  !EN(:,:) = 0.d0
  !VP(:,:) = 0.d0
  !VN(:,:) = 0.d0

  !do i = 1, nbasis
    !VP(i, i) = 1.d0
    !EP(i, :) = VP(i, :)

    !do j = 1, i-1
      !do k = 1, j
        !EP(i, k) = EP(i, k) - EP(j, k) *projector(EP(j,:), VP(i, :), OmegaMin, OmegaMax)
      !enddo
    !enddo

    !norm = dsqrt(projector(EP(i,:), EP(i,:), OmegaMin, OmegaMax))
    !EP(i, :) = EP(i, :)/norm
  !enddo

  !do i = 1, nbasis
    !VN(i, i) = 1.d0
    !EN(i, :) = VN(i, :)

    !do k = 1, j
      !do j = 1, i-1
        !EN(i, k) = EN(i, k) - EN(j, k) *projector(EN(j,:), VN(i, :), -OmegaMax, -OmegaMin)
      !enddo
    !enddo
    !norm = dsqrt(projector(EN(i,:), EN(i,:), -OmegaMax, -OmegaMin))
    !EN(i, :) = EN(i, :)/norm
  !enddo

!END SUBROUTINE calculate_basis



!SUBROUTINE calculate_Gamma_basis(nbasisGamma, OmegaMin, OmegaMax, EPN, ENP, &
  !&  EPPR, ENNR, EPPL, ENNL)
  !implicit none
  !integer, intent(in):: nbasisGamma, OmegaMin, OmegaMax
  !double precision :: ENP(nbasisGamma, nbasisGamma),  EPN(nbasisGamma, nbasisGamma)
  !double precision :: EPPR(nbasisGamma, nbasisGamma), ENNR(nbasisGamma, nbasisGamma)
  !double precision :: EPPL(nbasisGamma, nbasisGamma), ENNL(nbasisGamma, nbasisGamma)
  !double precision :: VNP(nbasisGamma, nbasisGamma),  VPN(nbasisGamma, nbasisGamma)
  !double precision :: VPPR(nbasisGamma, nbasisGamma), VNNR(nbasisGamma, nbasisGamma)
  !double precision :: VPPL(nbasisGamma, nbasisGamma), VNNL(nbasisGamma, nbasisGamma)
  !double precision :: norm
  !integer :: i, j, k
  !VPN(:,:) = 0.d0
  !VNP(:,:) = 0.d0
  !VPPR(:,:) = 0.d0
  !VNNR(:,:) = 0.d0
  !VPPL(:,:) = 0.d0
  !VNNL(:,:) = 0.d0

  !EPN(:,:) = 0.d0
  !ENP(:,:) = 0.d0
  !EPPR(:,:) = 0.d0
  !ENNR(:,:) = 0.d0
  !EPPL(:,:) = 0.d0
  !ENNL(:,:) = 0.d0

  !do i = 1, nbasisGamma
    !VPN(i, i) = 1.d0
    !EPN(i, :) = VPN(i, :)

    !do j = 1, i-1
      !do k = 1, j
        !EPN(i, k) = EPN(i, k) - EPN(j, k) *projector_Gamma(EPN(j,:), VPN(i, :), OmegaMin, OmegaMax, &
          !& -OmegaMax, -OmegaMin)
      !enddo
    !enddo

    !norm = dsqrt(projector_Gamma(EPN(i, :), EPN(i,:), OmegaMin, OmegaMax, &
      !& -OmegaMax, -OmegaMin))
    !EPN(i, :) = EPN(i, :)/norm 
  !enddo

  !do i = 1, nbasisGamma
    !VNP(i, i) = 1.d0
    !ENP(i, :) = VNP(i, :)

    !do j = 1, i-1
      !do k = 1, j
        !ENP(i, k) = ENP(i, k) - ENP(j, k) *projector_Gamma(ENP(j,:), VNP(i, :), -OmegaMax, -OmegaMin, &
          !& OmegaMin, OmegaMax)
      !enddo
    !enddo

    !norm = dsqrt(projector_Gamma(ENP(i, :), ENP(i,:), -OmegaMax, -OmegaMin, &
      !& OmegaMin, OmegaMax))
    !ENP(i, :) =  ENP(i, :)/norm
  !enddo

  !do i = 1, nbasisGamma
    !VPPR(i, i) = 1.d0
    !EPPR(i, :) = VPPR(i, :)

    !do j = 1, i-1
      !do k = 1, j
        !EPPR(i, k) = EPPR(i, k) - EPPR(j, k) *projector_Gamma(EPPR(j,:), VPPR(i, :), OmegaMin, OmegaMax, &
          !& OmegaMin, 0)
      !enddo
    !enddo

    !norm = dsqrt(projector_Gamma(EPPR(i, :), EPPR(i,:), OmegaMin, OmegaMax, &
      !& OmegaMin, 0))
    !EPPR(i, :) = EPPR(i, :)/norm 
  !enddo

  !do i = 1, nbasisGamma
    !VPPL(i, i) = 1.d0
    !EPPL(i, :) = VPPL(i, :)

    !do j = 1, i-1
      !do k = 1, j
        !EPPL(i, k) = EPPL(i, k) - EPPL(j, k) *projector_Gamma(EPPL(j,:), VPPL(i, :), OmegaMin, OmegaMax, &
          !& 0, OmegaMax)
      !enddo
    !enddo

    !norm = dsqrt(projector_Gamma(EPPL(i, :), EPPL(i,:), OmegaMin, OmegaMax, &
      !& 0, OmegaMax))
    !EPPL(i, :) = EPPL(i, :)/norm 
  !enddo

  !do i = 1, nbasisGamma
    !VNNR(i, i) = 1.d0
    !ENNR(i, :) = VNNR(i, :)

    !do j = 1, i-1
      !do k = 1, j
        !ENNR(i, k) = ENNR(i, k) - ENNR(j, k) *projector_Gamma(ENNR(j,:), VNNR(i, :), -OmegaMax, -OmegaMin, &
          !& -OmegaMax, 0)
      !enddo
    !enddo

    !norm = dsqrt(projector_Gamma(ENNR(i, :), ENNR(i,:), -OmegaMax, -OmegaMin, &
      !& -OmegaMax, 0))
    !ENNR(i, :) = ENNR(i, :)/norm
  !enddo

  !do i = 1, nbasisGamma
    !VNNL(i, i) = 1.d0
    !ENNL(i, :) = VNNL(i, :)

    !do j = 1, i-1
      !do k = 1, j
        !ENNL(i, k) = ENNL(i, k) - ENNL(j, k) *projector_Gamma(ENNL(j,:), VNNL(i, :), -OmegaMax, -OmegaMin, &
          !& 0, -OmegaMin)
      !enddo
    !enddo

    !norm = dsqrt(projector_Gamma(ENNL(i, :), ENNL(i,:), -OmegaMax, -OmegaMin, &
      !& 0, -OmegaMin))
    !ENNL(i, :) = ENNL(i, :)/norm
  !enddo
  !return
!END SUBROUTINE calculate_Gamma_basis




!!---------- weight calculate f(omega) ------------------
!DOUBLE PRECISION FUNCTION weight_basis(Coef, n)
  !implicit none
  !double precision :: Coef(nbasis)
  !integer :: j, n
  !weight_basis = 0.d0
  !do j = 1, nbasis
    !weight_basis = weight_basis + Coef(j)*OriginalBasis(j, n)
  !enddo
  !return 
!END FUNCTION weight_basis


!!---------- projector: \int_{omega} f1(omega)*f2(omega)  --------
!DOUBLE PRECISION FUNCTION projector(Coef1, Coef2, OmegaMin, OmegaMax)
  !implicit none
  !double precision :: Coef1(nbasis), Coef2(nbasis)
  !integer :: j, k, n, OmegaMin, OmegaMax

  !projector  = 0.d0
  !do j = 1, nbasis 
    !do k = 1, nbasis
      !do n = OmegaMin, OmegaMax
        !projector = projector + Coef1(j)*OriginalBasis(j, n)*Coef2(k)*OriginalBasis(k, n)
      !enddo
    !enddo
  !enddo
  !return 
!END FUNCTION projector



!!---------- weight calculate f(omega) ------------------
!DOUBLE PRECISION FUNCTION weight_basis_Gamma(Coef, n1, n2)
  !implicit none
  !double precision :: Coef(nbasisGamma)
  !integer :: j, n1, n2
  !weight_basis_Gamma = 0.d0
  !do j = 1, nbasisGamma
    !weight_basis_Gamma = weight_basis_Gamma+Coef(j)*OriginalBasisGamma(j, n1, n2)
  !enddo
  !return 
!END FUNCTION weight_basis_Gamma


!!---------- projector: \int_{omega1, omega2} f1(omega1,omega2)*f2(omega1,omega2)  --------
!DOUBLE PRECISION FUNCTION projector_Gamma(Coef1, Coef2, OmegaMin1, OmegaMax1, OmegaMin2, OmegaMax2)
  !implicit none
  !double precision :: Coef1(nbasisGamma), Coef2(nbasisGamma)
  !integer :: j, k, n1, n2, OmegaMin1, OmegaMax1, OmegaMin2, OmegaMax2

  !projector_Gamma  = 0.d0

  !do j = 1, nbasisGamma
    !do k = 1, nbasisGamma
      !do n1 = OmegaMin1, OmegaMax1
        !if(OmegaMin2==0) then
          !if(n1<OmegaMax1) then
            !do n2 = n1+1, OmegaMax1
              !projector_Gamma = projector_Gamma + Coef1(j)*OriginalBasisGamma(j, n1, n2)*Coef2(k)* &
                !& OriginalBasisGamma(k, n1, n2)*ReweightBasis(n1, n2)
            !enddo
          !endif
        !else if(OmegaMax2==0) then
          !if(n1>OmegaMin1) then
            !do n2 = OmegaMin2, n1-1
              !projector_Gamma = projector_Gamma + Coef1(j)*OriginalBasisGamma(j, n1, n2)*Coef2(k)* &
                !& OriginalBasisGamma(k, n1, n2)*ReweightBasis(n1, n2)
            !enddo
          !endif
        !else
          !do n2 = OmegaMin2, OmegaMax2
            !projector_Gamma = projector_Gamma + Coef1(j)*OriginalBasisGamma(j, n1, n2)*Coef2(k)* &
              !& OriginalBasisGamma(k, n1, n2)*ReweightBasis(n1, n2)
          !enddo
        !endif
      !enddo
    !enddo
  !enddo
  !return 
!END FUNCTION projector_Gamma




!!--------- test subroutine for the Gram-Schmidt basis ------
!SUBROUTINE test_basis(nbasis, OmegaMin, OmegaMax, CoefP, CoefN)
  !implicit none
  !integer, intent(in) :: nbasis, OmegaMin, OmegaMax
  !double precision :: CoefP(nbasis, nbasis), CoefN(nbasis, nbasis)
  !integer :: i, j, k, l, n
  !double precision :: omega, x, y

  !do i = 1, nbasis

    !y = 0.d0
    !do  n = OmegaMin, OmegaMax
      !y = y + weight_basis(CoefP(i,:), n)**2.d0
    !enddo
    !if(dabs(y-1.d0)>1.d-10) write(*, *) i, y, "Positive Basis error!!"

    !y = 0.d0
    !do  n = -OmegaMax, -OmegaMin
      !y = y + weight_basis(CoefN(i,:), n)**2.d0
    !enddo
    !if(dabs(y-1.d0)>1.d-10) write(*, *) i, y, "Negative Basis error!!"
  !enddo

  !do i = 1, nbasis
    !do j = i+1, nbasis
      !y = projector(CoefP(i, :), CoefP(j, :), OmegaMin, OmegaMax)
      !if(dabs(y)>1.d-10) write(*, *) i, j, y, "Positive Basis error!!"
      !y = projector(CoefN(i, :), CoefN(j, :), -OmegaMax, -OmegaMin)
      !if(dabs(y)>1.d-10) write(*, *) i, j, y, "Negative Basis error!!"
    !enddo
  !enddo

!END SUBROUTINE test_basis




!!--------- test subroutine for the Gram-Schmidt basis ------
!SUBROUTINE test_basis_Gamma(nbasisGamma, OmegaMin, OmegaMax, EPN, ENP, &
    !& EPPR, ENNR, EPPL, ENNL)
  !implicit none
  !integer, intent(in) :: nbasisGamma, OmegaMin, OmegaMax
  !double precision :: ENP(nbasisGamma, nbasisGamma), EPN(nbasisGamma,nbasisGamma)
  !double precision :: EPPR(nbasisGamma, nbasisGamma),ENNR(nbasisGamma,nbasisGamma)
  !double precision :: EPPL(nbasisGamma, nbasisGamma),ENNL(nbasisGamma,nbasisGamma)
  !integer :: i, j, k, l, n1, n2
  !double precision :: omega, x, y

  !do i = 1, nbasisGamma
    !y = 0.d0
    !do  n1 = OmegaMin+1, OmegaMax
      !do n2 = OmegaMin, n1-1
        !y = y + weight_basis_Gamma(EPPR(i,:), n1, n2)**2.d0*ReweightBasis(n1, n2)
      !enddo
    !enddo
    !if(dabs(y-1.d0)>1.d-8) write(*, *) i, y-1.d0, "++R Gamma Basis error!!"

    !y = 0.d0
    !do  n1 = OmegaMin, OmegaMax-1
      !do n2 = n1+1, OmegaMax
        !y = y + weight_basis_Gamma(EPPL(i,:), n1, n2)**2.d0*ReweightBasis(n1, n2)
      !enddo
    !enddo
    !if(dabs(y-1.d0)>1.d-8) write(*, *) i, y-1.d0, "++L Gamma Basis error!!"

    !y = 0.d0
    !do  n1 = OmegaMin, OmegaMax
      !do n2 = -OmegaMax, -OmegaMin
        !y = y + weight_basis_Gamma(EPN(i,:), n1, n2)**2.d0*ReweightBasis(n1, n2)
      !enddo
    !enddo
    !if(dabs(y-1.d0)>1.d-8) write(*, *) i, y-1.d0, "+- Gamma Basis error!!"

    !y = 0.d0
    !do  n1 = -OmegaMax, -OmegaMin
      !do n2 = OmegaMin, OmegaMax
        !y = y + weight_basis_Gamma(ENP(i,:), n1, n2)**2.d0*ReweightBasis(n1, n2)
      !enddo
    !enddo
    !if(dabs(y-1.d0)>1.d-8) write(*, *) i, y-1.d0, "-+ Gamma Basis error!!"

    !y = 0.d0
    !do  n1 = -OmegaMax+1, -OmegaMin
      !do n2 = -OmegaMax, n1-1
        !y = y + weight_basis_Gamma(ENNR(i,:), n1, n2)**2.d0*ReweightBasis(n1, n2)
      !enddo
    !enddo
    !if(dabs(y-1.d0)>1.d-8) write(*, *) i, y-1.d0, "--R Gamma Basis error!!"

    !y = 0.d0
    !do  n1 = -OmegaMax, -OmegaMin-1
      !do n2 = n1+1, -OmegaMin
        !y = y + weight_basis_Gamma(ENNL(i,:), n1, n2)**2.d0*ReweightBasis(n1, n2)
      !enddo
    !enddo
    !if(dabs(y-1.d0)>1.d-8) write(*, *) i, y-1.d0, "--L Gamma Basis error!!"
  !enddo


  !do i = 1, nbasisGamma
    !do j = i+1, nbasisGamma
      !y = projector_Gamma(EPN(i, :), EPN(j, :), OmegaMin, OmegaMax,-OmegaMax,-OmegaMin)
      !if(dabs(y)>1.d-8) write(*, *) i, j, y, "+- Gamma Basis error!!"
      !y = projector_Gamma(ENP(i, :), ENP(j, :),-OmegaMax,-OmegaMin, OmegaMin, OmegaMax)
      !if(dabs(y)>1.d-8) write(*, *) i, j, y, "-+ Gamma Basis error!!"
      !y = projector_Gamma(EPPR(i, :), EPPR(j, :), OmegaMin, OmegaMax, OmegaMin, 0)
      !if(dabs(y)>1.d-7) write(*, *) i, j, y, "++R Gamma Basis error!!"
      !y = projector_Gamma(EPPL(i, :), EPPL(j, :), OmegaMin, OmegaMax, 0, OmegaMax)
      !if(dabs(y)>1.d-7) write(*, *) i, j, y, "++L Gamma Basis error!!"
      !y = projector_Gamma(ENNR(i, :), ENNR(j, :),-OmegaMax,-OmegaMin,-OmegaMax, 0)
      !if(dabs(y)>1.d-7) write(*, *) i, j, y, "--R Gamma Basis error!!"
      !y = projector_Gamma(ENNL(i, :), ENNL(j, :),-OmegaMax,-OmegaMin, 0, -OmegaMin)
      !if(dabs(y)>1.d-7) write(*, *) i, j, y, "--L Gamma Basis error!!"
    !enddo
  !enddo

!END SUBROUTINE test_basis_Gamma
!====================================================================

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
    !integer,parameter :: Lx=4,Ly=4
    !integer,parameter :: Omega=2
    !integer,parameter :: NtypeW=7
    !double precision :: WR(NtypeW,0:Lx-1,0:Ly-1,-Omega:Omega,-Omega:Omega)
!end module

!program main
    !use test
    !integer :: ix,iy
    !print *, "Start"
    !WR(:,:,:,:,:)=1
    !WR(1,0,0,:,:)=3
    !write(*,*) "Before"
    !do iy=0,Ly-1
      !write(*,*) WR(1,:,iy,2,2)
    !enddo
    !call FFT(WR,NtypW,(2*Omega+1)**2,1)
    !print *, "First FFT"
    !do iy=0,Ly-1
      !write(*,*) WR(1,:,iy,2,2)
    !enddo
    !call FFT(WR,NtypW,(2*Omega+1)**2,-1)
    !print *, "Second FFT"
    !do iy=0,Ly-1
      !write(*,*) WR(1,:,iy,2,2)
    !enddo

!CONTAINS

SUBROUTINE transfer_W0(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    call FFT_r(W0PF,1,MxT,BackForth)
    call FFT_tau_single(W0PF,1,Lx*Ly,BackForth)
END SUBROUTINE

SUBROUTINE transfer_Gam0(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer :: it1, it2

    call FFT_r(Gam0PF,1,MxT**2,BackForth)

    if(BackForth/=-1) then
      do it2 = 0, MxT-1
        do it1 = 0, MxT-1
          Gam0PF(:,:,it1,it2) = Gam0PF(:,:,it1,it2)* cdexp(dcmplx(0.d0,-Pi*real(it1)/real(MxT)))
          Gam0PF(:,:,it1,it2) = Gam0PF(:,:,it1,it2)* cdexp(dcmplx(0.d0,-Pi*real(it2)/real(MxT)))
        enddo
      enddo

      call FFT_tau_double(Gam0PF,1,Lx*Ly,BackForth)
    else if(BackForth ==-1) then
      call FFT_tau_double(Gam0PF,1,Lx*Ly,BackForth)

      do it2 = 0, MxT-1
        do it1 = 0, MxT-1
          Gam0PF(:,:,it1,it2) = Gam0PF(:,:,it1,it2)* cdexp(dcmplx(0.d0,Pi*real(it1)/real(MxT)))
          Gam0PF(:,:,it1,it2) = Gam0PF(:,:,it1,it2)* cdexp(dcmplx(0.d0,Pi*real(it1)/real(MxT)))
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
    call FFT_tau_single(W,NtypeW,Lx*Ly,BackForth)
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

      call FFT_tau_double(Gam,NtypeGam,Lx*Ly,BackForth)
    else if(BackForth ==-1) then
      call FFT_tau_double(Gam,NtypeGam,Lx*Ly,BackForth)

      do it2 = 0, MxT-1
        do it1 = 0, MxT-1
          Gam(:,:,:,it1,it2) = Gam(:,:,:,it1,it2)* cdexp(dcmplx(0.d0,Pi*real(it1)/real(MxT)))
          Gam(:,:,:,it1,it2) = Gam(:,:,:,it1,it2)* cdexp(dcmplx(0.d0,Pi*real(it1)/real(MxT)))
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
    call FFT_tau_single(Chi,1,Lx*Ly,BackForth)
END SUBROUTINE

!SUBROUTINE transfer_Polar_r(BackForth)
    !implicit none
    !integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    !call FFT_r(Polar,1,MxT,BackForth)
!END SUBROUTINE

!SUBROUTINE transfer_Polar_t(BackForth)
    !implicit none
    !integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    !call FFT_tau_single(Polar,1,Lx*Ly,BackForth)
!END SUBROUTINE

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


SUBROUTINE FFT_r(XR,Ntype,Nz,BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer,intent(in) :: Nz
    integer,intent(in) :: Ntype
    complex*16    :: XR(Ntype,0:Lx-1,0:Ly-1,0:Nz-1)
    integer :: Power,Noma
    integer :: ix,iy,iz,it
    double precision,allocatable ::   Real1(:), Im1(:)

    do iy=0,Ly-1
      do ix=0,Lx-1
        if(ix>Lx/2 .and. iy<=Ly/2) then
          XR(:,ix,iy,:)=XR(:,Lx-ix,iy,:)
        else if(ix<=Lx/2 .and. iy>Ly/2) then
          XR(:,ix,iy,:)=XR(:,ix,Ly-iy,:)
        else if(ix>Lx/2 .and. iy>Ly/2)  then
          XR(:,ix,iy,:)=XR(:,Lx-ix,Ly-iy,:)
        endif
      enddo
    enddo

    !do FFT in x direction
    Power=log(Lx*1.d0)/log(2.d0)+1.d-14  
    Noma=2**power
    allocate(Real1(0:Noma-1),Im1(0:Noma-1))
    do iz=0,Nz-1
      do iy=0,Ly-1
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
    power=log(Ly*1.d0)/log(2.d0)+1.d-14  
    Noma=2**Power
    allocate(Real1(0:Noma-1),Im1(0:Noma-1))
    do iz=0,Nz-1
      do ix=0,Lx-1
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




!=================== VISUALIZATION  ==================================
SUBROUTINE DRAW
    IMPLICIT NONE
    DOUBLE PRECISION :: x1,y1, x2,y2, y3, scx, scy, sgn   
    DOUBLE PRECISION :: scydash, ca1,ca2,ra,a1,a2, radian
    DOUBLE PRECISION :: phi1, phi2, pi2, pi4
    INTEGER :: i
    integer :: iWLn,iGLn,Vertex1,Vertex2,Vertex3
    
    pi2=dasin(1.d0)
    pi4=pi2/2.d0
    radian=90.d0/pi2
    scx=500/beta
    scy=400./Lx/Ly
    x1=scx*beta
    y1=scy*Vol
    scydash=scy/40.

    open(11, file='graph.eps')
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
    write(11,*) '/C12 {/Times-Roman findfont 12 scalefont setfont show} def'
    write(11,*) '% Put an arrowhead at point x2 y2,'
    write(11,*) '% pointing away from x1 y1'
    write(11,*) '% Replace x2 y2 with coordinates of arrowbase:'
    write(11,*) '% the point to connect lines to'
    write(11,*) '% ArrowHeadSize gives the size of the arrow'
    write(11,*) '/ArrowHeadSize 20 def'
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
    write(11,*) '2 setlinewidth'
    write(11,*) '5 140 translate'
    write(11,*) '1 1 scale'
    write(11,*) '%%EndSetup'

    do i=1,NWLn;
      iWLn=WLnKey2Value(i)
      Vertex1=NeighLn(1,iWLn)
      x1=scx*T3Vertex(Vertex1)
      y1=scy*site_num(GXVertex(Vertex1),GYVertex(Vertex1))
      !write(*,*) "W"
      !write(*,*) T3Vertex(Vertex1),Vertex1
      Vertex2=NeighLn(2,iWLn)
      !write(*,*) T3Vertex(Vertex2),Vertex2
      x2=scx*T3Vertex(Vertex2)
      y2=scy*site_num(GXVertex(Vertex2),GYVertex(Vertex2))

      if(Vertex1 .ne. Ira .and. Vertex2 .ne. Masha &
          & .and. Vertex1 .ne. Masha .and. Vertex2 .ne. Ira) then
        if(TypeLn(iWLn)<=4) then
          write(11,*) '0 0 0 setrgbcolor'
        else
          write(11,*) '0 1 0 setrgbcolor'
        endif
      else
        if(Vertex1==Ira .or. Vertex2==Ira) write(11,*)'0 1 1 setrgbcolor'
        if(Vertex1==Masha .or. Vertex2==Masha) write(11,*)'0 1 1 setrgbcolor'
      endif
      write(11,791) x1,y1,x2,y2 
    enddo
  
    do i=1,NGLn;
      iGLn=GLnKey2Value(i)
      Vertex1=NeighLn(1,iGLn)
      x1=scx*T3Vertex(Vertex1)
      y1=scy*site_num(GXVertex(Vertex1),GYVertex(Vertex1))
      !write(*,*) "G"
      !write(*,*) "V1",Vertex1
      !write(*,*) T3Vertex(Vertex1),GXVertex(Vertex1),GYVertex(Vertex1)
      Vertex2=NeighLn(2,iGLn)
      !write(*,*) "V2",Vertex2
      x2=scx*T3Vertex(Vertex2)
      y2=scy*site_num(GXVertex(Vertex2),GYVertex(Vertex2))
      !write(*,*) T3Vertex(Vertex2),GXVertex(Vertex2),GYVertex(Vertex2)

      if(TypeLn(iGLn)==1) then
        write(11,*) '1 0 0 setrgbcolor'
      else
        write(11,*) '0 0 1 setrgbcolor'
      endif


      if(Vertex1/=Vertex2) then
        if(x1==x2) x2=500 !unphysical situation

        ra=dsqrt((x2-x1)**2+(y2-y1)**2)/2.d0
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
        y3=scy*site_num(GXVertex(Vertex3),GYVertex(Vertex3))
        sgn=1.d0
        IF(y3>y1) sgn=-1.d0   
        write(11,780) x1, y1+sgn*scy/5. , scy/5. 
      endif
    enddo

    do i=1,NVertex
      Vertex1=VertexKey2Value(i)
      x1=scx*T3Vertex(Vertex1)
      y1=scy*site_num(GXVertex(Vertex1),GYVertex(Vertex1))
      if(Vertex1==MeasureGam) then
        write(11,*) '0 1 0 setrgbcolor'
        write(11,777) x1, y1, scy/15.
        if(Vertex1==Ira .or. Vertex1==Masha) then
          write(11,*) '1 0 0 setrgbcolor'
          write(11,777) x1, y1, scy/20.
        endif
      else
        if(Vertex1==Ira .or. Vertex1==Masha) then
          write(11,*) '1 0 0 setrgbcolor'
          write(11,777) x1, y1, scy/20.
        else
          write(11,*) '0 0 0 setrgbcolor'
          write(11,777) x1, y1, scy/20.
        endif
      endif
      write(11,801) x1-5., y1+7., i
      write(11,802) 500.0,650.0-i*10,i,T3Vertex(Vertex1), &
         & GXVertex(Vertex1),GYVertex(Vertex2)
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
 791  format ('Ndashed ',f6.1,x,f6.1,' M',f6.1,x,f6.1,' L')
 801  format (f6.1,x,f6.1,x,' M (',i2,') C12')
 802  format (f6.1,x,f6.1,x,' M (',i2,':',f6.3,'; (',i3,',',i3,')) C12')
 
END SUBROUTINE DRAW

!__________________________________________ 
!   	'0 0 0 setrgbcolor'    black
!     '0 0 1 setrgbcolor'    blue
!     '1 0 0 setrgbcolor'    red
!     '0 1 0 setrgbcolor'    green

INTEGER FUNCTION site_num(X,Y)
  implicit none
  integer :: x,y
  site_num=y*Ly+x+1
  return
END FUNCTION site_num




!============   initialize diagrams =================================
SUBROUTINE first_order_diagram
  implicit none
  integer :: i
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
  Hash4G(kLn(1)) = 1
  kLn(2) = 32
  Hash4G(kLn(2)) = 1
  kLn(3) = add_k(kLn(2), -kLn(1))
  Hash4W(abs(kLn(3))) = 1
  kLn(4) = 23
  Hash4G(kLn(4)) = 1
  kLn(5) = add_k(kLn(3), kLn(4))
  Hash4G(kLn(5)) = 1
  kLn(6) = add_k(kLn(1), -kLn(4))
  Hash4W(abs(kLn(6))) = 1


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

  ! the site variables of Gamma 1, 2, 3 and 4
  GXVertex(1:4) = 0;            GYVertex(1:4) = 0
  WXVertex(1:4) = 0;            WYVertex(1:4) = 0

  ! the time variables of Gamma 1, 2, 3 and 4
  T1Vertex(1)=0.d0; T2Vertex(1)=0.d0; T3Vertex(1)=0.d0
  T1Vertex(2)=0.1; T2Vertex(2)=0.1; T3Vertex(2)=0.1
  T1Vertex(3)=0.2; T2Vertex(3)=0.2; T3Vertex(3)=0.2
  T1Vertex(4)=0.3; T2Vertex(4)=0.3; T3Vertex(4)=0.3

  ! Direction of Gamma: 1: left of W;  2: right of W
  DirecVertex(1) = 1;                DirecVertex(2) = 2
  DirecVertex(3) = 1;                DirecVertex(4) = 2

  ! NeighVertex(i, j): the ith neighbor line of gamma j
  NeighVertex(1,1) = 2;        NeighVertex(2,1) = 1;        NeighVertex(3,1) = 3
  NeighVertex(1,2) = 4;        NeighVertex(2,2) = 5;        NeighVertex(3,2) = 3
  NeighVertex(1,3) = 1;        NeighVertex(2,3) = 4;        NeighVertex(3,3) = 6
  NeighVertex(1,4) = 5;        NeighVertex(2,4) = 2;        NeighVertex(3,4) = 6

  ! weights for lines and vertexes
  WeightLn(1) = weight_line(StatusLn(1),-1, 1,0,0,0.d0,TypeLn(1))
  WeightLn(2) = weight_line(StatusLn(2),-1, 1,0,0,0.d0,TypeLn(2))
  WeightLn(3) = weight_line(StatusLn(3), 0, 2,0,0,0.d0,TypeLn(3))
  WeightLn(4) = weight_line(StatusLn(4),-1, 1,0,0,0.d0,TypeLn(4))
  WeightLn(5) = weight_line(StatusLn(5),-1, 1,0,0,0.d0,TypeLn(5))
  WeightLn(6) = weight_line(StatusLn(6), 0, 2,0,0,0.d0,TypeLn(6))

  WeightVertex(1) = weight_vertex(StatusVertex(1), 1, 0, 0, 0.d0, 0.d0, TypeVertex(1))
  WeightVertex(2) = weight_vertex(StatusVertex(2), 1, 0, 0, 0.d0, 0.d0, TypeVertex(2))
  WeightVertex(3) = weight_vertex(StatusVertex(3), 1, 0, 0, 0.d0, 0.d0, TypeVertex(3))
  WeightVertex(4) = weight_vertex(StatusVertex(4), 1, 0, 0, 0.d0, 0.d0, TypeVertex(4))


  ratio = CoefOfWeight(1)*(1.d0/Beta)**Order *SignFermiLoop
  Anew = d_times_cd(ratio, WeightLn(1)*WeightLn(2)*WeightLn(3)*WeightLn(4)* &
    & WeightLn(5)*WeightLn(6)*WeightVertex(1)*WeightVertex(2)*WeightVertex(3)* &
    & WeightVertex(4))

  WeightCurrent = abs(Anew)
  Phase = Anew/WeightCurrent
  !-------------------------------------------------------
end SUBROUTINE first_order_diagram




SUBROUTINE first_order_diagram_with_bubble
  implicit none
  integer :: i
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
  kLn(1) = 100
  Hash4G(kLn(1)) = 1
  kLn(2) = 32
  Hash4G(kLn(2)) = 1
  kLn(3) = add_k(kLn(1), -kLn(1))
  Hash4W(abs(kLn(3))) = 1
  kLn(4) = 23
  Hash4G(kLn(4)) = 1
  kLn(5) = add_k(kLn(4), kLn(3))
  Hash4G(kLn(5)) = 1
  kLn(6) = add_k(kLn(4), -kLn(2))
  Hash4W(abs(kLn(6))) = 1


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
  GXVertex(1:4) = 0;            GYVertex(1:4) = 0
  WXVertex(1:4) = 0;            WYVertex(1:4) = 0

  ! the time variables of Gamma 1, 2, 3 and 4
  T1Vertex(1)=0.d0; T2Vertex(1)=0.d0; T3Vertex(1)=0.d0
  T1Vertex(2)=0.1; T2Vertex(2)=0.1; T3Vertex(2)=0.1
  T1Vertex(3)=0.2; T2Vertex(3)=0.2; T3Vertex(3)=0.2
  T1Vertex(4)=0.3; T2Vertex(4)=0.3; T3Vertex(4)=0.3

  ! Direction of Gamma: 1: left of W;  2: right of W
  DirecVertex(1) = 1;                DirecVertex(2) = 2
  DirecVertex(3) = 1;                DirecVertex(4) = 2

  ! NeighVertex(i, j): the ith neighbor line of gamma j
  NeighVertex(1,1) = 1;        NeighVertex(2,1) = 1;        NeighVertex(3,1) = 3
  NeighVertex(1,2) = 5;        NeighVertex(2,2) = 4;        NeighVertex(3,2) = 3
  NeighVertex(1,3) = 4;        NeighVertex(2,3) = 2;        NeighVertex(3,3) = 6
  NeighVertex(1,4) = 2;        NeighVertex(2,4) = 5;        NeighVertex(3,4) = 6

  ! weights for lines and vertexes
  WeightLn(1) = weight_line(StatusLn(1),-1, 1,0,0,0.d0,TypeLn(1))
  WeightLn(2) = weight_line(StatusLn(2),-1, 1,0,0,0.d0,TypeLn(2))
  WeightLn(3) = weight_line(StatusLn(3), 0, 2,0,0,0.d0,TypeLn(3))
  WeightLn(4) = weight_line(StatusLn(4),-1, 1,0,0,0.d0,TypeLn(4))
  WeightLn(5) = weight_line(StatusLn(5),-1, 1,0,0,0.d0,TypeLn(5))
  WeightLn(6) = weight_line(StatusLn(6), 0, 2,0,0,0.d0,TypeLn(6))

  WeightVertex(1) = weight_vertex(StatusVertex(1), 1, 0, 0, 0.d0, 0.d0, TypeVertex(1))
  WeightVertex(2) = weight_vertex(StatusVertex(2), 1, 0, 0, 0.d0, 0.d0, TypeVertex(2))
  WeightVertex(3) = weight_vertex(StatusVertex(3), 1, 0, 0, 0.d0, 0.d0, TypeVertex(3))
  WeightVertex(4) = weight_vertex(StatusVertex(4), 1, 0, 0, 0.d0, 0.d0, TypeVertex(4))


  ratio = CoefOfWeight(1)*(1.d0/Beta)**Order *SignFermiLoop
  Anew = d_times_cd(ratio, WeightLn(1)*WeightLn(2)*WeightLn(3)*WeightLn(4)* &
    & WeightLn(5)*WeightLn(6)*WeightVertex(1)*WeightVertex(2)*WeightVertex(3)* &
    & WeightVertex(4))

  WeightCurrent = abs(Anew)
  Phase = Anew/WeightCurrent
  !-------------------------------------------------------
end SUBROUTINE first_order_diagram_with_bubble




SUBROUTINE second_order_diagram
  implicit none
  integer :: i
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
  kLn(1) = 10
  Hash4G(kLn(1)) = 1
  kLn(2) = 21
  Hash4G(kLn(2)) = 1
  kLn(3) = add_k(kLn(1), -kLn(2))
  Hash4W(abs(kLn(3))) = 1
  kLn(4) = 33
  Hash4G(kLn(4)) = 1
  kLn(5) = add_k(kLn(4), -kLn(3))
  Hash4G(kLn(5)) = 1
  kLn(6) = 46
  Hash4W(abs(kLn(6))) = 1
  kLn(7) = add_k(kLn(4), -kLn(6))
  Hash4G(kLn(7)) = 1
  kLn(8) = add_k(kLn(2), kLn(6))
  Hash4G(kLn(8)) = 1
  kLn(9) = add_k(kLn(1), -kLn(8))
  Hash4W(abs(kLn(9))) = 1


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
  GXVertex(1:6) = 0;            GYVertex(1:6) = 0
  WXVertex(1:6) = 0;            WYVertex(1:6) = 0

  ! the time variables of Gamma 1, 2, 3 and 4
  T1Vertex(1)=0.d0;  T2Vertex(1)=0.d0;  T3Vertex(1)=0.d0
  T1Vertex(2)=0.1d0; T2Vertex(2)=0.1d0; T3Vertex(2)=0.1d0
  T1Vertex(3)=0.2d0; T2Vertex(3)=0.2d0; T3Vertex(3)=0.2d0
  T1Vertex(4)=0.3d0; T2Vertex(4)=0.3d0; T3Vertex(4)=0.3d0
  T1Vertex(5)=0.4d0; T2Vertex(5)=0.4d0; T3Vertex(5)=0.4d0
  T1Vertex(6)=0.5d0; T2Vertex(6)=0.5d0; T3Vertex(6)=0.5d0

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
  WeightLn(1) = weight_line(StatusLn(1),-1, 1,0,0,0.d0,TypeLn(1))
  WeightLn(2) = weight_line(StatusLn(2),-1, 1,0,0,0.d0,TypeLn(2))
  WeightLn(3) = weight_line(StatusLn(3), 0, 2,0,0,0.d0,TypeLn(3))
  WeightLn(4) = weight_line(StatusLn(4),-1, 1,0,0,0.d0,TypeLn(4))
  WeightLn(5) = weight_line(StatusLn(5),-1, 1,0,0,0.d0,TypeLn(5))
  WeightLn(6) = weight_line(StatusLn(6), 0, 2,0,0,0.d0,TypeLn(6))
  WeightLn(7) = weight_line(StatusLn(7),-1, 1,0,0,0.d0,TypeLn(7))
  WeightLn(8) = weight_line(StatusLn(8),-1, 1,0,0,0.d0,TypeLn(8))
  WeightLn(9) = weight_line(StatusLn(9), 0, 2,0,0,0.d0,TypeLn(9))

  WeightVertex(1) = weight_vertex(StatusVertex(1), 1, 0, 0, 0.d0, 0.d0, TypeVertex(1))
  WeightVertex(2) = weight_vertex(StatusVertex(2), 1, 0, 0, 0.d0, 0.d0, TypeVertex(2))
  WeightVertex(3) = weight_vertex(StatusVertex(3), 1, 0, 0, 0.d0, 0.d0, TypeVertex(3))
  WeightVertex(4) = weight_vertex(StatusVertex(4), 1, 0, 0, 0.d0, 0.d0, TypeVertex(4))
  WeightVertex(5) = weight_vertex(StatusVertex(5), 1, 0, 0, 0.d0, 0.d0, TypeVertex(5))
  WeightVertex(6) = weight_vertex(StatusVertex(6), 1, 0, 0, 0.d0, 0.d0, TypeVertex(6))


  ratio = CoefOfWeight(2)*(1.d0/Beta)**Order *SignFermiLoop
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
