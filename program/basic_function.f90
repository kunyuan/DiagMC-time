!============== BASIC PROCEDURES IN UPDATE ==========================

!--------- exchange the location of Ira and Masha  -----
SUBROUTINE switch_ira_and_masha
  implicit none
  integer  :: k
  if(IsWormPresent .eqv. .false.)  return
  k=Ira;   Ira=Masha; Masha=k
  kMasha=-kMasha;     OmegaMasha=-OmegaMasha;  SpinMasha=-SpinMasha
END SUBROUTINE switch_ira_and_masha

!====================================================================


!=====================================================================
!==================FAKE WEIGHT FOR MEASURING AND WORM=================
!=====================================================================


SUBROUTINE calculate_GamNormWeight
  implicit none
  integer :: omega1, omegaW, omega2, ityp
  double precision :: Gam0

  GamNormWeight = 0.d0
  do omega1 = -MxOmegaDiag, MxOmegaDiag
    do omega2 = -MxOmegaDiag, MxOmegaDiag
      omegaW = omega1-omega2
      do ityp = 1, 6
        if(ityp==3 .or. ityp==4) cycle
        Gam0 = weight_meas_W(0,0,omegaW,ityp)
        Gam0 = Gam0 *weight_meas_Gamma(0,0,0,omega1,omega2,1)
        Gam0 = Gam0 *weight_Gamma0(0,0,omega1,omega2,1)
        GamNormWeight = GamNormWeight + Gam0
      enddo
    enddo
  enddo
  return
END SUBROUTINE calculate_GamNormWeight

!--------- worm weight function  ---------
DOUBLE PRECISION FUNCTION weight_worm(dxg, dyg, dxw, dyw, domega)
  implicit none
  integer, intent(in)  :: dxg, dyg, dxw, dyw, domega

  !weight_worm = 1.d0
  weight_worm = 1.d0/(1.d0+abs(domega)**2.d0)
  weight_worm = weight_worm*exp(-0.5d0*(abs(dxg)+abs(dyg)))
  weight_worm = weight_worm*exp(-0.5d0*(abs(dxw)+abs(dyw)))

  return
END FUNCTION weight_worm


!--------- worm weight function for W ---------
DOUBLE PRECISION FUNCTION weight_worm_W(dx, dy, omega1, ityp)
  implicit none
  integer, intent(in)  :: dx, dy, omega1, ityp

  weight_worm_W = weight_W(dx, dy, omega1, 1)
  return
END FUNCTION weight_worm_W

!--------- worm weight function for Gamma ---------
DOUBLE PRECISION FUNCTION weight_worm_Gamma(dx, dy, omega1, omega2, ityp)
  implicit none
  integer, intent(in)  :: dx, dy, omega1, omega2, ityp

  weight_worm_Gamma = weight_Gamma0(dx, dy, omega1, omega2, 1)
  return
END FUNCTION weight_worm_Gamma

!--------- measure weight function for Gamma ---------
DOUBLE PRECISION FUNCTION weight_meas_G(omega1, ityp)
  implicit none
  integer, intent(in)  :: omega1, ityp
  weight_meas_G = 1.d0
  return
END FUNCTION weight_meas_G

DOUBLE PRECISION FUNCTION weight_meas_W(dx, dy, omega1, ityp)
  implicit none
  integer, intent(in)  :: dx, dy, omega1, ityp

  weight_meas_W = 1.d0
  return
END FUNCTION weight_meas_W

DOUBLE PRECISION FUNCTION weight_meas_Gamma(iorder, dx, dy, omega1, omega2, ityp)
  implicit none
  integer, intent(in)  :: iorder, dx, dy, omega1, omega2, ityp

  weight_meas_Gamma = 0.d0
  if(dx==0 .and. dy==0) then
    if(ityp ==1 .or. ityp == 2) then
      weight_meas_Gamma = 1.d0
    else if(ityp == 3 .or. ityp == 4) then
      weight_meas_Gamma = 0.d0
    else if(ityp == 5 .or. ityp == 6) then
      weight_meas_Gamma = 1.d0
    endif
  endif

  if(iorder==0) then
    weight_meas_Gamma = weight_meas_Gamma *1.d0/((omega1**2.d0+1.d0)*(omega2**2.d0 +1.d0))
  else if(iorder>=1) then
    weight_meas_Gamma = weight_meas_Gamma *1.d0/((abs(omega1)+1.d0)*(abs(omega2)+1.d0))
  endif
  return
END FUNCTION weight_meas_Gamma



!------------- definition of the system symmetry ----------------
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

SUBROUTINE update_WeightCurrent
  implicit none
  integer :: i, Gam1, Gam2, G, W
  double precision :: weight
  double precision :: wln(MxNLn), wgam(MxNVertex)

  weight = 1.d0
  do i = 1, MxNLn
    if(StatusLn(i)<0)  cycle
    if(KindLn(i)==1) then
      wln(i) = weight_line(StatusLn(i),1,0,0,OmegaLn(i),TypeLn(i))
    else
      Gam1 = NeighLn(1,i);       Gam2 = NeighLn(2,i)
      wln(i) = weight_line(StatusLn(i),2,WXVertex(Gam1)-WXVertex(Gam2), &
        & WYVertex(Gam1)-WYVertex(Gam2), OmegaLn(i), TypeLn(i))
    endif
    WeightLn(i) = wln(i)
    weight = weight *wln(i)
  enddo

  do i = 1, MxNVertex
    if(StatusVertex(i)<0)  cycle
    G = NeighVertex(1, i);           W = NeighVertex(3, i)
    wgam(i) = weight_vertex(Order, StatusVertex(i),GXVertex(i)-WXVertex(i), &
      & GYVertex(i)-WYVertex(i), OmegaLn(G), OmegaLn(W), TypeVertex(i))
    WeightVertex(i) = wgam(i)
    weight = weight *wgam(i)
  enddo

  weight = weight*CoefOfWeight(Order)*(1.d0/Beta)**Order *(-1.d0)**NFermiLoop

  WeightCurrent = weight

  return
END SUBROUTINE update_WeightCurrent
!====================================================================
!============== ELEMENTS GENERATORS =================================
!====================================================================

!------------- definition of the probabilities ----------------
SUBROUTINE def_prob
  implicit none
  integer  :: j, k
  double precision :: Cupdate, CPresent, CAbsent

  !-------- probability for updates --------
  do k = 1, Nupdate
    Fupdate(k) = 0.d0
    do j = 1, k
      Fupdate(k)=Fupdate(k)+Pupdate(j)
    enddo
  enddo

  Cupdate = Fupdate(Nupdate)
  do k = 1, Nupdate
    Pupdate(k) = Pupdate(k)/Cupdate
    Fupdate(k) = Fupdate(k)/Cupdate
  enddo

  !--------- probability for omega ----------------------
  do j = -4, 4
    !Pomega(j) = exp(-1.d0*abs(j))
    Pomega(j) = 1.d0/(abs(j)-0.1)**2.d0

    Fomega(j) = 0.d0
    do k = -4, j
      Fomega(j) = Fomega(j) + Pomega(k)
    enddo
  enddo

  do j = -4, 4
    Pomega(j) = Pomega(j)/Fomega(4)
    Fomega(j) = Fomega(j)/Fomega(4)
  enddo


  !--------- probability for x ----------------------
  Px(-1) =1.d0/3.d0
  Px( 0) =1.d0/3.d0
  Px( 1) =1.d0/3.d0

  !--------- probability for y ----------------------
  Py(-1) =1.d0/3.d0
  Py( 0) =1.d0/3.d0
  Py( 1) =1.d0/3.d0
END SUBROUTINE def_prob


!---------- int omega -------------------------
INTEGER FUNCTION generate_omega()
  implicit none
  integer :: nr
  double precision :: rand

  rand = rn()
  generate_omega = -4
  do while(rand>=Fomega(generate_omega))
    generate_omega=generate_omega+1
  enddo

  return
END FUNCTION generate_omega

DOUBLE PRECISION FUNCTION prob_omega(omega)
  implicit none 
  integer, intent(in) :: omega
  prob_omega = Pomega(omega)
  return
END FUNCTION prob_omega

LOGICAL FUNCTION Is_delta_omega_not_valid(omega)
  implicit none
  integer,intent(in) :: omega
  if(abs(omega)<=4)  then
    Is_delta_omega_not_valid = .false.
  else
    Is_delta_omega_not_Valid = .true.
  endif
END FUNCTION Is_delta_omega_not_valid

LOGICAL FUNCTION Is_omega_not_valid(omega)
  implicit none
  integer,intent(in) :: omega
  if(abs(omega)<=MxOmega)  then
    Is_omega_not_valid = .false.
  else
    Is_Omega_not_Valid = .true.
  endif
END FUNCTION Is_omega_not_valid

!---------- int x y -------------------------
INTEGER FUNCTION generate_x()
  implicit none
  integer :: nr
  nr = Floor(rn()*3.d0)-1
  generate_x = nr
  return
END FUNCTION generate_x

INTEGER FUNCTION generate_y()
  implicit none
  integer :: nr
  nr = Floor(rn()*3.d0)-1
  generate_y = nr
  return
END FUNCTION generate_y

DOUBLE PRECISION FUNCTION prob_x(x)
  implicit none 
  integer, intent(in) :: x
  prob_x = Px(x)
  return
END FUNCTION prob_x

DOUBLE PRECISION FUNCTION prob_y(y)
  implicit none 
  integer, intent(in) :: y
  prob_y = Py(y)
  return
END FUNCTION prob_y

INTEGER FUNCTION find_neigh_x(x, dx)
  implicit none
  integer, intent(in) :: x, dx
  find_neigh_x = x+dx
  if(find_neigh_x<0) then
    find_neigh_x = find_neigh_x + Lx
  else if(find_neigh_x>=Lx) then
    find_neigh_x = find_neigh_x - Lx
  endif
END FUNCTION find_neigh_x

INTEGER FUNCTION find_neigh_y(y, dy)
  implicit none
  integer, intent(in) :: y, dy
  find_neigh_y = y+dy
  if(find_neigh_y<0) then
    find_neigh_y = find_neigh_y + Ly
  else if(find_neigh_y>=Ly) then
    find_neigh_y = find_neigh_y - Ly
  endif
END FUNCTION find_neigh_y

INTEGER FUNCTION diff_x(x1, x2)
  implicit none
  integer, intent(in) :: x1, x2
  diff_x = x1 - x2
  if(diff_x<0)     diff_x = Lx+diff_x
END FUNCTION diff_x

INTEGER FUNCTION diff_y(y1, y2)
  implicit none
  integer, intent(in) :: y1, y2
  diff_y = y1 - y2
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


!---------- int k -------------------------
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
  generate_gline = Ln4GList(rand)
  return
END FUNCTION generate_gline

!----------- randomly pick a wline -------------
INTEGER FUNCTION generate_wline()
  implicit none
  integer :: rand
  rand = Floor(rn()*NWLn)+1
  !if(NWLn/=Order+1) then
    !write(*, *) Order, NWLn
  !endif
  generate_wline = Ln4WList(rand)
  return
END FUNCTION generate_wline

!----------- randomly pick a gamma -------------
INTEGER FUNCTION generate_gamma()
  implicit none
  integer :: rand
  rand = Floor(rn()*NGam)+1
  generate_gamma = Vertex4GamList(rand)
  return
END FUNCTION generate_gamma
!!=======================================================================
!!=======================================================================
!!=======================================================================



!============== GRAM-SCHMIDT BASIS ==================================

!-------------- generate the Gram-Schmidt basis --------------
SUBROUTINE calculate_basis(nbasis, OmegaMin, OmegaMax, EP, EN)
  implicit none
  integer, intent(in):: nbasis, OmegaMin, OmegaMax
  double precision :: EP(nbasis, nbasis), EN(nbasis, nbasis)
  double precision :: VP(nbasis, nbasis), VN(nbasis, nbasis)
  double precision :: norm
  integer :: i, j, k

  EP(:,:) = 0.d0
  EN(:,:) = 0.d0
  VP(:,:) = 0.d0
  VN(:,:) = 0.d0

  do i = 1, nbasis
    VP(i, i) = 1.d0
    EP(i, :) = VP(i, :)

    do j = 1, i-1
      do k = 1, j
        EP(i, k) = EP(i, k) - EP(j, k) *projector(EP(j,:), VP(i, :), OmegaMin, OmegaMax)
      enddo
    enddo

    norm = dsqrt(projector(EP(i,:), EP(i,:), OmegaMin, OmegaMax))
    EP(i, :) = EP(i, :)/norm
  enddo

  do i = 1, nbasis
    VN(i, i) = 1.d0
    EN(i, :) = VN(i, :)

    do k = 1, j
      do j = 1, i-1
        EN(i, k) = EN(i, k) - EN(j, k) *projector(EN(j,:), VN(i, :), -OmegaMax, -OmegaMin)
      enddo
    enddo
    norm = dsqrt(projector(EN(i,:), EN(i,:), -OmegaMax, -OmegaMin))
    EN(i, :) = EN(i, :)/norm
  enddo

END SUBROUTINE calculate_basis



SUBROUTINE calculate_Gamma_basis(nbasisGamma, OmegaMin, OmegaMax, EPN, ENP, &
  &  EPPR, ENNR, EPPL, ENNL)
  implicit none
  integer, intent(in):: nbasisGamma, OmegaMin, OmegaMax
  double precision :: ENP(nbasisGamma, nbasisGamma),  EPN(nbasisGamma, nbasisGamma)
  double precision :: EPPR(nbasisGamma, nbasisGamma), ENNR(nbasisGamma, nbasisGamma)
  double precision :: EPPL(nbasisGamma, nbasisGamma), ENNL(nbasisGamma, nbasisGamma)
  double precision :: VNP(nbasisGamma, nbasisGamma),  VPN(nbasisGamma, nbasisGamma)
  double precision :: VPPR(nbasisGamma, nbasisGamma), VNNR(nbasisGamma, nbasisGamma)
  double precision :: VPPL(nbasisGamma, nbasisGamma), VNNL(nbasisGamma, nbasisGamma)
  double precision :: norm
  integer :: i, j, k
  VPN(:,:) = 0.d0
  VNP(:,:) = 0.d0
  VPPR(:,:) = 0.d0
  VNNR(:,:) = 0.d0
  VPPL(:,:) = 0.d0
  VNNL(:,:) = 0.d0

  EPN(:,:) = 0.d0
  ENP(:,:) = 0.d0
  EPPR(:,:) = 0.d0
  ENNR(:,:) = 0.d0
  EPPL(:,:) = 0.d0
  ENNL(:,:) = 0.d0

  do i = 1, nbasisGamma
    VPN(i, i) = 1.d0
    EPN(i, :) = VPN(i, :)

    do j = 1, i-1
      do k = 1, j
        EPN(i, k) = EPN(i, k) - EPN(j, k) *projector_Gamma(EPN(j,:), VPN(i, :), OmegaMin, OmegaMax, &
          & -OmegaMax, -OmegaMin)
      enddo
    enddo

    norm = dsqrt(projector_Gamma(EPN(i, :), EPN(i,:), OmegaMin, OmegaMax, &
      & -OmegaMax, -OmegaMin))
    EPN(i, :) = EPN(i, :)/norm 
  enddo

  do i = 1, nbasisGamma
    VNP(i, i) = 1.d0
    ENP(i, :) = VNP(i, :)

    do j = 1, i-1
      do k = 1, j
        ENP(i, k) = ENP(i, k) - ENP(j, k) *projector_Gamma(ENP(j,:), VNP(i, :), -OmegaMax, -OmegaMin, &
          & OmegaMin, OmegaMax)
      enddo
    enddo

    norm = dsqrt(projector_Gamma(ENP(i, :), ENP(i,:), -OmegaMax, -OmegaMin, &
      & OmegaMin, OmegaMax))
    ENP(i, :) =  ENP(i, :)/norm
  enddo

  do i = 1, nbasisGamma
    VPPR(i, i) = 1.d0
    EPPR(i, :) = VPPR(i, :)

    do j = 1, i-1
      do k = 1, j
        EPPR(i, k) = EPPR(i, k) - EPPR(j, k) *projector_Gamma(EPPR(j,:), VPPR(i, :), OmegaMin, OmegaMax, &
          & OmegaMin, 0)
      enddo
    enddo

    norm = dsqrt(projector_Gamma(EPPR(i, :), EPPR(i,:), OmegaMin, OmegaMax, &
      & OmegaMin, 0))
    EPPR(i, :) = EPPR(i, :)/norm 
  enddo

  do i = 1, nbasisGamma
    VPPL(i, i) = 1.d0
    EPPL(i, :) = VPPL(i, :)

    do j = 1, i-1
      do k = 1, j
        EPPL(i, k) = EPPL(i, k) - EPPL(j, k) *projector_Gamma(EPPL(j,:), VPPL(i, :), OmegaMin, OmegaMax, &
          & 0, OmegaMax)
      enddo
    enddo

    norm = dsqrt(projector_Gamma(EPPL(i, :), EPPL(i,:), OmegaMin, OmegaMax, &
      & 0, OmegaMax))
    EPPL(i, :) = EPPL(i, :)/norm 
  enddo

  do i = 1, nbasisGamma
    VNNR(i, i) = 1.d0
    ENNR(i, :) = VNNR(i, :)

    do j = 1, i-1
      do k = 1, j
        ENNR(i, k) = ENNR(i, k) - ENNR(j, k) *projector_Gamma(ENNR(j,:), VNNR(i, :), -OmegaMax, -OmegaMin, &
          & -OmegaMax, 0)
      enddo
    enddo

    norm = dsqrt(projector_Gamma(ENNR(i, :), ENNR(i,:), -OmegaMax, -OmegaMin, &
      & -OmegaMax, 0))
    ENNR(i, :) = ENNR(i, :)/norm
  enddo

  do i = 1, nbasisGamma
    VNNL(i, i) = 1.d0
    ENNL(i, :) = VNNL(i, :)

    do j = 1, i-1
      do k = 1, j
        ENNL(i, k) = ENNL(i, k) - ENNL(j, k) *projector_Gamma(ENNL(j,:), VNNL(i, :), -OmegaMax, -OmegaMin, &
          & 0, -OmegaMin)
      enddo
    enddo

    norm = dsqrt(projector_Gamma(ENNL(i, :), ENNL(i,:), -OmegaMax, -OmegaMin, &
      & 0, -OmegaMin))
    ENNL(i, :) = ENNL(i, :)/norm
  enddo
  return
END SUBROUTINE calculate_Gamma_basis




!---------- weight calculate f(omega) ------------------
DOUBLE PRECISION FUNCTION weight_basis(Coef, n)
  implicit none
  double precision :: Coef(nbasis)
  integer :: j, n
  weight_basis = 0.d0
  do j = 1, nbasis
    weight_basis = weight_basis + Coef(j)*OriginalBasis(j, n)
  enddo
  return 
END FUNCTION weight_basis


!---------- projector: \int_{omega} f1(omega)*f2(omega)  --------
DOUBLE PRECISION FUNCTION projector(Coef1, Coef2, OmegaMin, OmegaMax)
  implicit none
  double precision :: Coef1(nbasis), Coef2(nbasis)
  integer :: j, k, n, OmegaMin, OmegaMax

  projector  = 0.d0
  do j = 1, nbasis 
    do k = 1, nbasis
      do n = OmegaMin, OmegaMax
        projector = projector + Coef1(j)*OriginalBasis(j, n)*Coef2(k)*OriginalBasis(k, n)
      enddo
    enddo
  enddo
  return 
END FUNCTION projector



!---------- weight calculate f(omega) ------------------
DOUBLE PRECISION FUNCTION weight_basis_Gamma(Coef, n1, n2)
  implicit none
  double precision :: Coef(nbasisGamma)
  integer :: j, n1, n2
  weight_basis_Gamma = 0.d0
  do j = 1, nbasisGamma
    weight_basis_Gamma = weight_basis_Gamma+Coef(j)*OriginalBasisGamma(j, n1, n2)
  enddo
  return 
END FUNCTION weight_basis_Gamma


!---------- projector: \int_{omega1, omega2} f1(omega1,omega2)*f2(omega1,omega2)  --------
DOUBLE PRECISION FUNCTION projector_Gamma(Coef1, Coef2, OmegaMin1, OmegaMax1, OmegaMin2, OmegaMax2)
  implicit none
  double precision :: Coef1(nbasisGamma), Coef2(nbasisGamma)
  integer :: j, k, n1, n2, OmegaMin1, OmegaMax1, OmegaMin2, OmegaMax2

  projector_Gamma  = 0.d0

  do j = 1, nbasisGamma
    do k = 1, nbasisGamma
      do n1 = OmegaMin1, OmegaMax1
        if(OmegaMin2==0) then
          if(n1<OmegaMax1) then
            do n2 = n1+1, OmegaMax1
              projector_Gamma = projector_Gamma + Coef1(j)*OriginalBasisGamma(j, n1, n2)*Coef2(k)* &
                & OriginalBasisGamma(k, n1, n2)*ReweightBasis(n1, n2)
            enddo
          endif
        else if(OmegaMax2==0) then
          if(n1>OmegaMin1) then
            do n2 = OmegaMin2, n1-1
              projector_Gamma = projector_Gamma + Coef1(j)*OriginalBasisGamma(j, n1, n2)*Coef2(k)* &
                & OriginalBasisGamma(k, n1, n2)*ReweightBasis(n1, n2)
            enddo
          endif
        else
          do n2 = OmegaMin2, OmegaMax2
            projector_Gamma = projector_Gamma + Coef1(j)*OriginalBasisGamma(j, n1, n2)*Coef2(k)* &
              & OriginalBasisGamma(k, n1, n2)*ReweightBasis(n1, n2)
          enddo
        endif
      enddo
    enddo
  enddo
  return 
END FUNCTION projector_Gamma




!--------- test subroutine for the Gram-Schmidt basis ------
SUBROUTINE test_basis(nbasis, OmegaMin, OmegaMax, CoefP, CoefN)
  implicit none
  integer, intent(in) :: nbasis, OmegaMin, OmegaMax
  double precision :: CoefP(nbasis, nbasis), CoefN(nbasis, nbasis)
  integer :: i, j, k, l, n
  double precision :: omega, x, y

  do i = 1, nbasis

    y = 0.d0
    do  n = OmegaMin, OmegaMax
      y = y + weight_basis(CoefP(i,:), n)**2.d0
    enddo
    if(dabs(y-1.d0)>1.d-10) write(*, *) i, y, "Positive Basis error!!"

    y = 0.d0
    do  n = -OmegaMax, -OmegaMin
      y = y + weight_basis(CoefN(i,:), n)**2.d0
    enddo
    if(dabs(y-1.d0)>1.d-10) write(*, *) i, y, "Negative Basis error!!"
  enddo

  do i = 1, nbasis
    do j = i+1, nbasis
      y = projector(CoefP(i, :), CoefP(j, :), OmegaMin, OmegaMax)
      if(dabs(y)>1.d-10) write(*, *) i, j, y, "Positive Basis error!!"
      y = projector(CoefN(i, :), CoefN(j, :), -OmegaMax, -OmegaMin)
      if(dabs(y)>1.d-10) write(*, *) i, j, y, "Negative Basis error!!"
    enddo
  enddo

END SUBROUTINE test_basis




!--------- test subroutine for the Gram-Schmidt basis ------
SUBROUTINE test_basis_Gamma(nbasisGamma, OmegaMin, OmegaMax, EPN, ENP, &
    & EPPR, ENNR, EPPL, ENNL)
  implicit none
  integer, intent(in) :: nbasisGamma, OmegaMin, OmegaMax
  double precision :: ENP(nbasisGamma, nbasisGamma), EPN(nbasisGamma,nbasisGamma)
  double precision :: EPPR(nbasisGamma, nbasisGamma),ENNR(nbasisGamma,nbasisGamma)
  double precision :: EPPL(nbasisGamma, nbasisGamma),ENNL(nbasisGamma,nbasisGamma)
  integer :: i, j, k, l, n1, n2
  double precision :: omega, x, y

  do i = 1, nbasisGamma
    y = 0.d0
    do  n1 = OmegaMin+1, OmegaMax
      do n2 = OmegaMin, n1-1
        y = y + weight_basis_Gamma(EPPR(i,:), n1, n2)**2.d0*ReweightBasis(n1, n2)
      enddo
    enddo
    if(dabs(y-1.d0)>1.d-8) write(*, *) i, y-1.d0, "++R Gamma Basis error!!"

    y = 0.d0
    do  n1 = OmegaMin, OmegaMax-1
      do n2 = n1+1, OmegaMax
        y = y + weight_basis_Gamma(EPPL(i,:), n1, n2)**2.d0*ReweightBasis(n1, n2)
      enddo
    enddo
    if(dabs(y-1.d0)>1.d-8) write(*, *) i, y-1.d0, "++L Gamma Basis error!!"

    y = 0.d0
    do  n1 = OmegaMin, OmegaMax
      do n2 = -OmegaMax, -OmegaMin
        y = y + weight_basis_Gamma(EPN(i,:), n1, n2)**2.d0*ReweightBasis(n1, n2)
      enddo
    enddo
    if(dabs(y-1.d0)>1.d-8) write(*, *) i, y-1.d0, "+- Gamma Basis error!!"

    y = 0.d0
    do  n1 = -OmegaMax, -OmegaMin
      do n2 = OmegaMin, OmegaMax
        y = y + weight_basis_Gamma(ENP(i,:), n1, n2)**2.d0*ReweightBasis(n1, n2)
      enddo
    enddo
    if(dabs(y-1.d0)>1.d-8) write(*, *) i, y-1.d0, "-+ Gamma Basis error!!"

    y = 0.d0
    do  n1 = -OmegaMax+1, -OmegaMin
      do n2 = -OmegaMax, n1-1
        y = y + weight_basis_Gamma(ENNR(i,:), n1, n2)**2.d0*ReweightBasis(n1, n2)
      enddo
    enddo
    if(dabs(y-1.d0)>1.d-8) write(*, *) i, y-1.d0, "--R Gamma Basis error!!"

    y = 0.d0
    do  n1 = -OmegaMax, -OmegaMin-1
      do n2 = n1+1, -OmegaMin
        y = y + weight_basis_Gamma(ENNL(i,:), n1, n2)**2.d0*ReweightBasis(n1, n2)
      enddo
    enddo
    if(dabs(y-1.d0)>1.d-8) write(*, *) i, y-1.d0, "--L Gamma Basis error!!"
  enddo


  do i = 1, nbasisGamma
    do j = i+1, nbasisGamma
      y = projector_Gamma(EPN(i, :), EPN(j, :), OmegaMin, OmegaMax,-OmegaMax,-OmegaMin)
      if(dabs(y)>1.d-8) write(*, *) i, j, y, "+- Gamma Basis error!!"
      y = projector_Gamma(ENP(i, :), ENP(j, :),-OmegaMax,-OmegaMin, OmegaMin, OmegaMax)
      if(dabs(y)>1.d-8) write(*, *) i, j, y, "-+ Gamma Basis error!!"
      y = projector_Gamma(EPPR(i, :), EPPR(j, :), OmegaMin, OmegaMax, OmegaMin, 0)
      if(dabs(y)>1.d-7) write(*, *) i, j, y, "++R Gamma Basis error!!"
      y = projector_Gamma(EPPL(i, :), EPPL(j, :), OmegaMin, OmegaMax, 0, OmegaMax)
      if(dabs(y)>1.d-7) write(*, *) i, j, y, "++L Gamma Basis error!!"
      y = projector_Gamma(ENNR(i, :), ENNR(j, :),-OmegaMax,-OmegaMin,-OmegaMax, 0)
      if(dabs(y)>1.d-7) write(*, *) i, j, y, "--R Gamma Basis error!!"
      y = projector_Gamma(ENNL(i, :), ENNL(j, :),-OmegaMax,-OmegaMin, 0, -OmegaMin)
      if(dabs(y)>1.d-7) write(*, *) i, j, y, "--L Gamma Basis error!!"
    enddo
  enddo

END SUBROUTINE test_basis_Gamma
!====================================================================


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
  call FFT(W0InMoment, NtypW,1, BackForth)
END SUBROUTINE

SUBROUTINE transfer_W(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    call FFT(WR,NtypW,MxOmegaW1+1,BackForth)
    call FFT(WRTailP,NtypW,Nbasis,BackForth)
    call FFT(WRTailC,NtypW,1,BackForth)
END SUBROUTINE

SUBROUTINE transfer_Gamma(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation

    call FFT(GamR,NtypGa,(2*MxOmegaGamG1+1)*(2*MxOmegaGamG1+1),BackForth)

    call FFT(GamRTailNP,NtypGa,NbasisGamma,BackForth)
    call FFT(GamRTailPN,NtypGa,NbasisGamma,BackForth)
    call FFT(GamRTailPPR,NtypGa,NbasisGamma,BackForth)
    call FFT(GamRTailNNR,NtypGa,NbasisGamma,BackForth)
    call FFT(GamRTailPPL,NtypGa,NbasisGamma,BackForth)
    call FFT(GamRTailNNL,NtypGa,NbasisGamma,BackForth)

    call FFT(GamRTailDiagP,NtypGa,Nbasis,BackForth)
    call FFT(GamRTailDiagN,NtypGa,Nbasis,BackForth)

    call FFT(GamRTailMP,NtypGa,Nbasis*(2*MxOmegaGamG1+1),BackForth)
    call FFT(GamRTailMN,NtypGa,Nbasis*(2*MxOmegaGamG1+1),BackForth)
    call FFT(GamRTailPM,NtypGa,Nbasis*(2*MxOmegaGamG1+1),BackForth)
    call FFT(GamRTailNM,NtypGa,Nbasis*(2*MxOmegaGamG1+1),BackForth)
    call FFT(GamRTailC,NtypGa,1,BackForth)
END SUBROUTINE

SUBROUTINE transfer_Pi(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    call FFT(PiR,ntypPi,2*MxOmegaChi+1,BackForth)
END SUBROUTINE

SUBROUTINE transfer_Chi(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    call FFT(ChiR,ntypCh,2*MxOmegaChi+1,BackForth)
END SUBROUTINE


SUBROUTINE FFT(XR,Ntype,Nz,BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer,intent(in) :: Nz
    integer,intent(in) :: Ntype
    double precision :: XR(Ntype,0:Lx-1,0:Ly-1,0:Nz-1)
    integer :: Power,Noma
    integer :: ix,iy,iz,ixt,iyt,it
    double precision,allocatable ::   Real1(:), Im1(:)

    do ix=0,Lx-1
      do iy=0,Ly-1
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
    do it=1,Ntype
      do iy=0,Ly-1
         do iz=0,Nz-1
             Real1(:)=XR(it,:,iy,iz)
             Im1(:)=0.0
             call sffteu(Real1, Im1, Noma, Power, BackForth)
             XR(it,:,iy,iz)=Real1(:)
             if(sum(abs(Im1(:)))>1e-10) then
               write(*,*) sum(abs(Im1(:))), "Where is the imaginary part comes from?"
             endif
         end do
      end do
    enddo
    deallocate(Real1,Im1)

    !do FFT in y direction

    power=log(Ly*1.d0)/log(2.d0)+1.d-14  
    Noma=2**Power
    allocate(Real1(0:Noma-1),Im1(0:Noma-1))
    do it=1,Ntype
      do ix=0,Lx-1
         do iz=0,Nz-1
             Real1(:)=XR(it,ix,:,iz)
             Im1(:)=0.0
             call sffteu(Real1, Im1, Noma, Power, BackForth)
             XR(it,ix,:,iz)=Real1(:)
             if(sum(abs(Im1(:)))>1e-10) then
               write(*,*) sum(abs(Im1(:))), "Where is the imaginary part comes from?"
             endif
         end do
      end do
    enddo
    deallocate(Real1,Im1)

end SUBROUTINE FFT

    
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
!====================================================================
!====================================================================
!====================================================================
