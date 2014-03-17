!!====================== initialization =================================
!!=======================================================================

SUBROUTINE initialize_markov
    implicit none
    integer :: i

    !--------------- initialize variables ---------------
    GLnKey2Value(:) = 0
    WLnKey2Value(:) = 0
    VertexKey2Value(:) = 0

    StatusLn(:) = -1
    StatusVertex(:)= -1

    WeightLn(:) = 1.d0
    WeightVertex(:) = 1.d0

    Hash4G(:) = 0
    Hash4W(:) = 0

    do i = 1, MxNLn
      NextLn(i) = i+1
      if(i==MxNLn)  NextLn(i) = -1
    enddo

    do i = 1, MxNVertex
      NextVertex(i) = i+1
      if(i==MxNVertex)  NextVertex(i) = -1
    enddo

    call def_prob
    call def_spatial_weight
    call def_spin
    call def_diagram
    return
end SUBROUTINE initialize_markov

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
  return
END SUBROUTINE def_prob

SUBROUTINE def_spatial_weight
  implicit none
  double precision :: tt,logL(2)
  integer :: ix,iy,i
  integer,dimension(2) :: L
! initializing logscale sampling of space
! for seeding rx coordinate displacement do this
!     x=rndm(); rx=0.5d0*dexp(x*logL); ix=rx
!     IF(rndm()>0.5d0) ix=L(1)-1-ix
! this gives you displacement ix in the affine basis called with probability
! SpatialWeight(i,ix)
  L(1)=Lx
  L(2)=Ly
  logL(1)=logLx
  logL(2)=logLy

  do i=1,2
    ix=0;
    SpatialWeight(i,ix)=dlog(2.d0)/(2.d0*logL(i))
    SpatialWeight(i,L(i)-1-ix)=SpatialWeight(i,ix)
    do ix=1,L(i)/2-1
      tt=ix*1.d0;
      SpatialWeight(i,ix)=dlog((tt+1.d0)/tt)/(2.d0*logL(i)) 
      SpatialWeight(i,L(i)-1-ix)=SpatialWeight(i,ix)
    enddo
    tt=0.0
    do ix=0,L(i)-1
      tt=tt+SpatialWeight(i,ix)
    enddo 
    do ix=0,L(i)-1
      SpatialWeight(i,ix)=SpatialWeight(i,ix)/tt
    enddo 
  enddo

END SUBROUTINE def_spatial_weight


SUBROUTINE def_spin
  implicit none

  TypeGam2W(:,:) = 0
  TypeSp2Gam(:,:,:,:) = 0

  !-- calculate the type of W according to the neighbor vertexes --
  !---- TypeGam2W(Type(Gamleft), Type(Gamright)) = Type(W)
  TypeGam2W(1,1) = 1;      TypeGam2W(1,4) = 1
  TypeGam2W(4,1) = 1;      TypeGam2W(4,4) = 1

  TypeGam2W(2,2) = 2;      TypeGam2W(2,3) = 2
  TypeGam2W(3,2) = 2;      TypeGam2W(3,3) = 2

  TypeGam2W(1,2) = 3;      TypeGam2W(1,3) = 3
  TypeGam2W(4,2) = 3;      TypeGam2W(4,3) = 3

  TypeGam2W(2,1) = 4;      TypeGam2W(2,4) = 4
  TypeGam2W(3,1) = 4;      TypeGam2W(3,4) = 4

  TypeGam2W(5,6) = 5;      TypeGam2W(6,5) = 6

  !-- calculate the type of Gam according to the neighbor G, W --
  !---- TypeSp2Gam(Type(Gin), Type(Gout), Type(Gamin), Type(Gamout)) = Type(Gam)
  TypeSp2Gam(1,1,1,1) = 1
  TypeSp2Gam(2,2,2,2) = 2
  TypeSp2Gam(1,1,2,2) = 3
  TypeSp2Gam(2,2,1,1) = 4
  TypeSp2Gam(1,2,1,2) = 5
  TypeSp2Gam(2,1,2,1) = 6

END SUBROUTINE def_spin

!------------- definition of the config of diagram ----------------
SUBROUTINE def_diagram
  implicit none
  call first_order_diagram
  return
END SUBROUTINE def_diagram

!=======================================================================
!=======================================================================
!=======================================================================




!=======================================================================
!========================= UPDATES =====================================
!=======================================================================

SUBROUTINE markov
  implicit none
  integer :: iflag, istep,isamp
  double precision :: nr

  do isamp = 1, Nsamp
    istep = 0
    do while(istep < NStep)
      nr=rn()
      iupdate = 1
      do while(nr>Fupdate(iupdate))
        iupdate = iupdate + 1
      enddo

      select case(iupdate)
        case( 1) 
          call create_worm_along_wline          
        case( 2)      
          call delete_worm_along_wline                       
        !case( 3) 
          !call create_worm_along_gline  
        !case( 4) 
          !call delete_worm_along_gline  
        case( 5) 
          call move_worm_along_wline             
        case( 6) 
          call move_worm_along_gline              
        !case( 7) 
          !call add_interaction                  
        !case( 8) 
          !call remove_interaction             
        !case( 9) 
          !call add_interaction_cross              
        !case(10)
          !call remove_interaction_cross                  
        !case(11) 
          !call reconnect                      
        !case(12) 
          !call change_gline_space          
        case(13)  
          call change_wline_space         
        case(14) 
          call change_Gamma_type     
        case(15) 
          call move_measuring_index       
        !case(16)  
          !call change_gline_time        
        !case(17)  
          !call change_wline_time        
        case(18)  
          call change_wline_isdelta       
        case(19)  
          call change_Gamma_isdelta       
      end select

      imc = imc + 1

      if(IsWormPresent) then
        GamWormOrder(Order) = GamWormOrder(Order) + 1
        if(rn()<=0.5d0) call switch_ira_and_masha
      else
        istep = istep + 1
        GamOrder(Order) = GamOrder(Order) + 1
      endif
      if(IsToss==0) call measure
    enddo
    !if(IsToss==0) call measure
  enddo

END SUBROUTINE markov


!====================================================================
!============================== UPDATES =============================
!====================================================================


!---------------- create worm on wline : Pupdate(1) ------------
!------       \iGout    jGout/           \           /     -----
!------   iGam ======>======= jGam  => Ira====>======Masha -----
!------       /iGin      jGin\           /           \     -----
!---------------------------------------------------------------
SUBROUTINE create_worm_along_wline
  implicit none
  integer :: iWLn, iGam, jGam, iGin, jGin, iGout, jGout
  integer :: statiGam, statjGam, statW
  integer :: statiGin, statiGout, statjGin, statjGout
  integer :: k, kiW
  integer :: kiWold, flag
  integer :: spin
  integer :: tempGam
  double precision :: tau, tau1, tau2, Pacc
  double precision :: WWorm
  complex*16 :: WiGin, WiGout, WjGin, WjGout
  complex*16 :: WIra, WMasha, WWNew, Anew, Aold, sgn

  !------------ step1 : check if worm is present ------------------
  if(IsWormPresent .eqv. .true.)   return
  !ProbProp(iupdate) = ProbProp(iupdate) + 1

  !------------ step2 : propose the new config ------------------
  iWLn = generate_wline()

  iGam = NeighLn(1, iWLn);              jGam = NeighLn(2, iWLn)
  iGin = NeighVertex(1, iGam);          jGin = NeighVertex(1, jGam)
  iGout = NeighVertex(2, iGam);         jGout = NeighVertex(2, jGam)

  spin = 4*(Floor(rn()*2.d0))-2         !spin = +/- 2
  !if(spin-2*(TypeLn(jGin)-TypeLn(jGout))==0)  return  !in this case, I/M cannot move

  k = generate_k()
  kiWold = kLn(iWLn)
  kLn(iWLn) = add_k(kiWold, k)

  !------------ step3 : configuration check -------------------
  flag=0
  if(flag==0) then
    if(Is_reducible_W(iWLn))                 flag=1
  endif
  if(flag==0) then
    if(Is_reducible_W_Gam(iWLn))           flag=1
  endif

  if(flag==1) then
    kLn(iWLn) = kiWold
    return
  endif

  statiGam = add_ira_stat(StatusVertex(iGam))
  statjGam = add_ira_stat(StatusVertex(jGam))
  statW    = line_stat(statiGam, statjGam)

  !------------ step4 : weight calculation ----------------------
  tau = T2Vertex(NeighLn(2, iGin))- T1Vertex(NeighLn(1, iGin))
  WiGin  = weight_line(StatusLn(iGin),  IsDeltaLn(iGin), 1, 0, 0, tau, TypeLn(iGin))

  tau = T2Vertex(NeighLn(2, iGout))- T1Vertex(NeighLn(1, iGout))
  WiGout = weight_line(StatusLn(iGout), IsDeltaLn(iGin), 1, 0, 0, tau, TypeLn(iGout))

  tau = T2Vertex(NeighLn(2, jGin))- T1Vertex(NeighLn(1, jGin))
  WjGin  = weight_line(StatusLn(jGin),  IsDeltaLn(iGin), 1, 0, 0, tau, TypeLn(jGin))

  tau = T2Vertex(NeighLn(2, jGout))- T1Vertex(NeighLn(1, jGout))
  WjGout = weight_line(StatusLn(jGout), IsDeltaLn(iGin), 1, 0, 0, tau, TypeLn(jGout))

  tau1 = T3Vertex(iGam)- T2Vertex(iGam)
  tau2 = T1Vertex(iGam)- T3Vertex(iGam)
  WIra = weight_vertex(statiGam, IsDeltaVertex(iGam), GXVertex(iGam)-WXVertex(iGam),  &
    & GYVertex(iGam)-WYVertex(iGam), tau1, tau2, TypeVertex(iGam))

  tau1 = T3Vertex(jGam)- T2Vertex(jGam)
  tau2 = T1Vertex(jGam)- T3Vertex(jGam)
  WMasha = weight_vertex(statjGam, IsDeltaVertex(jGam), GXVertex(jGam)-WXVertex(jGam),  &
    & GYVertex(jGam)-WYVertex(jGam), tau1, tau2, TypeVertex(jGam))

  tau = T3Vertex(jGam) - T3Vertex(iGam)
  WWNew = weight_line(statW, IsDeltaLn(iWLn), 2, WXVertex(jGam)-WXVertex(iGam), &
    & WYVertex(jGam)-WYVertex(iGam), tau, TypeLn(iWLn))

  tau = T3Vertex(jGam) - T3Vertex(iGam)
  WWorm = weight_worm(GXVertex(iGam)-GXVertex(jGam), GYVertex(iGam)-GYVertex(jGam),&
    & WXVertex(iGam)-WXVertex(jGam), WYVertex(iGam)-WYVertex(jGam), tau) 

  Anew = WIra *WMasha *WWNew *WiGin *WiGout *WjGin *WjGout
  Aold = WeightLn(iWLn) *WeightVertex(iGam) *WeightVertex(jGam) *WeightLn(iGin) &
  & *WeightLn(iGout) *WeightLn(jGin) *WeightLn(jGout)

  if(abs(Anew)==0.d0)   return

  if(iGin==jGout)  then
    Anew = Anew/WjGout
    Aold = Aold/WeightLn(jGout)
  endif
  if(jGin==iGout)  then
    Anew = Anew/WjGin
    Aold = Aold/WeightLn(jGin)
  endif
  if(iGin==iGout) then
    Anew = Anew/WiGout
    Aold = Aold/WeightLn(iGout)
  endif
  if(jGin==jGout) then
    Anew = Anew/WjGout
    Aold = Aold/WeightLn(jGout)
  endif

  call weight_ratio(Pacc, sgn, Anew, Aold)
  Pacc = Pacc *Wworm *(Order+1.d0)/(CoefOfWorm*0.5d0) 
  Pacc = Pacc *Pupdate(2)/Pupdate(1)

  !------------ step5 : accept the update -----------------------
  ProbProp(Order, iupdate) = ProbProp(Order, iupdate) + 1

  if(rn()<=Pacc) then

    !-------------- update the diagram info --------------------
    Phase = Phase*sgn
    IsWormPresent = .true.

    !---------- update Ira and Masha -------------------
    Ira = iGam
    Masha = jGam
    kMasha = k
    SpinMasha = spin

    !----------- update k and omega of elements --------
    call update_Hash4W(kiWold, kLn(iWLn))

    !----------- update status of elements -------------
    StatusVertex(Ira)  = statiGam
    StatusVertex(Masha)= statjGam
    StatusLn(iWLn)  = statW

    !----------- update weight of elements -------------
    WeightWorm = WWorm
    WeightVertex(Ira) = WIra
    WeightVertex(Masha) = WMasha
    WeightLn(iWLn) = WWNew

    WeightLn(iGin)  = WiGin
    WeightLn(iGout) = WiGout
    WeightLn(jGin)  = WjGin
    WeightLn(jGout) = WjGout

    call update_weight(Anew, Aold)

    ProbAcc(Order, iupdate) = ProbAcc(Order, iupdate) + 1
  else
    kLn(iWLn) = kiWold
  endif
  return
END SUBROUTINE create_worm_along_wline



!!--------------- delete worm along wline: Pupdate(2) ----------------
!------       \iGout    jGout/           \           /     -----
!------   Ira  ======>======= Masha =>iGam====>======jGam  -----
!------       /iGin      jGin\           /           \     -----
!---------------------------------------------------------------
SUBROUTINE delete_worm_along_wline
  implicit none
  integer :: dir
  integer :: iWLn, iGam, jGam, iGin, jGin, iGout, jGout
  integer :: typIra, typMasha, typW
  integer :: statIra, statMasha, statW
  integer :: statiGin, statiGout, statjGin, statjGout
  integer :: k, kiW, kiWold, flag
  double precision :: Pacc, tau, tau1, tau2
  complex*16 :: WiGin, WiGout, WjGin, WjGout
  complex*16 :: Wi, Wj, WWNew, Anew, Aold, sgn

  !------------ step1 : check if worm is present ------------------
  if(IsWormPresent .eqv. .false.)  return
  !ProbProp(iupdate) = ProbProp(iupdate) + 1

  if(NeighVertex(3, Ira)/= NeighVertex(3, Masha))   return

  !------------ step2 : propose the new config ------------------
  iWLn = NeighVertex(3, Ira)
  dir = DirecVertex(Ira)

  iGin  = NeighVertex(1, Ira);          jGin  = NeighVertex(1, Masha)
  iGout = NeighVertex(2, Ira);          jGout = NeighVertex(2, Masha)

  k = kMasha
  kiWold = kLn(iWLn)
  kLn(iWLn) = add_k(kiWold, (-1)**dir*k)

  !------------ step3 : configuration check ---------------------
  flag=0
  if(flag==0) then
    if(Is_reducible_W(iWLn))   flag=1
  endif
  if(flag==0) then
    if(Is_reducible_W_Gam(iWLn))  flag=1
  endif

  if(flag==1) then
    kLn(iWLn) = kiWold
    return
  endif

  !------------ status, type for the new configuration ------------
  statIra   = delete_ira_stat(StatusVertex(Ira))
  statMasha = delete_ira_stat(StatusVertex(Masha))
  statW     = line_stat(statIra, statMasha)

  typIra = TypeSp2Gam(TypeLn(iGin),TypeLn(iGout),SpInVertex(1, Ira), SpInVertex(2, Ira))
  typMasha = TypeSp2Gam(TypeLn(jGin),TypeLn(jGout),SpInVertex(1,Masha),SpInVertex(2, Masha))
  if(dir==1) then
    typW = TypeGam2W(typIra, typMasha)
  else
    typW = TypeGam2W(typMasha, typIra)
  endif


  !------------ step4 : weight calculation ---------------------
  tau = T2Vertex(NeighLn(2, iGin))- T1Vertex(NeighLn(1, iGin))
  WiGin  = weight_line(StatusLn(iGin),  IsDeltaLn(iGin), 1, 0, 0, tau, TypeLn(iGin))

  tau = T2Vertex(NeighLn(2, iGout))- T1Vertex(NeighLn(1, iGout))
  WiGout = weight_line(StatusLn(iGout), IsDeltaLn(iGin), 1, 0, 0, tau, TypeLn(iGout))

  tau = T2Vertex(NeighLn(2, jGin))- T1Vertex(NeighLn(1, jGin))
  WjGin  = weight_line(StatusLn(jGin),  IsDeltaLn(iGin), 1, 0, 0, tau, TypeLn(jGin))

  tau = T2Vertex(NeighLn(2, jGout))- T1Vertex(NeighLn(1, jGout))
  WjGout = weight_line(StatusLn(jGout), IsDeltaLn(iGin), 1, 0, 0, tau, TypeLn(jGout))

  tau1 = T3Vertex(Ira)- T2Vertex(Ira)
  tau2 = T1Vertex(Ira)- T3Vertex(Ira)
  Wi = weight_vertex(statIra, IsDeltaVertex(Ira), GXVertex(Ira)-WXVertex(Ira),  &
    & GYVertex(Ira)-WYVertex(Ira), tau1, tau2, typIra)

  tau1 = T3Vertex(Masha)- T2Vertex(Masha)
  tau2 = T1Vertex(Masha)- T3Vertex(Masha)
  Wj = weight_vertex(statMasha, IsDeltaVertex(Masha), GXVertex(Masha)-WXVertex(Masha),  &
    & GYVertex(Masha)-WYVertex(Masha), tau1, tau2, typMasha)

  tau = T3Vertex(Masha) - T3Vertex(Ira)
  WWNew = weight_line(statW, IsDeltaLn(iWLn), 2, WXVertex(Masha)-WXVertex(Ira), WYVertex(Masha)- &
    & WYVertex(Ira), tau, typW)

  Anew = Wi *Wj *WWNew *WiGin *WiGout *WjGin *WjGout
  Aold = WeightLn(iWLn) *WeightVertex(Ira) *WeightVertex(Masha) *WeightLn(iGin) &
  & *WeightLn(iGout) *WeightLn(jGin) *WeightLn(jGout)

  if(abs(Anew)==0.d0)   return

  if(iGin==jGout)  then
    Anew = Anew/WjGout
    Aold = Aold/WeightLn(jGout)
  endif

  if(jGin==iGout)  then
    Anew = Anew/WjGin
    Aold = Aold/WeightLn(jGin)
  endif

  if(iGin==iGout) then
    Anew = Anew/WiGout
    Aold = Aold/WeightLn(iGout)
  endif

  if(jGin==jGout) then
    Anew = Anew/WjGout
    Aold = Aold/WeightLn(jGout)
  endif

  call weight_ratio(Pacc, sgn, Anew, Aold)
  Pacc = Pacc *(CoefOfWorm*0.5d0)/WeightWorm/(Order+1.d0) 
  Pacc = Pacc *Pupdate(1)/Pupdate(2)

  ProbProp(Order, iupdate) = ProbProp(Order, iupdate) + 1
  !------------ step5 : accept the update -----------------------
  if(rn()<=Pacc) then

    !---------- update the diagram info ----------------
    Phase = Phase*sgn
    IsWormPresent = .false.

    !----------- update k and omega of elements --------
    call update_Hash4W(kiWold, kLn(iWLn))

    !----------- update type ---------------------------
    TypeVertex(Ira) = typIra
    TypeVertex(Masha) = typMasha
    TypeLn(iWLn) = typW

    !----------- update status of elements -------------
    StatusVertex(Ira)  = statIra
    StatusVertex(Masha)= statMasha
    StatusLn(iWLn)  = statW

    !----------- update weight of elements -------------
    WeightVertex(Ira) = Wi
    WeightVertex(Masha) = Wj

    WeightLn(iWLn) = WWNew
    WeightLn(iGin)  = WiGin
    WeightLn(iGout) = WiGout
    WeightLn(jGin)  = WjGin
    WeightLn(jGout) = WjGout

    call update_weight(Anew, Aold)

    ProbAcc(Order, 2) = ProbAcc(Order, 2) + 1
  else
    kLn(iWLn) = kiWold
  endif
  return
END SUBROUTINE delete_worm_along_wline



!----- move worm along wline : Pupdate(5) --------------
!-----   \2       4/                  \2        4/    ------
!-----Ira ========= jGam     ==>  iGam ========== Ira ------
!-----   /1       3\                  /1        3\    ------
!-----------------------------------------------------------
SUBROUTINE move_worm_along_wline
  implicit none
  integer :: dir, kW, kWold, flag
  integer :: iGam, jGam, WLn, GLn1, GLn2
  integer :: typiGam, typW
  integer :: statiGam, statjGam, statW
  double precision :: tau1, tau2, tau
  double precision :: WWorm, Pacc
  complex*16 :: WiGam, WjGam, WW
  complex*16 :: Anew, Aold, sgn

  !------- step1 : check if worm is present -------------
  if(IsWormPresent .eqv. .false.)  return

  !------- step2 : propose a new config -----------------
  !propose to move Ira from iGam to jGam
  iGam = Ira;  dir = DirecVertex(iGam)
  WLn = NeighVertex(3, Ira)
  jGam = NeighLn(3-dir, WLn)
  if(jGam == Masha .or. jGam==Ira)  return

  GLn1 = NeighVertex(1, iGam);     GLn2 = NeighVertex(2, iGam)

  !------- step3 : configuration check ------------------
  !-------- irreducibility check ------------------------
  kWold = kLn(WLn)
  kLn(WLn) = add_k(kWold, (-1)**dir *kMasha)

  flag = 0
  if(Is_reducible_W(WLn))   flag=1
  if(flag==0) then
    if(Is_reducible_W_Gam(WLn))  flag=1
  endif

  if(flag==1) then
    kLn(WLn) = kWold
    return
  endif

  !-------- the stat, type of the new config ----
  statiGam = delete_ira_stat(StatusVertex(iGam))
  statjGam = add_ira_stat(StatusVertex(jGam))
  statW = line_stat(statiGam, statjGam)

  typW = TypeLn(WLn)
  typiGam = TypeSp2Gam(TypeLn(GLn1), TypeLn(GLn2), SpInVertex(1, iGam), SpInVertex(2, iGam))

  !------- step4 : weight calculation -------------------
  WWorm = weight_worm(GXVertex(jGam)-GXVertex(Masha), GYVertex(jGam)-GYVertex(Masha),  &
    & WXVertex(jGam)-WXVertex(Masha), WYVertex(jGam)-WYVertex(Masha),T3Vertex(jGam)-T3Vertex(Masha))

  tau1 = T3Vertex(iGam)-T2Vertex(iGam)
  tau2 = T1Vertex(iGam)-T3Vertex(iGam)
  WiGam = weight_vertex(statiGam, IsDeltaVertex(iGam), GXVertex(iGam)-WXVertex(iGam), &
    GYVertex(iGam)-WYVertex(iGam), tau1, tau2, typiGam)

  tau1 = T3Vertex(jGam)-T2Vertex(jGam)
  tau2 = T1Vertex(jGam)-T3Vertex(jGam)
  WjGam = weight_vertex(statjGam, IsDeltaVertex(jGam), GXVertex(jGam)-WXVertex(jGam), &
    GYVertex(jGam)-WYVertex(jGam), tau1, tau2, TypeVertex(jGam))

  tau = (-1.d0)**dir *(T3Vertex(iGam)-T3Vertex(jGam))
  WW = weight_line(statW, IsDeltaLn(WLn), 2, WXVertex(iGam)-WXVertex(jGam), WYVertex(iGam) &
    & -WYVertex(jGam), tau, typW)

  Anew = WiGam *WjGam *WW 
  Aold = WeightVertex(iGam) *WeightVertex(jGam) *WeightLn(WLn) 

  call weight_ratio(Pacc, sgn, Anew, Aold)

  Pacc = Pacc *WWorm/WeightWorm

  !------- step5 : accept the update --------------------
  ProbProp(Order, iupdate) = ProbProp(Order, iupdate) + 1
  if(rn()<=Pacc) then

    !------ update diagram info ---------------
    Phase = Phase *sgn

    !------ update Ira ------------------------
    Ira = jGam

    !------ update k --------------------------
    call update_Hash4W(kWold, kLn(WLn)) 

    !------ update type of elements -----------
    TypeVertex(iGam) = typiGam
    TypeLn(WLn) = typW

    !------ update status of elements ---------
    StatusVertex(iGam) = statiGam
    StatusVertex(jGam) = statjGam

    StatusLn(WLn) = statW

    !------ update weight of elements ---------
    WeightWorm = WWorm
    WeightLn(WLn) = WW
    WeightVertex(iGam) = WiGam
    WeightVertex(jGam) = WjGam

    call update_weight(Anew, Aold)

    ProbAcc(Order, 5) = ProbAcc(Order, 5) + 1
  else
    kLn(WLn) = kWold
  endif
  return
END SUBROUTINE move_worm_along_wline



!----- move worm along gline : Pupdate(6) --------------
!---------------------- dir = 1 ----------------------------
!-----      Ira   jGam              iGam    Ira
!----- --<--||--<--||--<--  =>  --<--||--<--||--<--
!-----      ||     ||                ||     ||  
!-----      ||     ||                ||     ||
!-----------------------------------------------------------
!---------------------- dir = 2 ----------------------------
!-----      Ira   jGam              iGam    Ira
!----- -->--||-->--||-->--  =>  -->--||-->--||-->--
!-----      ||     ||                ||     ||  
!-----      ||     ||                ||     ||
!-----------------------------------------------------------
SUBROUTINE move_worm_along_gline
  implicit none
  integer :: dir, kGold, flag, dx, dy
  integer :: iGam, jGam, iW, jW, GLn
  integer :: vertex1, vertex2
  integer :: typG, typiGam, typjGam, typiW
  integer, dimension(2) :: spiGam, spjGam
  integer :: statiGam, statjGam, statiW, statjW, stat1, stat2
  double precision :: tau1, tau2, tau
  double precision :: WWorm, Pacc
  complex*16 :: WiW, WjW, WiGam, WjGam, WG, Anew, Aold, sgn

  !------- step1 : check if worm is present -------------
  if(IsWormPresent .eqv. .false.)   return
  !ProbProp(iupdate) = ProbProp(iupdate) + 1

  !------- step2 : propose a new config -----------------
  iGam = Ira;  iW = NeighVertex(3, iGam)
  dir = Floor(rn()*2.d0) + 1
  GLn = NeighVertex(dir, iGam)
  jGam = NeighLn(dir, GLn);    jW = NeighVertex(3, jGam)
  if(jGam==Masha .or. jGam==Ira)    return

  if(dir==1) then
    if(SpinMasha+TypeLn(GLn)==0 .or. SpinMasha+TypeLn(GLn)==3)  return
  else 
    if(SpinMasha+TypeLn(GLn)==-1 .or. SpinMasha+TypeLn(GLn)==4) return
  endif
  
  kGold = kLn(GLn)
  kLn(GLn) = add_k(kLn(GLn), (-1)**(dir+1) *kMasha)

  !------- step3 : configuration check ------------------
  !-------- irreducibility  check -----------------------
  flag = 0
  if(Is_reducible_G(GLn))          flag = 1
  if(flag==0) then
    if(Is_reducible_G_Gam(GLn))    flag = 1
  endif

  if(flag == 1) then
    kLn(GLn) = kGold
    return
  endif

  !----- type and status of the new config ---
  statiGam = delete_ira_stat(StatusVertex(iGam))
  statjGam = add_ira_stat(StatusVertex(jGam))

  stat1 = StatusVertex(NeighLn(3-DirecVertex(iGam), iW))
  statiW = line_stat(statiGam,  stat1)

  stat2 = StatusVertex(NeighLn(3-DirecVertex(jGam), jW))
  statjW = line_stat(statjGam,  stat2)

  if(iW == jW) then
    statiW = line_stat(statiGam, statjGam)
    statjW = statiW
  endif

  typG = 3 - TypeLn(GLn)

  if(SpInVertex(dir, iGam)==TypeLn(GLn)) then
    !spiGam(dir) = 3- SpInVertex(dir, iGam)
    spiGam(dir) = typG
    spiGam(3-dir) = SpInVertex(3-dir, iGam)
  else
    spiGam(dir) = SpInVertex(dir, iGam)
    spiGam(3-dir) = 3 - SpInVertex(3-dir, iGam)
  endif

  if(SpInVertex(3-dir, jGam)==TypeLn(GLn)) then
    !spjGam(3-dir) = 3- SpInVertex(3-dir, jGam)
    spjGam(3-dir) = typG
    spjGam(dir) = SpInVertex(dir, jGam)
  else
    spjGam(3-dir) = SpInVertex(3-dir, jGam)
    spjGam(dir) = 3 - SpInVertex(dir, jGam)
  endif

  if(dir == 1) then
    typiGam = TypeSp2Gam(typG, TypeLn(NeighVertex(2, iGam)), spiGam(1), spiGam(2)) 
    typjGam = TypeSp2Gam(TypeLn(NeighVertex(1, jGam)), typG, spjGam(1), spjGam(2))
  else
    typiGam = TypeSp2Gam(TypeLn(NeighVertex(1, iGam)), typG, spiGam(1), spiGam(2)) 
    typjGam = TypeSp2Gam(typG, TypeLn(NeighVertex(2, jGam)), spjGam(1), spjGam(2))
  endif

  if(statiW==2 .or. statiW==3) then
    typiW = TypeLn(iW)
  else
    if(DirecVertex(iGam)==1) then
      typiW = TypeGam2W(typiGam, TypeVertex(NeighLn(2, iW)))
    else 
      typiW = TypeGam2W(TypeVertex(NeighLn(1, iW)), typiGam)
    endif
  endif

  !------- step4 : weight calculation -------------------
  WWorm = weight_worm(GXVertex(jGam)-GXVertex(Masha), GYVertex(jGam)-GYVertex(Masha),  &
    & WXVertex(jGam)-WXVertex(Masha), WYVertex(jGam)-WYVertex(Masha),T3Vertex(jGam)-T3Vertex(Masha))

  tau1 = T3Vertex(iGam)-T2Vertex(iGam)
  tau2 = T1Vertex(iGam)-T1Vertex(iGam)
  WiGam = weight_vertex(statiGam, IsDeltaVertex(iGam),GXVertex(iGam)-WXVertex(iGam),   &
    & GYVertex(iGam)-WYVertex(iGam), tau1, tau2, typiGam)

  tau1 = T3Vertex(jGam)-T2Vertex(jGam)
  tau2 = T1Vertex(jGam)-T1Vertex(jGam)
  WjGam = weight_vertex(statjGam, IsDeltaVertex(jGam),GXVertex(jGam)-WXVertex(jGam),   &
    & GYVertex(jGam)-WYVertex(jGam), tau1, tau2, typjGam)

  vertex1 = NeighLn(1, iW)
  vertex2 = NeighLn(2, iW)
  tau = T3Vertex(vertex2) - T3Vertex(vertex1)
  dx  = WXVertex(vertex2) - WXVertex(vertex1)
  dy  = WYVertex(vertex2) - WYVertex(vertex1)
  WiW = weight_line(statiW, IsDeltaLn(iW), 2, dx, dy, tau, typiW)

  vertex1 = NeighLn(1, jW)
  vertex2 = NeighLn(2, jW)
  tau = T3Vertex(vertex2) - T3Vertex(vertex1)
  dx  = WXVertex(vertex2) - WXVertex(vertex1)
  dy  = WYVertex(vertex2) - WYVertex(vertex1)
  WjW = weight_line(statjW, IsDeltaLn(jW), 2, dx, dy, tau, TypeLn(jW))

  tau = T2Vertex(NeighLn(2, GLn)) - T1Vertex(NeighLn(1, GLn))
  WG = weight_line(StatusLn(GLn), IsDeltaLn(GLn), 1, 0, 0, tau, typG)

  Anew = WiGam *WjGam *WiW *WjW *WG
  Aold = WeightVertex(iGam) *WeightVertex(jGam) *WeightLn(iW) *WeightLn(jW) *WeightLn(GLn)

  if(Anew==0.d0)  return

  if(iW==jW) then
    Anew = Anew/WjW
    Aold = Aold/WeightLn(jW)
  endif

  call weight_ratio(Pacc, sgn, Anew, Aold)

  Pacc = Pacc *WWorm/WeightWorm

  !------- step5 : accept the update --------------------
  ProbProp(Order, iupdate) = ProbProp(Order, iupdate) + 1
  if(rn()<=Pacc) then

    !----- update the diagram info -------------- 
    Phase = Phase *sgn

    !---- update Ira and Masha ------------------
    Ira = jGam

    !---- update omega and k --------------------
    call update_Hash4G(kGold, kLn(GLn))
    
    !----- update type configuration ------------
    TypeLn(GLn) = typG
    TypeLn(iW)  = typiW

    TypeVertex(iGam) = typiGam
    TypeVertex(jGam) = typjGam

    SpInVertex(:, iGam) = spiGam(:)
    SpInVertex(:, jGam) = spjGam(:)

    !---- update the status of elements ---------
    StatusVertex(iGam) = statiGam
    StatusVertex(jGam) = statjGam
    StatusLn(iW) = statiW
    StatusLn(jW) = statjW

    !---- update the weight of elements ---------
    WeightWorm = WWorm

    WeightLn(GLn) = WG
    WeightLn(iW) = WiW
    WeightLn(jW) = WjW
    WeightVertex(iGam) = WiGam
    WeightVertex(jGam) = WjGam

    call update_weight(Anew, Aold)

    ProbAcc(Order, 6) = ProbAcc(Order, 6) + 1
  else
    kLn(GLn) = kGold
  endif
END SUBROUTINE move_worm_along_gline
!-----------------------------------------------------------






!------------- add interaction : Pupdate(7) -----------------
!---------------------- dir = 1, dirW = 1 ------------------------
!------               GamC                     GamA     GamC -----
!------    Ira ----<----            Ira ----<---||---<-----  -----
!------                     ==>                 \/           -----
!------  Masha ----<----          Masha ----<---||---<-----  -----
!------               GamD                     GamB     GamD -----
!-----------------------------------------------------------------
!---------------------- dir = 1, dirW = 2 ------------------------
!------               GamC                     GamA     GamC -----
!------    Ira ----<----            Ira ----<---||---<-----  -----
!------                     ==>                 /\           -----
!------  Masha ----<----          Masha ----<---||---<-----  -----
!------               GamD                     GamB     GamD -----
!-----------------------------------------------------------------
!---------------------- dir = 2, dirW = 1 ------------------------
!------               GamC                     GamA     GamC -----
!------    Ira ---->----            Ira ---->---||--->-----  -----
!------                     ==>                 \/           -----
!------  Masha ---->----          Masha ---->---||--->-----  -----
!------               GamD                     GamB     GamD -----
!-----------------------------------------------------------------
!---------------------- dir = 2, dirW = 2 ------------------------
!------               GamC                     GamA     GamC -----
!------    Ira ---->----            Ira ---->---||--->-----  -----
!------                     ==>                 /\           -----
!------  Masha ---->----          Masha ---->---||--->-----  -----
!------               GamD                     GamB     GamD -----
!-----------------------------------------------------------------
!SUBROUTINE add_interaction
  !implicit none
  !integer :: i, dir, dirW, flag
  !integer :: WAB, GIC, GMD, GIA, GMB
  !integer :: GamA, GamB, GamC, GamD
  !integer :: sp, typGamA, typGamB, typAB
  !integer :: statIA, statAC, statMB, statBD
  !integer :: kIA, kMB, kM, q
  !integer :: dxwA, dywA, dxwB, dywB, xwA, ywA, xwB, ywB
  !double precision :: WWorm, Pacc
  !complex*16 :: WA, WB, WGIA, WGAC, WGMB, WGBD, WIra, WMasha, WWAB, WMeasureGam
  !double precision :: Anew, Aold, sgn

  !return
!END SUBROUTINE add_interaction


SUBROUTINE change_gline_time
    implicit none
    
end SUBROUTINE

SUBROUTINE change_gline_space
    implicit none
    
end SUBROUTINE

SUBROUTINE change_wline_time
    implicit none
    !integer :: iWLn,iGam
    !double precision :: Pacc,WNewTau
    !complex*16 :: WOldWeight,WNewWeight,GamOldWeight,GamNewWeight,Anew,Aold,sgn
    !!------- step1 : check if worm is present -------------
    !if(IsWormPresent .eqv. .true.)    return
    !!ProbProp(iupdate) = ProbProp(iupdate) + 1

    !!------- step2 : propose a new config -----------------
    !iGam=generate_gamma()
    !if(IsDeltaVertex(iGam))return
    !WNewTau=generate_tau()

    
end SUBROUTINE

SUBROUTINE change_wline_space
    implicit none
    integer :: iGam, iWLn, jGam, dxw, dyw, xwi, ywi, xwj, ywj, dir
    double precision :: Pacc,T1,T2,T3,T4,T5,T6,WeightIX,WeightIY,WeightJX,WeightJY
    complex*16  ::  WiGam,WjGam, WW, Anew, Aold, sgn
    !------- step1 : check if worm is present -------------
    if(IsWormPresent .eqv. .true.)    return
    !ProbProp(iupdate) = ProbProp(iupdate) + 1

    !------- step2 : propose a new config -----------------
    iWLn=generate_wline()
    iGam = NeighLn(1, iWLn)
    jGam = NeighLn(2, iWLn)

    call generate_x(WXVertex(iGam),xwi,WeightIX);
    call generate_y(WYVertex(iGam),ywi,WeightIY)

    if(IsDeltaLn(iWLn)==0) then
      call generate_x(WXVertex(jGam),xwj,WeightJX);
      call generate_y(WYVertex(jGam),ywj,WeightJY)
    else
      xwj=xwi
      ywj=ywi
    endif

    !------- step3 : configuration check ------------------

    !------- step4 : weight calculation -------------------
    T1=T1Vertex(iGam);
    T2=T2Vertex(iGam);
    T3=T3Vertex(iGam);
    WiGam = weight_vertex(StatusVertex(iGam),IsDeltaVertex(iGam),GXVertex(iGam)-xwi,GYVertex(iGam)-ywi, &
      & T3-T2, T1-T3, TypeVertex(iGam))
    T4=T1Vertex(jGam);
    T5=T2Vertex(jGam);
    T6=T3Vertex(jGam);
    WjGam = weight_vertex(StatusVertex(jGam),IsDeltaVertex(jGam),GXVertex(jGam)-xwj,GYVertex(jGam)-ywj, &
      & T6-T5, T4-T6, TypeVertex(jGam))
    WW = weight_line(StatusLn(iWLn),IsDeltaLn(iWLn),2, xwi-xwj, ywi-ywj,T3-T6,TypeLn(iWLn))


    Anew = WiGam*WjGam*WW
    Aold = WeightLn(iWLn)*WeightVertex(iGam)*WeightVertex(jGam)

    call weight_ratio(Pacc, sgn, Anew, Aold)

    !------- step5 : accept the update --------------------
    ProbProp(Order, iupdate) = ProbProp(Order, iupdate) + 1
    if(rn()<=Pacc) then

      !------ update the diagram info -------------------
      Phase = Phase *sgn

      !------ update the site of elements --------------
      WXVertex(iGam) = xwi
      WYVertex(iGam) = ywi

      WXVertex(jGam) = xwj
      WYVertex(jGam) = ywj
      !------ update the weight of elements ------------
      WeightLn(iWLn) = WW
      WeightVertex(iGam) = WiGam
      WeightVertex(iGam) = WjGam

      call update_weight(Anew, Aold)

      ProbAcc(Order, 13) = ProbAcc(Order, 13) + 1
    endif
    return
end SUBROUTINE



!------------ change Gamma type: Pupdate(14) ------------
SUBROUTINE change_Gamma_type
  implicit none
  integer :: iGam, iWLn, jGam, dir, typGam, typW
  double precision :: tau1, tau2, tau, Pacc
  complex*16 :: WGam, WW, Anew, Aold, sgn

  !------- step1 : check if worm is present -------------
  if(IsWormPresent .eqv. .true.)    return
  !ProbProp(iupdate) = ProbProp(iupdate) + 1

  !------- step2 : propose a new config -----------------
  iGam = generate_gamma()
  if(IsDeltaVertex(iGam)==1)  return

  dir  = DirecVertex(iGam)
  iWLn = NeighVertex(3, iGam)
  jGam = NeighLn(3-dir, iWLn)

  if(TypeVertex(iGam)==1 .or. TypeVertex(iGam)==2) then
    typGam = TypeVertex(iGam)+2
  else if(TypeVertex(iGam)==3 .or. TypeVertex(iGam)==4) then
    typGam = TypeVertex(iGam)-2
  else 
    return
  endif

  if(dir==1) then
    typW = 5-TypeLn(iWLn)
  else
    if(TypeLn(iWLn)<=2) then
      typW = TypeLn(iWLn)+2
    else if(TypeLn(iWLn)<=4) then
      typW = TypeLn(iWLn)-2
    else
      write(*, *) "change Gamma type error!"
      stop
    endif
  endif

  !------- step4 : weight calculation -------------------
  tau1 = T3Vertex(iGam)-T2Vertex(iGam)
  tau2 = T1Vertex(iGam)-T3Vertex(iGam)
  WGam = weight_vertex(StatusVertex(iGam), IsDeltaVertex(iGam), GXVertex(iGam)-WXVertex(iGam), &
    & GYVertex(iGam)-WYVertex(iGam), tau1, tau2, typGam)

  tau = T3Vertex(NeighLn(2, iWLn))-T3Vertex(NeighLn(1, iWLn))
  WW  = weight_line(StatusLn(iWLn), IsDeltaLn(iWLn), 2, WXVertex(iGam)-WXVertex(jGam), &
    & WYVertex(iGam)-WYVertex(jGam), tau, typW) 

  Anew = WW *WGam
  Aold = WeightLn(iWLn)*WeightVertex(iGam)

  call weight_ratio(Pacc, sgn, Anew, Aold)

  !------- step5 : accept the update --------------------
  ProbProp(Order, iupdate) = ProbProp(Order, iupdate) + 1

  if(rn()<=Pacc) then

    !------ update the diagram info -------------------
    Phase = Phase *sgn

    !------ update the site of elements --------------
    TypeVertex(iGam)    = typGam
    SpInVertex(:, iGam)  = 3-SpInVertex(:, iGam)
    TypeLn(iWLn)     = typW

    !------ update the weight of elements ------------
    WeightLn(iWLn) = WW
    WeightVertex(iGam) = WGam

    call update_weight(Anew, Aold)

    ProbAcc(Order, 14) = ProbAcc(Order, 14) + 1
  endif
  return
END SUBROUTINE change_Gamma_type


!------------ move measuring index: Pupdate(15) ------------
SUBROUTINE move_measuring_index
  implicit none
  integer :: iW, iGin, iGout, iGam
  integer :: jW, jGin, jGout, jGam
  integer :: statiGam, statiGin, statiGout, statiW
  integer :: statjGam, statjGin, statjGout, statjW
  double precision :: tau, tau1, tau2, Pacc
  complex*16 :: WiGam, WiGin, WiGout, WiW
  complex*16 :: WjGam, WjGin, WjGout, WjW, Anew, Aold, sgn

  !------- step1 : check if worm is present -------------
  if(IsWormPresent .eqv. .true.)    return
  !ProbProp(iupdate) = ProbProp(iupdate) + 1

  !------- step2 : propose a new config -----------------
  iGam = generate_gamma()
  if(IsDeltaVertex(iGam)/=1)   return

  iGin  = NeighVertex(1, iGam);      iGout = NeighVertex(2, iGam)
  iW    = NeighVertex(3, iGam)
  if(IsDeltaLn(iW)/=0)         return

  jGam  = MeasureGam
  jGin  = NeighVertex(1, jGam);      jGout = NeighVertex(2, jGam)
  jW    = NeighVertex(3, jGam)
  
  if(iGam==jGam)            return

  !----- the status for the new config ------------------
  statiGam  = add_mea_stat(StatusVertex(iGam))
  statjGam  = delete_mea_stat(StatusVertex(jGam))

  if(jW/=iW) then
    statiW    = line_stat(statiGam, StatusVertex(NeighLn(3-DirecVertex(iGam), iW)))
    statjW    = line_stat(statjGam, StatusVertex(NeighLn(3-DirecVertex(jGam), jW)))
  else
    statiW    = line_stat(statiGam, statjGam)
    statjW    = line_stat(statjGam, statiGam)
  endif

  statiGin    = 1
  statiGout   = 1
  if(jGout/=iGin) then
    statjGout   = 0
  else
    statjGout   = 1
  endif

  if(iGout/=jGin) then
    statjGin    = 0
  else
    statjGin    = 1
  endif

  !------- step4 : weight calculation -------------------
  tau1 = T3Vertex(iGam)-T2Vertex(iGam)
  tau2 = T1Vertex(iGam)-T3Vertex(iGam)
  WiGam = weight_vertex(statiGam, IsDeltaVertex(iGam), GXVertex(iGam)-WXVertex(iGam), &
    & GYVertex(iGam)-WYVertex(iGam), tau1, tau2, TypeVertex(iGam))

  tau = T3Vertex(NeighLn(2, iW)) - T3Vertex(NeighLn(1, iW))
  WiW = weight_line(statiW, IsDeltaLn(iW), 2, WXVertex(NeighLn(2, iW))-WXVertex(NeighLn(1,iW)),&
    & WYVertex(NeighLn(2,iW))-WYVertex(NeighLn(1,iW)), tau, TypeLn(iW))

  tau = T2Vertex(NeighLn(2, iGin)) - T1Vertex(NeighLn(1, iGin))
  WiGin = weight_line(statiGin, IsDeltaLn(iGin), 1, 0, 0, tau, TypeLn(iGin))

  tau = T2Vertex(NeighLn(2, iGout)) - T1Vertex(NeighLn(1, iGout))
  WiGout = weight_line(statiGout, IsDeltaLn(iGout), 1, 0, 0, tau, TypeLn(iGout))

  tau1 = T3Vertex(jGam)-T2Vertex(jGam)
  tau2 = T1Vertex(jGam)-T3Vertex(jGam)
  WjGam = weight_vertex(statjGam, 1, GXVertex(jGam)-WXVertex(jGam), GYVertex(jGam)-&
    & WYVertex(jGam), tau1, tau2, TypeVertex(jGam))

  tau = T3Vertex(NeighLn(2, jW)) - T3Vertex(NeighLn(1, jW))
  WjW = weight_line(statjW, 0, 2, WXVertex(NeighLn(2, jW))-WXVertex(NeighLn(1,jW)),&
    & WYVertex(NeighLn(2,jW))-WYVertex(NeighLn(1,jW)), tau, TypeLn(jW))

  tau = T2Vertex(NeighLn(2, jGin)) - T1Vertex(NeighLn(1, jGin))
  WjGin = weight_line(statjGin, IsDeltaLn(jGin), 1, 0, 0, tau, TypeLn(jGin))

  tau = T2Vertex(NeighLn(2, jGout)) - T1Vertex(NeighLn(1, jGout))
  WjGout = weight_line(statjGout, IsDeltaLn(jGout), 1, 0, 0, tau, TypeLn(jGout))

  Anew = WjGam *WjW *WjGin *WjGout *WiGam *WiW *WiGin *WiGout
  Aold = WeightVertex(iGam) *WeightLn(iW) *WeightLn(iGin) *WeightLn(iGout) &
    & *WeightVertex(jGam) *WeightLn(jW) *WeightLn(jGin) *WeightLn(jGout)

  if(abs(Anew)==0.d0)   return

  if(iGin==jGout) then
    Anew = Anew/WjGout
    Aold = Aold/WeightLn(jGout)
  endif

  if(jGin==iGout) then
    Anew = Anew/WjGin
    Aold = Aold/WeightLn(jGin)
  endif

  if(iW==jW) then
    Anew = Anew/WjW
    Aold = Aold/WeightLn(jW)
  endif

  call weight_ratio(Pacc, sgn, Anew, Aold)

  !------- step5 : accept the update --------------------
  ProbProp(Order, iupdate) = ProbProp(Order, iupdate) + 1

  if(rn()<=Pacc) then

    !-------- update the diagram info ---------------
    Phase = Phase *sgn
    MeasureGam = iGam

    !-------- update the status of elements ---------
    StatusVertex(iGam) = statiGam
    StatusLn(iW)    = statiW
    StatusLn(iGin)  = statiGin
    StatusLn(iGout) = statiGout

    StatusVertex(jGam) = statjGam
    StatusLn(jW)    = statjW
    StatusLn(jGin)  = statjGin
    StatusLn(jGout) = statjGout

    !-------- update the weight of elements ---------
    WeightVertex(iGam) = WiGam
    WeightLn(iW) = WiW
    WeightLn(iGin) = WiGin
    WeightLn(iGout) = WiGout

    WeightVertex(jGam) = WjGam
    WeightLn(jW) = WjW
    WeightLn(jGin) = WjGin
    WeightLn(jGout) = WjGout

    call update_weight(Anew, Aold)

    ProbAcc(Order, 15) = ProbAcc(Order, 15) + 1

  endif
END SUBROUTINE move_measuring_index


!------------ change wline isdelta: Pupdate(18) ------------
SUBROUTINE change_wline_isdelta
  implicit none
  integer :: iWLn, iGam, jGam, jGin, jGout, backforth
  !backforth = 1: delta -> normal; backforth = 0: normal -> delta
  double precision :: t3iGam, t3jGam, t1jGam, t2jGam, tau1, tau2, dtau, Pacc
  complex*16 :: WW, WjGin, WjGout, WjGam, Anew, Aold, sgn

  !------- step1 : check if worm is present -------------
  if(IsWormPresent .eqv. .true.)    return
  !ProbProp(iupdate) = ProbProp(iupdate) + 1

  !------- step2 : propose a new config -----------------
  iWLn = generate_wline()
  if(StatusLn(iWLn)==1 .or. StatusLn(iWLn)==3)  return

  iGam = NeighLn(1, iWLn)
  jGam = NeighLn(2, iWLn)
  jGin = NeighVertex(1, jGam)
  jGout = NeighVertex(2, jGam)

  t3iGam = T3Vertex(iGam)
  t1jGam = T1Vertex(jGam)
  t2jGam = T2Vertex(jGam)
  t3jGam = T3Vertex(jGam)

  backforth = IsDeltaLn(iWLn)

  if(backforth==0) then
    t3jGam = t3iGam
  else 
    if(t3jGam/=t3iGam) then
      write(*, *) "There is a bug in delta-W!"
      call print_config
      stop
    endif
    t3jGam = generate_tau()
  endif

  if(IsDeltaVertex(jGam)==1) then
    t1jGam = t3jGam
    t2jGam = t3jGam
  endif
  
  !------- step4 : weight calculation -------------------
  if(backforth==0) then
    WW = weight_line(StatusLn(iWLn), 1, 2, WXVertex(jGam)-WXVertex(iGam), &
      & WYVertex(jGam)-WYVertex(iGam), 0.d0, TypeLn(iWLn)) 
  else
    dtau = t3jGam - t3iGam
    WW = weight_line(StatusLn(iWLn), 0, 2, WXVertex(jGam)-WXVertex(iGam), &
      & WYVertex(jGam)-WYVertex(iGam), dtau, TypeLn(iWLn)) 
  endif

  tau1 = t3jGam-t2jGam
  tau2 = t1jGam-t3jGam
  WjGam = weight_vertex(StatusVertex(jGam), IsDeltaVertex(jGam), GXVertex(jGam)-WXVertex(jGam), &
    & GYVertex(jGam)-WYVertex(jGam), tau1, tau2, TypeVertex(jGam))

  dtau = t2jGam-T1Vertex(NeighLn(1, jGin))
  WjGin = weight_line(StatusLn(jGin), IsDeltaLn(jGin), 1, 0, 0, dtau, TypeLn(jGin)) 

  dtau = T2Vertex(NeighLn(2, jGout))-t1jGam
  WjGout = weight_line(StatusLn(jGout), IsDeltaLn(jGout), 1, 0, 0, dtau, TypeLn(jGout)) 

  Anew = WjGam *WW *WjGin *WjGout
  Aold = WeightLn(iWLn) *WeightVertex(jGam) *WeightLn(jGin) *WeightLn(jGout)

  call weight_ratio(Pacc, sgn, Anew, Aold)
  if(backforth==0) then
    Pacc = Pacc*prob_tau(T3Vertex(jGam))
  else 
    Pacc = Pacc/prob_tau(t3jGam)
  endif

  !------- step5 : accept the update --------------------
  ProbProp(Order, iupdate) = ProbProp(Order, iupdate) + 1

  if(rn()<=Pacc) then

    !------ update the diagram info -------------------
    Phase = Phase *sgn

    !------- update the weight type of elements -------
    IsDeltaLn(iWLn) = 1-backforth

    !------ update the time of elements ---------------
    T3Vertex(iGam) = t3iGam
    T1Vertex(jGam) = t1jGam
    T2Vertex(jGam) = t2jGam
    T3Vertex(jGam) = t3jGam

    !------ update the weight of elements -------------
    WeightLn(iWLn) = WW
    WeightVertex(jGam) = WjGam
    WeightLn(jGin) = WjGin
    WeightLn(jGout) = WjGout

    call update_weight(Anew, Aold)

    ProbAcc(Order, 18) = ProbAcc(Order, 18) + 1
  endif
  return
END SUBROUTINE change_wline_isdelta
  


!------------ change gamma isdelta: Pupdate(19) ------------
SUBROUTINE change_gamma_isdelta
  implicit none
  integer :: iGam, iGin, iGout, iW, backforth
  !backforth = 1: delta -> normal; backforth = 0: normal -> delta
  double precision :: t1iGam, t2iGam, t3iGam, tau1, tau2, dtau, Pacc
  complex*16 :: WiGin, WiGout, WiGam, Anew, Aold, sgn

  !------- step1 : check if worm is present -------------
  if(IsWormPresent .eqv. .true.)    return
  !ProbProp(iupdate) = ProbProp(iupdate) + 1

  !------- step2 : propose a new config -----------------
  iGam = generate_gamma()
  if(iGam==MeasureGam)  return

  iGin = NeighVertex(1, iGam)
  iGout = NeighVertex(2, iGam)

  t1iGam = T1Vertex(iGam)
  t2iGam = T2Vertex(iGam)
  t3iGam = T3Vertex(iGam)

  backforth = IsDeltaVertex(iGam)

  if(backforth==0) then
    t1iGam = t3iGam
    t2iGam = t3iGam
  else 
    if(t1iGam/=t3iGam) then
      write(*, *) "There is a bug in delta-Gam!"
      call print_config
      stop
    endif

    if(t2iGam/=t3iGam) then
      write(*, *) "There is a bug in delta-Gam!"
      call print_config
      stop
    endif

    t1iGam = generate_tau()
    t2iGam = generate_tau()
  endif

  !------- step4 : weight calculation -------------------
  tau1 = t3iGam-t2iGam
  tau2 = t1iGam-t3iGam
  WiGam = weight_vertex(StatusVertex(iGam), 1-backforth, GXVertex(iGam)-WXVertex(iGam), &
    & GYVertex(iGam)-WYVertex(iGam), tau1, tau2, TypeVertex(iGam))

  dtau = t2iGam-T1Vertex(NeighLn(1, iGin))
  WiGin = weight_line(StatusLn(iGin), IsDeltaLn(iGin), 1, 0, 0, dtau, TypeLn(iGin)) 

  dtau = T2Vertex(NeighLn(2, iGout))-t1iGam
  WiGout = weight_line(StatusLn(iGout), IsDeltaLn(iGout), 1, 0, 0, dtau, TypeLn(iGout)) 

  Anew = WiGam *WiGin *WiGout
  Aold = WeightVertex(iGam) *WeightLn(iGin) *WeightLn(iGout)

  call weight_ratio(Pacc, sgn, Anew, Aold)

  if(backforth==0) then
    Pacc = Pacc*prob_tau(T1Vertex(iGam))*prob_tau(T2Vertex(iGam))
  else 
    Pacc = Pacc/(prob_tau(t1iGam)*prob_tau(t2iGam))
  endif

  !------- step5 : accept the update --------------------
  ProbProp(Order, iupdate) = ProbProp(Order, iupdate) + 1

  if(rn()<=Pacc) then

    !------ update the diagram info -------------------
    Phase = Phase *sgn

    !------- update the weight type of elements -------
    IsDeltaVertex(iGam) = 1-backforth

    !------ update the time of elements ---------------
    T1Vertex(iGam) = t1iGam
    T2Vertex(iGam) = t2iGam
    T3Vertex(iGam) = t3iGam

    !------ update the weight of elements -------------
    WeightVertex(iGam) = WiGam
    WeightLn(iGin) = WiGin
    WeightLn(iGout) = WiGout

    call update_weight(Anew, Aold)

    ProbAcc(Order, 19) = ProbAcc(Order, 19) + 1
  endif
  return
END SUBROUTINE change_Gamma_isdelta
  





!--------- exchange the location of Ira and Masha  -----
SUBROUTINE switch_ira_and_masha
  implicit none
  integer  :: k
  if(IsWormPresent .eqv. .false.)  return
  k=Ira;   Ira=Masha; Masha=k
  kMasha=-kMasha;     SpinMasha=-SpinMasha
END SUBROUTINE switch_ira_and_masha

!====================================================================
!====================================================================
!====================================================================


!====================================================================
!============================== WEIGHTS =============================
!====================================================================

!-------- the weight of a line -------------------------
! dx = x2-x1;  dy = y2- y1; tau = tau2-tau1
COMPLEX*16 FUNCTION weight_line(stat, isdelta, knd, dx0, dy0, tau, typ)
  implicit none
  integer :: stat, isdelta, knd, dx, dy, typ
  integer :: dx0, dy0
  double precision :: tau
  integer :: t

  t = Floor(tau*MxT/Beta)

  dx = diff_x(dx0)
  dy = diff_y(dy0)

  !if(stat == 0) then
    !if(knd==1) weight_line = weight_G(typ, t)
    !if(knd==2 .and. isdelta==0) weight_line = weight_W(typ, dx, dy, t)
    !if(knd==2 .and. isdelta==1) weight_line = weight_W0(typ, dx, dy)
  !else if(stat == 1 .or. stat==3) then
    !!  Have measuring vertex around
    !if(knd==1) weight_line = weight_meas_G(typ, t)
    !if(knd==2) weight_line = weight_meas_W(typ, dx, dy, t)
  !else if(stat == 2) then
    !! Have Ira or Masha around (no influence on G)
    !if(knd==2) weight_line = weight_worm_W(typ, dx, dy, t)
    !if(knd==1) then
      !write(*, *) IsWormPresent, iupdate, "gline status == 2 or 3! Error!" 
      !call print_config
      !stop
    !endif
  !else if(stat==-1) then
    !write(*, *) IsWormPresent, iupdate, "line status == -1! There is no weight!" 
    !call print_config
    !stop
  !else
    !write(*, *) IsWormPresent, iupdate, "line status error!", stat
    !call print_config
    !stop
  !endif

  !---------------------- for test --------------------------------------
  if(stat==0 .or. stat==1 .or. stat==3) then
    if(knd==1) weight_line = weight_meas_G(typ, t)
    if(knd==2) weight_line = weight_meas_W(typ, dx, dy, t)
  else if(stat==2) then
    if(knd==2) weight_line = weight_worm_W(typ, dx, dy, t)
    if(knd==1) then
      write(*, *) IsWormPresent, iupdate, "gline status == 2 or 3! Error!" 
      call print_config
      stop
    endif
  else if(stat==-1) then
    write(*, *) IsWormPresent, iupdate, "line status == -1! There is no weight!" 
    stop
  else
    write(*, *) IsWormPresent, iupdate, "line status error!", stat
    stop
  endif
  !------------------------ end -----------------------------------------

  return
END FUNCTION weight_line

!-------- the weight of a vertex -------------------------
!dx = xg-xw;  dy = yg-yw; dtau1 = tau3-tau2; dtau2 = tau1-tau3
COMPLEX*16 FUNCTION weight_vertex(stat, isdelta, dx0, dy0, dtau1, dtau2, typ)
  implicit none
  integer :: stat, dx, dy, t1, t2, typ, isdelta
  integer :: dx0, dy0
  double precision :: weight
  double precision :: dtau1, dtau2

  t1 = Floor(dtau1*MxT/Beta)
  t2 = Floor(dtau2*MxT/Beta)

  dx = diff_x(dx0)
  dy = diff_y(dy0)

  !if(stat==0) then
      !!if(isdelta==0) weight_vertex = weight_Gam(typ, dx, dy, t1, t2)
      !!if(isdelta==1) weight_vertex = weight_Gam0(typ, dx, dy)
      !!----------------- for bare Gamma ------------------------------
      !if(isdelta==0) weight_vertex = 0.d0
      !if(isdelta==1) weight_vertex = weight_Gam0(typ, dx, dy)
      !!------------------------ end ----------------------------------
  !else if(stat==2) then
    !weight_vertex = weight_worm_Gam(typ, dx, dy, t1, t2)
  !else if(stat==1 .or. stat==3) then
    !weight_vertex = weight_meas_Gam(typ, dx, dy, t1, t2)
  !else if(stat==-1) then 
    !write(*, *) IsWormPresent, iupdate, "vertex status == -1! There is no weight!" 
    !stop
  !else
    !write(*, *) IsWormPresent, iupdate, "vertex status error!", stat
    !stop
  !endif

  !---------------------- for test --------------------------------------
  if(stat==0 .or. stat==1 .or. stat==3) then
    weight_vertex = weight_meas_Gam(typ, dx, dy, t1, t2)
  else if(stat==2) then
    weight_vertex = weight_worm_Gam(1, dx, dy, t1, t2)
  else if(stat==-1) then
    write(*, *) IsWormPresent, iupdate, "vertex status == -1! There is no weight!" 
    stop
  else
    write(*, *) IsWormPresent, iupdate, "vertex status error!", stat
    stop
  endif
  !------------------------ end -----------------------------------------
  return
END FUNCTION weight_vertex



!-------- the weight ratio of new/old config -------------------------
SUBROUTINE weight_ratio(Pacc, sgn, Anew, Aold)
  implicit none
  complex*16,intent(in) ::Anew, Aold
  double precision,intent(out) :: Pacc
  complex*16, intent(out) :: sgn
  Pacc = abs(Anew/Aold)
  sgn = (Anew/Aold)/Pacc
  return
END SUBROUTINE weight_ratio

SUBROUTINE update_weight(Anew, Aold)
  implicit none 
  complex*16 :: Anew, Aold

  WeightCurrent = WeightCurrent *abs(Anew/Aold)
  return
END SUBROUTINE update_weight
!====================================================================
!====================================================================
!====================================================================




!====================================================================
!===================== PRINT CONFIGURATION ==========================
!====================================================================


SUBROUTINE print_config
  implicit none
  integer :: i, iln, iv
  
  open(8, access='append', file=trim(title1)//"_mc.conf")
  
  write(8, *) "============================================================"
  write(8, *) imc, IsWormPresent, iupdate

  if(IsWormPresent .eqv. .true.) then
    write(8, *) "Ira", Ira, "Masha", Masha, "SpinMasha", SpinMasha
    write(8, *) "kMasha", kMasha
  endif

  write(8, *) "Order", Order
  write(8, *) "SignFermiLoop", SignFermiLoop

  write(8, *) "Measuring Gamma", MeasureGam
  write(8, *) "Phase", Phase
  write(8, *) "Weight", WeightCurrent

  do i = 1, NGLn
    iln = GLnKey2Value(i)
    if(StatusLn(iln) <0) cycle
    write(8, 10) iln, KindLn(iln), IsDeltaLn(iln), TypeLn(iln), kLn(iln), StatusLn(iln), NeighLn(1:2,iln)
  enddo

  do i = 1, NWLn
    iln = WLnKey2Value(i)
    if(StatusLn(iln) <0) cycle
    write(8, 10) iln, KindLn(iln), IsDeltaLn(iln), TypeLn(iln), kLn(iln), StatusLn(iln), NeighLn(1:2,iln)
  enddo

  do i = 1, NVertex
    iv = VertexKey2Value(i)
    if(StatusVertex(iv) <0) cycle
    write(8, 12) iv,IsDeltaVertex(iv), TypeVertex(iv),SpInVertex(1, iv),SpInVertex(2, iv), &
      & GXVertex(iv),GYVertex(iv),WXVertex(iv),WYVertex(iv), T1Vertex(iv), T2Vertex(iv),  &
      & T3Vertex(iv), DirecVertex(iv), StatusVertex(iv), NeighVertex(:,iv)
  enddo
  write(8, *) "============================================================"

  10 format(' Line:',i2,2x,'kind:',i2,2x,'isdelta:',i2,2x,'type:',i2,2x,'k:',i8,2x,'stat:',i2, 2x,&
    & 'neigh:',i6,i6)
  12 format('Gamma:',i2,2x,'isdelta:',i2,2x,'type:',i2,2x,'typein:',i2,2x,'typeout:',i2,2x,&
    & 'gr:(',i4,i4,'), wr:(',i4,i4,')', 't:(', f7.4, f7.4, f7.4, ')',2x, &
    & 'direction:', i2,2x, 'stat:',i2, 2x,'neigh:', i6,i6,i6)

  close(8)
END SUBROUTINE print_config
!====================================================================
!====================================================================




!!=======================================================================
!!============================= MEASURE =================================
!!=======================================================================
SUBROUTINE measure
  implicit none
  integer :: i, iln,iGam,jGam
  integer :: flag, it
  integer :: spg, spw
  integer :: ityp, nloop
  integer :: MeaGin, MeaGout, MeaW, xg, yg, xw, yw, dir, typ
  integer :: dx, dy, dt1, dt2
  double precision :: factorM
  double precision :: tau1, tau2, tau3

  imeasure = imeasure + 1

  factorM = 1.d0

  !-------- find out the variables for Gamma ----------------
  MeaGin = NeighVertex(2, MeasureGam)
  MeaGout = NeighVertex(1, MeasureGam)
  MeaW = NeighVertex(3, MeasureGam)

  dir = DirecVertex(MeasureGam)

  !----- find the type for Gamma ----------------------
  if(dir==1) then
    if(TypeLn(MeaW)==1 .or. TypeLn(MeaW)==4) then
      spw = 1
    else if(TypeLn(MeaW)==2 .or. TypeLn(MeaW)==3) then
      spw = 2
    endif
  else
    if(TypeLn(MeaW)==1 .or. TypeLn(MeaW)==3) then
      spw = 1
    else if(TypeLn(MeaW)==2 .or. TypeLn(MeaW)==4) then
      spw = 2
    endif
  endif

  if(TypeVertex(MeasureGam)==5 .or. TypeVertex(MeasureGam)==6) then
    typ = 11-TypeVertex(MeasureGam)
  else if(TypeVertex(MeasureGam)==1 .or. TypeVertex(MeasureGam)==3) then
    typ = TypeSp2Gam(1,1,spw,spw)
  else if(TypeVertex(MeasureGam)==2 .or. TypeVertex(MeasureGam)==4) then
    typ = TypeSp2Gam(2,2,spw,spw)
  endif
  ityp = (typ+1)/2

  nloop = Mod(Floor(SignFermiloop/2.d0+0.5d0), 2)


  !----- find the space and time variables for Gamma -------
  xg = GXVertex(MeasureGam)
  yg = GYVertex(MeasureGam)
  xw = WXVertex(NeighLn(3-dir, MeaW))
  yw = WYVertex(NeighLn(3-dir, MeaW))
  dx = diff_x(xg-xw)
  dy = diff_y(yg-yw)

  tau1 = T1Vertex(NeighLn(1, MeaGin))
  tau2 = T2Vertex(NeighLn(2, MeaGout))
  tau3 = T3Vertex(NeighLn(3-dir, MeaW))

  dt1 = Floor((tau3-tau2)*MxT/Beta)
  if(dt1<0) then
    dt1 = dt1 + MxT
    factorM = factorM * (-1.d0)
  endif

  dt2 = Floor((tau1-tau3)*MxT/Beta)
  if(dt2<0) then
    dt2 = dt2 + MxT
    factorM = factorM * (-1.d0)
  endif
  
  factorM = factorM *CoefOfSymmetry(dx, dy)* CoefOfWeight(Order)*WeightLn(MeaW) &
    & *WeightVertex(MeasureGam)* WeightLn(MeaGin)*WeightLn(MeaGout)

  !------------------- accumulation -------------------------------------------------------
  GamMC(Order, nloop, ityp, dx, dy, dt1, dt2) = GamMC(Order, nloop, ityp, dx, dy, dt1, dt2) &
    & + Phase/factorM
  GamSqMC(Order,nloop, ityp, dx, dy, dt1, dt2 ) = GamSqMC(Order,nloop, ityp, dx, dy, dt1, dt2) &
    & + (Phase/factorM)**2.d0

  if(Order==0) then
    GamNorm = GamNorm + Phase
  endif

  !===============  test variables =================================
  !iGam=NeighLn(3,1)
  !if(IsDeltaLn(3)==0 .and. StatusLn(3)==0) TestData(1) = TestData(1) +1.d0
  if(IsWormPresent) then
    if(NeighVertex(3, Ira)==NeighVertex(3, Masha)) then
      TestData(1) = TestData(1) +1.d0
    else
      TestData(2) = TestData(2) +1.d0
    endif
  else
    if(TypeLn(3)==1 .and. TypeLn(6)==1) TestData(3) = TestData(3) + 1.d0
    TestData(0)=TestData(0)+1.d0
  endif
  !================================================================
  
        
END SUBROUTINE measure
!!=======================================================================
!!=======================================================================
!!=======================================================================


!====================================================================
