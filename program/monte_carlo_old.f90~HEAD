
!!====================== initialization =================================
!!=======================================================================

!------------- definition of the config of diagram ----------------
subroutine def_conf
  implicit none
  integer :: i

  CheckG = .true.
  CheckW = .true.
  CheckGamma = .false.

  call def_prob

  !--------------- initialization for weight ---------------
  Ln4GList(:) = 0
  Ln4WList(:) = 0
  Vertex4GamList(:) = 0
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

  call def_spin
  call def_diagram

  return
end subroutine def_conf


SUBROUTINE def_spin
  implicit none

  TypeGam2W(:,:) = 0
  TypeGW2Gam(:,:,:,:) = 0

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
  !---- TypeGW2Gam(Type(Gin), Type(Gout), Type(Gamin), Type(Gamout)) = Type(Gam)
  TypeGW2Gam(1,1,1,1) = 1
  TypeGW2Gam(2,2,2,2) = 2
  TypeGW2Gam(1,1,2,2) = 3
  TypeGW2Gam(2,2,1,1) = 4
  TypeGW2Gam(1,2,1,2) = 5
  TypeGW2Gam(2,1,2,1) = 6

  return
END SUBROUTINE def_spin


SUBROUTINE def_diagram
  implicit none

  !-------------- 1-order diagram ------------------------
  Order = 1
  ! the index of measuring gamma
  NGLn = 4;  NWLn = 2;  NGam = 4
  ! the number of glines, wlines, gamma
  MeasGamma = 1
  ! number of fermi loops
  NFermiLoop = 0
  ! the phase of the diagram
  Phase = 1.d0
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

  Ln4GList(1)= 1;     Ln4GList(2)= 2;     Ln4GList(3)= 4
  Ln4GList(4)= 5
  Ln4WList(1)= 3;     Ln4WList(2)= 6

  List4Ln(1)  = 1;     List4Ln(2)  = 2;     List4Ln(3)  = 1
  List4Ln(4)  = 3;     List4Ln(5)  = 4;     List4Ln(6)  = 2


  !kind of lines: 1: gline;  2: wline
  KindLn(1) = 1;     KindLn(2) = 1;     KindLn(3) = 2
  KindLn(4) = 1;     KindLn(5) = 1;     KindLn(6) = 2

  !type of glines: 1: spin up; 2: spin down
  TypeLn(1) = 1;     TypeLn(2) = 1
  TypeLn(4) = 1;     TypeLn(5) = 1

  !type of gamma inside spins: 1: spin up; 2: spin down
  TypeVertexIn(:) = 1;   TypeVertexOut(:) = 1

  !type of Gamma: 1: gin, 2: gout, 3: win, 4: wout
  TypeVertex(1) = 1
  TypeVertex(2) = 1
  TypeVertex(3) = 1
  TypeVertex(4) = 1

  !type of wlines
  TypeLn(3) = TypeGam2W(TypeVertex(1), TypeVertex(2))
  TypeLn(6) = TypeGam2W(TypeVertex(3), TypeVertex(4))

  ! the omega on lines
  OmegaLn(1) = 0
  OmegaLn(2) = 0
  OmegaLn(3) = OmegaLn(2)-OmegaLn(1)
  OmegaLn(4) = 0
  OmegaLn(5) = OmegaLn(3)+OmegaLn(4)
  OmegaLn(6) = OmegaLn(1)-OmegaLn(4)


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
    Vertex4GamList(i) = i
    List4Vertex(i)  = i
  enddo
  TailGam = 5

  ! the site variables of Gamma 1 and 2
  GXVertex(1) = 0;            GYVertex(1) = 0
  GXVertex(2) = GXVertex(1);     GYVertex(2) = GYVertex(1)
  GXVertex(3) = GXVertex(1);     GYVertex(3) = GYVertex(1)
  GXVertex(4) = GXVertex(1);     GYVertex(4) = GYVertex(1)
  WXVertex(1) = 0;            WYVertex(1) = 0
  WXVertex(2) = 0;            WYVertex(2) = 0
  WXVertex(3) = 0;            WYVertex(3) = 0
  WXVertex(4) = 0;            WYVertex(4) = 0

  ! Direction of Gamma: 1: left of W;  2: right of W
  DirecVertex(1) = 1;                DirecVertex(2) = 2
  DirecVertex(3) = 1;                DirecVertex(4) = 2

  ! NeighVertex(i, j): the ith neighbor line of gamma j
  NeighVertex(1,1) = 2;        NeighVertex(2,1) = 1;        NeighVertex(3,1) = 3
  NeighVertex(1,2) = 4;        NeighVertex(2,2) = 5;        NeighVertex(3,2) = 3
  NeighVertex(1,3) = 1;        NeighVertex(2,3) = 4;        NeighVertex(3,3) = 6
  NeighVertex(1,4) = 5;        NeighVertex(2,4) = 2;        NeighVertex(3,4) = 6

  ! weights for lines and vertexes
  WeightLn(1) = weight_line(StatusLn(1),1,0,0,OmegaLn(1),TypeLn(1))
  WeightLn(2) = weight_line(StatusLn(2),1,0,0,OmegaLn(2),TypeLn(2))
  WeightLn(3) = weight_line(StatusLn(3),2,0,0,OmegaLn(3),TypeLn(3))
  WeightLn(4) = weight_line(StatusLn(4),1,0,0,OmegaLn(4),TypeLn(4))
  WeightLn(5) = weight_line(StatusLn(5),1,0,0,OmegaLn(5),TypeLn(5))
  WeightLn(6) = weight_line(StatusLn(6),2,0,0,OmegaLn(6),TypeLn(6))

  WeightVertex(1) = weight_vertex(1, StatusVertex(1), 0, 0, OmegaLn(2), OmegaLn(1), TypeVertex(1))
  WeightVertex(2) = weight_vertex(1, StatusVertex(2), 0, 0, OmegaLn(4), OmegaLn(5), TypeVertex(2))
  WeightVertex(3) = weight_vertex(1, StatusVertex(3), 0, 0, OmegaLn(1), OmegaLn(4), TypeVertex(3))
  WeightVertex(4) = weight_vertex(1, StatusVertex(4), 0, 0, OmegaLn(5), OmegaLn(2), TypeVertex(4))


  WeightCurrent = CoefOfWeight(1)*WeightLn(1)*WeightLn(2)*WeightLn(3)*WeightLn(4)* &
    & WeightLn(5)*WeightLn(6)*WeightVertex(1)*WeightVertex(2)*WeightVertex(3)* &
    & WeightVertex(4)*(1.d0/Beta)**Order *(-1.d0)**NFermiLoop

  Phase = Phase *WeightCurrent/abs(WeightCurrent)

  !-------------------------------------------------------

  !call print_config

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
  integer :: iflag, imeasure
  double precision :: nr

  imeasure = 0
  do while(imeasure < nw)
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
      case( 3) 
        call create_worm_along_gline  
      case( 4) 
        call delete_worm_along_gline  
      case( 5) 
        call move_worm_along_wline             
      case( 6) 
        call move_worm_along_gline              
      case( 7) 
        call add_interaction                  
      case( 8) 
        call remove_interaction             
      case( 9) 
        call add_interaction_cross              
      case(10)
        call remove_interaction_cross                  
      case(11) 
        call reconnect                      
      case(12) 
        call shift_gline_in_space          
      case(13)  
        call shift_wline_in_space         
      case(14) 
        call change_Gamma_type     
      case(15) 
        call move_measuring_index       
    end select

    imc = imc + 1


    if(IsWormPresent ) then
      GamWormOrder(Order) = GamWormOrder(Order) + 1

      if(rn()<=0.5d0) call switch_ira_and_masha

      !if(Order==0) then
        !HistoOmegaW(OmegaMasha) = HistoOmegaW(OmegaMasha) + 1
        !HistoOmegaW(OmegaLn(NeighVertex(1, MeasGamma))) = HistoOmegaW(OmegaLn(NeighVertex(1, &
          !& MeasGamma))) + 1
        !HistoOmegaW(OmegaLn(NeighVertex(2, MeasGamma))) = HistoOmegaW(OmegaLn(NeighVertex(2, &
          !& MeasGamma))) + 1
        !HistoOmegaW(OmegaLn(NeighVertex(3, MeasGamma))) = HistoOmegaW(OmegaLn(NeighVertex(3, &
          !& MeasGamma))) + 1
        !if(NeighVertex(1, Ira)==NeighVertex(2, Ira) .or. NeighVertex(1, Masha) &
          !& ==NeighVertex(2, Masha)) then
          !Gam0Bubble = Gam0Bubble +1
          !if(Ira==MeasGamma) then
            !HistoOmegaW(OmegaLn(NeighVertex(1, Masha))) = HistoOmegaW(OmegaLn(NeighVertex(1, &
              !& Masha))) + 1
          !else if(Masha==MeasGamma) then
            !HistoOmegaW(OmegaLn(NeighVertex(1, Ira))) = HistoOmegaW(OmegaLn(NeighVertex(1, &
              !& Ira))) + 1
          !endif
        !endif
      !endif

    else
      GamOrder(Order) = GamOrder(Order) + 1
      imeasure = imeasure + 1
    endif
  enddo

END SUBROUTINE markov



!===============================================================
!========================= updates =============================
!===============================================================


!---------------- create worm : Pupdate(1) -----------------
!------       \iGout    jGout/           \           /     -----
!------   iGam ======>======= jGam  => Ira====>======Masha -----
!------       /iGin      jGin\           /           \     -----
!---------------------------------------------------------------
SUBROUTINE create_worm_along_wline
  implicit none
  integer :: NWLn, iWLn, iGam, jGam, iGin, jGin, iGout, jGout
  integer :: statiGam, statjGam, statW
  integer :: statiGin, statiGout, statjGin, statjGout
  integer :: k, omega, kiW, omegaiW
  integer :: spin
  integer :: tempGam
  double precision :: WiGin, WiGout, WjGin, WjGout
  double precision :: WWorm, WIra, WMasha, WWNew, Anew, Aold, Pacc, sgn
  double precision :: rand
  integer :: kiWold, flag

  !------------ step1 : check if worm is present ------------------
  if(IsWormPresent .eqv. .true.)   return
  !ProbProp(iupdate) = ProbProp(iupdate) + 1

  !------------ step2 : propose the new config ------------------
  iWLn = generate_wline()
  omega = generate_omega()
  omegaiW = omega + OmegaLn(iWLn)
  if(Is_omega_not_valid(OmegaiW))         return

  spin = 4*(Floor(rn()*2.d0))-2

  iGam = NeighLn(1, iWLn);              jGam = NeighLn(2, iWLn)
  iGin = NeighVertex(1, iGam);          jGin = NeighVertex(1, jGam)
  iGout = NeighVertex(2, iGam);         jGout = NeighVertex(2, jGam)
  if(spin+2*(TypeLn(jGout)-TypeLn(jGin))==0)  return

  k = generate_k()
  kiWold = kLn(iWLn)
  kLn(iWLn) = add_k(kiWold, k)

  !------------ step3 : configuration check -------------------
  flag=0
  if(flag==0) then
    if(Is_reducible_W(iWLn))                 flag=1
  endif
  if(flag==0) then
    if(Is_reducible_W_Gamma(iWLn))           flag=1
  endif

  if(flag==1) then
    kLn(iWLn) = kiWold
    return
  endif

  statiGam = add_ira_stat(StatusVertex(iGam))
  statjGam = add_ira_stat(StatusVertex(jGam))
  statW    = line_stat(statiGam, statjGam)

  !------------ step4 : weight calculation ----------------------
  WiGin  = weight_line(StatusLn(iGin),  1, 0, 0, OmegaLn(iGin),  TypeLn(iGin))
  WiGout = weight_line(StatusLn(iGout), 1, 0, 0, OmegaLn(iGout), TypeLn(iGout))
  WjGin  = weight_line(StatusLn(jGin),  1, 0, 0, OmegaLn(jGin),  TypeLn(jGin))
  WjGout = weight_line(StatusLn(jGout), 1, 0, 0, OmegaLn(jGout), TypeLn(jGout))

  WIra = weight_vertex(Order, statiGam, diff_x(GXVertex(iGam),WXVertex(iGam)),  &
    & diff_y(GYVertex(iGam),WYVertex(iGam)), OmegaLn(iGin), OmegaLn(iGout), TypeVertex(iGam))

  WMasha = weight_vertex(Order, statjGam, diff_x(GXVertex(jGam),WXVertex(jGam)), &
    & diff_y(GYVertex(jGam),WYVertex(jGam)), OmegaLn(jGin), OmegaLn(jGout), TypeVertex(jGam))

  WWNew = weight_line(statW, 2, diff_x(WXVertex(jGam),WXVertex(iGam)), &
    & diff_y(WYVertex(jGam),WYVertex(iGam)), omegaiW, TypeLn(iWLn))

  WWorm = weight_worm(diff_x(GXVertex(iGam), GXVertex(jGam)), diff_y(GYVertex(iGam),GYVertex(jGam)),&
    & diff_x(WXVertex(iGam),WXVertex(jGam)), diff_y(WYVertex(iGam),WYVertex(jGam)), omega)

  Anew = WWorm* WIra *WMasha *WWNew *WiGin *WiGout *WjGin *WjGout
  Aold = WeightLn(iWLn) *WeightVertex(iGam) *WeightVertex(jGam) *WeightLn(iGin) &
  & *WeightLn(iGout) *WeightLn(jGin) *WeightLn(jGout)

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
  Pacc = Pacc *(Order+1.d0)/(CoefOfWorm*0.5d0*prob_omega(omega)) 
  Pacc = Pacc *Pupdate(2)/Pupdate(1)

  !------------ step5 : accept the update -----------------------
  ProbProp(Order, iupdate) = ProbProp(Order, iupdate) + 1
  rand = rn()
  !write(*, *) "create_worm", Pacc, rand
  !call print_markov_before(Pacc, rand)
  if(rand<=Pacc) then

    !-------------- update the diagram info --------------------
    Phase = Phase*sgn
    IsWormPresent = .true.

    !---------- update Ira and Masha -------------------
    Ira = iGam
    Masha = jGam
    kMasha = k
    OmegaMasha = omega
    SpinMasha = spin

    !----------- update k and omega of elements --------
    call update_Hash4W(kiWold, kLn(iWLn))
    OmegaLn(iWLn) = omegaiW

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
!!-----------------------------------------------------------------




!!--------------- delete worm along wline: Pupdate(2) ----------------
!------       \iGout    jGout/           \           /     -----
!------   Ira  ======>======= Masha =>iGam====>======jGam  -----
!------       /iGin      jGin\           /           \     -----
!---------------------------------------------------------------
SUBROUTINE delete_worm_along_wline
  implicit none
  integer :: dir
  integer :: NWLn, iWLn, iGam, jGam, iGin, jGin, iGout, jGout
  integer :: typIra, typMasha, typW
  integer :: statIra, statMasha, statW
  integer :: statiGin, statiGout, statjGin, statjGout
  integer :: k, omega, kiW, omegaiW
  double precision :: WiGin, WiGout, WjGin, WjGout
  double precision :: Wi, Wj, WWNew, Anew, Aold, Pacc, sgn
  double precision :: rand
  integer :: kiWold
  integer :: flag

  !------------ step1 : check if worm is present ------------------
  if(IsWormPresent .eqv. .false.)  return
  !ProbProp(iupdate) = ProbProp(iupdate) + 1
  if(NeighVertex(3, Ira)/= NeighVertex(3, Masha))   return

  !------------ step2 : propose the new config ------------------
  iWLn = NeighVertex(3, Ira)
  dir = DirecVertex(Ira)

  omega = (-1)**dir*OmegaMasha
  if(Is_delta_omega_not_valid(omega))  return
  omegaiW = OmegaLn(iWLn) + omega
  if(Is_omega_not_valid(OmegaiW))  return

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
    if(Is_reducible_W_Gamma(iWLn))  flag=1
  endif

  if(flag==1) then
    kLn(iWLn) = kiWold
    return
  endif

  !------------ status, type for the new configuration ------------
  statIra   = delete_ira_stat(StatusVertex(Ira))
  statMasha = delete_ira_stat(StatusVertex(Masha))
  statW     = line_stat(statIra, statMasha)

  typIra = TypeGW2Gam(TypeLn(iGin),TypeLn(iGout),TypeVertexIn(Ira), TypeVertexOut(Ira))
  typMasha = TypeGW2Gam(TypeLn(jGin),TypeLn(jGout),TypeVertexIn(Masha),TypeVertexOut(Masha))
  if(dir==1) then
    typW = TypeGam2W(typIra, typMasha)
  else
    typW = TypeGam2W(typMasha, typIra)
  endif


  !------------ step4 : weight calculation ---------------------
  WiGin  = weight_line(StatusLn(iGin),  1, 0, 0, OmegaLn(iGin),  TypeLn(iGin))
  WiGout = weight_line(StatusLn(iGout), 1, 0, 0, OmegaLn(iGout), TypeLn(iGout))
  WjGin  = weight_line(StatusLn(jGin),  1, 0, 0, OmegaLn(jGin),  TypeLn(jGin))
  WjGout = weight_line(StatusLn(jGout), 1, 0, 0, OmegaLn(jGout), TypeLn(jGout))

  Wi = weight_vertex(Order, statIra, diff_x(GXVertex(Ira),WXVertex(Ira)), diff_y(GYVertex(Ira),&
    & WYVertex(Ira)), OmegaLn(iGin), OmegaLn(iGout), typIra)
  Wj = weight_vertex(Order, statMasha, diff_x(GXVertex(Masha),WXVertex(Masha)),  &
    & diff_y(GYVertex(Masha),WYVertex(Masha)), OmegaLn(jGin), OmegaLn(jGout), typMasha)
  WWNew = weight_line(statW, 2, diff_x(WXVertex(Masha),WXVertex(Ira)), diff_y(WYVertex(Masha), &
    & WYVertex(Ira)), omegaiW, typW)

  Anew = Wi *Wj *WWNew *WiGin *WiGout *WjGin *WjGout
  Aold = WeightWorm *WeightLn(iWLn) *WeightVertex(Ira) *WeightVertex(Masha) *WeightLn(iGin) &
  & *WeightLn(iGout) *WeightLn(jGin) *WeightLn(jGout)

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
  AveWeightRatio(Order, 1) = AveWeightRatio(Order, 1) + Pacc
  ProbProp(Order, iupdate) = ProbProp(Order, iupdate) + 1

  Pacc = Pacc *(CoefOfWorm*0.5d0*prob_omega(omega))/(Order+1.d0) 
  Pacc = Pacc *Pupdate(1)/Pupdate(2)

  !------------ step5 : accept the update -----------------------
  rand = rn()
  !write(*, *) "delete worm", Pacc, rand
  !call print_markov_before(Pacc, rand)
  if(rand<=Pacc) then

    !---------- update the diagram info ----------------
    Phase = Phase*sgn
    IsWormPresent = .false.

    !----------- update k and omega of elements --------
    call update_Hash4W(kiWold, kLn(iWLn))
    OmegaLn(iWLn) = omegaiW

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
  !call print_markov_after
  return
END SUBROUTINE delete_worm_along_wline
!-----------------------------------------------------------------





!------------- create worm along gline: Pupdate(3) -----------------
!-----     iGam   jGam               Ira   Masha
!----- -->--||-->--||-->--  ==> -->--||-->--||-->--
!-----      ||     ||                ||     ||  
!-----      ||     ||                ||     ||
!-----------------------------------------------------------
SUBROUTINE create_worm_along_gline
  implicit none
  integer :: GLn1, GLn2, GLn3, iGam, jGam, iW, jW
  integer :: statiGam, statjGam, statiW, statjW
  integer :: statGLn1, statGLn2, statGLn3
  integer :: dxiW, dyiW, dxjW, dyjW
  integer :: k, omega, kGLn2, omegaGLn2
  integer :: spin
  integer :: typiGamin, typiGamout, typjGamin, typjGamout
  double precision :: WGLn1, WGLn2, WGLn3, WiW, WjW
  double precision :: WWorm, WiGam, WjGam, Anew, Aold, Pacc, sgn
  double precision :: rand
  integer :: kGLn2old, flag

  !------------ step1 : check if worm is present ------------------
  if(IsWormPresent .eqv. .true.)   return
  !ProbProp(iupdate) = ProbProp(iupdate) + 1

  !------------ step2 : propose the new config ------------------
  GLn2 = generate_gline()

  iGam = NeighLn(1, GLn2);          jGam = NeighLn(2, GLn2)
  GLn1 = NeighVertex(1, iGam);      GLn3 = NeighVertex(2, jGam)
  iW   = NeighVertex(3, iGam);      jW   = NeighVertex(3, jGam)

  omega = generate_omega()
  omegaGLn2 = omega + OmegaLn(GLn2)
  if(Is_omega_not_valid(OmegaGLn2)) return 

  spin = 4*TypeLn(GLn2)-6     !s==1, ds==-2; s==2, ds==2
  k = generate_k()

  !------ update the frequency --------------------------------
  kGLn2old = kLn(GLn2)
  kLn(GLn2) = add_k(kGLn2old, k)

  !------------ step3 : configuration check -------------------
  flag=0
  if(flag==0) then
    if(Is_reducible_G(GLn2))                 flag=1
  endif
  if(flag==0) then
    if(Is_reducible_G_Gamma(GLn2))           flag=1
  endif

  if(flag==1) then
    kLn(GLn2) = kGLn2old
    return
  endif

  !----- site, type and status of the new config ---

  statiGam = add_ira_stat(StatusVertex(iGam))
  statjGam = add_ira_stat(StatusVertex(jGam))

  if(iW/=jW) then
    statiW   = line_stat(statiGam, StatusVertex(NeighLn(3-DirecVertex(iGam), iW)))
    statjW   = line_stat(statjGam, StatusVertex(NeighLn(3-DirecVertex(jGam), jW)))
  else
    statiW = line_stat(statiGam, statjGam)
    statjW = line_stat(statiGam, statjGam)
  endif

  typiGamin  = TypeVertexIn(iGam)
  typiGamout = TypeVertexOut(iGam)
  typjGamin  = TypeVertexIn(jGam)
  typjGamout = TypeVertexOut(jGam)

  if(typiGamin==TypeLn(GLn1)) then
    typiGamout = 3 - typiGamout
  else
    typiGamin = 3- typiGamin
  endif
  
  if(typjGamout==TypeLn(GLn3)) then
    typjGamin = 3 - typjGamin
  else
    typjGamout = 3- typjGamout
  endif

  dxiW  = diff_x(WXVertex(iGam), WXVertex(NeighLn(3-DirecVertex(iGam),iW)))
  dyiW  = diff_y(WYVertex(iGam), WYVertex(NeighLn(3-DirecVertex(iGam),iW)))
 
  dxjW  = diff_x(WXVertex(NeighLn(3-DirecVertex(jGam),jW)), WXVertex(jGam))
  dyjW  = diff_y(WYVertex(NeighLn(3-DirecVertex(jGam),jW)), WYVertex(jGam))

  !------- step4 : weight calculation -------------------
  WiW = weight_line(statiW, 2, dxiW, dyiW, OmegaLn(iW), TypeLn(iW))
  WjW = weight_line(statjW, 2, dxjW, dyjW, OmegaLn(jW), TypeLn(jW))
  WGLn1 = weight_line(StatusLn(GLn1), 1, 0, 0, OmegaLn(GLn1), TypeLn(GLn1))
  WGLn3 = weight_line(StatusLn(GLn3), 1, 0, 0, OmegaLn(GLn3), TypeLn(GLn3))
  WGLn2 = weight_line(StatusLn(GLn2), 1, 0, 0, omegaGLn2, 3-TypeLn(GLn2))

  WiGam = weight_vertex(Order, statiGam, diff_x(GXVertex(iGam),WXVertex(iGam)), diff_y(GYVertex(iGam) &
    & ,WYVertex(iGam)), OmegaLn(GLn1), omegaGLn2, TypeVertex(iGam))
  WjGam = weight_vertex(Order, statjGam, diff_x(GXVertex(jGam),WXVertex(jGam)), diff_y(GYVertex(jGam) &
    & ,WYVertex(jGam)), omegaGLn2, OmegaLn(GLn3), TypeVertex(jGam))

  WWorm = weight_worm(diff_x(GXVertex(iGam),GXVertex(jGam)),diff_y(GYVertex(iGam),GYVertex(jGam)), &
    & diff_x(WXVertex(jGam),WXVertex(iGam)),diff_y(WYVertex(jGam),WYVertex(iGam)), omega)

  Anew  = WWorm *WiW *WjW *WiGam *WjGam *WGLn1 *WGLn2 *WGLn3
  Aold  = WeightLn(iW) *WeightLn(jW) *WeightLn(GLn1) *WeightLn(GLn2) *WeightLn(GLn3) &
    & *WeightVertex(iGam) *WeightVertex(jGam)

  if(GLn1==GLn3) then
    Anew = Anew/WGLn3
    Aold = Aold/WeightLn(GLn3)
  endif

  if(iW==jW) then
    Anew = Anew/WjW
    Aold = Aold/WeightLn(jW)
  endif


  call weight_ratio(Pacc, sgn, Anew, Aold)
  if(Order>0) then
    Pacc = Pacc *(2.d0*Order+2.d0)/(CoefOfWorm*prob_omega(omega)) 
  else
    Pacc = Pacc/(CoefOfWorm*prob_omega(omega))
  endif
  Pacc = Pacc *Pupdate(4)/Pupdate(3)

  !------------ step5 : accept the update -----------------------
  ProbProp(Order, iupdate) = ProbProp(Order, iupdate) + 1
  rand = rn()
  if(rand<=Pacc) then

    !-------------- update the diagram info --------------------
    Phase = Phase*sgn
    IsWormPresent = .true.

    !---------- update Ira and Masha -------------------
    Ira = iGam
    Masha = jGam
    kMasha = k
    OmegaMasha = omega
    SpinMasha = spin

    !----------- update k and omega of elements --------
    call update_Hash4G(kGLn2old, kLn(GLn2))
    OmegaLn(GLn2) = omegaGLn2

    !----------- update status of elements -------------
    StatusVertex(Ira)  = statiGam
    StatusVertex(Masha)= statjGam
    StatusLn(iW)  = statiW
    StatusLn(jW)  = statjW

    !----------- update the type of elements -----------
    TypeLn(GLn2) = 3-TypeLn(GLn2)

    TypeVertexIn(iGam)  = typiGamin
    TypeVertexOut(iGam) = typiGamout
    TypeVertexIn(jGam)  = typjGamin
    TypeVertexOut(jGam) = typjGamout

    !----------- update weight of elements -------------
    WeightWorm = WWorm
    WeightVertex(Ira) = WiGam
    WeightVertex(Masha) = WjGam

    WeightLn(GLn1) = WGLn1
    WeightLn(GLn2) = WGLn2
    WeightLn(GLn3) = WGLn3

    WeightLn(iW)  = WiW
    WeightLn(jW)  = WjW

    call update_weight(Anew, Aold)

    ProbAcc(Order, 3) = ProbAcc(Order, 3) + 1
  else
    kLn(GLn2) = kGLn2old
  endif
  return
END SUBROUTINE create_worm_along_gline
!!-----------------------------------------------------------------




!!---------- delete worm along gline: Pupdate(4) ----------------
!-----    dir==1  ----------------------------------------------
!-----     iGam   jGam               Ira   Masha
!----- --<--||--<--||--<-- <==  --<--||--<--||--<--
!-----      ||     ||                ||     ||  
!-----      ||     ||                ||     ||
!---------------------------------------------------------------
!-----    dir==2  ----------------------------------------------
!-----     iGam   jGam               Ira   Masha
!----- -->--||-->--||-->-- <==  -->--||-->--||-->--
!-----      ||     ||                ||     ||  
!-----      ||     ||                ||     ||
!---------------------------------------------------------------
SUBROUTINE delete_worm_along_gline
  implicit none
  integer :: dir
  integer :: iGam, jGam, GLn1, GLn2, GLn3, iW, jW
  integer :: typIra, typMasha, typiGamin, typiGamout, typjGamin, typjGamout 
  integer :: typiW, typjW, typGLn2
  integer :: dxiW, dyiW, dxjW, dyjW
  integer :: statIra, statMasha, statiW, statjW
  integer :: statGLn1, statGLn2, statGLn3
  integer :: k, omega, kGLn2, omegaGLn2
  double precision :: WGLn1, WGLn2, WGLn3
  double precision :: WiW, WjW, WiGam, WjGam, Anew, Aold, Pacc, sgn
  double precision :: rand
  integer :: kGLn2old
  integer :: flag

  !------------ step1 : check if worm is present ------------------
  if(IsWormPresent .eqv. .false.)  return
  !ProbProp(iupdate) = ProbProp(iupdate) + 1
  if(NeighVertex(1,Ira)== NeighVertex(2,Masha)) then
    dir = 1
  else if (NeighVertex(2,Ira)==NeighVertex(1,Masha)) then
    dir = 2
  else
    return
  endif

  !------------ step2 : propose the new config ------------------
  GLn2  = NeighVertex(dir, Ira)
  if((-1)**dir*SpinMasha/= 6-4*TypeLn(GLn2)) return     !s==1, (-1)**dir*ds==2; s==2, (-1)**dir*ds==-2

  GLn1  = NeighVertex(3-dir, Ira);      GLn3 = NeighVertex(dir, Masha)
  iW    = NeighVertex(3, Ira);          jW   = NeighVertex(3, Masha)

  omega = OmegaMasha
  if(Is_delta_omega_not_valid(omega))  return

  omegaGLn2 = OmegaLn(GLn2) -(-1)**dir*omega
  if(flag==0 .and. Is_omega_not_valid(omegaGLn2))  return 

  k = kMasha
  kGLn2old = kLn(GLn2)
  kLn(GLn2) = add_k(kGLn2old, -(-1)**dir*k)

  !------------ step3 : configuration check ---------------------
  flag=0
  if(flag==0) then
    if(Is_reducible_G(GLn2))   flag=1
  endif
  if(flag==0) then
    if(Is_reducible_G_Gamma(GLn2))  flag=1
  endif

  if(flag==1) then
    kLn(GLn2) = kGLn2old
    return
  endif

  !------------ stat, type for the new configuration ------------
  statIra   = delete_ira_stat(StatusVertex(Ira))
  statMasha = delete_ira_stat(StatusVertex(Masha))

  if(iW/=jW) then
    statiW   = line_stat(statIra, StatusVertex(NeighLn(3-DirecVertex(Ira), iW)))
    statjW   = line_stat(statMasha, StatusVertex(NeighLn(3-DirecVertex(Masha), jW)))
  else
    statiW = line_stat(statIra, statMasha)
    statjW = line_stat(statIra, statMasha)
  endif

  typGLn2 = 3-TypeLn(GLn2)

  typiGamin  = TypeVertexIn(Ira)
  typiGamout = TypeVertexOut(Ira)
  typjGamin  = TypeVertexIn(Masha)
  typjGamout = TypeVertexOut(Masha)

  if(dir==1) then
    if(typiGamout==TypeLn(GLn1)) then
      typiGamin = 3 - typiGamin
    else
      typiGamout = 3 - typiGamout
    endif

    if(typjGamin==TypeLn(GLn3)) then
      typjGamout = 3 - typjGamout
    else
      typjGamin = 3 - typjGamin
    endif

    typIra = TypeGW2Gam(typGLn2, TypeLn(GLn1), typiGamin, typiGamout)
    typMasha = TypeGW2Gam(TypeLn(GLn3), typGLn2, typjGamin, typjGamout)
  else
    if(typiGamin==TypeLn(GLn1)) then
      typiGamout = 3 - typiGamout
    else
      typiGamin = 3 - typiGamin
    endif

    if(typjGamout==TypeLn(GLn3)) then
      typjGamin = 3 - typjGamin
    else
      typjGamout = 3 - typjGamout
    endif

    typIra = TypeGW2Gam(TypeLn(GLn1), typGLn2, typiGamin, typiGamout)
    typMasha = TypeGW2Gam(typGLn2, TypeLn(GLn3), typjGamin, typjGamout)
  endif

  if(iW/=jW) then
    if(DirecVertex(Ira)==1) then
      typiW = TypeGam2W(typIra, TypeVertex(NeighLn(2, iW)))
    else
      typiW = TypeGam2W(TypeVertex(NeighLn(1, iW)), typIra)
    endif

    if(DirecVertex(Masha)==1) then
      typjW = TypeGam2W(typMasha, TypeVertex(NeighLn(2, jW)))
    else
      typjW = TypeGam2W(TypeVertex(NeighLn(1, jW)), typMasha)
    endif
  else
    if(DirecVertex(Ira)==1) then
      typiW = TypeGam2W(typIra, typMasha)
    else if(DirecVertex(Masha)==1) then
      typiW = TypeGam2W(typMasha, typIra)
    endif
    typjW = typiW
  endif

  dxiW  = diff_x(WXVertex(Ira), WXVertex(NeighLn(3-DirecVertex(Ira),iW)))
  dyiW  = diff_y(WYVertex(Ira), WYVertex(NeighLn(3-DirecVertex(Ira),iW)))
 
  dxjW  = diff_x(WXVertex(NeighLn(3-DirecVertex(Masha),jW)), WXVertex(Masha))
  dyjW  = diff_y(WYVertex(NeighLn(3-DirecVertex(Masha),jW)), WYVertex(Masha))

  !------------ step4 : weight calculation ---------------------
  WiW = weight_line(statiW, 2, dxiW, dyiW, OmegaLn(iW), typiW)
  WjW = weight_line(statjW, 2, dxjW, dyjW, OmegaLn(jW), typjW)
  WGLn1 = weight_line(StatusLn(GLn1), 1, 0, 0, OmegaLn(GLn1), TypeLn(GLn1))
  WGLn3 = weight_line(StatusLn(GLn3), 1, 0, 0, OmegaLn(GLn3), TypeLn(GLn3))
  WGLn2 = weight_line(StatusLn(GLn2), 1, 0, 0, omegaGLn2, 3-TypeLn(GLn2))

  if(dir ==1) then
    WiGam = weight_vertex(Order, statIra, diff_x(GXVertex(Ira),WXVertex(Ira)), diff_y(GYVertex(Ira) &
      & ,WYVertex(Ira)), omegaGLn2, OmegaLn(GLn1), typIra)
    WjGam = weight_vertex(Order, statMasha, diff_x(GXVertex(Masha),WXVertex(Masha)), diff_y(GYVertex(Masha) &
      & ,WYVertex(Masha)), OmegaLn(GLn3), OmegaGLn2, typMasha)
  else
    WiGam = weight_vertex(Order, statIra, diff_x(GXVertex(Ira),WXVertex(Ira)), diff_y(GYVertex(Ira) &
      & ,WYVertex(Ira)), OmegaLn(GLn1), OmegaGLn2, typIra)
    WjGam = weight_vertex(Order, statMasha, diff_x(GXVertex(Masha),WXVertex(Masha)), diff_y(GYVertex(Masha) &
      & ,WYVertex(Masha)), omegaGLn2, OmegaLn(GLn3), typMasha)
  endif

  Anew = WiW *WjW *WiGam *WjGam *WGLn1 *WGLn2 *WGLn3
  Aold = WeightWorm *WeightLn(iW) *WeightLn(jW) *WeightLn(GLn1) *WeightLn(GLn2) &
    & *WeightLn(GLn3) *WeightVertex(Ira) *WeightVertex(Masha)

  if(GLn1==GLn3) then
    Anew = Anew/WGLn3
    Aold = Aold/WeightLn(GLn3)
  endif

  if(iW==jW) then
    Anew = Anew/WjW
    Aold = Aold/WeightLn(jW)
  endif

  call weight_ratio(Pacc, sgn, Anew, Aold)
  AveWeightRatio(Order, 2) = AveWeightRatio(Order, 2) + Pacc
  ProbProp(Order, iupdate) = ProbProp(Order, iupdate) + 1
  if(Order>0) then
    Pacc = Pacc *(CoefOfWorm*prob_omega(omega))/(2.d0*Order+2.d0) 
  else
    Pacc = Pacc *(CoefOfWorm*prob_omega(omega)) 
  endif
  Pacc = Pacc *Pupdate(3)/Pupdate(4)

  !------------ step5 : accept the update -----------------------
  rand = rn()
  if(rand<=Pacc) then

    !---------- update the diagram info ----------------
    Phase = Phase*sgn
    IsWormPresent = .false.

    !----------- update k and omega of elements --------
    call update_Hash4G(kGLn2old, kLn(GLn2))
    OmegaLn(GLn2) = omegaGLn2

    !----------- update type ---------------------------
    TypeLn(GLn2) = typGLn2

    TypeVertex(Ira) = typIra
    TypeVertex(Masha) = typMasha

    TypeVertexIn(Ira)    = typiGamin
    TypeVertexOut(Ira)   = typiGamout
    TypeVertexIn(Masha)  = typjGamin
    TypeVertexOut(Masha) = typjGamout

    TypeLn(iW) = typiW
    TypeLn(jW) = typjW

    !----------- update status of elements -------------
    StatusVertex(Ira)  = statIra
    StatusVertex(Masha)= statMasha
    StatusLn(iW)  = statiW
    StatusLn(jW)  = statjW

    !----------- update weight of elements -------------
    WeightVertex(Ira) = WiGam
    WeightVertex(Masha) = WjGam

    WeightLn(iW) = WiW
    WeightLn(jW) = WjW
    WeightLn(GLn1)  = WGLn1
    WeightLn(GLn2)  = WGLn2
    WeightLn(GLn3)  = WGLn3

    call update_weight(Anew, Aold)

    ProbAcc(Order, 4) = ProbAcc(Order, 4) + 1
  else
    kLn(GLn2) = kGLn2old
  endif

  return
END SUBROUTINE delete_worm_along_gline



!----- move worm along wline : Pupdate(5) --------------
!-----   \2       4/                  \2        4/    ------
!-----Ira ========= jGam     ==>  iGam ========== Ira ------
!-----   /1       3\                  /1        3\    ------
!-----------------------------------------------------------
SUBROUTINE move_worm_along_wline
  implicit none
  integer :: iGam, jGam, WLn, GLn1, GLn2, GLn3, GLn4
  integer :: dir, omegaW, kW
  integer :: dxW, dyW
  integer :: typiGam, typW
  integer :: statiGam, statjGam, statW, statG1, statG2, statG3, statG4
  integer :: tmpstat1, tmpstat2, tmpstat3, tmpstat4
  double precision :: WWorm, WiGam, WjGam, WW, WG1, WG2, WG3, WG4
  double precision :: Anew, Aold, Pacc, sgn
  integer :: kWold
  integer :: flag

  !------- step1 : check if worm is present -------------
  if(IsWormPresent .eqv. .false.)  return
  !ProbProp(iupdate) = ProbProp(iupdate) + 1

  !------- step2 : propose a new config -----------------
  dir  = DirecVertex(Ira)
  WLn  = NeighVertex(3, Ira)
  jGam = NeighLn(3-dir, WLn)
  if(jGam==Masha .or. jGam ==Ira)      return

  omegaW = OmegaLn(WLn) +(-1)**dir*OmegaMasha
  if(Is_omega_not_valid(OmegaW))  return

  GLn1 = NeighVertex(1, Ira)
  GLn2 = NeighVertex(2, Ira)
  GLn3 = NeighVertex(1, jGam)
  GLn4 = NeighVertex(2, jGam)

  kWold = kLn(WLn)
  kLn(WLn) = add_k(kWold, (-1)**dir*kMasha)

  !------- step3 : configuration check ------------------
  flag = 0
  if(flag==0) then
    if(Is_reducible_W(WLn))         flag = 1
  endif
  if(flag==0) then
    if(Is_reducible_W_Gamma(WLn))   flag = 1
  endif

  if(flag == 1) then
    kLn(WLn) = kWold
    return
  endif

  !-------- the stat, type, dx, dy of the new config ----
  dxW = diff_x(WXVertex(jGam), WXVertex(Ira))
  dyW = diff_y(WYVertex(jGam), WYVertex(Ira))


  statW = StatusLn(WLn)
  statiGam = delete_ira_stat(StatusVertex(Ira))
  statjGam = add_ira_stat(StatusVertex(jGam))

  typiGam = TypeGW2Gam(TypeLn(GLn1), TypeLn(GLn2), TypeVertexIn(Ira), &
    & TypeVertexOut(Ira))
  typW =  TypeLn(WLn)

  !------- step4 : weight calculation -------------------
  WW  = weight_line(statW, 2, dxW, dyW, omegaW, typW)
  WG1 = weight_line(StatusLn(GLn1), 1, 0, 0, OmegaLn(GLn1), TypeLn(GLn1))
  WG2 = weight_line(StatusLn(GLn2), 1, 0, 0, OmegaLn(GLn2), TypeLn(GLn2))
  WG3 = weight_line(StatusLn(GLn3), 1, 0, 0, OmegaLn(GLn3), TypeLn(GLn3))
  WG4 = weight_line(StatusLn(GLn4), 1, 0, 0, OmegaLn(GLn4), TypeLn(GLn4))

  WiGam = weight_vertex(Order, statiGam, diff_x(GXVertex(Ira),WXVertex(Ira)), diff_y(GYVertex(Ira) &
    & ,WYVertex(Ira)), OmegaLn(GLn1), OmegaLn(GLn2), typiGam)
  WjGam = weight_vertex(Order, statjGam, diff_x(GXVertex(jGam),WXVertex(jGam)), diff_y(GYVertex(jGam) &
    & ,WYVertex(jGam)), OmegaLn(GLn3), OmegaLn(GLn4), TypeVertex(jGam))

  WWorm= weight_worm(diff_x(GXVertex(jGam),GXVertex(Masha)),diff_y(GYVertex(jGam),GYVertex(Masha)), &
    & diff_x(WXVertex(jGam),WXVertex(Masha)),diff_y(WYVertex(jGam),WYVertex(Masha)), OmegaMasha)

  Anew = WWorm *WW *WG1 *WG2 *WG3 *WG4 *WiGam *WjGam
  Aold = WeightWorm *WeightLn(WLn) *WeightLn(GLn1) *WeightLn(GLn2) *WeightLn(GLn3) &
    & *WeightLn(GLn4) *WeightVertex(Ira) *WeightVertex(jGam)

  if(GLn1==GLn2) then
    Anew = Anew/WG2
    Aold = Aold/WeightLn(GLn2)
  endif

  if(GLn3==GLn4) then
    Anew = Anew/WG4
    Aold = Aold/WeightLn(GLn4)
  endif

  if(GLn1==GLn4) then
    Anew = Anew/WG4
    Aold = Aold/WeightLn(GLn4)
  endif

  if(GLn2==GLn3) then
    Anew = Anew/WG3
    Aold = Aold/WeightLn(Gln3)
  endif

  call weight_ratio(Pacc, sgn, Anew, Aold)

  !------- step5 : accept the update --------------------
  ProbProp(Order, iupdate) = ProbProp(Order, iupdate) + 1
  if(rn()<=Pacc) then

    !------ update diagram info ---------------
    Phase = Phase *sgn

    !------ update Ira and Masha --------------
    iGam = Ira
    Ira = jGam

    !------ update k and omega ----------------
    call update_Hash4W(kWold, kLn(WLn))
    OmegaLn(WLn) = omegaW

    !------ update type of elements -----------
    TypeLn(WLn) = typW
    TypeVertex(iGam) = typiGam

    !------ update status of elements ---------
    StatusLn(WLn) = statW
    StatusVertex(iGam) = statiGam
    StatusVertex(jGam) = statjGam

    !------ update weight of elements ---------
    WeightWorm = WWorm
    WeightLn(WLn) = WW
    WeightLn(GLn1) = WG1
    WeightLn(GLn2) = WG2
    WeightLn(GLn3) = WG3
    WeightLn(GLn4) = WG4

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
  integer :: dir 
  integer :: iGam, jGam, iW, jW, GLn1, GLn2, GLn3
  integer :: typiGam, typiGamin, typiGamout, typjGamin, typjGamout, typiW
  integer :: omegaG, kG
  integer :: dxiW, dyiW, dxjW, dyjW
  integer :: statiGam, statjGam, statG1, statG2, statG3, statiW, statjW
  double precision :: WWorm, WiW, WjW, WiGam, WjGam, WG1, WG2, WG3
  double precision :: Anew, Aold, Pacc, sgn
  integer :: kGold
  integer :: flag

  !------- step1 : check if worm is present -------------
  if(IsWormPresent .eqv. .false.)   return
  !ProbProp(iupdate) = ProbProp(iupdate) + 1

  !------- step2 : propose a new config -----------------
  iGam = Ira
  dir = Floor(rn()*2.d0)+1
  GLn2 = NeighVertex(dir, iGam)
  jGam = NeighLn(dir, GLn2)
  if(jGam==Masha .or. jGam==Ira)    return     

  if((dir+SpinMasha==3 .or. dir+SpinMasha==0)) then      !! dir=1 spin=2 or dir=2 spin=-2
    if(TypeLn(GLn2)==1)      return
  else                                                 !! dir=1 spin=-2 or dir=2 spin=2
    if(TypeLn(GLn2)==2)      return
  endif

  iW = NeighVertex(3, iGam)
  jW = NeighVertex(3, jGam)

  GLn1 = NeighVertex(3-dir, iGam)
  GLn3 = NeighVertex(dir, jGam)

  omegaG = OmegaLn(GLn2)-(-1)**dir *OmegaMasha
  if(Is_omega_not_valid(omegaG))  return


  kGold = kLn(GLn2)
  kLn(GLn2) = add_k(kLn(GLn2), -(-1)**dir *kMasha)

  !------- step3 : configuration check ------------------
  flag = 0
  if(flag==0) then
    if(Is_reducible_G(GLn2))          flag = 1
  endif
  if(flag==0) then
    if(Is_reducible_G_Gamma(GLn2))    flag = 1
  endif

  if(flag == 1) then
    kLn(GLn2) = kGold
    return
  endif

  !----- site, type and status of the new config ---

  statiGam = delete_ira_stat(StatusVertex(iGam))
  statjGam = add_ira_stat(StatusVertex(jGam))
  if(iW/=jW) then
    statiW   = line_stat(statiGam, StatusVertex(NeighLn(3-DirecVertex(iGam), iW)))
    statjW   = line_stat(statjGam, StatusVertex(NeighLn(3-DirecVertex(jGam), jW)))
  else
    statiW = line_stat(statiGam, statjGam)
    statjW = line_stat(statiGam, statjGam)
  endif

  typiGamin  = TypeVertexIn(iGam)
  typiGamout = TypeVertexOut(iGam)
  typjGamin  = TypeVertexIn(jGam)
  typjGamout = TypeVertexOut(jGam)

  if(dir==1) then
    if(typiGamout==TypeLn(GLn1)) then
      typiGamin = 3 - typiGamin
    else
      typiGamout = 3 - typiGamout
    endif

    if(typjGamin==TypeLn(GLn3)) then
      typjGamout = 3 - typjGamout
    else
      typjGamin = 3 - typjGamin
    endif

    typiGam = TypeGW2Gam(3-TypeLn(GLn2), TypeLn(GLn1), typiGamin, &
      & typiGamout)
  else
    if(typiGamin==TypeLn(GLn1)) then
      typiGamout = 3 - typiGamout
    else
      typiGamin = 3 - typiGamin
    endif

    if(typjGamout==TypeLn(GLn3)) then
      typjGamin = 3 - typjGamin
    else
      typjGamout = 3 - typjGamout
    endif

    typiGam = TypeGW2Gam(TypeLn(GLn1), 3-TypeLn(GLn2), typiGamin, &
      & typiGamout)
  endif

  if(DirecVertex(iGam)==1) then
    if(statiW<=1) then
      typiW = TypeGam2W(typiGam, TypeVertex(NeighLn(2, iW)))
    else
      typiW = TypeLn(iW)
    endif
  else
    if(statiW<=1) then
      typiW = TypeGam2W(TypeVertex(NeighLn(1, iW)), typiGam)
    else
      typiW = TypeLn(iW)
    endif
  endif

  dxiW  = diff_x(WXVertex(iGam), WXVertex(NeighLn(3-DirecVertex(iGam),iW)))
  dyiW  = diff_y(WYVertex(iGam), WYVertex(NeighLn(3-DirecVertex(iGam),iW)))
 
  dxjW  = diff_x(WXVertex(NeighLn(3-DirecVertex(jGam),jW)), WXVertex(jGam))
  dyjW  = diff_y(WYVertex(NeighLn(3-DirecVertex(jGam),jW)), WYVertex(jGam))

  !------- step4 : weight calculation -------------------
  WiW = weight_line(statiW, 2, dxiW, dyiW, OmegaLn(iW), typiW)
  WjW = weight_line(statjW, 2, dxjW, dyjW, OmegaLn(jW), TypeLn(jW))
  WG1 = weight_line(StatusLn(GLn1), 1, 0, 0, OmegaLn(GLn1), TypeLn(GLn1))
  WG3 = weight_line(StatusLn(GLn3), 1, 0, 0, OmegaLn(GLn3), TypeLn(GLn3))
  WG2 = weight_line(StatusLn(GLn2), 1, 0, 0, omegaG, 3-TypeLn(GLn2))

  if(dir==1) then
    WiGam = weight_vertex(Order, statiGam, diff_x(GXVertex(iGam),WXVertex(iGam)), diff_y(GYVertex(iGam) &
      & ,WYVertex(iGam)), omegaG, OmegaLn(GLn1), typiGam)
    WjGam = weight_vertex(Order, statjGam, diff_x(GXVertex(jGam),WXVertex(jGam)), diff_y(GYVertex(jGam) &
      & ,WYVertex(jGam)), OmegaLn(GLn3), omegaG, TypeVertex(jGam))
  else
    WiGam = weight_vertex(Order, statiGam, diff_x(GXVertex(iGam),WXVertex(iGam)), diff_y(GYVertex(iGam) &
      & ,WYVertex(iGam)), OmegaLn(GLn1), omegaG, typiGam)
    WjGam = weight_vertex(Order, statjGam, diff_x(GXVertex(jGam),WXVertex(jGam)), diff_y(GYVertex(jGam) &
      & ,WYVertex(jGam)), omegaG, OmegaLn(GLn3), TypeVertex(jGam))
  endif
  WWorm= weight_worm(diff_x(GXVertex(jGam),GXVertex(Masha)),diff_y(GYVertex(jGam),GYVertex(Masha)), &
    & diff_x(WXVertex(jGam),WXVertex(Masha)),diff_y(WYVertex(jGam),WYVertex(Masha)), omegaMasha)

  Anew = WWorm *WiW *WjW *WiGam *WjGam *WG1 *WG2 *WG3
  Aold = WeightWorm *WeightLn(iW) *WeightLn(jW) *WeightLn(GLn1) *WeightLn(GLn2) &
    & *WeightLn(GLn3) *WeightVertex(iGam) *WeightVertex(jGam)

  if(GLn1==GLn3) then
    Anew = Anew/WG3
    Aold = Aold/WeightLn(GLn3)
  endif

  if(iW==jW) then
    Anew = Anew/WjW
    Aold = Aold/WeightLn(jW)
  endif

  call weight_ratio(Pacc, sgn, Anew, Aold)

  !------- step5 : accept the update --------------------
  ProbProp(Order, iupdate) = ProbProp(Order, iupdate) + 1
  if(rn()<=Pacc) then

    !----- update the diagram info -------------- 
    Phase = Phase *sgn

    !---- update Ira and Masha ------------------
    Ira = jGam

    !---- update omega and k --------------------
    call update_Hash4G(kGold, kLn(GLn2))
    OmegaLn(GLn2) = omegaG
    
    !----- update type configuration ------------
    TypeLn(GLn2) = 3-TypeLn(GLn2)

    TypeVertexIn(iGam)  = typiGamin
    TypeVertexOut(iGam) = typiGamout
    TypeVertexIn(jGam)  = typjGamin
    TypeVertexOut(jGam) = typjGamout

    TypeLn(iW) = typiW
    TypeVertex(iGam) = typiGam

    !---- update the status of elements ---------
    StatusVertex(iGam) = statiGam
    StatusVertex(jGam) = statjGam
    StatusLn(iW) = statiW
    StatusLn(jW) = statjW

    !---- update the weight of elements ---------
    WeightWorm = WWorm
    WeightVertex(iGam) = WiGam
    WeightVertex(jGam) = WjGam
    WeightLn(iW) = WiW
    WeightLn(jW) = WjW
    WeightLn(GLn1) = WG1
    WeightLn(GLn2) = WG2
    WeightLn(GLn3) = WG3

    call update_weight(Anew, Aold)

    ProbAcc(Order, 6) = ProbAcc(Order, 6) + 1
  else
    kLn(GLn2) = kGold
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
SUBROUTINE add_interaction
  implicit none
  integer :: i, dir, dirW, q, omega
  integer :: kIA, omegaIA, kMB, omegaMB, kM, omegaM
  integer :: dxwA, dywA, dxwB, dywB
  integer :: xwA, ywA, xwB, ywB
  integer :: WAB, GIC, GMD, GIA, GMB
  integer :: GamA, GamB, GamC, GamD
  integer :: sp, typGamA, typGamB, typAB, statIA, statAC, statMB, statBD
  double precision :: Pacc, Anew, Aold, sgn
  double precision :: WWorm, WA, WB, WGIA, WGAC, WGMB, WGBD, WIra, WMasha, WWAB
  double precision :: WMeasGamma
  double precision :: rand
  integer :: confold, confnew
  integer :: flag

  !------------ step1 : check if worm is present ------------------
  if(IsWormPresent .eqv. .false.)    return
  !ProbProp(iupdate) = ProbProp(iupdate) + 1
  if(Order >= MCOrder) return    


  !------------ step2 : propose the new config ------------------
  dir = Floor(rn()*2.d0)+1;       dirW = Floor(rn()*2.d0)+1

  xwA = GXVertex(Ira);            ywA = GYVertex(Ira)
  xwB = GXVertex(Masha);          ywB = GYVertex(Masha)

  q    = generate_k()
  omega= generate_omega()
  kM = add_k(kMasha, -(-1)**dirW *q)
  OmegaM = OmegaMasha - (-1)**dirW *omega

  GIC  = NeighVertex(dir, Ira);  GMD  = NeighVertex(dir, Masha)
  GamC = NeighLn(dir, GIC);      GamD = NeighLn(dir, GMD)

  if(dir+dirW == 3) then
    kIA = add_k(kLn(GIC), q);          kMB = add_k(kLn(GMD), -q)
    omegaIA = OmegaLn(GIC) + omega;    omegaMB = OmegaLn(GMD) - omega
  else
    kIA = add_k(kLn(GIC),-q);          kMB = add_k(kLn(GMD),  q)
    omegaIA = OmegaLn(GIC) - omega;    omegaMB = OmegaLn(GMD) + omega
  endif
  if(Is_omega_not_valid(omegaIA)) return 
  if(Is_omega_not_valid(omegaMB)) return
  if(Is_omega_not_valid(omegaM))  return

  if(Is_k_valid(kM) .eqv. .false.)  return 
  if(kIA==kMB)  return       
  if(Is_k_valid(add_k(kIA, -kMB)).eqv. .false.)  return
  if(Hash4W(abs(add_k(kIA,-kMB)))==1)  return
  if(Is_k_valid(add_k(kIA, -kLn(GMD))).eqv. .false.)  return
  if(Hash4W(abs(add_k(kIA,-kLn(GMD))))==1)  return
  if(Is_k_valid(add_k(kLn(GIC), -kMB)).eqv. .false.)  return
  if(Hash4W(abs(add_k(kLn(GIC),-kMB)))==1)  return

  !-------- the new spin, type and status for the new config --
  typGamA = TypeGW2Gam(TypeLn(GIC),TypeLn(GIC),TypeLn(GIC),TypeLn(GIC))
  typGamB = TypeGW2Gam(TypeLn(GMD),TypeLn(GMD),TypeLn(GMD),TypeLn(GMD))
   
  if(dirW==1) then
    typAB = TypeGam2W(typGamA, typGamB)
  else
    typAB = TypeGam2W(typGamB, typGamA)
  endif

  statIA = Mod(StatusVertex(Ira),2)
  statMB = Mod(StatusVertex(Masha),2)
  statAC = Mod(StatusVertex(GamC),2)
  statBD = Mod(StatusVertex(GamD),2)

  !----------- step4 : weight calculation --------------------
  if(dir==1) then
    WA = weight_vertex(Order+1, 0, 0, 0, OmegaLn(GIC), omegaIA, TypeLn(GIC))
    WB = weight_vertex(Order+1, 0, 0, 0, OmegaLn(GMD), omegaMB, TypeLn(GMD))
    WIra = weight_vertex(Order+1, StatusVertex(Ira), diff_x(GXVertex(Ira),WXVertex(Ira)), &
      & diff_y(GYVertex(Ira),WYVertex(Ira)), omegaIA, &
      & OmegaLn(NeighVertex(2, Ira)), TypeVertex(Ira))
    WMasha = weight_vertex(Order+1, StatusVertex(Masha), diff_x(GXVertex(Masha),WXVertex(Masha)),  &
      & diff_y(GYVertex(Masha),WYVertex(Masha)), &
      & omegaMB, OmegaLn(NeighVertex(2, Masha)), TypeVertex(Masha))
  else
    WA = weight_vertex(Order+1, 0, 0, 0, omegaIA, OmegaLn(GIC), TypeLn(GIC))
    WB = weight_vertex(Order+1, 0, 0, 0, omegaMB, OmegaLn(GMD), TypeLn(GMD))
    WIra = weight_vertex(Order+1, StatusVertex(Ira), diff_x(GXVertex(Ira),WXVertex(Ira)), &
      & diff_y(GYVertex(Ira),WYVertex(Ira)), OmegaLn(NeighVertex(1, Ira)), &
      & omegaIA, TypeVertex(Ira))
    WMasha = weight_vertex(Order+1, StatusVertex(Masha), diff_x(GXVertex(Masha),WXVertex(Masha)),  &
      & diff_y(GYVertex(Masha),WYVertex(Masha)), OmegaLn(NeighVertex(1, Masha)), &
      & omegaMB, TypeVertex(Masha))
  endif


  WWAB = weight_line(0, 2, diff_x(xwA,xwB), diff_y(ywA,ywB), omega, typAB)
  WGIA = weight_line(statIA, 1, 0, 0, omegaIA, TypeLn(GIC)) 
  WGAC = weight_line(statAC, 1, 0, 0, OmegaLn(GIC), TypeLn(GIC)) 
  WGMB = weight_line(statMB, 1, 0, 0, omegaMB, TypeLn(GMD)) 
  WGBD = weight_line(statBD, 1, 0, 0, OmegaLn(GMD), TypeLn(GMD)) 
  WWorm= weight_worm(diff_x(GXVertex(Ira),GXVertex(Masha)),diff_y(GYVertex(Ira),GYVertex(Masha)), &
    & diff_x(WXVertex(Ira),WXVertex(Masha)),diff_y(WYVertex(Ira),WYVertex(Masha)), omegaM)
    
  !----------  change the topology for the configuration after update --
  call insert_gamma(GamA, GXVertex(Ira), GYVertex(Ira), xwA, ywA, dirW, typGamA, 0, WA)
  call insert_gamma(GamB, GXVertex(Masha), GYVertex(Masha), xwB, ywB, 3-dirW, typGamB, 0, WB)
  call insert_line(WAB, omega, q, 2, typAB, 0, WWAB)
  call insert_line(GIA, omegaIA, kIA, 1, TypeLn(GIC), statIA, WGIA)
  call insert_line(GMB, omegaMB, kMB, 1, TypeLn(GMD), statMB, WGMB)

  !------------- update the topology -----------------------------
  NeighVertex(dir, GamA) = GIC;     NeighVertex(3-dir, GamA) = GIA
  NeighVertex(dir, GamB) = GMD;     NeighVertex(3-dir, GamB) = GMB
  NeighVertex(3, GamA)   = WAB;     NeighVertex(3, GamB)     = WAB

  NeighLn(dir, GIA)   = GamA;    NeighLn(3-dir, GIA)   = Ira
  NeighLn(dir, GMB)   = GamB;    NeighLn(3-dir, GMB)   = Masha

  NeighLn(dirW, WAB)  = GamA;    NeighLn(3-dirW, WAB)  = GamB

  NeighVertex(dir, Ira)  = GIA;     NeighVertex(dir, Masha)  = GMB
  NeighLn(3-dir, GIC) = GamA;    NeighLn(3-dir, GMD)   = GamB 


  !----------- step3 : configuration check --------------------
  flag = 0
  
  if(flag==0 ) then
    if(Is_reducible_G(GIA))       flag = 1
  endif
  if(flag==0) then
    if(Is_reducible_G(GMB))       flag = 2
  endif
  if(flag==0) then
    if(Is_reducible_W(WAB))       flag = 3
  endif
  if(flag==0) then
    if(Is_reducible_G_Gamma(GIA))       flag = 4
  endif
  if(flag==0) then
    if(Is_reducible_G_Gamma(GMB))       flag = 5
  endif
  if(flag==0) then
    if(Is_reducible_G_Gamma(GIC))       flag = 6
  endif
  if(flag==0) then
    if(Is_reducible_G_Gamma(GMD))       flag = 7
  endif
  if(flag==0) then
    if(Is_reducible_W_Gamma(WAB))       flag = 8
  endif

  if(flag >= 1) then
    !------------- delete line and vertexes --------------------
    call undo_insert_line(GIA,1)
    call undo_insert_line(GMB,1)
    call undo_insert_line(WAB,2)
    call delete_gamma(GamA)
    call delete_gamma(GamB)

    NeighVertex(dir, Ira)  = GIC;     NeighVertex(dir, Masha)  = GMD
    NeighLn(3-dir, GIC) = Ira;     NeighLn(3-dir, GMD)   = Masha


    return
  endif

  !------------ weight calculation ----------------------------
  if(MeasGamma/=Ira .and. MeasGamma/=Masha) then
    WMeasGamma = weight_vertex(Order+1, StatusVertex(MeasGamma), &
      & diff_x(GXVertex(MeasGamma), WXVertex(MeasGamma)), diff_y(GYVertex(MeasGamma), &
      & WYVertex(MeasGamma)), OmegaLn(NeighVertex(1, MeasGamma)), &
      & OmegaLn(NeighVertex(2, MeasGamma)), TypeVertex(MeasGamma)) 
    Anew = CoefOfWeight(Order+1)*WWorm *WMeasGamma *WA *WB *WWAB *WIra *WMasha *WGIA *WGAC &
      & *WGMB *WGBD*(1.d0)/Beta
    Aold = CoefOfWeight(Order)*WeightWorm *WeightVertex(MeasGamma) *WeightLn(GIC) *WeightLn(GMD)&
      & *WeightVertex(Ira) *WeightVertex(Masha)
  else
    Anew = CoefOfWeight(Order+1)*WWorm *WA *WB *WWAB *WIra *WMasha *WGIA *WGAC &
      & *WGMB *WGBD*(1.d0)/Beta
    Aold = CoefOfWeight(Order)*WeightWorm *WeightLn(GIC) *WeightLn(GMD)&
      & *WeightVertex(Ira) *WeightVertex(Masha)
  endif

  call weight_ratio(Pacc, sgn, Anew, Aold)

  Pacc = Pacc/prob_omega(omega)
  Pacc = Pacc *Pupdate(8)/Pupdate(7)

  !write(*, *) iupdate, Pacc
  !------------ step5 : accept the update -----------------------
  ProbProp(Order, iupdate) = ProbProp(Order, iupdate) + 1
  rand = rn()
  if(rand<=Pacc) then

    !-------------- update the diagram info --------------------
    Order = Order + 1
    Phase = Phase *sgn

    !------------- update k and omega -------------------------
    kMasha = kM
    OmegaMasha = OmegaM

    call add_Hash4G(kIA)
    call add_Hash4G(kMB)
    call add_Hash4W(q)

    !------------- update the status of elements --------------
    StatusLn(GMD) = statBD
    StatusLn(GIC) = statAC

    !------------- update weight of elements ------------------
    WeightWorm = WWorm
    WeightVertex(Ira) = WIra
    WeightVertex(Masha) = WMasha
    WeightLn(GIC) = WGAC
    WeightLn(GMD) = WGBD
    if(MeasGamma/=Ira .and. MeasGamma/=Masha)  WeightVertex(MeasGamma) = WMeasGamma

    call update_weight(Anew, Aold)

    ProbAcc(Order-1, 7) = ProbAcc(Order-1, 7) + 1
  else
    
    !------------- delete line and vertexes --------------------
    call undo_insert_line(GIA,1)
    call undo_insert_line(GMB,1)
    call undo_insert_line(WAB,2)
    call delete_gamma(GamA)
    call delete_gamma(GamB)

    NeighVertex(dir, Ira)  = GIC;     NeighVertex(dir, Masha)  = GMD
    NeighLn(3-dir, GIC) = Ira;     NeighLn(3-dir, GMD)   = Masha

  endif
  !call print_markov_after
  return
END SUBROUTINE add_interaction
!-----------------------------------------------------------------




!------------ remove interaction : Pupdate(8) ---------------
!---------------------- dir = 1, dirW = 1 ------------------------
!------               GamC                     GamA     GamC -----
!------    Ira ----<----            Ira ----<---||---<-----  -----
!------                     <==                 \/           -----
!------  Masha ----<----          Masha ----<---||---<-----  -----
!------               GamD                     GamB     GamD -----
!-----------------------------------------------------------------
!---------------------- dir = 1, dirW = 2 ------------------------
!------               GamC                     GamA     GamC -----
!------    Ira ----<----            Ira ----<---||---<-----  -----
!------                     <==                 /\           -----
!------  Masha ----<----          Masha ----<---||---<-----  -----
!------               GamD                     GamB     GamD -----
!-----------------------------------------------------------------
!---------------------- dir = 2, dirW = 1 ------------------------
!------               GamC                     GamA     GamC -----
!------    Ira ---->----            Ira ---->---||--->-----  -----
!------                     <==                 \/           -----
!------  Masha ---->----          Masha ---->---||--->-----  -----
!------               GamD                     GamB     GamD -----
!-----------------------------------------------------------------
!---------------------- dir = 2, dirW = 2 ------------------------
!------               GamC                     GamA     GamC -----
!------    Ira ---->----            Ira ---->---||--->-----  -----
!------                     <==                 /\           -----
!------  Masha ---->----          Masha ---->---||--->-----  -----
!------               GamD                     GamB     GamD -----
!-----------------------------------------------------------------
SUBROUTINE remove_interaction
  implicit none
  integer :: i, dir, dirW, omega
  integer :: omegaM, kM
  integer :: xwA, ywA, xwB, ywB
  integer :: WAB, GIA, GAC, GMB, GBD
  integer :: GamA, GamB, GamC, GamD
  integer :: statIC, statMD
  double precision :: WWorm, WGIC, WGMD, WIra, WMasha
  double precision :: WMeasGamma
  double precision :: Pacc, Anew, Aold, sgn
  double precision :: rand
  integer :: flag

  !------------ step1 : check if worm is present ------------------
  if(IsWormPresent .eqv. .false.)    return
  !ProbProp(iupdate) = ProbProp(iupdate) + 1
  if(Order <= 0) return          !!! just for 1st and higher order

  !------------ step2 : propose the new config ------------------
  dir  = Floor(rn()*2.d0)+1
  GIA  = NeighVertex(dir, Ira);          GMB  = NeighVertex(dir, Masha)
  GamA = NeighLn(dir, GIA);           GamB = NeighLn(dir, GMB)
  if(TypeVertex(GamA)>2) return 
  if(TypeVertex(GamB)>2) return 

  WAB  = NeighVertex(3, GamA)
  if(StatusLn(WAB)/=0)            return
  if(WAB /= NeighVertex(3, GamB))    return 

  GAC  = NeighVertex(dir, GamA);         GBD  = NeighVertex(dir, GamB)
  GamC = NeighLn(dir, GAC);              GamD = NeighLn(dir, GBD)

  dirW  = DirecVertex(GamA)
  xwA   = WXVertex(GamA);                ywA  = WYVertex(GamA)
  xwB   = WXVertex(GamB);                ywB  = WYVertex(GamB)
  if(xwA/=GXVertex(Ira))   return 
  if(ywA/=GYVertex(Ira))   return 
  if(xwB/=GXVertex(Masha))   return 
  if(ywB/=GYVertex(Masha))   return 

  omega = OmegaLn(WAB)
  if(Is_delta_omega_not_valid(omega)) return

  kM = add_k(kMasha, (-1)**dirW *kLn(WAB))
  if(Is_k_valid(kM) .eqv. .false.)    return

  OmegaM = OmegaMasha +(-1)**dirW *omega
  if(Is_omega_not_valid(OmegaM))      return

  if(Ira==MeasGamma .or. GamC==MeasGamma) then
    statIC = 1
  else
    statIC = 0
  endif
  if(Masha==MeasGamma .or. GamD==MeasGamma) then
    statMD = 1
  else
    statMD = 0
  endif
  !------------- change the topology -------------------------
  !------------- delete line and vertexes --------------------
  call delete_line(GIA,1)
  call delete_line(GMB,1)
  call delete_line(WAB,2)
  call delete_gamma(GamA)
  call delete_gamma(GamB)

  !------------- update the topology ------------------------
  NeighVertex(dir, Ira) = GAC;      NeighVertex(dir, Masha) = GBD
  NeighLn(3-dir, GAC) = Ira;     NeighLn(3-dir, GBD) = Masha

  !----------- step3 : configuration check ----------------------

  flag = 0
  if(flag==0) then
    if(NeighVertex(1, Ira)==NeighVertex(2, Ira))  flag = 1
  endif

  if(flag==0) then
    if(NeighVertex(1, Masha)==NeighVertex(2, Masha))  flag = 1
  endif

  if(flag == 1) then
    !----------  change the topology for the configuration after update --
    call undo_delete_line(WAB,2, 0)
    call undo_delete_line(GMB,1, Mod(StatusVertex(Masha),2))
    call undo_delete_line(GIA,1, Mod(StatusVertex(Ira),2))

    call undo_delete_gamma(GamB)
    call undo_delete_gamma(GamA)

    !------------- update the topology ------------------------
    NeighVertex(dir, Ira) = GIA;      NeighVertex(dir, Masha) = GMB
    NeighLn(3-dir, GAC) = GamA;    NeighLn(3-dir, GBD) = GamB
    return
  endif

  !----------- step4 : weight calculation -----------------------
  WGIC = weight_line(statIC, 1, 0, 0, OmegaLn(GAC), TypeLn(GAC))
  WGMD = weight_line(statMD, 1, 0, 0, OmegaLn(GBD), TypeLn(GBD))
  if(dir==1) then
    WIra = weight_vertex(Order-1, StatusVertex(Ira), diff_x(GXVertex(Ira),WXVertex(Ira)), &
      & diff_y(GYVertex(Ira),WYVertex(Ira)),  &
      & OmegaLn(GAC), OmegaLn(NeighVertex(2, Ira)), TypeVertex(Ira))
    WMasha = weight_vertex(Order-1, StatusVertex(Masha), diff_x(GXVertex(Masha),WXVertex(Masha)), &
      & diff_y(GYVertex(Masha),WYVertex(Masha)), OmegaLn(GBD),  &
      & OmegaLn(NeighVertex(2,Masha)), TypeVertex(Masha))
  else
    WIra = weight_vertex(Order-1, StatusVertex(Ira), diff_x(GXVertex(Ira),WXVertex(Ira)), &
      & diff_y(GYVertex(Ira),WYVertex(Ira)),  &
      & OmegaLn(NeighVertex(1, Ira)), OmegaLn(GAC), TypeVertex(Ira))
    WMasha = weight_vertex(Order-1, StatusVertex(Masha), diff_x(GXVertex(Masha),WXVertex(Masha)), &
      & diff_y(GYVertex(Masha),WYVertex(Masha)), OmegaLn(NeighVertex(1, Masha)),  &
      & OmegaLn(GBD), TypeVertex(Masha))
  endif

  WWorm= weight_worm(diff_x(GXVertex(Ira),GXVertex(Masha)),diff_y(GYVertex(Ira),GYVertex(Masha)), &
    & diff_x(WXVertex(Ira),WXVertex(Masha)),diff_y(WYVertex(Ira),WYVertex(Masha)), omegaM)

  if(MeasGamma/=Ira .and. MeasGamma/=Masha) then
    WMeasGamma = weight_vertex(Order-1, StatusVertex(MeasGamma), &
      & diff_x(GXVertex(MeasGamma), WXVertex(MeasGamma)), diff_y(GYVertex(MeasGamma), &
      & WYVertex(MeasGamma)), OmegaLn(NeighVertex(1, MeasGamma)), &
      & OmegaLn(NeighVertex(2, MeasGamma)), TypeVertex(MeasGamma)) 

    Anew = CoefOfWeight(Order-1)*WWorm*WMeasGamma*WGIC *WGMD *WIra *WMasha
    Aold = CoefOfWeight(Order)*WeightWorm *WeightVertex(MeasGamma) *WeightLn(GIA) &
      & *WeightLn(GAC)*WeightLn(GMB)* WeightLn(GBD)*WeightVertex(Ira)*WeightVertex(Masha) &
      & *WeightVertex(GamA)*WeightVertex(GamB)*WeightLn(WAB)*(1.d0/Beta)
  else
    Anew = CoefOfWeight(Order-1)*WWorm*WGIC *WGMD *WIra *WMasha
    Aold = CoefOfWeight(Order)*WeightWorm *WeightLn(GIA) &
      & *WeightLn(GAC)*WeightLn(GMB)* WeightLn(GBD)*WeightVertex(Ira)*WeightVertex(Masha) &
      & *WeightVertex(GamA)*WeightVertex(GamB)*WeightLn(WAB)*(1.d0/Beta)
  endif

  call weight_ratio(Pacc, sgn, Anew, Aold)

  Pacc = Pacc *prob_omega(omega)
  Pacc = Pacc *Pupdate(7)/Pupdate(8)
  
  !write(*, *) iupdate, Pacc
  !------------ step5 : accept the update -----------------------
  ProbProp(Order, iupdate) = ProbProp(Order, iupdate) + 1
  rand = rn()
  !call print_markov_before(Pacc, rand)
  if(rand<=Pacc) then

    !-------------- update the diagram info --------------------
    Order = Order - 1
    Phase = Phase *sgn


    !------------- update k and omega of elements ------------
    kMasha = kM
    OmegaMasha = OmegaM
    call delete_Hash4G(kLn(GIA))
    call delete_Hash4G(kLn(GMB))

    !------------ update the status of elements -------------
    StatusLn(GAC) = statIC
    StatusLn(GBD) = statMD

    !------------ update the weight of elements -------------
    WeightWorm = WWorm
    WeightVertex(Ira) = WIra
    WeightVertex(Masha) = WMasha
    if(MeasGamma/=Ira .and. MeasGamma/=Masha)  WeightVertex(MeasGamma)=WMeasGamma
    WeightLn(GAC) = WGIC
    WeightLn(GBD) = WGMD

    call update_weight(Anew, Aold)

    ProbAcc(Order+1, 8) = ProbAcc(Order+1, 8) + 1
  else

    !----------  change the topology for the configuration after update --
    call undo_delete_line(WAB,2, 0)
    call undo_delete_line(GMB,1, Mod(StatusVertex(Masha),2))
    call undo_delete_line(GIA,1, Mod(StatusVertex(Ira),2))

    call undo_delete_gamma(GamB)
    call undo_delete_gamma(GamA)

    !------------- update the topology ------------------------
    NeighVertex(dir, Ira) = GIA;      NeighVertex(dir, Masha) = GMB
    NeighLn(3-dir, GAC) = GamA;    NeighLn(3-dir, GBD) = GamB

  endif
  !call print_markov_after
  return
END SUBROUTINE remove_interaction
!-----------------------------------------------------------



!------------- add interaction cross : Pupdate(9) -----------------------
!---------------------- dir = 1, dirW = 1 ------------------------
!------               GamC                     GamA     GamC -----
!------    Ira ----<----            Ira ----<---||---<-----  -----
!------                     ==>                 \/           -----
!------  Masha ---->----          Masha ---->---||--->-----  -----
!------               GamD                     GamB     GamD -----
!-----------------------------------------------------------------
!---------------------- dir = 1, dirW = 2 ------------------------
!------               GamC                     GamA     GamC -----
!------    Ira ----<----            Ira ----<---||---<-----  -----
!------                     ==>                 /\           -----
!------  Masha ---->----          Masha ---->---||--->-----  -----
!------               GamD                     GamB     GamD -----
!-----------------------------------------------------------------
!---------------------- dir = 2, dirW = 1 ------------------------
!------               GamC                     GamA     GamC -----
!------    Ira ---->----            Ira ---->---||--->-----  -----
!------                     ==>                 \/           -----
!------  Masha ----<----          Masha ----<---||---<-----  -----
!------               GamD                     GamB     GamD -----
!-----------------------------------------------------------------
!---------------------- dir = 2, dirW = 2 ------------------------
!------               GamC                     GamA     GamC -----
!------    Ira ---->----            Ira ---->---||--->-----  -----
!------                     ==>                 /\           -----
!------  Masha ----<----          Masha ----<---||---<-----  -----
!------               GamD                     GamB     GamD -----
!-----------------------------------------------------------------
!-----------------------------------------------------------------
SUBROUTINE add_interaction_cross
  implicit none
  integer :: i, dir, dirW, q, omega
  integer :: kIA, omegaIA, kMB, omegaMB, kM, omegaM
  integer :: dxwA, dywA, dxwB, dywB
  integer :: xwA, ywA, xwB, ywB
  integer :: WAB, GIC, GMD, GIA, GMB
  integer :: GamA, GamB, GamC, GamD
  integer :: sp, typGamA, typGamB, typAB, statIA, statAC, statMB, statBD
  double precision :: Pacc, Anew, Aold, sgn
  double precision :: WWorm, WA, WB, WGIA, WGAC, WGMB, WGBD, WIra, WMasha, WWAB
  double precision :: WMeasGamma
  double precision :: rand
  integer :: flag

  !------------ step1 : check if worm is present ------------------
  if(IsWormPresent .eqv. .false.)    return
  !ProbProp(iupdate) = ProbProp(iupdate) + 1
  if(Order >= MCOrder) return    


  !------------ step2 : propose the new config ------------------
  dir = Floor(rn()*2.d0)+1;    dirW = Floor(rn()*2.d0)+1

  GIC  = NeighVertex(dir, Ira);  GMD  = NeighVertex(3-dir, Masha)
  GamC = NeighLn(dir, GIC);      GamD = NeighLn(3-dir, GMD)
  if(GIC==GMD .or. GamC==GamD)  return

  omega= generate_omega()
  OmegaM = OmegaMasha - (-1)**dirW *omega
  q    = generate_k()
  kM = add_k(kMasha, -(-1)**dirW *q)
  if(dir+dirW == 3) then
    kIA = add_k(kLn(GIC), q);          kMB = add_k(kLn(GMD),  q)
    omegaIA = OmegaLn(GIC) + omega;    omegaMB = OmegaLn(GMD) + omega
  else
    kIA = add_k(kLn(GIC),-q);          kMB = add_k(kLn(GMD), -q)
    omegaIA = OmegaLn(GIC) - omega;    omegaMB = OmegaLn(GMD) - omega
  endif

  if(Is_omega_not_valid(omegaIA))  return
  if(Is_omega_not_valid(omegaMB))  return
  if(Is_omega_not_valid(omegaM))  return

  if((Is_k_valid(kM) .eqv. .false.))  return
  if(kIA==kMB)    return    
  if((Is_k_valid(add_k(kIA, -kMB)).eqv. .false.)) return
  if(Hash4W(abs(add_k(kIA,-kMB)))==1) return
  if((Is_k_valid(add_k(kIA, -kLn(GMD))).eqv. .false.)) return
  if(Hash4W(abs(add_k(kIA,-kLn(GMD))))==1) return
  if((Is_k_valid(add_k(kLn(GIC), -kMB)).eqv. .false.)) return
  if(Hash4W(abs(add_k(kLn(GIC),-kMB)))==1) return
  

  xwA = GXVertex(Ira);            ywA = GYVertex(Ira)
  xwB = GXVertex(Masha);          ywB = GYVertex(Masha)

  !-------- the new spin, type and status for the new config --
  typGamA = TypeGW2Gam(TypeLn(GIC),TypeLn(GIC),TypeLn(GIC),TypeLn(GIC))
  typGamB = TypeGW2Gam(TypeLn(GMD),TypeLn(GMD),TypeLn(GMD),TypeLn(GMD))
   
  if(dirW==1) then
    typAB = TypeGam2W(typGamA, typGamB)
  else
    typAB = TypeGam2W(typGamB, typGamA)
  endif

  statIA = Mod(StatusVertex(Ira),2);    statMB = Mod(StatusVertex(Masha),2)
  statAC = Mod(StatusVertex(GamC),2);   statBD = Mod(StatusVertex(GamD),2)

  !----------- step4 : weight calculation --------------------
  if(dir==1) then
    WA = weight_vertex(Order+1, 0, 0, 0, OmegaLn(GIC), omegaIA, TypeLn(GIC))
    WB = weight_vertex(Order+1, 0, 0, 0, OmegaLn(GMD), omegaMB, TypeLn(GMD))
    WIra = weight_vertex(Order+1, StatusVertex(Ira), diff_x(GXVertex(Ira),WXVertex(Ira)), &
      & diff_y(GYVertex(Ira),WYVertex(Ira)), omegaIA, OmegaLn(NeighVertex(2, Ira)), TypeVertex(Ira))
    WMasha = weight_vertex(Order+1, StatusVertex(Masha), diff_x(GXVertex(Masha),WXVertex(Masha)),  &
      & diff_y(GYVertex(Masha),WYVertex(Masha)), OmegaLn(NeighVertex(1,Masha)), &
      & omegaMB,  TypeVertex(Masha))
  else
    WA = weight_vertex(Order+1, 0, 0, 0, omegaIA, OmegaLn(GIC), TypeLn(GIC))
    WB = weight_vertex(Order+1, 0, 0, 0, omegaMB, OmegaLn(GMD), TypeLn(GMD))
    WIra = weight_vertex(Order+1, StatusVertex(Ira), diff_x(GXVertex(Ira),WXVertex(Ira)), &
      & diff_y(GYVertex(Ira),WYVertex(Ira)), OmegaLn(NeighVertex(1, Ira)), &
      & omegaIA, TypeVertex(Ira))
    WMasha = weight_vertex(Order+1, StatusVertex(Masha), diff_x(GXVertex(Masha),WXVertex(Masha)),  &
      & diff_y(GYVertex(Masha),WYVertex(Masha)), omegaMB, OmegaLn(NeighVertex(2, Masha)), &
      & TypeVertex(Masha))
  endif

  WWAB = weight_line(0, 2, diff_x(xwA,xwB), diff_y(ywA,ywB), omega, typAB)
  WGIA = weight_line(statIA, 1, 0, 0, omegaIA, TypeLn(GIC)) 
  WGAC = weight_line(statAC, 1, 0, 0, OmegaLn(GIC), TypeLn(GIC)) 
  WGMB = weight_line(statMB, 1, 0, 0, omegaMB, TypeLn(GMD)) 
  WGBD = weight_line(statBD, 1, 0, 0, OmegaLn(GMD), TypeLn(GMD)) 
  WWorm= weight_worm(diff_x(GXVertex(Ira),GXVertex(Masha)),diff_y(GYVertex(Ira),GYVertex(Masha)), &
    & diff_x(WXVertex(Ira),WXVertex(Masha)),diff_y(WYVertex(Ira),WYVertex(Masha)), omegaM)
    
  !----------  change the topology for the configuration after update --
  call insert_gamma(GamA, GXVertex(Ira), GYVertex(Ira), xwA, ywA, dirW, typGamA, 0, WA)
  call insert_gamma(GamB, GXVertex(Masha), GYVertex(Masha), xwB, ywB, 3-dirW, typGamB, 0, WB)
  call insert_line(WAB, omega, q, 2, typAB, 0, WWAB)
  call insert_line(GIA, omegaIA, kIA, 1, TypeLn(GIC), statIA, WGIA)
  call insert_line(GMB, omegaMB, kMB, 1, TypeLn(GMD), statMB, WGMB)

  !------------- update the topology -----------------------------
  NeighVertex(dir, GamA) = GIC;     NeighVertex(3-dir, GamA) = GIA
  NeighVertex(dir, GamB) = GMB;     NeighVertex(3-dir, GamB) = GMD
  NeighVertex(3, GamA)   = WAB;     NeighVertex(3, GamB)     = WAB

  NeighLn(dir, GIA)   = GamA;       NeighLn(3-dir, GIA)   = Ira
  NeighLn(dir, GMB)   = Masha;      NeighLn(3-dir, GMB)   = GamB

  NeighLn(dirW, WAB)  = GamA;       NeighLn(3-dirW, WAB)  = GamB

  NeighVertex(dir, Ira)  = GIA;     NeighVertex(3-dir, Masha)  = GMB
  NeighLn(3-dir, GIC) = GamA;       NeighLn(dir, GMD)   = GamB 


  !----------- step3 : configuration check --------------------
  flag = 0
  if(flag==0) then
    if(Is_reducible_G(GIA))       flag = 1
  endif
  if(flag==0) then
    if(Is_reducible_G(GMB))       flag = 2
  endif
  if(flag==0) then
    if(Is_reducible_W(WAB))       flag = 3
  endif
  if(flag==0) then
    if(Is_reducible_G_Gamma(GIA))       flag = 4
  endif
  if(flag==0) then
    if(Is_reducible_G_Gamma(GMB))       flag = 5
  endif
  if(flag==0) then
    if(Is_reducible_G_Gamma(GIC))       flag = 6
  endif
  if(flag==0) then
    if(Is_reducible_G_Gamma(GMD))       flag = 7
  endif
  if(flag==0) then
    if(Is_reducible_W_Gamma(WAB))       flag = 8
  endif

  if(flag >= 1) then
    !------------- delete line and vertexes --------------------
    call undo_insert_line(GIA,1)
    call undo_insert_line(GMB,1)
    call undo_insert_line(WAB,2)
    call delete_gamma(GamA)
    call delete_gamma(GamB)

    NeighVertex(dir, Ira)  = GIC;     NeighVertex(3-dir, Masha)  = GMD
    NeighLn(3-dir, GIC) = Ira;        NeighLn(dir, GMD)   = Masha

    return
  endif

  !------------ weight calculation ----------------------------
  if(MeasGamma/=Ira .and. MeasGamma/=Masha) then
    WMeasGamma = weight_vertex(Order+1, StatusVertex(MeasGamma), &
      & diff_x(GXVertex(MeasGamma), WXVertex(MeasGamma)), diff_y(GYVertex(MeasGamma), &
      & WYVertex(MeasGamma)), OmegaLn(NeighVertex(1, MeasGamma)), &
      & OmegaLn(NeighVertex(2, MeasGamma)), TypeVertex(MeasGamma)) 
    Anew = CoefOfWeight(Order+1)*WWorm *WMeasGamma *WA *WB *WWAB *WIra *WMasha *WGIA *WGAC &
      & *WGMB *WGBD*(1.d0)/Beta
    Aold = CoefOfWeight(Order)*WeightWorm *WeightVertex(MeasGamma) *WeightLn(GIC) *WeightLn(GMD)&
      & *WeightVertex(Ira) *WeightVertex(Masha)
  else
    Anew = CoefOfWeight(Order+1)*WWorm *WA *WB *WWAB *WIra *WMasha *WGIA *WGAC &
      & *WGMB *WGBD*(1.d0)/Beta
    Aold = CoefOfWeight(Order)*WeightWorm *WeightLn(GIC) *WeightLn(GMD)&
      & *WeightVertex(Ira) *WeightVertex(Masha)
  endif

  call weight_ratio(Pacc, sgn, Anew, Aold)

  Pacc = Pacc/prob_omega(omega)
  Pacc = Pacc *Pupdate(10)/Pupdate(9)

  !write(*, *) iupdate, Pacc
  !------------ step5 : accept the update -----------------------
  ProbProp(Order, iupdate) = ProbProp(Order, iupdate) + 1
  rand = rn()
  if(rand<=Pacc) then

    !-------------- update the diagram info --------------------
    Order = Order + 1
    Phase = Phase *sgn

    !------------- update k and omega -------------------------
    kMasha = kM
    OmegaMasha = OmegaM

    call add_Hash4G(kIA)
    call add_Hash4G(kMB)
    call add_Hash4W(q)

    !------------- update the status of elements --------------
    StatusLn(GMD) = statBD
    StatusLn(GIC) = statAC

    !------------- update weight of elements ------------------
    WeightWorm = WWorm
    WeightVertex(Ira) = WIra
    WeightVertex(Masha) = WMasha
    if(MeasGamma/=Ira .and. MeasGamma/=Masha)  WeightVertex(MeasGamma) = WMeasGamma
    WeightLn(GIC) = WGAC
    WeightLn(GMD) = WGBD

    call update_weight(Anew, Aold)

    ProbAcc(Order-1, 9) = ProbAcc(Order-1, 9) + 1
  else
    
    !------------- delete line and vertexes --------------------
    call undo_insert_line(GIA,1)
    call undo_insert_line(GMB,1)
    call undo_insert_line(WAB,2)
    call delete_gamma(GamA)
    call delete_gamma(GamB)

    NeighVertex(dir, Ira)  = GIC;    NeighVertex(3-dir, Masha)  = GMD
    NeighLn(3-dir, GIC) = Ira;       NeighLn(dir, GMD)   = Masha

  endif

END SUBROUTINE add_interaction_cross







!------------- remove interaction cross : Pupdate(10) -----------------------
!---------------------- dir = 1, dirW = 1 ------------------------
!------               GamC                     GamA     GamC -----
!------    Ira ----<----            Ira ----<---||---<-----  -----
!------                     <==                 \/           -----
!------  Masha ---->----          Masha ---->---||--->-----  -----
!------               GamD                     GamB     GamD -----
!-----------------------------------------------------------------
!---------------------- dir = 1, dirW = 2 ------------------------
!------               GamC                     GamA     GamC -----
!------    Ira ----<----            Ira ----<---||---<-----  -----
!------                     <==                 /\           -----
!------  Masha ---->----          Masha ---->---||--->-----  -----
!------               GamD                     GamB     GamD -----
!-----------------------------------------------------------------
!---------------------- dir = 2, dirW = 1 ------------------------
!------               GamC                     GamA     GamC -----
!------    Ira ---->----            Ira ---->---||--->-----  -----
!------                     <==                 \/           -----
!------  Masha ----<----          Masha ----<---||---<-----  -----
!------               GamD                     GamB     GamD -----
!-----------------------------------------------------------------
!---------------------- dir = 2, dirW = 2 ------------------------
!------               GamC                     GamA     GamC -----
!------    Ira ---->----            Ira ---->---||--->-----  -----
!------                     <==                 /\           -----
!------  Masha ----<----          Masha ----<---||---<-----  -----
!------               GamD                     GamB     GamD -----
!-----------------------------------------------------------------
!-----------------------------------------------------------------
SUBROUTINE remove_interaction_cross
  implicit none
  integer :: i, dir, dirW, omega
  integer :: omegaM, kM
  integer :: xwA, ywA, xwB, ywB
  integer :: WAB, GIA, GAC, GMB, GBD
  integer :: GamA, GamB, GamC, GamD
  integer :: statIC, statMD
  double precision :: WWorm, WGIC, WGMD, WIra, WMasha
  double precision :: WMeasGamma
  double precision :: Pacc, Anew, Aold, sgn
  double precision :: rand
  integer :: flag

  !------------ step1 : check if worm is present ------------------
  if(IsWormPresent .eqv. .false.)    return
  !ProbProp(iupdate) = ProbProp(iupdate) + 1
  if(Order <= 0) return          !!! just for 1st and higher order

  !------------ step2 : propose the new config ------------------
  dir  = Floor(rn()*2.d0)+1
  GIA  = NeighVertex(dir, Ira);          GMB  = NeighVertex(3-dir, Masha)
  GamA = NeighLn(dir, GIA);              GamB = NeighLn(3-dir, GMB)
  WAB  = NeighVertex(3, GamA)
  if(StatusLn(WAB)/=0)               return
  if(WAB /= NeighVertex(3, GamB))    return 
  if(TypeVertex(GamA)>2) return 
  if(TypeVertex(GamB)>2) return 

  GAC  = NeighVertex(dir, GamA);         GBD  = NeighVertex(3-dir, GamB)
  GamC = NeighLn(dir, GAC);              GamD = NeighLn(3-dir, GBD)
  if(GIA==GMB .or. GamA==GamB .or. GamC==GamD)  return

  dirW  = DirecVertex(GamA)
  xwA   = WXVertex(GamA);                ywA  = WYVertex(GamA)
  xwB   = WXVertex(GamB);                ywB  = WYVertex(GamB)
  if(xwA/=GXVertex(Ira))   return 
  if(ywA/=GYVertex(Ira))   return 
  if(xwB/=GXVertex(Masha))   return 
  if(ywB/=GYVertex(Masha))   return 

  omega = OmegaLn(WAB)
  if(Is_delta_omega_not_valid(omega))  return

  kM = add_k(kMasha, (-1)**dirW *kLn(WAB))
  if(Is_k_valid(kM) .eqv. .false.)    return

  OmegaM = OmegaMasha +(-1)**dirW *omega
  if(Is_omega_not_valid(OmegaM))   return

  if(Ira==MeasGamma .or. GamC==MeasGamma) then
    statIC = 1
  else
    statIC = 0
  endif
  if(Masha==MeasGamma .or. GamD==MeasGamma) then
    statMD = 1
  else
    statMD = 0
  endif

  !------------- change the topology -------------------------
  !------------- delete line and vertexes --------------------
  call delete_line(GIA,1)
  call delete_line(GMB,1)
  call delete_line(WAB,2)
  call delete_gamma(GamA)
  call delete_gamma(GamB)

  !------------- update the topology ------------------------
  NeighVertex(dir, Ira) = GAC;      NeighVertex(3-dir, Masha) = GBD
  NeighLn(3-dir, GAC) = Ira;        NeighLn(dir, GBD) = Masha

  !----------- step3 : configuration check ----------------------
  flag = 0
  if(flag==0) then
    if(NeighVertex(1, Ira)==NeighVertex(2, Ira))  flag = 1
  endif

  if(flag==0) then
    if(NeighVertex(1, Masha)==NeighVertex(2, Masha))  flag = 1
  endif

  if(flag == 1) then

    !----------  change the topology for the configuration after update --
    call undo_delete_line(WAB,2, 0)
    call undo_delete_line(GMB,1, Mod(StatusVertex(Masha),2))
    call undo_delete_line(GIA,1, Mod(StatusVertex(Ira),2))

    call undo_delete_gamma(GamB)
    call undo_delete_gamma(GamA)

    !------------- update the topology ------------------------
    NeighVertex(dir, Ira) = GIA;      NeighVertex(3-dir, Masha) = GMB
    NeighLn(3-dir, GAC) = GamA;       NeighLn(dir, GBD) = GamB
    return
  endif

  !----------- step4 : weight calculation -----------------------
  WGIC = weight_line(statIC, 1, 0, 0, OmegaLn(GAC), TypeLn(GAC))
  WGMD = weight_line(statMD, 1, 0, 0, OmegaLn(GBD), TypeLn(GBD))
  if(dir==1) then
    WIra = weight_vertex(Order-1, StatusVertex(Ira), diff_x(GXVertex(Ira),WXVertex(Ira)), &
      & diff_y(GYVertex(Ira),WYVertex(Ira)),  &
      & OmegaLn(GAC), OmegaLn(NeighVertex(2, Ira)), TypeVertex(Ira))
    WMasha = weight_vertex(Order-1, StatusVertex(Masha), diff_x(GXVertex(Masha),WXVertex(Masha)), &
      & diff_y(GYVertex(Masha),WYVertex(Masha)), OmegaLn(NeighVertex(1, Masha)),  &
      & OmegaLn(GBD), TypeVertex(Masha))
  else
    WIra = weight_vertex(Order-1, StatusVertex(Ira), diff_x(GXVertex(Ira),WXVertex(Ira)), &
      & diff_y(GYVertex(Ira),WYVertex(Ira)),OmegaLn(NeighVertex(1, Ira)), &
      & OmegaLn(GAC), TypeVertex(Ira))
    WMasha = weight_vertex(Order-1, StatusVertex(Masha), diff_x(GXVertex(Masha),WXVertex(Masha)), &
      & diff_y(GYVertex(Masha),WYVertex(Masha)), OmegaLn(GBD),  &
      & OmegaLn(NeighVertex(2,Masha)), TypeVertex(Masha))
  endif
  WWorm= weight_worm(diff_x(GXVertex(Ira),GXVertex(Masha)),diff_y(GYVertex(Ira),GYVertex(Masha)), &
    & diff_x(WXVertex(Ira),WXVertex(Masha)),diff_y(WYVertex(Ira),WYVertex(Masha)), omegaM)

  if(MeasGamma/=Ira .and. MeasGamma/=Masha) then
    WMeasGamma = weight_vertex(Order-1, StatusVertex(MeasGamma), &
      & diff_x(GXVertex(MeasGamma), WXVertex(MeasGamma)), diff_y(GYVertex(MeasGamma), &
      & WYVertex(MeasGamma)), OmegaLn(NeighVertex(1, MeasGamma)), &
      & OmegaLn(NeighVertex(2, MeasGamma)), TypeVertex(MeasGamma)) 

    Anew = CoefOfWeight(Order-1)*WWorm*WMeasGamma*WGIC *WGMD *WIra *WMasha
    Aold = CoefOfWeight(Order)*WeightWorm *WeightVertex(MeasGamma) *WeightLn(GIA) &
      & *WeightLn(GAC)*WeightLn(GMB)* WeightLn(GBD)*WeightVertex(Ira)*WeightVertex(Masha) &
      & *WeightVertex(GamA)*WeightVertex(GamB)*WeightLn(WAB)*(1.d0/Beta)
  else
    Anew = CoefOfWeight(Order-1)*WWorm*WGIC *WGMD *WIra *WMasha
    Aold = CoefOfWeight(Order)*WeightWorm *WeightLn(GIA) &
      & *WeightLn(GAC)*WeightLn(GMB)* WeightLn(GBD)*WeightVertex(Ira)*WeightVertex(Masha) &
      & *WeightVertex(GamA)*WeightVertex(GamB)*WeightLn(WAB)*(1.d0/Beta)
  endif

  call weight_ratio(Pacc, sgn, Anew, Aold)

  Pacc = Pacc *prob_omega(omega)
  Pacc = Pacc *Pupdate(9)/Pupdate(10)
  
  !write(*, *) iupdate, Pacc
  !------------ step5 : accept the update -----------------------
  ProbProp(Order, iupdate) = ProbProp(Order, iupdate) + 1
  rand = rn()
  if(rand<=Pacc) then

    !-------------- update the diagram info --------------------
    Order = Order - 1
    Phase = Phase *sgn

    !------------- update k and omega of elements ------------
    kMasha = kM
    OmegaMasha = OmegaM
    call delete_Hash4G(kLn(GIA))
    call delete_Hash4G(kLn(GMB))

    !------------ update the status of elements -------------
    StatusLn(GAC) = statIC
    StatusLn(GBD) = statMD

    !------------ update the weight of elements -------------
    WeightWorm = WWorm
    WeightVertex(Ira) = WIra
    WeightVertex(Masha) = WMasha
    if(MeasGamma/=Ira .and. MeasGamma/=Masha)  WeightVertex(MeasGamma)=WMeasGamma
    WeightLn(GAC) = WGIC
    WeightLn(GBD) = WGMD

    call update_weight(Anew, Aold)

    ProbAcc(Order+1, 10) = ProbAcc(Order+1, 10) + 1
  else

    !----------  change the topology for the configuration after update --
    call undo_delete_line(WAB,2, 0)
    call undo_delete_line(GMB,1, Mod(StatusVertex(Masha),2))
    call undo_delete_line(GIA,1, Mod(StatusVertex(Ira),2))

    call undo_delete_gamma(GamB)
    call undo_delete_gamma(GamA)

    !------------- update the topology ------------------------
    NeighVertex(dir, Ira) = GIA;      NeighVertex(3-dir, Masha) = GMB
    NeighLn(3-dir, GAC) = GamA;       NeighLn(dir, GBD) = GamB

  endif

END SUBROUTINE remove_interaction_cross




!------------- reconnect : Pupdate(11) -----------------------
!------------------------- dirW = 1 ------------------------------
!------  Ira           GamA            Ira            GamA   -----
!------    ------<-------                  --<-\    /-<--    -----
!------                                         \  /         -----
!------                           ==>            \/          -----
!------    ------<-------                  --<---/\---<---   -----
!------  Masha         GamB            Masha          GamB   -----
!-----------------------------------------------------------------
!------------------------- dirW = 2 ------------------------------
!------  Ira           GamA            Ira            GamA   -----
!------    ------>-------                  -->-\    /->--    -----
!------                                         \  /         -----
!------                           ==>            \/          -----
!------    ------>-------                  -->---/\--->---   -----
!------  Masha         GamB            Masha          GamB   -----
!-----------------------------------------------------------------
SUBROUTINE reconnect
  implicit none
  integer :: GamA, GamB, GIA, GMB, WLnIra, WLnMasha
  integer :: GIC, GMD
  integer :: dir, omegaM, kM
  integer :: statIA, statMB
  double precision :: WWorm, WGIA, WGMB, WIra, WMasha
  double precision :: Pacc, sgn, Anew, Aold
  integer :: flag
  
  !------- step1 : check if worm is present -------------
  If(IsWormPresent .eqv. .false.)  return
  !ProbProp(iupdate) = ProbProp(iupdate) + 1
  !If(Order<=1)     return

  !------- step2 : propose a new config -----------------
  if((GXVertex(Ira)/=GXVertex(Masha) .or. GYVertex(Ira)/=GYVertex(Masha)))  return

  dir = Floor(rn()*2.d0)+1
  GIA = NeighVertex(dir, Ira);       GMB  = NeighVertex(dir, Masha)
  if(TypeLn(GIA)/=TypeLn(GMB))            return

  GIC = NeighVertex(3-dir, Ira);     GMD  = NeighVertex(3-dir, Masha)
  WLnIra = NeighVertex(3, Ira);      WLnMasha = NeighVertex(3, Masha)
  GamA = NeighLn(dir, GIA);       GamB = NeighLn(dir, GMB)

  omegaM = OmegaMasha + (-1)**dir*(OmegaLn(GMB) -OmegaLn(GIA))
  kM = add_k(kMasha, (-1)**dir*(add_k(kLn(GMB), -kLn(GIA))))
  if(Is_omega_not_valid(omegaM)) return
  if(Is_k_valid(kM) .eqv. .false.) return


  !------ update the topology ---------------
  NeighVertex(dir, Ira)  = GMB
  NeighLn(3-dir, GMB) = Ira

  NeighVertex(dir, Masha) = GIA
  NeighLn(3-dir, GIA) = Masha 

  !----- the status for the new config ------------------
  if(GamA==MeasGamma .or. Masha==MeasGamma) then
    statIA = 1
  else
    statIA = 0
  endif
  if(GamB==MeasGamma .or. Ira==MeasGamma) then
    statMB = 1
  else
    statMB = 0
  endif


  !------- step3 : configuration check ------------------
  flag = 0

  if(flag==0) then
    if(Is_reducible_G_Gamma(GIA))  flag = 1
  endif
  if(flag==0) then
    if(Is_reducible_G_Gamma(GMB))  flag = 1
  endif

  if(flag==0) then
    if(NeighVertex(1, Ira)==NeighVertex(2, Ira)) flag=1
  endif
  if(flag==0) then
    if(NeighVertex(1, Masha)==NeighVertex(2, Masha)) flag=1
  endif

  if(flag == 1) then
    NeighVertex(dir, Ira)  = GIA
    NeighLn(3-dir, GIA) = Ira

    NeighVertex(dir, Masha) = GMB
    NeighLn(3-dir, GMB) = Masha

    return
  endif


  !------- step4 : weight calculation -------------------
  WGIA = weight_line(statIA, 1, 0, 0, OmegaLn(GIA), TypeLn(GIA))
  WGMB = weight_line(statMB, 1, 0, 0, OmegaLn(GMB), TypeLn(GMB))

  if(dir==1) then
    WIra = weight_vertex(Order, StatusVertex(Ira), diff_x(GXVertex(Ira),WXVertex(Ira)), &
      diff_y(GYVertex(Ira),WYVertex(Ira)), OmegaLn(GMB), OmegaLn(GIC), TypeVertex(Ira))
    WMasha = weight_vertex(Order, StatusVertex(Masha), diff_x(GXVertex(Masha),WXVertex(Masha)), &
      diff_y(GYVertex(Masha),WYVertex(Masha)), OmegaLn(GIA), OmegaLn(GMD), TypeVertex(Masha))
  else
    WIra = weight_vertex(Order, StatusVertex(Ira), diff_x(GXVertex(Ira),WXVertex(Ira)), &
      diff_y(GYVertex(Ira),WYVertex(Ira)), OmegaLn(GIC), OmegaLn(GMB), TypeVertex(Ira))
    WMasha = weight_vertex(Order, StatusVertex(Masha), diff_x(GXVertex(Masha),WXVertex(Masha)), &
      diff_y(GYVertex(Masha),WYVertex(Masha)), OmegaLn(GMD), OmegaLn(GIA), TypeVertex(Masha))
  endif
  WWorm= weight_worm(diff_x(GXVertex(Ira),GXVertex(Masha)),diff_y(GYVertex(Ira),GYVertex(Masha)), &
    & diff_x(WXVertex(Ira),WXVertex(Masha)),diff_y(WYVertex(Ira),WYVertex(Masha)), omegaM)

  Anew = (-1.d0) *WWorm *WIra *WMasha *WGIA *WGMB
  Aold = WeightWorm *WeightVertex(Ira) *WeightVertex(Masha) *WeightLn(GIA) *WeightLn(GMB)

  call weight_ratio(Pacc, sgn, Anew, Aold)

  !------- step5 : accept the update --------------------
  ProbProp(Order, iupdate) = ProbProp(Order, iupdate) + 1
  if(rn()<=Pacc) then

    !------ update diagram info ---------------
    NFermiLoop = NFermiLoop+1      ! 1 for odd; 2 for even
    Phase = Phase *sgn

    !------ update Ira and Masha --------------
    OmegaMasha = omegaM
    kMasha = kM

    !------ update the status -----------------
    StatusLn(GIA) = statIA
    StatusLn(GMB) = statMB

    !------ update weight of elements ---------
    WeightWorm = WWorm
    WeightLn(GIA) = WGIA
    WeightLn(GMB) = WGMB
    WeightVertex(Ira) = WIra
    WeightVertex(Masha) = WMasha

    call update_weight(Anew, Aold)

    ProbAcc(Order, 11) = ProbAcc(Order, 11) + 1

    !call print_config
  else

    !------ update the topology ---------------
    NeighVertex(dir, Ira)  = GIA
    NeighLn(3-dir, GIA) = Ira

    NeighVertex(dir, Masha) = GMB
    NeighLn(3-dir, GMB) = Masha

  endif
END SUBROUTINE reconnect







!------------ shift gline in space: Pupdate(12) ------------
SUBROUTINE shift_gline_in_space
  implicit none
  integer :: iGam, Gam(MxNVertex), i, dxg, dyg, xg, yg, xw, yw, iG
  integer :: dir, iWLn, WLn(MxNVertex), nGam, flag, j, xwj, ywj
  double precision :: WGam(MxNVertex), WW(MxNVertex), Anew, Aold, Pacc, sgn

  !------- step1 : check if worm is present -------------
  if(IsWormPresent .eqv. .true.)    return
  !ProbProp(iupdate) = ProbProp(iupdate) + 1

  !------- step2 : propose a new config -----------------
  Gam(:) = -1
  WLn(:) = -1
  i=1
  iGam = generate_gamma()
  Gam(1) = iGam
  WLn(1) = NeighVertex(3, iGam)
  iG = NeighVertex(2, iGam)
  do while(NeighLn(2, iG)/=Gam(1))
    i = i+1
    iGam = NeighLn(2, iG) 
    iWLn = NeighVertex(3, iGam)
    Gam(i) = iGam
    WLn(i) = iWLn
    iG = NeighVertex(2, iGam)
  enddo
  nGam = i

  dxg = generate_x();             dyg = generate_y()
  xg = find_neigh_x(GXVertex(Gam(1)), dxg)
  yg = find_neigh_y(GYVertex(Gam(1)), dyg)

  !------- step4 : weight calculation -------------------
  Anew = 1.d0;                   Aold = 1.d0
  do i = 1, MxNVertex
    if(Gam(i)==-1) cycle

    xw = find_neigh_x(WXVertex(Gam(i)), dxg)
    yw = find_neigh_y(WYVertex(Gam(i)), dyg)
    dir = DirecVertex(Gam(i))

    WGam(i) = weight_vertex(Order, StatusVertex(Gam(i)), diff_x(xg,xw),  &
      & diff_y(yg,yw), OmegaLn(NeighVertex(1, &
      & Gam(i))), OmegaLn(NeighVertex(2, Gam(i))), TypeVertex(Gam(i)))

    xwj = WXVertex(NeighLn(3-dir, WLn(i)))
    ywj = WYVertex(NeighLn(3-dir, WLn(i)))

    do j = 1, MxNVertex
      if(Gam(j)==-1) cycle
      if(NeighLn(3-dir, WLn(i))==Gam(j)) then
        xwj = find_neigh_x(xwj, dxg)
        ywj = find_neigh_y(ywj, dyg)
      endif
    enddo

    WW(i) = weight_line(StatusLn(WLn(i)), 2, diff_x(xw,xwj), diff_y(yw,ywj), &
      & OmegaLn(WLn(i)), TypeLn(WLn(i)))

    Anew = Anew *WGam(i) *WW(i)
    Aold = Aold *WeightVertex(Gam(i)) *WeightLn(WLn(i))

    do j = 1, MxNVertex
      if(Gam(j)==-1) cycle
      if(NeighLn(3-dir, WLn(i))==Gam(j)) then
        Anew = Anew/sqrt(abs(WW(i)))
        Aold = Aold/sqrt(abs(WeightLn(WLn(i))))
      endif
    enddo
  enddo

  call weight_ratio(Pacc, sgn, Anew, Aold)

  !------- step5 : accept the update --------------------
  ProbProp(Order, iupdate) = ProbProp(Order, iupdate) + 1
  if(rn()<=Pacc) then

    !------ update the diagram info -----------------
    Phase = Phase *sgn

    do i = 1, MxNVertex
      if(Gam(i)==-1) cycle
      !------ update the site of elements -------------
      GXVertex(Gam(i)) = xg
      GYVertex(Gam(i)) = yg
      WXVertex(Gam(i)) = find_neigh_x(WXVertex(Gam(i)), dxg)
      WYVertex(Gam(i)) = find_neigh_y(WYVertex(Gam(i)), dyg)

      !------ update the weight of elements -------------
      WeightVertex(Gam(i)) = WGam(i)
      WeightLn(WLn(i))  = WW(i)
    enddo

    call update_weight(Anew, Aold)

    ProbAcc(Order, 12) = ProbAcc(Order, 12) + 1
  endif
END SUBROUTINE shift_gline_in_space




!------------ shift wline in space: Pupdate(13) ------------
SUBROUTINE shift_wline_in_space
  implicit none
  integer :: iGam, iWLn, jGam, dxw, dyw, xw, yw, xwj, ywj, dir
  double precision :: WGam, WW, Anew, Aold, Pacc, sgn

  !------- step1 : check if worm is present -------------
  if(IsWormPresent .eqv. .true.)    return
  !ProbProp(iupdate) = ProbProp(iupdate) + 1

  !------- step2 : propose a new config -----------------
  iGam = generate_gamma()
  dir  = DirecVertex(iGam)
  iWLn = NeighVertex(3, iGam)
  jGam = NeighLn(3-dir, iWLn)

  dxw = generate_x();         dyw = generate_y()
  xw = find_neigh_x(WXVertex(iGam), dxw)
  yw = find_neigh_y(WYVertex(iGam), dyw)

  xwj = WXVertex(jGam);          ywj = WYVertex(jGam)
  
  !------- step3 : configuration check ------------------

  !------- step4 : weight calculation -------------------
  WGam = weight_vertex(Order, StatusVertex(iGam), diff_x(GXVertex(iGam),xw), diff_y(GYVertex(iGam),yw), &
    & OmegaLn(NeighVertex(1, iGam)), OmegaLn(NeighVertex(2, iGam)), TypeVertex(iGam))
  WW = weight_line(StatusLn(iWLn), 2, diff_x(xw,xwj), diff_y(yw,ywj), &
    & OmegaLn(iWLn), TypeLn(iWLn)) 


  Anew = WGam *WW
  Aold = WeightLn(iWLn)*WeightVertex(iGam)

  call weight_ratio(Pacc, sgn, Anew, Aold)

  !Pacc = Pacc*prob_x(WXVertex(iGam))*prob_y(WYVertex(iGam))/(prob_x(dxw)*prob_y(dyw))

  !------- step5 : accept the update --------------------
  ProbProp(Order, iupdate) = ProbProp(Order, iupdate) + 1
  if(rn()<=Pacc) then

    !------ update the diagram info -------------------
    Phase = Phase *sgn

    !------ update the site of elements --------------
    WXVertex(iGam) = xw
    WYVertex(iGam) = yw

    !------ update the weight of elements ------------
    WeightLn(iWLn) = WW
    WeightVertex(iGam) = WGam

    call update_weight(Anew, Aold)

    ProbAcc(Order, 13) = ProbAcc(Order, 13) + 1
  endif
  return
END SUBROUTINE shift_wline_in_space
!-----------------------------------------------------------


!------------ change Gamma type: Pupdate(14) ------------
SUBROUTINE change_Gamma_type
  implicit none
  integer :: iGam, iWLn, jGam, dir, typGam, typW
  double precision :: WGam, WW, Anew, Aold, Pacc, sgn

  !------- step1 : check if worm is present -------------
  if(IsWormPresent .eqv. .true.)    return
  !ProbProp(iupdate) = ProbProp(iupdate) + 1

  !------- step2 : propose a new config -----------------
  iGam = generate_gamma()
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

  !------- step3 : configuration check ------------------

  !------- step4 : weight calculation -------------------
  WGam = weight_vertex(Order, StatusVertex(iGam), diff_x(GXVertex(iGam),WXVertex(iGam)), &
    & diff_y(GYVertex(iGam),WYVertex(iGam)), &
    & OmegaLn(NeighVertex(1, iGam)), OmegaLn(NeighVertex(2, iGam)), typGam)


  WW = weight_line(StatusLn(iWLn), 2, diff_x(WXVertex(iGam),WXVertex(jGam)), &
    & diff_y(WYVertex(iGam),WYVertex(jGam)), OmegaLn(iWLn), typW) 

  Anew = WGam *WW
  Aold = WeightLn(iWLn)*WeightVertex(iGam)

  call weight_ratio(Pacc, sgn, Anew, Aold)

  !------- step5 : accept the update --------------------
  ProbProp(Order, iupdate) = ProbProp(Order, iupdate) + 1
  if(rn()<=Pacc) then

    !------ update the diagram info -------------------
    Phase = Phase *sgn

    !------ update the site of elements --------------
    TypeVertex(iGam)    = typGam
    TypeVertexIn(iGam)  = 3-TypeVertexIn(iGam)
    TypeVertexOut(iGam) = 3-TypeVertexOut(iGam)

    TypeLn(iWLn)     = typW

    !------ update the weight of elements ------------
    WeightLn(iWLn) = WW
    WeightVertex(iGam) = WGam

    call update_weight(Anew, Aold)

    ProbAcc(Order, 14) = ProbAcc(Order, 14) + 1
  endif
  return
END SUBROUTINE change_Gamma_type
!-----------------------------------------------------------



!------------ move measuring index: Pupdate(15) ------------
SUBROUTINE move_measuring_index
  implicit none
  integer :: iW, iGin, iGout, iGam
  integer :: jW, jGin, jGout, jGam
  integer :: statiGam, statiGin, statiGout, statiW
  integer :: statjGam, statjGin, statjGout, statjW
  double precision :: WiGam, WiGin, WiGout, WiW
  double precision :: WjGam, WjGin, WjGout, WjW, Anew, Aold, Pacc, sgn

  !------- step1 : check if worm is present -------------
  if(IsWormPresent .eqv. .true.)    return
  !ProbProp(iupdate) = ProbProp(iupdate) + 1

  !------- step2 : propose a new config -----------------
  iGam = generate_gamma()

  iGin  = NeighVertex(1, iGam);      iGout = NeighVertex(2, iGam)
  iW    = NeighVertex(3, iGam)

  jGam  = MeasGamma
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
  WiGam = weight_vertex(Order, statiGam, diff_x(GXVertex(iGam),WXVertex(iGam)), &
    & diff_y(GYVertex(iGam),WYVertex(iGam)), OmegaLn(iGin), OmegaLn(iGout), TypeVertex(iGam))

  WiW = weight_line(statiW, 2, diff_x(WXVertex(NeighLn(2, iW)),WXVertex(NeighLn(1,iW))),&
    & diff_y(WYVertex(NeighLn(2,iW)),WYVertex(NeighLn(1,iW))), OmegaLn(iW), TypeLn(iW))

  WiGin = weight_line(statiGin, 1, 0, 0, OmegaLn(iGin), TypeLn(iGin))
  WiGout = weight_line(statiGout, 1, 0, 0, OmegaLn(iGout), TypeLn(iGout))


  WjGam = weight_vertex(Order, statjGam, diff_x(GXVertex(jGam),WXVertex(jGam)), diff_y(GYVertex(jGam) &
    & ,WYVertex(jGam)), OmegaLn(jGin), OmegaLn(jGout), TypeVertex(jGam))

  WjW = weight_line(statjW, 2, diff_x(WXVertex(NeighLn(2, jW)),WXVertex(NeighLn(1,jW))),&
    & diff_y(WYVertex(NeighLn(2,jW)),WYVertex(NeighLn(1,jW))), OmegaLn(jW), TypeLn(jW))

  WjGin = weight_line(statjGin, 1, 0, 0, OmegaLn(jGin), TypeLn(jGin))
  WjGout = weight_line(statjGout, 1, 0, 0, OmegaLn(jGout), TypeLn(jGout))

  Anew = WjGam *WjW *WjGin *WjGout *WiGam *WiW *WiGin *WiGout
  Aold = WeightVertex(iGam) *WeightLn(iW) *WeightLn(iGin) *WeightLn(iGout) &
    & *WeightVertex(jGam) *WeightLn(jW) *WeightLn(jGin) *WeightLn(jGout)

  if(Anew==0.d0)   return

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
    MeasGamma = iGam

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

!!=======================================================================
!!=======================================================================
!!=======================================================================





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


LOGICAL FUNCTION Is_reducible_G_Gamma(GLn)
  implicit none
  integer, intent(in) :: GLn
  integer :: nG, Gam1, Gam2, W1, W2, nW
  integer :: newk, kG 
  integer :: i, nnk1, nnk2

  Is_reducible_G_Gamma = .false.
  newk = kLn(GLn)

  Gam1 = NeighLn(1, GLn)
  Gam2 = NeighLn(2, GLn)
  W1 = NeighVertex(3, Gam1)
  W2 = NeighVertex(3, Gam2)

  if(CheckGamma) then
    do i = 1, NGLn
      nG = Ln4GList(i)
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
        Is_reducible_G_Gamma = .true.
      endif
    enddo
  endif
    
  return
END FUNCTION Is_reducible_G_Gamma




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
LOGICAL FUNCTION Is_reducible_W_Gamma(WLn)
  implicit none
  integer, intent(in) :: WLn
  integer :: absk, i, Gam1, Gam2, G1, G2, G3, G4, kG5, kG6
  integer :: nG, kG, pkGG, nkGG

  absk = abs(kLn(WLn))
  Is_reducible_W_Gamma = .false.
  Gam1 = NeighLn(1, WLn)
  Gam2 = NeighLn(2, WLn)
  G1 = NeighVertex(1, Gam1)
  G2 = NeighVertex(2, Gam1)
  G3 = NeighVertex(1, Gam2)
  G4 = NeighVertex(2, Gam2)

  if(CheckGamma) then
    do i = 1, NGLn
      nG = Ln4GList(i)
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
        Is_reducible_W_Gamma = .true.
      endif
    enddo
  endif
  return
END FUNCTION Is_reducible_W_Gamma

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
  if(delete_ira_stat == -2) delete_ira_stat = 0
  if(delete_ira_stat == -1) delete_ira_stat = 1
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
SUBROUTINE insert_line(newline, omega, k, knd, typ, stat, weigh)
  implicit none
  integer, intent(out) :: newline
  integer, intent(in) :: omega, k, knd, typ, stat
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
    Ln4GList(NGLn) = newline
    List4Ln(newline) = NGLn
  else
    NWLn = NWLn + 1
    Ln4WList(NWLn) = newline
    List4Ln(newline) = NWLn
  endif

  OmegaLn(newline) = omega
  kLn(newline) = k
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

    tmp = Ln4GList(NGLn)
    Ln4GList(NGLn) = 0
    Ln4GList(List4Ln(occline)) = tmp
    List4Ln(tmp) = List4Ln(occline)
    NGLn = NGLn -1
  else

    tmp = Ln4WList(NWLn)
    Ln4WList(NWLn) = 0
    Ln4WList(List4Ln(occline)) = tmp
    List4Ln(tmp) = List4Ln(occline)
    NWLn = NWLn -1
  endif

  if(TailLn == -1) then
    write(*, *) "Tail=-1! Too many lines!"
    stop
  endif
  return
END SUBROUTINE undo_insert_line


!------------- insert a gamma to the link -------------------
SUBROUTINE insert_gamma(newgamma, gx, gy, wx, wy, dir, typ, stat, weigh)
  implicit none
  integer, intent(out) :: newgamma
  integer, intent(in) :: gx, gy, wx, wy, dir, typ, stat
  double precision, intent(in) :: weigh

  newgamma = TailGam
  TailGam = NextVertex(TailGam)
  if(StatusVertex(TailGam)>=0) then
    write(*, *) IsWormPresent, iupdate, "insert_gamma error!!!"
    call print_config
    stop
  endif

  if(TailGam == -1) then
    write(*, *) "Tail=-1! Too many gammas!"
    stop
  endif
   
  NGam = NGam + 1
  Vertex4GamList(NGam) = newgamma
  List4Vertex(newgamma) = NGam

  GXVertex(newgamma) = gx
  GYVertex(newgamma) = gy
  WXVertex(newgamma) = wx
  WYVertex(newgamma) = wy
  DirecVertex(newgamma) = dir
  TypeVertex(newgamma) = typ
  TypeVertexIn(newgamma) = typ
  TypeVertexOut(newgamma) = typ
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

    tmp = Ln4GList(NGLn)
    Ln4GList(NGLn) = 0
    Ln4GList(List4Ln(occline)) = tmp
    List4Ln(tmp) = List4Ln(occline)
    NGLn = NGLn -1
  else

    tmp = Ln4WList(NWLn)
    Ln4WList(NWLn) = 0
    Ln4WList(List4Ln(occline)) = tmp
    List4Ln(tmp) = List4Ln(occline)
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
    Ln4GList(NGLn) = newline
    List4Ln(newline) = NGLn
    !call add_Hash4G(kLn(newline))
  else
    NWLn = NWLn + 1
    Ln4WList(NWLn) = newline
    List4Ln(newline) = NWLn
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
  NextVertex(occgamma) = TailGam
  StatusVertex(occgamma) = -1
  TailGam = occgamma

  tmp = Vertex4GamList(NGam)
  Vertex4GamList(NGam) = 0
  Vertex4GamList(List4Vertex(occgamma)) = tmp
  List4Vertex(tmp) = List4Vertex(occgamma)
  NGam = NGam -1

  if(TailGam == -1) then
    write(*, *) "Tail=-1! Too many vertexes!"
    stop
  endif
  return
END SUBROUTINE delete_gamma


SUBROUTINE undo_delete_gamma(newgamma)
  implicit none
  integer, intent(in) :: newgamma

  if(TailGam/=newgamma)    then
    write(*, *) "undo_delete_gamma error!"
    stop
  endif
  StatusVertex(newgamma) = 0

  TailGam = NextVertex(newgamma)
  if(TailGam == -1) then
    write(*, *) "Tail=-1! Too many gammas!"
    stop
  endif
   
  NGam = NGam + 1
  Vertex4GamList(NGam) = newgamma
  List4Vertex(newgamma) = NGam

  return
END SUBROUTINE undo_delete_gamma
!!=======================================================================
!!=======================================================================
!!=======================================================================






!!=======================================================================
!!================= DIAGRAM WEIGHT CALCULATION ==========================
!!=======================================================================

!-------- the weight of a vertex -------------------------
!dx = xg-xw;  dy = yg-yw
DOUBLE PRECISION FUNCTION weight_vertex(iorder, stat, dx, dy, omega1, omega2, typ)
  implicit none
  integer :: stat, dx, dy, omega1, omega2, typ, iorder
  double precision :: weight

  !---------------------- for test --------------------------------------
  !if(stat>=0 .and. stat<=1) then
    !weight_vertex = weight_meas_gamma(dx, dy, omega1, omega2, typ)
  !else if(stat>=2 .and. stat<=3) then
    !weight_vertex = weight_meas_gamma(dx, dy, omega1, omega2, 1)
  !else if(stat==-1) then
    !write(*, *) IsWormPresent, iupdate, "vertex status == -1! There is no weight!" 
    !stop
  !else
    !write(*, *) IsWormPresent, iupdate, "vertex status error!", stat
    !stop
  !endif
  !------------------------ end -----------------------------------------
  if(stat==0) then
    weight_vertex = weight_Gamma0(dx, dy, omega1, omega2, typ)
    !if(iorder>=1) then
      !weight_vertex = weight_gamma(dx, dy, omega1, omega2, typ)
    !else if(iorder==0) then
      !weight_vertex = weight_gamma0(dx, dy, omega1, omega2, typ)
    !endif
  else if(stat==2) then
    weight_vertex = weight_worm_gamma(dx, dy, omega1, omega2, typ)
  else if(stat==1 .or. stat==3) then
    weight_vertex = weight_meas_gamma(iorder, dx, dy, omega1, omega2, typ)
  else if(stat==-1) then 
    write(*, *) IsWormPresent, iupdate, "vertex status == -1! There is no weight!" 
    stop
  else
    write(*, *) IsWormPresent, iupdate, "vertex status error!", stat
    stop
  endif
  return
END FUNCTION weight_vertex

!-------- the weight of a line -------------------------
! dx = x2-x1;  dy = y2- y1
DOUBLE PRECISION FUNCTION weight_line(stat, knd, dx, dy, omega, typ)
  implicit none
  integer :: stat, knd, dx, dy, omega, typ

  !---------------------- for test --------------------------------------
  !if(stat >= 0 .and. stat<=1) then
    !if(knd==1) weight_line = weight_meas_G(omega, typ)
    !if(knd==2) weight_line = weight_meas_W(dx, dy, omega, typ)
  !else if(stat >= 2 .and. stat<=3) then
    !if(knd==1) weight_line = weight_meas_G(omega, 1)
    !if(knd==2) weight_line = weight_meas_W(dx, dy, omega, 1)
  !else if(stat==-1) then
    !write(*, *) IsWormPresent, iupdate, "line status == -1! There is no weight!" 
    !stop
  !else
    !write(*, *) IsWormPresent, iupdate, "line status error!", stat
    !stop
  !endif
  !------------------------ end -----------------------------------------

  if(stat == 0) then
    if(knd==1) weight_line = weight_G(omega, typ)
    if(knd==2) weight_line = weight_W(dx, dy, omega, typ)
  else if(stat == 1 .or. stat==3) then
    if(knd==1) weight_line = weight_meas_G(omega, typ)
    if(knd==2) weight_line = weight_meas_W(dx, dy, omega, typ)
  else if(stat == 2) then
    if(knd==2) weight_line = weight_worm_W(dx, dy, omega, typ)
    if(knd==1) then
      write(*, *) IsWormPresent, iupdate, "gline status == 2 or 3! Error!" 
      call print_config
      stop
    endif
  else if(stat==-1) then
    write(*, *) IsWormPresent, iupdate, "line status == -1! There is no weight!" 
    call print_config
    stop
  else
    write(*, *) IsWormPresent, iupdate, "line status error!", stat
    call print_config
    stop
  endif

  return
END FUNCTION weight_line

!-------- the weight ratio of new/old config -------------------------
SUBROUTINE weight_ratio(Pacc, sgn, Anew, Aold)
  implicit none
  double precision,intent(in) ::Anew, Aold
  double precision,intent(out) :: Pacc, sgn
  Pacc = Anew/Aold
  if(Pacc>1.d-20) then
    sgn = 1.d0
  else if(Pacc<-1.d-20) then
    sgn = -1.d0
    Pacc = -1.d0*Pacc
  else 
    sgn = 1.d0
    Pacc = 0.d0
  endif
  return
END SUBROUTINE weight_ratio

SUBROUTINE update_weight(Anew, Aold)
  implicit none 
  double precision :: Anew, Aold

  WeightCurrent = WeightCurrent *Anew/Aold
  return
END SUBROUTINE update_weight
!!=======================================================================
!!=======================================================================
!!=======================================================================










!!=======================================================================
!!============================ PRINT CONFIGURATION ======================
!!=======================================================================
SUBROUTINE print_updates
  implicit none
  write(*, *) "================================================="
  write(*, *) "IsWormPresent == .true.:" 
  write(*, *) "1: delete worm"
  write(*, *) "2: move worm along gline"
  write(*, *) "3: move worm along wline"
  write(*, *) "4: remove bubble"
  write(*, *) "5: add interaction"
  write(*, *) "6: remove interaction"
  write(*, *) "7: reconnect"

  write(*, *) "IsWormPresent == .false.:" 
  write(*, *) "1: create worm"
  write(*, *) "2: shift gline in space"
  write(*, *) "3: shift wline in space"
  write(*, *) "4: add bubble"
  write(*, *) "5: move measuring index"
  write(*, *) "6: change Gamma type"

  write(*, *) "================================================="
  return
END SUBROUTINE print_updates


SUBROUTINE print_config
  implicit none
  integer :: i, iln, iv
  
  open(8, access='append', file=trim(title3)//"_mc.conf")
  
  write(8, *) "============================================================"
  write(8, *) imc, IsWormPresent, iupdate

  if(IsWormPresent .eqv. .true.) then
    write(8, *) "Ira", Ira, "Masha", Masha, "SpinMasha", SpinMasha, "OmegaMasha",OmegaMasha
    write(8, *) "kMasha", kMasha
  endif

  write(8, *) "Order", Order
  write(8, *) "NFermiLoop", NFermiLoop

  write(8, *) "Measuring Gamma", MeasGamma
  write(8, *) "Phase", Phase
  write(8, *) "Weight", WeightCurrent

  do i = 1, NGLn
    iln = Ln4GList(i)
    if(StatusLn(iln) <0) cycle
    write(8, 10) iln, KindLn(iln), TypeLn(iln), kLn(iln), OmegaLn(iln), StatusLn(iln), NeighLn(1:2,iln)
  enddo

  do i = 1, NWLn
    iln = Ln4WList(i)
    if(StatusLn(iln) <0) cycle
    write(8, 10) iln, KindLn(iln), TypeLn(iln), kLn(iln), OmegaLn(iln), StatusLn(iln), NeighLn(1:2,iln)
  enddo

  do i = 1, NGam
    iv = Vertex4GamList(i)
    if(StatusVertex(iv) <0) cycle
    write(8, 12) iv,TypeVertex(iv),TypeVertexIn(iv),TypeVertexOut(iv),GXVertex(iv),GYVertex(iv),&
      & WXVertex(iv),WYVertex(iv), DirecVertex(iv), StatusVertex(iv), NeighVertex(:,iv)
  enddo
  write(8, *) "============================================================"

  10 format(' Line:',i2,2x,'kind:',i2,2x,'type:',i2,2x,'k:',i8,2x,'omega:',i6,2x,'stat:',i2, 2x,&
    & 'neigh:',i6,i6)
  12 format('Gamma:',i2,2x,'type:',i2,2x,'typein:',i2,2x,'typeout:',i2,2x,'gr:(',i4,i4,'), wr:(',i4,i4,')',&
  & 'direction:', i2,2x, 'stat:',i2, 2x,'neigh:', i6,i6,i6)

  close(8)
END SUBROUTINE print_config


INTEGER FUNCTION find_config()
  implicit none
  integer :: iconf, itopo
  if(Order==0) then
    if(IsWormPresent) then
      iconf = 2
    else
      iconf = 1
    endif
  else if(Order == 1) then
    if(IsWormPresent) then
      if(Ira==MeasGamma) then
        itopo = 0
      else 
        itopo = 1
      endif
      iconf = 4 + itopo
    else
      iconf = 3 
    endif
  else
    write(*, *) "Order/=0 or 1"
    stop
  endif
  find_config = iconf
END FUNCTION find_config
!!=======================================================================
!!=======================================================================
!!=======================================================================


!!=======================================================================
!!========================= FIND TOPOLOGY FOR ORDER 2 ===================
!!=======================================================================

LOGICAL FUNCTION is_topo_1()
  implicit none
  integer :: i, iGLn
  integer :: MeaGin, MeaGout, MeaW
  integer :: G1, G2, G3, G4, W1, W2, Gam1, Gam2, Gam3, Gam4

  MeaGin = NeighVertex(2, MeasGamma)
  MeaGout = NeighVertex(1, MeasGamma)
  MeaW = NeighVertex(3, MeasGamma)

  if(Order/=2) then
    is_topo_1 = .false.
    return
  endif

  if(OmegaLn(MeaGin)/=OmegaLn(MeaGout)) then
    is_topo_1 = .false.
    return
  else if(OmegaLn(MeaW)/=0) then
    write(*, *) "Omega Error when find topo 1"
    call print_config
    stop
  endif

  if(Mod(NFermiloop,2)==0) then
    is_topo_1 = .false.
    return
  endif

  G1 = NeighVertex(2, MeasGamma)
  Gam1 = NeighLn(2, G1)
  G2 = NeighVertex(2, Gam1)
  Gam2 = NeighLn(2, G2)
  W1 = NeighVertex(3, Gam1)
  W2 = NeighVertex(3, Gam2)
  Gam3 = NeighLn(3-DirecVertex(Gam1), W1)
  Gam4 = NeighLn(3-DirecVertex(Gam2), W2)
  if(NeighVertex(1, Gam3)==NeighVertex(2, Gam4)) then
    is_topo_1 = .true.
  else if(NeighVertex(2, Gam3)==NeighVertex(1, Gam4)) then
    is_topo_1 = .false.
  else
    write(*, *) "Error when find topo 1"
    call print_config
    stop
  endif
  return
END FUNCTION is_topo_1


LOGICAL FUNCTION is_topo_2()
  implicit none
  integer :: i, iGLn
  integer :: MeaGin, MeaGout, MeaW
  integer :: G1, G2, G3, G4, W1, W2, Gam1, Gam2, Gam3, Gam4

  MeaGin = NeighVertex(2, MeasGamma)
  MeaGout = NeighVertex(1, MeasGamma)
  MeaW = NeighVertex(3, MeasGamma)

  if(Order/=2) then
    is_topo_2 = .false.
    return
  endif

  if(OmegaLn(MeaGin)/=OmegaLn(MeaGout)) then
    is_topo_2 = .false.
    return
  else if(OmegaLn(MeaW)/=0) then
    write(*, *) "Omega Error when find topo 2"
    call print_config
    stop
  endif

  if(Mod(NFermiloop,2)==1) then
    is_topo_2 = .false.
    return
  endif

  G1 = NeighVertex(2, MeasGamma)
  Gam1 = NeighLn(2, G1)
  G2 = NeighVertex(2, Gam1)
  Gam2 = NeighLn(2, G2)
  W1 = NeighVertex(3, Gam2)
  if(W1==NeighVertex(3, MeasGamma)) then
    is_topo_2 = .true.
  else
    is_topo_2 = .false.
  endif
  return
END FUNCTION is_topo_2


!!=======================================================================
!!=======================================================================
!!=======================================================================




!!=======================================================================
!!============================= MEASURE =================================
!!=======================================================================
SUBROUTINE measure
  implicit none
  integer :: i, iln, iGamma
  integer :: flag, it
  integer :: spg, spw
  integer :: iconf
  integer :: ityp, nloop
  integer :: MeaGin, MeaGout, MeaW, xg, yg, xw, yw, dir, typ
  integer :: NormW, Gam1, Gam2, Gam3
  double precision :: factorM


  !-------- find out the variables for Gamma ----------------
  MeaGin = NeighVertex(2, MeasGamma)
  MeaGout = NeighVertex(1, MeasGamma)
  MeaW = NeighVertex(3, MeasGamma)

  dir = DirecVertex(MeasGamma)

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

  if(TypeVertex(MeasGamma)==5 .or. TypeVertex(MeasGamma)==6) then
    typ = 11-TypeVertex(MeasGamma)
  else if(TypeVertex(MeasGamma)==1 .or. TypeVertex(MeasGamma)==3) then
    typ = TypeGW2Gam(1,1,spw,spw)
  else if(TypeVertex(MeasGamma)==2 .or. TypeVertex(MeasGamma)==4) then
    typ = TypeGW2Gam(2,2,spw,spw)
  endif

  xg = GXVertex(MeasGamma);        yg = GYVertex(MeasGamma)
  xw = WXVertex(NeighLn(3-dir, MeaW))
  yw = WYVertex(NeighLn(3-dir, MeaW))


  if(abs(OmegaLn(MeaGin))<=MxOmegaDiag .and. abs(OmegaLn(MeaGout))<=MxOmegaDiag) then

    factorM = CoefOfSymmetry(diff_x(xg,xw), diff_y(yg,yw))* &
      & CoefOfWeight(Order)*WeightLn(MeaW)*WeightVertex(MeasGamma)* &
      & WeightLn(MeaGin)*WeightLn(MeaGout)

    ime = ime + 1
    !write(175, *) imc, ime, Order, WeightCurrent
  
    nloop = Mod(NFermiloop, 2)

    GamMC(Order, nloop,(typ+1)/2, diff_x(xg,xw), diff_y(yg,yw), OmegaLn(MeaGin), OmegaLn(MeaGout)) = &
      & GamMC(Order, nloop,(typ+1)/2, diff_x(xg,xw),diff_y(yg,yw),OmegaLn(MeaGin), OmegaLn(MeaGout)) + &
      & Phase/factorM

    GamSqMC(Order,nloop,(typ+1)/2, diff_x(xg,xw), diff_y(yg,yw), OmegaLn(MeaGin), OmegaLn(MeaGout)) = &
      & GamSqMC(Order,nloop,(typ+1)/2,diff_x(xg,xw),diff_y(yg,yw), OmegaLn(MeaGin), OmegaLn(MeaGout)) + &
      & (Phase/factorM)**2.d0

    if(Order==2 .and. typ==1) then
      factorM = CoefOfWeight(Order)*WeightLn(MeaW)*WeightVertex(MeasGamma)* &
        & WeightLn(MeaGin)*WeightLn(MeaGout)
      if(is_topo_1()) then
        gam2topo(1, omegaln(meagin)) = gam2topo(1, omegaln(meagin)) + phase/factorm
      endif
      if(is_topo_2()) then
        gam2topo(2, omegaln(meagin)) = gam2topo(2, omegaln(meagin)) + phase/factorm
      endif
    endif


    if(Order==0) then
      GamNorm = GamNorm + Phase
    endif
  endif
        
END SUBROUTINE measure
!!=======================================================================
!!=======================================================================
!!=======================================================================


