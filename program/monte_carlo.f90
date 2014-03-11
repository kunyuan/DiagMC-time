!!====================== initialization =================================
!!=======================================================================

SUBROUTINE initialize_markov
    implicit none
    integer :: i

    Pupdate(:)  = 1.d0
    Pupdate(2)  = 5.d0
    Pupdate(3)  = 0.d0
    Pupdate(4)  = 0.d0

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
    call def_spin
    call def_diagram
    return
end SUBROUTINE
!------------- definition of the config of diagram ----------------

SUBROUTINE def_spin
  implicit none

END SUBROUTINE def_spin

SUBROUTINE def_diagram
  implicit none
  integer :: i, deltat

  !-------------- 1-order diagram ------------------------
  Order = 1
  ! the index of measuring gamma
  NGLn = 4;  NWLn = 2;  NVertex = 4
  ! the number of glines, wlines, gamma
  MeasureGam = 1
  ! number of fermi loops
  SignFermiLoop = 1.d0
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

  ! the site variables of Gamma 1, 2, 3 and 4
  GXVertex(1:4) = 0;            GYVertex(1:4) = 0
  WXVertex(1:4) = 0;            WYVertex(1:4) = 0

  ! the time variables of Gamma 1, 2, 3 and 4
  T1Vertex(1:4) = 0;            T2Vertex(1:4) = 0
  T3Vertex(1:4) = 0           

  ! Direction of Gamma: 1: left of W;  2: right of W
  DirecVertex(1) = 1;                DirecVertex(2) = 2
  DirecVertex(3) = 1;                DirecVertex(4) = 2

  ! NeighVertex(i, j): the ith neighbor line of gamma j
  NeighVertex(1,1) = 2;        NeighVertex(2,1) = 1;        NeighVertex(3,1) = 3
  NeighVertex(1,2) = 4;        NeighVertex(2,2) = 5;        NeighVertex(3,2) = 3
  NeighVertex(1,3) = 1;        NeighVertex(2,3) = 4;        NeighVertex(3,3) = 6
  NeighVertex(1,4) = 5;        NeighVertex(2,4) = 2;        NeighVertex(3,4) = 6

  ! weights for lines and vertexes
  !deltat = T2Vertex(NeighLn(2, 1))-T1Vertex(NeighLn(1, 1))
  !WeightLn(1) = weight_line(StatusLn(1),1,0,0,deltat,TypeLn(1))
  !WeightLn(2) = weight_line(StatusLn(2),1,0,0,deltat,TypeLn(2))
  !WeightLn(3) = weight_line(StatusLn(3),2,0,0,deltat,TypeLn(3))
  !WeightLn(4) = weight_line(StatusLn(4),1,0,0,deltat,TypeLn(4))
  !WeightLn(5) = weight_line(StatusLn(5),1,0,0,deltat,TypeLn(5))
  !WeightLn(6) = weight_line(StatusLn(6),2,0,0,deltat,TypeLn(6))

  !WeightVertex(1) = weight_vertex(1, StatusVertex(1), 0, 0, OmegaLn(2), OmegaLn(1), TypeVertex(1))
  !WeightVertex(2) = weight_vertex(1, StatusVertex(2), 0, 0, OmegaLn(4), OmegaLn(5), TypeVertex(2))
  !WeightVertex(3) = weight_vertex(1, StatusVertex(3), 0, 0, OmegaLn(1), OmegaLn(4), TypeVertex(3))
  !WeightVertex(4) = weight_vertex(1, StatusVertex(4), 0, 0, OmegaLn(5), OmegaLn(2), TypeVertex(4))


  WeightCurrent = CoefOfWeight(1)*WeightLn(1)*WeightLn(2)*WeightLn(3)*WeightLn(4)* &
    & WeightLn(5)*WeightLn(6)*WeightVertex(1)*WeightVertex(2)*WeightVertex(3)* &
    & WeightVertex(4)*(1.d0/Beta)**Order *SignFermiLoop

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

  call def_prob
  imeasure = 0
  do while(imeasure < NStep)
    nr=rn()
    iupdate = 1
    do while(nr>Fupdate(iupdate))
      iupdate = iupdate + 1
    enddo

    select case(iupdate)
      !case( 1) 
        !call create_worm_along_wline          
      !case( 2)      
        !call delete_worm_along_wline                       
      !case( 3) 
        !call create_worm_along_gline  
      !case( 4) 
        !call delete_worm_along_gline  
      !case( 5) 
        !call move_worm_along_wline             
      !case( 6) 
        !call move_worm_along_gline              
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
      case(12) 
        call change_gline_space          
      case(13)  
        call change_wline_space         
      !case(14) 
        !call change_Gamma_type     
      !case(15) 
        !call move_measuring_index       
      case(16)  
        call change_gline_time        
      case(17)  
        call change_wline_time        
    end select

    imc = imc + 1


    if(IsWormPresent ) then
      GamWormOrder(Order) = GamWormOrder(Order) + 1
      if(rn()<=0.5d0) call switch_ira_and_masha
    else
      GamOrder(Order) = GamOrder(Order) + 1
      imeasure = imeasure + 1
    endif
  enddo

END SUBROUTINE markov


!===============================================================
!========================= updates =============================
!===============================================================
SUBROUTINE change_gline_time
    implicit none
    
end SUBROUTINE

SUBROUTINE change_wline_time
    implicit none
    integer :: iGam, iWLn, jGam, dxw, dyw, xw, yw, xwj, ywj, dir
    double precision :: WGam, WW, Anew, Aold, Pacc, sgn

    !------- step1 : check if worm is present -------------
    if(IsWormPresent .eqv. .true.)    return
    !ProbProp(iupdate) = ProbProp(iupdate) + 1

    !------- step2 : propose a new config -----------------
    iGam = generate_gamma()
    if(IsDeltaVertex(iGam)==1) return
    dir  = DirecVertex(iGam)
    iWLn = NeighVertex(3, iGam)
    jGam = NeighLn(3-dir, iWLn)

    dxw = generate_x();         dyw = generate_y()
    xw = find_neigh_x(WXVertex(iGam), dxw)
    yw = find_neigh_y(WYVertex(iGam), dyw)

    xwj = WXVertex(jGam);       ywj = WYVertex(jGam)
    
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
end SUBROUTINE

