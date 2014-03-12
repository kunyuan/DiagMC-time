!!====================== initialization =================================
!!=======================================================================

SUBROUTINE initialize_markov
    implicit none
    integer :: i

    CoefOfWeight(:) = 1.d0

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

END SUBROUTINE def_spin

!------------- definition of the config of diagram ----------------
SUBROUTINE def_diagram
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

  !IsDeltaVertex: 1: delta function(Gam); 0: normal function(Gam)
  IsDeltaVertex(1:4) = 1

  ! the site variables of Gamma 1, 2, 3 and 4
  GXVertex(1:4) = 0;            GYVertex(1:4) = 0
  WXVertex(1:4) = 0;            WYVertex(1:4) = 0

  ! the time variables of Gamma 1, 2, 3 and 4
  T1Vertex(1:4) = 0.d0;            T2Vertex(1:4) = 0.d0
  T3Vertex(1:4) = 0.d0           

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

  call print_config

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
  integer :: iflag, istep
  double precision :: nr

  istep = 0
  do while(istep < NStep)
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
      !case(12) 
        !call shift_gline_in_space          
      !case(13)  
        !call shift_wline_in_space         
      !case(14) 
        !call change_Gamma_type     
      !case(15) 
        !call move_measuring_index       
      case(16)  
        call change_gline_time        
      !case(17)  
        !call shift_wline_in_space         
    end select

    imc = imc + 1


    if(IsWormPresent ) then
      GamWormOrder(Order) = GamWormOrder(Order) + 1
      if(rn()<=0.5d0) call switch_ira_and_masha
    else
      GamOrder(Order) = GamOrder(Order) + 1
      istep = istep + 1
    endif
  enddo

END SUBROUTINE markov


!====================================================================
!============================== UPDATES =============================
!====================================================================
SUBROUTINE change_gline_time
    implicit none
    
end SUBROUTINE

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






!====================================================================
!============================== WEIGHTS =============================
!====================================================================

!-------- the weight of a line -------------------------
! dx = x2-x1;  dy = y2- y1; tau = tau2-tau1
COMPLEX*16 FUNCTION weight_line(stat, isdelta, knd, dx, dy, tau, typ)
  implicit none
  integer :: stat, isdelta, knd, dx, dy, typ
  double precision :: tau
  integer :: t

  t = Floor(tau*MxT/Beta)

  dx = diff_x(dx)
  dy = diff_y(dy)

  !---------------------- for test --------------------------------------
  !if(stat >= 0 .and. stat<=3) then
    !if(knd==1) weight_line = weight_meas_G(1, t)
    !if(knd==2) weight_line = weight_meas_W(1, dx, dy, t)
  !else if(stat==-1) then
    !write(*, *) IsWormPresent, iupdate, "line status == -1! There is no weight!" 
    !stop
  !else
    !write(*, *) IsWormPresent, iupdate, "line status error!", stat
    !stop
  !endif
  !------------------------ end -----------------------------------------

  if(stat == 0) then
    if(knd==1) weight_line = weight_G(typ, t)
    if(knd==2 .and. isdelta==0) weight_line = weight_W(typ, dx, dy, t)
    if(knd==2 .and. isdelta==1) weight_line = weight_W0(typ, dx, dy)
  else if(stat == 1 .or. stat==3) then
    if(knd==1) weight_line = weight_meas_G(typ, t)
    if(knd==2) weight_line = weight_meas_W(typ, dx, dy, t)
  else if(stat == 2) then
    if(knd==2) weight_line = weight_worm_W(typ, dx, dy, t)
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

!-------- the weight of a vertex -------------------------
!dx = xg-xw;  dy = yg-yw; dtau1 = tau3-tau2; dtau2 = tau1-tau3
COMPLEX*16 FUNCTION weight_vertex(stat, isdelta, dx, dy, dtau1, dtau2, typ)
  implicit none
  integer :: stat, dx, dy, t1, t2, typ, isdelta
  double precision :: weight
  double precision :: dtau1, dtau2

  t1 = Floor(dtau1*MxT/Beta)
  t2 = Floor(dtau2*MxT/Beta)

  dx = diff_x(dx)
  dy = diff_y(dy)

  !---------------------- for test --------------------------------------
  !if(stat>=0 .and. stat<=3) then
    !weight_vertex = weight_meas_gam(1, dx, dy, t1, t2)
  !else if(stat==-1) then
    !write(*, *) IsWormPresent, iupdate, "vertex status == -1! There is no weight!" 
    !stop
  !else
    !write(*, *) IsWormPresent, iupdate, "vertex status error!", stat
    !stop
  !endif
  !------------------------ end -----------------------------------------

  if(stat==0) then
      !if(isdelta==0) weight_vertex = weight_Gam(typ, dx, dy, t1, t2)
      !if(isdelta==1) weight_vertex = weight_Gam0(typ, dx, dy)
      !----------------- for bare Gamma ------------------------------
      if(isdelta==0) weight_vertex = 0.d0
      if(isdelta==1) weight_vertex = weight_Gam0(typ, dx, dy)
      !------------------------ end ----------------------------------
  else if(stat==2) then
    weight_vertex = weight_worm_Gam(typ, dx, dy, t1, t2)
  else if(stat==1 .or. stat==3) then
    weight_vertex = weight_meas_Gam(typ, dx, dy, t1, t2)
  else if(stat==-1) then 
    write(*, *) IsWormPresent, iupdate, "vertex status == -1! There is no weight!" 
    stop
  else
    write(*, *) IsWormPresent, iupdate, "vertex status error!", stat
    stop
  endif
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
    write(8, 12) iv,IsDeltaVertex(iv), TypeVertex(iv),TypeVertexIn(iv),TypeVertexOut(iv), &
      & GXVertex(iv),GYVertex(iv),WXVertex(iv),WYVertex(iv), T1Vertex(iv), T2Vertex(iv),  &
      & T3Vertex(iv), DirecVertex(iv), StatusVertex(iv), NeighVertex(:,iv)
  enddo
  write(8, *) "============================================================"

  10 format(' Line:',i2,2x,'kind:',i2,2x,'isdelta:',i2,2x,'type:',i2,2x,'k:',i8,2x,'stat:',i2, 2x,&
    & 'neigh:',i6,i6)
  12 format('Gamma:',i2,2x,'isdelta:',i2,2x,'type:',i2,2x,'typein:',i2,2x,'typeout:',i2,2x,&
    & 'gr:(',i4,i4,'), wr:(',i4,i4,')', 't:(', f6.4, f6.4, f6.4, ')',2x, &
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
  integer :: i, iln
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
    typ = TypeGW2Gam(1,1,spw,spw)
  else if(TypeVertex(MeasureGam)==2 .or. TypeVertex(MeasureGam)==4) then
    typ = TypeGW2Gam(2,2,spw,spw)
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
        
END SUBROUTINE measure
!!=======================================================================
!!=======================================================================
!!=======================================================================


!====================================================================
