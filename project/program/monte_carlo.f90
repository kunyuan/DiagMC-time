!====================== initialization =================================
!!=======================================================================

SUBROUTINE initialize_markov
    implicit none
    integer :: i
    double precision :: ratioleft


    call LogFile%QuickLog("Initializing the configuration...")


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

    call LogFile%QuickLog("Initialization Done!")
    call LogFile%QuickLog("Weight of Initialization Diagram: "//trim(str(WeightCurrent)))

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
  double precision :: tt
  integer :: ix,iy,i

  ! initializing logscale sampling of space
  ! for seeding rx coordinate displacement do this
  !     x=rndm(); rx=0.5d0*dexp(x*logL); ix=rx
  !     IF(rndm()>0.5d0) ix=L(1)-1-ix
  ! this gives you displacement ix in the affine basis called with probability
  ! SpatialWeight(i,ix)

  do i=1,D
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
  TypeGam2W(1,2) = 3;      TypeGam2W(1,3) = 3

  TypeGam2W(3,1) = 4;      TypeGam2W(3,4) = 4
  TypeGam2W(3,2) = 2;      TypeGam2W(3,3) = 2

  TypeGam2W(2,1) = 4;      TypeGam2W(2,4) = 4
  TypeGam2W(2,2) = 2;      TypeGam2W(2,3) = 2

  TypeGam2W(4,1) = 1;      TypeGam2W(4,4) = 1
  TypeGam2W(4,2) = 3;      TypeGam2W(4,3) = 3


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

SUBROUTINE markov(IsToss)
  implicit none
  integer :: istep,isamp,i,MaxSamp,iblck
  double precision :: nr,x, ratioleft
  logical :: IsToss

  iblck=0
  mc_version = 0
  if(IsToss) then
    MaxSamp = Ntoss
  else
    MaxSamp = Nsamp
  endif

  !------- assign the config ratio for each order ------
  TimeRatio(0) = 0.25d0
  ratioleft = 1.d0 - TimeRatio(0)
  do i = MCOrder, 1, -1
    TimeRatio(i) = 0.5d0*(ratioleft)
    ratioleft = ratioleft - TimeRatio(i)
  enddo
  TimeRatio(1) = TimeRatio(1) + ratioleft

  call LogFile%QuickLog("Starting Markov...")
  call check_config

  do while(.true.)
    iblck = iblck +1
    do isamp = 1, MaxSamp
      do istep=1, Nstep
        imc = imc + 1.0

        if(IsWormPresent .and. rn()<0.50) call switch_ira_and_masha

        nr=rn()
        if(nr<Fupdate(1)) then
          iupdate = 1
          call create_worm_along_wline          
        else if(nr<Fupdate(2)) then
          iupdate = 2
          call delete_worm_along_wline                       
        else if(nr<Fupdate(3)) then
          iupdate = 3
          !call create_worm_along_gline  
        else if(nr<Fupdate(4)) then
          iupdate = 4
          !call delete_worm_along_gline  
        else if(nr<Fupdate(5)) then
          iupdate = 5
          call move_worm_along_wline             
        else if(nr<Fupdate(6)) then
          iupdate = 6
          !call move_worm_along_gline              
          call move_worm_along_gline_test           
        else if(nr<Fupdate(7)) then
          iupdate = 7
          ProbCall(Order, iupdate) = ProbCall(Order, iupdate) + 1
          call add_interaction                  
        else if(nr<Fupdate(8)) then
          iupdate = 8
          ProbCall(Order, iupdate) = ProbCall(Order, iupdate) + 1
          call remove_interaction             
        else if(nr<Fupdate(9)) then
          iupdate = 9
          !call add_interaction_cross              
        else if(nr<Fupdate(10)) then
          iupdate = 10
          !call remove_interaction_cross                  
        else if(nr<Fupdate(11)) then
          iupdate = 11
          call reconnect                      
        else if(nr<Fupdate(12)) then
          iupdate = 12
          call change_gline_space          
        else if(nr<Fupdate(13)) then
          iupdate = 13
          call change_wline_space         
        else if(nr<Fupdate(14)) then
          iupdate = 14
          call change_Gamma_type     
        else if(nr<Fupdate(15)) then
          iupdate = 15
          call move_measuring_index       
        else if(nr<Fupdate(16)) then
          iupdate = 16
          call change_Gamma_time        
        else if(nr<Fupdate(17)) then
          iupdate = 17
          call change_wline_isdelta       
        else if(nr<Fupdate(18)) then
          iupdate = 18
          call change_Gamma_isdelta       
        endif

        if(iupdate/=7 .and. iupdate/=8) then
          ProbCall(Order, iupdate) = ProbCall(Order, iupdate) + 1
        endif
      enddo

      if( .not. IsToss) call measure
    enddo 

    if(IsToss .and. iblck==1) return

    !========================== REWEIGHTING =========================
    if(mod(iblck, 10)==0) then

      call check_config

      call statistics
      call output_GamMC

      call print_status
      call print_config

      call LogFile%QuickLog("Check if there is a new G,W data...")
      call read_flag

      if(mc_version/=file_version) then
        call LogFile%QuickLog(str(mc_version)+', '+str(file_version))
        call LogFile%QuickLog("Updating G, W, and Gamma...")

        call read_GWGamma
        call update_WeightCurrent
        call check_config
        call print_config
        mc_version = file_version

        call recalculate_Reweighting
      endif


      if(iblck <= 10) then
        GamWormOrder=0.d0
        GamOrder=0.d0
      endif
    endif
    !================================================================

    if(mod(iblck,30)==0) then

      call check_config

      call LogFile%QuickLog("Writing data and configuration...")

      call statistics
      call write_monte_carlo_conf
      call write_monte_carlo_data

      call LogFile%QuickLog("Writing data and configuration done!")
    endif

    if(iblck>Mxint) then
      call LogFile%QuickLog("Block number:"+str(Mxint)+" achieved! Code stops.")
      stop
    endif
  enddo
END SUBROUTINE markov

SUBROUTINE recalculate_Reweighting
  implicit none
  integer :: i
  double precision :: x
  
  call LogFile%WriteStamp()
  call LogFile%WriteLine("Reweighting order of diagrams...")

  !------ if the errorbar is too big, just don't reweight ----
  do i = 0, MCOrder
    if(Error(i)*Norm(i)/Quan(i)>MxError) then
      call LogFile%WriteLine("Error bar too large! Order "+str(i)+" err: "+ &
        & str(Error(i)*Norm(i)/Quan(i)))
      call LogFile%WriteLine("Reweighting is done!")
      return
    endif
  enddo

  x = SUM(GamWormOrder(:))
  CoefOfWeight(:)=TimeRatio(:)*CoefOfWeight(:)*x/(GamWormOrder(:)+50.d0)
  CoefOfWeight(1:MCOrder) = CoefOfWeight(1:MCOrder)/CoefOfWeight(0)
  CoefOfWeight(0) = 1.d0

  call LogFile%WriteLine("Reweight Ratios:")

  do i=0,MCOrder
    call LogFile%WriteLine('Order'+str(i)+' :'+str(CoefOfWeight(i)))
    call LogFile%WriteLine('     worm: '+str(GamWormOrder(i)))
    call LogFile%WriteLine('     phycical: '+str(GamOrder(i)))
  enddo

  call LogFile%WriteLine("Reweighting is done!")
  return
END SUBROUTINE

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
  integer :: k, kiW, ktemp
  integer :: kiWold
  integer :: spin
  integer :: tempGam
  double precision :: tau, tau1, tau2, Pacc
  double precision :: WWorm
  complex*16 :: WiGin, WiGout, WjGin, WjGout
  complex*16 :: WIra, WMasha, WWNew, Anew, Aold, sgn
  logical :: flag

  !------------ step1 : check if worm is present ------------------
  if(IsWormPresent .eqv. .true.)   return

  !------------ step2 : propose the new config ------------------
  iWLn = generate_wline()

  iGam = NeighLn(1, iWLn);              jGam = NeighLn(2, iWLn)
  iGin = NeighVertex(1, iGam);          jGin = NeighVertex(1, jGam)
  iGout = NeighVertex(2, iGam);         jGout = NeighVertex(2, jGam)

  spin = 4*(Floor(rn()*2.d0))-2         !spin = +/- 2
  if(spin-2*(TypeLn(jGin)-TypeLn(jGout))==0)  return  !in this case, I/M cannot move

  k = generate_k()
  kiWold = kLn(iWLn)
  ktemp=add_k(kiWold, k)
  if(Is_k_valid_for_W(ktemp)==.false.) return
  call update_line(iWLn, ktemp, 2)
   

  !------------ step3 : configuration check -------------------
  flag=Is_reducible_W_Gam_both_side(ktemp, kLn(iGin), kLn(iGout), kLn(jGin), kLn(jGout))

  if(flag) then
    call undo_update_line(iWLn, kiWold, 2)
    return
  endif

  statiGam = add_ira_stat(StatusVertex(iGam))
  statjGam = add_ira_stat(StatusVertex(jGam))
  statW    = wline_stat(statiGam, statjGam)

  !------------ step4 : weight calculation ----------------------
  tau = TVertex(2, NeighLn(2, iGin))- TVertex(1, NeighLn(1, iGin))
  WiGin  = weight_gline(StatusLn(iGin), tau, TypeLn(iGin))

  tau = TVertex(2, NeighLn(2, iGout))- TVertex(1, NeighLn(1, iGout))
  WiGout = weight_gline(StatusLn(iGout), tau, TypeLn(iGout))

  tau = TVertex(2, NeighLn(2, jGin))- TVertex(1, NeighLn(1, jGin))
  WjGin  = weight_gline(StatusLn(jGin), tau, TypeLn(jGin))

  tau = TVertex(2, NeighLn(2, jGout))- TVertex(1, NeighLn(1, jGout))
  WjGout = weight_gline(StatusLn(jGout), tau, TypeLn(jGout))

  tau1 = TVertex(3, iGam)- TVertex(2, iGam)
  tau2 = TVertex(1, iGam)- TVertex(3, iGam)
  WIra = weight_vertex(statiGam, IsDeltaVertex(iGam),diff_r(D,GRVertex(iGam),WRVertex(iGam)),  &
    & tau1, tau2, TypeVertex(iGam))

  tau1 = TVertex(3, jGam)- TVertex(2, jGam)
  tau2 = TVertex(1, jGam)- TVertex(3, jGam)
  WMasha = weight_vertex(statjGam, IsDeltaVertex(jGam), diff_r(D,GRVertex(jGam),WRVertex(jGam)),  &
    & tau1, tau2, TypeVertex(jGam))

  tau = TVertex(3, jGam) - TVertex(3, iGam)
  WWNew = weight_wline(statW, IsDeltaLn(iWLn), diff_r(D,WRVertex(jGam),WRVertex(iGam)), &
    & tau, TypeLn(iWLn))

  tau = TVertex(3, jGam) - TVertex(3, iGam)
  WWorm = weight_worm(diff_r(D,GRVertex(iGam),GRVertex(jGam)), diff_r(D, WRVertex(iGam), &
    & WRVertex(jGam)), tau) 

  Anew = WIra *WMasha *WWNew *WiGin *WiGout *WjGin *WjGout
  Aold = WeightLn(iWLn) *WeightVertex(iGam) *WeightVertex(jGam) *WeightLn(iGin) &
  & *WeightLn(iGout) *WeightLn(jGin) *WeightLn(jGout)

  if(abs(Anew)==0.d0) then
    call undo_update_line(iWLn, kiWold, 2)
    return
  endif

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
  Pacc = Pacc *CoefOfWorm *Wworm *(Order+1.d0)/0.5d0 
  Pacc = Pacc *Pupdate(2)/Pupdate(1)

  !------------ step5 : accept the update -----------------------
  ProbProp(Order, 1) = ProbProp(Order, 1) + 1

  if(rn()<Pacc) then

    !-------------- update the diagram info --------------------
    Phase = Phase*sgn
    IsWormPresent = .true.

    !---------- update Ira and Masha -------------------
    Ira = iGam
    Masha = jGam
    SpinMasha = spin
    kMasha = k

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

    ProbAcc(Order, 1) = ProbAcc(Order, 1) + 1
  else
    call undo_update_line(iWLn, kiWold, 2)
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
  integer :: k, kiW, kiWold,ktemp
  double precision :: Pacc, tau, tau1, tau2
  complex*16 :: WiGin, WiGout, WjGin, WjGout
  complex*16 :: Wi, Wj, WWNew, Anew, Aold, sgn
  logical :: flag

  !------------ step1 : check if worm is present ------------------
  if(IsWormPresent .eqv. .false.)  return
  if(NeighVertex(3, Ira)/= NeighVertex(3, Masha))   return

  !------------ step2 : propose the new config ------------------
  iWLn = NeighVertex(3, Ira)
  dir = DirecVertex(Ira)

  iGin  = NeighVertex(1, Ira);          jGin  = NeighVertex(1, Masha)
  iGout = NeighVertex(2, Ira);          jGout = NeighVertex(2, Masha)

  k = kMasha
  kiWold = kLn(iWLn)
  ktemp = add_k(kiWold, (-1)**dir*k)
  if(Is_k_valid_for_W(ktemp)==.false.) return
  call update_line(iWLn, ktemp, 2)

  !------------ step3 : configuration check ---------------------
  flag=Is_reducible_W_Gam_both_side(ktemp, kLn(iGin), kLn(iGout), kLn(jGin), kLn(jGout))

  if(flag) then
    call undo_update_line(iWLn, kiWold, 2)
    return
  endif

  !------------ status, type for the new configuration ------------
  statIra   = delete_ira_stat(StatusVertex(Ira))
  statMasha = delete_ira_stat(StatusVertex(Masha))
  statW     = wline_stat(statIra, statMasha)

  typIra = TypeSp2Gam(TypeLn(iGin),TypeLn(iGout),SpInVertex(1, Ira), SpInVertex(2, Ira))
  typMasha = TypeSp2Gam(TypeLn(jGin),TypeLn(jGout),SpInVertex(1,Masha),SpInVertex(2, Masha))
  if(dir==1) then
    typW = TypeGam2W(typIra, typMasha)
  else
    typW = TypeGam2W(typMasha, typIra)
  endif


  !------------ step4 : weight calculation ---------------------
  tau = TVertex(2, NeighLn(2, iGin))- TVertex(1, NeighLn(1, iGin))
  WiGin  = weight_gline(StatusLn(iGin), tau, TypeLn(iGin))

  tau = TVertex(2, NeighLn(2, iGout))- TVertex(1, NeighLn(1, iGout))
  WiGout = weight_gline(StatusLn(iGout), tau, TypeLn(iGout))

  tau = TVertex(2, NeighLn(2, jGin))- TVertex(1, NeighLn(1, jGin))
  WjGin  = weight_gline(StatusLn(jGin), tau, TypeLn(jGin))

  tau = TVertex(2, NeighLn(2, jGout))- TVertex(1, NeighLn(1, jGout))
  WjGout = weight_gline(StatusLn(jGout), tau, TypeLn(jGout))

  tau1 = TVertex(3, Ira)- TVertex(2, Ira)
  tau2 = TVertex(1, Ira)- TVertex(3, Ira)
  Wi = weight_vertex(statIra, IsDeltaVertex(Ira), diff_r(D, GRVertex(Ira),WRVertex(Ira)),  &
    & tau1, tau2, typIra)

  tau1 = TVertex(3, Masha)- TVertex(2, Masha)
  tau2 = TVertex(1, Masha)- TVertex(3, Masha)
  Wj = weight_vertex(statMasha, IsDeltaVertex(Masha), diff_r(D,GRVertex(Masha),WRVertex(Masha)),  &
    & tau1, tau2, typMasha)

  tau = (-1)**dir*(TVertex(3, Ira) - TVertex(3, Masha))
  WWNew = weight_wline(statW, IsDeltaLn(iWLn),diff_r(D,WRVertex(Masha),WRVertex(Ira)), &
    & tau, typW)

  Anew = Wi *Wj *WWNew *WiGin *WiGout *WjGin *WjGout
  Aold = WeightLn(iWLn) *WeightVertex(Ira) *WeightVertex(Masha) *WeightLn(iGin) &
  & *WeightLn(iGout) *WeightLn(jGin) *WeightLn(jGout)

  if(abs(Anew)==0.d0) then
    call undo_update_line(iWLn, kiWold, 2)
    return
  endif

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
  Pacc = Pacc *0.5d0/(CoefOfWorm*WeightWorm*(Order+1.d0)) 
  Pacc = Pacc *Pupdate(1)/Pupdate(2)

  ProbProp(Order, 2) = ProbProp(Order, 2) + 1
  !------------ step5 : accept the update -----------------------
  if(Pacc>1.d-12 .and. rn()<Pacc) then

    !---------- update the diagram info ----------------
    Phase = Phase*sgn
    IsWormPresent = .false.

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
    call undo_update_line(iWLn, kiWold, 2)
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
  integer :: dir, kW, kWold, ktemp
  integer :: iGam, jGam, WLn, GLn1, GLn2, GLn3, GLn4
  integer :: typiGam, typW
  integer :: statiGam, statjGam, statW
  double precision :: tau1, tau2, tau
  double precision :: WWorm, Pacc
  complex*16 :: WiGam, WjGam, WW
  complex*16 :: Anew, Aold, sgn
  logical :: flag

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
  ktemp = add_k(kWold, (-1)**dir *kMasha)

  !print *, Ira, kWold, kMasha, ktemp
  !call print_config
  !if(SignFermiloop==1.d0) stop

  if(Is_k_valid_for_W(ktemp)==.false.) return
  call update_line(WLn, ktemp, 2)

  flag=Is_reducible_W_Gam_one_side(ktemp, kLn(GLn1), kLn(GLn2))

  if(flag) then
    call undo_update_line(WLn, kWold, 2)
    return
  endif

  !-------- the stat, type of the new config ----
  statiGam = delete_ira_stat(StatusVertex(iGam))
  statjGam = add_ira_stat(StatusVertex(jGam))
  statW = wline_stat(statiGam, statjGam)

  typW = TypeLn(WLn)
  typiGam = TypeSp2Gam(TypeLn(GLn1), TypeLn(GLn2), SpInVertex(1, iGam), SpInVertex(2, iGam))

  !------- step4 : weight calculation -------------------
  WWorm = weight_worm(diff_r(D, GRVertex(jGam),GRVertex(Masha)), &
    & diff_r(D,WRVertex(jGam),WRVertex(Masha)), TVertex(3, jGam)-TVertex(3, Masha))

  tau1 = TVertex(3, iGam)-TVertex(2, iGam)
  tau2 = TVertex(1, iGam)-TVertex(3, iGam)
  WiGam = weight_vertex(statiGam, IsDeltaVertex(iGam), diff_r(D, GRVertex(iGam),WRVertex(iGam)), &
    & tau1, tau2, typiGam)

  tau1 = TVertex(3, jGam)-TVertex(2, jGam)
  tau2 = TVertex(1, jGam)-TVertex(3, jGam)
  WjGam = weight_vertex(statjGam, IsDeltaVertex(jGam), diff_r(D, GRVertex(jGam),WRVertex(jGam)), &
    & tau1, tau2, TypeVertex(jGam))

  tau = (-1.d0)**dir *(TVertex(3, iGam)-TVertex(3, jGam))
  WW = weight_wline(statW, IsDeltaLn(WLn), diff_r(D, WRVertex(iGam),WRVertex(jGam)), &
    & tau, typW)

  Anew = WiGam *WjGam *WW 
  Aold = WeightVertex(iGam) *WeightVertex(jGam) *WeightLn(WLn) 

  call weight_ratio(Pacc, sgn, Anew, Aold)

  Pacc = Pacc *WWorm/WeightWorm

  !------- step5 : accept the update --------------------
  ProbProp(Order, 5) = ProbProp(Order, 5) + 1
  if(Pacc>1.d-12 .and. rn()<Pacc) then

    !------ update diagram info ---------------
    Phase = Phase *sgn

    !------ update Ira ------------------------
    Ira = jGam

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
    call undo_update_line(WLn, kWold, 2)
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
  integer :: dir, kGold, ktemp, dr
  double precision :: Ppro
  integer :: iGam, jGam, iW, jW, GLn, kiGLn
  integer :: vertex1, vertex2
  integer :: typG, typiGam, typjGam, typiW
  integer, dimension(2) :: spiGam, spjGam
  integer :: statiGam, statjGam, statiW, statjW, stat1, stat2
  double precision :: tau1, tau2, tau
  double precision :: WWorm, Pacc
  complex*16 :: WiW, WjW, WiGam, WjGam, WG, Anew, Aold, sgn
  logical :: flag

  !------- step1 : check if worm is present -------------
  if(IsWormPresent .eqv. .false.)   return

  !------- step2 : propose a new config -----------------
  iGam = Ira
  iW = NeighVertex(3, iGam)
  dir = Floor(rn()*2.d0) + 1
  GLn = NeighVertex(dir, iGam)
  jGam = NeighLn(dir, GLn)
  jW = NeighVertex(3, jGam)
  if(jGam==Masha .or. jGam==Ira)    return

  if(dir==1) then
    if(SpinMasha+TypeLn(GLn)==0 .or. SpinMasha+TypeLn(GLn)==3)  return
  else 
    if(SpinMasha+TypeLn(GLn)==-1 .or. SpinMasha+TypeLn(GLn)==4) return
  endif
  
  kGold = kLn(GLn)
  ktemp = add_k(kGold, (-1)**(dir+1) *kMasha)
  
  if(Is_k_valid_for_G(ktemp)==.false.) return
  call update_line(GLn, ktemp, 1)

  !------- step3 : configuration check ------------------
  !-------- irreducibility  check -----------------------
  kiGLn = kLn(NeighVertex(3-dir, iGam))
  flag=Is_reducible_G_Gam_one_side(ktemp, kiGLn)

  if(flag) then
    call undo_update_line(GLn, kGold, 1)
    return
  endif

  !----- type and status of the new config ---
  statiGam = delete_ira_stat(StatusVertex(iGam))
  statjGam = add_ira_stat(StatusVertex(jGam))

  stat1 = StatusVertex(NeighLn(3-DirecVertex(iGam), iW))
  statiW = wline_stat(statiGam,  stat1)

  stat2 = StatusVertex(NeighLn(3-DirecVertex(jGam), jW))
  statjW = wline_stat(statjGam,  stat2)

  if(iW == jW) then
    statiW = wline_stat(statiGam, statjGam)
    statjW = statiW
  endif

  typG = 3 - TypeLn(GLn)

  if(TypeVertex(iGam)==1 .or. TypeVertex(iGam)==2) then
    Ppro = 2.d0

    spiGam(dir) = 3-SpInVertex(dir, iGam)
    spiGam(3-dir) = SpInVertex(3-dir, iGam)
  else if(TypeVertex(iGam)==3 .or. TypeVertex(iGam)==4) then
    Ppro = 2.d0

    spiGam(dir) = SpInVertex(dir, iGam)
    spiGam(3-dir) = 3-SpInVertex(3-dir, iGam)
  else if(TypeVertex(iGam)==5 .or. TypeVertex(iGam)==6) then
    Ppro = 0.5d0

    if(rn()<0.5d0) then
      spiGam(dir) = 3-SpInVertex(dir, iGam)
      spiGam(3-dir) = SpInVertex(3-dir, iGam)
    else
      spiGam(dir) = SpInVertex(dir, iGam)
      spiGam(3-dir) = 3-SpInVertex(3-dir, iGam)
    endif
  endif

  if(TypeVertex(jGam)==1 .or. TypeVertex(jGam)==2) then
    Ppro = Ppro*2.d0

    spjGam(3-dir) = 3-SpInVertex(3-dir, jGam)
    spjGam(dir) = SpInVertex(dir, jGam)
  else if(TypeVertex(jGam)==3 .or. TypeVertex(jGam)==4) then
    Ppro = Ppro*2.d0

    spjGam(3-dir) = SpInVertex(3-dir, jGam)
    spjGam(dir) = 3-SpInVertex(dir, jGam)
  else if(TypeVertex(jGam)==5 .or. TypeVertex(jGam)==6) then
    Ppro = Ppro*0.5d0

    if(rn()<0.5d0) then
      spjGam(3-dir) = 3-SpInVertex(3-dir, jGam)
      spjGam(dir) = SpInVertex(dir, jGam)
    else
      spjGam(3-dir) = SpInVertex(3-dir, jGam)
      spjGam(dir) = 3-SpInVertex(dir, jGam)
    endif
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
  WWorm = weight_worm(diff_r(D, GRVertex(jGam),GRVertex(Masha)), &
    & diff_r(D, WRVertex(jGam), WRVertex(Masha)), TVertex(3, jGam)-TVertex(3, Masha))

  tau1 = TVertex(3, iGam)-TVertex(2, iGam)
  tau2 = TVertex(1, iGam)-TVertex(3, iGam)
  WiGam = weight_vertex(statiGam, IsDeltaVertex(iGam),diff_r(D, GRVertex(iGam),WRVertex(iGam)),   &
    & tau1, tau2, typiGam)

  tau1 = TVertex(3, jGam)-TVertex(2, jGam)
  tau2 = TVertex(1, jGam)-TVertex(3, jGam)
  WjGam = weight_vertex(statjGam, IsDeltaVertex(jGam),diff_r(D, GRVertex(jGam),WRVertex(jGam)),   &
    & tau1, tau2, typjGam)

  vertex1 = NeighLn(1, iW)
  vertex2 = NeighLn(2, iW)
  tau = TVertex(3, vertex2) - TVertex(3, vertex1)
  dr  = diff_r(D, WRVertex(vertex2),WRVertex(vertex1))
  WiW = weight_wline(statiW, IsDeltaLn(iW), dr, tau, typiW)

  vertex1 = NeighLn(1, jW)
  vertex2 = NeighLn(2, jW)
  tau = TVertex(3, vertex2) - TVertex(3, vertex1)
  dr  = diff_r(D, WRVertex(vertex2),WRVertex(vertex1))
  WjW = weight_wline(statjW, IsDeltaLn(jW), dr, tau, TypeLn(jW))

  tau = TVertex(2, NeighLn(2, GLn)) - TVertex(1, NeighLn(1, GLn))
  WG = weight_gline(StatusLn(GLn), tau, typG)

  Anew = WiGam *WjGam *WiW *WjW *WG
  Aold = WeightVertex(iGam) *WeightVertex(jGam) *WeightLn(iW) *WeightLn(jW) *WeightLn(GLn)

  if(abs(Anew)==0.d0) then
    call undo_update_line(GLn, kGold, 1)
    return
  endif

  if(iW==jW) then
    Anew = Anew/WjW
    Aold = Aold/WeightLn(jW)
  endif

  call weight_ratio(Pacc, sgn, Anew, Aold)

  Pacc = Pacc *WWorm/WeightWorm/Ppro

  !------- step5 : accept the update --------------------
  ProbProp(Order, 6) = ProbProp(Order, 6) + 1
  if(Pacc>1.d-12 .and. rn()<Pacc) then

    !----- update the diagram info -------------- 
    Phase = Phase *sgn

    !---- update Ira and Masha ------------------
    Ira = jGam

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
    call undo_update_line(GLn, kGold, 1)
    return
  endif
END SUBROUTINE move_worm_along_gline
!-----------------------------------------------------------




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
SUBROUTINE move_worm_along_gline_test
  implicit none
  integer :: iGam, jGam, iGam2, jGam2, iW, jW, GLn,kiGLn, kjGLn
  integer :: vertex1, vertex2
  integer :: sG, typG, typiGam, typjGam, typiW
  integer, dimension(2) :: spiGam, spjGam
  integer :: dir, kGold, dr, ktemp
  double precision :: Ppro
  integer :: statiGam, statjGam, statiW, statjW, statG, stat1, stat2
  double precision :: tau1, tau2, tau
  double precision :: WWorm, Pacc
  complex*16 :: WiW, WjW, WiGam, WjGam, WG, Anew, Aold, sgn
  logical :: flag

  if(IsWormPresent .eqv. .false.)  return

  iGam = Ira
  dir = Floor(rn()*2.d0) + 1
  GLn = NeighVertex(dir, iGam)
  jGam = NeighLn(dir, GLn)
  if(jGam==Ira .or. jGam==Masha) return

  iW = NeighVertex(3, iGam)
  jW = NeighVertex(3, jGam)

  sG = get_spin_G(TypeLn(GLn))

  !!deltasG+spinMasha*(-1)**dir==0
  !!deltasG = -2*sG
  if(-2*sG+SpinMasha*(-1)**dir/=0)  return

  kGold = kLn(GLn)
  ktemp = add_k(kGold, -(-1)**dir*kMasha)

  if(Is_k_valid_for_G(ktemp)==.false.) return
  call update_line(GLn, ktemp, 1)

  !-------- irreducibility  check -----------------------
  kiGLn = kLn(NeighVertex(3-dir, iGam))
  flag=Is_reducible_G_Gam_one_side(ktemp, kiGLn)

  if(flag) then
    call undo_update_line(GLn, kGold, 1)
    return
  endif

  statiGam = delete_ira_stat(StatusVertex(iGam))
  statjGam = add_ira_stat(StatusVertex(jGam))

  iGam2 = NeighLn(3-DirecVertex(iGam), iW)
  statiW = wline_stat(statiGam, StatusVertex(iGam2))
  if(iGam2==jGam) then
    statiW = wline_stat(statiGam, statjGam)
  endif

  jGam2 = NeighLn(3-DirecVertex(jGam), jW)
  statjW = wline_stat(statjGam, StatusVertex(jGam2))
  if(jGam2==iGam) then
    statjW = wline_stat(statiGam, statjGam)
  endif

  typG = flip_typ_G(TypeLn(GLn))

  call flip_typ_Vertex(dir, TypeVertex(iGam), spiGam, typiGam)
  call flip_typ_Vertex(3-dir, TypeVertex(jGam), spjGam, typjGam)

  if(is_worm_nearby(statiW)) then
    typiW = TypeLn(iW)
  else
    typiW = get_typW_from_Vertex(DirecVertex(iGam), typiGam, TypeVertex(iGam2))
  endif

  WWorm = weight_worm(diff_r(D, GRVertex(jGam),GRVertex(Masha)), &
    & diff_r(D,WRVertex(jGam),WRVertex(Masha)), TVertex(3, jGam)-TVertex(3, Masha))

  tau1 = TVertex(3, iGam)-TVertex(2, iGam)
  tau2 = TVertex(1, iGam)-TVertex(3, iGam)
  WiGam = weight_vertex(statiGam, IsDeltaVertex(iGam),diff_r(D, GRVertex(iGam),WRVertex(iGam)),   &
    & tau1, tau2, typiGam)

  tau1 = TVertex(3, jGam)-TVertex(2, jGam)
  tau2 = TVertex(1, jGam)-TVertex(3, jGam)
  WjGam = weight_vertex(statjGam, IsDeltaVertex(jGam),diff_r(D, GRVertex(jGam),WRVertex(jGam)),   &
    & tau1, tau2, typjGam)

  vertex1 = NeighLn(1, iW)
  vertex2 = NeighLn(2, iW)
  tau = TVertex(3, vertex2) - TVertex(3, vertex1)
  dr  = diff_r(D, WRVertex(vertex2), WRVertex(vertex1))
  WiW = weight_wline(statiW, IsDeltaLn(iW), dr, tau, typiW)

  vertex1 = NeighLn(1, jW)
  vertex2 = NeighLn(2, jW)
  tau = TVertex(3, vertex2) - TVertex(3, vertex1)
  dr  = diff_r(D, WRVertex(vertex2), WRVertex(vertex1))
  WjW = weight_wline(statjW, IsDeltaLn(jW), dr, tau, TypeLn(jW))

  tau = TVertex(2, NeighLn(2, GLn)) - TVertex(1, NeighLn(1, GLn))
  WG = weight_gline(StatusLn(GLn), tau, typG)

  Anew = WiGam *WjGam *WiW *WjW *WG
  Aold = WeightVertex(iGam) *WeightVertex(jGam) *WeightLn(iW) *WeightLn(jW) *WeightLn(GLn)

  if(abs(Anew)==0.d0) then
    call undo_update_line(GLn, kGold, 1)
    return
  endif

  if(iW==jW) then
    Anew = Anew/WjW
    Aold = Aold/WeightLn(jW)
  endif

  call weight_ratio(Pacc, sgn, Anew, Aold)

  if(TypeVertex(iGam)==5 .or. TypeVertex(iGam)==6) then
    Ppro = 0.5d0
  else if(TypeVertex(iGam)<=4) then
    Ppro = 2.0d0
  endif

  if(TypeVertex(jGam)==5 .or. TypeVertex(jGam)==6) then
    Ppro = Ppro*0.5d0
  else if(TypeVertex(jGam)<=4) then
    Ppro = Ppro*2.0d0
  endif

  Pacc = Pacc *WWorm/WeightWorm/Ppro

  !------- step5 : accept the update --------------------
  ProbProp(Order, 6) = ProbProp(Order, 6) + 1
  if(Pacc>1.d-12 .and. rn()<Pacc) then

    !----- update the diagram info -------------- 
    Phase = Phase *sgn

    !---- update Ira and Masha ------------------
    Ira = jGam

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
    call undo_update_line(GLn, kGold, 1)
    return
  endif

  return
END SUBROUTINE move_worm_along_gline_test

INTEGER FUNCTION get_spin_G(typ)
  implicit none
  integer, intent(in) :: typ
  if(typ==1) then
    get_spin_G = 1
  else if(typ==2) then
    get_spin_G = -1
  endif
  return
END FUNCTION get_spin_G

INTEGER FUNCTION flip_typ_G(typ)
  implicit none
  integer, intent(in) :: typ
  if(typ==1) then
    flip_typ_G = 2
  elseif(typ==2) then
    flip_typ_G = 1
  endif
  return
END FUNCTION flip_typ_G


SUBROUTINE get_spin_Vertex(typ, sp)
  implicit none
  integer, intent(in) :: typ
  integer, intent(out) :: sp(4)

  !TypeSp2Gam(1,1,1,1) = 1
  !TypeSp2Gam(2,2,2,2) = 2
  !TypeSp2Gam(1,1,2,2) = 3
  !TypeSp2Gam(2,2,1,1) = 4
  !TypeSp2Gam(1,2,1,2) = 5
  !TypeSp2Gam(2,1,2,1) = 6

  if(typ==1) then
    sp(1)=1; sp(2)=1; sp(3)=1; sp(4)=1
  else if(typ==2) then
    sp(1)=2; sp(2)=2; sp(3)=2; sp(4)=2
  else if(typ==3) then
    sp(1)=1; sp(2)=1; sp(3)=2; sp(4)=2
  else if(typ==4) then
    sp(1)=2; sp(2)=2; sp(3)=1; sp(4)=1
  else if(typ==5) then
    sp(1)=1; sp(2)=2; sp(3)=1; sp(4)=2
  else if(typ==6) then
    sp(1)=2; sp(2)=1; sp(3)=2; sp(4)=1
  endif
  return
END SUBROUTINE get_spin_Vertex

SUBROUTINE flip_typ_Vertex(dir, typ, newspinside, newtyp)
  implicit none
  integer, intent(in) :: dir, typ
  integer, intent(out) :: newspinside(2), newtyp
  integer :: oldsp(4), newsp(4), oldspinside(2)

  call get_spin_Vertex(typ, oldsp)

  newsp(dir) = flip_typ_G(oldsp(dir))
  newsp(3-dir) = oldsp(3-dir)

  oldspinside(1) = oldsp(3)
  oldspinside(2) = oldsp(4)

  if(typ==1 .or. typ==2) then
    newspinside(dir) = flip_typ_G(oldspinside(dir))
    newspinside(3-dir) = oldspinside(3-dir)

  else if(typ==3 .or. typ==4) then
    newspinside(dir) = oldspinside(dir)
    newspinside(3-dir) = flip_typ_G(oldspinside(3-dir))

  else if(typ==5 .or. typ==6) then
    if(rn()<0.5d0) then
      newspinside(dir) = flip_typ_G(oldspinside(dir))
      newspinside(3-dir) = oldspinside(3-dir)
    else 
      newspinside(dir) = oldspinside(dir)
      newspinside(3-dir) = flip_typ_G(oldspinside(3-dir))
    endif
  endif

  newtyp = TypeSp2Gam(newsp(1),newsp(2),newspinside(1), newspinside(2))

  !call LogFile%QuickLog("flip_typ_Vertex")
  !call LogFile%QuickLog(str(dir)+str(typ))
  !call LogFile%QuickLog(str(newsp(1))+str(newsp(2)))
  !call LogFile%QuickLog(str(newspinside(1))+str(newspinside(2)))

  return
END SUBROUTINE flip_typ_Vertex

LOGICAL FUNCTION is_worm_nearby(statline)
  implicit none
  integer :: statline
  if(statline==2 .or. statline==3) then
    is_worm_nearby = .true.
  else
    is_worm_nearby = .false.
  endif
  return
END FUNCTION is_worm_nearby

INTEGER FUNCTION get_typW_from_Vertex(dir, typGam1, typGam2)
  implicit none
  integer, intent(in) :: dir, typGam2, typGam1
  if(dir==1) then
    get_typW_from_Vertex = TypeGam2W(typGam1, typGam2)
  elseif(dir==2) then
    get_typW_from_Vertex = TypeGam2W(typGam2, typGam1)
  endif
  return
END FUNCTION get_typW_from_Vertex


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
!------               GamD                     GamB     GamD ----
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
  integer :: isdelta, dir, dirW, kIA, kMB, kM, q, kIC, kMD
  integer :: GIC, GMD, GamC, GamD, kNeighG
  integer :: WAB, GIA, GMB, GamA, GamB
  integer :: typGamA, typGamB, typAB
  integer :: statA, statB, statIA, statAC, statMB, statBD
  double precision :: tauA, tauB
  double precision :: tau, tau1, tau2, Pacc
  complex*16 :: WA, WB, WGIA, WGAC, WGMB, WGBD, WWAB, WMeasureGam
  complex*16 :: Anew, Aold, sgn
  logical :: flag
  
  !---------- step1 : check if worm is present ------------------
  if(IsWormPresent .eqv. .false.)    return
  if(Order>=MCOrder) return

  !------------ step2 : propose the new config ------------------
  dir = Floor(rn()*2.d0) + 1
  GIC = NeighVertex(dir, Ira);      GamC = NeighLn(dir, GIC)
  GMD = NeighVertex(dir, Masha);    GamD = NeighLn(dir, GMD)

  dirW = Floor(rn()*2.d0) + 1
  !isdelta = Floor(rn()*2.d0)
  isdelta = 0

  q = generate_k()
  if(Is_k_valid_for_W(q)==.false.) return
  kM = add_k(kMasha, -(-1)**dirW*q)

  kIC = kLn(GIC)
  kIA = add_k(kIC, -(-1)**(dir+dirW)*q)

  kNeighG=kLn(NeighVertex(3-dir, Ira))
  if(Is_Bold .and. abs(add_k(KIA,-kNeighG))==abs(kLn(NeighVertex(3,Ira))))  return
  !This line is used to reject the configuration with |k_GAC-kNeighG|==|k of W attached to Ira|
  !Such configuration will always become reducibile once new interaction line is added
  !Please refer to the comments in the reducibility check codes block in basic_function.f90
  !for further information on the definition of reducibility

  if(Is_k_valid_for_G(kIA)==.false.) return


  kMD = kLn(GMD)
  kMB = add_k(kMD, (-1)**(dir+dirW)*q)

  kNeighG=kLn(NeighVertex(3-dir, Masha))
  if(Is_Bold .and. abs(add_k(KMB,-kNeighG))==abs(kLn(NeighVertex(3,Masha))))  return

  if(Is_k_valid_for_G(kMB)==.false.) return

  if(kIA==kMB)  return       

  !-------- the new time, spin, type and status for the new config ----
  tauA = generate_tau()
  if(isdelta==0) then
    tauB = generate_tau()
  else 
    tauB = tauA
  endif

  typGamA = TypeSp2Gam(TypeLn(GIC), TypeLn(GIC), TypeLn(GIC), TypeLn(GIC))
  typGamB = TypeSp2Gam(TypeLn(GMD), TypeLn(GMD), TypeLn(GMD), TypeLn(GMD))

  if(dirW==1) then
    typAB = TypeGam2W(typGamA, typGamB)
  else
    typAB = TypeGam2W(typGamB, typGamA)
  endif

  statA = 0
  statB = 0
  statIA = gline_stat(statA, StatusVertex(Ira)) 
  statMB = gline_stat(statB, StatusVertex(Masha)) 
  statAC = gline_stat(statA, StatusVertex(GamC))
  statBD = gline_stat(statB, StatusVertex(GamD))

  Order = Order + 1 !!!!! add Order to the next Order

  !-------------- step3 : weight calculation --------------------
  !! TVertex(:, GamA)=tauA
  !! TVertex(:, GamB)=tauB
  WA = weight_vertex(statA, 1, 0, 0.d0, 0.d0, typGamA)  
  WB = weight_vertex(statB, 1, 0, 0.d0, 0.d0, typGamB)  

  tau = (-1)**dir*(tauA - TVertex(3-dir, Ira))
  WGIA = weight_gline(statIA, tau, TypeLn(GIC))

  tau = (-1)**dir*(tauB - TVertex(3-dir, Masha))
  WGMB = weight_gline(statMB, tau, TypeLn(GMD))

  tau = (-1)**dirW*(tauA-tauB)
  WWAB = weight_wline(0, isdelta, diff_r(D, GRVertex(GamC),GRVertex(GamD)), &
    & tau, typAB)
  
  !---  change the topology for the configuration after update --
  call insert_gamma(GamA, 1, GRVertex(GamC), GRVertex(GamC), &
    &  tauA, tauA, tauA, dirW, typGamA, statA, WA)

  call insert_gamma(GamB, 1, GRVertex(GamD), GRVertex(GamD), &
    &  tauB, tauB, tauB, 3-dirW, typGamB, statB, WB)

  call insert_line(GIA, -1, kIA, 1, TypeLn(GIC), statIA, WGIA)
  call insert_line(GMB, -1, kMB, 1, TypeLn(GMD), statMB, WGMB)
  call insert_line(WAB,  isdelta,q, 2,    typAB,      0, WWAB) 
  
  !------------ update the topology -----------------------------
  NeighLn(3-dir, GIC) = GamA;        NeighLn(3-dir, GMD) = GamB
  NeighVertex(dir,Ira) = GIA;        NeighVertex(dir,Masha)=GMB

  NeighVertex(dir, GamA) = GIC;      NeighVertex(dir, GamB) = GMD
  NeighVertex(3-dir,GamA)= GIA;      NeighVertex(3-dir,GamB)= GMB
  NeighVertex(3,   GamA) = WAB;      NeighVertex(3,   GamB) = WAB

  NeighLn(dir, GIA) = GamA;          NeighLn(dir, GMB) = GamB
  NeighLn(3-dir,GIA)= Ira;           NeighLn(3-dir,GMB)= Masha

  NeighLn(dirW, WAB)= GamA;          NeighLn(3-dirW,WAB)= GamB

  !------------ step4 : configuration check ---------------------
  flag=Is_reducible_add_interaction(q, kIA, kIC, kMB, kMD)

  if(flag) then
    Order = Order - 1
    call delete_gamma(GamA)
    call delete_gamma(GamB)
    call undo_insert_line(GIA, 1)
    call undo_insert_line(GMB, 1)
    call undo_insert_line(WAB, 2)

    NeighLn(3-dir, GIC) = Ira;        NeighLn(3-dir, GMD) = Masha
    NeighVertex(dir,Ira) = GIC;        NeighVertex(dir,Masha)=GMD
    return
  endif

  !------------- weight calculation -----------------------------
  tau = (-1)**dir*(TVertex(dir, GamC)-tauA)
  WGAC = weight_gline(statAC, tau, TypeLn(GIC))

  tau = (-1)**dir*(TVertex(dir, GamD)-tauB)
  WGBD = weight_gline(statBD, tau, TypeLn(GMD))

  tau1 = TVertex(3,MeasureGam)-TVertex(2,MeasureGam)
  tau2 = TVertex(1,MeasureGam)-TVertex(3,MeasureGam)
  WMeasureGam = weight_vertex(StatusVertex(MeasureGam), IsDeltaVertex(MeasureGam),  &
    & diff_r(D,GRVertex(MeasureGam),WRVertex(MeasureGam)), tau1, tau2, TypeVertex(MeasureGam))

  Anew = WA *WB *WGIA *WGMB *WWAB *WGAC *WGBD *WMeasureGam
  Anew = -1.d0*Anew

  Aold = WeightLn(GIC) *WeightLn(GMD) *WeightVertex(MeasureGam)

  call weight_ratio(Pacc, sgn, Anew, Aold)
  !Pacc = Pacc *CoefOfWeight(Order)/(0.5d0*prob_tau(tauA)*prob_tau(tauB)*CoefOfWeight(Order-1))
  Pacc = Pacc *CoefOfWeight(Order)/(prob_tau(tauA)*prob_tau(tauB)*CoefOfWeight(Order-1))
  Pacc = Pacc *Pupdate(8)/Pupdate(7)

  !------------ step5 : accept the update -----------------------
  ProbProp(Order-1, 7) = ProbProp(Order-1, 7) + 1
  if(Pacc>1.d-12 .and. rn()<Pacc) then

    !--------------- update the diagram info --------------------
    Phase = Phase *sgn

    !--------------- update k and omega -------------------------
    kMasha = kM

    !--------------- update the status of elements --------------
    StatusLn(GIC) = statAC
    StatusLn(GMD) = statBD

    !--------------- update weight of elements ------------------
    WeightVertex(MeasureGam) = WMeasureGam
    WeightLn(GIC) = WGAC
    WeightLn(GMD) = WGBD

    call update_weight(Anew, Aold)

    ProbAcc(Order-1, 7) = ProbAcc(Order-1, 7) + 1
  else
    !-------------- delete line and vertexes --------------------
    Order = Order - 1
    call delete_gamma(GamA)
    call delete_gamma(GamB)
    call undo_insert_line(GIA, 1)
    call undo_insert_line(GMB, 1)
    call undo_insert_line(WAB, 2)

    NeighLn(3-dir, GIC) = Ira;        NeighLn(3-dir, GMD) = Masha
    NeighVertex(dir,Ira)= GIC;        NeighVertex(dir,Masha)= GMD

    return
  endif
  return
END SUBROUTINE add_interaction





!----------------- remove_interaction: Pupdate(8) ----------------
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
  integer :: dir, dirW, kM, kIA, kMB, kAB, kAC, kBD, kNeighG
  integer :: GAC, GBD, GamC, GamD
  integer :: WAB, GIA, GMB, GamA, GamB
  integer :: statIC, statMD, statIA, statMB, statAB
  double precision :: tau, tau1, tau2, Pacc
  complex*16 :: WGIC, WGMD, WMeasureGam
  complex*16 :: Anew, Aold, sgn
  logical :: flag
  
  !---------- step1 : check if worm is present ------------------
  if(IsWormPresent .eqv. .false.)    return
  !if(Order==1) return
  if(Order==0) return

  !------------ step2 : propose the new config ------------------
  dir = Floor(rn()*2.d0)+1
  GIA = NeighVertex(dir, Ira);          GMB = NeighVertex(dir, Masha)
  GamA = NeighLn(dir, GIA);             GamB = NeighLn(dir, GMB)

  if(NeighVertex(3, GamA)/=NeighVertex(3, GamB))  return
  if(StatusVertex(GamA)/=0 .or. StatusVertex(GamB)/=0) return
  if(TypeVertex(GamA)>2 .or. TypeVertex(GamB)>2) return
  if(IsDeltaVertex(GamA)==0 .or. IsDeltaVertex(GamB)==0)  return
  if(GRVertex(GamA)/=WRVertex(GamA)) return
  if(GRVertex(GamB)/=WRVertex(GamB)) return

  WAB = NeighVertex(3, GamA)
  if(IsDeltaLn(WAB)==1)  return
  GAC = NeighVertex(dir, GamA);         GBD = NeighVertex(dir, GamB)
  GamC = NeighLn(dir, GAC);             GamD = NeighLn(dir, GBD)

  dirW = DirecVertex(GamA)
  kAB=kLn(WAB)
  kM = add_k(kMasha, (-1)**dirW*kAB)

  kIA = kLn(GIA)
  kAC = kLn(GAC)

  kMB = kLn(GMB)
  kBD = kLn(GBD)

  statIA = StatusLn(GIA)
  statMB = StatusLn(GMB)
  statAB = StatusLn(WAB)

  !------------ update the topology -----------------------------
  Order = Order - 1
  call delete_gamma(GamA)
  call delete_gamma(GamB)
  call delete_line(GIA, 1)
  call delete_line(GMB, 1)
  call delete_line(WAB, 2)

  NeighLn(3-dir, GAC) = Ira;            NeighLn(3-dir, GBD) = Masha
  NeighVertex(dir, Ira) = GAC;          NeighVertex(dir, Masha)=GBD

  !-------- the new status for the new config ----
  statIC = gline_stat(StatusVertex(Ira), StatusVertex(GamC))
  statMD = gline_stat(StatusVertex(Masha),StatusVertex(GamD))

  !-------------- step4 : weight calculation --------------------
  tau = (-1)**dir *(TVertex(dir, GamC)-TVertex(3-dir, Ira))
  WGIC = weight_gline(statIC, tau, TypeLn(GAC))
  tau = (-1)**dir *(TVertex(dir, GamD)-TVertex(3-dir, Masha))
  WGMD = weight_gline(statMD, tau, TypeLn(GBD))

  tau1 = TVertex(3,MeasureGam)-TVertex(2,MeasureGam)
  tau2 = TVertex(1,MeasureGam)-TVertex(3,MeasureGam)
  WMeasureGam = weight_vertex(StatusVertex(MeasureGam), IsDeltaVertex(MeasureGam), &
    & diff_r(D,GRVertex(MeasureGam),WRVertex(MeasureGam)), &
    & tau1, tau2, TypeVertex(MeasureGam))

  Anew = WGIC *WGMD *WMeasureGam
  Aold = WeightVertex(GamA)*WeightVertex(GamB)*WeightLn(GIA)*WeightLn(GMB)*WeightLn(WAB)* &
    & WeightLn(GAC) *WeightLn(GBD) *WeightVertex(MeasureGam)
  Aold = (-1.d0)*Aold

  call weight_ratio(Pacc, sgn, Anew, Aold)

  !Pacc = Pacc *CoefOfWeight(Order)*0.5d0*prob_tau(TVertex(1,GamA))*prob_tau(TVertex(1,GamB))/ &
    !& CoefOfWeight(Order+1)
  Pacc = Pacc *CoefOfWeight(Order)*prob_tau(TVertex(1,GamA))*prob_tau(TVertex(1,GamB))/ &
    & CoefOfWeight(Order+1)
  Pacc = Pacc *Pupdate(7)/Pupdate(8)

  !------------ step5 : accept the update -----------------------
  ProbProp(Order+1, 8) = ProbProp(Order+1, 8) + 1
  if(Pacc>1.d-12 .and. rn()<Pacc) then

    !--------------- update the diagram info --------------------
    Phase = Phase *sgn

    !--------------- update k -----------------------------------
    kMasha = kM

    !--------------- update the status of elements --------------
    StatusLn(GAC) = statIC
    StatusLn(GBD) = statMD

    !--------------- update weight of elements ------------------
    WeightVertex(MeasureGam) = WMeasureGam
    WeightLn(GAC) = WGIC
    WeightLn(GBD) = WGMD

    call update_weight(Anew, Aold)

    ProbAcc(Order+1, 8) = ProbAcc(Order+1, 8) + 1
  else
    !-------------- undo delete line and vertexes --------------------
    Order = Order + 1
    call undo_delete_gamma(GamB)
    call undo_delete_gamma(GamA)
    call undo_delete_line(WAB, 2, statAB, kAB)
    call undo_delete_line(GMB, 1, statMB, kMB)
    call undo_delete_line(GIA, 1, statIA, kIA)

    NeighLn(3-dir, GAC) = GamA;        NeighLn(3-dir, GBD) = GamB
    NeighVertex(dir, Ira) = GIA;       NeighVertex(dir, Masha)=GMB
    return
  endif

END SUBROUTINE remove_interaction

!----------------- reconnect: Pupdate(11) ----------------
!------------------------- dir = 1 ------------------------------
!------  Ira           GamA            Ira            GamA   -----
!------    ------<-------                  --<-\    /-<--    -----
!------                                         \  /         -----
!------                           ==>            \/          -----
!------    ------<-------                  --<---/\---<---   -----
!------  Masha         GamB            Masha          GamB   -----
!-----------------------------------------------------------------
!------------------------- dir = 2 ------------------------------
!------  Ira           GamA            Ira            GamA   -----
!------    ------>-------                  -->-\    /->--    -----
!------                                         \  /         -----
!------                           ==>            \/          -----
!------    ------>-------                  -->---/\--->---   -----
!------  Masha         GamB            Masha          GamB   -----
!-----------------------------------------------------------------
SUBROUTINE reconnect
  implicit none
  integer :: GamA,GamB,dir
  integer :: GIA,GMB, kNeighG
  integer :: statIA,statMB
  complex*16 :: Anew, Aold, sgn
  complex*16 :: WGIA,WGMB
  double precision :: tau, Pacc
  logical :: flag
  !---------- step1 : check if worm is present ------------------
  if(IsWormPresent .eqv. .false.)    return
  if(GRVertex(Ira)/=GRVertex(Masha))  return 

  dir=Floor(rn()*2.d0)+1
  GIA=NeighVertex(dir,Ira)
  GMB=NeighVertex(dir,Masha)
  if(TypeLn(GIA)/=TypeLn(GMB))return

  GamA=NeighLn(dir,GIA)
  GamB=NeighLn(dir,GMB)

  !------------ step2 : propose the new config ------------------
  NeighLn(3-dir,GIA)=Masha
  NeighVertex(dir, Masha) = GIA
  NeighLn(3-dir,GMB)=Ira
  NeighVertex(dir, Ira) = GMB

  !Don't have changed delta_k of Ira and Masha yet here

  !------------ step4 : configuration check ---------------------
  ! do the step4 here so we can save some time

  !-------- the new spin, type and status for the new config ----
  statIA = gline_stat(StatusVertex(GamA), StatusVertex(Masha))
  statMB = gline_stat(StatusVertex(GamB), StatusVertex(Ira))

  !-------------- step3 : weight calculation --------------------
  tau = (-1)**dir *(TVertex(dir, GamA)-TVertex(3-dir, Masha))
  WGIA=weight_gline(statIA, tau,TypeLn(GIA))

  tau = (-1)**dir *(TVertex(dir, GamB)-TVertex(3-dir, Ira))
  WGMB=weight_gline(statMB, tau,TypeLn(GMB))

  Anew = d_times_cd(-1.d0, WGIA *WGMB)
  Aold = WeightLn(GIA)*WeightLn(GMB)

  call weight_ratio(Pacc, sgn, Anew, Aold)

  !------- step5 : accept the update --------------------
  ProbProp(Order, 11) = ProbProp(Order, 11) + 1
  if(Pacc>1.d-12 .and. rn()<Pacc) then

    !--------------- update the diagram info --------------------
    SignFermiloop=-SignFermiloop
    Phase = Phase *sgn

    !--------------- update k and omega -------------------------
    kMasha = add_k(kMasha, (-1)**dir*(add_k(kLn(GMB),-kLn(GIA))))

    !--------------- update the status of elements --------------
    StatusLn(GIA) = statIA
    StatusLn(GMB) = statMB

    !--------------- update weight of elements ------------------
    WeightLn(GIA)=WGIA
    WeightLn(GMB)=WGMB

    call update_weight(Anew, Aold)

    ProbAcc(Order, 11) = ProbAcc(Order, 11) + 1
  else
    !-------------- delete line and vertexes --------------------
    NeighLn(3-dir, GIA)=Ira
    NeighVertex(dir, Ira)=GIA
    NeighLn(3-dir, GMB)=Masha
    NeighVertex(dir, Masha)=GMB
  endif
  return

END SUBROUTINE reconnect




!------------ change gline space: Pupdate(12) ------------
SUBROUTINE change_gline_space
  implicit none
  integer :: Num,InitialGam,iGam,jGam,iWLn,NewRG,dR
  integer :: GamList(MxNVertex),NewRW(MxNVertex),i
  integer :: flagW(MxNLn)
  double precision :: Pacc, WeightR, WeightRW, tau
  complex*16  ::  Anew, Aold, sgn, WeightW(MxNVertex)

  !------- step1 : check if worm is present -------------
  if(IsWormPresent .eqv. .true.)    return

  !------- step2 : propose a new config -----------------
  iGam=generate_vertex()
  InitialGam=iGam

  WeightR=1.d0

  call generate_r(D,GRVertex(iGam),NewRG,dR,WeightR,.true.)

  flagW(:) = 0
  Num=0
  do while(.true.)
    Num=Num+1
    iWLn=NeighVertex(3,iGam)

    flagW(iWLn) = flagW(iWLn) + 1

    GamList(Num)=iGam
    iGam=NeighLn(2,NeighVertex(2,iGam))
    if(InitialGam==iGam) exit
  enddo

  AOld=(1.d0, 0.d0)
  ANew=(1.d0, 0.d0)
  do i = 1, Num
    iGam = GamList(i)
    iWLn=NeighVertex(3,iGam)
    WeightRW = 1.d0

    call generate_r(D,WRVertex(iGam),NewRW(i),dR,WeightRW,.false.)

    if(flagW(iWLn)==1) then
      jGam=NeighLn(3-DirecVertex(iGam),iWLn)
      tau = (-1)**DirecVertex(iGam)*(TVertex(3, iGam)-TVertex(3, jGam))
      WeightW(i)=weight_wline(StatusLn(iWLn),IsDeltaLn(iWLn),&
          & diff_r(D, NewRW(i),WRVertex(jGam)),tau,TypeLn(iWLn))

      Anew=Anew*WeightW(i)
      Aold=Aold*WeightLn(iWLn)
      if(abs(Anew)<macheps*abs(Aold)) return
    endif
  enddo

  call weight_ratio(Pacc, sgn, Anew, Aold)

  if(WeightR<1.d-12)  return
  Pacc=Pacc/WeightR

  !------- step5 : accept the update --------------------
  ProbProp(Order, 12) = ProbProp(Order, 12) + 1
  if(Pacc>1.d-12 .and. rn()<Pacc) then

    !------ update the diagram info -------------------
    Phase = Phase *sgn

    do i = 1, Num
      iGam=GamList(i)
      iWLn=NeighVertex(3,iGam)

      !------ update the site of elements -------------
      GRVertex(iGam) = NewRG
      WRVertex(iGam) = NewRW(i)

      !------ update the weight of elements -------------
      if(flagW(iWLn)==1) then
        WeightLn(iWLn)  = WeightW(i)
      endif
    enddo

    call update_weight(Anew, Aold)
    ProbAcc(Order, 12) = ProbAcc(Order, 12) + 1
  endif
  return
  
END SUBROUTINE


!------------ change wline space: Pupdate(13) ------------
SUBROUTINE change_wline_space
    implicit none
    integer :: iGam, iWLn, jGam, rwi, rwj, dir,dr
    double precision :: Pacc,T1,T2,T3,T4,T5,T6,WeightR
    complex*16  ::  WiGam,WjGam, WW, Anew, Aold, sgn

    !------- step1 : check if worm is present -------------
    if(IsWormPresent .eqv. .true.)    return

    !------- step2 : propose a new config -----------------
    iWLn=generate_wline()
    iGam = NeighLn(1, iWLn)
    jGam = NeighLn(2, iWLn)

    WeightR=1.0
    call generate_r(D,WRVertex(iGam),rwi,dr,WeightR,.true.)

    if(IsDeltaLn(iWLn)==0) then
      call generate_r(D,WRVertex(jGam),rwj,dr,WeightR,.true.)
    else
      rwj=rwi
    endif

    !------- step3 : configuration check ------------------

    !------- step4 : weight calculation -------------------
    T1=TVertex(1, iGam);
    T2=TVertex(2, iGam);
    T3=TVertex(3, iGam);
    WiGam = weight_vertex(StatusVertex(iGam),IsDeltaVertex(iGam), &
      & diff_r(D,GRVertex(iGam),rwi), T3-T2, T1-T3, TypeVertex(iGam))

    T4=TVertex(1, jGam);
    T5=TVertex(2, jGam);
    T6=TVertex(3, jGam);
    WjGam = weight_vertex(StatusVertex(jGam),IsDeltaVertex(jGam), &
      & diff_r(D,GRVertex(jGam),rwj), T6-T5, T4-T6, TypeVertex(jGam))

    WW = weight_wline(StatusLn(iWLn),IsDeltaLn(iWLn),diff_r(D,rwi,rwj),T6-T3,TypeLn(iWLn))

    Anew = WiGam*WjGam*WW
    Aold = WeightLn(iWLn)*WeightVertex(iGam)*WeightVertex(jGam)

    call weight_ratio(Pacc, sgn, Anew, Aold)
    if(WeightR<1.d-12)  return
    Pacc = Pacc/WeightR

    !------- step5 : accept the update --------------------
    ProbProp(Order, 13) = ProbProp(Order, 13) + 1
    if(Pacc>1.d-12 .and. rn()<Pacc) then

      !------ update the diagram info -------------------
      Phase = Phase *sgn

      !------ update the site of elements --------------
      WRVertex(iGam) = rwi
      WRVertex(jGam) = rwj

      !------ update the weight of elements ------------
      WeightLn(iWLn) = WW
      WeightVertex(iGam) = WiGam
      WeightVertex(jGam) = WjGam

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

  !------- step2 : propose a new config -----------------
  iGam = generate_vertex()
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
    typW = TypeGam2W(typGam, TypeVertex(jGam))
  else
    typW = TypeGam2W(TypeVertex(jGam), typGam)
  endif

  !------- step4 : weight calculation -------------------
  tau1 = TVertex(3, iGam)-TVertex(2, iGam)
  tau2 = TVertex(1, iGam)-TVertex(3, iGam)
  WGam = weight_vertex(StatusVertex(iGam), IsDeltaVertex(iGam), &
    & diff_r(D, GRVertex(iGam),WRVertex(iGam)), tau1, tau2, typGam)

  tau = TVertex(3, NeighLn(2, iWLn))-TVertex(3, NeighLn(1, iWLn))
  WW  = weight_wline(StatusLn(iWLn), IsDeltaLn(iWLn), diff_r(D,WRVertex(iGam),WRVertex(jGam)), &
    & tau, typW) 

  Anew = WW *WGam
  Aold = WeightLn(iWLn)*WeightVertex(iGam)

  call weight_ratio(Pacc, sgn, Anew, Aold)

  !------- step5 : accept the update --------------------
  ProbProp(Order, 14) = ProbProp(Order, 14) + 1

  if(Pacc>1.d-12 .and. rn()<Pacc) then

    if(TypeVertex(iGam)==1 .or. TypeVertex(iGam)==2) then
      BalenceCheck(Order,1,1)=BalenceCheck(Order,1,1)+1.d0
    else if(TypeVertex(iGam)==3 .or. TypeVertex(iGam)==4) then
      BalenceCheck(Order,1,2)=BalenceCheck(Order,1,2)+1.d0
    endif

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

  !------- step2 : propose a new config -----------------
  jGam  = MeasureGam
  jGin  = NeighVertex(1, jGam);      jGout = NeighVertex(2, jGam)
  jW    = NeighVertex(3, jGam)

  iGam = generate_vertex()
  if(iGam==MeasureGam)         return
  if(IsDeltaVertex(iGam)/=1)   return

  iGin  = NeighVertex(1, iGam);      iGout = NeighVertex(2, iGam)
  iW    = NeighVertex(3, iGam)

  if(IsDeltaLn(iW)/=0)         return

  !----- the status for the new config ------------------
  statiGam  = add_mea_stat(StatusVertex(iGam))
  statjGam  = delete_mea_stat(StatusVertex(jGam))

  if(jW/=iW) then
    statiW    = wline_stat(statiGam, StatusVertex(NeighLn(3-DirecVertex(iGam), iW)))
    statjW    = wline_stat(statjGam, StatusVertex(NeighLn(3-DirecVertex(jGam), jW)))
  else
    statiW    = wline_stat(statiGam, statjGam)
    statjW    = wline_stat(statjGam, statiGam)
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
 !dx = xg-xw;  dy = yg-yw; dtau1 = tau3-tau2; dtau2 = tau1-tau3
  tau1 = TVertex(3, iGam)-TVertex(2, iGam)
  tau2 = TVertex(1, iGam)-TVertex(3, iGam)
!   if(abs(tau1)>1e-8 .or. abs(tau2)>1e-8) then
!     print *," I am a Bug1!"
!     stop
!   endif
  WiGam = weight_vertex(statiGam, IsDeltaVertex(iGam), diff_r(D,GRVertex(iGam),WRVertex(iGam)), &
    & tau1, tau2, TypeVertex(iGam))

  tau = TVertex(3, NeighLn(2, iW)) - TVertex(3, NeighLn(1, iW))
  WiW = weight_wline(statiW, IsDeltaLn(iW), diff_r(D,WRVertex(NeighLn(2, iW)), &
    & WRVertex(NeighLn(1,iW))), tau, TypeLn(iW))

  tau = TVertex(2, NeighLn(2, iGin)) - TVertex(1, NeighLn(1, iGin))
  WiGin = weight_gline(statiGin, tau, TypeLn(iGin))

  tau = TVertex(2, NeighLn(2, iGout)) - TVertex(1, NeighLn(1, iGout))
  WiGout = weight_gline(statiGout, tau, TypeLn(iGout))

  tau1 = TVertex(3, jGam)-TVertex(2, jGam)
  tau2 = TVertex(1, jGam)-TVertex(3, jGam)
!   if(abs(tau1)>1e-8 .or. abs(tau2)>1e-8) then
!     print *," I am a Bug2!"
!     stop
!   endif
  WjGam = weight_vertex(statjGam, IsDeltaVertex(jGam), diff_r(D,GRVertex(jGam),WRVertex(jGam)),  &
    & tau1, tau2, TypeVertex(jGam))

  tau = TVertex(3, NeighLn(2, jW)) - TVertex(3, NeighLn(1, jW))
  WjW = weight_wline(statjW, IsDeltaLn(jW), diff_r(D,WRVertex(NeighLn(2, jW)),WRVertex(NeighLn(1,jW))),&
    &  tau, TypeLn(jW))

  tau = TVertex(2, NeighLn(2, jGin)) - TVertex(1, NeighLn(1, jGin))
  WjGin = weight_gline(statjGin, tau, TypeLn(jGin))

  tau = TVertex(2, NeighLn(2, jGout)) - TVertex(1, NeighLn(1, jGout))
  WjGout = weight_gline(statjGout, tau, TypeLn(jGout))

  Anew = WjGam*WiGam *WiW *WiGin *WiGout
  Aold = WeightVertex(iGam) *WeightLn(iW) *WeightLn(iGin) *WeightLn(iGout) &
    & *WeightVertex(jGam)

  if(abs(Anew)==0.d0)   return

!   if(iGin==iGout .or. jGin==jGout) then
!     call LogFile%QuickLog("a buble diagram!"+str(imc))
!     call print_config
!     stop
!   endif

  if(iGin/=jGout) then
    Anew = Anew*WjGout
    Aold = Aold*WeightLn(jGout)
  endif

  if(jGin/=iGout) then
    Anew = Anew*WjGin
    Aold = Aold*WeightLn(jGin)
  endif

  if(iW/=jW) then
    Anew = Anew*WjW
    Aold = Aold*WeightLn(jW)
  endif

  call weight_ratio(Pacc, sgn, Anew, Aold)

  !------- step5 : accept the update --------------------
  ProbProp(Order, 15) = ProbProp(Order, 15) + 1

  if(Pacc>1.d-12 .and. rn()<Pacc) then

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



!------------- change gamma time : Pupdate(16) -----------------
SUBROUTINE change_Gamma_time
  implicit none
  integer :: isdelta
  integer :: iLn, iLn1, iLn2, iLn3, iGam, jGam, dir, dirGam
  double precision :: Pacc, tau, tau1, tau2, newtau
  complex*16 :: WGam, WLn, WLn1, WLn2, WLn3, Anew, Aold, sgn

  !------- step1 : check if worm is present -------------
  if(IsWormPresent .eqv. .true.)    return

  !------- step2 : propose a new config -----------------
  iGam=generate_vertex()
  dir = Floor(rn()*3.d0) + 1  !dir=1: in;  dir=2: out; dir=3 : w
  iLn = NeighVertex(dir, iGam)

  isdelta=0
  if(IsDeltaVertex(iGam)==1) isdelta=1

  if(isdelta==0) then
    if(dir==3 .and. IsDeltaLn(iLn)==1)  return
  else
    if(IsDeltaLn(NeighVertex(3, iGam))==1) return
  endif

  newtau=generate_tau()

  !------- step3 : weight calculation -------------------
  if(isdelta==0) then
    if(dir==1) then
      tau1 = TVertex(3, iGam) - newtau
      tau2 = TVertex(1, iGam) - TVertex(3, iGam)
    else if(dir==2) then
      tau1 = TVertex(3, iGam) - TVertex(2, iGam)
      tau2 = newtau - TVertex(3, iGam)
    else if(dir==3) then
      tau1 = newtau - TVertex(2, iGam)
      tau2 = TVertex(1, iGam) - newtau
    endif
  else
    tau1 = 0
    tau2 = 0
  endif

  WGam = weight_vertex(StatusVertex(iGam),IsDeltaVertex(iGam),diff_r(D,GRVertex(iGam) &
    & ,WRVertex(iGam)), tau1, tau2, TypeVertex(iGam))

  if(isdelta==0) then
    if(dir==1 .or. dir==2) then
      jGam = NeighLn(dir, iLn)
      tau = (-1)**dir *(TVertex(dir, jGam) - newtau)
      WLn = weight_gline(StatusLn(iLn), tau, TypeLn(iLn))
    else if(dir==3) then
      dirGam = DirecVertex(iGam)
      jGam = NeighLn(3-dirGam, iLn)
      tau = (-1)**dirGam*(newtau -TVertex(3, jGam))
      WLn = weight_wline(StatusLn(iLn),IsDeltaLn(iLn), diff_r(D, WRVertex(iGam),WRVertex(jGam)), &
        & tau, TypeLn(iLn))
    endif
    Anew = WGam*WLn
    Aold = WeightLn(iLn)*WeightVertex(iGam)

  else if(isdelta==1) then
    iLn1 = NeighVertex(1, iGam)
    jGam = NeighLn(1, iLn1)
    tau = newtau - TVertex(1, jGam)
    WLn1 = weight_gline(StatusLn(iLn1), tau, TypeLn(iLn1))

    iLn2 = NeighVertex(2, iGam)
    jGam = NeighLn(2, iLn2)
    tau = TVertex(2, jGam) - newtau
    WLn2 = weight_gline(StatusLn(iLn2), tau, TypeLn(iLn2))

    iLn3 = NeighVertex(3, iGam)
    dirGam = DirecVertex(iGam)
    jGam = NeighLn(3-dirGam, iLn3)
    tau = (-1)**dirGam*(newtau -TVertex(3, jGam))
    WLn3 = weight_wline(StatusLn(iLn3),IsDeltaLn(iLn3), diff_r(D, WRVertex(iGam),WRVertex(jGam)), &
      & tau, TypeLn(iLn3))

    Anew = WGam *WLn1 *WLn2 *WLn3
    Aold = WeightVertex(iGam) *WeightLn(iLn1) *WeightLn(iLn2) *WeightLn(iLn3)
  endif

  call weight_ratio(Pacc, sgn, Anew, Aold)

  !------- step5 : accept the update --------------------
  ProbProp(Order, 16) = ProbProp(Order, 16) + 1

  if(Pacc>1.d-12 .and. rn()<Pacc) then

    !------ update the diagram info ---------------------
    Phase = Phase *sgn
    call update_weight(Anew, Aold)

    !------ update the time and weight of elements ------
    if(isdelta==0) then
      if(dir<3) then
        TVertex(3-dir, iGam) = newtau
      else 
        TVertex(3, iGam) = newtau
      endif

      WeightLn(iLn) = WLn
    else 
      TVertex(:, iGam) = newtau
      WeightLn(iLn1) = WLn1
      WeightLn(iLn2) = WLn2
      WeightLn(iLn3) = WLn3
    endif

    WeightVertex(iGam) = WGam

    ProbAcc(Order, 16) = ProbAcc(Order, 16) + 1
  endif
  return
END SUBROUTINE change_Gamma_time




!------------ change wline isdelta: Pupdate(17) ------------
SUBROUTINE change_wline_isdelta
  implicit none
  integer :: iWLn, iGam, jGam, jGin, jGout, backforth
  !backforth = 1: delta -> normal; backforth = 0: normal -> delta
  double precision :: t3iGam, t3jGam, t1jGam, t2jGam, tau1, tau2, dtau, Pacc
  complex*16 :: WW, WjGin, WjGout, WjGam, Anew, Aold, sgn

  !------- step1 : check if worm is present -------------
  if(IsWormPresent .eqv. .true.)    return

  !------- step2 : propose a new config -----------------
  iWLn = generate_wline()
  if(StatusLn(iWLn)==1)  return

  iGam = NeighLn(1, iWLn)
  jGam = NeighLn(2, iWLn)
  jGin = NeighVertex(1, jGam)
  jGout = NeighVertex(2, jGam)

  t3iGam = TVertex(3, iGam)

  t1jGam = TVertex(1, jGam)
  t2jGam = TVertex(2, jGam)
  t3jGam = TVertex(3, jGam)

  backforth = IsDeltaLn(iWLn)

  if(backforth==0) then
    t3jGam = t3iGam
  else 
    t3jGam = generate_tau()
  endif

  if(IsDeltaVertex(jGam)==1) then
    t1jGam = t3jGam
    t2jGam = t3jGam
  endif
  
  !------- step4 : weight calculation -------------------
  if(backforth==0) then
    WW = weight_wline(StatusLn(iWLn), 1, diff_r(D, WRVertex(jGam),WRVertex(iGam)), &
      & 0.d0, TypeLn(iWLn)) 
  else
    dtau = t3jGam - t3iGam
    WW = weight_wline(StatusLn(iWLn), 0, diff_r(D, WRVertex(jGam),WRVertex(iGam)), &
      & dtau, TypeLn(iWLn)) 
  endif

  tau1 = t3jGam-t2jGam
  tau2 = t1jGam-t3jGam
  WjGam = weight_vertex(StatusVertex(jGam), IsDeltaVertex(jGam), diff_r(D, GRVertex(jGam), &
    & WRVertex(jGam)), tau1, tau2, TypeVertex(jGam))

  dtau = t2jGam-TVertex(1, NeighLn(1, jGin))
  WjGin = weight_gline(StatusLn(jGin), dtau, TypeLn(jGin)) 

  dtau = TVertex(2, NeighLn(2, jGout))-t1jGam
  WjGout = weight_gline(StatusLn(jGout), dtau, TypeLn(jGout)) 

  Anew = WjGam *WW *WjGin *WjGout
  Aold = WeightLn(iWLn) *WeightVertex(jGam) *WeightLn(jGin) *WeightLn(jGout)

  call weight_ratio(Pacc, sgn, Anew, Aold)
  if(backforth==0) then
    Pacc = Pacc*prob_tau(TVertex(3, jGam))
  else 
    Pacc = Pacc/prob_tau(t3jGam)
  endif

  !------- step5 : accept the update --------------------
  ProbProp(Order, 17) = ProbProp(Order, 17) + 1

  if(Pacc>1.d-12 .and. rn()<Pacc) then

    !------ update the diagram info -------------------
    Phase = Phase *sgn

    !------- update the weight type of elements -------
    IsDeltaLn(iWLn) = 1-backforth

    !------ update the time of elements ---------------
    TVertex(3, iGam) = t3iGam
    TVertex(1, jGam) = t1jGam
    TVertex(2, jGam) = t2jGam
    TVertex(3, jGam) = t3jGam

    !------ update the weight of elements -------------
    WeightLn(iWLn) = WW
    WeightVertex(jGam) = WjGam
    WeightLn(jGin) = WjGin
    WeightLn(jGout) = WjGout

    call update_weight(Anew, Aold)

    ProbAcc(Order, 17) = ProbAcc(Order, 17) + 1
  endif
  return
END SUBROUTINE change_wline_isdelta
  


!------------ change gamma isdelta: Pupdate(18) ------------
SUBROUTINE change_gamma_isdelta
  implicit none
  integer :: iGam, iGin, iGout, iW, backforth
  !backforth = 1: delta -> normal; backforth = 0: normal -> delta
  double precision :: t1iGam, t2iGam, t3iGam, tau1, tau2, dtau, Pacc
  complex*16 :: WiGin, WiGout, WiGam, Anew, Aold, sgn

  !------- step1 : check if worm is present -------------
  if(IsWormPresent .eqv. .true.)   return

  !------- step2 : propose a new config -----------------
  iGam = generate_vertex()
  if(iGam==MeasureGam)  return
  if(TypeVertex(iGam)==3 .or. TypeVertex(iGam)==4)  return
   
  backforth = IsDeltaVertex(iGam)

  iGin = NeighVertex(1, iGam)
  iGout = NeighVertex(2, iGam)

  t1iGam = TVertex(1, iGam)
  t2iGam = TVertex(2, iGam)
  t3iGam = TVertex(3, iGam)


  if(backforth==0) then
    t1iGam = t3iGam
    t2iGam = t3iGam
  else 
    t1iGam = generate_tau()
    t2iGam = generate_tau()
  endif

  !------- step4 : weight calculation -------------------
  tau1 = t3iGam-t2iGam
  tau2 = t1iGam-t3iGam
  WiGam = weight_vertex(StatusVertex(iGam), 1-backforth, diff_r(D, GRVertex(iGam),WRVertex(iGam)), &
    & tau1, tau2, TypeVertex(iGam))

  dtau = t2iGam-TVertex(1, NeighLn(1, iGin))
  WiGin = weight_gline(StatusLn(iGin), dtau, TypeLn(iGin)) 

  dtau = TVertex(2, NeighLn(2, iGout))-t1iGam
  WiGout = weight_gline(StatusLn(iGout), dtau, TypeLn(iGout)) 

  Anew = WiGam *WiGin *WiGout
  Aold = WeightVertex(iGam) *WeightLn(iGin) *WeightLn(iGout)

  call weight_ratio(Pacc, sgn, Anew, Aold)

  if(backforth==0) then
    Pacc = Pacc*prob_tau(TVertex(1, iGam))*prob_tau(TVertex(2, iGam))
  else 
    Pacc = Pacc/(prob_tau(t1iGam)*prob_tau(t2iGam))
  endif

  !------- step5 : accept the update --------------------
  ProbProp(Order, 18) = ProbProp(Order, 18) + 1

  if(Pacc>1.d-12 .and. rn()<Pacc) then

    !------ update the diagram info -------------------
    Phase = Phase *sgn

    !------- update the weight type of elements -------
    IsDeltaVertex(iGam) = 1-backforth

    !------ update the time of elements ---------------
    TVertex(1, iGam) = t1iGam
    TVertex(2, iGam) = t2iGam
    TVertex(3, iGam) = t3iGam

    !------ update the weight of elements -------------
    WeightVertex(iGam) = WiGam
    WeightLn(iGin) = WiGin
    WeightLn(iGout) = WiGout

    call update_weight(Anew, Aold)

    ProbAcc(Order, 18) = ProbAcc(Order, 18) + 1
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

!-------- the weight ratio of new/old config -------------------------
SUBROUTINE weight_ratio(Pacc, sgn, Anew, Aold)
  implicit none
  complex*16,intent(in) ::Anew, Aold
  double precision,intent(out) :: Pacc
  complex*16, intent(out) :: sgn
  Pacc = abs(Anew/Aold)
  if(Pacc>=1.d-12) then
    sgn = (Anew/Aold)/Pacc
  else
    Pacc = 0.d0
    sgn = (1.d0, 0.d0)
  endif
  return
END SUBROUTINE weight_ratio

SUBROUTINE update_weight(Anew, Aold)
  implicit none 
  complex*16 :: Anew, Aold

  !if(abs(Anew)<=1.d-12)  then
    !call LogFile%WriteStamp('e')
    !call LogFile%WriteLine("the weight for the new conf is too small, should not accept!")
    !call LogFile%WriteLine("the new weight :"+str(Anew))
    !call LogFile%WriteLine("the previous weight :"+str(Aold))
    !call print_config
    !stop
  !endif

  WeightCurrent = WeightCurrent *abs(Anew/Aold)
  return
END SUBROUTINE update_weight
!====================================================================
!====================================================================
!====================================================================


!!=======================================================================
!!============================= MEASURE =================================
!!=======================================================================
SUBROUTINE measure
  implicit none
  integer :: i, iln,iGam,jGam, iW
  integer :: it
  integer :: ibasis, ibin
  integer :: spg, spw
  integer :: ityp, nloop
  integer :: MeaGin, MeaGout, MeaW, rg, rw, dir, typ
  integer :: dr, dt1, dt2
  integer :: dx, dy
  integer :: ikey, sumt, sumd
  double precision  :: factorM
  double precision :: tau1, tau2, tau3
  double precision :: dtau1, dtau2
  double precision :: wbasis
  logical :: flag

  GamWormOrder(Order) = GamWormOrder(Order) + 1.d0
  !===========  Measure in worm space  ==========================
  if(IsWormPresent) then
    Z_worm=Z_worm+1.d0
  endif


  !=============================================================
  if(.not. IsWormPresent) then
    !===========  Measure in normal space  ==========================
    Z_normal=Z_normal+1.d0
    GamOrder(Order) = GamOrder(Order) + 1.d0

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

    if(typ==1 .or. typ==2) then
      ityp = 1
    else if(typ==3 .or. typ==4) then
      ityp = 2
    else if(typ==5 .or. typ==6) then
      ityp = 3
    endif

    !----- find the space and time variables for Gamma -------
    rg = GRVertex(MeasureGam)
    rw = WRVertex(NeighLn(3-dir, MeaW))
    dr = diff_r(D, rg, rw)

    tau1 = TVertex(1, NeighLn(1, MeaGout))
    tau2 = TVertex(2, NeighLn(2, MeaGin))
    tau3 = TVertex(3, NeighLn(3-dir, MeaW))

    factorM = 1.d0
    dtau1 = tau3-tau2
    if(dtau1<0) then
      dtau1 = dtau1 + Beta
      factorM = factorM * (-1.d0)
    endif

    dtau2 = tau1-tau3
    if(dtau2<0) then
      dtau2 = dtau2 + Beta
      factorM = factorM * (-1.d0)
    endif

    dt1 = Floor(dtau1*MxT/Beta)
    dt2 = Floor(dtau2*MxT/Beta)
    factorM = factorM *CoefOfSymmetry(dr)* CoefOfWeight(Order) *abs(WeightVertex(MeasureGam))

    !================= accumulation ===================================
    if(Order==0 .and. IsDeltaVertex(NeighLn(3-dir, MeaW))==1) then
      GamNorm = GamNorm + Phase/CoefOfWeight(0)
    endif

    !============ save the diagonal Gamma ================================
    if(dt1==dt2 .and. ityp==1) then
      GamMC(Order, dr,  dt1) = GamMC(Order, dr, dt1) + (Phase/factorM)
      ReGamSqMC(Order, dr, dt1) = ReGamSqMC(Order, dr, dt1) + (Real(Phase/factorM))**2.d0
      ImGamSqMC(Order, dr, dt1) = ImGamSqMC(Order, dr, dt1) + (dimag(Phase/factorM))**2.d0
    endif

    !============ save Gamma in the fitting coeffecients =====================================
    !flag = .false.
    !if(dt1+dt2>MxT-1) then
      !flag= .true.
      !dt1 = MxT-1-dt1
      !dt2 = MxT-1-dt2
    !endif

    ibin = get_bin_Gam(dt1, dt2)

    if(IsBasis2D(ibin)) then
      do ibasis = 1, NBasisGam
        wbasis = (Beta/dble(MxT))**2.d0*weight_basis_Gam(CoefGam &
          & (0:BasisOrderGam,0:BasisOrderGam, ibasis,ibin), dt1*Beta/MxT, dt2*Beta/MxT)

        GamBasis(Order, ityp, dr, ibin, ibasis) = GamBasis(Order, ityp, dr, ibin, &
          & ibasis) + (Phase/factorM)*wbasis

        ReGamSqBasis(Order, ityp, dr, ibin, ibasis) = ReGamSqBasis(Order, ityp, dr, ibin, &
          & ibasis) + (real(Phase)/factorM)**2.d0*wbasis
        ImGamSqBasis(Order, ityp, dr, ibin, ibasis) = ImGamSqBasis(Order, ityp, dr, ibin, &
          & ibasis) + (dimag(Phase)/factorM)**2.d0*wbasis
      enddo
    else 
      do ibasis = 1, NBasis
        wbasis = (Beta/dble(MxT))*weight_basis(CoefGam(0:BasisOrder,0, &
          & ibasis,ibin), dt1*Beta/MxT)

        GamBasis(Order, ityp, dr, ibin, ibasis) = GamBasis(Order, ityp, dr, ibin, &
          & ibasis) + (Phase/factorM)*wbasis
        ReGamSqBasis(Order, ityp, dr, ibin, ibasis) = ReGamSqBasis(Order, ityp, dr, ibin, &
          & ibasis) + (real(Phase)/factorM)**2.d0*wbasis
        ImGamSqBasis(Order, ityp, dr, ibin, ibasis) = ImGamSqBasis(Order, ityp, dr, ibin, &
          & ibasis) + (dimag(Phase)/factorM)**2.d0*wbasis
      enddo
    endif

    Quan(Order) = Quan(Order) + real(Phase)/factorM
    Norm(Order) = Norm(Order) + 1.d0
    !print *, Order, Quan(Order), Norm(Order)

    !==============  Error renormalization =========================
    if(Order==1 .and. ityp==1 .and. dr==0) then
      if(dt1==0 .and. dt2==0) then
        Norm(MCOrder+1) = Norm(MCOrder+1)+1.d0
        Quan(MCOrder+1) = Quan(MCOrder+1) + real(Phase)/factorM
      endif
    endif
    !================================================================

    !===============  test variables =================================
    sumt = 0
    do ikey = 1, NWLn
      i = WLnKey2Value(ikey)
      sumt = sumt+ TypeLn(i)
    enddo

    sumd = 0
    do ikey = 1, NVertex
      i = VertexKey2Value(ikey)
      sumd = sumd+ IsDeltaVertex(i)
    enddo

    if(sumt==NWLn .and. sumd==NVertex) then
      Quan(MCOrder+Order+2) = Quan(MCOrder+Order+2) + 1.d0/CoefOfWeight(Order)
      Norm(MCOrder+Order+2) = Norm(MCOrder+Order+2) + 1.d0
    endif
    !=============================================================

  endif
END SUBROUTINE measure

!================   statistics =======================================

SUBROUTINE statistics
    implicit none
    integer :: i,iorder
    double precision :: x
    double precision,allocatable :: temp(:,:)
    StatNum=StatNum+1
    if(StatNum>=MaxStat) then
      allocate(temp(MaxStat,0:NObs-1))
      temp=ObsRecord
      MaxStat=MaxStat*2
      if(MaxStat>MxNblck) then
        call LogFile%QuickLog("Too many memory blocks, even bigger than "//trim(str(MxNblck)), 'e')
        stop
      endif
      deallocate(ObsRecord)
      allocate(ObsRecord(MaxStat,0:NObs-1))
      ObsRecord=0.0
      ObsRecord(1:MaxStat/2,:)=temp
      deallocate(temp)
    endif

    do i=0,NObs-1
      if(Norm(i)>1e-6) then
        x=Quan(i)/Norm(i)
        ObsRecord(StatNum,i)=x 
        call ERSTAT(ObsRecord(:,i),StatNum,amax,tmax,amin,tmin) 
        Error(i)=(amax-amin)/2.d0;
      else
        Error(i)=100.d0
      endif
    enddo
end SUBROUTINE

! error bar analysis from 3/4 of the file
subroutine ERSTAT(a,nre,amax,tmax,amin,tmin)
!   Analizing 3/4 print-out
   
  integer, intent(IN) :: nre
  integer :: i
  double precision :: a(nre), amax, tmax, amin, tmin, aa

  amax=-1.d200
  amin=1.d200
  DO i=nre/4+1, nre
     aa=a(i)
     if (aa > amax) then
        amax=aa
        tmax=i
     end if
     if (aa < amin) then
        amin=aa
        tmin=i
     end if
  END DO
  tmax=tmax/nre
  tmin=tmin/nre
end subroutine ERSTAT

!	call ERSTAT(record,prntout,amax,tmax,amin,tmin)
!      print 704, record(prntout),( amax-amin)/2.
! 704  format(6x,'record =',g12.5,4x,'+-',g10.3)  
!!=======================================================================
!!=======================================================================
!!=======================================================================


!====================================================================
