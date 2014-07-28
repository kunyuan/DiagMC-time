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
    wln = weight_wline(StatusLn(i),IsDeltaLn(i), diff_r(D, WRVertex(Gam1),WRVertex(Gam2)), &
      &  TVertex(3, Gam2)-TVertex(3,Gam1), TypeLn(i))
    WeightLn(i) = wln
    weight = weight *wln
  enddo


  do ikey = 1, NVertex
    i = VertexKey2Value(ikey)
    tau1 = TVertex(3, i)-TVertex(2, i)
    tau2 = TVertex(1, i)-TVertex(3, i)
    wgam = weight_vertex(StatusVertex(i), IsDeltaVertex(i), diff_r(D,GRVertex(i),WRVertex(i)), &
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

  if(HEAVY_DEBUG .and. CHECK_G) then
    if(Hash4G(newk)/=0 .or. Hash4G(oldk)/=GLn) then
      call LogFile%WriteStamp('e')
      call LogFile%WriteLine("Oops, update_Hash4G found a bug!")
      call LogFile%WriteLine("IsWormPresent:"+str(IsWormPresent)+", update number:"+str(imc))
      call LogFile%WriteLine("G Hash table for old k"+str(newk)+" is not 1!!"+str(Hash4G(newk)))
      call print_config
      stop
    endif
  endif

  Hash4G(newk)=Hash4G(oldk)
  Hash4G(oldk)=0
END SUBROUTINE update_Hash4G

SUBROUTINE add_Hash4G(newk, GLn)
  implicit none
  integer, intent(in) :: newk, GLn
  !MAKE SURE Hash4G(newk)==0 BEFORE CALL THIS SUBROUTINE
  if(HEAVY_DEBUG .and. CHECK_G .and. Hash4G(newk)/=0) then
    call LogFile%WriteStamp('e')
    call LogFile%WriteLine("Oops, add_Hash4G found a bug!")
    call LogFile%WriteLine("IsWormPresent:"+str(IsWormPresent)+", update number:"+str(imc))
    call LogFile%WriteLine("G Hash table for old k"+str(newk)+" is not 1!!"+str(Hash4G(newk)))
    call print_config
    stop
  endif
  Hash4G(newk)=GLn
END SUBROUTINE add_Hash4G

SUBROUTINE delete_Hash4G(oldk, GLn)
  implicit none
  integer, intent(in) :: oldk, GLn
  if(Hash4G(oldk)==GLn) then
    Hash4G(oldk)=0
  else
    if(HEAVY_DEBUG .and. CHECK_G) then
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
!reducibility check strategy
!for both worm and normal space, pick up two G lines and one W line, the following rule give the definition
!of reducibility:
!     |k_G1-kG2|==|k_W|
!unless: G1, G2 and W are connect to the same Gamma and this Gamma is not Ira or Masha
!You may refer check_conf.f90: check_irreducibility subroutine to find the implementation of the rule

!--------- check the irreducibility for G -----------------------
LOGICAL FUNCTION Is_k_valid_for_G(k)
  implicit none
  integer,intent(in) :: k
  if(CHECK_G .and. Hash4G(k)/=0) then
    Is_k_valid_for_G=.false.
  else
    Is_k_valid_for_G=.true.
  endif
END FUNCTION
!--------- check the irreducibility for W -----------------------

LOGICAL FUNCTION Is_k_valid_for_W(k)
  implicit none
  integer,intent(in) :: k
  integer :: ak
  ak=abs(k)
  if(CHECK_W .and. (ak==0 .or. Hash4W(ak)/=0)) then
    Is_k_valid_for_W=.false.
  else
    Is_k_valid_for_W=.true.
  endif
END FUNCTION
!-------------- check the ireeducibility of Gamma --------------
LOGICAL FUNCTION Is_reducible_G_Gam_one_side(kG, kNeighG)
  implicit none
  integer, intent(in) :: kG, kNeighG
  integer :: i, knG

  Is_reducible_G_Gam_one_side = .false.
  if(.not. CHECK_GAM) return

  do i = 1, NGLn
    knG = kLn(GLnKey2Value(i))
    if(knG/=kG .and. knG/=kNeighG) then
      if(Hash4W(abs(add_k(knG, -kG)))/=0) then
          Is_reducible_G_Gam_one_side = .true.
          return
      endif
    endif
  enddo

  return
END FUNCTION Is_reducible_G_Gam_one_side

LOGICAL FUNCTION Is_reducible_W_Gam(kW)
  implicit none
  integer, intent(in) :: kW
  integer :: i, ktemp, kG

  Is_reducible_W_Gam = .false.
  if(.not. CHECK_GAM) return

  do i = 1, NGLn
    kG = kLn(GLnKey2Value(i))
    ktemp = add_k(kG, kW)
    if(Hash4G(ktemp)/=0) then
      Is_reducible_W_Gam=.true.
      return
    endif

    ktemp = add_k(kG, -kW)
    if(Hash4G(ktemp)/=0) then
      Is_reducible_W_Gam=.true.
      return
    endif
  enddo
END FUNCTION Is_reducible_W_Gam

LOGICAL FUNCTION Is_reducible_W_Gam_one_side(kW, kG1, kG2)
  implicit none
  integer, intent(in) :: kW, kG1, kG2
  integer :: i, ktemp, kG

  Is_reducible_W_Gam_one_side = .false.
  if(.not. CHECK_GAM) return

  do i = 1, NGLn
    kG = kLn(GLnKey2Value(i))
    ktemp = add_k(kG, kW)
    if((kG/=kG1 .or. ktemp/=kG2) .and. (kG/=kG2 .or. ktemp/=kG1)) then
      if(Hash4G(ktemp)/=0) then
        Is_reducible_W_Gam_one_side=.true.
        return
      endif
    endif

    ktemp = add_k(kG, -kW)
    if((kG/=kG1 .or. ktemp/=kG2) .and. (kG/=kG2 .or. ktemp/=kG1)) then
      if(Hash4G(ktemp)/=0) then
        Is_reducible_W_Gam_one_side=.true.
        return
      endif
    endif
  enddo
END FUNCTION Is_reducible_W_Gam_one_side


LOGICAL FUNCTION Is_reducible_W_Gam_both_side(kW, kIA, kIC, kMB, kMD)
  implicit none
  integer,intent(in) :: kW, kIA, kMB, kIC, kMD
  integer :: kG
  integer :: ktemp

  Is_reducible_W_Gam_both_side = .false.
  if(.not. CHECK_GAM) return

  do i = 1, NGLn
    kG = kLn(GLnKey2Value(i))
    ktemp = add_k(kG, kW)
    !ktemp is G here
    if((kG/=kIA .or. ktemp/=kIC) .and. (kG/=kIC .or. ktemp/=kIA) &
      & .and. (kG/=kMB .or. ktemp/=kMD) .and. (kG/=kMD .or. ktemp/=kMB)) then
      if(Hash4G(ktemp)/=0) then
          Is_reducible_W_Gam_both_side=.true.
          return
      endif
    endif

    ktemp = add_k(kG, -kW)
    !ktemp is G here
    if((kG/=kIA .or. ktemp/=kIC) .and. (kG/=kIC .or. ktemp/=kIA) &
      & .and. (kG/=kMB .or. ktemp/=kMD) .and. (kG/=kMD .or. ktemp/=kMB)) then
      if(Hash4G(ktemp)/=0) then
          Is_reducible_W_Gam_both_side=.true.
          return
      endif
    endif
  enddo 
end FUNCTION

!------------- check the irreducibility for add interaction operation --------------------

LOGICAL FUNCTION Is_reducible_add_interaction(kW, kIA, kIC, kMB, kMD)
  implicit none
  integer,intent(in) :: kW, kIA, kMB, kIC, kMD
  integer :: kG
  integer :: ktemp

  Is_reducible_add_interaction = .false.
  if(.not. CHECK_GAM) return

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

    if(kG/=kIA .and. kG/=kIC) then
    !check GIA
      if(Hash4W(abs(add_k(kG, -kIA)))/=0) then
        Is_reducible_add_interaction=.true.
        return
      endif
    endif

    if(kG/=kMB .and. kG/=kMD) then
    !check GMB
      if(Hash4W(abs(add_k(kG, -kMB)))/=0) then
        Is_reducible_add_interaction=.true.
        return
      endif
    endif
  enddo
END FUNCTION Is_reducible_add_interaction

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
SUBROUTINE insert_gamma(newgamma, isdelta, gr, wr, t1, t2, t3, dir, typ, stat, weigh)
  implicit none
  integer, intent(out) :: newgamma
  integer, intent(in) :: gr, wr, isdelta, dir, typ, stat
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
  GRVertex(newgamma) = gr
  WRVertex(newgamma) = wr
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

!============   initialize diagrams =================================
SUBROUTINE first_order_diagram
  implicit none
  integer :: i, dr0
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
  GRVertex(1:4) = 0
  WRVertex(1:4) = 0

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

  dr0 = 0
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
  integer :: i, dr0
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
  GRVertex(1:4) = 0
  WRVertex(1:4) = 0

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

  dr0 = 0
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
  integer :: i, dr0
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
  GRVertex(1:6) = 0
  WRVertex(1:6) = 0

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

  dr0 = 0
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
