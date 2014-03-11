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
    if
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
end SUBROUTINE

