!!====================== initialization =================================
!!=======================================================================

SUBROUTINE initialize_markov
    implicit none

    Pupdate(:)  = 1.d0
    Pupdate(2)  = 5.d0

    Pupdate(3)  = 0.d0
    Pupdate(4)  = 0.d0

    !--------------- initialize variables ---------------
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
  do while(imeasure < nw)
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
