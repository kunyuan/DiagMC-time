
!============ Initialization of the self consistent loop ==========
SUBROUTINE initialize_self_consistent
  implicit none
  integer :: it, ityp, ix, iy
  integer :: it1, it2

  !------- Initialization of the self-consistent loop ------------
  call initialize_G
  call initialize_W
  call initialize_Gam

  call initialize_G0F
  call initialize_W0PF
  call initialize_Gam0PF
  
END SUBROUTINE initialize_self_consistent


!--------------------- Initialization of G ----------------------
SUBROUTINE initialize_G
  implicit none
  integer :: typ, t

  G(:,:) = (0.d0, 0.d0)

  do t = 0, MxT-1
    do typ = 1, NTypeG
      G(typ, t) = weight_G0(typ, t)
    enddo
  enddo
END SUBROUTINE initialize_G

!!------------------- Initialization of W ----------------------
SUBROUTINE initialize_W
  implicit none
  integer :: typ, t, dx, dy

  W(:,:,:,:) = (0.d0, 0.d0)
END SUBROUTINE initialize_W

!!--------------- Initialization of Gamma -----------------------
SUBROUTINE initialize_Gam
  implicit none
  integer :: ityp, it1, it2

  Gam(:,:,:,:,:) = (0.d0, 0.d0)
END SUBROUTINE initialize_Gam
 
!!------- Initialization of G0 in frequency ----------
SUBROUTINE initialize_G0F
  implicit none
  integer :: typ, t

  G0F = (0.d0, 0.d0)
  do t = 0, MxT-1
    G0F(t) = weight_G0(1, t)
  enddo
  call transfer_G0(1)
END SUBROUTINE initialize_G0F


!!------- Initialization of W0 in monmentum and frequency ----------
SUBROUTINE initialize_W0PF
  implicit none
  double precision :: ratio

  W0PF(:,:,:) = (0.d0, 0.d0)
  ratio = 0.25d0*Jcp*real(MxT)/Beta

  W0PF(1,    0, 0) = dcmplx(ratio, 0.d0)
  W0PF(Lx-1, 0, 0) = dcmplx(ratio, 0.d0)
  W0PF(0,    1, 0) = dcmplx(ratio, 0.d0)
  W0PF(0, Ly-1, 0) = dcmplx(ratio, 0.d0)

  call transfer_W0(1)
END SUBROUTINE initialize_W0PF

!!------- Initialization of Gam0 in monmentum and frequency ----------
SUBROUTINE initialize_Gam0PF
  implicit none
  double precision :: ratio

  ratio = (real(MxT)/Beta)**2.d0
  Gam0PF(:,:,:,:) = (0.d0, 0.d0)
  Gam0PF(0,0,0,0) = dcmplx(ratio, 0.d0)

  call transfer_Gam0(1)
END SUBROUTINE initialize_Gam0PF
 
!!============== WEIGHT CALCULATING ==================================

SUBROUTINE calculate_Polar
  implicit none
  integer :: px, py, omega
  integer :: omegaGin, omegaGout
  complex(kind=8) :: Gin, Gout, Gam1
  double precision :: ratio

  Polar(:,:,:) = (0.d0, 0.d0)

  ratio = 2.d0/real(MxT)*(Beta/real(MxT))**4.d0
  do omega = 0, MxT-1
    do omegaGin = 0, MxT-1
      do px = 0, Lx-1
        do py = 0, Ly-1
          omegaGout = omegaGin - omega

          Gin = weight_G(1, omegaGin)
          if(omegaGout>=0) then
            Gout = weight_G(1, omegaGout)
            Gam1 = weight_Gam(5, px, py, omegaGin, omegaGout)
          else
            Gout = weight_G(1, omegaGout+MxT)
            Gam1 = weight_Gam(5, px, py, omegaGin, omegaGout+MxT)
          endif
          
          Polar(px, py, omega) = Polar(px, py, omega)+d_times_cd(ratio, cdexp((0.d0, -1.d0) &
            & *2.d0*omegaGout*Pi/MxT)*Gin*Gout*Gam1)
        enddo
      enddo
    enddo
  enddo
  return
END SUBROUTINE calculate_Polar


SUBROUTINE calculate_Sigma
  implicit none
  integer :: px, py, omega
  integer :: omegaG, omegaW
  complex(kind=8) :: G1, W1, Gam1
  double precision :: ratio

  Sigma(:) = (0.d0, 0.d0)

  ratio = -3.d0/(real(Lx)*real(Ly)*real(MxT))*(Beta/real(MxT))**4.d0
  do omega = 0, MxT-1
    do omegaG = 0, MxT-1
      do py = 0, Ly-1
        do px = 0, Lx-1
          omegaW = omega-omegaG

          G1 = weight_G(1, omegaG)
          if(omegaW>=0) then
            W1 = weight_W(1, px, py, omegaW)
          else
            W1 = weight_W(1, px, py, omegaW+MxT)
          endif
          Gam1 = weight_Gam(5, px, py, omegaG, omega)
          
          Sigma(omega) = Sigma(omega)+d_times_cd(ratio, G1*W1*Gam1)
        enddo
      enddo
    enddo
  enddo
  return
END SUBROUTINE calculate_Sigma


!!--------- calculate weight for W matrix ---------
SUBROUTINE calculate_W
  implicit none
  integer :: omega, px, py

  !-------- calculate W = W0/(1-W0*G^2*Gamma) ----------------------------
  W(:,:,:,:) = (0.d0, 0.d0)
  W(1,:,:,:) = W0PF(:,:,:)/((1.d0, 0.d0)-W0PF(:,:,:)*Polar(:,:,:))

  do  omega = 0, MxT-1
    do py = 0, Ly-1
      do px = 0, Lx-1
        W(3,px,py,omega) = d_times_cd(-1.d0,W(1,px,py,omega))
        W(5,px,py,omega) = d_times_cd( 2.d0,W(1,px,py,omega))
      enddo
    enddo
  enddo

  W(2,:,:,:) = W(1,:,:,:)
  W(4,:,:,:) = W(3,:,:,:)
  W(6,:,:,:) = W(5,:,:,:)

  !!-------------- update the matrix and tail ------------
END SUBROUTINE calculate_W

!!--------- calculate weight for G matrix ---------
SUBROUTINE calculate_G
  implicit none
  complex(kind=8) :: G0

  G(:,:) = (0.d0, 0.d0)

  !--------- G = G0/(1-G0*Sigma)----------------------
  G(1, :) =  G0F(:)/((1.d0,0.d0)-G0F(:)*Sigma(:))
  G(2, :) =  G(1, :)

  !!-------------- update the matrix and tail ------------
END SUBROUTINE calculate_G




!!------------- calculate the susceptibility ---------------------------
SUBROUTINE calculate_Chi
  implicit none
  integer :: px, py, omega
  double precision :: ratio

  !!-------- calculate Chi = Pi/(1 - W0 * Pi) ------------------
  !!--------- already sum over the spins -----------------------
  Chi(:,:,:) = 0.d0
  ratio = -1.d0*(real(MxT)/Beta)**2.d0*0.75d0

  Chi(:,:,:) = Polar(:,:,:)/((1.d0,0.d0) -W0PF(:,:,:)*Polar(:,:,:))

  do omega = 0, MxT-1
    do py = 0, Ly-1
      do px = 0, Lx-1
        Chi(px, py, omega) = d_times_cd(ratio, Chi(px, py, omega))
      enddo
    enddo
  enddo

  return
END SUBROUTINE calculate_Chi

!!====================================================================

SUBROUTINE plus_minus_W0(Backforth)
  implicit none
  integer, intent(in) :: Backforth
  integer :: ityp, px, py, omega

  if(Backforth/=-1) then
    W(1,:,:,:) = W(1,:,:,:) + W0PF(:,:,:)
  else 
    W(1,:,:,:) = W(1,:,:,:) - W0PF(:,:,:)
  endif

  do  omega = 0, MxT-1
    do py = 0, Ly-1
      do px = 0, Lx-1
        W(3,px,py,omega) = d_times_cd(-1.d0,W(1,px,py,omega))
        W(5,px,py,omega) = d_times_cd( 2.d0,W(1,px,py,omega))
      enddo
    enddo
  enddo

  W(2,:,:,:) = W(1,:,:,:)
  W(4,:,:,:) = W(3,:,:,:)
  W(6,:,:,:) = W(5,:,:,:)

  return
END SUBROUTINE

SUBROUTINE plus_minus_Gam0(Backforth)
  implicit none
  integer, intent(in) :: Backforth
  integer :: ityp, px, py, omega1, omega2

  if(Backforth/=-1) then
    Gam(1,:,:,:,:) = Gam(1,:,:,:,:) + Gam0PF(:,:,:,:)
    Gam(2,:,:,:,:) = Gam(2,:,:,:,:) + Gam0PF(:,:,:,:)
    Gam(5,:,:,:,:) = Gam(5,:,:,:,:) + Gam0PF(:,:,:,:)
    Gam(6,:,:,:,:) = Gam(6,:,:,:,:) + Gam0PF(:,:,:,:)
  else 
    Gam(1,:,:,:,:) = Gam(1,:,:,:,:) - Gam0PF(:,:,:,:)
    Gam(2,:,:,:,:) = Gam(2,:,:,:,:) - Gam0PF(:,:,:,:)
    Gam(5,:,:,:,:) = Gam(5,:,:,:,:) - Gam0PF(:,:,:,:)
    Gam(6,:,:,:,:) = Gam(6,:,:,:,:) - Gam0PF(:,:,:,:)
  endif

  return
END SUBROUTINE


!!======================== WEIGHT EXTRACTING =========================

!--------- weight for bare propagator ----------------
Complex*16 FUNCTION weight_G0(typ, t)
  implicit none
  integer, intent(in)  :: typ, t  
  double precision     :: tau
  complex(kind=8)      :: muc  

  muc = dcmplx(0.d0, Mu(1)*pi/(2.d0*Beta))
  tau = real(t)*Beta/MxT
  if(tau>=0) then
    weight_G0 = cdexp(muc*tau)/(1.d0, 1.d0) 
  else
    weight_G0 = -cdexp(muc*(tau+Beta))/(1.d0, 1.d0) 
  endif
  return
END FUNCTION weight_G0


!!--------- calculate weight for bare interaction ----
Complex*16 FUNCTION weight_W0(typ, dx, dy)
  implicit none
  integer, intent(in) :: dx, dy, typ
  integer :: dx1, dy1
  double precision :: ratio

  ratio = Jcp

  dx1 = dx;       dy1 = dy
  if(dx1>=0  .and. dx1<Lx .and. dy1>=0 .and. dy1<Ly) then
    if(dx1>dLx)     dx1 = Lx-dx1
    if(dy1>dLy)     dy1 = Ly-dy1

    weight_W0 = (0.d0, 0.d0)

    if((dx1==1.and.dy1==0).or.(dx1==0.and.dy1==1)) then
      if(typ ==1 .or. typ == 2) then
        weight_W0 = dcmplx(0.25d0*ratio, 0.d0)
      else if(typ == 3 .or. typ == 4) then
        weight_W0 = dcmplx(-0.25d0*ratio, 0.d0)
      else if(typ == 5 .or. typ == 6) then
        weight_W0 = dcmplx(0.5d0*ratio, 0.d0)
      endif
    endif
  else
    write(*, *) dx1, dy1, "dx, dy bigger than system size!"
    stop
  endif
END FUNCTION weight_W0


!!--------- calculate weight for bare Gamma ---------
COMPLEX*16 FUNCTION weight_Gam0(typ, dx, dy)
  implicit none
  integer, intent(in)  :: dx, dy, typ
  double precision :: ratio

  if(dx>=0 .and. dx<Lx .and. dy>=0 .and. dy<Ly) then
    if(dx==0.and.dy==0) then
      if(typ==1 .or. typ==2 .or. typ==5 .or. typ==6) then
        weight_Gam0 = (1.d0, 0.d0)
      else
        weight_Gam0 = (0.d0, 0.d0)
      endif
    else
      weight_Gam0 = (0.d0, 0.d0)
    endif
  else
    write(*, *) dx, dy, "dx, dy bigger than system size!"
    stop
  endif
END FUNCTION weight_Gam0

!!--------- extract weight for G ---------
COMPLEX*16 FUNCTION weight_G(typ1, t1)
  implicit none
  integer, intent(in)  :: t1, typ1

  if(t1>=0) then
    weight_G = G(typ1, t1)
  else
    weight_G = -G(typ1, t1+MxT)
  endif
END FUNCTION weight_G

!!--------- extract weight for W ---------
COMPLEX*16 FUNCTION weight_W(typ1, dx1, dy1, t1)
  implicit none
  integer, intent(in)  :: dx1, dy1, t1, typ1

  if(t1>=0) then
    weight_W = W(typ1, dx1, dy1, t1)
  else
    weight_W = W(typ1, dx1, dy1, t1+MxT)
  endif
END FUNCTION weight_W

!!--------- extract weight for Gamma ---------
COMPLEX*16 FUNCTION weight_Gam(typ1, dx1, dy1, t1, t2)
  implicit none
  integer, intent(in)  :: dx1, dy1, t1, t2, typ1
  double precision :: GaR
  integer :: ib, jb, dx, dy
  
  if(t1>=0 .and. t2>=0) then
    weight_Gam = Gam(typ1, dx1, dy1, t1, t2)
  else if(t1<0 .and. t2>=0) then
    weight_Gam = -Gam(typ1, dx1, dy1, t1+MxT, t2)
  else if(t1>=0 .and. t2<0) then
    weight_Gam = -Gam(typ1, dx1, dy1, t1, t2+MxT)
  else
    weight_Gam = -Gam(typ1, dx1, dy1, t1+MxT, t2+MxT)
  endif
END FUNCTION weight_Gam

!!====================================================================



SUBROUTINE Gam_mc2matrix_mc
  implicit none
  integer :: iorder, dx, dy, ityp, iloop, t1, t2, typ
  complex*16 :: norm

  call initialize_Gam

  norm = GamNormWeight/GamNorm

  do ityp = 1, NTypeGam/2
    do iloop = 1, 2
      do iorder = 1, MCOrder
      typ = 2*(ityp-1)+1
      Gam(typ, :, :, :, :)=Gam(typ, :, :, :, :)+GamMC(iorder,iloop,ityp,:,:,:,:)*norm
      enddo
    enddo
  enddo

  Gam(2,:,:,:,:) = Gam(1,:,:,:,:)
  Gam(4,:,:,:,:) = Gam(3,:,:,:,:)
  Gam(6,:,:,:,:) = Gam(5,:,:,:,:)

END SUBROUTINE Gam_mc2matrix_mc
 



!SUBROUTINE update_Gamma_matrix(norder)
  !implicit none
  !integer, intent(in) :: norder
  !integer :: iorder, dx, dy, ityp, omega1, omega2, ib, jb
  !double precision :: gam1, gam2, perr, norm, nmc
  !double precision :: Ga0R, GaR1
  !double precision :: tempGamMC(-MxOmegaGamG2:MxOmegaGamG2, -MxOmegaGamG2:MxOmegaGamG2)

  !tempGamMC(:, :) = 0.d0


  !if(norder==0) then
    !call initialize_Gamma
    !return
  !endif

  !norm = GamNormWeight*ime/GamNorm

  !do ityp = 1, ntypGa
    !do dx = 0, dLx
      !do dy = 0, dLy


      !enddo
    !enddo
  !enddo

!END SUBROUTINE update_Gamma_matrix
!====================================================================
