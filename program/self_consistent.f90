
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
  
END SUBROUTINE initialize_self_consistent


!--------------------- Initialization of G ----------------------
SUBROUTINE initialize_G
  implicit none
  integer :: typ, t

  G = (0.d0, 0.d0)

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

  W = (0.d0, 0.d0)
  do t = 0, MxT-1
    do dy = 0, Ly-1
      do dx = 0, Lx-1
        do typ = 1, NtypeW
          W(typ, dx, dy, t) = weight_W0(typ, dx, dy, t)
        enddo
      enddo
    enddo
  enddo
END SUBROUTINE initialize_W

!!--------------- Initialization of Gamma -----------------------
SUBROUTINE initialize_Gam
  implicit none
  integer :: ityp, it1, it2

  Gam = (0.d0, 0.d0)
  do it2 = 0, MxT-1
    do it1 = 0, MxT-1
      do ityp = 1, NTypeGam
        Gam(ityp,0,0,it1,it2) = weight_Gam0(ityp, 0, 0, it1, it2)
      enddo
    enddo
  enddo
  return
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
  integer :: typ, t, dx, dy

  W0PF = (0.d0, 0.d0)
  do t = 0, MxT-1
    do dy = 0, Ly-1
      do dx = 0, Lx-1
        W0PF(dx, dy, t) = weight_W0(1, dx, dy, t)
      enddo
    enddo
  enddo
  call transfer_W0(1)
END SUBROUTINE initialize_W0PF

!!============== WEIGHT CALCULATING ==================================

SUBROUTINE calculate_Polar
  implicit none
  integer :: px, py, omega
  integer :: omegaGin, omegaGout
  complex(kind=8) :: Gin, Gout, Gam1
  double precision :: ratio

  ratio = 2.d0/Beta
  Polar(:,:,:) = (0.d0, 0.d0)
  do omega = 0, MxT-1
    do omegaGin = omega, MxT-1
      do px = 0, Lx-1
        do py = 0, Ly-1
          omegaGout = omegaGin - omega

          Gin = weight_G(1, omegaGin)
          Gout = weight_G(1, omegaGout)
          Gam1 = weight_Gam(5, px, py, omegaGin, omegaGout)
          
          Polar(px, py, omega) = Polar(px, py, omega)+d_times_cd(ratio, Gin*Gout*Gam1)
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

  ratio = -3.d0/(real(Lx)*real(Ly)*Beta)
  Sigma(:) = (0.d0, 0.d0)

  do omega = 0, MxT-1
    do omegaG = 0, omega
      do py = 0, Ly-1
        do px = 0, Lx-1
          omegaW = omega-omegaG

          G1 = weight_G0(1, omegaG)
          W1 = weight_W0(1, px, py, omegaW)
          Gam1 = weight_Gam0(5, px, py, omegaG, omega)
          
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
!SUBROUTINE calculate_Chi
  !implicit none
  !integer :: px, py, omega, ityp
  !double precision :: Pi1, Pi2, W1, W2
  !double precision :: PiW1, PiW2, temp
  !double precision :: coef1, coef2
  !integer :: dx, dy
  !double precision :: W0R, ChR, W0ChR
  !double precision :: trChR



  !!-------- calculate Chi = Pi/(1 - W0 * Pi) ------------------
  !ChiR(:,:,:,:) = 0.d0

  !call transfer_Pi(1)
  !call transfer_Chi(1)

  !do px = 0, dLx
    !do py = 0, dLy
      !do omega = -MxOmegaChi, MxOmegaChi

        !W1 = W0InMoment(1, px, py)
        !W2 = W0InMoment(3, px, py)
        !Pi1 = PiR(1, px, py, omega)
        !Pi2 = PiR(3, px, py, omega)

        !PiW1 = Pi1*W1 + Pi2*W2
        !PiW2 = Pi1*W2 + Pi2*W1
        !temp = (1.d0-PiW1)**2.d0-PiW2**2.d0

        !coef1 = (1.d0-PiW1)/temp
        !coef2 = PiW2/temp

        !ChiR(1, px, py, omega) = coef1*Pi1 + coef2*Pi2
        !ChiR(3, px, py, omega) = coef1*Pi2 + coef2*Pi1

        !do ityp = 2, 4, 2
          !ChiR(ityp, px, py, omega) = ChiR(ityp-1, px, py, omega)
        !enddo

      !enddo
    !enddo
  !enddo

  !call transfer_Pi(-1)
  !call transfer_Chi(-1)

  !trChiR(:,:,:) = 0.d0
  !do dx = 0, dLx
    !do dy = 0, dLy
      !do omega = -MxOmegaChi, MxOmegaChi

        !trChR = 0.d0
        !do ityp = 1, 4
          !if(ityp<=2) then
            !W0R = 0.25
          !else
            !W0R = -0.25
          !endif
          !ChR = ChiR(ityp, dx, dy, omega)
          !W0ChR = ChR *W0R

          !trChR = trChR + W0ChR
        !enddo
        !trChiR(dx, dy, omega) = trChR
      !enddo
    !enddo
  !enddo
!END SUBROUTINE calculate_Chi

!!====================================================================


!!======================== WEIGHT EXTRACTING =========================

!--------- weight for bare propagator ----------------
Complex*16 FUNCTION weight_G0(typ, t)
  implicit none
  integer, intent(in)  :: typ, t  
  double precision     :: tau
  complex(kind=8)      :: muc  

  muc = dcmplx(0.d0, Mu(1)*pi/(2.d0*Beta))
  tau = (t+0.5d0)*Beta/MxT
  if(tau>=0) then
    weight_G0 = cdexp(muc*tau)/(1.d0, 1.d0) 
  else if(tau>=-MxT) then
    weight_G0 = -cdexp(muc*(tau+Beta))/(1.d0, 1.d0) 
  endif
  return
END FUNCTION weight_G0


!!--------- calculate weight for bare interaction ----
Complex*16 FUNCTION weight_W0(typ, dx, dy, t)
  implicit none
  integer, intent(in) :: dx, dy, typ, t
  integer :: dx1, dy1

  dx1 = dx;       dy1 = dy
  if(dx1>=0  .and. dx1<Lx .and. dy1>=0 .and. dy1<Ly) then
    if(dx1>dLx)     dx1 = Lx-dx1
    if(dy1>dLy)     dy1 = Ly-dy1

    weight_W0 = (0.d0, 0.d0)

    if(t==0) then
      !if((dx1==1.and.dy1==0).or.(dx1==0.and.dy1==1)) then
      if(dx1==0.and.dy1==0) then
        if(typ ==1 .or. typ == 2) then
          weight_W0 = dcmplx(0.25d0*Jcp, 0.d0)
        else if(typ == 3 .or. typ == 4) then
          weight_W0 = dcmplx(-0.25d0*Jcp, 0.d0)
        else if(typ == 5 .or. typ == 6) then
          weight_W0 = dcmplx(0.5d0*Jcp, 0.d0)
        endif
      endif
    endif
  else
    write(*, *) dx1, dy1, "dx, dy bigger than system size!"
    stop
  endif

END FUNCTION weight_W0

!!--------- calculate weight for bare Gamma ---------
COMPLEX*16 FUNCTION weight_Gam0(typ, dx, dy, t1, t2)
  implicit none
  integer, intent(in)  :: dx, dy, t1, t2, typ

  if(dx>=0 .and. dx<Lx .and. dy>=0 .and. dy<Ly) then
    if(t1==0 .and. t2==0 .and. dx==0.and.dy==0) then
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
  double precision:: GGI
  integer :: ib

  weight_G = weight_G0(typ1, t1)
  !if(omega1>=-MxOmegaG1 .and. omega1<=MxOmegaG1) then
    !GGI = GI(typ1, omega1)
  !else if(omega1<-MxOmegaG1) then
    !GGI = 0.d0
    !do ib = 1, nbasis
      !GGI = GGI +GITailN(typ1, ib)*weight_basis(GCoefN(ib,:),omega1)
    !enddo
  !else if(omega1>MxOmegaG1) then
    !GGI = 0.d0
    !do ib = 1, nbasis
      !GGI = GGI +GITailP(typ1, ib)*weight_basis(GCoefP(ib,:),omega1)
    !enddo
  !endif
  !weight_G = GGI
END FUNCTION weight_G

!!--------- extract weight for W ---------
COMPLEX*16 FUNCTION weight_W(typ1, dx1, dy1, t1)
  implicit none
  integer, intent(in)  :: dx1, dy1, t1, typ1
  double precision :: WWR
  integer :: dx, dy, ib

  weight_W = weight_W0(typ1, dx1, dy1, t1)

  !dx = dx1;      dy = dy1
  !if(dx>=0 .and. dx<Lx .and. dy>=0 .and. dy<Ly) then
    !if(dx>dLx)     dx = Lx-dx
    !if(dy>dLy)     dy = Ly-dy

    !if(abs(omega1)<=MxOmegaW1) then
      !WWR = WR(typ1, dx, dy, abs(omega1))
    !else if(abs(omega1)>MxOmegaW1) then
      !WWR = WRTailC(typ1, dx, dy)
      !do ib = 1, nbasis
        !WWR = WWR +WRTailP(typ1,dx,dy,ib)*weight_basis(WCoefP(ib,:),abs(omega1))
      !enddo
      !!write(*, *) WWR
    !endif
  !else
    !write(*, *) dx, dy, "dx, dy bigger than system size!"
    !stop
  !endif
  !weight_W = WWR
END FUNCTION weight_W

!!--------- extract weight for Gamma ---------
COMPLEX*16 FUNCTION weight_Gam(typ1, dx1, dy1, t1, t2)
  implicit none
  integer, intent(in)  :: dx1, dy1, t1, t2, typ1
  double precision :: GaR
  integer :: ib, jb, dx, dy
  
  weight_Gam = weight_Gam0(typ1, dx1, dy1, t1, t2)

  
  !dx = dx1;      dy = dy1
  !if(dx>=0 .and. dx<Lx .and. dy>=0 .and. dy<Ly) then
    !if(dx>dLx)     dx = Lx-dx
    !if(dy>dLy)     dy = Ly-dy

    !if(abs(omega1)<=MxOmegaGamG1 .and. abs(omega2)<=MxOmegaGamG1) then
      !GaR = GamR(typ1, dx, dy, omega1, omega2)

    !else if(omega1>MxOmegaGamG1 .and. abs(omega2)<=MxOmegaGamG1) then
      !GaR = GamRTailC(typ1, dx, dy)
      !do ib = 1, nbasis
        !GaR = GaR + GamRTailPM(typ1, dx, dy, ib, omega2)*weight_basis(&
          !& GamGCoefP(ib,:), omega1)
      !enddo
    !else if(omega1<-MxOmegaGamG1 .and. abs(omega2)<=MxOmegaGamG1) then
      !GaR = GamRTailC(typ1, dx, dy)
      !do ib = 1, nbasis
        !GaR = GaR + GamRTailNM(typ1, dx, dy, ib, omega2)*weight_basis(&
          !& GamGCoefN(ib,:), omega1)
      !enddo
    !else if(omega2>MxOmegaGamG1 .and. abs(omega1)<=MxOmegaGamG1) then
      !GaR = GamRTailC(typ1, dx, dy)
      !do ib = 1, nbasis
        !GaR = GaR + GamRTailMP(typ1, dx, dy, omega1, ib)*weight_basis(&
          !& GamGCoefP(ib,:), omega2)
      !enddo
    !else if(omega2<-MxOmegaGamG1 .and. abs(omega1)<=MxOmegaGamG1) then
      !GaR = GamRTailC(typ1, dx, dy)
      !do ib = 1, nbasis
        !GaR = GaR + GamRTailMN(typ1, dx, dy, omega1, ib)*weight_basis(&
          !& GamGCoefN(ib,:), omega2)
      !enddo

    !else if(omega1>MxOmegaGamG1 .and. omega2>MxOmegaGamG1) then
      !if(omega1==omega2) then
        !GaR = GamRTailC(typ1, dx, dy)
        !do ib = 1, nbasis
          !GaR = GaR + GamRTailDiagP(typ1, dx, dy, ib)*weight_basis(&
            !& GamGCoefP(ib,:), omega1)
        !enddo
      !else if(omega1<omega2) then
        !GaR = GamRTailC(typ1, dx, dy)
        !do ib = 1, nbasisGamma
          !GaR = GaR + GamRTailPPL(typ1, dx, dy, ib)*weight_basis_Gamma(&
            !& GamCoefPPL(ib,:), omega1, omega2)
        !enddo
      !else if(omega1>omega2) then
        !GaR = GamRTailC(typ1, dx, dy)
        !do ib = 1, nbasisGamma
          !GaR = GaR + GamRTailPPR(typ1, dx, dy, ib)*weight_basis_Gamma(&
            !& GamCoefPPR(ib,:), omega1, omega2)
        !enddo
      !endif

    !else if(omega1<-MxOmegaGamG1 .and. omega2>MxOmegaGamG1) then
      !GaR = GamRTailC(typ1, dx, dy)
      !do ib = 1, nbasisGamma
        !GaR = GaR + GamRTailNP(typ1, dx, dy, ib)*weight_basis_Gamma(&
          !& GamCoefNP(ib,:), omega1, omega2)
      !enddo
    !else if(omega1>MxOmegaGamG1 .and. omega2<-MxOmegaGamG1) then
      !GaR = GamRTailC(typ1, dx, dy)
      !do ib = 1, nbasisGamma
        !GaR = GaR + GamRTailPN(typ1, dx, dy, ib)*weight_basis_Gamma(&
          !& GamCoefPN(ib,:), omega1, omega2)
      !enddo
    !else if(omega1<-MxOmegaGamG1 .and. omega2<-MxOmegaGamG1) then
      !if(omega1==omega2) then
        !GaR = GamRTailC(typ1, dx, dy)
        !do ib = 1, nbasis
          !GaR = GaR + GamRTailDiagN(typ1, dx, dy, ib)*weight_basis(&
            !& GamGCoefN(ib,:), omega1)
        !enddo
      !else if(omega1<omega2) then
        !GaR = GamRTailC(typ1, dx, dy)
        !do ib = 1, nbasisGamma
          !GaR = GaR + GamRTailNNL(typ1, dx, dy, ib)*weight_basis_Gamma(&
            !& GamCoefNNL(ib,:), omega1, omega2)
        !enddo
      !else if(omega1>omega2) then
        !GaR = GamRTailC(typ1, dx, dy)
        !do ib = 1, nbasisGamma
          !GaR = GaR + GamRTailNNR(typ1, dx, dy, ib)*weight_basis_Gamma(&
            !& GamCoefNNR(ib,:), omega1, omega2)
        !enddo
      !endif
    !endif
  !else
    !write(*, *) dx, dy, "dx, dy bigger than system size!"
    !stop
  !endif
  !weight_Gamma = GaR
END FUNCTION weight_Gam

!!====================================================================



!SUBROUTINE Gamma_mc2matrix_mc
  !implicit none
  !integer :: iorder, dx, dy, ityp, omega1, omega2, ib, jb
  !double precision :: gam1, gam2, perr, norm, nmc
  !double precision :: Ga0R, GaR1
  !double precision :: tempGamMC(-MxOmegaGamG2:MxOmegaGamG2, -MxOmegaGamG2:MxOmegaGamG2)

  !call initialize_Gamma
  !tempGamMC(:, :) = 0.d0

  !norm = GamNormWeight*ime/GamNorm

  !do ityp = 1, ntypGa
    !do dx = 0, dLx
      !do dy = 0, dLy

        !do omega1 = -MxOmegaDiag, MxOmegaDiag
          !do omega2 = -MxOmegaDiag, MxOmegaDiag
            !tempGamMC(omega1,omega2) = 0.d0
            !do iorder = 1, MCOrder
              !gam1 = SUM(GamMC(iorder, :,(ityp+1)/2, dx, dy, omega1,omega2))/ime
              !gam2 = SUM(GamSqMC(iorder, :,(ityp+1)/2, dx, dy, omega1, omega2))/ime
              !if(abs(gam1)>1.d-30) then
                !perr = sqrt(gam2-gam1**2.d0)/(sqrt(ime-1)*abs(gam1))
                !if(perr<=MxError) then
                  !tempGamMC(omega1, omega2) = tempGamMC(omega1, omega2) + gam1*norm
                !endif
              !endif
            !enddo
          !enddo
        !enddo

        !do omega1 = -MxOmegaGamG2, -MxOmegaGamG1-1
          !do omega2 = -MxOmegaGamG2, -MxOmegaGamG1-1
            !if(omega1==omega2) then
              !do ib = 1, nbasis
                !GamRTailDiagN(ityp,dx,dy,ib) = GamRTailDiagN(ityp,dx,dy,ib) &
                  !& + weight_basis(GamGCoefN(ib, :), omega1)* &
                  !& tempGamMC(omega1,omega2)
              !enddo
            !else if(omega1<omega2) then
              !do ib = 1, nbasisGamma
                !GamRTailNNL(ityp,dx,dy,ib) = GamRTailNNL(ityp,dx,dy,ib) &
                  !& + weight_basis_Gamma(GamCoefNNL(ib, :), omega1, omega2)* &
                  !& ReweightBasis(omega1, omega2) *tempGamMC(omega1,omega2)
              !enddo
            !else if(omega1>omega2) then
              !do ib = 1, nbasisGamma
                !GamRTailNNR(ityp,dx,dy,ib) = GamRTailNNR(ityp,dx,dy,ib) &
                  !& + weight_basis_Gamma(GamCoefNNR(ib, :), omega1, omega2)* &
                  !& ReweightBasis(omega1, omega2) *tempGamMC(omega1,omega2)
              !enddo
            !endif
          !enddo

          !do omega2 = -MxOmegaGamG1, MxOmegaGamG1
            !do ib = 1, nbasis
              !GamRTailNM(ityp,dx,dy,ib,omega2) = GamRTailNM(ityp,dx,dy,ib,omega2) &
                !& + weight_basis(GamGCoefN(ib, :), omega1)* &
                !& tempGamMC(omega1,omega2)
            !enddo
          !enddo

          !do omega2 = MxOmegaGamG1+1, MxOmegaGamG2
            !do ib = 1, nbasisGamma
              !GamRTailNP(ityp,dx,dy,ib) = GamRTailNP(ityp,dx,dy,ib) &
                !& + weight_basis_Gamma(GamCoefNP(ib, :), omega1, omega2)* &
                !& ReweightBasis(omega1, omega2) *tempGamMC(omega1,omega2)
            !enddo
          !enddo
        !enddo

        !do omega1 = -MxOmegaGamG1, MxOmegaGamG1
          !do omega2 = -MxOmegaGamG2, -MxOmegaGamG1-1
            !do ib = 1, nbasis
              !GamRTailMN(ityp,dx,dy,omega1,ib) = GamRTailMN(ityp,dx,dy,omega1,ib) &
                !& + weight_basis(GamGCoefN(ib,:), omega2)* &
                !& tempGamMC(omega1,omega2)
            !enddo
          !enddo

          !do omega2 = -MxOmegaGamG1, MxOmegaGamG1
            !GamR(ityp,dx, dy, omega1, omega2) = GamR(ityp,dx,dy,omega1,omega2)+tempGamMC(omega1, omega2)
          !enddo

          !do omega2 = MxOmegaGamG1+1, MxOmegaGamG2
            !do ib = 1, nbasis
              !GamRTailMP(ityp,dx,dy,omega1,ib) = GamRTailMP(ityp,dx,dy,omega1,ib) &
                !& + weight_basis(GamGCoefP(ib,:), omega2)* &
                !& tempGamMC(omega1,omega2)
            !enddo
          !enddo
        !enddo

        !do omega1 = MxOmegaGamG1+1, MxOmegaGamG2
          !do omega2 = -MxOmegaGamG2, -MxOmegaGamG1-1
            !do ib = 1, nbasisGamma
              !GamRTailPN(ityp,dx,dy,ib) = GamRTailPN(ityp,dx,dy,ib) &
                !& + weight_basis_Gamma(GamCoefPN(ib, :), omega1, omega2)* &
                !& ReweightBasis(omega1, omega2)* tempGamMC(omega1,omega2)
            !enddo
          !enddo

          !do omega2 = -MxOmegaGamG1, MxOmegaGamG1
            !do ib = 1, nbasis
              !GamRTailPM(ityp,dx,dy,ib,omega2) = GamRTailPM(ityp,dx,dy,ib,omega2) &
                !& + weight_basis(GamGCoefP(ib,:), omega1)* &
                !& tempGamMC(omega1,omega2)
            !enddo
          !enddo

          !do omega2 = MxOmegaGamG1+1, MxOmegaGamG2
            !if(omega1==omega2) then
              !do ib = 1, nbasis
                !GamRTailDiagP(ityp,dx,dy,ib) = GamRTailDiagP(ityp,dx,dy,ib) &
                  !& + weight_basis(GamGCoefP(ib, :), omega1)* &
                  !& tempGamMC(omega1,omega2)
              !enddo
            !else if(omega1<omega2) then
              !do ib = 1, nbasisGamma
                !GamRTailPPL(ityp,dx,dy,ib) = GamRTailPPL(ityp,dx,dy,ib) &
                  !& + weight_basis_Gamma(GamCoefPPL(ib, :), omega1, omega2)* &
                  !& ReweightBasis(omega1, omega2) *tempGamMC(omega1,omega2)
              !enddo
            !else if(omega1>omega2) then
              !do ib = 1, nbasisGamma
                !GamRTailPPR(ityp,dx,dy,ib) = GamRTailPPR(ityp,dx,dy,ib) &
                  !& + weight_basis_Gamma(GamCoefPPR(ib, :), omega1, omega2)* &
                  !& ReweightBasis(omega1, omega2) *tempGamMC(omega1,omega2)
              !enddo
            !endif
          !enddo
        !enddo
      !enddo
    !enddo
  !enddo

!END SUBROUTINE Gamma_mc2matrix_mc
 



!SUBROUTINE update_Gamma_matrix(norder)
  !implicit none
  !integer, intent(in) :: norder
  !integer :: iorder, dx, dy, ityp, omega1, omega2, ib, jb
  !double precision :: gam1, gam2, perr, norm, nmc
  !double precision :: Ga0R, GaR1
  !double precision :: tempGamMC(-MxOmegaGamG2:MxOmegaGamG2, -MxOmegaGamG2:MxOmegaGamG2)

  !tempGamMC(:, :) = 0.d0

  !GamR(:,:,:,:,:) = 0.d0
  !GamRTailC(:,:,:) = 0.d0
  !GamRTailNP(:,:,:,:) = 0.d0
  !GamRTailPN(:,:,:,:) = 0.d0
  !GamRTailPPR(:,:,:,:) = 0.d0
  !GamRTailNNR(:,:,:,:) = 0.d0
  !GamRTailPPL(:,:,:,:) = 0.d0
  !GamRTailNNL(:,:,:,:) = 0.d0

  !GamRTailPM(:,:,:,:,:) = 0.d0
  !GamRTailNM(:,:,:,:,:) = 0.d0
  !GamRTailMP(:,:,:,:,:) = 0.d0
  !GamRTailMN(:,:,:,:,:) = 0.d0
  !GamRTailDiagP(:,:,:,:) = 0.d0
  !GamRTailDiagN(:,:,:,:) = 0.d0

  !if(norder==0) then
    !call initialize_Gamma
    !return
  !endif

  !norm = GamNormWeight*ime/GamNorm

  !do ityp = 1, ntypGa
    !do dx = 0, dLx
      !do dy = 0, dLy

        !do omega1 = -MxOmegaDiag, MxOmegaDiag
          !do omega2 = -MxOmegaDiag, MxOmegaDiag
            !tempGamMC(omega1,omega2) = 0.d0
            !gam1 = SUM(GamMC(norder, :,(ityp+1)/2, dx, dy, omega1,omega2))/ime
            !gam2 = SUM(GamSqMC(norder, :,(ityp+1)/2, dx, dy, omega1, omega2))/ime
            !if(abs(gam1)>1.d-30) then
              !perr = sqrt(gam2-gam1**2.d0)/(sqrt(ime-1)*abs(gam1))
              !if(perr<=MxError) then
                !tempGamMC(omega1, omega2) = tempGamMC(omega1, omega2) + gam1*norm
              !endif
            !endif
          !enddo
        !enddo

        !do omega1 = -MxOmegaGamG2, -MxOmegaGamG1-1
          !do omega2 = -MxOmegaGamG2, -MxOmegaGamG1-1
            !if(omega1==omega2) then
              !do ib = 1, nbasis
                !GamRTailDiagN(ityp,dx,dy,ib) = GamRTailDiagN(ityp,dx,dy,ib) &
                  !& + weight_basis(GamGCoefN(ib, :), omega1)* &
                  !& tempGamMC(omega1,omega2)
              !enddo
            !else if(omega1<omega2) then
              !do ib = 1, nbasisGamma
                !GamRTailNNL(ityp,dx,dy,ib) = GamRTailNNL(ityp,dx,dy,ib) &
                  !& + weight_basis_Gamma(GamCoefNNL(ib, :), omega1, omega2)* &
                  !& ReweightBasis(omega1, omega2) *tempGamMC(omega1,omega2)
              !enddo
            !else if(omega1>omega2) then
              !do ib = 1, nbasisGamma
                !GamRTailNNR(ityp,dx,dy,ib) = GamRTailNNR(ityp,dx,dy,ib) &
                  !& + weight_basis_Gamma(GamCoefNNR(ib, :), omega1, omega2)* &
                  !& ReweightBasis(omega1, omega2) *tempGamMC(omega1,omega2)
              !enddo
            !endif
          !enddo

          !do omega2 = -MxOmegaGamG1, MxOmegaGamG1
            !do ib = 1, nbasis
              !GamRTailNM(ityp,dx,dy,ib,omega2) = GamRTailNM(ityp,dx,dy,ib,omega2) &
                !& + weight_basis(GamGCoefN(ib, :), omega1)* &
                !& tempGamMC(omega1,omega2)
            !enddo
          !enddo

          !do omega2 = MxOmegaGamG1+1, MxOmegaGamG2
            !do ib = 1, nbasisGamma
              !GamRTailNP(ityp,dx,dy,ib) = GamRTailNP(ityp,dx,dy,ib) &
                !& + weight_basis_Gamma(GamCoefNP(ib, :), omega1, omega2)* &
                !& ReweightBasis(omega1, omega2) *tempGamMC(omega1,omega2)
            !enddo
          !enddo
        !enddo

        !do omega1 = -MxOmegaGamG1, MxOmegaGamG1
          !do omega2 = -MxOmegaGamG2, -MxOmegaGamG1-1
            !do ib = 1, nbasis
              !GamRTailMN(ityp,dx,dy,omega1,ib) = GamRTailMN(ityp,dx,dy,omega1,ib) &
                !& + weight_basis(GamGCoefN(ib,:), omega2)* &
                !& tempGamMC(omega1,omega2)
            !enddo
          !enddo

          !do omega2 = -MxOmegaGamG1, MxOmegaGamG1
            !GamR(ityp,dx, dy, omega1, omega2) = GamR(ityp,dx,dy,omega1,omega2)+tempGamMC(omega1, omega2)
          !enddo

          !do omega2 = MxOmegaGamG1+1, MxOmegaGamG2
            !do ib = 1, nbasis
              !GamRTailMP(ityp,dx,dy,omega1,ib) = GamRTailMP(ityp,dx,dy,omega1,ib) &
                !& + weight_basis(GamGCoefP(ib,:), omega2)* &
                !& tempGamMC(omega1,omega2)
            !enddo
          !enddo
        !enddo

        !do omega1 = MxOmegaGamG1+1, MxOmegaGamG2
          !do omega2 = -MxOmegaGamG2, -MxOmegaGamG1-1
            !do ib = 1, nbasisGamma
              !GamRTailPN(ityp,dx,dy,ib) = GamRTailPN(ityp,dx,dy,ib) &
                !& + weight_basis_Gamma(GamCoefPN(ib, :), omega1, omega2)* &
                !& ReweightBasis(omega1, omega2)* tempGamMC(omega1,omega2)
            !enddo
          !enddo

          !do omega2 = -MxOmegaGamG1, MxOmegaGamG1
            !do ib = 1, nbasis
              !GamRTailPM(ityp,dx,dy,ib,omega2) = GamRTailPM(ityp,dx,dy,ib,omega2) &
                !& + weight_basis(GamGCoefP(ib,:), omega1)* &
                !& tempGamMC(omega1,omega2)
            !enddo
          !enddo

          !do omega2 = MxOmegaGamG1+1, MxOmegaGamG2
            !if(omega1==omega2) then
              !do ib = 1, nbasis
                !GamRTailDiagP(ityp,dx,dy,ib) = GamRTailDiagP(ityp,dx,dy,ib) &
                  !& + weight_basis(GamGCoefP(ib, :), omega1)* &
                  !& tempGamMC(omega1,omega2)
              !enddo
            !else if(omega1<omega2) then
              !do ib = 1, nbasisGamma
                !GamRTailPPL(ityp,dx,dy,ib) = GamRTailPPL(ityp,dx,dy,ib) &
                  !& + weight_basis_Gamma(GamCoefPPL(ib, :), omega1, omega2)* &
                  !& ReweightBasis(omega1, omega2) *tempGamMC(omega1,omega2)
              !enddo
            !else if(omega1>omega2) then
              !do ib = 1, nbasisGamma
                !GamRTailPPR(ityp,dx,dy,ib) = GamRTailPPR(ityp,dx,dy,ib) &
                  !& + weight_basis_Gamma(GamCoefPPR(ib, :), omega1, omega2)* &
                  !& ReweightBasis(omega1, omega2) *tempGamMC(omega1,omega2)
              !enddo
            !endif
          !enddo
        !enddo

      !enddo
    !enddo
  !enddo

!END SUBROUTINE update_Gamma_matrix
!====================================================================
