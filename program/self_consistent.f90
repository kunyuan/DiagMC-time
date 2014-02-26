
!============ Initialization of the self consistent loop ==========
SUBROUTINE initialize_self_consistent
  implicit none

  !------- Initialization of the self-consistent loop ------------
  call initialize_G
  call initialize_W
  !call initialize_Gam
  !call initialize_Pi

  !!------ Initialization of W0 in p domain ----------------------
  !call initialize_W0InMoment

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

  do t = 0, MxT-1
    do dy = 0, Ly-1
      do dx = 0, Lx-1
        do typ = 1, NtypeW
          WR(typ, dx, dy, t) = weight_W0(typ, dx, dy, t)
        enddo
      enddo
    enddo
  enddo
END SUBROUTINE initialize_W

!!--------------- Initialization of Gamma -----------------------
!SUBROUTINE initialize_Gamma
  !implicit none
  !integer :: typ, omega, omega2, dx, dy
  !double precision :: Ga0R

  !GamR(:,:,:,:,:) = 0.d0
  !GamRTailC(:,:,:) = 0.d0
  
  !GamRTailPPR(:,:,:,:) = 0.d0
  !GamRTailNNR(:,:,:,:) = 0.d0
  !GamRTailPPL(:,:,:,:) = 0.d0
  !GamRTailNNL(:,:,:,:) = 0.d0

  !GamRTailDiagP(:,:,:,:) = 0.d0
  !GamRTailDiagN(:,:,:,:) = 0.d0

  !GamRTailNP(:,:,:,:) = 0.d0
  !GamRTailPN(:,:,:,:) = 0.d0

  !GamRTailMP(:,:,:,:,:) = 0.d0
  !GamRTailMN(:,:,:,:,:) = 0.d0
  !GamRTailPM(:,:,:,:,:) = 0.d0
  !GamRTailNM(:,:,:,:,:) = 0.d0

  !do typ = 1, ntypGa
    !GamRTailC(typ, 0, 0) = weight_Gamma0(0, 0, 0, 0, typ)
    !do omega = -MxOmegaGamG1, MxOmegaGamG1
      !do omega2 = -MxOmegaGamG1, MxOmegaGamG1
        !GamR(typ,0,0,omega,omega2) = weight_Gamma0(0,0, omega, omega2, typ)
      !enddo
    !enddo
  !enddo
  !return
!END SUBROUTINE initialize_Gamma
 
!!--------------- Initialization of Pi -----------------------
!SUBROUTINE initialize_Pi
  !implicit none
  !integer :: typ, omega
  !double precision :: G2R

  !PiR(:,:,:,:) = 0.d0
  !do typ = 1, ntypPi
    !G2R = weight_G2(0)
    !PiR(typ, 0, 0, 0) = G2R/Beta
  !enddo
!END SUBROUTINE initialize_Pi

!!--------------- Initialization of Pi -----------------------
!SUBROUTINE initialize_W0InMoment
  !implicit none
  !integer :: typ, dx, dy, omega
  !double precision :: W0R

  !do typ = 1, ntypW
    !do dx = 0, Lx-1
      !do dy = 0, Ly-1
        !W0InMoment(typ, dx, dy) = weight_W0(dx, dy, 0, typ)
      !enddo
    !enddo
  !enddo

  !call transfer_W0(1)
!END SUBROUTINE initialize_W0InMoment


!!============== WEIGHT CALCULATING ==================================

!!--------- calculate weight for W matrix ---------
!SUBROUTINE calculate_W
  !implicit none
  !integer :: iloop, px, py, omega, ib, jb, ityp
  !integer :: omegaGin, omegaGout
  !double precision :: W0R, W01R, GinI, GoutI, GaR1, GaR2
  !double precision :: Gin0I, Gout0I, Ga0R, G2R
  !double precision :: WintR, GGGaR
  !double precision :: WMR(ntypW,0:dLx,0:dLy,-MxOmegaW2:MxOmegaW2)
  !double precision :: WOldR
  !double precision :: err

  !!-------- calculate W = W0/(1-W0*G^2*Gamma) ----------------------------
  !WMR(:,:,:,:) = 0.d0
  !do px = 0, dLx
    !do py = 0, dLy
      !do omega = 0, MxOmegaW2
      
        !WintR = 0.d0

        !GGGaR = 0.d0
        !do omegaGin = -MxOmegaWInt, MxOmegaWInt
          !omegaGout = omegaGin -omega

          !GinI = weight_G(omegaGin, 1)
          !GoutI = weight_G(omegaGout, 1)
          !GaR1 = weight_Gamma(px, py, omegaGin, omegaGout, 1)
          !GaR2 = weight_Gamma(px, py, omegaGin, omegaGout, 3)

          !Gin0I = weight_G0up(omegaGin, 1)
          !Gout0I = weight_G0up(omegaGout, 1)
          
          !GGGaR = GGGaR + (-GinI *GoutI *(GaR1+(-1.d0)*GaR2) +Gin0I *Gout0I)/Beta
        !enddo

        !W01R = W0InMoment(1, px, py)
        !G2R = weight_G2(omega)

        !GGGaR = GGGaR + G2R/Beta
        !WintR = W01R/(1.d0-2.d0*W01R*GGGaR)

        !WMR(1, px, py, omega) =       WintR
        !WMR(3, px, py, omega) = -1.d0*WintR
        !WMR(5, px, py, omega) =  2.d0*WintR

        !do ityp = 2, 6, 2
          !WMR(ityp, px, py, omega) = WMR(ityp-1, px, py, omega)
        !enddo
      !enddo
    !enddo
  !enddo

  !!-------------- update the matrix and tail ------------
  !WRTailP(:,:,:,:) = 0.d0

  !do ityp = 1, ntypW
    !do px = 0, dLx
      !do py = 0, dLy
        !do omega = 0, MxOmegaW1
          !WR(ityp, px, py, omega) = WMR(ityp, px, py, omega)
        !enddo

        !do omega = MxOmegaW1+1, MxOmegaW2
          !do ib = 1, nbasis
            !WRTailP(ityp,px,py,ib) = WRTailP(ityp,px,py,ib) + weight_basis(WCoefP(ib,:), omega)* &
              !& (WMR(ityp,px,py,omega)-WRTailC(ityp, px, py))
          !enddo
        !enddo
      !enddo
    !enddo
  !enddo

!END SUBROUTINE calculate_W




!!--------- calculate weight for G matrix ---------
!SUBROUTINE calculate_G
  !implicit none
  !integer :: ib, pxw, pyw, omega, omegaW, omegaG1, ityp, it
  !double precision :: GintI, G0I
  !double precision :: G1I, WWR, GaR1, GaR2, GWGaI
  !double precision :: GMI(ntypG, -MxOmegaG2:MxOmegaG2)
  !double precision :: err

  !GMI(:,:) = 0.d0

  !!--------- G = G0/(1+G0*G*W*Gamma)----------------------
  !do omega = -MxOmegaG2, MxOmegaG2
    !G0I = weight_G0up(omega, 1)

    !GintI = 0.d0
    !do pxw = 0, Lx-1
      !do pyw = 0, Ly-1
        !do omegaG1 = -MxOmegaGInt, MxOmegaGInt
          !omegaW = omega-omegaG1

          !G1I = weight_G(omegaG1, 1)
          !WWR = weight_W(pxw, pyw, omegaW, 1)
          !GaR1 = weight_Gamma(pxw, pyw, omegaG1, omega, 1)
          !GaR2 = weight_Gamma(pxw, pyw, omegaG1, omega, 3)

          !GWGaI = G1I *WWR *3.d0*(GaR1+(-1.d0)*GaR2)
          !GintI = GintI + GWGaI/(Lx*Ly*Beta)

        !enddo
      !enddo
    !enddo

    !GMI(1, omega) =  G0I/(1.d0 - GintI *G0I)

    !GMI(2, omega) =  GMI(1, omega)
  !enddo



  !!-------------- update the matrix and tail ------------
  !GITailP(:,:) = 0.d0
  !GITailN(:,:) = 0.d0

  !do ityp = 1, ntypG
    !do omega = -MxOmegaG1, MxOmegaG1
      !GI(ityp, omega) = GMI(ityp, omega)
    !enddo

    !do omega = MxOmegaG1+1, MxOmegaG2
      !do ib = 1, nbasis
        !GITailP(ityp,ib) = GITailP(ityp,ib) &
          !& + weight_basis(GCoefP(ib,:),omega)*GMI(ityp,omega)
      !enddo
    !enddo

    !do omega = -MxOmegaG2, -MxOmegaG1-1
      !do ib = 1, nbasis
        !GITailN(ityp,ib) = GITailN(ityp,ib) &
          !& + weight_basis(GCoefN(ib,:),omega)*GMI(ityp,omega)
      !enddo
    !enddo
  !enddo

  !!do omega = MxOmegaG1+1, MxOmegaG2
    !!err = GMI(1,omega)
    !!do ib = 1, nbasis
      !!err = err-(weight_basis(GCoefP(ib,:),omega) &
        !!& *GITailP(1,ib))
    !!enddo
    !!write(*, *) "G", omega, err
  !!enddo
!END SUBROUTINE calculate_G




!!--------- calculate Sigma ---------
!SUBROUTINE calculate_Sigma
  !implicit none
  !integer :: ib, xw, yw, omega, omega1, omega2, omegaW, omegaG1, ityp, it
  !integer :: omegaG2, omegaG3
  !integer :: iorder
  !double precision :: GintI, G0I
  !double precision :: G1I, WWR, GaR1, GaR2, GWGaI
  !double precision :: G2I, G3I, W1R, W2R

  !SigmaI(:,:,:) = 0.d0

  !!--------- Sigma = G*W*Gamma ----------------------

  !do iorder = 0, MCOrder
    !call update_Gamma_matrix(iorder)
    !do omega = -MxOmegaSigma, MxOmegaSigma

      !GintI = 0.d0
      !do xw = 0, Lx-1
        !do yw = 0, Ly-1
          !do omegaG1 = -MxOmegaSigmaInt, MxOmegaSigmaInt
            !omegaW = omega-omegaG1

            !G1I = weight_G(omegaG1, 1)
            !WWR = weight_W(xw, yw, omegaW, 1)
            !GaR1 = weight_Gamma(xw, yw, omegaG1, omega, 1)
            !GaR2 = weight_Gamma(xw, yw, omegaG1, omega, 3)

            !GWGaI = G1I *WWR *3.d0*(GaR1+(-1.d0)*GaR2)
            !GintI = GintI - GWGaI*(Lx*Ly)/Beta
          !enddo
        !enddo
      !enddo

      !SigmaI(iorder, 1, omega) =  GintI 
      !SigmaI(iorder, 2, omega) =  SigmaI(iorder, 1, omega)
    !enddo
  !enddo
  !call Gamma_mc2matrix_mc




  !!------- test code for Sigma up to order 1 --------------------
  !!do omega = -MxOmegaSigma, MxOmegaSigma

    !!GintI = 0.d0
    !!do omegaG1 = -MxOmegaSigmaInt, MxOmegaSigmaInt
      !!omegaW = omega-omegaG1

      !!G1I = weight_G(omegaG1, 1)
      !!WWR = weight_W(0, 0, omegaW, 1)

      !!GWGaI = G1I *WWR *3.d0
      !!GintI = GintI - GWGaI*(Lx*Ly)/Beta
    !!enddo
    !!SigmaI(0, 1, omega) =  GintI 
    !!SigmaI(0, 2, omega) =  SigmaI(0, 1, omega)
  !!enddo



  !!do omega = -MxOmegaSigma, MxOmegaSigma

    !!GintI = 0.d0
    !!do omega1 = -MxOmegaSigmaInt, MxOmegaSigmaInt
      !!do omega2 = -MxOmegaSigmaInt, MxOmegaSigmaInt
        !!omegaG1 = omega-omega1
        !!omegaG2 = omega-omega1-omega2
        !!omegaG3 = omega-omega2

        !!G1I = weight_G(omegaG1, 1)
        !!G2I = weight_G(omegaG2, 1)
        !!G3I = weight_G(omegaG3, 1)

        !!W1R = weight_W(0, 0, omega1, 1)
        !!W2R = weight_W(0, 0, omega2, 1)

        !!GWGaI = G1I *G2I *G3I *W1R *W2R *3.d0
        !!GintI = GintI + GWGaI*(Lx*Ly)/(Beta**2.d0)
      !!enddo
    !!enddo
    !!SigmaI(1, 1, omega) =  GintI 
    !!SigmaI(1, 2, omega) =  SigmaI(1, 1, omega)
  !!enddo


  !return
!END SUBROUTINE calculate_Sigma



!!------------- calculate Pi -------------------------------------------
!SUBROUTINE calculate_Pi
  !implicit none
  !integer :: dx, dy, ib, jb, ityp
  !integer :: omega, omegain, omegaout
  !double precision :: G01I, G02I, Ga0R, G1I, G2I, GaR1, GaR2, Ga0R2
  !double precision :: G2R, Pi1R 

  !!-------- calculate Pi = G^2 * Gamma --------------------
  !PiR(:,:,:,:) = 0.d0

  !do dx = 0, dLx
    !do dy = 0, dLy
      !do omega = -MxOmegaChi, MxOmegaChi

        !do omegain = -MxOmegaPiInt,  MxOmegaPiInt
          !omegaout = omegain-omega

          !G1I = weight_G(omegain, 1)
          !G2I = weight_G(omegaout, 1)

          !G01I = weight_G0up(omegain, 1)
          !G02I = weight_G0up(omegaout, 1)

          !GaR1 = weight_Gamma(dx, dy, omegain, omegaout, 1)
          !Ga0R = weight_Gamma0(dx, dy, omegain, omegaout, 1)

          !Pi1R = -(GaR1 *G1I *G2I - Ga0R *G01I *G02I)/Beta
          !PiR(1, dx, dy, omega) = PiR(1, dx, dy, omega) + Pi1R

          !GaR2 = weight_Gamma(dx, dy, omegain, omegaout, 3)
          !Ga0R2 = weight_Gamma0(dx, dy, omegain, omegaout, 3)

          !Pi1R = -(GaR2 *G1I *G2I - Ga0R2 *G01I *G02I)/Beta
          !PiR(3, dx, dy, omega) = PiR(3, dx, dy, omega) + Pi1R
        !enddo

        !G2R = weight_G2(omega)
        !if(dx==0 .and. dy==0) then
          !PiR(1, dx, dy, omega) = PiR(1, dx, dy, omega) + G2R/Beta
        !endif

        !do ityp = 2, 4, 2
          !PiR(ityp, dx, dy, omega) = PiR(ityp-1, dx, dy, omega)
        !enddo

      !enddo
    !enddo
  !enddo
!END SUBROUTINE calculate_Pi



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

!--------- weight for 1-order propagator ----------------
Complex FUNCTION weight_G0(typ, t)
  implicit none
  integer, intent(in)  :: typ, t  
  double precision     :: tau
  complex(kind=8)      :: muc  

  muc = cmplx(0.d0, Mu(1)*pi/(2.d0*Beta))
  tau = t*Beta/MxT
  if(tau>=0) then
    weight_G0 = exp(muc*tau)/(1.d0, 1.d0) 
  else if(tau>=-MxT) then
    weight_G0 = -exp(muc*(tau+Beta))/(1.d0, 1.d0) 
  endif
  write(*, *) weight_G0
  return
END FUNCTION weight_G0

!!---------- G^2 -----------------------------------------
!DOUBLE PRECISION FUNCTION weight_G2(omega1)
  !implicit none
  !integer, intent(in)  :: omega1  
  !if(omega1 == 0 ) then
    !weight_G2 = -0.5d0*Beta**2.d0
  !else
    !weight_G2 = 0.d0
  !endif
!END FUNCTION weight_G2

!!--------- calculate weight for 1-order interaction ----
Complex FUNCTION weight_W0(typ, dx, dy, t)
  implicit none
  integer, intent(in) :: dx, dy, typ, t
  integer :: dx1, dy1

  dx1 = dx;       dy1 = dy
  if(dx1>=0  .and. dx1<Lx .and. dy1>=0 .and. dy1<Ly) then
    if(dx1>dLx)     dx1 = Lx-dx1
    if(dy1>dLy)     dy1 = Ly-dy1

    weight_W0 = (0.d0, 0.d0)

    if(t==0) then
      if((dx1==1.and.dy1==0).or.(dx1==0.and.dy1==1)) then
        if(typ1 ==1 .or. typ1 == 2) then
          W0R = cmplx(0.25d0*Jcp, 0.d0)
        else if(typ1 == 3 .or. typ1 == 4) then
          W0R = cmplx(-0.25d0*Jcp, 0.d0)
        else if(typ1 == 5 .or. typ1 == 6) then
          W0R = cmplx(0.5d0*Jcp, 0.d0)
        endif
      endif
    endif
  else
    write(*, *) dx1, dy1, "dx, dy bigger than system size!"
    stop
  endif

  weight_W0 = W0R
END FUNCTION weight_W0

!!--------- calculate weight for 1-order Gamma ---------
!DOUBLE PRECISION FUNCTION weight_Gamma0(dx1, dy1, omega1, omega2, typ1)
  !implicit none
  !integer, intent(in)  :: dx1, dy1, omega1, omega2, typ1
  !double precision :: Gamma0R
  !integer :: dx, dy

  !dx = dx1;      dy = dy1
  !if(dx>=0 .and. dx<Lx .and. dy>=0 .and. dy<Ly) then
    !if(dx>dLx)     dx = Lx-dx
    !if(dy>dLy)     dy = Ly-dy

    !if(dx==0.and.dy==0) then
      !if(typ1==1 .or. typ1==2 .or. typ1==5 .or. typ1==6) then
        !Gamma0R = 1.d0
      !else
        !Gamma0R = 0.d0
      !endif
    !else
      !Gamma0R = 0.d0
    !endif
  !else
    !write(*, *) dx1, dy1, "dx, dy bigger than system size!"
    !stop
  !endif
  !weight_Gamma0 = Gamma0R
!END FUNCTION weight_Gamma0

!!--------- extract weight for G ---------
!DOUBLE PRECISION FUNCTION weight_G(omega1, typ1)
  !implicit none
  !integer, intent(in)  :: omega1, typ1
  !double precision :: muc, omegac
  !double precision:: GGI
  !integer :: ib
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
!END FUNCTION weight_G

!!--------- extract weight for W ---------
!DOUBLE PRECISION FUNCTION weight_W(dx1, dy1, omega1, typ1)
  !implicit none
  !integer, intent(in)  :: dx1, dy1, omega1, typ1
  !double precision :: WWR
  !integer :: dx, dy, ib

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
!END FUNCTION weight_W

!!--------- extract weight for Gamma ---------
!DOUBLE PRECISION FUNCTION weight_Gamma(dx1, dy1, omega1, omega2, typ1)
  !implicit none
  !integer, intent(in)  :: dx1, dy1, omega1, omega2, typ1
  !double precision :: GaR
  !integer :: ib, jb, dx, dy

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
!END FUNCTION weight_Gamma

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
