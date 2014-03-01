

!=================== The 1st order of W matrix ======================
!--------- calculate weight for W matrix ---------
SUBROUTINE calculate_W1
  implicit none
  integer :: dx, dy, omega
  integer :: ib, jb, ityp
  integer :: xg, yg, dxw0, dyw0, dxw, dyw
  double precision :: W0R, W0R1, G2R, G2W0R, WintR
  double precision :: WMR(ntypW,0:dLx,0:dLy,-MxOmegaW2:MxOmegaW2)


  !-------- W = W0/1-W0*G^2 ---------------------------------------
  WMR(:,:,:,:) = 0.d0
  G2R = weight_G2(0)

  do dx = 0, dLx
    do dy = 0, dLy
      W0R = weight_W0(dx, dy, 0, 1)
      do xg = 0, Lx-1
        do yg = 0, Ly-1
          dxw0 = abs(xg);       dyw0 = abs(yg)
          dxw  = abs(dx-xg);    dyw  = abs(dy-yg)
          W0R1 = weight_W0(dxw0, dyw0, 0, 1)
          G2W0R = G2R*W0R1/Beta
        enddo
      enddo
      WintR = W0R/(1.d0-2.d0*G2W0R)

      WMR(1, dx, dy, 0) =       WintR
      WMR(3, dx, dy, 0) = -1.d0*WintR
      WMR(5, dx, dy, 0) =  2.d0*WintR

      do ityp = 2, 6, 2
        WMR(ityp, dx, dy, 0) = WMR(ityp-1, dx, dy, 0)
      enddo
    enddo
  enddo

  !-------------- update the matrix and tail -----------------
  WRTailP(:,:,:,:) = 0.d0

  do ityp = 1, ntypW
    do dx = 0, dLx
      do dy = 0, dLy
        do omega = 0, MxOmegaW1
          W0R = weight_W0(dx, dy, omega, ityp)
          WR(ityp, dx, dy, omega) = WMR(ityp, dx, dy, omega)+W0R
        enddo
      enddo
    enddo
  enddo
  !write(*, *) WR(1, 0, 1, 0)

END SUBROUTINE calculate_W1
!====================================================================


!====================== the 1st order of Gamma ======================
!===================== need to change the type of Gamma =============
SUBROUTINE calculate_Gamma1
  implicit none
  integer :: dx, dy, omega1, omega2, dir, ityp
  integer :: xw1, yw1, xw2, yw2, omegaW, dirW
  integer :: ib,jb
  integer :: typ(3), gintyp(3), gouttyp(3), wtyp(3), ga1typ(3), ga2typ(3), ga3typ(3)
  integer :: omega3, omegain, omegaout
  double precision :: Ga0R, GinI, GoutI, WWR, GaR1, GaR2, GaR3
  double precision :: DeltaGaR
  double precision :: Gam0R, Gam1R
  double precision :: err
  double precision :: tempGamMC(-MxOmegaGamG2:MxOmegaGamG2, -MxOmegaGamG2:MxOmegaGamG2)

  Gam1MR(:,:,:,:,:) = 0.d0

  typ(1) = 1
  gintyp(1) = 1;      gouttyp(1) = 1
  ga1typ(1) = 1;      ga2typ(1)  = 1;     ga3typ(1) = 1
  wtyp(1)   = 1

  typ(2) = 3
  gintyp(2) = 2;      gouttyp(2) = 2
  ga1typ(2) = 2;      ga2typ(2)  = 5;     ga3typ(2) = 6
  wtyp(2)   = 5


  typ(3) = 5
  gintyp(3) = 1;      gouttyp(3) = 2
  ga1typ(3) = 5;      ga2typ(3)  = 1;     ga3typ(3) = 2
  wtyp(3)   = 3

  dir = 1

  do ityp = 1, 3
    do dx = 0, dLx
    if(dx/=0) cycle
    do dy = 0, dLx
    if(dy/=0) cycle
    do omega1 = -MxOmegaDiag, MxOmegaDiag
      do omega2 = -MxOmegaDiag, MxOmegaDiag

        DeltaGaR = 0.d0
        do omegaW = -MxOmegaGamInt, MxOmegaGamInt
          omegain = omega1 +omegaW
          omegaout = omega2 +omegaW
          GinI = weight_G(omegain, gintyp(ityp))
          GoutI = weight_G(omegaout, gouttyp(ityp))
          GaR1 = weight_Gamma0(dx,dy,omegain,omegaout,ga1typ(ityp))
          do xw1 = 0, Lx-1
            if(xw1/= 0) cycle
            do yw1 = 0, Ly-1
              if(yw1/= 0) cycle
              do xw2 = 0, Lx-1
                if(xw2/= 0) cycle
                do yw2 = 0, Ly-1
                  if(yw2/= 0) cycle

                  GaR2 = weight_Gamma0(xw1,yw1,omega1,omegain,ga2typ(ityp)) 
                  GaR3 = weight_Gamma0(xw2,yw2,omegaout,omega2,ga3typ(ityp)) 
                  WWR = weight_W(xw1-xw2,yw1-yw2,omegaW,wtyp(ityp))
                  DeltaGaR = DeltaGaR +GinI *GoutI *WWR *GaR1 *GaR2 *GaR3
                enddo
              enddo
            enddo
          enddo
        enddo

        Gam1R = DeltaGaR/Beta
        Gam1MR(dx,dy,typ(ityp),omega1,omega2) = Gam1R
        Gam1MR(dx,dy,typ(ityp)+1,omega1,omega2) = Gam1R
      enddo
    enddo
    enddo
    enddo
  enddo


  call initialize_Gamma

  do dx = 0, dLx
    if(dx/=0) cycle
    do dy = 0, dLy
      if(dy/=0) cycle
      do ityp = 1, ntypGa

        do omega1 = -MxOmegaGamG2, MxOmegaGamG2
          do omega2 = -MxOmegaGamG2, MxOmegaGamG2
            tempGamMC(omega1, omega2) = Gam1MR( dx, dy, ityp, omega1, omega2)
          enddo
        enddo

        do omega1 = -MxOmegaGamG2, -MxOmegaGamG1-1
          do omega2 = -MxOmegaGamG2, -MxOmegaGamG1-1
            if(omega1==omega2) then
              do ib = 1, nbasis
                GamRTailDiagN(ityp,dx,dy,ib) = GamRTailDiagN(ityp,dx,dy,ib) &
                  & + weight_basis(GamGCoefN(ib, :), omega1)* &
                  & tempGamMC(omega1,omega2)
              enddo
            else if(omega1<omega2) then
              do ib = 1, nbasisGamma
                GamRTailNNL(ityp,dx,dy,ib) = GamRTailNNL(ityp,dx,dy,ib) &
                  & + weight_basis_Gamma(GamCoefNNL(ib, :), omega1, omega2)* &
                  & ReweightBasis(omega1, omega2) *tempGamMC(omega1,omega2)
              enddo
            else if(omega1>omega2) then
              do ib = 1, nbasisGamma
                GamRTailNNR(ityp,dx,dy,ib) = GamRTailNNR(ityp,dx,dy,ib) &
                  & + weight_basis_Gamma(GamCoefNNR(ib, :), omega1, omega2)* &
                  & ReweightBasis(omega1, omega2) *tempGamMC(omega1,omega2)
              enddo
            endif
          enddo

          do omega2 = -MxOmegaGamG1, MxOmegaGamG1
            do ib = 1, nbasis
              GamRTailNM(ityp,dx,dy,ib,omega2) = GamRTailNM(ityp,dx,dy,ib,omega2) &
                & + weight_basis(GamGCoefN(ib, :), omega1)* &
                & tempGamMC(omega1,omega2)
            enddo
          enddo

          do omega2 = MxOmegaGamG1+1, MxOmegaGamG2
            do ib = 1, nbasisGamma
              GamRTailNP(ityp,dx,dy,ib) = GamRTailNP(ityp,dx,dy,ib) &
                & + weight_basis_Gamma(GamCoefNP(ib, :), omega1, omega2)* &
                & ReweightBasis(omega1, omega2) *tempGamMC(omega1,omega2)
            enddo
          enddo
        enddo

        do omega1 = -MxOmegaGamG1, MxOmegaGamG1
          do omega2 = -MxOmegaGamG2, -MxOmegaGamG1-1
            do ib = 1, nbasis
              GamRTailMN(ityp,dx,dy,omega1,ib) = GamRTailMN(ityp,dx,dy,omega1,ib) &
                & + weight_basis(GamGCoefN(ib,:), omega2)* &
                & tempGamMC(omega1,omega2)
            enddo
          enddo

          do omega2 = -MxOmegaGamG1, MxOmegaGamG1
            GamR(ityp,dx, dy, omega1, omega2) = GamR(ityp,dx,dy,omega1,omega2)+tempGamMC(omega1, omega2)
          enddo

          do omega2 = MxOmegaGamG1+1, MxOmegaGamG2
            do ib = 1, nbasis
              GamRTailMP(ityp,dx,dy,omega1,ib) = GamRTailMP(ityp,dx,dy,omega1,ib) &
                & + weight_basis(GamGCoefP(ib,:), omega2)* &
                & tempGamMC(omega1,omega2)
            enddo
          enddo
        enddo

        do omega1 = MxOmegaGamG1+1, MxOmegaGamG2
          do omega2 = -MxOmegaGamG2, -MxOmegaGamG1-1
            do ib = 1, nbasisGamma
              GamRTailPN(ityp,dx,dy,ib) = GamRTailPN(ityp,dx,dy,ib) &
                & + weight_basis_Gamma(GamCoefPN(ib, :), omega1, omega2)* &
                & ReweightBasis(omega1, omega2)* tempGamMC(omega1,omega2)
            enddo
          enddo

          do omega2 = -MxOmegaGamG1, MxOmegaGamG1
            do ib = 1, nbasis
              GamRTailPM(ityp,dx,dy,ib,omega2) = GamRTailPM(ityp,dx,dy,ib,omega2) &
                & + weight_basis(GamGCoefP(ib,:), omega1)* &
                & tempGamMC(omega1,omega2)
            enddo
          enddo

          do omega2 = MxOmegaGamG1+1, MxOmegaGamG2
            if(omega1==omega2) then
              do ib = 1, nbasis
                GamRTailDiagP(ityp,dx,dy,ib) = GamRTailDiagP(ityp,dx,dy,ib) &
                  & + weight_basis(GamGCoefP(ib, :), omega1)* &
                  & tempGamMC(omega1,omega2)
              enddo
            else if(omega1<omega2) then
              do ib = 1, nbasisGamma
                GamRTailPPL(ityp,dx,dy,ib) = GamRTailPPL(ityp,dx,dy,ib) &
                  & + weight_basis_Gamma(GamCoefPPL(ib, :), omega1, omega2)* &
                  & ReweightBasis(omega1, omega2) *tempGamMC(omega1,omega2)
              enddo
            else if(omega1>omega2) then
              do ib = 1, nbasisGamma
                GamRTailPPR(ityp,dx,dy,ib) = GamRTailPPR(ityp,dx,dy,ib) &
                  & + weight_basis_Gamma(GamCoefPPR(ib, :), omega1, omega2)* &
                  & ReweightBasis(omega1, omega2) *tempGamMC(omega1,omega2)
              enddo
            endif
          enddo
        enddo
      enddo
    enddo
  enddo

END SUBROUTINE calculate_Gamma1
!====================================================================




!====================== the 1st order of Gamma ======================
!===================== need to change the type of Gamma =============
SUBROUTINE calculate_Gamma2
  implicit none
  integer :: dx, dy, omega1, omega2, omegap,dir, ityp
  integer :: xw1, yw1, xw2, yw2, omegaW, dirW
  integer :: topo
  integer :: omega, omega3, omega4, omega5
  double precision :: G1I, G2I, G3I, G4I, W1R, W2R, GaR1, GaR2, GaR3
  double precision :: DeltaGaR
  double precision :: Gam0R, Gam1R

  Gam2MR(:,:,:) = 0.d0

  topo = 1
  do omega = 5, MxOmegaDiag
    do omegap = 5, MxOmegaDiag
      DeltaGaR = 0.d0
      do omega1 = -MxOmegaGamInt, MxOmegaGamInt
        do omega2 = -MxOmegaGamInt, MxOmegaGamInt

          omega3 = omega  - omega1
          omega4 = omega1 + omega2
          omega5 = omega2 + omega - omegap
          omegaW = omega - omegap - omega1

          G1I = weight_G(omega2, 1)
          G2I = weight_G(omega5, 1)
          G3I = weight_G(omega3, 1)
          G4I = weight_G(omega4, 1)

          W1R = weight_W(1, 0, omega1, 1)
          W2R = weight_W(1, 0, omegaW, 1)

          DeltaGaR = DeltaGaR + G1I *G2I *G3I *G4I *W1R *W2R

          omega3 = omega  + omega1
          omega4 = omega1 + omega2
          omega5 = omega2 + omegap - omega
          omegaW = omegap - omega1 - omega

          G1I = weight_G(omega2, 1)
          G2I = weight_G(omega5, 1)
          G3I = weight_G(omega3, 1)
          G4I = weight_G(omega4, 1)

          W1R = weight_W(1, 0, omega1,1)
          W2R = weight_W(1, 0, omegaW,1)

          DeltaGaR = DeltaGaR + 5.d0 *G1I *G2I *G3I *G4I *W1R *W2R
        enddo
      enddo

      Gam1R = -DeltaGaR/(Beta**2.d0)
      Gam2MR(topo, omega, omegap) = Gam1R
    enddo
  enddo

  !topo = 2
  !do omega = -MxOmegaDiag, MxOmegaDiag
    !DeltaGaR = 0.d0
    !do omega1 = -MxOmegaGamInt, MxOmegaGamInt
      !do omega2 = -MxOmegaGamInt, MxOmegaGamInt

        !omega3 = omega -omega1
        !omega4 = omega -omega1 -omega2
        !omega5 = omega -omega2

        !!G1I = weight_G(omega3, 1)
        !!G2I = weight_G(omega3, 1)
        !!G3I = weight_G(omega4, 1)
        !!G4I = weight_G(omega4, 1)

        !!W1R = weight_W(0,0,omega1,1)
        !!W2R = weight_W(0,0,omega2,1)

        !!DeltaGaR = DeltaGaR + 5.d0*G1I *G2I *G3I *G4I *W1R *W2R

        !G1I = weight_G(omega3, 1)
        !G2I = weight_G(omega4, 1)
        !G3I = weight_G(omega5, 1)
        !G4I = weight_G(omega5, 1)

        !W1R = weight_W(0,0,omega1,1)
        !W2R = weight_W(0,0,omega2,1)

        !DeltaGaR = DeltaGaR + (-1.d0)*G1I *G2I *G3I *G4I *W1R *W2R

        !!G1I = weight_G(omega3, 1)
        !!G2I = weight_G(omega3, 1)
        !!G3I = weight_G(omega4, 1)
        !!G4I = weight_G(omega5, 1)

        !!W1R = weight_W(0,0,omega1,1)
        !!W2R = weight_W(0,0,omega2,1)

        !!DeltaGaR = DeltaGaR + (-1.d0)*G1I *G2I *G3I *G4I *W1R *W2R

        !!G1I = weight_G(omega3, 1)
        !!G2I = weight_G(omega4, 1)
        !!G3I = weight_G(omega4, 1)
        !!G4I = weight_G(omega5, 1)

        !!W1R = weight_W(0,0,omega1,1)
        !!W2R = weight_W(0,0,omega2,1)

        !!DeltaGaR = DeltaGaR + 1.d0*G1I *G2I *G3I *G4I *W1R *W2R
      !enddo
    !enddo
    !Gam1R = DeltaGaR/(Beta**2.d0)
    !Gam2MR(topo,omega) = Gam1R
  !enddo
END SUBROUTINE calculate_Gamma2
!====================================================================

