
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
  W0PF(L(1)-1, 0, 0) = dcmplx(ratio, 0.d0)
  W0PF(0,    1, 0) = dcmplx(ratio, 0.d0)
  W0PF(0, L(2)-1, 0) = dcmplx(ratio, 0.d0)

  call transfer_W0_r(1)
  call transfer_W0_t(1)
END SUBROUTINE initialize_W0PF

!!------- Initialization of Gam0 in monmentum and frequency ----------
SUBROUTINE initialize_Gam0PF
  implicit none
  double precision :: ratio

  ratio = (real(MxT)/Beta)**2.d0
  Gam0PF(:,:,:,:) = (0.d0, 0.d0)
  Gam0PF(0,0,0,0) = dcmplx(ratio, 0.d0)

  call transfer_Gam0_r(1)
  call transfer_Gam0_t(1)
END SUBROUTINE initialize_Gam0PF
 
!!============== WEIGHT CALCULATING ==================================

SUBROUTINE calculate_Polar
  implicit none
  integer :: px, py, p(2), omega
  integer :: omegaGin, omegaGout
  complex(kind=8) :: Gin, Gout, Gam1
  double precision :: ratio

  Polar(:,:,:) = (0.d0, 0.d0)

  ratio = 2.d0/real(MxT)*(Beta/real(MxT))**4.d0
  do omega = 0, MxT-1
    do omegaGin = 0, MxT-1
      do py = 0, L(2)-1
        do px = 0, L(1)-1
          p = (/px, py/)
          omegaGout = omegaGin - omega

          Gin = weight_G(1, omegaGin)
          if(omegaGout>=0) then
            Gout = weight_G(1, omegaGout)
            Gam1 = weight_Gam(5, p, omegaGin, omegaGout)
          else
            Gout = weight_G(1, omegaGout+MxT)
            Gam1 = weight_Gam(5, p, omegaGin, omegaGout+MxT)
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
  integer :: px, py, p(2), omega
  integer :: omegaG, omegaW
  complex(kind=8) :: G1, W1, Gam1
  double precision :: ratio

  Sigma(:) = (0.d0, 0.d0)

  ratio = -3.d0/(real(L(1))*real(L(2))*real(MxT))*(Beta/real(MxT))**4.d0
  do omega = 0, MxT-1
    do omegaG = 0, MxT-1
      do py = 0, L(2)-1
        do px = 0, L(1)-1
          p = (/px, py/)
          omegaW = omega-omegaG

          G1 = weight_G(1, omegaG)
          if(omegaW>=0) then
            W1 = weight_W(1, p, omegaW)
          else
            W1 = weight_W(1, p, omegaW+MxT)
          endif
          Gam1 = weight_Gam(5, p, omegaG, omega)
          
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

  do  omega = 0, MxT-1
    do py = 0, L(2)-1
      do px = 0, L(1)-1
        Denom(px,py,omega) = (1.d0, 0.d0) -W0PF(px,py,omega)*Polar(px,py,omega)
        W(1,px,py,omega) = W0PF(px,py,omega)/Denom(px,py,omega)
        if(W(1,px,py,omega)/=W(1,px,py,omega)) then
          call LogFile%QuickLog("calculate_W NaN appears!")
          stop
        endif
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
  integer :: omega

  G(:,:) = (0.d0, 0.d0)

  !--------- G = G0/(1-G0*Sigma)----------------------
  do omega = 0, MxT-1

    G(1, omega) =  G0F(omega)/((1.d0,0.d0)-G0F(omega)*Sigma(omega))

    if(G(1,omega)/=G(1,omega)) then
      call LogFile%QuickLog("calculate_G NaN appears!")
      stop
    endif

    G(2, omega) =  G(1, omega)
  enddo

  !!-------------- update the matrix and tail ------------
END SUBROUTINE calculate_G


!!------------- calculate the denominator of the susceptibility ---------
SUBROUTINE calculate_Denom
  implicit none
  Denom(:,:,:) = (1.d0, 0.d0) -W0PF(:,:,:)*Polar(:,:,:)
  return
END SUBROUTINE calculate_Denom



!!------------- calculate the susceptibility ---------------------------
SUBROUTINE calculate_Chi
  implicit none
  integer :: px, py, omega
  double precision :: ratio

  !!-------- calculate Chi = Pi/(1 - W0 * Pi) ------------------
  !!--------- already sum over the spins -----------------------
  ratio = -1.d0*(real(MxT)/Beta)**2.d0*0.75d0

  do  omega = 0, MxT-1
    do py = 0, L(2)-1
      do px = 0, L(1)-1
        Chi(px,py,omega) = Polar(px,py,omega)/Denom(px,py,omega)
        Chi(px,py,omega) = d_times_cd(ratio, Chi(px,py,omega))
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
    do py = 0, L(2)-1
      do px = 0, L(1)-1
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





SUBROUTINE Gam_mc2matrix_mc
  implicit none
  integer :: iorder, dx, dy, ityp, iloop, it1, it2, typ
  complex*16 :: cgam, normal
  logical :: flag(MxOrder)
  double precision :: rgam2, igam2, rerr, ierr, rpercenterr, ipercenterr

  call initialize_Gam

  normal = GamNormWeight*Z_normal/GamNorm

  flag(:) = .true.
  do it2 = 0, MxT-1
    do it1 = 0, MxT-1
      do iorder = 1, MCOrder

        cgam = GamMC(iorder, 1, 0, 0,it1,it2) /Z_normal

        rgam2 = ReGamSqMC(iorder, 1, 0, 0, it1, it2)/Z_normal
        rerr = sqrt(abs(rgam2)-(real(cgam))**2.d0)/sqrt(Z_normal-1)
        rerr = rerr* ratioerr

        if(abs(real(cgam))<1.d-30) then
          rpercenterr = 0.d0
        else
          rpercenterr = rerr/abs(real(cgam))
        endif
        if(rpercenterr>MxError)  flag(iorder)=.false.

        igam2 = ImGamSqMC(iorder, 1, 0, 0, it1, it2)/Z_normal
        ierr = sqrt(abs(igam2)-(dimag(cgam))**2.d0)/sqrt(Z_normal-1)
        ierr = ierr* ratioerr

        if(abs(dimag(cgam))<1.d-30) then
          ipercenterr = 0.d0
        else
          ipercenterr = ierr/abs(dimag(cgam))
        endif
        if(ipercenterr>MxError) flag(iorder)=.false.
      enddo
    enddo
  enddo

  do iorder = 1, MCOrder
    if(flag(iorder)) then
      do it2 = 0, MxT-1
        do it1 = 0, MxT-1
          do dx = 0, L(1)-1
            do dy = 0, L(2)-1
              do ityp = 1, NTypeGam/2
                cgam = GamMC(iorder,ityp,dx,dy,it1,it2)/Z_normal
                typ = 2*(ityp-1)+1
                Gam(typ,dx,dy,it1,it2) = Gam(typ,dx,dy,it1,it2)+ cgam*normal
              enddo
            enddo
          enddo
        enddo
      enddo
    endif

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
    !do dx = 0, dL(1)
      !do dy = 0, dL(2)


      !enddo
    !enddo
  !enddo

!END SUBROUTINE update_Gamma_matrix
!====================================================================
