
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

  W(:,:,:) = (0.d0, 0.d0)
END SUBROUTINE initialize_W

!!--------------- Initialization of Gamma -----------------------
SUBROUTINE initialize_Gam
  implicit none
  integer :: ityp, it1, it2

  Gam(:,:,:,:) = (0.d0, 0.d0)
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
  integer :: site
  double precision :: ratio

  W0PF(:,:) = (0.d0, 0.d0)
  ratio = real(MxT)/Beta
  do site = 0, Vol-1
    W0PF(site, 0:MxT-1) = ratio *weight_W0(1, site)
  enddo
  call transfer_W0_r(1)
END SUBROUTINE initialize_W0PF

!!------- Initialization of Gam0 in monmentum and frequency ----------
SUBROUTINE initialize_Gam0PF
  implicit none
  double precision :: ratio

  ratio = (real(MxT)/Beta)**2.d0
  Gam0PF(0:Vol-1,0:MxT-1,0:MxT-1) = ratio*weight_Gam0(1, 0)
END SUBROUTINE initialize_Gam0PF
 
!!============== WEIGHT CALCULATING ==================================

SUBROUTINE calculate_Polar
  implicit none
  integer :: p, omega
  integer :: omegaGin, omegaGout
  complex(kind=8) :: Gin, Gout, Gam1
  double precision :: ratio

  Polar(:,:) = (0.d0, 0.d0)

  ratio = 2.d0/real(MxT)*(Beta/real(MxT))**4.d0
  do omega = 0, MxT-1
    do omegaGin = 0, MxT-1
      do p = 0, Vol-1
        omegaGout = omegaGin - omega

        Gin = G(1, omegaGin)
        if(omegaGout>=0) then
          Gout = G(1, omegaGout)
          Gam1 = Gam(5, p, omegaGin, omegaGout)
        else
          Gout = -1.d0*G(1, omegaGout+MxT)
          Gam1 = -1.d0*Gam(5, p, omegaGin, omegaGout+MxT)
        endif
        Polar(p, omega) = Polar(p, omega)+d_times_cd(ratio, cdexp(dcmplx(0.d0, -1.d0)* &
          & omegaGout*2.d0*Pi/MxT) *Gin*Gout*Gam1)
      enddo
    enddo
  enddo
  return
END SUBROUTINE calculate_Polar


SUBROUTINE calculate_Sigma
  implicit none
  integer :: p, omega
  integer :: omegaG, omegaW
  complex(kind=8) :: G1, W1, Gam1
  double precision :: ratio

  Sigma(:) = (0.d0, 0.d0)

  ratio = -3.d0/(Vol*real(MxT))*(Beta/real(MxT))**4.d0
  do omega = 0, MxT-1
    do omegaG = 0, MxT-1
      do p = 0, Vol-1
        omegaW = omega-omegaG

        G1 = G(1, omegaG)
        if(omegaW>=0) then
          W1 = W(1, p, omegaW)
        else
          W1 = W(1, p, omegaW+MxT)
        endif
        Gam1 = Gam(5, p, omegaG, omega)
        
        Sigma(omega) = Sigma(omega)+d_times_cd(ratio, G1*W1*Gam1)
      enddo
    enddo
  enddo
  return
END SUBROUTINE calculate_Sigma


!!--------- calculate weight for W matrix ---------
SUBROUTINE calculate_W
  implicit none
  integer :: omega, p

  !-------- calculate W = W0/(1-W0*G^2*Gamma) ----------------------------
  W(:,:,:) = (0.d0, 0.d0)

  do  omega = 0, MxT-1
    do p = 0, Vol-1
      Denom(p,omega) = (1.d0, 0.d0) -W0PF(p,omega)*Polar(p,omega)
      W(1,p,omega) = W0PF(p,omega)/Denom(p,omega)
      if(W(1,p,omega)/=W(1,p,omega)) then
        call LogFile%QuickLog("calculate_W NaN appears!")
        stop
      endif
      W(3,p,omega) = d_times_cd(-1.d0,W(1,p,omega))
      W(5,p,omega) = d_times_cd( 2.d0,W(1,p,omega))
    enddo
  enddo

  W(2,:,:) = W(1,:,:)
  W(4,:,:) = W(3,:,:)
  W(6,:,:) = W(5,:,:)

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
  Denom(:,:) = (1.d0, 0.d0) -W0PF(:,:)*Polar(:,:)
  return
END SUBROUTINE calculate_Denom



!!------------- calculate the susceptibility ---------------------------
SUBROUTINE calculate_Chi
  implicit none
  integer :: p, omega
  double precision :: ratio

  !!-------- calculate Chi = Pi/(1 - W0 * Pi) ------------------
  !!--------- already sum over the spins -----------------------
  ratio = -1.d0*(real(MxT)/Beta)**2.d0*0.75d0

  do  omega = 0, MxT-1
    do p = 0, Vol-1
      Chi(p,omega) = Polar(p,omega)/Denom(p,omega)
      Chi(p,omega) = d_times_cd(ratio, Chi(p,omega))
    enddo
  enddo
  return
END SUBROUTINE calculate_Chi

!!====================================================================

SUBROUTINE plus_minus_W0(Backforth)
  implicit none
  integer, intent(in) :: Backforth

  if(Backforth/=-1) then
    W(1,:,:) = W(1,:,:) + W0PF(:,:)
  else 
    W(1,:,:) = W(1,:,:) - W0PF(:,:)
  endif

  W(3,:,:) = -1.d0 * W(1,:,:)
  W(5,:,:) =  2.d0 * W(1,:,:)

  W(2,:,:) = W(1,:,:)
  W(4,:,:) = W(3,:,:)
  W(6,:,:) = W(5,:,:)

  return
END SUBROUTINE

SUBROUTINE plus_minus_Gam0(Backforth)
  implicit none
  integer, intent(in) :: Backforth
  integer :: ityp, px, py, omega1, omega2

  if(Backforth/=-1) then
    Gam(1,:,:,:) = Gam(1,:,:,:) + Gam0PF(:,:,:)
    Gam(2,:,:,:) = Gam(2,:,:,:) + Gam0PF(:,:,:)
    Gam(5,:,:,:) = Gam(5,:,:,:) + Gam0PF(:,:,:)
    Gam(6,:,:,:) = Gam(6,:,:,:) + Gam0PF(:,:,:)
  else 
    Gam(1,:,:,:) = Gam(1,:,:,:) - Gam0PF(:,:,:)
    Gam(2,:,:,:) = Gam(2,:,:,:) - Gam0PF(:,:,:)
    Gam(5,:,:,:) = Gam(5,:,:,:) - Gam0PF(:,:,:)
    Gam(6,:,:,:) = Gam(6,:,:,:) - Gam0PF(:,:,:)
  endif

  return
END SUBROUTINE plus_minus_Gam0

 
SUBROUTINE Gam_mc2matrix_mc
  implicit none
  integer :: iorder, dr, ityp, iloop, it1, it2, typ, drr, itt1, itt2
  integer :: ibin, ibasis
  complex*16 :: cgam, normal
  logical :: flag(MxOrder)
  double precision :: totrerr, totierr, tau1, tau2
  double precision :: rgam, igam, rgam2, igam2, rerr, ierr, rpercenterr, ipercenterr

  call initialize_Gam

  normal = GamNormWeight/GamNorm

  call LogFile%QuickLog("(ErrorRatio):"+str(ratioerr))

  flag(:) = .true.

  call LogFile%WriteStamp('i')

  looporder: do iorder = 1, MCOrder
    totrerr = 0.d0
    totierr = 0.d0
    do it1 = 0, MxT-1
      rgam = real(GamMC(iorder, 0, it1))/Z_normal
      rgam2 = ReGamSqMC(iorder, 0, it1)/Z_normal

      rerr = sqrt(abs(rgam2)-rgam**2.d0)/sqrt(Z_normal-1)
      rerr = rerr* ratioerr
      totrerr = totrerr + abs(rerr)

      igam = dimag(GamMC(iorder, 0, it1))/Z_normal
      igam2 = ImGamSqMC(iorder, 0, it1)/Z_normal

      ierr = sqrt(abs(igam2)-igam**2.d0)/sqrt(Z_normal-1)
      ierr = ierr* ratioerr
      totierr = totierr + abs(ierr)
    enddo

    rgam = SUM(abs(real(GamMC(iorder, 0, :))))/Z_normal
    if(rgam<1.d-30) then
      rpercenterr = 0.d0
    else
      rpercenterr = totrerr/rgam
    endif

    if(rpercenterr>MxError) then
      flag(iorder)=.false.
    endif

    igam = SUM(abs(dimag(GamMC(iorder, 0, :))))/Z_normal
    if(igam<1.d-30) then
      ipercenterr = 0.d0
    else
      ipercenterr = totierr/igam
    endif

    if(ipercenterr>MxError) then
      flag(iorder)=.false.
    endif

    call LogFile%WriteLine("Order "+str(iorder)+" relative error: "+str(rpercenterr,'(f6.3)')+","+str(ipercenterr,'(f6.3)')+", accept: "+str(flag(iorder)))


    if(flag(iorder)) then
      do it2 = 0, MxT-1
        do it1 = 0, MxT-1
          tau1 = dble(it1)*Beta/dble(MxT)
          tau2 = dble(it2)*Beta/dble(MxT)
          ibin = get_bin_Gam(it1, it2)

          do dr = 0, Vol-1
            do ityp = 1, NTypeGam/2
              typ = 2*(ityp-1)+1

              if(IsBasis2D(ibin)) then
                do ibasis = 1, NBasisGam
                  cgam = GamBasis(iorder,ityp,dr,ibin,ibasis)* weight_basis_Gam( &
                    & CoefGam(0:BasisOrderGam,0:BasisOrderGam,ibasis,ibin), tau1, tau2)
                  Gam(typ,dr,it1,it2) = Gam(typ,dr,it1,it2) + cgam*normal
                enddo
              else 
                do ibasis = 1, NBasis
                  cgam = GamBasis(iorder, ityp, dr, ibin, ibasis)* weight_basis( &
                    & CoefGam(0:BasisOrder,0,ibasis,ibin), tau1)
                  Gam(typ,dr,it1,it2) = Gam(typ,dr,it1,it2) + cgam*normal
                enddo
              endif
            enddo
          enddo
        enddo
      enddo
    endif
  enddo looporder

  Gam(2,:,:,:) = Gam(1,:,:,:)
  Gam(4,:,:,:) = Gam(3,:,:,:)
  Gam(6,:,:,:) = Gam(5,:,:,:)

END SUBROUTINE Gam_mc2matrix_mc
