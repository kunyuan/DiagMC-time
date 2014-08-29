
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
SUBROUTINE calculate_W(iloop)
  implicit none
  integer, intent(in) :: iloop
  integer :: omega, p
  double precision :: ratio

  !-------- calculate W = W0/(1-W0*G^2*Gamma) ----------------------------
  newW(:,:,:) = (0.d0, 0.d0)

  do  omega = 0, MxT-1
    do p = 0, Vol-1
      Denom(p,omega) = (1.d0, 0.d0) -W0PF(p,omega)*Polar(p,omega)
      newW(1,p,omega) = W0PF(p,omega)/Denom(p,omega)
      if(newW(1,p,omega)/=newW(1,p,omega)) then
        call LogFile%QuickLog("calculate_W NaN appears!")
        stop
      endif
      newW(3,p,omega) = d_times_cd(-1.d0,newW(1,p,omega))
      newW(5,p,omega) = d_times_cd( 2.d0,newW(1,p,omega))
    enddo
  enddo

  newW(2,:,:) = newW(1,:,:)
  newW(4,:,:) = newW(3,:,:)
  newW(6,:,:) = newW(5,:,:)

  ratio = dble(iloop)/(1.d0+dble(iloop))
  W(:,:,:) = ratio*W(:,:,:)+(1.d0-ratio)*newW(:,:,:)

END SUBROUTINE calculate_W

!!--------- calculate weight for G matrix ---------
SUBROUTINE calculate_G(iloop)
  implicit none
  integer, intent(in) :: iloop
  complex(kind=8) :: G0
  integer :: omega
  double precision :: ratio

  newG(:,:) = (0.d0, 0.d0)

  !--------- G = G0/(1-G0*Sigma)----------------------
  do omega = 0, MxT-1

    newG(1, omega) =  G0F(omega)/((1.d0,0.d0)-G0F(omega)*Sigma(omega))

    if(newG(1,omega)/=newG(1,omega)) then
      call LogFile%QuickLog("calculate_G NaN appears!")
      stop
    endif

    newG(2, omega) =  newG(1, omega)
  enddo

  ratio = dble(iloop)/(1.d0+dble(iloop))
  G(:, :) = ratio*G(:, :) + (1.d0-ratio)*newG(:, :)
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

 
SUBROUTINE Gam_mc2matrix_mc(changeBeta)
  implicit none
  logical, intent(out) :: changeBeta
  integer :: iorder, ir, dr, ityp, iloop, it1, it2, typ, drr, itt1, itt2
  integer :: ibin, ibasis
  complex*16 :: cgam, normal
  logical :: flag(MxOrder)
  double precision :: totrerr, totierr, tau1, tau2
  double precision :: rgam, igam, rgam2, igam2, rerr, ierr, rpercenterr, ipercenterr

  call LogFile%QuickLog("(ErrorRatio):"+str(ratioerr))

  flag(:) = .false.

  call LogFile%WriteStamp('i')

  normal = GamNormWeight/GamNorm

  GamBasis = (0.d0, 0.d0)

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

    if(rpercenterr<=MxError) then
      flag(iorder)=.true.
    endif

    igam = SUM(abs(dimag(GamMC(iorder, 0, :))))/Z_normal
    if(igam<1.d-30) then
      ipercenterr = 0.d0
    else
      ipercenterr = totierr/igam
    endif

    if(ipercenterr<=MxError) then
      flag(iorder)=.true.
    endif

    call LogFile%WriteLine("Order "+str(iorder)+" relative error: "+str(rpercenterr,'(f6.3)')+","+str(ipercenterr,'(f6.3)')+", accept: "+str(flag(iorder)))


    if(flag(iorder)) then
      do ityp = 1, NTypeGam/2
        typ = 2*(ityp-1) + 1
        do dr = 0, Vol-1
          ir = diff_r(D, dr, 0)
          GamBasis(typ,dr,:,:) = GamBasis(typ,dr,:,:) + normal*GamMCBasis(iorder,ityp,ir,:,:)
        enddo
      enddo
    else 
      exit looporder
    endif
  enddo looporder

  GamBasis(2,:,:,:) = GamBasis(1,:,:,:)
  GamBasis(4,:,:,:) = GamBasis(3,:,:,:)
  GamBasis(6,:,:,:) = GamBasis(5,:,:,:)

  changeBeta = .false.
  if((MCOrder<3 .and. flag(MCOrder)) .or.(MCOrder>=3 .and. flag(3))) then
    changeBeta = .true.
  endif

  call Gam_basis2matrix

END SUBROUTINE Gam_mc2matrix_mc

SUBROUTINE Gam_mc2matrix_mc_by_order(iorder)
  implicit none
  integer, intent(in) :: iorder
  integer :: ir, dr, ityp, typ
  complex*16 :: normal

  normal = GamNormWeight/GamNorm

  GamBasis = (0.d0, 0.d0)

  do ityp = 1, NTypeGam/2
    typ = 2*(ityp-1) + 1
    do dr = 0, Vol-1
      ir = diff_r(D, dr, 0)
      GamBasis(typ,dr,:,:) = normal*GamMCBasis(iorder,ityp,ir,:,:)
    enddo
  enddo

  GamBasis(2,:,:,:) = GamBasis(1,:,:,:)
  GamBasis(4,:,:,:) = GamBasis(3,:,:,:)
  GamBasis(6,:,:,:) = GamBasis(5,:,:,:)

  call Gam_basis2matrix
END SUBROUTINE Gam_mc2matrix_mc_by_order

SUBROUTINE Gam_mc2matrix_mc_up2order(iorder)
  implicit none
  integer, intent(in) :: iorder
  integer :: ir, dr, ityp, typ, i
  complex*16 :: normal

  normal = GamNormWeight/GamNorm

  GamBasis = (0.d0, 0.d0)

  do ityp = 1, NTypeGam/2
    typ = 2*(ityp-1) + 1
    do dr = 0, Vol-1
      ir = diff_r(D, dr, 0)
      do i = 1, iorder
        GamBasis(typ,dr,:,:) = GamBasis(typ,dr,:,:) + normal*GamMCBasis(i,ityp,ir,:,:)
      enddo
    enddo
  enddo

  GamBasis(2,:,:,:) = GamBasis(1,:,:,:)
  GamBasis(4,:,:,:) = GamBasis(3,:,:,:)
  GamBasis(6,:,:,:) = GamBasis(5,:,:,:)

  call Gam_basis2matrix
END SUBROUTINE Gam_mc2matrix_mc_up2order

SUBROUTINE Gam_basis2matrix
  implicit none
  integer :: ityp, ir, it1, it2

  call initialize_Gam

  do it2 = 0, MxT-1
    do it1 = 0, MxT-1
      do ir = 0, Vol-1
        do ityp = 1, NtypeGam
          Gam(ityp, ir, it1, it2) = Gam_basis(it1, it2, GamBasis(ityp, ir, :,:))
        enddo
      enddo
    enddo
  enddo
  return
END SUBROUTINE Gam_basis2matrix

COMPLEX*16 FUNCTION Gam_basis(it1, it2, GammaBasis)
  implicit none
  integer, intent(in) :: it1, it2
  complex*16, intent(in) :: GammaBasis(1:NbinGam, 1:NBasisGam)
  double precision :: tau1, tau2
  integer :: ibin, ibasis, jt1, jt2
  complex*16 :: cgam

  ibin = get_bin_Gam(it1, it2)

  tau1 = (dble(it1)+0.5d0)*Beta/dble(MxT)
  tau2 = (dble(it2)+0.5d0)*Beta/dble(MxT)

  if(ibin==1) then
    cgam = (0.d0, 0.d0)
    do ibasis = 1, NBasisGam
      cgam = cgam + GammaBasis(ibin,ibasis)* weight_basis_Gam( &
        & CoefGam(0:BasisOrderGam,0:BasisOrderGam,ibasis,ibin), it1, it2)
    enddo
    Gam_basis = cgam
  else if(ibin==2) then
    jt1 = MxT-1-it1
    jt2 = MxT-1-it2
    cgam = (0.d0, 0.d0)
    do ibasis = 1, NBasisGam
      cgam = cgam + GammaBasis(1,ibasis)* weight_basis_Gam( &
        & CoefGam(0:BasisOrderGam,0:BasisOrderGam,ibasis,1), jt1, jt2)
    enddo
    Gam_basis = dcmplx(real(cgam), -dimag(cgam))
  endif
  return
END FUNCTION Gam_basis

