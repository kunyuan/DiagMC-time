
!============  the diagonal Gamma ================================
SUBROUTINE accumulate_diagonal_Gamma(dt1, val, dGam, RedGam, ImdGam)
  implicit none
  integer, intent(in) :: dt1
  complex*16, intent(in) :: val
  complex*16 :: dGam(0:MxT-1)
  double precision :: RedGam(0:MxT-1), ImdGam(0:MxT-1)
  integer :: dt

  dt = dt1
  dGam(dt) = dGam(dt) + dcmplx(real(val),  dimag(val))
  RedGam(dt) = RedGam(dt) + (Real(val))**2.d0
  ImdGam(dt) = ImdGam(dt) + (dimag(val))**2.d0

  dt = Mxt-1-dt1
  dGam(dt) = dGam(dt) + dcmplx(real(val),  -1.d0*dimag(val))
  RedGam(dt) = RedGam(dt) + (Real(val))**2.d0
  ImdGam(dt) = ImdGam(dt) + (dimag(val))**2.d0

  return
END SUBROUTINE accumulate_diagonal_Gamma




SUBROUTINE accumulate_Gamma_basis(dt1, dt2, val, basisGam, RebasisGam, ImbasisGam)
  implicit none
  integer, intent(in) :: dt1, dt2
  complex*16, intent(in) :: val
  complex*16 :: basisGam(1:NbinGam, 1:NbasisGam)
  double precision :: RebasisGam(1:NbinGam, 1:NbasisGam)
  double precision :: ImbasisGam(1:NbinGam, 1:NbasisGam)
  complex*16 :: wbasis
  integer :: ibin, ibasis, dtp1, dtp2

  !============  Gamma in the fitting coeffecients =====================================
  ibin = get_bin_Gam(dt1, dt2)

  if(ibin==1) then
    do ibasis = 1, NBasisGam
      wbasis = (Beta/dble(MxT))**2.d0*weight_basis_Gam(CoefGam &
        & (0:BasisOrderGam,0:BasisOrderGam, ibasis,ibin), (real(dt1)+0.5d0)*Beta/MxT, &
        & (real(dt2)+0.5d0)*Beta/MxT)

      basisGam(ibin, ibasis) = basisGam(ibin, ibasis) + dcmplx(real(val), dimag(val))*wbasis

      RebasisGam(ibin, ibasis) = RebasisGam(ibin, ibasis) + real(val)**2.d0*wbasis
      ImbasisGam(ibin, ibasis) = ImbasisGam(ibin, ibasis) + dimag(val)**2.d0*wbasis
    enddo
  endif

  dtp1 = MxT-1 - dt1
  dtp2 = MxT-1 - dt2

  ibin = get_bin_Gam(dtp1, dtp2)

  if(ibin==1) then
    do ibasis = 1, NBasisGam
      wbasis = (Beta/dble(MxT))**2.d0*weight_basis_Gam(CoefGam &
        & (0:BasisOrderGam,0:BasisOrderGam, ibasis,ibin), (real(dtp1)+0.5d0)*Beta/MxT, &
        & (real(dtp2)+0.5d0)*Beta/MxT)

      basisGam(ibin, ibasis) = basisGam(ibin, ibasis) + dcmplx(real(val), -1.d0*dimag(val))*wbasis

      RebasisGam(ibin, ibasis) = RebasisGam(ibin, ibasis) + real(val)**2.d0*wbasis
      ImbasisGam(ibin, ibasis) = ImbasisGam(ibin, ibasis) + dimag(val)**2.d0*wbasis
    enddo
  endif

END SUBROUTINE accumulate_Gamma_basis



COMPLEX*16 FUNCTION Gam_from_Basis(it1, it2, basisGam)
  implicit none
  integer, intent(in) :: it1, it2
  complex*16, intent(in) :: basisGam(1:NbinGam, 1:NbasisGam)
  integer :: ibin, ibasis
  double precision :: tau1, tau2
  complex*16 :: cgam

  tau1 = (dble(it1)+0.5d0)*Beta/dble(MxT)
  tau2 = (dble(it2)+0.5d0)*Beta/dble(MxT)

  ibin = get_bin_Gam(it1, it2)

  if(ibin==1 .or. ibin==3) then
    cgam = (0.d0, 0.d0)
    do ibasis = 1, NBasisGam
      cgam = cgam + basisGam(1,ibasis)* weight_basis_Gam( &
        & CoefGam(0:BasisOrderGam,0:BasisOrderGam,ibasis,1), tau1, tau2)
    enddo
    Gam_from_Basis = cgam
  else if(ibin ==2) then
    tau1 = Beta - tau1
    tau2 = Beta - tau2

    cgam = (0.d0, 0.d0)
    do ibasis = 1, NBasisGam
      cgam = cgam + basisGam(1,ibasis)* weight_basis_Gam( &
        & CoefGam(0:BasisOrderGam,0:BasisOrderGam,ibasis,1), tau1, tau2)
    enddo
    Gam_from_Basis = dcmplx(real(cgam), -1.d0*dimag(cgam))
  endif

  return
END FUNCTION Gam_from_Basis

