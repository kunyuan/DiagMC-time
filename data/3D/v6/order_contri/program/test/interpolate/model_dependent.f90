
!!--------- extract weight for G ---------
COMPLEX*16 FUNCTION weight_G(typ1, tau1)
  implicit none
  double precision, intent(in)  :: tau1
  integer, intent(in) :: typ1
  call interpolation_GW(tau1, G(typ1, :), weight_G)
  if(tau1<=1.d-14)  weight_G = -1.d0* weight_G
END FUNCTION weight_G

!!--------- extract weight for W ---------
COMPLEX*16 FUNCTION weight_W(typ1, site, tau1)
  implicit none
  integer, intent(in)  :: site, typ1
  double precision, intent(in) :: tau1
  call interpolation_GW(tau1, W(typ1, site, :), weight_W)
END FUNCTION weight_W

!!--------- extract weight for Gamma ---------
COMPLEX*16 FUNCTION weight_Gam(typ1, site, tau1, tau2)
  implicit none
  integer, intent(in)  :: site, typ1
  double precision, intent(in) :: tau1, tau2
  call interpolation_Gam(tau1, tau2, Gam(typ1,site,:,:), weight_Gam)
  if(tau1<=1.d-14 .and. tau2>1.d-14) then
    weight_Gam = -1.d0 *weight_Gam
  else if(tau1>1.d-14 .and. tau2<=1.d-14) then
    weight_Gam = -1.d0 *weight_Gam
  endif
END FUNCTION weight_Gam



SUBROUTINE interpolation_GW(tau, Table, Val)
  implicit none
  double precision, intent(in) :: tau
  complex*16, dimension(0:MxT-1), intent(in) :: Table
  complex*16, intent(out) :: Val
  integer :: ts, tb
  double precision :: ptau, taus, taub

  ptau = shift_tau(tau)

  ts = Floor((ptau-0.5d0*Beta/MxT)*MxT/Beta)
  tb = ts + 1

  if(ts<0) then
    ts=0
    tb=1
  else if(tb>=MxT) then
    ts=MxT-2
    tb=MxT-1
  endif

  taus = (real(ts)+0.5d0)*Beta/MxT
  taub = (real(tb)+0.5d0)*Beta/MxT

  Val = interpolate(ptau, taus, Table(ts), taub, Table(tb))
  return
END SUBROUTINE interpolation_GW



SUBROUTINE interpolation_Gam(tau1, tau2, Table, Val)
  implicit none
  double precision, intent(in) :: tau1, tau2
  complex*16, dimension(0:MxT-1, 0:MxT-1), intent(in) :: Table
  complex*16, intent(out) :: Val
  integer :: i, ts(1:2), tb(1:2)
  double precision :: ptau(1:2), taus(1:2), taub(1:2)
  complex*16 :: Val1, Val2

  ptau(1) = shift_tau(tau1)
  ptau(2) = shift_tau(tau2)

  do i = 1, 2
    ts(i) = Floor(ptau(i)*MxT/Beta)
    tb(i) = ts(i) + 1
    if(ts(i)<0 .or. ts(i)>=MxT) call LogFile%QuickLog("ts("+str(i)+") error!"+str(ts(i)),'e')

    taus(i) = real(ts(i))*Beta/MxT
    taub(i) = real(tb(i))*Beta/MxT
  enddo

  Val1 = interpolate(ptau(1), taus(1), Table(ts(1), ts(2)), taub(1), Table(tb(1), ts(2)))
  Val2 = interpolate(ptau(1), taus(1), Table(ts(1), tb(2)), taub(1), Table(tb(1), tb(2)))

  Val  = interpolate(ptau(2), taus(2), Val1, taub(2), Val2)
  return
END SUBROUTINE interpolation_Gam

complex*16 FUNCTION interpolate(tau, taus, Fs, taub, Fb)
  double precision, intent(in) :: tau
  double precision, intent(in) :: taus, taub
  complex*16, intent(in) :: Fs, Fb
  interpolate = ((tau-taus)*Fb+(taub-tau)*Fs)*MxT/Beta
  return
END FUNCTION interpolate

DOUBLE PRECISION FUNCTION shift_tau(tau)
  implicit none
  double precision, intent(in) :: tau
  shift_tau = tau
  if(tau<=1.d-14) then
    if(abs(tau)<=1.d-14) then
      shift_tau = 0.d0
    else
      shift_tau = tau + Beta
    endif
  endif
  return
END FUNCTION shift_tau
