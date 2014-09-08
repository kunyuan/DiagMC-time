
PROGRAM MAIN
  implicit none
  integer :: i, dt1, dt2, N
  double precision :: dtau1, dtau2, tau1, tau2
  integer, parameter :: MxT=100
  double precision, parameter :: Beta=1.d0

  N = 200
  do i = -N, N-1
    tau1 = i*Beta/N+1.d-7
    dtau1 = shift_tau(tau1)
    dt1 = discrete_tau(dtau1)
    print *, tau1, dtau1, dt1
  enddo

CONTAINS

DOUBLE PRECISION FUNCTION shift_tau(tau)
  implicit none
  double precision, intent(in) :: tau

  shift_tau = tau
  if(tau<=1.d-14) then
    if(abs(tau)<=1.d-14)  shift_tau = 0.d0
    shift_tau = shift_tau + Beta
  endif
  return
END FUNCTION shift_tau


INTEGER FUNCTION discrete_tau(tau)
  implicit none
  double precision, intent(in) :: tau
  discrete_tau = Floor(tau*MxT/Beta)
  if(discrete_tau==MxT)  discrete_tau = MxT-1
  return
END FUNCTION discrete_tau

END PROGRAM
