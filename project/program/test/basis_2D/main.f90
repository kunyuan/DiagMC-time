INCLUDE "mylib/mylib.f90"
INCLUDE "vrbls_mc.f90"
PROGRAM MAIN
  USE string_basic
  USE logging_module
  USE vrbls_mc
  implicit none
  complex*16 :: Gamma(0:MxT-1, 0:MxT-1)
  double precision :: ReGamma(0:MxT-1, 0:MxT-1)
  double precision :: ImGamma(0:MxT-1, 0:MxT-1)

  complex*16 :: GamDiag(0:MxT-1)
  double precision :: ReGamDiag(0:MxT-1)
  double precision :: ImGamDiag(0:MxT-1)

  complex*16 :: GammaBasis(1:NbinGam, 1:NBasisGam)
  double precision :: ReGammaBasis(1:NbinGam, 1:NBasisGam)
  double precision :: ImGammaBasis(1:NbinGam, 1:NBasisGam)
  complex*16 :: gam1
  integer :: it1, it2


  Beta = 1.d0


  Gamma = (0.d0, 0.d0)
  ReGamma = 0.d0
  ImGamma = 0.d0

  GamDiag = (0.d0, 0.d0)
  ReGamDiag = 0.d0
  ImGamDiag = 0.d0

  GammaBasis = (0.d0, 0.d0)
  ReGammaBasis = 0.d0
  ImGammaBasis = 0.d0

  call initialize_polynomials
  call initialize_bins
  call calculate_basis_GWGam

  call initial_Gamma(Gamma)

  open(10, file='Gamma.dat')
  write(10, *) "##########Gamma"
  write(10, *) "#tau1:", MxT, " tau2:", MxT
  do it1 = 0, MxT-1
    do it2 = 0, MxT-1
      gam1 = Gamma(it1, it2)
      write(10, *) real(gam1), dimag(gam1)
    enddo
  enddo
  close(10)

  do it1 = 0, MxT-1
    !call accumulate_diagonal_Gamma(it1, Gamma(it1, MxT-1-it1)/2.d0, GamDiag(:), &
      !& ReGamDiag(:), ImGamDiag(:))
    do it2 = 0, MxT-1
      call accumulate_Gamma_basis(it1, it2, Gamma(it1,it2)/2.d0, GammaBasis(:,:), &
        & ReGammaBasis(:,:), ImGammaBasis(:,:))
    enddo
  enddo

  open(11, file='Gamma_basis.dat')
  write(11, *) "##########GammaBasis"
  write(11, *) "#tau1:", MxT, " tau2:", MxT
  do it1 = 0, MxT-1
    do it2 = 0, MxT-1
      gam1 = Gam_from_Basis(it1, it2, GammaBasis)
      write(11, *) real(gam1), dimag(gam1)
    enddo
  enddo
  close(11)

  !open(12, file='Gamma_diag.dat')
  !do it1 = 0, MxT-1
    !gam1 = Gamma(it1, MxT-1-it1)
    !write(12, *) 1, it1, real(gam1), dimag(gam1)
    !gam1 = GamDiag(it1)
    !write(12, *) 2, it1, real(gam1), dimag(gam1)
  !enddo
  !close(12)


CONTAINS
INCLUDE "basis.f90"
INCLUDE "fitting.f90"

SUBROUTINE initial_Gamma(Gam)
  implicit none
  integer :: it1, it2
  double precision :: tau1, tau2
  complex*16 :: Gam(0:MxT-1, 0:MxT-1)
  do it1 = 0, MxT-1
    do it2 = 0, MxT-1
      tau1 = (real(it1)+0.5d0)*Beta/MxT
      tau2 = (real(it2)+0.5d0)*Beta/MxT
      Gam(it1, it2) = (tau1-Beta/2.d0)**2.d0+(tau2-Beta/2.d0)**2.d0
    enddo
  enddo
  return
END SUBROUTINE initial_Gamma

END PROGRAM
