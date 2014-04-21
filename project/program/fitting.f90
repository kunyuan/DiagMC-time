INCLUDE "mylib/mylib.f90"
INCLUDE "vrbls_mc.f90"
PROGRAM MAIN
  USE string_basic
  USE logging_module
  USE vrbls_mc


CONTAINS

!!============== GRAM-SCHMIDT BASIS ==================================
!!-------------- generate the Gram-Schmidt basis --------------
SUBROUTINE calculate_basis(nbasis, OmegaMin, OmegaMax, EP, EN)
  implicit none
  integer, intent(in):: nbasis, OmegaMin, OmegaMax
  double precision :: EP(nbasis, nbasis), EN(nbasis, nbasis)
  double precision :: VP(nbasis, nbasis), VN(nbasis, nbasis)
  double precision :: norm
  integer :: i, j, k

  EP(:,:) = 0.d0
  EN(:,:) = 0.d0
  VP(:,:) = 0.d0
  VN(:,:) = 0.d0

  do i = 1, nbasis
    VP(i, i) = 1.d0
    EP(i, :) = VP(i, :)

    do j = 1, i-1
      do k = 1, j
        EP(i, k) = EP(i, k) - EP(j, k) *projector(EP(j,:), VP(i, :), OmegaMin, OmegaMax)
      enddo
    enddo

    norm = dsqrt(projector(EP(i,:), EP(i,:), OmegaMin, OmegaMax))
    EP(i, :) = EP(i, :)/norm
  enddo

  do i = 1, nbasis
    VN(i, i) = 1.d0
    EN(i, :) = VN(i, :)

    do k = 1, j
      do j = 1, i-1
        EN(i, k) = EN(i, k) - EN(j, k) *projector(EN(j,:), VN(i, :), -OmegaMax, -OmegaMin)
      enddo
    enddo
    norm = dsqrt(projector(EN(i,:), EN(i,:), -OmegaMax, -OmegaMin))
    EN(i, :) = EN(i, :)/norm
  enddo

END SUBROUTINE calculate_basis



!---------- weight calculate f(omega) ------------------
DOUBLE PRECISION FUNCTION weight_basis(Coef, n)
  implicit none
  double precision :: Coef(nbasis)
  integer :: j, n
  weight_basis = 0.d0
  do j = 1, nbasis
    weight_basis = weight_basis + Coef(j)*OriginalBasis(j, n)
  enddo
  return 
END FUNCTION weight_basis


!---------- projector: \int_{omega} f1(omega)*f2(omega)  --------
DOUBLE PRECISION FUNCTION projector(Coef1, Coef2, OmegaMin, OmegaMax)
  implicit none
  double precision :: Coef1(nbasis), Coef2(nbasis)
  integer :: j, k, n, OmegaMin, OmegaMax

  projector  = 0.d0
  do j = 1, nbasis 
    do k = 1, nbasis
      do n = OmegaMin, OmegaMax
        projector = projector + Coef1(j)*OriginalBasis(j, n)*Coef2(k)*OriginalBasis(k, n)
      enddo
    enddo
  enddo
  return 
END FUNCTION projector




!--------- test subroutine for the Gram-Schmidt basis ------
SUBROUTINE test_basis(nbasis, OmegaMin, OmegaMax, CoefP, CoefN)
  implicit none
  integer, intent(in) :: nbasis, OmegaMin, OmegaMax
  double precision :: CoefP(nbasis, nbasis), CoefN(nbasis, nbasis)
  integer :: i, j, k, l, n
  double precision :: omega, x, y

  do i = 1, nbasis

    y = 0.d0
    do  n = OmegaMin, OmegaMax
      y = y + weight_basis(CoefP(i,:), n)**2.d0
    enddo
    if(dabs(y-1.d0)>1.d-10) then
      write(logstr, *) i, y, "Positive Basis error!!"
      call write_log
    endif

    y = 0.d0
    do  n = -OmegaMax, -OmegaMin
      y = y + weight_basis(CoefN(i,:), n)**2.d0
    enddo
    if(dabs(y-1.d0)>1.d-10) then
      write(logstr, *) i, y, "Negative Basis error!!"
      call write_log
    endif
  enddo

  do i = 1, nbasis
    do j = i+1, nbasis
      y = projector(CoefP(i, :), CoefP(j, :), OmegaMin, OmegaMax)
      if(dabs(y)>1.d-10) then
        write(logstr, *) i, j, y, "Positive Basis error!!"
        call write_log
      endif
      y = projector(CoefN(i, :), CoefN(j, :), -OmegaMax, -OmegaMin)
      if(dabs(y)>1.d-10) then
        write(logstr, *) i, j, y, "Negative Basis error!!"
        call write_log
      endif
    enddo
  enddo

END SUBROUTINE test_basis





!SUBROUTINE calculate_Gamma_basis(nbasisGamma, OmegaMin, OmegaMax, EPN, ENP, &
  !&  EPPR, ENNR, EPPL, ENNL)
  !implicit none
  !integer, intent(in):: nbasisGamma, OmegaMin, OmegaMax
  !double precision :: ENP(nbasisGamma, nbasisGamma),  EPN(nbasisGamma, nbasisGamma)
  !double precision :: EPPR(nbasisGamma, nbasisGamma), ENNR(nbasisGamma, nbasisGamma)
  !double precision :: EPPL(nbasisGamma, nbasisGamma), ENNL(nbasisGamma, nbasisGamma)
  !double precision :: VNP(nbasisGamma, nbasisGamma),  VPN(nbasisGamma, nbasisGamma)
  !double precision :: VPPR(nbasisGamma, nbasisGamma), VNNR(nbasisGamma, nbasisGamma)
  !double precision :: VPPL(nbasisGamma, nbasisGamma), VNNL(nbasisGamma, nbasisGamma)
  !double precision :: norm
  !integer :: i, j, k
  !VPN(:,:) = 0.d0
  !VNP(:,:) = 0.d0
  !VPPR(:,:) = 0.d0
  !VNNR(:,:) = 0.d0
  !VPPL(:,:) = 0.d0
  !VNNL(:,:) = 0.d0

  !EPN(:,:) = 0.d0
  !ENP(:,:) = 0.d0
  !EPPR(:,:) = 0.d0
  !ENNR(:,:) = 0.d0
  !EPPL(:,:) = 0.d0
  !ENNL(:,:) = 0.d0

  !do i = 1, nbasisGamma
    !VPN(i, i) = 1.d0
    !EPN(i, :) = VPN(i, :)

    !do j = 1, i-1
      !do k = 1, j
        !EPN(i, k) = EPN(i, k) - EPN(j, k) *projector_Gamma(EPN(j,:), VPN(i, :), OmegaMin, OmegaMax, &
          !& -OmegaMax, -OmegaMin)
      !enddo
    !enddo

    !norm = dsqrt(projector_Gamma(EPN(i, :), EPN(i,:), OmegaMin, OmegaMax, &
      !& -OmegaMax, -OmegaMin))
    !EPN(i, :) = EPN(i, :)/norm 
  !enddo

  !do i = 1, nbasisGamma
    !VNP(i, i) = 1.d0
    !ENP(i, :) = VNP(i, :)

    !do j = 1, i-1
      !do k = 1, j
        !ENP(i, k) = ENP(i, k) - ENP(j, k) *projector_Gamma(ENP(j,:), VNP(i, :), -OmegaMax, -OmegaMin, &
          !& OmegaMin, OmegaMax)
      !enddo
    !enddo

    !norm = dsqrt(projector_Gamma(ENP(i, :), ENP(i,:), -OmegaMax, -OmegaMin, &
      !& OmegaMin, OmegaMax))
    !ENP(i, :) =  ENP(i, :)/norm
  !enddo

  !do i = 1, nbasisGamma
    !VPPR(i, i) = 1.d0
    !EPPR(i, :) = VPPR(i, :)

    !do j = 1, i-1
      !do k = 1, j
        !EPPR(i, k) = EPPR(i, k) - EPPR(j, k) *projector_Gamma(EPPR(j,:), VPPR(i, :), OmegaMin, OmegaMax, &
          !& OmegaMin, 0)
      !enddo
    !enddo

    !norm = dsqrt(projector_Gamma(EPPR(i, :), EPPR(i,:), OmegaMin, OmegaMax, &
      !& OmegaMin, 0))
    !EPPR(i, :) = EPPR(i, :)/norm 
  !enddo

  !do i = 1, nbasisGamma
    !VPPL(i, i) = 1.d0
    !EPPL(i, :) = VPPL(i, :)

    !do j = 1, i-1
      !do k = 1, j
        !EPPL(i, k) = EPPL(i, k) - EPPL(j, k) *projector_Gamma(EPPL(j,:), VPPL(i, :), OmegaMin, OmegaMax, &
          !& 0, OmegaMax)
      !enddo
    !enddo

    !norm = dsqrt(projector_Gamma(EPPL(i, :), EPPL(i,:), OmegaMin, OmegaMax, &
      !& 0, OmegaMax))
    !EPPL(i, :) = EPPL(i, :)/norm 
  !enddo

  !do i = 1, nbasisGamma
    !VNNR(i, i) = 1.d0
    !ENNR(i, :) = VNNR(i, :)

    !do j = 1, i-1
      !do k = 1, j
        !ENNR(i, k) = ENNR(i, k) - ENNR(j, k) *projector_Gamma(ENNR(j,:), VNNR(i, :), -OmegaMax, -OmegaMin, &
          !& -OmegaMax, 0)
      !enddo
    !enddo

    !norm = dsqrt(projector_Gamma(ENNR(i, :), ENNR(i,:), -OmegaMax, -OmegaMin, &
      !& -OmegaMax, 0))
    !ENNR(i, :) = ENNR(i, :)/norm
  !enddo

  !do i = 1, nbasisGamma
    !VNNL(i, i) = 1.d0
    !ENNL(i, :) = VNNL(i, :)

    !do j = 1, i-1
      !do k = 1, j
        !ENNL(i, k) = ENNL(i, k) - ENNL(j, k) *projector_Gamma(ENNL(j,:), VNNL(i, :), -OmegaMax, -OmegaMin, &
          !& 0, -OmegaMin)
      !enddo
    !enddo

    !norm = dsqrt(projector_Gamma(ENNL(i, :), ENNL(i,:), -OmegaMax, -OmegaMin, &
      !& 0, -OmegaMin))
    !ENNL(i, :) = ENNL(i, :)/norm
  !enddo
  !return
!END SUBROUTINE calculate_Gamma_basis




!!---------- weight calculate f(omega) ------------------
!DOUBLE PRECISION FUNCTION weight_basis_Gamma(Coef, n1, n2)
  !implicit none
  !double precision :: Coef(nbasisGamma)
  !integer :: j, n1, n2
  !weight_basis_Gamma = 0.d0
  !do j = 1, nbasisGamma
    !weight_basis_Gamma = weight_basis_Gamma+Coef(j)*OriginalBasisGamma(j, n1, n2)
  !enddo
  !return 
!END FUNCTION weight_basis_Gamma


!!---------- projector: \int_{omega1, omega2} f1(omega1,omega2)*f2(omega1,omega2)  --------
!DOUBLE PRECISION FUNCTION projector_Gamma(Coef1, Coef2, OmegaMin1, OmegaMax1, OmegaMin2, OmegaMax2)
  !implicit none
  !double precision :: Coef1(nbasisGamma), Coef2(nbasisGamma)
  !integer :: j, k, n1, n2, OmegaMin1, OmegaMax1, OmegaMin2, OmegaMax2

  !projector_Gamma  = 0.d0

  !do j = 1, nbasisGamma
    !do k = 1, nbasisGamma
      !do n1 = OmegaMin1, OmegaMax1
        !if(OmegaMin2==0) then
          !if(n1<OmegaMax1) then
            !do n2 = n1+1, OmegaMax1
              !projector_Gamma = projector_Gamma + Coef1(j)*OriginalBasisGamma(j, n1, n2)*Coef2(k)* &
                !& OriginalBasisGamma(k, n1, n2)*ReweightBasis(n1, n2)
            !enddo
          !endif
        !else if(OmegaMax2==0) then
          !if(n1>OmegaMin1) then
            !do n2 = OmegaMin2, n1-1
              !projector_Gamma = projector_Gamma + Coef1(j)*OriginalBasisGamma(j, n1, n2)*Coef2(k)* &
                !& OriginalBasisGamma(k, n1, n2)*ReweightBasis(n1, n2)
            !enddo
          !endif
        !else
          !do n2 = OmegaMin2, OmegaMax2
            !projector_Gamma = projector_Gamma + Coef1(j)*OriginalBasisGamma(j, n1, n2)*Coef2(k)* &
              !& OriginalBasisGamma(k, n1, n2)*ReweightBasis(n1, n2)
          !enddo
        !endif
      !enddo
    !enddo
  !enddo
  !return 
!END FUNCTION projector_Gamma


!!--------- test subroutine for the Gram-Schmidt basis ------
!SUBROUTINE test_basis_Gamma(nbasisGamma, OmegaMin, OmegaMax, EPN, ENP, &
    !& EPPR, ENNR, EPPL, ENNL)
  !implicit none
  !integer, intent(in) :: nbasisGamma, OmegaMin, OmegaMax
  !double precision :: ENP(nbasisGamma, nbasisGamma), EPN(nbasisGamma,nbasisGamma)
  !double precision :: EPPR(nbasisGamma, nbasisGamma),ENNR(nbasisGamma,nbasisGamma)
  !double precision :: EPPL(nbasisGamma, nbasisGamma),ENNL(nbasisGamma,nbasisGamma)
  !integer :: i, j, k, l, n1, n2
  !double precision :: omega, x, y

  !do i = 1, nbasisGamma
    !y = 0.d0
    !do  n1 = OmegaMin+1, OmegaMax
      !do n2 = OmegaMin, n1-1
        !y = y + weight_basis_Gamma(EPPR(i,:), n1, n2)**2.d0*ReweightBasis(n1, n2)
      !enddo
    !enddo
    !if(dabs(y-1.d0)>1.d-8) then
      !write(logstr, *) i, y-1.d0, "++R Gamma Basis error!!"
      !call write_log
    !endif

    !y = 0.d0
    !do  n1 = OmegaMin, OmegaMax-1
      !do n2 = n1+1, OmegaMax
        !y = y + weight_basis_Gamma(EPPL(i,:), n1, n2)**2.d0*ReweightBasis(n1, n2)
      !enddo
    !enddo
    !if(dabs(y-1.d0)>1.d-8) then
      !write(logstr, *) i, y-1.d0, "++L Gamma Basis error!!"
      !call write_log
    !endif

    !y = 0.d0
    !do  n1 = OmegaMin, OmegaMax
      !do n2 = -OmegaMax, -OmegaMin
        !y = y + weight_basis_Gamma(EPN(i,:), n1, n2)**2.d0*ReweightBasis(n1, n2)
      !enddo
    !enddo
    !if(dabs(y-1.d0)>1.d-8) then
      !write(logstr, *) i, y-1.d0, "+- Gamma Basis error!!"
      !call write_log
    !endif

    !y = 0.d0
    !do  n1 = -OmegaMax, -OmegaMin
      !do n2 = OmegaMin, OmegaMax
        !y = y + weight_basis_Gamma(ENP(i,:), n1, n2)**2.d0*ReweightBasis(n1, n2)
      !enddo
    !enddo
    !if(dabs(y-1.d0)>1.d-8) then
      !write(logstr, *) i, y-1.d0, "-+ Gamma Basis error!!"
      !call write_log
    !endif

    !y = 0.d0
    !do  n1 = -OmegaMax+1, -OmegaMin
      !do n2 = -OmegaMax, n1-1
        !y = y + weight_basis_Gamma(ENNR(i,:), n1, n2)**2.d0*ReweightBasis(n1, n2)
      !enddo
    !enddo
    !if(dabs(y-1.d0)>1.d-8) then
      !write(logstr, *) i, y-1.d0, "--R Gamma Basis error!!"
      !call write_log
    !endif

    !y = 0.d0
    !do  n1 = -OmegaMax, -OmegaMin-1
      !do n2 = n1+1, -OmegaMin
        !y = y + weight_basis_Gamma(ENNL(i,:), n1, n2)**2.d0*ReweightBasis(n1, n2)
      !enddo
    !enddo
    !if(dabs(y-1.d0)>1.d-8) then
      !write(logstr, *) i, y-1.d0, "--L Gamma Basis error!!"
      !call write_log
    !endif
  !enddo


  !do i = 1, nbasisGamma
    !do j = i+1, nbasisGamma
      !y = projector_Gamma(EPN(i, :), EPN(j, :), OmegaMin, OmegaMax,-OmegaMax,-OmegaMin)
      !if(dabs(y)>1.d-8) then
        !write(logstr, *) i, j, y, "+- Gamma Basis error!!"
        !call write_log
      !endif
      !y = projector_Gamma(ENP(i, :), ENP(j, :),-OmegaMax,-OmegaMin, OmegaMin, OmegaMax)
      !if(dabs(y)>1.d-8) then
        !write(logstr, *) i, j, y, "-+ Gamma Basis error!!"
        !call write_log
      !endif
      !y = projector_Gamma(EPPR(i, :), EPPR(j, :), OmegaMin, OmegaMax, OmegaMin, 0)
      !if(dabs(y)>1.d-7) then
        !write(logstr, *) i, j, y, "++R Gamma Basis error!!"
        !call write_log
      !endif
      !y = projector_Gamma(EPPL(i, :), EPPL(j, :), OmegaMin, OmegaMax, 0, OmegaMax)
      !if(dabs(y)>1.d-7) then
        !write(logstr, *) i, j, y, "++L Gamma Basis error!!"
        !call write_log
      !endif
      !y = projector_Gamma(ENNR(i, :), ENNR(j, :),-OmegaMax,-OmegaMin,-OmegaMax, 0)
      !if(dabs(y)>1.d-7) then
        !write(logstr, *) i, j, y, "--R Gamma Basis error!!"
        !call write_log
      !endif
      !y = projector_Gamma(ENNL(i, :), ENNL(j, :),-OmegaMax,-OmegaMin, 0, -OmegaMin)
      !if(dabs(y)>1.d-7) then
        !write(logstr, *) i, j, y, "--L Gamma Basis error!!"
        !call write_log
      !endif
    !enddo
  !enddo

!END SUBROUTINE test_basis_Gamma
END PROGRAM MAIN
