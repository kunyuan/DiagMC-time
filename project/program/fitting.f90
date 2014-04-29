!include "mylib/mylib.f90"
!INCLUDE "vrbls_mc.f90"
!PROGRAM MAIN
  !USE string_basic
  !USE logging_module
  !USE vrbls_mc
  !implicit none
  !integer :: i, j, t1, t2, t
  !double precision :: tau, tau1, tau2
  !complex*16 :: muc, FG(0:MxT-1)
  !complex*16 :: GBasis(1:NbinG, 1:NBasis)
  !complex*16 :: Gfit
  !complex*16 :: FGamma(0:MxT-1, 0:MxT-1)
  !complex*16 :: GamBasis(1:NbinGam, 1:NBasisGam)
  !complex*16 :: Gamfit

  !Beta = 1.d0

  !call initialize_polynomials
  !call initialize_bins
  !call calculate_basis_GWGam 

  !!=========fack G ==============================
  !muc = dcmplx(0.d0, pi/(2.d0*Beta))
  !do t = 0, MxT-1
    !tau = real(t)*Beta/MxT
    !FG(t) = cdexp(muc*tau)/(1.d0, 1.d0) 
  !enddo

  !!======================== calculate the fitting coeffcients for G ==========================
  !call LogFile%QuickLog("test for G fitting ...")
  !do j = 1, NBasis
    !do i = 1, NbinG
      !GBasis(i, j) = (0.d0, 0.d0)
      !do t = FromG(i), ToG(i)
        !GBasis(i, j) = GBasis(i, j) + (Beta/dble(MxT))*FG(t)*weight_basis(CoefG(:, j, i), t)
      !enddo
    !enddo
  !enddo

  !!========================= test the fitting function for G ================================
  !do i = 1, NbinG
    !do t = FromG(i), ToG(i)
      !Gfit = (0.d0, 0.d0)
      !do j = 1, NBasis
        !Gfit = Gfit + weight_basis(CoefG(:, j, i), t)*GBasis(i, j)
      !enddo
      !if(abs(Gfit-FG(t))>1.d-6)  then
        !call LogFile%QuickLog(str(t)+str(Gfit-FG(t))+"G fitting error!")
      !endif
    !enddo
  !enddo
  !call LogFile%QuickLog("test for G fitting done~")



  !!=========fack Gamma ==============================
  !do t1 = 0, MxT-1
    !do t2 = 0, MxT-1
      !tau1 = real(t1)*Beta/MxT
      !tau2 = real(t2)*Beta/MxT
      !FGamma(t1, t2) = dcmplx((tau1+tau2)**2.d0, 0.d0)
    !enddo
  !enddo

  !!======================== calculate the fitting coeffcients for Gamma ==========================
  !do i = 1, NbinGam
    !if(IsBasis2D(i)) then
      !do j = 1, NBasisGam
        !GamBasis(i, j) = (0.d0, 0.d0)
        !do t1 = FromGamT1(i), ToGamT1(i)
          !do t2 = FromGamT2(t1, i), ToGamT2(t1, i)
            !GamBasis(i, j) = GamBasis(i, j) + (Beta/dble(MxT))**2.d0*FGamma(t1,t2)* &
              !& weight_basis_Gam(CoefGam(0:BasisOrderGam, 0:BasisOrderGam, j, i), t1, t2)
          !enddo
        !enddo
      !enddo
    !else
      !do j = 1, NBasis
        !GamBasis(i, j) = (0.d0, 0.d0)
        !do t1 = FromGamT1(i), ToGamT1(i)
          !GamBasis(i, j) = GamBasis(i, j) + (Beta/dble(MxT))*FGamma(t1, FromGamT2(t1, i))* &
            !& weight_basis(CoefGam(0:BasisOrder, 0, j, i), t1)
        !enddo
      !enddo
    !endif
  !enddo

  !!========================= test the fitting function for Gamma ================================
  !call LogFile%QuickLog("test for Gamma fitting ...")
  !do i = 1, NbinGam
    !if(IsBasis2D(i)) then
      !do t1 = FromGamT1(i), ToGamT1(i)
        !do t2 = FromGamT2(t1, i), ToGamT2(t1, i)
          !Gamfit = (0.d0, 0.d0)
          !do j = 1, NBasisGam
            !Gamfit = Gamfit + GamBasis(i, j) * weight_basis_Gam(CoefGam(0:BasisOrderGam, &
              !& 0:BasisOrderGam, j, i), t1, t2)
          !enddo
          !if(abs(Gamfit-FGamma(t1,t2))>1.d-9)  then
            !call LogFile%QuickLog(str(t1)+str(t2)+str(Gamfit-FGamma(t1,t2))+"Gamma fitting error!")
          !endif
        !enddo
      !enddo
    !else 
      !do t1 = FromGamT1(i), ToGamT1(i)
        !Gamfit = (0.d0, 0.d0)
        !do j = 1, NBasis
          !Gamfit = Gamfit + GamBasis(i, j) * weight_basis(CoefGam(0:BasisOrder, 0, j, i), t1)
        !enddo
        !t2 = FromGamT2(t1, i)
        !if(abs(Gamfit-FGamma(t1,t2))>1.d-9)  then
          !call LogFile%QuickLog(str(t1)+str(t2)+str(Gamfit-FGamma(t1,t2))+"Gamma fitting error!")
        !endif
      !enddo
    !endif
  !enddo
  !call LogFile%QuickLog("test for Gamma fitting done~")

!CONTAINS




!============== divide the space into N bins for G,W, Gamma =============================
SUBROUTINE initialize_bins
  implicit none
  integer :: t1, t2

  !================== define the bins for G,W, Gamma ====================================
  FromG(1) = 0
  ToG(1)   = MxT-1

  FromW(1) = 0
  ToW(1)   = MxT-1

  IsBasis2D(1) = .true.
  FromGamT1(1) = 0
  ToGamT1(1) = MxT-2
  do t1 = 0, MxT-2
    FromGamT2(t1, 1) = 0
    ToGamT2(t1, 1)   = MxT-2-t1
  enddo

  IsBasis2D(2) = .false.
  FromGamT1(2) = 0
  ToGamT1(2) = MxT-1
  do t1 = 0, MxT-1
    FromGamT2(t1, 2) = MxT-1-t1
    ToGamT2(t1, 2) = MxT-1-t1
  enddo

  IsBasis2D(3) = .true.
  FromGamT1(3) = 1
  ToGamT1(3) = MxT-1
  do t1 = 1, MxT-1
    FromGamT2(t1, 3) = MxT-t1
    ToGamT2(t1, 3)   = MxT-1
  enddo

  !!=============== test the get_bin_Gam ======================
  !do t1 = 0, MxT-1
    !t2 = MxT-1-t1-2
    !if(t2>=0) then
      !if(get_bin_Gam(t1, t2)/=1) print *, "bin 1 error!"
    !endif

    !t2 = MxT-1-t1
    !if(get_bin_Gam(t1, t2)/=2) then
      !print *, "bin 2 error!"
    !endif

    !t2 = MxT-1-t1+2
    !if(t2<MxT) then
      !if(get_bin_Gam(t1, t2)/=3) print *, "bin 3 error!"
    !endif
  !enddo
  return
END SUBROUTINE initialize_bins


INTEGER FUNCTION get_bin_Gam(it1, it2)
  implicit none
  integer :: it1, it2
  if(it1+it2<MxT-1)  get_bin_Gam = 1
  if(it1+it2==MxT-1) get_bin_Gam = 2
  if(it1+it2>MxT-1)  get_bin_Gam = 3
  return
END FUNCTION get_bin_Gam

INTEGER FUNCTION get_bin_W(it1)
  implicit none
  integer :: it1
  get_bin_W = 1
  return
END FUNCTION get_bin_W

INTEGER FUNCTION get_bin_G(it1)
  implicit none
  integer :: it1
  get_bin_G = 1
  return
END FUNCTION get_bin_G

!============== set the polynomial vectors ==============================================
!---------- you can set any polynomial functions as you want -----------------
SUBROUTINE initialize_polynomials
  implicit none
  integer :: i, it1, it2, ibasis

  Polynomial(:, :) = 0.d0
  do i = 1, BasisOrder
    Polynomial(i, i) = 1.d0
  enddo
  Polynomial(0, NBasis) = 1.d0

  PolynomialGam(:, :, :) = 0.d0
  ibasis = 0
  do it1 = 0, BasisOrderGam
    do it2 = 0, BasisOrderGam
      if(it1==0 .and. it2==0) cycle
      ibasis = ibasis + 1
      PolynomialGam(it1, it2, ibasis) = 1.d0
    enddo
  enddo

  PolynomialGam(0, 0, NBasisGam) = 1.d0

END SUBROUTINE initialize_polynomials


DOUBLE PRECISION FUNCTION get_polynomial(iorder, tau)
  implicit none
  integer, intent(in) :: iorder
  double precision, intent(in) :: tau
  get_polynomial = (tau)**iorder
  return
END FUNCTION get_polynomial

DOUBLE PRECISION FUNCTION get_polynomial_Gam(iorder1, iorder2, tau1, tau2)
  implicit none
  integer, intent(in) :: iorder1, iorder2
  double precision, intent(in) :: tau1, tau2
  get_polynomial_Gam = (tau1)**iorder1 *(tau2)**iorder2
  return
END FUNCTION get_polynomial_Gam


SUBROUTINE calculate_basis_GWGam
  implicit none
  integer :: i

  do i = 1, NbinG
    call calculate_basis(FromG(i), ToG(i), CoefG(:, :, i))
  enddo

  do i = 1, NbinW
    call calculate_basis(FromW(i), ToW(i), CoefW(:, :, i))
  enddo

  do i = 1, NbinGam
    if(IsBasis2D(i)) then
      call calculate_basis_Gamma_2D(FromGamT1(i), ToGamT1(i), FromGamT2(:,i), ToGamT2(:,i), &
        & CoefGam(0:BasisOrderGam, 0:BasisOrderGam, 1:NBasisGam, i))
    else
      call calculate_basis(FromGamT1(i), ToGamT1(i), CoefGam(0:BasisOrder, &
        & 0, 1:NBasis, i))
    endif
  enddo
  return
END SUBROUTINE calculate_basis_GWGam


!!============== GRAM-SCHMIDT BASIS ==================================

!============== basis for G and W ====================================
!!-------------- generate the Gram-Schmidt basis --------------
SUBROUTINE calculate_basis(tmin, tmax, Coef)
  implicit none
  integer, intent(in):: tmin, tmax
  double precision, dimension(0:BasisOrder, NBasis) :: Coef
  double precision :: norm, projtmp
  integer :: i, j, k

  Coef(:, :) = 0.d0

  !========== projecting and normalizing the basis ===========================
  do i = 1, NBasis

    Coef(:, i) = Polynomial(:, i)

    do j = 1, i-1
      projtmp = projector(tmin, tmax, Coef(:,j), Polynomial(:,i))
      Coef(:, i) = Coef(:, i) - Coef(:, j) *projtmp
    enddo

    norm = dsqrt(projector(tmin, tmax, Coef(:, i), Coef(:, i)))
    Coef(:, i) = Coef(:, i)/norm
  enddo

  call test_basis(tmin, tmax, Coef)

END SUBROUTINE calculate_basis



!=============================================================================
!=================== basis for Gamma==========================================
!=============================================================================



SUBROUTINE calculate_basis_Gamma_2D(t1min, t1max, t2min, t2max, Coef)
  implicit none
  integer, intent(in):: t1min, t1max
  integer, intent(in):: t2min(0:MxT-1), t2max(0:MxT-1)
  double precision, dimension(0:BasisOrderGam, 0:BasisOrderGam, NBasisGam) :: Coef
  double precision :: norm, projtmp
  integer :: i, j, k

  Coef(:, :, :) = 0.d0

  !========== projecting and normalizing the basis ===========================
  do i = 1, NBasisGam

    Coef(:, :, i) = PolynomialGam(:, :, i)

    do j = 1, i-1
      projtmp = projector_Gam(t1min, t1max, t2min, t2max, Coef(:,:,j), PolynomialGam(:,:,i))
      Coef(:, :, i) = Coef(:,:,i) - Coef(:,:,j) *projtmp
    enddo

    norm = dsqrt(projector_Gam(t1min, t1max, t2min, t2max, Coef(:, :, i), Coef(:,:,i)))
    Coef(:, :, i) = Coef(:, :, i)/norm
  enddo

  call test_basis_Gam(t1min, t1max, t2min, t2max, Coef)

END SUBROUTINE calculate_basis_Gamma_2D


DOUBLE PRECISION FUNCTION projector_Gam(t1min, t1max, t2min, t2max, Coef1, Coef2)
  implicit none
  integer, intent(in) :: t1min, t1max
  integer, intent(in) :: t2min(0:MxT-1), t2max(0:MxT-1)
  double precision, intent(in) :: Coef1(0:BasisOrderGam, 0:BasisOrderGam), Coef2(0:BasisOrderGam, 0:BasisOrderGam)
  integer :: it1, jt1, t1, it2, jt2, t2
  double precision :: tau1, tau2

  projector_Gam  = 0.d0
  do t1 = t1min, t1max
    do t2 = t2min(t1), t2max(t1)
      tau1 = dble(t1)*Beta/dble(MxT)
      tau2 = dble(t2)*Beta/dble(MxT)
      do it1 = 0, BasisOrderGam 
        do jt1 = 0, BasisOrderGam
          do it2 = 0, BasisOrderGam 
            do jt2 = 0, BasisOrderGam
              projector_Gam = projector_Gam + (Beta/dble(MxT))**2.d0*Coef1(it1,it2)*Coef2(jt1, &
                &  jt2)*tau1**(dble(it1+jt1))*tau2**(dble(it2+jt2))
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
  return 
END FUNCTION projector_Gam

DOUBLE PRECISION FUNCTION weight_basis_Gam(Coef, t1, t2)
  implicit none
  double precision :: Coef(0:BasisOrderGam, 0:BasisOrderGam)
  integer :: j1, t1, j2, t2
  double precision :: tau1, tau2
  tau1 = dble(t1)*Beta/dble(MxT)
  tau2 = dble(t2)*Beta/dble(MxT)

  weight_basis_Gam = 0.d0
  do j1 = 0, BasisOrderGam
    do j2 = 0, BasisOrderGam
      weight_basis_Gam = weight_basis_Gam + Coef(j1, j2)*tau1**dble(j1)*tau2**dble(j2)
    enddo
  enddo
  return 
END FUNCTION weight_basis_Gam


!--------- test subroutine for the Gram-Schmidt basis ------
SUBROUTINE test_basis_Gam(t1min, t1max, t2min, t2max, Coef)
  implicit none
  integer, intent(in) :: t1min, t1max
  integer, intent(in) :: t2min(0:MxT-1), t2max(0:MxT-1)
  double precision, intent(in) :: Coef(0:BasisOrderGam, 0:BasisOrderGam, NbasisGam)
  integer :: i, j, k, l, n, t1, t2
  double precision :: omega, x, y

  do i = 1, NBasisGam
    y = 0.d0
    do  t1 = t1min, t1max
      do  t2 = t2min(t1), t2max(t1)
        y = y + (Beta/dble(MxT))**2.d0*weight_basis_Gam(Coef(:, :, i), t1, t2)**2.d0
      enddo
    enddo
    if(dabs(y-1.d0)>1.d-8) then
      call LogFile%QuickLog(str(i)+str(y)+"Gamma Basis error!! Not normal!")
    endif
  enddo

  do i = 1, NBasisGam
    do j = i+1, NBasisGam
      y = projector_Gam(t1min, t1max, t2min, t2max, Coef(:, :, i), Coef(:, :, j))
      if(dabs(y)>5.d-7) then
        call LogFile%QuickLog(str(i)+str(j)+str(y)+"Gamma Basis error!! Not orthogonal!")
      endif
    enddo
  enddo

END SUBROUTINE test_basis_Gam

<<<<<<< HEAD
DOUBLE PRECISION FUNCTION projector(tmin, tmax, Coef1, Coef2)
  implicit none
  integer, intent(in) :: tmin, tmax
  double precision, intent(in) :: Coef1(0:BasisOrder), Coef2(0:BasisOrder)
  integer :: j, k, t
  double precision :: tau

  projector  = 0.d0
  do t = tmin, tmax
    tau = dble(t)*Beta/dble(MxT)
    do j = 0, BasisOrder 
      do k = 0, BasisOrder
        projector = projector + Beta/dble(MxT)*Coef1(j)*Coef2(k)*tau**(dble(j+k))
      enddo
    enddo
  enddo
  return 
END FUNCTION projector

DOUBLE PRECISION FUNCTION weight_basis(Coef, tau)
  implicit none
  double precision :: Coef(0:BasisOrder)
  integer :: j
  double precision :: tau

  weight_basis = 0.d0
  do j = 0, BasisOrder
    weight_basis = weight_basis + Coef(j)*tau**dble(j)
  enddo
  return 
END FUNCTION weight_basis


!--------- test subroutine for the Gram-Schmidt basis ------
SUBROUTINE test_basis(tmin, tmax, Coef)
  implicit none
  integer, intent(in) :: tmin, tmax
  double precision, intent(in) :: Coef(0:BasisOrder, Nbasis)
  integer :: i, j, k, l, n, t
  double precision :: x, y
  double precision :: tau, tau1, tau2

  do i = 1, NBasis
    y = 0.d0
    do  t = tmin, tmax
      tau = dble(t)*Beta/dble(MxT)
      y = y + Beta/dble(MxT)*weight_basis(Coef(:, i), tau)**2.d0
    enddo
    if(dabs(y-1.d0)>1.d-8) then
      call LogFile%QuickLog(str(i)+str(y)+"Basis error!! Not normal!")
    endif
  enddo

  do i = 1, Nbasis
    do j = i+1, Nbasis
      y = projector(tmin, tmax, Coef(:, i), Coef(:, j))
      if(dabs(y)>1.d-8) then
        call LogFile%QuickLog(str(i)+str(j)+str(y)+"Basis error!! Not orthogonal!")
      endif
    enddo
  enddo

END SUBROUTINE test_basis




!=============================================================================
!=================== basis for Gamma==========================================
!=============================================================================



SUBROUTINE calculate_basis_Gamma_2D(t1min, t1max, t2min, t2max, Coef)
  implicit none
  integer, intent(in):: t1min, t1max
  integer, intent(in):: t2min(0:MxT-1), t2max(0:MxT-1)
  double precision, dimension(0:BasisOrderGam, 0:BasisOrderGam, NBasisGam) :: Coef
  double precision :: norm, projtmp
  integer :: i, j, k

  Coef(:, :, :) = 0.d0

  !========== projecting and normalizing the basis ===========================
  do i = 1, NBasisGam

    Coef(:, :, i) = PolynomialGam(:, :, i)

    do j = 1, i-1
      projtmp = projector_Gam(t1min, t1max, t2min, t2max, Coef(:,:,j), PolynomialGam(:,:,i))
      Coef(:, :, i) = Coef(:,:,i) - Coef(:,:,j) *projtmp
    enddo

    norm = dsqrt(projector_Gam(t1min, t1max, t2min, t2max, Coef(:, :, i), Coef(:,:,i)))
    Coef(:, :, i) = Coef(:, :, i)/norm
  enddo

  call test_basis_Gam(t1min, t1max, t2min, t2max, Coef)

END SUBROUTINE calculate_basis_Gamma_2D


DOUBLE PRECISION FUNCTION projector_Gam(t1min, t1max, t2min, t2max, Coef1, Coef2)
  implicit none
  integer, intent(in) :: t1min, t1max
  integer, intent(in) :: t2min(0:MxT-1), t2max(0:MxT-1)
  double precision, intent(in) :: Coef1(0:BasisOrderGam, 0:BasisOrderGam), Coef2(0:BasisOrderGam, 0:BasisOrderGam)
  integer :: it1, jt1, t1, it2, jt2, t2
  double precision :: tau1, tau2

  projector_Gam  = 0.d0
  do t1 = t1min, t1max
    do t2 = t2min(t1), t2max(t1)
      tau1 = dble(t1)*Beta/dble(MxT)
      tau2 = dble(t2)*Beta/dble(MxT)
      do it1 = 0, BasisOrderGam 
        do jt1 = 0, BasisOrderGam
          do it2 = 0, BasisOrderGam 
            do jt2 = 0, BasisOrderGam
              projector_Gam = projector_Gam + (Beta/dble(MxT))**2.d0*Coef1(it1,it2)*Coef2(jt1, &
                &  jt2)*tau1**(dble(it1+jt1))*tau2**(dble(it2+jt2))
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
  return 
END FUNCTION projector_Gam

DOUBLE PRECISION FUNCTION weight_basis_Gam(Coef, tau1, tau2)
  implicit none
  double precision :: Coef(0:BasisOrderGam, 0:BasisOrderGam)
  integer :: j1, t1, j2, t2
  double precision :: tau1, tau2
  weight_basis_Gam = 0.d0
  do j1 = 0, BasisOrderGam
    do j2 = 0, BasisOrderGam
      weight_basis_Gam = weight_basis_Gam + Coef(j1, j2)*tau1**dble(j1)*tau2**dble(j2)
    enddo
  enddo
  return 
END FUNCTION weight_basis_Gam


!--------- test subroutine for the Gram-Schmidt basis ------
SUBROUTINE test_basis_Gam(t1min, t1max, t2min, t2max, Coef)
  implicit none
  integer, intent(in) :: t1min, t1max
  integer, intent(in) :: t2min(0:MxT-1), t2max(0:MxT-1)
  double precision, intent(in) :: Coef(0:BasisOrderGam, 0:BasisOrderGam, NbasisGam)
  integer :: i, j, k, l, n, t1, t2
  double precision :: tau, tau1, tau2
  double precision :: omega, x, y

  do i = 1, NBasisGam
    y = 0.d0
    do  t1 = t1min, t1max
      tau1 = dble(t1)*Beta/dble(MxT)
      do  t2 = t2min(t1), t2max(t1)
        tau2 = dble(t2)*Beta/dble(MxT)

        y = y + (Beta/dble(MxT))**2.d0*weight_basis_Gam(Coef(:, :, i), tau1, tau2)**2.d0
      enddo
    enddo
    if(dabs(y-1.d0)>1.d-8) then
      call LogFile%QuickLog(str(i)+str(y)+"Gamma Basis error!! Not normal!")
    endif
  enddo

  do i = 1, NBasisGam
    do j = i+1, NBasisGam
      y = projector_Gam(t1min, t1max, t2min, t2max, Coef(:, :, i), Coef(:, :, j))
      if(dabs(y)>5.d-7) then
        call LogFile%QuickLog(str(i)+str(j)+str(y)+"Gamma Basis error!! Not orthogonal!")
      endif
    enddo
  enddo

END SUBROUTINE test_basis_Gam

!END PROGRAM MAIN
