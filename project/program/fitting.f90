
!============== divide the space into N bins for G,W, Gamma =============================
SUBROUTINE initialize_bins
  implicit none
  integer :: t1, t2

  !================== define the bins for G,W, Gamma ====================================
  FromG(1) = 0
  ToG(1)   = MxT-1

  FromW(1) = 0
  ToW(1)   = MxT-1

  !============ Gamma with 3 bins ================================
  IsBasis2D(1) = .true.
  FromGamT1(1) = 0
  ToGamT1(1) = MxT-1
  do t1 = FromGamT1(1), ToGamT1(1)
    FromGamT2(t1, 1) = 0
    ToGamT2(t1, 1)   = MxT-2-t1
  enddo
  

  IsBasis2D(2) = .true.
  FromGamT1(2) = 1
  ToGamT1(2) = MxT-1
  do t1 = FromGamT1(2), ToGamT1(2)
    FromGamT2(t1, 2) = MxT-t1
    ToGamT2(t1, 2)   = MxT-1
  enddo

  return
END SUBROUTINE initialize_bins

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



INTEGER FUNCTION get_bin_Gam(it1, it2)
  implicit none
  integer :: it1, it2
  if(it1+it2<MxT-1)  get_bin_Gam = 1
  if(it1+it2>MxT-1)  get_bin_Gam = 2
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



DOUBLE PRECISION FUNCTION projector(tmin, tmax, Coef1, Coef2)
  implicit none
  integer, intent(in) :: tmin, tmax
  double precision, intent(in) :: Coef1(0:BasisOrder), Coef2(0:BasisOrder)
  integer :: j, k, t
  double precision :: tau

  projector  = 0.d0
  do t = tmin, tmax
    tau = (dble(t)+0.5d0)/dble(MxT)
    do j = 0, BasisOrder 
      do k = 0, BasisOrder
        projector = projector + 1.d0/dble(MxT)*Coef1(j)*Coef2(k)*tau**(dble(j+k))
      enddo
    enddo
  enddo
  return 
END FUNCTION projector

DOUBLE PRECISION FUNCTION weight_basis(Coef, t)
  implicit none
  double precision :: Coef(0:BasisOrder)
  integer :: j, t
  double precision :: tau

  tau = (dble(t)+0.5d0)/dble(MxT)
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
      y = y + 1.d0/dble(MxT)*weight_basis(Coef(:, i), t)**2.d0
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
      tau1 = (dble(t1)+0.5d0)/dble(MxT)
      tau2 = (dble(t2)+0.5d0)/dble(MxT)
      do it1 = 0, BasisOrderGam 
        do jt1 = 0, BasisOrderGam
          do it2 = 0, BasisOrderGam 
            do jt2 = 0, BasisOrderGam
              projector_Gam = projector_Gam + (1.d0/dble(MxT))**2.d0*Coef1(it1,it2)*Coef2(jt1, &
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
  integer :: j1, t1, j2, t2, it1, it2
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
      do  t2 = t2min(t1), t2max(t1)
        y = y + (1.d0/dble(MxT))**2.d0*weight_basis_Gam(Coef(:, :, i),  &
          &  (dble(t1)+0.5d0)/dble(MxT), (dble(t2)+0.5d0)/dble(MxT))**2.d0
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

