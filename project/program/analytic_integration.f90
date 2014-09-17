
!=============== Notice ======================
! in tau discrete integral
! t1-t2 should be represented as t1-t2-1
! -t should be -t-1
!=============================================

!====================== the order-0 of Gamma ======================

SUBROUTINE calculate_GamNormWeight
  implicit none
  integer :: t, t1, t2, t3, ityp
  complex*16 :: Gam0

  GamNormWeight = (0.d0, 0.d0)

  !--------- bare Gamma --------------------
  do ityp = 1, 6
    Gam0 = weight_Gam0(ityp, 0)
    Gam0 = Gam0 *weight_Gam0(ityp, 0)
    GamNormWeight = GamNormWeight + Gam0*(real(MxT)/Beta)**2.d0
  enddo

  return
END SUBROUTINE calculate_GamNormWeight


!====================== the 1st order of Gamma ======================
SUBROUTINE calculate_Gam1
  implicit none
  integer :: t1, t2, ityp
  integer :: omega1, omega2, omega3
  integer :: omegaW, omegaout
  integer :: tg(5), tgam1(5), tgam2(5), tgam3(5), tw(5)
  integer :: x1, y1, x3, y3, isite1, isite3, isiteW
  complex*16 :: Gin, Gout, iW, Gam1, Gam2, Gam3
  complex*16 :: weight
  complex*16 :: FGam(0:MxT-1, 0:MxT-1)
  DOUble precision :: ratio

  !!============ fake G function =================================
  !do t1 = 0, MxT-1
    !!G(:, t1) = cdexp((0.d0, 1.d0)*pi*real(t1)/(2.d0*MxT))/(1.d0, 1.d0)
    !G(:, t1) = (1.d0, 0.d0)
  !enddo

  !!============ fake W function =================================
  !W = (0.d0, 0.d0)
  !W(1,0,0,:) = (1.d0, 0.d0)
  !W(2,0,0,:) = (1.d0, 0.d0)
  !W(3,0,0,:) = (-1.d0, 0.d0)
  !W(4,0,0,:) = (-1.d0, 0.d0)
  !W(5,0,0,:) = (2.d0, 0.d0)
  !W(6,0,0,:) = (2.d0, 0.d0)

  !!============ fake Gamma function =================================
  !Gam = (0.d0, 0.d0)
  !do t1 = 0, MxT-1
    !do t2 = 0, MxT-1
      !!FGam(t1, t2) = dcmplx((t1*Beta/MxT)**2.d0+(t2*Beta/MxT)**2.d0+1.d0, 0.d0)
      !FGam(t1, t2) = (1.d0, 0.d0)
    !enddo
  !enddo
  !Gam(1,0,0,:,:) = FGam(:, :)
  !Gam(2,0,0,:,:) = FGam(:, :)
  !Gam(3,0,0,:,:) = FGam(:, :)
  !Gam(4,0,0,:,:) = FGam(:, :)
  !Gam(5,0,0,:,:) = FGam(:, :)
  !Gam(6,0,0,:,:) = FGam(:, :)
  !================== bold gamma ===============================

  call transfer_G_t(1)
  call transfer_W_t(1)
  call transfer_Gam_t(1)

  tg(1:4)  = 1;   tgam2(1:4) = 1
  tgam1(1) = 1;   tgam3(1)   = 1;        tw(1) = 1
  tgam1(4) = 3;   tgam3(4)   = 3;        tw(4) = 2
  tgam1(3) = 1;   tgam3(3)   = 3;        tw(3) = 3
  tgam1(2) = 3;   tgam3(2)   = 1;        tw(2) = 4

  tg(5)   = 2;    tgam2(5)   = 4
  tgam1(5) = 5;   tgam3(5)   = 6;        tw(5) = 5

  ratio = -1.d0/MxT *(Beta/MxT)**6.d0

  !GamOrder1(:,:,:) = (0.d0, 0.d0)
  !do omega2 = 0, MxT-1
    !do omega1 = 0, MxT-1

      !do ityp = 1, 5
        !do omega3 = 0, MxT-1

          !omegaW = omega1-omega3
          !if(omegaW<0)  omegaW = omegaW+MxT

          !omegaout = omega2-omegaW
          !if(omegaout<0)  omegaout = omegaout+MxT

          !Gin = G(tg(ityp), omega3)
          !Gout = G(tg(ityp), omegaout)
          !Gam2 = Gam(tgam2(ityp), 0, omega3, omegaout)+weight_Gam0(tgam2(ityp), 0)

          !do isite1 = 0, Vol-1
            !do isite3 = 0, Vol-1

              !isiteW = diff_r(D, isite1, isite3)
              !Gam1 = Gam(tgam1(ityp), isite1, omega1, omega3)+weight_Gam0(tgam1(ityp), 0)
              !Gam3 = Gam(tgam3(ityp), isite3, omegaout, omega2)+weight_Gam0(tgam3(ityp), 0)
              !iW = W(tw(ityp), isiteW, omegaW)

              !weight = Gin *Gout *iW *Gam1 *Gam2 *Gam3
              !GamOrder1(1,omega1, omega2)=GamOrder1(1, omega1, omega2)+ratio*weight
            !enddo
          !enddo

        !enddo
      !enddo

    !enddo
  !enddo

  call transfer_G_t(-1)
  call transfer_W_t(-1)
  call transfer_Gam_t(-1)

  call transfer_GamOrder1_t(-1)


  !================== bare Gamma ===============================
  !GamOrder1(:,:,:) = (0.d0, 0.d0)
  !do t2 = 0, MxT-1
    !do t1 = 0, MxT-1

        !Gin = weight_G(1, t1)
        !Gout = weight_G(1, t2)
        !Gam1 = weight_Gam0(1, 0)
        !Gam2 = weight_Gam0(1, 0) 
        !Gam3 = weight_Gam0(1, 0) 
        !if(t1+t2<MxT) then
          !iW = weight_W(1, 0, t1+t2)
        !else 
          !iW = weight_W(1, 0, t1+t2-MxT)
        !endif

        !weight = Gin *Gout *iW *Gam1 *Gam2 *Gam3

        !GamOrder1(1,t1,t2) = -1.d0* weight
    !enddo
  !enddo
END SUBROUTINE calculate_Gam1
!====================================================================




