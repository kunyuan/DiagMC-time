INCLUDE "../mylib/mylib.f90"

MODULE test
    implicit none
    integer, parameter,dimension(2) :: L =(/4,4/)   ! the system size
    integer,parameter :: MxT=64
    integer,parameter :: NtypeG=2
    integer,parameter :: NtypeW=6
    integer,parameter :: NtypeGam=6
    double precision            :: Beta
    double precision, parameter :: Pi    = 3.14159265358979323846d0
    double precision, parameter :: Pi2   = 6.2831853071795865d0
    double precision, parameter :: Pi4   = 2.0d0*Pi2
    complex*16 :: G(NtypeG,0:MxT-1)
    complex*16 :: W(NtypeW,0:L(1)-1,0:L(2)-1,0:MxT-1)
    complex*16 :: W0PF(0:L(1)-1,0:L(2)-1,0:MxT-1)
    complex*16 :: Gam(NtypeGam,0:L(1)-1,0:L(2)-1,0:MxT-1, 0:MxT-1)
    complex*16 :: Gam0PF(0:L(1)-1,0:L(2)-1,0:MxT-1, 0:MxT-1)
end module

PROGRAM MAIN
  USE test
  USE string_basic
  implicit none
  integer :: it, it1, it2, ityp
  double precision :: ratio

  Beta = 0.9d0
  !------- test G0 ----------------------
  !do it = 0, MxT-1
    !do ityp = 1, NtypeG
      !G(ityp, it) = weight_G0(ityp, it)
    !enddo
  !enddo

  !call output_G(0)

  !call transfer_G_t(1)
  !call output_G(1)

  !call transfer_G_t(-1)
  !call output_G(2)

  !----- test W0PF -----------------------
  !W0PF(:,:,:) = (0.d0, 0.d0)
  !ratio = 0.25d0*Jcp*real(MxT)/Beta

  !W0PF(1,      0, 0) = dcmplx(ratio, 0.d0)
  !W0PF(L(1)-1, 0, 0) = dcmplx(ratio, 0.d0)
  !W0PF(0,      1, 0) = dcmplx(ratio, 0.d0)
  !W0PF(0, L(2)-1, 0) = dcmplx(ratio, 0.d0)

  !call transfer_W0_r(1)
  !call transfer_W0_t(1)

  !------ test Gam0PF-------------------
  ratio = (real(MxT)/Beta)**2.d0
  Gam0PF(:,:,:,:) = (0.d0, 0.d0)
  Gam0PF(0,0,0,0) = dcmplx(ratio, 0.d0)

  call output_Gam0(0)

  call transfer_Gam0_r(1)
  call transfer_Gam0_t(1)
  call output_Gam0(1)

  call transfer_Gam0_r(-1)
  call transfer_Gam0_t(-1)
  call output_Gam0(2)

  
CONTAINS

SUBROUTINE output_G(i)
  implicit none
  integer :: i, ityp, it1

  open(104, access='append', file='test_fourier.dat')
  write(104, *) "##################################G",trim(adjustl(str(i)))
  write(104, *) "#tau:", MxT
  do it1 = 0, MxT-1
    write(104, *)  real(G(1, it1)), dimag(G(1, it1))
  enddo
  write(104, *)
  close(104)
END SUBROUTINE output_G

SUBROUTINE output_Gam0(i)
  implicit none
  integer :: i, ityp, it1, it2

  open(104, access='append', file='test_fourier.dat')
  write(104, *) "##################################Gam0PF",trim(adjustl(str(i)))
  write(104, *) "#tau1:", MxT, ",tau2:", MxT
  do it2 = 0, MxT-1
    do it1 = 0, MxT-1
      write(104, *)  real(Gam0PF(0, 0, it1, it2)), dimag(Gam0PF(0, 0, it1, it2))
    enddo
  enddo
  write(104, *)
  close(104)
END SUBROUTINE output_Gam0


Complex*16 FUNCTION weight_G0(typ, t)
  implicit none
  integer, intent(in)  :: typ, t  
  double precision     :: tau
  complex(kind=8)      :: muc  

  muc = dcmplx(0.d0, pi/(2.d0*Beta))
  tau = (real(t)+0.5d0)*Beta/MxT
  if(tau>=0) then
    weight_G0 = cdexp(muc*tau)/(1.d0, 1.d0) 
  else
    weight_G0 = -cdexp(muc*(tau+Beta))/(1.d0, 1.d0) 
  endif
  return
END FUNCTION weight_G0


!!================================================================================
!!============ transfer quantities in (r,k) domain ===========================
!!================================================================================
SUBROUTINE transfer_r(Backforth)
  implicit none
  integer,intent(in) :: Backforth

  call transfer_W_r(Backforth)
  call transfer_Gam_r(Backforth)
  return
END SUBROUTINE

SUBROUTINE transfer_W0_r(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    call FFT_r(W0PF,1,MxT,BackForth)
END SUBROUTINE

SUBROUTINE transfer_Gam0_r(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer :: it1, it2
    call FFT_r(Gam0PF,1,MxT**2,BackForth)
end SUBROUTINE

SUBROUTINE transfer_W_r(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer :: ix, iy
    call FFT_r(W,NtypeW,MxT,BackForth)
END SUBROUTINE transfer_W_r

SUBROUTINE transfer_Gam_r(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    call FFT_r(Gam,NtypeGam,MxT**2,BackForth)
END SUBROUTINE

!!================================================================================
!!============ transfer quantities in (t,omega) domain ===========================
!!================================================================================
SUBROUTINE transfer_t(Backforth)
  implicit none
  integer,intent(in) :: Backforth

  call transfer_G_t(Backforth)
  call transfer_W_t(Backforth)
  call transfer_Gam_t(Backforth)
END SUBROUTINE



SUBROUTINE transfer_W0_t(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    call FFT_tau_single(W0PF,1,L(1)*L(2),BackForth)
END SUBROUTINE

SUBROUTINE transfer_Gam0_t(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer :: it1, it2

    call FFT_shift_t_double(Gam0PF,1,L(1)*L(2),BackForth)
    call FFT_tau_double(Gam0PF,1,L(1)*L(2),BackForth)
    call FFT_shift_omega_double(Gam0PF,1,L(1)*L(2),BackForth)
END SUBROUTINE

SUBROUTINE transfer_G_t(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer            :: it
    call FFT_shift_t_single(G,NtypeG,1,BackForth)
    call FFT_tau_single(G,NtypeG,1,BackForth)
    call FFT_shift_omega_single(G,NtypeG,1,BackForth)
END SUBROUTINE transfer_G_t

SUBROUTINE transfer_W_t(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    call FFT_tau_single(W,NtypeW,L(1)*L(2),BackForth)
END SUBROUTINE transfer_W_t


SUBROUTINE transfer_Gam_t(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer :: it1, it2
    call FFT_shift_t_double(Gam,NtypeGam,L(1)*L(2),BackForth)
    call FFT_tau_double(Gam,NtypeGam,L(1)*L(2),BackForth)
    call FFT_shift_omega_double(Gam,NtypeGam,L(1)*L(2),BackForth)
END SUBROUTINE


SUBROUTINE FFT_shift_t_single(XR, Ntype, Nz, BackForth)
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer,intent(in) :: Nz
    integer,intent(in) :: Ntype
    complex(kind=8)    :: XR(Ntype,0:Nz-1,0:MxT-1)
    integer :: it

    if(BackForth/=-1) then
      do it = 0, MxT-1
        XR(:,:,it) = XR(:,:,it)* cdexp(dcmplx(0.d0,-Pi*real(it)/real(MxT)))
      enddo
    else
      do it = 0, MxT-1
        XR(:,:,it) = XR(:,:,it)* cdexp(dcmplx(0.d0,Pi*real(it)/real(MxT)))
      enddo
    endif
    return
END SUBROUTINE FFT_shift_t_single

SUBROUTINE FFT_shift_omega_single(XR, Ntype, Nz, BackForth)
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer,intent(in) :: Nz
    integer,intent(in) :: Ntype
    complex(kind=8)    :: XR(Ntype,0:Nz-1,0:MxT-1)
    integer :: it

    if(BackForth/=-1) then
      do it = 0, MxT-1
        XR(:,:,it) = XR(:,:,it)* cdexp(dcmplx(0.d0,-Pi*(it+0.5d0)/real(MxT)))
      enddo
    else
      do it = 0, MxT-1
        XR(:,:,it) = XR(:,:,it)* cdexp(dcmplx(0.d0,Pi*(it+0.5d0)/real(MxT)))
      enddo
    endif
    return
END SUBROUTINE FFT_shift_omega_single


SUBROUTINE FFT_shift_t_double(XR, Ntype, Nz, BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer,intent(in) :: Nz
    integer,intent(in) :: Ntype
    complex(kind=8)    :: XR(Ntype,0:Nz-1,0:MxT-1,0:MxT-1)
    integer :: it1, it2

    if(BackForth/=-1) then
      do it2 = 0, MxT-1
        do it1 = 0, MxT-1
          XR(:,:,it1, it2) = XR(:,:,it1, it2)* cdexp(dcmplx(0.d0,-Pi*real(it1+it2)/real(MxT)))
        enddo
      enddo
    else
      do it2 = 0, MxT-1
        do it1 = 0, MxT-1
          XR(:,:,it1, it2) = XR(:,:,it1, it2)* cdexp(dcmplx(0.d0,Pi*real(it1+it2)/real(MxT)))
        enddo
      enddo
    endif
    return
END SUBROUTINE FFT_shift_t_double

SUBROUTINE FFT_shift_omega_double(XR, Ntype, Nz, BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer,intent(in) :: Nz
    integer,intent(in) :: Ntype
    complex(kind=8)    :: XR(Ntype,0:Nz-1,0:MxT-1,0:MxT-1)
    integer :: it1, it2

    if(BackForth/=-1) then
      do it2 = 0, MxT-1
        do it1 = 0, MxT-1
          XR(:,:,it1, it2) = XR(:,:,it1, it2)* cdexp(dcmplx(0.d0,-Pi*(it1+it2+1.d0)/real(MxT)))
        enddo
      enddo
    else
      do it2 = 0, MxT-1
        do it1 = 0, MxT-1
          XR(:,:,it1, it2) = XR(:,:,it1, it2)* cdexp(dcmplx(0.d0,Pi*(it1+it2+1.d0)/real(MxT)))
        enddo
      enddo
    endif
    return
END SUBROUTINE FFT_shift_omega_double

SUBROUTINE FFT_r(XR,Ntype,Nz,BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer,intent(in) :: Nz
    integer,intent(in) :: Ntype
    complex*16    :: XR(Ntype,0:L(1)-1,0:L(2)-1,0:Nz-1)
    integer :: Power,Noma
    integer :: ix,iy,iz,it
    double precision,allocatable ::   Real1(:), Im1(:)

    do iy=0,L(2)-1
      do ix=0,L(1)-1
        if(ix>L(1)/2 .and. iy<=L(2)/2) then
          XR(:,ix,iy,:)=XR(:,L(1)-ix,iy,:)
        else if(ix<=L(1)/2 .and. iy>L(2)/2) then
          XR(:,ix,iy,:)=XR(:,ix,L(2)-iy,:)
        else if(ix>L(1)/2 .and. iy>L(2)/2)  then
          XR(:,ix,iy,:)=XR(:,L(1)-ix,L(2)-iy,:)
        endif
      enddo
    enddo

    !do FFT in x direction
    Power=log(L(1)*1.d0)/log(2.d0)+1.d-14  
    Noma=2**power
    allocate(Real1(0:Noma-1),Im1(0:Noma-1))
    do iz=0,Nz-1
      do iy=0,L(2)-1
        do it=1,Ntype
           Real1(:)=real(XR(it,:,iy,iz))
           Im1(:)=dimag(XR(it,:,iy,iz))
           call sffteu(Real1, Im1, Noma, Power, BackForth)
           XR(it,:,iy,iz)=dcmplx(Real1(:), Im1(:))
         enddo
      enddo
    enddo
    deallocate(Real1,Im1)

    !do FFT in y direction
    power=log(L(2)*1.d0)/log(2.d0)+1.d-14  
    Noma=2**Power
    allocate(Real1(0:Noma-1),Im1(0:Noma-1))
    do iz=0,Nz-1
      do ix=0,L(1)-1
        do it=1,Ntype
           Real1(:)=real(XR(it,ix,:,iz))
           Im1(:)=dimag(XR(it,ix,:,iz))
           call sffteu(Real1, Im1, Noma, Power, BackForth)
           XR(it,ix,:,iz)=dcmplx(Real1(:), Im1(:))
         enddo
      enddo
    enddo
    deallocate(Real1,Im1)

end SUBROUTINE FFT_r



SUBROUTINE FFT_tau_single(XR,Ntype,Nz,BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer,intent(in) :: Nz
    integer,intent(in) :: Ntype
    complex(kind=8)    :: XR(Ntype,0:Nz-1,0:MxT-1)
    integer :: Power,Noma
    integer :: itau,iomega,iz,it
    double precision,allocatable ::   Real1(:), Im1(:)

    !do FFT in tau 
    Power=log(MxT*1.d0)/log(2.d0)+1.d-14  
    Noma=2**power
    allocate(Real1(0:Noma-1),Im1(0:Noma-1))
    do iz=0,Nz-1
      do it=1,Ntype
        Real1(:)=real(XR(it,iz,:))
        Im1(:)=dimag(XR(it,iz,:))
        call sffteu(Real1, Im1, Noma, Power, BackForth)
        XR(it,iz,:)=dcmplx(Real1(:), Im1(:))
       enddo
    enddo
    deallocate(Real1,Im1)

end SUBROUTINE FFT_tau_single
    
SUBROUTINE FFT_tau_double(XR,Ntype,Nz,BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer,intent(in) :: Nz
    integer,intent(in) :: Ntype
    complex(kind=8)    :: XR(Ntype,0:Nz-1,0:MxT-1,0:MxT-1)
    integer :: Power,Noma
    integer :: itau1,itau2,iz,it
    double precision,allocatable ::   Real1(:), Im1(:)

    !do FFT in tau1 direction
    Power=log(MxT*1.d0)/log(2.d0)+1.d-14  
    Noma=2**power
    allocate(Real1(0:Noma-1),Im1(0:Noma-1))
    do itau2=0,MxT-1
      do iz=0,Nz-1
        do it=1,Ntype
          Real1(:)=real(XR(it,iz,:,itau2))
          Im1(:)=dimag(XR(it,iz,:,itau2))
          call sffteu(Real1, Im1, Noma, Power, BackForth)
          XR(it,iz,:,itau2)=dcmplx(Real1(:),Im1(:))
        enddo
      enddo
    enddo
    deallocate(Real1,Im1)

    !do FFT in tau2 direction
    Power=log(MxT*1.d0)/log(2.d0)+1.d-14  
    Noma=2**power
    allocate(Real1(0:Noma-1),Im1(0:Noma-1))
    do itau1=0,MxT-1
      do iz=0,Nz-1
        do it=1,Ntype
            Real1(:)=real(XR(it,iz,itau1,:))
            Im1(:)=dimag(XR(it,iz,itau1,:))
            call sffteu(Real1, Im1, Noma, Power, BackForth)
            XR(it,iz,itau1,:)=dcmplx(Real1(:),Im1(:))
        enddo
      enddo
    enddo
    deallocate(Real1,Im1)

end SUBROUTINE FFT_tau_double
    
!-------------------------------------------------------------c
!                                                             c
!  Subroutine sffteu( x, y, n, m, itype )                     c
!                                                             c
!  This routine is a slight modification of a complex split   c
!  radix FFT routine presented by C.S. Burrus.  The original  c
!  program header is shown below.                             c
!                                                             c
!  Arguments:                                                 c
!     x - real array containing real parts of transform       c
!              sequence (in/out)                              c
!     y - real array containing imag parts of transform       c
!              sequence (in/out)                              c
!     n - integer length of transform (in)                    c
!     m - integer such that n = 2**m  (in)                    c
!     itype - integer job specifier (in)                      c
!              itype .ne. -1 --> foward transform             c
!              itype .eq. -1 --> backward transform           c
!                                                             c
!  The forward transform computes                             c
!     X(k) = sum_{j=0}^{N-1} x(j)*exp(-2ijk*pi/N)             c
!                                                             c
!  The backward transform computes                            c
!     x(j) = (1/N) * sum_{k=0}^{N-1} X(k)*exp(2ijk*pi/N)      c
!                                                             c
!  Requires standard FORTRAN functions - sin, cos             c
!                                                             c
!  Steve Kifowit, 9 July 1997                                 c
!                                                             c
!-------------------------------------------------------------C
!  A Duhamel-Hollman Split-Radix DIF FFT                      C
!  Reference:  Electronics Letters, January 5, 1984           C
!  Complex input and output in data arrays X and Y            C
!  Length is N = 2**M                                         C
!                                                             C
!  C.S. Burrus          Rice University         Dec 1984      C
!-------------------------------------------------------------C

subroutine SFFTEU( X, Y, N, M, ITYPE )
integer  N, M, ITYPE
real(kind=8) ::  X(*), Y(*)
integer  I, J, K, N1, N2, N4, IS, ID, I0, I1, I2, I3
real(kind=8) ::  E, A, A3, CC1, SS1, CC3, SS3
real(kind=8) ::  R1, R2, S1, S2, S3, XT
!      INTRINSIC  SIN, COS
real(kind=8), parameter :: TWOPI = 6.2831853071795864769 

if ( N .eq. 1 ) return

if ( ITYPE .eq. -1 ) then
   do I = 1, N
      Y(I) = - Y(I)
   end do
end if

N2 = 2 * N
do K = 1, M-1
   N2 = N2 / 2
   N4 = N2 / 4
   E = TWOPI / N2
   A = 0.0
   do J = 1, N4
      A3 = 3 * A
      CC1 = DCOS( A )
      SS1 = DSIN( A )
      CC3 = DCOS( A3 )
      SS3 = DSIN( A3 )
      A = J * E
      IS = J
      ID = 2 * N2
40        do I0 = IS, N-1, ID
         I1 = I0 + N4
         I2 = I1 + N4
         I3 = I2 + N4
         R1 = X(I0) - X(I2)
         X(I0) = X(I0) + X(I2)
         R2 = X(I1) - X(I3)
         X(I1) = X(I1) + X(I3)
         S1 = Y(I0) - Y(I2)
         Y(I0) = Y(I0) + Y(I2)
         S2 = Y(I1) - Y(I3)
         Y(I1) = Y(I1) + Y(I3)
         S3 = R1 - S2
         R1 = R1 + S2
         S2 = R2 - S1
         R2 = R2 + S1
         X(I2) = R1 * CC1 - S2 * SS1
         Y(I2) = - S2 * CC1 - R1 * SS1
         X(I3) = S3 * CC3 + R2 * SS3
         Y(I3) = R2 * CC3 - S3 * SS3
      end do
      IS = 2 * ID - N2 + J
      ID = 4 * ID
      if ( IS .lt. N ) goto 40
   end do
end do

!--------LAST STAGE, LENGTH-2 BUTTERFLY ----------------------C

IS = 1
ID = 4
50  do I0 = IS, N, ID
   I1 = I0 + 1
   R1 = X(I0)
   X(I0) = R1 + X(I1)
   X(I1) = R1 - X(I1)
   R1 = Y(I0)
   Y(I0) = R1 + Y(I1)
   Y(I1) = R1 - Y(I1)
end do
IS = 2 * ID - 1
ID = 4 * ID
if ( IS .lt. N ) goto 50

!-------BIT REVERSE COUNTER-----------------------------------C

100 J = 1
N1 = N - 1
do I = 1, N1
   if ( I .ge. J ) goto 101
   XT = X(J)
   X(J) = X(I)
   X(I) = XT
   XT = Y(J)
   Y(J) = Y(I)
   Y(I) = XT
101    K = N / 2
102    if ( K .ge. J ) goto 103
   J = J - K
   K = K / 2
   goto 102
103    J = J + K
end do

if ( ITYPE .eq. -1 ) then
   do I = 1, N
      X(I) = X(I) / N
      Y(I) = - Y(I) / N
   end do
end if

return

end  subroutine SFFTEU
END PROGRAM


