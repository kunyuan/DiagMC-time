!MODULE test
    !implicit none
    !integer,parameter :: Lx=4,Ly=4
    !integer,parameter :: Omega=2
    !integer,parameter :: NtypeW=7
    !double precision :: WR(NtypeW,0:Lx-1,0:Ly-1,-Omega:Omega,-Omega:Omega)
!end module

!program main
    !use test
    !integer :: ix,iy
    !print *, "Start"
    !WR(:,:,:,:,:)=1
    !WR(1,0,0,:,:)=3
    !write(*,*) "Before"
    !do iy=0,Ly-1
      !write(*,*) WR(1,:,iy,2,2)
    !enddo
    !call FFT(WR,NtypW,(2*Omega+1)**2,1)
    !print *, "First FFT"
    !do iy=0,Ly-1
      !write(*,*) WR(1,:,iy,2,2)
    !enddo
    !call FFT(WR,NtypW,(2*Omega+1)**2,-1)
    !print *, "Second FFT"
    !do iy=0,Ly-1
      !write(*,*) WR(1,:,iy,2,2)
    !enddo

!CONTAINS

SUBROUTINE transfer_W0(BackForth)
  implicit none
  integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
  call FFT(W0InMoment, NtypW,1, BackForth)
END SUBROUTINE

SUBROUTINE transfer_W(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    call FFT(WR,NtypW,MxOmegaW1+1,BackForth)
    call FFT(WRTailP,NtypW,Nbasis,BackForth)
    call FFT(WRTailC,NtypW,1,BackForth)
END SUBROUTINE

SUBROUTINE transfer_Gamma(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation

    call FFT(GamR,NtypGa,(2*MxOmegaGamG1+1)*(2*MxOmegaGamG1+1),BackForth)

    call FFT(GamRTailNP,NtypGa,NbasisGamma,BackForth)
    call FFT(GamRTailPN,NtypGa,NbasisGamma,BackForth)
    call FFT(GamRTailPPR,NtypGa,NbasisGamma,BackForth)
    call FFT(GamRTailNNR,NtypGa,NbasisGamma,BackForth)
    call FFT(GamRTailPPL,NtypGa,NbasisGamma,BackForth)
    call FFT(GamRTailNNL,NtypGa,NbasisGamma,BackForth)

    call FFT(GamRTailDiagP,NtypGa,Nbasis,BackForth)
    call FFT(GamRTailDiagN,NtypGa,Nbasis,BackForth)

    call FFT(GamRTailMP,NtypGa,Nbasis*(2*MxOmegaGamG1+1),BackForth)
    call FFT(GamRTailMN,NtypGa,Nbasis*(2*MxOmegaGamG1+1),BackForth)
    call FFT(GamRTailPM,NtypGa,Nbasis*(2*MxOmegaGamG1+1),BackForth)
    call FFT(GamRTailNM,NtypGa,Nbasis*(2*MxOmegaGamG1+1),BackForth)
    call FFT(GamRTailC,NtypGa,1,BackForth)
END SUBROUTINE

SUBROUTINE transfer_Pi(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    call FFT(PiR,ntypPi,2*MxOmegaChi+1,BackForth)
END SUBROUTINE

SUBROUTINE transfer_Chi(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    call FFT(ChiR,ntypCh,2*MxOmegaChi+1,BackForth)
END SUBROUTINE


SUBROUTINE FFT(XR,Ntype,Nz,BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer,intent(in) :: Nz
    integer,intent(in) :: Ntype
    double precision :: XR(Ntype,0:Lx-1,0:Ly-1,0:Nz-1)
    integer :: Power,Noma
    integer :: ix,iy,iz,ixt,iyt,it
    double precision,allocatable ::   Real1(:), Im1(:)

    do ix=0,Lx-1
      do iy=0,Ly-1
        if(ix>Lx/2 .and. iy<=Ly/2) then
          XR(:,ix,iy,:)=XR(:,Lx-ix,iy,:)
        else if(ix<=Lx/2 .and. iy>Ly/2) then
          XR(:,ix,iy,:)=XR(:,ix,Ly-iy,:)
        else if(ix>Lx/2 .and. iy>Ly/2)  then
          XR(:,ix,iy,:)=XR(:,Lx-ix,Ly-iy,:)
        endif
      enddo
    enddo

    !do FFT in x direction
    Power=log(Lx*1.d0)/log(2.d0)+1.d-14  
    Noma=2**power
    allocate(Real1(0:Noma-1),Im1(0:Noma-1))
    do it=1,Ntype
      do iy=0,Ly-1
         do iz=0,Nz-1
             Real1(:)=XR(it,:,iy,iz)
             Im1(:)=0.0
             call sffteu(Real1, Im1, Noma, Power, BackForth)
             XR(it,:,iy,iz)=Real1(:)
             if(sum(abs(Im1(:)))>1e-10) then
               write(*,*) sum(abs(Im1(:))), "Where is the imaginary part comes from?"
             endif
         end do
      end do
    enddo
    deallocate(Real1,Im1)

    !do FFT in y direction

    power=log(Ly*1.d0)/log(2.d0)+1.d-14  
    Noma=2**Power
    allocate(Real1(0:Noma-1),Im1(0:Noma-1))
    do it=1,Ntype
      do ix=0,Lx-1
         do iz=0,Nz-1
             Real1(:)=XR(it,ix,:,iz)
             Im1(:)=0.0
             call sffteu(Real1, Im1, Noma, Power, BackForth)
             XR(it,ix,:,iz)=Real1(:)
             if(sum(abs(Im1(:)))>1e-10) then
               write(*,*) sum(abs(Im1(:))), "Where is the imaginary part comes from?"
             endif
         end do
      end do
    enddo
    deallocate(Real1,Im1)

end SUBROUTINE FFT

    
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
!end program
