
!!====== transfer to a D-dimensional coordinates from a 1-d array =============
function get_cord_from_site(dims, site)
  implicit none
  integer, intent(in) :: dims, site
  integer, dimension(dims) :: get_cord_from_site
  integer :: i, tmp

  tmp = site
  do i = dims, 1, -1
    get_cord_from_site(i) =  tmp/dVol(i)
    tmp = tmp - get_cord_from_site(i)*dVol(i)
  enddo
  return
END FUNCTION get_cord_from_site

!!====== transfer a D-dimensional space into a 1-d array =============
integer function get_site_from_cord(dims, cord)
  implicit none
  integer :: dims
  integer, dimension(dims) :: cord
  integer :: i, j

  get_site_from_cord = 0
  do i = 1, dims
    get_site_from_cord = get_site_from_cord + cord(i)*dVol(i)
  enddo
  return
END FUNCTION get_site_from_cord


!======= decompose a matrix from 1d array to D-dimensional matrix ==============
SUBROUTINE decompose_matrix(Ntyp, Lx, Ly, Lz,Nleft, Mat1D, Mat)
  implicit none
  integer, intent(in) :: Ntyp, Lx, Ly, Lz, Nleft
  complex*16, intent(in) :: Mat1D(Ntyp,0:Lx*Ly*Lz-1,0:Nleft-1)
  complex*16, intent(out):: Mat(Ntyp,0:Lx-1,0:Ly-1,0:Lz-1,0:Nleft-1)
  integer :: it, il, ix, iy, iz, site
  do it = 1, Ntyp
    do il = 0, Nleft-1
      do ix = 0, Lx-1
        do iy = 0, Ly-1
          do iz = 0, Lz-1
            site = get_site_from_cord(3, (/ix, iy, iz/))
            Mat(it, ix, iy, iz, il) = Mat1D(it, site, il)
          enddo
        enddo
      enddo
    enddo
  enddo
  return
END SUBROUTINE decompose_matrix



!======= combine a matrix from D-dimensional matrix to 1-d matrix ==============
SUBROUTINE combine_matrix(Ntyp, Lx, Ly, Lz, Nleft, Mat, Mat1D)
  implicit none
  integer, intent(in) :: Ntyp, Lx, Ly, Lz, Nleft
  complex*16, intent(in):: Mat(Ntyp,0:Lx-1,0:Ly-1,0:Lz-1,0:Nleft-1)
  complex*16, intent(out) :: Mat1D(Ntyp,0:Lx*Ly*Lz-1,0:Nleft-1)
  integer :: it, il, ir(3), site
  do it = 1, Ntyp
    do il = 0, Nleft-1
      do site = 0, Lx*Ly*Lz-1
        ir = get_cord_from_site(3, site)
        Mat1D(it, site, il) = Mat(it, ir(1), ir(2), ir(3), il) 
      enddo
    enddo
  enddo
  return
END SUBROUTINE combine_matrix



!!------------- definition of the lattice and type symmetry ----------------
SUBROUTINE def_symmetry
  implicit none
  integer :: site, dt1, dt2, i
  integer :: dr(1:D)

  CoefOfSymmetry(:, :, :) = 2.d0   !typ 1,2; 3,4; 5,6
  
  do site = 0, Vol-1
    dr = get_cord_from_site(D, site)
    do i = 1, D
      if(dr(i)/=0 .and. dr(i)/=L(i)/2) then
        CoefOfSymmetry(site, :, :) = CoefOfSymmetry(site, :, :)*2.d0
      endif
    enddo
  enddo

  do dt1 = 0, MxT-1
    do dt2 = 0, MxT-1
      CoefOfSymmetry(:, dt1, dt2) = CoefOfSymmetry(:, dt1, dt2)*2.d0
    enddo
  enddo

  return
END SUBROUTINE def_symmetry


Integer FUNCTION diff_r(dims, site1, site2)
  implicit none
  integer, intent(in) :: site1, site2
  integer, intent(in) :: dims
  integer :: r1(dims), r2(dims), dr(dims)

  r1 = get_cord_from_site(dims, site1)
  r2 = get_cord_from_site(dims, site2)

  dr = r1 - r2

  do i = 1, dims
    if(dr(i)<0)  dr(i) = dr(i)+L(i)
    if(dr(i)>dL(i))  dr(i) = L(i)-dr(i)
  enddo

  diff_r = get_site_from_cord(dims, dr)
  return
END FUNCTION diff_r

!!======================== WEIGHT EXTRACTING =========================
!! most basic interface to the matrix element

!--------- weight for bare propagator ----------------
Complex*16 FUNCTION weight_G0(typ, t)
  implicit none
  integer, intent(in)  :: typ, t  
  double precision     :: tau
  complex(kind=8)      :: muc  

  muc = dcmplx(0.d0, Mu(1)*pi/(2.d0*Beta))
  tau = (real(t)+0.5d0)*Beta/MxT
  if(tau>=0) then
    weight_G0 = cdexp(muc*tau)/(1.d0, 1.d0) 
  else
    weight_G0 = -cdexp(muc*(tau+Beta))/(1.d0, 1.d0) 
  endif
  return
END FUNCTION weight_G0


!!--------- calculate weight for bare interaction ----
Complex*16 FUNCTION weight_W0(typ, site)
  implicit none
  integer, intent(in) :: site, typ
  double precision :: ratio

  weight_W0 = (0.d0, 0.d0)

  if(is_W0_nonzero(D, site)==1) then
    if(typ ==1 .or. typ == 2) then
      weight_W0 = dcmplx(0.25d0, 0.d0)
    else if(typ == 3 .or. typ == 4) then
      weight_W0 = dcmplx(-0.25d0, 0.d0)
    else if(typ == 5 .or. typ == 6) then
      weight_W0 = dcmplx(0.5d0, 0.d0)
    endif
  endif

  if(Is_J1J2 .and. is_W0_nonzero(D, site)==2) then
    if(typ ==1 .or. typ == 2) then
      weight_W0 = dcmplx(0.25d0*Jcp, 0.d0)
    else if(typ == 3 .or. typ == 4) then
      weight_W0 = dcmplx(-0.25d0*Jcp, 0.d0)
    else if(typ == 5 .or. typ == 6) then
      weight_W0 = dcmplx(0.5d0*Jcp, 0.d0)
    endif
  endif

END FUNCTION weight_W0

!!--------- calculate weight for bare gamma ---------
COMPLEX*16 FUNCTION weight_Gam0(typ, site)
  implicit none
  integer, intent(in)  :: site, typ
  double precision :: ratio

  weight_Gam0 = (0.d0, 0.d0)

  if(is_Gam0_nonzero(D, site)) then
    if(typ==1 .or. typ==2 .or. typ==5 .or. typ==6) then
      weight_Gam0 = (1.d0, 0.d0)
    endif
  endif
END FUNCTION weight_Gam0


!======================= Complex Operations =============================
COMPLEX*16 FUNCTION d_times_cd(nd, ncd)
  implicit none 
  double precision, intent(in) :: nd
  complex*16, intent(in) :: ncd
  d_times_cd = dcmplx(nd*real(ncd), nd*dimag(ncd))
  return
END FUNCTION d_times_cd


!================ find the lowest denominator in brillouin zone ===========
COMPLEX*16 FUNCTION find_lowest_in_bz(Matrix, klow)
  implicit none
  complex*16, intent(in) :: Matrix(0:Vol-1)
  integer, intent(out) :: klow(1:D)
  integer :: ik
  find_lowest_in_bz = Matrix(0)
  do  ik = 0, Vol-1
    if(real(Matrix(ik))-real(find_lowest_in_bz)<-1.d-12) then
      find_lowest_in_bz = Matrix(ik)
      klow = get_cord_from_site(D, ik)
    endif
  enddo
  return 
END FUNCTION find_lowest_in_bz




!!=======================================================================
!!======================= Fourier Transformation ========================
!!=======================================================================

SUBROUTINE transfer_r(Backforth)
  implicit none
  integer,intent(in) :: Backforth

  call transfer_W_r(Backforth)
  call transfer_Gam_r(Backforth)
  return
END SUBROUTINE



SUBROUTINE transfer_t(Backforth)
  implicit none
  integer,intent(in) :: Backforth

  call transfer_G_t(Backforth)
  call transfer_W_t(Backforth)
  call transfer_Gam_t(Backforth)
END SUBROUTINE

!================================================================================
!============ transfer quantities in (r,k) domain ===========================
!================================================================================

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

SUBROUTINE transfer_Chi_r(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    call FFT_r(Chi,1,MxT,BackForth)
    return
END SUBROUTINE

SUBROUTINE transfer_Polar_r(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    call FFT_r(Polar,1,MxT,BackForth)
END SUBROUTINE

!!================================================================================
!!============ transfer quantities in (t,omega) domain ===========================
!!================================================================================

!!!TODO: need to be update!!!!
!========W0 and Gam0 are just constants in omega domain ====================
SUBROUTINE transfer_W0_t(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    call FFT_tau_single(W0PF,1,Vol,.false.,BackForth)
END SUBROUTINE

SUBROUTINE transfer_Gam0_t(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer :: it1, it2

    call FFT_tau_double(Gam0PF,1,Vol,.false.,BackForth)
END SUBROUTINE
!===========================================================================

SUBROUTINE transfer_G0(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer :: it
    call FFT_tau_single(G0F,1,1,.true.,BackForth)
END SUBROUTINE

SUBROUTINE transfer_GamOrder1_t(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer ::it1, it2
    call FFT_tau_double(GamOrder1,NtypeGam,1,.true.,BackForth)
END SUBROUTINE

SUBROUTINE transfer_G_t(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer            :: it
    call FFT_tau_single(G,NtypeG,1,.true.,BackForth)
END SUBROUTINE transfer_G_t

SUBROUTINE transfer_Gam_t(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer :: it1, it2
    call FFT_tau_double(Gam,NtypeGam,Vol,.true.,BackForth)
END SUBROUTINE

SUBROUTINE transfer_Sigma_t(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer       :: it
    call FFT_tau_single(Sigma,1,1,.true.,BackForth)
END SUBROUTINE

SUBROUTINE transfer_Polar_t(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    call FFT_tau_single(Polar,1,Vol,.false.,BackForth)
END SUBROUTINE

SUBROUTINE transfer_W_t(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    call FFT_tau_single(W,NtypeW,Vol,.false.,BackForth)
END SUBROUTINE transfer_W_t

SUBROUTINE transfer_Chi_t(BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    call FFT_tau_single(Chi,1,Vol,.false.,BackForth)
END SUBROUTINE


SUBROUTINE FFT_r(XR,Ntype,Nz,BackForth)
    implicit none
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer,intent(in) :: Nz
    integer,intent(in) :: Ntype
    complex*16 :: XR(Ntype,0:Vol-1,0:Nz-1)
    complex*16,allocatable :: XR3D(:,:,:,:,:)
    integer :: Power,Noma
    integer :: ix, iy, iz, iother,it
    integer :: FL(3)
    double precision,allocatable ::   Real1(:), Im1(:)
    
    do i =  1, 3
      if(D>=i) then
        FL(i) = L(i)
      else
        FL(i) = 1
      endif
    enddo

    allocate(XR3D(Ntype,0:FL(1)-1, 0:FL(2)-1, 0:FL(3)-1, 0:Nz-1))
    
    call decompose_matrix(Ntype, FL(1), FL(2), FL(3), Nz, XR, XR3D)


    !do FFT in x direction
    if(D>=1) then
      Power=log(FL(1)*1.d0)/log(2.d0)+1.d-14  
      Noma=2**power
      allocate(Real1(0:Noma-1),Im1(0:Noma-1))

      do iother=0,Nz-1
        do iy=0,FL(2)-1
          do iz=0,FL(3)-1
            do it=1,Ntype
              Real1(:)=real(XR3D(it,:,iy,iz,iother))
              Im1(:)=dimag(XR3D(it,:,iy,iz,iother))
              call sffteu(Real1, Im1, Noma, Power, BackForth)
              XR3D(it,:,iy,iz,iother)=dcmplx(Real1(:), Im1(:))
            enddo
          enddo
        enddo
      enddo
      deallocate(Real1,Im1)
    endif

    !do FFT in y direction
    if(D>=2) then
      power=log(FL(2)*1.d0)/log(2.d0)+1.d-14  
      Noma=2**Power
      allocate(Real1(0:Noma-1),Im1(0:Noma-1))
      do iother=0,Nz-1
        do ix=0,FL(1)-1
          do iz=0,FL(3)-1
            do it=1,Ntype
              Real1(:)=real(XR3D(it,ix,:,iz,iother))
              Im1(:)=dimag(XR3D(it,ix,:,iz,iother))
              call sffteu(Real1, Im1, Noma, Power, BackForth)
              XR3D(it,ix,:,iz,iother)=dcmplx(Real1(:), Im1(:))
            enddo
          enddo
        enddo
      enddo
      deallocate(Real1,Im1)
    endif

    !do FFT in z direction
    if(D>=3) then
      power=log(FL(3)*1.d0)/log(2.d0)+1.d-14  
      Noma=2**Power
      allocate(Real1(0:Noma-1),Im1(0:Noma-1))
      do iother=0,Nz-1
        do ix=0,FL(1)-1
          do iy=0,FL(2)-1
            do it=1,Ntype
              Real1(:)=real(XR3D(it,ix,iy,:,iother))
              Im1(:)=dimag(XR3D(it,ix,iy,:,iother))
              call sffteu(Real1, Im1, Noma, Power, BackForth)
              XR3D(it,ix,iy,:,iother)=dcmplx(Real1(:), Im1(:))
            enddo
          enddo
        enddo
      enddo
      deallocate(Real1,Im1)
    endif

    call combine_matrix(Ntype, FL(1), FL(2), FL(3), Nz, XR3D, XR)
    deallocate(XR3D)

end SUBROUTINE FFT_r


SUBROUTINE FFT_tau_single(XR,Ntype,Nz,IsAntiSym,BackForth)
    implicit none
    logical,intent(in) :: IsAntiSym    !IsAntiSym=.true. antisymmetric function
                                       !IsAntiSym=.false. symmetric function
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer,intent(in) :: Nz
    integer,intent(in) :: Ntype
    complex(kind=8)    :: XR(Ntype,0:Nz-1,0:MxT-1)
    integer :: Power,Noma
    integer :: itau,iomega,iz,it
    double precision,allocatable ::   Real1(:), Im1(:)


    if(BackForth/=-1) then
      !multiply exp(-i*Pi*(n+1/2)/MxT) term to anti-symmetric tau function
      if(IsAntiSym) then
        do it = 0, MxT-1
          XR(:,:,it) = XR(:,:,it)* cdexp(dcmplx(0.d0,-Pi*(real(it)+0.5d0)/real(MxT)))
        enddo
      endif
    else
      !shift tau with 1/2
      !multiply  exp(i*Pi*m/MxT) term to omega function
      !do it = 0, MxT-1
        !XR(:,:,it) = XR(:,:,it)* cdexp(dcmplx(0.d0,Pi*real(it)/real(MxT)))
      !enddo
    endif

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

    if(BackForth/=-1) then
      !shift tau with 1/2
      !multiply  exp(-i*Pi*m/MxT) term to omega function
      !do it = 0, MxT-1
        !XR(:,:,it) = XR(:,:,it)* cdexp(dcmplx(0.d0,-Pi*real(it)/real(MxT)))
      !enddo
    else
      !multiply exp(i*Pi*(n+1/2)/MxT) term to anti-symmetric tau function
      if(IsAntiSym) then
        do it = 0, MxT-1
          XR(:,:,it) = XR(:,:,it)* cdexp(dcmplx(0.d0,Pi*(real(it)+0.5d0)/real(MxT)))
        enddo
      endif
    endif

end SUBROUTINE FFT_tau_single
    
SUBROUTINE FFT_tau_double(XR,Ntype,Nz,IsAntiSym,BackForth)
    implicit none
    logical,intent(in) :: IsAntiSym    !IsAntiSym=.true. antisymmetric function
                                       !IsAntiSym=.false. symmetric function
    integer,intent(in) :: BackForth    !Backforth=-1 reverse tranformation
    integer,intent(in) :: Nz
    integer,intent(in) :: Ntype
    complex(kind=8)    :: XR(Ntype,0:Nz-1,0:MxT-1,0:MxT-1)
    integer :: Power,Noma
    integer :: itau1,itau2,iz,it,it1,it2
    double precision,allocatable ::   Real1(:), Im1(:)


    if(BackForth/=-1) then
      !multiply exp(-i*Pi*(n1+1/2+n2+1/2)/MxT) term to anti-symmetric tau function
      if(IsAntiSym) then
        do it2 = 0, MxT-1
          do it1 = 0, MxT-1
            XR(:,:,it1, it2) = XR(:,:,it1, it2)* cdexp(dcmplx(0.d0,-Pi*real(it1+it2+1.d0)/real(MxT)))
          enddo
        enddo
      endif
    else
      !shift tau with 1/2
      !multiply  exp(i*Pi*(m1+m2)/MxT) term to omega function
      !do it2 = 0, MxT-1
        !do it1 = 0, MxT-1
          !XR(:,:,it1, it2) = XR(:,:,it1, it2)* cdexp(dcmplx(0.d0,Pi*real(it1+it2)/real(MxT)))
        !enddo
      !enddo
    endif

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

    if(BackForth/=-1) then
      !shift tau with 1/2
      !multiply  exp(-i*Pi*(m1+m2)/MxT) term to omega function
      !do it2 = 0, MxT-1
        !do it1 = 0, MxT-1
          !XR(:,:,it1, it2) = XR(:,:,it1, it2)* cdexp(dcmplx(0.d0,-Pi*real(it1+it2)/real(MxT)))
        !enddo
      !enddo
    else
      !multiply exp(i*Pi*(n1+1/2+n2+1/2)/MxT) term to anti-symmetric tau function
      if(IsAntiSym) then
        do it2 = 0, MxT-1
          do it1 = 0, MxT-1
            XR(:,:,it1, it2) = XR(:,:,it1, it2)* cdexp(dcmplx(0.d0,Pi*real(it1+it2+1.d0)/real(MxT)))
          enddo
        enddo
      endif
    endif
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
!end program
!!=======================================================================
!!=======================================================================
!!=======================================================================



!!=======================================================================
!!================ RANDOM NUMBER GENERATORS =============================
!!=======================================================================

!---------- shift register random number generator ------
SUBROUTINE set_RNG
  implicit none
  integer :: i_r, k_r, k1_r, iseed

  nrannr = mxrn
  iseed  = iabs(Seed)+1
  k_r    = 3**18+2*iseed
  k1_r   = 1313131*iseed
  k1_r   = k1_r-(k1_r/mod2)*mod2

  do i_r = 1, len1
    k_r  = k_r *mult
    k1_r = k1_r*mul2
    k1_r = k1_r-(k1_r/mod2)*mod2
    ir1(i_r) = k_r+k1_r*8193
  enddo

  do i_r = 1, len2
    k_r  = k_r *mult
    k1_r = k1_r*mul2
    k1_r = k1_r-(k1_r/mod2)*mod2
    ir2(i_r) = k_r+k1_r*4099
  enddo

  do i_r = 1, len1
    inxt1(i_r) = i_r+1
  enddo
  inxt1(len1) = 1
  ipnt1 = 1
  ipnf1 = ifd1+1

  do i_r = 1, len2
    inxt2(i_r) = i_r+1
  enddo
  inxt2(len2) = 1
  ipnt2 = 1
  ipnf2 = ifd2 + 1
END SUBROUTINE set_RNG 

!---------- calculate next random number --------------
DOUBLE PRECISION FUNCTION rn()
  implicit none
  integer   :: i_r, l_r, k_r
  nrannr = nrannr +1
  if(nrannr>=mxrn) then
    nrannr = 1
    do i_r= 1, mxrn
      l_r = ieor(ir1(ipnt1),ir1(ipnf1))
      k_r = ieor(ir2(ipnt2),ir2(ipnf2))
      irn(i_r) = ieor(k_r,l_r)
      ir1(ipnt1)=l_r
      ipnt1 = inxt1(ipnt1)
      ipnf1 = inxt1(ipnf1)
      ir2(ipnt2) = k_r
      ipnt2 = inxt2(ipnt2)
      ipnf2 = inxt2(ipnf2)
    enddo
  endif 
  rn = irn(nrannr)*tm32+0.5d0
END FUNCTION rn

!===================================================================



!==============Trace elapsed time ==================================
!! THIS IS PROJECT-INDEPENDENT 
SUBROUTINE set_time_elapse
  implicit none
  !-- read and calculate time (in seconds) -------------------------
  call date_and_time(date, time, zone, tval)
  t_curr = tval(5)*3600.d0+tval(6)*60.d0+tval(7)+tval(8)*0.001d0 
  h_curr = tval(5)
  t_prev = t_curr
  h_prev = h_curr
  return
END SUBROUTINE set_time_elapse
  

!==============Trace elapsed time ==================================
!! THIS IS PROJECT-INDEPENDENT 
SUBROUTINE time_elapse
  implicit none
  
  !-- read and calculate time (in seconds) -------------------------
  call date_and_time(date, time, zone, tval)
  t_curr = tval(5)*3600.d0+tval(6)*60.d0+tval(7)+tval(8)*0.001d0 
  h_curr = tval(5)

  t_elap = t_curr-t_prev
  if(h_curr<h_prev) t_elap = t_elap+24*3600.d0
  t_prev = t_curr
  h_prev = h_curr 
  return
END SUBROUTINE time_elapse

