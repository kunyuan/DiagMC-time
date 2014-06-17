Program main
  implicit none
  !integer, parameter :: D=1
  integer, parameter :: D=2
  !integer, parameter :: D=3
  integer, parameter :: L=4
  integer, parameter :: Vol=L**D
  integer, parameter :: Ntyp = 6
  integer, parameter :: Nleft = 4
  integer :: ix,iy,iz,it,il, ir(D)
  integer :: site
  complex*16 :: F1D(Ntyp, 0:Vol-1, 0:Nleft-1)
  complex*16 :: F(Ntyp, 0:L-1, 0:L-1, 0:L-1, 0:Nleft-1)

  !============ test example for 1D ================

  !============ test example for 2D ================
  do site = 0, Vol-1
    do it = 1, Ntyp
      do il = 0, Nleft-1 
        ir = get_matrix_from_site(D, site)
        F1D(it, site, il) = cmplx(ir(1)**2.d0+ir(2)**2.d0, 0.0)
      enddo
    enddo
    print *, 1, site, 0, F1D(1, site, 0)
  enddo

  call decompose_matrix(Ntyp, L, L, 1, Nleft, F1D, F)

  do ix = 0, L-1
    do iy = 0, L-1
      iz = 0
      site = get_site_from_matrix(D, (/ix,iy,iz/))
      print *, ix, iy, iz, site, F(1, ix, iy, iz, 0)
    enddo
  enddo

  !============ test example for 3D ================
  !do site = 0, Vol-1
    !do it = 1, Ntyp
      !do il = 0, Nleft-1 
        !F1D(it, site, il) = cmplx(site, 0.0)
      !enddo
    !enddo
    !print *, 1, site, 0, F1D(1, site, 0)
  !enddo

  !call decompose_matrix(Ntyp, L, L, L, Nleft, F1D, F)

  !do ix = 0, L-1
    !do iy = 0, L-1
      !do iz = 0, L-1
        !site = get_site_from_matrix(D, (/ix,iy,iz/))
        !print *, ix, iy, iz, site, F(1, ix, iy, iz, 0)
      !enddo
    !enddo
  !enddo

CONTAINS

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
              site = get_site_from_matrix(3, (/ix, iy, iz/))
              Mat(it, ix, iy, iz, il) = Mat1D(it, site, il)
            enddo
          enddo
        enddo
      enddo
    enddo
    return
  END SUBROUTINE decompose_matrix

  integer function get_site_from_matrix(dims, matrix)
    implicit none
    integer :: dims
    integer, dimension(dims) :: matrix
    integer :: i

    get_site_from_matrix = 0
    do i = 1, dims
      get_site_from_matrix = get_site_from_matrix + matrix(i)*L**(i-1)
    enddo
    return
  END FUNCTION get_site_from_matrix

  function get_matrix_from_site(dims, site)
    implicit none
    integer, intent(in) :: dims, site
    integer, dimension(dims) :: get_matrix_from_site
    integer :: i, tmp

    tmp = site
    do i = dims, 1, -1
      get_matrix_from_site(i) =  tmp/L**(i-1)
      tmp = tmp - get_matrix_from_site(i)*L**(i-1)
    enddo
    return
  END FUNCTION get_matrix_from_site

END PROGRAM
