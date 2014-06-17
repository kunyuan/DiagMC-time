Program main
  implicit none
  !integer, parameter :: D=1
  !integer, parameter :: D=2
  integer, parameter :: D=3
  integer, parameter :: L=8
  integer :: i,j,k
  integer :: site
  integer, dimension(D) :: mt, mt_new

  !============ test example for 1D ================
  !do i = 1, L
    !mt = (/i/)
    !site = get_site_from_matrix(D, mt)
    !if(site/=mt(1)) then
      !print *, '1d transfer error: get_site_from_matrix'
    !endif

    !mt_new = get_matrix_from_site(D, site)
    !if(mt_new(1)/=mt(1) .or. mt_new(1)/=i) then
      !print *, '1d transfer error: get_matrix_from_site'
    !endif
  !enddo

  !============ test example for 2D ================
  !do i = 0, L-1
    !do j = 0, L-1
      !mt = (/i, j/)
      !site = get_site_from_matrix(D, mt)
      !mt_new = get_matrix_from_site(D, site)
      !if(mt_new(1)/=i .or. mt_new(2)/=j) then
        !print *, '2d transfer error: get_matrix_from_site'
      !endif
    !enddo
  !enddo

  !============ test example for 3D ================
  do i = 0, L-1
    do j = 0, L-1
      do k = 0, L-1
        mt = (/i, j, k/)
        site = get_site_from_matrix(D, mt)
        mt_new = get_matrix_from_site(D, site)
        if(mt_new(1)/=i .or. mt_new(2)/=j .or. mt_new(3)/=k) then
          print *, '3d transfer error: get_matrix_from_site'
        endif
      enddo
    enddo
  enddo

  do i = 0, L**3-1
    mt = get_matrix_from_site(D, i)
    print *, i, ': ', mt(1), mt(2), mt(3)
  enddo

CONTAINS
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
