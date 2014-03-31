PROGRAM MAIN
  implicit none
  integer, dimension(2) :: dr
  double precision :: W

  dr(:) = 1
  W = weight_W(dr)
  write(*, *) W

contains

double precision function weight_W(r)
  implicit none
  integer, intent(in) :: r(2)

  weight_W = r(1)**2.d0 + r(2)**2.d0
  return
end function

END PROGRAM

