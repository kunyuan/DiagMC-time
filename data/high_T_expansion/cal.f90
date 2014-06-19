PROGRAM Main
  implicit none
  integer, parameter :: Order=8
  double precision :: a(0:Order)
  double precision :: beta
  double precision :: stChi
  integer :: i, tmpOrder

  print *, "Input beta: "
  read *, beta
  a(0) = 0.25
  a(1) = 0.375
  a(2) = 0.4375
  a(3) = 0.46875
  a(4) = 0.483203
  a(5) = 0.501367
  a(6) = 0.512950
  a(7) = 0.512135
  a(8) = 0.505525

  do tmpOrder = 0, Order
    stChi = 0.d0
    do i = 0, tmpOrder
      stChi = stChi + a(i)*beta**i
    enddo
    print *, "Order: ", tmpOrder, "Staggered Susceptibility: ", stChi
  enddo


END PROGRAM
