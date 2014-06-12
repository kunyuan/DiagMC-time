
SUBROUTINE test_subroutine
    implicit none
    !integer :: isamp
    !!======== test x,y distribution =========================
    !integer :: i,nr(2),cr(2),dr(2),N,x,y
    !double precision :: hist(2,0:MxL(1)-1),weight
    !call initialize_markov
    !print *,"Testing..."
    !hist(:,:)=0.d0
    !N=100000
    !weight=0.0
    !cr(:)=0
    !do i=1,N
      !call generate_xy(cr,nr,dr,weight,.true.)
      !hist(1,nr(1))=hist(1,nr(1))+1
      !hist(2,nr(2))=hist(2,nr(2))+1
    !enddo
    !open(11,file="testx.dat")
    !write(11,*) "X:",L(1),logL(1)
    !do x=0,L(1)
      !write(11,*) x, hist(1,x)/N, SpatialWeight(1,x)
    !enddo
    !close(11)
    !open(11,file="testy.dat")
    !write(11,*) "Y:",L(2),logL(2)
    !do y=0,L(2)
      !write(11,*) y, hist(2,y)/N, SpatialWeight(2,y)
    !enddo
    !close(11)

    !========  test drawing subroutine =====================
    !call initialize_markov
    !call print_config

    !======== analytic_integration =========================
    !call read_GWGamma
    !call calculate_Gam1
    !call output_Gam1

    return
END SUBROUTINE
