
!!=======================================================================
!!====================== CONFIGURATION CHECK ============================
!!=======================================================================
SUBROUTINE check_config
  implicit none

  call check_topo
  call check_stat
  call check_irreducibility
  call check_omega
  call check_type
  call check_site
  call check_weight
  return
END SUBROUTINE check_config



SUBROUTINE check_topo
  implicit none
  integer :: i, j, k
  integer :: nextLn, nLn
  
  if(NGLn/=2*(Order+1) .or. NWLn/=Order+1 .or. NGam/=2*(Order+1)) then
    write(*, *) "================================================="
    write(*, *) "Oops, check_topo found a bug!"
    write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
    write(*, *) "Diagram order", Order
    write(*, *) "number of glines", NGLn, "number of wlines", NWLn, &
      & "number of gamma", NGam
    write(*, *) "================================================="
    call print_config
    stop
  endif

  do i = 1, NGLn
    if(List4Ln(Ln4GList(i))/=i) then
      write(*, *) "================================================="
      write(*, *) "Oops, check_topo found a bug!"
      write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      write(*, *) "gline's location is wrong!"
      write(*, *) "real location", i, "gline's number", Ln4GList(i), "gline's &
        & location", List4Ln(Ln4GList(i))
      write(*, *) "================================================="
      call print_config
      stop
    endif
  enddo

  if(IsWormPresent .eqv. .false.) then
    nLn = 1
    nextLn = NeighVertex(2, NeighLn(2, Ln4GList(1)))
    do while(nextLn/=Ln4GList(1))
      nextLn = NeighVertex(2, NeighLn(2, nextLn))
      nLn = nLn + 1
      if(nLn>NGLn) then
        write(*, *) "================================================="
        write(*, *) "Oops, check_topo found a bug!"
        write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
        write(*, *) "The Fermi loop is wrong!"
        write(*, *) "The gline number in a loop:", nLn
        write(*, *) "================================================="
        call print_config
        stop
      endif
    enddo

    !if(nLn/=NGLn) then
      !write(*, *) "================================================="
      !write(*, *) "Oops, check_topo found a bug!"
      !write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      !write(*, *) "The Fermi loop is wrong!"
      !write(*, *) "The gline number in a loop:", nLn
      !write(*, *) "================================================="
      !call print_config
      !stop
    !endif
  endif

  do i = 1, NWLn
    if(List4Ln(Ln4WList(i))/=i) then
      write(*, *) "================================================="
      write(*, *) "Oops, check_topo found a bug!"
      write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      write(*, *) "wline's location is wrong!"
      write(*, *) "real location", i, "wline's number", Ln4WList(i), "wline's &
        & location", List4Ln(Ln4WList(i))
      write(*, *) "================================================="
      call print_config
      stop
    endif

    if(IsWormPresent .eqv. .false.) then
      if(NeighLn(1, Ln4WList(i))==NeighLn(2, Ln4WList(i))) then
        write(*, *) "================================================="
        write(*, *) "Oops, check_topo found a bug!"
        write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
        write(*, *) "wline's topology is wrong!"
        write(*, *) "wline's number", Ln4WList(i), "wline's left Gam", NeighLn(1, &
          &  Ln4WList(i)), "right Gam", NeighLn(2, Ln4WList(i))
        write(*, *) "================================================="
        call print_config
        stop
      endif
    endif
  enddo

  do i = 1, NGam
    if(List4Vertex(Vertex4GamList(i))/=i) then
      write(*, *) "================================================="
      write(*, *) "Oops, check_topo found a bug!"
      write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      write(*, *) "gamma's location is wrong!"
      write(*, *) "real location", i, "gamma's number", Vertex4GamList(i), "gamma's &
        & location", List4Ln(Vertex4GamList(i))
      write(*, *) "================================================="
      call print_config
      stop
    endif
  enddo

END SUBROUTINE check_topo



SUBROUTINE check_stat
  implicit none
  integer :: i, j, k
  integer :: stat

  if(IsWormPresent) then
    if(StatusVertex(Ira)<=1)  then
      write(*, *) "================================================="
      write(*, *) "Oops, check_stat found a bug!"
      write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      write(*, *) "Ira's status is wrong!"
      write(*, *) "Ira's number", Ira, "Ira's status",StatusVertex(Ira)
      write(*, *) "================================================="
      call print_config
      stop
    endif

    if(StatusVertex(Masha)<=1)  then
      write(*, *) "================================================="
      write(*, *) "Oops, check_stat found a bug!"
      write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      write(*, *) "Masha's status is wrong!"
      write(*, *) "Masha's number", Masha, "Masha's status",StatusVertex(Masha)
      write(*, *) "================================================="
      call print_config
      stop
    endif
  endif

  do i = 1, MxNVertex
    if(StatusVertex(i)==-1)  cycle
    stat = 0
    if(MeasGamma==i)         stat = stat+1
    if(IsWormPresent) then
      if(Ira==i .or. Masha==i) stat = stat+2
    endif
    if(StatusVertex(i)/=stat) then

      write(*, *) "================================================="
      write(*, *) "Oops, check_stat found a bug!"
      write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      write(*, *) "Gamma's status is wrong!"
      write(*, *) "Gamma's number", i, "real status", stat, "Gamma's status",StatusVertex(i)
      write(*, *) "================================================="

      call print_config
      stop
    endif
  enddo

  do i = 1, MxNLn
    if(StatusLn(i)==-1)  cycle
    stat = 0
    if(NeighLn(1,i)==MeasGamma .or. NeighLn(2,i)==MeasGamma) stat = stat+1
    if(KindLn(i)==2) then
      if(IsWormPresent) then
        if(NeighLn(1,i)==Ira .or. NeighLn(2,i)==Ira .or. NeighLn(1,i)==Masha &
          & .or. NeighLn(2,i)==Masha)                            stat = stat+2
      endif
    endif
    if(StatusLn(i)/=stat) then
      write(*, *) "================================================="
      write(*, *) "Oops, check_stat found a bug!"
      write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      write(*, *) "line's status is wrong!"
      write(*, *) "line's number", i, "real status", stat, "line's status",StatusLn(i)
      write(*, *) "================================================="

      call print_config
      stop
    endif
  enddo
END SUBROUTINE check_stat
  

SUBROUTINE check_site
  implicit none
  integer :: i, j, k
  integer :: omega
  do i = 1, MxNVertex
    if(StatusVertex(i)<0)  cycle
    if(GXVertex(i)<0 .or. GXVertex(i)>Lx-1) then
      write(*, *) "================================================="
      write(*, *) "Oops, check_site found a bug!"
      write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      write(*, *) "GX is wrong!"
      write(*, *) "gamma's number", i, "GX", GXVertex(i)
      write(*, *) "================================================="
      call print_config
      stop
    endif
    if(GYVertex(i)<0 .or. GYVertex(i)>Ly-1) then
      write(*, *) "================================================="
      write(*, *) "Oops, check_site found a bug!"
      write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      write(*, *) "GY is wrong!"
      write(*, *) "gamma's number", i, "GY", GYVertex(i)
      write(*, *) "================================================="
      call print_config
      stop
    endif
    if(WXVertex(i)<0 .or. WXVertex(i)>Lx-1) then
      write(*, *) "================================================="
      write(*, *) "Oops, check_site found a bug!"
      write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      write(*, *) "WX is wrong!"
      write(*, *) "gamma's number", i, "WX", WXVertex(i)
      write(*, *) "================================================="
      call print_config
      stop
    endif
    if(WYVertex(i)<0 .or. WYVertex(i)>Ly-1) then
      write(*, *) "================================================="
      write(*, *) "Oops, check_site found a bug!"
      write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      write(*, *) "WY is wrong!"
      write(*, *) "gamma's number", i, "WY", WYVertex(i)
      write(*, *) "================================================="
      call print_config
      stop
    endif
  enddo
END SUBROUTINE check_site
  



SUBROUTINE check_omega
  implicit none
  integer :: i, j, k
  integer :: omega
  
  do i = 1, MxNVertex
    if(StatusVertex(i)==-1) cycle
    if(DirecVertex(i)==1) then
      omega = OmegaLn(NeighVertex(1,i))-OmegaLn(NeighVertex(2,i)) &
        & -OmegaLn(NeighVertex(3,i))
    else
      omega = OmegaLn(NeighVertex(1,i))-OmegaLn(NeighVertex(2,i)) &
        & +OmegaLn(NeighVertex(3,i))
    endif
    if(IsWormPresent .eqv. .true.) then
      if(i==Ira)   omega = omega + OmegaMasha
      if(i==Masha) omega = omega - OmegaMasha
    endif
    if(omega/=0) then
      write(*, *) "================================================="
      write(*, *) "Oops, check_omega found a bug!"
      write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      write(*, *) "Omega on gamma is not conserved!"
      write(*, *) "gamma's number", i, "omega", omega
      write(*, *) "================================================="
      call print_config
      stop
    endif
  enddo

  do i = 1, MxNVertex
    if(StatusVertex(i)==-1) cycle
    if(DirecVertex(i)==1) then
      k = add_k(add_k(kLn(NeighVertex(1,i)), -kLn(NeighVertex(2,i))), -kLn(NeighVertex(3,i)))
    else
      k = add_k(add_k(kLn(NeighVertex(1,i)), -kLn(NeighVertex(2,i))), kLn(NeighVertex(3,i)))
    endif
    if(IsWormPresent .eqv. .true.) then
      if(i==Ira)   k = k + kMasha
      if(i==Masha) k = k - kMasha
    endif
    if(k/=0) then
      write(*, *) "================================================="
      write(*, *) "Oops, check_omega found a bug!"
      write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      write(*, *) "k on gamma is not conserved!"
      write(*, *) "gamma's number", i, "k", k
      write(*, *) "================================================="
      call print_config
      stop
    endif
  enddo
  return
END SUBROUTINE check_omega


SUBROUTINE check_type
  implicit none
  integer :: i
  integer :: Gam1, Gam2, G1, G2, G3, G4
  integer :: flag
  integer :: sum1

  !========== Only for the first round in self-consistent loop =============
  !do i = 1, MxNVertex
    !if(StatusVertex(i)==-1) cycle
    !if(TypeVertex(i)==3 .or. TypeVertex(i)==4)  then
      !write(*, *) "================================================="
      !write(*, *) "Oops, check_type found a bug!"
      !write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      !write(*, *) "In the first round, the weight of Gamma with type 3 and 4 is not zero!"
      !write(*, *) "Gamma's number", i, "type of Gamma", TypeVertex(i), "status", &
        !& StatusVertex(i), "weight", WeightVertex(i)
      !write(*, *) "================================================="
      !call print_config
      !stop
    !endif
    !if(TypeVertexIn(i)>2 .or. TypeVertexIn(i)<1) then
      !write(*, *) "================================================="
      !write(*, *) "Oops, check_type found a bug!"
      !write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      !write(*, *) "The type of Gamma inside spin is not up nor down!"
      !write(*, *) "Gamma's number", i, "type of Gamma inside spin in", TypeVertexIn(i)
      !write(*, *) "================================================="
      !call print_config
      !stop
    !endif
    !if(TypeVertexOut(i)>2 .or. TypeVertexOut(i)<1) then
      !write(*, *) "================================================="
      !write(*, *) "Oops, check_type found a bug!"
      !write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      !write(*, *) "The type of Gamma inside spin is not up nor down!"
      !write(*, *) "Gamma's number", i, "type of Gamma inside spin out", TypeVertexOut(i)
      !write(*, *) "================================================="
      !call print_config
      !stop
    !endif
  !enddo
  !==========================================================================

  do i = 1, MxNLn
    if(StatusLn(i)==-1) cycle
    if(KindLn(i)==1)    cycle

    Gam1 = NeighLn(1, i)
    Gam2 = NeighLn(2, i)

    sum1 = 0
    if(IsWormPresent) then
      if(Gam1==Ira) then
        sum1 = sum1 -SpinMasha
      else if(Gam1==Masha) then
        sum1 = sum1 + SpinMasha
      endif
      if(Gam2==Ira) then
        sum1 = sum1 -SpinMasha
      else if(Gam2==Masha) then
        sum1 = sum1 + SpinMasha
      endif
    endif

    G1 = NeighVertex(1, Gam1);      G2 = NeighVertex(2, Gam1)
    G3 = NeighVertex(1, Gam2);      G4 = NeighVertex(2, Gam2)

    flag = 0
    if(StatusLn(i)<=1) then
      if(TypeVertex(Gam1)/=TypeGW2Gam(TypeLn(G1),TypeLn(G2),&
        & TypeVertexIn(Gam1),TypeVertexOut(Gam1))) flag=1
      if(TypeVertex(Gam2)/=TypeGW2Gam(TypeLn(G3),TypeLn(G4),&
        & TypeVertexIn(Gam2),TypeVertexOut(Gam2))) flag=2
      if(TypeLn(i)/=TypeGam2W(TypeVertex(Gam1), TypeVertex(Gam2)))  flag=3
    endif
    if(sum1+2*(TypeLn(G1)+TypeLn(G3)-TypeLn(G2)-TypeLn(G4))/=0) flag=4
    if(sum1+2*(TypeVertexIn(Gam1)+TypeVertexIn(Gam2)-TypeVertexOut(Gam1)-TypeVertexOut(Gam2))/=0) flag=5

    if(flag==1 .or. flag==2) then
      write(*, *) "================================================="
      write(*, *) "Oops, check_type found a bug!"
      write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      write(*, *) "The type of Gamma is wrong!"
      write(*, *) "wline's number", i, "Gamma", Gam1, Gam2
      write(*, *) "================================================="
      call print_config
      stop
    else if(flag==3) then
      write(*, *) "================================================="
      write(*, *) "Oops, check_type found a bug!"
      write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      write(*, *) "The type of wline is wrong!"
      write(*, *) "wline's number", i, "Gamma", Gam1, Gam2
      write(*, *) TypeLn(i), TypeVertex(Gam1), TypeVertex(Gam2)
      write(*, *) "================================================="
      call print_config
      stop
    else if(flag==4) then
      write(*, *) "================================================="
      write(*, *) "Oops, check_type found a bug!"
      write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      write(*, *) "The type of glines are wrong!"
      write(*, *) "wline's number", i, "Gamma", Gam1, Gam2
      write(*, *) sum1, TypeLn(G1), TypeLn(G3), TypeLn(G2), TypeLn(G4)
      write(*, *) "================================================="
      call print_config
      stop
    else if(flag==5) then
      write(*, *) "================================================="
      write(*, *) "Oops, check_type found a bug!"
      write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      write(*, *) "The type of gamma inside lines are wrong!"
      write(*, *) "wline's number", i, "Gamma", Gam1, Gam2
      write(*, *) sum1, TypeVertexIn(Gam1), TypeVertexIn(Gam2), TypeVertexOut(Gam2), &
        & TypeVertexOut(Gam2)
      write(*, *) "================================================="
      call print_config
      stop
    endif
  enddo

  return
END SUBROUTINE check_type


SUBROUTINE check_irreducibility
  implicit none
  integer :: i, j, k
  integer :: Gi, Wi
  integer :: Gj, Wj
  integer :: Gk, Wk
  
  if(CheckG) then
    do i = 1, NGLn
      Gi = Ln4GList(i)
      do j = i+1, NGLn
        Gj = Ln4GList(j)
        if(kLn(Gj)==kLn(Gi)) then
          write(*, *) "================================================="
          write(*, *) "Oops, check_irreducibility found a bug!"
          write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
          write(*, *) "Gline is not irreducible!"
          write(*, *) "Gline's number", Gi, "k", kLn(Gi)
          write(*, *) "================================================="
          call print_config
          stop
        endif
      enddo
    enddo
  endif

  if(CheckW) then
    do i = 1, NWLn
      Wi = Ln4WList(i)
      if(kLn(Wi)==0) then
        write(*, *) "================================================="
        write(*, *) "Oops, check_irreducibility found a bug!"
        write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
        write(*, *) "Wline is not irreducible!"
        write(*, *) "Wline's number", Wi, "k", kLn(Wi)
        write(*, *) "================================================="
        call print_config
        stop
      endif

      do j = i+1, NWLn
        Wj = Ln4WList(j)
        if(abs(kLn(Wj))==abs(kLn(Wi))) then
          write(*, *) "================================================="
          write(*, *) "Oops, check_irreducibility found a bug!"
          write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
          write(*, *) "Wline is not irreducible!"
          write(*, *) "Wline's number", Wi, "k", kLn(Wi)
          write(*, *) "================================================="
          call print_config
          stop
        endif
      enddo
    enddo
  endif


  if(CheckGamma) then
    do i = 1, NGLn
      Gi = Ln4GList(i)
      do j = i+1, NGLn
        Gj = Ln4GList(j)
        if(NeighLn(1,Gi)==NeighLn(2,Gj) .or. NeighLn(2,Gi)==NeighLn(1,Gj)) cycle
        do k = 1, NWLn
          Wk = Ln4WList(k)
          if(NeighLn(1,Wk)==NeighLn(1,Gi) .or. NeighLn(1,Wk)==NeighLn(1,Gj)) cycle
          if(NeighLn(1,Wk)==NeighLn(2,Gi) .or. NeighLn(1,Wk)==NeighLn(2,Gj)) cycle
          if(NeighLn(2,Wk)==NeighLn(1,Gi) .or. NeighLn(2,Wk)==NeighLn(1,Gj)) cycle
          if(NeighLn(2,Wk)==NeighLn(2,Gi) .or. NeighLn(2,Wk)==NeighLn(2,Gj)) cycle
          if(abs(add_k(kLn(Gi), -kLn(Gj)))==abs(kLn(Wk))) then
            write(*, *) "================================================="
            write(*, *) "Oops, check_irreducibility found a bug!"
            write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
            write(*, *) "Gamma is not irreducible!"
            write(*, *) "Gline's number", Gi, Gj, "Wline", Wk, "k", kLn(Gi), kLn(Gj), kLn(Wk)
            write(*, *) "================================================="
            call print_config
            stop
          endif
        enddo
      enddo
    enddo
  endif

  return
END SUBROUTINE check_irreducibility



SUBROUTINE check_weight
  implicit none
  integer :: i, Gam1, Gam2, G, G1, G2, W
  double precision :: weight
  double precision :: wln(MxNLn), wgam(MxNVertex)

  weight = 1.d0
  do i = 1, MxNLn
    if(StatusLn(i)<0)  cycle
    if(KindLn(i)==1) then
      wln(i) = weight_line(StatusLn(i),1,0,0,OmegaLn(i),TypeLn(i))
    else
      Gam1 = NeighLn(1,i);       Gam2 = NeighLn(2,i)
      wln(i) = weight_line(StatusLn(i),2,diff_x(WXVertex(Gam1),WXVertex(Gam2)), &
        & diff_y(WYVertex(Gam1),WYVertex(Gam2)), OmegaLn(i), TypeLn(i))
    endif
    weight = weight *wln(i)
  enddo

  do i = 1, MxNVertex
    if(StatusVertex(i)<0)  cycle
    G1 = NeighVertex(1, i);           G2 = NeighVertex(2, i)
    W = NeighVertex(3, i)
    wgam(i) = weight_vertex(Order, StatusVertex(i), diff_x(GXVertex(i),WXVertex(i)), &
      & diff_y(GYVertex(i),WYVertex(i)), OmegaLn(G1), OmegaLn(G2), TypeVertex(i))
    weight = weight *wgam(i)
  enddo

  weight = weight*CoefOfWeight(Order)*(1.d0/Beta)**Order *(-1.d0)**NFermiLoop
  if(IsWormPresent) weight = weight*WeightWorm

  if(abs(WeightCurrent-weight)>1.d-5) then
    write(*, *) "================================================="
    write(*, *) "Oops, check_weight found a bug!"
    write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate

    write(*, *) "real weight", weight, "current weight", WeightCurrent
    write(*, *) Order, Beta

    do i = 1, MxNLn
      if(StatusLn(i)<0)  cycle
      write(*, *) i, wln(i), WeightLn(i)
    enddo

    do i = 1, MxNVertex
      if(StatusVertex(i)<0)  cycle
      write(*, *) i, wgam(i), WeightVertex(i)
    enddo
    write(*, *) "================================================="

    call print_config
    stop
  endif
  return
END SUBROUTINE check_weight

!!=======================================================================
!!=======================================================================
!!=======================================================================


