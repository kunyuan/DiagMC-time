
!!=======================================================================
!!====================== CONFIGURATION CHECK ============================
!!=======================================================================
SUBROUTINE check_config
  implicit none

  call check_topo
  call check_stat
  call check_irreducibility
  call check_k_conserve
  call check_type
  call check_time
  call check_site
  call check_weight
  return
END SUBROUTINE check_config



SUBROUTINE check_topo
  implicit none
  integer :: i, j, k
  integer :: nextLn, nLn
  
  if(NGLn/=2*(Order+1) .or. NWLn/=Order+1 .or. NVertex/=2*(Order+1)) then
    write(*, *) "================================================="
    write(*, *) "Oops, check_topo found a bug!"
    write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
    write(*, *) "Diagram order", Order
    write(*, *) "number of glines", NGLn, "number of wlines", NWLn, &
      & "number of gamma", NVertex
    write(*, *) "================================================="
    call print_config
    stop
  endif

  do i = 1, NGLn
    if(LnValue2Key(GLnKey2Value(i))/=i) then
      write(*, *) "================================================="
      write(*, *) "Oops, check_topo found a bug!"
      write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      write(*, *) "gline's location is wrong!"
      write(*, *) "real location", i, "gline's number", GLnKey2Value(i), "gline's &
        & location", LnValue2Key(GLnKey2Value(i))
      write(*, *) "================================================="
      call print_config
      stop
    endif
  enddo

  if(IsWormPresent .eqv. .false.) then
    nLn = 1
    nextLn = NeighVertex(2, NeighLn(2, GLnKey2Value(1)))
    do while(nextLn/=GLnKey2Value(1))
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
    if(LnValue2Key(WLnKey2Value(i))/=i) then
      write(*, *) "================================================="
      write(*, *) "Oops, check_topo found a bug!"
      write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      write(*, *) "wline's location is wrong!"
      write(*, *) "real location", i, "wline's number", WLnKey2Value(i), "wline's &
        & location", LnValue2Key(WLnKey2Value(i))
      write(*, *) "================================================="
      call print_config
      stop
    endif

    if(IsWormPresent .eqv. .false.) then
      if(NeighLn(1, WLnKey2Value(i))==NeighLn(2, WLnKey2Value(i))) then
        write(*, *) "================================================="
        write(*, *) "Oops, check_topo found a bug!"
        write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
        write(*, *) "wline's topology is wrong!"
        write(*, *) "wline's number", WLnKey2Value(i), "wline's left Gam", NeighLn(1, &
          &  WLnKey2Value(i)), "right Gam", NeighLn(2, WLnKey2Value(i))
        write(*, *) "================================================="
        call print_config
        stop
      endif
    endif
  enddo

  do i = 1, NVertex
    if(VertexValue2Key(VertexKey2Value(i))/=i) then
      write(*, *) "================================================="
      write(*, *) "Oops, check_topo found a bug!"
      write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      write(*, *) "gamma's location is wrong!"
      write(*, *) "real location", i, "gamma's number", VertexKey2Value(i), "gamma's &
        & location", VertexValue2Key(VertexKey2Value(i))
      write(*, *) "================================================="
      call print_config
      stop
    endif
  enddo

END SUBROUTINE check_topo



SUBROUTINE check_stat
  implicit none
  integer :: i, ikey, j, k
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

  do ikey = 1, NVertex
    i = VertexKey2Value(ikey)

    stat = 0
    if(MeasureGam==i)         stat = stat+1
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

  do ikey = 1, NGLn
    i = GLnKey2Value(ikey)
    stat = 0
    if(NeighLn(1,i)==MeasureGam .or. NeighLn(2,i)==MeasureGam) stat = stat+1
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

  do ikey = 1, NWLn
    i = WLnKey2Value(ikey)
    stat = 0
    if(NeighLn(1,i)==MeasureGam .or. NeighLn(2,i)==MeasureGam) stat = stat+1
    if(IsWormPresent) then
      if(NeighLn(1,i)==Ira .or. NeighLn(2,i)==Ira .or. NeighLn(1,i)==Masha &
        & .or. NeighLn(2,i)==Masha)                            stat = stat+2
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
  
SUBROUTINE check_time
  implicit none
  integer :: ikey, i

  do ikey = 1, NVertex
    i = VertexKey2Value(ikey)
    if(IsDeltaVertex(i)==1) then
      if(T1Vertex(i)/=T2Vertex(i) .or. T1Vertex(i)/=T3Vertex(i) .or. T2Vertex(i) &
        & /=T3Vertex(i)) then
        write(*, *) "================================================="
        write(*, *) "Oops, check_time found a bug!"
        write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
        write(*, *) "delta Gamma is wrong!"
        write(*, *) "gamma's number", i, "time", T1Vertex(i), T2Vertex(i), T3Vertex(i)
        write(*, *) "================================================="
        call print_config
        stop
      endif
    endif
  enddo

  do ikey = 1, NWLn
    i = WLnKey2Value(ikey)
    if(IsDeltaLn(i)==1) then
      if(T3Vertex(NeighLn(1, i))/=T3Vertex(NeighLn(2, i))) then
        write(*, *) "================================================="
        write(*, *) "Oops, check_time found a bug!"
        write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
        write(*, *) "delta W is wrong!"
        write(*, *) "W's number", i, "time", T3Vertex(NeighLn(1,i)),T3Vertex(NeighLn(2,i))
        write(*, *) "================================================="
        call print_config
        stop
      endif
    endif
  enddo
  return
END SUBROUTINE check_time

SUBROUTINE check_site
  implicit none
  integer :: ikey, i, j, k

  do ikey = 1, NVertex
    i = VertexKey2Value(ikey)
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
  



SUBROUTINE check_k_conserve
  implicit none
  integer :: ikey
  integer :: i, j, k
  integer :: omega
  
  do ikey = 1, NVertex
    i = VertexKey2Value(ikey)
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
      write(*, *) "Oops, check_k_conserve found a bug!"
      write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate
      write(*, *) "k on gamma is not conserved!"
      write(*, *) "gamma's number", i, "k", k
      write(*, *) "================================================="
      call print_config
      stop
    endif
  enddo
  return
END SUBROUTINE check_k_conserve


SUBROUTINE check_type
  implicit none
  integer :: i, ikey
  integer :: Gam1, Gam2, G1, G2, G3, G4
  integer :: flag
  integer :: sum1

  do ikey = 1, NWLn
    i = WLnKey2Value(ikey)

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
      if(TypeVertex(Gam1)/=TypeSp2Gam(TypeLn(G1),TypeLn(G2),&
        & SpInVertex(1, Gam1),SpInVertex(2, Gam1))) flag=1
      if(TypeVertex(Gam2)/=TypeSp2Gam(TypeLn(G3),TypeLn(G4),&
        & SpInVertex(1, Gam2),SpInVertex(2, Gam2))) flag=2
      if(TypeLn(i)/=TypeGam2W(TypeVertex(Gam1), TypeVertex(Gam2)))  flag=3
    endif
    if(sum1+2*(TypeLn(G1)+TypeLn(G3)-TypeLn(G2)-TypeLn(G4))/=0) flag=4
    if(sum1+2*(SpInVertex(1, Gam1)+SpInVertex(1, Gam2)-SpInVertex(2, Gam1)-SpInVertex(2, Gam2))/=0) flag=5

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
      write(*, *) sum1, SpInVertex(:, Gam1), SpInVertex(:, Gam2)
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
      Gi = GLnKey2Value(i)
      do j = i+1, NGLn
        Gj = GLnKey2Value(j)
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
      Wi = WLnKey2Value(i)
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
        Wj = WLnKey2Value(j)
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


  if(CheckGam) then
    do i = 1, NGLn
      Gi = GLnKey2Value(i)
      do j = i+1, NGLn
        Gj = GLnKey2Value(j)
        if(NeighLn(1,Gi)==NeighLn(2,Gj) .or. NeighLn(2,Gi)==NeighLn(1,Gj)) cycle
        do k = 1, NWLn
          Wk = WLnKey2Value(k)
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
  double precision :: tau1, tau2, tau
  complex*16 :: weight
  complex*16 :: wln(MxNLn), wgam(MxNVertex)

  weight = 1.d0
  do i = 1, MxNLn
    if(StatusLn(i)<0)  cycle
    if(KindLn(i)==1) then
      tau = T2Vertex(NeighLn(2, i))-T1Vertex(NeighLn(1, i))
      wln(i) = weight_line(StatusLn(i),IsDeltaLn(i),1,0,0, tau, TypeLn(i))
    else
      Gam1 = NeighLn(1,i);       Gam2 = NeighLn(2,i)
      tau = T3Vertex(NeighLn(2, i))-T3Vertex(NeighLn(1, i))
      wln(i) = weight_line(StatusLn(i),IsDeltaLn(i),2,WXVertex(Gam1)-WXVertex(Gam2), &
        & WYVertex(Gam1)-WYVertex(Gam2), tau, TypeLn(i))
    endif
    weight = weight *wln(i)
  enddo

  do i = 1, MxNVertex
    if(StatusVertex(i)<0)  cycle
    G1 = NeighVertex(1, i);           G2 = NeighVertex(2, i)
    W = NeighVertex(3, i)
    tau1 = T3Vertex(i)-T2Vertex(i)
    tau2 = T1Vertex(i)-T3Vertex(i)
    wgam(i) = weight_vertex(StatusVertex(i), IsDeltaVertex(i), GXVertex(i)-WXVertex(i), &
      & GYVertex(i)-WYVertex(i), tau1, tau2, TypeVertex(i))
    weight = weight *wgam(i)
  enddo

  weight = weight *CoefOfWeight(Order) *(1.d0/Beta)**Order *SignFermiLoop

  if(IsWormPresent) weight = weight*WeightWorm

  if(abs(Phase*WeightCurrent - weight)>1.d-5) then
    write(*, *) "================================================="
    write(*, *) "Oops, check_weight found a bug!"
    write(*, *) "IsWormPresent", IsWormPresent, "update number", iupdate

    write(*, *) "real weight", weight
    write(*, *) "current weight", Phase*WeightCurrent
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


