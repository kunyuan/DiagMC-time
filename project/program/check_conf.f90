
!!=======================================================================
!!====================== CONFIGURATION CHECK ============================
!!=======================================================================
SUBROUTINE check_config
  implicit none

  call LogFile%QuickLog("Checking the configuration...")
  call check_topo
  call check_stat
  call check_irreducibility
  call check_isdelta
  call check_k_conserve
  call check_type
  call check_time
  call check_site
  call check_weight
  call LogFile%QuickLog("Checking configuration is done!")
  return
END SUBROUTINE check_config


SUBROUTINE check_isdelta
  implicit none
  
  if(IsDeltaVertex(MeasureGam)/=1) then
    call LogFile%WriteStamp('e')
    call LogFile%WriteLine("Oops, check_isdelta found a bug!")
    call LogFile%WriteLine("IsWormPresent:"+str(IsWormPresent)+", imc:"+str(imc))
    call LogFile%WriteLine("Diagramorder:"+str(Order)+", MeasureGam:"+str(MeasureGam))
    stop
  endif

  if(IsDeltaLn(NeighVertex(3, MeasureGam))/=0) then
    call LogFile%WriteStamp('e')
    call LogFile%WriteLine("Oops, check_isdelta found a bug!")
    call LogFile%WriteLine("IsWormPresent:"+str(IsWormPresent)+", imc:"+str(imc))
    call LogFile%WriteLine("Diagramorder:"+str(Order)+", MeasureGam:"+str(NeighVertex(3,MeasureGam)))
    stop
  endif
END SUBROUTINE


SUBROUTINE check_topo
  implicit none
  integer :: i, j, k
  integer :: nextLn, nLn, curG
  logical :: flag(MxNGLn)
  
  if(NGLn/=2*(Order+1) .or. NWLn/=Order+1 .or. NVertex/=2*(Order+1)) then
    call LogFile%WriteStamp('e')
    call LogFile%WriteLine("Oops, check_topo found a bug!")
    call LogFile%WriteLine("IsWormPresent"+str(IsWormPresent)+", imc"+str(imc))
    call LogFile%WriteLine("Diagram order"+str(Order))
    call LogFile%WriteLine("number of glines"+str(NGLn)+",number of wlines"+str(NWLn)+ &
      & ",number of gamma"+str(NVertex))
    call print_config
    stop
  endif

  do i = 1, NGLn
    if(LnValue2Key(GLnKey2Value(i))/=i) then
      call LogFile%WriteStamp('e')
      call LogFile%WriteLine("Oops, check_topo found a bug!")
      call LogFile%WriteLine("IsWormPresent"+str(IsWormPresent)+", imc"+str(imc))
      call LogFile%WriteLine("gline's location is wrong!")
      call LogFile%WriteLine("real location"+str(i)+",gline's number"+str(GLnKey2Value(i))+",gline's &
        & location"+str(LnValue2Key(GLnKey2Value(i))))
      call print_config
      stop
    endif
  enddo

  flag(:)=.false.
  do i = 1, NGLn
    curG = GLnKey2Value(i)
    if(flag(curG)==.true.) cycle
    flag(curG)=.true.
    NextLn = NeighVertex(2, NeighLn(2, curG))
    nLn=1
    do while(NextLn/=curG)
      nLn=nLn+1
      if(flag(NextLn)==.true.) then
        !print *,curG, NextLn
        call LogFile%WriteStamp('e')
        call LogFile%WriteLine("Oops, check_topo found a bug!")
        call LogFile%WriteLine("IsWormPresent"+str(IsWormPresent)+", imc"+str(imc))
        call LogFile%WriteLine("The Fermi loop has mulitiple branches!")
        call LogFile%WriteLine("The gline number in a loop:"+str(NextLn))
        call print_config
        stop
      endif
      flag(NextLn)=.true.
      if(nLn>NGLn) then
        call LogFile%WriteStamp('e')
        call LogFile%WriteLine("Oops, check_topo found a bug!")
        call LogFile%WriteLine("IsWormPresent"+str(IsWormPresent)+", imc"+str(imc))
        call LogFile%WriteLine("The Fermi loop is wrong!")
        call LogFile%WriteLine("The gline number in a loop:"+str(nLn))
        call print_config
        stop
      endif
      NextLn = NeighVertex(2, NeighLn(2, NextLn))
    enddo
  enddo

  do i = 1, NWLn
    if(LnValue2Key(WLnKey2Value(i))/=i) then
      call LogFile%WriteStamp('e')
      call LogFile%WriteLine("Oops, check_topo found a bug!")
      call LogFile%WriteLine("IsWormPresent"+str(IsWormPresent)+", imc"+str(imc))
      call LogFile%WriteLine("wline's location is wrong!")
      call LogFile%WriteLine("real location"+str(i)+",wline's number"+str(WLnKey2Value(i))+",wline's &
        & location"+str(LnValue2Key(WLnKey2Value(i))))
      call print_config
      stop
    endif

    if(IsWormPresent .eqv. .false.) then
      if(NeighLn(1, WLnKey2Value(i))==NeighLn(2, WLnKey2Value(i))) then
        call LogFile%WriteStamp('e')
        call LogFile%WriteLine("Oops, check_topo found a bug!")
        call LogFile%WriteLine("IsWormPresent"+str(IsWormPresent)+", imc"+str(imc))
        call LogFile%WriteLine("wline's topology is wrong!")
        call LogFile%WriteLine("wline's number"+str(WLnKey2Value(i))+",wline's left Gam"+str(NeighLn(1, &
          &  WLnKey2Value(i)))+",right Gam"+str(NeighLn(2, WLnKey2Value(i))))
        call print_config
        stop
      endif
    endif
  enddo

  do i = 1, NVertex
    if(VertexValue2Key(VertexKey2Value(i))/=i) then
      call LogFile%WriteStamp('e')
      call LogFile%WriteLine("Oops, check_topo found a bug!")
      call LogFile%WriteLine("IsWormPresent"+str(IsWormPresent)+", imc"+str(imc))
      call LogFile%WriteLine("gamma's location is wrong!")
      call LogFile%WriteLine("real location"+str(i)+",gamma's number"+str(VertexKey2Value(i))+",gamma's &
        & location"+str(VertexValue2Key(VertexKey2Value(i))))
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
      call LogFile%WriteStamp('e')
      call LogFile%WriteLine("Oops, check_stat found a bug!")
      call LogFile%WriteLine("IsWormPresent"+str(IsWormPresent)+",imc"+str(imc))
      call LogFile%WriteLine("Ira's status is wrong!")
      call LogFile%WriteLine("Ira's number"+str(Ira)+",Ira's status"+str(StatusVertex(Ira)))
      call print_config
      stop
    endif

    if(StatusVertex(Masha)<=1)  then
      call LogFile%WriteStamp('e')
      call LogFile%WriteLine("Oops, check_stat found a bug!")
      call LogFile%WriteLine("IsWormPresent"+str(IsWormPresent)+",imc"+str(imc))
      call LogFile%WriteLine("Masha's status is wrong!")
      call LogFile%WriteLine("Masha's number"+str(Masha)+",Masha's status"+str(StatusVertex(Masha)))
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
      call LogFile%WriteStamp('e')
      call LogFile%WriteLine("Oops, check_stat found a bug!")
      call LogFile%WriteLine("IsWormPresent"+str(IsWormPresent)+",imc"+str(imc))
      call LogFile%WriteLine("Gamma's status is wrong!")
      call LogFile%WriteLine("Gamma's number"+str(i)+",real status"+str(stat)+",Gamma's status"+str(StatusVertex(i)))
      call print_config
      stop
    endif
  enddo

  do ikey = 1, NGLn
    i = GLnKey2Value(ikey)
    stat = 0
    if(NeighLn(1,i)==MeasureGam .or. NeighLn(2,i)==MeasureGam) stat = stat+1
    if(StatusLn(i)/=stat) then
      call LogFile%WriteStamp('e')
      call LogFile%WriteLine("Oops, check_stat found a bug!")
      call LogFile%WriteLine("IsWormPresent"+str(IsWormPresent)+",imc"+str(imc))
      call LogFile%WriteLine("line's status is wrong!")
      call LogFile%WriteLine("line's number"+str(i)+",real status"+str(stat)+",line's status"+str(StatusLn(i)))
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
      call LogFile%WriteStamp('e')
      call LogFile%WriteLine("Oops, check_stat found a bug!")
      call LogFile%WriteLine("IsWormPresent"+str(IsWormPresent)+",imc"+str(imc))
      call LogFile%WriteLine("line's status is wrong!")
      call LogFile%WriteLine("line's number"+str(i)+",real status"+str(stat)+",line's status"+str(StatusLn(i)))
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
      if(TVertex(1, i)/=TVertex(2, i) .or. TVertex(1, i)/=TVertex(3, i) .or. TVertex(2, i) &
        & /=TVertex(3, i)) then
        call LogFile%WriteStamp('e')
        call LogFile%WriteLine("Oops, check_time found a bug!")
        call LogFile%WriteLine("IsWormPresent"+str(IsWormPresent)+",imc"+str(imc))
        call LogFile%WriteLine("delta Gamma is wrong!")
        call LogFile%WriteLine("gamma's number"+str(i)+",time"+str(TVertex(1, i))+str(TVertex(2, i))+str(TVertex(3, i)))
        call print_config
        stop
      endif
    endif
  enddo

  do ikey = 1, NWLn
    i = WLnKey2Value(ikey)
    if(IsDeltaLn(i)==1) then
      if(TVertex(3, NeighLn(1, i))/=TVertex(3, NeighLn(2, i))) then
        call LogFile%WriteStamp('e')
        call LogFile%WriteLine("Oops, check_time found a bug!")
        call LogFile%WriteLine("IsWormPresent"+str(IsWormPresent)+",imc"+str(imc))
        call LogFile%WriteLine("delta W is wrong!")
        call LogFile%WriteLine("W's number"+str(i)+",time"+str(TVertex(3,NeighLn(1,i)))+str(TVertex(3,NeighLn(2,i))))
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
    if(GRVertex(1, i)<0 .or. GRVertex(1, i)>L(1)-1) then
      call LogFile%WriteStamp('e')
      call LogFile%WriteLine("Oops, check_site found a bug!")
      call LogFile%WriteLine("IsWormPresent"+str(IsWormPresent)+",imc"+str(imc))
      call LogFile%WriteLine("GX is wrong!")
      call LogFile%WriteLine("gamma's number"+str(i)+",GX"+str(GRVertex(1, i)))
      call print_config
      stop
    endif
    if(GRVertex(2, i)<0 .or. GRVertex(2, i)>L(2)-1) then
      call LogFile%WriteStamp('e')
      call LogFile%WriteLine("Oops, check_site found a bug!")
      call LogFile%WriteLine("IsWormPresent"+str(IsWormPresent)+",imc"+str(imc))
      call LogFile%WriteLine("GY is wrong!")
      call LogFile%WriteLine("gamma's number"+str(i)+",GY"+str(GRVertex(2, i)))
      call print_config
      stop
    endif
    if(WRVertex(1, i)<0 .or. WRVertex(1, i)>L(1)-1) then
      call LogFile%WriteStamp('e')
      call LogFile%WriteLine("Oops, check_site found a bug!")
      call LogFile%WriteLine("IsWormPresent"+str(IsWormPresent)+",imc"+str(imc))
      call LogFile%WriteLine("GY is wrong!")
      call LogFile%WriteLine("gamma's number"+str(i)+",WX"+str(WRVertex(1, i)))
      call print_config
      stop
    endif
    if(WRVertex(2, i)<0 .or. WRVertex(2, i)>L(2)-1) then
      call LogFile%WriteStamp('e')
      call LogFile%WriteLine("Oops, check_site found a bug!")
      call LogFile%WriteLine("IsWormPresent"+str(IsWormPresent)+",imc"+str(imc))
      call LogFile%WriteLine("GY is wrong!")
      call LogFile%WriteLine("gamma's number"+str(i)+",WY"+str(WRVertex(2, i)))
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
      if(i==Ira)   k = add_k(k, kMasha)
      if(i==Masha) k = add_k(k, -kMasha)
    endif
    if(imc==172436574 .and. i==1) then
      print *, kMasha, kLn(NeighVertex(1,i)), -kLn(NeighVertex(2,i)), kLn(NeighVertex(3,i))
    endif
    if(k/=0) then
      call LogFile%WriteStamp('e')
      call LogFile%WriteLine("Oops, check_k_conserve found a bug!")
      call LogFile%WriteLine("IsWormPresent"+str(IsWormPresent)+",imc"+str(imc))
      call LogFile%WriteLine("k on gamma is not conserved!")
      call LogFile%WriteLine("gamma's number"+str(i)+",k"+str(k))
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
      call LogFile%WriteStamp('e')
      call LogFile%WriteLine("Oops, check_type found a bug!")
      call LogFile%WriteLine("IsWormPresent"+str(IsWormPresent)+",imc"+str(imc))
      call LogFile%WriteLine("The type of Gamma is wrong!")
      call LogFile%WriteLine("wline's number"+str(i)+",Gamma"+str(Gam1)+str(Gam2))
      call print_config
      stop
    else if(flag==3) then
      call LogFile%WriteStamp('e')
      call LogFile%WriteLine("Oops, check_type found a bug!")
      call LogFile%WriteLine("IsWormPresent"+str(IsWormPresent)+",imc"+str(imc))
      call LogFile%WriteLine("The type of wline is wrong!")
      call LogFile%WriteLine("wline's number"+str(i)+",Gamma"+str(Gam1)+str(Gam2))
      call LogFile%WriteLine(str(TypeLn(i))+str(TypeVertex(Gam1))+str(TypeVertex(Gam2)))
      call print_config
      stop
    else if(flag==4) then
      call LogFile%WriteStamp('e')
      call LogFile%WriteLine("Oops, check_type found a bug!")
      call LogFile%WriteLine("IsWormPresent"+str(IsWormPresent)+",imc"+str(imc))
      call LogFile%WriteLine("The type of glines are wrong!")
      call LogFile%WriteLine("wline's number"+str(i)+",Gamma"+str(Gam1)+str(Gam2))
      call LogFile%WriteLine(str(sum1)+str(TypeLn(G1))+str(TypeLn(G3))+str(TypeLn(G2))+str(TypeLn(G4)))
      call print_config
      stop
    else if(flag==5) then
      call LogFile%WriteStamp('e')
      call LogFile%WriteLine("Oops, check_type found a bug!")
      call LogFile%WriteLine("IsWormPresent"+str(IsWormPresent)+",imc"+str(imc))
      call LogFile%WriteLine("The type of gamma inside lines are wrong!")
      call LogFile%WriteLine("wline's number"+str(i)+",Gamma"+str(Gam1)+str(Gam2))
      call LogFile%WriteLine(str(sum1)+str(SpInVertex(1, &
          &  Gam1))+str(SpInVertex(2,Gam1))+str(SpInVertex(1,Gam2))+str(SpInVertex(2,Gam2)))
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
  integer :: Gam1, Gam2
  
  if(CHECK_G) then
    do i = 1, NGLn
      Gi = GLnKey2Value(i)
      do j = i+1, NGLn
        Gj = GLnKey2Value(j)
        if(kLn(Gj)==kLn(Gi)) then
          call LogFile%WriteStamp('e')
          call LogFile%WriteLine('imc:'+str(imc))
          call LogFile%WriteLine("Oops, check_irreducibility found a bug!")
          call LogFile%WriteLine("IsWormPresent"+str(IsWormPresent)+",imc"+str(imc))
          call LogFile%WriteLine("Gline is not irreducible!")
          call LogFile%WriteLine("Gline's number"+str(Gi)+",k"+str(kLn(Gi)))
          call print_config
          stop
        endif
      enddo
    enddo
  endif

  if(CHECK_W) then
    do i = 1, NWLn
      Wi = WLnKey2Value(i)
      if(kLn(Wi)==0) then
        call LogFile%WriteStamp('e')
        call LogFile%WriteLine('imc:'+str(imc))
        call LogFile%WriteLine("Oops, check_irreducibility found a bug!")
        call LogFile%WriteLine("IsWormPresent"+str(IsWormPresent)+",imc"+str(imc))
        call LogFile%WriteLine("Wline is not irreducible!")
        call LogFile%WriteLine("Wline's number"+str(Wi)+",k"+str(kLn(Wi)))
        call print_config
        stop
      endif

      do j = i+1, NWLn
        Wj = WLnKey2Value(j)
        if(abs(kLn(Wj))==abs(kLn(Wi))) then
          call LogFile%WriteStamp('e')
          call LogFile%WriteLine('imc:'+str(imc))
          call LogFile%WriteLine("Oops, check_irreducibility found a bug!")
          call LogFile%WriteLine("IsWormPresent"+str(IsWormPresent)+",imc"+str(imc))
          call LogFile%WriteLine("Wline is not irreducible!")
          call LogFile%WriteLine("Wline's number"+str(Wi)+",k"+str(kLn(Wi)))
          call print_config
          stop
        endif
      enddo
    enddo
  endif


  if(CHECK_GAM) then
    do i = 1, NGLn
      Gi = GLnKey2Value(i)
      do j = i+1, NGLn
        Gj = GLnKey2Value(j)
        do k = 1, NWLn
          Wk = WLnKey2Value(k)
          if(abs(add_k(kLn(Gi), -kLn(Gj)))==abs(kLn(Wk))) then
            Gam1=NeighLn(1,Wk)
            Gam2=NeighLn(2,Wk)
            if(is_connected(Gam1, Gi, Gj)) cycle
            if(is_connected(Gam2, Gi, Gj)) cycle
            call LogFile%WriteStamp('e')
            call LogFile%WriteLine('imc:'+str(imc))
            call LogFile%WriteLine("Oops, check_irreducibility found a bug!")
            call LogFile%WriteLine("IsWormPresent"+str(IsWormPresent)+",iupdate"+str(iupdate))
            call LogFile%WriteLine("Gamma is not irreducible!")
            call LogFile%WriteLine("Gline's number"+str(Gi)+str(Gj)+",Wline"+str(Wk)+ &
                & ",k"+str(kLn(Gi))+str(kLn(Gj))+str(kLn(Wk)))
            call print_config
            stop
          endif
        enddo
      enddo
    enddo
  endif

  return
END SUBROUTINE check_irreducibility

logical FUNCTION is_connected(Gam, G1, G2)
  implicit none
  integer :: Gam, G1, G2
  if((NeighVertex(1,Gam)==G1 .and. NeighVertex(2,Gam)==G2) .or. &
   & (NeighVertex(2,Gam)==G1 .and. NeighVertex(1,Gam)==G2)) then
    if(IsWormPresent .and. (Gam==Ira .or. Gam==Masha)) then
      is_connected=.true.
    else
      is_connected=.true.
    endif
  else
    is_connected=.false.
  endif
  return
end FUNCTION



SUBROUTINE check_weight
  implicit none
  integer :: i, ikey, Gam1, Gam2, G, G1, G2, W
  double precision :: tau1, tau2, tau
  complex*16 :: weight
  complex*16 :: gln(NGLn), wln(NWLn), gam(NVertex)

  weight = (1.d0, 0.d0)
  do ikey = 1, NGLn
    i = GLnKey2Value(ikey)
    Gam1 = NeighLn(1,i);       Gam2 = NeighLn(2,i)
    tau = TVertex(2, Gam2)-TVertex(1,Gam1)
    gln(ikey) = weight_gline(StatusLn(i),tau,TypeLn(i))
    weight = weight *gln(ikey)
  enddo

  do ikey = 1, NWLn
    i = WLnKey2Value(ikey)
    Gam1 = NeighLn(1,i);       Gam2 = NeighLn(2,i)
    wln(ikey) = weight_wline(StatusLn(i),IsDeltaLn(i), WRVertex(:, Gam1)-WRVertex(:, Gam2), &
      & TVertex(3, Gam2)-TVertex(3,Gam1), TypeLn(i))
    weight = weight *wln(ikey)
  enddo

  do ikey = 1, NVertex
    i = VertexKey2Value(ikey)
    tau1 = TVertex(3, i)-TVertex(2, i)
    tau2 = TVertex(1, i)-TVertex(3, i)
    gam(ikey) = weight_vertex(StatusVertex(i), IsDeltaVertex(i), GRVertex(:, i)-WRVertex(:, i), &
      & tau1, tau2, TypeVertex(i))
    weight = weight *gam(ikey)
  enddo

  !weight = weight*(1.d0/Beta)**Order *SignFermiLoop
  weight = weight *(-1.d0)**Order*SignFermiLoop

  if(abs(real(Phase*WeightCurrent-weight))>1.d-6.or.abs(dimag(Phase*WeightCurrent-weight))>1.d-6) then
    call LogFile%WriteStamp('e')
    call LogFile%WriteLine("Oops, check_weight found a bug!")
    call LogFile%WriteLine("IsWormPresent"+str(IsWormPresent)+",imc"+str(imc))
    call LogFile%WriteLine("real weight"+str(weight))
    call LogFile%WriteLine("current weight"+str(Phase*WeightCurrent))
    call LogFile%WriteLine(str(Order)+str(Beta))

    do ikey = 1, NGLn
      i = GLnKey2Value(ikey)
      call LogFile%WriteLine("G"+str(i)+str(gln(ikey))+str(WeightLn(i)))
    enddo

    do ikey = 1, NWLn
      i = WLnKey2Value(ikey)
      call LogFile%WriteLine("W"+str(i)+str(wln(ikey))+str(WeightLn(i)))
    enddo

    do ikey = 1, NVertex
      i = VertexKey2Value(ikey)
      call LogFile%WriteLine("Gam"+str(i)+str(gam(ikey))+str(WeightVertex(i)))
    enddo
    call print_config
    stop
  endif
  return
END SUBROUTINE check_weight

!!=======================================================================
!!=======================================================================
!!=======================================================================


