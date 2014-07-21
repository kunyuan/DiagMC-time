!================================================================
!================= print/draw information    ====================
!================================================================

!==================================================================
!===================== PRINT STATUS ================================
!==================================================================
!If you want to log staff when runing markov, log here!
SUBROUTINE print_status
    implicit none
    integer :: iorder,i
    character*30 :: updatename(Nupdate)

    call LogFile%WriteStamp()
    call LogFile%WriteLine("MC steps:"+str(imc))
    call time_elapse
    call LogFile%WriteLine("Printing interval:"+str(t_elap,'(f12.3)')+'s')
    call LogFile%WriteLine("Efficiency: "+str(imc/t_elap,'(f12.0)')+"steps per second.")
    call LogFile%WriteLine('Statistics Number ='+str(StatNum))

    do i=0,MCOrder
      if(Norm(i)>1e-6) then
        write(logstr,"(i2, A,f15.6,'+/-',f15.6)") i,QuanName(i),Quan(i)/Norm(i),Error(i) 
        call LogFile%WriteLine(logstr)
      endif
    enddo

    call LogFile%WriteLine("------------------------------------------------")

    updatename(1)= " 1: create worm along wline"
    updatename(2)= " 2: delete worm along wline"
    updatename(3)= " 3: create worm along gline"
    updatename(4)= " 4: delete worm along gline"
    updatename(5)= " 5: move worm along wline"
    updatename(6)= " 6: move worm along gline"
    updatename(7)= " 7: add interaction"
    updatename(8)= " 8: remove interaction"
    updatename(9)= " 9: add interaction cross"
    updatename(10)= "10: remove interaction cross"
    updatename(11)= "11: reconnect"
    updatename(12)= "12: shift gline in space"
    updatename(13)= "13: shift wline in space"
    updatename(14)= "14: change Gamma type"
    updatename(15)= "15: move measuring index"
    updatename(16)= "16: change Gamma time"
    updatename(17)= "17: change wline isdelta"
    updatename(18)= "18: change Gamma isdelta"

    do iorder = 0, MCOrder
      call LogFile%WriteLine("Order"+str(iorder))
      do i = 1, Nupdate
        if(ProbProp(iorder, i)/=0.d0) then
          write(logstr, '(A,2f13.1,f13.5,f13.1,f13.5)') updatename(i), ProbCall(iorder, i), ProbProp(iorder, i), &
            & ProbProp(iorder, i)/ProbCall(iorder, i), ProbAcc(iorder, i), &
            & ProbAcc(iorder, i)/ProbProp(iorder, i)
          call LogFile%WriteLine(logstr)
        endif
      enddo
    enddo

    call LogFile%WriteLine("------------------------------------------------")
    do iorder = 0, MCOrder
      call LogFile%WriteLine("Order"+str(iorder))
      BalenceCheck(iorder,1,3)=(BalenceCheck(iorder,1,1)-BalenceCheck(iorder,1,2)) &
        & /sqrt(BalenceCheck(iorder,1,1))
      write(logstr, '(A,3f17.5)') "Gamma Type 1,2 <==> 3,4:",BalenceCheck(iorder,1,:)
      call LogFile%WriteLine(logstr)
    enddo
    call LogFile%WriteLine("------------------------------------------------")
    call LogFile%WriteLine("Reducibility ratio "+str(TestData(2)/TestData(1))+" : "+str(TestData(3)/TestData(1)))

END SUBROUTINE print_status

!====================================================================
!===================== PRINT CONFIGURATION ==========================
!====================================================================


SUBROUTINE print_config
  implicit none
  integer :: i, iln, iv
  
  open(108, access='append', file=trim(title_mc)//"_mc.conf")
  
  write(108, *) "============================================================"
  write(108, *) imc, IsWormPresent, iupdate

  if(IsWormPresent .eqv. .true.) then
    write(108, *) "Ira", Ira, "Masha", Masha, "SpinMasha", SpinMasha
    write(108, *) "kMasha", kMasha
  endif

  write(108, *) "Order", Order
  write(108, *) "SignFermiLoop", SignFermiLoop

  write(108, *) "Measuring Gamma", MeasureGam
  write(108, *) "Phase", Phase
  write(108, *) "Weight", WeightCurrent

  do i = 1, NGLn
    iln = GLnKey2Value(i)
    if(StatusLn(iln) <0) cycle
    write(108, 10) iln, KindLn(iln), IsDeltaLn(iln), TypeLn(iln), kLn(iln), StatusLn(iln), NeighLn(1:2,iln)
  enddo

  do i = 1, NWLn
    iln = WLnKey2Value(i)
    if(StatusLn(iln) <0) cycle
    write(108, 10) iln, KindLn(iln), IsDeltaLn(iln), TypeLn(iln), kLn(iln), StatusLn(iln), NeighLn(1:2,iln)
  enddo

  do i = 1, NVertex
    iv = VertexKey2Value(i)
    if(StatusVertex(iv) <0) cycle
    write(108, 12) iv,IsDeltaVertex(iv), TypeVertex(iv), SpInVertex(1, iv),SpInVertex(2, iv), &
      & GRVertex(iv),WRVertex(iv), TVertex(1, iv), TVertex(2, iv),  &
      & TVertex(3, iv), DirecVertex(iv), StatusVertex(iv), NeighVertex(:,iv)
  enddo
  write(108, *) "============================================================"

  10 format(' Line:',i2,2x,'kind:',i2,2x,'isdelta:',i2,2x,'type:',i2,2x,'k:',i8,2x,'stat:',i2, 2x,&
    & 'neigh:',i6,i6)
  12 format('Gamma:',i2,2x,'isdelta:',i2,2x,'type:',i2,2x,'typein:',i2,2x,'typeout:',i2,2x,&
    & 'gr:(',i4,'), wr:(',i4,')', 't:(', f7.4, f7.4, f7.4, ')',2x, &
    & 'direction:', i2,2x, 'stat:',i2, 2x,'neigh:', i6,i6,i6)

  close(108)
  call DRAW
END SUBROUTINE print_config

!=================== VISUALIZATION  ==================================
SUBROUTINE DRAW
    IMPLICIT NONE
    DOUBLE PRECISION :: x1,y1, x2,y2, y3, scx, scy, sgn   
    DOUBLE PRECISION :: scydash, ca1,ca2,ra,a1,a2, radian
    DOUBLE PRECISION :: phi1, phi2, pi2, pi4,ini,seg,theta2
    INTEGER :: i
    INTEGER :: FEXIST, RES
    integer :: iWLn,iGLn,Vertex1,Vertex2,Vertex3
    character(10) :: imcstr
    character(100) :: tempstr
    
    pi2=dasin(1.d0)
    pi4=pi2/2.d0
    theta2=dasin(0.2d0)
    !theta2=pi4
    radian=90.d0/pi2
    scx=500/beta
    scy=400./Vol
    x1=scx*beta
    y1=scy*Vol
    scydash=scy/40.

    INQUIRE(DIRECTORY="graph",EXIST=FEXIST)
    IF(.not.FEXIST) THEN
      call system("mkdir graph")
    ENDIF

    write(imcstr,'(i10)') int(imc)
    open(11, file='graph/graph_'//trim(adjustl(title2))//'_'//trim(adjustl(imcstr))//'.eps')
    write(11,*) '%!'
    write(11,*) '%BoundingBox: 0 0 ', x1, y1
    write(11,*) '%%EndComments'
    write(11,*) '%%BeginProlog'
    write(11,*) '/L { lineto stroke} def'
    write(11,*) '/M { moveto } def'
    write(11,*) '/N {newpath } def'
    write(11,*) '/Ndashed {[5 5] 0 setdash newpath } def'
    write(11,*) '/Nsolid {[] 0 setdash newpath } def'
    write(11,*) '/Y { 0 360 arc closepath gsave fill'
    write(11,*) '     grestore stroke } def'
    write(11,*) '/YL { 0 360 arc stroke} def' 
    write(11,*) '/Y45 { arc stroke} def'
    write(11,*) '/C {/Times-Roman findfont 15 scalefont setfont show} def'
    write(11,*) '% Put an arrowhead at point x2 y2,'
    write(11,*) '% pointing away from x1 y1'
    write(11,*) '% Replace x2 y2 with coordinates of arrowbase:'
    write(11,*) '% the point to connect lines to'
    write(11,*) '% ArrowHeadSize gives the size of the arrow'
    write(11,*) '/ArrowHeadSize 10 def'
    write(11,*) '/ahead {'
    write(11,*) '    1 index 4 index sub'
    write(11,*) '    1 index 4 index sub'
    write(11,*) '    exch atan'
    write(11,*) '    ArrowHeadSize -.8 mul'
    write(11,*) '    dup'
    write(11,*) '    2 index cos mul 4 index add'
    write(11,*) '    exch'
    write(11,*) '    2 index sin mul 3 index add'
    write(11,*) '    5 2 roll'
    write(11,*) '    gsave'
    write(11,*) '        3 1 roll'
    write(11,*) '        translate'
    write(11,*) '        rotate'
    write(11,*) '        newpath'
    write(11,*) '        0 0 moveto'
    write(11,*) '        ArrowHeadSize dup neg exch .25 mul'
    write(11,*) '        2 copy lineto'
    write(11,*) '        ArrowHeadSize -.8 mul 0'
    write(11,*) '        2 copy'
    write(11,*) '        6 4 roll'
    write(11,*) '        neg curveto'
    write(11,*) '        closepath fill'
    write(11,*) '    grestore'
    write(11,*) '} bind def'
    write(11,*) ''
    write(11,*) '%%EndProlog'
    write(11,*) '%%BeginSetup'
    write(11,*) '1 setlinewidth'
    write(11,*) '5 140 translate'
    write(11,*) '1 1 scale'
    write(11,*) '%%EndSetup'

    ini=640.0
    seg=15.0
    write(11,*) '0 0 0 setrgbcolor'
    write(tempstr, *) "Beta: ",beta,"    L: ",L(1:D)
    write(11,803) 0.,ini,tempstr
    write(tempstr,*)  "Jcp: ",Jcp," Seed: ",Seed
    ini=ini-seg
    write(11,803) 0.,ini,tempstr
    write(tempstr, *) "imc: ",imc," Is Worm Here:", IsWormPresent
    ini=ini-seg
    write(11,803) 0.,ini,tempstr
    if(IsWormPresent .eqv. .true.) then
      write(tempstr, *) "Ira: ", Ira, "Masha: ", Masha
      ini=ini-seg
      write(11,803) 0.,ini,tempstr
      write(tempstr, *) "Spin of Masha: ", SpinMasha, "k of Masha: ", kMasha
      ini=ini-seg
      write(11,803) 0.,ini,tempstr
    endif

    write(tempstr, *) "Order: ", Order
    ini=ini-seg
    write(11,803) 0.,ini,tempstr
    write(tempstr, *) "Fermi Loop Sign: ", SignFermiLoop
    ini=ini-seg
    write(11,803) 0.,ini,tempstr

    write(tempstr, *) "Measuring Gamma: ", MeasureGam
    ini=ini-seg
    write(11,803) 0.,ini,tempstr
    write(tempstr, *) "Phase: ", Phase
    ini=ini-seg
    write(11,803) 0.,ini,tempstr
    write(tempstr, *) "Weight: ", WeightCurrent
    ini=ini-seg
    write(11,803) 0.,ini,tempstr

    ini=455.0
    write(11,*) '0 1 0 setrgbcolor'
    write(11,777) 500.0, ini, scy/20.
    write(11,803) 520.0, ini-5.0, "MeasureGam"
    ini=ini-seg
    write(11,*) '1 0 0 setrgbcolor'
    write(11,777) 500.0, ini, scy/30.
    write(11,803) 520.0, ini-5.0, "Ira"
    ini=ini-seg
    write(11,*) '0 0 1 setrgbcolor'
    write(11,777) 500.0, ini, scy/30.
    write(11,803) 520.0, ini-5.0, "Masha"

    write(11,*) '0 0 0 setrgbcolor'
    do i=1,NWLn;
      iWLn=WLnKey2Value(i)
      Vertex1=NeighLn(1,iWLn)
      x1=scx*TVertex(3, Vertex1)
      y1=scy*GRVertex(Vertex1)
      Vertex2=NeighLn(2,iWLn)
      x2=scx*TVertex(3, Vertex2)
      y2=scy*GRVertex(Vertex2)

      ra=dsqrt((x2-x1)**2+(y2-y1)**2)/2.d0
      if(dabs(ra)<1e-6) then
        write(11,792) x1, y1+scy/5. , scy/5. 
        cycle
      endif
      phi1=(y2-y1)/(2.d0*ra)
      phi2=(x2-x1)/(2.d0*ra)
      ca1=(x1+x2)/2.d0
      ca2=(y1+y2)/2.d0;
      ra=ra/sin(theta2)*cos(theta2)
      ca1=ca1+phi1*ra
      ca2=ca2-phi2*ra; 
      ra=ra/cos(theta2)
      a1=pi2-theta2+dasin(phi1)
      a2=a1+2.0*theta2
      IF(phi2.lt.0) then
        a2=4*pi2-pi2+theta2-dasin(phi1)
        a1=a2-2.0*theta2
      endif
      a1=a1*radian
      a2=a2*radian 

      write(11,791)  ca1, ca2, ra, a1, a2
    enddo
  
    do i=1,NGLn;
      iGLn=GLnKey2Value(i)
      Vertex1=NeighLn(1,iGLn)
      x1=scx*TVertex(3, Vertex1)
      y1=scy*GRVertex(Vertex1)
      Vertex2=NeighLn(2,iGLn)
      x2=scx*TVertex(3, Vertex2)
      y2=scy*GRVertex(Vertex2)

      if(TypeLn(iGLn)==1) then
        write(11,*) '1 0 0 setrgbcolor'
      else
        write(11,*) '0 0 1 setrgbcolor'
      endif


      if(Vertex1/=Vertex2) then

        ra=dsqrt((x2-x1)**2+(y2-y1)**2)/2.d0
        if(dabs(ra)<1e-6) then
          iWLn=NeighVertex(3,Vertex1)
          Vertex3=NeighLn(3-DirecVertex(Vertex1),iWLn)
          y3=scy*GRVertex(Vertex3)
          sgn=1.d0
          IF(y3>y1) sgn=-1.d0   
          write(11,780) x1, y1+sgn*scy/5. , scy/5. 
          cycle
        endif
        phi1=(y2-y1)/(2.d0*ra)
        phi2=(x2-x1)/(2.d0*ra)
        ca1=(x1+x2)/2.d0
        ca2=(y1+y2)/2.d0;
        ca1=ca1+phi1*ra
        ca2=ca2-phi2*ra; 
        ra=ra*dsqrt(2.d0)
        a1=pi4+dasin(phi1)
        a2=a1+pi2
        IF(phi2.lt.0) then
          a2=4*pi2-pi4-dasin(phi1)
          a1=a2-pi2
        endif
        a1=a1*radian
        a2=a2*radian 
        write(11,781)  ca1, ca2, ra, a1, a2  ! propagator lines - arcs 

        ca1=x1+(phi2-phi1)*scy/2.;     
        ca2=y1+(phi2+phi1)*scy/2. 

        write(11,790) ca1, ca2, x1, y1       ! arrows
      else
        iWLn=NeighVertex(3,Vertex1)
        Vertex3=NeighLn(3-DirecVertex(Vertex1),iWLn)
        y3=scy*GRVertex(Vertex3)
        sgn=1.d0
        IF(y3>y1) sgn=-1.d0   
        write(11,780) x1, y1+sgn*scy/5. , scy/5. 
      endif
    enddo

    write(11,*) '0 0 0 setrgbcolor'
    write(11,"(f6.1,x,f6.1,x,' M (Gamma info) C')") 520.0,350.0
    ini=350.0-seg
    ! write Gamma info
    do i=1,NVertex
      Vertex1=VertexKey2Value(i)
      x1=scx*TVertex(3, Vertex1)
      y1=scy*GRVertex(Vertex1)
      if(Vertex1==MeasureGam) then
        write(11,*) '0 1 0 setrgbcolor'
        write(11,777) x1, y1, scy/20.
      endif
      if(Vertex1==Ira .and. IsWormPresent) then
        write(11,*) '1 0 0 setrgbcolor'
      elseif(Vertex1==Masha .and. IsWormPresent) then
        write(11,*) '0 0 1 setrgbcolor'
      else
        write(11,*) '0 0 0 setrgbcolor'
      endif
      write(11,777) x1, y1, scy/30.
      write(11,801) x1-5., y1+7., Vertex1
      if(Vertex1==MeasureGam) then
        write(11,*) '0 1 0 setrgbcolor'
      endif
      write(11,802) 500.0,ini,Vertex1,TVertex(3, Vertex1), &
         & GRVertex(Vertex1)
      ini=ini-seg
    enddo

    write(11,*) '0 0 0 setrgbcolor'
    write(11,"(f6.1,x,f6.1,x,' M (G info) C')") 530.0,ini
    ini=ini-seg
    do i=1,NGLn
      iGLn=GLnKey2Value(i)
      Vertex1=NeighLn(1,iGLn)
      Vertex2=NeighLn(2,iGLn)
      if(TypeLn(iGLn)==1) then
        write(11,*) '1 0 0 setrgbcolor'
      else
        write(11,*) '0 0 1 setrgbcolor'
      endif
      write(11,804) 500.0, ini, iGLn, Vertex2,Vertex1
      ini=ini-seg
    enddo
    write(11,*) '0 0 0 setrgbcolor'
    write(11,"(f6.1,x,f6.1,x,' M (W info) C')") 530.0,ini
    ini=ini-seg
    do i=1,NWLn;
      iWLn=WLnKey2Value(i)
      Vertex1=NeighLn(1,iWLn)
      Vertex2=NeighLn(2,iWLn)
      write(11,805) 500.0, ini, iWLn, Vertex1,Vertex2
      ini=ini-seg
    enddo


    write(11,*) ''
    write(11,*) 'stroke showpage'
    write(11,*) '%%Trailer'

    close (11)
  return
!******************************************************************
 777  format ('Nsolid ',f6.1,x,f6.1,x,f9.3,' Y')
 778  format ('N ',f6.1,x,f6.1,x,' M')
 779  format (     f6.1,x,f6.1,x,' L')
 780  format ('Nsolid ',f6.1,x,f6.1,x,f9.3,' YL')
 781  format ('Nsolid ',f6.1,x,f6.1,x,f6.1,x,f6.1,x,f6.1,' Y45')
 790  format ('Nsolid ',f6.1,x,f6.1,x,f6.1,x,f6.1,x,' ahead')
 !791  format ('Ndashed ',f6.1,x,f6.1,' M',f6.1,x,f6.1,' L')
 791  format ('Ndashed ',f6.1,x,f6.1,x,f6.1,x,f6.1,x,f6.1,' Y45')
 792  format ('Ndashed ',f6.1,x,f6.1,x,f9.3,' YL')
 801  format (f6.1,x,f6.1,x,' M (',i2,') C')
 802  format (f6.1,x,f6.1,x,' M (',i2,':',f6.3,'; (',i3,')) C')
 803  format (f6.1,x,f6.1,x,' M (',A,') C')
 804  format (f6.1,x,f6.1,x,' M (G ',i2,' : ',i2,' --->',i2') C')
 805  format (f6.1,x,f6.1,x,' M (W ',i2,' : ',i2,' <=>',i2') C')
 
END SUBROUTINE DRAW

!__________________________________________ 
!   	'0 0 0 setrgbcolor'    black
!     '0 0 1 setrgbcolor'    blue
!     '1 0 0 setrgbcolor'    red
!     '0 1 0 setrgbcolor'    green

!====================================================================
!====================================================================

!!================================================================
!!================= READ/WRITE CONFIGURATIONS ====================
!!================================================================
SUBROUTINE write_monte_carlo_data
  implicit none
  integer :: iorder, itopo, ir, irr, ityp, it1, it2
  double precision :: rgam2, rerr
  integer :: ibin, ibasis
  complex*16 :: gam1

  ratioerr = 1.d0

  !=========== write into files =========================================
  open(104, status="replace", &
    & file=trim(title_mc)//"_monte_carlo_data.bin.dat",form="binary")

  write(104) finalBeta, Beta, MCOrder, L(1:D)
  write(104) imc, GamNorm, GamNormWeight
  write(104) Z_normal, ratioerr
  do it1 = 0, MxT-1
    do ir = 0, Vol-1
      do iorder = 0, MCOrder
        write(104)  GamMC(iorder, ir, it1)
        write(104)  ReGamSqMC(iorder, ir, it1)
        write(104)  ImGamSqMC(iorder, ir, it1)
      enddo
    enddo
  enddo

  do ibasis = 1, NBasisGam
    do ibin = 1, NbinGam
      do ir = 0, Vol-1
        do ityp = 1, NtypeGam/2
          do iorder = 0, MCOrder
            write(104) GamBasis(iorder, ityp, ir, ibin, ibasis)
            write(104) ReGamSqBasis(iorder, ityp, ir, ibin, ibasis)
            write(104) ImGamSqBasis(iorder, ityp, ir, ibin, ibasis)
          enddo
        enddo
      enddo
    enddo
  enddo
  !=========  write on the screen ========================================

  
  close(104)
END SUBROUTINE write_monte_carlo_data


SUBROUTINE read_monte_carlo_data
  implicit none
  integer :: iorder, ir, ityp, it1, it2, itopo,ios
  integer :: ibin, ibasis
  integer :: iL(1:D)
  double precision :: fBeta
  logical :: alive
  

  inquire(file=trim(title)//"_monte_carlo_data.bin.dat",exist=alive)
  if(.not. alive) then
    call LogFile%QuickLog("There is no monte carlo binary data yet!",'e')
    stop -1
  endif

  open(105, status="old", file=trim(title)//"_monte_carlo_data.bin.dat",form="binary")
  read(105,iostat=ios) fBeta, iorder, iL(1:D)
  read(105,iostat=ios) imc, GamNorm, GamNormWeight
  read(105,iostat=ios) Z_normal, ratioerr

  do it1 = 0, MxT-1
    do ir = 0, Vol-1
      do iorder = 0, MCOrder
        read(105,iostat=ios)  GamMC(iorder, ir, it1)
        read(105,iostat=ios)  ReGamSqMC(iorder,  ir, it1)
        read(105,iostat=ios)  ImGamSqMC(iorder,  ir, it1)
      enddo
    enddo
  enddo

  do ibasis = 1, NBasisGam
    do ibin = 1, NbinGam
      do ir = 0, Vol-1
        do ityp = 1, NtypeGam/2
          do iorder = 0, MCOrder
            read(105) GamBasis(iorder, ityp, ir, ibin, ibasis)
            read(105) ReGamSqBasis(iorder, ityp, ir, ibin, ibasis)
            read(105) ImGamSqBasis(iorder, ityp, ir, ibin, ibasis)
          enddo
        enddo
      enddo
    enddo
  enddo

  close(105)
  if(ios/=0) then
    call LogFile%QuickLog("The monte carlo binary data is broken?",'e')
    stop -1
  endif
END SUBROUTINE read_monte_carlo_data


SUBROUTINE read_Gamma
  implicit none
  logical :: alive, changeBeta

  inquire(file=trim(title)//"_monte_carlo_data.bin.dat",exist=alive)
  if(.not. alive) then
    call LogFile%QuickLog("There is no monte carlo binary data yet!")

    call initialize_Gam
    call LogFile%QuickLog("Initialize the Gamma Done!")
  else

    call LogFile%QuickLog("Reading MC data...")
    call read_monte_carlo_data

    call LogFile%QuickLog("Reading Done!...")

    !!-------- update the Gamma matrix with MC data -------
    call Gam_mc2matrix_mc(changeBeta)
    call LogFile%QuickLog("Tabulate the MC data Done!...")

  endif

END SUBROUTINE read_Gamma



SUBROUTINE write_monte_carlo_conf
  implicit none
  integer :: i

  open(103, status="replace", &
   & file=trim(title_mc)//"_monte_carlo_conf.bin.dat",form="binary")

  write(103)  mc_version, Z_normal, Z_worm, StatNum
  write(103)  ProbCall(:,:), ProbProp(:,:), ProbAcc(:,:)
  write(103)  GamOrder(:), GamWormOrder(:)
  write(103)  CoefOfWorm, CoefOfWeight(:)
  write(103)  TestData(:)
  write(103)  TailLn, TailVertex
  write(103)  NextLn(:)
  write(103)  NextVertex(:)
  write(103)  Order, MeasureGam, SignFermiLoop, IsWormPresent 
  write(103)  Ira, Masha, SpinMasha, kMasha
  write(103)  StatusVertex(:)
  write(103)  TypeVertex(:)
  write(103)  GRVertex(:)
  write(103)  WRVertex(:)
  write(103)  TVertex(1, :), TVertex(2, :), TVertex(3, :)
  write(103)  DirecVertex(:)             
  write(103)  IsDeltaVertex(:)
  write(103)  IsDeltaLn(:)
  write(103)  kLn(:)
  write(103)  KindLn(:)

  do i = 1, 2
    write(103)  NeighLn(i,:)
  enddo

  do i = 1, 3
    write(103)  NeighVertex(i,:) 
  enddo

  close(103)
END SUBROUTINE write_monte_carlo_conf



SUBROUTINE read_monte_carlo_conf
  implicit none
  integer :: i,ios
  integer :: ikey, ikeyG, ikeyW
  integer :: iGin, iGout, iW
  integer :: iGam, jGam
  complex*16 :: ComCurrent
  double precision :: tau
  logical :: alive

  inquire(file=trim(title)//"_monte_carlo_conf.bin.dat",exist=alive)
  if(.not. alive) then
    call LogFile%QuickLog("There is no monte carlo binary configuration yet!",'e')
    stop -1
  endif
  open(106, status="old", file=trim(title)//"_monte_carlo_conf.bin.dat",form="binary")

  read(106,iostat=ios)  mc_version, Z_normal, Z_worm, StatNum
  read(106,iostat=ios)  ProbCall(:,:), ProbProp(:,:), ProbAcc(:,:)
  read(106,iostat=ios)  GamOrder(:), GamWormOrder(:)
  read(106,iostat=ios)  CoefOfWorm, CoefOfWeight(:)
  read(106,iostat=ios)  TestData(:)
  read(106,iostat=ios)  TailLn, TailVertex
  read(106,iostat=ios)  NextLn(:)
  read(106,iostat=ios)  NextVertex(:)
  read(106,iostat=ios)  Order, MeasureGam, SignFermiLoop, IsWormPresent 
  read(106,iostat=ios)  Ira, Masha, SpinMasha, kMasha
  read(106,iostat=ios)  StatusVertex(:)
  read(106,iostat=ios)  TypeVertex(:)
  read(106,iostat=ios)  GRVertex(:)
  read(106,iostat=ios)  WRVertex(:)
  read(106,iostat=ios)  TVertex(1, :), TVertex(2, :), TVertex(3, :)
  read(106,iostat=ios)  DirecVertex(:)             
  read(106,iostat=ios)  IsDeltaVertex(:)
  read(106,iostat=ios)  IsDeltaLn(:)
  read(106,iostat=ios)  kLn(:)
  read(106,iostat=ios)  KindLn(:)

  do i = 1, 2
    read(106,iostat=ios)  NeighLn(i,:)
  enddo

  do i = 1, 3
    read(106,iostat=ios)  NeighVertex(i,:) 
  enddo
  close(106)

  if(ios/=0) then
    call LogFile%QuickLog("The monte carlo binary configuration is broken?",'e')
    stop -1
  endif

  NVertex = 2*(Order+1)
  NGLn = NVertex
  NWLn = Order+1
  ComCurrent = (1.d0, 0.d0)

  ikey = 1
  do i = 1, MxNVertex
    if(StatusVertex(i)==-1) cycle
    VertexValue2Key(i) = ikey
    VertexKey2Value(ikey) = i
    ikey = ikey + 1
  enddo

  if(ikey/=NVertex+1) then
    call LogFile%QuickLog("read_monte_carlo_conf: Number of Vertex Error!",'e')
    stop
  endif

  StatusLn(:)=-1
  do ikey = 1, NVertex
    i = VertexKey2Value(ikey)
    if(StatusVertex(i)==-1) then
      call LogFile%QuickLog("read_monte_carlo_conf: Status of Vertex Error!",'e')
      stop
    endif

    iGin = NeighVertex(1, i)
    iGout = NeighVertex(2, i)
    iW = NeighVertex(3, i)

    if(StatusLn(iGin)==-1) then
      StatusLn(iGin) = gline_stat(StatusVertex(NeighLn(1, iGin)), StatusVertex(NeighLn(2, iGin)))
      KindLn(iGin) = 1
    endif
    if(StatusLn(iGout)==-1) then
      StatusLn(iGout) = gline_stat(StatusVertex(NeighLn(1, iGout)), StatusVertex(NeighLn(2, iGout)))
      KindLn(iGout) = 1
    endif
    if(StatusLn(iW)==-1) then
      StatusLn(iW) = wline_stat(StatusVertex(NeighLn(1, iW)), StatusVertex(NeighLn(2, iW)))
      KindLn(iW) = 2
    endif

    if(TypeVertex(i)==1 .or. TypeVertex(i)==4) then
      SpInVertex(:, i) = 1
    else if(TypeVertex(i)==2 .or. TypeVertex(i)==3) then
      SpInVertex(:, i) = 2
    else if(TypeVertex(i)==5 .or. TypeVertex(i)==6) then
      SpInVertex(2, i) = Mod(TypeVertex(i), 2)+1
      SpInVertex(1, i) = 3 - SpInVertex(2, i)
    endif

    WeightVertex(i) = weight_vertex(StatusVertex(i), IsDeltaVertex(i), &
      & diff_r(D, GRVertex(i),WRVertex(i)), &
      & TVertex(3,i)-TVertex(2,i), TVertex(1,i)-TVertex(3,i), TypeVertex(i))

    ComCurrent = ComCurrent* WeightVertex(i)
  enddo

  ikeyG = 1
  ikeyW = 1
  do i = 1, MxNLn
    if(StatusLn(i)==-1) cycle
    if(KindLn(i)==1) then
      LnValue2Key(i) = ikeyG
      GLnKey2Value(ikeyG) = i
      ikeyG = ikeyG + 1
    else if(KindLn(i)==2) then
      LnValue2Key(i) = ikeyW
      WLnKey2Value(ikeyW) = i
      ikeyW = ikeyW + 1
    endif
  enddo

  if(ikeyG/=NGLn+1) then
    call LogFile%QuickLog("read_monte_carlo_conf: Number of Glines Error!",'e')
    stop
  endif

  if(ikeyW/=NWLn+1) then
    call LogFile%QuickLog("read_monte_carlo_conf: Number of Glines Error!",'e')
    stop
  endif

  do ikeyG = 1, NGLn
    i = GLnKey2Value(ikeyG)
    if(Is_k_valid_for_G(kLn(i))) then
      call add_Hash4G(kLn(i),i)
    else
      call LogFile%QuickLog("read_monte_carlo_conf: k of G Error!",'e')
      stop
    endif
    TypeLn(i) = mod(TypeVertex(NeighLn(2,i)), 2)
    if(TypeLn(i)==0)  TypeLn(i) = 2
    tau = TVertex(2, NeighLn(2, i))-TVertex(1,NeighLn(1,i))
    WeightLn(i) = weight_gline(StatusLn(i), tau, TypeLn(i))
    ComCurrent = ComCurrent* WeightLn(i)
  enddo

  do ikeyW = 1, NWLn
    i = WLnKey2Value(ikeyW)
    if(Is_k_valid_for_W(kLn(i))) then
      call add_Hash4W(kLn(i),i)
    else
      call LogFile%QuickLog("read_monte_carlo_conf: k of W Error!",'e')
      stop
    endif
    iGam = NeighLn(1, i)
    jGam = NeighLn(2, i)
    if(StatusLn(i)<=1) then
      TypeLn(i) = TypeGam2W(TypeVertex(iGam), TypeVertex(jGam))
    else
      TypeLn(i) = 1
    endif

    tau = TVertex(3, iGam)-TVertex(3, jGam)
    WeightLn(i) = weight_wline(StatusLn(i), IsDeltaLn(i), diff_r(D, WRVertex(jGam),WRVertex(iGam)), &
      & tau, TypeLn(i))
    ComCurrent = ComCurrent* WeightLn(i)
  enddo

  WeightCurrent = abs(ComCurrent) *(1.d0/Beta)**Order *SignFermiLoop
  Phase = ComCurrent/abs(ComCurrent) *SignFermiLoop

  return
END SUBROUTINE read_monte_carlo_conf
!!================================================================
!!================================================================
!!================================================================





  
!!================================================================
!!========== PRINT OUT THE DATA FILES ============================
!!================================================================

SUBROUTINE output_GamMC
  implicit none
  integer :: i, j, iorder, it1, it2
  double precision :: rerr, ierr, rpercenterr, ipercenterr
  complex*16 :: gam, gamn, normal 
  double precision :: rgam2, igam2

  open(35, access="append", file=trim(title_mc)//"_Gam_MC.dat")

  normal = GamNormWeight*Z_normal/GamNorm

  ratioerr = 1.d0

  write(35, *) "============================================"
  write(35, *) "Beta", Beta, "L", L(1), "Order", MCOrder, "Seed",Seed
  write(35, *) imc, Z_normal, GamNormWeight, GamNorm, normal, ratioerr

  do iorder = 1, MCOrder
    write(35, *) "Order", iorder
    write(35, *) "dx = 0, dy = 0"
    do it1 = 0, MxT-1
      gam = GamMC(iorder, 0, it1)/Z_normal
      gamn = gam*normal

      rgam2 = ReGamSqMC(iorder, 0, it1)/Z_normal
      rerr = sqrt(abs(rgam2)-(real(gam))**2.d0)/sqrt(Z_normal-1)
      rerr = rerr* ratioerr

      if(abs(real(gam))<1.d-30) then
        rpercenterr = 0.d0
      else
        rpercenterr = rerr/abs(real(gam))
      endif

      igam2 = ImGamSqMC(iorder, 0, it1)/Z_normal
      ierr = sqrt(abs(igam2)-(dimag(gam))**2.d0)/sqrt(Z_normal-1)
      ierr = ierr* ratioerr

      if(abs(dimag(gam))<1.d-30) then
        ipercenterr = 0.d0
      else
        ipercenterr = ierr/abs(dimag(gam))
      endif

      write(35, '(i3,2x, E20.10E3,"+/-",f13.6,"%    +i",E20.10E3,"+/-",f13.6,"%")') it1, &
        & real(gamn), rpercenterr, dimag(gamn), ipercenterr
    enddo
    write(35, *)
  enddo

  close(35)
END SUBROUTINE output_GamMC

SUBROUTINE output_test
  implicit none
  integer :: iorder
  !open(104, status='replace', file=trim(title_mc)//"_Gam_Order_test.dat")
  open(104, access='append', file=trim(title_mc)//"_Gam_Order_test.dat")
  do iorder = 1, MCOrder
    write(104, *) iorder, Quan(iorder+2)/Quan(2), Error(iorder+2)/Quan(2)
  enddo
  write(104, *)

  write(104, *) "Type 1,2 <==> 3,4"
  do iorder = 0, MCOrder
    BalenceCheck(iorder,1,3)=(BalenceCheck(iorder,1,1)-BalenceCheck(iorder,1,2)) &
      & /sqrt(BalenceCheck(iorder,1,1))
    write(104, '(i3,3f17.5)') iorder, BalenceCheck(iorder,1,:)
  enddo
  write(104, *)
  close(104)
END SUBROUTINE
!!================================================================
!!================================================================
!!================================================================

