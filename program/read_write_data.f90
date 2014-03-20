!================================================================
!================= print/draw information    ====================
!================================================================

!==================================================================
!===================== PRINT STATUS ================================
!==================================================================
!If you want to log staff when runing markov, log here!
SUBROUTINE print_status
    implicit none
    integer :: iorder

    open(36, access="append", file=trim(title4)//".log")
    write(36,*) "MC steps: ",imc
    call time_elapse
    write(36,251) t_elap
  251 format(' Printing interval:',f16.7,2x,'s')
    write(36,*) "Efficiency: ",imc/t_elap," per second."

    write(36, *) " 1: create worm along wline"
    write(36, *) " 2: delete worm along wline"
    write(36, *) " 3: create worm along gline"
    write(36, *) " 4: delete worm along gline"
    write(36, *) " 5: move worm along wline"
    write(36, *) " 6: move worm along gline"
    write(36, *) " 7: add interaction"
    write(36, *) " 8: remove interaction"
    write(36, *) " 9: add interaction cross"
    write(36, *) "10: remove interaction cross"
    write(36, *) "11: reconnect"
    write(36, *) "12: shift gline in space"
    write(36, *) "13: shift wline in space"
    write(36, *) "14: change Gamma type"
    write(36, *) "15: move measuring index"
    write(36, *) "16: change Gamma time"
    write(36, *) "17: change wline isdelta"
    write(36, *) "18: change Gamma isdelta"
    do iorder = 0, MCOrder
      write(36, *) "Order", iorder
      do i = 1, Nupdate
        if(ProbProp(iorder, i)/=0.d0) then
          write(36, '(i3,3f17.5)') i, ProbProp(iorder, i), ProbAcc(iorder, i), ProbAcc(iorder, i)/ProbProp(iorder, i)
        endif
      enddo
      write(36, *)
    enddo

    close(36)
END SUBROUTINE print_status


SUBROUTINE write_log
  implicit none
  open(36, access="append", file=trim(title4)//".log")
  write(36,'(A)', advance='no') trim(logstr)
  close(36)
end SUBROUTINE

!====================================================================
!===================== PRINT CONFIGURATION ==========================
!====================================================================


SUBROUTINE print_config
  implicit none
  integer :: i, iln, iv
  
  open(8, access='append', file=trim(title1)//"_mc.conf")
  
  write(8, *) "============================================================"
  write(8, *) imc, IsWormPresent, iupdate

  if(IsWormPresent .eqv. .true.) then
    write(8, *) "Ira", Ira, "Masha", Masha, "SpinMasha", SpinMasha
    write(8, *) "kMasha", kMasha
  endif

  write(8, *) "Order", Order
  write(8, *) "SignFermiLoop", SignFermiLoop

  write(8, *) "Measuring Gamma", MeasureGam
  write(8, *) "Phase", Phase
  write(8, *) "Weight", WeightCurrent

  do i = 1, NGLn
    iln = GLnKey2Value(i)
    if(StatusLn(iln) <0) cycle
    write(8, 10) iln, KindLn(iln), IsDeltaLn(iln), TypeLn(iln), kLn(iln), StatusLn(iln), NeighLn(1:2,iln)
  enddo

  do i = 1, NWLn
    iln = WLnKey2Value(i)
    if(StatusLn(iln) <0) cycle
    write(8, 10) iln, KindLn(iln), IsDeltaLn(iln), TypeLn(iln), kLn(iln), StatusLn(iln), NeighLn(1:2,iln)
  enddo

  do i = 1, NVertex
    iv = VertexKey2Value(i)
    if(StatusVertex(iv) <0) cycle
    write(8, 12) iv,IsDeltaVertex(iv), TypeVertex(iv), SpInVertex(1, iv),SpInVertex(2, iv), &
      & GRVertex(1, iv),GRVertex(2, iv),WRVertex(1, iv),WRVertex(2, iv), TVertex(1, iv), TVertex(2, iv),  &
      & TVertex(3, iv), DirecVertex(iv), StatusVertex(iv), NeighVertex(:,iv)
  enddo
  write(8, *) "============================================================"

  10 format(' Line:',i2,2x,'kind:',i2,2x,'isdelta:',i2,2x,'type:',i2,2x,'k:',i8,2x,'stat:',i2, 2x,&
    & 'neigh:',i6,i6)
  12 format('Gamma:',i2,2x,'isdelta:',i2,2x,'type:',i2,2x,'typein:',i2,2x,'typeout:',i2,2x,&
    & 'gr:(',i4,i4,'), wr:(',i4,i4,')', 't:(', f7.4, f7.4, f7.4, ')',2x, &
    & 'direction:', i2,2x, 'stat:',i2, 2x,'neigh:', i6,i6,i6)

  close(8)
  call DRAW()
END SUBROUTINE print_config

!=================== VISUALIZATION  ==================================
SUBROUTINE DRAW()
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
    scy=400./Lx/Ly
    x1=scx*beta
    y1=scy*Vol
    scydash=scy/40.

    INQUIRE(DIRECTORY="graph",EXIST=FEXIST)
    IF(.not.FEXIST) THEN
      call system("mkdir graph")
    ENDIF

    write(imcstr,'(i10)') int(imc)
    open(11, file='graph/graph_'//trim(adjustl(imcstr))//'.eps')
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
    write(tempstr, *) "Beta: ",beta,"    Lx: ",Lx,"       Ly: ",Ly
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
      y1=scy*site_num(GRVertex(1, Vertex1),GRVertex(2, Vertex1))
      Vertex2=NeighLn(2,iWLn)
      x2=scx*TVertex(3, Vertex2)
      y2=scy*site_num(GRVertex(1, Vertex2),GRVertex(2, Vertex2))

      !if(TypeLn(iWLn)<=4) then
        !write(11,*) '0 0 0 setrgbcolor'
      !else
        !write(11,*) '1 0 1 setrgbcolor'
      !endif
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
      !write(11,791) x1,y1,x2,y2 
    enddo
  
    do i=1,NGLn;
      iGLn=GLnKey2Value(i)
      Vertex1=NeighLn(1,iGLn)
      x1=scx*TVertex(3, Vertex1)
      y1=scy*site_num(GRVertex(1, Vertex1),GRVertex(2, Vertex1))
      Vertex2=NeighLn(2,iGLn)
      x2=scx*TVertex(3, Vertex2)
      y2=scy*site_num(GRVertex(1, Vertex2),GRVertex(2, Vertex2))

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
          y3=scy*site_num(GRVertex(1, Vertex3),GRVertex(2, Vertex3))
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
        y3=scy*site_num(GRVertex(1, Vertex3),GRVertex(2, Vertex3))
        sgn=1.d0
        IF(y3>y1) sgn=-1.d0   
        write(11,780) x1, y1+sgn*scy/5. , scy/5. 
      endif
    enddo

    write(11,*) '0 0 0 setrgbcolor'
    write(11,"(f6.1,x,f6.1,x,' M (Gamma info) C')") 520.0,350.0
    ini=350.0-seg
    do i=1,NVertex
      Vertex1=VertexKey2Value(i)
      x1=scx*TVertex(3, Vertex1)
      y1=scy*site_num(GRVertex(1, Vertex1),GRVertex(2, Vertex1))
      if(Vertex1==MeasureGam) then
        write(11,*) '0 1 0 setrgbcolor'
        write(11,777) x1, y1, scy/20.
      endif
      if(Vertex1==Ira) then
        write(11,*) '1 0 0 setrgbcolor'
      elseif(Vertex1==Masha) then
        write(11,*) '0 0 1 setrgbcolor'
      else
        write(11,*) '0 0 0 setrgbcolor'
      endif
      write(11,777) x1, y1, scy/30.
      write(11,801) x1-5., y1+7., i
      write(11,802) 500.0,ini,i,TVertex(3, Vertex1), &
         & GRVertex(1, Vertex1),GRVertex(2, Vertex2)
      ini=ini-seg
    enddo
    write(11,*) '0 0 0 setrgbcolor'
    write(11,"(f6.1,x,f6.1,x,' M (G info) C')") 530.0,ini
    ini=ini-seg
    do i=1,NGLn
      iGLn=GLnKey2Value(i)
      Vertex1=NeighLn(1,iGLn)
      Vertex2=NeighLn(2,iGLn)
      if(TypeLn(i)==1) then
        write(11,*) '1 0 0 setrgbcolor'
      else
        write(11,*) '0 0 1 setrgbcolor'
      endif
      write(11,804) 500.0, ini, Vertex2,Vertex1
      ini=ini-seg
    enddo
    write(11,*) '0 0 0 setrgbcolor'
    write(11,"(f6.1,x,f6.1,x,' M (W info) C')") 530.0,ini
    ini=ini-seg
    do i=1,NWLn;
      iWLn=WLnKey2Value(i)
      Vertex1=NeighLn(1,iWLn)
      Vertex2=NeighLn(2,iWLn)
      write(11,805) 500.0, ini, Vertex1,Vertex2
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
 802  format (f6.1,x,f6.1,x,' M (',i2,':',f6.3,'; (',i3,',',i3,')) C')
 803  format (f6.1,x,f6.1,x,' M (',A,') C')
 804  format (f6.1,x,f6.1,x,' M (G: ',i2,' ------>',i2') C')
 805  format (f6.1,x,f6.1,x,' M (W: ',i2,' <====>',i2') C')
 
END SUBROUTINE DRAW

!__________________________________________ 
!   	'0 0 0 setrgbcolor'    black
!     '0 0 1 setrgbcolor'    blue
!     '1 0 0 setrgbcolor'    red
!     '0 1 0 setrgbcolor'    green

INTEGER FUNCTION site_num(X,Y)
  implicit none
  integer :: x,y
  site_num=y*Ly+x+1
  return
END FUNCTION site_num
!====================================================================
!====================================================================

!!================================================================
!!================= READ/WRITE CONFIGURATIONS ====================
!!================================================================

!SUBROUTINE read_first_Gamma
  !implicit none
  !integer :: ix, iy, ityp, iomega1, iomega2

  !open(104, status='old', file=trim(title1)//"_first_order_Gamma.dat")
  !read(104, *)  GamNormWeight
  !do ix = 0, dLx
  !do iy = 0, dLy
    !do ityp = 1, ntypGa
      !do iomega1 = -MxOmegaDiag, MxOmegaDiag
        !do iomega2 = 0, MxOmegaDiag
          !read(104, *) Gam1MR(ix, iy, ityp, iomega1, iomega2) 
        !enddo
      !enddo
    !enddo
  !enddo
  !enddo
  !close(104)
!END SUBROUTINE read_first_Gamma



!SUBROUTINE write_first_Gamma
  !implicit none
  !integer :: ix, iy, ityp, iomega1, iomega2

  !open(104, status='replace', file=trim(title1)//"_first_order_Gamma.dat")
  !write(104, *)  GamNormWeight
  !do ix = 0, dLx
  !do iy = 0, dLy
    !do ityp = 1, ntypGa
      !do iomega1 = -MxOmegaDiag, MxOmegaDiag
        !do iomega2 = 0, MxOmegaDiag
          !write(104, *) Gam1MR(ix, iy, ityp, iomega1, iomega2) 
        !enddo
      !enddo
    !enddo
  !enddo
  !enddo
  !close(104)
!END SUBROUTINE write_first_Gamma


SUBROUTINE read_GWGamma
  implicit none
  integer :: ix, iy, ityp, it1, it2, ios

  open(100, status="old", file=trim(title1)//"_G_file.dat",iostat=ios)
  open(101, status="old", file=trim(title1)//"_W_file.dat",iostat=ios)
  open(102, status="old", file=trim(title1)//"_Gamma_file.dat",iostat=ios)

  if(.not. (ios==0)) then
    write(*,*) "You have to run self consistent loop first!"
    stop
  endif

  do it1 = 0, MxT-1
    do ityp = 1, NTypeG
      read(100, *) G(ityp, it1)
    enddo
  enddo

  do it1 = 0, MxT-1
    do iy = 0, Ly-1
      do ix = 0, Lx-1
        do ityp = 1, NTypeW
          read(101, *) W(ityp, ix, iy, it1)
        enddo
      enddo
    enddo
  enddo

  do it2 = 0, MxT-1
    do it1 = 0, MxT-1
      do iy = 0, Ly-1
        do ix = 0, Lx-1
          do ityp = 1, NTypeG
            read(102, *) Gam(ityp, ix, iy, it1, it2)
          enddo
        enddo
      enddo
    enddo
  enddo



  close(100)
  close(101)
  close(102)
  return
END SUBROUTINE read_GWGamma

SUBROUTINE write_GWGamma
  implicit none
  integer :: ix, iy, ityp
  integer :: it1, it2
  character*26 ich

  open(100, status="replace", file=trim(title1)//"_G_file.dat")
  open(101, status="replace", file=trim(title1)//"_W_file.dat")
  open(102, status="replace", file=trim(title1)//"_Gamma_file.dat")

  do it1 = 0, MxT-1
    do ityp = 1, NTypeG
      write(100, *) G(ityp, it1)
    enddo
  enddo

  do it1 = 0, MxT-1
    do iy = 0, Ly-1
      do ix = 0, Lx-1
        do ityp = 1, NTypeW
          write(101, *) W(ityp, ix, iy, it1)
        enddo
      enddo
    enddo
  enddo

  do it2 = 0, MxT-1
    do it1 = 0, MxT-1
      do iy = 0, Ly-1
        do ix = 0, Lx-1
          do ityp = 1, NTypeG
            write(102, *) Gam(ityp, ix, iy, it1, it2)
          enddo
        enddo
      enddo
    enddo
  enddo

  close(100)
  close(101)
  close(102)
  return
END SUBROUTINE write_GWGamma




SUBROUTINE write_monte_carlo_data
  implicit none
  integer :: iorder, itopo, ix, iy, ityp, it1, it2

  !=========== write into files =========================================
  open(104, status="replace", &
    & file=trim(title3)//"_monte_carlo_data.bin.dat",form="binary")

  write(104) imc, GamNorm, GamNormWeight
  do it2 = 0, MxT-1
    do it1 = 0, MxT-1
      do iy = 0, Ly-1
        do ix = 0, Lx-1
          do ityp = 1, NtypeGam/2
            do itopo = 0, 1
              do iorder = 0, MCOrder
                write(104)  GamMC(iorder, itopo, ityp, ix, iy, it1, it2)
                write(104)  GamSqMC(iorder, itopo, ityp, ix, iy, it1, it2)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
  !=========  write on the screen ========================================

  
  close(104)
END SUBROUTINE write_monte_carlo_data

SUBROUTINE write_monte_carlo_conf
  implicit none
  integer :: i

  open(103, status="replace", &
   & file=trim(title3)//"_monte_carlo_conf.bin.dat",form="binary")

  write(103)  Order, NGLn, NWLn, NVertex, MeasureGam, SignFermiLoop, IsWormPresent 
  write(103)  Ira, Masha, SpinMasha, kMasha
  write(103)  TailLn, TailVertex
  write(103)  NextLn(:)
  write(103)  NextVertex(:)
  write(103)  WeightCurrent, Phase
  write(103)  StatusLn(:)
  write(103)  GLnKey2Value(1:NGLn), WLnKey2Value(1:NWLn)
  write(103)  LnValue2Key(:)
  write(103)  StatusVertex(:)
  write(103)  VertexKey2Value(:)
  write(103)  VertexValue2Key(:)
  write(103)  KindLn(:)
  write(103)  TypeLn(:)
  write(103)  SpInVertex(:, :)
  write(103)  TypeVertex(:)
  write(103)  kLn(:)
  write(103)  Hash4G(:) 
  write(103)  Hash4W(:)
  write(103)  GRVertex(1, :), GRVertex(2, :)
  write(103)  WRVertex(1, :), WRVertex(2, :)
  write(103)  TVertex(1, :), TVertex(2, :), TVertex(3, :)
  write(103)  DirecVertex(:)             
  write(103)  WeightLn(:) 
  write(103)  WeightVertex(:) 

  do i = 1, 2
    write(103)  NeighLn(i,:)
  enddo

  do i = 1, 3
    write(103)  NeighVertex(i,:) 
  enddo

  close(103)
END SUBROUTINE write_monte_carlo_conf


SUBROUTINE read_monte_carlo_data
  implicit none
  integer :: iorder, ix, iy, ityp, it1, it2, itopo

  open(105, status="old", file=trim(title)//"_monte_carlo_data.bin.dat",form="binary")

  read(105) imc, GamNorm, GamNormWeight
  do it2 = 0, MxT-1
    do it1 = 0, MxT-1
      do iy = 0, Ly-1
        do ix = 0, Lx-1
          do ityp = 1, NtypeGam/2
            do itopo = 0, 1
              do iorder = 0, MCOrder
                read(105)  GamMC(iorder, itopo, ityp, ix, iy, it1, it2)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
  

  close(105)
END SUBROUTINE read_monte_carlo_data



SUBROUTINE read_monte_carlo_conf
  implicit none
  integer :: i

  open(106, status="old", file=trim(title)//"_monte_carlo_conf.bin.dat",form="binary")

  read(106)  Order, NGLn, NWLn, NVertex, MeasureGam, SignFermiLoop, IsWormPresent 
  read(106)  Ira, Masha, SpinMasha, kMasha
  read(106)  TailLn, TailVertex
  read(106)  NextLn(:)
  read(106)  NextVertex(:)
  read(106)  WeightCurrent, Phase
  read(106)  StatusLn(:)
  read(106)  GLnKey2Value(1:NGLn), WLnKey2Value(1:NWLn)
  read(106)  LnValue2Key(:)
  read(106)  StatusVertex(:)
  read(106)  VertexKey2Value(:)
  read(106)  VertexValue2Key(:)
  read(106)  KindLn(:)
  read(106)  TypeLn(:)
  read(106)  SpInVertex(:, :)
  read(106)  TypeVertex(:)
  read(106)  kLn(:)
  read(106)  Hash4G(:) 
  read(106)  Hash4W(:)
  read(106)  GRVertex(1, :), GRVertex(2, :)
  read(106)  WRVertex(1, :), WRVertex(2, :)
  read(106)  TVertex(1, :), TVertex(2, :), TVertex(3, :)
  read(106)  DirecVertex(:)             
  read(106)  WeightLn(:) 
  read(106)  WeightVertex(:) 

  do i = 1, 2
    read(106)  NeighLn(i,:)
  enddo

  do i = 1, 3
    read(106)  NeighVertex(i,:) 
  enddo

  close(106)
END SUBROUTINE read_monte_carlo_conf
!!================================================================
!!================================================================
!!================================================================





  
!!================================================================
!!========== PRINT OUT THE DATA FILES ============================
!!================================================================

SUBROUTINE output_Quantities
  implicit none
  integer :: ifile, iorder
  integer :: omega, omega1, i, ityp, dx, dy
  integer :: omega2
  integer :: it1, it2, it
  double precision :: GGI
  double precision :: GaR1, GaR2, GaR3
  double precision :: WWR1, WWR2, WWR3, WWR4

  open(15, file=trim(title1)//"_Chi.dat")
  !open(16, file=trim(title1)//"_Sigma.dat")

  !open(11, access="append", file=trim(title1)//"_G.dat")
  !open(12, access="append", file=trim(title1)//"_G_omega_0.dat")
  !open(13, access="append", file=trim(title1)//"_W.dat")
  !open(14, access="append", file=trim(title1)//"_W_omega_0.dat")
  !open(15, access="append", file=trim(title1)//"_Chi.dat")
  !open(16, access="append", file=trim(title1)//"_Chi_omega_0.dat")
  !open(17, access="append", file=trim(title1)//"_Gamma.dat")
  !open(18, access="append", file=trim(title1)//"_Gamma_omega_0.dat")
  !open(19, access="append", file=trim(title1)//"_Sigma.dat")
  !open(20, access="append", file=trim(title1)//"_Sigma_omega_0.dat")
  !open(21, access="append", file=trim(title1)//"_Pi.dat")
  !open(22, access="append", file=trim(title1)//"_Pi_omega_0.dat")
  !open(23, access="append", file=trim(title1)//"_Gamma_matrix.dat")

  do it = 0, MxT-1
    do dy = 0, dLy
      do dx = 0, dLx
        write(15, *) dx, dy, it, Chi(dx, dy, it)
      enddo
    enddo
  enddo

  close(15)

  !do it = 0, MxT-1
    !write(16, *) it, real(Sigma(it)), dimag(Sigma(it))
  !enddo

  !close(16)

  !do ifile = 11, 22
    !write(ifile, *) "=============================================================="
    !write(ifile, *) "Beta", Beta, "Lx, Ly", Lx, Ly, "Order", MCOrder, "Seed",Seed
    !write(ifile, *) ime, GamNorm, GamNormWeight
  !enddo

  !do omega = -MxOmegaG2, MxOmegaG2
    !GGI = weight_G(omega, 2)
    !write(11, *) omega, GGI
    !if(omega==0)      write(12, *) GGI
  !enddo

  !close(11)
  !close(12)

  !do ityp = 1, 2, 5
    !do omega = -MxOmegaW2, MxOmegaW2
      !WWR1 = weight_W( 0, 0, omega, 1)
      !WWR2 = weight_W( 0, 1, omega, 1)
      !WWR3 = weight_W( 1, 1, omega, 1)
      !write(13, *) ityp, omega, WWR1, WWR2, WWR3

      !if(omega==0)      write(14, *) ityp, WWR1, WWR2, WWR3
    !enddo
  !enddo

  !close(13)
  !close(14)

  !do dx = 0, dLx
    !do dy = 0, dLy
      !do omega = -MxOmegaChi, MxOmegaChi
        !write(15, *) Beta, dx, dy, omega, -3.d0*trChiR(dx, dy, omega)/Beta
        !if(omega==0)      write(16, *) Beta, dx, dy, -3.d0*trChiR(dx, dy, 0)/Beta
      !enddo
    !enddo
  !enddo

  !close(15)
  !close(16)

  !do omega1 = -MxOmegaDiag, MxOmegaDiag
    !do omega2 = -MxOmegaDiag, MxOmegaDiag
      !GaR1 = weight_Gamma(1, 0, omega1, omega2, 1)
      !write(23, '(f14.8)', advance='no') GaR1 
    !enddo

    !write(23, *)
    !GaR1 = weight_Gamma(0, 0, omega1, omega1, 1)
    !GaR2 = weight_Gamma(1, 0, omega1, omega1, 1)
    !GaR3 = weight_Gamma(1, 1, omega1, omega1, 1)
    !write(17, *) GaR1-1.d0, GaR2, GaR3
    !if(omega1==0)      write(18, *) GaR1-1.d0, GaR2, GaR3
  !enddo

  !close(17)
  !close(18)
  !close(23)

  !do iorder = 0, MCOrder
    !write(19, *) "Order:", iorder
    !do omega = -MxOmegaSigma, MxOmegaSigma
      !write(19, *) Beta, omega, SigmaI(iorder, 1, omega)
      !if(omega==0)      write(20, *) Beta, iorder, SigmaI(iorder, 1, 0)
    !enddo
  !enddo

  !close(19)
  !close(20)

  !do ityp = 1, ntypPi
    !do dx = 0, dLx
      !do dy = 0, dLy
        !do omega = -MxOmegaChi, MxOmegaChi
          !write(21, *) Beta, dx, dy, omega, PiR(ityp, dx, dy, omega)
          !if(omega==0)      write(22, *) Beta, dx, dy, PiR(ityp, dx, dy, 0)
        !enddo
      !enddo
    !enddo
  !enddo

  !close(21)
  !close(22)
END SUBROUTINE output_Quantities



!SUBROUTINE output_prob_MC
  !implicit none
  !integer :: omega, omega1, omega2, ityp, iorder, iloop
  !integer :: iconf, i, j, iobs, iomega
  !double precision :: Ga0R,GaMR1, GaCR1, GaMR2, GaMR3
  !double precision :: Norm

  !open(33, access="append", file=trim(title3)//"_prob_MC.dat")
  !!open(34, access="append", file=trim(title3)//"_over_complete_1_2.dat")
  !!open(35, access="append", file=trim(title3)//"_over_complete_3_4.dat")
  !!open(36, access="append", file=trim(title3)//"_over_complete_7_8.dat")
  !!open(37, access="append", file=trim(title3)//"_over_complete_9_10.dat")

  !write(33, *) "============================================"
  !write(33, *) "Beta", Beta, "Lx, Ly", Lx, Ly, "Order", MCOrder, "Seed",Seed
  !write(33, *) " 1: create worm along wline"
  !write(33, *) " 2: delete worm along wline"
  !write(33, *) " 3: create worm along gline"
  !write(33, *) " 4: delete worm along gline"
  !write(33, *) " 5: move worm along wline"
  !write(33, *) " 6: move worm along gline"
  !write(33, *) " 7: add interaction"
  !write(33, *) " 8: remove interaction"
  !write(33, *) " 9: add interaction cross"
  !write(33, *) "10: remove interaction cross"
  !write(33, *) "11: reconnect"
  !write(33, *) "12: shift gline in space"
  !write(33, *) "13: shift wline in space"
  !write(33, *) "14: change Gamma type"
  !write(33, *) "15: move measuring index"

  !write(33, *) "MC", imc, "Measure", ime
  !do iorder = 0, MCOrder
    !write(33, *) "Order", iorder
    !do i = 1, Nupdate
      !if(ProbProp(iorder, i)/=0.d0) then
        !write(33, *) i, ProbProp(iorder, i), ProbAcc(iorder, i), ProbAcc(iorder, i)/ProbProp(iorder, i)
      !endif
    !enddo
    !write(33, *)
  !enddo
  !write(33, *)

  !write(33, *)  "Physical configurations for different orders"
  !do iorder = 0, MCOrder
    !!if(iorder ==0) then
      !!write(33, *) iorder, GamOrder(iorder), GamWormOrder(iorder), Gam0Bubble
    !!else
      !!write(33, *) iorder, GamOrder(iorder), GamWormOrder(iorder)
    !!endif
    !write(33, *) iorder, GamOrder(iorder), GamWormOrder(iorder)
  !enddo
  !write(33, *) 

  !write(33, *)  "Average Weight Ratios"
  !do iorder = 0, MCOrder
    !write(33, *) iorder, 2, AveWeightRatio(iorder, 1)/ProbProp(iorder,2)
    !write(33, *) iorder, 4, AveWeightRatio(iorder, 2)/ProbProp(iorder,4)
  !enddo
  !write(33, *) 

  !!write(33, *)  "Histogram for omega on W"
  !!do omega = -MxOmega, MxOmega
    !!if(HistoOmegaW(omega)/=0.d0) then
      !!write(33, *) omega, HistoOmegaW(omega)
    !!endif
  !!enddo

  !!if(ProbAcc(1)/=0.d0 .and. ProbAcc(2)/=0.d0) then
    !!write(34, *) ime, imc, (ProbAcc(1)-ProbAcc(2))/sqrt(ProbAcc(1)+ProbAcc(2))
  !!endif

  !!if(ProbAcc(3)/=0.d0 .and. ProbAcc(4)/=0.d0) then
    !!write(35, *) ime, imc, (ProbAcc(3)-ProbAcc(4))/sqrt(ProbAcc(3)+ProbAcc(4))
  !!endif

  !!if(ProbAcc(7)/=0.d0 .and. ProbAcc(8)/=0.d0) then
    !!write(36, *) ime, imc, (ProbAcc(7)-ProbAcc(8))/sqrt(ProbAcc(7)+ProbAcc(8))
  !!endif

  !!if(ProbAcc(9)/=0.d0 .and. ProbAcc(10)/=0.d0) then
    !!write(37, *) ime, imc, (ProbAcc(9)-ProbAcc(10))/sqrt(ProbAcc(10)+ProbAcc(9))
  !!endif

  !close(33)
  !!close(34)
  !!close(35)
  !!close(36)
  !!close(37)
!END SUBROUTINE output_prob_MC


!SUBROUTINE output_GamMC
  !implicit none
  !integer :: iconf, i, j, iorder, iomega, iomega2
  !double precision :: err, percenterr, norm
  !double precision :: gam1, gam2, gamn

  !open(34, access="append", file=trim(title3)//"_Gamma_MC_matrix.dat")
  !open(35, access="append", file=trim(title3)//"_Gamma_MC.dat")
  !open(36, access="append", file=trim(title3)//"_Gamma_order2_MC.dat")

  !write(34, *) "============================================"
  !write(34, *) "Beta", Beta, "Lx, Ly", Lx, Ly, "Order", MCOrder, "Seed",Seed
  !write(34, *) imc, ime, GamNormWeight, GamNorm

  !write(35, *) "============================================"
  !write(35, *) "Beta", Beta, "Lx, Ly", Lx, Ly, "Order", MCOrder, "Seed",Seed
  !write(35, *) imc, ime, GamNormWeight, GamNorm

  !write(36, *) "============================================"
  !write(36, *) "Beta", Beta, "Lx, Ly", Lx, Ly, "Order", MCOrder, "Seed",Seed
  !write(36, *) imc, ime, GamNormWeight, GamNorm

  !norm = GamNormWeight*ime/GamNorm

  !write(34, *) "Order 1, dx=0, dy=0"
  !do iomega = -MxOmegaDiag, MxOmegaDiag
    !do iomega2 = -MxOmegaDiag, MxOmegaDiag
      !gam1 = GamMC(1,0,1, 0, 0, iomega, iomega2)/ime
      !gamn = gam1*norm
      !write(34, '(f14.8)', advance='no') gamn
    !enddo
    !write(34, *) 
  !enddo
  !write(34, *) 

  !write(34, *) "Order 2, dx=0, dy=0"
  !do iomega = -MxOmegaDiag, MxOmegaDiag
    !do iomega2 = -MxOmegaDiag, MxOmegaDiag
      !gam1 = GamMC(2, 0, 1, 0, 0, iomega, iomega2)/ime
      !gamn = gam1*norm
      !write(34, '(f14.8)', advance='no') gamn
    !enddo
    !write(34, *) 
  !enddo
  !write(34, *) 

  !write(34, *) "Order 2, dx=1, dy=0"
  !do iomega = -MxOmegaDiag, MxOmegaDiag
    !do iomega2 = -MxOmegaDiag, MxOmegaDiag
      !gam1 = GamMC(2, 1, 1, 1, 0, iomega, iomega2)/ime
      !gamn = gam1*norm
      !write(34, '(f14.8)', advance='no') gamn
    !enddo
    !write(34, *) 
  !enddo
  !write(34, *) 

  !write(34, *) "Order 3, dx=0, dy=0"
  !do iomega = -MxOmegaDiag, MxOmegaDiag
    !do iomega2 = -MxOmegaDiag, MxOmegaDiag
      !gam1 = GamMC(3, 0, 1, 0, 0, iomega, iomega2)/ime
      !gamn = gam1*norm
      !write(34, '(f14.8)', advance='no') gamn
    !enddo
    !write(34, *) 
  !enddo
  !write(34, *) 

  !write(34, *) "Order 3, dx=1, dy=0"
  !do iomega = -MxOmegaDiag, MxOmegaDiag
    !do iomega2 = -MxOmegaDiag, MxOmegaDiag
      !gam1 = GamMC(3, 1, 1, 1, 0, iomega, iomega2)/ime
      !gamn = gam1*norm
      !write(34, '(f14.8)', advance='no') gamn
    !enddo
    !write(34, *) 
  !enddo
  !write(34, *) 

  !write(34, *) "total, dx=1, dy=0"
  !do iomega = -MxOmegaDiag, MxOmegaDiag
    !do iomega2 = -MxOmegaDiag, MxOmegaDiag
      !gam1 = SUM(GamMC(0:MCOrder, 1, 1, 1, 0, iomega, iomega2))/ime
      !gamn = gam1*norm
      !write(34, '(f14.8)', advance='no') gamn
    !enddo
    !write(34, *) 
  !enddo
  !write(34, *) 

  !do iorder = 1, MCOrder
    !write(35, *) "Order", iorder
    !write(35, *) "dx = 0, dy = 0"
    !do iomega = -MxOmegaDiag, MxOmegaDiag
      !gam1 = GamMC(iorder, 0, 1, 0, 0, iomega, iomega)/ime
      !gam2 = GamSqMC(iorder,0, 1, 0, 0, iomega, iomega)/ime
      !err = sqrt(gam2-gam1**2.d0)/sqrt(ime-1)
      !if(abs(gam1)<1.d-30) then
        !percenterr = 0.d0
      !else
        !percenterr = err/abs(gam1)
      !endif
      !gamn = gam1*norm
      !write(35, *) iomega, gamn, percenterr*gamn, percenterr
    !enddo

    !write(35, *) "dx = 1, dy = 0, omega1=omega2"
    !do iomega = -MxOmegaDiag, MxOmegaDiag
      !gam1 = GamMC(iorder, 1, 1, 1, 0, iomega, iomega)/ime
      !gam2 = GamSqMC(iorder,1, 1, 1, 0, iomega, iomega)/ime
      !err = sqrt(gam2-gam1**2.d0)/sqrt(ime-1)
      !if(abs(gam1)<1.d-30) then
        !percenterr = 0.d0
      !else
        !percenterr = err/abs(gam1)
      !endif
      !gamn = gam1*norm
      !write(35, *) iomega, gamn, percenterr*gamn, percenterr
    !enddo

    !write(35, *) "dx = 1, dy = 0, omega2=-1"
    !do iomega = -MxOmegaDiag, MxOmegaDiag
      !gam1 = GamMC(iorder, 1, 1, 1, 0, iomega, -1)/ime
      !gam2 = GamSqMC(iorder,1, 1, 1, 0, iomega, -1)/ime
      !err = sqrt(gam2-gam1**2.d0)/sqrt(ime-1)
      !if(abs(gam1)<1.d-30) then
        !percenterr = 0.d0
      !else
        !percenterr = err/abs(gam1)
      !endif
      !gamn = gam1*norm
      !write(35, *) iomega, gamn, percenterr*gamn, percenterr
    !enddo
    !write(35, *) 
    !write(35, *) 
  !enddo

  !write(36, *) "topo 1:"
  !do iomega = -MxOmegaDiag, MxOmegaDiag
    !gam1 = Gam2Topo(1,iomega)/ime
    !gamn = gam1*norm
    !write(36, *) iomega, gamn
  !enddo

  !write(36, *) "topo 2:"
  !do iomega = -MxOmegaDiag, MxOmegaDiag
    !gam1 = Gam2Topo(2,iomega)/ime
    !gamn = gam1*norm
    !write(36, *) iomega, gamn
  !enddo


  !!write(36, *) "Order 2, typ=1, iloop = 2"
  !!do iomega = -MxOmegaDiag, MxOmegaDiag
    !!do iomega2 = -MxOmegaDiag, MxOmegaDiag
      !!gamn = Gam2MR(1, iomega, iomega2)
      !!write(36, '(f14.8)', advance='no') gamn
    !!enddo
    !!write(36, *) 
  !!enddo


  !close(34)
  !close(35)
  !close(36)
!END SUBROUTINE output_GamMC

SUBROUTINE write_monte_carlo_test
  implicit none
  integer :: iorder
  !open(11, access="append", file=trim(title3)//"_test.dat")
  !do iorder = 0, MCOrder
    !write(11, *) iorder, "conf(total)/conf(all spin up)", TestData(iorder)/TestData(MCOrder+iorder+1)
  !enddo
  !write(11, *)

  !do iorder = 1, MCOrder
    !write(11, *) iorder, iorder-1, TestData(MCOrder+iorder+1)/TestData(MCOrder+iorder)
  !enddo
  !write(11, *) 
  !close(11)
  return
END SUBROUTINE

!!================================================================
!!================================================================
!!================================================================

