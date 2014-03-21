

!====================== the 1st order of Gamma ======================
!===================== need to change the type of Gamma =============
SUBROUTINE calculate_Gam1
  implicit none
  integer :: dx, dy, t1, t2, dir, ityp
  integer :: typ(3), gintyp(3), gouttyp(3), wtyp(3), ga1typ(3), ga2typ(3), ga3typ(3)
  complex*16 :: Gin, Gout, iW, Gam1, Gam2, Gam3
  complex*16 :: weight

  GamOrder1(:,:,:) = (0.d0, 0.d0)

  !================== bare Gamma ===============================
  typ(1) = 1
  gintyp(1) = 1;      gouttyp(1) = 1
  ga1typ(1) = 1;      ga2typ(1)  = 1;     ga3typ(1) = 1
  wtyp(1)   = 1

  typ(2) = 3
  gintyp(2) = 2;      gouttyp(2) = 2
  ga1typ(2) = 2;      ga2typ(2)  = 5;     ga3typ(2) = 6
  wtyp(2)   = 5


  typ(3) = 5
  gintyp(3) = 1;      gouttyp(3) = 2
  ga1typ(3) = 5;      ga2typ(3)  = 1;     ga3typ(3) = 2
  wtyp(3)   = 3

  do t2 = 0, MxT-1
    do t1 = 0, MxT-1
      do ityp = 1, 3

        Gin = weight_G(gintyp(ityp), -t1)
        Gout = weight_G(gouttyp(ityp), t2)
        Gam1 = weight_Gam0(ga1typ(ityp), 0, 0)
        Gam2 = weight_Gam0(ga2typ(ityp), 0, 0) 
        Gam3 = weight_Gam0(ga3typ(ityp), 0, 0) 
        iW = weight_W(wtyp(ityp), 0, 0, t2-t1)

        weight = Gin *Gout *iW *Gam1 *Gam2 *Gam3

        GamOrder1(typ(ityp),t1,t2) = d_times_cd(Beta/(dble(MxT)**2.d0), weight)
        GamOrder1(typ(ityp)+1,t1,t2) = GamOrder1(typ(ityp), t1, t2)
      enddo
    enddo
  enddo


END SUBROUTINE calculate_Gam1
!====================================================================




