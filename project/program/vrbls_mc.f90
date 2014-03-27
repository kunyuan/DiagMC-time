!********************************************************************
MODULE vrbls_mc
  IMPLICIT NONE

  !======================== Parameters ====================================
  double precision, parameter :: Pi    = 3.14159265358979323846d0
  double precision, parameter :: Pi2   = 6.2831853071795865d0
  double precision, parameter :: Pi4   = 2.0d0*Pi2
  double precision, parameter :: tm32  = 1.d0/(2.d0**32)
  double precision, parameter :: macheps = epsilon(0d0)
  integer, parameter          :: Mxint = 2147483647
  integer, parameter          :: Mnint =-2147483647

  integer, parameter :: D = 2                       ! 2-dimensional system
  integer, parameter :: MxLx  = 8, MxLy = 8       ! the largest system
  !                                                   size,MxLx>MxLy
  integer, parameter :: MxVol = MxLx**D             ! the maximum system volume
  integer, parameter :: MxT   = 128                  ! the maximum number of time segments
  integer, parameter :: MxK   = 1000000             ! the maximum momentum
  !integer, parameter :: MxOmegaBasis = 2048         ! the maximum omega used in basis

  double precision, parameter :: MxError = 0.25     ! the maximum error for MC
  integer, parameter          :: MxNblck = 1000000    ! the maximum memory blocks in MC simulations

  integer, parameter :: MxOrder =  10               ! the maximum order of the diagram
  integer, parameter :: MxNLn  = 3*MxOrder+3        ! the maximum number of lines
  integer, parameter :: MxNGLn = 2*MxOrder+2        ! the maximum number of glines
  integer, parameter :: MxNWLn = MxOrder+1          ! the maximum number of wlines
  integer, parameter :: MxNVertex = 2*MxOrder+2     ! the maximum number of vertexes

  !=======================================================================

  character*100 :: title1
  character*100 :: title2
  character*100 :: title_loop
  character*100 :: title_mc
  character*100 :: title
  character*100 :: logstr

  integer:: file_version
  integer:: mc_version

  !======================== Input parameter ==============================
  integer          ::  Lx, Ly, L(2), Vol                ! System size
  double precision  ::  logL(2)
  integer          ::  dLx, dLy
  double precision ::  Jcp                        ! interaction
  double precision ::  Mu(2)                      ! Chem. potential for spin down & up
  double precision ::  Beta                       ! inverse temperature
  integer          ::  MCOrder                    ! the max order for Gamma in MC
  logical          ::  CheckG, CheckW, CheckGam ! if turn on the irreducibility check
  !=======================================================================



  !======================== Spin type for different quantities ===============
  integer, parameter :: NTypeG = 2                   ! types of G  
  !----------------------------------------------------------------------------
  !----- type=1: up;   type=2: down -------------------------------------------
  !----------------------------------------------------------------------------
  integer, parameter :: NTypeW = 6                   ! types of W  
  !----------------------------------------------------------------------------
  !               (out)           (out)
  !               b \      W       / d
  !                  ==============
  !               a /              \ c
  !               (in)            (in)
  !----- type = 1:  a = up;   c = up;   b = up;    d = up  -------------------
  !----- type = 2:  a = down; c = down; b = down;  d = down  -----------------
  !----- type = 3:  a = up;   c = down; b = up;    d = down  -----------------
  !----- type = 4:  a = down; c = up;   b = down;  d = up  -------------------
  !----- type = 5:  a = up;   c = down; b = down;  d = up  -------------------
  !----- type = 6:  a = down; c = up;   b = up;    d = down  -----------------
  !----------------------------------------------------------------------------
  integer, parameter :: NTypeGam = 6                  ! types of Gamma  
  !----------------------------------------------------------------------------
  !                       (out)
  !                 (out) / b
  !                   d /||
  !                ====  ||
  !                   c \||
  !                  (in) \ a
  !                       (in)
  !----- type = 1:  a = up;   d = up;   b = up;    c = up  -------------------
  !----- type = 2:  a = down; d = down; b = down;  c = down  -----------------
  !----- type = 3:  a = up;   d = down; b = up;    c = down  -----------------
  !----- type = 4:  a = down; d = up;   b = down;  c = up  -------------------
  !----- type = 5:  a = up;   d = down; b = down;  c = up  -------------------
  !----- type = 6:  a = down; d = up;   b = up;    c = down  -----------------
  !----------------------------------------------------------------------------
  integer, parameter :: NTypeChi = 4                  ! types of Chi 
  !----------------------------------------------------------------------------
  !          a-----         ----b
  !          |                  |
  !          |                  |
  !          a-----         ----b
  !----- type = 1:  a = up;   b = up   ----------------------------------------
  !----- type = 2:  a = down; b = down   --------------------------------------
  !----- type = 3:  a = up;   b = down  ---------------------------------------
  !----- type = 4:  a = down; b = up  -----------------------------------------
  !----------------------------------------------------------------------------
  !=======================================================================


  !========================= Self-consistent loop ========================
  !============== unfinished =============================================
  complex(kind=8) :: G(NtypeG, 0:MxT-1)
  complex(kind=8) :: G0F(0:MxT-1)
  complex(kind=8) :: Sigma(0:MxT-1)

  complex(kind=8), allocatable :: W(:,:,:,:)
  complex(kind=8), allocatable :: Gam(:,:,:,:,:)

  complex(kind=8), allocatable :: W0PF(:,:,:)
  complex(kind=8), allocatable :: Gam0PF(:,:,:,:)

  complex(kind=8), allocatable :: Polar(:,:,:)
  complex(kind=8), allocatable :: Chi(:,:,:)
  !=======================================================================

  !======================== analytic integration =========================
  complex(kind=8) :: GamOrder1(NTypeGam, 0:MxT-1, 0:MxT-1) 

  !====================== MC Simulation ==================================
  complex*16 :: GamNorm, GamNormWeight           ! the weight of the normalization diagram
  complex*16, allocatable :: GamMC(:,:,:,:,:,:)      ! the measurement of Gamma in MC
  double precision, allocatable :: ReGamSqMC(:,:,:,:,:,:)      ! the measurement of Gamma in MC
  double precision, allocatable :: ImGamSqMC(:,:,:,:,:,:)      ! the measurement of Gamma in MC

  double precision :: GamOrder(0:MxOrder)              ! the configuration number of different orders
  double precision :: GamWormOrder(0:MxOrder)          ! the configuration number in whole section

  double precision :: WeightCurrent            ! the current weight of the configuration
  double precision :: CoefOfWorm
  double precision :: CoefOfWeight(0:MxOrder)  ! the coeffecients for different orders and worm section
  double precision :: CoefOfSymmetry(0:MxLx, 0:MxLy)
  double precision :: SpatialWeight(1:2,0:MxLx-1)  ! 1, fox x direction 2, for y direction;
                                                 !!! attention: dx,dy=0~L!!!!

  !------------- MC steps -----------------------------------
  integer          :: iupdate           ! the update number
  double precision :: imc               ! imc: the MC step
  double precision :: imeasure
  integer          :: NSamp             ! # total MC steps
  integer          :: NToss             ! # MC steps for toss
  integer          :: NStep             ! # MC steps for one measurement
  logical          :: IsToss
  logical          :: IsForever

  !------------ basic variables for a diagram --------------------------
  integer          :: Order             ! order of the simulating diagram 
  COMPLEX*16       :: Phase             ! phase of the present diagram
  integer          :: NGLn, NWLn, NVertex ! number of glines, wlines, gamma
  integer          :: MeasureGam        ! measuring vertex
  integer          :: SignFermiLoop     ! the sign of Fermi loops(1: even; -1: odd)
  logical          :: IsWormPresent     ! to represent if worm is present

  !------------ Ira and Masha variables for a diagram -------------------
  integer          :: Ira, Masha        ! vertex number of Ira and Masha
  integer          :: kMasha            ! excess momentum for Masha (-kMasha for Ira)
  integer          :: SpinMasha         ! excess spin for Masha (-SpinMasha for Ira)
  double precision :: WeightWorm        ! weight for worm sector

  !----------- data structure for lines -------------------------------------
  integer, dimension(MxNLn) :: kLn                  ! momentum of a line
  integer, dimension(MxNLn) :: KindLn               ! kind of a line: 1 Glines; 2 Wlines
  integer, dimension(MxNLn) :: TypeLn               ! type of a line: 1-2 Glines; 1-6 Wlines
  integer, dimension(2,MxNLn) :: NeighLn            ! 1: begin gamma; 2: end gamma 
  integer, dimension(MxNLn) :: StatusLn             ! 2:I&M; 0:normal; 1: measure; -1: empty
  integer, dimension(MxNLn) :: IsDeltaLn            ! 0: W, not delta; 1: W, delta; -1: G
  COMPLEX*16, dimension(MxNLn) :: WeightLn          ! weight of glines and wlines
  integer, dimension(MxNLn)  :: NextLn
  integer                    :: TailLn      ! for add or delete use

  integer, dimension(MxNGLn) :: GLnKey2Value        
  integer, dimension(MxNWLn) :: WLnKey2Value       
  integer, dimension(MxNLn) ::  LnValue2Key       
                                                   

  !----------- data structure for vertexes -------------------------------------
  integer, dimension(2, MxNVertex) :: GRVertex      ! gline sites for Gamma
  integer, dimension(2, MxNVertex) :: WRVertex      ! wline sites for Gamma
  double precision, dimension(3, MxNVertex) :: TVertex! times for Gamma
  !  1: outgoing G, 2: ingoing G, 3: W attached
  integer, dimension(MxNVertex) :: DirecVertex              ! Gamma is 1: left, 2: right
  integer, dimension(MxNVertex) :: TypeVertex               ! type of Gamma: 1-6
  integer, dimension(2,MxNVertex) :: SpInVertex             ! type of Gamma inside spin: 1-2
  integer, dimension(3,MxNVertex) :: NeighVertex            ! 1: in gline; 2: out gline; 3: wline
  COMPLEX*16, dimension(MxNVertex) :: WeightVertex          ! weight of  gamma
  integer, dimension(MxNLn) :: IsDeltaVertex                ! 0: not delta; 1: delta
  integer, dimension(MxNVertex):: StatusVertex              ! 2:I&M; 0:normal; 1: measure; -1: empty
  integer, dimension(MxNVertex):: NextVertex
  integer                      :: TailVertex    ! for add or delete use

  integer, dimension(MxNVertex):: VertexKey2Value        
  integer, dimension(MxNVertex):: VertexValue2Key      

  !----------- spin for W and Gamma---------------------------------------------
  integer, dimension(6,6)   :: TypeGam2W
  integer, dimension(2,2,2,2) :: TypeSp2Gam

  !------------ probabilities functions for updates -----------------------------
  integer, parameter :: Nupdate = 18
  double precision :: Pupdate(Nupdate)  ! the probabilities to call different subroutines       
  double precision :: Fupdate(Nupdate)  ! a function of the summation of Pcall

  double precision :: ProbProp(0:MxOrder, Nupdate)
  double precision :: ProbAcc(0:MxOrder, Nupdate)

  !------------ hashtable for irreducibility ------------------------------------
  integer, dimension(-MxK:MxK) ::  Hash4G     ! hashtable for k in G
  integer, dimension(0:MxK) ::     Hash4W     ! hashtable for k in W

  !=======================================================================
  !------------ statistics ------------------------------------------------------
  !-- Observables --------------------------------------------------
  !! THIS IS PROJECT-DEPENDENT 
  logical :: prt
  integer          :: MaxStat
  integer          :: StatNum
  integer, parameter :: NObs = 20              ! Total # observables
  double precision   :: Quan(NObs)                 ! 1st--#quan.  2nd--#block
  double precision   :: Norm(NObs)
  double precision   :: Error(NObs)
  character(len=30),dimension(NObs) :: QuanName

  double precision, allocatable :: ObsRecord(:,:)                 ! 1st--#quan.  2nd--#block
  double precision :: Z_normal
  double precision :: Z_worm
	DOUBLE PRECISION :: amax, tmax, amin, tmin
  !================= Random-number generator =============================
  integer                      :: Seed              ! random-number seed
  integer, parameter           :: mult=32781
  integer, parameter           :: mod2=2796203, mul2=125
  integer, parameter           :: len1=9689,    ifd1=471
  integer, parameter           :: len2=127,     ifd2=30
  integer, dimension(1:len1)   :: inxt1
  integer, dimension(1:len2)   :: inxt2
  integer, dimension(1:len1)   :: ir1
  integer, dimension(1:len2)   :: ir2
  integer                      :: ipnt1, ipnf1
  integer                      :: ipnt2, ipnf2
  integer, parameter           :: mxrn = 10000
  integer, dimension(1:mxrn)   :: irn(mxrn)
  integer                      :: nrannr                 ! random-number counter
  !=======================================================================

  !==================Time-checking =======================================
  character( 8)         :: date
  character(10)         :: time
  character(5 )         :: zone
  integer, dimension(8) :: tval
  double precision      :: t_prev, t_curr, t_elap
  integer               :: h_prev, h_curr
  double precision      :: t_init, t_simu, t_meas, t_toss
  !=======================================================================

  !================ Test =================================================
  double precision,dimension(0:20)  :: TestData
  !=======================================================================

end module vrbls_mc