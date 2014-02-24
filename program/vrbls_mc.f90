!********************************************************************
MODULE vrbls_mc
  IMPLICIT NONE

  !======================== Parameters ====================================
  double precision, parameter :: Pi    = 3.14159265358979323846d0
  double precision, parameter :: Pi2   = 6.2831853071795865d0
  double precision, parameter :: Pi4   = 2.0d0*Pi2
  integer, parameter          :: Mxint = 2147483647
  integer, parameter          :: Mnint =-2147483647

  integer, parameter :: D = 2                       ! 2-dimensional system
  integer, parameter :: MxLx  = 64, MxLy = 64       ! the largest system size
  integer, parameter :: MxVol = MxLx**D             ! the maximum system volume
  integer, parameter :: MxT   = 10000               ! the maximum number of time segments
  integer, parameter :: MxK   = 1000000             ! the maximum momentum

  !integer, parameter :: MxOmega =  600000           ! the maximum omega used in MC
  !integer, parameter :: MxOmegaBasis = 2048         ! the maximum omega used in basis
  !integer, parameter :: MxOmegaDiag = 32            ! the maximum omega of the measured diagram

  double precision, parameter :: MxError = 0.25     ! the maximum error for MC
  integer, parameter          :: MxNblck = 10240    ! the maximum blocks in MC simulations

  integer, parameter :: MxOrder =  10               ! the maximum order of the diagram
  integer, parameter :: MxNLn  = 3*MxOrder+3        ! the maximum number of lines
  integer, parameter :: MxNGLn = 2*MxOrder+2        ! the maximum number of glines
  integer, parameter :: MxNWLn = MxOrder+1          ! the maximum number of wlines
  integer, parameter :: MxNVertex = 2*MxOrder+2     ! the maximum number of vertexes
  !=======================================================================

  character*100 :: title1
  character*100 :: title2
  character*100 :: title3
  character*100 :: title

  integer:: file_version

  !======================== Input parameter ==============================
  integer          ::  Lx, Ly, Vol                ! System size
  double precision ::  Jcp                        ! interaction
  double precision ::  Mu(2)                      ! Chem. potential for spin down & up
  double precision ::  Beta                       ! inverse temperature
  integer          ::  MCOrder                    ! the max order for Gamma in MC
  logical          ::  CheckG, CheckW, CheckGam ! if turn on the irreducibility check
  !=======================================================================

  !=================== Schmidt basis for G and W =========================
  !integer, parameter :: nbasis = 3                  ! number of basises
  !integer, parameter :: nbasisGamma = nbasis**2+2   ! number of basises

  !double precision :: OriginalBasis(nbasis, -MxOmegaBasis:MxOmegaBasis)
  !double precision :: OriginalBasisGamma(nbasisGamma, -MxOmegaBasis:MxOmegaBasis, &
    !& -MxOmegaBasis:MxOmegaBasis)
  !double precision :: ReweightBasis(-MxOmegaBasis:MxOmegaBasis, -MxOmegaBasis:MxOmegaBasis)

  !double precision :: GCoefP(nbasis, nbasis)        ! coeffecients for positive tail of G
  !double precision :: GCoefN(nbasis, nbasis)        ! coeffecients for negative tail of G

  !double precision :: WCoefP(nbasis, nbasis)        ! coeffecients for positive tail of W
  !double precision :: WCoefN(nbasis, nbasis)        ! coeffecients for negative tail of W


  !double precision :: GamGCoefP(nbasis, nbasis)  ! coeffecients for the positive omegaG tail for Gamma
  !double precision :: GamGCoefN(nbasis, nbasis)  ! coeffecients for the negative omegaG tail for Gamma

  !double precision :: GamCoefPPR(nbasisGamma, nbasisGamma)  !coeffecients for the tail of Gamma
  !double precision :: GamCoefPPL(nbasisGamma, nbasisGamma)  !coeffecients for the tail of Gamma
  !double precision :: GamCoefNNR(nbasisGamma, nbasisGamma)  ! coeffecients for tail of Gamma
  !double precision :: GamCoefNNL(nbasisGamma, nbasisGamma)  ! coeffecients for tail of Gamma

  !double precision :: GamCoefPN(nbasisGamma, nbasisGamma)  ! coeffecients for tail of Gamma
  !double precision :: GamCoefNP(nbasisGamma, nbasisGamma)  !coeffecients for the tail of Gamma

  !=======================================================================

  !========================= Self-consistent loop ========================
  !integer, parameter :: MxOmegaG1 = 4               ! maximum blocks of the tabulate
  !integer, parameter :: MxOmegaG2 = 32              ! maximum blocks of the tail

  !integer, parameter :: MxOmegaW1 = 4               ! maximum blocks of the tabulate
  !integer, parameter :: MxOmegaW2 = 32              ! maximum blocks of the tail

  !integer, parameter :: MxOmegaGamG1 = 8            ! maximum blocks of the tabulate
  !integer, parameter :: MxOmegaGamG2 = 32            ! maximum blocks of the tail


  !integer, parameter :: MxOmegaGInt = 128           ! maximum blocks of the tail in the integration
  !integer, parameter :: MxOmegaWInt = 128           ! maximum blocks of the tail in the integration
  !integer, parameter :: MxOmegaGamInt = 128           ! maximum blocks of the tail in the integration
  !integer, parameter :: MxOmegaPiInt = 128          ! maximum blocks of the tail in the integration
  !integer, parameter :: MxOmegaSigmaInt = 2048       ! maximum blocks of the tail in the integration

  !integer, parameter :: MxOmegaChi = 16             ! maximum blocks of the tail
  !integer, parameter :: MxOmegaSigma = 1024         ! maximum blocks of the tail


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
  integer, parameter :: NTypePi = 4                  ! types of Pi 
  !----------------------------------------------------------------------------
  !          a------------------b
  !          |                  |
  !          |                  |
  !          a------------------b
  !----- type = 1:  a = up;   b = up   ----------------------------------------
  !----- type = 2:  a = down; b = down   --------------------------------------
  !----- type = 3:  a = up;   b = down  ---------------------------------------
  !----- type = 4:  a = down; b = up  -----------------------------------------
  !----------------------------------------------------------------------------
  integer, parameter :: NtypeSigma = 2               ! types of Sigma


  !double precision   :: GI(NtypG, -MxOmegaG1:MxOmegaG1)      ! matrix of G
  !double precision   :: GITailP(NtypG,Nbasis),GITailN(NtypG,Nbasis)   ! positive and negative tails of G         
  !double precision, allocatable   :: WR(:, :, :, :)         ! matrix of W        
  !double precision, allocatable   :: WRTailP(:, :, :, :)    ! positive tail of W
  !double precision, allocatable   :: WRTailC(:, :, :)       ! constant for the tail of W
  !double precision, allocatable   :: W0InMoment(:, :, :)    ! matrix of W in momentum space

  !double precision, allocatable   :: GamR(:, :, :, :, :)    ! matrix of Gamma 
  !double precision, allocatable   :: GamRTailNP(:, :, :, :) ! negative Gin-tail & positive Gout-tail 
  !double precision, allocatable   :: GamRTailPN(:, :, :, :) ! positive Gin-tail & positive Gout-tail 

  !double precision, allocatable   :: GamRTailPPR(:, :, :, :) ! positive Gin-tail & positive Gout-tail 
  !double precision, allocatable   :: GamRTailPPL(:, :, :, :) ! positive Gin-tail & positive Gout-tail 
  !double precision, allocatable   :: GamRTailNNR(:, :, :, :) ! negative Gin-tail & positive Gout-tail 
  !double precision, allocatable   :: GamRTailNNL(:, :, :, :) ! negative Gin-tail & positive Gout-tail 

  !double precision, allocatable   :: GamRTailPM(:, :, :, :, :) ! positive Gin-tail & Gout-matrix
  !double precision, allocatable   :: GamRTailNM(:, :, :, :, :) ! negative Gin-tail & Gout-matrix
  !double precision, allocatable   :: GamRTailMP(:, :, :, :, :) ! positive Gout-tail & Gin-matrix
  !double precision, allocatable   :: GamRTailMN(:, :, :, :, :) ! negative Gout-tail & Gin-matrix
  !double precision, allocatable   :: GamRTailDiagP(:, :, :, :)  ! Gout==Gin tail
  !double precision, allocatable   :: GamRTailDiagN(:, :, :, :)  ! Gout==Gin tail

  !double precision, allocatable   :: GamRTailC(:, :, :)        ! constant for the tails

  !double precision, allocatable   :: ChiR( :, :, :,  :)        ! matrix of chi        
  !double precision, allocatable   :: trChiR( :, :, :)          ! matrix of the trace of chi      
  !double precision, allocatable   :: PiR(:, :, :, :)           ! matrix of Pi
  !double precision, allocatable   :: SigmaI(:, :, :)           ! matrix of sigma

  !=======================================================================

  !====================== MC Simulation ==================================
  double precision, allocatable   :: Gam1MR(:, :, :, :, :)       ! matrix of Gamma 
  double precision :: Gam2MR(2, -MxOmegaDiag:MxOmegaDiag, -MxOmegaDiag:MxOmegaDiag)  ! matrix of Gamma 

  double precision :: GamNorm, GamNormWeight            ! the weight of the normalization diagram
  double precision, allocatable   :: GamMC(:,:,:,:,:,:,:) ! the measurement of Gamma in MC
  double precision, allocatable   :: GamSqMC(:,:,:,:,:,:,:) ! the measurement of Gamma in MC
  double precision :: GamOrder(0:MxOrder)  ! the configuration number of different orders
  double precision :: GamWormOrder(0:MxOrder)  ! the configuration number of different orders in worm section
  !double precision :: Gam0Bubble
  !double precision :: Gam2Topo(2, -MxOmegaDiag:MxOmegaDiag)

  !double precision :: AveWeightRatio(0:MxOrder, 2)
  !double precision :: HistoOmegaW(-MxOmega: MxOmega)

  double precision :: WeightCurrent                  ! the current weight of the configuration
  double precision, allocatable :: CoefOfSymmetry(:,:)
  double precision :: CoefOfWorm
  double precision :: CoefOfWeight(0:MxOrder)  ! the coeffecients for different orders and worm section

  !------------- MC steps -----------------------------------
  integer          :: iupdate           ! the update
  double precision :: imc, ime          ! imc: the MC step; ime: the Measurement step
  integer          :: NSamp             ! # total MC steps
  integer          :: NBlck 
  integer          :: NToss             ! # MC steps for themeralization

  !------------ basic variables for a diagram --------------------------
  integer          :: Order             ! order of the simulating diagram 
  double precision :: Phase             ! phase of the present diagram
  integer          :: NGLn, NWLn, NVertex ! number of glines, wlines, gamma
  integer          :: MeasureGam         ! measuring vertex
  integer          :: NFermiLoop        ! number of Fermi loops
  logical          :: IsWormPresent     ! to represent if worm is present

  !------------ Ira and Masha variables for a diagram -------------------
  integer          :: Ira, Masha        ! vertex number of Ira and Masha
  integer          :: kMasha            ! excess momentum for Masha (-kMasha for Ira)
  integer          :: OmegaMasha        ! excess frequency for Masha (-OmegaMasha for Ira)
  integer          :: SpinMasha         ! excess spin for Masha (-SpinMasha for Ira)
  double precision :: WeightWorm        ! weight for worm sector

  !----------- data structure for lines -------------------------------------
  !!! the properties of Lines !!!!!!!!!!!!!!!!!!!!!!
  integer, dimension(MxNLn) :: OmegaLn              ! frequency of a line
  integer, dimension(MxNLn) :: kLn                  ! momentum of a line
  integer, dimension(MxNLn) :: KindLn               ! kind of a line: 1 Glines; 2 Wlines
  integer, dimension(MxNLn) :: TypeLn               ! type of a line: 1-2 Glines; 1-6 Wlines
  integer, dimension(2,MxNLn) :: NeighLn            ! 1: begin gamma; 2: end gamma 
  double precision, dimension(MxNLn) :: WeightLn    ! weight of glines and wlines
  integer, dimension(MxNLn) :: StatusLn             ! 2:I&M; 0:normal; 1: measure; -1: empty
  integer, dimension(MxNLn) :: NextLn, List4Ln        ! the next empty block;  the
                                                    ! location in OccuGLn and OccuWLn
  !!! the properties of Lines !!!!!!!!!!!!!!!!!!!!!!
  integer                   :: TailLn               ! the tail of the link
  integer, dimension(MxNGLn) :: Ln4GList             ! the occupied lines
  integer, dimension(MxNWLn) :: Ln4WList             ! the occupied lines

  !----------- data structure for vertexes -------------------------------------
  !!! the properties of Gamma(Vertex) !!!!!!!!!!!!!!!!!!!!!!
  integer, dimension(MxNVertex) :: GXVertex, GYVertex       ! gline sites for Gamma
  integer, dimension(MxNVertex) :: WXVertex, WYVertex       ! wline sites for Gamma
  integer, dimension(MxNVertex) :: DirecVertex              ! Gamma is 1: left, 2: right
  integer, dimension(MxNVertex) :: TypeVertex               ! type of Gamma: 1-6
  integer, dimension(MxNVertex) :: TypeVertexIn, TypeVertexOut ! type of Gamma inside spin: 1-2
  integer, dimension(3,MxNVertex) :: NeighVertex            ! 1: in gline; 2: out gline; 3: wline
  double precision, dimension(MxNVertex) :: WeightVertex    ! weight of  gamma
  integer, dimension(MxNVertex):: StatusVertex              ! 2:I&M; 0:normal; 1: measure; -1: empty
  integer, dimension(MxNVertex) :: NextVertex, List4Vertex  ! NextVertex:  the next empty block
                                                    ! List4Vertex: the location of gammas in List 
  !!! the properties of Gamma(Vertex) !!!!!!!!!!!!!!!!!!!!!!
  integer                   :: TailGam                ! the tail of the link
  integer, dimension(MxNVertex):: Vertex4GamList         ! the occupied gammas

  !----------- spin for W and Gamma---------------------------------------------
  integer, dimension(6,6)   :: TypeGam2W
  integer, dimension(2,2,2,2) :: TypeGW2Gam

  !------------ probabilities functions for updates -----------------------------
  integer, parameter :: Nupdate = 15
  double precision :: Pupdate(Nupdate)  ! the probabilities to call different subroutines       
  double precision :: Fupdate(Nupdate)  ! a function of the summation of Pcall

  double precision :: ProbProp(0:MxOrder, Nupdate)
  double precision :: ProbAcc(0:MxOrder, Nupdate)

  !------------ hashtable for irreducibility --------------------------------------
  integer, dimension(-MxK:MxK) ::  Hash4G     ! hashtable for k in G
  integer, dimension(0:MxK) ::     Hash4W     ! hashtable for k in W

  !------------ probabilities functions for variables -----------------------------
  !---omega: P(omega) = exp(-lamomega * omega)
  double precision, dimension(-MxOmega:MxOmega) :: Pomega, Fomega

  !!---x, y: P(x) = exp(-lamx * x);  P(y) = exp(-lamy * y)
  double precision,parameter :: lamx = 1.d0
  double precision, dimension(-MxLx+1:MxLx-1) :: Px, Fx
  double precision,parameter :: lamy = 1.d0
  double precision, dimension(-MxLy+1:MxLy-1) :: Py, Fy

  integer, parameter :: Nobs_b = 67
  logical :: prt
  double precision :: Obs(Nobs_b, MxNblck)
  double precision , dimension(1:Nobs_b) :: Quan, Ave, Dev, Cor
  
  !=======================================================================

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

end module vrbls_mc
