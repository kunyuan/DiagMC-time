INCLUDE "mylib/mylib.f90"
INCLUDE "vrbls_mc.f90"
PROGRAM MAIN
  USE string_basic
  USE logging_module
  USE vrbls_mc
  implicit none
  integer :: it, i, isub,ID, ios

  call get_command_argument(1,infile)
  if(len_trim(infile)==0) then
    call LogFile%QuickLog("No input file is specified!",'e')
    stop
  endif

  call read_infile(id, isub)

  MCOrder = iniMCOrder
  print *, iniBeta, dBeta, finalBeta, MCOrder, isub

  !============ assign all the titles ============================
  write(title_loop, '(f5.2)') finalBeta  !the title for G,W files
  write(title1, '(i2)')  MCOrder
  write(title2,'(i4)') ID

  ! the title for different MC jobs
  title_mc = trim(adjustl(title2))//'_'//trim(adjustl(title_loop))//'_'//trim(adjustl(title1))  

  if(isub==1) then
    call LogFile%Initial("project.log","loop")
    call LogTerm%Initial('*','loop')
  else if(isub==2) then
    call LogFile%Initial(trim(title_mc)//".log","job."//trim(adjustl(title2)))
    call LogTerm%Initial('*',"job."//trim(adjustl(title2)))
  else if(isub==3) then
    call LogFile%Initial("project.log","integral")
    call LogTerm%Initial('*','integral')
  else if(isub==4) then
    call LogFile%Initial("project.log","output")
    call LogTerm%Initial('*','output')
  else if(isub==5) then
    call LogFile%Initial('*','test')
    call LogTerm%Initial('*','test')
  endif

  call init_basic

  call LogFile%QuickLog("initialize basic properties done!")

  !!=====================================================================
  if(ISub==1) then
    call self_consistent
  else if(ISub==2) then
    call monte_carlo
  else if(ISub==3) then
    call numerical_integeration
  else if(ISub==4) then
    call just_output
  else if(ISub==5) then
    call test_subroutine
  endif

CONTAINS

INCLUDE "analytic_integration.f90"
INCLUDE "fitting.f90"

INCLUDE "self_consistent.f90"
INCLUDE "monte_carlo.f90"

INCLUDE "sc_functions.f90"
INCLUDE "mc_functions.f90"

INCLUDE "model_dependent.f90"
INCLUDE "read_write.f90"

INCLUDE "check_conf.f90"


subroutine numerical_integeration
  implicit none
  logical :: ifchange
  double precision :: mcBeta

  call read_input(.true.)
  call update_T_dependent

  call read_GWGamma

  call calculate_Gam1
  call output_Gam1
  return
end subroutine

SUBROUTINE just_output
  implicit none
  logical :: flag, ifchange
  double precision :: mcBeta

  call LogFile%QuickLog("Just output something!") 
  call read_input(.true.)
  call update_T_dependent

  call LogFile%QuickLog("Reading G, W...")
  if(read_GW()) call LogFile%QuickLog("Read G, W done!")

  call read_Gamma_MC(ifchange, mcBeta)
  !if(abs(mcBeta-Beta)>1.d-5)  then
    !call LogFile%QuickLog("Beta for Gamma is not the same with beta in input file!",'e')
    !stop -1
  !endif
  call LogFile%QuickLog("Reading Gamma done!")

  flag = self_consistent_GW(.false.)

  call output_Quantities
end SUBROUTINE just_output

SUBROUTINE self_consistent
  implicit none
  logical :: flag, ifchange
  double precision :: mcBeta

  !------- read the G, W, and Gamma  -------------------
  if(IsLoad==.false.) then

    Beta = iniBeta
    MCOrder = iniMCOrder
    file_version = 0 

    call update_T_dependent

    flag = self_consistent_GW(.true.)

    call output_Quantities

    call write_GWGamma
    call write_input(.false., .true.)

    !!!======================================================================
  else if(IsLoad) then

    call read_input(.true.)

    call update_T_dependent

    call LogFile%QuickLog("Reading G, W...")
    if(read_GW()) call LogFile%QuickLog("Read G, W done!")

    call read_Gamma_MC(ifchange, mcBeta)
    !if(abs(mcBeta-Beta)>1.d-5)  then
      !call LogFile%QuickLog("Beta for Gamma is not the same with beta in input file!",'e')
      !stop -1
    !endif
    call LogFile%QuickLog("Reading Gamma done!")

    flag = self_consistent_GW(.false.)

    call output_Quantities

    call write_GWGamma
    call write_input(ifchange, .false.)
  endif
  return
END SUBROUTINE self_consistent

LOGICAL FUNCTION self_consistent_GW(isloop)
  implicit none
  logical, intent(in) :: isloop
  integer :: iloop, istag
  integer :: px, py
  complex*16 :: WOld, WNow, denominator

  call transfer_r(1)
  call transfer_t(1)
  call plus_minus_W0(1)
  call plus_minus_Gam0(1)

  !!------ calculate G, W in momentum domain --------------
  istag = get_site_from_cord(D, L(1:D)/2)
  self_consistent_GW = .true.

  if(isloop) then
    WOld = (10.d0, 0.d0)
    WNow = W(1, istag, 0)
    call calculate_Polar
    call calculate_W

    do iloop = 1, 20
      WOld = WNow

      call calculate_Sigma
      call calculate_Polar

      call calculate_G
      call calculate_W

      call calculate_Denom
      call calculate_Chi

      WNow = W(1, istag, 0)
      call LogFile%QuickLog("G-W loop:"//str(iloop)//str(WNow/W0PF(istag, 0)))

      denominator = Denom(istag, 0)
      call LogFile%QuickLog("denominator: "+str(denominator), 'i')

      !!for test
      !call output_denominator

      if(real(denominator)<1.d-14)  then
        self_consistent_GW = .false.
        call LogFile%QuickLog("denominator touches zero! "+str(denominator), 'e')
        stop -1
      endif
    enddo

  else
    call calculate_Sigma
    call calculate_Polar

    call calculate_G
    call calculate_W

    call calculate_Denom
    call calculate_Chi

    denominator = Denom(istag, 0)
    call LogFile%QuickLog("denominator: "+str(denominator), 'i')

    if(real(denominator)<1.d-14)  then
      self_consistent_GW = .false.
      call LogFile%QuickLog("denominator touches zero! "+str(denominator), 'e')
      stop -1
    endif
  endif

  !!-------------------------------------------------------
  call plus_minus_W0(-1)
  call plus_minus_Gam0(-1)

  call transfer_r(-1)
  call transfer_t(-1)

  call transfer_Polar_r(-1)
  call transfer_Polar_t(-1)
  call transfer_Sigma_t(-1)

  call transfer_Chi_r(-1)
  call transfer_Chi_t(-1)
  return
END FUNCTION self_consistent_GW

SUBROUTINE monte_carlo
  implicit none
  integer :: i, mc_version
  logical :: ifchange
  double precision :: WR, GamR, mcBeta

  call write_readlist

  !===== read the starting beta and L and MCOrder ================
  call read_input(.true.)
  mc_version = file_version

  call update_T_dependent !initialize all the properties dependent with Beta
  call LogFile%QuickLog("updating T dependent properties done!")

  call read_GWGamma

  if(IsLoad==.false.) then

    !==============THERMOLIZATION=====================================
    call LogFile%QuickLog("Start Thermalization ...")

    ProbCall(:,:) = 0.d0
    ProbProp(:,:) = 0.d0
    ProbAcc(:,:) = 0.d0
    BalenceCheck(:,:,:)=0.d0

    !-------- throw away some configurations to thermalize -----------
    call initialize_markov
    call update_WeightCurrent
    call check_config

    call markov(.true.)

    call LogFile%QuickLog("Thermalization done!")
    call time_elapse
    t_simu = t_elap
    call LogFile%QuickLog('Thermalization time: '//trim(str(t_simu,'(f12.2)'))//'s')
    !=================================================================

    !!================ MC SIMULATION FOR GAMMA =============================
    imc = 0.d0
    Z_normal=0.d0
    Z_worm=0.d0
    StatNum=0
    GamNorm = (0.d0, 0.d0)

    ProbCall(:,:) = 0.d0
    ProbProp(:,:) = 0.d0
    ProbAcc(:,:) = 0.d0

    GamOrder(:) = 0.d0
    GamWormOrder(:) = 0.d0

    GamMC(:,:,:) = (0.d0, 0.d0)
    ReGamSqMC(:,:,:) = 0.d0
    ImGamSqMC(:,:,:) = 0.d0

    GamMCBasis(:,:,:,:,:) = (0.d0, 0.d0)
    ReGamSqBasis(:,:,:,:,:) = 0.d0
    ImGamSqBasis(:,:,:,:,:) = 0.d0

    TestData(:)=0.d0

  else
    call read_monte_carlo_data(mcBeta)
    !if(abs(mcBeta-Beta)>1.d-5)  then
      !call LogFile%QuickLog("Beta for Gamma is not the same with beta in input file!",'e')
      !stop -1
    !endif

    call initialize_markov
    call update_WeightCurrent
    call check_config
  endif

  LoopTimes = 1.d0
  call LogFile%QuickLog("Running MC Simulations...")
  call markov(.false.)

  call time_elapse
  t_simu = t_elap
  call LogFile%QuickLog('Simulation time: '//trim(str(t_simu,'(f12.2)'))//'s')

  return
END SUBROUTINE monte_carlo

SUBROUTINE init_space
  implicit none
  integer :: i

  !================== space variables ==================================
  Vol = 1.d0

  do i = 1, D
    dVol(i) = Vol

    Vol = Vol *L(i)
    dL(i) = Floor(L(i)/2.d0)
  enddo

  do i = D+1, 3
    dVol(i) = Vol
  enddo

  logL(:)=dlog(L(:)*1.d0)
  SpatialWeight(:,:)=0.d0

  call def_symmetry

  return
END SUBROUTINE init_space

SUBROUTINE init_matrix
  implicit none

  allocate(W(NTypeW, 0:Vol-1, 0:MxT-1))
  allocate(Gam(NTypeGam, 0:Vol-1, 0:MxT-1, 0:MxT-1))
  allocate(GamBasis(NTypeGam, 0:Vol-1, 1:NbinGam, 1:NBasisGam))

  allocate(W0PF(0:Vol-1, 0:MxT-1))
  allocate(Gam0PF(0:Vol-1, 0:MxT-1, 0:MxT-1))
  allocate(Polar(0:Vol-1, 0:MxT-1))
  allocate(Denom(0:Vol-1, 0:MxT-1))
  allocate(Chi(0:Vol-1, 0:MxT-1))

  allocate(GamMC(0:MCOrder, 0:Vol-1, 0:MxT-1))
  allocate(ReGamSqMC(0:MCOrder, 0:Vol-1, 0:MxT-1))
  allocate(ImGamSqMC(0:MCOrder, 0:Vol-1, 0:MxT-1))

  allocate(GamMCBasis(0:MCOrder,1:NTypeGam/2, 0:Vol-1, 1:NbinGam, 1:NBasisGam))
  allocate(ReGamSqBasis(0:MCOrder,1:NTypeGam/2, 0:Vol-1, 1:NbinGam, 1:NBasisGam))
  allocate(ImGamSqBasis(0:MCOrder,1:NTypeGam/2, 0:Vol-1, 1:NbinGam, 1:NBasisGam))
END SUBROUTINE init_matrix


SUBROUTINE init_basic
  implicit none

  call LogFile%QuickLog("Initializing basic properties...")

  !!================= INITIALIZATION =======================================
  Mu(1)  = 1.d0
  Mu(2)  = 1.d0

  !!================ init basis properties ==============================
  call initialize_polynomials
  call initialize_bins
  call calculate_basis_GWGam

  !================ updates frequency   ================================
  Pupdate( :)  = 0.d0
  Pupdate( 1)  = 1.d0      ! create_worm_along_wline
  Pupdate( 2)  = 1.d0      ! delete_worm_along_wline
  Pupdate( 5)  = 1.d0      ! move_worm_along_wline
  Pupdate( 6)  = 1.d0      ! move_worm_along_gline
  Pupdate( 7)  = 1.d0      ! add_interaction
  Pupdate( 8)  = 1.d0      ! remove_interaction
  Pupdate(11)  = 1.d0      ! reconnect
  Pupdate(12)  = 1.d0      ! change_gline_space
  Pupdate(13)  = 1.d0      ! change_wline_space
  Pupdate(14)  = 1.d0      ! change_Gamma_type
  Pupdate(15)  = 1.d0      ! move_measuring_index
  Pupdate(16)  = 1.d0      ! change_Gamma_time
  Pupdate(17)  = 1.d0      ! change_wline_isdelta
  Pupdate(18)  = 1.d0      ! change_gamma_isdelta
  !===================================================================

  !===============  Test variables ==================================
  TestData(:)=0.d0
  !===================================================================

  !==============  Statistics for MC  ===============================
  Quan(:)=0.d0
  Norm(:)=0.d0
  Error(:)=0.d0
  QuanName(:)="Undefined"

  do i= 0, MCOrder
    QuanName(i) = "(Order "+str(i)+"Gamma)"
  enddo

  MaxStat=1024
  allocate(ObsRecord(1:MaxStat,0:NObs-1))

  !--------------- initialize topology variables for MC ---------------
  GLnKey2Value(:) = 0
  WLnKey2Value(:) = 0
  VertexKey2Value(:) = 0

  StatusLn(:) = -1
  StatusVertex(:)= -1

  WeightLn(:) = 1.d0
  WeightVertex(:) = 1.d0

  Hash4G(:) = 0
  Hash4W(:) = 0

  do i = 1, MxNLn
    NextLn(i) = i+1
    if(i==MxNLn)  NextLn(i) = -1
  enddo

  do i = 1, MxNVertex
    NextVertex(i) = i+1
    if(i==MxNVertex)  NextVertex(i) = -1
  enddo

  call def_prob
  call def_spin
  !=================================================================

  call init_space
  call init_matrix

  call def_spatial_weight

  call set_time_elapse
  call set_RNG

  call LogFile%QuickLog("Initialization Done!...")
  return
END SUBROUTINE init_basic

SUBROUTINE update_T_dependent
  implicit none

  call initialize_self_consistent
  call initialize_Gam
  GamBasis = (0.d0, 0.d0)

  call calculate_GamNormWeight !dependent with beta
  return
END SUBROUTINE update_T_dependent

SUBROUTINE test_subroutine
  implicit none
  integer :: ir, i, it1, it2
  integer :: isamp
  logical :: ifchange
  double precision :: mcBeta

  call read_input(.true.)
  call update_T_dependent

  call read_GWGamma

  call output_Quantities

  return
END SUBROUTINE test_subroutine

END PROGRAM MAIN

