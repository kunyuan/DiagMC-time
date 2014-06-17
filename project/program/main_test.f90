INCLUDE "mylib/mylib.f90"
INCLUDE "vrbls_mc.f90"
PROGRAM MAIN
  USE string_basic
  USE logging_module
  USE vrbls_mc
  implicit none
  integer :: it, i, ISub,ID, ios
  logical :: IsLoad
  character(len=128) :: infile

  call get_command_argument(1,infile)
  if(len_trim(infile)==0) then
    call LogFile%QuickLog("No input file is specified!",'e')
    stop
  endif
  open(unit=11, file=trim(infile), action='read')
    read(11,*) ID
    read(11,*) L(1:D)
    read(11,*) Jcp
    read(11,*) Beta
    read(11,*) MCOrder
    read(11,*) IsLoad
    read(11,*) ISub
    if(ISub==2) then
      read(11,*) Ntoss
      read(11,*) Nsamp
      read(11,*) Nstep
      read(11,*) Seed
      read(11,*) title
      read(11,*) CoefOfWorm
      CoefOfWeight(0) = 1.d0
      read(11,*) CoefOfWeight(1:MCOrder)
    elseif(ISub==1 .or. ISub==4) then
      read(11,*) title
    endif
  close(11)

  write(*,*) Beta, MCOrder
  write(title_loop, '(f5.2)') beta
  write(title1, '(i2)')  MCOrder
  write(title2,'(i4)') ID

  title1 = trim(adjustl(title_loop))//'_'//trim(adjustl(title1))
  title_mc = trim(adjustl(title2))//'_'//trim(adjustl(title1))

  if(ISub==1) then
    call LogFile%Initial("project.log","loop")
    call LogTerm%Initial('*','loop')
  else if(ISub==2) then
    call LogFile%Initial(trim(title_mc)//".log","job."//trim(adjustl(title2)))
    call LogTerm%Initial('*',"job."//trim(adjustl(title2)))
  else if(ISub==3) then
    call LogFile%Initial("project.log","integral")
    call LogTerm%Initial('*','integral')
  else if(ISub==4) then
    call LogFile%Initial("project.log","output")
    call LogTerm%Initial('*','output')
  else if(ISub==5) then
    call LogFile%Initial('*','test')
    call LogTerm%Initial('*','test')
  endif

  call LogFile%QuickLog("Initializing basic properties...")

  !!================= INITIALIZATION =======================================
  Mu(1)  = 1.d0
  Mu(2)  = 1.d0

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

  !===============  Test variables ==================================
  TestData(:)=0.d0
  !===================================================================
  !==============   Statistics ======================================
  Quan(:)=0.d0
  Norm(:)=0.d0
  Error(:)=0.d0
  QuanName(:)="Undefined"

  allocate(W(NTypeW, 0:Vol-1, 0:MxT-1))
  allocate(Gam(NTypeGam, 0:Vol-1, 0:MxT-1, 0:MxT-1))

  allocate(W0PF(0:Vol-1, 0:MxT-1))
  allocate(Gam0PF(0:Vol-1, 0:MxT-1, 0:MxT-1))
  allocate(Polar(0:Vol-1, 0:MxT-1))
  allocate(Denom(0:Vol-1, 0:MxT-1))
  allocate(Chi(0:Vol-1, 0:MxT-1))

  allocate(GamMC(0:MCOrder, 0:Vol-1, 0:MxT-1))
  allocate(ReGamSqMC(0:MCOrder, 0:Vol-1, 0:MxT-1))
  allocate(ImGamSqMC(0:MCOrder, 0:Vol-1, 0:MxT-1))

  allocate(GamBasis(0:MCOrder,1:NTypeGam/2, 0:Vol-1, 1:NbinGam, 1:NBasisGam))
  allocate(ReGamSqBasis(0:MCOrder,1:NTypeGam/2, 0:Vol-1, 1:NbinGam, 1:NBasisGam))
  allocate(ImGamSqBasis(0:MCOrder,1:NTypeGam/2, 0:Vol-1, 1:NbinGam, 1:NBasisGam))

  MaxStat=1024
  allocate(ObsRecord(1:MaxStat,0:NObs-1))

  call LogFile%QuickLog("Initializing time and RNG...")

  call set_time_elapse
  call set_RNG

  !=========== initialization of basis =======================
  call initialize_polynomials
  call initialize_bins
  call calculate_basis_GWGam

  call initialize_self_consistent
  call def_symmetry

  call LogFile%QuickLog("Initialization Done!...")
  !!=====================================================================

  call self_consistent

  if(ISub==1) then
    call self_consistent
  !else if(ISub==2) then
    !ios=1
    !do while(ios/=0) 
      !open(10,access="append", iostat=ios, file="read_list.dat")
      !write(10, *) trim(adjustl(title_mc))//"_monte_carlo_data.bin.dat"
      !close(10)
    !enddo
    !call monte_carlo
  else if(ISub==3) then
    call numerical_integeration
  !else if(ISub==4) then
    !call just_output
  !else if(ISub==5) then
    !call test_subroutine
  endif

CONTAINS

INCLUDE "self_consistent.f90"
INCLUDE "sc_functions.f90"
INCLUDE "sc_read_write.f90"
INCLUDE "analytic_integration.f90"
INCLUDE "fitting.f90"

!INCLUDE "monte_carlo.f90"
!INCLUDE "mc_functions.f90"
!INCLUDE "mc_read_write.f90"
!INCLUDE "check_conf.f90"


subroutine numerical_integeration
  implicit none
  call LogFile%QuickLog("Reading G,W, and Gamma...")
  call read_GWGamma

  call calculate_Gam1
end subroutine

!SUBROUTINE just_output
  !implicit none
  !logical :: flag
  !call LogFile%QuickLog("Just output something!")
  !call LogFile%QuickLog("Reading G,W, and Gamma...")
  !call read_GWGamma


  !call LogFile%QuickLog("Reading MC data...")
  !call read_monte_carlo_data

  !call LogFile%QuickLog("Reading Done!...")

  !!!-------- update the Gamma matrix with MC data -------
  !call Gam_mc2matrix_mc

  !flag = self_consistent_GW(1.d-8)

  !call calculate_Chi
  !call transfer_Chi_r(-1)
  !call transfer_Chi_t(-1)
  !call transfer_Sigma_t(-1)

  !call output_Quantities
!end SUBROUTINE just_output

SUBROUTINE self_consistent
  implicit none
  integer :: iloop
  logical :: flag

  !call output_Quantities

  !------- read the G, W, and Gamma  -------------------
  if(IsLoad==.false.) then

    flag = self_consistent_GW(1.d-8)

    call calculate_Chi
    call transfer_Chi_r(-1)
    call transfer_Chi_t(-1)
    call transfer_Sigma_t(-1)

    call output_Quantities

    call write_GWGamma
    !!!======================================================================
  !else if(IsLoad) then

    !call LogFile%QuickLog("Reading G,W, and Gamma...")
    !call read_GWGamma


    !call LogFile%QuickLog("Reading MC data...")
    !call read_monte_carlo_data

    !call LogFile%QuickLog("Reading Done!...")

    !!!-------- update the Gamma matrix with MC data -------
    !call Gam_mc2matrix_mc

    !flag = self_consistent_GW(1.d-8)

    !call calculate_Chi
    !call transfer_Chi_r(-1)
    !call transfer_Chi_t(-1)
    !call transfer_Sigma_t(-1)

    !call output_Quantities

    !call update_flag
    !call write_GWGamma
  endif
  return
END SUBROUTINE self_consistent

LOGICAL FUNCTION self_consistent_GW(err)
  implicit none
  double precision, intent(in) :: err
  integer :: iloop
  integer :: px, py
  complex*16 :: WOld, WNow

  call transfer_r(1)
  call transfer_t(1)

  call plus_minus_W0(1)
  call plus_minus_Gam0(1)

  call output_Quantities

  !!------ calculate G, W in momentum domain --------------
  WOld = (10.d0, 0.d0)
  WNow = weight_W(1, 0, 0)
  self_consistent_GW = .true.

  iloop = 0

  call calculate_Polar
  call calculate_W

  do while(abs(real(WNow)-real(WOld))>err) 
    WOld = WNow
    iloop = iloop + 1

    call calculate_Sigma
    call calculate_Polar

    call calculate_G
    call calculate_W

    WNow = weight_W(1, 0, 0)

    call LogFile%QuickLog("G-W loop:"//str(iloop)//str(real(WOld))//str(real(WNOw)))
  enddo
  call calculate_Sigma
  call calculate_Polar
  call calculate_Denom
  call calculate_Chi

  !!-------------------------------------------------------
  call plus_minus_W0(-1)
  call plus_minus_Gam0(-1)

  call transfer_r(-1)
  call transfer_t(-1)
  return
END FUNCTION self_consistent_GW

!SUBROUTINE monte_carlo
  !implicit none
  !integer :: i, mc_version
  !double precision :: WR, GamR

  !call LogFile%QuickLog("Initializing monte carlo...")
  !call read_GWGamma

  !call calculate_GamNormWeight

  !call initialize_markov

  !call LogFile%QuickLog("Initializing monte carlo done!")

  !do i= 0, MCOrder
    !QuanName(i) = "(Order "+str(i)+"Gamma)"
  !enddo

  !if(IsLoad==.false.) then

    !call LogFile%QuickLog("Start Thermalization ...")

    !ProbCall(:,:) = 0.d0
    !ProbProp(:,:) = 0.d0
    !ProbAcc(:,:) = 0.d0
    !BalenceCheck(:,:,:)=0.d0

    !!-------- throw away some configurations to thermalize -----------
    !call markov(.true.)

    !call LogFile%QuickLog("Thermalization done!")

    !call time_elapse
    !t_simu = t_elap
    !call LogFile%QuickLog('Thermalization time: '//trim(str(t_simu,'(f12.2)'))//'s')
    

    !!!================ MC SIMULATION FOR GAMMA =============================
    !imc = 0.d0
    !Z_normal=0.0
    !Z_worm=0.0
    !StatNum=0

    !ProbCall(:,:) = 0.d0
    !ProbProp(:,:) = 0.d0
    !ProbAcc(:,:) = 0.d0

    !GamOrder(:) = 0.d0
    !GamWormOrder(:) = 0.d0

    !GamMC(:,:,:,:) = (0.d0, 0.d0)
    !ReGamSqMC(:,:,:,:) = 0.d0
    !ImGamSqMC(:,:,:,:) = 0.d0

    !GamBasis(:,:,:,:,:,:) = (0.d0, 0.d0)
    !ReGamSqBasis(:,:,:,:,:,:) = 0.d0
    !ImGamSqBasis(:,:,:,:,:,:) = 0.d0

    !GamNorm = (0.d0, 0.d0)
    !TestData(:)=0.d0

    !call read_flag
    !mc_version = file_version

  !else if(IsLoad) then

    !!------- read the configuration and MC data from previous simulation --
    !call LogFile%QuickLog("Reading the previous MC ...")

    !call read_monte_carlo_conf
    !call LogFile%QuickLog("Read the previous MC conf Done!...")
    !call read_monte_carlo_data
    !call LogFile%QuickLog("Read the previous MC data Done!...")

    !call LogFile%QuickLog(str(mc_version)+', '+str(file_version))
    !call LogFile%QuickLog("Updating G, W, and Gamma...")

    !call read_GWGamma
    !call output_Quantities
    !call update_WeightCurrent

    !call check_config
    !call print_config
  !endif

  !call LogFile%QuickLog("Running MC Simulations...")

  !call markov(.false.)

  !call time_elapse
  !t_simu = t_elap
  !call LogFile%QuickLog('Simulation time: '//trim(str(t_simu,'(f12.2)'))//'s')

  !return
!END SUBROUTINE monte_carlo

SUBROUTINE read_flag
  implicit none
  integer :: ios

  open(11, status="old", iostat=ios, file="loop.inp")
  read(11, *) file_version
  close(11)

  if(ios/=0) then
    call LogFile%QuickLog(str(ios)+"Fail to read the loop number, continue to MC!")
    return
  endif

  return
END SUBROUTINE read_flag


SUBROUTINE update_flag
  implicit none
  integer :: ios 

  open(11, status="old", iostat=ios, file="loop.inp")
  read(11, *) file_version
  close(11)

  if(ios/=0) then
    call LogFile%QuickLog(str(ios)+"Fail to read the loop number in loop.inp!")
    return
  endif

  open(12, status="replace", iostat=ios, file="loop.inp")
  write(12, *) file_version+1
  close(12)

  if(ios/=0) then
    call LogFile%QuickLog(str(ios)+"Fail to write the loop number to loop.inp!")
    return
  endif
  return
END SUBROUTINE update_flag

END PROGRAM MAIN
