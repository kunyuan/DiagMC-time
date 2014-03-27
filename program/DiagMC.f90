INCLUDE "vrbls_mc.f90"
PROGRAM MAIN
  USE vrbls_mc
  implicit none
  integer :: InpMC, it, i, ISub,ID

  print *, 'Lx, Ly, Ntoss, Nsamp, IsForever, NStep, Jcp, beta, MCOrder, Seed, ISub, InpMC, ID, title'
  read  *,  Lx, Ly, Ntoss, Nsamp, IsForever, NStep, Jcp, beta, MCOrder, Seed, ISub, InpMC, ID, title

  L(1)=Lx
  L(2)=Ly
  logL(:)=dlog(L(:)*1.d0)
  SpatialWeight(:,:)=0.d0

  CoefOfWeight(0) = 1.d0
  do i = 1, MCOrder
    read *, CoefOfWeight(i)
  enddo
  read *, CoefOfWorm

  write(title1, '(f5.2)') beta
  write(title2, '(i2)')  MCOrder
  write(title3, '(i14)') Seed
  write(title4,'(i4)') ID

  title2 = trim(adjustl(title1))//'_'//trim(adjustl(title2))
  title3 = trim(adjustl(title2))//'_'//trim(adjustl(title3))
  title4 = trim(adjustl(title4))//'_'//trim(adjustl(title3))
  title3 = title4

  write(logstr,*) "Initializing..."
  call write_log
  !!================= INITIALIZATION =======================================
  Mu(1)  = 1.d0
  Mu(2)  = 1.d0

  !================== space variables ==================================
  Vol = Lx*Ly
  dLx = Floor(Lx/2.d0)
  dLy = Floor(Ly/2.d0)

  !================ irreducibility check ===============================
  CheckG = .true.
  CheckW = .true.
  CheckGam = .false.
  !CheckGam = .true.

  !================ updates frequency   ================================
  Pupdate( :)  = 0.d0
  Pupdate( 1)  = 1.d0
  Pupdate( 2)  = 1.d0
  Pupdate( 5)  = 1.d0
  Pupdate( 6)  = 1.d0
  Pupdate( 7)  = 1.d0
  Pupdate( 8)  = 1.d0
  Pupdate(11)  = 1.d0
  Pupdate(12)  = 1.d0
  Pupdate(13)  = 1.d0
  Pupdate(14)  = 1.d0
  Pupdate(15)  = 1.d0
  Pupdate(16)  = 1.d0
  Pupdate(17)  = 1.d0
  !Pupdate(18)  = 1.d0

  !===============  Test variables ==================================
  TestData(:)=0.d0
  !===================================================================
  !==============   Statistics ======================================
  Quan(:)=0.d0
  Norm(:)=0.d0
  Error(:)=0.d0
  QuanName(:)="Undefined"

  allocate(W(NTypeW, 0:Lx-1, 0:Ly-1, 0:MxT-1))
  allocate(Gam(NTypeGam, 0:Lx-1, 0:Ly-1, 0:MxT-1, 0:MxT-1))

  allocate(W0PF(0:Lx-1, 0:Ly-1, 0:MxT-1))
  allocate(Gam0PF(0:Lx-1, 0:Ly-1, 0:MxT-1, 0:MxT-1))
  allocate(Polar(0:Lx-1, 0:Ly-1, 0:MxT-1))
  allocate(Chi(0:Lx-1, 0:Ly-1, 0:MxT-1))

  allocate(GamMC(0:MCOrder,1:NTypeGam/2, 0:Lx-1, 0:Ly-1, 0:MxT-1, 0:MxT-1))
  allocate(ReGamSqMC(0:MCOrder,1:NTypeGam/2, 0:Lx-1, 0:Ly-1, 0:MxT-1, 0:MxT-1))
  allocate(ImGamSqMC(0:MCOrder,1:NTypeGam/2, 0:Lx-1, 0:Ly-1, 0:MxT-1, 0:MxT-1))

  MaxStat=1024
  allocate(ObsRecord(1:MaxStat,1:NObs))

  write(logstr,*) "Initializing more..."
  call write_log

  call set_time_elapse
  call set_RNG
  call initialize_self_consistent
  call def_symmetry
  write(logstr,*) "Initializing done!" 
  call write_log

  !!=====================================================================

  if(ISub==1) then
    call self_consistent
  else if(ISub==2) then
    call monte_carlo
  else if(ISub==3) then
    call test_subroutine
  endif

CONTAINS
INCLUDE "basic_function.f90"
INCLUDE "self_consistent.f90"
INCLUDE "monte_carlo.f90"
INCLUDE "check_conf.f90"
INCLUDE "analytic_integration.f90"
INCLUDE "read_write_data.f90"



SUBROUTINE self_consistent
  implicit none
  integer :: iloop
  logical :: flag

  !------- read the G, W, and Gamma  -------------------
  if(InpMC==0) then

    flag = self_consistent_GW(1.d-8)

    call calculate_Chi
    call transfer_Chi_r(-1)
    call transfer_Chi_t(-1)

    call transfer_Sigma_t(-1)

    call output_Quantities

    call write_GWGamma
    !!!======================================================================
  else if(InpMC==1) then

    call read_GWGamma

    call read_monte_carlo_data
    !!-------- update the Gamma matrix with MC data -------
    call Gam_mc2matrix_mc

    !call calculate_Gam1
    !call output_Gam1

    flag = self_consistent_GW(1.d-8)

    call calculate_Chi
    call transfer_Chi_r(-1)
    call transfer_Chi_t(-1)

    call transfer_Sigma_t(-1)
    call output_Quantities

    call write_GWGamma
    call update_flag
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

  !!------ calculate G, W in momentum domain --------------
  WOld = (10.d0, 0.d0)
  WNow = weight_W(1, 0, 0, 0)
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

    WNow = weight_W(1, 0, 0, 0)

    write(logstr, *) "G-W loop:", iloop, real(WOld), real(WNow)
    call write_log
  enddo
  call calculate_Sigma
  call calculate_Polar

  !!-------------------------------------------------------
  call plus_minus_W0(-1)
  call plus_minus_Gam0(-1)

  call transfer_r(-1)
  call transfer_t(-1)
  return
END FUNCTION self_consistent_GW

SUBROUTINE monte_carlo
  implicit none
  integer :: isamp, iblck, mc_version
  double precision :: WR, GamR

  write(logstr,*) "Initializing monte carlo..."
  call write_log

  call read_GWGamma
  call calculate_GamNormWeight  

  call initialize_markov

  write(logstr,*) "Initializing monte carlo done!"
  call write_log

  if(InpMC==0) then

    write(logstr,*) "Start Thermalization..."
    call write_log

    ProbProp(:,:) = 0.d0
    ProbAcc(:,:) = 0.d0

    !-------- throw away some configurations to thermalize -----------
    IsToss=.true.
    call markov(Ntoss)

    write(logstr,*) "Thermalization done!"
    call write_log

    call time_elapse
    t_simu = t_elap
    write(logstr,52) t_simu
    call write_log
    52 format(' Thermalization time:',f16.7,2x,'s')

  else if(InpMC==1) then

    !------- read the configuration and MC data from previous simulation --
    call read_monte_carlo_conf
    call read_monte_carlo_data

    call print_config
    call check_config

  endif

  !!================ MC SIMULATION FOR GAMMA =============================
  imc = 0.d0
  Z_normal=0.0
  Z_worm=0.0
  StatNum=0

  ProbProp(:,:) = 0.d0
  ProbAcc(:,:) = 0.d0

  GamOrder(:) = 0.d0
  GamWormOrder(:) = 0.d0

  GamMC(:,:,:,:,:,:) = (0.d0, 0.d0)
  ReGamSqMC(:,:,:,:,:,:) = 0.d0
  ImGamSqMC(:,:,:,:,:,:) = 0.d0
  GamNorm = (0.d0, 0.d0)
  TestData(:)=0.d0

  do i = 1, MCOrder+1
    QuanName(i)="(total conf)"
    Norm(i) = 1.d0
  enddo

  mc_version = 0

  write(logstr,*) "Simulation started!"
  call write_log

  IsToss=.false.
  call markov(Nsamp)

  call time_elapse
  t_simu = t_elap
  write(logstr,51) t_simu
  call write_log
  51 format(' Simulation time:',f16.7,2x,'s')

  return
END SUBROUTINE monte_carlo

SUBROUTINE read_flag
  implicit none
  open(11, status="old", file="loop.inp")
  read(11, *) file_version
  close(11)
  return
END SUBROUTINE read_flag


SUBROUTINE update_flag
  implicit none

  open(11, status="old", file="loop.inp")
  read(11, *) file_version
  close(11)

  open(12, status="replace", file="loop.inp")
  write(12, *) file_version+1
  close(12)
  return
END SUBROUTINE update_flag

SUBROUTINE test_subroutine
    implicit none
    !integer :: isamp
    !!======== test x,y distribution =========================
    !integer :: i,nr(2),cr(2),dr(2),N,x,y
    !double precision :: hist(2,0:MxLx-1),weight
    !call initialize_markov
    !print *,"Testing..."
    !hist(:,:)=0.d0
    !N=100000
    !weight=0.0
    !cr(:)=0
    !do i=1,N
      !call generate_xy(cr,nr,dr,weight,.true.)
      !hist(1,nr(1))=hist(1,nr(1))+1
      !hist(2,nr(2))=hist(2,nr(2))+1
    !enddo
    !open(11,file="testx.dat")
    !write(11,*) "X:",L(1),logL(1)
    !do x=0,Lx
      !write(11,*) x, hist(1,x)/N, SpatialWeight(1,x)
    !enddo
    !close(11)
    !open(11,file="testy.dat")
    !write(11,*) "Y:",L(2),logL(2)
    !do y=0,Ly
      !write(11,*) y, hist(2,y)/N, SpatialWeight(2,y)
    !enddo
    !close(11)

    !========  test drawing subroutine =====================
    !call initialize_markov
    !call print_config

    !======== analytic_integration =========================
    !call read_GWGamma
    !call calculate_Gam1
    !call output_Gam1

    return
END SUBROUTINE


END PROGRAM MAIN

