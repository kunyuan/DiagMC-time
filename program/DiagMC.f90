INCLUDE "vrbls_mc.f90"
PROGRAM MAIN
  USE vrbls_mc
  implicit none
  integer :: InpMC, it, i, ISub

  print *, 'Lx, Ly, Ntoss, Nsamp, Nblck, NStep, Jcp, beta, MCOrder, Seed, ISub, InpMC, title'
  read  *,  Lx, Ly, Ntoss, Nsamp, Nblck, NStep, Jcp, beta, MCOrder, Seed, ISub, InpMC, title

  logLx=dlog(Lx*1.d0)
  logLy=dlog(Ly*1.d0)
  SpatialWeight(:,:)=0.d0

  do i = 1, MCOrder
    read *, CoefOfWeight(i)
  enddo
  read *, CoefOfWorm

  write(title1, '(f5.2)') beta
  write(title2, '(i2)')  MCOrder
  write(title3, '(i14)') Seed

  title2 = trim(adjustl(title1))//'_'//trim(adjustl(title2))
  title3 = trim(adjustl(title2))//'_'//trim(adjustl(title3))

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

  !================ updates frequency   ================================
  Pupdate( :)  = 0.d0
  Pupdate( 1)  = 1.d0
  Pupdate( 2)  = 1.d0
  Pupdate(13)  = 1.d0
  Pupdate(14)  = 1.d0
  Pupdate(15)  = 1.d0
  Pupdate(18)  = 1.d0
  Pupdate(19)  = 1.d0

  !===============  Test variables ==================================
  TestData(:)=0.d0
  !===================================================================


  allocate(W(NTypeW, 0:Lx-1, 0:Ly-1, 0:MxT-1))
  allocate(Gam(NTypeGam, 0:Lx-1, 0:Ly-1, 0:MxT-1, 0:MxT-1))

  allocate(W0PF(0:Lx-1, 0:Ly-1, 0:MxT-1))
  allocate(Gam0PF(0:Lx-1, 0:Ly-1, 0:MxT-1, 0:MxT-1))
  allocate(Polar(0:Lx-1, 0:Ly-1, 0:MxT-1))
  allocate(Chi(0:Lx-1, 0:Ly-1, 0:MxT-1))

  allocate(GamMC(0:MCOrder,0:1,1:NTypeGam/2, 0:Lx-1, 0:Ly-1, 0:MxT-1, 0:MxT-1))
  allocate(GamSqMC(0:MCOrder,0:1,1:NTypeGam/2, 0:Lx-1, 0:Ly-1, 0:MxT-1, 0:MxT-1))

  call set_time_elapse
  call set_RNG
  call initialize_self_consistent
  call def_symmetry

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
!INCLUDE "analytic_integration.f90"
INCLUDE "read_write_data.f90"
!INCLUDE "statistics.f90"



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

    call output_Quantities

    call write_GWGamma

    !!===== G, W with 1-order Gamma =============================
    !call calculate_Gamma1   

    !call calculate_Pi
    !call calculate_Chi
    !call calculate_Sigma
    !call output_Quantities

    !!===== G, W with 2-order Gamma =============================
    !call calculate_Gamma2   
    !call output_GamMC

    !!!======= self_consistent loop of G, W, Gamma up to 1st order =======
    !WOldR = 10.d0
    !WWR = weight_W(1, 0, 0, 1)
    !iloop = 0

    !do while(abs(WWR-WOldR)>1.d-7) 
      !WOldR = WWR
      !iloop = iloop + 1 
      !call calculate_Gamma1   
      !flag = self_consistent_GW(1.d-7)
      !WWR = weight_W(1, 0, 0, 1)
      !write(*, *) "first order self-consistent loop:", iloop, WOldR, WWR
    !enddo

    !call calculate_Gamma1

    !call calculate_Pi
    !call calculate_Chi
    !call calculate_Sigma
    !call output_Quantities


    !!!======================================================================
  !else if(InpMC==1) then
    !call read_GWGamma
    !call read_monte_carlo_data
    !call output_GamMC

    !!-------- update the Gamma matrix with MC data -------
    !call Gamma_mc2matrix_mc

    !if(self_consistent_GW(1.d-7)) then

      !!--- calculation of Pi and Chi -------------------------
      !call calculate_Pi
      !call calculate_Chi
      !call calculate_Sigma
      !call output_Quantities

      !!----- update the G, W -------------------------------
      !call write_GWGamma
      !call update_flag

    !endif
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

  if(InpMC==0) then
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

      write(*, *) "G-W loop:", iloop, real(WOld), real(WNow)
    enddo
  !else
    !do iloop = 1, 10 
      !WOldR = WWR

      !call calculate_G
      !call calculate_W
      !WWR = weight_W(0, 0, 1, 1)

      !write(*, *) "G-W loop:", iloop, WOldR, WWR
    !enddo
  endif

  !!-------------------------------------------------------
  call plus_minus_W0(-1)
  call plus_minus_Gam0(-1)

  call transfer_r(-1)
  call transfer_t(-1)
  return
END FUNCTION self_consistent_GW

SUBROUTINE transfer_r(Backforth)
  implicit none
  integer,intent(in) :: Backforth

  call transfer_W_r(Backforth)
  call transfer_Gam_r(Backforth)
  return
END SUBROUTINE

SUBROUTINE transfer_t(Backforth)
  implicit none
  integer,intent(in) :: Backforth

  call transfer_G_t(Backforth)
  call transfer_W_t(Backforth)
  call transfer_Gam_t(Backforth)
  return
END SUBROUTINE

SUBROUTINE monte_carlo
  implicit none
  integer :: isamp, iblck, mc_version
  double precision :: WR, GamR

  call read_GWGamma
  call calculate_GamNormWeight   ! need to be updated
  
  call initialize_markov

  if(InpMC==0) then

    ProbProp(:,:) = 0.d0
    ProbAcc(:,:) = 0.d0

    !-------- throw away some configurations to thermalize -----------
    IsToss=1
    do isamp = 1, Ntoss
      call markov
    enddo
    !call print_config

    call time_elapse
    t_simu = t_elap
    write(*,52) t_simu
    52 format(/'thermalization time:',f16.7,2x,'s')

  else if(InpMC==1) then

    !------- read the configuration and MC data from previous simulation --
    call read_monte_carlo_conf
    call read_monte_carlo_data

    call print_config
    call check_config

  endif

  !!================ MC SIMULATION FOR GAMMA =============================
  imc = 0.d0
  imeasure = 0.d0

  ProbProp(:,:) = 0.d0
  ProbAcc(:,:) = 0.d0

  GamOrder(:) = 0.d0
  GamWormOrder(:) = 0.d0

  GamMC(:,:,:,:,:,:,:) = 0.d0
  GamSqMC(:,:,:,:,:,:,:) = 0.d0
  GamNorm = 0.d0

  mc_version = 0
  IsToss=0

  do iblck = 1, Nblck
    write(*,*) "Block",iblck," Started!"
    call markov

    !call output_GamMC
    !call output_prob_MC

    !call read_flag
    !if(mc_version/=file_version) then
      !call read_GWGamma
      !call update_WeightCurrent
      !mc_version = file_version
    !endif

    call print_config
    !call write_monte_carlo_conf
    !call write_monte_carlo_data
    !call write_monte_carlo_test
  enddo

  call time_elapse
  t_simu = t_elap
  write(*,51) t_simu
  51 format(/'simulation time:',f16.7,2x,'s')

  return
END SUBROUTINE monte_carlo

!SUBROUTINE read_flag
  !implicit none
  !open(11, status="old", file="selfconsist_loop")
  !read(11, *) file_version
  !close(11)
  !return
!END SUBROUTINE read_flag


!SUBROUTINE update_flag
  !implicit none
  !open(11, status="old", file="selfconsist_loop")
  !read(11, *) file_version
  !close(11)

  !open(12, status="replace", file="selfconsist_loop")
  !write(12, *) file_version+1
  !close(12)
  !return
!END SUBROUTINE update_flag

SUBROUTINE test_subroutine
    implicit none
    integer :: isamp
    !======== test x,y distribution =========================
    !integer :: i,x,y,N
    !double precision :: histx(0:MxLx-1),xweight
    !double precision :: histy(0:MxLy-1),yweight
    !call initialize_markov
    !histx(:)=0.d0
    !histy(:)=0.d0
    !N=10000000
    !xweight=0.0
    !yweight=0.0
    !do i=1,N
      !call generate_x(0,x,xweight)
      !histx(x)=histx(x)+1
      !call generate_y(0,x,yweight)
      !histy(y)=histy(y)+1
    !enddo
    !open(11,file="testx.dat")
    !write(11,*) "X:",Lx,logLx
    !do x=0,Lx-1
      !write(11,*) x, histx(x)/N, SpatialWeight(1,x)
    !enddo
    !close(11)
    !open(11,file="testy.dat")
    !write(11,*) "Y:",Ly,logLy
    !do y=0,Ly-1
      !write(11,*) y, histy(y)/N, SpatialWeight(2,y)
    !enddo
    !close(11)
    !========  test drawing subroutine =====================
    call initialize_markov
    do isamp = 1, Ntoss
      call create_worm_along_wline
      if(IsWormPresent) call DRAW
    enddo
    !call DRAW
END SUBROUTINE


END PROGRAM MAIN

