INCLUDE "vrbls_mc.f90"
PROGRAM MAIN
  USE vrbls_mc
  implicit none
  integer ::  i, j, iloop, sciloop, iblck, istep
  double precision :: scWOldR, scWWR
  integer :: flag
  integer :: ISub, InpMC
  integer :: isamp, imeas

  print *, 'Lx, Ly, Ntoss, Nsamp, nw, Nblck, Jcp, beta, MCOrder, Seed, ISub, InpMC, title'
  read  *,  Lx, Ly, Ntoss, Nsamp, nw, Nblck, Jcp, beta, MCOrder, Seed, ISub, InpMC, title

  do i = 1, MCOrder
    read *, CoefOfWeight(i)
  enddo
  read *, CoefOfWorm

  write(title1, '(f4.2)') beta
  write(title2, '(i2)')  MCOrder
  write(title3, '(i14)') Seed

  title2 = trim(adjustl(title1))//'_'//trim(adjustl(title2))
  title3 = trim(adjustl(title2))//'_'//trim(adjustl(title3))

  !!================= INITIALIZATION =======================================
  Mu(1)  = 1.d0
  Mu(2)  = 1.d0

  Vol = Lx*Ly
  dLx = Floor(Lx/2.d0)
  dLy = Floor(Ly/2.d0)

  CoefOfWeight(0) = 1.d0

  !!--- allocate the matrix for W, Gamma, Pi, Chi, Sigma -----
  allocate(W0InMoment(ntypW, 0:Lx-1, 0:Ly-1))

  allocate(WRTailC(ntypW, 0:Lx-1, 0:Ly-1))
  allocate(WR(ntypW, 0:Lx-1, 0:Ly-1, 0:MxOmegaW1))
  allocate(WRTailP(ntypW, 0:Lx-1, 0:Ly-1, nbasis))


  allocate(GamRTailC(ntypGa, 0:Lx-1, 0:Ly-1))
  allocate(GamR(ntypGa, 0:Lx-1, 0:Ly-1, -MxOmegaGamG1:MxOmegaGamG1, &
    & -MxOmegaGamG1:MxOmegaGamG1))

  allocate(GamRTailMP(ntypGa, 0:Lx-1, 0:Ly-1, -MxOmegaGamG1:MxOmegaGamG1, nbasis))
  allocate(GamRTailMN(ntypGa, 0:Lx-1, 0:Ly-1, -MxOmegaGamG1:MxOmegaGamG1, nbasis))
  allocate(GamRTailPM(ntypGa, 0:Lx-1, 0:Ly-1, nbasis, -MxOmegaGamG1:MxOmegaGamG1))
  allocate(GamRTailNM(ntypGa, 0:Lx-1, 0:Ly-1, nbasis, -MxOmegaGamG1:MxOmegaGamG1))

  allocate(GamRTailDiagP(ntypGa, 0:Lx-1, 0:Ly-1, nbasis))
  allocate(GamRTailDiagN(ntypGa, 0:Lx-1, 0:Ly-1, nbasis))

  allocate(GamRTailNP(ntypGa, 0:Lx-1, 0:Ly-1, nbasisGamma))
  allocate(GamRTailPN(ntypGa, 0:Lx-1, 0:Ly-1, nbasisGamma))
  allocate(GamRTailPPR(ntypGa, 0:Lx-1, 0:Ly-1, nbasisGamma))
  allocate(GamRTailNNR(ntypGa, 0:Lx-1, 0:Ly-1, nbasisGamma))
  allocate(GamRTailPPL(ntypGa, 0:Lx-1, 0:Ly-1, nbasisGamma))
  allocate(GamRTailNNL(ntypGa, 0:Lx-1, 0:Ly-1, nbasisGamma))

  allocate(ChiR(ntypCh, 0:Lx-1, 0:Ly-1, -MxOmegaChi:MxOmegaChi))
  allocate(trChiR(0:Lx-1, 0:Ly-1, -MxOmegaChi:MxOmegaChi))
  allocate(PiR(ntypPi, 0:Lx-1, 0:Ly-1, -MxOmegaChi:MxOmegaChi))
  allocate(SigmaI(0:MCOrder, ntypSigma, -MxOmegaSigma:MxOmegaSigma))

  allocate(GamMC(0:MCOrder, 0:1, ntypGa/2, 0:Lx-1, 0:Ly-1, -MxOmegaDiag:MxOmegaDiag, &
    & -MxOmegaDiag:MxOmegaDiag))

  !----------------------------------------------------------

  call set_time_elapse
  call set_RNG
  call initialize_self_consistent
  call def_symmetry
  !!=====================================================================

  if(ISub==1) then
    call monte_carlo
  else if(ISub==2) then
    call self_consistent
  endif



CONTAINS
INCLUDE "basic_function.f90"
INCLUDE "self_consistent.f90"
INCLUDE "monte_carlo.f90"
INCLUDE "check_conf.f90"
INCLUDE "first_order_integration.f90"
INCLUDE "read_write_data.f90"
INCLUDE "statistics.f90"



SUBROUTINE self_consistent
  implicit none
  double precision :: WOldR, WWR
  integer :: iloop
  logical :: flag

  !------- read the G, W, and Gamma  -------------------
  if(InpMC==1) then
    call read_GWGamma
    call read_monte_carlo_data
    call output_GamMC

    !-------- update the Gamma matrix with MC data -------
    call Gamma_mc2matrix_mc

    if(self_consistent_GW(1.d-7)) then

      !--- calculation of Pi and Chi -------------------------
      call calculate_Pi
      call calculate_Chi
      call calculate_Sigma
      call output_Quantities

      !----- update the G, W -------------------------------
      call write_GWGamma
      call update_flag

    endif
  else if(InpMC==0) then

    flag = self_consistent_GW(1.d-8)

    call calculate_Pi
    call calculate_Chi
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
  endif


  return
END SUBROUTINE self_consistent


SUBROUTINE monte_carlo
  implicit none
  integer :: istep, iblck, mc_version
  double precision :: WR, GamR

  call read_GWGamma
  call calculate_GamNormWeight
  call def_conf

  if(InpMC==0) then

    ProbProp(:,:) = 0.d0
    ProbAcc(:,:) = 0.d0

    !-------- throw away some configurations to thermalize -----------
    do istep = 1, Ntoss
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
    call print_config
    call check_config

  endif

  !!================ MC SIMULATION FOR GAMMA =============================
  imc = 0.d0
  ime = 0.d0
  ProbProp(:,:) = 0.d0
  ProbAcc(:,:) = 0.d0

  GamOrder(:) = 0.d0
  GamWormOrder(:) = 0.d0

  GamMC(:,:,:,:,:,:,:) = 0.d0
  GamSqMC(:,:,:,:,:,:,:) = 0.d0
  GamNorm = 0.d0

  Gam2Topo(:,:) = 0.d0

  AveWeightRatio(:,:) = 0.d0


  mc_version = 0


  do iblck = 1, Nblck
    do isamp = 1, Nsamp
      call markov
      call measure

      if(Mod(isamp,1000000)==0) then
        call check_config
        call print_config

        call output_GamMC
        call output_prob_MC

        call read_flag
        if(mc_version/=file_version) then
          call read_GWGamma
          call update_WeightCurrent
          mc_version = file_version
        endif
      endif

      if(Mod(isamp,5000000)==0) then
        call write_monte_carlo_conf
        call write_monte_carlo_data
      endif
    enddo
  enddo

  call write_monte_carlo_conf
  call write_monte_carlo_data

  call time_elapse
  t_simu = t_elap
  write(*,51) t_simu
  51 format(/'simulation time:',f16.7,2x,'s')

  return
END SUBROUTINE monte_carlo

SUBROUTINE read_flag
  implicit none
  open(11, status="old", file="selfconsist_loop")
  read(11, *) file_version
  close(11)
  return
END SUBROUTINE read_flag


SUBROUTINE update_flag
  implicit none
  open(11, status="old", file="selfconsist_loop")
  read(11, *) file_version
  close(11)

  open(12, status="replace", file="selfconsist_loop")
  write(12, *) file_version+1
  close(12)
  return
END SUBROUTINE update_flag


LOGICAL FUNCTION self_consistent_GW(err)
  implicit none
  double precision, intent(in) :: err
  integer :: iloop
  double precision :: WOldR, WWR

  call transfer_W(1)
  call transfer_Gamma(1)

  !!------ calculate G, W in momentum domain --------------
  WOldR = 10.d0
  WWR = weight_W(0, 0, 1, 1)
  self_consistent_GW = .true.

  if(InpMC==0) then
    iloop = 0
    do while(abs(WWR-WOldR)>err) 
      WOldR = WWR
      iloop = iloop + 1

      call calculate_G
      call calculate_W

      WWR = weight_W(0, 0, 1, 1)
      write(*, *) "G-W loop:", iloop, WOldR, WWR
    enddo
  else
    do iloop = 1, 10 
      WOldR = WWR

      call calculate_G
      call calculate_W
      WWR = weight_W(0, 0, 1, 1)

      write(*, *) "G-W loop:", iloop, WOldR, WWR
    enddo
  endif

  write(*, *) "G-W loop:", iloop, WOldR, WWR
  !!-------------------------------------------------------

  call transfer_W(-1)
  call transfer_Gamma(-1)
  return
END FUNCTION self_consistent_GW

END PROGRAM MAIN

