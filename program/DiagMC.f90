INCLUDE "vrbls_mc.f90"
PROGRAM MAIN
  USE vrbls_mc
  implicit none
  integer :: InpMC, it

  !print *, 'Lx, Ly, Ntoss, Nsamp, nw, Nblck, Jcp, beta, MCOrder, Seed, ISub, InpMC, title'
  !read  *,  Lx, Ly, Ntoss, Nsamp, nw, Nblck, Jcp, beta, MCOrder, Seed, ISub, InpMC, title

  !do i = 1, MCOrder
    !read *, CoefOfWeight(i)
  !enddo
  !read *, CoefOfWorm

  !write(title1, '(f4.2)') beta
  !write(title2, '(i2)')  MCOrder
  !write(title3, '(i14)') Seed

  !title2 = trim(adjustl(title1))//'_'//trim(adjustl(title2))
  !title3 = trim(adjustl(title2))//'_'//trim(adjustl(title3))

  !!================= INITIALIZATION =======================================
  Mu(1)  = 1.d0
  Mu(2)  = 1.d0

  Jcp = 1.d0
  Beta = 0.5d0
  Lx = 4
  Ly = 4
  Vol = Lx*Ly
  dLx = Floor(Lx/2.d0)
  dLy = Floor(Ly/2.d0)

  InpMC = 0

  allocate(W(NTypeW, 0:Lx-1, 0:Ly-1, 0:MxT-1))
  allocate(Gam(NTypeGam, 0:Lx-1, 0:Ly-1, 0:MxT-1, 0:MxT-1))

  allocate(W0PF(0:Lx-1, 0:Ly-1, 0:MxT-1))
  allocate(Polar(0:Lx-1, 0:Ly-1, 0:MxT-1))
  allocate(Chi(0:Lx-1, 0:Ly-1, 0:MxT-1))

  call set_time_elapse
  call set_RNG
  call initialize_self_consistent

  call self_consistent

  open(15, file="G_t.dat")
  
  do it = 0, MxT-1
    write(15, *) (real(it)+0.5d0)*Beta/real(MxT), real(G(1, it)), dimag(G(1, it))
  enddo
  close(15)


  !call def_symmetry

  !!=====================================================================

  !if(ISub==1) then
    !call monte_carlo
  !else if(ISub==2) then
    !call self_consistent
  !endif



CONTAINS
INCLUDE "basic_function.f90"
INCLUDE "self_consistent.f90"
!INCLUDE "monte_carlo.f90"
!INCLUDE "check_conf.f90"
!INCLUDE "analytic_integration.f90"
!INCLUDE "read_write_data.f90"
!INCLUDE "statistics.f90"



SUBROUTINE self_consistent
  implicit none
  integer :: iloop
  logical :: flag

  !------- read the G, W, and Gamma  -------------------
  if(InpMC==0) then

    flag = self_consistent_GW(1.d-8)

    !call calculate_Chi
    !call transfer_Chi(-1)

    !call output_Quantities

    !call write_GWGamma

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

  !!------ calculate G, W in momentum domain --------------
  WOld = (10.d0, 0.d0)
  WNow = weight_W(1, 0, 0, 0)
  self_consistent_GW = .true.

  if(InpMC==0) then
    iloop = 0
    write(*, *) "G-W loop:", iloop, real(Wold), real(WNow)

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

!SUBROUTINE monte_carlo
  !implicit none
  !integer :: istep, iblck, mc_version
  !double precision :: WR, GamR

  !call read_GWGamma
  !call calculate_GamNormWeight
  !call def_conf

  !if(InpMC==0) then

    !ProbProp(:,:) = 0.d0
    !ProbAcc(:,:) = 0.d0

    !!-------- throw away some configurations to thermalize -----------
    !do istep = 1, Ntoss
      !call markov
    !enddo
    !!call print_config

    !call time_elapse
    !t_simu = t_elap
    !write(*,52) t_simu
    !52 format(/'thermalization time:',f16.7,2x,'s')

  !else if(InpMC==1) then

    !!------- read the configuration and MC data from previous simulation --
    !call read_monte_carlo_conf
    !call print_config
    !call check_config

  !endif

  !!!================ MC SIMULATION FOR GAMMA =============================
  !imc = 0.d0
  !ime = 0.d0
  !ProbProp(:,:) = 0.d0
  !ProbAcc(:,:) = 0.d0

  !GamOrder(:) = 0.d0
  !GamWormOrder(:) = 0.d0

  !GamMC(:,:,:,:,:,:,:) = 0.d0
  !GamSqMC(:,:,:,:,:,:,:) = 0.d0
  !GamNorm = 0.d0

  !Gam2Topo(:,:) = 0.d0

  !AveWeightRatio(:,:) = 0.d0


  !mc_version = 0


  !do iblck = 1, Nblck
    !do isamp = 1, Nsamp
      !call markov
      !call measure

      !if(Mod(isamp,1000000)==0) then
        !call check_config
        !call print_config

        !call output_GamMC
        !call output_prob_MC

        !call read_flag
        !if(mc_version/=file_version) then
          !call read_GWGamma
          !call update_WeightCurrent
          !mc_version = file_version
        !endif
      !endif

      !if(Mod(isamp,5000000)==0) then
        !call write_monte_carlo_conf
        !call write_monte_carlo_data
      !endif
    !enddo
  !enddo

  !call write_monte_carlo_conf
  !call write_monte_carlo_data

  !call time_elapse
  !t_simu = t_elap
  !write(*,51) t_simu
  !51 format(/'simulation time:',f16.7,2x,'s')

  !return
!END SUBROUTINE monte_carlo

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


END PROGRAM MAIN

