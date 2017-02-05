!========================================================================
!==== Fore Ecological strategy simulator (ForestESS) ====================
!============   Main program   ==========================================
!=============== 12-25-2016 =============================================
!========================================================================
!
! This work was financially supported by US Forest Service and Princeton
! Environment Institute. The technical details of this model are in:
!
! Weng, E. S., Farrior, C. E., Dybzinski, R., Pacala, S. W., 2016. 
! Predicting vegetation type through physiological and environmental 
! interactions with leaf traits: evergreen and deciduous forests in an 
! earth system modeling framework. Global Change Biology, 
! doi: 10.1111/gcb.13542.
!
! Weng, E. S., Malyshev, S., Lichstein, J. W., Farrior, C. E., 
! Dybzinski, R., Zhang, T., Shevliakova, E., Pacala, S. W., 2015. 
! Scaling from individual trees to forests in an Earth system modeling 
! framework using a mathematically tractable model of height-structured 
! competition. Biogeosciences, 12: 2655–2694, doi:10.5194/bg-12-2655-2015.
!
!
! Contact Ensheng Weng (wengensheng@gmail.com) for qeustions.
!                      (02/03/2017)
!
!------------------------------------------------------------------------
! This simulator uses the subroutines from esdvm.F90.
! Included in this simulator are:
!     Phenology
!     Individual-level carbon budget (gain and respiration)
!     Plant growth: Allometry and allocation
!     Reproduction
!     Mortality
!     Population dynamics
!     Soil C-N dynamics
!
!
!----------------------------------------
! Subroutine call structure:
! ! Initilization:
!     ! Read in forcing data
!     call read_forcingdata(forcingData,datalines,days_data,yr_data,timestep)
!     ! Parameter initialization
!     call initialize_PFT_data()
!     ! Initialize vegetation tile and plant cohorts
!     call initialize_vegn_tile(vegn,nCohorts)
!
! ! Daily time step:
!        call vegn_phenology(vegn,j)
!        call vegn_C_N_budget(vegn, tsoil, soil_theta)
!        call vegn_starvation(vegn)
!        call vegn_growth_EW(vegn)
!
! ! Yearly time step:
!            call vegn_reproduction(vegn)
!            call vegn_nat_mortality(vegn, real(seconds_per_year))
!            call vegn_starvation(vegn)
!            ! Re-organize cohorts
!            call relayer_cohorts(vegn)
!            call vegn_mergecohorts(vegn)
!            call kill_lowdensity_cohorts(vegn)
!            ! update LAImax for each PFT according to available N
!            call vegn_annualLAImax_update(vegn)
!            ! set annual variables zero
!            call vegn_annual_diagnostics_zero(vegn)
!----- END -----------------------------------------------------------
!

program BiomeESS
   use esdvm_mod
   implicit none
   type(tile_type),  pointer :: vegn
   type(cohort_type),pointer :: cp,cc

   integer,parameter :: rand_seed = 86456
   integer,parameter :: totalyears = 10
   integer,parameter :: nCohorts = 1
   integer :: datalines ! the total lines in forcing data file
   integer :: yr_data   ! Years of the forcing data
   integer :: days_data ! days of the forcing data
   integer :: steps_per_day ! 24 or 48
   real    :: timestep  ! hour, Time step of forcing data, usually hourly (1.0)
   real    :: tsoil, soil_theta
   real    :: NPPtree,fseed, fleaf, froot, fwood ! for output
   real    :: dDBH ! yearly growth of DBH, mm
   character(len=50) :: plantcohorts,plantCNpools,soilCNpools,allpools       ! output file names
   logical :: new_annual_cycle = .False.
   integer :: istat1,istat2,istat3
   integer :: year0, year1, iyears
   integer :: totyears, totdays
   integer :: i, j, k, idays, idoy
   integer :: simu_steps,idata


   ! create output files
   plantcohorts = 'Annual_cohorts.txt'
   plantCNpools = 'Plant_C_N_pools.csv'  ! daily
   soilCNpools  = 'Soil_C_N_pools.csv'
   allpools     = 'AnnualEcosystemDynamics.csv'

   open(101,file=plantcohorts, ACTION='write', IOSTAT=istat1)
   open(102,file=plantCNpools, ACTION='write', IOSTAT=istat2)
   open(103,file=soilCNpools,  ACTION='write', IOSTAT=istat3)
   open(104,file=allpools,     ACTION='write', IOSTAT=istat3)
   ! head
   write(101,'(3(a5,","),25(a9,","))')               &
        'cID','PFT','layer','density', 'f_layer',    &
        'dDBH','dbh','height','Acrown','totwood',    &
        'nsc', 'NSN',                                &
        'NPPtr','f_seed','f_leaf','f_root','f_wood', &
        'GPP-yr','NPP-yr','Ra_yr','N_uptk','maxLAI'


   write(102,'(5(a5,","),25(a8,","))')               &
        'year','doy','cID','PFT',                    &
        'layer','density', 'f_layer', 'LAI',         &
        'NSC','seedC','leafC','rootC','SW-C','HW-C', &
        'NSN','seedN','leafN','rootN','SW-N','HW-N'

   write(103,'(2(a5,","),25(a8,","))')  'year','doy',         &
        'GPP', 'NPP', 'Rh',   &
        'McrbC', 'fineL', 'struL', 'McrbN', 'fineN', 'struN', &
        'mineralN', 'N_uptk'

   write(104,'(1(a5,","),25(a8,","))')  'year',               &
        'GPP', 'NPP', 'Rh',                                   &
        'NSC','SeedC','leafC','rootC','SW-C', 'HW-C','LAI',   &
        'McrbC', 'fineL', 'struL', 'McrbN', 'fineN', 'struN', &
        'mineralN', 'N_uptk'
   ! Read in forcing data
   call read_forcingdata(forcingData,datalines,days_data,yr_data,timestep)
   steps_per_day = int(24.0/timestep)

   ! Parameter initialization: Initialize PFT parameters
   call initialize_PFT_data()

   ! Initialize vegetation tile and plant cohorts
   allocate(vegn)
   call initialize_vegn_tile(vegn,nCohorts)
   ! total years of model run
   totyears = model_run_years
   totdays  = INT(totyears/yr_data+1)*days_data

   ! Sort and relayer cohorts
   call relayer_cohorts(vegn)

   ! ----- model run ----------
   year0 = forcingData(1)%year
   iyears = 1
   idoy   = 0
   simu_steps = 0
   call vegn_annual_diagnostics_zero(vegn)
   ! output initial cohorts
   write(*,*)"Initial cohorts:", vegn%n_cohorts
   write(*,'(3(I6,","))')0, vegn%n_cohorts
   do i=1,vegn%n_cohorts
      cc => vegn%cohorts(i)
      write(*,*)i,vegn%cohorts(i)%species,vegn%cohorts(i)%height,vegn%cohorts(i)%crownarea
      dDBH = 0.0  ! mm in diameter
      write(*,'(3(I5,","),1(F9.1,","),25(F9.3,","))')          &
            cc%ccID,cc%species,cc%layer,                       &
            cc%nindivs*10000, cc%layerfrac, dDBH,              &
            cc%dbh,cc%height,cc%crownarea,                     &
            cc%bsw+cc%bHW,cc%nsc,cc%NSN,                       &
            NPPtree,fseed, fleaf, froot, fwood,                &
            cc%annualGPP,cc%annualNPP, cc%annualResp,          &
            cc%N_up_yr*1000,spdata(cc%species)%laimax
   enddo

   ! Model run starts here !!
   do idays =1, totdays ! 1*days_data ! days for the model run
        idoy = idoy + 1
        ! get daily mean temperature
        vegn%Tc_daily = 0.0
        tsoil         = 0.0
        do i=1,steps_per_day
             idata = MOD(simu_steps, datalines)+1
             year0 = forcingData(idata)%year  ! Current year
             vegn%Tc_daily = vegn%Tc_daily + forcingData(idata)%Tair
             tsoil         = tsoil + forcingData(idata)%tsoil
             simu_steps = simu_steps + 1

             ! fast-step calls
             ! ***** none ******
             ! leaf photosynthesis
             ! transpiration
             ! plant respiration
             ! soil water dynamics
        enddo
        vegn%Tc_daily = vegn%Tc_daily/steps_per_day
        tsoil         = tsoil/steps_per_day
        soil_theta    = vegn%soil_theta

        ! daily calls
        call vegn_phenology(vegn,j)
        call vegn_C_N_budget(vegn, tsoil, soil_theta)
        call vegn_starvation(vegn)
        call vegn_growth_EW(vegn)

        ! daily output
        do i = 1, vegn%n_cohorts
            cc => vegn%cohorts(i)
            ! cohorts
            !write(102,'(5(I5,","),1(F8.1,","),6(F8.3,","),2(F8.2,","),25(F8.2,","))')  &
            !    iyears,idoy,cc%ccID,cc%species,cc%layer,   &
            !    cc%nindivs*10000, cc%layerfrac, cc%LAI, &
            !    cc%NSC, cc%seedC, cc%bl, cc%br, cc%bsw, cc%bHW, &
            !    cc%NSN*1000, cc%seedN*1000, cc%leafN*1000, &
            !    cc%rootN*1000,cc%sapwdN*1000,cc%woodN*1000
        enddo
        ! Tile level, daily
        write(103,'(2(I5,","),15(F8.4,","))') iyears, idoy,  &
                vegn%GPP, vegn%NPP, vegn%Rh, &
                vegn%MicrobialC, vegn%metabolicL, vegn%structuralL, &
                vegn%MicrobialN*1000, vegn%metabolicN*1000, vegn%structuralN*1000, &
                vegn%mineralN*1000,   vegn%N_uptake*1000

        ! annual calls
        idata = MOD(simu_steps+1, datalines)+1 !
        year1 = forcingData(idata)%year  ! Check if it is the last day of a year
        new_annual_cycle = ((year0 /= year1).OR. &
                           (idata == steps_per_day .and. simu_steps > datalines))

        if(new_annual_cycle)then
            idoy = 0

            write(101,'(2(I6,","),1(F9.2,","))')iyears, vegn%n_cohorts,vegn%annualN*1000
            write(*,  '(2(I6,","),1(F9.2,","))')iyears, vegn%n_cohorts,vegn%annualN*1000
            ! output yearly variables
            write(*,'(3(a5,","),25(a9,","))') &
                'chtID','PFT','layer','density', 'f_layer',  &
                'dDBH','dbh','height','Acrown', &
                'wood','nsc', 'NSN','NPPtr',     &
                'NPPL','NPPR','NPPW','GPP-yr','NPP-yr','N_uptk','spLAI'

            do i = 1, vegn%n_cohorts
               cc => vegn%cohorts(i)
               NPPtree = cc%seedC + cc%NPPleaf + cc%NPProot + cc%NPPwood
               fseed = cc%seedC/NPPtree
               fleaf = cc%NPPleaf/NPPtree
               froot = cc%NPProot/NPPtree
               fwood = cc%NPPwood/NPPtree

               dDBH = (cc%DBH   - cc%DBH_ys)*1000.
               write(101,'(3(I5,","),1(F9.1,","),5(F9.3,","),1(F9.2,","),20(F9.3,","))')       &
                    cc%ccID,cc%species,cc%layer,                       &
                    cc%nindivs*10000, cc%layerfrac, dDBH,              &
                    cc%dbh,cc%height,cc%crownarea,                     &
                    cc%bsw+cc%bHW,cc%nsc,cc%NSN,                       &
                    NPPtree,fseed, fleaf, froot, fwood,                &
                    cc%annualGPP,cc%annualNPP, cc%annualResp,          &
                    cc%N_up_yr*1000,spdata(cc%species)%laimax

               ! Screen output
               write(*,'(3(I5,","),1(F9.1,","),5(F9.3,","),1(F9.2,","),20(F9.3,","))')       &
                    cc%ccID,cc%species,cc%layer,                        &
                    cc%nindivs*10000, cc%layerfrac,dDBH,                &
                    cc%dbh,cc%height,cc%crownarea,                      &
                    cc%bsw+cc%bHW,cc%nsc,cc%NSN,                        &
                    fseed, fleaf, froot, fwood,                         &
                    cc%annualGPP,cc%annualNPP,                          &
                    cc%N_up_yr*1000,      &
                    spdata(cc%species)%laimax

                    ! Vegn pools:
                    vegn%maxNSC     = vegn%maxNSC   + cc%NSC * cc%nindivs
                    vegn%maxSeedC   = vegn%maxSeedC + cc%seedC * cc%nindivs
                    vegn%maxleafC   = vegn%maxleafC + cc%bl * cc%nindivs
                    vegn%maxrootC   = vegn%maxrootC + cc%br * cc%nindivs
                    vegn%SapwoodC   = vegn%SapwoodC + cc%bsw * cc%nindivs
                    vegn%WoodC      = vegn%WoodC    + cc%bHW * cc%nindivs
                    vegn%maxLAI     = vegn%maxLAI   + cc%leafarea * cc%nindivs
            enddo
            ! tile pools
            write(104,'(1(I5,","),13(F8.4,","),6(F8.2,","),2(F8.3,","))') &
                  iyears,       &
                  vegn%annualGPP, vegn%annualNPP, vegn%annualRh, &
                  vegn%maxNSC, vegn%maxSeedC, vegn%maxleafC, vegn%maxrootC,  &
                  vegn%SapwoodC, vegn%WoodC, vegn%maxLAI,                    &
                  vegn%MicrobialC, vegn%metabolicL, vegn%structuralL, &
                  vegn%MicrobialN*1000, vegn%metabolicN*1000, vegn%structuralN*1000, &
                  vegn%mineralN*1000,   vegn%accu_Nup*1000
            ! 'year',                 &
            !'NSC','SeedC','leafC', 'rootC', 'SW-C',  'HW-C', 'LAI', &
            !'McrbC', 'fineL', 'struL', 'McrbN', 'fineN', 'struN',   &
            !'mineralN', 'N_uptake'

            ! ---------- annual call -------------
            ! update the LAImax of each PFT according to available N for next year
            call vegn_annualLAImax_update(vegn)

            ! Reproduction and mortality
            call vegn_reproduction(vegn)
            call vegn_nat_mortality(vegn, real(seconds_per_year))
            call vegn_starvation(vegn)

            ! Re-organize cohorts
            call relayer_cohorts(vegn)
            call vegn_mergecohorts(vegn)
            call kill_lowdensity_cohorts(vegn)

            ! set annual variables zero
            call vegn_annual_diagnostics_zero(vegn)
            iyears = iyears + 1
        endif

   enddo

   !deallocate(cc)
   close(101)
   close(102)
   close(103)
   close(104)
   deallocate(vegn%cohorts)
   deallocate(forcingData)

  contains
!========================================================================

! read in forcing data (Users need to write their own data input procedure)
subroutine read_forcingdata(forcingData,datalines,days_data,yr_data,timestep)
  type(climate_data_type),pointer,intent(inout) :: forcingData(:)
  integer,intent(inout) :: datalines,days_data,yr_data
  real, intent(inout)   :: timestep
  !------------local var -------------------
  type(climate_data_type), pointer :: climateData(:)
  character(len=50)  filepath_in
  character(len=50)  climfile
  character(len=50)  parafile       ! parameter file
  character(len=80)  commts
  integer, parameter :: iiterms=9        ! MDK data for Oak Ridge input
  integer, parameter :: ilines=12*366*24 ! the maxmum records of Oak Ridge FACE, 1999~2007
  integer,dimension(ilines) :: year_data
  real,   dimension(ilines) :: doy_data,hour_data
  real input_data(iiterms,ilines)
  real inputstep
  integer :: istat1,istat2,istat3
  integer :: doy,idays
  integer :: i,j,k
  integer :: m,n

  filepath_in = 'input/'
  climfile    = 'Temperate_forcing.txt'
  climfile=trim(filepath_in)//trim(climfile)

! open forcing data
  open(11,file=climfile,status='old',ACTION='read',IOSTAT=istat2)
  write(*,*)istat2
! skip 2 lines of input met data file
  read(11,'(a160)') commts
! read(11,'(a160)') commts ! MDK data only has one line comments
  m       = 0  ! to record the lines in a file
  idays   = 1  ! the total days in a data file
  yr_data = 0 ! to record years of a dataset
  do    ! read forcing files
      m=m+1
      read(11,*,IOSTAT=istat3)year_data(m),doy_data(m),hour_data(m),   &
                              (input_data(n,m),n=1,iiterms)
      if(istat3<0)exit
      if(m == 1) then
          doy = doy_data(m)
      else
          doy = doy_data(m-1)
      endif
      if(doy /= doy_data(m)) idays = idays + 1

  enddo ! end of reading the forcing file

  timestep = hour_data(2) - hour_data(1)
  !write(*,*)"forcing",datalines,yr_data,timestep
  if (timestep==1.0)then
      write(*,*)"the data freqency is hourly"
  elseif(timestep==0.5)then
      write(*,*)"the data freqency is half hourly"
  else
      write(*,*)"Please check time step!"
      stop
  endif
  close(11)    ! close forcing file
! Put the data into forcing 
  datalines = m - 1
  days_data = idays
  yr_data  = year_data(datalines-1) - year_data(1) + 1

  allocate(climateData(datalines))
  do i=1,datalines
     climateData(i)%year      = year_data(i)          ! Year
     climateData(i)%doy       = doy_data(i)           ! day of the year
     climateData(i)%hod       = hour_data(i)          ! hour of the day
     climateData(i)%PAR       = input_data(1,i)       ! umol/m2/s
     climateData(i)%radiation = input_data(2,i)       ! W/m2
     climateData(i)%Tair      = input_data(3,i) + 273.16  ! air temperature, K
     climateData(i)%Tsoil     = input_data(4,i) + 273.16  ! soil temperature, K
     climateData(i)%RH        = input_data(5,i)           ! relative humidity
     climateData(i)%rain      = input_data(6,i)/(timestep * 3600)! ! kgH2O m-2 s-1
     climateData(i)%windU     = input_data(7,i)        ! wind velocity (m s-1)
     climateData(i)%pressure  = input_data(8,i)        ! pa
     climateData(i)%soilwater = 0.8    ! soil moisture, vol/vol
  enddo
  forcingData => climateData
  write(*,*)"forcing", datalines,days_data,yr_data
end subroutine read_forcingdata

!=====================================================
end program BiomeESS



