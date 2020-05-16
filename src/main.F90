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
!
!  (Updated on 05/16/2020)
!
!------------------------------------------------------------------------
! This simulator can simulate evolutionarily stable strategy (ESS) of LMA
! and reproduce the forest succession patterns. But, since it
! does not include the models of photosynthesis, leaf stomatal
! conductance, transpiration, soil water dynamics, and energy balance, it
! cannot simulate the ESS of allocation as reported in Weng et al. 2015
! Biogeosciences.
!
! Processes included in this simulator are:
!     Phenology
!     Individual-level carbon budget (gain and respiration)
!     Plant growth: Allometry and allocation
!     Reproduction
!     Mortality
!     Population dynamics
!     Soil C-N dynamics
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

! for Vegetation-Remotesensing data assimilation
#define DATA_ASSIMILATION

program ForestESS
   use esdvm
   implicit none
   type(tile_type),  pointer :: vegn
   type(cohort_type),pointer :: cp,cc

   character(len=50),parameter :: namelistfile = 'parameters_CN.nml'
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
   real    :: plantC,plantN, soilC, soilN
   character(len=50) :: plantcohorts,plantCNpools,soilCNpools,allpools       ! output file names
   logical :: new_annual_cycle = .False.
   logical :: switch = .True. ! for invasion by another PFT
   integer :: istat1,istat2,istat3
   integer :: year0, year1, iyears
   integer :: totyears, totdays
   integer :: i, j, k, idays, idoy
   integer :: simu_steps,idata
   character(len=50) filepath_out

   filepath_out='output/'
   ! create output files
   plantcohorts = trim(filepath_out)//'Annual_cohorts.csv'
   plantCNpools = trim(filepath_out)//'Plant_C_N_pools.csv'  ! daily
   soilCNpools  = trim(filepath_out)//'Soil_C_N_pools.csv'   ! daily
   allpools     = trim(filepath_out)//'EcosystemDynamics.csv'

   open(101,file=trim(plantcohorts), ACTION='write', IOSTAT=istat1)
   open(102,file=trim(plantCNpools), ACTION='write', IOSTAT=istat2)
   open(103,file=trim(soilCNpools),  ACTION='write', IOSTAT=istat3)
   open(104,file=trim(allpools),     ACTION='write', IOSTAT=istat3)
   ! head
   if(outputdaily)then
     write(101,'(5(a5,","),25(a8,","))')               &
        'year','doy','cID','PFT',                    &
        'layer','density', 'f_layer', 'LAI',         &
        'NSC','seedC','leafC','rootC','SW-C','HW-C', &
        'NSN','seedN','leafN','rootN','SW-N','HW-N'

     write(102,'(2(a5,","),25(a8,","))')  'year','doy',         &
        'GPP', 'NPP', 'Rh',   &
        'McrbC', 'fineL', 'struL', 'McrbN', 'fineN', 'struN', &
        'mineralN', 'N_uptk'
   endif
   ! Yearly output
   write(103,'(3(a5,","),25(a9,","))')             &
        'cID','PFT','layer','density', 'f_layer',  &
        'dDBH','dbh','height','Acrown',            &
        'wood','nsc', 'NSN','NPPtr','seed',        &
        'NPPL','NPPR','NPPW','GPP-yr','NPP-yr',    &
        'N_uptk','N_fix','maxLAI'
   write(104,'(1(a5,","),80(a8,","))')  'year',              &
        'CAI','LAI','GPP', 'NPP',   'Rh',                    &
        'plantC','soilC',    'plantN', 'soilN',              &
        'NSC', 'SeedC', 'leafC', 'rootC', 'SapwoodC', 'WoodC',    &
        'NSN', 'SeedN', 'leafN', 'rootN', 'SapwoodN', 'WoodN', &
        'McrbC','fineL','struL', 'McrbN', 'fineN',    'struN', &
        'mineralN', 'N_fxed','N_uptk','N_yrMin','N_P2S','N_loss', &
        'reproC','reproN','NewC-C','NewC-N'

   ! Parameter initialization: Initialize PFT parameters
   call initialize_PFT_data(namelistfile)

   ! Initialize vegetation tile and plant cohorts
   allocate(vegn)
   call initialize_vegn_tile(vegn,nCohorts,namelistfile)

   ! Read in forcing data
   call read_forcingdata(forcingData,datalines,days_data,yr_data,timestep)
   !call read_NACPforcing(forcingData,datalines,days_data,yr_data,timestep)
   steps_per_day = int(24.0/timestep)

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
   do i=1,vegn%n_cohorts
      cc => vegn%cohorts(i)
      write(*,*)i,vegn%cohorts(i)%species,vegn%cohorts(i)%height,vegn%cohorts(i)%crownarea

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
             ! ***** not included processes ******
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
        call vegn_growth_EW(vegn, tsoil, soil_theta)

        !!! daily output
        if(outputdaily)then
          do i = 1, vegn%n_cohorts
            cc => vegn%cohorts(i)
            !! cohorts
            write(101,'(5(I5,","),1(F8.1,","),6(F8.3,","),2(F8.2,","),25(F8.2,","))')  &
                iyears,idoy,cc%ccID,cc%species,cc%layer,   &
                cc%nindivs*10000, cc%layerfrac, cc%LAI, &
                cc%NSC, cc%seedC, cc%bl, cc%br, cc%bsw, cc%bHW, &
                cc%NSN*1000, cc%seedN*1000, cc%leafN*1000, &
                cc%rootN*1000,cc%sapwN*1000,cc%woodN*1000
          enddo
          !! Tile level, daily
          write(102,'(2(I5,","),15(F8.4,","))') iyears, idoy,  &
                vegn%GPP, vegn%NPP, vegn%Rh, &
                vegn%MicrobialC, vegn%metabolicL, vegn%structuralL, &
                vegn%MicrobialN*1000, vegn%metabolicN*1000, vegn%structuralN*1000, &
                vegn%mineralN*1000,   vegn%N_uptake*1000
        endif ! daily output
        ! annual calls
        idata = MOD(simu_steps+1, datalines)+1 !
        year1 = forcingData(idata)%year  ! Check if it is the last day of a year
        new_annual_cycle = ((year0 /= year1).OR. &
                           (idata == steps_per_day .and. simu_steps > datalines))

        if(new_annual_cycle)then
            idoy = 0

            write(103,'(2(I6,","),1(F9.2,","))')iyears, vegn%n_cohorts,vegn%annualN*1000
            write(*,  '(2(I6,","),1(F9.2,","))')iyears, vegn%n_cohorts,vegn%annualN*1000
            ! output yearly variables
            write(*,'(3(a5,","),25(a9,","))') &
                'chtID','PFT','layer','density', 'f_layer',  &
                'dDBH','dbh','height','Acrown', &
                'wood','nsc', 'NSN','fNPPseed',     &
                'fNPPL','fNPPR','fNPPW','GPP-yr','NPP-yr', &
                'N_uptk','N_fix','spLAI'

            ! Cohotrs ouput
            do i = 1, vegn%n_cohorts
               cc => vegn%cohorts(i)
               NPPtree = cc%seedC + cc%NPPleaf + cc%NPProot + cc%NPPwood
               fseed = cc%seedC/NPPtree
               fleaf = cc%NPPleaf/NPPtree
               froot = cc%NPProot/NPPtree
               fwood = cc%NPPwood/NPPtree

               dDBH = (cc%DBH   - cc%DBH_ys)*1000.
               write(103,'(1(I7,","),2(I4,","),1(F9.1,","),25(F9.3,","))')        &
                    cc%ccID,cc%species,cc%layer,                        &
                    cc%nindivs*10000, cc%layerfrac,dDBH,                &
                    cc%dbh,cc%height,cc%crownarea,                      &
                    cc%bsw+cc%bHW,cc%nsc,cc%NSN,                        &
                    NPPtree,fseed, fleaf, froot, fwood,                 &
                    cc%annualGPP,cc%annualNPP,                          &
                    cc%N_up_yr*1000,cc%fixedN_yr*1000,                  &
                    spdata(cc%species)%laimax

               ! Screen output
               write(*,'(1(I7,","),2(I4,","),1(F9.1,","),25(F9.3,","))')          &
                    cc%ccID,cc%species,cc%layer,                        &
                    cc%nindivs*10000, cc%layerfrac,dDBH,                &
                    cc%dbh,cc%height,cc%crownarea,                      &
                    cc%bsw+cc%bHW,cc%nsc,cc%NSN,                        &
                    fseed, fleaf, froot, fwood,                         &
                    cc%annualGPP,cc%annualNPP,                          &
                    cc%N_up_yr*1000,cc%fixedN_yr*1000,                  &
                    spdata(cc%species)%laimax
            enddo

            ! ---------- annual call -------------
            ! update the LAImax of each PFT according to available N for next year ! HuangXin: comment this if N cycle is turned off
            if(update_annualLAImax) call vegn_annualLAImax_update(vegn)

            ! Reproduction and mortality
            call vegn_reproduction(vegn)
            call vegn_nat_mortality(vegn, real(seconds_per_year))
            call vegn_starvation(vegn)  ! called daily

            ! Re-organize cohorts
            call relayer_cohorts(vegn)
            call kill_lowdensity_cohorts(vegn)
            call vegn_mergecohorts(vegn)

            ! Summarize tile-based variables for output
            do i = 1, vegn%n_cohorts
               cc => vegn%cohorts(i)
                    ! Vegn C pools:
                    vegn%NSC     = vegn%NSC   + cc%NSC      * cc%nindivs
                    vegn%SeedC   = vegn%SeedC + cc%seedC    * cc%nindivs
                    vegn%leafC   = vegn%leafC + cc%bl       * cc%nindivs
                    vegn%rootC   = vegn%rootC + cc%br       * cc%nindivs
                    vegn%SapwoodC= vegn%SapwoodC + cc%bsw   * cc%nindivs
                    vegn%woodC   = vegn%woodC    + cc%bHW   * cc%nindivs
                    vegn%LAI     = vegn%LAI   + cc%leafarea * cc%nindivs
                    ! Vegn N pools
                    vegn%NSN     = vegn%NSN   + cc%NSN      * cc%nindivs
                    vegn%SeedN   = vegn%SeedN + cc%seedN    * cc%nindivs
                    vegn%leafN   = vegn%leafN + cc%leafN    * cc%nindivs
                    vegn%rootN   = vegn%rootN + cc%rootN    * cc%nindivs
                    vegn%SapwoodN= vegn%SapwoodN + cc%sapwN * cc%nindivs
                    vegn%woodN   = vegn%woodN    + cc%woodN * cc%nindivs
                    vegn%fixedN_yr  = vegn%fixedN_yr  + cc%fixedN_yr * cc%nindivs
            enddo

            plantC = vegn%NSC + vegn%SeedC + vegn%leafC + vegn%rootC + vegn%SapwoodC + vegn%woodC
            plantN = vegn%NSN + vegn%SeedN + vegn%leafN + vegn%rootN + vegn%SapwoodN + vegn%woodN
            soilC  = vegn%MicrobialC + vegn%metabolicL + vegn%structuralL
            soilN  = vegn%MicrobialN + vegn%metabolicN + vegn%structuralN + vegn%mineralN
            write(104,'(1(I5,","),22(F9.4,","),6(F9.3,","),18(F10.4,","))') &
                  iyears,       &
                  vegn%CAI,vegn%LAI, &
                  vegn%annualGPP, vegn%annualNPP, vegn%annualRh, &
                  plantC,soilC,plantN *1000, soilN * 1000, &
                  vegn%NSC, vegn%SeedC, vegn%leafC, vegn%rootC,  &
                  vegn%SapwoodC, vegn%woodC,                     &
                  vegn%NSN*1000, vegn%SeedN*1000, vegn%leafN*1000, vegn%rootN*1000, &
                  vegn%SapwoodN *1000,  vegn%WoodN *1000,  &
                  vegn%MicrobialC, vegn%metabolicL, vegn%structuralL, &
                  vegn%MicrobialN*1000, vegn%metabolicN*1000, vegn%structuralN*1000, &
                  vegn%mineralN*1000,   vegn%fixedN_yr*1000, vegn%accu_Nup*1000, &
                  vegn%annualN*1000,vegn%N_P2S_yr*1000, vegn%Nloss_yr*1000, &
                  vegn%totseedC*1000,vegn%totseedN*1000,vegn%totNewCC*1000,vegn%totNewCN*1000


            ! set annual variables zero
            call vegn_annual_diagnostics_zero(vegn)

            !! invasion test
            if(do_migration)then
              if(iyears>=600 .and. mod(iyears,300)==0 .and. switch)then
                !call vegn_species_switch(vegn,2,iyears,200)
                vegn%cohorts(1)%species = 5 - vegn%cohorts(1)%species ! mod(iyears/200,2) + 2

                !! move foliage carbon and nitrogen to soil
                if(vegn%cohorts(1)%bl > 0.0) then
                    vegn%structuralL = vegn%structuralL + vegn%cohorts(1)%bl * vegn%cohorts(1)%nindivs
                    vegn%structuralN = vegn%structuralN + vegn%cohorts(1)%leafN * vegn%cohorts(1)%nindivs
                    vegn%cohorts(1)%bl = 0.0
                    vegn%cohorts(1)%leafN = 0.0
                endif

                !switch = .false.
            endif
          endif ! do_migration
            ! update the years of model run
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
  !character(len=50)  filepath_in
  !character(len=120)  climfile
  !character(len=50)  parafile       ! parameter file
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


  !filepath_in = '/Users/eweng/Documents/BiomeESS/forcingData/'
  !climfile = 'ORNL_forcing.csv' ! 'US-WCrforcing.txt'
  climfile=trim(filepath_in)//trim(climfile)
  write(*,*)climfile
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
     climateData(i)%P_air     = input_data(8,i)        ! pa
     climateData(i)%soilwater = 0.8    ! soil moisture, vol/vol
  enddo
  forcingData => climateData
  write(*,*)"forcing", datalines,days_data,yr_data
end subroutine read_forcingdata

!=============================================================
! for reading in NACP site synthesis forcing
subroutine read_NACPforcing(forcingData,datalines,days_data,yr_data,timestep)
  type(climate_data_type),pointer,intent(inout) :: forcingData(:)
  integer,intent(inout) :: datalines,days_data,yr_data
  real, intent(inout)   :: timestep
  !------------local var -------------------
  type(climate_data_type), pointer :: climateData(:)
  !character(len=50)  filepath_in
!  character(len=50)  climfile
  character(len=80)  commts
  integer, parameter :: niterms=15       ! NACP site forcing
  integer, parameter :: ilines=22*366*48 ! the maxmum records
  integer,dimension(ilines) :: year_data
  real,   dimension(ilines) :: doy_data,hour_data
  real input_data(niterms,ilines)
  real inputstep
  integer :: istat1,istat2,istat3
  integer :: doy,idays
  integer :: i,j,k
  integer :: m,n

  !filepath_in = 'input/'
  climfile=trim(filepath_in)//trim(climfile)
  write(*,*)'inputfile: ',climfile
! open forcing data
  open(11,file=climfile,status='old',ACTION='read',IOSTAT=istat2)
  write(*,*)istat2
! skip 2 lines of input met data file
  read(11,'(a160)') commts
  read(11,'(a160)') commts
  m       = 0  ! to record the lines in a file
  idays   = 1  ! the total days in a data file
  yr_data = 0 ! to record years of a dataset
  do    ! read forcing files
      m=m+1
      read(11,*,IOSTAT=istat3)year_data(m),doy_data(m),hour_data(m),   &
                              (input_data(n,m),n=1,niterms)

      if(istat3<0)exit
      if(m == 1) then
          doy = doy_data(m)
      else
          doy = doy_data(m-1)
      endif
      if(doy /= doy_data(m)) idays = idays + 1
      !write(*,*)year_data(m),doy_data(m),hour_data(m)
      ! discard one line
      !read(11,*,IOSTAT=istat3)year_data(m),doy_data(m),hour_data(m),   &
      !                        (input_data(n,m),n=1,niterms)
  enddo ! end of reading the forcing file

  timestep = hour_data(2) - hour_data(1)
  write(*,*)"forcing",datalines,yr_data,timestep,dt_fast_yr
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
     climateData(i)%PAR       = input_data(11,i)*2.0  ! umol/m2/s
     climateData(i)%radiation = input_data(11,i)      ! W/m2
     climateData(i)%Tair      = input_data(1,i)       ! air temperature, K
     climateData(i)%Tsoil     = input_data(1,i)       ! soil temperature, K
     climateData(i)%rain      = input_data(7,i)       ! kgH2O m-2 s-1
     climateData(i)%windU     = input_data(5,i)        ! wind velocity (m s-1)
     climateData(i)%P_air     = input_data(9,i)        ! pa
     climateData(i)%RH        = input_data(3,i) !/mol_h2o*mol_air* & ! relative humidity (0.xx)
                                !climateData(i)%P_air/esat(climateData(i)%Tair-273.16)
     climateData(i)%CO2       = input_data(15,i) * 1.0e-6       ! mol/mol
     climateData(i)%soilwater = 0.8    ! soil moisture, vol/vol
  enddo
  forcingData => climateData
  write(*,*)"forcing", datalines,days_data,yr_data

end subroutine read_NACPforcing

!=====================================================
!=====================================================
end program ForestESS
