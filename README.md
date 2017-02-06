Fore Ecological strategy simulator (ForestESS)

These codes are extracted from the version of LM3-PPA used to simulate the forest successional dynamics in the paper of Weng, E., Farrior, CE, Dybzinski, R, and Pacala, SW, et al. 2016 Global Change Biology. Only the subroutines related to plant growth, allometry, reproduction, mortality, and population dynamics are included. It can be used as a forest simulator. 

To run this model, you may just type in “./modelrun.x” if you have a gfortran compiler.
Or, you can type in:

fortran src/esdvm.f90 src/main.F90 -o ess

./ess

fortran src/BasalAreaAnalysis.f90 -o ana

./ana

where, “fortran” is the name of your fortran compiler.

This work was financially supported by US Forest Service and Princeton Environment Institute. The technical details of this model are in: 

Weng, E. S., Farrior, C. E., Dybzinski, R., Pacala, S. W., 2016. Predicting vegetation type through physiological and environmental interactions with leaf traits: evergreen and deciduous forests in an earth system modeling framework. Global Change Biology,  doi: 10.1111/gcb.13542.

Weng, E. S., Malyshev, S., Lichstein, J. W., Farrior, C. E., Dybzinski, R., Zhang, T., Shevliakova, E., Pacala, S. W., 2015. Scaling from individual trees to forests in an Earth system modeling framework using a mathematically tractable model of height-structured competition. Biogeosciences, 12: 2655–2694, doi:10.5194/bg-12-2655-2015.

Please contact Ensheng Weng (wengensheng@gmail.com) for any qeustions. (02/03/2017)

The processes included in this simulator are:

Phenology

Individual-level carbon budget (gain and respiration)

Plant growth: Allometry and allocation

Reproduction

Mortality

Population dynamics
Soil C-N dynamics


----------------------------------------

! Subroutine call structure:

! ! Initilization:

!     ! Read in forcing data

!     call read_forcingdata(forcingData,datalines,days_data,yr_data,timestep)

!     ! Parameter initialization

!     call initialize_PFT_data()

!     ! Initialize vegetation tile and plant cohorts

!     call initialize_vegn_tile(vegn,nCohorts)

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


