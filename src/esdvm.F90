! The subroutines are from LM3PPA, the version used in Weng et al. 2016.
! This simulator can simulate evolutionarily stable strategy (ESS) of LMA
! and reproduce the forest succession patterns shown in Weng et al.,
! 2016 Global Change Biology along the graidient of temperature. But, it
! does not include the models of photosynthesis, leaf stomatal
! conductance, transpiration, soil water dynamics, and energy balance in LM3PPA.
! So, for simulating the ESS of allocation and the ESS of wood density over
! the gradient of precipitation and NPP, we need to include:
!     1. Plant hydraulics
!     2. Soil water dynamics:
!          dWs/dt = Prcp - Transp - Evap - Runoff


module esdvm
 implicit none
! ---- public types -------
 public :: spec_data_type, cohort_type, tile_type

 !===============constants===============
 logical, public, parameter :: read_from_parameter_file = .TRUE.
 integer, public, parameter :: days_per_year  = 365
 integer, public, parameter :: hours_per_year = 365 * 24  ! 8760
 integer, public, parameter :: seconds_per_year = 365 * 24 * 3600
 real,    public, parameter :: dt_fast_yr = 1.0/365.0 ! daily

 real,    public, parameter :: PI = 3.1415926
 integer, public, parameter :: MSPECIES = 15
 integer, public, parameter :: max_lev  = 30 ! Soil layers, for soil water dynamics
 integer, public, parameter :: LEAF_ON  = 1
 integer, public, parameter :: LEAF_OFF = 0
 ! Soil SOM reference C/N ratios
 real, public, parameter :: CN0metabolicL  = 15.0 ! 25.0 ! 15.0
 real, public, parameter :: CN0structuralL = 35.0 ! 55.0 ! 35.0

! ---- public variables ---------
 public :: forcingData,spdata, MaxCohortID, &
    K1, K2, K_nitrogen, etaN, MLmixRatio, &
    fsc_fine, fsc_wood,  &
    GR_factor,  l_fract, &
    DBH_mort, A_mort, B_mort

! ------ public subroutines ---------
public :: initialize_PFT_data, initialize_cohort_from_biomass, &
          initialize_vegn_tile
public :: vegn_phenology,vegn_C_N_budget, vegn_growth_EW
public :: vegn_reproduction, vegn_annualLAImax_update
public :: vegn_starvation, vegn_nat_mortality, vegn_species_switch
public :: relayer_cohorts, vegn_mergecohorts, kill_lowdensity_cohorts
public :: vegn_annual_diagnostics_zero


!===============data types ==============================
!-----------PFT data type----------------
type spec_data_type
  integer :: lifeform     ! 0 for grasses, 1 for trees
  integer :: phenotype    ! phenology type: 0 for deciduous, 1 for evergreen
  integer :: pt           ! photosynthetic physiology of species
  ! leaf traits
  real    :: LMA          ! leaf mass per unit area, kg C/m2
  real    :: leafLS       ! leaf life span
  real    :: alpha_L      ! leaf turn over rate
  real    :: LNA          ! leaf Nitrogen per unit area, kg N/m2
  real    :: LNbase       ! basal leaf Nitrogen per unit area, kg N/m2, (Rubisco)
  real    :: CNleafsupport! leaf structural tissues, 175
  real    :: leaf_size    ! characteristic leaf size
  real    :: alpha_phot   ! photosynthesis efficiency
  real    :: m_cond       ! factor of stomatal conductance
  real    :: Vmax         ! max rubisco rate, mol m-2 s-1
  real    :: Vannual      ! annual productivity per unit area at full fun (kgC m-2 yr-1)
  real    :: gamma_L      ! leaf respiration coeficient
  real    :: gamma_LN     ! leaf respiration coeficient per unit N
  ! root traits
  real    :: rho_FR       ! material density of fine roots (kgC m-3)
  real    :: root_r       ! radius of the fine roots, m
  real    :: SRA          ! speific fine root area, m2/kg C
  real    :: gamma_FR     ! Fine root respiration rate, kgC kgC-1 yr-1
  real    :: alpha_FR     ! Turnover rate of Fine roots, fraction yr-1
  real    :: root_perm    ! fine root membrane permeability per unit area, kg/(m3 s)
!  real    :: rho_N_up0   ! maximum N uptake rate
!  real    :: N_roots0    ! root biomass at half of max. N-uptake rate
  real    :: NfixRate0    ! Reference N fixation rate (kgN kgC-1 root)
  real    :: NfixCost0    ! Carbon cost of N fixation (kgC kgN-1)
  ! wood traits
  real    :: rho_wood     ! woody density, kg C m-3 wood
  real    :: gamma_SW     ! sapwood respiration rate, kgC m-2 Acambium yr-1
  real    :: taperfactor

  ! Allometry
  real    :: alphaHT, thetaHT ! height = alphaHT * DBH ** thetaHT
  real    :: alphaCA, thetaCA ! crown area = alphaCA * DBH ** thetaCA
  real    :: alphaBM, thetaBM ! biomass = alphaBM * DBH ** thetaBM
  real    :: phiRL            ! ratio of fine root to leaf area
  real    :: phiCSA           ! ratio of sapwood CSA to target leaf area
  real    :: tauNSC           ! residence time of C in NSC (to define storage capacity)
  ! Default C/N ratios
  real    :: CNleaf0
  real    :: CNroot0
  real    :: CNsw0
  real    :: CNwood0
  real    :: CNseed0
  ! phenology
  real    :: tc_crit         ! K, for turning OFF a growth season
  real    :: tc_crit_on      ! K, for turning ON a growth season
  real    :: gdd_crit        ! K, critical value of GDD5 for turning ON growth season
  !  vital rates
  real    :: maturalage       ! the age that can reproduce
  real    :: v_seed           ! fracton of G_SF to G_F
  real    :: seedlingsize     ! size of the seedlings, kgC/indiv
  real    :: prob_g,prob_e    ! germination and establishment probabilities
  real    :: mortrate_d_c     ! daily mortality rate in canopy
  real    :: mortrate_d_u     ! daily mortality rate in understory
  ! Population level variables
  real    :: LAImax,underLAImax ! max. LAI
  real    :: LAI_light        ! light controlled maximum LAI
  integer :: n_cc             ! for calculating LAImax via cc%LAImax derived from cc%NSN
  real    :: layerfrac        ! species layer fraction
  real    :: internal_gap_frac ! fraction of internal gaps in the canopy
  ! "internal" gaps are the gaps that are created within the canopy by the
  ! branch fall processes.
end type
! PFT-specific parameters
type(spec_data_type), save :: spdata(0:MSPECIES) ! define PFTs

!----------cohort-----------------
type :: cohort_type
  ! ---- biological prognostic variables
  integer :: ccID    = 0   ! cohort ID
  integer :: species = 0   ! vegetation species
  real    :: gdd     = 0.0   ! for phenology
  integer :: status  = 0   ! growth status of plant: 1 for ON, 0 for OFF
  integer :: layer   = 1   ! the layer of this cohort (numbered from top, top layer=1)
  integer :: firstlayer = 0 ! 0 = never been in the first layer; 1 = at least one year in first layer
  real    :: layerfrac  = 0.0 ! fraction of layer area occupied by this cohort

! for populatin structure
  real    :: nindivs   = 1.0 ! density of vegetation, individuals/m2
  real    :: age       = 0.0 ! age of cohort, years
  real    :: dbh       = 0.0 ! diameter at breast height, m
  real    :: height    = 0.0 ! vegetation height, m
  real    :: crownarea = 1.0 ! crown area, m2/individual
  real    :: leafarea  = 0.0 ! total area of leaves, m2/individual
  real    :: lai       = 0.0 ! leaf area index, m2/m2
! carbon pools
  real    :: bl      = 0.0 ! biomass of leaves, kg C/individual
  real    :: br      = 0.0 ! biomass of fine roots, kg C/individual
  real    :: bsw     = 0.0 ! biomass of sapwood, kg C/individual
  real    :: bHW     = 0.0 ! biomass of heartwood, kg C/individual
  real    :: seedC   = 0.0 ! biomass put aside for future progeny, kg C/individual
  real    :: nsc     = 0.0 ! non-structural carbon, kg C/individual

! ----- carbon fluxes
  real :: gpp  = 0.0 ! gross primary productivity kg C/timestep
  real :: npp  = 0.0 ! net primary productivity kg C/timestep
  real :: resp = 0.0 ! plant respiration
  real :: resl = 0.0 ! leaf respiration
  real :: resr = 0.0 ! root respiration
  real :: resg = 0.0 ! growth respiration
  real :: NPPleaf,NPProot,NPPwood ! to record C allocated to leaf, root, and wood
  real :: annualGPP
  real :: annualNPP
  real :: annualResp

! ---- Nitrogen model related parameters
  real    :: NSNmax = 0.
  real    :: NSN = 0.    ! non-structural N pool
  real    :: leafN = 0.
  real    :: sapwN= 0.
  real    :: woodN = 0. ! N of heart wood
  real    :: rootN = 0. ! N of fine roots
  real    :: seedN = 0. !
  real    :: N_uptake = 0.
  real    :: N_up_yr  = 0.0
  real    :: fixedN ! fixed N at each stem per tree
  real    :: fixedN_yr = 0.0 ! annual N fixation per unit crown area

  ! TODO: see if we can make bl_max, br_max local variables
  real    :: bl_max  = 0.0 ! Max. leaf biomass, kg C/individual
  real    :: br_max  = 0.0 ! Max. fine root biomass, kg C/individual
  real    :: CSAsw   = 0.0
  real    :: topyear = 0.0 ! the years that a plant in top layer
  real    :: DBH_ys             ! DBH at the begining of a year (growing season)

! ---- water uptake-related variables
  real    :: root_length(max_lev) = 0.0 ! individual's root length per unit depth, m of root/m
  real    :: K_r = 0.0 ! root membrane permeability per unit area, kg/(m3 s)
  real    :: r_r = 0.0 ! radius of fine roots, m
  real    :: uptake_frac(max_lev) = 0.0 ! normalized vertical distribution of uptake

! for photosynthesis
!  real :: An_op = 0.0 ! mol C/(m2 of leaf per year)
!  real :: An_cl = 0.0 ! mol C/(m2 of leaf per year)
!  real :: w_scale =-9999
  real :: carbon_gain = 0.0 ! carbon gain since last growth, kg C/individual
  real :: extinct = 0.5     ! light extinction coefficient in the canopy for photosynthesis

end type cohort_type

!---------------------------
type :: tile_type
   integer :: n_cohorts = 0
   integer :: n_years   = 0
   integer :: n_canopycc = 0
   type(cohort_type), pointer :: cohorts(:)=>NULL()
   real :: area  ! m2
   real :: age=0 ! tile age
   ! leaf area index
   real :: LAI  ! leaf area index
   real :: CAI  ! crown area index
   real :: LAIlayer(0:9) = 0.0 ! LAI of each crown layer, max. 9
   ! uptake-related variables
   real :: root_distance(max_lev) ! characteristic half-distance between fine roots, m
   ! averaged quantities for PPA phenology
   real :: tc_daily = 0.0
   real :: gdd      = 0.0 ! growing degree-days
   real :: tc_pheno = 0.0 ! smoothed canopy air temperature for phenology

   ! litter and soil carbon pools
   real :: litter = 0.0 ! litter flux
   real :: MicrobialC  = 0  ! Microbes (kg C/m2)
   real :: metabolicL  = 0  ! fast soil carbon pool, (kg C/m2)
   real :: structuralL = 0  ! slow soil carbon pool, (kg C/m2)

!!  Nitrogen pools, Weng 2014-08-08
   real :: MicrobialN= 0
   real :: metabolicN = 0  ! fast soil nitrogen pool, (kg N/m2)
   real :: structuralN = 0  ! slow soil nitrogen pool, (kg N/m2)
   real :: mineralN = 0  ! Mineral nitrogen pool, (kg N/m2)
   real :: N_input        ! annual N input (kgN m-2 yr-1)
   real :: N_uptake  = 0  ! kg N m-2 yr-1
   real :: accu_Nup       ! accumulated N uptake kgN m-2
   real :: fixedN_yr = 0.  ! fixe N in a tile
   real :: annualN = 0.0  ! annual available N in a year
   real :: Nloss_yr= 0.0  ! annual N loss
   real :: N_P2S_yr= 0.0  ! annual N from plants to soil
   real :: previousN      ! an weighted annual available N
   real :: Wrunoff        ! Water runoff of the veg tile, unit?
   real :: soil_theta     ! volumetric soil water content vol/vol
!  Carbon fluxes
   real :: gpp =0 ! gross primary production, kgC m-2 yr-1
   real :: npp =0 ! net primary productivity
   real :: resp = 0 ! auto-respiration of plants
   real :: nep =0 ! net ecosystem productivity
   real :: rh  =0 ! soil carbon lost to the atmosphere
   ! for annual diagnostics
   real :: annualGPP = 0.0 ! kgC m-2 yr-1
   real :: annualNPP = 0.0
   real :: annualResp = 0.0
   real :: annualRh   = 0.0
   ! for annual reporting at tile level
   real :: NSC, SeedC, leafC, rootC, SapwoodC, WoodC
   real :: NSN, SeedN, leafN, rootN, SapwoodN, WoodN
   real :: totSeedC,totSeedN,totNewCC, totNewCN
end type tile_type

!---------------------------

type :: climate_data_type
   integer :: year          ! Year
   integer :: doy           ! day of the year
   real    :: hod           ! hour of the day
   real    :: PAR           ! check uit
   real    :: radiation     ! W/m2
   real    :: Tair          ! K
   real    :: Tsoil         ! soil temperature, K
   real    :: RH            ! relative humidity
   real    :: rain          ! kgH2O m-2 s-1
   real    :: windU         ! wind velocity (m s-1)
   real    :: P_air         ! pa
   real    :: CO2           ! ppm
   real    :: soilwater     ! soil moisture, vol/vol
end type climate_data_type
!---------------------------
type(climate_data_type),pointer, save :: forcingData(:)

integer :: MaxCohortID = 0
! Constants:
real :: K1 = 2 ! Fast soil C decomposition rate (yr-1)
real :: K2 = 0.05 ! slow soil C decomposition rate (yr-1)
real :: K_nitrogen = 8.0     ! mineral Nitrogen turnover rate
real :: MLmixRatio = 0.8     ! the ratio of C and N returned to litters from microbes
real :: etaN       = 0.025    ! N loss through runoff (organic and mineral)
real :: LMAmin     = 0.02    ! minimum LMA, boundary condition
real :: fsc_fine   = 1.0     ! fraction of fast turnover carbon in fine biomass
real :: fsc_wood   = 0.2     ! fraction of fast turnover carbon in wood biomass
real :: GR_factor  = 0.33 ! growth respiration factor
real :: l_fract    = 0.0 ! 0.25  ! 0.5 ! fraction of the carbon retained after leaf drop

! Ensheng's growth parameters:
real :: wood_fract_min = 0.33333
! for understory mortality rate is calculated as:
! deathrate = mortrate_d_u * ( 1 + A * exp(B*(DBH_mort-DBH))/(1 + exp(B*(DBH_mort-DBH)))
real :: DBH_mort   = 0.025 ! characteristic DBH for mortality
real :: A_mort     = 4.0   ! A coefficient in understory mortality rate correction, 1/year
real :: B_mort     = 30.0  ! B coefficient in understory mortality rate correction, 1/m
! for leaf life span and LMA (leafLS = c_LLS * LMA
real :: c_LLS  = 28.57143 ! yr/ (kg C m-2), 1/LMAs, where LMAs = 0.035

! reduction of bl_max and br_max for the understory vegetation, unitless
real :: understory_lai_factor = 0.25

! -------- PFT-specific parameters ----------
! c4grass  c3grass  temp-decid  tropical  evergreen  BE  BD  BN  NE  ND  G  D  T  A
integer :: pt(0:MSPECIES) = 0
!(/1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0/) ! 0 for C3, 1 for C4
integer :: phenotype(0:MSPECIES)= 0
! (/0,  0,  0,  0,  1,  1,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0 /) ! 0 for Deciduous, 1 for evergreen
integer :: lifeform(0:MSPECIES) = 1 ! life form of PFTs: 0 for grasses, 1 for trees
real :: alpha_FR(0:MSPECIES) = 0.5 ! 1.2 ! Fine root turnover rate yr-1
!(/0.8, 0.8,0.8, 0.8, 0.8,0.8,0.8,0.8,1.0,1.0,0.6, 1.0, 0.55, 0.9, 0.55, 0.55/)

! root parameters
real :: rho_FR(0:MSPECIES) = 200 ! woody density, kgC m-3
real :: root_r(0:MSPECIES) = 2.9E-4
!(/1.1e-4, 1.1e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 1.1e-4, 1.1e-4, 2.2e-4, 2.2e-4/)
real :: root_perm(0:MSPECIES)= 1.0E-5
!(/1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5/)
   ! fine root membrane permeability per unit membrane area, kg/(m3 s).
   ! Root membrane permeability is "high" for the value from Siqueira et al., 2008,
! Water Resource Research Vol. 44, W01432, converted to mass units
!real :: rho_N_up0(0:MSPECIES) = 0.5 ! fraction of mineral N per hour
!real :: N_roots0(0:MSPECIES) = 0.3 ! kgC m-2

real :: leaf_size(0:MSPECIES)= 0.04 !
! photosynthesis parameters
real :: Vmax(0:MSPECIES)= 70.0E-6 !
real :: Vannual(0:MSPECIES) = 1.2 ! kgC m-2 yr-1
!(/1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2/)
real :: m_cond(0:MSPECIES)= 9.0 !
real :: alpha_phot(0:MSPECIES)= 0.06 !
real :: gamma_L(0:MSPECIES)= 0.02 !
real :: gamma_LN(0:MSPECIES)= 25.0 ! kgC kgN-1 yr-1
real :: gamma_SW(0:MSPECIES)= 0.0025 ! kgC m-2 Acambium yr-1
real :: gamma_FR(0:MSPECIES)= 12.0  !kgC kgN-1 yr-1 ! 0.6  ! kgC kgC-1 yr-1
real :: tc_crit(0:MSPECIES)= 283.16
real :: tc_crit_on(0:MSPECIES)= 280.16 !
real :: gdd_crit(0:MSPECIES)= 300.0 !

! Allometry parameters
real :: alphaHT(0:MSPECIES)      = 36.01
real :: thetaHT(0:MSPECIES)      = 0.5 !
real :: alphaCA(0:MSPECIES)      = 200.0
real :: thetaCA(0:MSPECIES)      = 1.5
real :: alphaBM(0:MSPECIES)      = 5200.0
real :: thetaBM(0:MSPECIES)      = 2.5

! Reproduction prarameters
real :: maturalage(0:MSPECIES)   = 5.0  ! year
real :: v_seed(0:MSPECIES)       = 0.1  ! fraction of allocation to wood+seeds
real :: seedlingsize(0:MSPECIES) = 0.05 ! kgC
real :: prob_g(0:MSPECIES)       = 1.0
real :: prob_e(0:MSPECIES)       = 1.0

! Mortality
real :: mortrate_d_c(0:MSPECIES) = 0.01
real :: mortrate_d_u(0:MSPECIES) = 0.05

! Leaf parameters
real :: LMA(0:MSPECIES)          = 0.035  !  leaf mass per unit area, kg C/m2
!(/0.04,    0.04,    0.035,   0.035,   0.140,  0.032, 0.032,  0.036,   0.036,   0.036,   0.036,   0.036,   0.036,   0.036,   0.036,   0.036  /)
real :: leafLS(0:MSPECIES) = 1.0
!(/1., 1., 1., 1., 3., 3., 1., 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 /)
real :: LNbase(0:MSPECIES)        = 0.6E-3 !& !  basal leaf Nitrogen per unit area, kg N/m2
real :: CNleafsupport(0:MSPECIES) = 60.0 ! CN ratio of leaf supporting tissues
real :: rho_wood(0:MSPECIES)      = 265.0 ! kgC m-3
real :: taperfactor(0:MSPECIES)   = 0.65 ! taper factor, from a cylinder to a tree
real :: LAImax(0:MSPECIES)        = 4.0 ! maximum LAI for a tree
real :: LAI_light(0:MSPECIES)     = 4.0 ! maximum LAI limited by light
real :: tauNSC(0:MSPECIES)        = 3.0 ! NSC residence time,years
real :: phiRL(0:MSPECIES)         = 1.2 ! ratio of fine root area to leaf area
real :: phiCSA(0:MSPECIES)        = 1.25E-4 ! ratio of sapwood area to leaf area
! C/N ratios for plant pools
real :: CNleaf0(0:MSPECIES)   = 50. ! C/N ratios for leaves
real :: CNsw0(0:MSPECIES)     = 350.0 ! C/N ratios for woody biomass
real :: CNwood0(0:MSPECIES)   = 350.0 ! C/N ratios for woody biomass
real :: CNroot0(0:MSPECIES)   = 60.0 ! C/N ratios for leaves
real :: CNseed0(0:MSPECIES)   = 20.0 ! C/N ratios for leaves
real :: NfixRate0(0:MSPECIES) = 0.0 !Reference N fixation rate (0.03 kgN kgC-1 root yr-1)
real :: NfixCost0(0:MSPECIES) = 12.0 ! FUN model, Fisher et al. 2010, GBC
real :: internal_gap_frac(0:MSPECIES)= 0.1 ! The gaps between trees

namelist /vegn_parameters_nml/  &
  pt, phenotype, lifeform, &
  Vmax, Vannual,   &
  gamma_L, gamma_LN, gamma_SW, gamma_FR,  &
  rho_FR, root_r, root_perm, &
  !rho_N_up0, N_roots0, &
  leaf_size, leafLS, LAImax, LAI_light,   &
  LMA, LNbase, CNleafsupport, c_LLS,      &
  K1,K2, K_nitrogen, etaN, MLmixRatio,    &
  LMAmin, fsc_fine, fsc_wood, &
  GR_factor, l_fract, wood_fract_min,  &
  gdd_crit,tc_crit, tc_crit_on, &
  alphaHT, thetaHT, alphaCA, thetaCA, alphaBM, thetaBM, &
  maturalage, v_seed, seedlingsize, prob_g,prob_e,      &
  mortrate_d_c, mortrate_d_u,                           &
  DBH_mort, A_mort, B_mort,                       &
  phiRL, phiCSA, rho_wood, taperfactor, &
  tauNSC, understory_lai_factor, &
  CNleaf0,CNsw0,CNwood0,CNroot0,CNseed0, &
  NfixRate0, NfixCost0,  &
  internal_gap_frac

!----- Initial conditions -------------
integer, parameter :: MAX_INIT_COHORTS = 10 ! Weng, 2014-10-01
integer :: init_n_cohorts                        = MAX_INIT_COHORTS
integer :: init_cohort_species(MAX_INIT_COHORTS) = 2
real    :: init_cohort_nindivs(MAX_INIT_COHORTS) = 1.0  ! initial individual density, individual/m2
real    :: init_cohort_bl(MAX_INIT_COHORTS)      = 0.0  ! initial biomass of leaves, kg C/individual
real    :: init_cohort_br(MAX_INIT_COHORTS)      = 0.0  ! initial biomass of fine roots, kg C/individual
real    :: init_cohort_bsw(MAX_INIT_COHORTS)     = 0.05 ! initial biomass of sapwood, kg C/individual
real    :: init_cohort_bHW(MAX_INIT_COHORTS)     = 0.0  ! initial biomass of heartwood, kg C/tree
real    :: init_cohort_seedC(MAX_INIT_COHORTS)   = 0.0  ! initial biomass of seeds, kg C/individual
real    :: init_cohort_nsc(MAX_INIT_COHORTS)     = 0.05 ! initial non-structural biomass, kg C/
!  initial soil Carbon and Nitrogen for a vegn tile, Weng 2012-10-24
real   :: init_fast_soil_C  = 0.0  ! initial fast soil C, kg C/m2
real   :: init_slow_soil_C  = 0.0  ! initial slow soil C, kg C/m2
real   :: init_Nmineral = 0.015  ! Mineral nitrogen pool, (kg N/m2)
real   :: N_input    = 0.0008 ! annual N input to soil N pool, kgN m-2 yr-1

character(len=80) :: filepath_in = '/Users/eweng/Documents/BiomeESS/forcingData/'
character(len=160) :: climfile = 'US-WCrforcing.txt'
integer:: model_run_years = 100
integer   :: equi_days       = 0 ! 100 * 365
logical   :: outputhourly = .False.
logical   :: outputdaily  = .False.
logical   :: do_U_shaped_mortality = .False.
logical   :: update_annualLAImax = .False.
logical   :: do_migration = .False.
logical   :: do_fire = .False.
logical   :: do_closedN_run = .True. !.False.

namelist /initial_state_nml/ &
    init_n_cohorts, init_cohort_species, init_cohort_nindivs, &
    init_cohort_bl, init_cohort_br, init_cohort_bsw, &
    init_cohort_bHW, init_cohort_seedC, init_cohort_nsc, &
    init_fast_soil_C, init_slow_soil_C,    &
    init_Nmineral, N_input, &
    filepath_in,climfile, model_run_years, &
    outputhourly, outputdaily, equi_days, &
    do_U_shaped_mortality,update_annualLAImax, &
    do_fire, do_migration, &
    do_closedN_run
!---------------------------------

 contains
! =============== ESS subroutines ============================================

! Net carbon gain
subroutine vegn_C_N_budget(vegn, tsoil, theta)
! hourly ! Changed to daily, Weng 2016-11-25
! include Nitrogen uptake and carbon budget
! carbon_gain is calculated here to drive plant growth and reproduciton
  type(tile_type), intent(inout) :: vegn
  real, intent(in) :: tsoil ! average temperature of soil, deg K
  real, intent(in) :: theta ! average soil wetness, unitless

  !-------local var
  type(cohort_type), pointer :: cc    ! current cohort
  real :: C_input  ! carbon assimilated per tree per fast time step
  !real :: dBL, dBR, dBSW ! leaf and fine root carbon tendencies
  !real :: turnoverC  ! temporary var for several calculations
  integer :: i

  real :: NSC_supply,LR_demand,LR_deficit
  real :: LeafGrowthMin, RootGrowthMin,NSCtarget,v
  real :: LR_growth,WS_growth
  real :: R_days,fNSC,fLFR,fStem

  ! Carbon gain trhough photosynthesis
  call vegn_C_gain(vegn,forcingData)

  ! Nitrogen uptake
  call vegn_N_uptake(vegn, tsoil, theta)

  ! update soil carbon
  call SOMdecomposition(vegn, tsoil, theta)


end subroutine vegn_C_N_budget

! ============================================================================
subroutine vegn_growth_EW(vegn, tsoil, theta)
! updates cohort biomass pools, LAI, and height using accumulated
! carbon_gain and bHW_gain
  type(tile_type), intent(inout) :: vegn
  real, intent(in) :: tsoil ! average temperature of soil, deg K
  real, intent(in) :: theta ! average soil wetness, unitless
  
  ! ---- local vars
  type(cohort_type), pointer :: cc    ! current cohort
  real :: CSAtot ! total cross section area, m2
  real :: CSAsw  ! Sapwood cross sectional area, m2
  real :: CSAwd  ! Heartwood cross sectional area, m2
  real :: DBHwd  ! diameter of heartwood at breast height, m
  real :: BSWmax ! max sapwood biomass, kg C/individual
  real :: dB_LRS, G_LFR  ! amount of carbon spent on leaf and root growth
  real :: dSeed ! allocation to seeds, Weng, 2016-11-26
  real :: dBL, dBR ! tendencies of leaf and root biomass, kgC/individual
  real :: dBSW ! tendency of sapwood biomass, kgC/individual
  real :: dBHW ! tendency of wood biomass, kgC/individual
  real :: dDBH ! tendency of breast height diameter, m
  real :: dCA ! tendency of crown area, m2/individual
  real :: dHeight ! tendency of vegetation height
  real :: dNsw    ! Nitrogen from SW to HW
  real :: sw2nsc = 0.0 ! conversion of sapwood to non-structural carbon
  real :: b,BL_u,BL_c
  real :: alphaBL, alphaBR
  real :: DBHtp
  real :: N_supply, N_demand,fNr,Nsupplyratio,extrasapwN
  integer :: i

  ! Available carbon and nitrogen for growth, and calculate fluxes
  vegn%gpp = 0.
  vegn%npp = 0.
  vegn%Resp = 0.
  ! Respiration and allocation for growth
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     associate ( sp => spdata(cc%species) )
     ! Maintenance respiration
     call plant_respiration(cc,tsoil) ! get resp per tree per time step

     ! Fetch carbon for growth
     call carbon_for_growth(cc)  ! put carbon into carbon_gain for growth

     cc%resp = cc%resp + cc%resg ! put growth respiration into total resp.
     cc%npp  = cc%gpp - cc%resp ! kgC tree-1 time step-1

!    Weng 2015-09-18
     cc%annualGPP  = cc%annualGPP  + cc%gpp/cc%crownarea ! * dt_fast_yr
     cc%annualNPP  = cc%annualNPP  + cc%npp/cc%crownarea ! * dt_fast_yr
     cc%annualResp = cc%annualResp + cc%resp/cc%crownarea ! * dt_fast_yr
     cc%fixedN_yr  = cc%fixedN_yr  + cc%fixedN/cc%crownarea
     ! accumulate tile-level GPP and NPP
     vegn%gpp = vegn%gpp + cc%gpp * cc%nindivs /dt_fast_yr ! kgC m-2 yr-1
     vegn%npp = vegn%npp + cc%npp * cc%nindivs /dt_fast_yr ! kgC m-2 yr-1
     vegn%resp= vegn%resp+ cc%resp* cc%nindivs /dt_fast_yr ! kgC m-2 yr-1

     END ASSOCIATE
  enddo
  cc => null()
  ! NEP is equal to NNP minus soil respiration
  vegn%nep = vegn%npp - vegn%rh ! kgC m-2 yr-1, though time step is daily
  ! Annual summary:
  vegn%annualGPP  = vegn%annualGPP  + vegn%gpp * dt_fast_yr
  vegn%annualNPP  = vegn%annualNPP  + vegn%npp * dt_fast_yr
  vegn%annualResp = vegn%annualResp + vegn%resp * dt_fast_yr
  vegn%annualRh   = vegn%annualRh   + vegn%rh   * dt_fast_yr ! annual Rh

  ! Growth/Allocation
  DBHtp = 0.8
  fNr   = 0.25
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
 !    call biomass_allocation(cc)
     associate (sp => spdata(cc%species)) ! F2003
     if (cc%status == LEAF_ON) then
        ! calculate the carbon spent on growth of leaves and roots
        G_LFR = max(0.0, min(cc%bl_max+cc%br_max-cc%bl-cc%br,  &
                            (1.-Wood_fract_min)*cc%carbon_gain))
        ! and distribute it between roots and leaves
        dBL  = min(G_LFR, max(0.0, &
          (G_LFR*cc%bl_max + cc%bl_max*cc%br - cc%br_max*cc%bl)/(cc%bl_max + cc%br_max) &
          ))
        dBR  = G_LFR - dBL
        ! calculate carbon spent on growth of sapwood growth
        if(cc%layer == 1 .AND. cc%age > sp%maturalage)then
            dSeed=      sp%v_seed * (cc%carbon_gain - G_LFR)
            dBSW = (1.0-sp%v_seed)* (cc%carbon_gain - G_LFR)
        else
            dSeed= 0.0
            dBSW = cc%carbon_gain - G_LFR
        endif
!       Specially for grasses, temporary
        if(sp%lifeform ==0 )then
            dSeed = dSeed + 0.15*G_LFR
            G_LFR = 0.85 * G_LFR
            dBR   = 0.85 * dBR
            dBL   = 0.85 * dBL
        endif
!!       Nitrogen effects on allocations between wood and leaves+roots
!
!!       Nitrogen demand by leaves, roots, and seeds (Their C/N ratios are fixed.)
        N_demand = dBL/sp%CNleaf0 + dBR/sp%CNroot0 + dSeed/sp%CNseed0
!!       Nitrogen available for all tisues, including wood
        N_supply= MAX(0.0, fNr*cc%NSN)
!!       same ratio reduction for leaf, root, and seed if(N_supply < N_demand)
        IF(N_demand > N_supply )then
            dB_LRS = dBL+dBR+dSeed
            Nsupplyratio = N_supply / N_demand
            dBR  =  Nsupplyratio * dBR
            dBL  =  Nsupplyratio * dBL
            dSeed=  Nsupplyratio * dSeed
            dBSW =  dBSW + (1.0 - Nsupplyratio) * dB_LRS ! dB_LRS: = dBL+dBR+dSeed, original guess)
            N_demand = N_supply ! dBL/sp%CNleaf0 + dBR/sp%CNroot0 + dSeed/sp%CNseed0
        ENDIF
!       update biomass pools
        cc%bl     = cc%bl    + dBL
        cc%br     = cc%br    + dBR
        cc%bsw    = cc%bsw   + dBSW
        cc%seedC  = cc%seedC + dSeed
!!      update nitrogen pools, Nitrogen allocation
        cc%NSN   = cc%NSN   - N_supply
        cc%leafN = cc%leafN + dBL   /sp%CNleaf0
        cc%rootN = cc%rootN + dBR   /sp%CNroot0
        cc%seedN = cc%seedN + dSeed /sp%CNseed0
        cc%sapwN = cc%sapwN + MAX((N_supply - N_demand),0.0)
!       Return excessiive Nitrogen in SW back to NSN
        if(cc%sapwN > cc%bsw/sp%CNsw0)then
           extrasapwN = cc%sapwN - cc%bsw/sp%CNsw0
           cc%NSN     = cc%NSN   + extrasapwN ! MAX(0.0, cc%sapwN - cc%bsw/sp%CNsw0)
           cc%sapwN   = cc%sapwN - extrasapwN ! MAX(0.0, cc%sapwN - cc%bsw/sp%CNsw0)
        endif

!       accumulated C allocated to leaf, root, and wood
        cc%NPPleaf = cc%NPPleaf + dBL
        cc%NPProot = cc%NPProot + dBR
        cc%NPPwood = cc%NPPwood + dBSW

!       update breast height diameter given increase of bsw
        dDBH   = dBSW / (sp%thetaBM * sp%alphaBM * cc%DBH**(sp%thetaBM-1))
        dHeight= sp%thetaHT * sp%alphaHT * cc%DBH**(sp%thetaHT-1) * dDBH
        dCA    = sp%thetaCA * sp%alphaCA * cc%DBH**(sp%thetaCA-1) * dDBH
!       update plant architecture
        cc%DBH       = cc%DBH       + dDBH
        cc%height    = cc%height    + dHeight
        cc%crownarea = cc%crownarea + dCA
        cc%leafarea  = leaf_area_from_biomass(cc%bl,cc%species,cc%layer,cc%firstlayer)
        cc%lai       = cc%leafarea/(cc%crownarea *(1.0-sp%internal_gap_frac))
!       conversion of sapwood to heartwood ! Nitrogen from sapwood to heart wood
        if(sp%lifeform>0)then
           CSAsw  = cc%bl_max/sp%LMA * sp%phiCSA * cc%height ! with Plant hydraulics, Weng, 2016-11-30
           CSAtot = 0.25 * PI * cc%DBH**2
           CSAwd  = max(0.0, CSAtot - CSAsw)
           DBHwd  = 2*sqrt(CSAwd/PI)
           BSWmax = sp%alphaBM * (cc%DBH**sp%thetaBM - DBHwd**sp%thetaBM)
           dBHW   = max(cc%bsw - BSWmax, 0.0)
           dNsw   = dBHW/cc%bsw *cc%sapwN
           ! update C and N of sapwood and wood
           cc%bHW   = cc%bHW   + dBHW
           cc%bsw   = cc%bsw   - dBHW
           cc%sapwN = cc%sapwN - dNsw
           cc%woodN = cc%woodN + dNsw
        endif

!       update bl_max and br_max daily
        BL_c = sp%LMA * sp%LAImax * cc%crownarea * &
               (1.0-sp%internal_gap_frac)
        BL_u = sp%LMA*cc%crownarea*(1.0-sp%internal_gap_frac)* &
                    sp%underLAImax
        if (cc%layer == 1) cc%topyear = cc%topyear + 1.0 /365.0
        if (cc%layer > 1 .and. cc%firstlayer == 0) then ! changed back, Weng 2014-01-23
            cc%bl_max = BL_u
!           Keep understory tree's root low and constant
            cc%br_max = 0.8*cc%bl_max/(sp%LMA*sp%SRA) ! sp%phiRL
        else
            cc%bl_max = BL_u + min(cc%topyear/5.0,1.0)*(BL_c - BL_u)
            cc%br_max = sp%phiRL*cc%bl_max/(sp%LMA*sp%SRA)
        endif
     else  ! if(cc%status == LEAF_OFf)then
        cc%nsc = cc%nsc + cc%carbon_gain
     endif ! "cc%status == LEAF_ON"
     ! reset carbon acculmulation terms
     cc%carbon_gain = 0
    end associate ! F2003
  enddo
  cc => null()

  ! Leaf and fine root turnover
  call vegn_leaf_fine_root_turnover(vegn, tsoil, theta)
  ! update tile and cohort ages
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     cc%age = cc%age + dt_fast_yr
  enddo
  vegn%age = vegn%age + dt_fast_yr

end subroutine vegn_growth_EW ! daily

!============================================================================
subroutine vegn_phenology(vegn,doy) ! daily step
  type(tile_type), intent(inout) :: vegn
  integer, intent(in) :: doy

  ! ---- local vars
  type(cohort_type), pointer :: cc
  integer :: i
  real    :: loss_coarse, loss_fine, lossN_coarse, lossN_fine
  real    :: stem_fall, stem_litter, grassdensity   ! for grasses only
  real    :: dAleaf, dBL, dBR, dNL, dNR       ! per day
  real    :: leaf_fall_rate, root_mort_rate      ! per day
  real    :: BL_u,BL_c
  real    :: retransN  ! retranslocation coefficient of Nitrogen
  logical :: cc_firstday = .false.
  logical :: growingseason
  logical :: TURN_ON_life, TURN_OFF_life, do_fake_phenology

  do_fake_phenology = .FALSE. ! .TRUE.

  retransN = 0.5
  leaf_fall_rate = 0.075; root_mort_rate = 0.0

  vegn%litter = 0   ! daily litter

  ! update vegn GDD and tc_pheno
  vegn%gdd      = vegn%gdd + max(0.0, vegn%tc_daily - 278.15)
  vegn%tc_pheno = vegn%tc_pheno * 0.85 + vegn%Tc_daily * 0.15

! ON and OFF of phenology: change the indicator of growing season for deciduous
  do i = 1,vegn%n_cohorts
     cc => vegn%cohorts(i)
     ! update GDD for each cohort
     cc%gdd = cc%gdd + max(0.0, vegn%tc_daily - 278.15) ! GDD5

     associate (sp => spdata(cc%species) )
!    for evergreen
     if(sp%phenotype==1 .and. cc%status==LEAF_OFF) cc%status=LEAF_ON
!    for deciduous and grasses
     if (do_fake_phenology)then
         TURN_ON_life  = (doy == 90)
         TURN_OFF_life = (doy == 300)
     else
         TURN_ON_life = (sp%phenotype == 0            .and. &
                         cc%status    == LEAF_OFF     .and. &
                         cc%gdd        > sp%gdd_crit  .and. &
                         vegn%tc_pheno > sp%tc_crit_on) !.and. &
             !(sp%lifeform /=0 .OR.(sp%lifeform ==0 .and.cc%layer==1))

         TURN_OFF_life = (sp%phenotype  == 0 .and.        &
                          cc%status     == LEAF_ON  .and. &
                          cc%gdd > sp%gdd_crit+600. .and. &
                          vegn%tc_pheno < sp%tc_crit)
     endif

     cc_firstday = .false.
     if(TURN_ON_life)then
         cc%status = LEAF_ON ! Turn on a growing season
         cc_firstday = .True.

     elseif(TURN_OFF_life )then
        cc%status = LEAF_OFF  ! Turn off a growing season
        cc%gdd   = 0.0        ! Start to counting a new cycle of GDD
        vegn%gdd = 0.0
     endif


!    End a growing season: leaves fall for deciduous
     if(cc%status == LEAF_OFF .AND. cc%bl > 0.0)then
        dBL = cc%bl ! min(leaf_fall_rate * cc%bl_max, cc%bl)
        dBR = 0.0   ! min( root_mort_rate * cc%br_max, cc%br)  ! Just for test: keep roots
        dNL = dBL/sp%CNleaf0
        dNR = dBR/sp%CNroot0

        dAleaf = leaf_area_from_biomass(dBL,cc%species,cc%layer,cc%firstlayer)
        stem_fall = 0.0 ! trees
        if(sp%lifeform==0)then  ! grasses
            stem_fall = MIN(1.0,dBL/cc%bl) * cc%bsw
        endif

!       Retranslocation to NSC and NSN
        cc%nsc = cc%nsc + l_fract  * (dBL + dBR + stem_fall)
        cc%NSN = cc%NSN + retransN * (dNL + dNR)
!       update plant pools
        cc%bl    = cc%bl  - dBL
        cc%br    = cc%br  - dBR
        !cc%bsw   = cc%bsw - stem_fall ! for grasses
        cc%leafN = cc%leafN - dNL ! cc%leafN * (1. - dBL    / cc%bl)
        cc%rootN = cc%rootN - dNR ! cc%rootN * (1.- dBR / cc%br)
!       update NPP for leaves, fine roots, and wood

        cc%NPPleaf = cc%NPPleaf - l_fract * dBL
        cc%NPProot = cc%NPProot - l_fract * dBR
        ! cc%NPPwood = cc%NPPwood - l_fract * stem_fall  ! for grasses
        cc%leafarea= leaf_area_from_biomass(cc%bl,cc%species,cc%layer,cc%firstlayer)
        cc%lai     = cc%leafarea/(cc%crownarea *(1.0-sp%internal_gap_frac))

!       put C and N into soil pools:  Substraction of C and N from leaf and root pools
        loss_coarse  = cc%nindivs * (dBL - dAleaf * LMAmin)
        loss_fine    = cc%nindivs * (dBR + dAleaf * LMAmin) ! + stem_fall ! for grasses
        lossN_coarse = cc%nindivs * (dNL - dAleaf * sp%LNbase)
        lossN_fine   = cc%nindivs * (dNR + dAleaf * sp%LNbase)  !  + stem_fall/sp%CNwood
        vegn%metabolicL = vegn%metabolicL + (1.-l_fract) *  &
                         (fsc_fine * loss_fine + fsc_wood * loss_coarse)
        vegn%structuralL = vegn%structuralL + (1.-l_fract) *     &
                         ((1.-fsc_fine)*loss_fine + (1.-fsc_wood)*loss_coarse)

!       Nitrogen to soil SOMs
        vegn%metabolicN  = vegn%metabolicN + (1.-retransN) *    &
                          (fsc_fine * lossN_fine + fsc_wood * lossN_coarse)
        vegn%structuralN = vegn%structuralN + (1.-retransN) *    &
                          ((1.-fsc_fine) * lossN_fine + (1.-fsc_wood) * lossN_coarse)

!       annual N from plants to soil
        vegn%N_P2S_yr = vegn%N_P2S_yr + lossN_fine + lossN_coarse
     endif
     end associate
     cc => null()
  enddo

end subroutine vegn_phenology

!=======================================================================
! the reproduction of each canopy cohort, yearly time step
! calculate the new cohorts added in this step and states:
! tree density, DBH, woddy and fine biomass
subroutine vegn_reproduction (vegn)
  type(tile_type), intent(inout) :: vegn

! ---- local vars
  type(cohort_type), pointer :: cc ! parent and child cohort pointers
  type(cohort_type), dimension(:),pointer :: ccold, ccnew   ! pointer to old cohort array
  integer,dimension(16) :: reproPFTs
  real,   dimension(16) :: seedC, seedN ! seed pool of productible PFTs
  real :: failed_seeds, N_failedseed !, prob_g, prob_e
  integer :: newcohorts, matchflag, nPFTs ! number of new cohorts to be created
  integer :: nCohorts, istat
  integer :: i, j, k ! cohort indices

! Looping through all reproductable cohorts and Check if reproduction happens
  reproPFTs = -999 ! the code of reproductive PFT
  vegn%totseedC = 0.0
  vegn%totseedN = 0.0
  vegn%totNewCC = 0.0
  vegn%totNewCN = 0.0
  seedC = 0.0
  seedN = 0.0
  nPFTs = 0
  do k=1, vegn%n_cohorts
     cc => vegn%cohorts(k)
     if(cohort_can_reproduce(cc))then
        matchflag = 0
        do i=1,nPFTs
           if(cc%species == reproPFTs(i))then
               seedC(i) = seedC(i) + cc%seedC  * cc%nindivs
               seedN(i) = seedN(i) + cc%seedN  * cc%nindivs
               ! reset parent's seed C and N
               vegn%totSeedC = vegn%totSeedC + cc%seedC  * cc%nindivs
               vegn%totSeedN = vegn%totSeedN + cc%seedN  * cc%nindivs
               cc%seedC = 0.0
               cc%seedN = 0.0

               matchflag = 1
               exit
           endif
        enddo
        if(matchflag==0)then ! when it is a new PFT, put it to the next place
            nPFTs            = nPFTs + 1 ! update the number os reproducible PFTs
            reproPFTs(nPFTs) = cc%species ! PFT number
            seedC(nPFTs)     = cc%seedC * cc%nindivs ! seed carbon
            seedN(nPFTs)     = cc%seedN * cc%nindivs ! seed nitrogen
            vegn%totSeedC = vegn%totSeedC + cc%seedC  * cc%nindivs
            vegn%totSeedN = vegn%totSeedN + cc%seedN  * cc%nindivs
            ! reset parent's seed C and N
            cc%seedC = 0.0
            cc%seedN = 0.0
        endif
     endif ! cohort_can_reproduce
  enddo ! k, vegn%n_cohorts

  ! Generate new cohorts
  newcohorts = nPFTs
  if (newcohorts >= 1) then   ! build new cohorts for seedlings
     ccold => vegn%cohorts ! keep old cohort information
     nCohorts = vegn%n_cohorts + newcohorts
     allocate(ccnew(1:nCohorts), STAT = istat)
     ccnew(1:vegn%n_cohorts) = ccold(1:vegn%n_cohorts) ! copy old cohort information
     vegn%cohorts => ccnew

     deallocate (ccold)

     ! set up new cohorts
     k = vegn%n_cohorts
     do i = 1, newcohorts
        k = k+1 ! increment new cohort index
        cc => vegn%cohorts(k)
        ! Give the new cohort an ID
        cc%ccID = MaxCohortID + i
        ! update child cohort parameters
        associate (sp => spdata(reproPFTs(i))) ! F2003
        ! density
        cc%nindivs = seedC(i)/sp%seedlingsize

        cc%species = reproPFTs(i)
        cc%status  = LEAF_OFF
        cc%firstlayer = 0
        cc%topyear = 0.0
        cc%age     = 0.0

        ! Carbon pools
        cc%bl      = 0.0 * sp%seedlingsize
        cc%br      = 0.0 * sp%seedlingsize
        cc%bsw     = 0.5 * sp%seedlingsize
        cc%bHW     = 0.0 * sp%seedlingsize
        cc%seedC   = 0.0
        cc%nsc     = sp%seedlingsize - cc%bsw !

!!      Nitrogen pools
        cc%leafN  = cc%bl/sp%CNleaf0
        cc%rootN  = cc%br/sp%CNroot0
        cc%sapwN  = cc%bsw/sp%CNsw0
        cc%woodN  = cc%bHW/sp%CNwood0
        cc%seedN  = 0.0
        cc%NSN    = sp%seedlingsize * seedN(i) / seedC(i) -  &
                    (cc%leafN + cc%rootN + cc%sapwN + cc%woodN)

        vegn%totNewCC = vegn%totNewCC + cc%nindivs*(cc%bl + cc%br + cc%bsw + cc%bHW + cc%nsc)
        vegn%totNewCN = vegn%totNewCN + cc%nindivs*(cc%leafN + cc%rootN + cc%sapwN + cc%woodN + cc%NSN)

        call init_cohort_allometry(cc)
        cc%carbon_gain = 0.0
        cc%gpp          = 0.0
        cc%npp          = 0.0
        cc%resp         = 0.0
        cc%resl         = 0.0
        cc%resr         = 0.0
        cc%resg         = 0.0
        cc%annualGPP = 0.0
        cc%annualNPP = 0.0
        cc%NPPleaf = 0.0
        cc%NPProot = 0.0
        cc%NPPwood = 0.0
        cc%N_up_yr = 0.0 ! annual cohort N uptake

!!        !! seeds fail
        !cc%nindivs = cc%nindivs * sp%prob_g * sp%prob_e
!!       put failed seeds to soil carbon pools
!        failed_seeds = 0.0 ! (1. - sp%prob_g*sp%prob_e) * seedC(i)!!

!        vegn%litter = vegn%litter + failed_seeds
!        vegn%metabolicL = vegn%metabolicL +        fsc_fine *failed_seeds
!        vegn%structuralL = vegn%structuralL + (1.0 - fsc_fine)*failed_seeds

!!      Nitrogen of seeds to soil SOMs
!        N_failedseed= 0.0 ! (1.-sp%prob_g*sp%prob_e)   * seedN(i)
!        vegn%metabolicN  = vegn%metabolicN   +        fsc_fine * N_failedseed
!        vegn%structuralN = vegn%structuralN  + (1.0 - fsc_fine)* N_failedseed

!       annual N from plants to soil
 !   vegn%N_P2S_yr = vegn%N_P2S_yr + N_failedseed

        end associate   ! F2003
     enddo
     MaxCohortID = MaxCohortID + newcohorts
     vegn%n_cohorts = k
     ccnew => null()

  endif ! set up new born cohorts

end subroutine vegn_reproduction

!============================================================================
!------------------------Mortality------------------------------------
subroutine vegn_nat_mortality (vegn, deltat)
! TODO: update background mortality rate as a function of wood density (Weng, Jan. 07 2017)
  type(tile_type), intent(inout) :: vegn
  real, intent(in) :: deltat ! time since last mortality calculations, s

  ! ---- local vars
  type(cohort_type), pointer :: cc => null()
  type(spec_data_type),   pointer :: sp
  real :: loss_fine, loss_coarse
  real :: lossN_fine,lossN_coarse
  real :: deathrate ! mortality rate, 1/year
  real :: deadtrees ! number of trees that died over the time step
  real :: DBHtp, tmp
  integer :: i, k

  real, parameter :: min_nindivs = 1e-5 ! 2e-15 ! 1/m. If nindivs is less than this number,
  ! then the entire cohort is killed; 2e-15 is approximately 1 individual per Earth
  ! surface area
  !write(*,*)'total cohorts:', vegn%n_cohorts
  !write(*,'(a81)')'i,PFT,layer,density,layerfrac,dDBH,dbh,height,Acrown,wood,nsc,NPPL,NPPW,aGPP,aNPP'
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     associate ( sp => spdata(cc%species))

     ! mortality rate can be a function of growth rate, age, and environmental
     ! conditions. Here, we only used two constants for canopy layer and under-
     ! story layer (mortrate_d_c and mortrate_d_u)
     ! for trees
     if(cc%layer > 1) then
            tmp = (1 + A_mort*exp(B_mort*(DBH_mort-cc%dbh)) &
                       /(1.0 + exp(B_mort*(DBH_mort-cc%dbh))))
            deathrate = spdata(cc%species)%mortrate_d_u * tmp
     else
            deathrate = spdata(cc%species)%mortrate_d_c !sp%mortrate_d_c

     endif
     deadtrees = cc%nindivs * (1.0-exp(-deathrate*deltat/seconds_per_year)) ! individuals / m2
     ! add dead C from leaf and root pools to fast soil carbon
     loss_coarse  = deadtrees * (cc%bHW + cc%bsw   + cc%bl - cc%leafarea*LMAmin)
     loss_fine    = deadtrees * (cc%nsc + cc%seedC + cc%br + cc%leafarea*LMAmin)
     lossN_coarse = deadtrees * (cc%woodN + cc%sapwN + cc%leafN - sp%LNbase*cc%leafarea)
     lossN_fine   = deadtrees * (cc%rootN + cc%seedN + cc%NSN   + sp%LNbase*cc%leafarea)

     vegn%metabolicL  = vegn%metabolicL +    fsc_fine *loss_fine +    fsc_wood *loss_coarse
     vegn%structuralL = vegn%structuralL + (1.-fsc_fine)*loss_fine + (1.-fsc_wood)*loss_coarse
     ! Nitrogen to soil pools

     vegn%metabolicN = vegn%metabolicN +  fsc_fine * lossN_fine +   &
                                          fsc_wood * lossN_coarse
     vegn%structuralN = vegn%structuralN +(1.-fsc_fine) * lossN_fine +   &
                                          (1.-fsc_wood) * lossN_coarse

!    annual N from plants to soil
     vegn%N_P2S_yr = vegn%N_P2S_yr + lossN_fine + lossN_coarse
     ! Update plant density
     cc%nindivs = cc%nindivs-deadtrees
     end associate
  enddo
end subroutine vegn_nat_mortality

! ============================================================================
function cohort_can_reproduce(cc); logical cohort_can_reproduce
  type(cohort_type), intent(in) :: cc

  associate (sp => spdata(cc%species) )! F2003
  cohort_can_reproduce = (cc%layer == 1 .and. &
                          cc%age   > sp%maturalage.and. &
                          cc%seedC > sp%seedlingsize .and. &
                          cc%seedN > sp%seedlingsize/sp%CNseed0)
  end associate

end function

!========================================================================
! Starvation due to low NSC
subroutine vegn_starvation (vegn)
  type(tile_type), intent(inout) :: vegn

  ! ---- local vars --------
  real :: deathrate ! mortality rate, 1/year
  real :: deadtrees ! number of trees that died over the time step
  real :: loss_fine,loss_coarse
  real :: lossN_fine,lossN_coarse
  integer :: i, k
  type(cohort_type), pointer :: cc
  type(cohort_type), dimension(:),pointer :: ccold, ccnew

  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     associate ( sp => spdata(cc%species)  )

!   Mortality due to starvation
    deathrate = 0.0
!   if (cc%bsw<0 .or. cc%nsc < 0.00001*cc%bl_max .OR.(cc%layer >1 .and. sp%lifeform ==0)) then
    if (cc%bsw<0 .or. cc%nsc < 0.00001*cc%bl_max) then
         deathrate = 1.0
         deadtrees = cc%nindivs * deathrate !individuals / m2

         ! Carbon and Nitrogen from plants to soil pools
         loss_coarse  = deadtrees * (cc%bHW + cc%bsw   + cc%bl - cc%leafarea*LMAmin)
         loss_fine    = deadtrees * (cc%nsc + cc%seedC + cc%br + cc%leafarea*LMAmin)
         lossN_coarse = deadtrees * (cc%woodN + cc%sapwN + cc%leafN - sp%LNbase*cc%leafarea)
         lossN_fine   = deadtrees * (cc%rootN + cc%seedN + cc%NSN   + sp%LNbase*cc%leafarea)

         vegn%metabolicL  = vegn%metabolicL + fsc_fine *loss_fine + fsc_wood *loss_coarse
         vegn%structuralL = vegn%structuralL + (1.0-fsc_fine)*loss_fine + (1.0-fsc_wood)*loss_coarse

         vegn%metabolicN = vegn%metabolicN + &
                           fsc_fine *lossN_fine +    fsc_wood *lossN_coarse
         vegn%structuralN = vegn%structuralN + &
                       (1.-fsc_fine)*lossN_fine +(1.-fsc_wood)*lossN_coarse

         ! annual N from plants to soil
         vegn%N_P2S_yr = vegn%N_P2S_yr + lossN_fine + lossN_coarse

!        update cohort individuals
         cc%nindivs = cc%nindivs*(1.0-deathrate)
     else
         deathrate = 0.0
     endif

     end associate
  enddo
end subroutine vegn_starvation
!=======================================================================
! switch the species of the first cohort to another species
! bugs !!!!!!
 subroutine vegn_species_switch(vegn,N_SP,iyears,FREQ)
  type(tile_type), intent(inout) :: vegn
  integer, intent(in):: N_SP  ! total species in model run settings
  integer, intent(in):: iyears
  integer, intent(in):: FREQ  ! frequency of species switching

  ! ---- local vars --------
  real :: loss_fine,loss_coarse
  real :: lossN_fine,lossN_coarse
  integer :: i, k
  type(cohort_type), pointer :: cc

     cc => vegn%cohorts(1)
     associate (sp => spdata(cc%species)) ! F2003
     if(cc%bl > 0.0) then ! remove all leaves to keep mass balance
        loss_coarse  = cc%nindivs * (cc%bl - cc%leafarea*LMAmin)
        loss_fine    = cc%nindivs *  cc%leafarea*LMAmin
        lossN_coarse = cc%nindivs * (cc%leafN - sp%LNbase*cc%leafarea)
        lossN_fine   = cc%nindivs *  sp%LNbase*cc%leafarea
        ! Carbon to soil pools
        vegn%metabolicL  = vegn%metabolicL  + fsc_fine *loss_fine + &
                                              fsc_wood *loss_coarse
        vegn%structuralL = vegn%structuralL + (1.0-fsc_fine)*loss_fine + &
                                              (1.0-fsc_wood)*loss_coarse
        ! Nitrogen to soil pools
        vegn%metabolicN = vegn%metabolicN + fsc_fine  *lossN_fine +   &
                                        fsc_wood *lossN_coarse
        vegn%structuralN = vegn%structuralN +(1.-fsc_fine) *lossN_fine +   &
                                      (1.-fsc_wood)*lossN_coarse
        ! annual N from plants to soil
        vegn%N_P2S_yr = vegn%N_P2S_yr + lossN_fine + lossN_coarse
        ! remove leaves
        cc%bl = 0.0
     endif
     end associate
     ! Change species
     cc%species = mod(iyears/FREQ,N_SP)+2

 end subroutine vegn_species_switch

! ============================================================================
! Arrange crowns into canopy layers according to their height and crown areas.
subroutine relayer_cohorts (vegn)
  type(tile_type), intent(inout) :: vegn ! input cohorts

  ! ---- local constants
  real, parameter :: tolerance = 1e-4
  real, parameter :: layer_vegn_cover = 1.0
  ! ---- local vars
  integer :: idx(vegn%n_cohorts) ! indices of cohorts in decreasing height order
  integer :: i ! new cohort index
  integer :: k ! old cohort index
  integer :: L ! layer index (top-down)
  integer :: N0,N1 ! initial and final number of cohorts
  real    :: frac ! fraction of the layer covered so far by the canopies
  type(cohort_type), pointer :: cc(:),new(:)
  real    :: nindivs

!  rand_sorting = .TRUE. ! .False.

  ! rank cohorts in descending order by height. For now, assume that they are
  ! in order
  N0 = vegn%n_cohorts; cc=>vegn%cohorts
  call rank_descending(cc(1:N0)%height,idx)

  ! calculate max possible number of new cohorts : it is equal to the number of
  ! old cohorts, plus the number of layers -- since the number of full layers is
  ! equal to the maximum number of times an input cohort can be split by a layer
  ! boundary.
  N1 = vegn%n_cohorts + int(sum(cc(1:N0)%nindivs*cc(1:N0)%crownarea))
  allocate(new(N1))

  ! copy cohort information to the new cohorts, splitting the old cohorts that
  ! stride the layer boundaries
  i = 1 ; k = 1 ; L = 1 ; frac = 0.0 ; nindivs = cc(idx(k))%nindivs
  do
     new(i)         = cc(idx(k))
     new(i)%nindivs = min(nindivs,(layer_vegn_cover-frac)/cc(idx(k))%crownarea)
     new(i)%layer   = L
     if (L==1) new(i)%firstlayer = 1
!    if (L>1)  new(i)%firstlayer = 0  ! switch off "push-down effects"
     frac = frac+new(i)%nindivs*new(i)%crownarea
     nindivs = nindivs - new(i)%nindivs

     if (abs(nindivs*cc(idx(k))%crownarea)<tolerance) then
       new(i)%nindivs = new(i)%nindivs + nindivs ! allocate the remainder of individuals to the last cohort
       if (k==N0) exit ! end of loop
       k = k+1 ; nindivs = cc(idx(k))%nindivs  ! go to the next input cohort
     endif

     if (abs(layer_vegn_cover - frac)<tolerance) then
       L = L+1 ; frac = 0.0              ! start new layer
     endif
!     write(*,*)i, new(i)%layer
     i = i+1
  enddo

  ! replace the array of cohorts
  deallocate(vegn%cohorts)
  vegn%cohorts => new ; vegn%n_cohorts = i
  ! update layer fraction for each cohort
  do i=1, vegn%n_cohorts
     vegn%cohorts(i)%layerfrac = vegn%cohorts(i)%nindivs * vegn%cohorts(i)%crownarea
  enddo

end subroutine relayer_cohorts

! ============= Plant physiology =============================================

! ============================================================================
subroutine vegn_C_gain(vegn,forcing)
!@ Calculate daily carbon gain per tree based on V and self-shading of leaves
!@ It is used to generate daily GPP (photosynthesis)
!@ This subroutine can be replaced by a photosynthesis model working at hourly
!@ time scale

  type(tile_type), intent(inout) :: vegn
  type(climate_data_type),intent(in):: forcing(:)

  !-------local var
  type(cohort_type), pointer :: cc    ! current cohort
  logical :: extra_light_in_lower_layers
  real :: f_light(10)      ! light fraction of each layer
  real :: V_annual     ! max V for each layer
  real :: f_gap ! additional GPP for lower layer cohorts due to gaps
  integer :: i, layer

  f_gap = 0.2 ! 0.1
! update accumulative LAI for each corwn layer
  vegn%CAI      = 0.0
  vegn%LAI      = 0.0
  vegn%LAIlayer = 0.0
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     associate ( sp => spdata(cc%species) )
     cc%leafarea=leaf_area_from_biomass(cc%bl,cc%species,cc%layer,cc%firstlayer)
     cc%lai     = cc%leafarea/(cc%crownarea *(1.0-sp%internal_gap_frac))

     layer = Max (1, Min(cc%layer,9)) + 1 ! next layer
     ! LAI above this layer: Layer1: 0; Layer2: LAI of Layer1 cohorts; ...
     vegn%LAIlayer(layer) = vegn%LAIlayer(layer) + cc%leafarea * cc%nindivs
     vegn%LAI = vegn%LAI + cc%leafarea  * cc%nindivs
     vegn%CAI = vegn%CAI + cc%crownarea * cc%nindivs
     END associate
  enddo
  ! Light fraction
  f_light(1) = 1.0
  do i =2, layer !MIN(int(vegn%CAI+1.0),9)
      f_light(i) = f_light(i-1) * &
                  (exp(-0.5*vegn%LAIlayer(i))*(1.-f_gap) + f_gap)
  enddo

! Assumption: no gaps  --> GPP of understory trees is too low!
! Assimilation of carbon for each cohort considering their light envrionment
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)

     layer = Max (1, Min(cc%layer,9))
     ! Photosynthesis can be calculated by a photosynthesis model
     V_annual = f_light(layer) * spdata(cc%species)%Vannual

     if(cc%status == LEAF_ON) then
         ! Add temperature response function of photosynthesis
         cc%gpp = V_annual/0.5 * (1.0 - exp(-0.5 * cc%LAI))  &
                  * cc%crownarea * dt_fast_yr                &
                  * exp(9000.0 * (1./298.16 - 1./vegn%tc_daily)) ! temperature response function
                 ! =1.2/0.5/cc%layer**2 * (1.0 - exp(-0.5* cc%LAI)) & ! 0.5 & !
                 ! * cc%crownarea * dt_fast_yr

       ! kgC tree-1 time step-1
     else
         cc%gpp = 0.0
     endif
     ! Update NSC
     cc%nsc = cc%nsc + cc%gpp
  enddo

  cc => null()
end subroutine vegn_C_gain
!=========================================================================

subroutine carbon_for_growth(cc)
  type(cohort_type), intent(inout) :: cc

  !-------local var
  real :: NSC_supply,LR_demand,LR_deficit
  real :: NSCtarget
  !real :: LR_growth,WS_growth
  real :: R_days,fNSC,fLFR,fsup

! Grab carbon from NSC pool and put them into "carbon_gain"
!   modified 9/3/2013 based on Steve's suggestions
   associate ( sp => spdata(cc%species) )
   if(cc%status == LEAF_ON) then
      R_days = 5.0
      fNSC = 0.05 * days_per_year * dt_fast_yr ! 0.2(daily) -->0.2/24 (hourly) 2014-10-22
      fLFR = 0.2 * days_per_year * dt_fast_yr
      fsup = dt_fast_yr/spdata(cc%species)%tauNSC  ! 0.05
      NSCtarget = 3.0 * (cc%bl_max + cc%br_max)

      LR_demand = 0.0; NSC_supply= 0.0
      IF(cc%nsc > 0. .AND. cc%status == LEAF_ON)then
          LR_deficit = max(cc%bl_max+cc%br_max-cc%bl-cc%br,0.0)
          LR_demand  = min(fLFR*LR_deficit, fNSC*cc%nsc)
          NSC_supply = cc%nsc * fsup ! max((cc%nsc - NSCtarget)*fsup,0.0) ! Weng 2014-01-23 for smoothing dDBH
      ENDIF
      cc%nsc  = cc%nsc - (LR_demand + NSC_supply)

!     Deduct growth respirtion from (LR_demand + NSC_supply)
      cc%resg    = GR_factor/(1.+GR_factor) * (LR_demand + NSC_supply) ! kgC tree-1 step-1
      LR_demand  = LR_demand  /(1.+ GR_factor) ! for building up tissues
      NSC_supply = NSC_supply /(1.+ GR_factor)
      ! carbon_gain is used to drive plant growth and reproduction
      cc%carbon_gain = cc%carbon_gain + (LR_demand + NSC_supply) !
   endif
   end associate

end subroutine carbon_for_growth

! ============================================================================
subroutine plant_respiration(cc, tsoil)
  type(cohort_type), intent(inout) :: cc
  real, intent(in) :: tsoil

  real :: tf,tfs ! thermal inhibition factors for above- and below-ground biomass
  real :: r_leaf, r_stem, r_root
  real :: Acambium  ! cambium area, m2/tree
  ! real :: LeafN     ! leaf nitrogen, kgN/Tree
  real :: NSCtarget ! used to regulation respiration rate
  real :: r_Nfix    ! respiration due to N fixation

  integer :: sp ! shorthand for cohort species
  sp = cc%species
  ! temperature response function
  tf  = exp(9000.0*(1.0/298.16-1.0/tsoil))
!  tfs = thermal_inhibition(tsoil)  ! original
  tfs = tf ! Rm_T_response_function(tsoil) ! Weng 2014-01-14
! With nitrogen model, leaf respiration is a function of leaf nitrogen
  NSCtarget = 3.0 * (cc%bl_max + cc%br_max)
  Acambium = PI * cc%DBH * cc%height * 1.2

  ! Facultive Nitrogen fixation
  !if(cc%NSN < cc%NSNmax .and. cc%NSC > 0.5 * NSCtarget)then
  !   cc%fixedN = spdata(sp)%NfixRate0 * cc%br * tf * dt_fast_yr ! kgN tree-1 step-1
  !else
  !   cc%fixedN = 0.0 ! spdata(sp)%NfixRate0 * cc%br * tf * dt_fast_yr ! kgN tree-1 step-1
  !endif

  ! Obligate Nitrogen Fixation
  cc%fixedN = spdata(sp)%NfixRate0 * cc%br * tf * dt_fast_yr ! kgN tree-1 step-1

  ! LeafN    = spdata(sp)%LNA * cc%leafarea
  r_stem   = spdata(sp)%gamma_SW * Acambium * tf * dt_fast_yr ! kgC tree-1 step-1
  r_root   = spdata(sp)%gamma_FR * cc%rootN * tf * dt_fast_yr ! root respiration ~ root N
  r_leaf   = spdata(sp)%gamma_LN * cc%leafN * tf * dt_fast_yr
  r_Nfix   = spdata(sp)%NfixCost0 * cc%fixedN

  cc%resp = (r_leaf + r_stem + r_root + r_Nfix) !* max(0.0, cc%nsc/NSCtarget)
  cc%resl = r_leaf !* max(0.0, cc%nsc/NSCtarget)
  cc%resr = r_root + r_Nfix !* max(0.0, cc%nsc/NSCtarget)
  ! Update NSC
  cc%nsc = cc%nsc - cc%resp
  cc%NSN = cc%NSN + cc%fixedN

end subroutine plant_respiration

! ============================================================================
subroutine vegn_leaf_fine_root_turnover(vegn, tsoil, theta)
  type(tile_type), intent(inout) :: vegn
  real, intent(in) :: tsoil ! average temperature of soil, deg K
  real, intent(in) :: theta ! average soil wetness, unitless

  !-------local var
  type(cohort_type), pointer :: cc    ! current cohort
  real :: loss_coarse, loss_fine, lossN_coarse, lossN_fine
  real :: dBL, dBR  ! leaf and fine root carbon tendencies
  real :: dNL, dNR  ! leaf and fine root nitrogen tendencies
  real :: dAleaf ! leaf area decrease due to dBL
  real :: retransN = 0.5  ! retranslocation coefficient of Nitrogen
  integer :: i

  retransN = 0.5
  ! update plant carbon and nitrogen for all cohorts
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     associate ( sp => spdata(cc%species) )
!    Turnover of leaves and roots regardless of the STATUS of leaf
!    longevity. Deciduous: 0; Evergreen 0.035/LMa
!    root turnover
     dBL = cc%bl    * sp%alpha_L  * dt_fast_yr
     dBR = cc%br    * sp%alpha_FR * dt_fast_yr
     dNL = cc%leafN * sp%alpha_L  * dt_fast_yr
     dNR = cc%rootN * sp%alpha_FR * dt_fast_yr
     dAleaf = leaf_area_from_biomass(dBL,cc%species,cc%layer,cc%firstlayer)

!    Retranslocation to NSC and NSN
     cc%nsc = cc%nsc + l_fract  * (dBL + dBR)
     cc%NSN = cc%NSN + retransN * (dNL + dNR)
!    update plant pools
     cc%bl    = cc%bl    - dBL
     cc%br    = cc%br    - dBR
     cc%leafN = cc%leafN - dNL
     cc%rootN = cc%rootN - dNR

!    update leaf area and LAI
     cc%leafarea= leaf_area_from_biomass(cc%bl,cc%species,cc%layer,cc%firstlayer)
     cc%lai     = cc%leafarea/(cc%crownarea *(1.0-sp%internal_gap_frac))

!    update NPP for leaves, fine roots, and wood
     cc%NPPleaf = cc%NPPleaf - l_fract * dBL
     cc%NPProot = cc%NPProot - l_fract * dBR

!    put C and N into soil pools
     loss_coarse  = cc%nindivs * (dBL - dAleaf * LMAmin)
     loss_fine    = cc%nindivs * (dBR + dAleaf * LMAmin)
     lossN_coarse = cc%nindivs * (dNL - dAleaf * sp%LNbase)
     lossN_fine   = cc%nindivs * (dNR + dAleaf * sp%LNbase)

     vegn%metabolicL = vegn%metabolicL + (1.-l_fract) *  &
                        (fsc_fine * loss_fine + fsc_wood * loss_coarse)
     vegn%structuralL = vegn%structuralL + (1.-l_fract) *     &
                         ((1.-fsc_fine)*loss_fine + (1.-fsc_wood)*loss_coarse)

!    Nitrogen to soil SOMs
     vegn%metabolicN  = vegn%metabolicN + (1.-retransN) *    &
                          (fsc_fine * lossN_fine + fsc_wood * lossN_coarse)
     vegn%structuralN = vegn%structuralN + (1.-retransN) *    &
                          ((1.-fsc_fine) * lossN_fine + (1.-fsc_wood) * lossN_coarse)

!    annual N from plants to soil
     vegn%N_P2S_yr = vegn%N_P2S_yr + lossN_fine + lossN_coarse

    END ASSOCIATE
  enddo

end subroutine vegn_leaf_fine_root_turnover
!=====================================================
! Weng, 2016-11-28
subroutine vegn_N_uptake(vegn, tsoil, theta)
  type(tile_type), intent(inout) :: vegn
  real, intent(in) :: tsoil ! average temperature of soil, deg K
  real, intent(in) :: theta ! average soil wetness, unitless

  !-------local var
  type(cohort_type),pointer :: cc
  real    :: rho_N_up0 = 0.02 ! hourly N uptake rate, fraction of the total mineral N
  real    :: N_roots0  = 0.1  ! root biomass at half max N-uptake rate,kg C m-2
  real    :: totNup    ! kgN m-2
  real    :: avgNup
  real    :: rho_N_up,N_roots   ! actual N uptake rate
  logical :: NSN_not_full
  integer :: i

!! Nitrogen uptake parameter
! It considers competition here. How much N one can absorp depends on
! how many roots it has and how many roots other individuals have.
  N_Roots  = 0.0
  vegn%N_uptake = 0.0
  if(vegn%mineralN > 0.0)then
     do i = 1, vegn%n_cohorts
        cc => vegn%cohorts(i)
        associate (sp => spdata(cc%species))
!!       A scheme for deciduous to get enough N:
        cc%NSNmax = 0.2 * cc%crownarea  ! 5*(cc%bl_max/sp%CNleaf0 + cc%br_max/sp%CNroot0)) !
        NSN_not_full = (cc%NSN < cc%NSNmax) !
        if(NSN_not_full) N_Roots = N_Roots + cc%br * cc%nindivs

        end associate
     enddo
     ! M-M equation for Nitrogen absoption, McMurtrie et al. 2012, Ecology & Evolution
     ! rate at given root biomass and period of time
     if(N_roots>0.0)then
        ! Add a temperature response equation herefor rho_N_up0 (Zhu Qing 2016)
        rho_N_up = 1.-exp(-rho_N_up0 * N_roots/(N_roots0+N_roots) * hours_per_year * dt_fast_yr) ! rate at given root density and time period
        totNup = rho_N_up * vegn%mineralN  &
                * exp(9000.0 * (1./298.16 - 1./tsoil)) ! kgN m-2 time step-1
        avgNup = totNup / N_roots ! kgN time step-1 kg roots-1
        ! Nitrogen uptaken by each cohort, N_uptake
        vegn%N_uptake = 0.0
        do i = 1, vegn%n_cohorts
           cc => vegn%cohorts(i)
           associate ( sp => spdata(cc%species) )
           NSN_not_full = (cc%NSN < cc%NSNmax)
           if(NSN_not_full)then
               cc%N_uptake = cc%br  * avgNup  !
               cc%nsn      = cc%nsn + cc%N_uptake
               cc%N_up_yr  = cc%N_up_yr + cc%N_uptake/cc%crownarea

               ! subtract N from mineral N
               vegn%mineralN = vegn%mineralN - cc%N_uptake * cc%nindivs
               vegn%N_uptake = vegn%N_uptake + cc%N_uptake * cc%nindivs
               vegn%accu_Nup  = vegn%accu_Nup + cc%N_uptake * cc%nindivs

           endif
           end associate
        enddo
        cc =>null()
     endif ! N_roots>0

  endif

end subroutine vegn_N_uptake
! ============================================================================
! Nitrogen mineralization and immoblization with microbial C & N pools
! it's a new decomposition model with coupled C & N pools and variable
! carbon use efficiency
subroutine SOMdecomposition(vegn, soilt, theta)
  type(tile_type), intent(inout) :: vegn
  real                , intent(in)    :: soilt ! soil temperature, deg K
  real                , intent(in)    :: theta

  real :: CUE0=0.4  ! default microbial CUE
  real :: phoMicrobial = 2.5 ! turnover rate of microbes (yr-1)
  real :: CUEfast,CUEslow
  real :: CNm = 10.0  ! Microbial C/N ratio
  real :: NforM, fNM=0.0  ! mineral N available for microbes
  real :: micr_C_loss, fast_L_loss, slow_L_loss
  real :: runoff ! kg m-2 /step
  real :: N_loss
  real :: DON_fast,DON_slow,DON_loss ! Dissolved organic N loss, kg N m-2 step-1
  real :: fDON=0.0   ! 0.02     ! fractio of DON production in decomposition
  real :: fast_N_free
  real :: slow_N_free
  real :: CNfast, CNslow
  real :: A  ! decomp rate reduction due to moisture and temperature

!  runoff = vegn%Wrunoff * 365*24*3600 *dt_fast_yr !kgH2O m-2 s-1 ->kg m-2/time step
  runoff = vegn%Wrunoff * dt_fast_yr !kgH2O m-2 yr-1 ->kgH2O m-2/time step
! CN ratios of soil C pools
  CNfast = vegn%metabolicL/vegn%metabolicN
  CNslow = vegn%structuralL/vegn%structuralN

!! C decomposition
!  A=A_function(soilt,theta)
!  micr_C_loss = vegn%microbialC *A*phoMicrobial* dt_fast_yr
!  fast_L_loss = vegn%metabolicL*A*K1           * dt_fast_yr
!  slow_L_loss = vegn%structuralL*A*K2          * dt_fast_yr

! C decomposition
  A=A_function(soilt,theta)
  micr_C_loss = vegn%microbialC * (1.0 - exp(-A*phoMicrobial* dt_fast_yr))
  fast_L_loss = vegn%metabolicL * (1.0 - exp(-A*K1          * dt_fast_yr))
  slow_L_loss = vegn%structuralL* (1.0 - exp(-A*K2          * dt_fast_yr))

! Carbon use efficiencies of microbes
  NforM = fNM * vegn%mineralN
  CUEfast = MIN(CUE0,CNm*(fast_L_loss/CNfast + NforM)/fast_L_loss)
  CUEslow = MIN(CUE0,CNm*(slow_L_loss/CNslow + NforM)/slow_L_loss)

! update C and N pools
! Carbon pools
  vegn%microbialC  = vegn%microbialC - micr_C_loss &
                    + fast_L_loss * CUEfast &
                    + slow_L_loss * CUEslow
  vegn%metabolicL = vegn%metabolicL - fast_L_loss
  vegn%structuralL = vegn%structuralL - slow_L_loss

! Find papers about soil DON losses
! DON loss, revised by Weng. 2016-03-03  ??
  fDON        = 0.25 ! 0.25 ! * dt_fast_yr ! 0.05 !* dt_fast_yr
  runoff      = 0.2 ! 0.2 ! mm day-1
  ! Assume it is proportional to decomposition rates
  ! Find some papers!!
  DON_fast    = fDON * fast_L_loss/CNfast * (etaN*runoff)
  DON_slow    = fDON * slow_L_loss/CNslow * (etaN*runoff)
  DON_loss    = DON_fast + DON_slow

! Update Nitrogen pools
  vegn%microbialN= vegn%microbialC/CNm
  vegn%metabolicN  = vegn%metabolicN  - fast_L_loss/CNfast - DON_fast
  vegn%structuralN = vegn%structuralN - slow_L_loss/CNslow - DON_slow

! Mixing of microbes to litters
  vegn%metabolicL   = vegn%metabolicL + MLmixRatio*fast_L_loss * CUEfast
  vegn%metabolicN   = vegn%metabolicN + MLmixRatio*fast_L_loss * CUEfast/CNm

  vegn%structuralL = vegn%structuralL + MLmixRatio*slow_L_loss * CUEslow
  vegn%structuralN = vegn%structuralN + MLmixRatio*slow_L_loss * CUEslow/CNm

  vegn%microbialC  = vegn%microbialC  - MLmixRatio*(fast_L_loss*CUEfast+slow_L_loss*CUEslow)
  vegn%microbialN  = vegn%microbialC/CNm

! update mineral N pool (mineralN)
  fast_N_free = MAX(0.0, fast_L_loss*(1./CNfast - CUEfast/CNm))
  slow_N_free = MAX(0.0, slow_L_loss*(1./CNslow - CUEslow/CNm))

! N_loss = MAX(0.,vegn%mineralN)        * A * K_nitrogen * dt_fast_yr
  N_loss = MAX(0.,vegn%mineralN) * (1.0-exp(-etaN*runoff - A*K_nitrogen * dt_fast_yr))
  vegn%Nloss_yr = vegn%Nloss_yr + N_loss + DON_loss

  vegn%mineralN = vegn%mineralN - N_loss       &
                  + vegn%N_input * dt_fast_yr  &
                  + fast_N_free + slow_N_free  &
                  + micr_C_loss/CNm
  vegn%annualN   = vegn%annualN - N_loss       &
                  + vegn%N_input * dt_fast_yr  &
                  + fast_N_free + slow_N_free  &
                  + micr_C_loss/CNm

! Check if soil C/N is above CN0
  fast_N_free = MAX(0., vegn%metabolicN  - vegn%metabolicL/CN0metabolicL)
  slow_N_free = MAX(0., vegn%structuralN - vegn%structuralL/CN0structuralL)
  vegn%metabolicN  = vegn%metabolicN  - fast_N_free
  vegn%structuralN = vegn%structuralN - slow_N_free
  vegn%mineralN    = vegn%mineralN + fast_N_free + slow_N_free
  vegn%annualN     = vegn%annualN  + fast_N_free + slow_N_free

! Heterotrophic respiration: decomposition of litters and SOM, kgC m-2 yr-1
  vegn%rh =   (micr_C_loss + fast_L_loss*(1.-CUEfast)+ slow_L_loss*(1.-CUEslow)) &
              /dt_fast_yr

end subroutine SOMdecomposition

! ============================================================================
! The combined reduction in decomposition rate as a funciton of TEMP and MOIST
! Based on CENTURY Parton et al 1993 GBC 7(4):785-809 and Bolker's copy of
! CENTURY code
function A_function(soilt, theta) result(A)
  real :: A                 ! return value, resulting reduction in decomposition rate
  real, intent(in) :: soilt ! effective temperature for soil carbon decomposition
  real, intent(in) :: theta

  real :: soil_temp ! temperature of the soil, deg C
  real :: Td        ! rate multiplier due to temp
  real :: Wd        ! rate reduction due to mositure

  ! coefficeints and terms used in temperaturex term
  real :: Topt,Tmax,t1,t2,tshl,tshr

  soil_temp = soilt-273.16

  ! EFFECT OF TEMPERATURE , ! from Bolker's century code
  Tmax=45.0;
  if (soil_temp > Tmax) soil_temp = Tmax;
  Topt=35.0;
  tshr=0.2; tshl=2.63;
  t1=(Tmax-soil_temp)/(Tmax-Topt);
  t2=exp((tshr/tshl)*(1.-t1**tshl));
  Td=t1**tshr*t2;

  if (soil_temp > -10) Td=Td+0.05;
  if (Td > 1.) Td=1.;

  ! EFFECT OF MOISTURE
  ! Linn and Doran, 1984, Soil Sci. Amer. J. 48:1267-1272
  ! This differs from the Century Wd
  ! was modified by slm/ens based on the figures from the above paper
  !     (not the reported function)

  if(theta <= 0.3) then
     Wd = 0.2;
  else if(theta <= 0.6) then
     Wd = 0.2+0.8*(theta-0.3)/0.3
  else
     Wd = 1.0 ! exp(2.3*(0.6-theta)); ! Weng, 2016-11-26
  endif

  A = (Td*Wd); ! the combined (multiplicative) effect of temp and water
               ! on decomposition rates
end function A_function

! =================== Cohort management ================================
! ======================================================================

subroutine rank_descending(x,idx)
! ranks array x in descending order: on return, idx() contains indices
! of elements of array x in descending order of x values. These codes
! are from Sergey Malyshev (LM3PPA, Weng et al. 2015 Biogeosciences)
   real,    intent(in)  :: x(:)
   integer, intent(out) :: idx(:)
   integer :: i,n
   integer, allocatable :: t(:)

   n = size(x)
   do i = 1,n
      idx(i) = i
   enddo

   allocate(t((n+1)/2))
   call mergerank(x,idx,n,t)
   deallocate(t)
end subroutine

! =====================================================================
! based on:
! http://rosettacode.org/wiki/Sorting_algorithms/Merge_sort#Fortran
subroutine merge(x,a,na,b,nb,c,nc)
   integer, intent(in) :: na,nb,nc ! Normal usage: NA+NB = NC
   real, intent(in)       :: x(*)
   integer, intent(in)    :: a(na)    ! B overlays C(NA+1:NC)
   integer, intent(in)    :: b(nb)
   integer, intent(inout) :: c(nc)

   integer :: i,j,k

   i = 1; j = 1; k = 1;
   do while(i <= na .and. j <= nb)
      if (x(a(i)) >= x(b(j))) then
         c(k) = a(i) ; i = i+1
      else
         c(k) = b(j) ; j = j+1
      endif
      k = k + 1
   enddo
   do while (i <= na)
      c(k) = a(i) ; i = i + 1 ; k = k + 1
   enddo
end subroutine merge

recursive subroutine mergerank(x,a,n,t)
  integer, intent(in) :: n
  real,    intent(in) :: x(*)
  integer, dimension(n), intent(inout) :: a
  integer, dimension((n+1)/2), intent (out) :: t

  integer :: na,nb
  integer :: v

  if (n < 2) return
  if (n == 2) then
     if ( x(a(1)) < x(a(2)) ) then
        v = a(1) ; a(1) = a(2) ; a(2) = v
     endif
     return
  endif
  na=(n+1)/2
  nb=n-na

  call mergerank(x,a,na,t)
  call mergerank(x,a(na+1),nb,t)

  if (x(a(na)) < x(a(na+1))) then
     t(1:na)=a(1:na)
     call merge(x,t,na,a(na+1),nb,a,n)
  endif
end subroutine mergerank

!============================================================================
! Merge similar cohorts in a tile
subroutine vegn_mergecohorts(vegn)
  type(tile_type), intent(inout) :: vegn

! ---- local vars
  type(cohort_type), pointer :: cc(:) ! array to hold new cohorts
  logical :: merged(vegn%n_cohorts)        ! mask to skip cohorts that were already merged
  real, parameter :: mindensity = 1.0E-6
  integer :: i,j,k

  allocate(cc(vegn%n_cohorts))
  merged(:)=.FALSE. ; k = 0
  do i = 1, vegn%n_cohorts
     if(merged(i)) cycle ! skip cohorts that were already merged
     k = k+1
     cc(k) = vegn%cohorts(i)
     ! try merging the rest of the cohorts into current one
     do j = i+1, vegn%n_cohorts
        if (merged(j)) cycle ! skip cohorts that are already merged
        if (cohorts_can_be_merged(vegn%cohorts(j),cc(k))) then
           call merge_cohorts(vegn%cohorts(j),cc(k))
           merged(j) = .TRUE.
        endif
     enddo
  enddo
  ! at this point, k is the number of new cohorts
  vegn%n_cohorts = k
  deallocate(vegn%cohorts)
  vegn%cohorts=>cc

end subroutine vegn_mergecohorts

! ============================================================================
! kill low density cohorts, a new function seperated from vegn_mergecohorts
! Weng, 2014-07-22
subroutine kill_lowdensity_cohorts(vegn)
  type(tile_type), intent(inout) :: vegn

! ---- local vars
  type(cohort_type), pointer :: cp, cc(:) ! array to hold new cohorts
  logical :: merged(vegn%n_cohorts)        ! mask to skip cohorts that were already merged
  real, parameter :: mindensity = 0.25E-4
  real :: loss_fine,loss_coarse
  real :: lossN_fine,lossN_coarse
  integer :: i,j,k


 ! calculate the number of cohorts with indivs>mindensity
  k = 0
  do i = 1, vegn%n_cohorts
     if (vegn%cohorts(i)%nindivs > mindensity) k=k+1
  enddo
  if (k==0) write(*,*)'kill_lowdensity_cohorts: ','All cohorts have died'

  ! exclude cohorts that have low individuals
  if (k < vegn%n_cohorts) then
     allocate(cc(k))
     k=0
     do i = 1,vegn%n_cohorts
        cp =>vegn%cohorts(i)
        associate(sp=>spdata(cp%species))
        if (cp%nindivs > mindensity) then
           k=k+1
           cc(k) = cp
        else
           ! Carbon and Nitrogen from plants to soil pools
           loss_coarse  = cp%nindivs * (cp%bHW + cp%bsw   + cp%bl - cp%leafarea*LMAmin)
           loss_fine    = cp%nindivs * (cp%nsc + cp%seedC + cp%br + cp%leafarea*LMAmin)
           lossN_coarse = cp%nindivs * (cp%woodN + cp%sapwN + cp%leafN - sp%LNbase*cp%leafarea)
           lossN_fine   = cp%nindivs * (cp%rootN + cp%seedN + cp%NSN   + sp%LNbase*cp%leafarea)

           vegn%metabolicL  = vegn%metabolicL + fsc_fine *loss_fine + fsc_wood *loss_coarse
           vegn%structuralL = vegn%structuralL + (1.0-fsc_fine)*loss_fine + (1.0-fsc_wood)*loss_coarse

           vegn%metabolicN = vegn%metabolicN + &
                           fsc_fine *lossN_fine +    fsc_wood *lossN_coarse
           vegn%structuralN = vegn%structuralN + &
                       (1.-fsc_fine)*lossN_fine +(1.-fsc_wood)*lossN_coarse
           ! annual N from plants to soil
           vegn%N_P2S_yr = vegn%N_P2S_yr + lossN_fine + lossN_coarse
        endif
        end associate
     enddo
     vegn%n_cohorts = k
     deallocate (vegn%cohorts)
     vegn%cohorts=>cc
  endif
end subroutine kill_lowdensity_cohorts


! ============================================================================
subroutine merge_cohorts(c1,c2)
  type(cohort_type), intent(in) :: c1
  type(cohort_type), intent(inout) :: c2

  real :: x1, x2 ! normalized relative weights

  if(c1%nindivs > 0.0 .or. c2%nindivs > 0.0)then
     x1 = c1%nindivs/(c1%nindivs+c2%nindivs)
     x2 = c2%nindivs/(c1%nindivs+c2%nindivs)
  !else
  !   x1 = 0.5
  !   x2 = 0.5
  !endif
  ! update number of individuals in merged cohort
     c2%nindivs = c1%nindivs + c2%nindivs
  !  Carbon
     c2%bl  = x1*c1%bl  + x2*c2%bl
     c2%br  = x1*c1%br  + x2*c2%br
     c2%bsw = x1*c1%bsw + x2*c2%bsw
     c2%bHW = x1*c1%bHW + x2*c2%bHW
     c2%seedC = x1*c1%seedC + x2*c2%seedC
     c2%nsc = x1*c1%nsc + x2*c2%nsc
  !   Allometry
     c2%dbh = x1*c1%dbh + x2*c2%dbh
     c2%height = x1*c1%height + x2*c2%height
     c2%crownarea = x1*c1%crownarea + x2*c2%crownarea
     c2%age = x1*c1%age + x2*c2%age
     c2%carbon_gain = x1*c1%carbon_gain + x2*c2%carbon_gain
     c2%topyear = x1*c1%topyear + x2*c2%topyear

  !  Nitrogen
     c2%leafN = x1*c1%leafN + x2*c2%leafN
     c2%rootN = x1*c1%rootN + x2*c2%rootN
     c2%sapwN = x1*c1%sapwN + x2*c2%sapwN
     c2%woodN = x1*c1%woodN + x2*c2%woodN
     c2%seedN = x1*c1%seedN + x2*c2%seedN
     c2%NSN   = x1*c1%NSN   + x2*c2%NSN

  !  calculate the resulting dry heat capacity
     c2%leafarea = leaf_area_from_biomass(c2%bl, c2%species, c2%layer, c2%firstlayer)
  endif
end subroutine merge_cohorts

! ============================================================================
function cohorts_can_be_merged(c1,c2); logical cohorts_can_be_merged
   type(cohort_type), intent(in) :: c1,c2

   real, parameter :: mindensity = 1.0E-4
   logical :: sameSpecies, sameLayer, sameSize, sameSizeTree, sameSizeGrass, lowDensity

   sameSpecies  = c1%species == c2%species
   sameLayer    = (c1%layer == c2%layer) ! .and. (c1%firstlayer == c2%firstlayer)
   sameSizeTree = (spdata(c1%species)%lifeform > 0).and.  &
                  (spdata(c2%species)%lifeform > 0).and.  &
                 ((abs(c1%DBH - c2%DBH)/c2%DBH < 0.2 ) .or.  &
                  (abs(c1%DBH - c2%DBH)        < 0.001))  ! it'll be always true for grasses
   sameSizeGrass= (spdata(c1%species)%lifeform ==0) .and. &
                  (spdata(c2%species)%lifeform ==0) .and. &
                 (((c1%DBH == c2%DBH) .And.(c1%nsc == c2%nsc)) .OR. &
                    c1%topyear==c2%topyear)  ! it'll be always true for grasses
   sameSize = sameSizeTree .OR. sameSizeGrass
   lowDensity  = .FALSE. ! c1%nindivs < mindensity
                         ! Weng, 2014-01-27, turned off
   cohorts_can_be_merged = sameSpecies .and. sameLayer .and. sameSize
end function

! ============================================================================
! calculate tree height, DBH, height, and crown area by initial biomass
! The allometry equations are from Ray Dybzinski et al. 2011 and Forrior et al. in review
!         HT = alphaHT * DBH ** (gamma-1)   ! DBH --> Height
!         CA = alphaCA * DBH ** gamma       ! DBH --> Crown Area
!         BM = alphaBM * DBH ** (gamma + 1) ! DBH --> tree biomass
subroutine initialize_cohort_from_biomass(cc,btot)
  type(cohort_type), intent(inout) :: cc
  real,intent(in)    :: btot ! total biomass per individual, kg C

  associate(sp=>spdata(cc%species))
     cc%DBH        = (btot / sp%alphaBM) ** ( 1.0/sp%thetaBM )
     cc%height     = sp%alphaHT * cc%dbh ** sp%thetaHT
     cc%crownarea  = sp%alphaCA * cc%dbh ** sp%thetaCA

     cc%bl_max = sp%LMA   * sp%LAImax        * cc%crownarea
     cc%br_max = sp%phiRL * sp%LAImax/sp%SRA * cc%crownarea
     cc%NSNmax = 0.2 * cc%crownarea ! 5.0*(cc%bl_max/(sp%CNleaf0*sp%leafLS)+cc%br_max/sp%CNroot0)
     cc%nsc    = 2.0 * (cc%bl_max + cc%br_max)
!    N pools
     cc%NSN    = 5.0*(cc%bl_max/sp%CNleaf0 + cc%br_max/sp%CNroot0)
     cc%leafN  = cc%bl/sp%CNleaf0
     cc%rootN  = cc%br/sp%CNroot0
     cc%sapwN  = cc%bsw/sp%CNsw0
     cc%woodN  = cc%bHW/sp%CNwood0
  end associate
end subroutine initialize_cohort_from_biomass

! ============================================================================

subroutine init_cohort_allometry(cc)
  type(cohort_type), intent(inout) :: cc
  ! ----- local var -----------
  real    :: btot ! total biomass per individual, kg C

  btot = max(0.0001,cc%bHW+cc%bsw)
  associate(sp=>spdata(cc%species))
     cc%DBH        = (btot / sp%alphaBM) ** ( 1.0/sp%thetaBM )
!    cc%treeBM     = sp%alphaBM * cc%dbh ** sp%thetaBM
     cc%height     = sp%alphaHT * cc%dbh ** sp%thetaHT
     cc%crownarea  = sp%alphaCA * cc%dbh ** sp%thetaCA

     ! calculations of bl_max and br_max are here only for the sake of the
     ! diagnostics, because otherwise those fields are inherited from the
     ! parent cohort and produce spike in the output, even though these spurious
     ! values are not used by the model
     cc%bl_max = sp%LMA   * sp%LAImax        * cc%crownarea
     cc%br_max = sp%phiRL * sp%LAImax/sp%SRA * cc%crownarea
     cc%NSNmax = 0.2 * cc%crownarea ! 5.0*(cc%bl_max/sp%CNleaf0 + cc%br_max/sp%CNroot0)
  end associate
end subroutine

! =============================================================================
! Added by Weng 2015-02-29
subroutine vegn_annualLAImax_update(vegn)
! used for updating LAImax according to mineral N in soil
! Potential problems:
!   1. All species LAImax are updated
!   2. For evergreen, LAImax can be less than current LAI.
!  Weng, 2017-08-02
  type(tile_type), intent(inout) :: vegn

  ! ---- local vars
  type(cohort_type), pointer :: cc
  real   :: LAImin, LAIfixedN, LAImineralN
  real   :: LAI_Nitrogen
  real   :: fixedN
  logical:: fixedN_based
  integer :: i
  ! Calculating LAI max based on mineral N or mineralN + fixed N
  fixedN_based = .True. ! .False. !
  LAImin       = 0.5

  !fixedN = 0.0
  !do i = 1,vegn%n_cohorts
  !      cc => vegn%cohorts(i)
  !      fixedN = fixedN + cc%fixedN_yr * cc%crownarea * cc%nindivs
  !enddo

 ! Mineral+fixed N-based LAImax
 ! LAI_fixedN = sp%Nfixrate0 * sp%LMA * sp%CNleaf0 * sp%leafLS / sp%LMA
 ! cc%br_max = sp%phiRL*cc%bl_max/(sp%LMA*sp%SRA)
  vegn%previousN = 0.8 * vegn%previousN + 0.2 * vegn%annualN
  do i=0,MSPECIES
      associate (sp => spdata(i) )

      LAIfixedN  = 0.5 * sp%Nfixrate0 * sp%CNleaf0 * sp%leafLS
      LAImineralN = 0.5*vegn%previousN*sp%CNleaf0*sp%leafLS/sp%LMA
      LAI_nitrogen = LAIfixedN + LAImineralN

      spdata(i)%LAImax = MAX(LAImin, MIN(LAI_nitrogen,sp%LAI_light))
      spdata(i)%underLAImax = MIN(sp%LAImax,1.2)
      end associate
  enddo

!  ! update the PFTs in the first layer based on fixed N
!  if(fixedN_based)then ! based on "cc%fixedN_yr + vegn%previousN"
!!    Reset sp%LAImax
!     do i = 1,vegn%n_cohorts
!        cc => vegn%cohorts(i)
!        associate (sp => spdata(cc%species) )
!        sp%LAImax    = 0.0  ! max(sp%LAImax,ccLAImax)
!        sp%layerfrac = 0.0
!        sp%n_cc      = 0
!        end associate
!     enddo
!!   Sum ccLAImax in the first layer
!     do i = 1,vegn%n_cohorts
!        cc => vegn%cohorts(i)
!        associate ( sp => spdata(cc%species) )
!        if(sp%LAImax < LAImin)then
!           LAI_nitrogen = 0.5*(vegn%previousN+cc%fixedN_yr)*sp%CNleaf0*sp%leafLS/sp%LMA
!           if(sp%Nfixrate0 > 0.0)
!           sp%LAImax    = MAX(LAImin, MIN(LAI_nitrogen,sp%LAI_light))
!        endif
!        end associate
!     enddo
!  endif
end subroutine vegn_annualLAImax_update


! ============================================================================
function leaf_area_from_biomass(bl,species,layer,firstlayer) result (area)
  real :: area ! returned value
  real,    intent(in) :: bl      ! biomass of leaves, kg C/individual
  integer, intent(in) :: species ! species
  integer, intent(in) :: layer, firstlayer

! modified by Weng (2014-01-09), 07-18-2017
  area = bl/spdata(species)%LMA
  !if(layer > 1.AND. firstlayer == 0)then
  !   area = bl/(0.5*spdata(species)%LMA) ! half thickness for leaves in understory
  !else
  !   area = bl/spdata(species)%LMA
  !endif
end function
! ============================================================================

! Weng, 2016-11-28
subroutine vegn_annual_diagnostics_zero(vegn)
! for annual update
  type(tile_type), intent(inout) :: vegn
  !-------local var
  type(cohort_type),pointer :: cc
  integer :: i

  vegn%N_P2S_yr  = 0.
  vegn%annualN   = 0.
  vegn%Nloss_yr  = 0.
  vegn%accu_Nup  = 0.0
  vegn%fixedN_yr = 0.
  vegn%annualGPP = 0.0
  vegn%annualNPP = 0.0
  vegn%annualResp = 0.0
  vegn%annualRh   = 0.0

  vegn%LAI     = 0.0

  vegn%NSC     = 0.0
  vegn%SeedC   = 0.0
  vegn%leafC   = 0.0
  vegn%rootC   = 0.0
  vegn%SapwoodC= 0.0
  vegn%WoodC   = 0.0

  vegn%NSN     = 0.0
  vegn%SeedN   = 0.0
  vegn%leafN   = 0.0
  vegn%rootN   = 0.0
  vegn%SapwoodN= 0.0
  vegn%WoodN   = 0.0

  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     cc%annualGPP = 0.0
     cc%annualNPP = 0.0
     cc%annualResp= 0.0
     cc%N_up_yr   = 0.0
     cc%fixedN_yr = 0.0
     cc%NPPleaf   = 0.0
     cc%NPProot   = 0.0
     cc%NPPwood   = 0.0
     cc%DBH_ys    = cc%DBH
  enddo
end subroutine vegn_annual_diagnostics_zero

!============= Parameter and Model initializations =====================
subroutine initialize_PFT_data(namelistfile)
   character(len=50),intent(in) :: namelistfile
  ! ---- local vars
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  integer :: i
  integer :: nml_unit

!  Read parameters from the parameter file (namelist)
  if(read_from_parameter_file)then
      nml_unit = 999
      open(nml_unit, file=namelistfile, form='formatted', action='read', status='old')
      read (nml_unit, nml=vegn_parameters_nml, iostat=io, end=10)
10    close (nml_unit)
   endif
      write(*,nml=vegn_parameters_nml)
  ! initialize vegetation data structure
  spdata%pt         = pt
  spdata%phenotype  = phenotype
  spdata%Vmax       = Vmax
  spdata%Vannual    = Vannual
  spdata%m_cond     = m_cond
  spdata%alpha_phot = alpha_phot
  spdata%gamma_L  = gamma_L
  spdata%gamma_LN = gamma_LN
  spdata%gamma_SW = gamma_SW
  spdata%gamma_FR = gamma_FR

  spdata%rho_FR    = rho_FR
  spdata%root_r    = root_r
  spdata%root_perm = root_perm
!  spdata%rho_N_up0 = rho_N_up0
!  spdata%N_roots0  = N_roots0

  spdata%leaf_size = leaf_size
  spdata%tc_crit   = tc_crit
  spdata%gdd_crit  = gdd_crit

! Plant traits
  spdata%LMA            = LMA      ! leaf mass per unit area, kg C/m2
  spdata%LNbase         = LNbase   ! Basal leaf nitrogen per unit area, kg N/m2
  spdata%CNleafsupport  = CNleafsupport
  spdata%lifeform     = lifeform
  spdata%alphaHT      = alphaHT
  spdata%thetaHT      = thetaHT
  spdata%alphaCA      = alphaCA
  spdata%thetaCA      = thetaCA
  spdata%alphaBM      = alphaBM
  spdata%thetaBM      = thetaBM

  spdata%maturalage   = maturalage
  spdata%v_seed       = v_seed
  spdata%seedlingsize = seedlingsize
  spdata%prob_g       = prob_g
  spdata%prob_e       = prob_e
  spdata%mortrate_d_c = mortrate_d_c
  spdata%mortrate_d_u = mortrate_d_u
  spdata%rho_wood     = rho_wood
  spdata%taperfactor  = taperfactor
  spdata%LAImax       = LAImax
  spdata%underLAImax  = LAImax
  spdata%LAI_light    = LAI_light
  spdata%tauNSC       = tauNSC
  spdata%phiRL        = phiRL
  spdata%phiCSA       = phiCSA
  ! root urnover rate
  spdata%alpha_FR = alpha_FR


!! Nitrogen Weng 2012-10-24
! spdata%CNleaf0 = CNleaf0
  spdata%CNsw0   = CNsw0
  spdata%CNwood0 = CNwood0
  spdata%CNroot0 = CNroot0
  spdata%CNseed0 = CNseed0
  spdata%NfixRate0 = NfixRate0
  spdata%NfixCost0 = NfixCost0

  spdata%internal_gap_frac = internal_gap_frac
  do i = 0, MSPECIES
     call init_derived_species_data(spdata(i))
  enddo
  end subroutine initialize_pft_data
!------------------------------------------
 subroutine init_derived_species_data(sp)
   type(spec_data_type), intent(inout) :: sp
   ! ---- local vars ------
   integer :: j
   real :: specific_leaf_area  ! m2/kgC
   real :: leaf_life_span      ! months

   sp%SRA = 2.0/(sp%root_r*sp%rho_FR)

   ! calculate alphaBM parameter of allometry. note that rho_wood was re-introduced for this calculation
   sp%alphaBM    = sp%rho_wood * sp%taperfactor * PI/4. * sp%alphaHT ! 5200

!  Vmax as a function of LNbase
   sp%Vmax = 0.025 * sp%LNbase ! Vmax/LNbase= 25E-6/0.8E-3 = 0.03125
!  CN0 of leaves
   sp%LNA     = sp%LNbase +  sp%LMA/sp%CNleafsupport
   sp%CNleaf0 = sp%LMA/sp%LNA
!  Leaf life span as a function of LMA
   sp%leafLS = MAX(c_LLS * sp%LMA,1.0)

!  Leaf turnover rate, (leaf longevity as a function of LMA)
   sp%alpha_L = 1.0/sp%leafLS * sp%phenotype

 end subroutine init_derived_species_data

!========================================================
subroutine initialize_vegn_tile(vegn,nCohorts,namelistfile)
   type(tile_type),intent(inout),pointer :: vegn
   integer,intent(in) :: nCohorts
   character(len=50),intent(in) :: namelistfile
!--------local vars -------

   type(cohort_type),dimension(:), pointer :: cc
   type(cohort_type),pointer :: cp
   integer,parameter :: rand_seed = 86456
   real    :: r
   real    :: btotal
   integer :: i, istat
   integer :: io           ! i/o status for the namelist
   integer :: ierr         ! error code, returned by i/o routines
   integer :: nml_unit

!  Read parameters from the parameter file (namelist)
   if(read_from_parameter_file)then
      ! --- Generate cohorts according to "initial_state_nml" ---
      nml_unit = 999
      open(nml_unit, file=namelistfile, form='formatted', action='read', status='old')
      read (nml_unit, nml=initial_state_nml, iostat=io, end=20)
20    close (nml_unit)
      write(*,nml=initial_state_nml)
      ! Initialize plant cohorts
      allocate(cc(1:init_n_cohorts), STAT = istat)
      vegn%cohorts => cc
      vegn%n_cohorts = init_n_cohorts
      cc => null()

      do i=1,init_n_cohorts
         cp => vegn%cohorts(i)
         cp%status  = LEAF_OFF ! ON=1, OFF=0 ! ON
         cp%layer   = 1
         cp%species = init_cohort_species(i)
         cp%ccID =  i
         cp%nsc     = init_cohort_nsc(i)
         cp%nindivs = init_cohort_nindivs(i) ! trees/m2
         cp%bsw     = init_cohort_bsw(i)
         cp%bHW   = init_cohort_bHW(i)
         btotal     = cp%bsw + cp%bHW  ! kgC /tree
         call initialize_cohort_from_biomass(cp,btotal)
      enddo
      MaxCohortID = cp%ccID
      ! Sorting these cohorts
      call relayer_cohorts(vegn)
      ! Initial Soil pools and environmental conditions
      vegn%metabolicL   = init_fast_soil_C ! kgC m-2
      vegn%structuralL  = init_slow_soil_C ! slow soil carbon pool, (kg C/m2)
      vegn%metabolicN   = vegn%metabolicL/CN0metabolicL  ! fast soil nitrogen pool, (kg N/m2)
      vegn%structuralN  = vegn%structuralL/CN0structuralL  ! slow soil nitrogen pool, (kg N/m2)
      vegn%N_input      = N_input  ! kgN m-2 yr-1, N input to soil
      vegn%mineralN     = init_Nmineral  ! Mineral nitrogen pool, (kg N/m2)
      vegn%previousN    = vegn%mineralN

   else
     ! ------- Generate cohorts randomly --------
     ! Initialize plant cohorts
      allocate(cc(1:nCohorts), STAT = istat)
      vegn%cohorts => cc
      vegn%n_cohorts = nCohorts
      cc => null()
      r = rand(rand_seed)
      do i=1,nCohorts
         cp => vegn%cohorts(i)
         cp%status  = LEAF_OFF ! ON=1, OFF=0 ! ON
         cp%layer   = 1
         cp%species = INT(rand()*5)+1
         cp%nindivs = rand()/10. ! trees/m2
         btotal     = rand()*100.0  ! kgC /tree
         call initialize_cohort_from_biomass(cp,btotal)
      enddo
      ! Sorting these cohorts
      call relayer_cohorts(vegn)
      ! ID each cohort
      do i=1,nCohorts
         cp => vegn%cohorts(i)
         cp%ccID = MaxCohortID + i
      enddo
      MaxCohortID = cp%ccID
      ! Initial Soil pools and environmental conditions
      vegn%metabolicL  = 0.2 ! kgC m-2
      vegn%structuralL = 7.0 ! slow soil carbon pool, (kg C/m2)
      vegn%metabolicN  = vegn%metabolicL/CN0metabolicL  ! fast soil nitrogen pool, (kg N/m2)
      vegn%structuralN = vegn%structuralL/CN0structuralL  ! slow soil nitrogen pool, (kg N/m2)
      vegn%N_input     = N_input  ! kgN m-2 yr-1, N input to soil
      vegn%mineralN    = 0.005  ! Mineral nitrogen pool, (kg N/m2)
      vegn%previousN   = vegn%mineralN

   endif  ! initialization: random or pre-described
end subroutine initialize_vegn_tile
!=========================================================================

end module esdvm
