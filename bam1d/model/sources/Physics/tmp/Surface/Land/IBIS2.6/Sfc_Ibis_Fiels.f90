MODULE Sfc_Ibis_Fiels
  USE Parallelism, ONLY:   &
       MsgOne
  USE Constants, ONLY :     &
       r4,i4,r8,i8,tice

  USE Options, ONLY: &
       reducedGrid,initlz,monl, &
       ifalb  , ifsst , &
       ifslm  , ifsnw , &
       ifozone, sstlag, &
       intsst , fint  , &
       iglsm_w, mxiter,nftgz0,nfzol,&
       nfsoiltp,nfvegtp,nfslmtp,nfsibi,isimp,intndvi,&
       nfprt  , nfctrl, nfsibd, nfalb,filta,ifndvi ,ifslmSib2,iftracer,&
       epsflt,istrt,Model1D,schemes,fNameTg3zrl,fNameRouLen,omlmodel,oml_hml0

  USE InputOutput, ONLY: &
       getsbc

  USE IOLowLevel, ONLY: &
       ReadGetNFTGZ

  USE FieldsPhysics, ONLY: &
       gtsea,soilm,o3mix,tg1,tg2,tg3,ssib,wsib3d,gl0,Mmlen,tg0,tgm,&
       tseam,z0,zorl,AlbVisDiff,sheleg,rVisDiff,gndvi,imask,MskAnt,tc0,&
       ppli,ppci,capac0,capacm,td0,w0,wm,tdm,tcm,qsfc0,tsfc0,qsfcm,tsfcm,tkemyj,tracermix,&
       HML,HUML,HVML,TSK,z0sea


  USE Utils, ONLY: &
       IJtoIBJB, &
       LinearIJtoIBJB, &
       NearestIJtoIBJB,&
       FreqBoxIJtoIBJB, &
       lati ,&
       vfirec

  USE Parallelism, ONLY: &
       MsgOne, FatalError

  IMPLICIT NONE
  PRIVATE
  !INTEGER, PARAMETER :: r4 = SELECTED_REAL_KIND(6)  ! Kind for 32-bits Real Numbers
  !INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)   ! Kind for 32-bits Integer Numbers
  !INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers 
  !INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers
  LOGICAL, PUBLIC :: doalb  =.TRUE.! true if surface albedo calculation time step

  ! -------------------------------
  ! state description configuration
  ! -------------------------------
  !
  !
  ! --------------------------
  ! typical ibis configuration
  ! --------------------------
  REAL(KIND=r8), PUBLIC    :: dtime ! model timestep (seconds)
  INTEGER, PUBLIC :: ibMax ! longitude dimension of domain
  INTEGER, PUBLIC :: jbMax ! latitude dimension of domain
  INTEGER, PUBLIC :: iMax  ! longitude dimension of domain
  INTEGER, PUBLIC :: jMax  ! latitude dimension of domain
  INTEGER, PUBLIC :: kMax  ! latitude dimension of domain
  REAL(KIND=r8), PUBLIC    :: xres  ! longitude resolution (in degrees)
  REAL(KIND=r8), PUBLIC    :: yres  ! latitude resolution (in degrees)
  INTEGER, PUBLIC :: idateprev (4)
  INTEGER, PUBLIC :: iyear0

  ! --------------
  ! some constants
  ! --------------

  REAL(KIND=r8), PUBLIC   , PARAMETER :: pi   = 3.1415927_r8   ! you know, that constant thingy

  !
  ! Arguments (input)     
  !
  INTEGER, PUBLIC, PARAMETER :: isimfire=0
  INTEGER, PUBLIC, PARAMETER :: isimco2 =1
  INTEGER, PUBLIC, PARAMETER :: spinmax =14
  INTEGER, PUBLIC, PARAMETER :: isimveg =0         ! 0 = static veg, 
  ! 1 = dynamic veg, 
  ! 2 = dynamic veg with cold start

  INTEGER, PUBLIC            :: irestart=0       ! 0 = initial run, 1 = restart run
  REAL(KIND=r8)              :: yrl22=365.2500_r8

  ! -------------------------------
  ! state description configuration
  ! -------------------------------
  !
  INTEGER, PUBLIC, PARAMETER :: nsoilay=6        ! number of soil layers
  INTEGER, PUBLIC, PARAMETER :: nsnolay=3        ! number of snow layers
  INTEGER, PUBLIC, PARAMETER :: nband=2          ! number of solar radiation wavebands : vis, nir

  INTEGER, PUBLIC, PARAMETER :: npft=12          ! number of plant functional types



  ! ----------------------------------------
  ! Soil texture-related parameters for ibis
  ! ----------------------------------------

  INTEGER, PUBLIC, PARAMETER :: ndat = 11      ! number of soil types, excludes organics for now.
  INTEGER, PUBLIC, PARAMETER :: npftu= 8       ! number of upper canopy pfts (used in modified 
  ! version of subroutine DYNAVEG) 
  !---------------------------------------------------------------------------------
  !  minimum density of woody biomass required for upper canopy closure (kg C m-2)
  !---------------------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER , PUBLIC :: woodnorm=7.5_r8 !  woodnorm  ! value of woody biomass for upper canopy closure (ie when 
  !wood = woodnorm fu = 1.0) (kg_C m-2)
  ! leaf optical properties from Sellers et al., 1996 and Bonan, 1995
  !-----------------------------------------------------------------------------------------
  ! leaf reflectance (rhoveg) and transmittance (tauveg), visible and NIR, for each canopy
  !-----------------------------------------------------------------------------------------
  REAL(KIND=r8),ALLOCATABLE, PUBLIC :: tauwood  (:,:,:)      ! wood biomass turnover time constant (years)

  REAL(KIND=r8), PUBLIC :: rhoveg  (nband,2)   ! reflectance of an average leaf/stem
  REAL(KIND=r8), PUBLIC :: tauveg  (nband,2)   ! transmittance of an average leaf/stem
  !-----------------------------------------------------------------------
  ! linear dimensions for aerodynamic flux parameterization: dleaf, dstem
  !-----------------------------------------------------------------------
  REAL(KIND=r8), PARAMETER   , PUBLIC :: dleaf(1:2)=(/ 0.100_r8 ,0.100_r8/)! typical linear leaf dimension in aerodynamic transfer coefficient (m)
  REAL(KIND=r8), PARAMETER   , PUBLIC :: dstem(1:2)=(/ 0.100_r8 ,0.100_r8/)! typical linear stem dimension in aerodynamic transfer coefficient (m)
  !--------------------------------------------------------------------
  ! normalization constants for canopy drag coefficients (m2 m-2)
  !---------------------------------------------------------------------
  REAL(KIND=r8), PARAMETER   , PUBLIC :: alaiml=8.0_r8! lower canopy leaf & stem maximum area (2 sided) for normalization of
  ! drag coefficient (m2 m-2)
  REAL(KIND=r8), PARAMETER   , PUBLIC :: alaimu=8.0_r8! upper canopy leaf & stem area (2 sided) for normalization of drag
  ! coefficient (m2 m-2)
  !----------------------------------------------------------------------------
  ! empirical coefficients for aerodynamic transfer parameterization (m s-0.5) 
  ! From Pollard & Thompson (1995, eq. A39a)
  !----------------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER   , PUBLIC :: cleaf= 0.01_r8  ! cleaf  : upper canopy leaf-air          ! empirical constant in upper canopy leaf-air aerodynamic transfer
  ! coefficient (m s-0.5) (A39a Pollard & Thompson 95) 
  REAL(KIND=r8), PARAMETER   , PUBLIC :: cstem  = 0.01_r8! cstem  : upper canopy stem-air          ! empirical constant in upper canopy stem-air aerodynamic transfer
  ! coefficient (m s-0.5) (A39a Pollard & Thompson 95) 
  REAL(KIND=r8), PARAMETER   , PUBLIC :: cgrass = 0.01_r8! cgrass : lower canopy-air                  ! empirical constant in lower canopy-air aerodynamic transfer
  ! coefficient (m s-0.5) (A39a Pollard & Thompson 95)
  !----------------------------------------------------------------------------
  ! heat capacities of leaves and stems  (J kg-1 C-1 m-2)
  ! derived from specific heat of liquid water (ch2o = 4.218 J g-1)
  !----------------------------------------------------------------------------
  REAL(KIND=r8), PARAMETER   , PUBLIC :: chs = 0.125e+05_r8!2.109e+05_r8! heat capacity of upper canopy stems per unit stem area (J kg-1 m-2)  
  REAL(KIND=r8), PARAMETER   , PUBLIC :: chu = 2.109e+03_r8!8.436e+03_r8! heat capacity of upper canopy leaves per unit leaf area (J kg-1 m-2) 
  REAL(KIND=r8), PARAMETER   , PUBLIC :: chl = 8.436e+03_r8!8.436e+03_r8! heat capacity of lower canopy leaves & stems per unit leaf/stem area (J kg-1 m-2)
  !-----------------------------------------------------------------------
  ! intercepted water capacity (mm h2o per unit leaf area == kg m-2)
  !-----------------------------------------------------------------------
  REAL(KIND=r8), PARAMETER   , PUBLIC :: wliqumax = 0.20_r8 ! maximum intercepted water on a unit upper canopy leaf area (kg m-2)
  REAL(KIND=r8), PARAMETER   , PUBLIC :: wliqsmax = 0.40_r8 ! maximum intercepted water on a unit upper canopy stem area (kg m-2)
  REAL(KIND=r8), PARAMETER   , PUBLIC :: wliqlmax = 0.20_r8 ! maximum intercepted water on a unit lower canopy stem & leaf area (kg m-2)
  !-----------------------------------------------------------------------
  ! intercepted snow capacity (mm h2o per unit leaf area == kg m-2)
  !-----------------------------------------------------------------------
  REAL(KIND=r8), PARAMETER   , PUBLIC :: wsnoumax = 2.00_r8! intercepted snow capacity for upper canopy leaves (kg m-2)
  REAL(KIND=r8), PARAMETER   , PUBLIC :: wsnosmax = 4.00_r8! intercepted snow capacity for upper canopy stems (kg m-2)
  REAL(KIND=r8), PARAMETER   , PUBLIC :: wsnolmax = 2.00_r8! intercepted snow capacity for lower canopy leaves & stems (kg m-2

  REAL(KIND=r8)   , PUBLIC :: wliqmin                 ! minimum intercepted water on unit vegetated area (kg m-2)
  REAL(KIND=r8)   , PUBLIC :: wsnomin                 ! minimum intercepted snow on unit vegetated area (kg m-2)
  !------------------------------------------------------------
  ! decay time for intercepted liquid dripoff (sec)
  !------------------------------------------------------------
  REAL(KIND=r8), PARAMETER   , PUBLIC :: tdripu = 7200.0_r8! decay time for dripoff of liquid intercepted by upper canopy leaves (sec)
  REAL(KIND=r8), PARAMETER   , PUBLIC :: tdrips = 7200.0_r8! decay time for dripoff of liquid intercepted by upper canopy stems (sec) 
  REAL(KIND=r8), PARAMETER   , PUBLIC :: tdripl = 7200.0_r8! decay time for dripoff of liquid intercepted by lower canopy leaves & stem (sec)
  !------------------------------------------------------------
  ! decay time for snow blowoff (sec)
  !------------------------------------------------------------
  REAL(KIND=r8), PARAMETER   , PUBLIC :: tblowu = 43200.0_r8! decay time for blowoff of snow intercepted by upper canopy leaves (sec)
  REAL(KIND=r8), PARAMETER   , PUBLIC :: tblows = 43200.0_r8! decay time for blowoff of snow intercepted by upper canopy stems (sec)
  REAL(KIND=r8), PARAMETER   , PUBLIC :: tblowl = 43200.0_r8! decay time for blowoff of snow intercepted by lower canopy leaves & stems (sec)

  !
  !      INCLUDE 'compft.h'
  !
  !--------------------------------------------------
  ! Other properties of vegetation -- for 12 PFTs
  !--------------------------------------------------
  ! vmax_pft : max Rubisco activity at 15 C, at top of canopy (mol[CO2] m-2 s-1) 
  ! specla   : specific leaf area (m2 kg-1)
  ! tauleaf  : foliar biomass turnover time constant (years)
  ! tauroot  : root biomass turnover time constant (years)
  ! tauwood0 : wood biomass turnover time constant (years)
  ! aleaf    : foliar allocation coefficient (fraction)
  ! aroot    : root allocation coefficient (fraction)
  ! awood    : wood allocation coefficient (fraction, = 1 - aleaf - aroot)
  !=================================================
  !--------------------------------------------------
  ! Other properties of vegetation -- for 12 PFTs
  !--------------------------------------------------
  ! vmax_pft : max Rubisco activity at 15 C, at top of canopy (mol[CO2] m-2 s-1) 
  ! specla   : specific leaf area (m2 kg-1)
  ! tauleaf  : foliar biomass turnover time constant (years)
  ! tauroot  : root biomass turnover time constant (years)
  ! tauwood  : wood biomass turnover time constant (years)
  ! aleaf    : foliar allocation coefficient (fraction)
  ! aroot    : root allocation coefficient (fraction)
  ! awood    : wood allocation coefficient (fraction, = 1 - aleaf - aroot)
  !      dummyvarpk(1:96)=(/&
  !------------------------------------------------------------------------
  !         vmax_pft  specla   tauleaf  tauroot tauwood  aleaf  aroot  awood    PFT
  !------------------------------------------------------------------------
  !      65.0e-06_r8, 25.0_r8, 1.01_r8, 1.0_r8,  25.0_r8, 0.30_r8, 0.20_r8, 0.50_r8,& !   1 
  !      65.0e-06_r8, 25.0_r8, 1.00_r8, 1.0_r8,  25.0_r8, 0.30_r8, 0.20_r8, 0.50_r8,& !   2 
  !      40.0e-06_r8, 25.0_r8, 1.00_r8, 1.0_r8,  25.0_r8, 0.30_r8, 0.20_r8, 0.50_r8,& !   3 
  !      30.0e-06_r8, 12.5_r8, 2.00_r8, 1.0_r8,  50.0_r8, 0.30_r8, 0.40_r8, 0.30_r8,& !   4 
  !      30.0e-06_r8, 25.0_r8, 1.00_r8, 1.0_r8,  50.0_r8, 0.30_r8, 0.20_r8, 0.50_r8,& !   5 
  !      25.0e-06_r8, 12.5_r8, 2.50_r8, 1.0_r8, 100.0_r8, 0.30_r8, 0.40_r8, 0.30_r8,& !   6 
  !      30.0e-06_r8, 25.0_r8, 1.00_r8, 1.0_r8, 100.0_r8, 0.30_r8, 0.20_r8, 0.50_r8,& !   7 
  !      30.0e-06_r8, 25.0_r8, 1.00_r8, 1.0_r8, 100.0_r8, 0.30_r8, 0.20_r8, 0.50_r8,& !   8 
  !      27.5e-06_r8, 12.5_r8, 1.50_r8, 1.0_r8,   5.0_r8, 0.45_r8, 0.40_r8, 0.15_r8,& !   9 
  !      27.5e-06_r8, 25.0_r8, 1.00_r8, 1.0_r8,   5.0_r8, 0.45_r8, 0.35_r8, 0.20_r8,& !  10 
  !      15.0e-06_r8, 20.0_r8, 1.25_r8, 1.0_r8, 999.0_r8, 0.45_r8, 0.55_r8, 0.00_r8,& !  11 
  !      25.0e-06_r8, 20.0_r8, 1.50_r8, 1.0_r8, 999.0_r8, 0.45_r8, 0.55_r8, 0.00_r8/)  !  12 
  !========================================================================
  REAL(KIND=r8)   , PUBLIC , PARAMETER :: vmax_pft(1:npft) =(/&
       !    vmax_pft        ! PFT
       120.0e-06_r8,&    !   1   65.0e-06_r8   120.0e-06
       65.0e-06_r8,&     !   2   65.0e-06_r8    65.0e-06
       40.0e-06_r8,&     !   3   40.0e-06_r8    40.0e-06
       30.0e-06_r8,&     !   4   30.0e-06_r8    30.0e-06
       30.0e-06_r8,&     !   5   30.0e-06_r8    30.0e-06
       25.0e-06_r8,&     !   6   25.0e-06_r8    25.0e-06
       30.0e-06_r8,&     !   7   30.0e-06_r8    30.0e-06
       30.0e-06_r8,&     !   8   30.0e-06_r8    30.0e-06
       27.5e-06_r8,&     !   9   27.5e-06_r8    27.5e-06
       27.5e-06_r8,&     !  10   27.5e-06_r8    27.5e-06
       15.0e-06_r8,&     !  11   15.0e-06_r8    15.0e-06
       25.0e-06_r8/)     !  12   25.0e-06_r8    25.0e-06
  REAL(KIND=r8)   , PUBLIC  , PARAMETER :: specla  (1:npft) =(/&     ! specific leaf area (m**2/kg) 
       !  specla        ! PFT
       25.0_r8,&         !   1   25.0_r8  25.0        !  1: tropical broadleaf evergreen trees
       25.0_r8,&         !   2   25.0_r8  25.0        !  2: tropical broadleaf drought-deciduous trees
       25.0_r8,&         !   3   25.0_r8  25.0        !  3: warm-temperate broadleaf evergreen trees
       12.5_r8,&         !   4   12.5_r8  12.5        !  4: temperate conifer evergreen trees
       25.0_r8,&         !   5   25.0_r8  25.0        !  5: temperate broadleaf cold-deciduous trees
       12.5_r8,&         !   6   12.5_r8  12.5        !  6: boreal conifer evergreen trees
       25.0_r8,&         !   7   25.0_r8  25.0        !  7: boreal broadleaf cold-deciduous trees
       25.0_r8,&         !   8   25.0_r8  25.0        !  8: boreal conifer cold-deciduous trees
       12.5_r8,&         !   9   12.5_r8  12.5        !  9: evergreen shrubs
       25.0_r8,&         !  10   25.0_r8  25.0        ! 10: cold-deciduous shrubs
       20.0_r8,&         !  11   20.0_r8  20.0        ! 11: warm (c4) grasses
       20.0_r8/)         !  12         20.0_r8  20.0        ! 12: cool (c3) grasses
  REAL(KIND=r8)   , PUBLIC  , PARAMETER :: tauleaf (1:npft)=(/&       ! foliar biomass turnover time constant (years)
       !   tauleaf      ! PFT
       1.01_r8,&         !   1  1.01_r8  1.01
       1.00_r8,&         !   2  1.00_r8  1.00
       1.00_r8,&         !   3  1.00_r8  1.00
       2.00_r8,&         !   4  2.00_r8  2.00
       1.00_r8,&         !   5  1.00_r8  1.00
       2.50_r8,&         !   6  2.50_r8  2.50
       1.00_r8,&         !   7  1.00_r8  1.00
       1.00_r8,&         !   8  1.00_r8  1.00
       1.50_r8,&         !   9  1.50_r8  1.50
       1.00_r8,&         !  10  1.00_r8  1.00
       1.25_r8,&         !  11  1.25_r8  1.25
       1.50_r8/)         !  12        1.50_r8  1.50
  REAL(KIND=r8)   , PUBLIC  , PARAMETER :: tauroot(1:npft)=(/& ! fine root biomass turnover time constant (years)
       !   tauroot      ! PFT        
       1.0_r8,& !        1  1.0_r8
       1.0_r8,& !        2  1.0_r8
       1.0_r8,& !        3  1.0_r8
       1.0_r8,& !        4  1.0_r8
       1.0_r8,& !        5  1.0_r8
       1.0_r8,& !        6  1.0_r8
       1.0_r8,& !        7  1.0_r8
       1.0_r8,& !        8  1.0_r8
       1.0_r8,& !        9  1.0_r8
       1.0_r8,& !  10  1.0_r8
       1.0_r8,& !  11  1.0_r8
       1.0_r8/) !  12  1.0_r8
  REAL(KIND=r8), PUBLIC  , PARAMETER :: tauwood0(1:npft)=(/& ! normal (unstressed) turnover time for wood biomass (years)
       !   tauwood0      ! PFT    
       25.0_r8,& !   1    25.0_r8
       25.0_r8,& !   2    25.0_r8
       25.0_r8,& !   3    25.0_r8
       50.0_r8,& !   4    50.0_r8
       50.0_r8,& !   5    50.0_r8
       100.0_r8,& !   6   100.0_r8
       100.0_r8,& !   7   100.0_r8
       100.0_r8,& !   8   100.0_r8
       5.0_r8,& !   9     5.0_r8
       5.0_r8,& !  10     5.0_r8
       999.0_r8,& !  11   999.0_r8
       999.0_r8/) !  12   999.0_r8
  REAL(KIND=r8)   , PUBLIC  , PARAMETER :: aleaf (1:npft) =(/&      ! ! carbon allocation fraction to leaves
       !   aleaf      ! PFT    
       0.30_r8,& !   1   0.30
       0.30_r8,& !   2   0.30
       0.30_r8,& !   3   0.30
       0.30_r8,& !   4   0.30
       0.30_r8,& !   5   0.30
       0.30_r8,& !   6   0.30
       0.30_r8,& !   7   0.30
       0.30_r8,& !   8   0.30
       0.45_r8,& !   9   0.45
       0.45_r8,& !  10   0.45
       0.45_r8,& !  11   0.45
       0.45_r8/) !  12   0.45
  REAL(KIND=r8)   , PUBLIC  , PARAMETER :: aroot   (1:npft) =(/&        ! carbon allocation fraction to fine roots
       !   aroot      ! PFT    
       0.20_r8,& !   1 0.20
       0.20_r8,& !   2 0.20
       0.20_r8,& !   3 0.20
       0.40_r8,& !   4 0.40
       0.20_r8,& !   5 0.20
       0.40_r8,& !   6 0.40
       0.20_r8,& !   7 0.20
       0.20_r8,& !   8 0.20
       0.40_r8,& !   9 0.40
       0.35_r8,& !  10 0.35
       0.55_r8,& !  11 0.55
       0.55_r8/) !  12 0.55
  REAL(KIND=r8)   , PUBLIC  , PARAMETER :: awood   (1:npft)=(/&         ! carbon allocation fraction to wood 
       !   awood      ! PFT    
       0.50_r8,& !   1 0.50
       0.50_r8,& !   2 0.50
       0.50_r8,& !   3 0.50
       0.30_r8,& !   4 0.30
       0.50_r8,& !   5 0.50
       0.30_r8,& !   6 0.30
       0.50_r8,& !   7 0.50
       0.50_r8,& !   8 0.50
       0.15_r8,& !   9 0.15
       0.20_r8,& !  10 0.20
       0.00_r8,& !  11 0.00
       0.00_r8/) !  12 0.00
  !-!------------------------------------------------------------------------
  !  Fundamental plant physiological parameters, name definitions
  ! ------------------------------------------------------------------------
  !      tau15 = 4500.0_r8    ! tau15 : co2/o2 specificity ratio at 15 degrees C (dimensionless)    
  !      kc15  = 1.5e-04_r8   ! kc15  : co2 kinetic parameter at 15 C (mol/mol)
  !      ko15  = 2.5e-01_r8   ! ko15  : o2 kinetic parameter at 15 C (mol/mol)
  !      cimax = 2000.e-06_r8 ! cimax : maximum value for ci (for model stability)

  REAL(KIND=r8), PARAMETER   , PUBLIC :: tau15   = 4500.0_r8              ! co2/o2 specificity ratio at 15 degrees C (dimensionless)
  REAL(KIND=r8), PARAMETER   , PUBLIC :: kc15    = 1.5e-04_r8              ! co2 kinetic parameter (mol/mol)
  REAL(KIND=r8), PARAMETER   , PUBLIC :: ko15    = 2.5e-01_r8             ! o2 kinetic parameter (mol/mol) 
  REAL(KIND=r8), PARAMETER   , PUBLIC :: cimax   = 2000.0e-06_r8           ! maximum value for ci (needed for model stability)
  !------------------------------------------------------------------------
  ! PFTs (top to bottom)
  !------------------------------------------------------------------------
  !  1: tropical broadleaf evergreen trees
  !  2: tropical broadleaf drought-deciduous trees
  !  3: warm-temperate broadleaf evergreen trees
  !  4: temperate conifer evergreen trees
  !  5: temperate broadleaf cold-deciduous trees
  !  6: boreal conifer evergreen trees
  !  7: boreal broadleaf cold-deciduous trees
  !  8: boreal conifer cold-deciduous trees
  !  9: evergreen shrubs
  ! 10: cold-deciduous shrubs
  ! 11: warm (C4) grasses
  ! 12: cool (C3) grasses
  !========================================================================

  !--------------------------------------------------------------------
  ! C3 and C4 physiology-specific parameters
  !--------------------------------------------------------------------
  REAL(KIND=r8)   ,PUBLIC, PARAMETER :: alpha3 =0.080_r8!0.040_r8 ! alpha3 - C3 intrinsic quantum efficiency (dimensionless)
  REAL(KIND=r8)   ,PUBLIC, PARAMETER :: theta3 =0.950_r8!0.550_r8 ! theta3 - C3 photosynthesis coupling coefficient
  REAL(KIND=r8)   ,PUBLIC, PARAMETER :: beta3  =0.990_r8!0.590_r8 ! beta3  - C3 photosynthesis coupling coefficient
  REAL(KIND=r8)   ,PUBLIC, PARAMETER :: alpha4 =0.050_r8!0.010_r8 ! alpha4 - C4 intrinsic quantum efficiency (dimensionless)
  REAL(KIND=r8)   ,PUBLIC, PARAMETER :: theta4 =0.970_r8!0.570_r8 ! theta4 - C4 photosynthesis coupling coefficient
  REAL(KIND=r8)   ,PUBLIC, PARAMETER :: beta4  =0.800_r8!0.400_r8 ! beta4  - C4 photosynthesis coupling coefficient 
  !--------------------------------------------------------------------
  ! Plant physiological properties - 5 classes
  !--------------------------------------------------------------------
  ! gamma    : leaf respiration coefficients 
  ! coefm    : 'm' coefficients for stomatal conductance relationship
  ! coefb    : 'b' coefficients for stomatal conductance relationship
  ! gsmin    : absolute minimum stomatal conductances
  !-------------------------------------------------
  ! gamma  coefm  coefb    gsmin   Physiol. Class
  !-------------------------------------------------
  REAL(KIND=r8)   ,PUBLIC, PARAMETER :: gammaub = 0.015_r8   ! Broadleaf trees
  REAL(KIND=r8)   ,PUBLIC, PARAMETER :: coefmub = 11.0_r8  !10.0_r8    ! Broadleaf trees 'm' coefficient for stomatal conductance relationship
  REAL(KIND=r8)   ,PUBLIC, PARAMETER :: coefbub = 0.01_r8    ! Broadleaf trees 'b' coefficient for stomatal conductance relationship
  REAL(KIND=r8)   ,PUBLIC, PARAMETER :: gsubmin = 0.00001_r8 ! Broadleaf trees  absolute minimum stomatal conductance

  REAL(KIND=r8)   ,PUBLIC, PARAMETER :: gammauc = 0.015_r8   ! leaf respiration coefficient
  REAL(KIND=r8)   ,PUBLIC, PARAMETER :: coefmuc = 6.0_r8     !'m' coefficient for stomatal conductance relationship 
  REAL(KIND=r8)   ,PUBLIC, PARAMETER :: coefbuc = 0.01_r8    !'b' coefficient for stomatal conductance relationship  
  REAL(KIND=r8)   ,PUBLIC, PARAMETER :: gsucmin = 0.00001_r8 ! absolute minimum stomatal conductance

  REAL(KIND=r8)   ,PUBLIC, PARAMETER :: gammals = 0.015_r8  ! Shrubs ! leaf respiration coefficient
  REAL(KIND=r8)   ,PUBLIC, PARAMETER :: coefmls = 9.0_r8    ! Shrubs ! 'm' coefficient for stomatal conductance relationship 
  REAL(KIND=r8)   ,PUBLIC, PARAMETER :: coefbls = 0.01_r8   ! Shrubs ! 'b' coefficient for stomatal conductance relationship 
  REAL(KIND=r8)   ,PUBLIC, PARAMETER :: gslsmin = 0.00001_r8! Shrubs ! absolute minimum stomatal conductance 

  REAL(KIND=r8)   ,PUBLIC, PARAMETER :: gammal4 = 0.030_r8  ! C4 grasses ! leaf respiration coefficient
  REAL(KIND=r8)   ,PUBLIC, PARAMETER :: coefml4 = 4.0_r8    ! C4 grasses ! 'm' coefficient for stomatal conductance relationship
  REAL(KIND=r8)   ,PUBLIC, PARAMETER :: coefbl4 = 0.04_r8   ! C4 grasses ! 'b' coefficient for stomatal conductance relationship
  REAL(KIND=r8)   ,PUBLIC, PARAMETER :: gsl4min = 0.00001_r8! C4 grasses ! absolute minimum stomatal conductance

  REAL(KIND=r8)   ,PUBLIC, PARAMETER :: gammal3 = 0.015_r8  ! C3 grasses ! leaf respiration coefficient
  REAL(KIND=r8)   ,PUBLIC, PARAMETER :: coefml3 = 9.0_r8    ! C3 grasses ! 'm' coefficient for stomatal conductance relationship 
  REAL(KIND=r8)   ,PUBLIC, PARAMETER :: coefbl3 = 0.01_r8   ! C3 grasses ! 'b' coefficient for stomatal conductance relationship 
  REAL(KIND=r8)   ,PUBLIC, PARAMETER :: gsl3min = 0.00001_r8! C3 grasses ! absolute minimum stomatal conductance
  !-----------------------------------------------------------------
  ! leaf orientation factors (-1 vertical, 0 random, 1 horizontal)
  !-----------------------------------------------------------------
  ! chifuz : upper canopy leaf orientation
  ! chiflz : lower canopy leaf orientation
  !-----------------------------------------
  REAL(KIND=r8)   ,PARAMETER, PUBLIC :: chifuz=   0.0_r8    ! chifuz              ! upper canopy leaf orientation factor
  REAL(KIND=r8)   ,PARAMETER, PUBLIC :: chiflz=  -0.5_r8    ! chiflz              ! lower canopy leaf orientation factor
  !------------------------------------------------------------------------
  ! PFTs (top to bottom)
  !------------------------------------------------------------------------
  !  1: tropical broadleaf evergreen trees
  !  2: tropical broadleaf drought-deciduous trees
  !  3: warm-temperate broadleaf evergreen trees
  !  4: temperate conifer evergreen trees
  !  5: temperate broadleaf cold-deciduous trees
  !  6: boreal conifer evergreen trees
  !  7: boreal broadleaf cold-deciduous trees
  !  8: boreal conifer cold-deciduous trees
  !  9: evergreen shrubs
  ! 10: cold-deciduous shrubs
  ! 11: warm (c4) grasses
  ! 12: cool (c3) grasses

  !--------------------------------------------------------------------------
  ! PFT climatic constraint definitions (left to right)
  !--------------------------------------------------------------------------
  ! TminL  : absolute minimum temperature (lower limit, C) 
  ! TminU  : absolute minimum temperature (upper limit, C) 
  ! Twarm  : temperature of the warmest month (mean??, C) [C4 only]
  ! GDD    : min growing degree days above 5 C threshold [upper canopy], or
  !          min growing degree days above 0 C threshold [lower canopy]

  ! DTP 2001/06/07: Changed this after studying code in climate.f. 
  !      Values of 9999 indicate this constraint is not used to 
  !      determine existence of the PFT.  
  REAL(KIND=r8), PARAMETER   , PUBLIC :: TminL(1:npft)=(/&! Absolute minimum temperature -- lower limit (upper canopy PFTs)
       !        TminL     PFT
       0.0_r8,&  !   1
       0.0_r8,&  !   2
       -10.0_r8,&  !   3
       -45.0_r8,&  !   4
       -45.0_r8,&  !   5
       -57.5_r8,&  !   6
       -57.5_r8,&  !   7
       9999.0_r8,&  !   8
       9999.0_r8,&  !   9
       9999.0_r8,&  !  10
       9999.0_r8,&  !  11
       9999.0_r8/)  !  12

  REAL(KIND=r8), PARAMETER   , PUBLIC :: TminU(1:npft)=(/&! Absolute minimum temperature -- upper limit (upper canopy PFTs)
       !        TminU     PFT
       9999.0_r8,&  !   1
       9999.0_r8,&  !   2
       0.0_r8,&  !   3
       0.0_r8,&  !   4
       0.0_r8,&  !   5
       -45.0_r8,&  !   6
       -45.0_r8,&  !   7
       -45.0_r8,&  !   8
       9999.0_r8,&  !   9
       9999.0_r8,&  !  10
       9999.0_r8,&  !  11
       9999.0_r8/)  !  12
  REAL(KIND=r8), PARAMETER   , PUBLIC :: Twarm(1:npft)=(/&! Temperature of warmest month (lower canopy PFTs)
       9999.0_r8,&  !   1
       9999.0_r8,&  !   2
       9999.0_r8,&  !   3
       9999.0_r8,&  !   4
       9999.0_r8,&  !   5
       9999.0_r8,&  !   6
       9999.0_r8,&  !   7
       9999.0_r8,&  !   8
       9999.0_r8,&  !   9
       9999.0_r8,&  !  10
       22.0_r8,&  !  11
       9999.0_r8/)  !  12
  REAL(KIND=r8), PARAMETER   , PUBLIC :: GDD(1:npft)=(/&! minimum GDD needed (base 5 C for upper canopy PFTs, 
                                ! base 0 C for lower canopy PFTs)
       !           GDD    PFT
       9999.0_r8,&  !   1
       9999.0_r8,&  !   2
       9999.0_r8,&  !   3
       1200.0_r8,&  !   4
       1200.0_r8,&  !   5
       350.0_r8,&  !   6
       350.0_r8,&  !   7
       350.0_r8,&  !   8
       100.0_r8,&  !   9
       100.0_r8,&  !  10
       100.0_r8,&  !  11
       100.0_r8/)  !  12


  REAL(KIND=r8)   , PUBLIC :: plai_init(4,15)     ! initial total LAI for each vegtype (used in iniveg)

  !------------------------------------------------------------------------
  ! Other miscellaneous variables needed for initializing plant LAI.
  !------------------------------------------------------------------------
  ! plaiupper    : Potental LAI of upper canopy (uniform initial vegetation) 
  ! plailower    : Potental LAI of lower canopy (uniform initial vegetation) 
  ! xminlai      : Minimum LAI for each existing PFT
  ! sapfrac_init : Initial value of sapwood fraction used for all woody PFTs
  !------------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER   , PUBLIC :: plaiupper    =0.5_r8          ! Potental LAI of upper canopy (uniform initial vegetation)
  REAL(KIND=r8), PARAMETER   , PUBLIC :: plailower    =0.5_r8          ! Potental LAI of lower canopy (uniform initial vegetation)
  REAL(KIND=r8), PARAMETER   , PUBLIC :: xminlai      =0.01_r8         ! Minimum LAI for each existing PFT
  REAL(KIND=r8), PARAMETER   , PUBLIC :: sapfrac_init =0.1_r8      ! Initial value of sapwood fraction used for all woody PFTs
  ! ************************************************************************
  ! define rooting profiles
  ! ************************************************************************
  !
  ! define rooting profiles based upon data published in:
  !
  ! Jackson et al., 1996:  A global analysis of root distributions
  ! for terrestrial biomes, Oecologia, 108, 389-411.
  !
  ! and
  !
  ! Jackson et al., 1997:  A global budget for fine root biomass, 
  ! surface area, and nutrient contents, Proceedings of the National
  ! Academy of Sciences, 94, 7362-7366.
  !
  ! rooting profiles are defined by the "beta" parameter
  !
  ! beta1 is assigned to the lower vegetation layer (grasses and shrubs)
  ! beta2 is assigned to the upper vegetation layer (trees)
  !
  ! according to Jackson et al. (1996, 1997), the values of beta
  ! typically fall in the following range
  !
  ! note that the 1997 paper specifically discusses the distribution
  ! of *fine roots* (instead of total root biomass), which may be more
  ! important for water and nutrient uptake
  !
  ! --------------                 ------------   ------------
  ! forest systems                 beta2 (1996)   beta2 (1997)
  ! --------------                 ------------   ------------
  ! tropical evergreen forest:        0.962          0.972
  ! tropical deciduous forest:        0.961          0.982
  ! temperate conifer forest:         0.976          0.980
  ! temperate broadleaf forest:       0.966          0.967
  ! all tropical/temperate forest:    0.970  
  ! boreal forest:                    0.943          0.943
  ! all trees:                                       0.976
  !
  ! -------------------------      ------------   ------------
  ! grassland / shrub systems      beta1 (1996)   beta1 (1997)
  ! -------------------------      ------------   ------------
  ! tropical grassland / savanna:     0.972          0.972
  ! temperate grassland:              0.943          0.943
  ! all grasses:                      0.952          0.952
  ! schlerophyllous shrubs:           0.964          0.950
  ! all shrubs:                       0.978          0.975
  ! crops:                            0.961
  ! desert:                           0.975          0.970
  ! tundra:                           0.914
  !
  ! --------------                 ------------
  ! all ecosystems                 beta  (1996)
  ! --------------                 ------------
  ! all ecosystems:                   0.966
  !
  ! for global simulations, we typically assign the following
  ! values to the beta parameters
  !
  ! beta1 = 0.950, which is typical for tropical/temperate grasslands
  ! beta2 = 0.970, which is typical for tropical/temperate forests
  !
  ! however, these values could be (and should be) further refined
  ! when using the model for specific regions
  ! 
  ! beta1: for lower layer herbaceous plants
  ! beta2: for upper layer trees      
  REAL(KIND=r8), PARAMETER   , PUBLIC :: beta1 = 0.970_r8! parameter for Jackson rooting profile, lower canopy
  REAL(KIND=r8), PARAMETER   , PUBLIC :: beta2 = 0.985_r8! parameter for Jackson rooting profile, upper canopy
  !
  !========================================================================
  ! params.soi : soil parameters....
  !========================================================================
  !
  !   6    ! nsoilay : number of soil layers (actually a constant in comsoi.h)
  !        !           could be read in here, if a new constant, maxsoilay,
  !        !           was introduced to dimension the soil layer arrays  
  !------------------------------------------------------
  ! Soil layer thicknesses (m) 
  ! N.B. Number of layers must equal nsoilay!!!
  !------------------------------------------------------

  REAL(KIND=r8), PARAMETER   , PUBLIC :: hsoi(1:nsoilay+1) =(/&    ! soil layer thickness (m)           
       0.10_r8,&  ! hsoi(1)
       0.15_r8,&  ! hsoi(2)
       0.25_r8,&  ! hsoi(3)
       0.50_r8,&  ! hsoi(4)
       5.00_r8,&  ! hsoi(5)
       10.00_r8,&  ! hsoi(6)
       16.00_r8/)  ! hsoi(6)
  !------------------------------------------------------
  ! Other miscellaneous soil parameters
  !------------------------------------------------------

  REAL(KIND=r8), PARAMETER   , PUBLIC :: bperm        = 1.00_r8 ! lower b.c. for soil profile drainage 
  ! (0.0 = impermeable; 1.0 = fully permeable)   
  REAL(KIND=r8), PARAMETER   , PUBLIC :: wpudmax = 4.5_r8! normalization constant for puddles (kg m-2)  
  REAL(KIND=r8), PARAMETER   , PUBLIC :: zwpmax  = 0.5_r8! assumed maximum fraction of soil surface 
  ! covered by puddles (dimensionless)
  !------------------------------------------------------
  ! Soil properties data from Rawls et al. (1992)
  ! Organic properties data compiled by Mustapha El Maayar (2000)
  ! Organic FC and WP taken from Nijssen et al., 1997 (JGR; table 3 OBS-top)

  !------------------------------------------------------
  ! Variable column header definitions
  !------------------------------------------------------
  ! Sand     : sand fraction
  ! Silt     : silt fraction
  ! Clay     : clay fraction
  ! Porosity : porosity (volume fraction)
  ! FC       : field capacity (volume fraction)
  ! WP       : wilting point (volume fraction)
  ! bexp     : Campbell's 'b' exponent
  ! AEP      : air entry potential (m-H20)
  ! SHC      : saturated hydraulic conductivity (m s-1)

  REAL(KIND=r8)   , PUBLIC :: texdat    (3,ndat)  ! sand/silt/clay fractions
  REAL(KIND=r8)   , PUBLIC :: porosdat  (ndat)    ! porosity volume fraction
  REAL(KIND=r8)   , PUBLIC :: sfielddat (ndat)    ! field capacity volume fraction
  REAL(KIND=r8)   , PUBLIC :: swiltdat  (ndat)    ! wilting point volume fraction
  REAL(KIND=r8)   , PUBLIC :: bexdat    (ndat)    ! Campbell moisture-release b exponent
  REAL(KIND=r8)   , PUBLIC :: suctiondat(ndat)    ! Air entry potential (m-H20)
  REAL(KIND=r8)   , PUBLIC :: hydrauldat(ndat)    ! saturated hydraulic conductivity (m s-1)
  !----------------------------------------------------------------------
  ! Decomposition pool/transformation parameters (see also Kucharik 
  ! et al. 2000)

  !----------------------------------------------------------------------
  ! lig_frac: split of lignified litter material between protected and
  ! non-protected slow OM pools
  !----------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER   , PUBLIC :: lig_frac =  0.50_r8 ! split of lignified litter material between protected/non-protected slow OM pools

  !----------------------------------------------------------------------
  ! fbsom: protected biomass as a fraction of total soil organic carbon
  ! from Verberne et al., 1990
  !----------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER   , PUBLIC :: fbsom = 0.017_r8! protected biomass as a fraction of total soil organic C from Verberne et al., 1990
  !----------------------------------------------------------------------
  ! effac: efficiency of microbial biomass reincorporated into biomass
  ! pool (from NCSOIL parameterizations; Molina et al., 1983)
  !----------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER   , PUBLIC :: effac = 0.40_r8 ! efficiency of microbial biomass reincorporated into biomass pool. 
  ! (From NCSOIL parameterizations; Molina et al., 1983)

  !======================================================================
  ! Define C:N ratios of substrate pools and biomass: 
  ! Values for metabolic and structural plant material and for lignin are 
  ! from Parton et al., 1987 and Whitmore and Parry, 1988, indexed as follows:
  !      cnr(1): c:n ratio of microbial biomass
  !      cnr(2): c:n ratio of passive soil carbon
  !      cnr(3): c:n ratio of protected slow soil carbon
  !      cnr(4): c:n ratio of non-protected slow soil C
  !      cnr(5): c:n ratio of resistant litter lignin
  !      cnr(6): c:n ratio of structural plant (leaf and root) litter
  !      cnr(7): c:n ratio of metabolic (plant and root) litter
  !      cnr(8): c:n ratio of woody biomass components

  REAL(KIND=r8), PARAMETER   , PUBLIC :: cnr(1:10)=(/&! C:N ratios of substrate pools and biomass for leaves and roots.
                                ! Values from Parton et al., 1987 and Whitmore and Parry, 1988
       !---------------------------------------------------------------------
       !   cnr(1)  cnr(2)  cnr(3)  cnr(4)  cnr(5)  cnr(6)  cnr(7)  cnr(8)  cnr(9)  cnr(10)
       !---------------------------------------------------------------------
       8.0_r8,   15.0_r8,   10.0_r8,   15.0_r8,  100.0_r8,  150.0_r8,    6.0_r8,   250.0_r8 ,0.0_r8 ,0.0_r8/)

  ! Miscellaneous other C:N factors...
  !      fmax  : maximum fraction allowed in resistant fraction
  !      rconst: rconst is a constant defined as 1200 [Huh?]
  !      cnleaf: average c:n ratio for leaf litterfall 
  !      cnroot: average c:n ratio for root turnover
  !      cnwood: average c:n ratio for woody debris

  REAL(KIND=r8), PARAMETER   , PUBLIC :: fmax        =    0.45_r8        ! maximum fraction allowed in resistant fraction (Verbene 1997)
  REAL(KIND=r8), PARAMETER   , PUBLIC :: rconst= 1200.0_r8                  ! constant defined as 1200 (from Verbene 1997 equations)
  REAL(KIND=r8), PARAMETER   , PUBLIC :: cnleaf=   40.0_r8                  ! average c:n ratio for leaf litterfall
  REAL(KIND=r8), PARAMETER   , PUBLIC :: cnroot=   60.0_r8                  ! average c:n ratio for root turnover
  REAL(KIND=r8), PARAMETER   , PUBLIC :: cnwood=  200.0_r8                  ! average c:n ratio for woody debris 
  !--------------------------------------------------------------------------------------

  ! Specific maximum decay rate or growth constants; rates are per day.
  ! Constants are taken from Parton et al., 1987 and Verberne et al., 1990
  ! and special issue of Geoderma (comparison of 9 organic matter models) in Dec. 1997

  ! Leaching parameterization was changed to agree with field data, and led to 
  ! changes in the values of the constants given below.  

  ! Approximate factors for Verberne et al. model where efficiencies are 100%
  ! for some of the transformations: one problem was that their rate constants were
  ! based on 25C, and our modifying functions are based on 15 C...thus the rate constants
  ! are somewhat smaller compared to the Verberne et al. (1990) model parameters.
  ! Rates are based on a daily decomposition timestep (per day)
  !      klm: dpm leaf litter --> microbial biomass
  !      kls: spm leaf litter --> microbial biomass
  !      kll: rpm leaf litter --> non or protected om
  !      krm: dpm root litter --> microbial biomass
  !      krs: spm root litter --> microbial biomass
  !      krl: rpm root litter --> non or protected om 
  !      kwm: dpm woody litter --> microbial biomass
  !      kws: spm woody litter --> microbial biomass
  !      kwl: rpm woody litter --> non or protected om 

  REAL(KIND=r8), PARAMETER   , PUBLIC :: klm = 0.150_r8  ! leaf metabolic litter 
  REAL(KIND=r8), PARAMETER   , PUBLIC :: kls = 0.010_r8  ! leaf structural litter
  REAL(KIND=r8), PARAMETER   , PUBLIC :: kll = 0.010_r8  ! leaf lignin
  REAL(KIND=r8), PARAMETER   , PUBLIC :: krm = 0.100_r8  ! root metabolic litter
  REAL(KIND=r8), PARAMETER   , PUBLIC :: krs = 0.005_r8  ! root structural litter
  REAL(KIND=r8), PARAMETER   , PUBLIC :: krl = 0.005_r8  ! root lignin
  REAL(KIND=r8), PARAMETER   , PUBLIC :: kwm = 0.001_r8  ! woody metabolic litter
  REAL(KIND=r8), PARAMETER   , PUBLIC :: kws = 0.001_r8  ! woody structural litter
  REAL(KIND=r8), PARAMETER   , PUBLIC :: kwl = 0.001_r8  ! wood  lignin
  !----------------------------------------------------------------------
  !
  !      kbn: biomass --> non protected organic matter 
  !      kbp: biomass --> protected organic matter
  !      knb: non-protected om --> biomass
  !      kns: non-protected om --> stabilized om
  !      kpb: protected om --> biomass
  !      kps: protected om --> stabilized om
  !      ksb: stabilized om --> biomass


  REAL(KIND=r8), PARAMETER   , PUBLIC :: kbn = 0.045_r8         ! microbial biomass --> nonprotected om 
  REAL(KIND=r8), PARAMETER   , PUBLIC :: kbp = 0.005_r8         ! microbial biomass --> protected om
  REAL(KIND=r8), PARAMETER   , PUBLIC :: knb = 0.001_r8         ! nonprotected om   --> biomass
  REAL(KIND=r8), PARAMETER   , PUBLIC :: kns = 1.0e-06_r8         ! nonprotected om   --> passive c 
  REAL(KIND=r8), PARAMETER   , PUBLIC :: kpb = 0.0001_r8         ! protected om      --> biomass
  REAL(KIND=r8), PARAMETER   , PUBLIC :: kps = 1.0e-06_r8         ! protected om      --> passive c
  REAL(KIND=r8), PARAMETER   , PUBLIC :: ksb = 8.0e-07_r8         ! passive c             --> biomass
  !
  ! Yields (efficiencies) with which microbes gain biomass from C 
  ! source; the rest is driven off as CO2 (microbial respiration). All 
  ! microbial CO2 is assumed to leave the soil profile over the course 
  ! of a year. Values are taken primarily from the models of Verberne 
  ! and from CENTURY.
  !----------------------------------------------------------------------
  !      ylm: efficiency for metabolic plant material - leaf matter
  !      yrm: efficiency for metabolic plant material - root matter
  !      ywm: efficiency for metabolic plant material - woody matter
  !      yls: efficiency for structural plant material - leaf matter
  !      yrs: efficiency for structural plant material - root matter
  !      yws: efficiency for structural plant material - woody matter
  !      yll: plant material resistant fraction - leaf matter
  !      yrl: plant material resistant fraction - root matter
  !      ywl: plant material resistant fraction - woody matter

  REAL(KIND=r8), PARAMETER   , PUBLIC :: ylm = 0.40_r8 ! leaf metabolic litter decomposition 
  REAL(KIND=r8), PARAMETER   , PUBLIC :: yrm = 0.40_r8 ! root metabolic litter decomposition
  REAL(KIND=r8), PARAMETER   , PUBLIC :: ywm = 0.40_r8 ! woody metabolic litter decomposition
  REAL(KIND=r8), PARAMETER   , PUBLIC :: yls = 0.30_r8 ! leaf structural litter decomposition
  REAL(KIND=r8), PARAMETER   , PUBLIC :: yrs = 0.30_r8 ! root structural litter decomposition
  REAL(KIND=r8), PARAMETER   , PUBLIC :: yws = 0.30_r8 ! woody structural litter decomposition
  REAL(KIND=r8), PARAMETER   , PUBLIC :: yll = 1.00_r8 ! leaf lignin
  REAL(KIND=r8), PARAMETER   , PUBLIC :: yrl = 1.00_r8 ! root lignin
  REAL(KIND=r8), PARAMETER   , PUBLIC :: ywl = 1.00_r8 ! wood lignin


  !      ybn: biomass       --> non-protected pool
  !      ybp: biomass       --> protected pool
  !      yps: protected     --> passive
  !      yns: non-protected --> passive
  !      ysb: passive pool  --> biomass
  !      ypb: protected     --> biomass
  !      ynb: non-protected --> biomass

  REAL(KIND=r8), PARAMETER   , PUBLIC :: ybn  = 1.00_r8! microbial biomass to nonprotected om
  REAL(KIND=r8), PARAMETER   , PUBLIC :: ybp  = 1.00_r8! microbial biomass to protected om
  REAL(KIND=r8), PARAMETER   , PUBLIC :: yps  = 1.00_r8! protected om to passive c
  REAL(KIND=r8), PARAMETER   , PUBLIC :: yns  = 1.00_r8! nonprotected om to passive  c
  REAL(KIND=r8), PARAMETER   , PUBLIC :: ysb  = 0.20_r8! passive c to biomass 
  REAL(KIND=r8), PARAMETER   , PUBLIC :: ypb  = 0.20_r8! protected om to biomass
  REAL(KIND=r8), PARAMETER   , PUBLIC :: ynb  = 0.25_r8! nonprotected om to biomass 
  !
  !    #    #    #     #     #####     #      ##    #
  !    #    ##   #     #       #       #     #  #   #
  !    #    # #  #     #       #       #    #    #  #
  !    #    #  # #     #       #       #    ######  #
  !    #    #   ##     #       #       #    #    #  #
  !    #    #    #     #       #       #    #    #  ######
  !


  REAL(KIND=r8)   , PUBLIC :: stef   ! stefan-boltzmann constant (W m-2 K-4)
  REAL(KIND=r8)   , PUBLIC :: vonk   ! von karman constant (dimensionless)
  REAL(KIND=r8)   , PUBLIC :: grav   ! gravitational acceleration (m s-2)
  REAL(KIND=r8)   , PUBLIC :: tmelt  ! freezing point of water (K)
  REAL(KIND=r8)   , PUBLIC :: hvap   ! latent heat of vaporization of water (J kg-1)
  REAL(KIND=r8)   , PUBLIC :: hfus   ! latent heat of fusion of water (J kg-1)
  REAL(KIND=r8)   , PUBLIC :: hsub   ! latent heat of sublimation of ice (J kg-1)
  REAL(KIND=r8)   , PUBLIC :: ch2o   ! specific heat of liquid water (J deg-1 kg-1)
  REAL(KIND=r8)   , PUBLIC :: cice   ! specific heat of ice (J deg-1 kg-1)
  REAL(KIND=r8)   , PUBLIC :: cair   ! specific heat of dry air at constant pressure (J deg-1 kg-1)
  REAL(KIND=r8)   , PUBLIC :: cvap   ! specific heat of water vapor at constant pressure (J deg-1 kg-1)
  REAL(KIND=r8)   , PUBLIC :: rair   ! gas constant for dry air (J deg-1 kg-1)
  REAL(KIND=r8)   , PUBLIC :: rvap   ! gas constant for water vapor (J deg-1 kg-1)
  REAL(KIND=r8)   , PUBLIC :: cappa  ! rair/cair
  REAL(KIND=r8)   , PUBLIC :: rhow   ! density of liquid water (all types) (kg m-3)
  REAL(KIND=r8)   , ALLOCATABLE, PUBLIC :: vzero  (:,:)! a REAL(KIND=r8) array of zeros, of length npoi
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: td       (:,:)    ! global! daily average temperature (K)

  REAL(KIND=r8)   , PUBLIC :: epsilon! small quantity to avoid zero-divides and other
  ! truncation or machine-limit troubles with small
  ! values. should be slightly greater than o(1)
  ! machine precision QSfc0SSiB qsfc0tsfcm qsfcm
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: Tsfc0IBIS(:,:) 
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: Qsfc0IBIS(:,:) 
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: TsfcmIBIS(:,:) 
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: QsfcmIBIS(:,:) 

  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: tsoi0  (:,:,:)! soil temperature for each layer (K)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: tsoi   (:,:,:)! soil temperature for each layer (K)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: tsoim  (:,:,:)! soil temperature for each layer (K)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: hvasug(:,:)  ! latent heat of vap/subl, for soil surface (J kg-1)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: hvasui(:,:)  ! latent heat of vap/subl, for snow surface (J kg-1)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: wsoim (:,:,:)! fraction of soil pore space containing liquid water
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: wsoi  (:,:,:)! fraction of soil pore space containing liquid water
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: wsoi0 (:,:,:)! fraction of soil pore space containing liquid water
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: wisoi(:,:,:)! fraction of soil pore space containing ice
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: consoi(:,:,:)! thermal conductivity of each soil layer (W m-1 K-1)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: wpud (:,:)  ! liquid content of puddles per soil area (kg m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: wipud(:,:)  ! ice content of puddles per soil area (kg m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: qglif(:,:,:)! 1: fraction of soil evap (fvapg) from soil liquid
  ! 2: fraction of soil evap (fvapg) from soil ice
  ! 3: fraction of soil evap (fvapg) from puddle liquid
  ! 4: fraction of soil evap (fvapg) from puddle ice

  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: z0soi(:,:)  ! roughness length of soil surface (m)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: tg    (:,:)  ! soil skin temperature (K)

  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: ti   (:,:)  ! snow skin temperature (K)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: albsav(:,:) ! saturated soil surface albedo (visible waveband)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: albsan(:,:) ! saturated soil surface albedo (near-ir waveband)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: rhosoi(:,:,:)! soil density (without pores, not bulk) (kg m-3)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: csoi  (:,:,:)! specific heat of soil, no pore spaces (J kg-1 deg-1)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: poros (:,:,:)! porosity (mass of h2o per unit vol at sat / rhow)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: sfield(:,:,:)! field capacity soil moisture value (fraction of pore space)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: swilt (:,:,:)! wilting soil moisture value (fraction of pore space)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: bex   (:,:,:)! exponent "b" in soil water potential
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: upsoiu(:,:,:)! soil water uptake from transpiration (kg_h2o m-2 s-1)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: upsoil(:,:,:)! soil water uptake from transpiration (kg_h2o m-2 s-1)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: heatg (:,:)! net heat flux into soil surface (W m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: heati (:,:)! net heat flux into snow surface (W m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: porosflo(:,:,:)! porosity after reduction by ice content
  INTEGER      , PUBLIC  , ALLOCATABLE :: ibex  (:,:,:)! nint(bex), used for cpu speed
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: hflo  (:,:,:)! downward heat transport through soil layers (W m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: suction(:,:,:)! saturated matric potential (m-h2o)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: hydraul(:,:,:)! saturated hydraulic conductivity (m/s)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: stressl(:,:,:)! soil moisture stress factor for the lower canopy (dimensionless)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: stressu(:,:,:)! soil moisture stress factor for the upper canopy (dimensionless)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: stresstu(:,:)! sum of stressu over all 6 soil layers (dimensionless)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: stresstl (:,:)! sum of stressl over all 6 soil layers (dimensionless)

  REAL(KIND=r8), PUBLIC   :: z0sno                    ! roughness length of snow surface (m)
  REAL(KIND=r8), PUBLIC   :: rhos                    ! density of snow (kg m-3)
  REAL(KIND=r8), PUBLIC   :: consno                    ! thermal conductivity of snow (W m-1 K-1)
  REAL(KIND=r8), PUBLIC   :: hsnotop                    ! thickness of top snow layer (m)
  REAL(KIND=r8), PUBLIC   :: hsnomin                    ! minimum total thickness of snow (m)
  REAL(KIND=r8), PUBLIC   :: fimin                    ! minimum fractional snow cover
  REAL(KIND=r8), PUBLIC   :: fimax                    ! maximum fractional snow cover
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: fi          (:,:)  ! fractional snow cover
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: tsno   (:,:,:) ! temperature of snow layers (K)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: hsno   (:,:,:)! thickness of snow layers (m)

  INTEGER, PUBLIC, ALLOCATABLE :: iwet   (:,:) ! wet day / dry day flag
  INTEGER, PUBLIC, ALLOCATABLE :: iwetday(:,:,:) 
  REAL(KIND=r8), PUBLIC         , ALLOCATABLE :: precipday(:,:,:)       
  REAL(KIND=r8), PUBLIC         , ALLOCATABLE :: asurd    (:,:,:) ! direct albedo of surface system
  REAL(KIND=r8), PUBLIC         , ALLOCATABLE :: asuri    (:,:,:)! diffuse albedo of surface system 
  REAL(KIND=r8), PUBLIC         , ALLOCATABLE :: xstore  (:,:,:) ! weather generator 'memory' matrix
  !
  !      INCLUDE 'comhyd.h'
  !
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: ginvap (:,:)! total evaporation rate from all intercepted h2o (kg_h2o m-2 s-1)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: gsuvap (:,:)! total evaporation rate from surface (snow/soil) (kg_h2o m-2 s-1)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: gtrans (:,:)! total transpiration rate from all vegetation canopies (kg_h2o m-2 s-1)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: gtransu(:,:)! transpiration from upper canopy (kg_h2o m-2 s-1)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: gtransl(:,:)! transpiration from lower canopy (kg_h2o m-2 s-1)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: grunof (:,:)! surface runoff rate (kg_h2o m-2 s-1)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: gdrain (:,:)! drainage rate out of bottom of lowest soil layer (kg_h2o m-2 s-1)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: gadjust(:,:)! h2o flux due to adjustments in subroutine wadjust (kg_h2o m-2 s-1)

  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: wtot   (:,:)! total amount of water stored in snow, soil, puddels, and on vegetation (kg_h2o)
  !
  !      INCLUDE 'comsum.h'
  !
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: a10td               (:,:)! 10-day average daily air temperature (K)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: a10ancub     (:,:)! 10-day average canopy photosynthesis rate - broadleaf (mol_co2 m-2 s-1)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: a10ancuc     (:,:)! 10-day average canopy photosynthesis rate - conifer (mol_co2 m-2 s-1)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: a10ancls     (:,:)! 10-day average canopy photosynthesis rate - shrubs (mol_co2 m-2 s-1)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: a10ancl4     (:,:)! 10-day average canopy photosynthesis rate - c4 grasses (mol_co2 m-2 s-1)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: a10ancl3     (:,:)! 10-day average canopy photosynthesis rate - c3 grasses (mol_co2 m-2 s-1)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: a10scalparamu(:,:)! 10-day average day-time scaling parameter - upper canopy (dimensionless)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: a10daylightu (:,:)! 10-day average day-time PAR - upper canopy (micro-Ein m-2 s-1)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: a10scalparaml(:,:)! 10-day average day-time scaling parameter - lower canopy (dimensionless)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: a10daylightl (:,:)! 10-day average day-time PAR - lower canopy (micro-Ein m-2 s-1)

  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: adrain    (:,:)! daily average rainfall rate (mm/day)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: adsnow    (:,:)! daily average snowfall rate (mm/day)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: adaet            (:,:)! daily average aet (mm/day)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: adtrunoff (:,:)! daily average total runoff (mm/day)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: adsrunoff (:,:)! daily average surface runoff (mm/day)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: addrainage(:,:)! daily average drainage (mm/day)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: adrh            (:,:)! daily average rh (percent)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: adsnod    (:,:)! daily average snow depth (m)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: adsnof    (:,:)! daily average snow fraction (fraction)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: adwsoi    (:,:)! daily average soil moisture (fraction)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: adtsoi    (:,:)! daily average soil temperature (c)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: adwisoi   (:,:)! daily average soil ice (fraction)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: adtlaysoi (:,:)! daily average soil temperature (c) of top layer
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: adwlaysoi (:,:)! daily average soil moisture of top layer(fraction)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: adwsoic   (:,:)! daily average soil moisture using root profile weighting (fraction)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: adtsoic   (:,:)! daily average soil temperature (c) using profile weighting
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: adco2mic  (:,:)! daily accumulated co2 respiration from microbes (kg_C m-2 /day)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: adco2root (:,:)! daily accumulated co2 respiration from roots (kg_C m-2 /day)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: adco2soi  (:,:)! daily accumulated co2 respiration from soil(total) (kg_C m-2 /day)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: adco2ratio(:,:)! ratio of root to total co2 respiration

  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: adnmintot (:,:)! daily accumulated net nitrogen mineralization (kg_N m-2 /day)

  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amtemp     (:,:)! monthly average air temperature (C)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amrain     (:,:)! monthly average rainfall rate (mm/day)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amsnow     (:,:)! monthly average snowfall rate (mm/day)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amaet      (:,:)! monthly average aet (mm/day)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amtrunoff  (:,:)! monthly average total runoff (mm/day)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amsrunoff  (:,:)! monthly average surface runoff (mm/day)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amdrainage (:,:)! monthly average drainage (mm/day)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amqa       (:,:)! monthly average specific humidity (kg-h2o/kg-air)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amsolar    (:,:)! monthly average incident solar radiation (W/m**2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amirup     (:,:)! monthly average upward ir radiation (W/m**2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amirdown   (:,:)! monthly average downward ir radiation (W/m**2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amsens     (:,:)! monthly average sensible heat flux (W/m**2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amlatent   (:,:)! monthly average latent heat flux (W/m**2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amlaiu     (:,:)! monthly average lai for upper canopy (m**2/m**2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amlail     (:,:)! monthly average lai for lower canopy (m**2/m**2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amtsoi     (:,:)! monthly average 1m soil temperature (C)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amwsoi     (:,:)! monthly average 1m soil moisture (fraction)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amwisoi    (:,:)! monthly average 1m soil ice (fraction)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amvwc      (:,:)! monthly average 1m volumetric water content (fraction)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amawc      (:,:)! monthly average 1m plant-available water content (fraction)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amsnod     (:,:)! monthly average snow depth (m)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amsnof     (:,:)! monthly average snow fraction (fraction)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amco2mic   (:,:)! monthly total CO2 flux from microbial respiration (kg-C/m**2/month)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amco2root  (:,:)! monthly total CO2 flux from soil due to root respiration (kg-C/m**2/month)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amnmintot  (:,:)! monthly total N mineralization from microbes (kg-N/m**2/month)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amnpp      (:,:,:)! monthly total npp for each plant type (kg-C/m**2/month)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amts2             (:,:)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amtransu   (:,:)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amtransl   (:,:)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amsuvap    (:,:)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: aminvap    (:,:)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amneetot   (:,:)     ! monthly total net ecosystem exchange of CO2 (kg-C/m**2/month)

  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amco2ratio(:,:)! monthly ratio of root to total co2 flux
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amnpptot  (:,:)! monthly total npp for ecosystem (kg-C/m**2/month)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amco2soi  (:,:)! monthly total soil CO2 flux from microbial
  ! and root respiration (kg-C/m**2/month)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amalbedo(:,:)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amtsoil (:,:,:)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amwsoil (:,:,:)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: amwisoil(:,:,:)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: aysolar    (:,:)! annual average incident solar radiation (w/m**2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: ayirup     (:,:)! annual average upward ir radiation (w/m**2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: ayirdown   (:,:)! annual average downward ir radiation (w/m**2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: aysens     (:,:)! annual average sensible heat flux (w/m**2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: aylatent   (:,:)! annual average latent heat flux (w/m**2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: ayprcp     (:,:)! annual average precipitation (mm/yr)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: ayaet             (:,:)! annual average aet (mm/yr)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: aytrans    (:,:)! annual average transpiration (mm/yr)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: aytrunoff  (:,:)! annual average total runoff (mm/yr)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: aysrunoff  (:,:)! annual average surface runoff (mm/yr)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: aydrainage (:,:)! annual average drainage (mm/yr)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: aydwtot    (:,:)! annual average soil+vegetation+snow water recharge (mm/yr or kg_h2o/m**2/yr)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: aygpptot   (:,:)! annual total gpp for ecosystem (kg-c/m**2/yr)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: aynpp      (:,:,:)! annual total npp for each plant type(kg-c/m**2/yr)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: aynpptot   (:,:)! annual total npp for ecosystem (kg-c/m**2/yr)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: ayco2soi   (:,:)! annual total soil CO2 flux from microbial and root respiration (kg-C/m**2/yr)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: ayneetot   (:,:)! annual total NEE for ecosystem (kg-C/m**2/yr)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: totnsoi    (:,:)! total nitrogen in soil (kg_N m-2)

  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: aywsoi     (:,:)! annual average 1m soil moisture (fraction)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: aywisoi    (:,:)! annual average 1m soil ice (fraction)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: aytsoi     (:,:)! annual average 1m soil temperature (C)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: ayvwc             (:,:)! annual average 1m volumetric water content (fraction)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: ayawc             (:,:)! annual average 1m plant-available water content (fraction)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: aystresstu (:,:)! annual average soil moisture stress parameter for upper canopy (dimensionless)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: aystresstl (:,:)! annual average soil moisture stress parameter for lower canopy (dimensionless)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: ayco2mic   (:,:)! annual total CO2 flux from microbial respiration (kg-C/m**2/yr)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: ayco2root  (:,:)! annual total CO2 flux from soil due to root respiration (kg-C/m**2/yr)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: ayrootbio  (:,:)! annual average live root biomass (kg-C / m**2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: aynmintot  (:,:)! annual total nitrogen mineralization (kg-N/m**2/yr)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: ayalit     (:,:)! aboveground litter (kg-c/m**2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: ayblit     (:,:)! belowground litter (kg-c/m**2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: aycsoi     (:,:)! total soil carbon (kg-c/m**2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: aycmic     (:,:)! total soil carbon in microbial biomass (kg-c/m**2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: ayanlit    (:,:)! aboveground litter nitrogen (kg-N/m**2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: aybnlit    (:,:)! belowground litter nitrogen (kg-N/m**2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: aynsoi     (:,:)! total soil nitrogen (kg-N/m**2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: ayalbedo   (:,:)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: aygpp             (:,:,:)! annual gross npp for each plant type(kg-c/m**2/yr)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: ayanpp     (:,:,:)! annual above-ground npp for each plant type(kg-c/m**2/yr)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: ayanpptot  (:,:)! annual above-ground npp for ecosystem (kg-c/m**2/yr)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: ayrratio   (:,:)! annual average runoff ratio (fraction)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: aytratio   (:,:)! annual average transpiration ratio (fraction)

  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: totcondub(:,:)    ! 
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: totconduc(:,:)    !
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: totcondls(:,:)    ! 
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: totcondl3(:,:)    !
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: totcondl4(:,:) !
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: frac     (:,:,:)! fraction of canopy occupied by each plant functional type
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: tu0      (:,:)! temperature of upper canopy leaves (K)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: tu       (:,:)! temperature of upper canopy leaves (K)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: tum      (:,:)! temperature of upper canopy leaves (K)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: ts       (:,:)! temperature of upper canopy stems (K)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: tl       (:,:)! temperature of lower canopy leaves & stems(K)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: topparu  (:,:)! total photosynthetically active raditaion absorbed by top leaves of upper canopy (W m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: topparl  (:,:)! total photosynthetically active raditaion absorbed by top leaves of lower canopy (W m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: agcub    (:,:)! canopy average gross photosynthesis rate - broadleaf  (mol_co2 m-2 s-1)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: agcuc    (:,:)! canopy average gross photosynthesis rate - conifer        (mol_co2 m-2 s-1)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: ancub    (:,:)! canopy average net photosynthesis rate - broadleaf        (mol_co2 m-2 s-1)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: ancuc    (:,:)! canopy average net photosynthesis rate - conifer        (mol_co2 m-2 s-1)

  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: tlsub    (:,:)! temperature of lower canopy vegetation buried by snow (K)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: t12      (:,:)! air temperature at z12 (K)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: t34      (:,:)! air temperature at z34 (K)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: q12      (:,:)! specific humidity of air at z12
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: q34      (:,:)! specific humidity of air at z34
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: ciub     (:,:)! intercellular co2 concentration - broadleaf (mol_co2/mol_air)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: ciuc     (:,:)! intercellular co2 concentration - conifer   (mol_co2/mol_air)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: cils     (:,:)! intercellular co2 concentration - shrubs    (mol_co2/mol_air)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: cil3     (:,:)! intercellular co2 concentration - c3 plants (mol_co2/mol_air)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: cil4     (:,:)! intercellular co2 concentration - c4 plants (mol_co2/mol_air)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: csub     (:,:)! leaf boundary layer co2 concentration - broadleaf (mol_co2/mol_air)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: csuc     (:,:)! leaf boundary layer co2 concentration - conifer   (mol_co2/mol_air)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: csls     (:,:)! leaf boundary layer co2 concentration - shrubs    (mol_co2/mol_air)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: csl3     (:,:)! leaf boundary layer co2 concentration - c3 plants (mol_co2/mol_air)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: csl4     (:,:)! leaf boundary layer co2 concentration - c4 plants (mol_co2/mol_air)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: gsub     (:,:)! upper canopy stomatal conductance - broadleaf  (mol_co2 m-2 s-1)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: gsuc     (:,:)! upper canopy stomatal conductance - conifer             (mol_co2 m-2 s-1)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: gsls     (:,:)! lower canopy stomatal conductance - shrubs             (mol_co2 m-2 s-1)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: gsl3     (:,:)! lower canopy stomatal conductance - c3 grasses (mol_co2 m-2 s-1)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: gsl4     (:,:)! lower canopy stomatal conductance - c4 grasses (mol_co2 m-2 s-1)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: agcls    (:,:)! canopy average gross photosynthesis rate - shrubs        (mol_co2 m-2 s-1)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: agcl4    (:,:)! canopy average gross photosynthesis rate - c4 grasses (mol_co2 m-2 s-1)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: agcl3    (:,:)! canopy average gross photosynthesis rate - c3 grasses (mol_co2 m-2 s-1)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: ancls    (:,:)! canopy average net photosynthesis rate - shrubs        (mol_co2 m-2 s-1)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: ancl4    (:,:)! canopy average net photosynthesis rate - c4 grasses        (mol_co2 m-2 s-1)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: ancl3    (:,:)! canopy average net photosynthesis rate - c3 grasses        (mol_co2 m-2 s-1)

  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: clitlm   (:,:)! carbon in leaf litter pool - metabolic       (kg_C m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: clitls   (:,:)! carbon in leaf litter pool - structural      (kg_C m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: clitll   (:,:)! carbon in leaf litter pool - lignin                      (kg_C m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: clitrm   (:,:)! carbon in fine root litter pool - metabolic  (kg_C m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: clitrs   (:,:)! carbon in fine root litter pool - structural (kg_C m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: clitrl   (:,:)! carbon in fine root litter pool - lignin     (kg_C m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: clitwm   (:,:)! carbon in woody litter pool - metabolic      (kg_C m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: clitws   (:,:)! carbon in woody litter pool - structural     (kg_C m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: clitwl   (:,:)! carbon in woody litter pool - lignin                (kg_C m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: totcmic  (:,:)! total carbon residing in microbial pools (kg_C m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: csoislop (:,:)! carbon in soil - slow protected humus         (kg_C m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: csoislon (:,:)! carbon in soil - slow nonprotected humus     (kg_C m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: csoipas  (:,:)! carbon in soil - passive humus             (kg_C m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: totlit   (:,:)! total carbon in all litter pools (kg_C m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: totnlit  (:,:)! total nitrogen in all litter pools (kg_N m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: totfall  (:,:)! total litterfall and root turnover (kg_C m-2/year)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: totalit  (:,:)! total standing aboveground litter (kg_C m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: totrlit  (:,:)! total root litter carbon belowground (kg_C m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: totanlit (:,:)! total standing aboveground nitrogen in litter (kg_N m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: totrnlit (:,:)! total root litter nitrogen belowground (kg_N m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: totcsoi  (:,:)! total carbon in all soil pools (kg_C m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: totnmic  (:,:)! total nitrogen residing in microbial pool (kg_N m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: tco2mic  (:,:)! instantaneous microbial co2 flux from soil (mol-CO2 / m-2 / second)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: tnpptot  (:,:)! instantaneous npp (mol-CO2 / m-2 / second)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: tneetot  (:,:)! instantaneous net ecosystem exchange of co2 per timestep (kg_C m-2/timestep)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: tnmin    (:,:)! instantaneous nitrogen mineralization (kg_N m-2/timestep)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: cdisturb (:,:)! annual amount of vegetation carbon lost 
  ! to atmosphere due to fire  (biomass burning) (kg_C m-2/year)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE  :: tempu   (:,:)! cold-phenology trigger for trees (non-dimensional)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE  :: templ   (:,:)! cold-phenology trigger for grasses/shrubs (non-dimensional)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE  :: dropu   (:,:)! drought-phenology trigger for trees (non-dimensional)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE  :: dropls  (:,:)! drought-phenology trigger for shrubs (non-dimensional)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE  :: dropl4  (:,:)! drought-phenology trigger for c4 grasses (non-dimensional)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE  :: dropl3  (:,:)! drought-phenology trigger for c3 grasses (non-dimensional)
  !
  !      include 'comveg.h'
  !
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE  :: wliqu   (:,:)! intercepted liquid h2o on upper canopy leaf area (kg m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE  :: wliqs   (:,:)! intercepted liquid h2o on upper canopy stem area (kg m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE  :: wliql   (:,:)! intercepted liquid h2o on lower canopy leaf and stem area (kg m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE  :: wsnou   (:,:)! intercepted frozen h2o (snow) on upper canopy leaf area (kg m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE  :: wsnos   (:,:)! intercepted frozen h2o (snow) on upper canopy stem area (kg m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE  :: wsnol   (:,:)! intercepted frozen h2o (snow) on lower canopy leaf & stem area (kg m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE  :: su      (:,:)! air-vegetation transfer coefficients (*rhoa) for upper canopy
  ! leaves (m s-1 * kg m-3) (A39a Pollard & Thompson 1995)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: ss     (:,:)  ! air-vegetation transfer coefficients (*rhoa) for upper canopy
  ! stems (m s-1 * kg m-3) (A39a Pollard & Thompson 1995)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: sl     (:,:)  ! air-vegetation transfer coefficients (*rhoa) for lower canopy
  ! leaves & stems (m s-1*kg m-3) (A39a Pollard & Thompson 1995)

  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: agddu  (:,:)  ! annual accumulated growing degree days for bud burst,
  ! upper canopy (day-degrees)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: agddl  (:,:)  ! annual accumulated growing degree days for bud burst,
  ! lower canopy (day-degrees)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: storedn (:,:) ! total storage of N in soil profile (kg_N m-2) 
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: yrleach (:,:) ! annual total amount C leached from soil profile (kg_C m-2/yr)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: ynleach (:,:)

  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: falll    (:,:)  ! annual leaf litter fall (kg_C m-2/year)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: fallr    (:,:)  ! annual root litter input                   (kg_C m-2/year)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: fallw    (:,:)  ! annual wood litter fall                   (kg_C m-2/year)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: exist    (:,:,:)! probability of existence of each plant functional type in a gridcell
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: vegtype0 (:,:)  ! annual vegetation type - ibis classification
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: plai     (:,:,:)! total leaf area index of each plant functional type
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: sapfrac  (:,:)  ! fraction of woody biomass that is in sapwood
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: cbiol    (:,:,:)! carbon in leaf biomass pool (kg_C m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: cbior    (:,:,:)! carbon in fine root biomass pool (kg_C m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: cbiow    (:,:,:)! carbon in woody biomass pool (kg_C m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: biomass  (:,:,:)! total biomass of each plant functional type  (kg_C m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: totlaiu  (:,:)  ! total leaf area index for the upper canopy
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: totlail  (:,:)  ! total leaf area index for the lower canopy

  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: totbiou  (:,:)  ! total biomass in the upper canopy (kg_C m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: totbiol  (:,:)  ! total biomass in the lower canopy (kg_C m-2)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: sai      (:,:,:)! current single-sided stem area index
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: fu       (:,:)  ! fraction of overall area covered by upper canopy
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: fl       (:,:)  ! fraction of snow-free area covered by lower  canopy
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: lai      (:,:,:)! canopy single-sided leaf area index (area leaf/area veg)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: zbot     (:,:,:)! height of lowest branches above ground (m)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE :: ztop     (:,:,:)! height of plant top above ground (m)
  REAL(KIND=r8), PUBLIC   :: oriev    (2)            ! fraction of leaf/stems with vertical
  REAL(KIND=r8), PUBLIC   :: orieh    (2)            ! fraction of leaf/stems with horizontal orientation
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE  :: froot    (:,:)   ! fraction of root in soil layer 
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE  :: decompl  (:,:) ! litter decomposition factor                  (dimensionless)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE  :: decomps  (:,:) ! soil organic matter decomposition factor        (dimensionless)
  REAL(KIND=r8), PUBLIC  , ALLOCATABLE  :: firefac  (:,:) ! factor that respresents the annual average fuel
  ! dryness of a grid cell, and hence characterizes the readiness to burn


  REAL(KIND=r8), PUBLIC   :: o2conc              ! o2 concentration (mol/mol)
  REAL(KIND=r8), PUBLIC   :: co2conc                      ! co2 concentration (mol/mol)
  REAL(KIND=r8), PARAMETER, PUBLIC    ::  co2init =0.000350_r8  ! initial co2 concentration in mol/mol
  REAL(KIND=r8), PARAMETER, PUBLIC    ::  o2init  =0.209000_r8  ! initial o2 concentration in mol/mol
  INTEGER(KIND=i8), PUBLIC, ALLOCATABLE :: iMaskIBIS (:,:)
  INTEGER(KIND=i8), PUBLIC, ALLOCATABLE :: MskAntIBIS(:,:)
  INTEGER, PUBLIC , ALLOCATABLE :: nlpoints(:)    ! number of land points

  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: brf    (:,:)     ! baffer 
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: lonscale(:,:)   ! longitude of nth point in degrees east
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: latscale(:,:)   ! latitude of nth point in degrees morth
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: xintopo(:,:)     ! topography (m)
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: xinveg (:,:)     ! fixed vegetation map
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: deltat (:,:)     ! absolute minimum temperature -

  ! temp on average of coldest month (C)

  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: garea  (:,:) 
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: sand   (:,:,:) ! percent sand of soil
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: clay   (:,:,:) ! percent clay of soil
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: clmwet (:,:,:) ! climatological wet days (days/month)
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: clmt   (:,:,:) ! climatological temperature (C)
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: clmtrng(:,:,:) ! climatological temp range (C)
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: clmprec(:,:,:) ! climatological precipitation (mm/day)
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: xinwind(:,:,:) ! climatological wind speed + anomaly (m s-1)
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: clmcld (:,:,:) ! climatological cloudiness (%)
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: clmq   (:,:,:) ! climatological relative humidity (%)
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: xint   (:,:,:) ! climatological temp + anomaly (C)
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: xintrng(:,:,:) ! climatological temp range + anomaly(C)
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: xinprec(:,:,:) ! climatological precipition + anomaly (mm/day)
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: xincld (:,:,:) ! climatological cloudiness + anomaly(%)
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: xinq   (:,:,:) ! climatological relative humidity + anomaly (%)
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: xinwet (:,:,:) ! climatological wet days + anomaly (days/month)
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: gdd0   (:,:)   ! growing degree days > 0C
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: gdd5   (:,:)   ! growing degree days > 5C   
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: gdd0this(:,:) ! annual total growing degree days for current year
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: tcthis  (:,:) ! coldest monthly temperature of current year (C)
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: twthis  (:,:) ! warmest monthly temperature of current year (C)
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: gdd5this(:,:) ! annual total growing degree days for current year

  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: tc          (:,:)   !  tc     = coldest monthly temperature
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: tw          (:,:)   !  tw     = warmest monthly temperature
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: tcmin  (:,:)   ! coldest daily temperature of current year (C)

  INTEGER, PUBLIC , ALLOCATABLE  :: ndtimes(:,:) ! counter for daily average calculations
  INTEGER, PUBLIC , ALLOCATABLE  :: nmtimes(:,:) ! counter for monthly average calculations
  INTEGER, PUBLIC , ALLOCATABLE  :: nytimes(:,:) ! counter for yearly average calculations
  REAL(KIND=r8), PUBLIC , ALLOCATABLE         :: nppdummy(:,:,:)! local ! canopy NPP before accounting for stem and root respiration
  REAL(KIND=r8), PUBLIC , ALLOCATABLE         :: tco2root(:,:)    ! local instantaneous fine co2 flux from soil (mol-CO2 / m-2 / second)

  INTEGER(KIND=i8), ALLOCATABLE :: iMaskSSiB     (:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: zdepth        (:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: poros_sib     (:)

  REAL(KIND=r8)   , ALLOCATABLE :: bee           (:)
  REAL(KIND=r8)   , ALLOCATABLE :: phsat         (:)
  REAL(KIND=r8)   , ALLOCATABLE :: zlt_fixed     (:,:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: xcover_fixed  (:,:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: wsib(:,:)


  INTEGER, PARAMETER, PUBLIC :: ndaypm (1:12)=(/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

  !  REAL(KIND=r8), PUBLIC, PARAMETER :: lati2(96)=(/ &
  !       -88.57217_r8, -86.72253_r8, -84.86197_r8, -82.99894_r8, -81.13498_r8,-79.27056_r8,-77.40589_r8,-75.54106_r8,&
  !       -73.67613_r8, -71.81113_r8, -69.94608_r8, -68.08099_r8, -66.21587_r8,-64.35073_r8,-62.48557_r8,-60.62040_r8,&
  !       -58.75521_r8, -56.89001_r8, -55.02481_r8, -53.15960_r8, -51.29438_r8,-49.42915_r8,-47.56393_r8,-45.69869_r8,&
  !       -43.83346_r8, -41.96822_r8, -40.10298_r8, -38.23774_r8, -36.37249_r8,-34.50724_r8,-32.64199_r8,-30.77674_r8,&
  !       -28.91149_r8, -27.04624_r8, -25.18099_r8, -23.31573_r8, -21.45048_r8,-19.58522_r8,-17.71996_r8,-15.85470_r8,&
  !       -13.98945_r8, -12.12419_r8, -10.25893_r8,  -8.39367_r8,  -6.52841_r8, -4.66315_r8, -2.79789_r8, -0.93263_r8,&
  !         0.93263_r8,   2.79789_r8,   4.66315_r8,   6.52841_r8,   8.39367_r8, 10.25893_r8, 12.12419_r8, 13.98945_r8,&
  !        15.85470_r8,  17.71996_r8,  19.58522_r8,  21.45048_r8,  23.31573_r8, 25.18099_r8, 27.04624_r8, 28.91149_r8,&
  !        30.77674_r8,  32.64199_r8,  34.50724_r8,  36.37249_r8,  38.23774_r8, 40.10298_r8, 41.96822_r8, 43.83346_r8,&
  !        45.69869_r8,  47.56393_r8,  49.42915_r8,  51.29438_r8,  53.15960_r8, 55.02481_r8, 56.89001_r8, 58.75521_r8,&
  !        60.62040_r8,  62.48557_r8,  64.35073_r8,  66.21587_r8,  68.08099_r8, 69.94608_r8, 71.81113_r8, 73.67613_r8,&
  !        75.54106_r8,  77.40589_r8,  79.27056_r8,  81.13498_r8,  82.99894_r8, 84.86197_r8, 86.72253_r8, 88.57217_r8/)

  PUBLIC :: InitFieldsIbis
  PUBLIC :: ReStartIBIS
  PUBLIC :: Finalize_IBIS
CONTAINS
  SUBROUTINE InitFieldsIbis (iMax_in,jMax_in,kMax_in,ibMax_in,jbMax_in,dtime_in,xres_in,&
       yres_in,idate,idatec,nfsibd,nfprt,nfsibt,fNameSibVeg,&
       fNameSibmsk,ifday ,ibMaxPerJB,tod ,fNameIBISMask,&
       fNameIBISDeltaTemp,fNameSandMask,fNameClayMask,fNameClimaTemp,RESTART,fgtmp,fgq     )
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: ibMax_in
    INTEGER, INTENT(IN   ) :: jbMax_in
    INTEGER, INTENT(IN   ) :: iMax_in
    INTEGER, INTENT(IN   ) :: jMax_in
    INTEGER, INTENT(IN   ) :: kMax_in
    REAL(KIND=r8)   , INTENT(IN   ) :: dtime_in
    REAL(KIND=r8)   , INTENT(IN   ) :: xres_in
    REAL(KIND=r8)   , INTENT(IN   ) :: yres_in
    INTEGER, INTENT(IN   ) :: idate(4) 
    INTEGER, INTENT(IN   ) :: idatec(4)
    INTEGER, INTENT(IN   ) :: nfsibd
    INTEGER, INTENT(IN   ) :: nfprt
    INTEGER, INTENT(IN   ) :: nfsibt
    CHARACTER(LEN=*), INTENT(IN   ) ::  fNameSibVeg
    CHARACTER(LEN=*), INTENT(IN   ) ::  fNameSibmsk
    INTEGER         , INTENT(IN   ) :: ifday 
    INTEGER         , INTENT(IN   ) :: ibMaxPerJB(jbMax_in)
    REAL(KIND=r8)   , INTENT(IN   ) :: tod  
    CHARACTER(LEN=*), INTENT(IN   ) ::  fNameIBISMask
    CHARACTER(LEN=*), INTENT(IN   ) ::  fNameIBISDeltaTemp
    CHARACTER(LEN=*), INTENT(IN   ) ::  fNameSandMask
    CHARACTER(LEN=*), INTENT(IN   ) ::  fNameClayMask
    CHARACTER(LEN=*), INTENT(IN   ) ::  fNameClimaTemp
    LOGICAL        , INTENT(IN   ) ::  RESTART    
    REAL(KIND=r8)   , INTENT(IN   ) :: fgtmp(ibMax_in,kMax_in,jbMax_in)
    REAL(KIND=r8)   , INTENT(IN   ) :: fgq  (ibMax_in,kMax_in,jbMax_in)
    IF (RESTART)THEN
       irestart=1       ! 1 = restart run
    ELSE
       irestart=0       ! 0 = initial run
    END IF
    xres  = xres_in
    yres  = yres_in
    ibMax = ibMax_in
    jbMax = jbMax_in
    iMax  = iMax_in
    jMax  = jMax_in
    kMax  = kMax_in
    dtime = dtime_in
    idateprev=idatec
    iyear0   =idate(4) 
    ALLOCATE(iMaskIBIS(ibMax,       jbMax));iMaskIBIS=0_i8
    ALLOCATE(MskAntIBIS(ibMax,      jbMax));MskAntIBIS=0
    ALLOCATE(nlpoints(              jbMax));nlpoints=0
    ALLOCATE(brf     (ibMax,        jbMax));brf=0.0_r8
    ALLOCATE(lonscale(ibMax,        jbMax));lonscale=0.0_r8
    ALLOCATE(latscale(ibMax,        jbMax));latscale=0.0_r8
    ALLOCATE(xintopo (ibMax,        jbMax));xintopo=0.0_r8
    ALLOCATE(garea   (ibMax,        jbMax));garea=0.0_r8
    ALLOCATE(xinveg  (ibMax,        jbMax));xinveg=0.0_r8
    ALLOCATE(deltat  (ibMax,        jbMax));deltat=0.0_r8
    ALLOCATE(sand    (ibMax,nsoilay,jbMax));sand=0.0_r8
    ALLOCATE(clay    (ibMax,nsoilay,jbMax));clay=0.0_r8
    ALLOCATE(clmwet  (ibMax,12     ,jbMax));clmwet=0.0_r8
    ALLOCATE(clmt    (ibMax,12     ,jbMax));clmt=0.0_r8
    ALLOCATE(clmtrng (ibMax,12     ,jbMax));clmtrng=0.0_r8
    ALLOCATE(clmprec (ibMax,12     ,jbMax));clmprec=0.0_r8
    ALLOCATE(xinwind (ibMax,12     ,jbMax));xinwind=0.0_r8
    ALLOCATE(clmcld  (ibMax,12     ,jbMax));clmcld=0.0_r8
    ALLOCATE(clmq    (ibMax,12     ,jbMax));clmq=0.0_r8
    ALLOCATE(xint    (ibMax,12     ,jbMax));xint=0.0_r8  
    ALLOCATE(xintrng (ibMax,12     ,jbMax));xintrng=0.0_r8
    ALLOCATE(xinprec (ibMax,12     ,jbMax));xinprec=0.0_r8
    ALLOCATE(xincld  (ibMax,12     ,jbMax));xincld=0.0_r8
    ALLOCATE(xinq    (ibMax,12     ,jbMax));xinq=0.0_r8  
    ALLOCATE(xinwet  (ibMax,12     ,jbMax));xinwet=0.0_r8 
    ALLOCATE(gdd0    (ibMax,        jbMax));gdd0=0.0_r8
    ALLOCATE(gdd5    (ibMax,        jbMax));gdd5=0.0_r8
    ALLOCATE(gdd0this(ibMax,        jbMax));gdd0this=0.0_r8
    ALLOCATE(tcthis  (ibMax,        jbMax));tcthis  =0.0_r8
    ALLOCATE(twthis  (ibMax,        jbMax));twthis  =0.0_r8
    ALLOCATE(gdd5this(ibMax,        jbMax));gdd5this=0.0_r8
    ALLOCATE(tc      (ibMax,        jbMax));tc=0.0_r8
    ALLOCATE(tw      (ibMax,        jbMax));tw=0.0_r8
    ALLOCATE(tcmin   (ibMax,        jbMax));tcmin=0.0_r8
    ALLOCATE(ndtimes (ibMax,        jbMax));ndtimes=0
    ALLOCATE(nmtimes (ibMax,        jbMax));nmtimes=0
    ALLOCATE(nytimes (ibMax,        jbMax));nytimes=0
    ALLOCATE(nppdummy(ibMax,npft   ,jbMax));nppdummy=0.0_r8! local ! canopy NPP before accounting for stem and root respiration
    ALLOCATE(tco2root(ibMax        ,jbMax));tco2root=0.0_r8 ! local instantaneous fine co2 flux from soil (mol-CO2 / m-2 / second)

    ALLOCATE(wsib    (ibMax,        jbMax));wsib=0.0_r8
    ALLOCATE(iMaskSSiB(ibMax,       jbMax));iMaskSSiB=0_i8
    ALLOCATE(exist   (ibMax,npft   ,jbMax));exist=0.0_r8
    ALLOCATE(fi      (ibMax,        jbMax));fi=0.0_r8
    ALLOCATE(Tsfc0IBIS(ibMax,       jbMax));Tsfc0IBIS   =0.0_r8
    ALLOCATE(Qsfc0IBIS(ibMax,       jbMax));Qsfc0IBIS   =0.0_r8
    ALLOCATE(TsfcmIBIS(ibMax,       jbMax));TsfcmIBIS   =0.0_r8
    ALLOCATE(QsfcmIBIS(ibMax,       jbMax));QsfcmIBIS   =0.0_r8
    ALLOCATE(tsoi0   (ibMax,nsoilay,jbMax));tsoi0=0.0_r8
    ALLOCATE(tsoi    (ibMax,nsoilay,jbMax));tsoi=0.0_r8
    ALLOCATE(tsoim   (ibMax,nsoilay,jbMax));tsoim=0.0_r8
    ALLOCATE(hvasug  (ibMax,        jbMax));hvasug=0.0_r8
    ALLOCATE(hvasui  (ibMax,        jbMax));hvasui=0.0_r8
    ALLOCATE(wsoim   (ibMax,nsoilay,jbMax));wsoim=0.0_r8
    ALLOCATE(wsoi    (ibMax,nsoilay,jbMax));wsoi=0.0_r8
    ALLOCATE(wsoi0   (ibMax,nsoilay,jbMax));wsoi0=0.0_r8
    ALLOCATE(wisoi   (ibMax,nsoilay,jbMax));wisoi=0.0_r8
    ALLOCATE(consoi  (ibMax,nsoilay,jbMax));consoi=0.0_r8
    ALLOCATE(wpud    (ibMax,        jbMax));wpud=0.0_r8
    ALLOCATE(wipud   (ibMax,        jbMax));wipud=0.0_r8
    ALLOCATE(qglif   (ibMax,4      ,jbMax));qglif=0.0_r8
    ALLOCATE(z0soi   (ibMax,        jbMax));z0soi=0.0_r8
    ALLOCATE(tg       (ibMax,       jbMax));tg=0.0_r8 
    ALLOCATE(ti      (ibMax,        jbMax));ti=0.0_r8
    ALLOCATE(albsav  (ibMax,        jbMax));albsav=0.0_r8
    ALLOCATE(albsan  (ibMax,        jbMax));albsan=0.0_r8
    ALLOCATE(rhosoi  (ibMax,nsoilay,jbMax));rhosoi=0.0_r8
    ALLOCATE(csoi    (ibMax,nsoilay,jbMax));csoi=0.0_r8
    ALLOCATE(poros   (ibMax,nsoilay,jbMax));poros=0.0_r8
    ALLOCATE(sfield  (ibMax,nsoilay,jbMax));sfield=0.0_r8
    ALLOCATE(swilt   (ibMax,nsoilay,jbMax));swilt=0.0_r8 
    ALLOCATE(bex     (ibMax,nsoilay,jbMax));bex=0.0_r8
    ALLOCATE(upsoiu  (ibMax,nsoilay,jbMax));upsoiu=0.0_r8
    ALLOCATE(upsoil  (ibMax,nsoilay,jbMax));upsoil=0.0_r8
    ALLOCATE(heatg   (ibMax,        jbMax));heatg=0.0_r8
    ALLOCATE(heati   (ibMax,        jbMax));heati=0.0_r8
    ALLOCATE(ibex    (ibMax,nsoilay,jbMax));ibex=0   
    ALLOCATE(hflo    (ibMax,nsoilay+1,jbMax));hflo=0.0_r8
    ALLOCATE(suction (ibMax,nsoilay,jbMax));suction=0.0_r8 
    ALLOCATE(hydraul (ibMax,nsoilay,jbMax));hydraul=0.0_r8
    ALLOCATE(porosflo(ibMax,nsoilay,jbMax));porosflo=0.0_r8
    ALLOCATE(stressl (ibMax,nsoilay,jbMax));stressl=0.0_r8
    ALLOCATE(stressu (ibMax,nsoilay,jbMax));stressu=0.0_r8
    ALLOCATE(stresstu(ibMax        ,jbMax));stresstu=0.0_r8 
    ALLOCATE(stresstl(ibMax        ,jbMax));stresstl=0.0_r8 
    ALLOCATE(tsno    (ibMax,nsnolay,jbMax));tsno=0.0_r8
    ALLOCATE(hsno    (ibMax,nsnolay,jbMax));hsno=0.0_r8
    ALLOCATE(iwet    (ibMax        ,jbMax));iwet=0
    ALLOCATE(iwetday (ibMax,31     ,jbMax));iwetday=0
    ALLOCATE(precipday(ibMax,31    ,jbMax));precipday=0.0_r8
    ALLOCATE(asurd   (ibMax,nband  ,jbMax));asurd=0.0_r8
    ALLOCATE(asuri   (ibMax,nband  ,jbMax));asuri=0.0_r8
    ALLOCATE(xstore  (ibMax,3      ,jbMax));xstore=0.0_r8
    ALLOCATE(ginvap  (ibMax        ,jbMax));ginvap=0.0_r8 
    ALLOCATE(gsuvap  (ibMax        ,jbMax));gsuvap=0.0_r8 
    ALLOCATE(gtrans  (ibMax        ,jbMax));gtrans=0.0_r8 
    ALLOCATE(gtransu (ibMax        ,jbMax));gtransu=0.0_r8 
    ALLOCATE(gtransl (ibMax        ,jbMax));gtransl=0.0_r8 
    ALLOCATE(grunof  (ibMax        ,jbMax));grunof=0.0_r8 
    ALLOCATE(gdrain  (ibMax        ,jbMax));gdrain=0.0_r8 
    ALLOCATE(gadjust (ibMax        ,jbMax));gadjust=0.0_r8 
    ALLOCATE(wtot    (ibMax        ,jbMax));wtot=0.0_r8         
    ALLOCATE(a10td   (ibMax        ,jbMax));a10td        =0.0_r8
    ALLOCATE(a10ancub     (ibMax        ,jbMax));a10ancub     =0.0_r8
    ALLOCATE(a10ancuc     (ibMax        ,jbMax));a10ancuc     =0.0_r8
    ALLOCATE(a10ancls     (ibMax        ,jbMax));a10ancls     =0.0_r8
    ALLOCATE(a10ancl4     (ibMax        ,jbMax));a10ancl4     =0.0_r8
    ALLOCATE(a10ancl3     (ibMax        ,jbMax));a10ancl3     =0.0_r8
    ALLOCATE(a10scalparamu(ibMax        ,jbMax));a10scalparamu=0.0_r8
    ALLOCATE(a10scalparaml(ibMax        ,jbMax));a10scalparaml=0.0_r8
    ALLOCATE(a10daylightu (ibMax        ,jbMax));a10daylightu =0.0_r8
    ALLOCATE(a10daylightl (ibMax        ,jbMax));a10daylightl =0.0_r8
    ALLOCATE(adrain       (ibMax ,jbMax));adrain     =0.0_r8
    ALLOCATE(adsnow       (ibMax ,jbMax));adsnow     =0.0_r8
    ALLOCATE(adaet        (ibMax ,jbMax));adaet      =0.0_r8
    ALLOCATE(adtrunoff    (ibMax ,jbMax));adtrunoff  =0.0_r8
    ALLOCATE(adsrunoff    (ibMax ,jbMax));adsrunoff  =0.0_r8
    ALLOCATE(addrainage   (ibMax ,jbMax));addrainage =0.0_r8
    ALLOCATE(adrh         (ibMax ,jbMax));adrh       =0.0_r8
    ALLOCATE(adsnod       (ibMax ,jbMax));adsnod     =0.0_r8
    ALLOCATE(adsnof       (ibMax ,jbMax));adsnof     =0.0_r8
    ALLOCATE(adwsoi       (ibMax ,jbMax));adwsoi     =0.0_r8
    ALLOCATE(adtsoi       (ibMax ,jbMax));adtsoi     =0.0_r8
    ALLOCATE(adwisoi      (ibMax ,jbMax));adwisoi    =0.0_r8
    ALLOCATE(adtlaysoi (ibMax ,jbMax));adtlaysoi  =0.0_r8
    ALLOCATE(adwlaysoi (ibMax ,jbMax));adwlaysoi  =0.0_r8
    ALLOCATE(adwsoic   (ibMax ,jbMax));adwsoic    =0.0_r8
    ALLOCATE(adtsoic   (ibMax ,jbMax));adtsoic    =0.0_r8
    ALLOCATE(adco2mic  (ibMax ,jbMax));adco2mic   =0.0_r8
    ALLOCATE(adco2root (ibMax ,jbMax));adco2root  =0.0_r8
    ALLOCATE(adco2soi  (ibMax ,jbMax));adco2soi    =0.0_r8
    ALLOCATE(adco2ratio(ibMax ,jbMax));adco2ratio  =0.0_r8
    ALLOCATE(adnmintot (ibMax ,jbMax));adnmintot  =0.0_r8

    ALLOCATE(amtemp    (ibMax ,jbMax));amtemp     =0.0_r8
    ALLOCATE(amrain    (ibMax ,jbMax));amrain     =0.0_r8
    ALLOCATE(amsnow    (ibMax ,jbMax));amsnow     =0.0_r8
    ALLOCATE(amaet     (ibMax ,jbMax));amaet      =0.0_r8
    ALLOCATE(amtrunoff (ibMax ,jbMax));amtrunoff  =0.0_r8
    ALLOCATE(amsrunoff (ibMax ,jbMax));amsrunoff  =0.0_r8
    ALLOCATE(amdrainage(ibMax ,jbMax));amdrainage =0.0_r8
    ALLOCATE(amqa      (ibMax ,jbMax));amqa       =0.0_r8
    ALLOCATE(amsolar   (ibMax ,jbMax));amsolar    =0.0_r8
    ALLOCATE(amirup    (ibMax ,jbMax));amirup     =0.0_r8
    ALLOCATE(amirdown  (ibMax ,jbMax));amirdown   =0.0_r8
    ALLOCATE(amsens    (ibMax ,jbMax));amsens     =0.0_r8
    ALLOCATE(amlatent  (ibMax ,jbMax));amlatent   =0.0_r8
    ALLOCATE(amlaiu    (ibMax ,jbMax));amlaiu     =0.0_r8
    ALLOCATE(amlail    (ibMax ,jbMax));amlail     =0.0_r8
    ALLOCATE(amtsoi    (ibMax ,jbMax));amtsoi     =0.0_r8
    ALLOCATE(amwsoi    (ibMax ,jbMax));amwsoi     =0.0_r8
    ALLOCATE(amwisoi   (ibMax ,jbMax));amwisoi    =0.0_r8
    ALLOCATE(amvwc     (ibMax ,jbMax));amvwc      =0.0_r8
    ALLOCATE(amawc     (ibMax ,jbMax));amawc      =0.0_r8
    ALLOCATE(amsnod    (ibMax ,jbMax));amsnod     =0.0_r8
    ALLOCATE(amsnof    (ibMax ,jbMax));amsnof     =0.0_r8
    ALLOCATE(amco2mic  (ibMax ,jbMax));amco2mic   =0.0_r8
    ALLOCATE(amco2root (ibMax ,jbMax));amco2root  =0.0_r8
    ALLOCATE(amnmintot (ibMax ,jbMax));amnmintot  =0.0_r8 
    ALLOCATE(tauwood   (ibMax, npft,jbMax));tauwood =0.0_r8    ! wood biomass turnover time constant (years)
    ALLOCATE(amnpp     (ibMax, npft,jbMax));amnpp =0.0_r8 
    ALLOCATE(amts2     (ibMax ,jbMax));amts2      =0.0_r8 
    ALLOCATE(amtransu  (ibMax ,jbMax));amtransu   =0.0_r8 
    ALLOCATE(amtransl  (ibMax ,jbMax));amtransl   =0.0_r8 
    ALLOCATE(amsuvap   (ibMax ,jbMax));amsuvap    =0.0_r8 
    ALLOCATE(aminvap   (ibMax ,jbMax));aminvap    =0.0_r8 
    ALLOCATE(amneetot  (ibMax ,jbMax));amneetot    =0.0_r8 
    ALLOCATE(amco2ratio(ibMax ,jbMax));amco2ratio  =0.0_r8 
    ALLOCATE(amnpptot  (ibMax ,jbMax));amnpptot    =0.0_r8 
    ALLOCATE(amco2soi  (ibMax ,jbMax));amco2soi    =0.0_r8 
    ALLOCATE(amalbedo  (ibMax ,jbMax));amalbedo=0.0_r8
    ALLOCATE(amtsoil   (ibMax,nsoilay,jbMax));amtsoil =0.0_r8
    ALLOCATE(amwsoil   (ibMax,nsoilay,jbMax));amwsoil =0.0_r8
    ALLOCATE(amwisoil  (ibMax,nsoilay,jbMax));amwisoil=0.0_r8
    ALLOCATE(aysolar   (ibMax ,jbMax));aysolar    =0.0_r8 
    ALLOCATE(ayirup    (ibMax ,jbMax));ayirup     =0.0_r8 
    ALLOCATE(ayirdown  (ibMax ,jbMax));ayirdown   =0.0_r8 
    ALLOCATE(aysens    (ibMax ,jbMax));aysens     =0.0_r8 
    ALLOCATE(aylatent  (ibMax ,jbMax));aylatent   =0.0_r8 
    ALLOCATE(ayprcp    (ibMax ,jbMax));ayprcp     =0.0_r8 
    ALLOCATE(ayaet     (ibMax ,jbMax));ayaet      =0.0_r8 
    ALLOCATE(aytrans   (ibMax ,jbMax));aytrans    =0.0_r8 
    ALLOCATE(aytrunoff (ibMax ,jbMax));aytrunoff  =0.0_r8 
    ALLOCATE(aysrunoff (ibMax ,jbMax));aysrunoff  =0.0_r8 
    ALLOCATE(aydrainage(ibMax ,jbMax));aydrainage =0.0_r8 
    ALLOCATE(aydwtot   (ibMax ,jbMax));aydwtot =0.0_r8
    ALLOCATE(aygpptot  (ibMax ,jbMax));aydrainage =0.0_r8
    ALLOCATE(aynpp     (ibMax, npft,jbMax));aynpp  =0.0_r8 
    ALLOCATE(aynpptot  (ibMax ,jbMax));aydrainage =0.0_r8
    ALLOCATE(ayco2soi  (ibMax ,jbMax));aydrainage =0.0_r8
    ALLOCATE(ayneetot  (ibMax ,jbMax));aydrainage =0.0_r8
    ALLOCATE(totnsoi   (ibMax ,jbMax));totnsoi =0.0_r8
    ALLOCATE(aywsoi    (ibMax ,jbMax));aywsoi     =0.0_r8 
    ALLOCATE(aywisoi   (ibMax ,jbMax));aywisoi    =0.0_r8 
    ALLOCATE(aytsoi    (ibMax ,jbMax));aytsoi     =0.0_r8 
    ALLOCATE(ayvwc     (ibMax ,jbMax));ayvwc      =0.0_r8 
    ALLOCATE(ayawc     (ibMax ,jbMax));ayawc      =0.0_r8 
    ALLOCATE(aystresstu(ibMax ,jbMax));aystresstu =0.0_r8 
    ALLOCATE(aystresstl(ibMax ,jbMax));aystresstl =0.0_r8 
    ALLOCATE(ayco2mic  (ibMax ,jbMax));ayco2mic   =0.0_r8 
    ALLOCATE(ayco2root (ibMax ,jbMax));ayco2root  =0.0_r8 
    ALLOCATE(ayrootbio (ibMax ,jbMax));ayrootbio  =0.0_r8 
    ALLOCATE(aynmintot (ibMax ,jbMax));aynmintot  =0.0_r8 
    ALLOCATE(ayalit    (ibMax ,jbMax));ayalit     =0.0_r8 
    ALLOCATE(ayblit    (ibMax ,jbMax));ayblit     =0.0_r8 
    ALLOCATE(aycsoi    (ibMax ,jbMax));aycsoi     =0.0_r8 
    ALLOCATE(aycmic    (ibMax ,jbMax));aycmic     =0.0_r8 
    ALLOCATE(ayanlit   (ibMax ,jbMax));ayanlit    =0.0_r8 
    ALLOCATE(aybnlit   (ibMax ,jbMax));aybnlit    =0.0_r8 
    ALLOCATE(aynsoi    (ibMax ,jbMax));aynsoi     =0.0_r8 
    ALLOCATE(ayalbedo  (ibMax ,jbMax));ayalbedo     =0.0_r8 
    ALLOCATE(aygpp     (ibMax, npft,jbMax));aygpp =0.0_r8 
    ALLOCATE(ayanpp    (ibMax, npft,jbMax));ayanpp =0.0_r8 
    ALLOCATE(ayanpptot (ibMax ,jbMax));ayanpptot =0.0_r8 
    ALLOCATE(ayrratio  (ibMax ,jbMax));ayrratio =0.0_r8 
    ALLOCATE(aytratio  (ibMax ,jbMax));aytratio =0.0_r8 
    ALLOCATE(totcondub (ibMax ,jbMax));totcondub =0.0_r8 
    ALLOCATE(totconduc (ibMax ,jbMax));totconduc =0.0_r8 
    ALLOCATE(totcondls (ibMax ,jbMax));totcondls =0.0_r8 
    ALLOCATE(totcondl3 (ibMax ,jbMax));totcondl3 =0.0_r8 
    ALLOCATE(totcondl4 (ibMax ,jbMax));totcondl4 =0.0_r8

    ALLOCATE(frac      (ibMax, npft,jbMax));frac =0.0_r8 
    ALLOCATE(tum       (ibMax ,jbMax));tum =0.0_r8
    ALLOCATE(tu        (ibMax ,jbMax));tu =0.0_r8
    ALLOCATE(tu0       (ibMax ,jbMax));tu0 =0.0_r8
    ALLOCATE(ts        (ibMax ,jbMax));ts =0.0_r8
    ALLOCATE(tl        (ibMax ,jbMax));tl =0.0_r8
    ALLOCATE(topparu   (ibMax ,jbMax));topparu  =0.0_r8
    ALLOCATE(topparl   (ibMax ,jbMax));topparl  =0.0_r8
    ALLOCATE(agcub     (ibMax ,jbMax));agcub  =0.0_r8
    ALLOCATE(agcuc     (ibMax ,jbMax));agcuc  =0.0_r8
    ALLOCATE(ancub     (ibMax ,jbMax));ancub  =0.0_r8
    ALLOCATE(ancuc     (ibMax ,jbMax));ancuc  =0.0_r8
    ALLOCATE(tlsub     (ibMax ,jbMax));tlsub     =0.0_r8
    ALLOCATE(t12       (ibMax ,jbMax));t12       =0.0_r8
    ALLOCATE(t34       (ibMax ,jbMax));t34       =0.0_r8
    ALLOCATE(q12       (ibMax ,jbMax));q12       =0.0_r8
    ALLOCATE(q34       (ibMax ,jbMax));q34       =0.0_r8
    ALLOCATE(ciub      (ibMax ,jbMax));ciub      =0.0_r8
    ALLOCATE(ciuc      (ibMax ,jbMax));ciuc      =0.0_r8
    ALLOCATE(cils      (ibMax ,jbMax));cils      =0.0_r8
    ALLOCATE(cil3      (ibMax ,jbMax));cil3      =0.0_r8
    ALLOCATE(cil4      (ibMax ,jbMax));cil4      =0.0_r8
    ALLOCATE(csub      (ibMax ,jbMax));csub      =0.0_r8
    ALLOCATE(csuc      (ibMax ,jbMax));csuc      =0.0_r8
    ALLOCATE(csls      (ibMax ,jbMax));csls      =0.0_r8
    ALLOCATE(csl3      (ibMax ,jbMax));csl3      =0.0_r8
    ALLOCATE(csl4      (ibMax ,jbMax));csl4      =0.0_r8
    ALLOCATE(gsub      (ibMax ,jbMax));gsub      =0.0_r8
    ALLOCATE(gsuc      (ibMax ,jbMax));gsuc      =0.0_r8
    ALLOCATE(gsls      (ibMax ,jbMax));gsls      =0.0_r8
    ALLOCATE(gsl3      (ibMax ,jbMax));gsl3      =0.0_r8
    ALLOCATE(gsl4      (ibMax ,jbMax));gsl4      =0.0_r8
    ALLOCATE(agcls     (ibMax ,jbMax));agcls    =0.0_r8
    ALLOCATE(agcl4     (ibMax ,jbMax));agcl4    =0.0_r8
    ALLOCATE(agcl3     (ibMax ,jbMax));agcl3    =0.0_r8
    ALLOCATE(ancls     (ibMax ,jbMax));ancls    =0.0_r8
    ALLOCATE(ancl4     (ibMax ,jbMax));ancl4    =0.0_r8
    ALLOCATE(ancl3     (ibMax ,jbMax));ancl3    =0.0_r8
    ALLOCATE(clitlm    (ibMax ,jbMax));clitlm    =0.0_r8
    ALLOCATE(clitls    (ibMax ,jbMax));clitls    =0.0_r8
    ALLOCATE(clitll    (ibMax ,jbMax));clitll    =0.0_r8
    ALLOCATE(clitrm    (ibMax ,jbMax));clitrm    =0.0_r8
    ALLOCATE(clitrs    (ibMax ,jbMax));clitrs    =0.0_r8
    ALLOCATE(clitrl    (ibMax ,jbMax));clitrl    =0.0_r8
    ALLOCATE(clitwm    (ibMax ,jbMax));clitwm    =0.0_r8
    ALLOCATE(clitws    (ibMax ,jbMax));clitws    =0.0_r8
    ALLOCATE(clitwl    (ibMax ,jbMax));clitwl    =0.0_r8
    ALLOCATE(totcmic   (ibMax ,jbMax));totcmic   =0.0_r8
    ALLOCATE(csoislop  (ibMax ,jbMax));csoislop  =0.0_r8
    ALLOCATE(csoislon  (ibMax ,jbMax));csoislon  =0.0_r8
    ALLOCATE(csoipas   (ibMax ,jbMax));csoipas   =0.0_r8
    ALLOCATE(totlit    (ibMax ,jbMax));totlit    =0.0_r8
    ALLOCATE(totnlit   (ibMax ,jbMax));totnlit   =0.0_r8
    ALLOCATE(totfall   (ibMax ,jbMax));totfall   =0.0_r8
    ALLOCATE(totalit   (ibMax ,jbMax));totalit   =0.0_r8
    ALLOCATE(totrlit   (ibMax ,jbMax));totrlit   =0.0_r8
    ALLOCATE(totanlit  (ibMax ,jbMax));totanlit  =0.0_r8
    ALLOCATE(totrnlit  (ibMax ,jbMax));totrnlit  =0.0_r8
    ALLOCATE(totcsoi   (ibMax ,jbMax));totcsoi   =0.0_r8
    ALLOCATE(totnmic   (ibMax ,jbMax));totnmic   =0.0_r8
    ALLOCATE(tco2mic   (ibMax ,jbMax));tco2mic   =0.0_r8
    ALLOCATE(tnpptot   (ibMax ,jbMax));tnpptot   =0.0_r8
    ALLOCATE(tneetot   (ibMax ,jbMax));tneetot   =0.0_r8
    ALLOCATE(tnmin     (ibMax ,jbMax));tnmin     =0.0_r8
    ALLOCATE(cdisturb  (ibMax ,jbMax));cdisturb  =0.0_r8
    ALLOCATE(tempu     (ibMax ,jbMax));tempu     =0.0_r8
    ALLOCATE(templ          (ibMax ,jbMax));templ     =0.0_r8
    ALLOCATE(dropu          (ibMax ,jbMax));dropu     =0.0_r8
    ALLOCATE(dropls  (ibMax ,jbMax));dropls    =0.0_r8
    ALLOCATE(dropl4  (ibMax ,jbMax));dropl4    =0.0_r8
    ALLOCATE(dropl3  (ibMax ,jbMax));dropl3    =0.0_r8
    ALLOCATE(wliqu          (ibMax ,jbMax));wliqu     =0.0_r8
    ALLOCATE(wliqs          (ibMax ,jbMax));wliqs     =0.0_r8
    ALLOCATE(wliql          (ibMax ,jbMax));wliql     =0.0_r8
    ALLOCATE(wsnou          (ibMax ,jbMax));wsnou     =0.0_r8
    ALLOCATE(wsnos          (ibMax ,jbMax));wsnos     =0.0_r8
    ALLOCATE(wsnol          (ibMax ,jbMax));wsnol     =0.0_r8
    ALLOCATE(su          (ibMax ,jbMax));su            =0.0_r8
    ALLOCATE(ss      (ibMax ,jbMax));ss        =0.0_r8  
    ALLOCATE(sl      (ibMax ,jbMax));sl        =0.0_r8
    ALLOCATE(agddu   (ibMax ,jbMax));agddu     =0.0_r8    
    ALLOCATE(agddl   (ibMax ,jbMax));agddl     =0.0_r8
    ALLOCATE(storedn (ibMax ,jbMax));storedn   =0.0_r8
    ALLOCATE(yrleach (ibMax ,jbMax));yrleach   =0.0_r8
    ALLOCATE(ynleach (ibMax ,jbMax));ynleach   =0.0_r8
    ALLOCATE(falll           (ibMax ,     jbMax));falll         =0.0_r8        
    ALLOCATE(fallr           (ibMax ,     jbMax));fallr         =0.0_r8        
    ALLOCATE(fallw           (ibMax ,     jbMax));fallw         =0.0_r8        
    ALLOCATE(vegtype0 (ibMax ,     jbMax));vegtype0 =0.0_r8         
    ALLOCATE(plai           (ibMax ,npft,jbMax));plai         =0.0_r8  
    ALLOCATE(sapfrac  (ibMax ,     jbMax));sapfrac  =0.0_r8          
    ALLOCATE(cbiol           (ibMax ,npft,jbMax));cbiol         =0.0_r8  
    ALLOCATE(cbior           (ibMax ,npft,jbMax));cbior         =0.0_r8  
    ALLOCATE(cbiow           (ibMax ,npft,jbMax));cbiow         =0.0_r8  
    ALLOCATE(biomass  (ibMax ,npft,jbMax));biomass  =0.0_r8  
    ALLOCATE(totlaiu  (ibMax ,     jbMax));totlaiu  =0.0_r8          
    ALLOCATE(totlail  (ibMax ,     jbMax));totlail  =0.0_r8          
    ALLOCATE(totbiou  (ibMax ,     jbMax));totbiou =0.0_r8        
    ALLOCATE(totbiol  (ibMax ,     jbMax));totbiol =0.0_r8        
    ALLOCATE(sai           (ibMax ,2   ,jbMax));sai        =0.0_r8        
    ALLOCATE(fu           (ibMax ,     jbMax));fu        =0.0_r8        
    ALLOCATE(fl           (ibMax ,     jbMax));fl        =0.0_r8        
    ALLOCATE(lai           (ibMax ,2   ,jbMax));lai=0.0_r8                
    ALLOCATE(zbot           (ibMax ,2   ,jbMax));zbot=0.0_r8                
    ALLOCATE(ztop           (ibMax ,2   ,jbMax));ztop=0.0_r8                
    ALLOCATE(decompl  (ibMax ,        jbMax));decompl =0.0_r8        
    ALLOCATE(decomps  (ibMax ,        jbMax));decomps =0.0_r8        
    ALLOCATE(firefac  (ibMax ,        jbMax));firefac =0.0_r8        
    ALLOCATE(froot    (nsoilay,2)); froot=0.0_r8        
    ALLOCATE(vzero    (ibMax ,        jbMax));vzero =0.0_r8
    ALLOCATE(td           (ibMax ,        jbMax));td =0.0_r8


    CALL vegin(iMax,jMax, nfsibd,nfprt,nfsibt,fNameSibVeg,fNameSibmsk)

    CALL RD_PARAM()

    co2conc = co2init
    o2conc  = o2init
    CALL MsgOne("** readit **" ," Init: IBIS")

    CALL readit(fNameIBISMask,fNameIBISDeltaTemp,fNameSandMask,fNameClayMask,fNameClimaTemp)

    CALL MsgOne("** readit2 **"," Init: IBIS")

    !
    ! check if diagnostic output is requested, if so read info from 'diag.infile'
    !

    !IF (idiag .ne. 0) CALL inidiag(idiag)

    !
    ! preliminary analysis of climate data
    !
    IF (irestart == 0) THEN! 0: not a restart run 1: restart run
       CALL climanl(TminL   , &
            TminU   , &
            Twarm   , &
            GDD     , &
            gdd0    , &
            gdd5    , &
            tc      , &  
            tw      , &
            tcmin   , &
            exist   , &
            xint    , &
            deltat  , &
            npft    , &
            ndaypm    )
    END IF
    CALL MsgOne('**(InitFieldsIbis)**','Start initialization of the IBIS')

    CALL initial(iMax   , &
         jMax   , &
         kMax   , &
         ibMax  , &
         jbMax  , &
         ifday  , &
         ibMaxPerJB, &
         tod    , &
         idate  , &
         idatec , & ! INTENT(IN   ))    
         fgtmp  , &
         fgq      )

    CALL MsgOne('**(InitFieldsIbis)**','End initialization of the IBIS')
  END SUBROUTINE InitFieldsIbis


  SUBROUTINE climanl(TminL   , &! INTENT(IN   )
       TminU   , &! INTENT(IN   )
       Twarm   , &! INTENT(IN   )
       GDD     , &! INTENT(IN   )
       gdd0    , &! INTENT(INOUT)
       gdd5    , &! INTENT(INOUT)
       tc      , &! INTENT(OUT  )
       tw      , &! INTENT(OUT  )
       tcmin   , &! INTENT(INOUT)local
       exist   , &! INTENT(OUT  )
       xint    , &! INTENT(IN   )
       deltat  , &! INTENT(IN   )
       npft    , &! INTENT(IN   )
       ndaypm    )! INTENT(IN   )
    ! ---------------------------------------------------------------------
    !
    ! this subsroutine is only used to initialize growing degree days,
    ! coldest temp, and warmest temp at very beginning - provides a
    ! climate 'history' based on monthly mean values
    !
    ! common blocks
    !
    IMPLICIT NONE
    !
    INTEGER      , INTENT(IN   ) :: npft             ! number of plant functional types
    INTEGER      , INTENT(IN   ) :: ndaypm(12)       ! number of days per month      
    REAL(KIND=r8), INTENT(IN   ) :: xint  (ibMax,12,jbMax)  ! climatological temp + anomaly (C)
    REAL(KIND=r8), INTENT(IN   ) :: deltat(ibMax,jbMax)    ! absolute minimum temperature - temp on average of coldest month (C)

    REAL(KIND=r8), INTENT(INOUT) :: gdd0  (ibMax,jbMax)     ! growing degree days > 0C 
    REAL(KIND=r8), INTENT(INOUT) :: gdd5  (ibMax,jbMax)  ! growing degree days > 5C
    REAL(KIND=r8), INTENT(OUT  ) :: tc    (ibMax,jbMax)   ! coldest monthly temperature (C)
    REAL(KIND=r8), INTENT(OUT  ) :: tw    (ibMax,jbMax)     ! warmest monthly temperature (C)
    REAL(KIND=r8), INTENT(INOUT) :: tcmin (ibMax,jbMax)    ! coldest daily temperature of current year (C)
    REAL(KIND=r8), INTENT(OUT  ) :: exist (ibMax,npft,jbMax)! probability of existence of each plant functional type in a gridcell
    REAL(KIND=r8), INTENT(IN   ) :: TminL (npft)     ! Absolute minimum temperature -- lower limit (upper canopy PFTs)
    REAL(KIND=r8), INTENT(IN   ) :: TminU (npft)     ! Absolute minimum temperature -- upper limit (upper canopy PFTs)
    REAL(KIND=r8), INTENT(IN   ) :: Twarm (npft)     ! Temperature of warmest month (lower canopy PFTs)
    REAL(KIND=r8), INTENT(IN   ) :: GDD   (npft)     ! minimum GDD needed (base 5 C for upper canopy PFTs, 
    ! base 0 C for lower canopy PFTs)

    !
    ! Local variables
    !
    INTEGER :: it1w      ! indice of previous month  (interpolation)
    INTEGER :: it2w  ! indice of following month (interpolation)
    INTEGER :: i           ! loop indices
    INTEGER :: j           ! loop indices
    INTEGER :: k           ! loop indices
    INTEGER :: lda           ! loop indices
    INTEGER :: nLndPts
    !
    REAL(KIND=r8)    :: rwork     ! work variable (1/ndaypm)
    REAL(KIND=r8)    :: dt           ! used for interpolation
    REAL(KIND=r8)    :: dtemp     ! interpolated temperature

    !
    ! initialize values
    !
    gdd0=0.0_r8 !CALL const2 (gdd0, npoi, 0.0)
    gdd5=0.0_r8 !CALL const2 (gdd5, npoi, 0.0)
    !
    DO j = 1,jbMax
       nLndPts=0
       DO i = 1, ibMax
          IF (iMaskIBIS(i,j) >= 1_i8) THEN
             nLndPts=nLndPts+1
             !
             ! coldest monthly temperature (year 0) in deg c
             !
             tc(nLndPts,j) = MIN (xint(nLndPts, 1,j),  xint(nLndPts, 2,j), xint(nLndPts, 3,j), &
                  xint(nLndPts, 4,j),  xint(nLndPts, 5,j), xint(nLndPts, 6,j), &
                  xint(nLndPts, 7,j),  xint(nLndPts, 8,j), xint(nLndPts, 9,j), &
                  xint(nLndPts,10,j),  xint(nLndPts,11,j), xint(nLndPts,12,j))
             !
             ! warmest monthly temperature (year 0) in deg c
             !
             tw(nLndPts,j) = MAX (xint(nLndPts, 1,j),  xint(nLndPts, 2,j), xint(nLndPts, 3,j), &
                  xint(nLndPts, 4,j),  xint(nLndPts, 5,j), xint(nLndPts, 6,j),  &
                  xint(nLndPts, 7,j),  xint(nLndPts, 8,j), xint(nLndPts, 9,j),  &
                  xint(nLndPts,10,j),  xint(nLndPts,11,j), xint(nLndPts,12,j))
             !
             tcmin(nLndPts,j) = tc(nLndPts,j) + deltat(nLndPts,j)
             !
          END IF
       END DO
    END DO
    !
    ! interpolating climatological monthly input values to daily
    !
    DO j = 1,jbMax
       DO i = 1, nlpoints(j)
          !
          DO  k = 1, 12
             !
             rwork = 1.0_r8 / REAL(ndaypm(k),kind=r8)
             !
             DO  lda = 1, ndaypm(k)
                !
                IF (REAL(lda,kind=r8).LT.REAL(ndaypm(k)+1,kind=r8)*0.5_r8) THEN
                   it1w = k - 1
                   it2w = k
                   dt   = (REAL(lda,kind=r8) - 0.5_r8) * rwork + 0.5_r8
                ELSE
                   it1w = k
                   it2w = k + 1
                   dt   = (REAL(lda,kind=r8) - 0.5_r8) * rwork - 0.5_r8
                END IF
                !
                IF (it1w.LT. 1) it1w = 12
                IF (it2w.GT.12) it2w = 1
                !
                dtemp = xint(i,it1w,j) + dt * (xint(i,it2w,j) - xint(i,it1w,j))
                !
                ! growing degree days, using deg c
                !
                gdd0(i,j) = gdd0(i,j) + MAX(0.0_r8, dtemp)
                gdd5(i,j) = gdd5(i,j) + MAX(0.0_r8, (dtemp - 5.0_r8))
                !
             END DO
             !
          END DO
          !
       END DO
    END DO
    !
    ! call routine to determine pft existence arrays
    !

    CALL existence(TminL  , &! INTENT(IN   )
         TminU  , &! INTENT(IN   )
         Twarm  , &! INTENT(IN   )
         GDD    , &! INTENT(IN   )
         exist  , &! INTENT(OUT  )
         tcmin  , &! INTENT(IN   )
         gdd5   , &! INTENT(IN   )
         gdd0   , &! INTENT(IN   )
         tw     , &! INTENT(IN   )
         npft     )! INTENT(IN   )

    RETURN
  END  SUBROUTINE climanl
  ! ---------------------------------------------------------------------
  SUBROUTINE readit ( fNameIBISMask,fNameIBISDeltaTemp,fNameSandMask,fNameClayMask,fNameClimaTemp)
    ! ---------------------------------------------------------------------
    !
    ! reads in initialization files and initializes some fields
    !
    IMPLICIT NONE
    !
    REAL(KIND=r8)    :: lonscale_in (iMax,jMax)   ! longitude of nth point in degrees east
    REAL(KIND=r8)    :: latscale_in (iMax,jMax)   ! latitude of nth point in degrees morth
    INTEGER          ::  ier  (iMax,jMax)

    INTEGER(KIND=i8) :: mskant_in(iMax,jMax)
    INTEGER(KIND=i8) :: imask_in(iMax,jMax)
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameIBISMask  
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameIBISDeltaTemp
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameSandMask
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameClayMask
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameClimaTemp
    !REAL(KIND=r4) ::  bfr(imax,jMax)

    !
    ! Local variables
    !

    INTEGER :: i ! loop indices
    INTEGER :: j ! loop indices
    INTEGER :: k ! loop indices
    INTEGER :: irec,LRecIn
    REAL(KIND=r8) :: xlat
    INTEGER       :: nLndPts,ierr
    REAL(KIND=r8) :: array2(iMax,jMax)
    REAL(KIND=r4) :: array(iMax,jMax)
    INTEGER(i4)   :: arrayi(iMax,jMax)
    INTEGER :: ibMask (ibMax,jbMax)     ! landmask 0=water, 1=land

    !
    !CHARACTER*80 :: filen
    !
    !
    ! 2-d surface and vegetation arrays
    !
    !icount(3) = 1
    !
    ! land mask, latitudes, and longitudes
    !    
    arrayi=0  
    INQUIRE (IOLENGTH=LRecIn)arrayi(1,1)
    PRINT*,TRIM(fNameIBISMask)
    OPEN (UNIT=2,FILE=TRIM(fNameIBISMask),FORM='FORMATTED', ACCESS='SEQUENTIAL', &
         ACTION='READ',STATUS='OLD', IOSTAT=ierr)

    IF (ierr /= 0) THEN
       WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(fNameIBISMask), ierr
       STOP "**(ERROR)**"
    END IF

    irec=0
    DO j=1,jMax
       DO i=1,iMax
          irec=irec+1
          READ(2,*, IOSTAT=ierr)arrayi(i,j)
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Read file ',a,' returned iostat=',i4)") &
                  TRIM(fNameIBISMask), ierr
             STOP "**(ERROR)**"
          END IF
       END DO
    END DO

    CLOSE(2)


    DO j = 1, jMax
       DO i = 1, iMax
          IF(arrayi(i,j) >=1 )THEN
             imask_in(i,j) = 1_i8
          ELSE
             imask_in(i,j) = 0_i8
          END IF
       END DO
    END DO

    IF (reducedGrid) THEN
       CALL FreqBoxIJtoIBJB(imask_in,iMaskIBIS)
    ELSE
       CALL IJtoIBJB( imask_in,iMaskIBIS)
    END IF
    ibMask = iMaskIBIS
    iMaskIBIS = 0.0_i8
    DO i=1,iMax
       lonscale_in(i,1) = REAL((i-1),kind=r8)/REAL(iMax,kind=r8)* 360.0_r8
    END DO

    DO j=2,jMax
       DO i=1,iMax
          lonscale_in(i,j) =lonscale_in(i,1)
       END DO
    END DO

    IF (reducedGrid) THEN
       CALL NearestIJtoIBJB(lonscale_in,lonscale)
    ELSE
       CALL IJtoIBJB(lonscale_in ,lonscale )
    END IF
    !  pi --- 180
    !   x          y
    !
    !    y = 180*x/pi
    !
    DO j=1,jMax
       latscale_in(1,j) =(lati(jMax+1-j)-(pi/2.0_r8))*(180.0_r8/pi)
    END DO
    DO j=1,jMax
       DO i=2,iMax
          latscale_in(i,j) =latscale_in(1,j)
       END DO
    END DO
    IF (reducedGrid) THEN
       CALL NearestIJtoIBJB(latscale_in,latscale)
    ELSE
       CALL IJtoIBJB(latscale_in ,latscale )
    END IF
    !
    !
    ! initialize lonindex, latindex for use in arr2vec, vec2arr, etc.
    ! and calculate the approximate the area of the land gridcells
    !
    DO   j = 1, jbMax
       DO   i = 1, ibMax
          IF (ibMask(i,j) == 1) THEN
             nlpoints(j)   = nlpoints(j) + 1
             xlat       = latscale(i,j) * pi / 180.0_r8
             garea(nlpoints(j),j) = yres * 111400.0_r8 * (xres * 111400.0_r8 *COS(xlat))
          END IF
       END DO
    END DO

    CALL ibismap(iMax           , &
         jMax           , &
         latscale_in    , &
         array2            )!

    IF (reducedGrid) THEN
       CALL LinearIJtoIBJB(array2,brf)
    ELSE
       CALL IJtoIBJB(array2 ,brf )
    END IF
    DO  j = 1, jbMax
       nLndPts=0
       DO  i = 1, ibMax
          IF (ibMask(i,j) == 1) THEN
             nLndPts=nLndPts+1
             garea(nLndPts,j) = brf(i,j)
          END IF
       END DO
    END DO
    !
    ! fixed vegetation map
    !
    INQUIRE (IOLENGTH=LRecIn)arrayi(1,1)
    OPEN (UNIT=2,FILE=TRIM(fNameIBISMask),FORM='FORMATTED', ACCESS='SEQUENTIAL', &
         ACTION='READ',STATUS='OLD', IOSTAT=ierr)

    IF (ierr /= 0) THEN
       WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(fNameIBISMask), ierr
       STOP "**(ERROR)**"
    END IF

    irec=0
    DO j=1,jMax
       DO i=1,iMax
          irec=irec+1
          READ(2,*, IOSTAT=ierr)arrayi(i,j)
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Read file ',a,' returned iostat=',i4)") &
                  TRIM(fNameIBISMask), ierr
             STOP "**(ERROR)**"
          END IF
       END DO
    END DO

    CLOSE(2)
    imask_in=INT(arrayi,kind=i8)
    IF (reducedGrid) THEN
       CALL FreqBoxIJtoIBJB(imask_in,iMaskIBIS)
    ELSE
       CALL IJtoIBJB( imask_in,iMaskIBIS)
    END IF

    iMask = iMaskIBIS

    DO j = 1, jbMax
       DO i = 1, ibMax
          brf(i,j) = REAL(iMaskIBIS(i,j),KIND=r8)
       END DO
    END DO

    IF (isimveg .LE. 1) THEN
       DO j = 1, jbMax
          nLndPts=0
          DO  i = 1, ibMax
             IF (iMaskIBIS(i,j) >= 1_i8) THEN
                nLndPts=nLndPts+1
                xinveg(nLndPts,j) = brf(i,j)
             END IF
          END DO
       END DO
    END IF

    DO j=1,jMax
       DO i=1,iMax        
          IF (imask_in(i,j) >= 1_r8) THEN
             ier(i,j) = 0
          ELSE
             ier(i,j) = 1
          END IF
       END DO
       IF (ANY( ier(1:iMax,j) /= 0)) THEN
          DO i=1,iMax              
             mskant_in(i,j) = 1_i8
          END DO
       ELSE
          DO i=1,iMax              
             mskant_in(i,j) = 0_i8
          END DO
       END IF
    END DO
    IF (reducedGrid) THEN
       CALL FreqBoxIJtoIBJB(mskant_in,MskAntIBIS)
    ELSE
       CALL IJtoIBJB( mskant_in,MskAntIBIS)
    END IF
    MskAnt=MskAntIBIS


    !
    ! 2-d soil array
    !
    !
    ! delta t
    !     
    array=0.0_r4 
    INQUIRE (IOLENGTH=LRecIn)array(1,1)
    OPEN (UNIT=2,FILE=TRIM(fNameIBISDeltaTemp),FORM='FORMATTED', ACCESS='SEQUENTIAL', &
         ACTION='READ',STATUS='OLD', IOSTAT=ierr)

    IF (ierr /= 0) THEN
       WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(fNameIBISDeltaTemp), ierr
       STOP "**(ERROR)**"
    END IF
    irec=0
    DO j=1,jMax
       DO i=1,iMax
          irec=irec+1
          READ(2,*, IOSTAT=ierr)array(i,j)
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Read file ',a,' returned iostat=',i4)") &
                  TRIM(fNameIBISDeltaTemp), ierr
             STOP "**(ERROR)**"
          END IF
          array2(i,j) = REAL(array(i,j),KIND=r8)
       END DO
    END DO
    CLOSE(2)

    IF (reducedGrid) THEN
       CALL LinearIJtoIBJB(array2,brf)
    ELSE
       CALL IJtoIBJB(array2 ,brf )
    END IF

    DO j = 1, jbMax
       nLndPts=0
       DO  i = 1, ibMax
          IF (iMaskIBIS(i,j) >= 1_i8) THEN
             nLndPts=nLndPts+1
             deltat(nLndPts,j) = brf(i,j)
          END IF
       END DO
    END DO

    !
    ! 3-d soil texture array
    !
    ! icount(3) is the 6 layers used in soita.sand.nc soita.clay.nc
    !
    INQUIRE (IOLENGTH=LRecIN) arrayi(1,1)
    OPEN(2, file=TRIM(fNameSandMask), form='formatted', &
         access='SEQUENTIAL', status='old', action='read', IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(fNameSandMask), ierr
       STOP "**(ERROR)**"
    END IF
    irec=0
    DO k=1,nsoilay
       DO j=1,jMax
          DO i=1,iMax
             irec=irec+1
             READ(2,*, IOSTAT=ierr)arrayi(i,j)
             IF (ierr /= 0) THEN
                WRITE(UNIT=nfprt,FMT="('**(ERROR)** Read file ',a,' returned iostat=',i4)") &
                     TRIM(fNameSandMask), ierr
                STOP "**(ERROR)**"
             END IF
             array2(i,j) = REAL(arrayi(i,j),KIND=r8)
          END DO
       END DO

       IF (reducedGrid) THEN
          CALL FreqBoxIJtoIBJB(array2,brf)
       ELSE
          CALL IJtoIBJB(array2 ,brf )
       END IF

       DO j = 1, jbMax
          nLndPts=0
          DO i = 1, ibMax
             IF (iMaskIBIS(i,j) >= 1_i8) THEN
                nLndPts=nLndPts+1
                sand(nLndPts,k,j) = brf(i,j)
             END IF
          END DO
       END DO
    END DO
    CLOSE(2)

    INQUIRE (IOLENGTH=LRecIN) arrayi(1,1)
    OPEN(2, file=TRIM(fNameClayMask), form='FORMATTED', &
         access='SEQUENTIAL', status='old', action='read', IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(fNameSandMask), ierr
       STOP "**(ERROR)**"
    END IF
    irec=0
    DO k=1,nsoilay
       DO j=1,jMax
          DO i=1,iMax
             irec=irec+1
             READ(2,*, IOSTAT=ierr)arrayi(i,j)
             IF (ierr /= 0) THEN
                WRITE(UNIT=nfprt,FMT="('**(ERROR)** Read file ',a,' returned iostat=',i4)") &
                     TRIM(fNameClayMask), ierr
                STOP "**(ERROR)**"
             END IF
             array2(i,j) = REAL(arrayi(i,j),KIND=r8)
          END DO
       END DO
       IF (reducedGrid) THEN
          CALL FreqBoxIJtoIBJB(array2,brf)
       ELSE
          CALL IJtoIBJB(array2 ,brf )
       END IF
       DO j = 1, jbMax
          nLndPts=0
          DO i = 1, ibMax
             IF (iMaskIBIS(i,j) >= 1_i8) THEN
                nLndPts=nLndPts+1
                clay(nLndPts,k,j) = brf(i,j)
             END IF
          END DO
       END DO
    END DO
    CLOSE(2)

    !
    ! 3-d climate arrays
    !
    !     filen = 'input/wetd.mon.nc'
    !     OPEN(2,file='input/wetd.mon.bin',ACCESS='DIRECT', &
    !    &  FORM='UNFORMATTED',RECL=1,ACTION='READ')
    !      irec=0
    !     DO k=1,12
    !        DO j=1,jMax
    !          DO i=1,iMax
    !              irec=irec+1
    !              READ(2,rec=irec)array(i,j)
    !          END DO
    !        END DO
    !        DO j = 1, jbMax
    !           DO i = 1, ibMax
    !               brf(i,j) = array(i,j)
    !           END DO
    !        END DO
    !        DO j = 1, jbMax
    !            nLndPts=0
    !           DO i = 1, ibMax
    !              IF (ibMask(i,j) == 1) THEN
    !                  nLndPts=nLndPts+1
    !                 clmwet(nLndPts,k,j) = brf(i,j)
    !              END IF
    !           END DO
    !        END DO
    !     END DO     
    !     CLOSE(2)
    !
    INQUIRE (IOLENGTH=LRecIN) array(1,1)
    OPEN(2, file=TRIM(fNameClimaTemp), form='FORMATTED', &
         access='SEQUENTIAL',status='old', action='read', IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(fNameClimaTemp), ierr
       STOP "**(ERROR)**"
    END IF
    irec=0
    DO k=1,12
       DO j=1,jMax
          DO i=1,iMax
             irec=irec+1
             READ(2,*, IOSTAT=ierr)array(i,j)
             IF (ierr /= 0) THEN
                WRITE(UNIT=nfprt,FMT="('**(ERROR)** Read file ',a,' returned iostat=',i4)") &
                     TRIM(fNameClimaTemp), ierr
                STOP "**(ERROR)**"
             END IF
             array2(i,j) = REAL(array(i,j),KIND=r8)
          END DO
       END DO
       IF (reducedGrid) THEN
          CALL LinearIJtoIBJB(array2,brf)
       ELSE
          CALL IJtoIBJB(array2 ,brf )
       END IF
       DO j = 1, jbMax
          nLndPts=0
          DO i = 1, ibMax
             IF (iMaskIBIS(i,j) >= 1_i8) THEN
                nLndPts=nLndPts+1
                clmt(nLndPts,k,j) = brf(i,j)
             END IF
          END DO
       END DO
    END DO
    CLOSE(2)

    !     filen = 'input/trange.mon.nc'
    !     OPEN(2,file='input/trange.mon.bin',ACCESS='DIRECT', &
    !    &  FORM='UNFORMATTED',RECL=1,ACTION='READ')
    !      irec=0
    !     DO k=1,12
    !        DO j=1,jMax
    !          DO i=1,iMax
    !              irec=irec+1
    !              READ(2,rec=irec)array(i,j)
    !          END DO
    !        END DO
    !        DO j = 1, jbMax
    !           DO i = 1, ibMax
    !               brf(i,j) = array(i,j)
    !           END DO
    !        END DO
    !        DO j = 1, jbMax
    !            nLndPts=0
    !           DO i = 1, ibMax
    !              IF (ibMask(i,j) == 1) THEN
    !                  nLndPts=nLndPts+1
    !                 clmtrng(nLndPts,k,j) = brf(i,j)
    !              END IF
    !           END DO
    !        END DO
    !     END DO     
    !     CLOSE(2)

    !     filen = 'input/prec.mon.nc'
    !     OPEN(2,file='input/prec.mon.bin',ACCESS='DIRECT', &
    !    &  FORM='UNFORMATTED',RECL=1,ACTION='READ')
    !      irec=0
    !     DO k=1,12
    !        DO j=1,jMax
    !          DO i=1,iMax
    !              irec=irec+1
    !              READ(2,rec=irec)array(i,j)
    !          END DO
    !        END DO
    !        DO j = 1, jbMax
    !           DO i = 1, ibMax
    !               brf(i,j) = array(i,j)
    !           END DO
    !        END DO
    !        DO j = 1, jbMax
    !            nLndPts=0
    !           DO i = 1, ibMax
    !              IF (ibMask(i,j) == 1) THEN
    !                  nLndPts=nLndPts+1
    !                 clmprec(nLndPts,k,j) = brf(i,j)
    !              END IF
    !           END DO
    !        END DO
    !     END DO     
    !     CLOSE(2)

    !     filen = 'input/wspd.mon.nc'
    !     OPEN(2,file='input/wspd.mon.bin',ACCESS='DIRECT', &
    !    &  FORM='UNFORMATTED',RECL=1,ACTION='READ')
    !      irec=0
    !     DO k=1,12
    !        DO j=1,jMax
    !          DO i=1,iMax
    !              irec=irec+1
    !              READ(2,rec=irec)array(i,j)
    !          END DO
    !        END DO
    !        DO j = 1, jbMax
    !           DO i = 1, ibMax
    !               brf(i,j) = array(i,j)
    !           END DO
    !        END DO
    !        DO j = 1, jbMax
    !            nLndPts=0
    !           DO i = 1, ibMax
    !              IF (ibMask(i,j) == 1) THEN
    !                  nLndPts=nLndPts+1
    !                 xinwind(nLndPts,k,j) = brf(i,j)
    !              END IF
    !           END DO
    !        END DO
    !     END DO     
    !     CLOSE(2)
    !
    !     filen = 'input/cld.mon.nc'
    !     OPEN(2,file='input/cld.mon.bin',ACCESS='DIRECT', &
    !    &  FORM='UNFORMATTED',RECL=1,ACTION='READ')
    !      irec=0
    !     DO k=1,12
    !        DO j=1,jMax
    !          DO i=1,iMax
    !              irec=irec+1
    !              READ(2,rec=irec)array(i,j)
    !          END DO
    !        END DO
    !        DO j = 1, jbMax
    !           DO i = 1, ibMax
    !               brf(i,j) = array(i,j)
    !           END DO
    !        END DO
    !        DO j = 1, jbMax
    !            nLndPts=0
    !           DO i = 1, ibMax
    !              IF (ibMask(i,j) == 1) THEN
    !                  nLndPts=nLndPts+1
    !                 clmcld(nLndPts,k,j) = brf(i,j)
    !              END IF
    !           END DO
    !        END DO
    !     END DO     
    !     CLOSE(2)
    !
    !     filen = 'input/rh.mon.nc'
    !     OPEN(2,file='input/rh.mon.bin',ACCESS='DIRECT', &
    !    &  FORM='UNFORMATTED',RECL=1,ACTION='READ')
    !      irec=0
    !     DO k=1,12
    !        DO j=1,jMax
    !          DO i=1,iMax
    !              irec=irec+1
    !              READ(2,rec=irec)array(i,j)
    !          END DO
    !        END DO
    !         
    !        DO j = 1, jbMax
    !           DO i = 1, ibMax
    !               brf(i,j) = array(i,j)
    !           END DO
    !        END DO
    !        DO j = 1, jbMax
    !            nLndPts=0
    !           DO i = 1, ibMax
    !              IF (ibMask(i,j) == 1) THEN
    !                  nLndPts=nLndPts+1
    !                 clmq(nLndPts,k,j) = brf(i,j)
    !              END IF
    !           END DO
    !        END DO
    !     END DO     
    !     CLOSE(2)

    DO k=1,12
       DO j = 1, jbMax
          DO i = 1, nlpoints(j)
             xint   (i,k,j) = clmt   (i,k,j)
             !xintrng(i,k,j) = clmtrng(i,k,j)
             !xinprec(i,k,j) = clmprec(i,k,j)
             !xincld (i,k,j) = clmcld (i,k,j)
             !xinq   (i,k,j) = clmq   (i,k,j)
             !xinwet (i,k,j) = clmwet (i,k,j)
          END DO
       END DO
    END DO
    !
!9000 FORMAT (1x,'ERROR in subroutine readit')
!9010 FORMAT (1x,' ')
!9020 FORMAT (1x,'number of land points: ', i10)
    !
    ! return to main program

    RETURN
  END SUBROUTINE readit
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE ibismap(iMax    , &
       jMax    , &
       clat    , &
       area      )
    ! ---------------------------------------------------------------------
    !
    ! The land surface model works by gathering all the land points on the
    ! CCM3 [iMax] x [jMax] grid into a vector of [lpt] land points.
    !
    ! This subroutine finds [ixy] and [jxy], which are indices for the 
    ! mapping: [iMax] x [jMax] grid <-> [lpt] vector of land points 
    !
    ! ---------------------------------------------------------------------
    !use ibispar
    !use comwork

    IMPLICIT NONE
    !-----------------------------Arguments---------------------------------
    ! Input arguments
    INTEGER, INTENT(IN   ) :: jMax             ! number of latitude points on grid
    INTEGER, INTENT(IN   ) :: iMax             ! number of longitude points on grid
    REAL(KIND=r8)   , INTENT(IN   ) :: clat(iMax,jMax)       ! latitude of model grid [radians]
    ! Output arguments
    REAL(KIND=r8)   , INTENT(OUT  ) :: area(iMax,jMax)          ! area of land gridcell centered on lati*loni [m^2]
    !
    !------------------------------Local variables--------------------------
    INTEGER :: i
    INTEGER :: j
    REAL(KIND=r8)    :: pi                 ! pi
    REAL(KIND=r8)    :: re                 ! radius of Earth (m)
    REAL(KIND=r8)    :: areagordon(iMax,jMax)
    REAL(KIND=r8)    :: lats      (jMax+1)
    REAL(KIND=r8)    :: lonw      (iMax+1,jMax)
    REAL(KIND=r8)    :: edgen
    REAL(KIND=r8)    :: edges
    REAL(KIND=r8)    :: edgew
    REAL(KIND=r8)    :: edgee
    INTEGER :: numlon   (jMax)
    REAL(KIND=r8)    :: lonscale (iMax)
    REAL(KIND=r8)    :: latscale (jMax)
    !
    !     180 --- pi
    !      y      x 
    !-----------------------------------------------------------------------
    !
    pi = 4.0_r8*ATAN(1.0_r8)
    re = 6.371227709e+6_r8
    !
    ! lonscale and latscale are used by the netcdf input/output subroutine
    !
    DO i = 1, iMax
       lonscale(i) = REAL((i-1),kind=r8)/REAL(iMax,kind=r8)* 360.0_r8
    END DO
    DO j = 1, jMax
       !         latscale(j) = clat(1,j) * 180._r8 /pi
       latscale(j) = clat(1,j) 
    END DO
    !
    !     the area of each grid cell is calculated following Gordon's
    !     subroutines

    edgen =  90.0_r8
    edges = -90.0_r8
    edgew = -REAL(iMax,kind=r8)/(2.0_r8*360.0_r8)
    edgee = 360.0_r8-REAL(iMax,kind=r8)/(2.0_r8*360.0_r8)
    DO j = 1, jMax
       numlon(j) = iMax
    END DO
    !
    CALL cellbox(jMax    , &
         iMax    , &
         edgen   , &
         edges   , &
         edgew   , &
         edgee   , &
         numlon  , &
         lonscale, &
         latscale, &
         lats    , &
         lonw      )

    CALL cellarea(jMax   , &
         iMax   , &
         edgen  , &
         edges  , &
         edgew  , &
         edgee  , &
         numlon , &
         lats   , &
         lonw   , &
         areagordon)
    DO j=1,jMax
       DO i=1,iMax
          area(i,j) = ABS(areagordon(i,j)*1e6)
       END DO
    END DO
    RETURN
  END SUBROUTINE ibismap

  ! --------------------------------------------------------------------
  SUBROUTINE cellbox (lat   , &
       lon   , &
       edgen , &
       edges , &
       edgew , &
       edgee , &
       numlon, &
       longxy, &
       latixy, &
       lats  , &
       lonw)
    ! --------------------------------------------------------------------

    ! ------------------------ code history ------------------------------
    ! source file       : cellbox.F
    ! purpose           : southern and western edges of grid cells
    ! date first created: March 1996 - lsm version 1 (dataset generation code)
    ! by whom           : Gordon Bonan
    ! date last revised : December 1998 - lsm version 2
    ! by whom           : Gordon Bonan
    ! --------------------------------------------------------------------

    IMPLICIT NONE

    ! ------------------------ input variables ---------------------------
    INTEGER, INTENT(IN   ) :: lat             !dimension: number of latitude points
    INTEGER, INTENT(IN   ) :: lon             !dimension: number of longitude points
    REAL(KIND=r8)   , INTENT(IN   ) :: edgen           !northern edge of grid (degrees)
    REAL(KIND=r8)   , INTENT(IN   ) :: edges           !southern edge of grid (degrees)
    REAL(KIND=r8)   , INTENT(IN   ) :: edgew           !western edge of grid (degrees)
    REAL(KIND=r8)   , INTENT(IN   ) :: edgee           !eastern edge of grid (degrees)
    INTEGER, INTENT(IN   ) :: numlon(lat)     !number of grid cells per latitude strip
    !      REAL(KIND=r8)*8 longxy(lon,lat)    !longitude at center of grid cell
    !      REAL(KIND=r8)*8 latixy(lon,lat)    !latitude at center of grid cell
    REAL(KIND=r8)   , INTENT(IN   ) :: longxy(lon)    !longitude at center of grid cell
    REAL(KIND=r8)   , INTENT(IN   ) :: latixy(lat)    !latitude at center of grid cell

    ! --------------------------------------------------------------------

    ! ------------------- output variables ----------------------------
    REAL(KIND=r8)   , INTENT(INOUT) :: lats(lat + 1)        !grid cell latitude, southern edge (degrees)
    REAL(KIND=r8)   , INTENT(INOUT) :: lonw(lon + 1,lat)    !grid cell longitude, western edge (degrees)
    ! --------------------------------------------------------------------

    ! ------------------- local variables -----------------------------
    INTEGER :: i               !longitude index
    INTEGER :: j               !latitude index
    REAL(KIND=r8)    :: dx                 !cell width
    ! --------------------------------------------------------------------

    ! --------------------------------------------------------------------
    ! Latitudes -- southern edges for each latitude strip. The southern
    ! and northern edges of latitude strip j are:
    !        southern = lats(j  )
    !        northern = lats(j+1)
    ! Hence, [lats] is dimensioned lats(lat+1)
    ! --------------------------------------------------------------------     
    lats(1) = edges 
    DO j = 2, lat
       !        lats(j) = (latixy(1,j-1) + latixy(1,j)) / 2.
       lats(j) = (latixy(j-1) + latixy(j)) / 2.0_r8
    END DO
    lats(lat+1) = edgen 

    ! --------------------------------------------------------------------
    ! Longitudes -- western edges. Longitudes for the western edge of the 
    ! cells must increase continuously and span 360 degrees. Three types of 
    ! grids are valid:
    !
    ! 1: grid starts at Dateline  with western edge on Dateline
    ! 2: grid starts at Greenwich with western edge on Greenwich
    ! 3: grid starts at Greenwich with center on Greenwich
    !
    ! For Grid 1 (Dateline)          , western edges range from:  -180 to 180
    ! For Grid 2 (Greenwich)         , western edges range from:     0 to 360
    ! For Grid 3 (Greenwich centered), western edges range from: -dx/2 to -dx/2 + 360 
    !
    ! Western edges correspond to [longxy] (longitude at center of cell) for
    ! Grid 1 only. In this case, western edges range from -180 to 180 with
    ! negative longitudes west of Greenwich. Hence, this is the preferred
    ! grid type. Grids 2 and 3 are supported because some data sets start
    ! at Greenwich rather than Dateline (Grid 2) and the NCAR CCM starts
    ! at Greenwich, centered on Greenwich (Grid 3). 
    !
    ! Partial grids that do not span 360 degrees are allowed so long as they
    ! have the convention of Grid 1 with 
    !      western edge of grid: >= -180 and < 180
    !      eastern edge of grid: > western edge  and <= 180
    !
    ! [lonw] must be dimensioned lonw(lon+1,lat) because each latitude
    ! strip can have variable longitudinal resolution
    ! --------------------------------------------------------------------

    DO j = 1, lat

       dx = (edgee - edgew) / numlon(j)

       ! western edge of first grid cell

       !         lonw(1,j) = longxy(1,j) - dx/2.
       lonw(1,j) = longxy(j) - dx/2.0_r8

       ! remaining grid cells

       DO i = 2, numlon(j)+1
          lonw(i,j) = lonw(1,j) + (i-1)*dx
       END DO

       ! set unused longitudes to non-valid number

       DO i = numlon(j)+2, lon
          lonw(i,j) = -999.0_r8
       END DO

    END DO

    RETURN
  END SUBROUTINE cellbox

  ! --------------------------------------------------------------------
  SUBROUTINE cellarea (lat   , &
       lon   , &
       edgen , &
       edges , &
       edgew , &
       edgee , &
       numlon, &
       lats  , &
       lonw  , &
       area    )
    ! --------------------------------------------------------------------

    ! ------------------------ code history ------------------------------
    ! source file       : cellarea.F
    ! purpose           : area of grid cells (square kilometers)
    ! date first created: March 1996 - lsm version 1 (dataset generation code)
    ! by whom           : Gordon Bonan
    ! date last revised : December 1998 - lsm version 2
    ! by whom           : Gordon Bonan
    ! --------------------------------------------------------------------

    IMPLICIT NONE

    ! ------------------------ input variables ---------------------------
    INTEGER, INTENT(IN   ) :: lat                 !dimension: number of latitude points
    INTEGER, INTENT(IN   ) :: lon                 !dimension: number of longitude points
    REAL(KIND=r8)   , INTENT(IN   ) :: edgen              !northern edge of grid (degrees)
    REAL(KIND=r8)   , INTENT(IN   ) :: edges              !southern edge of grid (degrees)
    REAL(KIND=r8)   , INTENT(IN   ) :: edgew              !western edge of grid (degrees)
    REAL(KIND=r8)   , INTENT(IN   ) :: edgee              !eastern edge of grid (degrees)
    INTEGER, INTENT(IN   ) :: numlon(lat)     !number of grid cells per latitude strip
    REAL(KIND=r8)   , INTENT(IN   ) :: lats  (lat + 1)        !grid cell latitude, southern edge (degrees)
    REAL(KIND=r8)   , INTENT(IN   ) :: lonw  (lon + 1,lat)    !grid cell longitude, western edge (degrees)
    ! --------------------------------------------------------------------

    ! ------------------------ output variables --------------------------
    REAL(KIND=r8)   , INTENT(INOUT) :: area(lon,lat)      !cell area (km**2)
    ! --------------------------------------------------------------------

    ! ------------------------ local variables ---------------------------
    REAL(KIND=r8)    :: deg2rad            !pi/180
    REAL(KIND=r8)    :: re                   !radius of earth (km)
    REAL(KIND=r8)    :: global             !summed area
    INTEGER :: j               !latitude index
    INTEGER :: i               !longitude index
    REAL(KIND=r8)    :: dx                   !cell width: E-W
    REAL(KIND=r8)    :: dy                   !cell width: N-S
    REAL(KIND=r8)    :: error                   !true area for error check
    ! --------------------------------------------------------------------

    deg2rad = (4.0_r8*ATAN(1.0_r8)) / 180.0_r8
    re = 6371.227709_r8
    global = 0.0_r8

    DO j = 1, lat
       DO i = 1, numlon(j)
          dx = (lonw(i+1,j) - lonw(i,j)) * deg2rad
          dy = SIN(lats(j+1)*deg2rad) - SIN(lats(j)*deg2rad)
          area(i,j) = dx*dy*re*re
          global = global + area(i,j)
       END DO
    END DO

    ! make sure total area from grid cells is same as area of grid
    ! as defined by its edges

    dx = (edgee - edgew) * deg2rad
    dy = SIN(edgen*deg2rad) - SIN(edges*deg2rad)
    error =  dx*dy*re*re
    !
    IF (ABS(global-error)/error .GT. 0.001_r8) THEN

       WRITE (nfprt,*) 'CELLAREA error: correct area is ',error, &
            ' but summed area of grid cells is ',global
       STOP
    END IF

    RETURN
  END  SUBROUTINE cellarea


  !
  !    #    #    #     #     #####     #      ##    #
  !    #    ##   #     #       #       #     #  #   #
  !    #    # #  #     #       #       #    #    #  #
  !    #    #  # #     #       #       #    ######  #
  !    #    #   ##     #       #       #    #    #  #
  !    #    #    #     #       #       #    #    #  ######
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE initial (iMax   , &
       jMax   , &
       kMax   , &
       ibMax  , &
       jbMax  , &
       ifday  , &
       ibMaxPerJB, &
       tod    , &
       idate  , &
       idatec , &
       fgtmp  , &
       fgq     )! INTENT(IN   ))
    ! ---------------------------------------------------------------------
    !
    IMPLICIT NONE    
    INTEGER         , INTENT(IN   ) :: iMax  
    INTEGER         , INTENT(IN   ) :: jMax  
    INTEGER         , INTENT(IN   ) :: kMax  
    INTEGER         , INTENT(IN   ) :: ibMax 
    INTEGER         , INTENT(IN   ) :: jbMax
    INTEGER         , INTENT(IN   ) :: ifday 
    INTEGER         , INTENT(IN   ) :: ibMaxPerJB(jbMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: tod    
    INTEGER         , INTENT(IN   ) :: idate(4)  
    INTEGER         , INTENT(IN   ) :: idatec (4)! INTENT(IN  )
    REAL(KIND=r8)   , INTENT(IN   ) :: fgtmp(ibMax,kMax,jbMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: fgq  (ibMax,kMax,jbMax)


    !
    !
    IF (irestart == 0) THEN
       CALL coldstart (iMax   , &
            jMax   , &
            kMax   , &
            ibMax  , &
            jbMax  , &
            ifday  , &
            ibMaxPerJB, &
            tod    , &
            idate  , &
            idatec , &! INTENT(IN        )
            nsoilay, &! INTENT(IN        )
            nsnolay, &! INTENT(IN        )
            fgtmp  , &! INTENT(IN        )
            fgq    , &! INTENT(IN        )
            hsno   , &! INTENT(OUT  )
            tsno   , &! INTENT(OUT  )
            tsoi   , &! INTENT(OUT  )
            tsoim  , &! INTENT(OUT  )
            wsoi   , &! INTENT(OUT  )
            wsoim  , &! INTENT(OUT  )
            wisoi    )! INTENT(OUT  )
    ELSE
       CALL restart (&! INTENT(IN   ) :: 
            npft           ,&! INTENT(IN   ) :: 
            fi           ,&! INTENT(OUT  ) :: 
            tsno           ,&! INTENT(OUT  ) :: 
            hsno           ,&! INTENT(OUT  ) :: 
            tsoi           ,&! INTENT(OUT  ) :: 
            tsoim           ,&! INTENT(OUT  ) ::
            wsoi           ,&! INTENT(OUT  ) :: 
            wsoim           ,&! INTENT(OUT  ) :: 
            wisoi           ,&! INTENT(OUT  ) :: 
            cbiol           ,&! INTENT(OUT  ) :: 
            cbiow           ,&! INTENT(OUT  ) :: 
            cbior           ,&! INTENT(OUT  ) :: 
            sapfrac           ,&! INTENT(OUT  ) :: 
            clitlm           ,&! INTENT(OUT  ) :: 
            clitls           ,&! INTENT(OUT  ) :: 
            clitll           ,&! INTENT(OUT  ) :: 
            clitrm           ,&! INTENT(OUT  ) :: 
            clitrs           ,&! INTENT(OUT  ) :: 
            clitrl           ,&! INTENT(OUT  ) :: 
            clitwm           ,&! INTENT(OUT  ) :: 
            clitws           ,&! INTENT(OUT  ) :: 
            clitwl           ,&! INTENT(OUT  ) :: 
            falll           ,&! INTENT(OUT  ) :: 
            fallr           ,&! INTENT(OUT  ) :: 
            fallw           ,&! INTENT(OUT  ) :: 
            totcmic           ,&! INTENT(OUT  ) :: 
            csoislop     ,&! INTENT(OUT  ) :: 
            csoislon     ,&! INTENT(OUT  ) :: 
            csoipas           ,&! INTENT(OUT  ) :: 
            gdd0           ,&! INTENT(INOUT) :: 
            gdd5           ,&! INTENT(INOUT) :: 
            tc           ,&! INTENT(OUT  ) :: 
            tw           ,&! INTENT(INOUT) :: 
            wipud           ,&! INTENT(OUT  ) :: 
            wpud           ,&! INTENT(OUT  ) :: 
            agddu           ,&! INTENT(OUT  ) :: 
            agddl           ,&! INTENT(OUT  ) :: 
            tempu           ,&! INTENT(OUT  ) :: 
            templ           ,&! INTENT(OUT  ) :: 
            a10td           ,&! INTENT(OUT  ) :: 
            a10ancub     ,&! INTENT(OUT  ) :: 
            a10ancls     ,&! INTENT(OUT  ) :: 
            a10ancl4     ,&! INTENT(OUT  ) :: 
            a10ancl3     ,&! INTENT(OUT  ) :: 
            a10scalparamu,&! INTENT(OUT  ) :: 
            a10scalparaml,&! INTENT(OUT  ) :: 
            a10daylightu ,&! INTENT(OUT  ) :: 
            a10daylightl ,&! INTENT(OUT  ) :: 
            dropu           ,&! INTENT(OUT  ) :: 
            dropls           ,&! INTENT(OUT  ) :: 
            dropl4           ,&! INTENT(OUT  ) :: 
            dropl3           ,&! INTENT(OUT  ) :: 
            tcmin           ,&! INTENT(INOUT) :: 
            deltat       ,&! INTENT(IN   ) :: 
            exist           ,&! INTENT(OUT  ) :: 
            TminL           ,&! INTENT(IN   ) :: 
            TminU           ,&! INTENT(IN   ) :: 
            Twarm           ,&! INTENT(IN   ) :: 
            GDD             ,&! INTENT(IN   ) :: 
            tcthis          ,&! INTENT(OUT  ) :: 
            twthis           )! INTENT(OUT  ) :: 
    END IF
    !
    ! initialize physical consts, dimensions, unit numbers, lsx model
    !
    CALL inisurf(irestart , &! INTENT(IN   )
         totcondub, &! INTENT(OUT  ) :: totcondub(npoi)    
         totconduc, &! INTENT(OUT  ) :: totconduc(npoi)    
         totcondls, &! INTENT(OUT  ) :: totcondls(npoi)    
         totcondl3, &! INTENT(OUT  ) :: totcondl3(npoi)    
         totcondl4, &! INTENT(OUT  ) :: totcondl4(npoi)    
         tu       , &! INTENT(OUT  ) :: tu      (npoi)         
         ts            , &! INTENT(OUT  ) :: ts      (npoi)         
         tl            , &! INTENT(OUT  ) :: tl      (npoi)         
         tlsub    , &! INTENT(OUT  ) :: tlsub   (npoi)         
         t12            , &! INTENT(OUT  ) :: t12     (npoi)         
         t34            , &! INTENT(OUT  ) :: t34     (npoi)         
         q12            , &! INTENT(OUT  ) :: q12     (npoi)         
         q34            , &! INTENT(OUT  ) :: q34     (npoi)         
         ciub            , &! INTENT(OUT  ) :: ciub    (npoi)         
         ciuc            , &! INTENT(OUT  ) :: ciuc    (npoi)         
         cils            , &! INTENT(OUT  ) :: cils    (npoi)         
         cil3            , &! INTENT(OUT  ) :: cil3    (npoi)         
         cil4            , &! INTENT(OUT  ) :: cil4    (npoi)         
         csub            , &! INTENT(OUT  ) :: csub    (npoi)         
         csuc            , &! INTENT(OUT  ) :: csuc    (npoi)         
         csls            , &! INTENT(OUT  ) :: csls    (npoi)         
         csl3            , &! INTENT(OUT  ) :: csl3    (npoi)         
         csl4            , &! INTENT(OUT  ) :: csl4    (npoi)         
         gsub            , &! INTENT(OUT  ) :: gsub    (npoi)         
         gsuc            , &! INTENT(OUT  ) :: gsuc    (npoi)         
         gsls            , &! INTENT(OUT  ) :: gsls    (npoi)         
         gsl3            , &! INTENT(OUT  ) :: gsl3    (npoi)         
         gsl4            , &! INTENT(OUT  ) :: gsl4    (npoi)         
         clitlm   , &! INTENT(OUT  ) :: clitlm        (npoi)   !
         clitls   , &! INTENT(OUT  ) :: clitls        (npoi)   !
         clitll   , &! INTENT(OUT  ) :: clitll        (npoi)   !
         clitrm   , &! INTENT(OUT  ) :: clitrm        (npoi)   !
         clitrs   , &! INTENT(OUT  ) :: clitrs        (npoi)   !
         clitrl   , &! INTENT(OUT  ) :: clitrl        (npoi)   !
         clitwm   , &! INTENT(OUT  ) :: clitwm        (npoi)   !
         clitws   , &! INTENT(OUT  ) :: clitws        (npoi)   !
         clitwl   , &! INTENT(OUT  ) :: clitwl        (npoi)   !
         totcmic  , &! INTENT(OUT  ) :: totcmic       (npoi)   !
         csoislop , &! INTENT(OUT  ) :: csoislop      (npoi)   !
         csoislon , &! INTENT(OUT  ) :: csoislon      (npoi)   !
         csoipas  , &! INTENT(OUT  ) :: csoipas       (npoi)   !
         totlit   , &! INTENT(OUT  ) :: totlit        (npoi)   !
         totnlit  , &! INTENT(OUT  ) :: totnlit  (npoi)   !
         totfall  , &! INTENT(OUT  ) :: totfall  (npoi)   !
         totalit  , &! INTENT(OUT  ) :: totalit  (npoi)   !
         totrlit  , &! INTENT(OUT  ) :: totrlit  (npoi)   !
         totanlit , &! INTENT(OUT  ) :: totanlit (npoi)   !
         totrnlit , &! INTENT(OUT  ) :: totrnlit (npoi)   !
         totcsoi  , &! INTENT(OUT  ) :: totcsoi  (npoi)   !
         totnmic  , &! INTENT(OUT  ) :: totnmic  (npoi)   !
         tco2mic  , &! INTENT(OUT  ) :: tco2mic  (npoi)   !
         tnpptot  , &! INTENT(OUT  ) :: tnpptot  (npoi)   !
         tneetot  , &! INTENT(OUT  ) :: tneetot  (npoi)   !
         tnmin    , &! INTENT(OUT  ) :: tnmin        (npoi)   !
         cdisturb , &! INTENT(OUT  ) :: cdisturb (npoi)   !
         tempu    , &! INTENT(OUT  ) :: tempu        (npoi)
         templ    , &! INTENT(OUT  ) :: templ        (npoi)
         dropu    , &! INTENT(OUT  ) :: dropu        (npoi)
         dropls   , &! INTENT(OUT  ) :: dropls        (npoi)
         dropl4   , &! INTENT(OUT  ) :: dropl4        (npoi)
         dropl3   , &! INTENT(OUT  ) :: dropl3        (npoi)
         wliqu    , &! INTENT(OUT  ) :: wliqu        (npoi)
         wliqs    , &! INTENT(OUT  ) :: wliqs        (npoi)
         wliql    , &! INTENT(OUT  ) :: wliql        (npoi)
         wsnou    , &! INTENT(OUT  ) :: wsnou        (npoi)
         wsnos    , &! INTENT(OUT  ) :: wsnos        (npoi)
         wsnol    , &! INTENT(OUT  ) :: wsnol        (npoi)
         su            , &! INTENT(OUT  ) :: su        (npoi)
         ss            , &! INTENT(OUT  ) :: ss        (npoi)
         sl            , &! INTENT(OUT  ) :: sl
         ginvap   , &! INTENT(OUT  ) :: ginvap        (npoi)
         gsuvap   , &! INTENT(OUT  ) :: gsuvap        (npoi)
         gtrans   , &! INTENT(OUT  ) :: gtrans        (npoi)
         grunof   , &! INTENT(OUT  ) :: grunof        (npoi)
         gdrain   , &! INTENT(OUT  ) :: gdrain        (npoi)
         iwet            , &! INTENT(OUT  ) :: iwet        (npoi)    
         iwetday  , &! INTENT(OUT  ) :: iwetday  (npoi,31)
         precipday, &! INTENT(OUT  ) :: precipday(npoi,31)
         asurd    , &! INTENT(OUT  ) :: asurd        (npoi,nband)
         asuri    , &! INTENT(OUT  ) :: asuri        (npoi,nband)
         xstore   , &! INTENT(OUT  ) :: xstore        (npoi,3)
         nband    , &! INTENT(IN   ) :: nband
         stef            , &! INTENT(OUT  ) :: stef 
         vonk            , &! INTENT(OUT  ) :: vonk 
         grav            , &! INTENT(OUT  ) :: grav 
         tmelt    , &! INTENT(OUT  ) :: tmelt
         hvap            , &! INTENT(OUT  ) :: hvap 
         hfus            , &! INTENT(OUT  ) :: hfus 
         hsub            , &! INTENT(OUT  ) :: hsub 
         ch2o            , &! INTENT(OUT  ) :: ch2o 
         cice            , &! INTENT(OUT  ) :: cice 
         cair            , &! INTENT(OUT  ) :: cair 
         cvap            , &! INTENT(OUT  ) :: cvap 
         rair            , &! INTENT(OUT  ) :: rair 
         rvap            , &! INTENT(OUT  ) :: rvap 
         cappa    , &! INTENT(OUT  ) :: cappa
         rhow            , &! INTENT(OUT  ) :: rhow 
         vzero    , &! INTENT(OUT  ) :: vzero        (npoi)
         epsilon    )! INTENT(OUT  ) :: epsilon
    !
    ! initialize snow model
    !
    CALL inisnow(rhow   , & ! INTENT(IN   ) :: rhow
         nsnolay, & ! INTENT(IN   ) :: nsnolay 
         dtime  , & ! INTENT(IN   ) :: dtime
         rhos   , & ! INTENT(OUT  ) :: rhos
         consno , & ! INTENT(OUT  ) :: consno
         hsnotop, & ! INTENT(OUT  ) :: hsnotop
         hsnomin, & ! INTENT(OUT  ) :: hsnomin
         fimin  , & ! INTENT(OUT  ) :: fimin
         fimax  , & ! INTENT(OUT  ) :: fimax
         z0sno    ) ! INTENT(OUT  ) :: z0sno
    !
    ! initialize soil model
    !
    CALL inisoil(irestart  , &! INTENT(IN   )
         texdat    , &! INTENT(IN   )
         porosdat  , &! INTENT(IN   )
         sfielddat , &! INTENT(IN   )
         swiltdat  , &! INTENT(IN   )
         bexdat    , &! INTENT(IN   )
         suctiondat, &! INTENT(IN   )
         hydrauldat, &! INTENT(IN   )
         ndat      , &! INTENT(IN   )
         wpud      , &! INTENT(OUT  )
         wipud     , &! INTENT(OUT  )
         z0soi     , &! INTENT(OUT  )
         wsoi      , &! INTENT(INOUT)
         wsoim     , &! INTENT(INOUT)
         wisoi     , &! INTENT(INOUT)
         tsoi      , &! INTENT(INOUT)
         tsoim     , &! INTENT(INOUT)
         tsno      , &! INTENT(INOUT)
         tg        , &! INTENT(OUT  )
         tgm       , &! INTENT(OUT  )
         ti        , &! INTENT(OUT  )
         sand      , &! INTENT(IN   )
         clay      , &! INTENT(IN   )
         albsav    , &! INTENT(OUT  )
         albsan    , &! INTENT(OUT  )
         rhosoi    , &! INTENT(OUT  )
         csoi      , &! INTENT(OUT  )
         poros     , &! INTENT(OUT  )
         sfield    , &! INTENT(OUT  )
         swilt     , &! INTENT(OUT  )
         bex       , &! INTENT(OUT  )
         ibex      , &! INTENT(OUT  )
         suction   , &! INTENT(OUT  )
         hydraul   , &! INTENT(OUT  )
         nband     , &! INTENT(IN   )
         nsoilay   , &! INTENT(IN   )
         nsnolay     )! INTENT(IN   )
    !
    ! initialize vegetation parameters
    !
    CALL MsgOne('**(initial)**','npft iniveg')

    CALL iniveg (isimveg        , &! INTENT(IN   )
         irestart        , &! INTENT(IN   )
         plai_init       , &! INTENT(IN   )
         plaiupper       , &! INTENT(IN   )
         plailower       , &! INTENT(IN   )
         xminlai         , &! INTENT(IN   )
         sapfrac_init    , &! INTENT(IN   )
         chiflz          , &! INTENT(IN   )
         chifuz          , &! INTENT(IN   )
         beta1           , &! INTENT(IN   )
         beta2           , &! INTENT(IN   )
         agddu           , &! INTENT(OUT  )
         agddl           , &! INTENT(OUT  )
         falll           , &! INTENT(OUT  )
         fallr           , &! INTENT(OUT  )
         fallw           , &! INTENT(OUT  )
         exist           , &! INTENT(IN   )
         vegtype0        , &! INTENT(OUT  )
         plai            , &! INTENT(OUT  )
         tw              , &! INTENT(IN   )
         sapfrac         , &! INTENT(OUT  )
         cbiol           , &! INTENT(INOUT)
         specla          , &! INTENT(IN   )
         cbior           , &! INTENT(INOUT)
         cbiow           , &! INTENT(INOUT)
         biomass         , &! INTENT(OUT  )
         totlaiu         , &! INTENT(OUT  )
         totlail         , &! INTENT(OUT  )
         totbiou         , &! INTENT(OUT  )
         totbiol         , &! INTENT(OUT  )
         sai             , &! INTENT(OUT  )
         fu              , &! INTENT(OUT  )
         woodnorm        , &! INTENT(IN   )
         fl              , &! INTENT(OUT  )
         lai             , &! INTENT(OUT  )
         zbot            , &! INTENT(OUT  )
         ztop            , &! INTENT(OUT  )
         oriev           , &! INTENT(OUT  )
         orieh           , &! INTENT(OUT  )
         froot           , &! INTENT(OUT  )
         a10td           , &! INTENT(OUT  )
         a10ancub        , &! INTENT(OUT  )
         a10ancuc        , &! INTENT(OUT  )
         a10ancls        , &! INTENT(OUT  )
         a10ancl4        , &! INTENT(OUT  )
         a10ancl3        , &! INTENT(OUT  )
         a10scalparamu   , &! INTENT(OUT  )
         a10scalparaml   , &! INTENT(OUT  )
         a10daylightu    , &! INTENT(OUT  )
         a10daylightl    , &! INTENT(OUT  )
         stresstu        , &! INTENT(OUT  )
         stresstl        , &! INTENT(OUT  )
         hsoi            , &! INTENT(IN   )
         xinveg          , &! INTENT(IN   )
         nsoilay         , &! INTENT(IN   )
         gdd0this        , &! INTENT(OUT  )
         gdd5this        , &! INTENT(OUT  )
         tcthis          , &! INTENT(OUT  )
         twthis          , &! INTENT(OUT  )
         npft              )! INTENT(IN   )

    !
    ! initialize variables for time averaging
    !
    IF (irestart == 0) &
         CALL inisum (wliqu          , &! INTENT(IN   ) :: wliqu  (npoi)        
         wsnou          , &! INTENT(IN   ) :: wsnou  (npoi)        
         fu             , &! INTENT(IN   ) :: fu     (npoi)        
         lai            , &! INTENT(IN   ) :: lai    (npoi,2) 
         wliqs          , &! INTENT(IN   ) :: wliqs  (npoi)        
         wsnos          , &! INTENT(IN   ) :: wsnos  (npoi)        
         sai            , &! INTENT(IN   ) :: sai    (npoi,2) 
         wliql          , &! INTENT(IN   ) :: wliql  (npoi)        
         wsnol          , &! INTENT(IN   ) :: wsnol  (npoi)        
         fl             , &! INTENT(IN   ) :: fl     (npoi)        
         decompl        , &! INTENT(OUT  ) :: decompl(npoi)        
         decomps        , &! INTENT(OUT  ) :: decomps(npoi)        
         firefac        , &! INTENT(OUT  ) :: firefac(npoi)        
         adrain         , &! INTENT(OUT  ) :: adrain    (npoi)     
         adsnow         , &! INTENT(OUT  ) :: adsnow    (npoi)     
         adaet          , &! INTENT(OUT  ) :: adaet     (npoi)     
         adtrunoff      , &! INTENT(OUT  ) :: adtrunoff (npoi)     
         adsrunoff      , &! INTENT(OUT  ) :: adsrunoff (npoi)     
         addrainage     , &! INTENT(OUT  ) :: addrainage(npoi)     
         adrh           , &! INTENT(OUT  ) :: adrh      (npoi)     
         adsnod         , &! INTENT(OUT  ) :: adsnod    (npoi)     
         adsnof         , &! INTENT(OUT  ) :: adsnof    (npoi)     
         adwsoi         , &! INTENT(OUT  ) :: adwsoi    (npoi)     
         adtsoi         , &! INTENT(OUT  ) :: adtsoi    (npoi)     
         adwisoi        , &! INTENT(OUT  ) :: adwisoi   (npoi)     
         adtlaysoi      , &! INTENT(OUT  ) :: adtlaysoi (npoi)     
         adwlaysoi      , &! INTENT(OUT  ) :: adwlaysoi (npoi)     
         adwsoic        , &! INTENT(OUT  ) :: adwsoic   (npoi)     
         adtsoic        , &! INTENT(OUT  ) :: adtsoic   (npoi)     
         adco2mic       , &! INTENT(OUT  ) :: adco2mic  (npoi)     
         adco2root      , &! INTENT(OUT  ) :: adco2root (npoi)     
         adnmintot      , &! INTENT(OUT  ) :: adnmintot (npoi)     
         amtemp         , &! INTENT(OUT  ) :: amtemp    (npoi)     
         amrain         , &! INTENT(OUT  ) :: amrain    (npoi)     
         amsnow         , &! INTENT(OUT  ) :: amsnow    (npoi)     
         amaet          , &! INTENT(OUT  ) :: amaet     (npoi)     
         amtrunoff      , &! INTENT(OUT  ) :: amtrunoff (npoi)     
         amsrunoff      , &! INTENT(OUT  ) :: amsrunoff (npoi)     
         amdrainage     , &! INTENT(OUT  ) :: amdrainage(npoi)     
         amqa           , &! INTENT(OUT  ) :: amqa      (npoi)     
         amsolar        , &! INTENT(OUT  ) :: amsolar   (npoi)     
         amirup         , &! INTENT(OUT  ) :: amirup    (npoi)     
         amirdown       , &! INTENT(OUT  ) :: amirdown  (npoi)     
         amsens         , &! INTENT(OUT  ) :: amsens    (npoi)     
         amlatent       , &! INTENT(OUT  ) :: amlatent  (npoi)     
         amlaiu         , &! INTENT(OUT  ) :: amlaiu    (npoi)     
         amlail         , &! INTENT(OUT  ) :: amlail    (npoi)     
         amtsoi         , &! INTENT(OUT  ) :: amtsoi    (npoi)     
         amwsoi         , &! INTENT(OUT  ) :: amwsoi    (npoi)     
         amwisoi        , &! INTENT(OUT  ) :: amwisoi   (npoi)     
         amvwc          , &! INTENT(OUT  ) :: amvwc     (npoi)     
         amawc          , &! INTENT(OUT  ) :: amawc     (npoi)     
         amsnod         , &! INTENT(OUT  ) :: amsnod    (npoi)     
         amsnof         , &! INTENT(OUT  ) :: amsnof    (npoi)     
         amco2mic       , &! INTENT(OUT  ) :: amco2mic  (npoi)     
         amco2root      , &! INTENT(OUT  ) :: amco2root (npoi)     
         amnmintot      , &! INTENT(OUT  ) :: amnmintot (npoi)     
         amnpp          , &! INTENT(OUT  ) :: amnpp     (npoi,npft)
         aysolar        , &! INTENT(OUT  ) :: aysolar   (npoi)     
         ayirup         , &! INTENT(OUT  ) :: ayirup    (npoi)     
         ayirdown       , &! INTENT(OUT  ) :: ayirdown  (npoi)     
         aysens         , &! INTENT(OUT  ) :: aysens    (npoi)     
         aylatent       , &! INTENT(OUT  ) :: aylatent  (npoi)     
         ayprcp         , &! INTENT(OUT  ) :: ayprcp    (npoi)     
         ayaet          , &! INTENT(OUT  ) :: ayaet     (npoi)     
         aytrans        , &! INTENT(OUT  ) :: aytrans   (npoi)     
         aytrunoff      , &! INTENT(OUT  ) :: aytrunoff (npoi)     
         aysrunoff      , &! INTENT(OUT  ) :: aysrunoff (npoi)     
         aydrainage     , &! INTENT(OUT  ) :: aydrainage(npoi)     
         aywsoi         , &! INTENT(OUT  ) :: aywsoi    (npoi)     
         aywisoi        , &! INTENT(OUT  ) :: aywisoi   (npoi)     
         aytsoi         , &! INTENT(OUT  ) :: aytsoi    (npoi)     
         ayvwc          , &! INTENT(OUT  ) :: ayvwc     (npoi)     
         ayawc          , &! INTENT(OUT  ) :: ayawc     (npoi)     
         aystresstu     , &! INTENT(OUT  ) :: aystresstu(npoi)     
         aystresstl     , &! INTENT(OUT  ) :: aystresstl(npoi)     
         ayco2mic       , &! INTENT(OUT  ) :: ayco2mic  (npoi)     
         ayco2root      , &! INTENT(OUT  ) :: ayco2root (npoi)     
         ayrootbio      , &! INTENT(OUT  ) :: ayrootbio (npoi)     
         aynmintot      , &! INTENT(OUT  ) :: aynmintot (npoi)     
         ayalit         , &! INTENT(OUT  ) :: ayalit    (npoi)     
         ayblit         , &! INTENT(OUT  ) :: ayblit    (npoi)     
         aycsoi         , &! INTENT(OUT  ) :: aycsoi    (npoi)     
         aycmic         , &! INTENT(OUT  ) :: aycmic    (npoi)     
         ayanlit        , &! INTENT(OUT  ) :: ayanlit   (npoi)     
         aybnlit        , &! INTENT(OUT  ) :: aybnlit   (npoi)     
         aynsoi         , &! INTENT(OUT  ) :: aynsoi    (npoi)     
         aygpp          , &! INTENT(OUT  ) :: aygpp     (npoi,npft)
         wpud           , &! INTENT(IN   )
         wipud          , &! INTENT(IN   )
         poros          , &! INTENT(IN   )
         wsoi           , &! INTENT(IN   )
         wisoi          , &! INTENT(IN   )
         hsoi           , &! INTENT(IN   )
         fi             , &! INTENT(IN   )
         rhos           , &! INTENT(IN   )
         hsno           , &! INTENT(IN   )
         wtot           , &! INTENT(OUT  )
         nsoilay        , &! INTENT(IN   )
         nsnolay        , &! INTENT(IN   )
         npft           , &! INTENT(IN   )
         rhow             )! INTENT(IN   )
    !
    ! return to main program
    !
    RETURN
  END SUBROUTINE initial
  ! ---------------------------------------------------------------------
  SUBROUTINE restart (npft         ,fi           ,  &
       tsno         ,hsno         ,tsoi         ,tsoim         , &
       wsoi         ,wsoim        ,&
       wisoi        ,cbiol        ,cbiow        ,cbior        , &
       sapfrac      ,clitlm       ,clitls       ,clitll       , &
       clitrm       ,clitrs       ,clitrl       ,clitwm       , &
       clitws       ,clitwl       ,falll        ,fallr        , &
       fallw        ,totcmic      ,csoislop     ,csoislon     , &
       csoipas      ,gdd0         ,gdd5         ,tc           , &
       tw           ,wipud        ,wpud         ,agddu        , &
       agddl        ,tempu        ,templ        ,a10td        , &
       a10ancub     ,a10ancls     ,a10ancl4     ,a10ancl3     , &
       a10scalparamu,a10scalparaml,a10daylightu ,a10daylightl , &
       dropu        ,dropls       ,dropl4       ,dropl3       , &
       tcmin        ,deltat       ,exist        ,TminL        , &
       TminU        ,Twarm        ,GDD          ,tcthis       , &
       twthis)
    ! ---------------------------------------------------------------------
    !
    ! reads in restart files, initializes some variables
    !
    ! this subroutine reads the restart values of:
    !
    !  fsnocov = fractional snow cover
    !  tsno    = temperature of snow
    !  hsno    = snow depth
    !  tsoi    = soil temperature
    !  wisoi   = soil ice content
    !  wsoi    = soil moisture content
    !  cbiol   = carbon in leaf biomass pool
    !  cbiow   = carbon in woody biomass pool
    !  cbior   = carbon in fine root biomass pool
    !  sapfrac = sapwood fraction
    !  clitlm  = leaf metabolic litter
    !  clitls  = leaf structural litter
    !  clitll  = leaf lignin litter
    !  clitrm  = root metabolic litter
    !  clitrs  = root structural litter
    !  clitrl  = root lignin litter
    !  clitwm  = woody metabolic litter
    !  clitws  = woody structural litter
    !  clitwl  = woody lignin litter
    !  falll   = annual leaf litterfall 
    !  fallr   = annual fine root turnover
    !  fallw   = annual woody turnover
    !  totcmic = total microbial carbon
    !  csoislop= slow soil carbon, protected humus
    !  csoislon= slow soil carbon, nonprotected humus
    !  csoipas = passive soil carbon
    !  gdd0    = growing degree days 0
    !  gdd5    = growing degree days 5
    !  tc      = coldest monthly temperature
    !  tw      = warmest monthly temperature
    !  wipud   = ice content of puddles per soil area
    !  wpud    = liquid content of puddles per soil area
    !  agddu   = annual accumulated growing degree days for bud burst, upper canopy
    !  agddl   = annual accumulated growing degree days for bud burst, lower canopy
    !  tempu   = cold-phenology trigger for trees
    !  templ   = cold-phenology trigger for grasses/shrubs
    !  a10td    = 10-day avg daily temp
    !  a10ancub = 10-day average canopy photosynthesis rate - broadleaf
    !  a10ancuc = 10-day average canopy photosynthesis rate - conifer
    !  a10ancls = 10-day average canopy photosynthesis rate - shrubs
    !  a10ancl4 = 10-day average canopy photosynthesis rate - c4 grasses
    !  a10ancl3 = 10-day average canopy photosynthesis rate - c3 grasses
    !  a10scalparamu = 10-day average canopy scaling parameter - upper canopy
    !  a10scalparaml = 10-day average canopy scaling parameter - lower canopy
    !  a10daylightu = 10-day average daylight - upper canopy
    !  a10daylightl = 10-day average daylight - lower canopy
    !  dropu   = drought-phenology trigger for trees
    !  dropls  = drought-phenology trigger for shrubs
    !  dropl4  = drought-phenology trigger for c4 grasses
    !  dropl3  = drought-phenology trigger for c3 grasses
    ! (NOTE: a10ancuc is not used at this point, so its restart entry 
    ! is commented out)
    !
    !
    ! Arguments
    !
    INTEGER, INTENT(IN   ) :: npft           
    REAL(KIND=r8)   , INTENT(OUT  ) :: fi           (ibMax,jbMax)          ! fraction of woody biomass that is in sapwood       
    REAL(KIND=r8)   , INTENT(OUT  ) :: tsno         (ibMax,nsnolay,jbMax)! temperature of snow layers (K)
    REAL(KIND=r8)   , INTENT(OUT  ) :: hsno         (ibMax,nsnolay,jbMax)! thickness of snow layers (m)
    REAL(KIND=r8)   , INTENT(OUT  ) :: tsoi         (ibMax,nsoilay,jbMax)! soil temperature for each layer (K)
    REAL(KIND=r8)   , INTENT(OUT  ) :: tsoim        (ibMax,nsoilay,jbMax)! soil temperature for each layer (K)
    REAL(KIND=r8)   , INTENT(OUT  ) :: wsoi         (ibMax,nsoilay,jbMax)! fraction of soil pore space containing liquid water
    REAL(KIND=r8)   , INTENT(OUT  ) :: wsoim        (ibMax,nsoilay,jbMax)! fraction of soil pore space containing liquid water
    REAL(KIND=r8)   , INTENT(OUT  ) :: wisoi        (ibMax,nsoilay,jbMax)! fraction of soil pore space containing ice
    REAL(KIND=r8)   , INTENT(OUT  ) :: cbiol        (ibMax,npft,jbMax)   ! carbon in leaf biomass pool (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: cbiow        (ibMax,npft,jbMax)   ! carbon in woody biomass pool (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: cbior        (ibMax,npft,jbMax)   ! carbon in fine root biomass pool (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: sapfrac      (ibMax,jbMax)             ! fraction of woody biomass that is in sapwood
    REAL(KIND=r8)   , INTENT(OUT  ) :: clitlm       (ibMax,jbMax)! carbon in leaf litter pool - metabolic         (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: clitls       (ibMax,jbMax)! carbon in leaf litter pool - structural         (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: clitll       (ibMax,jbMax)! carbon in leaf litter pool - lignin                (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: clitrm       (ibMax,jbMax)! carbon in fine root litter pool - metabolic  (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: clitrs       (ibMax,jbMax)! carbon in fine root litter pool - structural (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: clitrl       (ibMax,jbMax)! carbon in fine root litter pool - lignin         (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: clitwm       (ibMax,jbMax)! carbon in woody litter pool - metabolic         (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: clitws       (ibMax,jbMax)! carbon in woody litter pool - structural         (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: clitwl       (ibMax,jbMax)! carbon in woody litter pool - lignin            (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: falll        (ibMax,jbMax)! annual leaf litter fall (kg_C m-2/year)
    REAL(KIND=r8)   , INTENT(OUT  ) :: fallr        (ibMax,jbMax)! annual root litter input                         (kg_C m-2/year)
    REAL(KIND=r8)   , INTENT(OUT  ) :: fallw        (ibMax,jbMax)! annual wood litter fall                         (kg_C m-2/year)
    REAL(KIND=r8)   , INTENT(OUT  ) :: totcmic      (ibMax,jbMax)! total carbon residing in microbial pools (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: csoislop     (ibMax,jbMax)! carbon in soil - slow protected humus      (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: csoislon     (ibMax,jbMax)! carbon in soil - slow nonprotected humus         (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: csoipas      (ibMax,jbMax)! carbon in soil - passive humus                (kg_C m-2)
    REAL(KIND=r8)   , INTENT(INOUT) :: gdd0         (ibMax,jbMax)    
    REAL(KIND=r8)   , INTENT(INOUT) :: gdd5         (ibMax,jbMax)    
    REAL(KIND=r8)   , INTENT(OUT  ) :: tc           (ibMax,jbMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: tw           (ibMax,jbMax)    
    REAL(KIND=r8)   , INTENT(OUT  ) :: wipud        (ibMax,jbMax) ! ice content of puddles per soil area (kg m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: wpud         (ibMax,jbMax) ! liquid content of puddles per soil area (kg m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: agddu        (ibMax,jbMax)  ! annual accumulated growing degree days for bud burst,
    ! upper canopy (day-degrees)
    REAL(KIND=r8)   , INTENT(OUT  ) :: agddl        (ibMax,jbMax)! annual accumulated growing degree days for bud burst, 
    !lower canopy (day-degrees)
    REAL(KIND=r8)   , INTENT(OUT  ) :: tempu        (ibMax,jbMax)! cold-phenology trigger for trees (non-dimensional)
    REAL(KIND=r8)   , INTENT(OUT  ) :: templ        (ibMax,jbMax) ! cold-phenology trigger for grasses/shrubs (non-dimensional)
    REAL(KIND=r8)   , INTENT(OUT  ) :: a10td        (ibMax,jbMax) ! 10-day average daily air temperature (K)
    REAL(KIND=r8)   , INTENT(OUT  ) :: a10ancub     (ibMax,jbMax) ! 10-day average canopy photosynthesis rate 
    ! - broadleaf (mol_co2 m-2 s-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: a10ancls     (ibMax,jbMax) ! 10-day average canopy photosynthesis rate - 
    ! shrubs (mol_co2 m-2 s-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: a10ancl4     (ibMax,jbMax) ! 10-day average canopy photosynthesis rate - 
    ! c4 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: a10ancl3     (ibMax,jbMax) ! 10-day average canopy photosynthesis rate - 
    ! c3 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: a10scalparamu(ibMax,jbMax) ! 10-day average day-time scaling parameter - 
    ! upper canopy (dimensionless)
    REAL(KIND=r8)   , INTENT(OUT  ) :: a10scalparaml(ibMax,jbMax) ! 10-day average day-time scaling parameter - 
    ! lower canopy (dimensionless)
    REAL(KIND=r8)   , INTENT(OUT  ) :: a10daylightu (ibMax,jbMax) ! 10-day average day-time PAR - 
    ! upper canopy (micro-Ein m-2 s-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: a10daylightl (ibMax,jbMax) ! 10-day average day-time PAR - 
    ! lower canopy (micro-Ein m-2 s-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dropu        (ibMax,jbMax) ! drought-phenology trigger for trees (non-dimensional)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dropls       (ibMax,jbMax) ! drought-phenology trigger for shrubs (non-dimensional)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dropl4       (ibMax,jbMax) ! drought-phenology trigger for c4 grasses (non-dimensional)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dropl3       (ibMax,jbMax) ! drought-phenology trigger for c3 grasses (non-dimensional)
    REAL(KIND=r8)   , INTENT(INOUT) :: tcmin        (ibMax,jbMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: deltat       (ibMax,jbMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: exist        (ibMax,npft,jbMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: TminL        (npft)          
    REAL(KIND=r8)   , INTENT(IN   ) :: TminU        (npft)          
    REAL(KIND=r8)   , INTENT(IN   ) :: Twarm        (npft)          
    REAL(KIND=r8)   , INTENT(IN   ) :: GDD          (npft) 
    REAL(KIND=r8)   , INTENT(OUT  ) :: tcthis       (ibMax,jbMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: twthis       (ibMax,jbMax)    
    !
    ! Local variables
    !
    !
    !
    INTEGER  :: i,j,nLndPts
    !
    ! ---------------------------------------------------------------------
    !
    CALL MsgOne('**(ReStartIBIS)**','Reading physics state for restart')

    READ(UNIT=nfsibi)asurd   
    READ(UNIT=nfsibi)asuri   
    READ(UNIT=nfsibi)totcondub
    READ(UNIT=nfsibi)totconduc
    READ(UNIT=nfsibi)totcondls
    READ(UNIT=nfsibi)totcondl3
    READ(UNIT=nfsibi)totcondl4
    READ(UNIT=nfsibi)ginvap  
    READ(UNIT=nfsibi)gsuvap  
    READ(UNIT=nfsibi)gtrans  
    READ(UNIT=nfsibi)grunof  
    READ(UNIT=nfsibi)gdrain  
    READ(UNIT=nfsibi)totlit  
    READ(UNIT=nfsibi)totnlit 
    READ(UNIT=nfsibi)totfall 
    READ(UNIT=nfsibi)totalit 
    READ(UNIT=nfsibi)totrlit 
    READ(UNIT=nfsibi)totanlit
    READ(UNIT=nfsibi)totrnlit
    READ(UNIT=nfsibi)totcsoi 
    READ(UNIT=nfsibi)totnmic 
    READ(UNIT=nfsibi)tco2mic 
    READ(UNIT=nfsibi)totnsoi
    READ(UNIT=nfsibi)tnpptot 
    READ(UNIT=nfsibi)tneetot 
    READ(UNIT=nfsibi)tnmin   
    READ(UNIT=nfsibi)cdisturb
    READ(UNIT=nfsibi)su      
    READ(UNIT=nfsibi)ss      
    READ(UNIT=nfsibi)sl      
    READ(UNIT=nfsibi)stresstu
    READ(UNIT=nfsibi)stresstl
    READ(UNIT=nfsibi)stressl 
    READ(UNIT=nfsibi)stressu 
    READ(UNIT=nfsibi)tcthis
    READ(UNIT=nfsibi)twthis
    !
    !  daily average variables 
    !
    READ(UNIT=nfsibi)adrain          
    READ(UNIT=nfsibi)adsnow          
    READ(UNIT=nfsibi)adaet          
    READ(UNIT=nfsibi)adtrunoff 
    READ(UNIT=nfsibi)adsrunoff 
    READ(UNIT=nfsibi)addrainage
    READ(UNIT=nfsibi)adrh          
    READ(UNIT=nfsibi)adsnod          
    READ(UNIT=nfsibi)adsnof          
    READ(UNIT=nfsibi)adwsoi          
    READ(UNIT=nfsibi)adtsoi          
    READ(UNIT=nfsibi)adwisoi    
    READ(UNIT=nfsibi)adtlaysoi 
    READ(UNIT=nfsibi)adwlaysoi 
    READ(UNIT=nfsibi)adwsoic    
    READ(UNIT=nfsibi)adtsoic    
    READ(UNIT=nfsibi)adco2mic  
    READ(UNIT=nfsibi)adco2root 
    READ(UNIT=nfsibi)adco2soi  
    READ(UNIT=nfsibi)adco2ratio
    READ(UNIT=nfsibi)adnmintot 
    READ(UNIT=nfsibi)decompl
    READ(UNIT=nfsibi)decomps
    READ(UNIT=nfsibi)gdd0this
    READ(UNIT=nfsibi)gdd5this
    READ(UNIT=nfsibi)storedn
    READ(UNIT=nfsibi)yrleach
    READ(UNIT=nfsibi)ynleach
    !
    !  monthly average variables 
    !
    READ(UNIT=nfsibi)amrain    
    READ(UNIT=nfsibi)amsnow    
    READ(UNIT=nfsibi)amaet          
    READ(UNIT=nfsibi)amtrunoff 
    READ(UNIT=nfsibi)amsrunoff 
    READ(UNIT=nfsibi)amdrainage
    READ(UNIT=nfsibi)amtemp    
    READ(UNIT=nfsibi)amqa          
    READ(UNIT=nfsibi)amsolar   
    READ(UNIT=nfsibi)amirup    
    READ(UNIT=nfsibi)amirdown  
    READ(UNIT=nfsibi)amsens    
    READ(UNIT=nfsibi)amlatent  
    READ(UNIT=nfsibi)amlaiu    
    READ(UNIT=nfsibi)amlail    
    READ(UNIT=nfsibi)amtsoi    
    READ(UNIT=nfsibi)amwsoi    
    READ(UNIT=nfsibi)amwisoi   
    READ(UNIT=nfsibi)amvwc          
    READ(UNIT=nfsibi)amawc          
    READ(UNIT=nfsibi)amsnod    
    READ(UNIT=nfsibi)amsnof    
    READ(UNIT=nfsibi)amnpp          
    READ(UNIT=nfsibi)amnpptot  
    READ(UNIT=nfsibi)amco2mic  
    READ(UNIT=nfsibi)amco2root 
    READ(UNIT=nfsibi)amco2soi            
    READ(UNIT=nfsibi)amco2ratio
    READ(UNIT=nfsibi)amneetot  
    READ(UNIT=nfsibi)amnmintot 
    READ(UNIT=nfsibi)amts2          
    READ(UNIT=nfsibi)amtransu  
    READ(UNIT=nfsibi)amtransl  
    READ(UNIT=nfsibi)amsuvap   
    READ(UNIT=nfsibi)aminvap   
    READ(UNIT=nfsibi)amalbedo  
    READ(UNIT=nfsibi)amtsoil   
    READ(UNIT=nfsibi)amwsoil   
    READ(UNIT=nfsibi)amwisoil  
    !
    !  annual total variables 
    !
    READ(UNIT=nfsibi)aysolar   
    READ(UNIT=nfsibi)ayirup    
    READ(UNIT=nfsibi)ayirdown  
    READ(UNIT=nfsibi)aysens    
    READ(UNIT=nfsibi)aylatent  
    READ(UNIT=nfsibi)ayprcp    
    READ(UNIT=nfsibi)ayaet          
    READ(UNIT=nfsibi)aytrans   
    READ(UNIT=nfsibi)aytrunoff 
    READ(UNIT=nfsibi)aysrunoff 
    READ(UNIT=nfsibi)aydrainage
    READ(UNIT=nfsibi)aydwtot   
    READ(UNIT=nfsibi)aywsoi    
    READ(UNIT=nfsibi)aywisoi   
    READ(UNIT=nfsibi)aytsoi    
    READ(UNIT=nfsibi)ayvwc          
    READ(UNIT=nfsibi)ayawc          
    READ(UNIT=nfsibi)aystresstu
    READ(UNIT=nfsibi)aystresstl
    READ(UNIT=nfsibi)aygpp          
    READ(UNIT=nfsibi)aygpptot  
    READ(UNIT=nfsibi)aynpp          
    READ(UNIT=nfsibi)aynpptot  
    READ(UNIT=nfsibi)ayco2mic  
    READ(UNIT=nfsibi)ayco2root 
    READ(UNIT=nfsibi)ayco2soi  
    READ(UNIT=nfsibi)ayneetot  
    READ(UNIT=nfsibi)ayrootbio 
    READ(UNIT=nfsibi)aynmintot 
    READ(UNIT=nfsibi)ayalit    
    READ(UNIT=nfsibi)ayblit    
    READ(UNIT=nfsibi)aycsoi    
    READ(UNIT=nfsibi)aycmic    
    READ(UNIT=nfsibi)ayanlit   
    READ(UNIT=nfsibi)aybnlit   
    READ(UNIT=nfsibi)aynsoi    
    READ(UNIT=nfsibi)ayalbedo  
    READ(UNIT=nfsibi)firefac      
    READ(UNIT=nfsibi)wtot 

    READ(UNIT=nfsibi)tlsub
    READ(UNIT=nfsibi)t12  
    READ(UNIT=nfsibi)t34  
    READ(UNIT=nfsibi)q12  
    READ(UNIT=nfsibi)q34  
    READ(UNIT=nfsibi)ciub 
    READ(UNIT=nfsibi)ciuc 
    READ(UNIT=nfsibi)cils 
    READ(UNIT=nfsibi)cil3 
    READ(UNIT=nfsibi)cil4 
    READ(UNIT=nfsibi)csub 
    READ(UNIT=nfsibi)csuc 
    READ(UNIT=nfsibi)csls 
    READ(UNIT=nfsibi)csl3 
    READ(UNIT=nfsibi)csl4 
    READ(UNIT=nfsibi)gsub 
    READ(UNIT=nfsibi)gsuc 
    READ(UNIT=nfsibi)gsls 
    READ(UNIT=nfsibi)gsl3 
    READ(UNIT=nfsibi)gsl4 
    READ(UNIT=nfsibi)wliqu
    READ(UNIT=nfsibi)wliqs
    READ(UNIT=nfsibi)wliql
    READ(UNIT=nfsibi)wsnou
    READ(UNIT=nfsibi)wsnos
    READ(UNIT=nfsibi)wsnol

    READ(UNIT=nfsibi)fi,fu,fl,tu,ts,tl,tg,ti

    !
    ! nsnolay variables: tsno and hsno
    !

    !
    READ(UNIT=nfsibi)tsno

    READ(UNIT=nfsibi)hsno

    !
    ! nsoilay variables: tsoi, wisoi, wsoi
    !

    READ(UNIT=nfsibi)tsoim,tsoi

    READ(UNIT=nfsibi)wisoi

    READ(UNIT=nfsibi)wsoim,wsoi

    !
    ! npft variables
    !

    READ(UNIT=nfsibi)cbiol

    READ(UNIT=nfsibi)cbiow

    READ(UNIT=nfsibi)cbior
    READ(UNIT=nfsibi)ndtimes
    READ(UNIT=nfsibi)nmtimes
    READ(UNIT=nfsibi)nytimes
    !
    ! single level variables
    !

    READ(UNIT=nfsibi)sapfrac

    READ(UNIT=nfsibi)clitlm

    READ(UNIT=nfsibi)clitls

    READ(UNIT=nfsibi)clitll

    READ(UNIT=nfsibi)clitrm

    READ(UNIT=nfsibi)clitrs

    READ(UNIT=nfsibi)clitrl

    READ(UNIT=nfsibi)clitwm

    READ(UNIT=nfsibi)clitws

    READ(UNIT=nfsibi)clitwl

    READ(UNIT=nfsibi)falll

    READ(UNIT=nfsibi)fallr

    READ(UNIT=nfsibi)fallw

    READ(UNIT=nfsibi)totcmic

    READ(UNIT=nfsibi)csoislop

    READ(UNIT=nfsibi)csoislon

    READ(UNIT=nfsibi)csoipas

    READ(UNIT=nfsibi)gdd0

    READ(UNIT=nfsibi)gdd5

    READ(UNIT=nfsibi)tc

    READ(UNIT=nfsibi)tw

    READ(UNIT=nfsibi)wipud

    READ(UNIT=nfsibi)wpud

    READ(UNIT=nfsibi)agddu

    READ(UNIT=nfsibi)agddl

    READ(UNIT=nfsibi)tempu

    READ(UNIT=nfsibi)templ

    READ(UNIT=nfsibi)a10td

    READ(UNIT=nfsibi)a10ancub

    READ(UNIT=nfsibi)a10ancls

    READ(UNIT=nfsibi)a10ancl4

    READ(UNIT=nfsibi)a10ancl3

    READ(UNIT=nfsibi)a10scalparamu

    READ(UNIT=nfsibi)a10scalparaml

    READ(UNIT=nfsibi)a10daylightu

    READ(UNIT=nfsibi)a10daylightl

    READ(UNIT=nfsibi)dropu

    READ(UNIT=nfsibi)dropls

    READ(UNIT=nfsibi)dropl4

    READ(UNIT=nfsibi)dropl3

    READ(UNIT=nfsibi)lai

    READ(UNIT=nfsibi)zbot

    READ(UNIT=nfsibi)ztop 

    READ(UNIT=nfsibi)frac     

    READ(UNIT=nfsibi)td  

    READ(UNIT=nfsibi)  ppli,ppci

    READ(UNIT=nfsibi) gl0 ,zorl,gtsea,tseam,qsfc0,tsfc0,qsfcm,tsfcm,tkemyj,HML,HUML,HVML,TSK,z0sea

    READ(UNIT=nfsibi) w0,wm,capac0,capacm,td0,tdm,tcm,tc0,tgm,tg0,z0

    READ(UNIT=nfsibi) idateprev

    READ(UNIT=nfsibi) sai

    READ(UNIT=nfsibi) plai

    READ(UNIT=nfsibi) biomass

    READ(UNIT=nfsibi) totlaiu

    READ(UNIT=nfsibi) totlail

    READ(UNIT=nfsibi) totbiou

    READ(UNIT=nfsibi) totbiol

    READ(UNIT=nfsibi) exist

    READ(UNIT=nfsibi) vegtype0
    Mmlen=gl0

    RETURN
    !
    ! calculate tcmin
    !
    DO j = 1,jbMax
       nLndPts=0
       DO i = 1, ibMax
          IF (iMaskIBIS(i,j) >= 1_i8) THEN
             nLndPts=nLndPts+1
             tcmin(nLndPts,j) = tc(nLndPts,j) + deltat(nLndPts,j)
          END IF
       END DO
    END DO
    !
    CALL existence(TminL, &
         TminU, &
         Twarm, &
         GDD  , &
         exist, &
         tcmin, &
         gdd5 , & 
         gdd0 , &
         tw   , &
         npft   )  
    !
    RETURN
  END  SUBROUTINE restart
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE existence(TminL    , &! INTENT(IN   )
       TminU    , &! INTENT(IN   )
       Twarm    , &! INTENT(IN   )
       GDD      , &! INTENT(IN   )
       exist    , &! INTENT(OUT  )
       tcmin    , &! INTENT(IN   )
       gdd5     , &! INTENT(IN   ) 
       gdd0     , &! INTENT(IN   ) 
       tw       , &! INTENT(IN   )
       npft       )! INTENT(IN   )
    ! ---------------------------------------------------------------------
    !
    ! this routine determines which plant functional types (pft's) are allowed
    ! to exist in each gridcell, based on a simple set of climatic criteria
    !
    ! the logic here is based on the biome3 model of haxeltine and prentice
    !
    ! plant functional types:
    !
    ! 1)  tropical broadleaf evergreen trees
    ! 2)  tropical broadleaf drought-deciduous trees
    ! 3)  warm-temperate broadleaf evergreen trees
    ! 4)  temperate conifer evergreen trees
    ! 5)  temperate broadleaf cold-deciduous trees
    ! 6)  boREAL(KIND=r8) conifer evergreen trees
    ! 7)  boREAL(KIND=r8) broadleaf cold-deciduous trees
    ! 8)  boREAL(KIND=r8) conifer cold-deciduous trees
    ! 9)  evergreen shrubs
    ! 10) deciduous shrubs
    ! 11) warm (c4) grasses
    ! 12) cool (c3) grasses
    !
    !
    ! common blocks
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN   ) :: npft            ! number of plant functional types

    REAL(KIND=r8)   , INTENT(OUT  ) :: exist(ibMax,npft,jbMax)! probability of existence of each plant functional type in a gridcell
    REAL(KIND=r8)   , INTENT(IN   ) :: tcmin(ibMax,jbMax)     ! coldest daily temperature of current year (C)
    REAL(KIND=r8)   , INTENT(IN   ) :: gdd5 (ibMax,jbMax)     ! growing degree days > 5C
    REAL(KIND=r8)   , INTENT(IN   ) :: gdd0 (ibMax,jbMax)     ! growing degree days > 0C 
    REAL(KIND=r8)   , INTENT(IN   ) :: tw   (ibMax,jbMax)     ! warmest monthly temperature (C)
    REAL(KIND=r8)   , INTENT(IN   ) :: TminL(npft)     ! Absolute minimum temperature -- lower limit (upper canopy PFTs)
    REAL(KIND=r8)   , INTENT(IN   ) :: TminU(npft)     ! Absolute minimum temperature -- upper limit (upper canopy PFTs)
    REAL(KIND=r8)   , INTENT(IN   ) :: Twarm(npft)     ! Temperature of warmest month (lower canopy PFTs)
    REAL(KIND=r8)   , INTENT(IN   ) :: GDD  (npft)     ! minimum GDD needed (base 5 C for upper canopy PFTs, 
    ! base 0 C for lower canopy PFTs)

    !
    ! Local variables
    !
    INTEGER :: i   ,j  ,nLndPts ! loop indice
    !
    ! ---------------------------------------------------------------------
    !
    DO j = 1, jbMax
       nLndPts=0
       DO i = 1, ibMax
          IF (iMaskIBIS(i,j) >= 1_i8) THEN
             nLndPts=nLndPts+1
             !
             ! determine which plant types can exist in a given gridcell
             !
             exist(nLndPts,1,j)  = 0.0_r8
             exist(nLndPts,2,j)  = 0.0_r8
             exist(nLndPts,3,j)  = 0.0_r8
             exist(nLndPts,4,j)  = 0.0_r8
             exist(nLndPts,5,j)  = 0.0_r8
             exist(nLndPts,6,j)  = 0.0_r8
             exist(nLndPts,7,j)  = 0.0_r8
             exist(nLndPts,8,j)  = 0.0_r8
             exist(nLndPts,9,j)  = 0.0_r8
             exist(nLndPts,10,j) = 0.0_r8
             exist(nLndPts,11,j) = 0.0_r8
             exist(nLndPts,12,j) = 0.0_r8
             !
             ! 1) tropical broadleaf evergreen trees
             !
             !  - tcmin > 0.0
             !
             !        IF (tcmin(i).gt.0.0_r8)           exist(i,1) = 1.0_r8
             !
             ! 2) tropical broadleaf drought-deciduous trees
             !
             !  - tcmin > 0.0_r8
             !
             !        IF (tcmin(i).gt.0.0_r8)           exist(i,2) = 1.0_r8
             !
             ! 3) warm-temperate broadleaf evergreen trees
             !
             !  - tcmin <   0.0_r8 and
             !  - tcmin > -10.0_r8
             !
             !        IF ((tcmin(i).lt.0.0_r8).and. &
             !           (tcmin(i).gt.-10.0_r8))       exist(i,3) = 1.0_r8
             !
             ! 4) temperate conifer evergreen trees
             !
             !  - tcmin <    0.0_r8 and
             !  - tcmin >  -45.0_r8 and
             !  - gdd5  > 1200.0_r8
             !
             !        IF ((tcmin(i).lt.0.0_r8).and. &
             !           (tcmin(i).gt.-45.0_r8).and. &
             !           (gdd5(i).gt.1200.0_r8))       exist(i,4) = 1.0_r8
             !
             ! 5) temperate broadleaf cold-deciduous trees
             !
             !  - tcmin <    0.0 and
             !  - tcmin >  -45.0 and
             !  - gdd5  > 1200.0
             !
             !        IF ((tcmin(i).lt.0.0_r8).and.   &
             !           (tcmin(i).gt.-45.0_r8).and.  &
             !           (gdd5(i).gt.1200.0_r8))       exist(i,5) = 1.0_r8
             !
             ! 6) boREAL(KIND=r8) conifer evergreen trees
             !
             !  - tcmin <  -45.0_r8 or gdd5 < 1200.0_r8, and
             !  - tcmin >  -57.5_r8 and
             !  - gdd5  >  350.0_r8
             !
             !        IF (((tcmin(i).lt.-45.0_r8).or.(gdd5(i).lt.1200.0_r8)).and. &
             !            (tcmin(i).gt.-57.5_r8).and. &
             !            (gdd5(i).gt.350.0_r8))       exist(i,6) = 1.0_r8
             !
             ! 7) boREAL(KIND=r8) broadleaf cold-deciduous trees
             !
             !  - tcmin <  -45.0 or gdd5 < 1200.0, and
             !  - tcmin >  -57.5 and
             !  - gdd5  >  350.0
             !
             !        IF (((tcmin(i).lt.-45.0_r8).or.(gdd5(i).lt.1200.0_r8)).and. &
             !            (tcmin(i).gt.-57.5_r8).and. &
             !            (gdd5(i).gt.350.0_r8))       exist(i,7) = 1.0_r8
             !
             ! 8) boREAL(KIND=r8) conifer cold-deciduous trees
             !
             !  - tcmin <  -45.0 or gdd5 < 1200.0, and
             !  - gdd5  >  350.0
             !
             !        IF (((tcmin(i).lt.-45.0_r8).or.(gdd5(i).lt.1200.0_r8)).and. &
             !            (gdd5(i).gt.350.0_r8))       exist(i,8) = 1.0_r8
             !
             ! 9) evergreen shrubs
             !
             !  - gdd0 > 100.0
             !
             !        IF (gdd0(i).gt.100.0_r8)          exist(i,9) = 1.0_r8
             !
             ! 10) deciduous shrubs
             !
             !  - gdd0 > 100.0
             !
             !        IF (gdd0(i).gt.100.0_r8)          exist(i,10) = 1.0_r8
             !
             ! 11) warm (c4) grasses
             !
             !  - tw   >  22.0 and
             !  - gdd0 > 100.0
             !
             !        IF ((tw(i).gt.22.0_r8).and. &
             !            (gdd0(i).gt.100.0_r8))        exist(i,11) = 1.0_r8
             !
             ! 12) cool (c3) grasses
             !
             !  - gdd0 > 100.0
             !
             !        IF (gdd0(i).gt.100.0_r8)          exist(i,12) = 1.0_r8
             !
             !
             !*** DTP 2001/06/07: Modified version of above code reads in PFT
             !    existence criteria from external parameter file "params.veg"
             !    These are copied here for reference.... 
             !------------------------------------------------------------------
             !  TminL    TminU    Twarm    GDD    PFT
             !------------------------------------------------------------------
             !    0.0   9999.0   9999.0   9999  !   1
             !    0.0   9999.0   9999.0   9999  !   2
             !  -10.0      0.0   9999.0   9999  !   3
             !  -45.0      0.0   9999.0   1200  !   4
             !  -45.0      0.0   9999.0   1200  !   5
             !  -57.5    -45.0   9999.0    350  !   6
             !  -57.5    -45.0   9999.0    350  !   7
             ! 9999.0    -45.0   9999.0    350  !   8
             ! 9999.0   9999.0   9999.0    100  !   9
             ! 9999.0   9999.0   9999.0    100  !  10
             ! 9999.0   9999.0     22.0    100  !  11
             ! 9999.0   9999.0   9999.0    100  !  12
             !------------------------------------------------------------------

             ! 1) tropical broadleaf evergreen trees
             !
             !  - tcmin > 0.0
             !
             IF (tcmin(nLndPts,j).GT.TminL(1))      exist(nLndPts,1,j) = 1.0_r8
             !
             ! 2) tropical broadleaf drought-deciduous trees
             !
             !  - tcmin > 0.0
             !
             IF (tcmin(nLndPts,j).GT.TminL(2))      exist(nLndPts,2,j) = 1.0_r8
             !
             ! 3) warm-temperate broadleaf evergreen trees
             !
             !  - tcmin <   0.0 and
             !  - tcmin > -10.0
             !
             IF ((tcmin(nLndPts,j).LT.TminU(3)).AND.  &
                  (tcmin(nLndPts,j).GT.TminL(3)))    exist(nLndPts,3,j) = 1.0_r8
             !
             ! 4) temperate conifer evergreen trees
             !
             !  - tcmin <    0.0 and
             !  - tcmin >  -45.0 and
             !  - gdd5  > 1200.0
             !
             IF ((tcmin(nLndPts,j).LT.TminU(4)).AND.   &
                  (tcmin(nLndPts,j).GT.TminL(4)).AND.   &
                  (gdd5(nLndPts,j).GT.GDD(4)))       exist(nLndPts,4,j) = 1.0_r8
             !
             ! 5) temperate broadleaf cold-deciduous trees
             !
             !  - tcmin <    0.0 and
             !  - tcmin >  -45.0 and
             !  - gdd5  > 1200.0
             !
             IF ((tcmin(nLndPts,j).LT.TminU(5)).AND.     &
                  (tcmin(nLndPts,j).GT.TminL(5)).AND.     & 
                  (gdd5(nLndPts,j).GT.GDD(5)))       exist(nLndPts,5,j) = 1.0_r8
             !
             ! 6) boreal conifer evergreen trees
             !
             !  - tcmin <  -45.0 or gdd5 < 1200.0, and
             !  - tcmin >  -57.5 and
             !  - gdd5  >  350.0
             !
             IF (((tcmin(nLndPts,j).LT.TminU(6)).OR.   &
                  (gdd5(nLndPts,j).LT.GDD(4))).AND.     &
                  (tcmin(nLndPts,j).GT.TminL(6)).AND.   &
                  (gdd5(nLndPts,j).GT.GDD(6)))       exist(nLndPts,6,j) = 1.0_r8
             !
             ! 7) boreal broadleaf cold-deciduous trees
             !
             !  - tcmin <  -45.0 or gdd5 < 1200.0, and
             !  - tcmin >  -57.5 and
             !  - gdd5  >  350.0
             !
             IF (((tcmin(nLndPts,j).LT.TminU(7)).OR.  &
                  (gdd5(nLndPts,j).LT.GDD(5))).AND.    & 
                  (tcmin(nLndPts,j).GT.TminL(7)).AND.  &
                  (gdd5(nLndPts,j).GT.GDD(7)))       exist(nLndPts,7,j) = 1.0_r8
             !
             ! 8) boreal conifer cold-deciduous trees
             !
             !  - tcmin <  -45.0 or gdd5 < 1200.0, and
             !  - gdd5  >  350.0
             !
             IF (((tcmin(nLndPts,j).LT.TminU(8)).OR.  &
                  (gdd5(nLndPts,j).LT.TminL(4))).AND.  &
                  (gdd5(nLndPts,j).GT.GDD(8)))       exist(nLndPts,8,j) = 1.0_r8
             !
             ! 9) evergreen shrubs
             !
             !  - gdd0 > 100.0
             !
             IF (gdd0(nLndPts,j).GT.GDD(9))         exist(nLndPts,9,j) = 1.0_r8
             !
             ! 10) deciduous shrubs
             !
             !  - gdd0 > 100.0
             !
             IF (gdd0(nLndPts,j).GT.GDD(10))        exist(nLndPts,10,j) = 1.0_r8
             !
             ! 11) warm (c4) grasses
             !
             !  - tw   >  22.0 and
             !  - gdd0 > 100.0
             !
             IF ((tw(nLndPts,j).GT.Twarm(11)).AND.  &
                  (gdd0(nLndPts,j).GT.GDD(11)))      exist(nLndPts,11,j) = 1.0_r8
             !
             ! 12) cool (c3) grasses
             !
             !  - gdd0 > 100.0
             !
             IF (gdd0(nLndPts,j).GT.GDD(12))        exist(nLndPts,12,j) = 1.0_r8

          END IF
       END DO
    END DO
    !
    RETURN
  END SUBROUTINE existence

  ! ---------------------------------------------------------------------
  SUBROUTINE coldstart(iMax   , &
       jMax   , &
       kMax   , &
       ibMax  , &
       jbMax  , &
       ifday  , &
       ibMaxPerJB, &
       tod    , &
       idate  , &
       idatec , &! INTENT(IN   )
       nsoilay, &! INTENT(IN   )
       nsnolay, &! INTENT(IN   )
       fgtmp  , &! INTENT(IN        )
       fgq    , &! INTENT(IN        )
       hsno   , &! INTENT(OUT  )
       tsno   , &! INTENT(OUT  )
       tsoi   , &! INTENT(OUT  )
       tsoim  , &! INTENT(OUT  )
       wsoi   , &! INTENT(OUT  )
       wsoim  , &! INTENT(OUT  )
       wisoi    )! INTENT(OUT  )
    ! ---------------------------------------------------------------------
    !  
    IMPLICIT NONE
    !
    INTEGER         , INTENT(IN   ) :: iMax  
    INTEGER         , INTENT(IN   ) :: jMax  
    INTEGER         , INTENT(IN   ) :: kMax  
    INTEGER         , INTENT(IN   ) :: ibMax 
    INTEGER         , INTENT(IN   ) :: jbMax 
    INTEGER         , INTENT(IN   ) :: ifday 
    INTEGER         , INTENT(IN   ) :: ibMaxPerJB(jbMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: tod    
    INTEGER         , INTENT(IN   ) :: idate(4)  
    INTEGER         , INTENT(IN   ) :: idatec(4) ! INTENT(IN  )
    INTEGER         , INTENT(IN   ) :: nsoilay    ! number of soil layers
    INTEGER         , INTENT(IN   ) :: nsnolay    ! number of snow layers
    REAL(KIND=r8)   , INTENT(IN   ) :: fgtmp  (ibMax,kMax,jbMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: fgq    (ibMax,kMax,jbMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: hsno   (ibMax,nsnolay,jbMax)! thickness of snow layers (m)
    REAL(KIND=r8)   , INTENT(OUT  ) :: tsno   (ibMax,nsnolay,jbMax)! temperature of snow layers (K)
    REAL(KIND=r8)   , INTENT(OUT  ) :: tsoi   (ibMax,nsoilay,jbMax)! soil temperature for each layer (K)
    REAL(KIND=r8)   , INTENT(OUT  ) :: tsoim   (ibMax,nsoilay,jbMax)! soil temperature for each layer (K)
    REAL(KIND=r8)   , INTENT(OUT  ) :: wsoi   (ibMax,nsoilay,jbMax)! fraction of soil pore space containing liquid water
    REAL(KIND=r8)   , INTENT(OUT  ) :: wsoim   (ibMax,nsoilay,jbMax)! fraction of soil pore space containing liquid water
    REAL(KIND=r8)   , INTENT(OUT  ) :: wisoi  (ibMax,nsoilay,jbMax)! fraction of soil pore space containing ice
    REAL(KIND=r8)                :: buf      (iMax,jMax,4)
    REAL(KIND=r4) ::   brf (iMax,jMax)
    INTEGER :: i,j,k , nLndPts,irec,ierr,LRecIN,ncount
    REAL(KIND=r8)                :: sinmax
    REAL(KIND=r8), PARAMETER     :: xl0   =10.0_r8
    REAL(KIND=r8), PARAMETER     :: t0 =271.17_r8
    REAL(KIND=r8), PARAMETER     :: zero  =0.0e3_r8
    REAL(KIND=r8), PARAMETER     :: thousd=1.0e3_r8


    CLOSE(nftgz0)
    brf=0.0_r4
    INQUIRE (IOLENGTH=LRecIN) brf
    OPEN (UNIT=nftgz0,FILE=TRIM(fNameTg3zrl), FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIN, &
         ACTION='read', STATUS='OLD', IOSTAT=ierr) 
    IF (ierr /= 0) THEN
       WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(fNameTg3zrl), ierr
       STOP "**(ERROR)**"
    END IF
    CLOSE(nfzol)
    buf=0.0_r8
    INQUIRE (IOLENGTH=LRecIN) brf
    OPEN (UNIT=nfzol,FILE=TRIM(fNameRouLen),FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIN, &
         ACTION='READ', STATUS='OLD', IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(fNameRouLen), ierr
       STOP "**(ERROR)**"
    END IF

    CALL getsbc (iMax ,jMax  ,kMax, AlbVisDiff,gtsea,gndvi,soilm,sheleg,o3mix,tracermix,wsib3d,&
         ifday , tod  ,idate ,idatec, &
         ifalb,ifsst,ifndvi,ifslm ,ifslmSib2,ifsnw,ifozone,iftracer, &
         sstlag,intsst,intndvi,fint ,tice  , &
         yrl22  ,monl,ibMax,jbMax,ibMaxPerJB)

    irec=1
    CALL ReadGetNFTGZ(nftgz0,irec,buf(1:iMax,1:jMax,1),buf(1:iMax,1:jMax,2),buf(1:iMax,1:jMax,3))
    READ (UNIT=nfzol, REC=1) brf

    buf(1:iMax,1:jMax,4)=brf(1:iMax,1:jMax)

    IF (reducedGrid) THEN
             CALL LinearIJtoIBJB(buf(1:iMax,1:jMax,1),tg1)
            !CALL AveBoxIJtoIBJB(buf(1:iMax,1:jMax,1),tg1)
    ELSE
       CALL IJtoIBJB(buf(1:iMax,1:jMax,1) ,tg1 )
    END IF
    !tg1(1:ibMax,1:jbMax)=MAX(fgtmp(1:ibMax,1,1:jbMax),tg1(1:ibMax,1:jbMax))
    IF (reducedGrid) THEN
            CALL LinearIJtoIBJB(buf(1:iMax,1:jMax,2),tg2)
            !CALL AveBoxIJtoIBJB(buf(1:iMax,1:jMax,2),tg2)
    ELSE
       CALL IJtoIBJB(buf(1:iMax,1:jMax,2) ,tg2 )
    END IF
    !tg2(1:ibMax,1:jbMax)=MAX(fgtmp(1:ibMax,1,1:jbMax),tg2(1:ibMax,1:jbMax))
    IF (reducedGrid) THEN
             CALL LinearIJtoIBJB(buf(1:iMax,1:jMax,3),tg3)
            !CALL AveBoxIJtoIBJB(buf(1:iMax,1:jMax,3),tg3)
    ELSE
       CALL IJtoIBJB(buf(1:iMax,1:jMax,3) ,tg3 )
    END IF
    !tg3(1:ibMax,1:jbMax)=MAX(fgtmp(1:ibMax,1,1:jbMax),tg3(1:ibMax,1:jbMax))
    IF (reducedGrid) THEN
            CALL LinearIJtoIBJB(buf(1:iMax,1:jMax,4),zorl)
            !CALL AveBoxIJtoIBJB(buf(1:iMax,1:jMax,4),zorl)
    ELSE
       CALL IJtoIBJB(buf(1:iMax,1:jMax,4),zorl )
    END IF
    z0    =zorl
    sinmax=150.0_r8
    !
    !     use rvisd as temporary for abs(soilm)
    !
    DO j=1,jbMax
       DO i=1,ibMaxPerJB(j)
          rVisDiff(i,j)=ABS(soilm(i,j))
       END DO
    END DO


    CALL sibwet(ibMax,jbMax,rVisDiff,sinmax,iMaskSSiB,wsib,ssib, &
         mxiter,ibMaxPerJB)

    !-srf--------------------------------
    ppli  =0.0_r8
    ppci  =0.0_r8
    capac0=0.0_r8
    capacm=0.0_r8

    !
    !     td0 (deep soil temp) is temporarily defined as tg3
    !
    !$OMP PARALLEL DO PRIVATE(ncount,i)
    DO j=1,jbMax
       ncount=0
       DO i=1,ibMaxPerJB(j)
          gl0(i,j)=xl0
          Mmlen(i,j)=xl0
          IF(iMaskIBIS(i,j) .NE. 0_i8)gtsea(i,j)=290.0_r8
          tseam(i,j)=gtsea(i,j)
          TSK (I,J)=ABS(gtsea(i,j))
          IF (omlmodel) THEN
             HML  (i,j) = oml_hml0 - 13.5_r8*log(MAX(ABS(TSK(i,j))-tice+0.01_r8,1.0_r8))
             HUML (I,J)=0.0_r8
             HVML (I,J)=0.0_r8
          END IF

          IF(iMaskIBIS(i,j).EQ.0_i8) THEN
             IF(-gtsea(i,j).LT.t0) THEN
                iMaskIBIS(i,j)=-1_i8
             END IF
          ELSE
             ncount=ncount+1
             IF(iglsm_w == 0) THEN
                w0        (ncount,1,j)=wsib(i,j)
                w0        (ncount,2,j)=wsib(i,j)
                w0        (ncount,3,j)=wsib(i,j)
             ELSE
                !-srf--------------------------------
                w0        (ncount,1,j)=wsib3d(i,j,1)
                w0        (ncount,2,j)=wsib3d(i,j,2)
                w0        (ncount,3,j)=wsib3d(i,j,3)
                !-srf--------------------------------
             END IF

             td0   (ncount,  j)=tg3 (i,j)

             IF(iglsm_w == 0) THEN
                wm        (ncount,1,j)=wsib(i,j)
                wm        (ncount,2,j)=wsib(i,j)
                wm        (ncount,3,j)=wsib(i,j)
             ELSE
                !-srf--------------------------------
                wm        (ncount,1,j)=wsib3d(i,j,1)
                wm        (ncount,2,j)=wsib3d(i,j,2)
                wm        (ncount,3,j)=wsib3d(i,j,3)
                !-srf--------------------------------
             END IF

             tdm   (ncount,  j)=tg3 (i,j)
             tgm   (ncount,  j)=tg3 (i,j)
             tcm   (ncount,  j)=tg3 (i,j)
             ssib  (ncount,j  )=0.0_r8
             IF(soilm(i,j).LT.0.0_r8)ssib(ncount,j)=wsib(i,j)

             IF(sheleg(i,j).GT.zero) THEN
                capac0(ncount,2,j)=sheleg(i,j)/thousd
                capacm(ncount,2,j)=sheleg(i,j)/thousd
             END IF

          END IF
       END DO
    END DO
    !$OMP END PARALLEL DO
    !
    ! initialize some model variables for cold start conditions
    ! 
    DO j = 1,jbMax
       nLndPts=0
       DO i = 1, ibMax
          gl0(i,j)=xl0
          Mmlen(i,j)=xl0
          IF (iMaskIBIS(i,j) >= 1_i8) THEN
             nLndPts=nLndPts+1
             ! fi        (:,:)  ! fractional snow cover
             ! tsno   (:,:,:) ! temperature of snow layers (K)
             ! hsno   (:,:,:)! thickness of snow layers (m)
             IF(sheleg(i,j).GT.0.0_r8) THEN
                fi   (nLndPts,  j) = 1.0_r8
             END IF
          END IF
       END DO
    END DO
    DO k=1,nsnolay
       DO j = 1,jbMax
          nLndPts=0
          DO i = 1, ibMax
             IF (iMaskIBIS(i,j) >= 1_i8) THEN
                nLndPts=nLndPts+1
                IF(sheleg(i,j).GT.0.0_r8) THEN
                   hsno (nLndPts,k,j) = sheleg(i,j)/1000.0_r8
                ELSE
                   hsno (nLndPts,k,j) = 0.0_r8  
                END IF
                tsno (nLndPts,k,j) = 273.16_r8                     
             END IF
          END DO
       END DO
    END DO
    DO k=1,nsoilay
       DO j = 1,jbMax
          nLndPts=0
          DO i = 1, ibMax
             IF (iMaskIBIS(i,j) >= 1_i8) THEN
                nLndPts=nLndPts+1
                tg   (nLndPts,  j) = tg3 (i,j)
                tsoim(nLndPts,k,j) = tg3 (i,j)
                tsoi (nLndPts,k,j) = tg3 (i,j)
                wsoi (nLndPts,k,j) = wsib(i,j)
                wsoim(nLndPts,k,j) = wsib(i,j)
                IF (iMaskIBIS(i,j) >= 15_i8) wisoi(nLndPts,k,j) = 1.0_r8
             END IF
          END DO
       END DO
    END DO
    !     
    !     Initialize temperature and snow depths in Antarctica and Groenland
    !     
    !      IF (rdlsf) THEN 
    !lonscale
    !latscale
    DO j = 1,jbMax
       nLndPts=0
       DO i = 1, ibMax
          IF (iMaskIBIS(i,j) >= 1_i8) THEN
             nLndPts=nLndPts+1
             !
             !                 Antarctica
             !
             IF (latscale(i,j) .LE. -60.0_r8) THEN
                DO k = 1, nsoilay
                   tsoi (nLndPts,k,j) = tg3 (i,j)!xint(nLndPts,1,j) + 273.16_r8
                   tsoim(nLndPts,k,j) = tg3 (i,j)!xint(nLndPts,1,j) + 273.16_r8
                   wisoi(nLndPts,k,j) = 1.0_r8
                END DO
                tg  (nLndPts,j) = tg3 (i,j)!xint(nLndPts,1,j) + 273.16_r8
                hsno(nLndPts,1,j) = 0.05_r8
                hsno(nLndPts,2,j) = 2.0_r8
                hsno(nLndPts,3,j) = 2.0_r8
                fi  (nLndPts  ,j) = 1.0_r8
             END IF
             !
             !                 Greenland
             !
             IF (latscale(i,j) >= 60.0_r8  .AND. latscale(i,j) <= 85.0_r8  .AND. &
                  lonscale(i,j) <= 330.0_r8 .AND. lonscale(i,j) >  285.0_r8      ) THEN
                DO k = 1, nsoilay
                   tsoi (nLndPts,k,j) = tg3 (i,j)!xint(nLndPts,1,j) + 273.160_r8
                   tsoim(nLndPts,k,j) = tg3 (i,j)!xint(nLndPts,1,j) + 273.160_r8
                   wisoi(nLndPts,k,j) = 1.00_r8
                END DO
                tg  (nLndPts  ,j) = tg3 (i,j)!xint(nLndPts,1,j) + 273.160_r8
                hsno(nLndPts,1,j) = 0.050_r8
                hsno(nLndPts,2,j) = 2.0_r8
                hsno(nLndPts,3,j) = 2.0_r8
                fi  (nLndPts  ,j) = 1.0_r8
             END IF

             IF (latscale(i,j) >= 66.0_r8  .AND. latscale(i,j) <= 85.0_r8  .AND. &
                  lonscale(i,j) <= 345.0_r8 .AND. lonscale(i,j) > 330.0_r8      ) THEN
                DO k = 1, nsoilay
                   tsoi (nLndPts,k,j) = tg3 (i,j)!xint(nLndPts,1,j) + 273.160_r8
                   tsoim(nLndPts,k,j) = tg3 (i,j)!xint(nLndPts,1,j) + 273.160_r8
                   wisoi(nLndPts,k,j) = 1.00_r8
                END DO
                tg  (nLndPts  ,j) = tg3 (i,j)!xint(nLndPts,1,j) + 273.16_r8
                hsno(nLndPts,1,j) = 0.05_r8
                hsno(nLndPts,2,j) = 2.0_r8
                hsno(nLndPts,3,j) = 2.0_r8
                fi  (nLndPts  ,j) = 1.0_r8    
             END IF
          END IF
       END DO
    END DO
    !      END IF

    !         do k=lbeg, lend
    ! Antarctica
    !           if (lati(k)*180./pi .le. -60.) then
    !     
    !              do l = 1, nsoilay
    !                 tsoi(k,l) = xint(k,1) + 273.16
    !                 wisoi(k,l) = 1.0
    !              end do
    !              tg(k) = xint(k,1) + 273.16
    !              hsno(k,1) = 0.05
    !              hsno(k,2) = 2.
    !              hsno(k,3) = 2.
    !              fi(k) = 1.
    !
    !           end if
    ! Greenland
    !          if (lati(k)*180/pi .ge. 60  .and. 
    !    >         lati(k)*180/pi .le. 85  .and.
    !    >         loni(k)*180/pi .le. 330 .and. 
    !    >         loni(k)*180/pi .gt. 285      ) then
    !     
    !             do l = 1, nsoilay
    !                tsoi(k,l) = xint(k,1) + 273.16
    !                wisoi(k,l) = 1.0
    !             end do
    !             tg(k) = xint(k,1) + 273.16
    !             hsno(k,1) = 0.05
    !             hsno(k,2) = 2.
    !             hsno(k,3) = 2.
    !             fi(k) = 1.
    !     
    !          end if
    !     
    !          if (lati(k)*180/pi .ge. 66  .and. 
    !    >         lati(k)*180/pi .le. 85  .and.
    !    >         loni(k)*180/pi .le. 345 .and. 
    !    >         loni(k)*180/pi .gt. 330      ) then
    !     
    !             do l = 1, nsoilay
    !                tsoi(k,l) = xint(k,1) + 273.16
    !                wisoi(k,l) = 1.0
    !             end do
    !             tg(k) = xint(k,1) + 273.16
    !             hsno(k,1) = 0.05
    !             hsno(k,2) = 2.
    !             hsno(k,3) = 2.
    !              fi(k) = 1.
    !     
    !           end if
    !
    !        end do
    !     
    !      end if

    !
    ! return to main program
    !
    RETURN
  END SUBROUTINE coldstart


  ! vegin  :reads vegetation morphoLOGICAL and physioLOGICAL data.




  SUBROUTINE vegin (iMax,jMax, nfsibd,nfprt,nfsibt,fNameSibVeg,fNameSibmsk)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: iMax
    INTEGER, INTENT(IN) :: jMax
    INTEGER, INTENT(IN) :: nfsibd
    INTEGER, INTENT(IN) :: nfprt
    INTEGER, INTENT(IN) :: nfsibt
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameSibVeg
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameSibmsk
    INTEGER, PARAMETER ::  njj=6,nj=9, nk=3,ild=2
    INTEGER , PARAMETER :: ityp = 13 ! Number of Vegetation Types
    INTEGER , PARAMETER :: imon = 12 ! Number of Months with Defined Vegetation Types
    INTEGER , PARAMETER :: icg  = 2  ! Number of Vegetation Parameters
    INTEGER , PARAMETER :: iwv  = 3  ! Number of Radiation Wavelengths
    INTEGER , PARAMETER :: idp  = 3  ! Number of Soil Layer Parameters
    INTEGER , PARAMETER :: ibd  = 2  ! Number of Vegetation Stage
    INTEGER(KIND=i4) ::  ibuf (iMax,jMax)
    INTEGER ::  LRecIN
    INTEGER(KIND=i8) :: imask_in(iMax,jMax)

    ! Vegetation and Soil Parameters

    REAL (KIND=r4) rstpar_r4(ityp,icg,iwv), &
         chil_r4(ityp,icg), &
         topt_r4(ityp,icg), &
         tll_r4(ityp,icg), &
         tu_r4(ityp,icg), &
         defac_r4(ityp,icg), &
         ph1_r4(ityp,icg), &
         ph2_r4(ityp,icg), &
         rootd_r4(ityp,icg), &
         bee_r4(ityp), &
         phsat_r4(ityp), &
         satco_r4(ityp), &
         poros_r4(ityp), &
         zdepth_r4(ityp,idp), &
         green_r4(ityp,imon,icg), &
         xcover_r4(ityp,imon,icg), &
         zlt_r4(ityp,imon,icg), &
         x0x_r4(ityp,imon),&
         xd_r4(ityp,imon), &
         z2_r4   (ityp,imon), &
         z1_r4   (ityp,imon), &
         xdc_r4  (ityp,imon), &
         xbc_r4  (ityp,imon)


    !INTEGER :: jcg
    !INTEGER :: jmon
    !INTEGER :: jtyp
    !INTEGER :: iv
    !INTEGER :: im
    INTEGER :: i
    INTEGER :: ierr
    !
    !
    !     vegetation and soil parameters
    !
    ALLOCATE(bee   (ityp)          )
    ALLOCATE(phsat (ityp)          )
    ALLOCATE(poros_sib (ityp)          )
    ALLOCATE(zdepth(ityp,idp)          )
    ALLOCATE(xcover_fixed(ityp,imon,icg)  )
    ALLOCATE(zlt_fixed   (ityp,imon,icg)  )

    OPEN(UNIT=nfsibd, FILE=TRIM(fNameSibVeg),FORM='UNFORMATTED', ACCESS='SEQUENTIAL',&
         ACTION='READ',STATUS='OLD', IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(fNameSibVeg), ierr
       STOP "**(ERROR)**"
    END IF


    READ (UNIT=nfsibd) rstpar_r4, chil_r4, topt_r4, tll_r4, tu_r4, defac_r4, ph1_r4, ph2_r4, &
         rootd_r4, bee_r4, phsat_r4, satco_r4, poros_r4, zdepth_r4
    READ (UNIT=nfsibd) green_r4, xcover_r4, zlt_r4, x0x_r4, xd_r4, z2_r4, z1_r4, xdc_r4, xbc_r4

    bee          = bee_r4
    phsat        = phsat_r4
    poros_sib        = poros_r4
    zdepth       = zdepth_r4
    xcover_fixed = xcover_r4
    zlt_fixed    = zlt_r4

    CLOSE(nfsibd)

    bee(13) = 4.8_r8
    phsat(13) = -0.167_r8
    poros_sib(13) = 0.4352_r8
    DO i = 1, imon
       zlt_fixed(13,i,1)    = 0.0001_r8
       zlt_fixed(13,i,2)    = 0.0001_r8
       xcover_fixed(13,i,1) = 0.0001_r8
       xcover_fixed(13,i,2) = 0.0001_r8
    END DO

    zdepth(13,1) = 1.0_r8
    zdepth(13,2) = 1.0_r8
    zdepth(13,3) = 1.0_r8
    ibuf=0
    INQUIRE (IOLENGTH=LRecIN) ibuf

    OPEN (UNIT=nfsibt, FILE=TRIM(fNameSibmsk),FORM='FORMATTED', ACCESS='SEQUENTIAL', &
         ACTION='READ',STATUS='OLD', IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(fNameSibmsk), ierr
       STOP "**(ERROR)**"
    END IF

    READ(nfsibt, *) ibuf
    CLOSE(nfsibt)
    imask_in=ibuf
    IF (reducedGrid) THEN
       CALL FreqBoxIJtoIBJB(imask_in,iMaskSSiB)
    ELSE
       CALL IJtoIBJB( imask_in,iMaskSSiB)
    END IF
  END SUBROUTINE vegin


  ! sibwet :transform mintz-serafini and national meteoroLOGICAL center fields
  !         of soil moisture into sib compatible fields of soil moisture.




  SUBROUTINE sibwet &
       (ibmax,jbmax,sinp,sinmax,imask,wsib,ssib,mxiter,ibMaxPerJB)
    !
    !
    !     piers sellers : 29 april 1987
    !
    !
    !   input  :   sinp     = mintz-serafini or national meteoroLOGICAL
    !                         center soil moisture (mm)
    !              sinmax   = maximum value of sinp (mm)
    !              wsinp    = m-s or nmc fractional wetness
    !              ms       = 1, mintz-serafini
    !              nmc      = 1, national meteoroLOGICAL center
    !              bee      = sib : soil moisture potential factor
    !              phsat    = sib : soil potential at saturation (m)
    !              zdepth(3)= sib : depth of 3 soil layers (m)
    !              poros    = Porosidade do solo (m"3/m"3)
    !
    !   output :   wsibt    = sib : fractional wetness
    !              ssibt    = sib : soil moisture content (m)
    !              psit     = sib : soil moisture potential (m)
    !              factor   = sib : extraction factor
    !
    INTEGER, INTENT(in   ) :: ibmax
    INTEGER, INTENT(in   ) :: jbmax
    INTEGER, INTENT(in   ) :: mxiter
    REAL(KIND=r8)   , INTENT(in   ) :: sinp(ibmax,jbmax)
    REAL(KIND=r8)   , INTENT(in   ) :: sinmax
    !

    INTEGER(KIND=i8), INTENT(in   ) :: imask (ibmax,jbmax)
    REAL(KIND=r8)   , INTENT(inout  ) :: wsib  (ibmax,jbmax)
    REAL(KIND=r8)   , INTENT(inout  ) :: ssib  (ibmax,jbmax)
    INTEGER, INTENT(in   ) :: ibMaxPerJB(jbMax)

    INTEGER , PARAMETER :: ityp = 13 ! Number of Vegetation Types
    INTEGER , PARAMETER :: imon = 12 ! Number of Months with Defined Vegetation Types
    INTEGER , PARAMETER :: icg  = 2  ! Number of Vegetation Parameters
    INTEGER , PARAMETER :: iwv  = 3  ! Number of Radiation Wavelengths
    INTEGER , PARAMETER :: idp  = 3  ! Number of Soil Layer Parameters
    INTEGER , PARAMETER :: ibd  = 2  ! Number of Vegetation Stage
    REAL(KIND=r8) :: sm(ityp,mxiter)
    REAL(KIND=r8) :: time(ityp,mxiter)
    REAL(KIND=r8) :: fact(ityp,mxiter)

    REAL(KIND=r8), PARAMETER :: xph1(13,2) = RESHAPE( &
         (/-100.0_r8,-190.0_r8,-200.0_r8,-200.0_r8,-200.0_r8,-120.0_r8,-120.0_r8,-120.0_r8,-200.0_r8, &
         -200.0_r8, -10.0_r8,-190.0_r8, -10.0_r8,-100.0_r8,-190.0_r8,-200.0_r8,-200.0_r8,-200.0_r8, &
         -120.0_r8,-120.0_r8,-120.0_r8,-200.0_r8,-200.0_r8, -10.0_r8,-190.0_r8, -10.0_r8/), &
         (/13,2/))
    REAL(KIND=r8), PARAMETER :: xph2(13,2) = RESHAPE( &
         (/-500.0_r8,-250.0_r8,-250.0_r8,-250.0_r8,-250.0_r8,-230.0_r8,-230.0_r8,-280.0_r8,-400.0_r8, &
         -400.0_r8,-100.0_r8,-250.0_r8,-100.0_r8,-500.0_r8,-250.0_r8,-250.0_r8,-250.0_r8,-250.0_r8, &
         -230.0_r8,-230.0_r8,-280.0_r8,-400.0_r8,-400.0_r8,-100.0_r8,-250.0_r8,-100.0_r8/) , &
         (/13,2/))

    REAL(KIND=r8)    :: tzdep(3)
    REAL(KIND=r8)    :: tzltm(2)
    REAL(KIND=r8)    :: sibmax(ityp)
    REAL(KIND=r8)    :: tphsat
    REAL(KIND=r8)    :: tbee
    REAL(KIND=r8)    :: tporos
    INTEGER :: imm1
    INTEGER :: imm2
    INTEGER :: is
    INTEGER :: im
    INTEGER :: imm
    INTEGER :: ivegm
    REAL(KIND=r8)    :: cover
    REAL(KIND=r8)    :: tph1
    REAL(KIND=r8)    :: tph2
    REAL(KIND=r8)    :: sref
    REAL(KIND=r8)    :: smin
    REAL(KIND=r8)    :: dssib
    REAL(KIND=r8)    :: dw
    REAL(KIND=r8)    :: times
    REAL(KIND=r8)    :: soilmo
    REAL(KIND=r8)    :: w
    REAL(KIND=r8)    :: rsoilm
    INTEGER :: iter
    INTEGER :: latmax
    INTEGER :: lonmax
    INTEGER :: lat
    INTEGER :: lon
    REAL(KIND=r8)    :: tsinp
    REAL(KIND=r8)    :: etp
    REAL(KIND=r8)    :: facmod
    REAL(KIND=r8)    :: ssibt
    REAL(KIND=r8)    :: psit
    REAL(KIND=r8)    :: factor
    REAL(KIND=r8)    :: dt
    INTEGER :: itsoil
    INTEGER :: itfac

    sm  =0.0_r8
    time=0.0_r8
    fact=0.0_r8
    ssib=0.0_r8
    wsib=0.0_r8

    lonmax=ibmax
    latmax=jbmax

    DO is = 1,ityp
       !zdepth(3)= sib : depth of 3 soil layers (m)
       tzdep (1)= zdepth(is,1)
       tzdep (2)= zdepth(is,2)
       tzdep (3)= zdepth(is,3)
       tphsat   = phsat (is)
       tbee     = bee   (is)
       tporos   = poros_sib (is)
       imm1=1
       imm2=1
       tzltm(1)=zlt_fixed(is,1,1)
       tzltm(2)=zlt_fixed(is,1,2)
       DO im=2,12
          IF(tzltm(1).LE.zlt_fixed(is,im,1) ) THEN
             imm1=im
             tzltm(1)=zlt_fixed(is,im,1)
          END IF
          IF(tzltm(2).LE.zlt_fixed(is,im,2) )THEN
             imm2=im
             tzltm(2)=zlt_fixed(is,im,2)
          END IF
       END DO
       imm=imm1
       ivegm=1
       IF(tzltm(1).LE.tzltm(2)) THEN
          imm=imm2
          ivegm=2
       END IF
       !
       !     xcover......Fracao de cobertura vegetal icg=1 topo
       !     xcover......Fracao de cobertura vegetal icg=2 base
       !
       cover=xcover_fixed(is,imm,ivegm)
       tph1=xph1         (is,ivegm)
       tph2=xph2         (is,ivegm)
       !
       !                                                     m^3
       ! sibmax(is) =(Z1 + Z2 + Z3) * poros = [m + m + m] * ----- = m = Os
       !                                                     m^3
       !
       sibmax(is) = ( tzdep(1) + tzdep(2) + tzdep(3) ) * tporos
       !
       IF(nfctrl(83).GE.1)WRITE(UNIT=nfprt,FMT=999)is,sibmax(is),tzdep(1), &
            tzdep(2),tzdep(3),tporos
       !
       !            bee      = soil moisture potential factor
       !            phsat    = soil potential at saturation   (m)
       !
       !                   --              --
       !                  | log ( - tphsat/1)|
       !  O  = Os * EXP * | -----------------|
       !                  |        b         |
       !                   --              --
       !
       sref = sibmax(is) * EXP( LOG(tphsat /(-1.0e0_r8)) /tbee)
       !                   --                          --
       !                  | log ( - tphsat/(-1.0e10) )   |
       !Omin = Os * EXP * | -----------------------------|
       !                  |              b               |
       !                   --                          --
       !
       smin    = sibmax(is) * EXP( LOG(tphsat /(-1.0e10_r8)) / tbee)
       !
       !             O - Omin
       !dssib  = ------------------
       !              mxiter
       !
       dssib   = (sref - smin) / REAL(mxiter,r8)
       !
       !              O - Omin
       ! dw    =  ------------------
       !             mxiter*Os
       !
       dw      = dssib / sibmax(is)
       !
       times   = 0.0e0_r8
       soilmo  = sref
       !
       !       O
       ! w = -----
       !       Os
       !
       w = soilmo / sibmax(is)
       !
       !                      --             --
       !                     |       0.0027    |
       !rsoilm  = 101840.0 * |1.0 - w          |
       !                     |                 |
       !                      --             --
       !
       rsoilm  = 101840.0_r8 * (1.0_r8 - w**0.0027_r8)
       DO iter = 1, mxiter
          CALL extrak( w   ,dw  ,tbee,tphsat, rsoilm, cover, &
               tph1,tph2,psit,factor )
          !
          !       dssib
          !dt = ----------
          !       factor
          !
          dt            = dssib  / factor
          !
          soilmo        = soilmo - dssib
          !
          !       O
          ! w = -----
          !       Os
          !
          w             = soilmo / sibmax(is)
          times         = times  + dt
          sm  (is,iter) = soilmo
          time(is,iter) = times
          fact(is,iter) = factor
       END DO

    END DO
    !
    !     input soil moisture map is now transformed to sib fields.
    !
    DO lat = 1, latmax
       DO lon = 1, ibMaxPerJB(lat)
          is=imask(lon,lat)
          IF(is.NE.0)THEN
             tsinp = sinp(lon,lat)
             tsinp = MAX (sinmax/100.0e3_r8 , tsinp )
             tsinp = MIN (sinmax,tsinp)
             IF (tsinp .GT. 0.75e0_r8*sinmax ) etp = sinmax - tsinp
             facmod=MIN(1.0e0_r8,tsinp/(0.75e0_r8*sinmax) )
             IF (tsinp .LE. 0.75e0_r8*sinmax ) THEN
                etp = 0.75e0_r8*sinmax*LOG(0.75e0_r8*sinmax/tsinp ) + 0.25e0_r8*sinmax
             END IF
             etp = etp / 1000.0e0_r8
             DO iter = 1, mxiter
                itsoil=iter
                IF ( time(is,iter) - etp .GT. 0.0e0_r8  ) EXIT
             END DO
             DO iter=1,mxiter
                itfac=iter
                IF( fact(is,iter)-facmod-0.01e0_r8.LT.0.0e0_r8)EXIT
             END DO
             ssibt=MIN(sm(is,itsoil),sm(is,itfac))
             DO iter=1,mxiter
                IF(ssibt.GT.sm(is,iter))EXIT
             END DO
             ssib(lon,lat) = sm(is,iter)
             !
             !          O
             ! wsib = -----
             !         Os
             !
             wsib(lon,lat) = sm(is,iter) / sibmax(is)
          END IF
       END DO
    END DO
999 FORMAT(' IS,MAX,D1,D2,D3,POROS=',I2,1X,5E12.5)
  END SUBROUTINE sibwet



  !
  SUBROUTINE extrak( w, dw, tbee, tphsat, rsoilm, cover, tph1, tph2, &
       psit, factor )
    IMPLICIT NONE

    REAL(KIND=r8), INTENT(in   ) :: w
    REAL(KIND=r8), INTENT(in   ) :: dw
    REAL(KIND=r8), INTENT(in   ) :: tbee
    REAL(KIND=r8), INTENT(in   ) :: tphsat
    REAL(KIND=r8), INTENT(in   ) :: rsoilm
    REAL(KIND=r8), INTENT(in   ) :: cover
    REAL(KIND=r8), INTENT(in   ) :: tph1
    REAL(KIND=r8), INTENT(in   ) :: tph2
    REAL(KIND=r8), INTENT(inout  ) :: psit
    REAL(KIND=r8), INTENT(inout  ) :: factor
    REAL(KIND=r8) :: rsoil
    REAL(KIND=r8) :: argg
    REAL(KIND=r8) :: hr
    REAL(KIND=r8) :: rplant
    !                --     -- (-b)
    !               |      dw |                  0
    ! psit = PHYs * | w - --- |      where w = -----
    !               |      2  |                  0s
    !                --     --
    psit   = tphsat * ( w-dw/2.0e0_r8 ) ** (-tbee)
    !
    !                      --                        --
    !                     |       --     -- (0.0027)   |
    !                     |      |      dw |           |
    !rsoil   = 101840.0 * |1.0 - | w - --- |           |
    !                     |      |      2  |           |
    !                     |       --     --            |
    !                      --                        --
    !
    rsoil  = 101840.0_r8 * (1.0_r8-( w-dw/2.0_r8) ** 0.0027_r8)
    !
    !                9.81       1
    !argg = psit * -------- * -------
    !               461.50     310.0
    !
    argg   = MAX ( -10.0e0_r8 , ((psit * 9.81e0_r8 / 461.5e0_r8) / 310.e0_r8))
    !
    !            --                       --
    !           |         9.81       1      |
    !hr   = EXP |psit * -------- * -------  |
    !           |        461.50     310.0   |
    !            --                       --
    !
    hr     = EXP ( argg )
    !
    !         rsoilm
    ! rsoil =--------- * hr
    !         rsoil
    !
    rsoil  = rsoilm /rsoil * hr
    !
    !          ( psit - tph2 - 50.0)
    !rplant = -------------------------
    !             ( tph1 - tph2 )
    !
    rplant = ( psit - tph2 -50.0_r8) / ( tph1 - tph2 )
    rplant = MAX ( 0.0e0_r8, MIN ( 1.0e0_r8, rplant ) )
    !                                                                     --                   --
    !                  --                 --                             |     --     -- (0.0027)|
    !                 |( psit - tph2 - 50)  |                            |    |      dw |        |
    !factor = cover * |---------------------| + (1 - cover) * 101840.0 * |1 - | w - --- |        |
    !                 |   ( tph1 - tph2 )   |                            |    |      2  |        |
    !                  --                 --                             |     --     --         |
    !                                                                     --                   --
    factor = cover * rplant + ( 1.0e0_r8 - cover ) * rsoil
    factor = MAX ( 1.e-6_r8, factor )
  END SUBROUTINE extrak

  !
  ! ---------------------------------------------------------------------
  SUBROUTINE inisurf(irestart , &! INTENT(IN   )
       totcondub, &! INTENT(OUT  ) :: totcondub(npoi)    
       totconduc, &! INTENT(OUT  ) :: totconduc(npoi)    
       totcondls, &! INTENT(OUT  ) :: totcondls(npoi)    
       totcondl3, &! INTENT(OUT  ) :: totcondl3(npoi)    
       totcondl4, &! INTENT(OUT  ) :: totcondl4(npoi)    
       tu       , &! INTENT(OUT  ) :: tu        (npoi)    
       ts       , &! INTENT(OUT  ) :: ts        (npoi)    
       tl          , &! INTENT(OUT  ) :: tl        (npoi)    
       tlsub    , &! INTENT(OUT  ) :: tlsub        (npoi)    
       t12      , &! INTENT(OUT  ) :: t12        (npoi)    
       t34      , &! INTENT(OUT  ) :: t34        (npoi)    
       q12          , &! INTENT(OUT  ) :: q12        (npoi)    
       q34          , &! INTENT(OUT  ) :: q34        (npoi)    
       ciub     , &! INTENT(OUT  ) :: ciub        (npoi)    
       ciuc     , &! INTENT(OUT  ) :: ciuc        (npoi)    
       cils          , &! INTENT(OUT  ) :: cils        (npoi)    
       cil3     , &! INTENT(OUT  ) :: cil3        (npoi)    
       cil4     , &! INTENT(OUT  ) :: cil4        (npoi)    
       csub     , &! INTENT(OUT  ) :: csub        (npoi)    
       csuc          , &! INTENT(OUT  ) :: csuc        (npoi)    
       csls     , &! INTENT(OUT  ) :: csls        (npoi)    
       csl3     , &! INTENT(OUT  ) :: csl3        (npoi)    
       csl4     , &! INTENT(OUT  ) :: csl4        (npoi)    
       gsub          , &! INTENT(OUT  ) :: gsub        (npoi)    
       gsuc     , &! INTENT(OUT  ) :: gsuc        (npoi)    
       gsls     , &! INTENT(OUT  ) :: gsls        (npoi)    
       gsl3     , &! INTENT(OUT  ) :: gsl3        (npoi)    
       gsl4          , &! INTENT(OUT  ) :: gsl4        (npoi)    
       clitlm   , &! INTENT(OUT  ) :: clitlm   (npoi)   !
       clitls   , &! INTENT(OUT  ) :: clitls   (npoi)   !
       clitll   , &! INTENT(OUT  ) :: clitll   (npoi)   !
       clitrm   , &! INTENT(OUT  ) :: clitrm   (npoi)   !
       clitrs   , &! INTENT(OUT  ) :: clitrs   (npoi)   !
       clitrl   , &! INTENT(OUT  ) :: clitrl   (npoi)   !
       clitwm   , &! INTENT(OUT  ) :: clitwm   (npoi)   !
       clitws   , &! INTENT(OUT  ) :: clitws   (npoi)   !
       clitwl   , &! INTENT(OUT  ) :: clitwl   (npoi)   !
       totcmic  , &! INTENT(OUT  ) :: totcmic  (npoi)   !
       csoislop , &! INTENT(OUT  ) :: csoislop (npoi)   !
       csoislon , &! INTENT(OUT  ) :: csoislon (npoi)   !
       csoipas  , &! INTENT(OUT  ) :: csoipas  (npoi)   !
       totlit   , &! INTENT(OUT  ) :: totlit   (npoi)   !
       totnlit  , &! INTENT(OUT  ) :: totnlit  (npoi)   !
       totfall  , &! INTENT(OUT  ) :: totfall  (npoi)   !
       totalit  , &! INTENT(OUT  ) :: totalit  (npoi)   !
       totrlit  , &! INTENT(OUT  ) :: totrlit  (npoi)   !
       totanlit , &! INTENT(OUT  ) :: totanlit (npoi)   !
       totrnlit , &! INTENT(OUT  ) :: totrnlit (npoi)   !
       totcsoi  , &! INTENT(OUT  ) :: totcsoi  (npoi)   !
       totnmic  , &! INTENT(OUT  ) :: totnmic  (npoi)   !
       tco2mic  , &! INTENT(OUT  ) :: tco2mic  (npoi)   !
       tnpptot  , &! INTENT(OUT  ) :: tnpptot  (npoi)   !
       tneetot  , &! INTENT(OUT  ) :: tneetot  (npoi)   !
       tnmin    , &! INTENT(OUT  ) :: tnmin         (npoi)   !
       cdisturb , &! INTENT(OUT  ) :: cdisturb (npoi)   !
       tempu          , &! INTENT(OUT  ) :: tempu         (npoi)
       templ    , &! INTENT(OUT  ) :: templ         (npoi)
       dropu    , &! INTENT(OUT  ) :: dropu         (npoi)
       dropls   , &! INTENT(OUT  ) :: dropls   (npoi)
       dropl4   , &! INTENT(OUT  ) :: dropl4   (npoi)
       dropl3   , &! INTENT(OUT  ) :: dropl3   (npoi)
       wliqu    , &! INTENT(OUT  ) :: wliqu         (npoi)
       wliqs    , &! INTENT(OUT  ) :: wliqs         (npoi)
       wliql          , &! INTENT(OUT  ) :: wliql         (npoi)
       wsnou    , &! INTENT(OUT  ) :: wsnou         (npoi)
       wsnos    , &! INTENT(OUT  ) :: wsnos         (npoi)
       wsnol    , &! INTENT(OUT  ) :: wsnol         (npoi)
       su          , &! INTENT(OUT  ) :: su       (npoi)
       ss          , &! INTENT(OUT  ) :: ss       (npoi)
       sl       , &! INTENT(OUT  ) :: sl
       ginvap   , &! INTENT(OUT  ) :: ginvap   (npoi)
       gsuvap   , &! INTENT(OUT  ) :: gsuvap   (npoi)
       gtrans   , &! INTENT(OUT  ) :: gtrans   (npoi)
       grunof   , &! INTENT(OUT  ) :: grunof   (npoi)
       gdrain   , &! INTENT(OUT  ) :: gdrain   (npoi)
       iwet          , &! INTENT(OUT  ) :: iwet     (npoi)    
       iwetday  , &! INTENT(OUT  ) :: iwetday  (npoi,31)
       precipday, &! INTENT(OUT  ) :: precipday(npoi,31)
       asurd    , &! INTENT(OUT  ) :: asurd    (npoi,nband)
       asuri          , &! INTENT(OUT  ) :: asuri    (npoi,nband)
       xstore   , &! INTENT(OUT  ) :: xstore   (npoi,3)
       nband    , &! INTENT(IN   ) :: nband
       stef          , &! INTENT(OUT  ) :: stef 
       vonk     , &! INTENT(OUT  ) :: vonk 
       grav     , &! INTENT(OUT  ) :: grav 
       tmelt    , &! INTENT(OUT  ) :: tmelt
       hvap          , &! INTENT(OUT  ) :: hvap 
       hfus     , &! INTENT(OUT  ) :: hfus 
       hsub     , &! INTENT(OUT  ) :: hsub 
       ch2o     , &! INTENT(OUT  ) :: ch2o 
       cice          , &! INTENT(OUT  ) :: cice 
       cair     , &! INTENT(OUT  ) :: cair 
       cvap     , &! INTENT(OUT  ) :: cvap 
       rair     , &! INTENT(OUT  ) :: rair 
       rvap          , &! INTENT(OUT  ) :: rvap 
       cappa    , &! INTENT(OUT  ) :: cappa
       rhow     , &! INTENT(OUT  ) :: rhow 
       vzero    , &! INTENT(OUT  ) :: vzero    (npoi)
       epsilon    )! INTENT(OUT  ) :: epsilon
    ! ---------------------------------------------------------------------
    !
    ! does initialization for model
    !
    IMPLICIT NONE    
    !
    INTEGER, INTENT(IN   ) :: irestart     ! 0 = initial run, 1 = restart run
    INTEGER, INTENT(IN   ) :: nband        ! number of solar radiation wavebands
    REAL(KIND=r8)   , INTENT(OUT  ) :: stef         ! stefan-boltzmann constant (W m-2 K-4)
    REAL(KIND=r8)   , INTENT(OUT  ) :: vonk         ! von karman constant (dimensionless)
    REAL(KIND=r8)   , INTENT(OUT  ) :: grav         ! gravitational acceleration (m s-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: tmelt        ! freezing point of water (K)
    REAL(KIND=r8)   , INTENT(OUT  ) :: hvap         ! latent heat of vaporization of water (J kg-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: hfus         ! latent heat of fusion of water (J kg-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: hsub         ! latent heat of sublimation of ice (J kg-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ch2o         ! specific heat of liquid water (J deg-1 kg-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: cice         ! specific heat of ice (J deg-1 kg-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: cair         ! specific heat of dry air at constant pressure (J deg-1 kg-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: cvap         ! specific heat of water vapor at constant pressure (J deg-1 kg-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: rair         ! gas constant for dry air (J deg-1 kg-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: rvap         ! gas constant for water vapor (J deg-1 kg-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: cappa        ! rair/cair
    REAL(KIND=r8)   , INTENT(OUT  ) :: rhow         ! density of liquid water (all types) (kg m-3)
    REAL(KIND=r8)   , INTENT(OUT  ) :: vzero   (ibMax,jbMax)! a REAL(KIND=r8) array of zeros, of length npoi
    REAL(KIND=r8)   , INTENT(OUT  ) :: epsilon      ! small quantity to avoid zero-divides and other
    ! truncation or machine-limit troubles with small
    ! values. should be slightly greater than o(1)
    ! machine precision!      include 'comatm.h'
    INTEGER, INTENT(OUT  ) :: iwet     (ibMax,jbMax)             ! wet day / dry day flag
    INTEGER, INTENT(OUT  ) :: iwetday  (ibMax,31,jbMax)          ! wet day / dry day flag
    REAL(KIND=r8)   , INTENT(OUT  ) :: precipday(ibMax,31,jbMax)                 
    REAL(KIND=r8)   , INTENT(OUT  ) :: asurd    (ibMax,nband,jbMax) ! direct albedo of surface system
    REAL(KIND=r8)   , INTENT(OUT  ) :: asuri    (ibMax,nband,jbMax) ! diffuse albedo of surface system 
    REAL(KIND=r8)   , INTENT(OUT  ) :: xstore   (ibMax,3,jbMax)     ! weather generator 'memory' matrix
    REAL(KIND=r8)   , INTENT(OUT  ) :: ginvap   (ibMax,jbMax)   ! total evaporation rate from all intercepted h2o (kg_h2o m-2 s-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: gsuvap   (ibMax,jbMax)   ! total evaporation rate from surface (snow/soil) (kg_h2o m-2 s-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: gtrans   (ibMax,jbMax) ! total transpiration rate from all vegetation canopies (kg_h2o m-2 s-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: grunof   (ibMax,jbMax) ! surface runoff rate (kg_h2o m-2 s-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: gdrain   (ibMax,jbMax)  ! drainage rate out of bottom of lowest soil layer (kg_h2o m-2 s-1)

    REAL(KIND=r8)   , INTENT(OUT  ) :: totcondub(ibMax,jbMax)    ! 
    REAL(KIND=r8)   , INTENT(OUT  ) :: totconduc(ibMax,jbMax)    !
    REAL(KIND=r8)   , INTENT(OUT  ) :: totcondls(ibMax,jbMax)    ! 
    REAL(KIND=r8)   , INTENT(OUT  ) :: totcondl3(ibMax,jbMax)    !
    REAL(KIND=r8)   , INTENT(OUT  ) :: totcondl4(ibMax,jbMax)    !
    REAL(KIND=r8)   , INTENT(OUT  ) :: tu       (ibMax,jbMax)! temperature of upper canopy leaves (K)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ts       (ibMax,jbMax)! temperature of upper canopy stems (K)
    REAL(KIND=r8)   , INTENT(OUT  ) :: tl       (ibMax,jbMax)! temperature of lower canopy leaves & stems(K)
    REAL(KIND=r8)   , INTENT(OUT  ) :: tlsub    (ibMax,jbMax)! temperature of lower canopy vegetation buried by snow (K)
    REAL(KIND=r8)   , INTENT(OUT  ) :: t12      (ibMax,jbMax)! air temperature at z12 (K)
    REAL(KIND=r8)   , INTENT(OUT  ) :: t34      (ibMax,jbMax)! air temperature at z34 (K)
    REAL(KIND=r8)   , INTENT(OUT  ) :: q12      (ibMax,jbMax)! specific humidity of air at z12
    REAL(KIND=r8)   , INTENT(OUT  ) :: q34      (ibMax,jbMax)! specific humidity of air at z34
    REAL(KIND=r8)   , INTENT(OUT  ) :: ciub     (ibMax,jbMax)! intercellular co2 concentration - broadleaf (mol_co2/mol_air)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ciuc     (ibMax,jbMax)! intercellular co2 concentration - conifer   (mol_co2/mol_air)
    REAL(KIND=r8)   , INTENT(OUT  ) :: cils     (ibMax,jbMax)! intercellular co2 concentration - shrubs    (mol_co2/mol_air)
    REAL(KIND=r8)   , INTENT(OUT  ) :: cil3     (ibMax,jbMax)! intercellular co2 concentration - c3 plants (mol_co2/mol_air)
    REAL(KIND=r8)   , INTENT(OUT  ) :: cil4     (ibMax,jbMax)! intercellular co2 concentration - c4 plants (mol_co2/mol_air)
    REAL(KIND=r8)   , INTENT(OUT  ) :: csub     (ibMax,jbMax)! leaf boundary layer co2 concentration - broadleaf (mol_co2/mol_air)
    REAL(KIND=r8)   , INTENT(OUT  ) :: csuc     (ibMax,jbMax)! leaf boundary layer co2 concentration - conifer   (mol_co2/mol_air)
    REAL(KIND=r8)   , INTENT(OUT  ) :: csls     (ibMax,jbMax)! leaf boundary layer co2 concentration - shrubs    (mol_co2/mol_air)
    REAL(KIND=r8)   , INTENT(OUT  ) :: csl3     (ibMax,jbMax)! leaf boundary layer co2 concentration - c3 plants (mol_co2/mol_air)
    REAL(KIND=r8)   , INTENT(OUT  ) :: csl4     (ibMax,jbMax)! leaf boundary layer co2 concentration - c4 plants (mol_co2/mol_air)
    REAL(KIND=r8)   , INTENT(OUT  ) :: gsub     (ibMax,jbMax)! upper canopy stomatal conductance - broadleaf  (mol_co2 m-2 s-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: gsuc     (ibMax,jbMax)! upper canopy stomatal conductance - conifer  (mol_co2 m-2 s-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: gsls     (ibMax,jbMax)! lower canopy stomatal conductance - shrubs   (mol_co2 m-2 s-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: gsl3     (ibMax,jbMax)! lower canopy stomatal conductance - c3 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: gsl4     (ibMax,jbMax)! lower canopy stomatal conductance - c4 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: clitlm   (ibMax,jbMax)! carbon in leaf litter pool - metabolic       (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: clitls   (ibMax,jbMax)! carbon in leaf litter pool - structural      (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: clitll   (ibMax,jbMax)! carbon in leaf litter pool - lignin          (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: clitrm   (ibMax,jbMax)! carbon in fine root litter pool - metabolic  (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: clitrs   (ibMax,jbMax)! carbon in fine root litter pool - structural (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: clitrl   (ibMax,jbMax)! carbon in fine root litter pool - lignin     (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: clitwm   (ibMax,jbMax)! carbon in woody litter pool - metabolic      (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: clitws   (ibMax,jbMax)! carbon in woody litter pool - structural     (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: clitwl   (ibMax,jbMax)! carbon in woody litter pool - lignin          (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: totcmic  (ibMax,jbMax)! total carbon residing in microbial pools (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: csoislop (ibMax,jbMax)! carbon in soil - slow protected humus                (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: csoislon (ibMax,jbMax)! carbon in soil - slow nonprotected humus     (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: csoipas  (ibMax,jbMax)! carbon in soil - passive humus                   (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: totlit   (ibMax,jbMax)! total carbon in all litter pools (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: totnlit  (ibMax,jbMax)! total nitrogen in all litter pools (kg_N m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: totfall  (ibMax,jbMax)! total litterfall and root turnover (kg_C m-2/year)
    REAL(KIND=r8)   , INTENT(OUT  ) :: totalit  (ibMax,jbMax)! total standing aboveground litter (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: totrlit  (ibMax,jbMax)! total root litter carbon belowground (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: totanlit (ibMax,jbMax)! total standing aboveground nitrogen in litter (kg_N m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: totrnlit (ibMax,jbMax)! total root litter nitrogen belowground (kg_N m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: totcsoi  (ibMax,jbMax)! total carbon in all soil pools (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: totnmic  (ibMax,jbMax)! total nitrogen residing in microbial pool (kg_N m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: tco2mic  (ibMax,jbMax)! instantaneous microbial co2 flux from soil (mol-CO2 / m-2 / second)
    REAL(KIND=r8)   , INTENT(OUT  ) :: tnpptot  (ibMax,jbMax)! instantaneous npp (mol-CO2 / m-2 / second)
    REAL(KIND=r8)   , INTENT(OUT  ) :: tneetot  (ibMax,jbMax)! instantaneous net ecosystem exchange of co2 per timestep (kg_C m-2/timestep)
    REAL(KIND=r8)   , INTENT(OUT  ) :: tnmin    (ibMax,jbMax)! instantaneous nitrogen mineralization (kg_N m-2/timestep)
    REAL(KIND=r8)   , INTENT(OUT  ) :: cdisturb (ibMax,jbMax)! annual amount of vegetation carbon lost 
    ! to atmosphere due to fire  (biomass burning) (kg_C m-2/year)
    REAL(KIND=r8)   , INTENT(OUT  ) :: tempu    (ibMax,jbMax)! cold-phenology trigger for trees (non-dimensional)
    REAL(KIND=r8)   , INTENT(OUT  ) :: templ    (ibMax,jbMax)! cold-phenology trigger for grasses/shrubs (non-dimensional)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dropu    (ibMax,jbMax)! drought-phenology trigger for trees (non-dimensional)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dropls   (ibMax,jbMax)! drought-phenology trigger for shrubs (non-dimensional)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dropl4   (ibMax,jbMax)! drought-phenology trigger for c4 grasses (non-dimensional)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dropl3   (ibMax,jbMax)! drought-phenology trigger for c3 grasses (non-dimensional)
    REAL(KIND=r8)   , INTENT(OUT  ) :: wliqu    (ibMax,jbMax)! intercepted liquid h2o on upper canopy leaf area (kg m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: wliqs    (ibMax,jbMax)! intercepted liquid h2o on upper canopy stem area (kg m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: wliql    (ibMax,jbMax)! intercepted liquid h2o on lower canopy leaf and stem area (kg m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: wsnou    (ibMax,jbMax)! intercepted frozen h2o (snow) on upper canopy leaf area (kg m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: wsnos    (ibMax,jbMax)! intercepted frozen h2o (snow) on upper canopy stem area (kg m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: wsnol    (ibMax,jbMax)! intercepted frozen h2o (snow) on lower canopy leaf & stem area (kg m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: su       (ibMax,jbMax)! air-vegetation transfer coefficients (*rhoa) for upper
    ! canopy leaves (m s-1 * kg m-3) (A39a Pollard & Thompson 1995)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ss       (ibMax,jbMax)          ! air-vegetation transfer coefficients (*rhoa) for upper
    ! canopy stems (m s-1 * kg m-3) (A39a Pollard & Thompson 1995)
    REAL(KIND=r8)   , INTENT(OUT  ) :: sl       (ibMax,jbMax)          ! air-vegetation transfer coefficients (*rhoa) for lower
    ! canopy leaves & stems (m s-1*kg m-3) (A39a Pollard & Thompson 1995)
    !
    ! Arguments (input)     
    !
    !
    ! local variables
    !
    INTEGER :: i       ! loop indice
    INTEGER :: j       ! loop indice 
    INTEGER :: k       ! loop indice
    !
    ! set physical constants (mks)
    !
    stef  = 5.67e-8_r8 
    vonk  = 0.4_r8
    grav  = 9.80616_r8
    tmelt = 273.16_r8
    hvap  = 2.5104e+6_r8
    hfus  = 0.3336e+6_r8
    hsub  = hvap + hfus
    ch2o  = 4.218e+3_r8
    cice  = 2.106e+3_r8
    cair  = 1.00464e+3_r8
    cvap  = 1.81e+3_r8
    rair  = 287.04_r8
    rvap  = 461.0_r8
    cappa = rair / cair
    rhow  = 1.0e+3_r8
    !
    ! -----------------------------------------------------------------
    ! constant atmospheric co2 and o2
    ! -----------------------------------------------------------------
    o2conc  = 0.209000_r8
    co2conc = 0.000350_r8

    !
    !CALL const (vzero, npoi, 0.0)

    vzero=0.0_r8
    !
    ! specify the epsilon value for the model
    !
    epsilon = 1.0e-12_r8
    !
    ! initialize integer variables (can't use const for this)
    !
    ! wet day / dry day flag initialized to dry day (0)
    !      
    IF(irestart ==0)THEN

       DO j=1,jbMax 
          DO i = 1, nlpoints(j) 
             iwet(i,j) = 0
             DO k = 1,31
                iwetday  (i,k,j) = 0
                precipday(i,k,j) = 0
             END DO
          END DO
       END DO
       !
       ! zero flux arrays, and global diagnostic arrays
       !
       DO k = 1,nband
          DO j = 1,jbMax 
             DO i = 1, nlpoints(j) 
                !CALL const (asurd, npoi*nband, 0.0)
                asurd (i,k,j)= 0.0_r8
                !CALL const (asuri, npoi*nband, 0.0)
                asuri (i,k,j)=0.0_r8
             END DO
          END DO
       END DO
       DO j = 1,jbMax 
          DO i = 1, nlpoints(j) 
             !asuri
             !CALL const (totcondub, npoi, 0.0)
             totcondub(i,j) =   0.0_r8
             !CALL const (totconduc, npoi, 0.0)
             totconduc(i,j)=   0.0_r8
             !CALL const (totcondls, npoi, 0.0)
             totcondls(i,j)=   0.0_r8
             !CALL const (totcondl3, npoi, 0.0)
             totcondl3(i,j)=   0.0_r8
             !CALL const (totcondl4, npoi, 0.0)
             totcondl4(i,j)=   0.0_r8
             !ginvap
             !CALL const (ginvap, npoi, 0.0)
             ginvap(i,j)=0.0_r8
             !CALL const (gsuvap, npoi, 0.0)
             gsuvap(i,j)=0.0_r8
             !CALL const (gtrans, npoi, 0.0)
             gtrans(i,j)=0.0_r8
             !CALL const (grunof, npoi, 0.0)
             grunof(i,j)=0.0_r8
             !CALL const (gdrain, npoi, 0.0)
             gdrain(i,j)=0.0_r8
             !
             ! initialize vegetation prognostic variables
             !
             ! initialize all temperature fields to 10 degrees C
             !
             !CALL const (tu,    npoi, 283.16_r8)
             tu (i,j)= 283.16_r8
             !CALL const (ts,    npoi, 283.16_r8)
             ts (i,j)= 283.16_r8
             !CALL const (tl,    npoi, 283.16_r8)
             tl (i,j)= 283.16_r8
          END DO
       END DO
       !      END IF    

       !
       ! initialize weather generator 'memory'
       !
       !CALL const (xstore, npoi*3, 0.0_r8)
       !      IF (irestart == 0) THEN

       DO k = 1,nband
          DO j = 1,jbMax 
             DO i = 1, nlpoints(j) 
                xstore(i,k,j)=0.0_r8
             END DO
          END DO
       END DO
       DO j = 1,jbMax 
          DO i = 1, nlpoints(j) 
             !
             ! initialize temperature of lower canopy buried by
             ! snow to 0 degrees C
             !
             !CALL const (tlsub, npoi, 273.16)
             tlsub(i,j) =273.16_r8
             !
             ! initialize canopy air conditions (used in turvap)
             !
             !CALL const (t12, npoi, 283.16) 
             t12(i,j) = 283.16_r8
             !CALL const (t34, npoi, 283.16) 
             t34(i,j) = 283.16_r8
             !
             !CALL const (q12, npoi, 0.0)    
             q12(i,j) = 0.0_r8
             !CALL const (q34, npoi, 0.0)    
             q34(i,j) = 0.0_r8
             !
             ! initialize all co2 concentrations (mol/mol)
             !
             !CALL const (ciub, npoi, 350.0e-06)
             ciub(i,j) = 350.0e-06_r8
             !CALL const (ciuc, npoi, 350.0e-06)
             ciuc(i,j) = 350.0e-06_r8
             !CALL const (cils, npoi, 350.0e-06)
             cils(i,j) = 350.0e-06_r8
             !CALL const (cil3, npoi, 350.0e-06)
             cil3(i,j) = 350.0e-06_r8
             !CALL const (cil4, npoi, 350.0e-06)
             cil4(i,j) = 350.0e-06_r8
             !
             !CALL const (csub, npoi, 350.0e-06)
             csub(i,j) = 350.0e-06_r8
             !CALL const (csuc, npoi, 350.0e-06)
             csuc(i,j) = 350.0e-06_r8
             !CALL const (csls, npoi, 350.0e-06)
             csls(i,j) = 350.0e-06_r8
             !CALL const (csl3, npoi, 350.0e-06)
             csl3(i,j) = 350.0e-06_r8
             !CALL const (csl4, npoi, 350.0e-06)
             csl4(i,j) = 350.0e-06_r8
             !
             ! initialize stomatal conductance (mol-h2o/m**2/sec)
             !       
             !CALL const (gsub, npoi, 0.5)
             gsub(i,j) =0.5_r8
             !CALL const (gsuc, npoi, 0.5)
             gsuc (i,j)=0.5_r8
             !CALL const (gsls, npoi, 0.5)
             gsls (i,j)=0.5_r8
             !CALL const (gsl3, npoi, 0.5)
             gsl3 (i,j)=0.5_r8
             !CALL const (gsl4, npoi, 0.5)
             gsl4 (i,j)=0.5_r8
          END DO
       END DO
       !      END IF
       !
       ! initialize soil biogeochemistry variables
       !
       !      IF (irestart == 0) THEN
       DO j = 1,jbMax 
          DO i = 1, nlpoints(j) 
             !CALL const (clitlm,   npoi, 0.0)
             clitlm(i,j)=    0.0_r8
             !CALL const (clitls,   npoi, 0.0)
             clitls(i,j)=    0.0_r8
             !CALL const (clitll,   npoi, 0.0)
             clitll(i,j)=    0.0_r8
             !CALL const (clitrm,   npoi, 0.0)
             clitrm(i,j)=    0.0_r8
             !CALL const (clitrs,   npoi, 0.0)
             clitrs(i,j)=    0.0_r8
             !CALL const (clitrl,   npoi, 0.0)
             clitrl(i,j)=    0.0_r8
             !CALL const (clitwm,   npoi, 0.0)
             clitwm(i,j)=    0.0_r8
             !CALL const (clitws,   npoi, 0.0)
             clitws(i,j)=    0.0_r8
             !CALL const (clitwl,   npoi, 0.0)
             clitwl(i,j)=    0.0_r8
             !CALL const (totcmic,  npoi, 0.0)
             totcmic(i,j)=   0.0_r8
             !CALL const (csoislop, npoi, 0.0)
             csoislop(i,j)=  0.0_r8
             !CALL const (csoislon, npoi, 0.0)
             csoislon(i,j)=  0.0_r8
             !CALL const (csoipas,  npoi, 0.0)
             csoipas(i,j)=   0.0_r8
          END DO
       END DO
       DO j = 1,jbMax 
          DO i = 1, nlpoints(j) 

             !CALL const (totlit,   npoi, 0.0) 
             totlit(i,j) =     0.0_r8
             !CALL const (totnlit,  npoi, 0.0) 
             totnlit(i,j) =    0.0_r8
             !CALL const (totfall,  npoi, 0.0) 
             totfall(i,j) =    0.0_r8
             !CALL const (totalit,  npoi, 0.0) 
             totalit(i,j) =    0.0_r8
             !CALL const (totrlit,  npoi, 0.0) 
             totrlit(i,j) =    0.0_r8
             !CALL const (totanlit, npoi, 0.0) 
             totanlit(i,j) =   0.0_r8
             !CALL const (totrnlit, npoi, 0.0) 
             totrnlit(i,j) =   0.0_r8
             !CALL const (totcsoi,  npoi, 0.0) 
             totcsoi(i,j) =    0.0_r8
             !CALL const (totnmic,  npoi, 0.0) 
             totnmic(i,j) =    0.0_r8
             !CALL const (tco2mic,  npoi, 0.0) 
             tco2mic(i,j) =    0.0_r8
             !CALL const (tnpptot,  npoi, 0.0) 
             tnpptot(i,j) =    0.0_r8
             !CALL const (tneetot,  npoi, 0.0) 
             tneetot(i,j) =    0.0_r8
             !CALL const (tnmin,    npoi, 0.0) 
             tnmin (i,j)=           0.0_r8
             !
             ! initialize carbon lost to atmosphere due
             ! to biomass burning
             !
             !CALL const (cdisturb, npoi, 0.0)
             cdisturb(i,j) = 0.0_r8
          END DO
       END DO
       !      END IF

       !
       ! initialize phenology flags
       !
       !      IF (irestart == 0) THEN
       DO j = 1,jbMax 
          DO i = 1, nlpoints(j) 
             !CALL const (tempu,  npoi, 1.0) 
             tempu(i,j) = 1.0_r8
             !CALL const (templ,  npoi, 1.0) 
             templ(i,j) = 1.0_r8
             !CALL const (dropu,  npoi, 1.0) 
             dropu(i,j) = 1.0_r8
             !CALL const (dropls, npoi, 1.0) 
             dropls(i,j) =1.0_r8
             !CALL const (dropl4, npoi, 1.0) 
             dropl4(i,j) =1.0_r8
             !CALL const (dropl3, npoi, 1.0) 
             dropl3(i,j) =1.0_r8
          END DO
       END DO
       !      END IF
       !
       ! initialize water and snow interception fractions
       !          
       !      IF (irestart == 0) THEN

       DO j = 1,jbMax 
          DO i = 1, nlpoints(j) 

             !CALL const (wliqu, npoi, 0.0) 
             wliqu(i,j) = 0.0_r8
             !CALL const (wliqs, npoi, 0.0) 
             wliqs(i,j) = 0.0_r8
             !CALL const (wliql, npoi, 0.0) 
             wliql(i,j) = 0.0_r8
             !
             !CALL const (wsnou, npoi, 0.0) 
             wsnou(i,j) = 0.0_r8
             !CALL(i,j) const (wsnos, npoi, 0.0) 
             wsnos(i,j) = 0.0_r8
             !CALL const (wsnol, npoi, 0.0) 
             wsnol(i,j) = 0.0_r8
             !
             !CALL const (su, npoi, 0.0)
             su (i,j) = 0.0_r8
             !CALL const (ss, npoi, 0.0)
             ss (i,j) = 0.0_r8
             !CALL const (sl, npoi, 0.0)
             sl (i,j)= 0.0_r8

          END DO
       END DO
    END IF

    !
    ! return to main program
    !
    RETURN
  END SUBROUTINE inisurf
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE inisum (wliqu          , &! INTENT(IN   ) :: wliqu  (npoi)   
       wsnou          , &! INTENT(IN   ) :: wsnou  (npoi)   
       fu             , &! INTENT(IN   ) :: fu     (npoi)   
       lai            , &! INTENT(IN   ) :: lai    (npoi,2) 
       wliqs          , &! INTENT(IN   ) :: wliqs  (npoi)   
       wsnos          , &! INTENT(IN   ) :: wsnos  (npoi)   
       sai            , &! INTENT(IN   ) :: sai    (npoi,2) 
       wliql          , &! INTENT(IN   ) :: wliql  (npoi)   
       wsnol          , &! INTENT(IN   ) :: wsnol  (npoi)   
       fl             , &! INTENT(IN   ) :: fl     (npoi)   
       decompl        , &! INTENT(OUT  ) :: decompl(npoi)   
       decomps        , &! INTENT(OUT  ) :: decomps(npoi)   
       firefac        , &! INTENT(OUT  ) :: firefac(npoi)   
       adrain         , &! INTENT(OUT  ) :: adrain        (npoi)     
       adsnow         , &! INTENT(OUT  ) :: adsnow        (npoi)     
       adaet          , &! INTENT(OUT  ) :: adaet        (npoi)     
       adtrunoff      , &! INTENT(OUT  ) :: adtrunoff (npoi)     
       adsrunoff      , &! INTENT(OUT  ) :: adsrunoff (npoi)     
       addrainage     , &! INTENT(OUT  ) :: addrainage(npoi)     
       adrh           , &! INTENT(OUT  ) :: adrh        (npoi)     
       adsnod         , &! INTENT(OUT  ) :: adsnod        (npoi)     
       adsnof         , &! INTENT(OUT  ) :: adsnof        (npoi)     
       adwsoi         , &! INTENT(OUT  ) :: adwsoi        (npoi)     
       adtsoi         , &! INTENT(OUT  ) :: adtsoi        (npoi)     
       adwisoi        , &! INTENT(OUT  ) :: adwisoi        (npoi)     
       adtlaysoi      , &! INTENT(OUT  ) :: adtlaysoi (npoi)     
       adwlaysoi      , &! INTENT(OUT  ) :: adwlaysoi (npoi)     
       adwsoic        , &! INTENT(OUT  ) :: adwsoic        (npoi)     
       adtsoic        , &! INTENT(OUT  ) :: adtsoic        (npoi)     
       adco2mic       , &! INTENT(OUT  ) :: adco2mic  (npoi)     
       adco2root      , &! INTENT(OUT  ) :: adco2root (npoi)     
       adnmintot      , &! INTENT(OUT  ) :: adnmintot (npoi)     
       amtemp         , &! INTENT(OUT  ) :: amtemp        (npoi)     
       amrain         , &! INTENT(OUT  ) :: amrain        (npoi)     
       amsnow         , &! INTENT(OUT  ) :: amsnow        (npoi)     
       amaet          , &! INTENT(OUT  ) :: amaet        (npoi)     
       amtrunoff      , &! INTENT(OUT  ) :: amtrunoff (npoi)     
       amsrunoff      , &! INTENT(OUT  ) :: amsrunoff (npoi)     
       amdrainage     , &! INTENT(OUT  ) :: amdrainage(npoi)     
       amqa           , &! INTENT(OUT  ) :: amqa        (npoi)     
       amsolar        , &! INTENT(OUT  ) :: amsolar        (npoi)     
       amirup         , &! INTENT(OUT  ) :: amirup        (npoi)     
       amirdown       , &! INTENT(OUT  ) :: amirdown  (npoi)     
       amsens         , &! INTENT(OUT  ) :: amsens        (npoi)     
       amlatent       , &! INTENT(OUT  ) :: amlatent  (npoi)     
       amlaiu         , &! INTENT(OUT  ) :: amlaiu        (npoi)     
       amlail         , &! INTENT(OUT  ) :: amlail        (npoi)     
       amtsoi         , &! INTENT(OUT  ) :: amtsoi        (npoi)     
       amwsoi         , &! INTENT(OUT  ) :: amwsoi        (npoi)     
       amwisoi        , &! INTENT(OUT  ) :: amwisoi        (npoi)     
       amvwc          , &! INTENT(OUT  ) :: amvwc        (npoi)     
       amawc          , &! INTENT(OUT  ) :: amawc        (npoi)     
       amsnod         , &! INTENT(OUT  ) :: amsnod        (npoi)     
       amsnof         , &! INTENT(OUT  ) :: amsnof        (npoi)     
       amco2mic       , &! INTENT(OUT  ) :: amco2mic  (npoi)     
       amco2root      , &! INTENT(OUT  ) :: amco2root (npoi)     
       amnmintot      , &! INTENT(OUT  ) :: amnmintot (npoi)     
       amnpp          , &! INTENT(OUT  ) :: amnpp        (npoi,npft)
       aysolar        , &! INTENT(OUT  ) :: aysolar        (npoi)     
       ayirup         , &! INTENT(OUT  ) :: ayirup        (npoi)     
       ayirdown       , &! INTENT(OUT  ) :: ayirdown  (npoi)     
       aysens         , &! INTENT(OUT  ) :: aysens        (npoi)     
       aylatent       , &! INTENT(OUT  ) :: aylatent  (npoi)     
       ayprcp         , &! INTENT(OUT  ) :: ayprcp        (npoi)     
       ayaet          , &! INTENT(OUT  ) :: ayaet        (npoi)     
       aytrans        , &! INTENT(OUT  ) :: aytrans        (npoi)     
       aytrunoff      , &! INTENT(OUT  ) :: aytrunoff (npoi)     
       aysrunoff      , &! INTENT(OUT  ) :: aysrunoff (npoi)     
       aydrainage     , &! INTENT(OUT  ) :: aydrainage(npoi)     
       aywsoi         , &! INTENT(OUT  ) :: aywsoi        (npoi)     
       aywisoi        , &! INTENT(OUT  ) :: aywisoi        (npoi)     
       aytsoi         , &! INTENT(OUT  ) :: aytsoi        (npoi)     
       ayvwc          , &! INTENT(OUT  ) :: ayvwc        (npoi)     
       ayawc          , &! INTENT(OUT  ) :: ayawc        (npoi)     
       aystresstu     , &! INTENT(OUT  ) :: aystresstu(npoi)     
       aystresstl     , &! INTENT(OUT  ) :: aystresstl(npoi)     
       ayco2mic       , &! INTENT(OUT  ) :: ayco2mic  (npoi)     
       ayco2root      , &! INTENT(OUT  ) :: ayco2root (npoi)     
       ayrootbio      , &! INTENT(OUT  ) :: ayrootbio (npoi)     
       aynmintot      , &! INTENT(OUT  ) :: aynmintot (npoi)     
       ayalit         , &! INTENT(OUT  ) :: ayalit        (npoi)     
       ayblit         , &! INTENT(OUT  ) :: ayblit        (npoi)     
       aycsoi         , &! INTENT(OUT  ) :: aycsoi        (npoi)     
       aycmic         , &! INTENT(OUT  ) :: aycmic        (npoi)     
       ayanlit        , &! INTENT(OUT  ) :: ayanlit        (npoi)     
       aybnlit        , &! INTENT(OUT  ) :: aybnlit        (npoi)     
       aynsoi         , &! INTENT(OUT  ) :: aynsoi        (npoi)     
       aygpp          , &! INTENT(OUT  ) :: aygpp        (npoi,npft)
       wpud           , &! INTENT(IN   )
       wipud          , &! INTENT(IN   )
       poros                , &! INTENT(IN   )
       wsoi           , &! INTENT(IN   )
       wisoi          , &! INTENT(IN   )
       hsoi           , &! INTENT(IN   )
       fi             , &! INTENT(IN   )
       rhos           , &! INTENT(IN   )
       hsno           , &! INTENT(IN   )
       wtot           , &! INTENT(OUT  )
       nsoilay        , &! INTENT(IN   )
       nsnolay        , &! INTENT(IN   )
       npft           , &! INTENT(IN   )
       rhow             )! INTENT(IN   )
    ! ---------------------------------------------------------------------
    ! 
    ! does initialization for time averaging
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN   ) :: nsoilay ! number of soil layers
    INTEGER, INTENT(IN   ) :: nsnolay ! number of snow layers
    INTEGER, INTENT(IN   ) :: npft         ! number of plant functional types
    REAL(KIND=r8)   , INTENT(IN   ) :: rhow    ! density of liquid water (all types) (kg m-3)

    REAL(KIND=r8)   , INTENT(OUT  ) :: wtot  (ibMax,jbMax)     ! total amount of water stored in snow, soil,
    ! puddels, and on vegetation (kg_h2o)

    REAL(KIND=r8)   , INTENT(IN   ) :: fi    (ibMax,jbMax)       ! fractional snow cover
    REAL(KIND=r8)   , INTENT(IN   ) :: rhos                       ! density of snow (kg m-3)
    REAL(KIND=r8)   , INTENT(IN   ) :: hsno  (ibMax,nsnolay,jbMax)! thickness of snow layers (m)

    REAL(KIND=r8)   , INTENT(IN   ) :: wpud   (ibMax,jbMax)                  ! liquid content of puddles per soil area (kg m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: wipud  (ibMax,jbMax)                  ! ice content of puddles per soil area (kg m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: poros  (ibMax,nsoilay,jbMax)! porosity (mass of h2o per unit vol at sat / rhow)
    REAL(KIND=r8)   , INTENT(IN   ) :: wsoi   (ibMax,nsoilay,jbMax)! fraction of soil pore space containing liquid water
    REAL(KIND=r8)   , INTENT(IN   ) :: wisoi  (ibMax,nsoilay,jbMax)! fraction of soil pore space containing ice
    REAL(KIND=r8)   , INTENT(IN   ) :: hsoi   (nsoilay+1)          ! soil layer thickness (m)

    REAL(KIND=r8)   , INTENT(OUT  ) :: adrain   (ibMax,jbMax)! daily average rainfall rate (mm/day)
    REAL(KIND=r8)   , INTENT(OUT  ) :: adsnow   (ibMax,jbMax)! daily average snowfall rate (mm/day)
    REAL(KIND=r8)   , INTENT(OUT  ) :: adaet    (ibMax,jbMax)! daily average aet (mm/day)
    REAL(KIND=r8)   , INTENT(OUT  ) :: adtrunoff(ibMax,jbMax)! daily average total runoff (mm/day)
    REAL(KIND=r8)   , INTENT(OUT  ) :: adsrunoff(ibMax,jbMax)! daily average surface runoff (mm/day)
    REAL(KIND=r8)   , INTENT(OUT  ) :: addrainage(ibMax,jbMax)! daily average drainage (mm/day)
    REAL(KIND=r8)   , INTENT(OUT  ) :: adrh      (ibMax,jbMax)! daily average rh (percent)
    REAL(KIND=r8)   , INTENT(OUT  ) :: adsnod    (ibMax,jbMax)! daily average snow depth (m)
    REAL(KIND=r8)   , INTENT(OUT  ) :: adsnof    (ibMax,jbMax)! daily average snow fraction (fraction)
    REAL(KIND=r8)   , INTENT(OUT  ) :: adwsoi    (ibMax,jbMax)! daily average soil moisture (fraction)
    REAL(KIND=r8)   , INTENT(OUT  ) :: adtsoi    (ibMax,jbMax)! daily average soil temperature (c)
    REAL(KIND=r8)   , INTENT(OUT  ) :: adwisoi   (ibMax,jbMax)! daily average soil ice (fraction)
    REAL(KIND=r8)   , INTENT(OUT  ) :: adtlaysoi (ibMax,jbMax)! daily average soil temperature (c) of top layer
    REAL(KIND=r8)   , INTENT(OUT  ) :: adwlaysoi (ibMax,jbMax)! daily average soil moisture of top layer(fraction)
    REAL(KIND=r8)   , INTENT(OUT  ) :: adwsoic   (ibMax,jbMax)! daily average soil moisture using root profile weighting (fraction)
    REAL(KIND=r8)   , INTENT(OUT  ) :: adtsoic   (ibMax,jbMax)! daily average soil temperature (c) using profile weighting
    REAL(KIND=r8)   , INTENT(OUT  ) :: adco2mic  (ibMax,jbMax)! daily accumulated co2 respiration from microbes (kg_C m-2 /day)
    REAL(KIND=r8)   , INTENT(OUT  ) :: adco2root (ibMax,jbMax)! daily accumulated co2 respiration from roots (kg_C m-2 /day)
    REAL(KIND=r8)   , INTENT(OUT  ) :: adnmintot (ibMax,jbMax)! daily accumulated net nitrogen mineralization (kg_N m-2 /day)
    REAL(KIND=r8)   , INTENT(OUT  ) :: amtemp    (ibMax,jbMax)! monthly average air temperature (C)
    REAL(KIND=r8)   , INTENT(OUT  ) :: amrain    (ibMax,jbMax)! monthly average rainfall rate (mm/day)
    REAL(KIND=r8)   , INTENT(OUT  ) :: amsnow    (ibMax,jbMax)! monthly average snowfall rate (mm/day)
    REAL(KIND=r8)   , INTENT(OUT  ) :: amaet     (ibMax,jbMax)! monthly average aet (mm/day)
    REAL(KIND=r8)   , INTENT(OUT  ) :: amtrunoff (ibMax,jbMax)! monthly average total runoff (mm/day)
    REAL(KIND=r8)   , INTENT(OUT  ) :: amsrunoff (ibMax,jbMax)! monthly average surface runoff (mm/day)
    REAL(KIND=r8)   , INTENT(OUT  ) :: amdrainage(ibMax,jbMax)! monthly average drainage (mm/day)
    REAL(KIND=r8)   , INTENT(OUT  ) :: amqa      (ibMax,jbMax)! monthly average specific humidity (kg-h2o/kg-air)
    REAL(KIND=r8)   , INTENT(OUT  ) :: amsolar   (ibMax,jbMax)! monthly average incident solar radiation (W/m**2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: amirup    (ibMax,jbMax)! monthly average upward ir radiation (W/m**2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: amirdown  (ibMax,jbMax)! monthly average downward ir radiation (W/m**2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: amsens    (ibMax,jbMax)! monthly average sensible heat flux (W/m**2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: amlatent  (ibMax,jbMax)! monthly average latent heat flux (W/m**2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: amlaiu    (ibMax,jbMax)! monthly average lai for upper canopy (m**2/m**2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: amlail    (ibMax,jbMax)! monthly average lai for lower canopy (m**2/m**2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: amtsoi    (ibMax,jbMax)! monthly average 1m soil temperature (C)
    REAL(KIND=r8)   , INTENT(OUT  ) :: amwsoi    (ibMax,jbMax)! monthly average 1m soil moisture (fraction)
    REAL(KIND=r8)   , INTENT(OUT  ) :: amwisoi   (ibMax,jbMax)! monthly average 1m soil ice (fraction)
    REAL(KIND=r8)   , INTENT(OUT  ) :: amvwc     (ibMax,jbMax)! monthly average 1m volumetric water content (fraction)
    REAL(KIND=r8)   , INTENT(OUT  ) :: amawc     (ibMax,jbMax)! monthly average 1m plant-available water content (fraction)
    REAL(KIND=r8)   , INTENT(OUT  ) :: amsnod    (ibMax,jbMax)! monthly average snow depth (m)
    REAL(KIND=r8)   , INTENT(OUT  ) :: amsnof    (ibMax,jbMax)! monthly average snow fraction (fraction)
    REAL(KIND=r8)   , INTENT(OUT  ) :: amco2mic  (ibMax,jbMax)! monthly total CO2 flux from microbial respiration (kg-C/m**2/month)
    REAL(KIND=r8)   , INTENT(OUT  ) :: amco2root (ibMax,jbMax)! monthly total CO2 flux from soil due to root respiration (kg-C/m**2/month)
    REAL(KIND=r8)   , INTENT(OUT  ) :: amnmintot (ibMax,jbMax)! monthly total N mineralization from microbes (kg-N/m**2/month)
    REAL(KIND=r8)   , INTENT(OUT  ) :: amnpp     (ibMax,npft,jbMax)! monthly total npp for each plant type (kg-C/m**2/month)
    REAL(KIND=r8)   , INTENT(OUT  ) :: aysolar   (ibMax,jbMax)! annual average incident solar radiation (w/m**2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ayirup    (ibMax,jbMax)! annual average upward ir radiation (w/m**2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ayirdown  (ibMax,jbMax)! annual average downward ir radiation (w/m**2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: aysens    (ibMax,jbMax)! annual average sensible heat flux (w/m**2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: aylatent  (ibMax,jbMax)! annual average latent heat flux (w/m**2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ayprcp    (ibMax,jbMax)! annual average precipitation (mm/yr)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ayaet     (ibMax,jbMax)! annual average aet (mm/yr)
    REAL(KIND=r8)   , INTENT(OUT  ) :: aytrans   (ibMax,jbMax)! annual average transpiration (mm/yr)
    REAL(KIND=r8)   , INTENT(OUT  ) :: aytrunoff (ibMax,jbMax)! annual average total runoff (mm/yr)
    REAL(KIND=r8)   , INTENT(OUT  ) :: aysrunoff (ibMax,jbMax)! annual average surface runoff (mm/yr)
    REAL(KIND=r8)   , INTENT(OUT  ) :: aydrainage(ibMax,jbMax)! annual average drainage (mm/yr)
    REAL(KIND=r8)   , INTENT(OUT  ) :: aywsoi    (ibMax,jbMax)! annual average 1m soil moisture (fraction)
    REAL(KIND=r8)   , INTENT(OUT  ) :: aywisoi   (ibMax,jbMax)! annual average 1m soil ice (fraction)
    REAL(KIND=r8)   , INTENT(OUT  ) :: aytsoi    (ibMax,jbMax)! annual average 1m soil temperature (C)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ayvwc     (ibMax,jbMax)! annual average 1m volumetric water content (fraction)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ayawc     (ibMax,jbMax)! annual average 1m plant-available water content (fraction)
    REAL(KIND=r8)   , INTENT(OUT  ) :: aystresstu(ibMax,jbMax)! annual average soil moisture stress parameter for
    ! upper canopy (dimensionless)
    REAL(KIND=r8)   , INTENT(OUT  ) :: aystresstl(ibMax,jbMax)! annual average soil moisture stress parameter for]
    ! lower canopy (dimensionless)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ayco2mic  (ibMax,jbMax)! annual total CO2 flux from microbial respiration
    ! (kg-C/m**2/yr)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ayco2root (ibMax,jbMax)! annual total CO2 flux from soil due to root respiration (kg-C/m**2/yr)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ayrootbio (ibMax,jbMax)! annual average live root biomass (kg-C / m**2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: aynmintot (ibMax,jbMax)! annual total nitrogen mineralization (kg-N/m**2/yr)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ayalit    (ibMax,jbMax)! aboveground litter (kg-c/m**2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ayblit    (ibMax,jbMax)! belowground litter (kg-c/m**2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: aycsoi    (ibMax,jbMax)! total soil carbon (kg-c/m**2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: aycmic    (ibMax,jbMax)! total soil carbon in microbial biomass (kg-c/m**2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ayanlit   (ibMax,jbMax)! aboveground litter nitrogen (kg-N/m**2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: aybnlit   (ibMax,jbMax)! belowground litter nitrogen (kg-N/m**2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: aynsoi    (ibMax,jbMax)! total soil nitrogen (kg-N/m**2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: aygpp     (ibMax,npft,jbMax)! annual gross npp for each plant type(kg-c/m**2/yr)

    REAL(KIND=r8)   , INTENT(IN   ) :: wliqu     (ibMax,jbMax)! intercepted liquid h2o on upper canopy leaf area (kg m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: wsnou     (ibMax,jbMax)! intercepted frozen h2o (snow) on upper canopy leaf area (kg m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: fu     (ibMax,jbMax)! fraction of overall area covered by upper canopy
    REAL(KIND=r8)   , INTENT(IN   ) :: lai    (ibMax,2,jbMax) ! canopy single-sided leaf area index (area leaf/area veg)
    REAL(KIND=r8)   , INTENT(IN   ) :: wliqs  (ibMax,jbMax)! intercepted liquid h2o on upper canopy stem area (kg m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: wsnos  (ibMax,jbMax)! intercepted frozen h2o (snow) on upper canopy stem area (kg m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: sai    (ibMax,2,jbMax) ! current single-sided stem area index
    REAL(KIND=r8)   , INTENT(IN   ) :: wliql  (ibMax,jbMax)! intercepted liquid h2o on lower canopy leaf and stem area (kg m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: wsnol  (ibMax,jbMax)! intercepted frozen h2o (snow) on lower canopy leaf & stem area (kg m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: fl     (ibMax,jbMax)! fraction of snow-free area covered by lower  canopy
    REAL(KIND=r8)   , INTENT(OUT  ) :: decompl(ibMax,jbMax)           ! litter decomposition factor                  (dimensionless)
    REAL(KIND=r8)   , INTENT(OUT  ) :: decomps(ibMax,jbMax)           ! soil organic matter decomposition factor          (dimensionless)
    REAL(KIND=r8)   , INTENT(OUT  ) :: firefac(ibMax,jbMax)           ! factor that respresents the annual average fuel]
    ! dryness of a grid cell, and hence characterizes the readiness to burn

    !
    ! local variables
    !
    INTEGER :: i        ! loop indices
    INTEGER :: j        ! loop indices
    INTEGER :: k        ! loop indices

    !
    ! initialize total water content in soil+snow+vegetation (for mass conservation check)
    !
    DO j = 1,jbMax 
       DO i = 1, nlpoints(j) 

          wtot(i,j) = (wliqu(i,j)+wsnou(i,j)) * fu(i,j) * 2.0_r8 * lai(i,2,j) +  &
               (wliqs(i,j)+wsnos(i,j)) * fu(i,j) * 2.0_r8 * sai(i,2,j) +  &
               (wliql(i,j)+wsnol(i,j)) * fl(i,j) * 2.0_r8 *             &
               (lai(i,1,j) + sai(i,1,j)) * (1.0_r8 - fi(i,j))
          !
          wtot(i,j) = wtot(i,j) + wpud(i,j) + wipud(i,j)
          !
          DO k = 1, nsoilay
             wtot(i,j) = wtot(i,j) +                      &
                  poros(i,k,j)*wsoi(i,k,j)*(1.0_r8-wisoi(i,k,j))*hsoi(k)*rhow +  &
                  poros(i,k,j)*wisoi(i,k,j)*hsoi(k)*rhow
          END DO
          !
          DO k = 1, nsnolay
             wtot(i,j) = wtot(i,j) + fi(i,j)*rhos*hsno(i,k,j)
          END DO
       END DO
    END DO
    !
    ! Daily means
    !
    DO j = 1,jbMax 
       DO i = 1, nlpoints(j) 

          !CALL const (adrain,     npoi, 0.0) 
          adrain(i,j) = 0.0_r8
          !CALL const (adsnow,     npoi, 0.0) 
          adsnow(i,j)=0.0_r8
          !CALL const (adaet,      npoi, 0.0) 
          adaet(i,j)=0.0_r8
          !CALL const (adtrunoff,  npoi, 0.0)
          adtrunoff(i,j)=0.0_r8
          !CALL const (adsrunoff,  npoi, 0.0)
          adsrunoff(i,j)=0.0_r8
          !CALL const (addrainage, npoi, 0.0)
          addrainage(i,j)=0.0_r8
          !CALL const (adrh,       npoi, 0.0)
          adrh(i,j)=0.0_r8
          !CALL const (adsnod,     npoi, 0.0)
          adsnod(i,j)=0.0_r8
          !CALL const (adsnof,     npoi, 0.0)
          adsnof(i,j)=0.0_r8
          !CALL const (adwsoi,     npoi, 0.0)
          adwsoi(i,j)=0.0_r8
          !CALL const (adtsoi,     npoi, 0.0)
          adtsoi(i,j)=0.0_r8
          !CALL const (adwisoi,    npoi, 0.0)
          adwisoi(i,j)=0.0_r8
          !CALL const (adtlaysoi,  npoi, 0.0)
          adtlaysoi(i,j)=0.0_r8
          !CALL const (adwlaysoi,  npoi, 0.0)
          adwlaysoi(i,j)=0.0_r8
          !CALL const (adwsoic,    npoi, 0.0)
          adwsoic(i,j)=0.0_r8
          !CALL const (adtsoic,    npoi, 0.0)
          adtsoic(i,j)=0.0_r8
          !CALL const (adco2mic,   npoi, 0.0)
          adco2mic(i,j)=0.0_r8
          !CALL const (adco2root,  npoi, 0.0)
          adco2root(i,j)=0.0_r8
          !CALL const (decompl,    npoi, 0.0)
          decompl(i,j)=0.0_r8
          !CALL const (decomps,    npoi, 0.0)
          decomps(i,j)=0.0_r8
          !CALL const (adnmintot,  npoi, 0.0)
          adnmintot(i,j)=0.0_r8 
       END DO
    END DO

    !
    ! monthly mean quanties
    !
    DO j = 1,jbMax 
       DO i = 1, nlpoints(j) 

          !CALL const (amtemp,     npoi, 0.0)
          amtemp(i,j) = 0.0_r8
          !CALL const (amrain,     npoi, 0.0)
          amrain(i,j) =0.0_r8   
          !CALL const (amsnow,     npoi, 0.0_r8)
          amsnow(i,j) =0.0_r8      
          !CALL const (amaet,      npoi, 0.0_r8)
          amaet (i,j) =0.0_r8     
          !CALL const (amtrunoff,  npoi, 0.0_r8)
          amtrunoff(i,j) =0.0_r8   
          !CALL const (amsrunoff,  npoi, 0.0_r8)
          amsrunoff(i,j) =0.0_r8   
          !CALL const (amdrainage, npoi, 0.0_r8)
          amdrainage(i,j) =0.0_r8  
          !CALL const (amqa,       npoi, 0.0_r8)
          amqa(i,j) =0.0_r8     
          !CALL const (amsolar,    npoi, 0.0_r8)
          amsolar(i,j) =0.0_r8  
          !CALL const (amirup,     npoi, 0.0_r8)
          amirup(i,j) =0.0_r8    
          !CALL const (amirdown,   npoi, 0.0_r8)
          amirdown(i,j) =0.0_r8   
          !CALL const (amsens,     npoi, 0.0_r8)
          amsens(i,j) =0.0_r8     
          !CALL const (amlatent,   npoi, 0.0_r8)
          amlatent(i,j) =0.0_r8   
          !CALL const (amlaiu,     npoi, 0.0_r8)
          amlaiu(i,j) =0.0_r8  
          !CALL const (amlail,     npoi, 0.0_r8)
          amlail(i,j) =0.0_r8  
          !CALL const (amtsoi,     npoi, 0.0_r8)
          amtsoi(i,j) =0.0_r8  
          !CALL const (amwsoi,     npoi, 0.0_r8)
          amwsoi(i,j) =0.0_r8  
          !CALL const (amwisoi,    npoi, 0.0_r8)
          amwisoi(i,j) =0.0_r8  
          !CALL const (amvwc,      npoi, 0.0_r8)
          amvwc(i,j) =0.0_r8  
          !CALL const (amawc,      npoi, 0.0_r8)
          amawc(i,j) =0.0_r8    
          !CALL const (amsnod,     npoi, 0.0_r8)
          amsnod(i,j) =0.0_r8  
          !CALL const (amsnof,     npoi, 0.0_r8)
          amsnof(i,j) =0.0_r8  
          !CALL const (amco2mic,   npoi, 0.0_r8)
          amco2mic(i,j) =0.0_r8  
          !CALL const (amco2root,  npoi, 0.0_r8)
          amco2root(i,j) =0.0_r8  
          !CALL const (amnmintot,  npoi, 0.0_r8)
          amnmintot(i,j) =0.0_r8  
       END DO
    END DO
    DO k = 1,npft
       DO j = 1,jbMax 
          DO i = 1, nlpoints(j) 
             !CALL const (amnpp,      npoi*npft, 0.0_r8)
             amnpp(i,k,j)=0.0_r8    
          END DO
       END DO
    END DO
    !
    ! Annual mean quantities
    !
    DO j = 1,jbMax 
       DO i = 1, nlpoints(j) 
          !CALL const (aysolar,    npoi, 0.0_r8) 
          aysolar (i,j)     =  0.0_r8
          !CALL const (ayirup,     npoi, 0.0_r8) 
          ayirup  (i,j)     =  0.0_r8
          !CALL const (ayirdown,   npoi, 0.0_r8) 
          ayirdown(i,j)     =  0.0_r8
          !CALL const (aysens,     npoi, 0.0_r8) 
          aysens (i,j)      =  0.0_r8
          !CALL const (aylatent,   npoi, 0.0_r8) 
          aylatent (i,j)    =  0.0_r8
          !CALL const (ayprcp,     npoi, 0.0_r8) 
          ayprcp (i,j)      =  0.0_r8
          !CALL const (ayaet,      npoi, 0.0_r8) 
          ayaet  (i,j)      =  0.0_r8
          !CALL const (aytrans,    npoi, 0.0_r8) 
          aytrans (i,j)     =  0.0_r8
          !CALL const (aytrunoff,  npoi, 0.0_r8) 
          aytrunoff(i,j)    =  0.0_r8
          !CALL const (aysrunoff,  npoi, 0.0_r8) 
          aysrunoff(i,j)    =  0.0_r8
          !CALL const (aydrainage, npoi, 0.0_r8) 
          aydrainage(i,j)   =  0.0_r8
          !CALL const (aywsoi,     npoi, 0.0_r8) 
          aywsoi  (i,j)     =  0.0_r8
          !CALL const (aywisoi,    npoi, 0.0_r8) 
          aywisoi (i,j)     =  0.0_r8
          !CALL const (aytsoi,     npoi, 0.0_r8) 
          aytsoi (i,j)      =  0.0_r8
          !CALL const (ayvwc,      npoi, 0.0_r8) 
          ayvwc  (i,j)      =  0.0_r8
          !CALL const (ayawc,      npoi, 0.0_r8) 
          ayawc (i,j)       =  0.0_r8
          !CALL const (aystresstu, npoi, 0.0_r8) 
          aystresstu(i,j)   =  0.0_r8
          !CALL const (aystresstl, npoi, 0.0_r8) 
          aystresstl(i,j)   =  0.0_r8
          !CALL const (firefac,    npoi, 0.0_r8) 
          firefac  (i,j)    =  0.0_r8
          !CALL const (firefac,    npoi, 0.0_r8) 
          firefac(i,j)      =  0.0_r8
          !CALL const (ayco2mic,   npoi, 0.0_r8) 
          ayco2mic(i,j)     =  0.0_r8
          !CALL const (ayco2root,  npoi, 0.0_r8) 
          ayco2root(i,j)    =  0.0_r8
          !CALL const (ayrootbio,  npoi, 0.0_r8) 
          ayrootbio(i,j)    =  0.0_r8
          !CALL const (aynmintot,  npoi, 0.0_r8) 
          aynmintot(i,j)    =  0.0_r8
          !CALL const (ayalit,     npoi, 0.0_r8) 
          ayalit (i,j)      =  0.0_r8
          !CALL const (ayblit,     npoi, 0.0_r8) 
          ayblit (i,j)      =  0.0_r8
          !CALL const (aycsoi,     npoi, 0.0_r8) 
          aycsoi (i,j)      =  0.0_r8
          !CALL const (aycmic,     npoi, 0.0_r8) 
          aycmic(i,j)       =  0.0_r8
          !CALL const (ayanlit,    npoi, 0.0_r8) 
          ayanlit(i,j)      =  0.0_r8
          !CALL const (aybnlit,    npoi, 0.0_r8) 
          aybnlit(i,j)      =  0.0_r8
          !CALL const (aynsoi,     npoi, 0.0_r8) 
          aynsoi (i,j)     =  0.0_r8
       END DO
    END DO
    DO k = 1,npft
       DO j = 1,jbMax 
          DO i = 1, nlpoints(j) 
             !CALL const (amnpp,      npoi*npft, 0.0_r8)
             aygpp(i,k,j)=0.0_r8    
          END DO
       END DO
    END DO

    RETURN
  END SUBROUTINE inisum
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE inisnow(rhow   , & ! INTENT(IN   ) :: rhow    
       nsnolay, & ! INTENT(IN   ) :: nsnolay 
       dtime  , & ! INTENT(IN   ) :: dtime   
       rhos   , & ! INTENT(OUT  ) :: rhos    
       consno , & ! INTENT(OUT  ) :: consno  
       hsnotop, & ! INTENT(OUT  ) :: hsnotop 
       hsnomin, & ! INTENT(OUT  ) :: hsnomin 
       fimin  , & ! INTENT(OUT  ) :: fimin   
       fimax  , & ! INTENT(OUT  ) :: fimax   
       z0sno    ) ! INTENT(OUT  ) :: z0sno   
    ! ---------------------------------------------------------------------
    !
    ! does initialization for snow model
    !
    IMPLICIT NONE
    !
    REAL(KIND=r8)   , INTENT(IN   ) :: rhow    ! density of liquid water (all types) (kg m-3)
    INTEGER, INTENT(IN   ) :: nsnolay ! number of snow layers
    REAL(KIND=r8)   , INTENT(IN   ) :: dtime   ! model timestep (seconds)
    REAL(KIND=r8)   , INTENT(OUT  ) :: rhos         ! density of snow (kg m-3)
    REAL(KIND=r8)   , INTENT(OUT  ) :: consno              ! thermal conductivity of snow (W m-1 K-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: hsnotop      ! thickness of top snow layer (m)
    REAL(KIND=r8)   , INTENT(OUT  ) :: hsnomin      ! minimum total thickness of snow (m)
    REAL(KIND=r8)   , INTENT(OUT  ) :: fimin               ! minimum fractional snow cover
    REAL(KIND=r8)   , INTENT(OUT  ) :: fimax               ! maximum fractional snow cover
    REAL(KIND=r8)   , INTENT(OUT  ) :: z0sno               ! roughness length of snow surface (m)

    !
    ! rhos is density of snow
    !
    rhos = 0.25_r8 * rhow
    !
    ! consno is thermal conductivity of snow
    !
    consno = 0.20_r8
    !
    ! hsnotop is "adaptive-grid" thickness of top snow layer
    !
    hsnotop = 0.05_r8
    !
    ! hsnomin is minimum total snow thickness. total thickness
    ! is constrained to hsnomin for less than 100% cover. (hsnomin
    ! should be ge nsnolay*hsnotop for vadapt to work properly.)
    !
    hsnomin = MAX (0.15_r8, nsnolay * hsnotop)
    !
    ! fimin and fimax are minimum and maximum snowcover fractions
    !
    fimin = 0.00002_r8 * (2.0_r8*dtime / 1800.0_r8) * (0.15_r8 / hsnomin)
    fimax = 1.000_r8
    !
    ! z0sno is roughness lenth of snow cover
    !
    z0sno = 0.0005_r8
    !
    RETURN
  END SUBROUTINE inisnow
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE inisoil (irestart  , &! INTENT(IN   )
       texdat    , &! INTENT(IN   )
       porosdat  , &! INTENT(IN   )
       sfielddat , &! INTENT(IN   )
       swiltdat  , &! INTENT(IN   )
       bexdat    , &! INTENT(IN   )
       suctiondat, &! INTENT(IN   )
       hydrauldat, &! INTENT(IN   )
       ndat            , &! INTENT(IN   )
       wpud      , &! INTENT(OUT  )
       wipud     , &! INTENT(OUT  )
       z0soi     , &! INTENT(OUT  )
       wsoi            , &! INTENT(INOUT)
       wsoim            , &! INTENT(INOUT)
       wisoi     , &! INTENT(INOUT)
       tsoi      , &! INTENT(INOUT)
       tsoim     , &! INTENT(INOUT)
       tsno      , &! INTENT(INOUT)
       tg        , &! INTENT(OUT  )
       tgm       , &! INTENT(OUT  )
       ti            , &! INTENT(OUT  )
       sand      , &! INTENT(IN   )
       clay      , &! INTENT(IN   )
       albsav    , &! INTENT(OUT  )
       albsan    , &! INTENT(OUT  )
       rhosoi    , &! INTENT(OUT  )
       csoi      , &! INTENT(OUT  )
       poros     , &! INTENT(OUT  )
       sfield    , &! INTENT(OUT  )
       swilt     , &! INTENT(OUT  )
       bex            , &! INTENT(OUT  )
       ibex      , &! INTENT(OUT  )
       suction   , &! INTENT(OUT  )
       hydraul   , &! INTENT(OUT  )
       nband     , &! INTENT(IN   )
       nsoilay   , &! INTENT(IN   )
       nsnolay     )! INTENT(IN   )
    ! ---------------------------------------------------------------------
    !
    ! does initialization for soil database
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN   ) :: irestart! 0 = initial run, 1 = restart run
    INTEGER, INTENT(IN   ) :: nband   ! number of solar radiation wavebands
    INTEGER, INTENT(IN   ) :: nsoilay ! number of soil layers
    INTEGER, INTENT(IN   ) :: nsnolay ! number of snow layers

    REAL(KIND=r8)   , INTENT(OUT  ) :: wpud    (ibMax,jbMax) ! liquid content of puddles per soil area (kg m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: wipud   (ibMax,jbMax) ! ice content of puddles per soil area (kg m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: z0soi   (ibMax,jbMax) ! roughness length of soil surface (m)
    REAL(KIND=r8)   , INTENT(INOUT) :: wsoi    (ibMax,nsoilay,jbMax)! fraction of soil pore space containing liquid water
    REAL(KIND=r8)   , INTENT(INOUT) :: wsoim    (ibMax,nsoilay,jbMax)! fraction of soil pore space containing liquid water
    REAL(KIND=r8)   , INTENT(INOUT) :: wisoi   (ibMax,nsoilay,jbMax)! fraction of soil pore space containing ice
    REAL(KIND=r8)   , INTENT(INOUT) :: tsoi    (ibMax,nsoilay,jbMax)! soil temperature for each layer (K)
    REAL(KIND=r8)   , INTENT(INOUT) :: tsoim   (ibMax,nsoilay,jbMax)! soil temperature for each layer (K)    
    REAL(KIND=r8)   , INTENT(INOUT) :: tsno    (ibMax,nsnolay,jbMax) ! temperature of snow layers (K)
    REAL(KIND=r8)   , INTENT(OUT  ) :: tg      (ibMax,jbMax)  ! soil skin temperature (K)
    REAL(KIND=r8)   , INTENT(OUT  ) :: tgm     (ibMax,jbMax)  ! soil skin temperature (K)      
    REAL(KIND=r8)   , INTENT(OUT  ) :: ti      (ibMax,jbMax)  ! snow skin temperature (K)
    REAL(KIND=r8)   , INTENT(IN   ) :: sand    (ibMax,nsoilay,jbMax)  ! percent sand of soil
    REAL(KIND=r8)   , INTENT(IN   ) :: clay    (ibMax,nsoilay,jbMax)  ! percent clay of soil
    REAL(KIND=r8)   , INTENT(OUT  ) :: albsav  (ibMax,jbMax)    ! saturated soil surface albedo (visible waveband)
    REAL(KIND=r8)   , INTENT(OUT  ) :: albsan  (ibMax,jbMax)    ! saturated soil surface albedo (near-ir waveband)
    REAL(KIND=r8)   , INTENT(OUT  ) :: rhosoi  (ibMax,nsoilay,jbMax)! soil density (without pores, not bulk) (kg m-3)
    REAL(KIND=r8)   , INTENT(OUT  ) :: csoi    (ibMax,nsoilay,jbMax)! specific heat of soil, no pore spaces (J kg-1 deg-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: poros   (ibMax,nsoilay,jbMax)! porosity (mass of h2o per unit vol at sat / rhow)
    REAL(KIND=r8)   , INTENT(OUT  ) :: sfield  (ibMax,nsoilay,jbMax)! field capacity soil moisture value (fraction of pore space)
    REAL(KIND=r8)   , INTENT(OUT  ) :: swilt   (ibMax,nsoilay,jbMax)! wilting soil moisture value (fraction of pore space)
    REAL(KIND=r8)   , INTENT(OUT  ) :: bex     (ibMax,nsoilay,jbMax)! exponent "b" in soil water potential
    INTEGER, INTENT(OUT  ) :: ibex    (ibMax,nsoilay,jbMax)! nint(bex), used for cpu speed
    REAL(KIND=r8)   , INTENT(OUT  ) :: suction (ibMax,nsoilay,jbMax)! saturated matric potential (m-h2o)
    REAL(KIND=r8)   , INTENT(OUT  ) :: hydraul (ibMax,nsoilay,jbMax)! saturated hydraulic conductivity (m/s)

    INTEGER, INTENT(IN   ) :: ndat     ! number of soil types 
    ! excludes organics for now.
    REAL(KIND=r8)   , INTENT(IN   ) :: texdat    (3,ndat)! sand/silt/clay fractions
    REAL(KIND=r8)   , INTENT(IN   ) :: porosdat  (ndat)  ! porosity volume fraction
    REAL(KIND=r8)   , INTENT(IN   ) :: sfielddat (ndat)  ! field capacity volume fraction
    REAL(KIND=r8)   , INTENT(IN   ) :: swiltdat  (ndat)  ! wilting point volume fraction
    REAL(KIND=r8)   , INTENT(IN   ) :: bexdat    (ndat)  ! Campbell moisture-release b exponent
    REAL(KIND=r8)   , INTENT(IN   ) :: suctiondat(ndat)  ! Air entry potential (m-H20)
    REAL(KIND=r8)   , INTENT(IN   ) :: hydrauldat(ndat)  ! saturated hydraulic conductivity (m s-1)
    !
    ! Arguments (input)     
    !
    !
    ! local variables
    !
    INTEGER :: i
    INTEGER :: k
    INTEGER :: l
    INTEGER :: j
    INTEGER :: nLndPts
    !     INTEGER :: ndat,    ! number of textural classes
    INTEGER :: msand    ! % of sand in grid point
    INTEGER :: mclay    ! % of clay in grid point
    INTEGER :: lmin     ! closest textural class from texture of 
    ! grid point
    !      INTEGER :: textcls  ! textural class assignment (1..11)
    !
    REAL(KIND=r8)    :: fsand      ! fraction of sand in grid point
    REAL(KIND=r8)    :: fsilt             ! fraction of silt in grid point
    REAL(KIND=r8)    :: fclay             ! fraction of clay in grid point
    REAL(KIND=r8)    :: forganic   ! fraction of organic matter in soil
    !
    ! M. El Maayar modified this....
    !     >     dmin               ! small number

    !      parameter (ndat=11)
    !
    !     REAL(KIND=r8)    :: texdat(3,ndat)    ! % of sand, silt, clay in each textural class
    !     REAL(KIND=r8)    :: porosdat(ndat)    ! porosity in fraction of soil depth
    !     REAL(KIND=r8)    :: sfielddat(ndat)   ! field capacity in fraction of soil depth
    !     REAL(KIND=r8)    :: swiltdat(ndat)    ! wilting point in fraction of soil depth
    !     REAL(KIND=r8)    :: bexdat(ndat)      ! 'b' exponent in moisture release equation
    !     REAL(KIND=r8)    :: suctiondat(ndat)  ! saturated (air entry) potential (m-h2o)
    !     REAL(KIND=r8)    :: hydrauldat(ndat)  ! saturated hydraulic conductivity (m s-1)
    REAL(KIND=r8)    :: xdat(ndat)        ! % of sand in each textural class
    REAL(KIND=r8)    :: ydat(ndat)        ! % of silt in each textural class
    REAL(KIND=r8)    :: zdat(ndat)        ! % of clay in each textural class
    !
    !
    ! Rawls et al. (1992) soil properties data
    !
    !      ------------------
    !       sand  silt  clay
    !      ------------------
    !
    !      data texdat /
    !     >  0.92, 0.05, 0.03,  ! sand
    !     >  0.81, 0.12, 0.07,  ! loamy sand
    !     >  0.65, 0.25, 0.10,  ! sandy loam
    !     >  0.42, 0.40, 0.18,  ! loam
    !     >  0.20, 0.65, 0.15,  ! silt loam
    !     >  0.60, 0.13, 0.27,  ! sandy clay loam
    !     >  0.32, 0.34, 0.34,  ! clay loam
    !     >  0.09, 0.58, 0.33,  ! silty clay loam
    !     >  0.53, 0.07, 0.40,  ! sandy clay
    !     >  0.10, 0.45, 0.45,  ! silty clay
    !     >  0.20, 0.20, 0.60   ! clay
    !     >  /
    !
    ! porosity (fraction)
    !
    !      data porosdat /
    !     >  0.437,             ! sand
    !     >  0.437,             ! loamy sand
    !     >  0.453,             ! sandy loam
    !     >  0.463,             ! loam
    !     >  0.501,             ! silt loam
    !     >  0.398,             ! sandy clay loam
    !     >  0.464,             ! clay loam
    !     >  0.471,             ! silty clay loam
    !     >  0.430,             ! sandy clay
    !     >  0.479,             ! silty clay
    !     >  0.475              ! clay
    !     >  /
    !
    ! field capacity (fraction)
    !
    !      data sfielddat /
    !     >  0.091,             ! sand
    !     >  0.125,             ! loamy sand
    !     >  0.207,             ! sandy loam
    !     >  0.270,             ! loam
    !     >  0.330,             ! silt loam
    !     >  0.255,             ! sandy clay loam
    !     >  0.318,             ! clay loam
    !     >  0.366,             ! silty clay loam
    !     >  0.339,             ! sandy clay
    !     >  0.387,             ! silty clay
    !     >  0.396              ! clay
    !     >  /
    !
    ! wilting point (fraction)
    !
    !      data swiltdat /
    !     >  0.033,             ! sand
    !     >  0.055,             ! loamy sand
    !     >  0.095,             ! sandy loam
    !     >  0.117,             ! loam
    !     >  0.133,             ! silt loam
    !     >  0.148,             ! sandy clay loam
    !     >  0.197,             ! clay loam
    !     >  0.208,             ! silty clay loam
    !     >  0.239,             ! sandy clay
    !     >  0.250,             ! silty clay
    !     >  0.272              ! clay
    !     >  /
    !
    ! "b" exponent for the Campbell moisture-release equation
    !
    !      data bexdat /
    !     >  1.7,               ! sand
    !     >  2.1,               ! loamy sand
    !     >  3.1,               ! sandy loam
    !     >  4.5,               ! loam
    !     >  4.7,               ! silt loam
    !     >  4.0,               ! sandy clay loam
    !     >  5.2,               ! clay loam
    !     >  6.6,               ! silty clay loam
    !     >  6.0,               ! sandy clay
    !     >  7.9,               ! silty clay
    !     >  7.6                ! clay
    !     >  /
    !
    ! saturated (air entry) potential (m-h2o)
    !
    !      data suctiondat /
    !     >  0.070,             ! sand
    !     >  0.090,             ! loamy sand
    !     >  0.150,             ! sandy loam
    !     >  0.110,             ! loam
    !     >  0.210,             ! silt loam
    !     >  0.280,             ! sandy clay loam
    !     >  0.260,             ! clay loam
    !     >  0.330,             ! silty clay loam
    !     >  0.290,             ! sandy clay
    !     >  0.340,             ! silty clay
    !     >  0.370              ! clay
    !     >  /
    !
    ! saturated hydraulic conductivity (m s-1)
    !
    !      data hydrauldat /
    !     >  5.8330e-05,        ! sand
    !     >  1.6972e-05,        ! loamy sand
    !     >  7.1944e-06,        ! sandy loam
    !     >  3.6667e-06,        ! loam
    !     >  1.8889e-06,        ! silt loam
    !     >  1.1944e-06,        ! sandy clay loam
    !     >  6.3889e-07,        ! clay loam
    !     >  4.1667e-07,        ! silty clay loam
    !     >  3.3333e-07,        ! sandy clay
    !     >  2.5000e-07,        ! silty clay
    !     >  1.6667e-07         ! clay
    !     >  /
    !
    ! set sand/silt/clay vectors (xdat,ydat,zdat) for 11 data points
    !
    i=nband
    tsno=tsno
    
    DO l = 1, ndat
       xdat(l) = texdat(1,l)
       ydat(l) = texdat(2,l)
       zdat(l) = texdat(3,l)
    END DO
    !
    ! initialization and normalization constant for puddle model (kg m-2)
    !
    !      wpudmax = 4.5
    !
    IF (irestart == 0) THEN
       DO j = 1,jbMax 
          DO i = 1, nlpoints(j) 
             !CALL const (wpud,  npoi, 0.0)
             wpud (i,j)=0.0_r8
             !CALL const (wipud, npoi, 0.0)
             wipud(i,j)=0.0_r8
          END DO
       END DO
    END IF
    !
    ! set prescribed soil layer thicknesses
    !
    !      hsoi(1)  = 0.10
    !      hsoi(2)  = 0.15
    !      hsoi(3)  = 0.25
    !      hsoi(4)  = 0.50
    !      hsoi(5)  = 1.00
    !      hsoi(6)  = 2.00
    !
    ! set physical parameters of soil
    !
    !CALL const (z0soi, npoi, 0.005)
    DO j = 1,jbMax 
       DO i = 1, nlpoints(j) 
          z0soi(i,j) = 0.005_r8
       END DO
    END DO
    !
    ! initialize soil water and soil temperature fields
    !
    IF (irestart == 0) THEN
       !
       DO k=1,nsoilay
          DO j = 1,jbMax
             nLndPts=0
             DO i = 1, ibMax
                IF (iMaskIBIS(i,j) >= 1_i8) THEN
                   nLndPts=nLndPts+1
                   wsoi (nLndPts,k,j) = wsib(i,j)
                   wsoim(nLndPts,k,j) = wsib(i,j)
                   IF (iMaskIBIS(i,j) >= 15_i8) wisoi(nLndPts,k,j) = 1.00_r8
                   tsoi (nLndPts,k,j) = tg3 (i,j)
                   tsoim(nLndPts,k,j) = tg3 (i,j)
                END IF
             END DO
          END DO
       END DO

       DO j = 1,jbMax 
          nLndPts=0
          DO i = 1, ibMax
             IF (iMaskIBIS(i,j) >= 1_i8) THEN
                nLndPts=nLndPts+1
                tg  (nLndPts,j) = tg3 (i,j)
                tgm (nLndPts,j) = tg3 (i,j)
                !CALL const (ti, npoi, 273.13)
                ti (nLndPts,j) = 273.13_r8
             END IF
          END DO
       END DO
       !     ELSE
       DO j=1,jbMax
          DO i = 1, nlpoints(j) 
             ! tg (i,j) = tsoi (i,1,j)
             ! tgm(i,j) = tg (i,j)! tsoim(i,1,j)
             ! ti (i,j) = tsno (i,1,j)
          END DO
       END DO
    END IF
    !
    ! set soil surface parameters for the global domain
    !
    DO j = 1,jbMax
       DO i = 1, nlpoints(j) 
          !
          ! Convert input sand and clay percents to fractions
          !
          msand = NINT(sand(i,1,j))
          mclay = NINT(clay(i,1,j)) 
          !
          fsand = 0.01_r8 * msand
          fclay = 0.01_r8 * mclay
          fsilt = MAX(0.01_r8 * (100 - msand - mclay),0.0007_r8)
          ! M. El Maayar modified this.
          !        forganic = 1. - fsand - fclay - fsilt
          !
          ! soil surface albedo:
          !
          ! from bats table 3.ii assuming albedo depends on texture
          !
          albsav(i,j) = fsand * 0.120_r8 +  &
               fsilt * 0.085_r8 +  &
               fclay * 0.050_r8
          !
          !
          ! M. El Maayar modified this.
          !      if (nint(forganic).eq.1) then
          !        albsan(i) = 1.0
          !      else
          albsan(i,j) = 2.0_r8 * albsav(i,j)
          !      endif
          !
       END DO
    END DO
    !
    ! create soil properties look-up table
    !
    ! set soil parameters at each layer for the global domain
    ! soita.nc file is for layers only to 4 m; currently this
    ! is for a total of six layers. If there are any remaining  
    ! layers below that, set texture to be equal to that of the
    ! last layer (layer 6)
    ! analysis of the current WISE-IGBP soil textural dataset
    ! reveals very little information below 4 m.
    !
    DO k = 1, nsoilay 
       !DO i = 1, npoi       
       DO j = 1,jbMax
          DO i = 1, nlpoints(j) 

             !
             ! Convert input sand and clay percents to fractions
             !
             IF (k.LE.6) THEN
                msand = NINT(sand(i,k,j))
                mclay = NINT(clay(i,k,j)) 
             ELSE
                msand = NINT(sand(i,6,j)) 
                mclay = NINT(clay(i,6,j)) 
             END IF
             !
             ! M. El Maayar modified this
             !       if ((msand.ge.99).AND.(mclay.ge.99)) then
             !          fsand = 0.
             !          fclay = 0.
             !          fsilt = 0.
             !          forganic = 1.
             !       else
             fsand = 0.01_r8 * msand
             fclay = 0.01_r8 * mclay
             fsilt = MAX(0.01_r8 * (100 - msand - mclay),0.0007_r8)
             !
             ! for now, we assume that all soils have a 1% organic content -- 
             ! this is just a place holder until we couple the soil carbon
             ! dynamics to the soil physical properties
             !
             forganic = 0.010_r8
             !       endif
             !
             ! density of soil material (without pores, not bulk) (kg m-3)
             ! from Campbell and Norman, 1998
             !
             rhosoi(i,k,j) = 2650.0_r8 * (1.0_r8 - forganic) + &
                  1300.0_r8 * forganic 
             !
             ! specific heat of soil material (j kg-1 k-1):
             ! from Campbell and Norman, 1998
             !
             csoi(i,k,j) =  870.0_r8 * (1.0_r8 - forganic) +  &
                  1920.0_r8 * forganic 
             !
             ! cjk
             ! match textural fractions with soil textural class 
             ! calls two functions to match sand and clay fractions
             ! with proper soil textural class based on the usda
             ! classification system
             !
             lmin = textcls (msand,mclay)
             !
             ! porosity (fraction):
             ! 
             poros(i,k,j) = porosdat(lmin)
             !
             ! field capacity (defined relative to the porosity):
             !
             sfield(i,k,j) = 1.0_r8 / poros(i,k,j) * sfielddat(lmin)
             !
             ! wilting point (defined relative to the porosity):
             !
             swilt(i,k,j)  = 1.0_r8 / poros(i,k,j) * swiltdat(lmin)
             !
             ! "b" exponent for the Campbell moisture-release equation:
             !
             bex(i,k,j) = bexdat(lmin)
             !
             ! nearest integer of "b" exponent (for computational efficiency):
             !
             ibex(i,k,j) = NINT(bex(i,k,j))
             !
             ! saturated matric (air entry) potential (m-h2o):
             !
             suction(i,k,j) = suctiondat(lmin)
             !
             ! saturated hydraulic conductivity (m s-1):
             !
             hydraul(i,k,j) = hydrauldat(lmin)
             !
          END DO
       END DO     !DO 310 i = 1, npoi 
    END DO         !DO 300 k = 1, nsoilay 
    !
    ! return to main program
    !
    RETURN
  END SUBROUTINE inisoil

  !-------------------------------------------------------------------------
  INTEGER FUNCTION textcls (msand,mclay)
    !
    ! adapted for ibis by cjk 01/11/01
    !-------------------------------------------------------------------------
    ! |
    ! |                         T R I A N G L E
    ! | Main program that calls WHAT_TEXTURE, a function that classifies soil
    ! | in the USDA textural triangle using sand and clay %
    ! +-----------------------------------------------------------------------
    ! | Created by: aris gerakis, apr. 98 with help from brian baer
    ! | Modified by: aris gerakis, july 99: now all borderline cases are valid
    ! | Modified by: aris gerakis, 30 nov 99: moved polygon initialization to
    ! |              main program
    ! +-----------------------------------------------------------------------
    ! | COMMENTS
    ! | o Supply a data file with two columns, in free format:  1st column sand,
    ! |   2nd column clay %, no header.  The output is a file with the classes.
    ! +-----------------------------------------------------------------------
    ! | You may use, distribute and modify this code provided you maintain
    ! ! this header and give appropriate credit.
    ! +-----------------------------------------------------------------------
    !
    ! code adapted for IBIS by cjk 01-11-01
    !
    !
    INTEGER :: msand
    INTEGER :: mclay
    !
    !      LOGICAL :: inpoly
    !
    REAL(KIND=r8)    :: silty_loam     (1:7,1:2)
    REAL(KIND=r8)    :: sandy          (1:7,1:2)
    REAL(KIND=r8)    :: silty_clay_loam(1:7,1:2) 
    REAL(KIND=r8)    :: loam           (1:7,1:2)
    REAL(KIND=r8)    :: clay_loam      (1:7,1:2)
    REAL(KIND=r8)    :: sandy_loam     (1:7,1:2)
    REAL(KIND=r8)    :: silty_clay     (1:7,1:2)
    REAL(KIND=r8)    :: sandy_clay_loam(1:7,1:2) 
    REAL(KIND=r8)    :: loamy_sand     (1:7,1:2)
    REAL(KIND=r8)    :: clayey         (1:7,1:2)
    !     REAL(KIND=r8)    :: silt           (1:7,1:2) 
    REAL(KIND=r8)    :: sandy_clay     (1:7,1:2)
    !
    ! initalize polygon coordinates:
    ! each textural class reads in the sand coordinates (1,7) first, and
    ! then the corresponding clay coordinates (1,7)

    !     data silty_loam/0, 0, 23, 50, 20, 8, 0, 12, 27, 27, 0, 0, 12, 0/
    !
    ! because we do not have a separate silt category, have to redefine the
    ! polygon boundaries for the silt loam  
    !
    DATA sandy           /85.0_r8, 90.0_r8, 100.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8,  &
         10.0_r8,  0.0_r8,   0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8/
    DATA loamy_sand      /70.0_r8, 85.0_r8,  90.0_r8,85.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8,  & 
         15.0_r8, 10.0_r8,   0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8/
    DATA sandy_loam      /50.0_r8, 43.0_r8,  52.0_r8,52.0_r8,80.0_r8,85.0_r8, 70.0_r8,  &
         0.0_r8,  7.0_r8,   7.0_r8,20.0_r8,20.0_r8,15.0_r8,  0.0_r8/
    DATA loam            /43.0_r8, 23.0_r8,  45.0_r8,52.0_r8,52.0_r8, 0.0_r8,  0.0_r8,    &
         7.0_r8, 27.0_r8,  27.0_r8,20.0_r8, 7.0_r8, 0.0_r8,  0.0_r8/
    DATA silty_loam      / 0.0_r8,  0.0_r8,  23.0_r8,50.0_r8, 0.0_r8, 0.0_r8,  0.0_r8, 0.0_r8,    &
         27.0_r8, 27.0_r8,   0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8/ 
    !     DATA silt            /0, 0, 8, 20, 0, 0, 0, 0, 12, 12, 0, 0, 0, 0/
    DATA sandy_clay_loam /52.0_r8, 45.0_r8, 45.0_r8, 65.0_r8, 80.0_r8, 0.0_r8, 0.0_r8,    & 
         20.0_r8, 27.0_r8, 35.0_r8, 35.0_r8, 20.0_r8, 0.0_r8, 0.0_r8/
    DATA clay_loam       /20.0_r8, 20.0_r8, 45.0_r8, 45.0_r8, 0.0_r8, 0.0_r8, 0.0_r8,     &
         27.0_r8, 40.0_r8, 40.0_r8, 27.0_r8, 0.0_r8, 0.0_r8, 0.0_r8/
    DATA silty_clay_loam /0.0_r8, 0.0_r8, 20.0_r8, 20.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 27.0_r8,   &
         40.0_r8, 40.0_r8, 27.0_r8, 0.0_r8, 0.0_r8, 0.0_r8/
    DATA sandy_clay      /45.0_r8, 45.0_r8, 65.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8,      &
         35.0_r8, 55.0_r8, 35.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8/
    DATA silty_clay      /0.0_r8, 0.0_r8, 20.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 40.0_r8,    &
         60.0_r8, 40.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8/
    DATA clayey          /20.0_r8, 0.0_r8, 0.0_r8, 45.0_r8, 45.0_r8, 0.0_r8, 0.0_r8,      &
         40.0_r8, 60.0_r8, 100.0_r8, 55.0_r8, 40.0_r8, 0.0_r8, 0.0_r8/
    !
    ! polygon coordinates  
    !
    !     sand
    !
    !     >  85, 90, 100, 0, 0, 0, 0,       ! sand
    !     >  70, 85, 90, 85, 0, 0, 0,       ! loamy sand
    !     >  50, 43, 52, 52, 80, 85, 70,    ! sandy loam
    !     >  43, 23, 45, 52, 52, 0, 0,      ! loam
    !     >   0, 0, 23, 50, 0, 0, 0,        ! silt loam (combined with silt)
    !     >  52, 45, 45, 65, 80, 0, 0,      ! sandy clay loam
    !     >  20, 20, 45, 45, 0, 0, 0,       ! clay loam
    !     >   0, 0, 20, 20, 0, 0, 0,        ! silty clay loam
    !     >  45, 45, 65, 0, 0, 0, 0,        ! sandy clay
    !     >   0, 0, 20, 0, 0, 0, 0,         ! silty clay 
    !     >  20, 0, 0, 45, 45, 0, 0         ! clay
    !
    !      clay
    !
    !     > 0, 10, 0, 0, 0, 0, 0,           ! sand
    !     > 0, 15, 10, 0, 0, 0, 0,          ! loamy sand
    !     > 0, 7, 7, 20, 20, 15, 0,         ! sandy loam 
    !     > 7, 27, 27, 20, 7, 0, 0,         ! loam
    !     > 0, 27, 27, 0, 0, 0, 0,          ! silt loam (combined with silt)
    !     > 20, 27, 35, 35, 20, 0, 0,       ! sandy clay loam
    !     > 27, 40, 40, 27, 0, 0, 0,        ! clay loam
    !     > 27, 40, 40, 27, 0, 0, 0,        ! silty clay loam
    !     > 35, 55, 35, 0, 0, 0, 0,         ! sandy clay
    !     > 40, 60, 40, 0, 0, 0, 0,         ! silty clay
    !     > 40, 60, 100, 55, 40, 0, 0       ! clay
    !
    ! +-----------------------------------------------------------------------
    ! | figure out what texture grid cell and layer are part of  
    ! | classify a soil in the triangle based on sand and clay %
    ! +-----------------------------------------------------------------------
    ! | Created by: aris gerakis, apr. 98
    ! | Modified by: aris gerakis, june 99.  Now check all polygons instead of
    ! | stopping when a right solution is found.  This to cover all borderline 
    ! | cases.
    ! +-----------------------------------------------------------------------
    !
    ! find polygon(s) where the point is.  
    !
    textcls = 0 
    !
    IF (msand .GT. 0 .AND. mclay .GT. 0) THEN
       IF (inpoly(sandy, 3, msand, mclay)) THEN
          textcls = 1      ! sand
       END IF
       IF (inpoly(loamy_sand, 4, msand, mclay)) THEN
          textcls = 2      ! loamy sand
       END IF
       IF (inpoly(sandy_loam, 7, msand, mclay)) THEN
          textcls = 3      ! sandy loam
       END IF
       IF (inpoly(loam, 5, msand, mclay)) THEN
          textcls = 4      ! loam
       END IF
       IF (inpoly(silty_loam, 4, msand, mclay)) THEN
          textcls = 5      ! silt loam
       END IF
       IF (inpoly(sandy_clay_loam, 5, msand, mclay)) THEN
          textcls = 6      ! sandy clay loam
       END IF
       IF (inpoly(clay_loam, 4, msand, mclay)) THEN
          textcls = 7      ! clay loam
       END IF
       IF (inpoly(silty_clay_loam, 4, msand, mclay)) THEN
          textcls = 8      ! silty clay loam
       END IF
       IF (inpoly(sandy_clay, 3, msand, mclay)) THEN
          textcls = 9      ! sandy clay
       END IF
       IF (inpoly(silty_clay, 3, msand, mclay)) THEN
          textcls = 10     ! silty clay
       END IF
       IF (inpoly(clayey, 5, msand, mclay)) THEN
          textcls = 11     ! clay
       END IF
    END IF
    !
    IF (textcls .EQ. 0) THEN
       textcls = 5         ! silt loam
       !
       !        write (*, 1000) msand, mclay
       ! 1000   format (/, 1x, 'Texture not found for ', f5.1, ' sand and ', f5.1, ' clay')
    END IF
    !
    RETURN
  END FUNCTION textcls
  !
  !---------------------------------------------------------------------------
  LOGICAL FUNCTION inpoly (poly, npoints, xt, yt)
    !
    ! adapted for ibis by cjk 01/11/01
    !---------------------------------------------------------------------------
    !
    !                            INPOLY
    !   Function to tell if a point is inside a polygon or not.
    !--------------------------------------------------------------------------
    !   Copyright (c) 1995-1996 Galacticomm, Inc.  Freeware source code.
    !
    !   Please feel free to use this source code for any purpose, commercial
    !   or otherwise, as long as you don't restrict anyone else's use of
    !   this source code.  Please give credit where credit is due.
    !
    !   Point-in-polygon algorithm, created especially for World-Wide Web
    !   servers to process image maps with mouse-clickable regions.
    !
    !   Home for this file:  http://www.gcomm.com/develop/inpoly.c
    !
    !                                       6/19/95 - Bob Stein & Craig Yap
    !                                       stein@gcomm.com
    !                                       craig@cse.fau.edu
    !--------------------------------------------------------------------------
    !   Modified by:
    !   Aris Gerakis, apr. 1998: 1.  translated to Fortran
    !                            2.  made it work with REAL(KIND=r8) coordinates
    !                            3.  now resolves the case where point falls
    !                                on polygon border.
    !   Aris Gerakis, nov. 1998: Fixed error caused by hardware arithmetic
    !   Aris Gerakis, july 1999: Now all borderline cases are valid
    !--------------------------------------------------------------------------
    !   Glossary:
    !   function inpoly: true=inside, false=outside (is target point inside
    !                    a 2D polygon?)
    !   poly(*,2):  polygon points, [0]=x, [1]=y
    !   npoints: number of points in polygon
    !   xt: x (horizontal) of target point
    !   yt: y (vertical) of target point
    !--------------------------------------------------------------------------
    !
    ! declare arguments  
    !
    INTEGER :: npoints
    INTEGER :: xt
    INTEGER :: yt 
    !
    REAL(KIND=r8)    :: poly(7, 2)
    !
    ! local variables
    !
    REAL(KIND=r8)    :: xnew
    REAL(KIND=r8)    :: ynew
    REAL(KIND=r8)    :: xold
    REAL(KIND=r8)    :: yold
    REAL(KIND=r8)    :: x1
    REAL(KIND=r8)    :: y1
    REAL(KIND=r8)    :: x2
    REAL(KIND=r8)    :: y2
    !
    INTEGER ::  i
    !
    LOGICAL :: inside
    LOGICAL :: on_border

    inside = .FALSE.
    on_border = .FALSE.
    !
    IF (npoints .LT. 3)  THEN
       inpoly = .FALSE.
       RETURN
    END IF
    !
    xold = poly(npoints,1)
    yold = poly(npoints,2)

    DO i = 1 , npoints
       xnew = poly(i,1)
       ynew = poly(i,2)

       IF (xnew .GT. xold)  THEN
          x1 = xold
          x2 = xnew
          y1 = yold
          y2 = ynew
       ELSE
          x1 = xnew
          x2 = xold
          y1 = ynew
          y2 = yold
       END IF

       ! the outer IF is the 'straddle' test and the 'vertical border' test.
       ! the inner IF is the 'non-vertical border' test and the 'north' test.  

       ! the first statement checks whether a north pointing vector crosses  
       ! (stradles) the straight segment.  There are two possibilities, depe-
       ! nding on whether xnew < xold or xnew > xold.  The '<' is because edge 
       ! must be "open" at left, which is necessary to keep correct count when 
       ! vector 'licks' a vertix of a polygon.  

       IF ((xnew .LT. xt .AND. xt .LE. xold)   &
            .OR. (.NOT. xnew .LT. xt .AND.       &
            .NOT. xt .LE. xold)) THEN
          !
          ! the test point lies on a non-vertical border:
          !
          IF ((yt-y1)*(x2-x1) .EQ. (y2-y1)*(xt-x1)) THEN

             on_border = .TRUE. 
             !
             ! check if segment is north of test point.  If yes, reverse the 
             ! value of INSIDE.  The +0.001 was necessary to avoid errors due   
             ! arithmetic (e.g., when clay = 98.87 and sand = 1.13):   
             !
          ELSE IF ((yt-y1)*(x2-x1) .LT. (y2-y1)*(xt-x1) + 0.001) THEN

             inside = .NOT.inside ! cross a segment

          END IF
          !
          ! this is the rare case when test point falls on vertical border or  
          ! left edge of non-vertical border. The left x-coordinate must be  
          ! common.  The slope requirement must be met, but also point must be
          ! between the lower and upper y-coordinate of border segment.  There 
          ! are two possibilities,  depending on whether ynew < yold or ynew > 
          ! yold:
          !
       ELSE IF ((xnew .EQ. xt .OR. xold .EQ. xt)       &
            .AND. (yt-y1)*(x2-x1) .EQ.             &
            (y2-y1)*(xt-x1) .AND. ((ynew .LE. yt   &
            .AND. yt .LE. yold) .OR.               &
            (.NOT. ynew .LT. yt .AND. .NOT. yt .LT. yold))) THEN

          on_border = .TRUE. 

       END IF
       !
       xold = xnew
       yold = ynew
       !
    END DO!  DO i = 1 , npoints  
    !
    ! If test point is not on a border, the function result is the last state 
    ! of INSIDE variable.  Otherwise, INSIDE doesn't matter.  The point is
    ! inside the polygon if it falls on any of its borders:
    !
    IF (.NOT. on_border) THEN
       inpoly = inside
    ELSE
       inpoly = .TRUE.
    END IF
    !
    RETURN
  END FUNCTION inpoly
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE iniveg (isimveg      , &! INTENT(IN   )
       irestart     , &! INTENT(IN   )
       plai_init    , &! INTENT(IN   )
       plaiupper    , &! INTENT(IN   )
       plailower    , &! INTENT(IN   )
       xminlai      , &! INTENT(IN   )
       sapfrac_init , &! INTENT(IN   )
       chiflz       , &! INTENT(IN   )
       chifuz       , &! INTENT(IN   )
       beta1        , &! INTENT(IN   )
       beta2              , &! INTENT(IN   )
       agddu        , &! INTENT(OUT  )
       agddl              , &! INTENT(OUT  )
       falll        , &! INTENT(OUT  )
       fallr        , &! INTENT(OUT  )
       fallw        , &! INTENT(OUT  )
       exist              , &! INTENT(IN   )
       vegtype0     , &! INTENT(OUT  )
       plai         , &! INTENT(OUT  )
       tw           , &! INTENT(IN   )
       sapfrac      , &! INTENT(OUT  )
       cbiol        , &! INTENT(INOUT)
       specla              , &! INTENT(IN   )
       cbior        , &! INTENT(INOUT)
       cbiow              , &! INTENT(INOUT)
       biomass      , &! INTENT(OUT  )
       totlaiu      , &! INTENT(OUT  )
       totlail      , &! INTENT(OUT  )
       totbiou      , &! INTENT(OUT  )
       totbiol      , &! INTENT(OUT  )
       sai          , &! INTENT(OUT  )
       fu           , &! INTENT(OUT  )
       woodnorm     , &! INTENT(IN   )
       fl           , &! INTENT(OUT  )
       lai          , &! INTENT(OUT  )
       zbot         , &! INTENT(OUT  )
       ztop              , &! INTENT(OUT  )
       oriev        , &! INTENT(OUT  )
       orieh        , &! INTENT(OUT  )
       froot        , &! INTENT(OUT  )
       a10td              , &! INTENT(OUT  )
       a10ancub     , &! INTENT(OUT  )
       a10ancuc     , &! INTENT(OUT  )
       a10ancls     , &! INTENT(OUT  )
       a10ancl4     , &! INTENT(OUT  )
       a10ancl3     , &! INTENT(OUT  )
       a10scalparamu, &! INTENT(OUT  )
       a10scalparaml, &! INTENT(OUT  )
       a10daylightu , &! INTENT(OUT  )
       a10daylightl , &! INTENT(OUT  )
       stresstu     , &! INTENT(OUT  )
       stresstl     , &! INTENT(OUT  )
       hsoi              , &! INTENT(IN   )
       xinveg       , &! INTENT(IN   )
       nsoilay      , &! INTENT(IN   )
       gdd0this     , &! INTENT(OUT  )
       gdd5this     , &! INTENT(OUT  )
       tcthis       , &! INTENT(OUT  )
       twthis       , &! INTENT(OUT  )
       npft           )! INTENT(IN   )
    ! ---------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN   ) :: nsoilay                  ! number of soil layers
    INTEGER, INTENT(IN   ) :: npft                     ! number of plant functional types
    REAL(KIND=r8)   , INTENT(IN   ) :: xinveg        (ibMax,jbMax) ! fixed vegetation map
    REAL(KIND=r8)   , INTENT(OUT  ) :: stresstu      (ibMax,jbMax) ! sum of stressu over all 6 soil layers (dimensionless)
    REAL(KIND=r8)   , INTENT(OUT  ) :: stresstl      (ibMax,jbMax) ! sum of stressl over all 6 soil layers (dimensionless)
    REAL(KIND=r8)   , INTENT(IN   ) :: hsoi              (nsoilay+1)! soil layer thickness (m)
    REAL(KIND=r8)   , INTENT(OUT  ) :: a10td         (ibMax,jbMax)! 10-day average daily air temperature (K)
    REAL(KIND=r8)   , INTENT(OUT  ) :: a10ancub      (ibMax,jbMax)! 10-day average canopy photosynthesis rate -
    ! broadleaf (mol_co2 m-2 s-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: a10ancuc      (ibMax,jbMax)! 10-day average canopy photosynthesis rate - 
    ! conifer (mol_co2 m-2 s-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: a10ancls      (ibMax,jbMax)! 10-day average canopy photosynthesis rate - 
    ! shrubs (mol_co2 m-2 s-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: a10ancl4      (ibMax,jbMax)! 10-day average canopy photosynthesis rate - 
    ! c4 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: a10ancl3      (ibMax,jbMax)! 10-day average canopy photosynthesis rate -
    ! c3 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: a10scalparamu (ibMax,jbMax)! 10-day average day-time scaling parameter -
    ! upper canopy (dimensionless)
    REAL(KIND=r8)   , INTENT(OUT  ) :: a10scalparaml (ibMax,jbMax)! 10-day average day-time scaling parameter -
    ! lower canopy (dimensionless)
    REAL(KIND=r8)   , INTENT(OUT  ) :: a10daylightu  (ibMax,jbMax)! 10-day average day-time PAR - 
    ! upper canopy (micro-Ein m-2 s-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: a10daylightl  (ibMax,jbMax)! 10-day average day-time PAR - 
    ! lower canopy (micro-Ein m-2 s-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: agddu         (ibMax,jbMax)! annual accumulated growing degree days for bud 
    ! burst, upper canopy (day-degrees)
    REAL(KIND=r8)   , INTENT(OUT  ) :: agddl         (ibMax,jbMax)! annual accumulated growing degree days for bud 
    ! burst, lower canopy (day-degrees)
    REAL(KIND=r8)   , INTENT(OUT  ) :: falll         (ibMax,jbMax) ! annual leaf litter fall (kg_C m-2/year)
    REAL(KIND=r8)   , INTENT(OUT  ) :: fallr         (ibMax,jbMax) ! annual root litter input (kg_C m-2/year)
    REAL(KIND=r8)   , INTENT(OUT  ) :: fallw         (ibMax,jbMax) ! annual wood litter fall (kg_C m-2/year)
    REAL(KIND=r8)   , INTENT(IN   ) :: exist         (ibMax,npft,jbMax)! probability of existence of each plant
    ! functional type in a gridcell
    REAL(KIND=r8)   , INTENT(OUT  ) :: vegtype0      (ibMax,jbMax)! annual vegetation type - ibis classification
    REAL(KIND=r8)   , INTENT(OUT  ) :: plai          (ibMax,npft,jbMax)! total leaf area index of each plant functional type
    REAL(KIND=r8)   , INTENT(IN   ) :: tw            (ibMax,jbMax)    ! warmest monthly temperature (C)
    REAL(KIND=r8)   , INTENT(OUT  ) :: sapfrac       (ibMax,jbMax)    ! fraction of woody biomass that is in sapwood
    REAL(KIND=r8)   , INTENT(INOUT) :: cbiol         (ibMax,npft,jbMax)! carbon in leaf biomass pool (kg_C m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: specla        (npft)     ! specific leaf area (m**2/kg) 
    REAL(KIND=r8)   , INTENT(INOUT) :: cbior         (ibMax,npft,jbMax) ! carbon in fine root biomass pool (kg_C m-2)
    REAL(KIND=r8)   , INTENT(INOUT) :: cbiow         (ibMax,npft,jbMax)! carbon in woody biomass pool (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: biomass       (ibMax,npft,jbMax)! total biomass of each plant functional type  (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: totlaiu       (ibMax,jbMax)! total leaf area index for the upper canopy
    REAL(KIND=r8)   , INTENT(OUT  ) :: totlail       (ibMax,jbMax)! total leaf area index for the lower canopy
    REAL(KIND=r8)   , INTENT(OUT  ) :: totbiou       (ibMax,jbMax)! total biomass in the upper canopy (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: totbiol       (ibMax,jbMax)! total biomass in the lower canopy (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: sai           (ibMax,2,jbMax)   ! current single-sided stem area index
    REAL(KIND=r8)   , INTENT(OUT  ) :: fu            (ibMax,jbMax)! fraction of overall area covered by upper canopy
    REAL(KIND=r8)   , INTENT(OUT  ) :: gdd0this      (ibMax,jbMax)! INTENT(OUT  )
    REAL(KIND=r8)   , INTENT(OUT  ) :: gdd5this      (ibMax,jbMax)! INTENT(OUT  )
    REAL(KIND=r8)   , INTENT(OUT  ) :: tcthis        (ibMax,jbMax)! INTENT(OUT  )
    REAL(KIND=r8)   , INTENT(OUT  ) :: twthis        (ibMax,jbMax)! INTENT(OUT  )

    REAL(KIND=r8)   , INTENT(IN   ) :: woodnorm                 ! value of woody biomass for upper canopy closure 
    ! (ie when wood = woodnorm fu = 1.0) (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: fl            (ibMax,jbMax)! fraction of snow-free area covered by lower  canopy
    REAL(KIND=r8)   , INTENT(OUT  ) :: lai     (ibMax,2,jbMax)! canopy single-sided leaf area index (area leaf/area veg)
    REAL(KIND=r8)   , INTENT(OUT  ) :: zbot    (ibMax,2,jbMax) ! height of lowest branches above ground (m)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ztop    (ibMax,2,jbMax)! height of plant top above ground (m)
    REAL(KIND=r8)   , INTENT(OUT  ) :: oriev   (2)           ! fraction of leaf/stems with vertical
    REAL(KIND=r8)   , INTENT(OUT  ) :: orieh   (2)           ! fraction of leaf/stems with horizontal orientation
    REAL(KIND=r8)   , INTENT(OUT  ) :: froot   (nsoilay,2)   ! fraction of root in soil layer 
    REAL(KIND=r8)   , INTENT(IN   ) :: plai_init   (4,15)    ! initial total LAI for each vegtype (used in iniveg)
    REAL(KIND=r8)   , INTENT(IN   ) :: plaiupper             ! Potental LAI of upper canopy (uniform initial vegetation)
    REAL(KIND=r8)   , INTENT(IN   ) :: plailower             ! Potental LAI of lower canopy (uniform initial vegetation)
    REAL(KIND=r8)   , INTENT(IN   ) :: xminlai               ! Minimum LAI for each existing PFT
    REAL(KIND=r8)   , INTENT(IN   ) :: sapfrac_init          ! Initial value of sapwood fraction used for all woody PFTs
    REAL(KIND=r8)   , INTENT(IN   ) :: chiflz                ! lower canopy leaf orientation factor
    REAL(KIND=r8)   , INTENT(IN   ) :: chifuz                ! upper canopy leaf orientation factor 
    REAL(KIND=r8)   , INTENT(IN   ) :: beta1                 ! parameter for Jackson rooting profile, lower canopy
    REAL(KIND=r8)   , INTENT(IN   ) :: beta2                 ! parameter for Jackson rooting profile, upper canopy
    !
    ! Arguments (input)
    !
    INTEGER, INTENT(IN   ) :: irestart     ! 0: not a restart run 1: restart run
    INTEGER, INTENT(IN   ) :: isimveg      ! 0:static veg 1:dyn veg 2:dyn veg-cold start
    !
    ! local variables
    !
    INTEGER :: ideci        ! # deciduous plant functional types (pft)
    INTEGER :: ievgr        ! # evergreen pft 
    INTEGER :: ishrub       ! # shrub pft 
    INTEGER :: igrass       ! # herbaceous pft 
    INTEGER :: ilower       ! possible # pft for lower canopy
    INTEGER :: iupper       ! possible # pft for upper canopy
    INTEGER :: inveg        ! vegetation type
    INTEGER :: i            ! loop indices
    INTEGER :: j            ! loop indices
    INTEGER :: k            ! loop indices
    !
    REAL(KIND=r8)    :: plaievgr     ! potential lai of evergreen trees
    REAL(KIND=r8)    :: plaideci     ! potential lai of deciduous trees
    REAL(KIND=r8)    :: plaishrub    ! potential lai of shrubs
    REAL(KIND=r8)    :: plaigrass    ! potential lai of grasses
    !      REAL(KIND=r8)    :: plaiupper   ! potential lai of upper canopy (uniform initial veg)
    !      REAL(KIND=r8)    :: plailower   ! potential lai of lower canopy (uniform initial veg)
    !      REAL(KIND=r8)    :: xminlai     ! minimum lai for each existing plant type
    REAL(KIND=r8)    ::  wood        ! total wood biomas in grid cell
    !      REAL(KIND=r8)    :: chiflz      ! lower canopy leaf orientation factor 
    !    (-1 vertical, 0 random, 1 horizontal)
    !      REAL(KIND=r8)    :: chifuz      ! uppuer canopy leaf orientation factor 
    !      REAL(KIND=r8)    :: beta1       ! parameter for Jackson rooting profile, lower canopy
    !      REAL(KIND=r8)    :: beta2       ! parameter for Jackson rooting profile, upper canopy
    REAL(KIND=r8)    :: totdepth     ! total soil depth
    REAL(KIND=r8)    :: frootnorm1   ! normalization factor for Jackson rooting profile,low
    REAL(KIND=r8)    :: frootnorm2   ! normalization factor for Jackson rooting profile, up

    !
    REAL(KIND=r8)    :: depth(nsoilay)   ! soil layer depth (cm)
    !
    ! initialize specific leaf area values
    !
    !      data specla  / 25.0,  ! tropical broadleaf evergreen trees
    !     >               25.0,  ! tropical broadleaf drought-deciduous trees
    !     >               25.0,  ! warm-temperate broadleaf evergreen trees
    !     >               12.5,  ! temperate conifer evergreen trees
    !     >               25.0,  ! temperate broadleaf cold-deciduous trees
    !     >               12.5,  ! boREAL(KIND=r8) conifer evergreen trees
    !     >               25.0,  ! boreal broadleaf cold-deciduous trees  
    !     >               25.0,  ! boreal conifer cold-deciduous trees
    !     >               12.5,  ! evergreen shrubs 
    !     >               25.0,  ! deciduous shrubs 
    !     >               20.0,  ! warm (c4) grasses
    !     >               20.0 / ! cool (c3) grasses
    !
    !      woodnorm = 7.5
    !
    ! set c allocation coefficients for natural vegetation

    !      aleaf(1)  = 0.30
    !      aroot(1)  = 0.20
    !      awood(1)  = 1. - aleaf(1) - aroot(1)
    !
    !      aleaf(2)  = 0.30
    !      aroot(2)  = 0.20
    !      awood(2)  = 1. - aleaf(2) - aroot(2)
    !
    !      aleaf(3)  = 0.30
    !      aroot(3)  = 0.20
    !      awood(3)  = 1. - aleaf(3) - aroot(3)
    !
    !      aleaf(4)  = 0.30
    !      aroot(4)  = 0.40
    !      awood(4)  = 1. - aleaf(4) - aroot(4)
    !
    !      aleaf(5)  = 0.30
    !      aroot(5)  = 0.20
    !      awood(5)  = 1. - aleaf(5) - aroot(5)
    !
    !      aleaf(6)  = 0.30
    !      aroot(6)  = 0.40
    !      awood(6)  = 1. - aleaf(6) - aroot(6)
    !
    !      aleaf(7)  = 0.30
    !      aroot(7)  = 0.20
    !      awood(7)  = 1. - aleaf(7) - aroot(7)
    !
    !      aleaf(8)  = 0.30
    !      aroot(8)  = 0.20
    !      awood(8)  = 1. - aleaf(8) - aroot(8)
    !
    ! allocation coefficients for shrubs
    !
    !      aleaf(9)  = 0.45
    !      aroot(9)  = 0.40
    !      awood(9)  = 1. - aleaf(9) - aroot(9)
    !
    !      aleaf(10) = 0.45
    !      aroot(10) = 0.35
    !      awood(10) = 1. - aleaf(10) - aroot(10)
    !
    ! allocation coefficients for grasses
    !
    !      aleaf(11) = 0.45
    !      aroot(11) = 0.55
    !      awood(11) = 0.00
    !
    !      aleaf(12) = 0.45
    !      aroot(12) = 0.55
    !      awood(12) = 0.00
    !
    !DO i = 1, npoi
    DO j=1,jbMax
       DO i=1,nlpoints(j) 
          !
          ! initialize a few climatic variables needed for vegetation
          !
          IF (irestart ==  0) THEN
             agddu(i,j) = 1000.0_r8
             agddl(i,j) = 1000.0_r8
             !
             ! initialize the moisture stress factors
             !
             stresstu(i,j) = 1.0_r8
             stresstl(i,j) = 1.0_r8
             !
             DO k = 1, nsoilay
                stressl(i,k,j) = 1.0_r8
                stressu(i,k,j) = 1.0_r8
             END DO

             !
             ! initialize running-mean air temperature
             !
             !
             td   (i,j)    = 278.16_r8
             a10td(i,j)    = 273.16_r8
             !
             ! initialize running-mean values of canopy photosynthesis rates
             !
             a10ancub(i,j) = 10.0e-06_r8
             a10ancuc(i,j) = 10.0e-06_r8
             a10ancls(i,j) = 10.0e-06_r8
             a10ancl4(i,j) = 10.0e-06_r8
             a10ancl3(i,j) = 10.0e-06_r8
             ! 
             ! initialize running-mean values of the scaling parameter
             !
             a10scalparamu(i,j) = 0.5_r8 * 5.0_r8
             a10scalparaml(i,j) = 0.5_r8 * 5.0_r8
             a10daylightu(i,j) = 5.0_r8
             a10daylightl(i,j) = 5.0_r8
             !
             ! initialize litter fall
             !
             falll(i,j) = 0.0_r8
             fallr(i,j) = 0.0_r8
             fallw(i,j) = 0.0_r8
             !
             !     initialize this year's growing degree days and temperature of the
             !     warmest and coldest month for existence of pfts (in off-line IBIS,
             !     done in weather.f)
             ! 
             gdd0this(i,j) = 0.0_r8
             gdd5this(i,j) = 0.0_r8
             tcthis  (i,j) = 100.0_r8
             twthis  (i,j) = - 100.0_r8        

             !
             !END IF
             !
             ! reset counters
             !
             ievgr  = 0
             ideci  = 0
             ishrub = 0
             igrass = 0
             ilower = 0
             iupper = 0
             !
             ! determine number of evergreen plant functional types
             !
             IF (NINT(exist(i,1,j)).EQ.1) ievgr = ievgr + 1
             IF (NINT(exist(i,3,j)).EQ.1) ievgr = ievgr + 1
             IF (NINT(exist(i,4,j)).EQ.1) ievgr = ievgr + 1
             IF (NINT(exist(i,6,j)).EQ.1) ievgr = ievgr + 1
             !
             ! determine number of deciduous plant functional types
             !
             IF (NINT(exist(i,2,j)).EQ.1) ideci = ideci + 1
             IF (NINT(exist(i,5,j)).EQ.1) ideci = ideci + 1
             IF (NINT(exist(i,7,j)).EQ.1) ideci = ideci + 1
             IF (NINT(exist(i,8,j)).EQ.1) ideci = ideci + 1
             !
             ! make sure counter is at least 1 (to avoid division by zero)
             !
             ievgr = MAX (1, ievgr)
             ideci = MAX (1, ideci)
             !
             ! determine number of shrub functional types
             !
             IF (NINT(exist(i,9,j)).EQ.1)  ishrub = ishrub + 1
             IF (NINT(exist(i,10,j)).EQ.1) ishrub = ishrub + 1
             !
             ! determine number of herbaceous plant functional types
             !
             IF (NINT(exist(i,11,j)).EQ.1) igrass = igrass + 1
             IF (NINT(exist(i,12,j)).EQ.1) igrass = igrass + 1
             !
             ! make sure counter is at least 1 (to avoid division by zero)
             !
             ishrub = MAX (1, ishrub)
             igrass = MAX (1, igrass)
             !
             ! total number of possible pfts for each canopy
             !
             iupper = ievgr  + ideci
             ilower = ishrub + igrass 
             !
             ! make sure counter is at least 1 (to avoid division by zero)
             !
             iupper = MAX (1, iupper)
             ilower = MAX (1, ilower)
             !
             ! cold start of the vegetation
             !
             !IF (irestart == 0) THEN
             !
             !
             ! ************************************************************************
             ! case (0) assign vegetation characteristics for static vegetation
             ! ************************************************************************
             !
             ! and
             !
             ! ************************************************************************
             ! case (1) assign vegetation characteristics for dynamic vegtation
             !          that is initialized with fixed vegetation map
             ! ************************************************************************
             !
             IF (isimveg == 0 .OR. isimveg == 1) THEN
                !
                ! translate vegetation type (real) to nearest integer
                !
                inveg = NINT (xinveg(i,j))
                !
                ! for initialization purposes, set the predicted vegetation type
                ! to the initial vegetation type
                !
                vegtype0(i,j) = xinveg(i,j)
                !
                ! ---------------------------------------------------
                !  1: tropical evergreen forest / woodland
                !  2: tropical deciduous forest / woodland
                !  3: temperate evergreen broadleaf forest / woodland
                !  4: temperate evergreen conifer forest / woodland
                !  5: temperate deciduous forest / woodland
                !  6: boreal evergreen forest / woodland
                !  7: boreal deciduous forest / woodland
                !  8: mixed forest / woodland
                !  9: savanna
                ! 10: grassland / steppe
                ! 11: dense shrubland
                ! 12: open shrubland
                ! 13: tundra
                ! 14: desert
                ! 15: polar desert / rock / ice
                ! ---------------------------------------------------
                !
                ! these classes consist of some combination of 
                ! plant functional types:
                !
                ! ---------------------------------------------------
                !  1: tropical broadleaf evergreen trees
                !  2: tropical broadleaf drought-deciduous trees
                !  3: warm-temperate broadleaf evergreen trees
                !  4: temperate conifer evergreen trees
                !  5: temperate broadleaf cold-deciduous trees
                !  6: boreal conifer evergreen trees
                !  7: boreal broadleaf cold-deciduous trees
                !  8: boreal conifer cold-deciduous trees
                !  9: evergreen shrubs
                ! 10: cold-deciduous shrubs
                ! 11: warm (c4) grasses
                ! 12: cool (c3) grasses
                ! ---------------------------------------------------

                !*** DTP 2001/05/25. The following code replaces the 450+
                !    lines of stuff that follows it (hence the temporary goto
                !    statement). Note that values of plai_init are read in as
                !    parameters from params.veg. Note also that the declarations
                !    of the four local variables plaievgr, plaideci, plaishrub 
                !    and plaigrass can all be dropped.
                !
                plai(i,1,j)  = exist(i,1,j)  / REAL(ievgr,kind=r8)  * plai_init(1,inveg)
                plai(i,2,j)  = exist(i,2,j)  / REAL(ideci,kind=r8)  * plai_init(2,inveg)
                plai(i,3,j)  = exist(i,3,j)  / REAL(ievgr,kind=r8)  * plai_init(1,inveg)
                plai(i,4,j)  = exist(i,4,j)  / REAL(ievgr,kind=r8)  * plai_init(1,inveg)
                plai(i,5,j)  = exist(i,5,j)  / REAL(ideci,kind=r8)  * plai_init(2,inveg)
                plai(i,6,j)  = exist(i,6,j)  / REAL(ievgr,kind=r8)  * plai_init(1,inveg)
                plai(i,7,j)  = exist(i,7,j)  / REAL(ideci,kind=r8)  * plai_init(2,inveg)
                plai(i,8,j)  = exist(i,8,j)  / REAL(ideci,kind=r8)  * plai_init(2,inveg)
                plai(i,9,j)  = exist(i,9,j)  / REAL(ishrub,kind=r8) * plai_init(3,inveg)
                plai(i,10,j) = exist(i,10,j) / REAL(ishrub,kind=r8) * plai_init(3,inveg)
                plai(i,11,j) = exist(i,11,j) / REAL(igrass,kind=r8) * plai_init(4,inveg)
                plai(i,12,j) = exist(i,12,j) / REAL(igrass,kind=r8) * plai_init(4,inveg)
                IF ((inveg.EQ.9).OR.(inveg.EQ.10)) THEN
                   IF (tw(i,j).GT.22.0_r8) THEN
                      plai(i,11,j) = exist(i,11,j) * 0.80_r8 * plai_init(4,inveg)
                      plai(i,12,j) = exist(i,12,j) * 0.20_r8 * plai_init(4,inveg)
                   ELSE
                      plai(i,11,j) = exist(i,11,j) * 0.00_r8 * plai_init(4,inveg)
                      plai(i,12,j) = exist(i,12,j) * 1.00_r8 * plai_init(4,inveg)
                   END IF
                ELSE
                   plai(i,11,j) = exist(i,11,j) / REAL(igrass,kind=r8) * plai_init(4,inveg)
                   plai(i,12,j) = exist(i,12,j) / REAL(igrass,kind=r8) * plai_init(4,inveg)
                END IF

                GOTO 9999


                !
                ! initially all values are set to zero
                !
                plai(i,1,j)  = 0.0_r8
                plai(i,2,j)  = 0.0_r8
                plai(i,3,j)  = 0.0_r8
                plai(i,4,j)  = 0.0_r8
                plai(i,5,j)  = 0.0_r8
                plai(i,6,j)  = 0.0_r8
                plai(i,7,j)  = 0.0_r8
                plai(i,8,j)  = 0.0_r8
                plai(i,9,j)  = 0.0_r8
                plai(i,10,j) = 0.0_r8
                plai(i,11,j) = 0.0_r8
                plai(i,12,j) = 0.0_r8
                !
                ! ---------------------------------------------------
                !  1: tropical evergreen forest / woodland
                ! ---------------------------------------------------
                !
                IF (inveg.EQ.1) THEN
                   !
                   plaievgr  = 5.00_r8
                   plaideci  = 1.00_r8
                   plaishrub = 0.25_r8
                   plaigrass = 0.25_r8
                   !
                   plai(i,1,j)  = exist(i,1,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,2,j)  = exist(i,2,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,3,j)  = exist(i,3,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,4,j)  = exist(i,4,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,5,j)  = exist(i,5,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,6,j)  = exist(i,6,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,7,j)  = exist(i,7,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,8,j)  = exist(i,8,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,9,j)  = exist(i,9,j)  / REAL(ishrub,kind=r8) * plaishrub
                   plai(i,10,j) = exist(i,10,j) / REAL(ishrub,kind=r8) * plaishrub
                   plai(i,11,j) = exist(i,11,j) / REAL(igrass,kind=r8) * plaigrass
                   plai(i,12,j) = exist(i,12,j) / REAL(igrass,kind=r8) * plaigrass
                   !
                END IF
                !
                ! ---------------------------------------------------
                !  2: tropical deciduous forest / woodland
                ! ---------------------------------------------------
                !
                IF (inveg.EQ.2) THEN
                   !
                   plaievgr  = 1.00_r8
                   plaideci  = 5.00_r8
                   plaishrub = 0.25_r8
                   plaigrass = 0.25_r8
                   !
                   plai(i,1,j)  = exist(i,1,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,2,j)  = exist(i,2,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,3,j)  = exist(i,3,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,4,j)  = exist(i,4,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,5,j)  = exist(i,5,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,6,j)  = exist(i,6,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,7,j)  = exist(i,7,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,8,j)  = exist(i,8,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,9,j)  = exist(i,9,j)  / REAL(ishrub,kind=r8) * plaishrub
                   plai(i,10,j) = exist(i,10,j) / REAL(ishrub,kind=r8) * plaishrub
                   plai(i,11,j) = exist(i,11,j) / REAL(igrass,kind=r8) * plaigrass
                   plai(i,12,j) = exist(i,12,j) / REAL(igrass,kind=r8) * plaigrass
                   !
                END IF
                !
                ! ---------------------------------------------------
                !  3: temperate evergreen broadleaf forest / woodland
                ! ---------------------------------------------------
                !
                IF (inveg.EQ.3) THEN
                   !
                   plaievgr  = 4.00_r8
                   plaideci  = 1.00_r8
                   plaishrub = 0.25_r8
                   plaigrass = 0.25_r8
                   !
                   plai(i,1,j)  = exist(i,1,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,2,j)  = exist(i,2,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,3,j)  = exist(i,3,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,4,j)  = exist(i,4,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,5,j)  = exist(i,5,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,6,j)  = exist(i,6,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,7,j)  = exist(i,7,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,8,j)  = exist(i,8,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,9,j)  = exist(i,9,j)  / REAL(ishrub,kind=r8) * plaishrub
                   plai(i,10,j) = exist(i,10,j) / REAL(ishrub,kind=r8) * plaishrub
                   plai(i,11,j) = exist(i,11,j) / REAL(igrass,kind=r8) * plaigrass
                   plai(i,12,j) = exist(i,12,j) / REAL(igrass,kind=r8) * plaigrass
                   !
                END IF
                !
                ! ---------------------------------------------------
                !  4: temperate evergreen conifer forest / woodland
                ! ---------------------------------------------------
                !
                IF (inveg.EQ.4) THEN
                   !
                   plaievgr  = 3.00_r8
                   plaideci  = 1.00_r8
                   plaishrub = 0.25_r8
                   plaigrass = 0.25_r8
                   !
                   plai(i,1,j)  = exist(i,1,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,2,j)  = exist(i,2,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,3,j)  = exist(i,3,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,4,j)  = exist(i,4,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,5,j)  = exist(i,5,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,6,j)  = exist(i,6,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,7,j)  = exist(i,7,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,8,j)  = exist(i,8,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,9,j)  = exist(i,9,j)  / REAL(ishrub,kind=r8) * plaishrub
                   plai(i,10,j) = exist(i,10,j) / REAL(ishrub,kind=r8) * plaishrub
                   plai(i,11,j) = exist(i,11,j) / REAL(igrass,kind=r8) * plaigrass
                   plai(i,12,j) = exist(i,12,j) / REAL(igrass,kind=r8) * plaigrass
                   !
                END IF
                !
                ! ---------------------------------------------------
                !  5: temperate deciduous forest / woodland
                ! ---------------------------------------------------
                !
                IF (inveg.EQ.5) THEN
                   !
                   plaievgr  = 1.00_r8
                   plaideci  = 3.00_r8
                   plaishrub = 0.25_r8
                   plaigrass = 0.25_r8
                   !
                   plai(i,1,j)  = exist(i,1,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,2,j)  = exist(i,2,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,3,j)  = exist(i,3,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,4,j)  = exist(i,4,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,5,j)  = exist(i,5,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,6,j)  = exist(i,6,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,7,j)  = exist(i,7,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,8,j)  = exist(i,8,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,9,j)  = exist(i,9,j)  / REAL(ishrub,kind=r8) * plaishrub
                   plai(i,10,j) = exist(i,10,j) / REAL(ishrub,kind=r8) * plaishrub
                   plai(i,11,j) = exist(i,11,j) / REAL(igrass,kind=r8) * plaigrass
                   plai(i,12,j) = exist(i,12,j) / REAL(igrass,kind=r8) * plaigrass
                   !
                END IF
                !
                ! ---------------------------------------------------
                !  6: boreal evergreen forest / woodland
                ! ---------------------------------------------------
                !
                IF (inveg.EQ.6) THEN
                   !
                   plaievgr  = 3.00_r8
                   plaideci  = 1.00_r8
                   plaishrub = 0.25_r8
                   plaigrass = 0.25_r8
                   !
                   plai(i,1,j)  = exist(i,1,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,2,j)  = exist(i,2,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,3,j)  = exist(i,3,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,4,j)  = exist(i,4,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,5,j)  = exist(i,5,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,6,j)  = exist(i,6,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,7,j)  = exist(i,7,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,8,j)  = exist(i,8,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,9,j)  = exist(i,9,j)  / REAL(ishrub,kind=r8) * plaishrub
                   plai(i,10,j) = exist(i,10,j) / REAL(ishrub,kind=r8) * plaishrub
                   plai(i,11,j) = exist(i,11,j) / REAL(igrass,kind=r8) * plaigrass
                   plai(i,12,j) = exist(i,12,j) / REAL(igrass,kind=r8) * plaigrass
                   !
                END IF
                !
                ! ---------------------------------------------------
                !  7: boreal deciduous forest / woodland
                ! ---------------------------------------------------
                !
                IF (inveg.EQ.7) THEN
                   !
                   plaievgr  = 1.00_r8
                   plaideci  = 3.00_r8
                   plaishrub = 0.25_r8
                   plaigrass = 0.25_r8
                   !
                   plai(i,1,j)  = exist(i,1,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,2,j)  = exist(i,2,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,3,j)  = exist(i,3,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,4,j)  = exist(i,4,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,5,j)  = exist(i,5,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,6,j)  = exist(i,6,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,7,j)  = exist(i,7,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,8,j)  = exist(i,8,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,9,j)  = exist(i,9,j)  / REAL(ishrub,kind=r8) * plaishrub
                   plai(i,10,j) = exist(i,10,j) / REAL(ishrub,kind=r8) * plaishrub
                   plai(i,11,j) = exist(i,11,j) / REAL(igrass,kind=r8) * plaigrass
                   plai(i,12,j) = exist(i,12,j) / REAL(igrass,kind=r8) * plaigrass
                   !
                END IF
                !
                ! ---------------------------------------------------
                !  8: mixed forest / woodland
                ! ---------------------------------------------------
                !
                IF (inveg.EQ.8) THEN
                   !
                   plaievgr  = 2.00_r8
                   plaideci  = 2.00_r8
                   plaishrub = 0.25_r8
                   plaigrass = 0.25_r8
                   !
                   plai(i,1,j)  = exist(i,1,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,2,j)  = exist(i,2,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,3,j)  = exist(i,3,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,4,j)  = exist(i,4,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,5,j)  = exist(i,5,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,6,j)  = exist(i,6,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,7,j)  = exist(i,7,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,8,j)  = exist(i,8,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,9,j)  = exist(i,9,j)  / REAL(ishrub,kind=r8) * plaishrub
                   plai(i,10,j) = exist(i,10,j) / REAL(ishrub,kind=r8) * plaishrub
                   plai(i,11,j) = exist(i,11,j) / REAL(igrass,kind=r8) * plaigrass
                   plai(i,12,j) = exist(i,12,j) / REAL(igrass,kind=r8) * plaigrass
                   !
                END IF
                !
                ! ---------------------------------------------------
                !  9: savanna
                ! ---------------------------------------------------
                !
                IF (inveg.EQ.9) THEN
                   !
                   plaievgr  = 0.50_r8
                   plaideci  = 1.00_r8
                   plaishrub = 0.50_r8
                   plaigrass = 2.00_r8
                   !
                   plai(i,1,j)  = exist(i,1,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,2,j)  = exist(i,2,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,3,j)  = exist(i,3,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,4,j)  = exist(i,4,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,5,j)  = exist(i,5,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,6,j)  = exist(i,6,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,7,j)  = exist(i,7,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,8,j)  = exist(i,8,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,9,j)  = exist(i,9,j)  / REAL(ishrub,kind=r8) * plaishrub
                   plai(i,10,j) = exist(i,10,j) / REAL(ishrub,kind=r8) * plaishrub
                   !
                   ! determine if c4/c3 dominated grassland
                   !
                   IF (tw(i,j).GT.22.0) THEN
                      plai(i,11,j) = exist(i,11,j) * 0.80_r8 * plaigrass
                      plai(i,12,j) = exist(i,12,j) * 0.20_r8 * plaigrass
                   ELSE
                      plai(i,11,j) = exist(i,11,j) * 0.00_r8 * plaigrass
                      plai(i,12,j) = exist(i,12,j) * 1.00_r8 * plaigrass
                   END IF
                   !
                END IF
                !
                ! ---------------------------------------------------
                ! 10: grassland / steppe
                ! ---------------------------------------------------
                !
                IF (inveg.EQ.10) THEN
                   !
                   plaievgr  = 0.25_r8
                   plaideci  = 0.25_r8
                   plaishrub = 0.50_r8
                   plaigrass = 2.50_r8
                   !
                   plai(i,1,j)  = exist(i,1,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,2,j)  = exist(i,2,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,3,j)  = exist(i,3,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,4,j)  = exist(i,4,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,5,j)  = exist(i,5,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,6,j)  = exist(i,6,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,7,j)  = exist(i,7,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,8,j)  = exist(i,8,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,9,j)  = exist(i,9,j)  / REAL(ishrub,kind=r8) * plaishrub
                   plai(i,10,j) = exist(i,10,j) / REAL(ishrub,kind=r8) * plaishrub
                   !
                   ! determine if c4/c3 dominated savanna
                   !
                   IF (tw(i,j).GT.22.0_r8) THEN
                      plai(i,11,j) = exist(i,11,j) * 0.80_r8 * plaigrass
                      plai(i,12,j) = exist(i,12,j) * 0.20_r8 * plaigrass
                   ELSE
                      plai(i,11,j) = exist(i,11,j) * 0.00_r8 * plaigrass
                      plai(i,12,j) = exist(i,12,j) * 1.00_r8 * plaigrass
                   END IF
                   !
                END IF
                !
                ! ---------------------------------------------------
                ! 11: dense shrubland
                ! ---------------------------------------------------
                !
                IF (inveg.EQ.11) THEN
                   !
                   plaievgr  = 0.10_r8
                   plaideci  = 0.10_r8
                   plaishrub = 1.00_r8
                   plaigrass = 0.50_r8
                   !
                   plai(i,1,j)  = exist(i,1,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,2,j)  = exist(i,2,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,3,j)  = exist(i,3,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,4,j)  = exist(i,4,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,5,j)  = exist(i,5,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,6,j)  = exist(i,6,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,7,j)  = exist(i,7,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,8,j)  = exist(i,8,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,9,j)  = exist(i,9,j)  / REAL(ishrub,kind=r8) * plaishrub
                   plai(i,10,j) = exist(i,10,j) / REAL(ishrub,kind=r8) * plaishrub
                   plai(i,11,j) = exist(i,11,j) / REAL(igrass,kind=r8) * plaigrass
                   plai(i,12,j) = exist(i,12,j) / REAL(igrass,kind=r8) * plaigrass
                   !
                END IF
                !
                ! ---------------------------------------------------
                ! 12: open shrubland
                ! ---------------------------------------------------
                !
                IF (inveg.EQ.12) THEN
                   !
                   plaievgr  = 0.00_r8
                   plaideci  = 0.00_r8
                   plaishrub = 0.25_r8
                   plaigrass = 0.25_r8
                   !
                   plai(i,1,j)  = exist(i,1,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,2,j)  = exist(i,2,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,3,j)  = exist(i,3,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,4,j)  = exist(i,4,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,5,j)  = exist(i,5,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,6,j)  = exist(i,6,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,7,j)  = exist(i,7,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,8,j)  = exist(i,8,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,9,j)  = exist(i,9,j)  / REAL(ishrub,kind=r8) * plaishrub
                   plai(i,10,j) = exist(i,10,j) / REAL(ishrub,kind=r8) * plaishrub
                   plai(i,11,j) = exist(i,11,j) / REAL(igrass,kind=r8) * plaigrass
                   plai(i,12,j) = exist(i,12,j) / REAL(igrass,kind=r8) * plaigrass
                   !
                END IF
                !
                ! ---------------------------------------------------
                ! 13: tundra
                ! ---------------------------------------------------
                !
                IF (inveg.EQ.13) THEN
                   !
                   plaievgr  = 0.00_r8
                   plaideci  = 0.00_r8
                   plaishrub = 1.00_r8
                   plaigrass = 1.00_r8
                   !
                   plai(i,1,j)  = exist(i,1,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,2,j)  = exist(i,2,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,3,j)  = exist(i,3,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,4,j)  = exist(i,4,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,5,j)  = exist(i,5,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,6,j)  = exist(i,6,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,7,j)  = exist(i,7,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,8,j)  = exist(i,8,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,9,j)  = exist(i,9,j)  / REAL(ishrub,kind=r8) * plaishrub
                   plai(i,10,j) = exist(i,10,j) / REAL(ishrub,kind=r8) * plaishrub
                   plai(i,11,j) = exist(i,11,j) / REAL(igrass,kind=r8) * plaigrass
                   plai(i,12,j) = exist(i,12,j) / REAL(igrass,kind=r8) * plaigrass
                   !
                END IF
                !
                ! ---------------------------------------------------
                ! 14: desert
                ! ---------------------------------------------------
                !
                IF (inveg.EQ.14) THEN
                   !
                   plaievgr  = 0.00_r8
                   plaideci  = 0.00_r8
                   plaishrub = 0.05_r8
                   plaigrass = 0.05_r8
                   !
                   plai(i,1,j)  = exist(i,1,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,2,j)  = exist(i,2,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,3,j)  = exist(i,3,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,4,j)  = exist(i,4,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,5,j)  = exist(i,5,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,6,j)  = exist(i,6,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,7,j)  = exist(i,7,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,8,j)  = exist(i,8,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,9,j)  = exist(i,9,j)  / REAL(ishrub,kind=r8) * plaishrub
                   plai(i,10,j) = exist(i,10,j) / REAL(ishrub,kind=r8) * plaishrub
                   plai(i,11,j) = exist(i,11,j) / REAL(igrass,kind=r8) * plaigrass
                   plai(i,12,j) = exist(i,12,j) / REAL(igrass,kind=r8) * plaigrass
                   !
                END IF
                !
                ! ---------------------------------------------------
                ! 15: polar desert / rock / ice
                ! ---------------------------------------------------
                !
                IF (inveg.EQ.15) THEN
                   !
                   plaievgr  = 0.00_r8
                   plaideci  = 0.00_r8
                   plaishrub = 0.05_r8
                   plaigrass = 0.05_r8
                   !
                   plai(i,1,j)  = exist(i,1,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,2,j)  = exist(i,2,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,3,j)  = exist(i,3,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,4,j)  = exist(i,4,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,5,j)  = exist(i,5,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,6,j)  = exist(i,6,j)  / REAL(ievgr,kind=r8)  * plaievgr
                   plai(i,7,j)  = exist(i,7,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,8,j)  = exist(i,8,j)  / REAL(ideci,kind=r8)  * plaideci
                   plai(i,9,j)  = exist(i,9,j)  / REAL(ishrub,kind=r8) * plaishrub
                   plai(i,10,j) = exist(i,10,j) / REAL(ishrub,kind=r8) * plaishrub
                   plai(i,11,j) = exist(i,11,j) / REAL(igrass,kind=r8) * plaigrass
                   plai(i,12,j) = exist(i,12,j) / REAL(igrass,kind=r8) * plaigrass
                   !
                END IF
                !
             END IF
9999         CONTINUE
             !
             ! ************************************************************************
             ! case (2) assign vegetation characteristics for dynamic vegtation
             !          that is initialized with uniform vegetation conditions
             ! ************************************************************************
             !
             ! specify uniform initial conditions
             !
             IF (isimveg.EQ.2) THEN
                !
                !            plaiupper = 0.5
                !            plailower = 0.5
                !
                plai(i,1,j)  = exist(i,1,j)  / REAL(iupper,kind=r8) * plaiupper
                plai(i,2,j)  = exist(i,2,j)  / REAL(iupper,kind=r8) * plaiupper
                plai(i,3,j)  = exist(i,3,j)  / REAL(iupper,kind=r8) * plaiupper
                plai(i,4,j)  = exist(i,4,j)  / REAL(iupper,kind=r8) * plaiupper
                plai(i,5,j)  = exist(i,5,j)  / REAL(iupper,kind=r8) * plaiupper
                plai(i,6,j)  = exist(i,6,j)  / REAL(iupper,kind=r8) * plaiupper
                plai(i,7,j)  = exist(i,7,j)  / REAL(iupper,kind=r8) * plaiupper
                plai(i,8,j)  = exist(i,8,j)  / REAL(iupper,kind=r8) * plaiupper
                plai(i,9,j)  = exist(i,9,j)  / REAL(ilower,kind=r8) * plailower
                plai(i,10,j) = exist(i,10,j) / REAL(ilower,kind=r8) * plailower
                plai(i,11,j) = exist(i,11,j) / REAL(ilower,kind=r8) * plailower
                plai(i,12,j) = exist(i,12,j) / REAL(ilower,kind=r8) * plailower
                !
             END IF

             !
             ! ************************************************************************
             ! for both cases (1) and (2)
             ! ************************************************************************
             !
             ! set minimum lai for each existing plant type
             !
             !          xminlai = 0.010
             !
             plai(i,1,j)  = MAX (plai(i,1,j) , exist(i,1,j)  * xminlai)
             plai(i,2,j)  = MAX (plai(i,2,j) , exist(i,2,j)  * xminlai)
             plai(i,3,j)  = MAX (plai(i,3,j) , exist(i,3,j)  * xminlai)
             plai(i,4,j)  = MAX (plai(i,4,j) , exist(i,4,j)  * xminlai)
             plai(i,5,j)  = MAX (plai(i,5,j) , exist(i,5,j)  * xminlai)
             plai(i,6,j)  = MAX (plai(i,6,j) , exist(i,6,j)  * xminlai)
             plai(i,7,j)  = MAX (plai(i,7,j) , exist(i,7,j)  * xminlai)
             plai(i,8,j)  = MAX (plai(i,8,j) , exist(i,8,j)  * xminlai)
             plai(i,9,j)  = MAX (plai(i,9,j) , exist(i,9,j)  * xminlai)
             plai(i,10,j) = MAX (plai(i,10,j), exist(i,10,j) * xminlai)
             plai(i,11,j) = MAX (plai(i,11,j), exist(i,11,j) * xminlai)
             plai(i,12,j) = MAX (plai(i,12,j), exist(i,12,j) * xminlai)
             !
             !
             ! set sapwood fraction and biomass characteristics
             !
             !          sapfrac(i) = 0.1
             sapfrac(i,j) = sapfrac_init ! 0.1 from params.veg
             !
             wood = 0.0_r8
             !
             DO k = 1, npft
                !
                cbiol(i,k,j) = plai(i,k,j) / specla(k)
                cbior(i,k,j) = 0.5_r8 * cbiol(i,k,j)
                !
                cbiow(i,k,j) = 0.0_r8
                !
                IF (k.LT.9) cbiow(i,k,j) = plai(i,k,j) * 10.0_r8 / 6.0_r8
                !
                biomass(i,k,j) = cbiol(i,k,j) + cbiow(i,k,j) + cbior(i,k,j) 
                IF (K.LE.8) wood = wood + cbiow(i,k,j)
                !
             END DO
             !
             ! ************************************************************************
             ! determine basic vegetation structure characteristics
             ! ************************************************************************
             ! 
             ! total leaf area for upper and lower canopies
             !
             totlaiu(i,j)  =  plai(i,1,j) + plai(i,2,j) +  &
                  plai(i,3,j) + plai(i,4,j) +  &
                  plai(i,5,j) + plai(i,6,j) +  &
                  plai(i,7,j) + plai(i,8,j) 
             !
             totlail(i,j)  =  plai(i,9,j)  + plai(i,10,j) + &
                  plai(i,11,j) + plai(i,12,j) 
             !
             totbiou(i,j)  = biomass(i,1,j) + &
                  biomass(i,2,j) + &
                  biomass(i,3,j) + &
                  biomass(i,4,j) + &
                  biomass(i,5,j) + &
                  biomass(i,6,j) + &
                  biomass(i,7,j) + &
                  biomass(i,8,j) 
             !
             totbiol(i,j)  = biomass(i,9,j)  + &
                  biomass(i,10,j) + &
                  biomass(i,11,j) + &
                  biomass(i,12,j) 
             !
             ! fractional cover
             !
             !       fu(i) = wood / woodnorm
             !
             fu(i,j) = (1.0_r8 - EXP(-wood)) / (1.0_r8 - EXP(-woodnorm))
             !
             fl(i,j) = totlail(i,j) / 1.0_r8
             !
             fu(i,j) = MAX (0.25_r8, MIN (0.975_r8, fu(i,j)))
             fl(i,j) = MAX (0.25_r8, MIN (0.975_r8, fl(i,j)))
             !
             ! initial single-sided sai for upper and lower canopies
             !
             sai(i,1,j)  =  0.050_r8 * totlail(i,j)
             sai(i,2,j)  =  0.250_r8 * totlaiu(i,j)

             !
             ! initial lai for canopy physics
             !
             lai(i,1,j) = totlail(i,j) / fl(i,j)
             lai(i,2,j) = totlaiu(i,j) / fu(i,j)
             !
             ! specify canopy height parameters
             ! calculated as a function of only the vegetative fraction
             ! of each grid cell
             !
             zbot(i,1,j) =  0.05_r8
             !       ztop(i,1,j) =  max (0.25_r8, totlail(i) * 0.25_r8)
             ztop(i,1,j) =  MAX (0.25_r8, lai(i,1,j) * 0.25_r8)
             !
             !PK        zbot(i,2,j) =  ztop(i,1,j) + 1.0_r8 
             !       ztop(i,2,j) =  max (zbot(i,2,j) + 1.00_r8, 2.50_r8 * totbiou(i,j) * 0.75_r8)
             !PK        ztop(i,2,j) =  max (zbot(i,2,j) + 1.00_r8,    &
             !PK                          2.50_r8 * totbiou(i,j) / fu(i,j) * 0.75_r8)
             !
             zbot(i,2,j) =  3.0_r8
             ztop(i,2,j) =  MAX (zbot(i,2,j) + 1.0_r8, 2.5_r8*totbiou(i,j) / fu(i,j) * 0.75_r8)
             !
             ! constrain ztop of lower canopy to be at least 0.5 meter lower than
             ! zbot for upper canopy
             !
             ztop(i,1,j) = MIN (ztop(i,1,j), zbot(i,2,j) - 0.5_r8)

             !
             ! ************************************************************************
             ! case (3) assign vegetation characteristics for model restart
             ! ************************************************************************
             !
             !        ELSE    ! else for restart if loop
             !
             !          DO k = 1, npft
             !
             !            plai(i,k,j) = cbiol(i,k,j) * specla(k)
             !            biomass(i,k,j) = cbiol(i,k,j) + cbiow(i,k,j) + cbior(i,k,j)
             !
             !          END DO
             ! ************************************************************************
             ! determine basic vegetation structure characteristics
             ! ************************************************************************
             ! 
             ! total leaf area for upper and lower canopies
             !
             !          totlaiu(i,j)  =  plai(i,1,j) + plai(i,2,j) +  &
             !                           plai(i,3,j) + plai(i,4,j) +  &
             !                           plai(i,5,j) + plai(i,6,j) +  &
             !                           plai(i,7,j) + plai(i,8,j) 
             !
             !          totlail(i,j)  =  plai(i,9,j)  + plai(i,10,j) + &
             !                           plai(i,11,j) + plai(i,12,j) 
             !
             !          totbiou(i,j)  = biomass(i,1,j) + &
             !                               biomass(i,2,j) + &
             !                               biomass(i,3,j) + &
             !                               biomass(i,4,j) + &
             !                               biomass(i,5,j) + &
             !                               biomass(i,6,j) + &
             !                               biomass(i,7,j) + &
             !                               biomass(i,8,j) 
             !
             !          totbiol(i,j)  = biomass(i,9,j)  + &
             !                                biomass(i,10,j) + &
             !                               biomass(i,11,j) + &
             !                                biomass(i,12,j) 
             !!
             ! initial single-sided sai for upper and lower canopies
             !
             !          sai(i,1,j)  =  0.050_r8 * totlail(i,j)
             !          sai(i,2,j)  =  0.250_r8 * totlaiu(i,j)
             !
             ! Lai read from restart file
             !
             ! specify canopy height parameters
             ! calculated as a function of only the vegetative fraction
             ! of each grid cell
             !
             !        zbot(i,1,j) =  0.05_r8
             !       ztop(i,1,j) =  max (0.25_r8, totlail(i) * 0.25_r8)
             !        ztop(i,1,j) =  max (0.25_r8, lai(i,1,j) * 0.25_r8)
             !
             !        zbot(i,2,j) =  3.0_r8
             !        ztop(i,2,j) =  max (zbot(i,2,j) + 1.0_r8, 2.5_r8*totbiou(i,j) / fu(i,j) * 0.75_r8)
             !
             ! constrain ztop of lower canopy to be at least 0.5 meter lower than
             ! zbot for upper canopy
             !
             !        ztop(i,1,j) = min (ztop(i,1,j), zbot(i,2,j) - 0.5_r8)

          END IF  ! end restart if loop
          !
       END DO !DO i=1,npoi
    END DO

    !
    ! ************************************************************************
    ! assign some physical properties of vegetation
    ! ************************************************************************
    !
    ! leaf optical properties were taken from Sellers et al., 1996
    ! and Bonan, 1995
    !
    !      rhoveg(1,1) = 0.10     ! vis leaf reflectance, lower story
    !      rhoveg(1,2) = 0.10     ! vis leaf reflectance, upper story 
    !
    !      rhoveg(2,1) = 0.60     ! nir leaf reflectance, lower story
    !      rhoveg(2,2) = 0.40     ! nir leaf reflectance, upper story
    !
    !      tauveg(1,1) = 0.07     ! vis leaf transmittance, lower story
    !      tauveg(1,2) = 0.05     ! vis leaf transmittance, upper story
    !
    !      tauveg(2,1) = 0.25     ! nir leaf transmittance, lower story
    !      tauveg(2,2) = 0.20     ! nir leaf transmittance, upper story
    !
    !      chiflz = -0.5          ! leaf orientation factors (-1 vertical, 0 random, 1 horizontal)
    !      chifuz =  0.0          ! leaf orientation factors (-1 vertical, 0 random, 1 horizontal)
    !
    oriev(1) = MAX (-chiflz, 0.0_r8)
    oriev(2) = MAX (-chifuz, 0.0_r8)
    !
    orieh(1) = MAX ( chiflz, 0.0_r8)
    orieh(2) = MAX ( chifuz, 0.0_r8)
    !
    !      dleaf(1) = 0.10        ! linear dimensions for aerodynamic flux parameterization
    !      dstem(1) = 0.10        ! linear dimensions for aerodynamic flux parameterization
    !
    !      dleaf(2) = 0.10        ! linear dimensions for aerodynamic flux parameterization
    !      dstem(2) = 0.10        ! linear dimensions for aerodynamic flux parameterization
    !
    !      chu = ch2o *  2.0      ! heat capacity of upper leaves
    !      chl = ch2o *  2.0      ! heat capacity of lower leaves
    !      chs = ch2o * 50.0      ! heat capacity of stems
    !
    !      alaimu = 8.0           ! normalization constant for upper canopy aerodynamics
    !      alaiml = 8.0           ! normalization constant for lower canopy aerodynamics
    !
    !      cleaf  = 0.01          ! constant in leaf-air aero transfer parameterization
    !      cgrass = 0.01          ! constant in leaf-air aero transfer parameterization
    !      cstem  = 0.01          ! constant in leaf-air aero transfer parameterization
    !
    !      wliqumax = 0.20        ! intercepted water capacity (mm h2o per unit leaf area)
    !      wliqsmax = 0.40        ! intercepted water capacity (mm h2o per unit leaf area)
    !      wliqlmax = 0.20        ! intercepted water capacity (mm h2o per unit leaf area)
    !
    !      wsnoumax = 2.00        ! intercepted snow capacity (mm h2o per unit leaf area)
    !      wsnosmax = 4.00        ! intercepted snow capacity (mm h2o per unit leaf area)
    !      wsnolmax = 2.00        ! intercepted snow capacity (mm h2o per unit leaf area)
    !
    !      tdripu =  2.0 * 3600.0 ! decay time for intercepted liquid dripoff (sec)
    !      tdrips =  2.0 * 3600.0 ! decay time for intercepted liquid dripoff (sec)
    !      tdripl =  2.0 * 3600.0 ! decay time for intercepted liquid dripoff (sec)
    !
    !      tblowu = 12.0 * 3600.0 ! decay time for snow blowoff (sec)
    !      tblows = 12.0 * 3600.0 ! decay time for snow blowoff (sec)
    !      tblowl = 12.0 * 3600.0 ! decay time for snow blowoff (sec)
    !
    ! ************************************************************************
    ! define rooting profiles
    ! ************************************************************************
    !
    ! define rooting profiles based upon data published in:
    !
    ! Jackson et al., 1996:  A global analysis of root distributions
    ! for terrestrial biomes, Oecologia, 108, 389-411.
    !
    ! and
    !
    ! Jackson et al., 1997:  A global budget for fine root biomass, 
    ! surface area, and nutrient contents, Proceedings of the National
    ! Academy of Sciences, 94, 7362-7366.
    !
    ! rooting profiles are defined by the "beta" parameter
    !
    ! beta1 is assigned to the lower vegetation layer (grasses and shrubs)
    ! beta2 is assigned to the upper vegetation layer (trees)
    !
    ! according to Jackson et al. (1996, 1997), the values of beta
    ! typically fall in the following range
    !
    ! note that the 1997 paper specifically discusses the distribution
    ! of *fine roots* (instead of total root biomass), which may be more
    ! important for water and nutrient uptake
    !
    ! --------------                 ------------   ------------
    ! forest systems                 beta2 (1996)   beta2 (1997)
    ! --------------                 ------------   ------------
    ! tropical evergreen forest:        0.962          0.972
    ! tropical deciduous forest:        0.961          0.982
    ! temperate conifer forest:         0.976          0.980
    ! temperate broadleaf forest:       0.966          0.967
    ! all tropical/temperate forest:    0.970  
    ! boreal forest:                    0.943          0.943
    ! all trees:                                       0.976
    !
    ! -------------------------      ------------   ------------
    ! grassland / shrub systems      beta1 (1996)   beta1 (1997)
    ! -------------------------      ------------   ------------
    ! tropical grassland / savanna:     0.972          0.972
    ! temperate grassland:              0.943          0.943
    ! all grasses:                      0.952          0.952
    ! schlerophyllous shrubs:           0.964          0.950
    ! all shrubs:                       0.978          0.975
    ! crops:                            0.961
    ! desert:                           0.975          0.970
    ! tundra:                           0.914
    !
    ! --------------                 ------------
    ! all ecosystems                 beta  (1996)
    ! --------------                 ------------
    ! all ecosystems:                   0.966
    !
    ! for global simulations, we typically assign the following
    ! values to the beta parameters
    !
    ! beta1 = 0.950, which is typical for tropical/temperate grasslands
    ! beta2 = 0.970, which is typical for tropical/temperate forests
    !
    ! however, these values could be (and should be) further refined
    ! when using the model for specific regions
    ! 
    !      beta1 = 0.950  ! for lower layer herbaceous plants
    !      beta2 = 0.975  ! for upper layer trees
    !
    ! calculate total depth in centimeters
    !
    totdepth = 0.0_r8
    !
    DO k = 1, nsoilay
       totdepth = totdepth + hsoi(k) * 100.0_r8
    END DO
    !
    ! normalization factors
    !
    frootnorm1 = 1.0_r8 - beta1 ** totdepth
    frootnorm2 = 1.0_r8 - beta2 ** totdepth
    !
    ! calculate rooting profiles
    !
    DO k = 1, nsoilay
       !
       IF (k.EQ.1) THEN
          !
          depth(k) = hsoi(k) * 100.0_r8
          !
          froot(k,1) = 1.0_r8 - beta1 ** depth(k)
          froot(k,2) = 1.0_r8 - beta2 ** depth(k)
          !
       ELSE
          !
          depth(k) = depth(k-1) + hsoi(k) * 100.0_r8
          !
          froot(k,1) = (1.0_r8 - beta1 ** depth(k)) -  &
               (1.0_r8 - beta1 ** depth(k-1)) 
          !
          froot(k,2) = (1.0_r8 - beta2 ** depth(k)) -   & 
               (1.0_r8 - beta2 ** depth(k-1)) 
          !
       END IF
       !
       froot(k,1) = froot(k,1) / frootnorm1
       froot(k,2) = froot(k,2) / frootnorm2
       !
    END DO
    !
    ! return to main program
    !
    RETURN
  END SUBROUTINE iniveg
  !


  !
  ! #####   #####    ##    #####   #####     ##    #####    ####   
  ! #    #  #       #  #   #    #  #    #   #  #   #    #  #      
  ! #    #  #      #    #  #    #  #    #  #    #  #    #   ####       
  ! #####   #####  ######  #    #  #####   ######  #####        #  
  ! #   #   #      #    #  #    #  #       #    #  #   #   #    # 
  ! #    #  #####  #    #  #####   #       #    #  #    #   ####  
  !

  !----------------------------------------------------------
  SUBROUTINE RD_PARAM()
    !----------------------------------------------------------
    !
    !  Read various parameters from ibis.params
    ! 
    IMPLICIT NONE
    !  
    ! Local variables
    !      
    INTEGER, PARAMETER :: parm_unit =9  ! file unit assignment for input
    INTEGER :: npft2       ! number of PFTs reported in params.veg
    INTEGER :: npftu2      ! number of upper canopy PFTs reported in params.veg
    INTEGER :: nsoil2      ! number of soil texture classes reported in params.soi
    INTEGER :: nveg            ! number of vegetation classes reported in params.veg
    CHARACTER(LEN=20) :: parm_file
    !      co2init =0.000350_r8  ! co2init    initial co2 concentration in mol/mol (REAL(KIND=r8))
    !      o2init  =0.209000_r8  ! o2init     initial o2 concentration in mol/mol (real)

    parm_file='file_par'
    !------------------------------------------------------------------------
    !  Fundamental plant physiological parameters, name definitions
    ! ------------------------------------------------------------------------
    !      tau15 = 4500.0_r8    ! tau15 : co2/o2 specificity ratio at 15 degrees C (dimensionless)    
    !      kc15  = 1.5e-04_r8   ! kc15  : co2 kinetic parameter at 15 C (mol/mol)
    !      ko15  = 2.5e-01_r8   ! ko15  : o2 kinetic parameter at 15 C (mol/mol)
    !      cimax = 2000.e-06_r8 ! cimax : maximum value for ci (for model stability)

    !------------------------------------------------------------------------
    !========================================================================

    !------------------------------------------------------------------------
    npft2 =12    ! Number of PFTs in this parameter set   
    npftu2= 8    ! Number of upper canopy (tree) PFTs in this set

    IF (npft2 .NE. npft) THEN
       WRITE (nfprt, 9003) parm_file, npft2, npft
       GOTO 9006 ! In the circumstances this seems the best thing to do! 
    END IF
9003 FORMAT ('RD_PARAM Warning: Number of PFTs in ', A10, ' is: ', &
         I2, ' number in compar.h is: ', I2)  

    IF (npftu2 .NE. npftu) THEN
       WRITE (nfprt,9004) parm_file, npftu2, npftu
       GOTO 9006 ! In the circumstances this seems the best thing to do! 
    END IF
9004 FORMAT ('RD_PARAM Warning: Number of upper canopy (tree) PFTs ', &
         'in ', A10, ' is: ',  &
         I2, ' number in comage.h is: ', I2)  

    !------------------------------------------------------------------------
    ! PFTs (top to bottom)
    !------------------------------------------------------------------------
    !  1: tropical broadleaf evergreen trees
    !  2: tropical broadleaf drought-deciduous trees
    !  3: warm-temperate broadleaf evergreen trees
    !  4: temperate conifer evergreen trees
    !  5: temperate broadleaf cold-deciduous trees
    !  6: boreal conifer evergreen trees
    !  7: boreal broadleaf cold-deciduous trees
    !  8: boreal conifer cold-deciduous trees
    !  9: evergreen shrubs
    ! 10: cold-deciduous shrubs
    ! 11: warm (C4) grasses
    ! 12: cool (C3) grasses
    !========================================================================

    !--------------------------------------------------------------------
    ! C3 and C4 physiology-specific parameters
    !--------------------------------------------------------------------
    !      alpha3 =0.080_r8 ! alpha3 - C3 intrinsic quantum efficiency (dimensionless)
    !      theta3 =0.950_r8 ! theta3 - C3 photosynthesis coupling coefficient
    !      beta3  =0.990_r8 ! beta3  - C3 photosynthesis coupling coefficient
    !      alpha4 =0.050_r8 ! alpha4 - C4 intrinsic quantum efficiency (dimensionless)
    !      theta4 =0.970_r8 ! theta4 - C4 photosynthesis coupling coefficient
    !      beta4  =0.800_r8 ! beta4  - C4 photosynthesis coupling coefficient 
    !====================================================================
    !--------------------------------------------------------------------
    ! Plant physiological properties - 5 classes
    !--------------------------------------------------------------------
    ! gamma    : leaf respiration coefficients 
    ! coefm    : 'm' coefficients for stomatal conductance relationship
    ! coefb    : 'b' coefficients for stomatal conductance relationship
    ! gsmin    : absolute minimum stomatal conductances
    !-------------------------------------------------
    ! gamma  coefm  coefb    gsmin   Physiol. Class
    !-------------------------------------------------
    !      gammaub = 0.015_r8   ! Broadleaf trees
    !      coefmub = 10.0_r8    ! Broadleaf trees
    !      coefbub = 0.01_r8    ! Broadleaf trees
    !      gsubmin = 0.00001_r8 ! Broadleaf trees
    !      gammauc = 0.015_r8
    !      coefmuc = 6.0_r8
    !      coefbuc = 0.01_r8
    !      gsucmin = 0.00001_r8
    !      gammals = 0.015_r8  ! Shrubs
    !      coefmls = 9.0_r8    ! Shrubs
    !      coefbls = 0.01_r8   ! Shrubs
    !      gslsmin = 0.00001_r8! Shrubs
    !      gammal4 = 0.030_r8  ! C4 grasses
    !      coefml4 = 4.0_r8    ! C4 grasses
    !      coefbl4 = 0.04_r8   ! C4 grasses
    !      gsl4min = 0.00001_r8! C4 grasses
    !      gammal3 = 0.015_r8  ! C3 grasses
    !      coefml3 = 9.0_r8    ! C3 grasses
    !      coefbl3 = 0.01_r8   ! C3 grasses
    !      gsl3min = 0.00001_r8! C3 grasses

    !=================================================
    !--------------------------------------------------
    ! Other properties of vegetation -- for 12 PFTs
    !--------------------------------------------------
    ! vmax_pft : max Rubisco activity at 15 C, at top of canopy (mol[CO2] m-2 s-1) 
    ! specla   : specific leaf area (m2 kg-1)
    ! tauleaf  : foliar biomass turnover time constant (years)
    ! tauroot  : root biomass turnover time constant (years)
    ! tauwood  : wood biomass turnover time constant (years)
    ! aleaf    : foliar allocation coefficient (fraction)
    ! aroot    : root allocation coefficient (fraction)
    ! awood    : wood allocation coefficient (fraction, = 1 - aleaf - aroot)
    !      dummyvarpk(1:96)=(/&
    !------------------------------------------------------------------------
    ! vmax_pft  specla  tauleaf  tauroot tauwood  aleaf  aroot  awood    PFT
    !------------------------------------------------------------------------
    !      65.0e-06_r8, 25.0_r8, 1.01_r8, 1.0_r8,  25.0_r8, 0.30_r8, 0.20_r8, 0.50_r8,& !   1 
    !      65.0e-06_r8, 25.0_r8, 1.00_r8, 1.0_r8,  25.0_r8, 0.30_r8, 0.20_r8, 0.50_r8,& !   2 
    !      40.0e-06_r8, 25.0_r8, 1.00_r8, 1.0_r8,  25.0_r8, 0.30_r8, 0.20_r8, 0.50_r8,& !   3 
    !      30.0e-06_r8, 12.5_r8, 2.00_r8, 1.0_r8,  50.0_r8, 0.30_r8, 0.40_r8, 0.30_r8,& !   4 
    !      30.0e-06_r8, 25.0_r8, 1.00_r8, 1.0_r8,  50.0_r8, 0.30_r8, 0.20_r8, 0.50_r8,& !   5 
    !      25.0e-06_r8, 12.5_r8, 2.50_r8, 1.0_r8, 100.0_r8, 0.30_r8, 0.40_r8, 0.30_r8,& !   6 
    !      30.0e-06_r8, 25.0_r8, 1.00_r8, 1.0_r8, 100.0_r8, 0.30_r8, 0.20_r8, 0.50_r8,& !   7 
    !      30.0e-06_r8, 25.0_r8, 1.00_r8, 1.0_r8, 100.0_r8, 0.30_r8, 0.20_r8, 0.50_r8,& !   8 
    !      27.5e-06_r8, 12.5_r8, 1.50_r8, 1.0_r8,   5.0_r8, 0.45_r8, 0.40_r8, 0.15_r8,& !   9 
    !      27.5e-06_r8, 25.0_r8, 1.00_r8, 1.0_r8,   5.0_r8, 0.45_r8, 0.35_r8, 0.20_r8,& !  10 
    !      15.0e-06_r8, 20.0_r8, 1.25_r8, 1.0_r8, 999.0_r8, 0.45_r8, 0.55_r8, 0.00_r8,& !  11 
    !      25.0e-06_r8, 20.0_r8, 1.50_r8, 1.0_r8, 999.0_r8, 0.45_r8, 0.55_r8, 0.00_r8/)  !  12 
    !========================================================================
    !      i=0
    !      DO j = 1, npft 
    !      
    !         vmax_pft(j) = dummyvarpk(j+0+i)
    !         specla(j)   = dummyvarpk(j+1+i)
    !         tauleaf(j)  = dummyvarpk(j+2+i)
    !         tauroot(j)  = dummyvarpk(j+3+i)
    !         tauwood0(j) = dummyvarpk(j+4+i)
    !         aleaf(j)    = dummyvarpk(j+5+i)
    !         aroot(j)    = dummyvarpk(j+6+i)
    !         awood(j)    = dummyvarpk(j+7+i)
    !         i=i+7
    !      END DO
    !!---------------------------------------------------------------------------------
    !  minimum density of woody biomass required for upper canopy closure (kg C m-2)
    !---------------------------------------------------------------------------------
    !      woodnorm=7.5_r8 !  woodnorm 
    !================================================================================= 

    ! leaf optical properties from Sellers et al., 1996 and Bonan, 1995
    !      dummyvarpk(1:8)=(/ &
    !-----------------------------------------------------------------------------------------
    ! leaf reflectance (rhoveg) and transmittance (tauveg), visible and NIR, for each canopy
    !-----------------------------------------------------------------------------------------
    !     lower        upper  
    !---------------------------------------------------------
    !      0.10_r8,        0.10_r8,&    ! rhoveg(1,1); rhoveg(1,2) vis   
    !      0.60_r8,        0.40_r8,&    ! rhoveg(2,1); rhoveg(2,2) NIR
    !      0.07_r8,        0.05_r8,&    ! tauveg(1,1); tauveg(1,2) vis
    !      0.25_r8,        0.20_r8 /)   ! tauveg(2,1); tauveg(2,2) NIR  
    !---------------------------------------------------------
    !     0.10        0.60    ! rhoveg(1,1); rhoveg(1,2) vis   
    !     0.10        0.40    ! rhoveg(2,1); rhoveg(2,2) NIR
    !     0.07        0.25    ! tauveg(1,1); tauveg(1,2) vis
    !     0.05        0.20    ! tauveg(2,1); tauveg(2,2) NIR  
    !========================================================= 
    !
    !      i=0
    !      DO j = 1, nband 
    !         rhoveg(j,1) = dummyvarpk(j+0+i)
    !         rhoveg(j,2) = dummyvarpk(j+1+i)
    !         i=i+1
    !      END DO  
    !      i=i+2
    !      DO j = 1, nband 
    !         tauveg(j,1) = dummyvarpk(j+0+i)
    !         tauveg(j,2) = dummyvarpk(j+1+i)
    !         i=i+1
    !      END DO
    ! *********************************************
    ! assign some physical properties of vegetation
    ! *********************************************
    !
    ! leaf optical properties were taken from Sellers et al., 1996
    ! and Bonan, 1995
    !
    rhoveg(1,1) =  0.062_r8!0.10_r8     ! vis leaf reflectance, lower story
    rhoveg(1,2) =  0.062_r8 !0.10_r8     ! vis leaf reflectance, upper story 
    !
    rhoveg(2,1) = 0.60_r8     ! nir leaf reflectance, lower story
    rhoveg(2,2) = 0.40_r8     ! nir leaf reflectance, upper story
    !
    tauveg(1,1) = 0.07_r8     ! vis leaf transmittance, lower story
    tauveg(1,2) = 0.05_r8     ! vis leaf transmittance, upper story
    !
    tauveg(2,1) = 0.25_r8     ! nir leaf transmittance, lower story
    tauveg(2,2) = 0.20_r8     ! nir leaf transmittance, upper story
    !-----------------------------------------------------------------
    ! leaf orientation factors (-1 vertical, 0 random, 1 horizontal)
    !-----------------------------------------------------------------
    ! chifuz : upper canopy leaf orientation
    ! chiflz : lower canopy leaf orientation
    !-----------------------------------------
    !      chifuz=   0.0_r8           ! chifuz 
    !      chiflz=  -0.5_r8           ! chiflz              
    !=========================================
    !      dummyvarpk(1:4)=(/ &
    !-----------------------------------------------------------------------
    ! linear dimensions for aerodynamic flux parameterization: dleaf, dstem
    !-----------------------------------------------------------------------
    !    upper?     lower?
    !-----------------------------------------
    !      0.10_r8,         0.10_r8,        &  ! dleaf
    !      0.10_r8,   0.10_r8/)           ! dstem
    !=========================================
    !      i=0
    !      DO j = 1 ,1 
    !         dleaf(1) = dummyvarpk(j+0+i)
    !         dleaf(2) = dummyvarpk(j+1+i)
    !         i=i+1
    !      END DO  
    !      i=i+1
    !      DO j = 1, 1 
    !         dstem(1) = dummyvarpk(j+0+i)
    !         dstem(2) = dummyvarpk(j+1+i)
    !         i=i+1
    !      END DO      
    !
    !--------------------------------------------------------------------
    ! normalization constants for canopy drag coefficients (m2 m-2)
    !---------------------------------------------------------------------
    !      alaimu=8.0_r8   ! alaimu : upper canopy leaf & stem area (2 sided)
    !      alaiml=8.0_r8   ! alaiml : lower canopy leaf & stem maximum area (2 sided)
    !======================================================================

    !----------------------------------------------------------------------------
    ! empirical coefficients for aerodynamic transfer parameterization (m s-0.5) 
    ! From Pollard & Thompson (1995, eq. A39a)
    !----------------------------------------------------------------------------
    !      cleaf  = 0.01_r8! cleaf  : upper canopy leaf-air 
    !      cstem  = 0.01_r8! cstem  : upper canopy stem-air 
    !      cgrass = 0.01_r8! cgrass : lower canopy-air
    !===========================================================================

    !----------------------------------------------------------------------------
    ! heat capacities of leaves and stems  (J kg-1 C-1 m-2)
    ! derived from specific heat of liquid water (ch2o = 4.218 J g-1)
    !----------------------------------------------------------------------------
    !      chs= 2.109e+05_r8 ! chs : upper canopy stems per unit stem area
    !      chu= 8.436e+03_r8 ! chu : upper canopy leaves per unit leaf area
    !      chl= 8.436e+03_r8 ! chl : lower canopy leaves & stems per unit leaf/stem area
    !----------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! intercepted water capacity (mm h2o per unit leaf area == kg m-2)
    !-----------------------------------------------------------------------
    !           ! wliqmin  : minimum per unit vegetated area
    !      wliqumax = 0.20_r8  ! wliqumax : maximum per unit upper canopy leaf area 
    !      wliqsmax = 0.40_r8  ! wliqsmax : maximum per unit upper canopy stem area
    !      wliqlmax = 0.20_r8  ! wliqlmax : maximum per unit lower canopy stem & leaf area
    !=======================================================================

    !-----------------------------------------------------------------------
    ! intercepted snow capacity (mm h2o per unit leaf area == kg m-2)
    !-----------------------------------------------------------------------
    !           ! wsnomin  : minimum per unit vegetated area 
    !      wsnoumax = 2.00_r8  ! wsnoumax : maximum per unit upper canopy leaf area 
    !      wsnosmax = 4.00_r8  ! wsnosmax : maximum per unit upper canopy stem area
    !      wsnolmax = 2.00_r8  ! wsnolmax : maximum per unit lower canopy stem & leaf area
    !=======================================================================

    !------------------------------------------------------------
    ! decay time for intercepted liquid dripoff (sec)
    !------------------------------------------------------------
    !      tdripu = 7200.0_r8  ! tdripu : upper canopy leaves (2 hours)
    !      tdrips = 7200.0_r8  ! tdrips : upper canopy stems (2 hours)
    !      tdripl = 7200.0_r8  ! tdripl : lower canopy leaves & stems (2 hours)
    !============================================================

    !------------------------------------------------------------
    ! decay time for snow blowoff (sec)
    !------------------------------------------------------------
    !      tblowu = 43200.0_r8 ! tblowu : upper canopy leaves (12 hours)
    !      tblows = 43200.0_r8 ! tblows : upper canopy stems (12 hours)
    !      tblowl = 43200.0_r8 ! tblowl : lower canopy leaves & stems (12 hours)
    !============================================================
    !------------------------------------------------------------------------
    ! PFTs (top to bottom)
    !------------------------------------------------------------------------
    !  1: tropical broadleaf evergreen trees
    !  2: tropical broadleaf drought-deciduous trees
    !  3: warm-temperate broadleaf evergreen trees
    !  4: temperate conifer evergreen trees
    !  5: temperate broadleaf cold-deciduous trees
    !  6: boreal conifer evergreen trees
    !  7: boreal broadleaf cold-deciduous trees
    !  8: boreal conifer cold-deciduous trees
    !  9: evergreen shrubs
    ! 10: cold-deciduous shrubs
    ! 11: warm (c4) grasses
    ! 12: cool (c3) grasses

    !--------------------------------------------------------------------------
    ! PFT climatic constraint definitions (left to right)
    !--------------------------------------------------------------------------
    ! TminL  : absolute minimum temperature (lower limit, C) 
    ! TminU  : absolute minimum temperature (upper limit, C) 
    ! Twarm  : temperature of the warmest month (mean??, C) [C4 only]
    ! GDD    : min growing degree days above 5 C threshold [upper canopy], or
    !          min growing degree days above 0 C threshold [lower canopy]

    ! DTP 2001/06/07: Changed this after studying code in climate.f. 
    !      Values of 9999 indicate this constraint is not used to 
    !      determine existence of the PFT.  
    !      dummyvarpk(1:48)=(/ &
    !------------------------------------------------------------------------
    ! TminL    TminU    Twarm    GDD    PFT
    !------------------------------------------------------------------------
    !                 0.0_r8,  9999.0_r8,   9999.0_r8,   9999.0_r8,&  !   1
    !                0.0_r8,  9999.0_r8,   9999.0_r8,   9999.0_r8,&  !   2
    !        -10.0_r8,     0.0_r8,   9999.0_r8,   9999.0_r8,&  !   3
    !        -45.0_r8,     0.0_r8,   9999.0_r8,   1200.0_r8,&  !   4
    !        -45.0_r8,     0.0_r8,   9999.0_r8,   1200.0_r8,&  !   5
    !        -57.5_r8,   -45.0_r8,   9999.0_r8,    350.0_r8,&  !   6
    !        -57.5_r8,   -45.0_r8,   9999.0_r8,    350.0_r8,&  !   7
    !       9999.0_r8,   -45.0_r8,   9999.0_r8,    350.0_r8,&  !   8
    !       9999.0_r8,  9999.0_r8,   9999.0_r8,    100.0_r8,&  !   9
    !       9999.0_r8,  9999.0_r8,   9999.0_r8,    100.0_r8,&  !  10
    !       9999.0_r8,  9999.0_r8,     22.0_r8,    100.0_r8,&  !  11
    !       9999.0_r8,  9999.0_r8,   9999.0_r8,    100.0_r8/)   !  12
    !!========================================================================
    !      i=0
    !      DO j = 1, npft 
    !         TminL(j) =  dummyvarpk(j+0+i)
    !         TminU(j) =  dummyvarpk(j+1+i)
    !         Twarm(j) =  dummyvarpk(j+2+i)
    !         GDD(j)   =  dummyvarpk(j+3+i)
    !         i=i+3
    !      END DO

    !========================================================================

    nveg = 15  ! Number of allowed Vegetation Types in this parameter set   

    !------------------------------------------------------------------------
    ! Vegetation type classifications (used in subroutine iniveg)
    !------------------------------------------------------------------------
    !  1: tropical evergreen forest / woodland
    !  2: tropical deciduous forest / woodland
    !  3: temperate evergreen broadleaf forest / woodland
    !  4: temperate evergreen conifer forest / woodland
    !  5: temperate deciduous forest / woodland
    !  6: boreal evergreen forest / woodland
    !  7: boreal deciduous forest / woodland
    !  8: mixed forest / woodland
    !  9: savanna
    ! 10: grassland / steppe
    ! 11: dense shrubland
    ! 12: open shrubland
    ! 13: tundra
    ! 14: desert
    ! 15: polar desert / rock / ice

    !------------------------------------------------------------------------
    ! Variable name definitions (left to right)
    !------------------------------------------------------------------------
    ! plaievgr   : initial total LAI of evergreen tree (upper canopy) PFTs
    ! plaideci   : initial total LAI of deciduous tree (upper canopy) PFTs
    ! plaishrub  : initial total LAI of shrub (lower canopy) PFTs
    ! plaigrass  : initial total LAI of grass (lower canopy) PFTs
    ! NOTE: These data are read in to array plai_init. 
    !       The original variable names could be dropped. 
    !      dummyvarpk(1:60)=(/ &
    !------------------------------------------------------------------------
    !      plaievgr     plaideci    plaishrub    plaigrass      Veg Type
    !------------------------------------------------------------------------
    !       5.00_r8,     1.00_r8,     0.25_r8,     0.25_r8,&        !  1
    !       1.00_r8,     5.00_r8,     0.25_r8,     0.25_r8,&        !  2
    !       4.00_r8,     1.00_r8,     0.25_r8,     0.25_r8,&        !  3
    !       3.00_r8,     1.00_r8,     0.25_r8,     0.25_r8,&        !  4
    !       1.00_r8,     3.00_r8,     0.25_r8,     0.25_r8,&        !  5
    !       3.00_r8,     1.00_r8,     0.25_r8,     0.25_r8,&        !  6
    !       1.00_r8,     3.00_r8,     0.25_r8,     0.25_r8,&        !  7
    !       2.00_r8,     2.00_r8,     0.25_r8,     0.25_r8,&        !  8
    !       0.50_r8,     1.00_r8,     0.50_r8,     2.00_r8,&        !  9
    !       0.25_r8,     0.25_r8,     0.50_r8,     2.50_r8,&        ! 10
    !       0.10_r8,     0.10_r8,     1.00_r8,     0.50_r8,&        ! 11
    !       0.00_r8,     0.00_r8,     0.25_r8,     0.25_r8,&        ! 12
    !       0.00_r8,     0.00_r8,     1.00_r8,     1.00_r8,&        ! 13
    !       0.00_r8,     0.00_r8,     0.05_r8,     0.05_r8,&        ! 14
    !       0.00_r8,     0.00_r8,     0.05_r8,     0.05_r8/)        ! 15
    !========================================================================
    !      i=0
    !      DO j = 1, nveg 
    !          plai_init(1,j)=  dummyvarpk(j+0+i)
    !         plai_init(2,j)=  dummyvarpk(j+1+i)
    !          plai_init(3,j)=  dummyvarpk(j+2+i)
    !          plai_init(4,j)=  dummyvarpk(j+3+i)
    !        i=i+3
    !      END DO

    plai_init(1:4,1:15) =RESHAPE ((/&
         5.00_r8,     1.00_r8,        0.25_r8,     0.25_r8,&               !  1
         1.00_r8,     5.00_r8,        0.25_r8,     0.25_r8,&               !  2
         4.00_r8,     1.00_r8,        0.25_r8,     0.25_r8,&               !  3
         3.00_r8,     1.00_r8,        0.25_r8,     0.25_r8,&               !  4
         1.00_r8,     3.00_r8,        0.25_r8,     0.25_r8,&               !  5
         3.00_r8,     1.00_r8,        0.25_r8,     0.25_r8,&               !  6
         1.00_r8,     3.00_r8,        0.25_r8,     0.25_r8,&               !  7
         2.00_r8,     2.00_r8,        0.25_r8,     0.25_r8,&               !  8
         0.50_r8,     1.00_r8,        0.50_r8,     2.00_r8,&               !  9
         0.25_r8,     0.25_r8,        0.50_r8,     2.50_r8,&               ! 10
         0.10_r8,     0.10_r8,        1.00_r8,     0.50_r8,&               ! 11
         0.00_r8,     0.00_r8,        0.25_r8,     0.25_r8,&               ! 12
         0.00_r8,     0.00_r8,        1.00_r8,     1.00_r8,&               ! 13
         0.00_r8,     0.00_r8,        0.05_r8,     0.05_r8,&               ! 14
         0.00_r8,     0.00_r8,        0.05_r8,     0.05_r8/), (/4, nveg/) )! 15
    !
    !

    !------------------------------------------------------------------------
    ! Other miscellaneous variables needed for initializing plant LAI.
    !------------------------------------------------------------------------
    ! plaiupper    : Potental LAI of upper canopy (uniform initial vegetation) 
    ! plailower    : Potental LAI of lower canopy (uniform initial vegetation) 
    ! xminlai      : Minimum LAI for each existing PFT
    ! sapfrac_init : Initial value of sapwood fraction used for all woody PFTs
    !------------------------------------------------------------------------
    !      plaiupper    = 0.5_r8  ! plaiupper
    !      plailower    = 0.5_r8  ! plailower
    !      xminlai      = 0.01_r8 ! xminlai 
    !      sapfrac_init = 0.1_r8  ! sapfrac_init
    !========================================================================

    ! ************************************************************************
    ! define rooting profiles
    ! ************************************************************************
    !
    ! define rooting profiles based upon data published in:
    !
    ! Jackson et al., 1996:  A global analysis of root distributions
    ! for terrestrial biomes, Oecologia, 108, 389-411.
    !
    ! and
    !
    ! Jackson et al., 1997:  A global budget for fine root biomass, 
    ! surface area, and nutrient contents, Proceedings of the National
    ! Academy of Sciences, 94, 7362-7366.
    !
    ! rooting profiles are defined by the "beta" parameter
    !
    ! beta1 is assigned to the lower vegetation layer (grasses and shrubs)
    ! beta2 is assigned to the upper vegetation layer (trees)
    !
    ! according to Jackson et al. (1996, 1997), the values of beta
    ! typically fall in the following range
    !
    ! note that the 1997 paper specifically discusses the distribution
    ! of *fine roots* (instead of total root biomass), which may be more
    ! important for water and nutrient uptake
    !
    ! --------------                 ------------   ------------
    ! forest systems                 beta2 (1996)   beta2 (1997)
    ! --------------                 ------------   ------------
    ! tropical evergreen forest:        0.962          0.972
    ! tropical deciduous forest:        0.961          0.982
    ! temperate conifer forest:         0.976          0.980
    ! temperate broadleaf forest:       0.966          0.967
    ! all tropical/temperate forest:    0.970  
    ! boreal forest:                    0.943          0.943
    ! all trees:                                       0.976
    !
    ! -------------------------      ------------   ------------
    ! grassland / shrub systems      beta1 (1996)   beta1 (1997)
    ! -------------------------      ------------   ------------
    ! tropical grassland / savanna:     0.972          0.972
    ! temperate grassland:              0.943          0.943
    ! all grasses:                      0.952          0.952
    ! schlerophyllous shrubs:           0.964          0.950
    ! all shrubs:                       0.978          0.975
    ! crops:                            0.961
    ! desert:                           0.975          0.970
    ! tundra:                           0.914
    !
    ! --------------                 ------------
    ! all ecosystems                 beta  (1996)
    ! --------------                 ------------
    ! all ecosystems:                   0.966
    !
    ! for global simulations, we typically assign the following
    ! values to the beta parameters
    !
    ! beta1 = 0.950, which is typical for tropical/temperate grasslands
    ! beta2 = 0.970, which is typical for tropical/temperate forests
    !
    ! however, these values could be (and should be) further refined
    ! when using the model for specific regions
    ! 
    !      beta1 = 0.950_r8  ! beta1: for lower layer herbaceous plants
    !      beta2 = 0.975_r8  ! beta2: for upper layer trees      
    ! ******************************************************************************
    !
    !========================================================================
    ! params.soi : soil parameters....
    !========================================================================
    !
    !   6    ! nsoilay : number of soil layers (actually a constant in comsoi.h)
    !        !           could be read in here, if a new constant, maxsoilay,
    !        !           was introduced to dimension the soil layer arrays  
    !      dummyvarpk(1:6)=(/  &
    !------------------------------------------------------
    ! Soil layer thicknesses (m) 
    ! N.B. Number of layers must equal nsoilay!!!
    !------------------------------------------------------
    !         0.10_r8,&  ! hsoi(1)
    !         0.15_r8,&  ! hsoi(2)
    !         0.25_r8,&  ! hsoi(3)
    !         0.50_r8,&  ! hsoi(4)
    !         1.00_r8,&  ! hsoi(5)
    !         2.00_r8/)  ! hsoi(6)
    !======================================================
    !      DO j = 1, nsoilay
    !        hsoi(j) = dummyvarpk(j)
    !      END DO

    !------------------------------------------------------
    ! Other miscellaneous soil parameters
    !------------------------------------------------------

    !      bperm    = 0.10_r8  ! bperm   : soil hydraulic diffusivity lower b.c. (units???)
    !      wpudmax  = 4.5_r8   ! wpudmax : normalization constant for puddles (kg m-2) 
    !      zwpmax   = 0.5_r8   ! zwpmax  : maximum value of zwpud (currently assumed to be 0.5)

    nsoil2   = 11    ! nsoi         : number of soil texture classes

    !======================================================

    IF (nsoil2 .NE. ndat) THEN ! Is this needed????
       WRITE (nfprt, 9031) parm_file, nsoil2
       WRITE (nfprt, 9032) ndat
       GO TO 9006 ! In the circumstances this seems the best thing to do! 
    END IF


    !------------------------------------------------------
    ! Soil properties data from Rawls et al. (1992)
    ! Organic properties data compiled by Mustapha El Maayar (2000)
    ! Organic FC and WP taken from Nijssen et al., 1997 (JGR; table 3 OBS-top)

    !------------------------------------------------------
    ! Variable column header definitions
    !------------------------------------------------------
    ! Sand     : sand fraction
    ! Silt     : silt fraction
    ! Clay     : clay fraction
    ! Porosity : porosity (volume fraction)
    ! FC       : field capacity (volume fraction)
    ! WP       : wilting point (volume fraction)
    ! bexp     : Campbell's 'b' exponent
    ! AEP      : air entry potential (m-H20)
    ! SHC      : saturated hydraulic conductivity (m s-1)
    !      ! dummyvarpk=0.0_r8
    !      dummyvarpk(1:108)=(/ &
    !------------------------------------------------------------------------------------------------------
    !  Sand   Silt   Clay  Porosity    FC      WP    bexp    AEP       SHC        Texture class
    !------------------------------------------------------------------------------------------------------
    !      0.92_r8,0.05_r8,0.03_r8,0.437_r8,0.091_r8,0.033_r8,1.7_r8,0.07_r8,5.83300E-05_r8,&  ! Sand
    !      0.81_r8,0.12_r8,0.07_r8,0.437_r8,0.125_r8,0.055_r8,2.1_r8,0.09_r8,1.69720E-05_r8,&  ! Loamy Sand
    !      0.65_r8,0.25_r8,0.10_r8,0.453_r8,0.207_r8,0.095_r8,3.1_r8,0.15_r8,7.19440E-06_r8,&  ! Sandy Loam
    !      0.42_r8,0.40_r8,0.18_r8,0.463_r8,0.270_r8,0.117_r8,4.5_r8,0.11_r8,3.66670E-06_r8,&  ! Loam
    !      0.20_r8,0.65_r8,0.15_r8,0.501_r8,0.330_r8,0.133_r8,4.7_r8,0.21_r8,1.88890E-06_r8,&  ! Silty Loam
    !      0.60_r8,0.13_r8,0.27_r8,0.398_r8,0.255_r8,0.148_r8,4.0_r8,0.28_r8,1.19440E-06_r8,&  ! Sandy Clay Loam
    !      0.32_r8,0.34_r8,0.34_r8,0.464_r8,0.318_r8,0.197_r8,5.2_r8,0.26_r8,6.38890E-07_r8,&  ! Clay Loam
    !      0.09_r8,0.58_r8,0.33_r8,0.471_r8,0.366_r8,0.208_r8,6.6_r8,0.33_r8,4.16670E-07_r8,&  ! Silty Clay Loam
    !      0.53_r8,0.07_r8,0.40_r8,0.430_r8,0.339_r8,0.239_r8,6.0_r8,0.29_r8,3.33330E-07_r8,&  ! Sandy Clay
    !      0.10_r8,0.45_r8,0.45_r8,0.479_r8,0.387_r8,0.250_r8,7.9_r8,0.34_r8,2.50000E-07_r8,&  ! Silty Clay
    !      0.20_r8,0.20_r8,0.60_r8,0.475_r8,0.396_r8,0.272_r8,7.6_r8,0.37_r8,1.66670E-07_r8,&  ! Clay
    !      0.00_r8,0.00_r8,0.00_r8,0.800_r8,0.600_r8,0.200_r8,7.6_r8,0.37_r8,1.00000E-05_r8/)  ! Organic - TBA
    !====================================================================
    !      i=0
    !      DO j = 1, ndat
    !          texdat  (1,j)= dummyvarpk(j+0+i)
    !          texdat  (2,j)= dummyvarpk(j+1+i)
    !          texdat  (3,j)= dummyvarpk(j+2+i)
    !          porosdat  (j)= dummyvarpk(j+3+i)
    !          sfielddat (j)= dummyvarpk(j+4+i)
    !          swiltdat  (j)= dummyvarpk(j+5+i)
    !          bexdat    (j)= dummyvarpk(j+6+i)
    !          suctiondat(j)= dummyvarpk(j+7+i)
    !          hydrauldat(j)= dummyvarpk(j+8+i)
    !          i=i+8
    !      END DO

    !
    ! Rawls et al. (1992) soil properties data
    !
    !      ------------------
    !       sand  silt  clay
    !      ------------------
    texdat =RESHAPE ((/&
         0.920_r8, 0.050_r8, 0.030_r8,&  ! sand
         0.810_r8, 0.120_r8, 0.070_r8,&  ! loamy sand
         0.650_r8, 0.250_r8, 0.100_r8,&  ! sandy loam
         0.420_r8, 0.400_r8, 0.180_r8,&  ! loam
         0.200_r8, 0.650_r8, 0.150_r8,&  ! silt loam
         0.600_r8, 0.130_r8, 0.270_r8,&  ! sandy clay loam
         0.320_r8, 0.340_r8, 0.340_r8,&  ! clay loam
         0.090_r8, 0.580_r8, 0.330_r8,&  ! silty clay loam
         0.530_r8, 0.070_r8, 0.400_r8,&  ! sandy clay
         0.100_r8, 0.450_r8, 0.450_r8,&  ! silty clay
         0.200_r8, 0.200_r8, 0.600_r8/), (/3, ndat/) )   ! clay

    porosdat  (1:ndat) =RESHAPE ((/&! Porosity
         !------------------------------------------------------------------------------------------------------
         !      Porosity        Texture class                                                        SIB2     IBIS
         !------------------------------------------------------------------------------------------------------ 
         !                                                    Soil  Name            % clay  % sand    Poros     Poros
         0.373_r8,&  ! Sand                             1    sand                      3       92    0.373     0.437
         0.386_r8,&  ! Loamy Sand                       2    loamy sand                5       82    0.386     0.437
         0.407_r8,&  ! Sandy Loam                       3    sandy loam               10       65    0.407     0.453
         0.436_r8,&  ! Loam                             6    loam                     18       42    0.436     0.463
         0.461_r8,&  ! Silty Loam                       4    silt loam                13       22    0.461     0.501
         0.416_r8,&  ! Sandy Clay Loam                  7    sandy clay loam          28       58    0.416     0.398
         0.449_r8,&  ! Clay Loam                        9    clay loam                39       32    0.449     0.464
         0.476_r8,&  ! Silty Clay Loam                 10    silty clay loam          39       10    0.476     0.471
         0.423_r8,&  ! Sandy Clay                       8    sandy clay               40       52    0.423     0.430
         0.480_r8,&  ! Silty Clay                      11    silty clay               41        7    0.480     0.479
         0.465_r8,&  ! Clay                            12    clay                     65       19    0.465     0.475
         0.800_r8/), (/ndat/) )     ! Organic - TBA     5    silt                      7        7    0.436     0.800

    sfielddat  (1:ndat) =RESHAPE ((/&! field capacity (volume fraction)
         !------------------------------------------------------------------------------------------------------
         !      FC        Texture class                                                           SIB2       IBIS
         !------------------------------------------------------------------------------------------------------ 
         !                                                 Soil  Name               % clay  % sand        PhiSat       PhiSat
         0.047_r8,&  ! Sand                                1    sand                  3          92        0.047       0.091
         0.064_r8,&  ! Loamy Sand                          2    loamy sand            5          82        0.064       0.125
         0.107_r8,&  ! Sandy Loam                          3    sandy loam           10          65        0.107       0.207
         0.214_r8,&  ! Loam                                6    loam                 18          42        0.214       0.270
         0.391_r8,&  ! Silty Loam                          4    silt loam            13          22        0.391       0.330
         0.132_r8,&  ! Sandy Clay Loam                     7    sandy clay loam      28          58        0.132       0.255
         0.289_r8,&  ! Clay Loam                           9    clay loam            39          32        0.289       0.318
         0.561_r8,&  ! Silty Clay Loam                    10    silty clay loam      39          10        0.561       0.366
         0.158_r8,&  ! Sandy Clay                          8    sandy clay           40          52        0.158       0.339
         0.614_r8,&  ! Silty Clay                         11    silty clay           41           7        0.614       0.387
         0.428_r8,&  ! Clay                               12    clay                 65          19        0.428       0.396
         0.214_r8/), (/ndat/) )   ! Organic - TBA          5    silt                  7           7        0.214       0.600


    swiltdat  (1:ndat) =RESHAPE ((/&! wilting point (volume fraction)
         !------------------------------------------------------------------------------------------------------
         !      WP        Texture class
         !------------------------------------------------------------------------------------------------------ 
         0.033_r8,&  ! Sand
         0.055_r8,&  ! Loamy Sand
         0.095_r8,&  ! Sandy Loam
         0.117_r8,&  ! Loam
         0.133_r8,&  ! Silty Loam
         0.148_r8,&  ! Sandy Clay Loam
         0.197_r8,&  ! Clay Loam
         0.208_r8,&  ! Silty Clay Loam
         0.239_r8,&  ! Sandy Clay
         0.250_r8,&  ! Silty Clay
         0.272_r8,&  ! Clay
         0.200_r8/), (/ndat/) )   ! Organic - TBA

    bexdat  (1:ndat) =RESHAPE ((/&!  Campbell's 'b' exponent
         !------------------------------------------------------------------------------------------------------
         !      bexp        Texture class                                                         SIB2       IBIS
         !------------------------------------------------------------------------------------------------------ 
         !                                              Soil  Name            % clay  % sand         BEE       BEE 
         3.387_r8,&  ! Sand                              1    sand               3         92       3.387     1.7
         3.705_r8,&  ! Loamy Sand                        2    loamy sand         5         82       3.705     2.1
         4.500_r8,&  ! Sandy Loam                        3    sandy loam        10         65       4.500     3.1
         5.772_r8,&  ! Loam                              6    loam              18         42       5.772     4.5
         4.977_r8,&  ! Silty Loam                        4    silt loam         13         22       4.977     4.7
         7.362_r8,&  ! Sandy Clay Loam                   7    sandy clay loam   28         58       7.362     4.0
         9.111_r8,&  ! Clay Loam                    9    clay loam         39         32       9.111     5.2
         9.111_r8,&  ! Silty Clay Loam             10    silty clay loam   39         10       9.111     6.6
         9.270_r8,&  ! Sandy Clay                   8    sandy clay        40         52       9.270     6.0
         9.429_r8,&  ! Silty Clay                  11    silty clay        41          7       9.429     7.9
         13.245_r8,&  ! Clay                       12     clay              65         19      13.245     7.6
         7.6_r8/), (/ndat/) )  ! Organic - TBA      5    silt               7          7       5.772     7.6

    suctiondat  (1:ndat) =RESHAPE ((/&!  air entry potential (m-H20)
         !------------------------------------------------------------------------------------------------------
         !      AEP        Texture class
         !------------------------------------------------------------------------------------------------------ 
         0.07_r8,&  ! Sand
         0.09_r8,&  ! Loamy Sand
         0.15_r8,&  ! Sandy Loam
         0.11_r8,&  ! Loam
         0.21_r8,&  ! Silty Loam
         0.28_r8,&  ! Sandy Clay Loam
         0.26_r8,&  ! Clay Loam
         0.33_r8,&  ! Silty Clay Loam
         0.29_r8,&  ! Sandy Clay
         0.34_r8,&  ! Silty Clay
         0.37_r8,&  ! Clay
         0.37_r8/), (/ndat/) )  ! Organic - TBA

    hydrauldat  (1:ndat) =RESHAPE ((/&!saturated hydraulic conductivity (m s-1)
         !------------------------------------------------------------------------------------------------------
         !      SHC        Texture class           sib2            ibis
         !------------------------------------------------------------------------------------------------------ 
         !                                                                             Soil    Name            % clay  % sand      SatCo  
         0.236E-04_r8,&  ! Sand                0.236E-04     5.83300E-05_r8         1    sand               3       92    0.236E-04 
         0.166E-04_r8,&  ! Loamy Sand          0.166E-04     1.69720E-05_r8         2    loamy sand         5       82    0.166E-04 
         0.910E-05_r8,&  ! Sandy Loam          0.910E-05     7.19440E-06_r8         3    sandy loam         10      65    0.910E-05 
         0.405E-05_r8,&  ! Loam                0.405E-05     3.66670E-06_r8         6    loam               18      42    0.405E-05 
         0.200E-05_r8,&  ! Silty Loam          0.200E-05     1.88890E-06_r8         4    silt loam          13      22    0.200E-05 
         0.711E-05_r8,&  ! Sandy Clay Loam     0.711E-05     1.19440E-06_r8         7    sandy clay loam    28      58    0.711E-05 
         0.285E-05_r8,&  ! Clay Loam           0.285E-05     6.38890E-07_r8         9    clay loam          39      32    0.285E-05 
         0.131E-05_r8,&  ! Silty Clay Loam     0.131E-05     4.16670E-07_r8        10    silty clay loam    39      10    0.131E-05 
         0.576E-05_r8,&  ! Sandy Clay          0.576E-05     3.33330E-07_r8         8    sandy clay         40      52    0.576E-05 
         0.118E-05_r8,&  ! Silty Clay          0.118E-05     2.50000E-07_r8        11    silty clay         41      7     0.118E-05 
         0.180E-05_r8,&  ! Clay                0.180E-05     1.66670E-07_r8        12    clay               65      19    0.180E-05 
         1.00000E-05_r8/), (/ndat/) )   ! Organic - TBA 1.00000E-05                 5    silt               7       7     0.405E-05 

    !----------------------------------------------------------------------
    ! Decomposition pool/transformation parameters (see also Kucharik 
    ! et al. 2000)

    !----------------------------------------------------------------------
    ! lig_frac: split of lignified litter material between protected and
    ! non-protected slow OM pools
    !----------------------------------------------------------------------
    !      lig_frac =  0.50_r8 

    !----------------------------------------------------------------------
    ! fbsom: protected biomass as a fraction of total soil organic carbon
    ! from Verberne et al., 1990
    !----------------------------------------------------------------------
    !      fbsom   = 0.017_r8

    !----------------------------------------------------------------------
    ! effac: efficiency of microbial biomass reincorporated into biomass
    ! pool (from NCSOIL parameterizations; Molina et al., 1983)
    !----------------------------------------------------------------------
    !      effac =  0.40_r8 

    !======================================================================
    ! Define C:N ratios of substrate pools and biomass: 
    ! Values for metabolic and structural plant material and for lignin are 
    ! from Parton et al., 1987 and Whitmore and Parry, 1988, indexed as follows:
    !      cnr(1): c:n ratio of microbial biomass
    !      cnr(2): c:n ratio of passive soil carbon
    !      cnr(3): c:n ratio of protected slow soil carbon
    !      cnr(4): c:n ratio of non-protected slow soil C
    !      cnr(5): c:n ratio of resistant litter lignin
    !      cnr(6): c:n ratio of structural plant (leaf and root) litter
    !      cnr(7): c:n ratio of metabolic (plant and root) litter
    !      cnr(8): c:n ratio of woody biomass components
    !      dummyvarpk(1:8)=(/ &
    !---------------------------------------------------------------------
    !   cnr(1)  cnr(2)  cnr(3)  cnr(4)  cnr(5)  cnr(6)  cnr(7)  cnr(8)
    !---------------------------------------------------------------------
    !      8.0_r8,   15.0_r8,   10.0_r8,   15.0_r8,  100.0_r8,  150.0_r8,    6.0_r8,   250.0_r8/)
    !=====================================================================
    !      DO j=1,8
    !         cnr(j) = dummyvarpk(j)
    !      END DO

    ! Miscellaneous other C:N factors...
    !      fmax  : maximum fraction allowed in resistant fraction
    !      rconst: rconst is a constant defined as 1200 [Huh?]
    !      cnleaf: average c:n ratio for leaf litterfall 
    !      cnroot: average c:n ratio for root turnover
    !      cnwood: average c:n ratio for woody debris
    !      dummyvarpk(1:5)=(/ &
    !---------------------------------------------------------------------
    !      fmax    rconst    cnleaf    cnroot   cnwood
    !---------------------------------------------------------------------
    !       0.45_r8,   1200.0_r8,     40.0_r8,     60.0_r8,    200.0_r8/)
    !=====================================================================
    !      fmax   = dummyvarpk(1)
    !      rconst = dummyvarpk(2)
    !      cnleaf = dummyvarpk(3)
    !      cnroot = dummyvarpk(4)
    !      cnwood = dummyvarpk(5)


    !--------------------------------------------------------------------------------------

    ! Specific maximum decay rate or growth constants; rates are per day.
    ! Constants are taken from Parton et al., 1987 and Verberne et al., 1990
    ! and special issue of Geoderma (comparison of 9 organic matter models) in Dec. 1997

    ! Leaching parameterization was changed to agree with field data, and led to 
    ! changes in the values of the constants given below.  

    ! Approximate factors for Verberne et al. model where efficiencies are 100%
    ! for some of the transformations: one problem was that their rate constants were
    ! based on 25C, and our modifying functions are based on 15 C...thus the rate constants
    ! are somewhat smaller compared to the Verberne et al. (1990) model parameters.
    ! Rates are based on a daily decomposition timestep (per day)
    !--------------------------------------------------------------------------------------

    !      klm: dpm leaf litter --> microbial biomass
    !      kls: spm leaf litter --> microbial biomass
    !      kll: rpm leaf litter --> non or protected om
    !      krm: dpm root litter --> microbial biomass
    !      krs: spm root litter --> microbial biomass
    !      krl: rpm root litter --> non or protected om 
    !      kwm: dpm woody litter --> microbial biomass
    !      kws: spm woody litter --> microbial biomass
    !      kwl: rpm woody litter --> non or protected om 
    !      dummyvarpk(1:9)=(/ &
    !---------------------------------------------------------------------
    !     klm    kls    kll    krm    krs    krl    kwm    kws    kwl
    !---------------------------------------------------------------------
    !      0.150_r8, 0.010_r8,0.010_r8,0.100_r8,0.005_r8,0.005_r8,0.001_r8,0.001_r8,0.001_r8/)
    !=====================================================================
    !
    !      klm  = dummyvarpk(1)
    !      kls  = dummyvarpk(2)
    !      kll  = dummyvarpk(3)
    !      krm  = dummyvarpk(4)
    !      krs  = dummyvarpk(5)
    !      krl  = dummyvarpk(6)
    !      kwm  = dummyvarpk(7)
    !      kws  = dummyvarpk(8)
    !      kwl  = dummyvarpk(9)
    !----------------------------------------------------------------------
    !
    !      kbn: biomass --> non protected organic matter 
    !      kbp: biomass --> protected organic matter
    !      knb: non-protected om --> biomass
    !      kns: non-protected om --> stabilized om
    !      kpb: protected om --> biomass
    !      kps: protected om --> stabilized om
    !      ksb: stabilized om --> biomass
    !       dummyvarpk(1:7)=(/ &
    !---------------------------------------------------------------------
    !      kbn     kbp     knb     kns     kpb     kps     ksb
    !---------------------------------------------------------------------
    !      0.045_r8,   0.005_r8,   0.001_r8, 1.0e-06_r8,  0.0001_r8, 1.0e-06_r8, 8.0e-07_r8/)
    !=====================================================================
    !       kbn  = dummyvarpk(1)
    !       kbp  = dummyvarpk(2)
    !       knb  = dummyvarpk(3)
    !       kns  = dummyvarpk(4)
    !       kpb  = dummyvarpk(5)
    !       kps  = dummyvarpk(6)
    !       ksb  = dummyvarpk(7)
    !
    ! Yields (efficiencies) with which microbes gain biomass from C 
    ! source; the rest is driven off as CO2 (microbial respiration). All 
    ! microbial CO2 is assumed to leave the soil profile over the course 
    ! of a year. Values are taken primarily from the models of Verberne 
    ! and from CENTURY.
    !----------------------------------------------------------------------
    !      ylm: efficiency for metabolic plant material - leaf matter
    !      yrm: efficiency for metabolic plant material - root matter
    !      ywm: efficiency for metabolic plant material - woody matter
    !      yls: efficiency for structural plant material - leaf matter
    !      yrs: efficiency for structural plant material - root matter
    !      yws: efficiency for structural plant material - woody matter
    !      yll: plant material resistant fraction - leaf matter
    !      yrl: plant material resistant fraction - root matter
    !      ywl: plant material resistant fraction - woody matter
    !      dummyvarpk(1:9)=(/ &
    !---------------------------------------------------------------------
    !   ylm    yrm    ywm    yls    yrs    yws    yll    yrl    ywl
    !---------------------------------------------------------------------
    !      0.40_r8,0.40_r8,0.40_r8,0.30_r8,0.30_r8,0.30_r8,1.00_r8,1.00_r8,1.00_r8/)
    !=====================================================================       
    !      ylm = dummyvarpk(1)
    !      yrm = dummyvarpk(2)
    !      ywm = dummyvarpk(3)
    !      yls = dummyvarpk(4)
    !      yrs = dummyvarpk(5)
    !      yws = dummyvarpk(6)
    !      yll = dummyvarpk(7)
    !      yrl = dummyvarpk(8)
    !      ywl = dummyvarpk(9)


    !      ybn: biomass       --> non-protected pool
    !      ybp: biomass       --> protected pool
    !      yps: protected     --> passive
    !      yns: non-protected --> passive
    !      ysb: passive pool  --> biomass
    !      ypb: protected     --> biomass
    !      ynb: non-protected --> biomass
    !      dummyvarpk(1:7)=(/ &
    !---------------------------------------------------------------------
    !      ybn     ybp     yps     yns     ysb     ypb     ynb
    !---------------------------------------------------------------------
    !      1.00_r8,   1.00_r8,   1.00_r8,   1.00_r8,   0.20_r8,   0.20_r8,    0.25_r8/)
    !=====================================================================
    !
    !      ybn  = dummyvarpk(1)
    !      ybp  = dummyvarpk(2)
    !      yps  = dummyvarpk(3)
    !      yns  = dummyvarpk(4)
    !      ysb  = dummyvarpk(5)
    !      ypb  = dummyvarpk(6)
    !      ynb  = dummyvarpk(7)
    !
9031 FORMAT (' RD_PARAM Warning: Number of soil types in ', &
         A10, ' is: ', I2)
9032 FORMAT (' Number of soil types in comtex.h is: ', I2)  

    ! ******************************************************************************
    !open the parameter file 'params.hyd' for input; read in vegetation PFT parameters...
    !
    !      parm_file = 'params.hyd' ! What is this file supposed to contain???
    !
    !      open(UNIT=parm_unit, FILE=parm_file, STATUS='OLD', ERR=9001)
    !      close (parm_unit) 
    !
    !******************************************************************************

    CALL MsgOne('**(RD_PARAM)**','All data read in from parameter files successfully.')

    RETURN


    WRITE (nfprt,2003) parm_file, parm_unit
2003 FORMAT ('RD_PARAM: Unexpected EOF encountered in file ', &
         A12, ' on unit ', I2)
    STOP

    WRITE (nfprt,2001) parm_file, parm_unit 
2001 FORMAT ('RD_PARAM: Error opening parameter file ', A12, &
         ' on unit ', I2)      
    STOP

    WRITE (nfprt,2002) parm_file, parm_unit 
2002 FORMAT ('RD_PARAM: Error reading parameter file ', A12, &
         ' on unit ', I2)      
    STOP

9006 WRITE (nfprt,2006) parm_file, parm_unit
2006 FORMAT ('RD_PARAM: Data inconsistency in parameter file ', &
         A12, ' on unit ', I2) 
    STOP
  END  SUBROUTINE RD_PARAM


  SUBROUTINE ReStartIBIS (nfsibo)

    !INTEGER           ,INTENT(IN   ) :: jbMax
    !INTEGER           ,INTENT(IN   ) :: ifday
    !REAL(KIND=r8)     ,INTENT(IN   ) :: tod
    !INTEGER           ,INTENT(IN   ) :: idate(:)
    !INTEGER           ,INTENT(IN   ) :: idatec(:)
    INTEGER           ,INTENT(IN   ) :: nfsibo
    !INTEGER           ,INTENT(IN   ) :: ibMaxPerJB(:)
    !INTEGER                         :: i



    CALL MsgOne('**(ReStartIBIS)**','Saving physics state for restart')

    WRITE(UNIT=nfsibo)asurd   
    WRITE(UNIT=nfsibo)asuri   
    WRITE(UNIT=nfsibo)totcondub
    WRITE(UNIT=nfsibo)totconduc
    WRITE(UNIT=nfsibo)totcondls
    WRITE(UNIT=nfsibo)totcondl3
    WRITE(UNIT=nfsibo)totcondl4
    WRITE(UNIT=nfsibo)ginvap  
    WRITE(UNIT=nfsibo)gsuvap  
    WRITE(UNIT=nfsibo)gtrans  
    WRITE(UNIT=nfsibo)grunof  
    WRITE(UNIT=nfsibo)gdrain  
    WRITE(UNIT=nfsibo)totlit  
    WRITE(UNIT=nfsibo)totnlit 
    WRITE(UNIT=nfsibo)totfall 
    WRITE(UNIT=nfsibo)totalit 
    WRITE(UNIT=nfsibo)totrlit 
    WRITE(UNIT=nfsibo)totanlit
    WRITE(UNIT=nfsibo)totrnlit
    WRITE(UNIT=nfsibo)totcsoi 
    WRITE(UNIT=nfsibo)totnmic 
    WRITE(UNIT=nfsibo)tco2mic 
    WRITE(UNIT=nfsibo)totnsoi
    WRITE(UNIT=nfsibo)tnpptot 
    WRITE(UNIT=nfsibo)tneetot 
    WRITE(UNIT=nfsibo)tnmin   
    WRITE(UNIT=nfsibo)cdisturb
    WRITE(UNIT=nfsibo)su      
    WRITE(UNIT=nfsibo)ss      
    WRITE(UNIT=nfsibo)sl      
    WRITE(UNIT=nfsibo)stresstu
    WRITE(UNIT=nfsibo)stresstl
    WRITE(UNIT=nfsibo)stressl 
    WRITE(UNIT=nfsibo)stressu 
    WRITE(UNIT=nfsibo)tcthis
    WRITE(UNIT=nfsibo)twthis

    !
    !  daily average variables 
    !
    WRITE(UNIT=nfsibo)adrain          
    WRITE(UNIT=nfsibo)adsnow          
    WRITE(UNIT=nfsibo)adaet          
    WRITE(UNIT=nfsibo)adtrunoff 
    WRITE(UNIT=nfsibo)adsrunoff 
    WRITE(UNIT=nfsibo)addrainage
    WRITE(UNIT=nfsibo)adrh          
    WRITE(UNIT=nfsibo)adsnod          
    WRITE(UNIT=nfsibo)adsnof          
    WRITE(UNIT=nfsibo)adwsoi          
    WRITE(UNIT=nfsibo)adtsoi          
    WRITE(UNIT=nfsibo)adwisoi    
    WRITE(UNIT=nfsibo)adtlaysoi 
    WRITE(UNIT=nfsibo)adwlaysoi 
    WRITE(UNIT=nfsibo)adwsoic    
    WRITE(UNIT=nfsibo)adtsoic    
    WRITE(UNIT=nfsibo)adco2mic  
    WRITE(UNIT=nfsibo)adco2root 
    WRITE(UNIT=nfsibo)adco2soi  
    WRITE(UNIT=nfsibo)adco2ratio
    WRITE(UNIT=nfsibo)adnmintot 
    WRITE(UNIT=nfsibo)decompl
    WRITE(UNIT=nfsibo)decomps
    WRITE(UNIT=nfsibo)gdd0this
    WRITE(UNIT=nfsibo)gdd5this
    WRITE(UNIT=nfsibo)storedn
    WRITE(UNIT=nfsibo)yrleach
    WRITE(UNIT=nfsibo)ynleach
    !
    !  monthly average variables 
    !
    WRITE(UNIT=nfsibo)amrain    
    WRITE(UNIT=nfsibo)amsnow    
    WRITE(UNIT=nfsibo)amaet          
    WRITE(UNIT=nfsibo)amtrunoff 
    WRITE(UNIT=nfsibo)amsrunoff 
    WRITE(UNIT=nfsibo)amdrainage
    WRITE(UNIT=nfsibo)amtemp    
    WRITE(UNIT=nfsibo)amqa          
    WRITE(UNIT=nfsibo)amsolar   
    WRITE(UNIT=nfsibo)amirup    
    WRITE(UNIT=nfsibo)amirdown  
    WRITE(UNIT=nfsibo)amsens    
    WRITE(UNIT=nfsibo)amlatent  
    WRITE(UNIT=nfsibo)amlaiu    
    WRITE(UNIT=nfsibo)amlail    
    WRITE(UNIT=nfsibo)amtsoi    
    WRITE(UNIT=nfsibo)amwsoi    
    WRITE(UNIT=nfsibo)amwisoi   
    WRITE(UNIT=nfsibo)amvwc          
    WRITE(UNIT=nfsibo)amawc          
    WRITE(UNIT=nfsibo)amsnod    
    WRITE(UNIT=nfsibo)amsnof    
    WRITE(UNIT=nfsibo)amnpp          
    WRITE(UNIT=nfsibo)amnpptot  
    WRITE(UNIT=nfsibo)amco2mic  
    WRITE(UNIT=nfsibo)amco2root 
    WRITE(UNIT=nfsibo)amco2soi            
    WRITE(UNIT=nfsibo)amco2ratio
    WRITE(UNIT=nfsibo)amneetot  
    WRITE(UNIT=nfsibo)amnmintot 
    WRITE(UNIT=nfsibo)amts2          
    WRITE(UNIT=nfsibo)amtransu  
    WRITE(UNIT=nfsibo)amtransl  
    WRITE(UNIT=nfsibo)amsuvap   
    WRITE(UNIT=nfsibo)aminvap   
    WRITE(UNIT=nfsibo)amalbedo  
    WRITE(UNIT=nfsibo)amtsoil   
    WRITE(UNIT=nfsibo)amwsoil   
    WRITE(UNIT=nfsibo)amwisoil  
    !
    !  annual total variables 
    !
    WRITE(UNIT=nfsibo)aysolar   
    WRITE(UNIT=nfsibo)ayirup    
    WRITE(UNIT=nfsibo)ayirdown  
    WRITE(UNIT=nfsibo)aysens    
    WRITE(UNIT=nfsibo)aylatent  
    WRITE(UNIT=nfsibo)ayprcp    
    WRITE(UNIT=nfsibo)ayaet          
    WRITE(UNIT=nfsibo)aytrans   
    WRITE(UNIT=nfsibo)aytrunoff 
    WRITE(UNIT=nfsibo)aysrunoff 
    WRITE(UNIT=nfsibo)aydrainage
    WRITE(UNIT=nfsibo)aydwtot   
    WRITE(UNIT=nfsibo)aywsoi    
    WRITE(UNIT=nfsibo)aywisoi   
    WRITE(UNIT=nfsibo)aytsoi    
    WRITE(UNIT=nfsibo)ayvwc          
    WRITE(UNIT=nfsibo)ayawc          
    WRITE(UNIT=nfsibo)aystresstu
    WRITE(UNIT=nfsibo)aystresstl
    WRITE(UNIT=nfsibo)aygpp          
    WRITE(UNIT=nfsibo)aygpptot  
    WRITE(UNIT=nfsibo)aynpp          
    WRITE(UNIT=nfsibo)aynpptot  
    WRITE(UNIT=nfsibo)ayco2mic  
    WRITE(UNIT=nfsibo)ayco2root 
    WRITE(UNIT=nfsibo)ayco2soi  
    WRITE(UNIT=nfsibo)ayneetot  
    WRITE(UNIT=nfsibo)ayrootbio 
    WRITE(UNIT=nfsibo)aynmintot 
    WRITE(UNIT=nfsibo)ayalit    
    WRITE(UNIT=nfsibo)ayblit    
    WRITE(UNIT=nfsibo)aycsoi    
    WRITE(UNIT=nfsibo)aycmic    
    WRITE(UNIT=nfsibo)ayanlit   
    WRITE(UNIT=nfsibo)aybnlit   
    WRITE(UNIT=nfsibo)aynsoi    
    WRITE(UNIT=nfsibo)ayalbedo  
    WRITE(UNIT=nfsibo)firefac      
    WRITE(UNIT=nfsibo)wtot 

    WRITE(UNIT=nfsibo)tlsub
    WRITE(UNIT=nfsibo)t12  
    WRITE(UNIT=nfsibo)t34  
    WRITE(UNIT=nfsibo)q12  
    WRITE(UNIT=nfsibo)q34  
    WRITE(UNIT=nfsibo)ciub 
    WRITE(UNIT=nfsibo)ciuc 
    WRITE(UNIT=nfsibo)cils 
    WRITE(UNIT=nfsibo)cil3 
    WRITE(UNIT=nfsibo)cil4 
    WRITE(UNIT=nfsibo)csub 
    WRITE(UNIT=nfsibo)csuc 
    WRITE(UNIT=nfsibo)csls 
    WRITE(UNIT=nfsibo)csl3 
    WRITE(UNIT=nfsibo)csl4 
    WRITE(UNIT=nfsibo)gsub 
    WRITE(UNIT=nfsibo)gsuc 
    WRITE(UNIT=nfsibo)gsls 
    WRITE(UNIT=nfsibo)gsl3 
    WRITE(UNIT=nfsibo)gsl4 
    WRITE(UNIT=nfsibo)wliqu
    WRITE(UNIT=nfsibo)wliqs
    WRITE(UNIT=nfsibo)wliql
    WRITE(UNIT=nfsibo)wsnou
    WRITE(UNIT=nfsibo)wsnos
    WRITE(UNIT=nfsibo)wsnol

    WRITE(UNIT=nfsibo)fi,fu,fl,tu,ts,tl,tg,ti
    !
    ! nsnolay variables: tsno and hsno
    !

    !
    WRITE(UNIT=nfsibo)tsno

    WRITE(UNIT=nfsibo)hsno

    !
    ! nsoilay variables: tsoi, wisoi, wsoi
    !

    WRITE(UNIT=nfsibo)tsoim,tsoi

    WRITE(UNIT=nfsibo)wisoi

    WRITE(UNIT=nfsibo)wsoim,wsoi

    !
    ! npft variables
    !

    WRITE(UNIT=nfsibo)cbiol

    WRITE(UNIT=nfsibo)cbiow

    WRITE(UNIT=nfsibo)cbior
    WRITE(UNIT=nfsibo)ndtimes
    WRITE(UNIT=nfsibo)nmtimes
    WRITE(UNIT=nfsibo)nytimes
    !
    ! single level variables
    !

    WRITE(UNIT=nfsibo)sapfrac

    WRITE(UNIT=nfsibo)clitlm

    WRITE(UNIT=nfsibo)clitls

    WRITE(UNIT=nfsibo)clitll

    WRITE(UNIT=nfsibo)clitrm

    WRITE(UNIT=nfsibo)clitrs

    WRITE(UNIT=nfsibo)clitrl

    WRITE(UNIT=nfsibo)clitwm

    WRITE(UNIT=nfsibo)clitws

    WRITE(UNIT=nfsibo)clitwl

    WRITE(UNIT=nfsibo)falll

    WRITE(UNIT=nfsibo)fallr

    WRITE(UNIT=nfsibo)fallw

    WRITE(UNIT=nfsibo)totcmic

    WRITE(UNIT=nfsibo)csoislop

    WRITE(UNIT=nfsibo)csoislon

    WRITE(UNIT=nfsibo)csoipas

    WRITE(UNIT=nfsibo)gdd0

    WRITE(UNIT=nfsibo)gdd5

    WRITE(UNIT=nfsibo)tc

    WRITE(UNIT=nfsibo)tw

    WRITE(UNIT=nfsibo)wipud

    WRITE(UNIT=nfsibo)wpud

    WRITE(UNIT=nfsibo)agddu

    WRITE(UNIT=nfsibo)agddl

    WRITE(UNIT=nfsibo)tempu

    WRITE(UNIT=nfsibo)templ

    WRITE(UNIT=nfsibo)a10td

    WRITE(UNIT=nfsibo)a10ancub

    WRITE(UNIT=nfsibo)a10ancls

    WRITE(UNIT=nfsibo)a10ancl4

    WRITE(UNIT=nfsibo)a10ancl3

    WRITE(UNIT=nfsibo)a10scalparamu

    WRITE(UNIT=nfsibo)a10scalparaml

    WRITE(UNIT=nfsibo)a10daylightu

    WRITE(UNIT=nfsibo)a10daylightl

    WRITE(UNIT=nfsibo)dropu

    WRITE(UNIT=nfsibo)dropls

    WRITE(UNIT=nfsibo)dropl4

    WRITE(UNIT=nfsibo)dropl3

    WRITE(UNIT=nfsibo)lai

    WRITE(UNIT=nfsibo)zbot

    WRITE(UNIT=nfsibo)ztop 

    WRITE(UNIT=nfsibo)frac         

    WRITE(UNIT=nfsibo)td  

    WRITE(UNIT=nfsibo) ppli,ppci

    WRITE(UNIT=nfsibo) gl0 ,zorl,gtsea,tseam,qsfc0,tsfc0,qsfcm,tsfcm,tkemyj ,HML,HUML,HVML,TSK,z0sea

    WRITE(UNIT=nfsibo) w0,wm,capac0,capacm,td0,tdm,tcm,tc0,tgm,tg0,z0

    WRITE(UNIT=nfsibo) idateprev

    WRITE(UNIT=nfsibo) sai

    WRITE(UNIT=nfsibo) plai

    WRITE(UNIT=nfsibo) biomass

    WRITE(UNIT=nfsibo) totlaiu

    WRITE(UNIT=nfsibo) totlail

    WRITE(UNIT=nfsibo) totbiou

    WRITE(UNIT=nfsibo) totbiol

    WRITE(UNIT=nfsibo) exist

    WRITE(UNIT=nfsibo) vegtype0

  END SUBROUTINE ReStartIBIS


  SUBROUTINE Finalize_IBIS()
    DEALLOCATE(iMaskIBIS)
    DEALLOCATE(MskAntIBIS)
    DEALLOCATE(nlpoints)
    DEALLOCATE(brf     )
    DEALLOCATE(lonscale)
    DEALLOCATE(latscale)
    DEALLOCATE(xintopo )
    DEALLOCATE(garea   )
    DEALLOCATE(xinveg  )
    DEALLOCATE(deltat  )
    DEALLOCATE(sand    )
    DEALLOCATE(clay    )
    DEALLOCATE(clmwet  )
    DEALLOCATE(clmt    )
    DEALLOCATE(clmtrng )
    DEALLOCATE(clmprec )
    DEALLOCATE(xinwind )
    DEALLOCATE(clmcld  )
    DEALLOCATE(clmq    )
    DEALLOCATE(xint    )
    DEALLOCATE(xintrng )
    DEALLOCATE(xinprec )
    DEALLOCATE(xincld  )
    DEALLOCATE(xinq    )
    DEALLOCATE(xinwet  )
    DEALLOCATE(gdd0    )
    DEALLOCATE(gdd5    )
    DEALLOCATE(gdd0this)
    DEALLOCATE(tcthis  )
    DEALLOCATE(twthis  )
    DEALLOCATE(gdd5this)
    DEALLOCATE(tc      )
    DEALLOCATE(tw      )
    DEALLOCATE(tcmin   )
    DEALLOCATE(ndtimes )
    DEALLOCATE(nmtimes )
    DEALLOCATE(nytimes )
    DEALLOCATE(nppdummy)
    DEALLOCATE(tco2root)

    DEALLOCATE(wsib     )
    DEALLOCATE(iMaskSSiB)
    DEALLOCATE(exist    )
    DEALLOCATE(fi       )
    DEALLOCATE(Tsfc0IBIS)
    DEALLOCATE(Qsfc0IBIS)
    DEALLOCATE(TsfcmIBIS)
    DEALLOCATE(QsfcmIBIS)
    DEALLOCATE(tsoi0    )
    DEALLOCATE(tsoi     )
    DEALLOCATE(tsoim    )
    DEALLOCATE(hvasug   )
    DEALLOCATE(hvasui   )
    DEALLOCATE(wsoim    )
    DEALLOCATE(wsoi     )
    DEALLOCATE(wsoi0    )
    DEALLOCATE(wisoi    )
    DEALLOCATE(consoi   )
    DEALLOCATE(wpud     )
    DEALLOCATE(wipud    )
    DEALLOCATE(qglif    )
    DEALLOCATE(z0soi    )
    DEALLOCATE(tg      )
    DEALLOCATE(ti      )
    DEALLOCATE(albsav  )
    DEALLOCATE(albsan  )
    DEALLOCATE(rhosoi  )
    DEALLOCATE(csoi    )
    DEALLOCATE(poros   )
    DEALLOCATE(sfield  )
    DEALLOCATE(swilt   )
    DEALLOCATE(bex     )
    DEALLOCATE(upsoiu  )
    DEALLOCATE(upsoil  )
    DEALLOCATE(heatg   )
    DEALLOCATE(heati   )
    DEALLOCATE(ibex    )
    DEALLOCATE(hflo    )
    DEALLOCATE(suction )
    DEALLOCATE(hydraul )
    DEALLOCATE(porosflo)
    DEALLOCATE(stressl )
    DEALLOCATE(stressu )
    DEALLOCATE(stresstu)
    DEALLOCATE(stresstl)
    DEALLOCATE(tsno    )
    DEALLOCATE(hsno    )
    DEALLOCATE(iwet    )
    DEALLOCATE(iwetday )
    DEALLOCATE(precipday)
    DEALLOCATE(asurd   )
    DEALLOCATE(asuri   )
    DEALLOCATE(xstore  )
    DEALLOCATE(ginvap  )
    DEALLOCATE(gsuvap  )
    DEALLOCATE(gtrans  )
    DEALLOCATE(gtransu )
    DEALLOCATE(gtransl )
    DEALLOCATE(grunof  )
    DEALLOCATE(gdrain  )
    DEALLOCATE(gadjust )
    DEALLOCATE(wtot    )
    DEALLOCATE(a10td               )
    DEALLOCATE(a10ancub     )
    DEALLOCATE(a10ancuc     )
    DEALLOCATE(a10ancls     )
    DEALLOCATE(a10ancl4     )
    DEALLOCATE(a10ancl3     )
    DEALLOCATE(a10scalparamu)
    DEALLOCATE(a10scalparaml)
    DEALLOCATE(a10daylightu )
    DEALLOCATE(a10daylightl )
    DEALLOCATE(adrain    )
    DEALLOCATE(adsnow    )
    DEALLOCATE(adaet            )
    DEALLOCATE(adtrunoff )
    DEALLOCATE(adsrunoff )
    DEALLOCATE(addrainage)
    DEALLOCATE(adrh            )
    DEALLOCATE(adsnod    )
    DEALLOCATE(adsnof    )
    DEALLOCATE(adwsoi    )
    DEALLOCATE(adtsoi    )
    DEALLOCATE(adwisoi   )
    DEALLOCATE(adtlaysoi )
    DEALLOCATE(adwlaysoi )
    DEALLOCATE(adwsoic   )
    DEALLOCATE(adtsoic   )
    DEALLOCATE(adco2mic  )
    DEALLOCATE(adco2root )
    DEALLOCATE(adco2soi  )
    DEALLOCATE(adco2ratio)
    DEALLOCATE(adnmintot )

    DEALLOCATE(amtemp    )
    DEALLOCATE(amrain    )
    DEALLOCATE(amsnow    )
    DEALLOCATE(amaet            )
    DEALLOCATE(amtrunoff )
    DEALLOCATE(amsrunoff )
    DEALLOCATE(amdrainage)
    DEALLOCATE(amqa            )
    DEALLOCATE(amsolar   )
    DEALLOCATE(amirup    )
    DEALLOCATE(amirdown  )
    DEALLOCATE(amsens    )
    DEALLOCATE(amlatent  )
    DEALLOCATE(amlaiu    )
    DEALLOCATE(amlail    )
    DEALLOCATE(amtsoi    )
    DEALLOCATE(amwsoi    )
    DEALLOCATE(amwisoi   )
    DEALLOCATE(amvwc            )
    DEALLOCATE(amawc            )
    DEALLOCATE(amsnod    )
    DEALLOCATE(amsnof    )
    DEALLOCATE(amco2mic  )
    DEALLOCATE(amco2root )
    DEALLOCATE(amnmintot )
    DEALLOCATE(tauwood   )
    DEALLOCATE(amnpp     )
    DEALLOCATE(amts2            )
    DEALLOCATE(amtransu  )
    DEALLOCATE(amtransl  )
    DEALLOCATE(amsuvap   )
    DEALLOCATE(aminvap   )
    DEALLOCATE(amneetot  )
    DEALLOCATE(amco2ratio)
    DEALLOCATE(amnpptot  )
    DEALLOCATE(amco2soi  )
    DEALLOCATE(amalbedo  )
    DEALLOCATE(amtsoil   )
    DEALLOCATE(amwsoil   )
    DEALLOCATE(amwisoil  )
    DEALLOCATE(aysolar   )
    DEALLOCATE(ayirup    )
    DEALLOCATE(ayirdown  )
    DEALLOCATE(aysens    )
    DEALLOCATE(aylatent  )
    DEALLOCATE(ayprcp    )
    DEALLOCATE(ayaet            )
    DEALLOCATE(aytrans   )
    DEALLOCATE(aytrunoff )
    DEALLOCATE(aysrunoff )
    DEALLOCATE(aydrainage)
    DEALLOCATE(aydwtot   )
    DEALLOCATE(aygpptot  )
    DEALLOCATE(aynpp            )
    DEALLOCATE(aynpptot  )
    DEALLOCATE(ayco2soi  )
    DEALLOCATE(ayneetot  )
    DEALLOCATE(totnsoi   )
    DEALLOCATE(aywsoi    )
    DEALLOCATE(aywisoi   )
    DEALLOCATE(aytsoi    )
    DEALLOCATE(ayvwc            )
    DEALLOCATE(ayawc            )
    DEALLOCATE(aystresstu)
    DEALLOCATE(aystresstl)
    DEALLOCATE(ayco2mic  )
    DEALLOCATE(ayco2root )
    DEALLOCATE(ayrootbio )
    DEALLOCATE(aynmintot )
    DEALLOCATE(ayalit    )
    DEALLOCATE(ayblit    )
    DEALLOCATE(aycsoi    )
    DEALLOCATE(aycmic    )
    DEALLOCATE(ayanlit   )
    DEALLOCATE(aybnlit   )
    DEALLOCATE(aynsoi    )
    DEALLOCATE(ayalbedo  )
    DEALLOCATE(aygpp     )
    DEALLOCATE(ayanpp    )
    DEALLOCATE(ayanpptot )
    DEALLOCATE(ayrratio  )
    DEALLOCATE(aytratio  )
    DEALLOCATE(totcondub )
    DEALLOCATE(totconduc )
    DEALLOCATE(totcondls )
    DEALLOCATE(totcondl3 )
    DEALLOCATE(totcondl4 )

    DEALLOCATE(frac      )
    DEALLOCATE(tum       )
    DEALLOCATE(tu        )
    DEALLOCATE(tu0       )
    DEALLOCATE(ts            )
    DEALLOCATE(tl            )
    DEALLOCATE(topparu   )
    DEALLOCATE(topparl   )
    DEALLOCATE(agcub     )
    DEALLOCATE(agcuc     )
    DEALLOCATE(ancub     )
    DEALLOCATE(ancuc     )
    DEALLOCATE(tlsub           )
    DEALLOCATE(t12           )
    DEALLOCATE(t34           )
    DEALLOCATE(q12           )
    DEALLOCATE(q34           )
    DEALLOCATE(ciub           )
    DEALLOCATE(ciuc           )
    DEALLOCATE(cils           )
    DEALLOCATE(cil3           )
    DEALLOCATE(cil4           )
    DEALLOCATE(csub           )
    DEALLOCATE(csuc           )
    DEALLOCATE(csls           )
    DEALLOCATE(csl3           )
    DEALLOCATE(csl4           )
    DEALLOCATE(gsub           )
    DEALLOCATE(gsuc           )
    DEALLOCATE(gsls           )
    DEALLOCATE(gsl3           )
    DEALLOCATE(gsl4           )
    DEALLOCATE(agcls           )
    DEALLOCATE(agcl4           )
    DEALLOCATE(agcl3           )
    DEALLOCATE(ancls           )
    DEALLOCATE(ancl4           )
    DEALLOCATE(ancl3           )
    DEALLOCATE(clitlm   )
    DEALLOCATE(clitls   )
    DEALLOCATE(clitll   )
    DEALLOCATE(clitrm   )
    DEALLOCATE(clitrs   )
    DEALLOCATE(clitrl   )
    DEALLOCATE(clitwm   )
    DEALLOCATE(clitws   )
    DEALLOCATE(clitwl   )
    DEALLOCATE(totcmic  )
    DEALLOCATE(csoislop )
    DEALLOCATE(csoislon )
    DEALLOCATE(csoipas  )
    DEALLOCATE(totlit   )
    DEALLOCATE(totnlit  )
    DEALLOCATE(totfall  )
    DEALLOCATE(totalit  )
    DEALLOCATE(totrlit  )
    DEALLOCATE(totanlit )
    DEALLOCATE(totrnlit )
    DEALLOCATE(totcsoi  )
    DEALLOCATE(totnmic  )
    DEALLOCATE(tco2mic  )
    DEALLOCATE(tnpptot  )
    DEALLOCATE(tneetot  )
    DEALLOCATE(tnmin           )
    DEALLOCATE(cdisturb )
    DEALLOCATE(tempu           )
    DEALLOCATE(templ           )
    DEALLOCATE(dropu           )
    DEALLOCATE(dropls   )
    DEALLOCATE(dropl4   )
    DEALLOCATE(dropl3   )
    DEALLOCATE(wliqu           )
    DEALLOCATE(wliqs           )
    DEALLOCATE(wliql           )
    DEALLOCATE(wsnou           )
    DEALLOCATE(wsnos           )
    DEALLOCATE(wsnol           )
    DEALLOCATE(su           )
    DEALLOCATE(ss       )
    DEALLOCATE(sl       )
    DEALLOCATE(agddu    )
    DEALLOCATE(agddl    )
    DEALLOCATE(storedn  )
    DEALLOCATE(yrleach  )
    DEALLOCATE(ynleach  )
    DEALLOCATE(falll           )
    DEALLOCATE(fallr           )
    DEALLOCATE(fallw           )
    DEALLOCATE(vegtype0 )
    DEALLOCATE(plai           )
    DEALLOCATE(sapfrac  )
    DEALLOCATE(cbiol           )
    DEALLOCATE(cbior           )
    DEALLOCATE(cbiow           )
    DEALLOCATE(biomass  )
    DEALLOCATE(totlaiu  )
    DEALLOCATE(totlail  )
    DEALLOCATE(totbiou  )
    DEALLOCATE(totbiol  )
    DEALLOCATE(sai           )
    DEALLOCATE(fu           )
    DEALLOCATE(fl           )
    DEALLOCATE(lai           )
    DEALLOCATE(zbot           )
    DEALLOCATE(ztop           )
    DEALLOCATE(decompl  )
    DEALLOCATE(decomps  )
    DEALLOCATE(firefac  )
    DEALLOCATE(froot    )
    DEALLOCATE(vzero    )
    DEALLOCATE(td           )
  END SUBROUTINE Finalize_IBIS

  ! ---------------------------------------------------------------------
  SUBROUTINE const (arr   , & ! INTENT(OUT  )
       nar   , & ! INTENT(IN        )
       value   ) ! INTENT(IN        )
    ! ---------------------------------------------------------------------
    !
    ! sets all elements of real vector arr to value
    !
    IMPLICIT NONE
    !
    ! Arguments
    !
    INTEGER, INTENT(IN   ) :: nar
    REAL(KIND=r8)   , INTENT(IN   ) :: value
    REAL(KIND=r8)   , INTENT(OUT  ) :: arr(nar)
    !
    ! Local variables
    !
    INTEGER :: j
    !
    DO j = 1, nar
       arr(j) = value
    END DO
    !
    RETURN
  END SUBROUTINE const

END MODULE Sfc_Ibis_Fiels
