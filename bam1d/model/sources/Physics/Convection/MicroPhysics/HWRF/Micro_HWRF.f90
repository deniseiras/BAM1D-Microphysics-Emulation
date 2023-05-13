!WRF:MODEL_MP:PHYSICS
!
MODULE Micro_HWRF
  PRIVATE 
  ! Selecting Kinds
  INTEGER, PARAMETER         :: r4 = SELECTED_REAL_KIND(6)  ! Kind for 32-bits Real Numbers
  INTEGER, PARAMETER         :: i4 = SELECTED_INT_KIND(9)   ! Kind for 32-bits Integer Numbers
  INTEGER, PARAMETER         :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
  INTEGER, PARAMETER         :: i8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers
  INTEGER, PARAMETER         :: r16 = SELECTED_REAL_KIND(15)! Kind for 128-bits Real Numbers
  !-----------------------------------------------------------------------
  !-- The following changes were made on 24 July 2006.
  !   (1) All known version 2.1 dependencies were removed from the 
  !       operational WRF NMM model code (search for "!HWRF")
  !   (2) Incorporated code changes from the GFDL model (search for "!GFDL")
  !-----------------------------------------------------------------------
  REAL(KIND=r8)              :: ABFR
  REAL(KIND=r8)              :: CBFR
  REAL(KIND=r8)              :: CIACW
  REAL(KIND=r8)              :: CIACR
  REAL(KIND=r8)              :: C_N0r0
  REAL(KIND=r8)              :: CN0r0
  REAL(KIND=r8)              :: CN0r_DMRmin
  REAL(KIND=r8)              :: CN0r_DMRmax
  REAL(KIND=r8)              :: CRACW
  REAL(KIND=r8)              :: CRAUT
  REAL(KIND=r8)              :: ESW0
  REAL(KIND=r8)              :: RFmax
  REAL(KIND=r8)              :: RQR_DR1
  REAL(KIND=r8)              :: RQR_DR2
  REAL(KIND=r8)              :: RQR_DR3
  REAL(KIND=r8)              :: RQR_DRmin
  REAL(KIND=r8)              :: RQR_DRmax
  REAL(KIND=r8)              :: RR_DRmin
  REAL(KIND=r8)              :: RR_DR1
  REAL(KIND=r8)              :: RR_DR2
  REAL(KIND=r8)              :: RR_DR3
  REAL(KIND=r8)              :: RR_DRmax
  !
  INTEGER        , PARAMETER :: MY_T1=1
  INTEGER        , PARAMETER :: MY_T2=35
  REAL(KIND=r8)              :: MY_GROWTH(MY_T1:MY_T2)
  !
  REAL(KIND=r8)   ,PARAMETER :: DMImin=0.05e-3_r8
  REAL(KIND=r8)   ,PARAMETER :: DMImax=1.e-3_r8
  REAL(KIND=r8)   ,PARAMETER :: DelDMI=1.e-6_r8
  REAL(KIND=r8)   ,PARAMETER :: XMImin=1.e6_r8*DMImin
  REAL(KIND=r8)   ,PARAMETER :: XMImax=1.e6_r8*DMImax
  REAL(KIND=r8)   ,PARAMETER :: XMIexp=0.0536_r8
  INTEGER         ,PARAMETER :: MDImin=INT(XMImin)
  INTEGER         ,PARAMETER :: MDImax=INT(XMImax)

  REAL(KIND=r8)              :: ACCRI (MDImin:MDImax)
  REAL(KIND=r8)              :: SDENS (MDImin:MDImax)
  REAL(KIND=r8)              :: VSNOWI(MDImin:MDImax)
  REAL(KIND=r8)              :: VENTI1(MDImin:MDImax)
  REAL(KIND=r8)              :: VENTI2(MDImin:MDImax)
  REAL(KIND=r8)              :: MASSI (MDImin:MDImax)

  !
  REAL(KIND=r8)   ,PARAMETER :: DMRmin=0.05e-3_r8
  REAL(KIND=r8)   ,PARAMETER :: DMRmax=0.45e-3_r8
  REAL(KIND=r8)   ,PARAMETER :: DelDMR=1.e-6_r8
  REAL(KIND=r8)   ,PARAMETER :: XMRmin=1.e6_r8*DMRmin
  REAL(KIND=r8)   ,PARAMETER :: XMRmax=1.e6_r8*DMRmax
  INTEGER         ,PARAMETER :: MDRmin=INT(XMRmin)
  INTEGER         ,PARAMETER :: MDRmax=INT(XMRmax)
  REAL(KIND=r8)              :: ACCRR (MDRmin:MDRmax)
  REAL(KIND=r8)              :: MASSR (MDRmin:MDRmax)
  REAL(KIND=r8)              :: RRATE (MDRmin:MDRmax)
  REAL(KIND=r8)              :: VRAIN (MDRmin:MDRmax)
  REAL(KIND=r8)              :: VENTR1(MDRmin:MDRmax)
  REAL(KIND=r8)              :: VENTR2(MDRmin:MDRmax)
  !
  INTEGER         ,PARAMETER :: Nrime=40
  REAL(KIND=r8)              :: VEL_RF(2:9,0:Nrime)
  !
  INTEGER         ,PARAMETER :: NX=7501
  REAL(KIND=r8)   ,PARAMETER :: XMIN=180.0_r8
  REAL(KIND=r8)   ,PARAMETER :: XMAX=330.0_r8

  REAL(KIND=r8)              :: TBPVS (NX)
  REAL(KIND=r8)              :: TBPVS0(NX)
  REAL(KIND=r8)              :: C1XPVS0
  REAL(KIND=r8)              :: C2XPVS0
  REAL(KIND=r8)              :: C1XPVS
  REAL(KIND=r8)              :: C2XPVS
  !
  !--- Physical constants follow:
  REAL(KIND=r8)   ,PARAMETER ::  CP=1004.6_r8
  REAL(KIND=r8)   ,PARAMETER ::  EPSQ=1.E-12_r8
  REAL(KIND=r8)   ,PARAMETER ::  GRAV=9.806_r8
  REAL(KIND=r8)   ,PARAMETER ::  RHOL=1000.0_r8
  REAL(KIND=r8)   ,PARAMETER ::  RD=287.04_r8
  REAL(KIND=r8)   ,PARAMETER ::  RV=461.5_r8
  REAL(KIND=r8)   ,PARAMETER ::  T0C=273.15_r8
  REAL(KIND=r8)   ,PARAMETER ::  XLS=2.834E6_r8
  !--- Derived physical constants follow:
  REAL(KIND=r8)   ,PARAMETER ::  EPS=RD/RV
  REAL(KIND=r8)   ,PARAMETER ::  EPS1=RV/RD-1.0_r8
  REAL(KIND=r8)   ,PARAMETER ::  EPSQ1=1.001_r8*EPSQ
  REAL(KIND=r8)   ,PARAMETER ::  RCP=1.0_r8/CP
  REAL(KIND=r8)   ,PARAMETER ::  RCPRV=RCP/RV
  REAL(KIND=r8)   ,PARAMETER ::  RGRAV=1.0_r8/GRAV
  REAL(KIND=r8)   ,PARAMETER ::  RRHOL=1.0_r8/RHOL
  REAL(KIND=r8)   ,PARAMETER ::  XLS1=XLS*RCP
  REAL(KIND=r8)   ,PARAMETER ::  XLS2=XLS*XLS*RCPRV
  REAL(KIND=r8)   ,PARAMETER ::  XLS3=XLS*XLS/RV
  !--- Constants specific to the parameterization follow:
  !--- CLIMIT/CLIMIT1 are lower limits for treating accumulated precipitation
  REAL(KIND=r8)   ,PARAMETER ::  CLIMIT=10.0_r8*EPSQ
  REAL(KIND=r8)   ,PARAMETER ::  C1=1.0_r8/3.0_r8
  REAL(KIND=r8)   ,PARAMETER ::  DMR1=0.1E-3_r8
  REAL(KIND=r8)   ,PARAMETER ::  DMR2=0.2E-3_r8
  REAL(KIND=r8)   ,PARAMETER ::  DMR3=0.32E-3_r8
  REAL(KIND=r8)   ,PARAMETER ::  XMR1=1.e6_r8*DMR1
  REAL(KIND=r8)   ,PARAMETER ::  XMR2=1.e6_r8*DMR2
  REAL(KIND=r8)   ,PARAMETER ::  XMR3=1.e6_r8*DMR3
  INTEGER         ,PARAMETER ::  MDR1=INT(XMR1)
  INTEGER         ,PARAMETER ::  MDR2=INT(XMR2)
  INTEGER         ,PARAMETER ::  MDR3=INT(XMR3) 

  REAL(KIND=r8)              :: WC
  !
  ! ======================================================================
  !--- Important tunable parameters that are exported to other modules
  !GFDL * RHgrd - generic reference to the threshold relative humidity for 
  !GFDL           the onset of condensation
  !GFDL (new) * RHgrd_in  - "RHgrd" for the inner domain
  !GFDL (new) * RHgrd_out - "RHgrd" for the outer domain
  !HWRF 6/11/2010 mod - use lower RHgrd_out for p < 850 hPa
  !  * T_ICE - temperature (C) threshold at which all remaining liquid water
  !            is glaciated to ice
  !  * T_ICE_init - maximum temperature (C) at which ice nucleation occurs
  !  * NLImax - maximum number concentrations (m**-3) of large ice (snow/graupel/sleet) 
  !  * NLImin - minimum number concentrations (m**-3) of large ice (snow/graupel/sleet) 
  !  * N0r0 - assumed intercept (m**-4) of rain drops if drop diameters are between 0.2 and 0.45 mm
  !  * N0rmin - minimum intercept (m**-4) for rain drops 
  !  * NCW - number concentrations of cloud droplets (m**-3)
  !  * FLARGE1, FLARGE2 - number fraction of large ice to total (large+snow) ice 
  !          at T>0C and in presence of sublimation (FLARGE1), otherwise in
  !          presence of ice saturated/supersaturated conditions
  !  * PRINT_diag - for extended model diagnostics (code currently commented out)
  !  * PRINT_err - for run-time prints when water budgets are not conserved (for debugging)
  !      REAL, PRIVATE,  PARAMETER ::                                      &
  !    &  T_ICE=-10., T_ICE_init=-5.      !- Ver1
  !!!  &, T_ICE=-20.                      !- Ver2
  !    &  T_ICE=-40., T_ICE_init=-15.     !- Ver2
  !    &  T_ICE=-30., T_ICE_init=-5.      !- Ver2
  !
  ! ======================================================================
  REAL(KIND=r8)   ,PARAMETER :: RHgrd_in=1.0_r8         !GFDL
  REAL(KIND=r8)   ,PARAMETER :: RHgrd_out=0.975_r8     !GFDL
  REAL(KIND=r8)   ,PARAMETER :: P_RHgrd_out=850.E2_r8  !HWRF 6/11/2010
  REAL(KIND=r8)   ,PARAMETER :: T_ICE=-40.0_r8          !GFDL
  !REAL(KIND=r8)   ,PARAMETER :: T_ICE=-30.0_r8          !GFDL
  !GFDL     & ,T_ICE=-30.
  REAL(KIND=r8)   ,PARAMETER :: T_ICEK=T0C+T_ICE 
  REAL(KIND=r8)   ,PARAMETER :: T_ICE_init=-15.0_r8      !changed to 20E3 2-09-2012& ,NLImax=5.E3
  REAL(KIND=r8)   ,PARAMETER :: NLImax=20.E3_r8
  REAL(KIND=r8)   ,PARAMETER :: NLImin=0.1E3_r8
!  REAL(KIND=r8)   ,PARAMETER :: NLImin=1.E3_r8
  REAL(KIND=r8)   ,PARAMETER :: N0r0=8.E6_r8
  REAL(KIND=r8)   ,PARAMETER :: N0rmin=1.E4_r8
  !!2-09-2012     & ,NCW=60.E6 &  !GFDL
  !     REAL, PRIVATE,PARAMETER :: DMRmin=.05e-3,      DMRmax=.45e-3,     &
  !     &                           XMRmin=1.e6*DMRmin, XMRmax=1.e6*DMRmax,&
  !     &                           DelDMR=1.e-6,       NLImin=100.
  !    &,                          NLImin=100., NLImax=20.E3
  !      INTEGER, PRIVATE,PARAMETER :: MDRmin=XMRmin, MDRmax=XMRmax

  !! based on Aligo's email,NCW is changed to 250E6
  !REAL(KIND=r8)   ,PARAMETER :: NCW=250.E6_r8          !GFDL    !HWRF     & ,NCW=100.E6
  REAL(KIND=r8)   ,PARAMETER :: NCW=100.E6_r8           !GFDL    !HWRF     & ,NCW=100.E6

  REAL(KIND=r8)   ,PARAMETER :: FLARGE1=1.0_r8
  REAL(KIND=r8)   ,PARAMETER :: FLARGE2=0.2_r8
  !      LOGICAL, PARAMETER       :: PRINT_diag=.FALSE.  !GFDL
  LOGICAL         ,PARAMETER :: PRINT_err=.TRUE.    !GFDL** => eventually set this to .false.Micro_HWRF.f90
  !--- Other public variables passed to other routines:
  REAL(KIND=r8)              :: QAUT0
  !
  !
  PUBLIC :: Init_Micro_HWRF
  PUBLIC :: RunMicro_HWRF
CONTAINS

  !#######################################################################
  !------- Initialize constants & lookup tables for microphysics ---------
  !#######################################################################
  !

  ! SH 0211/2002

  !-----------------------------------------------------------------------
  SUBROUTINE Init_Micro_HWRF ( &
       DT               , &
       LOWLYR           , &
       restart          , &
       fNameMicro       , &
       F_ICE_PHY        , &
       F_RAIN_PHY       , &
       F_RIMEF_PHY      , &
       ibMax            , &
       jbMax            , &
       kMax               )
    !-----------------------------------------------------------------------
    !HWRF     &   MP_RESTART_STATE,TBPVS_STATE,TBPVS0_STATE,
    !-------------------------------------------------------------------------------
    !---  SUBPROGRAM DOCUMENTATION BLOCK
    !   PRGRMMR: Ferrier         ORG: W/NP22     DATE: February 2001
    !            Jin             ORG: W/NP22     DATE: January 2002 
    !        (Modification for WRF structure)
    !                                               
    !-------------------------------------------------------------------------------
    ! ABSTRACT:
    !   * Reads various microphysical lookup tables used in COLUMN_MICRO
    !   * Lookup tables were created "offline" and are read in during execution
    !   * Creates lookup tables for saturation vapor pressure w/r/t water & ice
    !-------------------------------------------------------------------------------
    !     
    ! USAGE: CALL etanewinit FROM SUBROUTINE GSMDRIVE AT MODEL START TIME
    !
    !   INPUT ARGUMENT LIST:
    !       DTPH - physics time step (s)
    !  
    !   OUTPUT ARGUMENT LIST: 
    !     NONE
    !     
    !   OUTPUT FILES:
    !     NONE
    !     
    !   SUBROUTINES:
    !     MY_GROWTH_RATES - lookup table for growth of nucleated ice
    !     GPVS            - lookup tables for saturation vapor pressure (water, ice)
    !
    !   UNIQUE: NONE
    !  
    !   LIBRARY: NONE
    !  
    !   COMMON BLOCKS:
    !     CMICRO_CONS - constants used in GSMCOLUMN
    !     CMY600       - lookup table for growth of ice crystals in 
    !                    water saturated conditions (Miller & Young, 1979)
    !     IVENT_TABLES - lookup tables for ventilation effects of ice
    !     IACCR_TABLES - lookup tables for accretion rates of ice
    !     IMASS_TABLES - lookup tables for mass content of ice
    !     IRATE_TABLES - lookup tables for precipitation rates of ice
    !     IRIME_TABLES - lookup tables for increase in fall speed of rimed ice
    !     MAPOT        - Need lat/lon grid resolution
    !     RVENT_TABLES - lookup tables for ventilation effects of rain
    !     RACCR_TABLES - lookup tables for accretion rates of rain
    !     RMASS_TABLES - lookup tables for mass content of rain
    !     RVELR_TABLES - lookup tables for fall speeds of rain
    !     RRATE_TABLES - lookup tables for precipitation rates of rain
    !   
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !   MACHINE : IBM SP
    !
    !-----------------------------------------------------------------------
    !
    !
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    !-----------------------------------------------------------------------
    !------------------------------------------------------------------------- 
    !-------------- Parameters & arrays for lookup tables -------------------- 
    !------------------------------------------------------------------------- 
    !
    !--- Common block of constants used in column microphysics
    !
    !WRF
    !     REAL(KIND=r8) DLMD,DPHD
    !WRF
    !
    !
    !     VARIABLES PASSED IN
    INTEGER,INTENT(IN) :: ibMax
    INTEGER,INTENT(IN) :: jbMax
    INTEGER,INTENT(IN) :: kMax
    LOGICAL,INTENT(IN) :: restart
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameMicro

    !WRF
    INTEGER,INTENT(INOUT) :: LOWLYR(1:ibMax,1:jbMax)
    !
    REAL(KIND=r8),INTENT(OUT) :: F_ICE_PHY  (1:ibMax, 1:kMax, 1:jbMax)
    REAL(KIND=r8),INTENT(OUT) :: F_RAIN_PHY (1:ibMax, 1:kMax, 1:jbMax)
    REAL(KIND=r8),INTENT(OUT) :: F_RIMEF_PHY(1:ibMax, 1:kMax, 1:jbMax)
    INTEGER,PARAMETER   :: ITLO=-60
    INTEGER,PARAMETER   :: ITHI=40
    REAL(KIND=r8)   ,PARAMETER   :: GSMDT=0.0_r8


    !HWRF      REAL(KIND=r8),DIMENSION(*), INTENT(INOUT) :: MP_RESTART_STATE
    !HWRF      REAL(KIND=r8),DIMENSION(NX), INTENT(INOUT) :: TBPVS_STATE,TBPVS0_STATE
    !     integer,DIMENSION(ITLO:ITHI,4),INTENT(INOUT) :: NSTATS
    !     REAL(KIND=r8),DIMENSION(ITLO:ITHI,5),INTENT(INOUT) :: QMAX
    !     REAL(KIND=r8),DIMENSION(ITLO:ITHI,22),INTENT(INOUT) :: QTOT
    !     REAL(KIND=r8),INTENT(INOUT) :: PRECtot(2),PRECmax(2)
    REAL(KIND=r8)   ,INTENT(IN) :: DT
    !
    !-----------------------------------------------------------------------
    !     LOCAL VARIABLES
    !-----------------------------------------------------------------------

    REAL(KIND=r8)    :: BBFR
    REAL(KIND=r8)    :: DTPH
    REAL(KIND=r8)    :: PI
    !REAL(KIND=r8)    :: DX
    REAL(KIND=r8)    :: Thour_print
    INTEGER :: I
    INTEGER :: J
    INTEGER :: K
    INTEGER :: etampnew_unit1
    LOGICAL, PARAMETER :: wrf_dm_on_monitor=.TRUE.
    !-----------------------------------------------------------------------
    !--- Parameters & data statement for local calculations
    !-----------------------------------------------------------------------
    !
    INTEGER, PARAMETER :: MDR1=INT(XMR1)
    INTEGER, PARAMETER :: MDR2=INT(XMR2)
    INTEGER, PARAMETER :: MDR3=INT(XMR3)
    REAL(KIND=r4)              :: VENTR1_R4(MDRmin:MDRmax)
    REAL(KIND=r4)              :: VENTR2_R4(MDRmin:MDRmax)
    REAL(KIND=r4)              :: ACCRR_R4 (MDRmin:MDRmax)
    REAL(KIND=r4)              :: MASSR_R4 (MDRmin:MDRmax)
    REAL(KIND=r4)              :: VRAIN_R4 (MDRmin:MDRmax)
    REAL(KIND=r4)              :: RRATE_R4 (MDRmin:MDRmax)

    REAL(KIND=r4)              :: VENTI1_R4(MDImin:MDImax)
    REAL(KIND=r4)              :: VENTI2_R4(MDImin:MDImax)
    REAL(KIND=r4)              :: ACCRI_R4 (MDImin:MDImax)
    REAL(KIND=r4)              :: MASSI_R4 (MDImin:MDImax)
    REAL(KIND=r4)              :: VSNOWI_R4(MDImin:MDImax)
    REAL(KIND=r4)              :: VEL_RF_R4(2:9,0:Nrime)



    !
    !-----------------------------------------------------------------------
    !
    !
    DO J=1,jbMax
       DO I=1,ibMax
          LOWLYR(I,J)=1
       ENDDO
    ENDDO
    !    
    IF(.NOT.RESTART) THEN    !HWRF
       DO J = 1,jbMax
          DO K = 1,kMax
             DO I= 1,ibMax
                F_ICE_PHY  (i,k,j)=0.0_r8
                F_RAIN_PHY (i,k,j)=0.0_r8
                F_RIMEF_PHY(i,k,j)=1.0_r8
             END DO
          END DO
       END DO
    END IF
    !    
    !-----------------------------------------------------------------------
    !IF(ALLOWED_TO_READ)THEN
       !-----------------------------------------------------------------------
       !
       !DX=SQRT((DELX)**2+(DELY)**2)/1000.0_r8    ! Model resolution at equator (km)  !GFDL
       !DX=MIN(100.0_r8, MAX(5.0_r8, DX) )
       !
       !-- Relative humidity threshold for the onset of grid-scale condensation
       !!-- 9/1/01:  Assume the following functional dependence for 5 - 100 km resolution:
       !!       RHgrd=0.90 for dx=100 km, 0.98 for dx=5 km, where
       !        RHgrd=0.90+.08*((100.-DX)/95.)**.5
       !
       DTPH=MAX (GSMDT*60.0_r8,DT)
       DTPH=NINT(DTPH/DT)*DT
       !
       !--- Create lookup tables for saturation vapor pressure w/r/t water & ice
       !
       CALL GPVS(NX,XMAX,XMIN,C1XPVS,C2XPVS,C1XPVS0,C2XPVS0,TBPVS,TBPVS0)
       !
       !--- Read in various lookup tables
       !
       !        IF ( wrf_dm_on_monitor ) THEN
       !          DO i = 31,99
       !            INQUIRE ( i , OPENED = opened )
       !            IF ( .NOT. opened ) THEN
       !              etampnew_unit1 = i
       !              GOTO 2061
       !            ENDIF
       !          ENDDO
       !          etampnew_unit1 = -1
       ! 2061     CONTINUE
       !        ENDIF
       !
       !        CALL wrf_dm_bcast_bytes ( etampnew_unit1 , IWORDSIZE )
       !!
       !        IF ( etampnew_unit1 < 0 ) THEN
       !          CALL wrf_error_fatal ( 'Micro_HWRF: etanewinit: Can not find unused fortran unit to read in lookup table.' )
       !        ENDIF
       !
       !        IF ( wrf_dm_on_monitor ) THEN
       !!was     OPEN (UNIT=1,FILE="eta_micro_lookup.dat",FORM="UNFORMATTED")
       etampnew_unit1=1
       print*, 'fNameMicro=',fNameMicro
       OPEN(UNIT=etampnew_unit1,FILE=TRIM(fNameMicro), &
            FORM="UNFORMATTED",STATUS="OLD",ERR=9061)

       !
       READ(etampnew_unit1) VENTR1_R4
       READ(etampnew_unit1) VENTR2_R4
       READ(etampnew_unit1) ACCRR_R4
       READ(etampnew_unit1) MASSR_R4
       READ(etampnew_unit1) VRAIN_R4
       READ(etampnew_unit1) RRATE_R4
       
       READ(etampnew_unit1) VENTI1_R4
       READ(etampnew_unit1) VENTI2_R4
       READ(etampnew_unit1) ACCRI_R4
       READ(etampnew_unit1) MASSI_R4
       READ(etampnew_unit1) VSNOWI_R4
       READ(etampnew_unit1) VEL_RF_R4
       VENTR1 =  VENTR1_R4
       VENTR2 =  VENTR2_R4
       ACCRR  =  ACCRR_R4
       MASSR  =  MASSR_R4
       VRAIN  =  VRAIN_R4
       RRATE  =  RRATE_R4

       VENTI1 =  VENTI1_R4
       VENTI2 =  VENTI2_R4
       ACCRI  =  ACCRI_R4
       MASSI  =  MASSI_R4
       VSNOWI =  VSNOWI_R4
       VEL_RF =  VEL_RF_R4
       !        read(etampnew_unit1) my_growth    ! Applicable only for DTPH=180 s
       CLOSE (etampnew_unit1)
       !        ENDIF
       !
       !
       !--- Calculates coefficients for growth rates of ice nucleated in water
       !    saturated conditions, scaled by physics time step (lookup table)
       !
       CALL MY_GROWTH_RATES (DTPH)
       !       CALL MY_GROWTH_RATES (DTPH,MY_GROWTH)
       !
       PI=ACOS(-1.0_r8)
       !
       !--- Constants associated with Biggs (1953) freezing of rain, as parameterized
       !    following Lin et al. (JCAM, 1983) & Reisner et al. (1998, QJRMS).
       !
       ABFR=-0.66_r8
       BBFR=100.0_r8
       CBFR=20.0_r8*PI*PI*BBFR*RHOL*1.E-21_r8
       !
       !--- CIACW is used in calculating riming rates
       !      The assumed effective collection efficiency of cloud water rimed onto
       !      ice is =0.5_r8 below:
       !
       CIACW=DTPH*0.25_r8*PI*0.5_r8*(1.E5_r8)**C1
       !
       !--- CIACR is used in calculating freezing of rain colliding with large ice
       !      The assumed collection efficiency is 1.0
       !
       CIACR=PI*DTPH
       !
       !--- Based on rain lookup tables for mean diameters from 0.05 to 0.45 mm
       !    * Four different functional relationships of mean drop diameter as 
       !      a function of rain rate (RR), derived based on simple fits to 
       !      mass-weighted fall speeds of rain as functions of mean diameter
       !      from the lookup tables.  
       !
       RR_DRmin=N0r0*RRATE(MDRmin)     ! RR for mean drop diameter of .05 mm
       RR_DR1=N0r0*RRATE(MDR1)         ! RR for mean drop diameter of .10 mm
       RR_DR2=N0r0*RRATE(MDR2)         ! RR for mean drop diameter of .20 mm
       RR_DR3=N0r0*RRATE(MDR3)         ! RR for mean drop diameter of .32 mm
       RR_DRmax=N0r0*RRATE(MDRmax)     ! RR for mean drop diameter of .45 mm
       !
       RQR_DRmin=N0r0*MASSR(MDRmin)    ! Rain content for mean drop diameter of .05 mm
       RQR_DR1=N0r0*MASSR(MDR1)        ! Rain content for mean drop diameter of .10 mm
       RQR_DR2=N0r0*MASSR(MDR2)        ! Rain content for mean drop diameter of .20 mm
       RQR_DR3=N0r0*MASSR(MDR3)        ! Rain content for mean drop diameter of .32 mm
       RQR_DRmax=N0r0*MASSR(MDRmax)    ! Rain content for mean drop diameter of .45 mm
       C_N0r0=PI*RHOL*N0r0
       CN0r0=1.E6_r8/C_N0r0**0.25_r8
       CN0r_DMRmin=1.0_r8/(PI*RHOL*DMRmin**4)
       CN0r_DMRmax=1.0_r8/(PI*RHOL*DMRmax**4)
       !
       !--- CRACW is used in calculating collection of cloud water by rain (an
       !      assumed collection efficiency of 1.0)
       !
       CRACW=DTPH*0.25_r8*PI*1.0_r8
       !
       ESW0=1000.0_r8*FPVS0(T0C)     ! Saturation vapor pressure at 0C
       RFmax=1.1_r8**Nrime          ! Maximum rime factor allowed
       !
       !------------------------------------------------------------------------
       !--------------- Constants passed through argument list -----------------
       !------------------------------------------------------------------------
       !
       !--- Important parameters for self collection (autoconversion) of 
       !    cloud water to rain. 
       !
       !--- CRAUT is proportional to the rate that cloud water is converted by
       !      self collection to rain (autoconversion rate)
       !
       CRAUT=1.0_r8-EXP(-1.E-3_r8*DTPH)
       !
       !--- QAUT0 is the threshold cloud content for autoconversion to rain 
       !      needed for droplets to reach a diameter of 20 microns (following
       !      Manton and Cotton, 1977; Banta and Hanson, 1987, JCAM)
       !--- QAUT0=1.2567, 0.8378, or 0.4189 g/m**3 for droplet number concentrations
       !          of 300, 200, and 100 cm**-3, respectively
       !
       QAUT0=PI*RHOL*NCW*(20.E-6_r8)**3/6.0_r8
       !
       !--- For calculating snow optical depths by considering bulk density of
       !      snow based on emails from Q. Fu (6/27-28/01), where optical 
       !      depth (T) = 1.5*SWP/(Reff*DENS), SWP is snow water path, Reff 
       !      is effective radius, and DENS is the bulk density of snow.
       !
       !    SWP (kg/m**2)=(1.E-3 kg/g)*SWPrad, SWPrad in g/m**2 used in radiation
       !    T = 1.5*1.E3*SWPrad/(Reff*DENS)
       !  
       !    See derivation for MASSI(INDEXS), note equal to RHO*QSNOW/NSNOW
       !
       !      SDENS=1.5e3/DENS, DENS=MASSI(INDEXS)/[PI*(1.E-6*INDEXS)**3]
       !
       DO I=MDImin,MDImax
          SDENS(I)=PI*1.5E-15_r8*FLOAT(I*I*I)/MASSI(I)
       ENDDO
       !
       Thour_print=-DTPH/3600.0_r8


    !ENDIF  ! Allowed_to_read

    RETURN
    !
    !-----------------------------------------------------------------------
    !
9061 CONTINUE
    STOP 'Micro_HWRF: error opening ETAMPNEW_DATA on unit '
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE Init_Micro_HWRF
  !

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE RunMicro_HWRF (&
       nCols       , &!INTEGER      , INTENT(IN)     :: nCols
       kMax        , &!INTEGER      , INTENT(IN)     :: kMax
       DT          , &!REAL(KIND=r8), INTENT(IN)     :: DT
       dz8w        , &!REAL(KIND=r8), INTENT(IN),    :: dz8w       (ims:ime, kms:kme, jms:jme)
       rho_phy     , &!REAL(KIND=r8), INTENT(IN),    :: rho_phy    (ims:ime, kms:kme, jms:jme)
       p_phy       , &!REAL(KIND=r8), INTENT(IN),    :: p_phy      (ims:ime, kms:kme, jms:jme)
       LOWLYR      , &!INTEGER      , INTENT(IN   )  :: LOWLYR     (ims:ime         , jms:jme)
       gt          , &!REAL(KIND=r8), INTENT(INOUT), :: gt         (ims:ime, kms:kme, jms:jme)
       qv          , &!REAL(KIND=r8), INTENT(INOUT), :: qv         (ims:ime, kms:kme, jms:jme)
       QC          , &!REAL(KIND=r8), INTENT(INOUT), :: qc         (ims:ime, kms:kme, jms:jme)
       QI          , &!REAL(KIND=r8), INTENT(INOUT), :: qi         (ims:ime, kms:kme, jms:jme)
       QR          , &!REAL(KIND=r8), INTENT(INOUT), :: qr         (ims:ime, kms:kme, jms:jme)
       F_ICE_PHY   , &!REAL(KIND=r8), INTENT(INOUT), :: F_ICE_PHY  (ims:ime, kms:kme, jms:jme)
       F_RAIN_PHY  , &!REAL(KIND=r8), INTENT(INOUT), :: F_RAIN_PHY (ims:ime, kms:kme, jms:jme)
       F_RIMEF_PHY , &!REAL(KIND=r8), INTENT(INOUT), :: F_RIMEF_PHY(ims:ime, kms:kme, jms:jme)
       RAINNCV     , &!REAL(KIND=r8), INTENT(INOUT), :: RAINNCV    (ims:ime,          jms:jme) !GID
       SNOWNCV       )!REAL(KIND=r8), INTENT(INOUT), :: SNOWNCV    (ims:ime,          jms:jme) !GID
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    !-----------------------------------------------------------------------

    INTEGER         ,INTENT(IN   ) :: nCols
    INTEGER         ,INTENT(IN   ) :: kMax
    REAL(KIND=r8)   ,INTENT(IN   ) :: DT
    REAL(KIND=r8)   ,INTENT(IN   ) :: dz8w       (1:nCols,1:kMax)
    REAL(KIND=r8)   ,INTENT(IN   ) :: p_phy      (1:nCols,1:kMax)
    REAL(KIND=r8)   ,INTENT(IN   ) :: rho_phy    (1:nCols,1:kMax)
    REAL(KIND=r8)   ,INTENT(INOUT) :: gt         (1:nCols,1:kMax)
    REAL(KIND=r8)   ,INTENT(INOUT) :: qv         (1:nCols,1:kMax)
    REAL(KIND=r8)   ,INTENT(INOUT) :: qc         (1:nCols,1:kMax)
    REAL(KIND=r8)   ,INTENT(INOUT) :: qi         (1:nCols,1:kMax)
    REAL(KIND=r8)   ,INTENT(INOUT) :: qr         (1:nCols,1:kMax)
    REAL(KIND=r8)   ,INTENT(INOUT) :: F_ICE_PHY  (1:nCols,1:kMax)
    REAL(KIND=r8)   ,INTENT(INOUT) :: F_RAIN_PHY (1:nCols,1:kMax)
    REAL(KIND=r8)   ,INTENT(INOUT) :: F_RIMEF_PHY(1:nCols,1:kMax)
!    REAL(KIND=r8)   ,INTENT(INOUT) :: RAINNC     (1:nCols)
    REAL(KIND=r8)   ,INTENT(OUT  ) :: RAINNCV    (1:nCols)
    REAL(KIND=r8)   ,INTENT(OUT  ) :: SNOWNCV    (1:nCols)    
    INTEGER         ,INTENT(IN   ) :: LOWLYR     (1:nCols)

    !-----------------------------------------------------------------------
    !     LOCAL VARS
    !-----------------------------------------------------------------------
    INTEGER, PARAMETER  :: GID=1

    !     NSTATS,QMAX,QTOT are diagnostic vars

    !     SOME VARS WILL BE USED FOR DATA ASSIMILATION (DON'T NEED THEM NOW). 
    !     THEY ARE TREATED AS LOCAL VARS, BUT WILL BECOME STATE VARS IN THE 
    !     FUTURE. SO, WE DECLARED THEM AS MEMORY SIZES FOR THE FUTURE USE

    !     TLATGS_PHY,TRAIN_PHY,APREC,PREC,ACPREC,SR are not directly related 
    !     the microphysics scheme. Instead, they will be used by Eta precip 
    !     assimilation.
    REAL(KIND=r8) :: SR         (1:nCols)
    REAL(KIND=r8) :: TLATGS_PHY(1:nCols,1:kMax)
    REAL(KIND=r8) :: TRAIN_PHY (1:nCols,1:kMax)
    REAL(KIND=r8) :: APREC     (1:nCols)
    REAL(KIND=r8) :: ASNOW     (1:nCols)

    REAL(KIND=r8) :: PREC      (1:nCols)
    REAL(KIND=r8) :: ACPREC    (1:nCols)
    REAL(KIND=r8) :: t_phy     (1:nCols,1:kMax)
    REAL(KIND=r8) :: CWM_PHY   (1:nCols,1:kMax)

    INTEGER :: I,K
    !
    !-----------------------------------------------------------------------
    !**********************************************************************
    !-----------------------------------------------------------------------
    !
    !
    DO k = 1,kMax
       DO i = 1,nCols
          t_phy(i,k) = gt(i,k)
!          qv(i,k)=qv(i,k)/(1.0_r8+qv(i,k)) !Convert to specific humidity
       ENDDO
    ENDDO

    !     initial diagnostic variables and data assimilation vars
    !     (will need to delete this part in the future)


    ! initial data assimilation vars (will need to delete this part in the future)

    DO k = 1,kMax
       DO i = 1,nCols
          TLATGS_PHY (i,k)=0.0_r8
          TRAIN_PHY  (i,k)=0.0_r8
       ENDDO
    ENDDO

    DO i = 1,nCols
       ACPREC(i)=0.0_r8
       APREC (i)=0.0_r8 
       ASNOW (i)=0.0_r8
       PREC  (i)=0.0_r8
       SR    (i)=0.0_r8
    ENDDO

    !-- 6/11/2010: Update CWM_PHY, F_ice, F_rain arrays
    DO k = 1,kMax
       DO i = 1,nCols
          CWM_PHY(I,K)=QC(I,K)+QR(I,K)+QI(I,K)
          IF (QI(I,K) <= EPSQ) THEN
             F_ICE_PHY(I,K)=0.0_r8
             IF (T_PHY(I,K) < T_ICEK) F_ICE_PHY(I,K)=1.0_r8
          ELSE
             F_ICE_PHY(I,K)=MAX( 0.0_r8, MIN(1.0_r8, QI(I,K)/CWM_PHY(I,K) ) )
          ENDIF
          IF (QR(I,K) <= EPSQ) THEN
             F_RAIN_PHY(I,K)=0.0_r8
          ELSE
             F_RAIN_PHY(I,K)=QR(I,K)/(QR(I,K)+QC(I,K))
          ENDIF
       ENDDO
    ENDDO

    !-----------------------------------------------------------------------

    CALL EGCP01DRV( &
         GID        , &!INTEGER,INTENT(IN   ) :: GID         ! grid%id gopal's doing
         DT         , &!REAL(KIND=r8)   ,INTENT(IN   ) :: DTPH
         LOWLYR     , &!INTEGER,INTENT(IN   ) :: LOWLYR      ( 1:nCols )
         APREC      , &!REAL(KIND=r8)   ,INTENT(INOUT) :: APREC       (1:nCols)
         ASNOW      , &!REAL(KIND=r8)   ,INTENT(INOUT) :: ASNOW       (1:nCols)
         PREC       , &!REAL(KIND=r8)   ,INTENT(INOUT) :: PREC        (1:nCols)
         ACPREC     , &!REAL(KIND=r8)   ,INTENT(INOUT) :: ACPREC      (1:nCols)
         SR         , &!REAL(KIND=r8)   ,INTENT(INOUT) :: SR          (1:nCols)
         dz8w       , &!REAL(KIND=r8)   ,INTENT(IN   ) :: dz8w        ( 1:nCols, 1:kMax )
         rho_phy    , &!REAL(KIND=r8)   ,INTENT(IN   ) :: RHO_PHY     ( 1:nCols, 1:kMax )
         CWM_PHY    , &!REAL(KIND=r8)   ,INTENT(INOUT) :: CWM_PHY     ( 1:nCols, 1:kMax )
         t_phy      , &!REAL(KIND=r8)   ,INTENT(INOUT) :: t_phy       ( 1:nCols, 1:kMax )
         qv         , &!REAL(KIND=r8)   ,INTENT(INOUT) :: Q_PHY       ( 1:nCols, 1:kMax )
         F_ICE_PHY  , &!REAL(KIND=r8)   ,INTENT(INOUT) :: F_ICE_PHY   ( 1:nCols, 1:kMax )
         P_PHY      , &!REAL(KIND=r8)   ,INTENT(IN   ) :: P_PHY       ( 1:nCols, 1:kMax )
         F_RAIN_PHY , &!REAL(KIND=r8)   ,INTENT(INOUT) :: F_RAIN_PHY  ( 1:nCols, 1:kMax )
         F_RIMEF_PHY, &!REAL(KIND=r8)   ,INTENT(INOUT) :: F_RIMEF_PHY ( 1:nCols, 1:kMax )
         TLATGS_PHY , &!REAL(KIND=r8)   ,INTENT(INOUT) :: TLATGS_PHY  ( 1:nCols, 1:kMax )
         TRAIN_PHY  , &!REAL(KIND=r8)   ,INTENT(INOUT) :: TRAIN_PHY   ( 1:nCols, 1:kMax )
         nCols        , &!INTEGER,INTENT(IN   ) :: nCols
         kMax          )!INTEGER,INTENT(IN   ) :: kMax
    !-----------------------------------------------------------------------

    DO k = 1,kMax
       DO i = 1,nCols
          gt(i,k) = t_phy(i,k)
          !qv(i,k)=qv(i,k)/(1.0_r8-qv(i,k))  !Convert to mixing ratio
          WC=CWM_PHY(I,K)
          QI(I,K)=0.0_r8
          QR(I,K)=0.0_r8
          QC(I,K)=0.0_r8
          IF(F_ICE_PHY(I,K)>=1.0_r8)THEN
             QI(I,K)=WC
          ELSEIF(F_ICE_PHY(I,K)<=0.0_r8)THEN
             QC(I,K)=WC
          ELSE
             QI(I,K)=F_ICE_PHY(I,K)*WC
             QC(I,K)=WC-QI(I,K)
          ENDIF
          !
          IF(QC(I,K)>0.0_r8.AND.F_RAIN_PHY(I,K)>0.0_r8)THEN
             IF(F_RAIN_PHY(I,K).GE.1.0_r8)THEN
                QR(I,K)=QC(I,K)
                QC(I,K)=0.0_r8
             ELSE
                QR(I,K)=F_RAIN_PHY(I,K)*QC(I,K)
                QC(I,K)=QC(I,K)-QR(I,K)
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    ! 
    ! update rain (from m to mm)

    DO i=1,nCols
    !   RAINNC (i)=APREC(i)*1000.0_r8+RAINNC(i)
       RAINNCV(i)=0.5_r8*APREC(i)!*1000.0_r8
       SNOWNCV(i)=0.5_r8*ASNOW(i)!*1000.0_r8

    ENDDO
    !
    !-----------------------------------------------------------------------

  END SUBROUTINE RunMicro_HWRF

  !-----------------------------------------------------------------------

  SUBROUTINE EGCP01DRV( &
       GID        , &!INTEGER,INTENT(IN   ) :: GID     ! grid%id gopal's doing
       DTPH       , &!REAL(KIND=r8)   ,INTENT(IN   ) :: DTPH
       LOWLYR     , &!INTEGER,INTENT(IN   ) :: LOWLYR( 1:nCols )
       APREC      , &!REAL(KIND=r8)   ,INTENT(INOUT) :: APREC (1:nCols)
       ASNOW      , &!REAL(KIND=r8)   ,INTENT(INOUT) :: ASNOW (1:nCols)
       PREC       , &!REAL(KIND=r8)   ,INTENT(INOUT) :: PREC  (1:nCols)
       ACPREC     , &!REAL(KIND=r8)   ,INTENT(INOUT) :: ACPREC(1:nCols)
       SR         , &!REAL(KIND=r8)   ,INTENT(INOUT) :: SR    (1:nCols)
       dz8w       , &!REAL(KIND=r8)   ,INTENT(IN   ) :: dz8w        ( 1:nCols, 1:kMax )
       RHO_PHY    , &!REAL(KIND=r8)   ,INTENT(IN   ) :: RHO_PHY     ( 1:nCols, 1:kMax )
       CWM_PHY    , &!REAL(KIND=r8)   ,INTENT(INOUT) :: CWM_PHY     ( 1:nCols, 1:kMax )
       T_PHY      , &!REAL(KIND=r8)   ,INTENT(INOUT) :: t_phy       ( 1:nCols, 1:kMax )
       Q_PHY      , &!REAL(KIND=r8)   ,INTENT(INOUT) :: Q_PHY       ( 1:nCols, 1:kMax )
       F_ICE_PHY  , &!REAL(KIND=r8)   ,INTENT(INOUT) :: F_ICE_PHY   ( 1:nCols, 1:kMax )
       P_PHY      , &!REAL(KIND=r8)   ,INTENT(IN   ) :: P_PHY       ( 1:nCols, 1:kMax )
       F_RAIN_PHY , &!REAL(KIND=r8)   ,INTENT(INOUT) :: F_RAIN_PHY  ( 1:nCols, 1:kMax )
       F_RIMEF_PHY, &!REAL(KIND=r8)   ,INTENT(INOUT) :: F_RIMEF_PHY ( 1:nCols, 1:kMax )
       TLATGS_PHY , &!REAL(KIND=r8)   ,INTENT(INOUT) :: TLATGS_PHY  ( 1:nCols, 1:kMax )
       TRAIN_PHY  , &!REAL(KIND=r8)   ,INTENT(INOUT) :: TRAIN_PHY   ( 1:nCols, 1:kMax )
       nCols      , &!INTEGER,INTENT(IN   ) :: nCols
       kMax         )!INTEGER,INTENT(IN   ) :: kMax
    !-----------------------------------------------------------------------
    ! DTPH           Physics time step (s)
    ! CWM_PHY ( )    Mixing ratio of the total condensate. kg/kg
    ! Q_PHY          Mixing ratio of water vapor. kg/kg
    ! F_RAIN_PHY     Fraction of rain. 
    ! F_ICE_PHY      Fraction of ice.
    ! F_RIMEF_PHY    Mass ratio of rimed ice (rime factor).
    !
    !TLATGS_PHY,TRAIN_PHY,APREC,PREC,ACPREC,SR are not directly related the
    !micrphysics sechme. Instead, they will be used by Eta precip assimilation.
    !
    !NSTATS,QMAX,QTOT are used for diagnosis purposes.
    !
    !-----------------------------------------------------------------------
    !--- Variables APREC,PREC,ACPREC,SR are calculated for precip assimilation
    !    and/or ZHAO's scheme in Eta and are not required by this microphysics 
    !    scheme itself.  
    !--- NSTATS,QMAX,QTOT are used for diagnosis purposes only.  They will be 
    !    printed out when PRINT_diag is true.
    !
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    !-----------------------------------------------------------------------
    !
    !     VARIABLES PASSED IN/OUT
    INTEGER,INTENT(IN   ) :: nCols
    INTEGER,INTENT(IN   ) :: kMax
    INTEGER,INTENT(IN   ) :: GID     ! grid%id gopal's doing
    REAL(KIND=r8)   ,INTENT(IN   ) :: DTPH
    INTEGER,INTENT(IN   ) :: LOWLYR( 1:nCols )
    REAL(KIND=r8)   ,INTENT(INOUT) :: APREC (1:nCols)
    REAL(KIND=r8)   ,INTENT(INOUT) :: ASNOW (1:nCols)

    REAL(KIND=r8)   ,INTENT(INOUT) :: PREC  (1:nCols)
    REAL(KIND=r8)   ,INTENT(INOUT) :: ACPREC(1:nCols)
    REAL(KIND=r8)   ,INTENT(INOUT) :: SR    (1:nCols)

    REAL(KIND=r8)   ,INTENT(INOUT) :: t_phy       ( 1:nCols, 1:kMax )
    REAL(KIND=r8)   ,INTENT(IN   ) :: dz8w        ( 1:nCols, 1:kMax )
    REAL(KIND=r8)   ,INTENT(IN   ) :: P_PHY       ( 1:nCols, 1:kMax )
    REAL(KIND=r8)   ,INTENT(IN   ) :: RHO_PHY     ( 1:nCols, 1:kMax )
    REAL(KIND=r8)   ,INTENT(INOUT) :: CWM_PHY     ( 1:nCols, 1:kMax )
    REAL(KIND=r8)   ,INTENT(INOUT) :: F_ICE_PHY   ( 1:nCols, 1:kMax )
    REAL(KIND=r8)   ,INTENT(INOUT) :: F_RAIN_PHY  ( 1:nCols, 1:kMax )
    REAL(KIND=r8)   ,INTENT(INOUT) :: F_RIMEF_PHY ( 1:nCols, 1:kMax )
    REAL(KIND=r8)   ,INTENT(INOUT) :: TLATGS_PHY  ( 1:nCols, 1:kMax )
    REAL(KIND=r8)   ,INTENT(INOUT) :: Q_PHY       ( 1:nCols, 1:kMax )
    REAL(KIND=r8)   ,INTENT(INOUT) :: TRAIN_PHY   ( 1:nCols, 1:kMax )
    !
    !-----------------------------------------------------------------------
    !LOCAL VARIABLES
    !-----------------------------------------------------------------------
    !
    !HWRF - Below are directives in the operational code that have been removed,
    !       where "TEMP_DEX" has been replaced with "I,J,L" and "TEMP_DIMS" has
    !       been replaced with "its:nCols,kMax"
    !HWRF#define CACHE_FRIENDLY_MP_ETANEW
    !HWRF#ifdef CACHE_FRIENDLY_MP_ETANEW
    !HWRF#  define TEMP_DIMS  kMax,its:nCols
    !HWRF#  define TEMP_DEX   L,I,J
    !HWRF#else
    !HWRF#  define TEMP_DIMS  its:nCols,kMax
    !HWRF#  define TEMP_DEX   I,J,L
    !HWRF#endif
    !HWRF!
    !HWRF      REAL(KIND=r8),DIMENSION(TEMP_DIMS) :: CWM,T,Q,TRAIN,TLATGS,P
    REAL(KIND=r8)    :: CWM      (1:nCols,1:kMax)
    REAL(KIND=r8)    :: T        (1:nCols,1:kMax)
    REAL(KIND=r8)    :: Q        (1:nCols,1:kMax)
    REAL(KIND=r8)    :: TRAIN    (1:nCols,1:kMax)
    REAL(KIND=r8)    :: TLATGS   (1:nCols,1:kMax)
    REAL(KIND=r8)    :: P        (1:nCols,1:kMax)
    REAL(KIND=r8)    :: F_ice    (1:kMax,1:nCols)
    REAL(KIND=r8)    :: F_rain   (1:kMax,1:nCols)
    REAL(KIND=r8)    :: F_RimeF  (1:kMax,1:nCols)
    INTEGER :: LMH      (1:nCols)
    REAL(KIND=r8)    :: P_col    (1:kMax)
    REAL(KIND=r8)    :: T_col    (1:kMax)
    REAL(KIND=r8)    :: QV_col   (1:kMax)
    REAL(KIND=r8)    :: WC_col   (1:kMax)
    REAL(KIND=r8)    :: RimeF_col(1:kMax)
    REAL(KIND=r8)    :: QI_col   (1:kMax)
    REAL(KIND=r8)    :: QR_col   (1:kMax)
    REAL(KIND=r8)    :: QW_col   (1:kMax)
    REAL(KIND=r8)    :: THICK_col(1:kMax)
    REAL(KIND=r8)    :: RHC_col  (1:kMax)
    REAL(KIND=r8)    :: DPCOL    (1:kMax)    !GFDL
    !REAL(KIND=r8)    :: PRECtot(2)
    !REAL(KIND=r8)    :: PRECmax(2)
    REAL(KIND=r8)    :: TC,WC,QI,QR,QW,Fice,Frain,DUM,ASNOW2,ARAIN
    INTEGER :: LSFC,I,L,K,KFLIP

    !-----------------------------------------------------------------------
    !
    DO I=1,nCols  
       LMH(I) = kMax-LOWLYR(I)+1
    END DO


    DO   I=1,nCols  
       DO L=1,kMax
          KFLIP=kMax+1-L
          CWM    (I,L)=CWM_PHY    (I,KFLIP)
          T      (I,L)=T_PHY      (I,KFLIP)
          Q      (I,L)=Q_PHY      (I,KFLIP)
          P      (I,L)=P_PHY      (I,KFLIP)
          TLATGS (I,L)=TLATGS_PHY (I,KFLIP)
          TRAIN  (I,L)=TRAIN_PHY  (I,KFLIP)
          F_ice  (L,I)=F_ice_PHY  (I,KFLIP)
          F_rain (L,I)=F_rain_PHY (I,KFLIP)
          F_RimeF(L,I)=F_RimeF_PHY(I,KFLIP)
       ENDDO
    END DO

    DO  I=1,nCols  
       LSFC=LMH(I)                      ! "L" of surface
       !
       DO K=1,kMax
          KFLIP=kMax+1-K
          DPCOL(K)=RHO_PHY(I,KFLIP)*GRAV*dz8w(I,KFLIP)
       ENDDO
       !   
       !
       !--- Initialize column data (1D arrays)
       !
       IF (CWM(I,1) .LE. EPSQ) CWM(I,1)=EPSQ
       F_ice  (1,I)=1.0_r8
       F_rain (1,I)=0.0_r8
       F_RimeF(1,I)=1.0_r8
       DO L=1,LSFC
          !
          !--- Pressure (Pa) = (Psfc-Ptop)*(ETA/ETA_sfc)+Ptop
          !
          P_col(L)=P(I,L)
          !
          !--- Layer thickness = RHO*DZ = -DP/G = (Psfc-Ptop)*D_ETA/(G*ETA_sfc)
          !
          THICK_col(L) =DPCOL(L)*RGRAV
          T_col    (L) =T(I,L)
          TC           =T_col(L)-T0C
          QV_col   (L) =MAX(EPSQ, Q(I,L))

          IF (CWM(I,L) .LE. EPSQ1) THEN
             WC_col(L)=0.
             IF (TC .LT. T_ICE) THEN
                F_ice(L,I)=1.0_r8
             ELSE
                F_ice(L,I)=0.0_r8
             ENDIF
             F_rain(L,I)=0.0_r8
             F_RimeF(L,I)=1.0_r8
          ELSE
             WC_col(L)=CWM(I,L)
          ENDIF
          !
          !--- Determine composition of condensate in terms of 
          !      cloud water, ice, & rain
          !
          WC=WC_col(L)
          QI=0.0_r8
          QR=0.0_r8
          QW=0.0_r8
          Fice =F_ice(L,I)
          Frain=F_rain(L,I)
!
!--- REAL*4 array storage
!
          IF (Fice .GE. 1.0_r8) THEN
             QI=WC
          ELSE IF (Fice .LE. 0.0_r8) THEN
             QW=WC
          ELSE
             QI=Fice*WC
             QW=WC-QI
          ENDIF
          IF (QW.GT.0.0_r8 .AND. Frain.GT.0.0_r8) THEN
             IF (Frain .GE. 1.0_r8) THEN
                QR=QW
                QW=0.0_r8
             ELSE
                QR=Frain*QW
                QW=QW-QR
             ENDIF
          ENDIF
          RimeF_col(L)=F_RimeF(L,I)
          QI_col(L)=QI
          QR_col(L)=QR
          QW_col(L)=QW
          !GFDL => New.  Added RHC_col to allow for height- and grid-dependent values for
          !GFDL          the relative humidity threshold for condensation ("RHgrd")
          !6/11/2010 mod - Use lower RHgrd_out threshold for < 850 hPa
          !------------------------------------------------------------
          IF(GID .EQ. 1 .AND. P_col(L)<P_RHgrd_out) THEN  ! gopal's doing based on GFDL
             RHC_col(L)=RHgrd_out        
          ELSE
             RHC_col(L)=RHgrd_in       
          ENDIF
          !------------------------------------------------------------
       ENDDO
       !
       !#######################################################################
       !
       !--- Perform the microphysical calculations in this column
       !
       CALL EGCP01COLUMN ( &
            ARAIN    , & !REAL(KIND=r8),INTENT(INOUT)    :: ARAIN
            ASNOW2    , & !REAL(KIND=r8),INTENT(INOUT)    :: ASNOW2
            DTPH     , & !REAL(KIND=r8),INTENT(IN   )    :: DTPH
            LSFC     , & !INTEGER,INTENT(IN)    :: LSFC
            P_col    , & !REAL(KIND=r8),INTENT(INOUT)    :: P_col    (1:kMax)
            QI_col   , & !REAL(KIND=r8),INTENT(INOUT)    :: QI_col   (1:kMax)
            QR_col   , & !REAL(KIND=r8),INTENT(INOUT)    :: QR_col   (1:kMax)
            QV_col   , & !REAL(KIND=r8),INTENT(INOUT)    :: QV_col   (1:kMax)
            QW_col   , & !REAL(KIND=r8),INTENT(INOUT)    :: QW_col   (1:kMax)
            RimeF_col, & !REAL(KIND=r8),INTENT(INOUT)    :: RimeF_col(1:kMax)
            T_col    , & !REAL(KIND=r8),INTENT(INOUT)    :: T_col    (1:kMax)
            THICK_col, & !REAL(KIND=r8),INTENT(INOUT)    :: THICK_col(1:kMax)
            WC_col   , & !REAL(KIND=r8),INTENT(INOUT)    :: WC_col   (1:kMax)
            RHC_col  , & !REAL(KIND=r8),INTENT(INOUT)    :: RHC_col  (1:kMax)!GFDL
            kMax       ) !INTEGER,INTENT(IN)    :: kMax


       !
       !#######################################################################
       !
       !
       !--- Update storage arrays
       !
       DO L=1,LSFC
          TRAIN (I,L)=(T_col(L)-T(I,L))/DTPH
          TLATGS(I,L)=T_col (L)-T(I,L)
          T     (I,L)=T_col (L)
          Q     (I,L)=QV_col(L)
          CWM   (I,L)=WC_col(L)
          !
          !--- REAL(KIND=r8)*4 array storage
          !
          F_RimeF(L,I)=MAX(1.0_r8, RimeF_col(L))
          IF (QI_col(L) .LE. EPSQ) THEN
             F_ice(L,I)=0.0_r8
             IF (T_col(L) .LT. T_ICEK) F_ice(L,I)=1.0_r8
          ELSE
             F_ice(L,I)=MAX( 0.0_r8, MIN(1.0_r8, QI_col(L)/WC_col(L)) )
          ENDIF
          IF (QR_col(L) .LE. EPSQ) THEN
             DUM=0
          ELSE
             DUM=QR_col(L)/(QR_col(L)+QW_col(L))
          ENDIF
          F_rain(L,I)=DUM
          !
       ENDDO
       !
       !--- Update accumulated precipitation statistics
       !
       !--- Surface precipitation statistics; SR is fraction of surface 
       !    precipitation (if >0) associated with snow
       !
       APREC(I)=(ARAIN )*RRHOL       ! Accumulated surface precip (depth in m)  !<--- Ying
       ASNOW(I)=(ASNOW2)*RRHOL       ! Accumulated surface snow (depth in m)  !<--- Ying

       PREC  (I)= PREC (I)  + (ARAIN+ASNOW2)*RRHOL
       ACPREC(I)=ACPREC(I)  + (ARAIN+ASNOW2)*RRHOL
       IF((ARAIN+ASNOW2)*RRHOL .LT. 1.E-8_r8) THEN
          SR(I)=0.0_r8
       ELSE
          SR(I)=RRHOL*ASNOW2/(ARAIN+ASNOW2)*RRHOL
       ENDIF
       !   !
       !   !--- Debug statistics 
       !   !
       !        IF (PRINT_diag) THEN
       !          PRECtot(1)=PRECtot(1)+ARAIN
       !          PRECtot(2)=PRECtot(2)+ASNOW2
       !          PRECmax(1)=MAX(PRECmax(1), ARAIN)
       !          PRECmax(2)=MAX(PRECmax(2), ASNOW2)
       !        ENDIF
       !#######################################################################
       !#######################################################################
       !
       !100       CONTINUE                          ! End "I" & "J" loops
    END DO

    DO  I=1,nCols  
       DO L=1,kMax
          KFLIP=kMax+1-L
          CWM_PHY    (I,KFLIP)=CWM    (I,L)
          T_PHY      (I,KFLIP)=T      (I,L)
          Q_PHY      (I,KFLIP)=Q      (I,L)
          TLATGS_PHY (I,KFLIP)=TLATGS (I,L)
          TRAIN_PHY  (I,KFLIP)=TRAIN  (I,L)
          F_ice_PHY  (I,KFLIP)=F_ice  (L,I)
          F_rain_PHY (I,KFLIP)=F_rain (L,I)
          F_RimeF_PHY(I,KFLIP)=F_RimeF(L,I)
       ENDDO
    END DO
  END SUBROUTINE EGCP01DRV
  !
  !
  !###############################################################################
  ! ***** VERSION OF MICROPHYSICS DESIGNED FOR HIGHER RESOLUTION MESO ETA MODEL
  !       (1) Represents sedimentation by preserving a portion of the precipitation
  !           through top-down integration from cloud-top.  Modified procedure to
  !           Zhao and Carr (1997).
  !       (2) Microphysical equations are modified to be less sensitive to time
  !           steps by use of Clausius-Clapeyron equation to account for changes in
  !           saturation mixing ratios in response to latent heating/cooling.  
  !       (3) Prevent spurious temperature oscillations across 0C due to 
  !           microphysics.
  !       (4) Uses lookup tables for: calculating two different ventilation 
  !           coefficients in condensation and deposition processes; accretion of
  !           cloud water by precipitation; precipitation mass; precipitation rate
  !           (and mass-weighted precipitation fall speeds).
  !       (5) Assumes temperature-dependent variation in mean diameter of large ice
  !           (Houze et al., 1979; Ryan et al., 1996).
  !        -> 8/22/01: This relationship has been extended to colder temperatures
  !           to parameterize smaller large-ice particles down to mean sizes of MDImin,
  !           which is 50 microns reached at -55.9C.
  !       (6) Attempts to differentiate growth of large and small ice, mainly for
  !           improved transition from thin cirrus to thick, precipitating ice
  !           anvils.
  !        -> 8/22/01: This feature has been diminished by effectively adjusting to
  !           ice saturation during depositional growth at temperatures colder than
  !           -10C.  Ice sublimation is calculated more explicitly.  The logic is
  !           that sources of are either poorly understood (e.g., nucleation for NWP) 
  !           or are not represented in the Eta model (e.g., detrainment of ice from 
  !           convection).  Otherwise the model is too wet compared to the radiosonde
  !           observations based on 1 Feb - 18 March 2001 retrospective runs.  
  !       (7) Top-down integration also attempts to treat mixed-phase processes,
  !           allowing a mixture of ice and water.  Based on numerous observational
  !           studies, ice growth is based on nucleation at cloud top &
  !           subsequent growth by vapor deposition and riming as the ice particles 
  !           fall through the cloud.  Effective nucleation rates are a function
  !           of ice supersaturation following Meyers et al. (JAM, 1992).  
  !        -> 8/22/01: The simulated relative humidities were far too moist compared 
  !           to the rawinsonde observations.  This feature has been substantially 
  !           diminished, limited to a much narrower temperature range of 0 to -10C.  
  !       (8) Depositional growth of newly nucleated ice is calculated for large time
  !           steps using Fig. 8 of Miller and Young (JAS, 1979), at 1 deg intervals
  !           using their ice crystal masses calculated after 600 s of growth in water
  !           saturated conditions.  The growth rates are normalized by time step
  !           assuming 3D growth with time**1.5 following eq. (6.3) in Young (1993).
  !        -> 8/22/01: This feature has been effectively limited to 0 to -10C.  
  !       (9) Ice precipitation rates can increase due to increase in response to
  !           cloud water riming due to (a) increased density & mass of the rimed
  !           ice, and (b) increased fall speeds of rimed ice.
  !        -> 8/22/01: This feature has been effectively limited to 0 to -10C.  
  !###############################################################################
  !###############################################################################
  !
  SUBROUTINE EGCP01COLUMN ( &
       ARAIN    , &!REAL(KIND=r8),INTENT(INOUT)    :: ARAIN
       ASNOW2   , &!REAL(KIND=r8),INTENT(INOUT)    :: ASNOW2
       DTPH     , &!REAL(KIND=r8),INTENT(IN   )    :: DTPH
       LSFC     , &!INTEGER,INTENT(IN)    :: LSFC
       P_col    , &!REAL(KIND=r8),INTENT(INOUT)    :: P_col    (1:kMax)
       QI_col   , &!REAL(KIND=r8),INTENT(INOUT)    :: QI_col   (1:kMax)
       QR_col   , &!REAL(KIND=r8),INTENT(INOUT)    :: QR_col   (1:kMax)
       QV_col   , &!REAL(KIND=r8),INTENT(INOUT)    :: QV_col   (1:kMax)
       QW_col   , &!REAL(KIND=r8),INTENT(INOUT)    :: QW_col   (1:kMax)
       RimeF_col, &!REAL(KIND=r8),INTENT(INOUT)    :: RimeF_col(1:kMax)
       T_col    , &!REAL(KIND=r8),INTENT(INOUT)    :: T_col    (1:kMax)
       THICK_col, &!REAL(KIND=r8),INTENT(INOUT)    :: THICK_col(1:kMax)
       WC_col   , &!REAL(KIND=r8),INTENT(INOUT)    :: WC_col   (1:kMax)
       RHC_col  , &!REAL(KIND=r8),INTENT(INOUT)    :: RHC_col  (1:kMax)!GFDL
       kMax       )!INTEGER,INTENT(IN)    :: kMax
    !
    !###############################################################################
    !###############################################################################
    !
    !-------------------------------------------------------------------------------
    !----- NOTE:  Code is currently set up w/o threading!  
    !-------------------------------------------------------------------------------
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !                .      .    .     
    ! SUBPROGRAM:  Grid-scale microphysical processes - condensation & precipitation
    !   PRGRMMR: Ferrier         ORG: W/NP22     DATE: 08-2001
    !   PRGRMMR: Jin  (Modification for WRF structure)
    !-------------------------------------------------------------------------------
    ! ABSTRACT:
    !   * Merges original GSCOND & PRECPD subroutines.   
    !   * Code has been substantially streamlined and restructured.
    !   * Exchange between water vapor & small cloud condensate is calculated using
    !     the original Asai (1965, J. Japan) algorithm.  See also references to
    !     Yau and Austin (1979, JAS), Rutledge and Hobbs (1983, JAS), and Tao et al.
    !     (1989, MWR).  This algorithm replaces the Sundqvist et al. (1989, MWR)
    !     parameterization.  
    !-------------------------------------------------------------------------------
    !     
    ! USAGE: 
    !   * CALL EGCP01COLUMN FROM SUBROUTINE EGCP01DRV
    !
    ! INPUT ARGUMENT LIST:
    !   DTPH       - physics time step (s)
    !   I_index    - I index
    !   J_index    - J index
    !   LSFC       - Eta level of level above surface, ground
    !   P_col      - vertical column of model pressure (Pa)
    !   QI_col     - vertical column of model ice mixing ratio (kg/kg)
    !   QR_col     - vertical column of model rain ratio (kg/kg)
    !   QV_col     - vertical column of model water vapor specific humidity (kg/kg)
    !   QW_col     - vertical column of model cloud water mixing ratio (kg/kg)
    !   RimeF_col  - vertical column of rime factor for ice in model (ratio, defined below)
    !   T_col      - vertical column of model temperature (deg K)
    !   THICK_col  - vertical column of model mass thickness (density*height increment)
    !   WC_col     - vertical column of model mixing ratio of total condensate (kg/kg)
    !   RHC_col    - vertical column of threshold relative humidity for onset of condensation (ratio)   !GFDL
    !   
    !
    ! OUTPUT ARGUMENT LIST: 
    !   ARAIN      - accumulated rainfall at the surface (kg)
    !   ASNOW2      - accumulated snowfall at the surface (kg)
    !   QV_col     - vertical column of model water vapor specific humidity (kg/kg)
    !   WC_col     - vertical column of model mixing ratio of total condensate (kg/kg)
    !   QW_col     - vertical column of model cloud water mixing ratio (kg/kg)
    !   QI_col     - vertical column of model ice mixing ratio (kg/kg)
    !   QR_col     - vertical column of model rain ratio (kg/kg)
    !   RimeF_col  - vertical column of rime factor for ice in model (ratio, defined below)
    !   T_col      - vertical column of model temperature (deg K)
    !     
    ! OUTPUT FILES:
    !     NONE
    !     
    ! Subprograms & Functions called:
    !   * REAL(KIND=r8) Function CONDENSE  - cloud water condensation
    !   * REAL(KIND=r8) Function DEPOSIT   - ice deposition (not sublimation)
    !
    ! UNIQUE: NONE
    !  
    ! LIBRARY: NONE
    !  
    ! COMMON BLOCKS:  
    !     CMICRO_CONS  - key constants initialized in GSMCONST
    !     CMICRO_STATS - accumulated and maximum statistics
    !     CMY_GROWTH   - lookup table for growth of ice crystals in 
    !                    water saturated conditions (Miller & Young, 1979)
    !     IVENT_TABLES - lookup tables for ventilation effects of ice
    !     IACCR_TABLES - lookup tables for accretion rates of ice
    !     IMASS_TABLES - lookup tables for mass content of ice
    !     IRATE_TABLES - lookup tables for precipitation rates of ice
    !     IRIME_TABLES - lookup tables for increase in fall speed of rimed ice
    !     RVENT_TABLES - lookup tables for ventilation effects of rain
    !     RACCR_TABLES - lookup tables for accretion rates of rain
    !     RMASS_TABLES - lookup tables for mass content of rain
    !     RVELR_TABLES - lookup tables for fall speeds of rain
    !     RRATE_TABLES - lookup tables for precipitation rates of rain
    !   
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !   MACHINE : IBM SP
    !
    !
    !------------------------------------------------------------------------- 
    !--------------- Arrays & constants in argument list --------------------- 
    !------------------------------------------------------------------------- 
    !
    IMPLICIT NONE
    !    
    INTEGER,INTENT(IN) :: kMax
    INTEGER,INTENT(IN) :: LSFC
    REAL(KIND=r8),INTENT(INOUT) :: ARAIN
    REAL(KIND=r8),INTENT(INOUT) :: ASNOW2
    REAL(KIND=r8),INTENT(IN   ) :: DTPH
    REAL(KIND=r8),INTENT(INOUT) ::  P_col    (1:kMax)
    REAL(KIND=r8),INTENT(INOUT) ::  QI_col   (1:kMax)
    REAL(KIND=r8),INTENT(INOUT) ::  QR_col   (1:kMax)
    REAL(KIND=r8),INTENT(INOUT) ::  QV_col   (1:kMax)
    REAL(KIND=r8),INTENT(INOUT) ::  QW_col   (1:kMax)
    REAL(KIND=r8),INTENT(INOUT) ::  RimeF_col(1:kMax)
    REAL(KIND=r8),INTENT(INOUT) ::  T_col    (1:kMax)
    REAL(KIND=r8),INTENT(INOUT) ::  THICK_col(1:kMax)
    REAL(KIND=r8),INTENT(INOUT) ::  WC_col   (1:kMax)
    REAL(KIND=r8),INTENT(INOUT) ::  RHC_col  (1:kMax)!GFDL
    !
    !------------------------------------------------------------------------- 
    !-------------- Common blocks for microphysical statistics ---------------
    !------------------------------------------------------------------------- 
    !
    !------------------------------------------------------------------------- 
    !--------- Common blocks for constants initialized in GSMCONST ----------
    !
    !INTEGER, PARAMETER    :: ITLO=-60, ITHI=40
    !INTEGER,INTENT(INOUT) :: NSTATS(ITLO:ITHI,4)
    !REAL(KIND=r8),INTENT(INOUT)    :: QMAX  (ITLO:ITHI,5)
    !REAL(KIND=r8),INTENT(INOUT)    :: QTOT  (ITLO:ITHI,22) 
    !
    !------------------------------------------------------------------------- 
    !--------------- Common blocks for various lookup tables -----------------
    !
    !--- Discretized growth rates of small ice crystals after their nucleation 
    !     at 1 C intervals from -1 C to -35 C, based on calculations by Miller 
    !     and Young (1979, JAS) after 600 s of growth.  Resultant growth rates
    !     are multiplied by physics time step in GSMCONST.
    !
    !------------------------------------------------------------------------- 
    !
    !--- Mean ice-particle diameters varying from 50 microns to 1000 microns
    !      (1 mm), assuming an exponential size distribution.  
    !
    !---- Meaning of the following arrays: 
    !        - mdiam - mean diameter (m)
    !        - VENTI1 - integrated quantity associated w/ ventilation effects 
    !                   (capacitance only) for calculating vapor deposition onto ice
    !        - VENTI2 - integrated quantity associated w/ ventilation effects 
    !                   (with fall speed) for calculating vapor deposition onto ice
    !        - ACCRI  - integrated quantity associated w/ cloud water collection by ice
    !        - MASSI  - integrated quantity associated w/ ice mass 
    !        - VSNOWI - mass-weighted fall speed of snow (large ice), used to calculate 
    !                   precipitation rates
    !
    !
    !------------------------------------------------------------------------- 
    !
    !--- VEL_RF - velocity increase of rimed particles as functions of crude
    !      particle size categories (at 0.1 mm intervals of mean ice particle
    !      sizes) and rime factor (different values of Rime Factor of 1.1**N, 
    !      where N=0 to Nrime).
    !
    !------------------------------------------------------------------------- 
    !
    !--- Mean rain drop diameters varying from 50 microns (0.05 mm) to 450 microns 
    !      (0.45 mm), assuming an exponential size distribution.  
    !
    !------------------------------------------------------------------------- 
    !------- Key parameters, local variables, & important comments ---------
    !-----------------------------------------------------------------------
    !
    !--- TOLER => Tolerance or precision for accumulated precipitation 
    !
    REAL(KIND=r8), PARAMETER :: TOLER=5.E-7_r8, C2=1.0_r8/6.0_r8, RHO0=1.194_r8, Xratio=.025_r8
    !
    !--- If BLEND=1:
    !      precipitation (large) ice amounts are estimated at each level as a 
    !      blend of ice falling from the grid point above and the precip ice
    !      present at the start of the time step (see TOT_ICE below).
    !--- If BLEND=0:
    !      precipitation (large) ice amounts are estimated to be the precip
    !      ice present at the start of the time step.
    !
    !--- Extended to include sedimentation of rain on 2/5/01 
    !
    REAL(KIND=r8), PARAMETER :: BLEND=1.0_r8
    !
    !-----------------------------------------------------------------------
    !--- Local variables
    !-----------------------------------------------------------------------
    !
    REAL(KIND=r8) EMAIRI, N0r, NLICE, NSmICE, RHgrd   !GFDL
    LOGICAL CLEAR, ICE_logical, DBG_logical, RAIN_logical
    INTEGER :: IDR,INDEX_MY,INDEXR,INDEXR1,INDEXS,IPASS,IXRF,    &
         &           IXS,LBEF,L
    !
    REAL(KIND=r8) :: ABI,ABW,AIEVP,ARAINnew,ASNOWnew,BLDTRH,BUDGET,            &
         &        CREVP,DELI,DELR,DELT,DELV,DELW,DENOMF,                    &
         &        DENOMI,DENOMW,DENOMWI,DIDEP,                              &
         &        DIEVP,DIFFUS,DLI,DTRHO,DUM,DUM1,                     &
         &        DUM2,DWV0,DWVI,DWVR,DYNVIS,ESI,ESW,FIR,FLARGE,FLIMASS,    &
         &        FSMALL,FWR,FWS,GAMMAR,GAMMAS,                             &
         &        PCOND,PIACR,PIACW,PIACWI,PIACWR,PICND,PIDEP,PIDEP_max,    &
         &        PIEVP,PILOSS,PIMLT,PP,PRACW,PRAUT,PREVP,PRLOSS,           &
         &        QI,QInew,QLICE,QR,QRnew,QSI,QSIgrd,QSInew,QSW,QSW0,       &
         &        QSWgrd,QSWnew,QT,QTICE,QTnew,QTRAIN,QV,QW,QW0,QWnew,      &
         &        RFACTOR,RHO,RIMEF,RIMEF1,RQR,RR,RRHO,SFACTOR,             &
         &        TC,TCC,TFACTOR,THERM_COND,THICK,TK,TK2,TNEW,              &
         &        TOT_ICE,TOT_ICEnew,TOT_RAIN,TOT_RAINnew,                  &
         &        VEL_INC,VENTR,VENTIL,VENTIS,VRAIN1,VRAIN2,VRIMEF,VSNOW,   &
         &        WC,WCnew,WSgrd,WS,WSnew,WV,WVnew,WVQW,                    &
                                !!     &        XLF,XLF1,XLI,XLV,XLV1,XLV2,XLIMASS,XRF,XSIMASS          
         &        XLF,XLF1,XLI,XLV,XLV1,XLV2,XLIMASS,XRF,XSIMASS,           &
         &        VRabove     !-- new variables          
    !
    !#######################################################################
    !########################## Begin Execution ############################
    !#######################################################################
    !
    !
    ARAIN=0.0_r8                ! Accumulated rainfall into grid box from above (kg/m**2)
    ASNOW2=0.0_r8                ! Accumulated snowfall into grid box from above (kg/m**2)
    !aligo
    VRabove=0.0_r8              ! Fall speed of rain into grid box from above (m/s)
    !aligo
    !
    !-----------------------------------------------------------------------
    !------------ Loop from top (L=1) to surface (L=LSFC) ------------------
    !-----------------------------------------------------------------------
    !

    DO L=1,LSFC

       !--- Skip this level and go to the next lower level if no condensate 
       !      and very low specific humidities
       !
       IF (QV_col(L).LE.EPSQ .AND. WC_col(L).LE.EPSQ) GO TO 10
       !
       !-----------------------------------------------------------------------
       !------------ Proceed with cloud microphysics calculations -------------
       !-----------------------------------------------------------------------
       !
       TK=T_col(L)         ! Temperature (deg K)
       TC=TK-T0C           ! Temperature (deg C)
       PP=P_col(L)         ! Pressure (Pa)
       QV=QV_col(L)        ! Specific humidity of water vapor (kg/kg)
       WV=QV/(1.0_r8-QV)   ! Water vapor mixing ratio (kg/kg)
       WC=WC_col(L)        ! Grid-scale mixing ratio of total condensate (water or ice; kg/kg)
       RHgrd=RHC_col(L)    ! Threshold relative humidity for the onset of condensation
       !
       !-----------------------------------------------------------------------
       !--- Moisture variables below are mixing ratios & not specifc humidities
       !-----------------------------------------------------------------------
       !
       CLEAR=.TRUE.
       !    
       !--- This check is to determine grid-scale saturation when no condensate is present
       !    
       ESW=1000.0_r8*FPVS0(TK)              ! Saturation vapor pressure w/r/t water
       QSW=EPS*ESW/(PP-ESW)             ! Saturation mixing ratio w/r/t water
       QSI = QSW                        ! Initialize variable
       WS=QSW                           ! General saturation mixing ratio (water/ice)
       IF (TC .LT. 0.0_r8) THEN
          ESI=1000.0_r8*FPVS(TK)             ! Saturation vapor pressure w/r/t ice
          QSI=EPS*ESI/(PP-ESI)           ! Saturation mixing ratio w/r/t water
          WS=QSI                         ! General saturation mixing ratio (water/ice)
       ENDIF
       !
       !--- Effective grid-scale Saturation mixing ratios
       !
       QSWgrd=RHgrd*QSW
       QSIgrd=RHgrd*QSI
       WSgrd=RHgrd*WS
       !
       !--- Check if air is subsaturated and w/o condensate
       !
       IF (WV.GT.WSgrd .OR. WC.GT.EPSQ) CLEAR=.FALSE.
       !
       !--- Check if any rain is falling into layer from above
       !
       IF (ARAIN .GT. CLIMIT) THEN
          CLEAR=.FALSE.
       ELSE
          ARAIN=0.0_r8
          !aligo
          VRabove=0.0_r8
          !aligo
       ENDIF
       !
       !--- Check if any ice is falling into layer from above
       !
       !--- NOTE that "SNOW" in variable names is synonomous with 
       !    large, precipitation ice particles
       !
       IF (ASNOW2 .GT. CLIMIT) THEN
          CLEAR=.FALSE.
       ELSE
          ASNOW2=0.0_r8
       ENDIF
       !
       !-----------------------------------------------------------------------
       !-- Loop to the end if in clear, subsaturated air free of condensate ---
       !-----------------------------------------------------------------------
       !
       IF (CLEAR) GO TO 10
       !
       !-----------------------------------------------------------------------
       !--------- Initialize RHO, THICK & microphysical processes -------------
       !-----------------------------------------------------------------------
       !
       !
       !--- Virtual temperature, TV=T*(1./EPS-1)*Q, Q is specific humidity;
       !    (see pp. 63-65 in Fleagle & Businger, 1963)
       !
       RHO=PP/(RD*TK*(1.0_r8+EPS1*QV))   ! Air density (kg/m**3)
       RRHO=1.0_r8/RHO                ! Reciprocal of air density
       DTRHO=DTPH*RHO             ! Time step * air density
       BLDTRH=BLEND*DTRHO         ! Blend parameter * time step * air density
       THICK=THICK_col(L)         ! Layer thickness = RHO*DZ = -DP/G = (Psfc-Ptop)*D_ETA/(G*ETA_sfc)
       !
       ARAINnew=0.0_r8                ! Updated accumulated rainfall
       ASNOWnew=0.0_r8                ! Updated accumulated snowfall
       QI=QI_col(L)               ! Ice mixing ratio
       QInew=0.0_r8                   ! Updated ice mixing ratio
       QR=QR_col(L)               ! Rain mixing ratio
       QRnew=0.0_r8                   ! Updated rain ratio
       QW=QW_col(L)               ! Cloud water mixing ratio
       QWnew=0.0_r8                   ! Updated cloud water ratio
       !
       PCOND=0.0_r8                   ! Condensation (>0) or evaporation (<0) of cloud water (kg/kg)
       PIDEP=0.0_r8                   ! Deposition (>0) or sublimation (<0) of ice crystals (kg/kg)
       PIACW=0.0_r8                   ! Cloud water collection (riming) by precipitation ice (kg/kg; >0)
       PIACWI=0.0_r8                  ! Growth of precip ice by riming (kg/kg; >0)
       PIACWR=0.0_r8                  ! Shedding of accreted cloud water to form rain (kg/kg; >0)
       PIACR=0.0_r8                   ! Freezing of rain onto large ice at supercooled temps (kg/kg; >0)
       PICND=0.0_r8                   ! Condensation (>0) onto wet, melting ice (kg/kg)
       PIEVP=0.0_r8                   ! Evaporation (<0) from wet, melting ice (kg/kg)
       PIMLT=0.0_r8                   ! Melting ice (kg/kg; >0)
       PRAUT=0.0_r8                   ! Cloud water autoconversion to rain (kg/kg; >0)
       PRACW=0.0_r8                   ! Cloud water collection (accretion) by rain (kg/kg; >0)
       PREVP=0.0_r8                   ! Rain evaporation (kg/kg; <0)
       !
       !--- Double check input hydrometeor mixing ratios
       !
       !          DUM=WC-(QI+QW+QR)
       !          DUM1=ABS(DUM)
       !          DUM2=TOLER*MIN(WC, QI+QW+QR)
       !          IF (DUM1 .GT. DUM2) THEN
       !            WRITE(6,"(/2(a,i4),a,i2)") '{@ i=',I_index,' j=',J_index,
       !     &                                 ' L=',L
       !            WRITE(6,"(4(a12,g11.4,1x))") 
       !     & '{@ TCold=',TC,'P=',.01*PP,'DIFF=',DUM,'WCold=',WC,
       !     & '{@ QIold=',QI,'QWold=',QW,'QRold=',QR
       !          ENDIF
       !
       !***********************************************************************
       !*********** MAIN MICROPHYSICS CALCULATIONS NOW FOLLOW! ****************
       !***********************************************************************
       !
       !--- Calculate a few variables, which are used more than once below
       !
       !--- Latent heat of vaporization as a function of temperature from
       !      Bolton (1980, JAS)
       !
       XLV=3.148E6_r8 -2370_r8*TK        ! Latent heat of vaporization (Lv)
       XLF=XLS-XLV                ! Latent heat of fusion (Lf)
       XLV1=XLV*RCP               ! Lv/Cp
       XLF1=XLF*RCP               ! Lf/Cp
       TK2=1.0_r8/(TK*TK)             ! 1./TK**2
       !GFDL          XLV2=XLV*XLV*QSW*TK2/RV    ! Lv**2*Qsw/(Rv*TK**2)
       XLV2=XLV*XLV*QSWgrd*TK2/RV    ! Lv**2*QSWgrd/(Rv*TK**2)   !GFDL
       DENOMW=1.0_r8+XLV2*RCP         ! Denominator term, Clausius-Clapeyron correction
       !
       !--- Basic thermodynamic quantities
       !      * DYNVIS - dynamic viscosity  [ kg/(m*s) ]
       !      * THERM_COND - thermal conductivity  [ J/(m*s*K) ]
       !      * DIFFUS - diffusivity of water vapor  [ m**2/s ]
       !
       TFACTOR=TK**1.5_r8/(TK+120.0_r8)
       DYNVIS=1.496E-6_r8*TFACTOR
       THERM_COND=2.116E-3_r8*TFACTOR
       DIFFUS=8.794E-5_r8*TK**1.81_r8/PP
       !
       !--- Air resistance term for the fall speed of ice following the
       !      basic research by Heymsfield, Kajikawa, others 
       !
       GAMMAS=(1.E5_r8/PP)**C1
       !
       !--- Air resistance for rain fall speed (Beard, 1985, JAS, p.470)
       !
       GAMMAR=(RHO0/RHO)**0.4_r8
       !
       !----------------------------------------------------------------------
       !-------------  IMPORTANT MICROPHYSICS DECISION TREE  -----------------
       !----------------------------------------------------------------------
       !
       !--- Determine if conditions supporting ice are present
       !
       IF (TC.LT.0.0_r8 .OR. QI.GT.EPSQ .OR. ASNOW2.GT.CLIMIT) THEN
          ICE_logical=.TRUE.
       ELSE
          ICE_logical=.FALSE.
          QLICE=0.0_r8
          QTICE=0.0_r8
       ENDIF
       !
       !--- Determine if rain is present
       !
       RAIN_logical=.FALSE.
       IF (ARAIN.GT.CLIMIT .OR. QR.GT.EPSQ) RAIN_logical=.TRUE.
       !
       IF (ICE_logical) THEN
          !
          !--- IMPORTANT:  Estimate time-averaged properties.
          !
          !---
          !  * FLARGE  - ratio of number of large ice to total (large & small) ice
          !  * FSMALL  - ratio of number of small ice crystals to large ice particles
          !  ->  Small ice particles are assumed to have a mean diameter of 50 microns.
          !  * XSIMASS - used for calculating small ice mixing ratio
          !---
          !  * TOT_ICE - total mass (small & large) ice before microphysics,
          !              which is the sum of the total mass of large ice in the 
          !              current layer and the input flux of ice from above
          !  * PILOSS  - greatest loss (<0) of total (small & large) ice by
          !              sublimation, removing all of the ice falling from above
          !              and the ice within the layer
          !  * RimeF1  - Rime Factor, which is the mass ratio of total (unrimed & rimed) 
          !              ice mass to the unrimed ice mass (>=1)
          !  * VrimeF  - the velocity increase due to rime factor or melting (ratio, >=1)
          !  * VSNOW   - Fall speed of rimed snow w/ air resistance correction
          !  * EMAIRI  - equivalent mass of air associated layer and with fall of snow into layer
          !  * XLIMASS - used for calculating large ice mixing ratio
          !  * FLIMASS - mass fraction of large ice
          !  * QTICE   - time-averaged mixing ratio of total ice
          !  * QLICE   - time-averaged mixing ratio of large ice
          !  * NLICE   - time-averaged number concentration of large ice
          !  * NSmICE  - number concentration of small ice crystals at current level
          !---
          !--- Assumed number fraction of large ice particles to total (large & small) 
          !    ice particles, which is based on a general impression of the literature.
          !
          WVQW=WV+QW                ! Water vapor & cloud water
          !


          IF (TC.GE.0.0_r8 .OR. WVQW.LT.QSIgrd) THEN
             !
             !--- Eliminate small ice particle contributions for melting & sublimation
             !
             FLARGE=FLARGE1
          ELSE
             !
             !--- Enhanced number of small ice particles during depositional growth
             !    (effective only when 0C > T >= T_ice [-10C] )
             !
             FLARGE=FLARGE2
             !
             !--- Larger number of small ice particles due to rime splintering
             !
             IF (TC.GE.-8.0_r8 .AND. TC.LE.-3.0_r8) FLARGE=0.5_r8*FLARGE
             !
          ENDIF            ! End IF (TC.GE.0. .OR. WVQW.LT.QSIgrd)
          !GFDL  => turned on in GFDL code, but not here =>           FLARGE=1.0
          FSMALL=(1.0_r8-FLARGE)/FLARGE
          XSIMASS=RRHO*MASSI(MDImin)*FSMALL
          IF (QI.LE.EPSQ .AND. ASNOW2.LE.CLIMIT) THEN
             INDEXS=MDImin
             TOT_ICE=0.0_r8
             PILOSS=0.0_r8
             RimeF1=1.0_r8
             VrimeF=1.0_r8
             VEL_INC=GAMMAS
             VSNOW=0.0_r8
             EMAIRI=THICK
             XLIMASS=RRHO*RimeF1*MASSI(INDEXS)
             FLIMASS=XLIMASS/(XLIMASS+XSIMASS)
             QLICE=0.0_r8
             QTICE=0.0_r8
             NLICE=0.0_r8
             NSmICE=0.0_r8
          ELSE
             !
             !--- For T<0C mean particle size follows Houze et al. (JAS, 1979, p. 160), 
             !    converted from Fig. 5 plot of LAMDAs.  Similar set of relationships 
             !    also shown in Fig. 8 of Ryan (BAMS, 1996, p. 66).
             !
             DUM=XMImax*EXP(0.0536_r8*TC)
             INDEXS=MIN(MDImax, MAX(MDImin, INT(DUM) ) )
             TOT_ICE=THICK*QI+BLEND*ASNOW2
             PILOSS=-TOT_ICE/THICK
             LBEF=MAX(1,L-1)
             DUM1=RimeF_col(LBEF)
             DUM2=RimeF_col(L)
             RimeF1=(DUM2*THICK*QI+DUM1*BLEND*ASNOW2)/TOT_ICE
             RimeF1=MIN(RimeF1, RFmax)
             DO IPASS=0,1
                IF (RimeF1 .LE. 1.0_r8) THEN
                   RimeF1=1.0_r8
                   VrimeF=1.0_r8
                ELSE
                   IXS=MAX(2, MIN(INDEXS/100, 9))
                   XRF=10.492_r8*LOG(RimeF1)
                   IXRF=MAX(0, MIN(INT(XRF), Nrime))
                   IF (IXRF .GE. Nrime) THEN
                      VrimeF=VEL_RF(IXS,Nrime)
                   ELSE
                      VrimeF=VEL_RF(IXS,IXRF)+(XRF-FLOAT(IXRF))*          &
                           &                    (VEL_RF(IXS,IXRF+1)-VEL_RF(IXS,IXRF))
                   ENDIF
                ENDIF            ! End IF (RimeF1 .LE. 1.)
                VEL_INC=GAMMAS*VrimeF
                VSNOW=VEL_INC*VSNOWI(INDEXS)
                !! Based on Alig's email, added on 2012-02-09
                !aligo
                IF (TC>0.0_r8) THEN
                   !-- Idea provided by Greg Thompson to smoothly increase the fall speed
                   !   of melting snow
                   DUM=MAX(VSNOW,VRabove)
                   VEL_INC=DUM/VSNOWI(INDEXS)
                   VSNOW=DUM
                ENDIF
                !aligo

                !!
                EMAIRI=THICK+BLDTRH*VSNOW
                XLIMASS=RRHO*RimeF1*MASSI(INDEXS)
                FLIMASS=XLIMASS/(XLIMASS+XSIMASS)
                QTICE=TOT_ICE/EMAIRI
                QLICE=FLIMASS*QTICE
                NLICE=QLICE/XLIMASS
                NSmICE=Fsmall*NLICE
                !
                IF ( (NLICE.GE.NLImin .AND. NLICE.LE.NLImax)            &
                     &                .OR. IPASS.EQ.1) THEN
                   EXIT
                ELSE
                   !
                   !--- Reduce excessive accumulation of ice at upper levels
                   !    associated with strong grid-resolved ascent
                   !
                   !--- Force NLICE to be between NLImin and NLImax
                   !
                   DUM=MAX(NLImin, MIN(NLImax, NLICE) )
                   !      XLI=RHO*(QTICE/DUM-XSIMASS)/RimeF1
!!! Bug fix 
!!! bug fix 2012-02-08, see Aligo's email
                   XLI=RHO*QLICE/(DUM*RimeF1)   !- QLICE is for large ice 
                   IF (XLI .LE. MASSI(MDImin) ) THEN
                      INDEXS=MDImin
                   ELSE IF (XLI .LE. MASSI(450) ) THEN
                      DLI=9.5885E5_r8*XLI**0.42066_r8         ! DLI in microns
                      INDEXS=MIN(MDImax, MAX(MDImin, INT(DLI) ) )
                   ELSE IF (XLI .LE. MASSI(MDImax) ) THEN
                      DLI=3.9751E6_r8*XLI**0.49870_r8         ! DLI in microns
                      INDEXS=MIN(MDImax, MAX(MDImin, INT(DLI) ) )
                   ELSE 
                      INDEXS=MDImax
                      !
                      !--- 8/22/01: Increase density of large ice if maximum limits 
                      !    are reached for number concentration (NLImax) and mean size 
                      !    (MDImax).  Done to increase fall out of ice.
                      !
                      IF (DUM .GE. NLImax)                                &
!!!     &                RimeF1=RHO*(QTICE/NLImax-XSIMASS)/MASSI(INDEXS)
!!! bug fix 2012-02-08, see Aligo's email
                           &                   RimeF1=RHO*QLICE/(NLImax*MASSI(INDEXS))    !- QLICE is for large ice 
                   ENDIF             ! End IF (XLI .LE. MASSI(MDImin) ) 
                   !            WRITE(6,"(4(a12,g11.4,1x))") 
                   !     & '{$ TC=',TC,'P=',.01*PP,'NLICE=',NLICE,'DUM=',DUM,
                   !     & '{$ XLI=',XLI,'INDEXS=',FLOAT(INDEXS),'RHO=',RHO,'QTICE=',QTICE,
                   !     & '{$ XSIMASS=',XSIMASS,'RimeF1=',RimeF1
                ENDIF                  ! End IF ( (NLICE.GE.NLImin .AND. NLICE.LE.NLImax) ...
             ENDDO                    ! End DO IPASS=0,1
          ENDIF                      ! End IF (QI.LE.EPSQ .AND. ASNOW2.LE.CLIMIT)
       ENDIF                        ! End IF (ICE_logical)
       !
       !----------------------------------------------------------------------
       !--------------- Calculate individual processes -----------------------
       !----------------------------------------------------------------------
       !
       !--- Cloud water autoconversion to rain and collection by rain
       !
       IF (QW.GT.EPSQ .AND. TC.GE.T_ICE) THEN
          !
          !--- QW0 could be modified based on land/sea properties, 
          !      presence of convection, etc.  This is why QAUT0 and CRAUT
          !      are passed into the subroutine as externally determined
          !      parameters.  Can be changed in the future if desired.
          !
          QW0=QAUT0*RRHO
          PRAUT=MAX(0.0_r8, QW-QW0)*CRAUT
          IF (QLICE .GT. EPSQ) THEN
             !
             !--- Collection of cloud water by large ice particles ("snow")
             !    PIACWI=PIACW for riming, PIACWI=0 for shedding
             !
             FWS=MIN(1.0_r8, CIACW*VEL_INC*NLICE*ACCRI(INDEXS)/PP**C1)
             PIACW=FWS*QW
             IF (TC .LT. 0.0_r8) PIACWI=PIACW    ! Large ice riming
          ENDIF           ! End IF (QLICE .GT. EPSQ)
       ENDIF             ! End IF (QW.GT.EPSQ .AND. TC.GE.T_ICE)
       !
       !----------------------------------------------------------------------
       !--- Loop around some of the ice-phase processes if no ice should be present
       !----------------------------------------------------------------------
       !
       IF (ICE_logical .EQV. .FALSE.) GO TO 20
       !
       !--- Now the pretzel logic of calculating ice deposition
       !
       IF (TC.LT.T_ICE .AND. (WV.GT.QSIgrd .OR. QW.GT.EPSQ)) THEN
          !
          !--- Adjust to ice saturation at T<T_ICE (-10C) if supersaturated.
          !    Sources of ice due to nucleation and convective detrainment are
          !    either poorly understood, poorly resolved at typical NWP 
          !    resolutions, or are not represented (e.g., no detrained 
          !    condensate in BMJ Cu scheme).
          !    
          PCOND=-QW
          DUM1=TK+XLV1*PCOND                 ! Updated (dummy) temperature (deg K)
          DUM2=WV+QW                         ! Updated (dummy) water vapor mixing ratio
          DUM=1000.0_r8*FPVS(DUM1)               ! Updated (dummy) saturation vapor pressure w/r/t ice
          DUM=RHgrd*EPS*DUM/(PP-DUM)         ! Updated (dummy) saturation mixing ratio w/r/t ice
          IF (DUM2 .GT. DUM) PIDEP=DEPOSIT (PP, DUM1, DUM2, RHgrd)  !GFDL
          DWVi=0.0_r8    ! Used only for debugging
          !
       ELSE IF (TC .LT. 0.0_r8) THEN
          !
          !--- These quantities are handy for ice deposition/sublimation
          !    PIDEP_max - max deposition or minimum sublimation to ice saturation
          !
          !GFDL            DENOMI=1.+XLS2*QSI*TK2
          !GFDL            DWVi=MIN(WVQW,QSW)-QSI
          DENOMI=1.0_r8+XLS2*QSIgrd*TK2   !GFDL
          DWVi=MIN(WVQW,QSWgrd)-QSIgrd   !GFDL
          PIDEP_max=MAX(PILOSS, DWVi/DENOMI)
          IF (QTICE .GT. 0.0_r8) THEN
             !
             !--- Calculate ice deposition/sublimation
             !      * SFACTOR - [VEL_INC**.5]*[Schmidt**(1./3.)]*[(RHO/DYNVIS)**.5],
             !        where Schmidt (Schmidt Number) =DYNVIS/(RHO*DIFFUS)
             !      * Units: SFACTOR - s**.5/m ;  ABI - m**2/s ;  NLICE - m**-3 ;
             !               VENTIL, VENTIS - m**-2 ;  VENTI1 - m ;  
             !               VENTI2 - m**2/s**.5 ; DIDEP - unitless
             !
             SFACTOR=SQRT(VEL_INC)*(RHO/(DIFFUS*DIFFUS*DYNVIS))**C2   !GFDL
             ABI=1.0_r8/(RHO*XLS3*QSI*TK2/THERM_COND+1.0_r8/DIFFUS)
             !
             !--- VENTIL - Number concentration * ventilation factors for large ice
             !--- VENTIS - Number concentration * ventilation factors for small ice
             !
             !--- Variation in the number concentration of ice with time is not
             !      accounted for in these calculations (could be in the future).
             !
             VENTIL=(VENTI1(INDEXS)+SFACTOR*VENTI2(INDEXS))*NLICE
             VENTIS=(VENTI1(MDImin)+SFACTOR*VENTI2(MDImin))*NSmICE
             DIDEP=ABI*(VENTIL+VENTIS)*DTPH
             !
             !--- Account for change in water vapor supply w/ time
             !
             IF (DIDEP .GE. Xratio)THEN
                DIDEP=(1.0_r8-EXP(-DIDEP*DENOMI))/DENOMI
             ENDIF
             IF (DWVi .GT. 0.0_r8) THEN
                PIDEP=MIN(DWVi*DIDEP, PIDEP_max)
             ELSE IF (DWVi .LT. 0.0_r8) THEN
                PIDEP=MAX(DWVi*DIDEP, PIDEP_max)
             ENDIF
             !
             !GFDL            ELSE IF (WVQW.GT.QSI .AND. TC.LE.T_ICE_init) THEN
          ELSE IF (WVQW.GT.QSIgrd .AND. TC.LE.T_ICE_init) THEN   !GFDL
             !
             !--- Ice nucleation in near water-saturated conditions.  Ice crystal
             !    growth during time step calculated using Miller & Young (1979, JAS).
             !--- These deposition rates could drive conditions below water saturation,
             !    which is the basis of these calculations.  Intended to approximate
             !    more complex & computationally intensive calculations.
             !
             INDEX_MY=MAX(MY_T1, MIN( INT(0.5_r8-TC), MY_T2 ) )
             !
             !--- DUM1 is the supersaturation w/r/t ice at water-saturated conditions
             !
             !--- DUM2 is the number of ice crystals nucleated at water-saturated 
             !    conditions based on Meyers et al. (JAM, 1992).
             !
             !--- Prevent unREAL(KIND=r8)istically large ice initiation (limited by PIDEP_max)
             !      if DUM2 values are increased in future experiments
             !
             DUM1=QSW/QSI-1.0_r8      
             DUM2=1.E3_r8*EXP(12.96_r8*DUM1-0.639_r8)
             PIDEP=MIN(PIDEP_max, DUM2*MY_GROWTH(INDEX_MY)*RRHO)
             !
          ENDIF       ! End IF (QTICE .GT. 0.)
          !
       ENDIF         ! End IF (TC.LT.T_ICE .AND. (WV.GT.QSIgrd .OR. QW.GT.EPSQ))
       !
       !------------------------------------------------------------------------
       !
20     CONTINUE     ! Jump here if conditions for ice are not present


       !
       !------------------------------------------------------------------------
       !
       !--- Cloud water condensation
       !
       IF (TC.GE.T_ICE .AND. (QW.GT.EPSQ .OR. WV.GT.QSWgrd)) THEN
          IF (PIACWI.EQ.0.0_r8 .AND. PIDEP.EQ.0.0_r8) THEN
             PCOND=CONDENSE (PP, QW, TK, WV, RHgrd)   !GFDL
          ELSE
             !
             !--- Modify cloud condensation in response to ice processes
             !
             DUM=XLV*QSWgrd*RCPRV*TK2
             DENOMWI=1.0_r8+XLS*DUM
             DENOMF=XLF*DUM
             DUM=MAX(0.0_r8, PIDEP)
             PCOND=(WV-QSWgrd-DENOMWI*DUM-DENOMF*PIACWI)/DENOMW
             DUM1=-QW
             DUM2=PCOND-PIACW
             IF (DUM2 .LT. DUM1) THEN
                !
                !--- Limit cloud water sinks
                !
                DUM=DUM1/DUM2
                PCOND=DUM*PCOND
                PIACW=DUM*PIACW
                PIACWI=DUM*PIACWI
             ENDIF        ! End IF (DUM2 .LT. DUM1)
          ENDIF          ! End IF (PIACWI.EQ.0. .AND. PIDEP.EQ.0.)
       ENDIF            ! End IF (TC.GE.T_ICE .AND. (QW.GT.EPSQ .OR. WV.GT.QSWgrd))
       !
       !--- Limit freezing of accreted rime to prevent temperature oscillations,
       !    a crude Schumann-Ludlam limit (p. 209 of Young, 1993). 
       !
       TCC=TC+XLV1*PCOND+XLS1*PIDEP+XLF1*PIACWI
       IF (TCC .GT. 0.0_r8) THEN
          PIACWI=0.0_r8
          TCC=TC+XLV1*PCOND+XLS1*PIDEP
       ENDIF
       IF (TC.GT.0.0_r8 .AND. TCC.GT.0.0_r8 .AND. ICE_logical) THEN
          !
          !--- Calculate melting and evaporation/condensation
          !      * Units: SFACTOR - s**.5/m ;  ABI - m**2/s ;  NLICE - m**-3 ;
          !               VENTIL - m**-2 ;  VENTI1 - m ;  
          !               VENTI2 - m**2/s**.5 ; CIEVP - /s
          !
          SFACTOR=SQRT(VEL_INC)*(RHO/(DIFFUS*DIFFUS*DYNVIS))**C2   !GFDL
          VENTIL=NLICE*(VENTI1(INDEXS)+SFACTOR*VENTI2(INDEXS))
          AIEVP=VENTIL*DIFFUS*DTPH
          IF (AIEVP .LT. Xratio) THEN
             DIEVP=AIEVP
          ELSE
             DIEVP=1.0_r8-EXP(-AIEVP)
          ENDIF
          QSW0=EPS*ESW0/(PP-ESW0)
          !GFDL            DWV0=MIN(WV,QSW)-QSW0
          DWV0=MIN(WV,QSWgrd)-QSW0*RHgrd   !GFDL
          DUM=QW+PCOND
          !GFDL            IF (WV.LT.QSW .AND. DUM.LE.EPSQ) THEN
          IF (WV.LT.QSWgrd .AND. DUM.LE.EPSQ) THEN   !GFDL
             !
             !--- Evaporation from melting snow (sink of snow) or shedding
             !    of water condensed onto melting snow (source of rain)
             !
             DUM=DWV0*DIEVP
             PIEVP=MAX( MIN(0.0_r8, DUM), PILOSS)
             PICND=MAX(0.0_r8, DUM)
          ENDIF            ! End IF (WV.LT.QSW .AND. DUM.LE.EPSQ)
          PIMLT=THERM_COND*TCC*VENTIL*RRHO*DTPH/XLF
          !
          !--- Limit melting to prevent temperature oscillations across 0C
          !
          DUM1=MAX( 0.0_r8, (TCC+XLV1*PIEVP)/XLF1 )
          PIMLT=MIN(PIMLT, DUM1)
          !
          !--- Limit loss of snow by melting (>0) and evaporation
          !
          DUM=PIEVP-PIMLT
          IF (DUM .LT. PILOSS) THEN
             DUM1=PILOSS/DUM
             PIMLT=PIMLT*DUM1
             PIEVP=PIEVP*DUM1
          ENDIF           ! End IF (DUM .GT. QTICE)
       ENDIF             ! End IF (TC.GT.0. .AND. TCC.GT.0. .AND. ICE_logical) 
       !
       !--- IMPORTANT:  Estimate time-averaged properties.
       !
       !  * TOT_RAIN - total mass of rain before microphysics, which is the sum of
       !               the total mass of rain in the current layer and the input 
       !               flux of rain from above
       !  * VRAIN1   - fall speed of rain into grid from above (with air resistance correction)
       !  * QTRAIN   - time-averaged mixing ratio of rain (kg/kg)
       !  * PRLOSS   - greatest loss (<0) of rain, removing all rain falling from
       !               above and the rain within the layer
       !  * RQR      - rain content (kg/m**3)
       !  * INDEXR   - mean size of rain drops to the nearest 1 micron in size
       !  * N0r      - intercept of rain size distribution (typically 10**6 m**-4)
       !
       TOT_RAIN=0.0_r8
       VRAIN1=0.0_r8
       QTRAIN=0.0_r8
       PRLOSS=0.0_r8
       RQR=0.0_r8
       N0r=0.0_r8
       INDEXR=MDRmin
       INDEXR1=INDEXR    ! For debugging only
       IF (RAIN_logical) THEN
          IF (ARAIN .LE. 0.0_r8) THEN
             INDEXR=MDRmin
             VRAIN1=0.0_r8
          ELSE
             !
             !--- INDEXR (related to mean diameter) & N0r could be modified 
             !      by land/sea properties, presence of convection, etc.
             !
             !--- Rain rate normalized to a density of 1.194 kg/m**3
             !
             RR=ARAIN/(DTPH*GAMMAR)
             !
             IF (RR .LE. RR_DRmin) THEN
                !
                !--- Assume fixed mean diameter of rain (0.2 mm) for low rain rates, 
                !      instead vary N0r with rain rate
                !
                INDEXR=MDRmin
             ELSE IF (RR .LE. RR_DR1) THEN
                !
                !--- Best fit to mass-weighted fall speeds (V) from rain lookup tables 
                !      for mean diameters (Dr) between 0.05 and 0.10 mm:
                !      V(Dr)=5.6023e4*Dr**1.136, V in m/s and Dr in m
                !      RR = PI*1000.*N0r0*5.6023e4*Dr**(4+1.136) = 1.408e15*Dr**5.136,
                !        RR in kg/(m**2*s)
                !      Dr (m) = 1.123e-3*RR**.1947 -> Dr (microns) = 1.123e3*RR**.1947
                !
                INDEXR=INT( 1.123E3_r8*RR**0.1947_r8 + 0.5_r8 )
                INDEXR=MAX( MDRmin, MIN(INDEXR, MDR1) )
             ELSE IF (RR .LE. RR_DR2) THEN
                !
                !--- Best fit to mass-weighted fall speeds (V) from rain lookup tables 
                !      for mean diameters (Dr) between 0.10 and 0.20 mm:
                !      V(Dr)=1.0867e4*Dr**.958, V in m/s and Dr in m
                !      RR = PI*1000.*N0r0*1.0867e4*Dr**(4+.958) = 2.731e14*Dr**4.958,
                !        RR in kg/(m**2*s)
                !      Dr (m) = 1.225e-3*RR**.2017 -> Dr (microns) = 1.225e3*RR**.2017
                !
                INDEXR=INT( 1.225E3_r8*RR**0.2017_r8 + 0.5_r8 )
                INDEXR=MAX( MDR1, MIN(INDEXR, MDR2) )
             ELSE IF (RR .LE. RR_DR3) THEN
                !
                !--- Best fit to mass-weighted fall speeds (V) from rain lookup tables 
                !      for mean diameters (Dr) between 0.20 and 0.32 mm:
                !      V(Dr)=2831.*Dr**.80, V in m/s and Dr in m
                !      RR = PI*1000.*N0r0*2831.*Dr**(4+.80) = 7.115e13*Dr**4.80, 
                !        RR in kg/(m**2*s)
                !      Dr (m) = 1.3006e-3*RR**.2083 -> Dr (microns) = 1.3006e3*RR**.2083
                !
                INDEXR=INT( 1.3006E3_r8*RR**0.2083_r8 + 0.5_r8 )
                INDEXR=MAX( MDR2, MIN(INDEXR, MDR3) )
             ELSE IF (RR .LE. RR_DRmax) THEN
                !
                !--- Best fit to mass-weighted fall speeds (V) from rain lookup tables 
                !      for mean diameters (Dr) between 0.32 and 0.45 mm:
                !      V(Dr)=944.8*Dr**.6636, V in m/s and Dr in m
                !      RR = PI*1000.*N0r0*944.8*Dr**(4+.6636) = 2.3745e13*Dr**4.6636,
                !        RR in kg/(m**2*s)
                !      Dr (m) = 1.355e-3*RR**.2144 -> Dr (microns) = 1.355e3*RR**.2144
                !
                INDEXR=INT( 1.355E3_r8*RR**0.2144_r8 + 0.5_r8 )
                INDEXR=MAX( MDR3, MIN(INDEXR, MDRmax) )
             ELSE 
                !
                !--- Assume fixed mean diameter of rain (0.45 mm) for high rain rates, 
                !      instead vary N0r with rain rate
                !
                INDEXR=MDRmax
             ENDIF              ! End IF (RR .LE. RR_DRmin) etc. 
             VRAIN1=GAMMAR*VRAIN(INDEXR)
          ENDIF              ! End IF (ARAIN .LE. 0.)
          INDEXR1=INDEXR     ! For debugging only
          TOT_RAIN=THICK*QR+BLEND*ARAIN
          QTRAIN=TOT_RAIN/(THICK+BLDTRH*VRAIN1)
          PRLOSS=-TOT_RAIN/THICK
          RQR=RHO*QTRAIN
          !
          !--- RQR - time-averaged rain content (kg/m**3)
          !
          IF (RQR .LE. RQR_DRmin) THEN
             N0r=MAX(N0rmin, CN0r_DMRmin*RQR)
             INDEXR=MDRmin
          ELSE IF (RQR .GE. RQR_DRmax) THEN
             N0r=CN0r_DMRmax*RQR
             INDEXR=MDRmax
          ELSE
             N0r=N0r0
             INDEXR=INT(MAX( XMRmin, MIN(CN0r0*RQR**0.25_r8, XMRmax) ))
          ENDIF
          !
          IF (TC .LT. T_ICE) THEN
             PIACR=-PRLOSS
          ELSE
             !GFDL              DWVr=WV-PCOND-QSW
             DWVr=WV-PCOND-QSWgrd   !GFDL
             DUM=QW+PCOND
             IF (DWVr.LT.0.0_r8 .AND. DUM.LE.EPSQ) THEN
                !
                !--- Rain evaporation
                !
                !    * RFACTOR - [GAMMAR**.5]*[Schmidt**(1./3.)]*[(RHO/DYNVIS)**.5],
                !        where Schmidt (Schmidt Number) =DYNVIS/(RHO*DIFFUS)
                !
                !    * Units: RFACTOR - s**.5/m ;  ABW - m**2/s ;  VENTR - m**-2 ;  
                !             N0r - m**-4 ;  VENTR1 - m**2 ;  VENTR2 - m**3/s**.5 ;
                !             CREVP - unitless
                !
                RFACTOR=SQRT(GAMMAR)*(RHO/(DIFFUS*DIFFUS*DYNVIS))**C2   !GFDL
                ABW=1.0_r8/(RHO*XLV2/THERM_COND+1.0_r8/DIFFUS)
                !
                !--- Note that VENTR1, VENTR2 lookup tables do not include the 
                !      1/Davg multiplier as in the ice tables
                !
                VENTR=N0r*(VENTR1(INDEXR)+RFACTOR*VENTR2(INDEXR))
                CREVP=ABW*VENTR*DTPH
                IF (CREVP .LT. Xratio) THEN
                   DUM=DWVr*CREVP
                ELSE
                   DUM=DWVr*(1.0_r8-EXP(-CREVP*DENOMW))/DENOMW
                ENDIF
                PREVP=MAX(DUM, PRLOSS)
             ELSE IF (QW .GT. EPSQ) THEN
                FWR=CRACW*GAMMAR*N0r*ACCRR(INDEXR)
                PRACW=MIN(1.0_r8,FWR)*QW
             ENDIF           ! End IF (DWVr.LT.0. .AND. DUM.LE.EPSQ)
             !
             IF (TC.LT.0.0_r8 .AND. TCC.LT.0.0_r8) THEN
                !
                !--- Biggs (1953) heteorogeneous freezing (e.g., Lin et al., 1983)
                !   - Rescaled mean drop diameter from microns (INDEXR) to mm (DUM) to prevent underflow
                !
                DUM=0.001_r8*FLOAT(INDEXR)
                DUM=(EXP(ABFR*TC)-1.0_r8)*DUM*DUM*DUM*DUM*DUM*DUM*DUM
                PIACR=MIN(CBFR*N0r*RRHO*DUM, QTRAIN)
                IF (QLICE .GT. EPSQ) THEN
                   !
                   !--- Freezing of rain by collisions w/ large ice
                   !
                   DUM=GAMMAR*VRAIN(INDEXR)
                   DUM1=DUM-VSNOW
                   !
                   !--- DUM2 - Difference in spectral fall speeds of rain and
                   !      large ice, parameterized following eq. (48) on p. 112 of 
                   !      Murakami (J. Meteor. Soc. Japan, 1990)
                   !
                   DUM2=SQRT(DUM1*DUM1+0.04_r8*DUM*VSNOW)    !GFDL
                   DUM1=5.E-12_r8*INDEXR*INDEXR+2.E-12_r8*INDEXR*INDEXS        &
                        &                 +0.5E-12_r8*INDEXS*INDEXS
                   FIR=MIN(1.0_r8, CIACR*NLICE*DUM1*DUM2)
                   !
                   !--- Future?  Should COLLECTION BY SMALL ICE SHOULD BE INCLUDED???
                   !
                   PIACR=MIN(PIACR+FIR*QTRAIN, QTRAIN)
                ENDIF        ! End IF (QLICE .GT. EPSQ)
                DUM=PREVP-PIACR
                IF (DUM .LT. PRLOSS) THEN
                   DUM1=PRLOSS/DUM
                   PREVP=DUM1*PREVP
                   PIACR=DUM1*PIACR
                ENDIF        ! End If (DUM .LT. PRLOSS)
             ENDIF          ! End IF (TC.LT.0. .AND. TCC.LT.0.)
          ENDIF            ! End IF (TC .LT. T_ICE)
       ENDIF              ! End IF (RAIN_logical) 
       !
       !----------------------------------------------------------------------
       !---------------------- Main Budget Equations -------------------------
       !----------------------------------------------------------------------
       !
       !
       !-----------------------------------------------------------------------
       !--- Update fields, determine characteristics for next lower layer ----
       !-----------------------------------------------------------------------
       !
       !--- Carefully limit sinks of cloud water
       !
       DUM1=PIACW+PRAUT+PRACW-MIN(0.,PCOND)
       IF (DUM1 .GT. QW) THEN
          DUM=QW/DUM1
          PIACW=DUM*PIACW
          PIACWI=DUM*PIACWI
          PRAUT=DUM*PRAUT
          PRACW=DUM*PRACW
          IF (PCOND .LT. 0.0_r8) PCOND=DUM*PCOND
       ENDIF
       PIACWR=PIACW-PIACWI          ! TC >= 0C
       !
       !--- QWnew - updated cloud water mixing ratio
       !
       DELW=PCOND-PIACW-PRAUT-PRACW
       QWnew=QW+DELW
       IF (QWnew .LE. EPSQ) QWnew=0.0_r8
       IF (QW.GT.0.0_r8 .AND. QWnew.NE.0.0_r8) THEN
          DUM=QWnew/QW
          IF (DUM .LT. TOLER) QWnew=0.0_r8
       ENDIF
       !
       !--- Update temperature and water vapor mixing ratios
       !
       DELT= XLV1*(PCOND+PIEVP+PICND+PREVP)                          &
            &         +XLS1*PIDEP+XLF1*(PIACWI+PIACR-PIMLT)
       Tnew=TK+DELT
       !
       DELV=-PCOND-PIDEP-PIEVP-PICND-PREVP
       WVnew=WV+DELV
       !
       !--- Update ice mixing ratios
       !
       !---
       !  * TOT_ICEnew - total mass (small & large) ice after microphysics,
       !                 which is the sum of the total mass of large ice in the 
       !                 current layer and the flux of ice out of the grid box below
       !  * RimeF      - Rime Factor, which is the mass ratio of total (unrimed & 
       !                 rimed) ice mass to the unrimed ice mass (>=1)
       !  * QInew      - updated mixing ratio of total (large & small) ice in layer
       !      -> TOT_ICEnew=QInew*THICK+BLDTRH*QLICEnew*VSNOW
       !        -> But QLICEnew=QInew*FLIMASS, so
       !      -> TOT_ICEnew=QInew*(THICK+BLDTRH*FLIMASS*VSNOW)
       !  * ASNOWnew   - updated accumulation of snow at bottom of grid cell
       !---
       !
       DELI=0.0_r8
       RimeF=1.0_r8
       IF (ICE_logical) THEN
          DELI=PIDEP+PIEVP+PIACWI+PIACR-PIMLT
          TOT_ICEnew=TOT_ICE+THICK*DELI
          IF (TOT_ICE.GT.0.0_r8 .AND. TOT_ICEnew.NE.0.0_r8) THEN
             DUM=TOT_ICEnew/TOT_ICE
             IF (DUM .LT. TOLER) TOT_ICEnew=0.0_r8
          ENDIF
          IF (TOT_ICEnew .LE. CLIMIT) THEN
             TOT_ICEnew=0.0_r8
             RimeF=1.0_r8
             QInew=0.0_r8
             ASNOWnew=0.0_r8
          ELSE
             !
             !--- Update rime factor if appropriate
             !
             DUM=PIACWI+PIACR
             IF (DUM.LE.EPSQ .AND. PIDEP.LE.EPSQ) THEN
                RimeF=RimeF1
             ELSE
                !
                !--- Rime Factor, RimeF = (Total ice mass)/(Total unrimed ice mass)
                !      DUM1 - Total ice mass, rimed & unrimed
                !      DUM2 - Estimated mass of *unrimed* ice
                !
                DUM1=TOT_ICE+THICK*(PIDEP+DUM)
                DUM2=TOT_ICE/RimeF1+THICK*PIDEP
                IF (DUM2 .LE. 0.0_r8) THEN
                   RimeF=RFmax
                ELSE
                   RimeF=MIN(RFmax, MAX(1.0_r8, DUM1/DUM2) )
                ENDIF
             ENDIF       ! End IF (DUM.LE.EPSQ .AND. PIDEP.LE.EPSQ)
             QInew=TOT_ICEnew/(THICK+BLDTRH*FLIMASS*VSNOW)
             IF (QInew .LE. EPSQ) QInew=0.0_r8
             IF (QI.GT.0.0_r8 .AND. QInew.NE.0.0_r8) THEN
                DUM=QInew/QI
                IF (DUM .LT. TOLER) QInew=0.0_r8
             ENDIF
             ASNOWnew=BLDTRH*FLIMASS*VSNOW*QInew
             IF (ASNOW2.GT.0.0_r8 .AND. ASNOWnew.NE.0.0_r8) THEN
                DUM=ASNOWnew/ASNOW2
                IF (DUM .LT. TOLER) ASNOWnew=0.0_r8
             ENDIF
          ENDIF         ! End IF (TOT_ICEnew .LE. CLIMIT)
       ENDIF           ! End IF (ICE_logical)


       !
       !--- Update rain mixing ratios
       !
       !---
       ! * TOT_RAINnew - total mass of rain after microphysics
       !                 current layer and the input flux of ice from above
       ! * VRAIN2      - time-averaged fall speed of rain in grid and below 
       !                 (with air resistance correction)
       ! * QRnew       - updated rain mixing ratio in layer
       !      -> TOT_RAINnew=QRnew*(THICK+BLDTRH*VRAIN2)
       !  * ARAINnew  - updated accumulation of rain at bottom of grid cell
       !---
       !
       DELR=PRAUT+PRACW+PIACWR-PIACR+PIMLT+PREVP+PICND
       TOT_RAINnew=TOT_RAIN+THICK*DELR
       IF (TOT_RAIN.GT.0.0_r8 .AND. TOT_RAINnew.NE.0.0_r8) THEN
          DUM=TOT_RAINnew/TOT_RAIN
          IF (DUM .LT. TOLER) TOT_RAINnew=0.0_r8
       ENDIF
       IF (TOT_RAINnew .LE. CLIMIT) THEN
          TOT_RAINnew=0.0_r8
          VRAIN2=0.0_r8
          QRnew=0.0_r8
          ARAINnew=0.0_r8
       ELSE
          !
          !--- 1st guess time-averaged rain rate at bottom of grid box
          !
          RR=TOT_RAINnew/(DTPH*GAMMAR)
          !
          !--- Use same algorithm as above for calculating mean drop diameter
          !      (IDR, in microns), which is used to estimate the time-averaged
          !      fall speed of rain drops at the bottom of the grid layer.  This
          !      isn't perfect, but the alternative is solving a transcendental 
          !      equation that is numerically inefficient and nasty to program
          !      (coded in earlier versions of GSMCOLUMN prior to 8-22-01).
          !
          IF (RR .LE. RR_DRmin) THEN
             IDR=MDRmin
          ELSE IF (RR .LE. RR_DR1) THEN
             IDR=INT( 1.123E3_r8*RR**0.1947_r8 + 0.5_r8 )
             IDR=MAX( MDRmin, MIN(IDR, MDR1) )
          ELSE IF (RR .LE. RR_DR2) THEN
             IDR=INT( 1.225E3_r8*RR**0.2017_r8 + 0.5_r8 )
             IDR=MAX( MDR1, MIN(IDR, MDR2) )
          ELSE IF (RR .LE. RR_DR3) THEN
             IDR=INT( 1.3006E3_r8*RR**0.2083_r8 + 0.5_r8 )
             IDR=MAX( MDR2, MIN(IDR, MDR3) )
          ELSE IF (RR .LE. RR_DRmax) THEN
             IDR=INT( 1.355E3_r8*RR**0.2144_r8 + 0.5_r8 )
             IDR=MAX( MDR3, MIN(IDR, MDRmax) )
          ELSE 
             IDR=MDRmax
          ENDIF              ! End IF (RR .LE. RR_DRmin)
          VRAIN2=GAMMAR*VRAIN(IDR)
          QRnew=TOT_RAINnew/(THICK+BLDTRH*VRAIN2)
          IF (QRnew .LE. EPSQ) QRnew=0.0_r8
          IF (QR.GT.0.0_r8 .AND. QRnew.NE.0.0_r8) THEN
             DUM=QRnew/QR
             IF (DUM .LT. TOLER) QRnew=0.0_r8
          ENDIF
          ARAINnew=BLDTRH*VRAIN2*QRnew
          IF (ARAIN.GT.0.0_r8 .AND. ARAINnew.NE.0.0_r8) THEN
             DUM=ARAINnew/ARAIN
             IF (DUM .LT. TOLER) ARAINnew=0.0_r8
          ENDIF
       ENDIF
       !
       WCnew=QWnew+QRnew+QInew
       !
       !----------------------------------------------------------------------
       !-------------- Begin debugging & verification ------------------------
       !----------------------------------------------------------------------
       !
       !--- QT, QTnew - total water (vapor & condensate) before & after microphysics, resp.
       !


       QT=THICK*(WV+WC)+ARAIN+ASNOW2
       QTnew=THICK*(WVnew+WCnew)+ARAINnew+ASNOWnew
       BUDGET=QT-QTnew
       !
       !--- Additional check on budget preservation, accounting for truncation effects
       !
       IF (PRINT_err) THEN
          DBG_logical=.FALSE.
          DUM=ABS(BUDGET)
          IF (DUM .GT. TOLER) THEN
             DUM=DUM/MIN(QT, QTnew)
             IF (DUM .GT. TOLER) DBG_logical=.TRUE.
          ENDIF
          !
          !          DUM=(RHgrd+.001)*QSInew
          !          IF ( (QWnew.GT.EPSQ) .OR. QRnew.GT.EPSQ .OR. WVnew.GT.DUM)    &
          !     &        .AND. TC.LT.T_ICE )  DBG_logical=.TRUE.
          !
          !          IF (TC.GT.5. .AND. QInew.GT.EPSQ) DBG_logical=.TRUE.
          !
          IF (WVnew.LT.EPSQ .OR. DBG_logical) THEN
             !
             !             WRITE(6,"(/2(a,i4),2(a,i2))") '{} i=',I_index,' j=',J_index, &
             !                  &                                    ' L=',L,' LSFC=',LSFC
             !
             ESW=1000.0_r8*FPVS0(Tnew)
             QSWnew=EPS*ESW/(PP-ESW)
             IF (TC.LT.0.0_r8 .OR. Tnew .LT. 0.0_r8) THEN
                ESI=1000.0_r8*FPVS(Tnew)
                QSInew=EPS*ESI/(PP-ESI)
             ELSE
                QSI=QSW
                QSInew=QSWnew
             ENDIF
             WSnew=QSInew
             !           WRITE(6,"(4(a12,g11.4,1x))")                                   &
             !    & '{} TCold=',TC,'TCnew=',Tnew-T0C,'P=',.01*PP,'RHO=',RHO,            &
             !    & '{} THICK=',THICK,'RHold=',WV/WS,'RHnew=',WVnew/WSnew,              &
             !    &   'RHgrd=',RHgrd,                                                   &
             !    & '{} RHWold=',WV/QSW,'RHWnew=',WVnew/QSWnew,'RHIold=',WV/QSI,        &
             !    &   'RHInew=',WVnew/QSInew,                                           &
             !    & '{} QSWold=',QSW,'QSWnew=',QSWnew,'QSIold=',QSI,'QSInew=',QSInew,   &
             !    & '{} WSold=',WS,'WSnew=',WSnew,'WVold=',WV,'WVnew=',WVnew,           &
             !    & '{} WCold=',WC,'WCnew=',WCnew,'QWold=',QW,'QWnew=',QWnew,           &
             !    & '{} QIold=',QI,'QInew=',QInew,'QRold=',QR,'QRnew=',QRnew,           &
             !    & '{} ARAINold=',ARAIN,'ARAINnew=',ARAINnew,'ASNOWold=',ASNOW,        &
             !    &   'ASNOWnew=',ASNOWnew,                                             &
             !    & '{} TOT_RAIN=',TOT_RAIN,'TOT_RAINnew=',TOT_RAINnew,                 &
             !    &   'TOT_ICE=',TOT_ICE,'TOT_ICEnew=',TOT_ICEnew,                      &
             !    & '{} BUDGET=',BUDGET,'QTold=',QT,'QTnew=',QTnew                       
             !  !
             !           WRITE(6,"(4(a12,g11.4,1x))")                                   &
             !    & '{} DELT=',DELT,'DELV=',DELV,'DELW=',DELW,'DELI=',DELI,             &
             !    & '{} DELR=',DELR,'PCOND=',PCOND,'PIDEP=',PIDEP,'PIEVP=',PIEVP,       &
             !    & '{} PICND=',PICND,'PREVP=',PREVP,'PRAUT=',PRAUT,'PRACW=',PRACW,     &
             !    & '{} PIACW=',PIACW,'PIACWI=',PIACWI,'PIACWR=',PIACWR,'PIMLT=',       &
             !    &    PIMLT,                                                           &
             !    & '{} PIACR=',PIACR                                                    
             !  !
             !           IF (ICE_logical) WRITE(6,"(4(a12,g11.4,1x))")                  &
             !    & '{} RimeF1=',RimeF1,'GAMMAS=',GAMMAS,'VrimeF=',VrimeF,              &
             !    &   'VSNOW=',VSNOW,                                                   &
             !    & '{} INDEXS=',FLOAT(INDEXS),'FLARGE=',FLARGE,'FSMALL=',FSMALL,       &
             !    &   'FLIMASS=',FLIMASS,                                               &
             !    & '{} XSIMASS=',XSIMASS,'XLIMASS=',XLIMASS,'QLICE=',QLICE,            &
             !    &   'QTICE=',QTICE,                                                   &
             !    & '{} NLICE=',NLICE,'NSmICE=',NSmICE,'PILOSS=',PILOSS,                &
             !    &   'EMAIRI=',EMAIRI,                                                 &
             !    & '{} RimeF=',RimeF                                                    
             !  !
             !           IF (TOT_RAIN.GT.0. .OR. TOT_RAINnew.GT.0.)                     &
             !    &        WRITE(6,"(4(a12,g11.4,1x))")                                 &
             !    & '{} INDEXR1=',FLOAT(INDEXR1),'INDEXR=',FLOAT(INDEXR),               &
             !    &   'GAMMAR=',GAMMAR,'N0r=',N0r,                                      &
             !    & '{} VRAIN1=',VRAIN1,'VRAIN2=',VRAIN2,'QTRAIN=',QTRAIN,'RQR=',RQR,   &
             !    & '{} PRLOSS=',PRLOSS,'VOLR1=',THICK+BLDTRH*VRAIN1,                   &
             !    &   'VOLR2=',THICK+BLDTRH*VRAIN2
             !  !
             !           IF (PRAUT .GT. 0.) WRITE(6,"(a12,g11.4,1x)") '{} QW0=',QW0
             !  !
             !           IF (PRACW .GT. 0.) WRITE(6,"(a12,g11.4,1x)") '{} FWR=',FWR
             !  !
             !           IF (PIACR .GT. 0.) WRITE(6,"(a12,g11.4,1x)") '{} FIR=',FIR
             !  !
             !            DUM=PIMLT+PICND-PREVP-PIEVP
             !            IF (DUM.GT.0. .or. DWVi.NE.0.)                                 &
             !     &        WRITE(6,"(4(a12,g11.4,1x))")                                 &
             !     & '{} TFACTOR=',TFACTOR,'DYNVIS=',DYNVIS,                             &
             !     &   'THERM_CON=',THERM_COND,'DIFFUS=',DIFFUS
             !   !
             !            IF (PREVP .LT. 0.) WRITE(6,"(4(a12,g11.4,1x))")                &
             !     & '{} RFACTOR=',RFACTOR,'ABW=',ABW,'VENTR=',VENTR,'CREVP=',CREVP,     &
             !     & '{} DWVr=',DWVr,'DENOMW=',DENOMW
             !   !
             !            IF (PIDEP.NE.0. .AND. DWVi.NE.0.)                              &
             !     &        WRITE(6,"(4(a12,g11.4,1x))")                                 &
             !     & '{} DWVi=',DWVi,'DENOMI=',DENOMI,'PIDEP_max=',PIDEP_max,            &
             !     &   'SFACTOR=',SFACTOR,                                               &
             !     & '{} ABI=',ABI,'VENTIL=',VENTIL,'VENTIL1=',VENTI1(INDEXS),           &
             !     &   'VENTIL2=',SFACTOR*VENTI2(INDEXS),                                &
             !     & '{} VENTIS=',VENTIS,'DIDEP=',DIDEP
             !   !
             !            IF (PIDEP.GT.0. .AND. PCOND.NE.0.)                             &
             !     &        WRITE(6,"(4(a12,g11.4,1x))")                                 &
             !     & '{} DENOMW=',DENOMW,'DENOMWI=',DENOMWI,'DENOMF=',DENOMF,            &
             !     &    'DUM2=',PCOND-PIACW
             !   !
             !            IF (FWS .GT. 0.) WRITE(6,"(4(a12,g11.4,1x))")                  &
             !     & '{} FWS=',FWS                     
             !   !
             !            DUM=PIMLT+PICND-PIEVP
             !            IF (DUM.GT. 0.) WRITE(6,"(4(a12,g11.4,1x))")                   &
             !     & '{} SFACTOR=',SFACTOR,'VENTIL=',VENTIL,'VENTIL1=',VENTI1(INDEXS),   &
             !     &   'VENTIL2=',SFACTOR*VENTI2(INDEXS),                                &
             !     & '{} AIEVP=',AIEVP,'DIEVP=',DIEVP,'QSW0=',QSW0,'DWV0=',DWV0       
             !
          ENDIF      !-- IF (WVnew.LT.EPSQ .OR. DBG_logical) THEN
       ENDIF      !-- IF (PRINT_err) THEN

       !
       !-----------------------------------------------------------------------
       !--------------- Water budget statistics & maximum values --------------
       !-----------------------------------------------------------------------
       !
       !          IF (PRINT_diag) THEN
       !            ITdx=MAX( ITLO, MIN( INT(Tnew-T0C), ITHI ) )
       !            IF (QInew .GT. EPSQ) NSTATS(ITdx,1)=NSTATS(ITdx,1)+1
       !            IF (QInew.GT.EPSQ  .AND.  QRnew+QWnew.GT.EPSQ)              &
       !     &        NSTATS(ITdx,2)=NSTATS(ITdx,2)+1
       !            IF (QWnew .GT. EPSQ) NSTATS(ITdx,3)=NSTATS(ITdx,3)+1 
       !            IF (QRnew .GT. EPSQ) NSTATS(ITdx,4)=NSTATS(ITdx,4)+1
       !  !
       !            QMAX(ITdx,1)=MAX(QMAX(ITdx,1), QInew)
       !            QMAX(ITdx,2)=MAX(QMAX(ITdx,2), QWnew)
       !            QMAX(ITdx,3)=MAX(QMAX(ITdx,3), QRnew)
       !            QMAX(ITdx,4)=MAX(QMAX(ITdx,4), ASNOWnew)
       !            QMAX(ITdx,5)=MAX(QMAX(ITdx,5), ARAINnew)
       !            QTOT(ITdx,1)=QTOT(ITdx,1)+QInew*THICK
       !            QTOT(ITdx,2)=QTOT(ITdx,2)+QWnew*THICK
       !            QTOT(ITdx,3)=QTOT(ITdx,3)+QRnew*THICK
       !  !
       !            QTOT(ITdx,4)=QTOT(ITdx,4)+PCOND*THICK
       !            QTOT(ITdx,5)=QTOT(ITdx,5)+PICND*THICK
       !            QTOT(ITdx,6)=QTOT(ITdx,6)+PIEVP*THICK
       !            QTOT(ITdx,7)=QTOT(ITdx,7)+PIDEP*THICK
       !            QTOT(ITdx,8)=QTOT(ITdx,8)+PREVP*THICK
       !            QTOT(ITdx,9)=QTOT(ITdx,9)+PRAUT*THICK
       !            QTOT(ITdx,10)=QTOT(ITdx,10)+PRACW*THICK
       !            QTOT(ITdx,11)=QTOT(ITdx,11)+PIMLT*THICK
       !            QTOT(ITdx,12)=QTOT(ITdx,12)+PIACW*THICK
       !            QTOT(ITdx,13)=QTOT(ITdx,13)+PIACWI*THICK
       !            QTOT(ITdx,14)=QTOT(ITdx,14)+PIACWR*THICK
       !            QTOT(ITdx,15)=QTOT(ITdx,15)+PIACR*THICK
       !  !
       !            QTOT(ITdx,16)=QTOT(ITdx,16)+(WVnew-WV)*THICK
       !            QTOT(ITdx,17)=QTOT(ITdx,17)+(QWnew-QW)*THICK
       !            QTOT(ITdx,18)=QTOT(ITdx,18)+(QInew-QI)*THICK
       !            QTOT(ITdx,19)=QTOT(ITdx,19)+(QRnew-QR)*THICK
       !            QTOT(ITdx,20)=QTOT(ITdx,20)+(ARAINnew-ARAIN)
       !            QTOT(ITdx,21)=QTOT(ITdx,21)+(ASNOWnew-ASNOW)
       !            IF (QInew .GT. 0.)                                          &
       !     &        QTOT(ITdx,22)=QTOT(ITdx,22)+QInew*THICK/RimeF
       !  !
       !          ENDIF
       !
       !----------------------------------------------------------------------
       !------------------------- Update arrays ------------------------------
       !----------------------------------------------------------------------
       !


       T_col(L)=Tnew                           ! Updated temperature
       !
       QV_col(L)=MAX(EPSQ, WVnew/(1.0_r8+WVnew))   ! Updated specific humidity
       WC_col(L)=MAX(EPSQ, WCnew)              ! Updated total condensate mixing ratio
       QI_col(L)=MAX(EPSQ, QInew)              ! Updated ice mixing ratio
       QR_col(L)=MAX(EPSQ, QRnew)              ! Updated rain mixing ratio
       QW_col(L)=MAX(EPSQ, QWnew)              ! Updated cloud water mixing ratio
       RimeF_col(L)=RimeF                      ! Updated rime factor
       ASNOW2=ASNOWnew                          ! Updated accumulated snow
       ARAIN=ARAINnew                          ! Updated accumulated rain
       !! 2012-02-09 Based on Aligo's email
       !aligo
       VRabove=VRAIN2                          ! Updated rain fall speed
       !aligo

       !!
       !
       !#######################################################################
       !
10     CONTINUE         ! ##### End "L" loop through model levels #####
    END DO ! DO 10 L=1,LSFC
!!!!!!!!        write(6,*)'2012-02-09,ncw,nlimax,vrabove=',ncw,nlimax,vrabove
    !
    !#######################################################################
    !
    !-----------------------------------------------------------------------
    !--------------------------- Return to GSMDRIVE -----------------------
    !-----------------------------------------------------------------------
    !
  CONTAINS
    !#######################################################################
    !--------- Produces accurate calculation of cloud condensation ---------
    !#######################################################################
    !
    REAL(KIND=r8) FUNCTION CONDENSE (PP, QW, TK, WV, RHgrd)   !GFDL
      !
      !---------------------------------------------------------------------------------
      !------   The Asai (1965) algorithm takes into consideration the release of ------
      !------   latent heat in increasing the temperature & in increasing the     ------
      !------   saturation mixing ratio (following the Clausius-Clapeyron eqn.).  ------
      !---------------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      INTEGER, PARAMETER :: HIGH_PRES=SELECTED_REAL_KIND(15)
      REAL (KIND=HIGH_PRES), PARAMETER :: RHLIMIT=0.001_r8
      REAL (KIND=HIGH_PRES), PARAMETER ::RHLIMIT1=-RHLIMIT
      REAL (KIND=HIGH_PRES) :: COND, SSAT, WCdum
      !
      REAL(KIND=r8),INTENT(IN) :: QW,PP,WV,TK,RHgrd   !GFDL
      REAL(KIND=r8)            ::  WVdum,Tdum,XLV2,DWV,WS,ESW,XLV1,XLV
      INTEGER :: nsteps
      !
      !-----------------------------------------------------------------------
      !
      !--- LV (T) is from Bolton (JAS, 1980)
      !
      XLV=3.148E6_r8-2370.0_r8*TK
      XLV1=XLV*RCP
      XLV2=XLV*XLV*RCPRV
      Tdum=TK
      WVdum=WV
      WCdum=QW
      ESW=1000.0_r8*FPVS0(Tdum)                     ! Saturation vapor press w/r/t water
      WS=RHgrd*EPS*ESW/(PP-ESW)                 ! Saturation mixing ratio w/r/t water
      DWV=WVdum-WS                              ! Deficit grid-scale water vapor mixing ratio
      SSAT=DWV/WS                               ! Supersaturation ratio
      CONDENSE=0.0_r8
      nsteps = 0
      DO WHILE ((SSAT.LT.RHLIMIT1 .AND. WCdum.GT.EPSQ).OR. SSAT.GT.RHLIMIT)
         nsteps = nsteps + 1
         COND=DWV/(1.0_r8+XLV2*WS/(Tdum*Tdum))       ! Asai (1965, J. Japan)
         COND=MAX(COND, -WCdum)                  ! Limit cloud water evaporation
         Tdum=Tdum+XLV1*COND                     ! Updated temperature
         WVdum=WVdum-COND                        ! Updated water vapor mixing ratio
         WCdum=WCdum+COND                        ! Updated cloud water mixing ratio
         CONDENSE=CONDENSE+COND                  ! Total cloud water condensation
         ESW=1000.0_r8*FPVS0(Tdum)                   ! Updated saturation vapor press w/r/t water
         WS=RHgrd*EPS*ESW/(PP-ESW)               ! Updated saturation mixing ratio w/r/t water
         DWV=WVdum-WS                            ! Deficit grid-scale water vapor mixing ratio
         SSAT=DWV/WS                             ! Grid-scale supersaturation ratio
      ENDDO
      !
    END FUNCTION CONDENSE
    !
    !#######################################################################
    !---------------- Calculate ice deposition at T<T_ICE ------------------
    !#######################################################################
    !
    REAL(KIND=r8) FUNCTION DEPOSIT (PP, Tdum, WVdum, RHgrd)   !GFDL
      !
      !--- Also uses the Asai (1965) algorithm, but uses a different target
      !      vapor pressure for the adjustment
      !
      IMPLICIT NONE      
      !
      INTEGER, PARAMETER :: HIGH_PRES=SELECTED_REAL_KIND(15)
      REAL (KIND=HIGH_PRES), PARAMETER :: RHLIMIT=0.001_r8
      REAL (KIND=HIGH_PRES), PARAMETER :: RHLIMIT1=-RHLIMIT
      REAL (KIND=HIGH_PRES) :: DEP, SSAT
      !    
      REAL(KIND=r8),INTENT(IN)    ::  PP,RHgrd   !GFDL
      REAL(KIND=r8),INTENT(INOUT) ::  WVdum,Tdum
      REAL(KIND=r8) :: ESI,WS,DWV
      !
      !-----------------------------------------------------------------------
      !
      ESI=1000.0_r8*FPVS(Tdum)                      ! Saturation vapor press w/r/t ice
      WS=RHgrd*EPS*ESI/(PP-ESI)                 ! Saturation mixing ratio
      DWV=WVdum-WS                              ! Deficit grid-scale water vapor mixing ratio
      SSAT=DWV/WS                               ! Supersaturation ratio
      DEPOSIT=0.0_r8
      DO WHILE (SSAT.GT.RHLIMIT .OR. SSAT.LT.RHLIMIT1)
         !
         !--- Note that XLVS2=LS*LV/(CP*RV)=LV*WS/(RV*T*T)*(LS/CP*DEP1), 
         !     where WS is the saturation mixing ratio following Clausius-
         !     Clapeyron (see Asai,1965; Young,1993,p.405) 
         !
         DEP=DWV/(1.0_r8+XLS2*WS/(Tdum*Tdum))        ! Asai (1965, J. Japan)
         Tdum=Tdum+XLS1*DEP                      ! Updated temperature
         WVdum=WVdum-DEP                         ! Updated ice mixing ratio
         DEPOSIT=DEPOSIT+DEP                     ! Total ice deposition
         ESI=1000.0_r8*FPVS(Tdum)                    ! Updated saturation vapor press w/r/t ice
         WS=RHgrd*EPS*ESI/(PP-ESI)               ! Updated saturation mixing ratio w/r/t ice
         DWV=WVdum-WS                            ! Deficit grid-scale water vapor mixing ratio
         SSAT=DWV/WS                             ! Grid-scale supersaturation ratio
      ENDDO
      !
    END FUNCTION DEPOSIT
    !
  END SUBROUTINE EGCP01COLUMN

  SUBROUTINE MY_GROWTH_RATES (DTPH)
    !     SUBROUTINE MY_GROWTH_RATES (DTPH,MY_GROWTH)
    !
    !--- Below are tabulated values for the predicted mass of ice crystals
    !    after 600 s of growth in water saturated conditions, based on 
    !    calculations from Miller and Young (JAS, 1979).  These values are
    !    crudely estimated from tabulated curves at 600 s from Fig. 6.9 of
    !    Young (1993).  Values at temperatures colder than -27C were 
    !    assumed to be invariant with temperature.  
    !
    !--- Used to normalize Miller & Young (1979) calculations of ice growth
    !    over large time steps using their tabulated values at 600 s.
    !    Assumes 3D growth with time**1.5 following eq. (6.3) in Young (1993).
    !
    IMPLICIT NONE
    !
    REAL(KIND=r8)  ,INTENT(IN ) :: DTPH
    !
    REAL(KIND=r8)              ::  DT_ICE
    REAL(KIND=r8),DIMENSION(35) :: MY_600
    !WRF
    !
    !-----------------------------------------------------------------------
    DATA MY_600 /                                                     &
         & 5.5e-8_r8 , 1.4E-7_r8 , 2.8E-7_r8, 6.E-7_r8  , 3.3E-6_r8 , & 
         & 2.E-6_r8  , 9.E-7_r8  , 8.8E-7_r8, 8.2E-7_r8 , 9.4e-7_r8 , & 
         & 1.2E-6_r8 , 1.85E-6_r8, 5.5E-6_r8, 1.5E-5_r8 , 1.7E-5_r8 , & 
         & 1.5E-5_r8 , 1.E-5_r8  , 3.4E-6_r8, 1.85E-6_r8, 1.35E-6_r8, & 
         & 1.05E-6_r8, 1.E-6_r8  , 9.5E-7_r8, 9.0E-7_r8 , 9.5E-7_r8 , & 
         & 9.5E-7_r8 , 9.E-7_r8  , 9.E-7_r8 , 9.E-7_r8  , 9.E-7_r8  , & 
         & 9.E-7_r8  , 9.E-7_r8  , 9.E-7_r8 , 9.E-7_r8  , 9.E-7_r8 /        ! -31 to -35 deg C
    !
    !-----------------------------------------------------------------------
    !
    DT_ICE   =(DTPH/600.0_r8)**1.5_r8
    MY_GROWTH=DT_ICE*MY_600
    !
    !-----------------------------------------------------------------------
    !
  END SUBROUTINE MY_GROWTH_RATES
  !
  !-----------------------------------------------------------------------
  !---------  Old GFS saturation vapor pressure lookup tables  -----------
  !-----------------------------------------------------------------------
  !
  SUBROUTINE GPVS(NX,XMAX,XMIN,C1XPVS,C2XPVS,C1XPVS0,C2XPVS0,TBPVS,TBPVS0)
    !     ******************************************************************
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !                .      .    .
    ! SUBPROGRAM:    GPVS        COMPUTE SATURATION VAPOR PRESSURE TABLE
    !   AUTHOR: N PHILLIPS       W/NP2      DATE: 30 DEC 82
    !
    ! ABSTRACT: COMPUTE SATURATION VAPOR PRESSURE TABLE AS A FUNCTION OF
    !   TEMPERATURE FOR THE TABLE LOOKUP FUNCTION FPVS.
    !   EXACT SATURATION VAPOR PRESSURES ARE CALCULATED IN SUBPROGRAM FPVSX.
    !   THE CURRENT IMPLEMENTATION COMPUTES A TABLE WITH A LENGTH
    !   OF 7501 FOR TEMPERATURES RANGING FROM 180.0 TO 330.0 KELVIN.
    !
    ! PROGRAM HISTORY LOG:
    !   91-05-07  IREDELL
    !   94-12-30  IREDELL             EXPAND TABLE
    !   96-02-19  HONG                ICE EFFECT
    !   01-11-29  JIN                 MODIFIED FOR WRF
    !
    ! USAGE:  CALL GPVS
    !
    ! SUBPROGRAMS CALLED:
    !   (FPVSX)  - INLINABLE FUNCTION TO COMPUTE SATURATION VAPOR PRESSURE
    !
    ! COMMON BLOCKS:
    !   COMPVS   - SCALING PARAMETERS AND TABLE FOR FUNCTION FPVS.
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: NX
    REAL(KIND=r8)   , INTENT(IN   ) :: XMAX
    REAL(KIND=r8)   , INTENT(IN   ) :: XMIN
    REAL(KIND=r8)   , INTENT(OUT  ) :: C1XPVS
    REAL(KIND=r8)   , INTENT(OUT  ) :: C2XPVS
    REAL(KIND=r8)   , INTENT(OUT  ) :: C1XPVS0
    REAL(KIND=r8)   , INTENT(OUT  ) :: C2XPVS0
    REAL(KIND=r8)   , INTENT(OUT  ) :: TBPVS (NX)
    REAL(KIND=r8)   , INTENT(OUT  ) :: TBPVS0(NX)

    REAL(KIND=r8)    :: X
    REAL(KIND=r8)    :: XINC
    REAL(KIND=r8)    :: T
    INTEGER :: JX
    !----------------------------------------------------------------------
    XINC=(XMAX-XMIN)/(NX-1)
    C1XPVS=1.0_r8-XMIN/XINC
    C2XPVS=1.0_r8/XINC
    C1XPVS0=1.0_r8-XMIN/XINC
    C2XPVS0=1.0_r8/XINC
    !
    DO JX=1,NX
       X=XMIN+(JX-1)*XINC
       T=X
       TBPVS (JX)=FPVSX (T)
       TBPVS0(JX)=FPVSX0(T)
    ENDDO
    ! 
  END SUBROUTINE GPVS
  !-----------------------------------------------------------------------
  !***********************************************************************
  !-----------------------------------------------------------------------
  REAL(KIND=r8)   FUNCTION FPVS(T)
    !-----------------------------------------------------------------------
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !                .      .    .
    ! SUBPROGRAM:    FPVS        COMPUTE SATURATION VAPOR PRESSURE
    !   AUTHOR: N PHILLIPS            W/NP2      DATE: 30 DEC 82
    !
    ! ABSTRACT: COMPUTE SATURATION VAPOR PRESSURE FROM THE TEMPERATURE.
    !   A LINEAR INTERPOLATION IS DONE BETWEEN VALUES IN A LOOKUP TABLE
    !   COMPUTED IN GPVS. SEE DOCUMENTATION FOR FPVSX FOR DETAILS.
    !   INPUT VALUES OUTSIDE TABLE RANGE ARE RESET TO TABLE EXTREMA.
    !   THE INTERPOLATION ACCURACY IS ALMOST 6 DECIMAL PLACES.
    !   ON THE CRAY, FPVS IS ABOUT 4 TIMES FASTER THAN EXACT CALCULATION.
    !   THIS FUNCTION SHOULD BE EXPANDED INLINE IN THE CALLING ROUTINE.
    !
    ! PROGRAM HISTORY LOG:
    !   91-05-07  IREDELL             MADE INTO INLINABLE FUNCTION
    !   94-12-30  IREDELL             EXPAND TABLE
    !   96-02-19  HONG                ICE EFFECT
    !   01-11-29  JIN                 MODIFIED FOR WRF
    !
    ! USAGE:   PVS=FPVS(T)
    !
    !   INPUT ARGUMENT LIST:
    !     T        - REAL TEMPERATURE IN KELVIN
    !
    !   OUTPUT ARGUMENT LIST:
    !     FPVS     - REAL SATURATION VAPOR PRESSURE IN KILOPASCALS (CB)
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !$$$
    IMPLICIT NONE
    REAL(KIND=r8),INTENT(IN) :: T
    REAL(KIND=r8) XJ
    INTEGER :: JX
    !-----------------------------------------------------------------------
    XJ=MIN(MAX(C1XPVS+C2XPVS*T,1.0_r8),FLOAT(NX))
    JX=INT(MIN(XJ,NX-1.0_r8))
    FPVS=TBPVS(JX)+(XJ-JX)*(TBPVS(JX+1)-TBPVS(JX))
    !
  END FUNCTION FPVS
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  REAL(KIND=r8) FUNCTION FPVS0(T)
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=r8),INTENT(IN) :: T
    REAL(KIND=r8)    :: XJ1
    INTEGER :: JX1
    !-----------------------------------------------------------------------
    XJ1=MIN(MAX(C1XPVS0+C2XPVS0*T,1.0_r8),FLOAT(NX))
    JX1=INT(MIN(XJ1,NX-1.0_r8))
    FPVS0=TBPVS0(JX1)+(XJ1-JX1)*(TBPVS0(JX1+1)-TBPVS0(JX1))
    !
  END FUNCTION FPVS0
  !-----------------------------------------------------------------------
  !***********************************************************************
  !-----------------------------------------------------------------------
  REAL(KIND=r8) FUNCTION FPVSX(T)
    !-----------------------------------------------------------------------
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !                .      .    .
    ! SUBPROGRAM:    FPVSX       COMPUTE SATURATION VAPOR PRESSURE
    !   AUTHOR: N PHILLIPS            W/NP2      DATE: 30 DEC 82
    !
    ! ABSTRACT: EXACTLY COMPUTE SATURATION VAPOR PRESSURE FROM TEMPERATURE.
    !   THE WATER MODEL ASSUMES A PERFECT GAS, CONSTANT SPECIFIC HEATS
    !   FOR GAS AND LIQUID, AND NEGLECTS THE VOLUME OF THE LIQUID.
    !   THE MODEL DOES ACCOUNT FOR THE VARIATION OF THE LATENT HEAT
    !   OF CONDENSATION WITH TEMPERATURE.  THE ICE OPTION IS NOT INCLUDED.
    !   THE CLAUSIUS-CLAPEYRON EQUATION IS INTEGRATED FROM THE TRIPLE POINT
    !   TO GET THE FORMULA
    !       PVS=PSATK*(TR**XA)*EXP(XB*(1.-TR))
    !   WHERE TR IS TTP/T AND OTHER VALUES ARE PHYSICAL CONSTANTS
    !   THIS FUNCTION SHOULD BE EXPANDED INLINE IN THE CALLING ROUTINE.
    !
    ! PROGRAM HISTORY LOG:
    !   91-05-07  IREDELL             MADE INTO INLINABLE FUNCTION
    !   94-12-30  IREDELL             EXACT COMPUTATION
    !   96-02-19  HONG                ICE EFFECT 
    !   01-11-29  JIN                 MODIFIED FOR WRF
    !
    ! USAGE:   PVS=FPVSX(T)
    ! REFERENCE:   EMANUEL(1994),116-117
    !
    !   INPUT ARGUMENT LIST:
    !     T        - REAL(KIND=r8) TEMPERATURE IN KELVIN
    !
    !   OUTPUT ARGUMENT LIST:
    !     FPVSX    - REAL(KIND=r8) SATURATION VAPOR PRESSURE IN KILOPASCALS (CB)
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !$$$
    IMPLICIT NONE
    REAL(KIND=r8) , INTENT(IN   ) :: T
    !-----------------------------------------------------------------------
    REAL(KIND=r8),PARAMETER :: TTP=2.7316E+2_r8
    REAL(KIND=r8),PARAMETER :: HVAP=2.5000E+6_r8
    REAL(KIND=r8),PARAMETER :: PSAT=6.1078E+2_r8
    REAL(KIND=r8),PARAMETER :: CLIQ=4.1855E+3_r8
    REAL(KIND=r8),PARAMETER :: CVAP= 1.8460E+3_r8
    REAL(KIND=r8),PARAMETER :: CICE=2.1060E+3_r8
    REAL(KIND=r8),PARAMETER :: HSUB=2.8340E+6_r8
    !
    REAL(KIND=r8),PARAMETER :: PSATK=PSAT*1.E-3_r8
    REAL(KIND=r8),PARAMETER :: DLDT=CVAP-CLIQ
    REAL(KIND=r8),PARAMETER :: XA=-DLDT/RV
    REAL(KIND=r8),PARAMETER :: XB=XA+HVAP/(RV*TTP)
    REAL(KIND=r8),PARAMETER :: DLDTI=CVAP-CICE
    REAL(KIND=r8),PARAMETER :: XAI=-DLDTI/RV
    REAL(KIND=r8),PARAMETER :: XBI=XAI+HSUB/(RV*TTP)
    REAL(KIND=r8) :: TR
    !-----------------------------------------------------------------------
    TR=TTP/T
    !
    IF(T.GE.TTP)THEN
       FPVSX=PSATK*(TR**XA)*EXP(XB*(1.0_r8-TR))
    ELSE
       FPVSX=PSATK*(TR**XAI)*EXP(XBI*(1.0_r8-TR))
    ENDIF
    ! 
  END FUNCTION FPVSX
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  REAL(KIND=r8)   FUNCTION FPVSX0(T)
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN   ) :: T

    REAL(KIND=r8),PARAMETER :: TTP=2.7316E+2_r8
    REAL(KIND=r8),PARAMETER :: HVAP=2.5000E+6_r8
    REAL(KIND=r8),PARAMETER :: PSAT=6.1078E+2_r8
    REAL(KIND=r8),PARAMETER :: CLIQ=4.1855E+3_r8
    REAL(KIND=r8),PARAMETER :: CVAP=1.8460E+3_r8
    REAL(KIND=r8),PARAMETER :: CICE=2.1060E+3_r8
    REAL(KIND=r8),PARAMETER :: HSUB=2.8340E+6_r8
    REAL(KIND=r8),PARAMETER :: PSATK=PSAT*1.E-3_r8
    REAL(KIND=r8),PARAMETER :: DLDT=CVAP-CLIQ
    REAL(KIND=r8),PARAMETER :: XA=-DLDT/RV
    REAL(KIND=r8),PARAMETER :: XB=XA+HVAP/(RV*TTP)
    REAL(KIND=r8),PARAMETER :: DLDTI=CVAP-CICE
    REAL(KIND=r8),PARAMETER :: XAI=-DLDT/RV
    REAL(KIND=r8),PARAMETER :: XBI=XA+HSUB/(RV*TTP)
    REAL(KIND=r8) :: TR
    !-----------------------------------------------------------------------
    TR=TTP/T
    FPVSX0=PSATK*(TR**XA)*EXP(XB*(1.0_r8-TR))
    !
  END FUNCTION FPVSX0
  !
END MODULE Micro_HWRF


!               PROGRAM Main
!                 USE Micro_HWRF
!                 IMPLICIT NONE
!               END PROGRAM Main
