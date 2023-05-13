MODULE MicroPhysics
  USE Constants, ONLY :  &
       delq             ,&
       r8,i8,qmin,grav,gasr

  USE Micro_Hack, ONLY:       &
      Init_Micro_Hack ,RunMicro_Hack   
 
  USE Micro_HWRF, ONLY:       &
      Init_Micro_HWRF ,RunMicro_HWRF   
 
  USE Micro_Ferrier, ONLY:       &
      Init_Micro_Ferrier ,RunMicro_FERRIER

  USE Micro_UKME, ONLY:       &
      Init_Micro_UKME ,RunMicro_UKME

  USE Micro_MORR, ONLY:       &
      Init_Micro_MORR ,RunMicro_MORR

  USE Micro_HugMorr, ONLY:       &
      Init_Micro_HugMorr ,RunMicro_HugMorr

  USE Micro_HugMorr_NN, ONLY:       &
      Init_Micro_HugMorr_NN ,RunMicro_HugMorr_NN, hugMorr_spinup_reached, fill_vars_timesteps

  USE Micro_LrgScl, ONLY:     &
      Init_Micro_LrgScl,RunMicro_LrgScl 

  USE FieldsPhysics, ONLY:  &
       LOWLYR             , &
       F_ICE_PHY          , &
       F_RAIN_PHY         , &
       F_RIMEF_PHY       


   ! TEMPORARY
   use csv_file


 IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER :: ppcnst=3
  REAL(KIND=r8), ALLOCATABLE  :: ql2(:,:,:)
  REAL(KIND=r8), ALLOCATABLE  :: qi2(:,:,:)
  REAL(KIND=r8), ALLOCATABLE  :: ql3(:,:,:)
  REAL(KIND=r8), ALLOCATABLE  :: qi3(:,:,:)
  REAL(KIND=r8), ALLOCATABLE  :: qr2(:,:,:)
  REAL(KIND=r8), ALLOCATABLE  :: qr3(:,:,:)

  REAL(KIND=r8), ALLOCATABLE  :: qs2(:,:,:)
  REAL(KIND=r8), ALLOCATABLE  :: qs3(:,:,:)

  REAL(KIND=r8), ALLOCATABLE  :: qg2(:,:,:)
  REAL(KIND=r8), ALLOCATABLE  :: qg3(:,:,:)

  REAL(KIND=r8), ALLOCATABLE  :: NI2(:,:,:)
  REAL(KIND=r8), ALLOCATABLE  :: NI3(:,:,:)

  REAL(KIND=r8), ALLOCATABLE  :: NS2(:,:,:)
  REAL(KIND=r8), ALLOCATABLE  :: NS3(:,:,:)
  
  REAL(KIND=r8), ALLOCATABLE  :: NR2(:,:,:)
  REAL(KIND=r8), ALLOCATABLE  :: NR3(:,:,:)

  REAL(KIND=r8), ALLOCATABLE  :: NG2(:,:,:)
  REAL(KIND=r8), ALLOCATABLE  :: NG3(:,:,:)

  REAL(KIND=r8), ALLOCATABLE  :: NC2(:,:,:)
  REAL(KIND=r8), ALLOCATABLE  :: NC3(:,:,:)


  PUBLIC InitMicroPhysics
  PUBLIC RunMicroPhysics
  PUBLIC FinalizeMicroPhysics
CONTAINS 
  SUBROUTINE InitMicroPhysics( &
                          dt      ,si        ,del       ,sl         , &
                          kMax    ,iMax      ,jMax      ,ibMax      , &
                          jbMax   ,fNameMicro, ILCON    ,microphys  , &
                          nClass  ,nAeros    ,EFFCS     ,EFFIS)
    IMPLICIT NONE
    REAL(KIND=r8)   , INTENT(IN   ) :: dt
    INTEGER         , INTENT(IN   ) :: kMax
    INTEGER         , INTENT(IN   ) :: iMax
    INTEGER         , INTENT(IN   ) :: jMax
    INTEGER         , INTENT(IN   ) :: ibMax
    INTEGER         , INTENT(IN   ) :: jbMax
    REAL(KIND=r8)   , INTENT(IN   ) :: si (kMax+1)
    REAL(KIND=r8)   , INTENT(IN   ) :: del(kMax  )
    REAL(KIND=r8)   , INTENT(IN   ) :: sl (kMax  )
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameMicro
    CHARACTER(LEN=*), INTENT(IN   ) :: ILCON
    LOGICAL         , INTENT(IN   ) :: microphys
    INTEGER         , INTENT(IN   ) :: nClass
    INTEGER         , INTENT(IN   ) :: nAeros
    REAL(KIND=r8)   , INTENT(OUT  ) :: EFFCS(ibMax,kMax,jbMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: EFFIS(ibMax,kMax,jbMax)
    CHARACTER(LEN=*), PARAMETER     :: h='**(InitMicroPhysics)**'
    LOGICAL                         :: restart

    IF(TRIM(ILCON).EQ.'YES'.or.TRIM(ILCON).EQ.'LSC') CALL Init_Micro_LrgScl()    
    IF(TRIM(ILCON).EQ.'MIC') CALL Init_Micro_Hack(kMax,jMax,ibMax,jbMax,ppcnst,si,sl,del)
    IF(TRIM(ILCON).EQ.'HWRF') CALL Init_Micro_HWRF( &
       dt               , &
       LOWLYR           , &
       restart          , &
       fNameMicro       , &
       F_ICE_PHY        , &
       F_RAIN_PHY       , &
       F_RIMEF_PHY      , &
       ibMax            , &
       jbMax            , &
       kMax               )
    IF(TRIM(ILCON).EQ.'HGFS') CALL Init_Micro_Ferrier (dt,iMax,jMax,ibmax,kMax,jbMax,F_ICE_PHY,F_RAIN_PHY ,F_RIMEF_PHY,si,sl)
    IF(TRIM(ILCON).EQ.'UKMO') CALL Init_Micro_UKME(iMax)
    IF(TRIM(ILCON).EQ.'MORR') CALL Init_Micro_MORR ()
    IF(TRIM(ILCON).EQ.'HUMO') then
        CALL Init_Micro_HugMorr (ibMax,kMax,jbMax,EFFCS,EFFIS)
    ELSE IF(TRIM(ILCON).EQ.'HUMN') then
      !   CALL Init_Micro_HugMorr (ibMax,kMax,jbMax,EFFCS,EFFIS)
        CALL Init_Micro_HugMorr_NN(dt)
    END IF

    IF (microphys) THEN
       ALLOCATE(ql3(ibMax,kMax,jbMax));ql3=0.00001e-12_r8
       ALLOCATE(qi3(ibMax,kMax,jbMax));qi3=0.00001e-12_r8

       ALLOCATE(ql2(ibMax,kMax,jbMax));ql2=0.00001e-12_r8
       ALLOCATE(qi2(ibMax,kMax,jbMax));qi2=0.00001e-12_r8
       IF((nClass+nAeros)>0)THEN
          ALLOCATE(qr3(ibMax,kMax,jbMax));qr3=0.00001e-12_r8
          IF( TRIM(ILCON).EQ.'MORR'.or.TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HUMN') THEN
             ALLOCATE(qs3(ibMax,kMax,jbMax));qs3=0.00001e-12_r8
             ALLOCATE(qg3(ibMax,kMax,jbMax));qg3=0.00001e-12_r8
             ALLOCATE(NI3(ibMax,kMax,jbMax));NI3=0.00001e-12_r8
             ALLOCATE(NS3(ibMax,kMax,jbMax));NS3=0.00001e-12_r8
             ALLOCATE(NR3(ibMax,kMax,jbMax));NR3=0.00001e-12_r8
             ALLOCATE(NG3(ibMax,kMax,jbMax));NG3=0.00001e-12_r8
             IF( TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HUMN')THEN
                ALLOCATE(NC3(ibMax,kMax,jbMax));NC3=0.00001e-12_r8
             END IF
          END IF
          ALLOCATE(qr2(ibMax,kMax,jbMax));qr2=0.00001e-12_r8
          IF( TRIM(ILCON).EQ.'MORR'.or.TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HUMN') THEN
             ALLOCATE(qs2(ibMax,kMax,jbMax));qs2=0.00001e-12_r8
             ALLOCATE(qg2(ibMax,kMax,jbMax));qg2=0.00001e-12_r8
             ALLOCATE(NI2(ibMax,kMax,jbMax));NI2=0.00001e-12_r8
             ALLOCATE(NS2(ibMax,kMax,jbMax));NS2=0.00001e-12_r8
             ALLOCATE(NR2(ibMax,kMax,jbMax));NR2=0.00001e-12_r8
             ALLOCATE(NG2(ibMax,kMax,jbMax));NG2=0.00001e-12_r8
             IF( TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HUMN') THEN 
                ALLOCATE(NC2(ibMax,kMax,jbMax));NC2=0.00001e-12_r8
             END IF
          END IF
       END IF 
    END IF

  END SUBROUTINE InitMicroPhysics
!**********************************************************************************
!**********************************************************************************


  SUBROUTINE RunMicroPhysics(&
      ! Run Flags
               mlrg        , & !INTEGER      , INTENT(in   ) :: mlrg
               ILCON       , & !CHARACTER(LEN=*), INTENT(IN   ) :: ILCON
               microphys   , & !LOGICAL      , INTENT(in   ) :: microphys
      ! Time info
               jdt         , & !INTEGER      , INTENT(in   ) :: jdt
               dt          , & !REAL(KIND=r8), INTENT(in   ) :: dt
      ! Model Geometry
               sl          , & !REAL(KIND=r8)   , INTENT(in   ) :: sl (kMax)
               si          , & !REAL(KIND=r8)   , INTENT(in   ) :: si (kMax+1)
               del         , & !REAL(KIND=r8)   , INTENT(in   ) :: del (kMax)
               colrad      , & !REAL(KIND=r8)   , INTENT(IN   ) :: colrad  (iMax)  
               LOWLYR      , & !INTEGER         , INTENT(IN   ) :: LOWLYR  (iMax)  
               terr        , & !REAL(KIND=r8)   , INTENT(IN   ) :: terr (iMax)  
      ! Model information
               mask        , & !INTEGER(KIND=i8), INTENT(IN   ) :: mask (iMax) 
               nClass      , & !INTEGER      , INTENT(in   ) :: nClass
               nAeros      , &
               iMax        , & !INTEGER      , INTENT(in   ) :: iMax
               kMax        , & !INTEGER      , INTENT(in   ) :: kMax
               latco       , & !INTEGER      , INTENT(in   ) :: latco
      ! Surface field
               tsfc        , & !REAL(KIND=r8)   , INTENT(IN   ) :: tsfc (iMax)
      ! PBL field
               PBL_CoefKh  , & !REAL(KIND=r8)   , INTENT(IN   ) :: PBL_CoefKh(1:iMax,1:kMax)         , &
               tke         , & !REAL(KIND=r8)   , INTENT(IN   ) :: tke(1:iMax,1:kMax)         , &
               pblh        , & !REAL(KIND=r8)   , INTENT(IN   ) :: pblh (iMax)  
      ! CONVECTION: Cloud field
               kuo         , & !INTEGER         , INTENT(in   ) :: kuo (iMax)  
               cmfmc       , & !REAL(KIND=r8)   , INTENT(IN   ) :: cmfmc   (iMax,kMax+1)   ! convective mass flux--m sub c
               cmfmc2      , & !REAL(KIND=r8)   , INTENT(IN   ) :: cmfmc2  (iMax,kMax+1)   ! shallow convective mass flux--m sub c
               dlf         , & !REAL(KIND=r8)   , INTENT(IN   ) :: dlf (iMax,kMax)! detrained water from ZM
               rliq        , & !REAL(KIND=r8)   , INTENT(IN   ) :: rliq (iMax)        ! vertical integral of liquid not yet in q(ixcldliq)
               concld      , & !REAL(KIND=r8)   , INTENT(INOUT) :: concld  (iMax,kMax)
               cld         , & !REAL(KIND=r8)   , INTENT(INOUT) :: cld (iMax,kMax)
               fdqn        , & !REAL(KIND=r8)   , INTENT(inout) :: fdqn (iMax,kMax)
               EFFCS       , & !REAL(KIND=r8)   , INTENT(INOUT) :: EFFCS(1:iMax,1:kMax)         , &
               EFFIS       , & !REAL(KIND=r8)   , INTENT(INOUT) :: EFFIS(1:iMax,1:kMax)         , &
               Total_Rain  , & !REAL(KIND=r8)   , INTENT(inout) :: Total_Rain  (iMax)
               Total_Snow  , & !REAL(KIND=r8)   , INTENT(inout) :: Total_Snow  (iMax)
               RAINNCV     , & !REAL(KIND=r8)   , INTENT(INOUT) :: RAINNCV     (iMax)
               SNOWNCV     , & !REAL(KIND=r8)   , INTENT(INOUT) :: SNOWNCV     (iMax)
               F_ICE_PHY   , & !REAL(KIND=r8)   , INTENT(INOUT) :: F_ICE_PHY   (1:iMax,1:kMax)
               F_RAIN_PHY  , & !REAL(KIND=r8)   , INTENT(INOUT) :: F_RAIN_PHY  (1:iMax,1:kMax)
               F_RIMEF_PHY , & !REAL(KIND=r8)   , INTENT(INOUT) :: F_RIMEF_PHY (1:iMax,1:kMax)
      ! Atmospheric fields
               ps2         , & !REAL(KIND=r8)   , INTENT(in   ) :: ps2 (iMax)
               t2          , & !REAL(KIND=r8)   , INTENT(IN   ) :: t2(iMax,kMax)
               t3          , & !REAL(KIND=r8)   , INTENT(inout) :: t3 (iMax,kMax) 
               q2          , & !REAL(KIND=r8)   , INTENT(IN   ) :: q2 (iMax,kMax)
               q3          , & !REAL(KIND=r8)   , INTENT(inout) :: q3 (iMax,kMax)
               ub          , & !REAL(KIND=r8)   , INTENT(IN   ) :: ub (iMax,kMax) ! (m/s) 
               vb          , & !REAL(KIND=r8)   , INTENT(IN   ) :: vb (iMax,kMax) ! (m/s)
               omgb        , & !REAL(KIND=r8)   , INTENT(IN   ) :: omgb (iMax,kMax) ! (Pa/s)
               dq          , & !REAL(KIND=r8)   , INTENT(inout) :: dq (iMax,kMax)
               tLrgs       , & !REAL(KIND=r8)   , INTENT(OUT  ) :: tLrgs       (1:iMax,1:kMax)
               qLrgs       , & !REAL(KIND=r8)   , INTENT(OUT  ) :: qLrgs       (1:iMax,1:kMax)
      ! Microphysics
               gicem       , & !REAL(KIND=r8)   , INTENT(INOUT) :: gicem       (iMax,kmax)
               gicep       , & !REAL(KIND=r8)   , INTENT(INOUT) :: gicep       (iMax,kmax)
               gliqm       , & !REAL(KIND=r8)   , INTENT(INOUT) :: gliqm       (iMax,kmax)
               gliqp       , & !REAL(KIND=r8)   , INTENT(INOUT) :: gliqp       (iMax,kmax)
               gvarm       , & !REAL(KIND=r8)   , OPTIONAL,   INTENT(INOUT) :: gvarm (iMax,kmax,nClass+nAeros)
               gvarp         ) !REAL(KIND=r8)   , OPTIONAL,   INTENT(INOUT) :: gvarp (iMax,kmax,nClass+nAeros)
    !************************************************************************
    !   The cumulus_driver subroutine calls deep and shallow cumulus
    !   parameterization schemes.
    !   more information nilo@cptec.inpe.br
    !   NOTE: This version is not official. You can use only for test.
    !************************************************************************
    !
    !  Definition/
    !---------------
    !             I N P U T  O U T P U T  F O R   G C M
    !             -------------------------------------
    ! INPUT
    !
    !** integer
    !    iMax                   ! end index for longitude domain
    !    kMax                   ! end index for u,v,t,p sigma levels
    !    jdt                    ! number of time step
    !    iccon                  ! cu schemes ex. KUO, ARA, GRE, RAS ..
    !   kuo                     ! convection yes(1) or not(0) for shallow convection
    !
    !** real
    !    dt                     ! time step (s)
    !    ta                     ! temperature (K) at time t-1
    !    tb                     ! temperature (K) at time t
    !    tc                     ! temperature (K) at time t+1
    !    qa                     ! water vapor mixing ratio (kg/kg) at time t-1
    !    qb                     ! water vapor mixing ratio (kg/kg) at time t
    !    qc                     ! water vapor mixing ratio (kg/kg) at time t+1
    !    psb                    ! surface pressure (cb)     at time t
    !    ub                     ! u-velocity (m/s) at time t
    !    vb                     ! v-velocity (m/s) at time t
    !    omgb                   ! vertical omega velocity (Pa/s) at time t
    !                           ! it is in half levels along with U,V,T,Q
    !    sl                     ! half sigma layers
    !    si                     ! full sigma layers
    !    del                    ! difference between full sigma layers
    !    xland                  ! land-sea mask (1 for land; 0 for water)
    !    zs                     ! topography (m)
    !    DX                     ! horizontal space interval (m)
    !    qrem,cldm              ! local variables for  RAS-Scheme
    !
    !    hrem,qrem              ! these arrays are needed for the heating 
    !                           ! and mostening from ras  scheme
    !
    !
    !    ktops, kbots           ! these arrays are needed for the new 
    !                           ! shallow convection scheme
    !    cldm                   ! needed for cloud fraction based on mass 
    !                           ! flux
    !    noshal1, kctop1, kcbot1! needed for cloud fraction based on mass 
    !                           ! flux new arrays needed for shallow 
    !                           ! convective clouds
    !     
    !
    !
    !   OUTPUT
    !**  integer
    !    kuo                    ! indix for shalow convection KUO,RAS,KUOG, GRELL
    !    ktop                   ! level of convective cloud top
    !    kbot                   ! level of convective cloud base
    !    plcl                   ! pressure at convection levl for shallow convection
    !                           ! in Kuo 
    !
    !** real
    !   RAINCV                  ! cumulus scheme precipitation (mm)
    !   tc                      ! new temperature (K) at time t+1  after CU precip
    !   qc                      ! new  water vapor mixing ratio (kg/kg) at time t+1.
    !
    !
    !*********************************************************************************
    IMPLICIT NONE
    !              I N P U T     O U T P U T    V A R I A B L E S
    !              ----------------------------------------------
    !              Xa at t-1   Xb at t    Xc at t+1

    INTEGER         , INTENT(in   ) :: jdt
    INTEGER         , INTENT(in   ) :: nClass
    INTEGER         , INTENT(IN   ) :: nAeros
    INTEGER         , INTENT(in   ) :: iMax
    INTEGER         , INTENT(in   ) :: kMax
    REAL(KIND=r8)   , INTENT(in   ) :: dt
    INTEGER         , INTENT(in   ) :: mlrg
    INTEGER         , INTENT(in   ) :: latco
    LOGICAL         , INTENT(in   ) :: microphys
    CHARACTER(LEN=*), INTENT(IN   ) :: ILCON
    REAL(KIND=r8)   , INTENT(in   ) :: sl      (kMax)
    REAL(KIND=r8)   , INTENT(in   ) :: si      (kMax+1)
    REAL(KIND=r8)   , INTENT(in   ) :: del     (kMax)
    INTEGER(KIND=i8), INTENT(IN   ) :: mask    (iMax) 
    REAL(KIND=r8)   , INTENT(IN   ) :: tsfc    (iMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: terr    (iMax)  
    REAL(KIND=r8)   , INTENT(IN   ) :: colrad  (iMax)  
    REAL(KIND=r8)   , INTENT(IN   ) :: pblh    (iMax)  
    INTEGER         , INTENT(in   ) :: kuo     (iMax)  
    INTEGER         , INTENT(IN   ) :: LOWLYR  (iMax)  
    REAL(KIND=r8)   , INTENT(in   ) :: ps2     (iMax)
    REAL(KIND=r8)   , INTENT(inout) :: t2      (iMax,kMax)
    REAL(KIND=r8)   , INTENT(inout) :: t3      (iMax,kMax) 
    REAL(KIND=r8)   , INTENT(inout) :: q2      (iMax,kMax)
    REAL(KIND=r8)   , INTENT(inout) :: q3      (iMax,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: ub      (iMax,kMax) ! (m/s) 
    REAL(KIND=r8)   , INTENT(IN   ) :: vb      (iMax,kMax) ! (m/s)
    REAL(KIND=r8)   , INTENT(IN   ) :: omgb    (iMax,kMax) ! (Pa/s)
    REAL(KIND=r8)   , INTENT(inout) :: dq      (iMax,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: cmfmc   (iMax,kMax+1)   ! convective mass flux--m sub c
    REAL(KIND=r8)   , INTENT(IN   ) :: cmfmc2  (iMax,kMax+1)   ! shallow convective mass flux--m sub c
    REAL(KIND=r8)   , INTENT(IN   ) :: dlf     (iMax,kMax)    ! detrained water from ZM
    REAL(KIND=r8)   , INTENT(IN   ) :: rliq    (iMax)        ! vertical integral of liquid not yet in q(ixcldliq)
    REAL(KIND=r8)   , INTENT(INOUT) :: concld  (iMax,kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: cld     (iMax,kMax)
    REAL(KIND=r8)   , INTENT(inout) :: fdqn    (iMax,kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: EFFCS      (1:iMax,1:kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: EFFIS      (1:iMax,1:kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: PBL_CoefKh (1:iMax,1:kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: tke(1:iMax,1:kMax) 
    REAL(KIND=r8)   , INTENT(inout) :: Total_Rain  (iMax)
    REAL(KIND=r8)   , INTENT(inout) :: Total_Snow  (iMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: RAINNCV     (iMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: SNOWNCV     (iMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: F_ICE_PHY   (1:iMax,1:kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: F_RAIN_PHY  (1:iMax,1:kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: F_RIMEF_PHY (1:iMax,1:kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: tLrgs       (1:iMax,1:kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: qLrgs       (1:iMax,1:kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: gicem       (iMax,kmax)
    REAL(KIND=r8)   , INTENT(INOUT) :: gicep       (iMax,kmax)
    REAL(KIND=r8)   , INTENT(INOUT) :: gliqm       (iMax,kmax)
    REAL(KIND=r8)   , INTENT(INOUT) :: gliqp       (iMax,kmax)
    REAL(KIND=r8)   , OPTIONAL,   INTENT(INOUT) :: gvarm (iMax,kmax,nClass+nAeros)
    REAL(KIND=r8)   , OPTIONAL,   INTENT(INOUT) :: gvarp (iMax,kmax,nClass+nAeros)

    ! WIND COMPONENTS FOR GRELL ENSEMBLE
    REAL(KIND=r8) :: u2(iMax,kMax)
    REAL(KIND=r8) :: v2(iMax,kMax)
    REAL(KIND=r8) :: w2(iMax,kMax)
    REAL(KIND=r8) :: prec_str(iMax)  ! [Total] sfc flux of precip from stratiform (m/s) 
    REAL(KIND=r8) :: snow_str(iMax)  ! [Total] sfc flux of snow from stratiform   (m/s)
    REAL(KIND=r8) :: prec_sed(iMax)  ! surface flux of total cloud water from sedimentation
    REAL(KIND=r8) :: snow_sed(iMax)  ! surface flux of cloud ice from sedimentation
    REAL(KIND=r8) :: prec_pcw(iMax)  ! sfc flux of precip from microphysics(m/s)
    REAL(KIND=r8) :: snow_pcw(iMax)  ! sfc flux of snow from microphysics (m/s)
    REAL(KIND=r8) :: icefrac (iMax)
    REAL(KIND=r8) :: landfrac(iMax)
    REAL(KIND=r8) :: ocnfrac (iMax)
    REAL(KIND=r8) :: landm   (iMax)  ! land fraction ramped over water
    REAL(KIND=r8) :: snowh   (iMax)  ! Snow depth over land, water equivalent (m)
    REAL(KIND=r8) :: ts      (iMax)      ! surface temperature
    REAL(KIND=r8) :: sst     (iMax)       !sea surface temperature
    REAL(KIND=r8) :: tv   (iMax,kMax)
    REAL(KIND=r8) :: press(iMax,kMax)
    REAL(KIND=r8) :: RHO  (iMax,kMax)
    REAL(KIND=r8) :: delz (iMax,kMax)
    REAL(KIND=r8) :: r1000
    REAL(KIND=r8) :: rbyg
    INTEGER       :: i
    INTEGER       :: k

    logical reached, passed

    IF (microphys) THEN
       DO k=1,kMax
          DO i=1,iMax
             ql3     (i,k,latco) = gliqp(i,k)
             qi3     (i,k,latco) = gicep(i,k)

             ql2     (i,k,latco) = gliqm(i,k) 
             qi2     (i,k,latco) = gicem(i,k)
             IF((nClass+nAeros)>0 .and. PRESENT(gvarm))THEN
                qr3     (i,k,latco) = gvarp(i,k,1)
                IF( TRIM(ILCON).EQ.'MORR'.or.TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HUMN') THEN
                   qs3     (i,k,latco) = gvarp(i,k,2)
                   qg3     (i,k,latco) = gvarp(i,k,3)
                   NI3     (i,k,latco) = gvarp(i,k,4)
                   NS3     (i,k,latco) = gvarp(i,k,5)
                   NR3     (i,k,latco) = gvarp(i,k,6)
                   NG3     (i,k,latco) = gvarp(i,k,7)
                   IF( TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HUMN') NC3     (i,k,latco) = gvarp(i,k,8)
                END IF
                qr2     (i,k,latco) = gvarm(i,k,1)
                IF( TRIM(ILCON).EQ.'MORR'.or.TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HUMN') THEN
                   qs2     (i,k,latco) = gvarm(i,k,2)
                   qg2     (i,k,latco) = gvarm(i,k,3)
                   NI2     (i,k,latco) = gvarm(i,k,4)
                   NS2     (i,k,latco) = gvarm(i,k,5)
                   NR2     (i,k,latco) = gvarm(i,k,6)
                   NG2     (i,k,latco) = gvarm(i,k,7)
                   IF( TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HUMN') NC2     (i,k,latco) = gvarm(i,k,8)
                END IF
             END IF 
          END DO
       END DO
    ELSE
       DO k=1,kMax
          DO i=1,iMax
             !ql3     (i,k,latco) = 0.0_r8
             !qi3     (i,k,latco) = 0.0_r8

             !ql2     (i,k,latco) = 0.0_r8
             !qi2     (i,k,latco) = 0.0_r8
             IF((nClass+nAeros)>0 )THEN
                qr3     (i,k,latco) = 0.0_r8 
                IF( TRIM(ILCON).EQ.'MORR'.or.TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HUMN') THEN
                   qs3     (i,k,latco) = 0.0_r8 
                   qg3     (i,k,latco) = 0.0_r8 
                   NI3     (i,k,latco) = 0.0_r8
                   NS3     (i,k,latco) = 0.0_r8
                   NR3     (i,k,latco) = 0.0_r8
                   NG3     (i,k,latco) = 0.0_r8
                   IF( TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HUMN')NC3     (i,k,latco) = 0.0_r8
                END IF
                qr2     (i,k,latco) = 0.0_r8
                IF( TRIM(ILCON).EQ.'MORR'.or.TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HUMN') THEN
                   qs2     (i,k,latco) = 0.0_r8
                   qg2     (i,k,latco) = 0.0_r8
                   NI2     (i,k,latco) = 0.0_r8
                   NS2     (i,k,latco) = 0.0_r8
                   NR2     (i,k,latco) = 0.0_r8
                   NG2     (i,k,latco) = 0.0_r8
                  IF( TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HUMN') NC2     (i,k,latco) = 0.0_r8
                END IF
             END IF
          END DO
       END DO
    END IF

    !-----------------------------------------------------------------
    ! Large Scale Precipitation
    !-----------------------------------------------------------------
    IF(TRIM(ILCON).EQ.'LSC' .OR. TRIM(ILCON).EQ.'YES' ) THEN
        ! print*, 'Enver before q3*del', iMax, kMax, size(del(:)), size(q3)
      CALL RunMicro_LrgScl(Total_Rain, t3, dq, q3, ps2, del, sl, dt, &
                mlrg, latco, iMax, kMax)
        !print*, 'Enver after q3*del',  q3, del
      CALL qnegat2 (q3, fdqn,(1.0_r8/dt), del, iMax, kMax)! time t+1
        snow_str=0.0_r8
        prec_sed=0.0_r8
        snow_sed=0.0_r8
        prec_pcw=0.0_r8
        snow_pcw=0.0_r8
    ELSE IF(TRIM(ILCON).EQ.'MIC') THEN
        snow_str=0.0_r8
        prec_sed=0.0_r8
        snow_sed=0.0_r8
        prec_pcw=0.0_r8
        snow_pcw=0.0_r8
        snowh   =0.0_r8
        icefrac =0.0_r8
        landfrac=0.0_r8
        ocnfrac =0.0_r8
        DO i=1,iMax
           IF(mask(i) == 13_i8 .or. mask(i) == 15_i8 )snowh(i)=5.0_r8
           IF(mask(i).GT.0_i8)THEN
              ! land
              icefrac(i)=0.0_r8
              landfrac(i)=1.0_r8
              ocnfrac(i)=0.0_r8
           ELSE
              ! water/ocean
              landfrac(i)=0.0_r8
              ocnfrac(i) =1.0_r8
              IF(ocnfrac(i).GT.0.01_r8.AND.tsfc(i).LT.260.0_r8) THEN
                 icefrac(i)=1.0_r8
                 ocnfrac(i) =0.0_r8
              ENDIF
           END IF
        END DO
        DO k=1,kMax
           DO i=1,iMax
              u2(i,k)=ub  (i,k)
              v2(i,k)=vb  (i,k)
              w2(i,k)=omgb(i,k)*1000.0_r8  ! (Pa/s)
           END DO
        END DO
!----------------------------
        DO k=1,kMax
           DO i=1,iMax
              press(i,k)=ps2(i)*1000_r8*sl(k)
              tv   (i,k)=t2 (i,k)*(1.0_r8+0.608_r8*q2(i,k))
           END DO
        END DO    
        !
        !  Calculate the distance between the surface and the first layer of the model
        !
        r1000=1000.0e0_r8 /gasr
        rbyg=gasr/grav*del(1)*0.5e0_r8
        DO i=1,iMax
           RHO     (i,1)=r1000*(ps2(i)*sl(1))/t2(i,1)
           delz    (i,1)=MAX((rbyg * tv(i,1)),0.5_r8)*0.75_r8
        END DO 

        DO k=2,kMax
           DO i=1,iMax
              RHO (i,k)=r1000*(ps2(i)*sl(k))/t2(i,k)
              delz(i,k)=0.5_r8*gasr*(tv(i,k-1)+tv(i,k))* &
                        LOG(press(i,k-1)/press(i,k))/grav
           END DO
        END DO
!--------------------------------------------
       landm          =0.0_r8  
       ts             =tsfc
       sst            =tsfc
       CALL RunMicro_Hack( &
       jdt                              , &! INTEGER , INTENT(in)  :: ibMax                     !number of columns (max)
       iMax                             , &! INTEGER , INTENT(in)  :: ibMax                     !number of columns (max)
       kMax                             , &! INTEGER , INTENT(in)  :: kMax                      !number of vertical levels
       kMax+1                           , &! INTEGER , INTENT(in)  :: kMax+1                    !number of vertical levels + 1
       ppcnst                           , &! INTEGER , INTENT(in)  :: ppcnst                    !number of constituent
       latco                            , &! INTEGER , INTENT(in)  :: latco                     !latitude
       dt                               , &! REAL(r8), INTENT(in)  :: dtime                     !timestep
       ps2         (1:iMax)*1000_r8     , &! REAL(r8), INTENT(in)  :: state_ps    (ibMax)       !(ibMax)     ! surface pressure(Pa)
       t2          (1:iMax,1:kMax)      , &! REAL(r8), INTENT(in)  :: state_t     (ibMax,kMax)  !(ibMax,kMax)! temperature (K)
       t3          (1:iMax,1:kMax)      , &! REAL(r8), INTENT(in)  :: state_t     (ibMax,kMax)  !(ibMax,kMax)! temperature (K)
       q2          (1:iMax,1:kMax)      , &! REAL(r8), INTENT(in)  :: state_qv    (ibMax,kMax)  !(ibMax,kMax,ppcnst)! vapor  mixing ratio (kg/kg moist or dry air depending on type)
       q3          (1:iMax,1:kMax)      , &! REAL(r8), INTENT(in)  :: state_qv    (ibMax,kMax)  !(ibMax,kMax,ppcnst)! vapor  mixing ratio (kg/kg moist or dry air depending on type)
       ql2         (1:iMax,1:kMax,latco), &! REAL(r8), INTENT(in)  :: state_ql    (pcols,pver)  !(pcols,pver,ppcnst)! liquid mixing ratio (kg/kg moist or dry air depending on type)
       ql3         (1:iMax,1:kMax,latco), &! REAL(r8), INTENT(in)  :: state_ql    (pcols,pver)  !(pcols,pver,ppcnst)! liquid mixing ratio (kg/kg moist or dry air depending on type)
       qi2         (1:iMax,1:kMax,latco), &! REAL(r8), INTENT(in)  :: state_qi    (pcols,pver)  !(pcols,pver,ppcnst)! ice    mixing ratio (kg/kg moist or dry air depending on type)
       qi3         (1:iMax,1:kMax,latco), &! REAL(r8), INTENT(in)  :: state_qi    (pcols,pver)  !(pcols,pver,ppcnst)! ice    mixing ratio (kg/kg moist or dry air depending on type)
       w2          (1:iMax,1:kMax)      , &! REAL(r8), INTENT(in)  :: state_omega (ibMax,kMax)  !(ibMax,kMax)! vertical pressure velocity (Pa/s) 
       icefrac     (1:iMax)             , &! REAL(r8), INTENT(in)  :: icefrac     (ibMax)       !sea ice fraction (fraction)
       landfrac    (1:iMax)             , &! REAL(r8), INTENT(in)  :: landfrac    (ibMax)       !land fraction (fraction)
       ocnfrac     (1:iMax)             , &! REAL(r8), INTENT(in)  :: ocnfrac     (ibMax)       !ocean fraction (fraction)
       landm       (1:iMax)             , &! REAL(r8), INTENT(in)  :: landm       (ibMax)       !land fraction ramped over water
       snowh       (1:iMax)             , &! REAL(r8), INTENT(in)  :: snowh       (ibMax)       !Snow depth over land, water equivalent (m)
       dlf         (1:iMax,1:kMax)      , &! REAL(r8), INTENT(in)  :: state_dlf   (ibMax,kMax)  !detrained water from ZM
       rliq        (1:iMax)             , &! REAL(r8), INTENT(in)  :: rliq        (ibMax)       !vertical integral of liquid not yet in q(ixcldliq)
       cmfmc       (1:iMax,1:kMax+1)    , &! REAL(r8), INTENT(in)  :: state_cmfmc (ibMax,kMax+1)!convective mass flux--m sub c
       cmfmc2      (1:iMax,1:kMax+1)    , &! REAL(r8), INTENT(in)  :: state_cmfmc2(ibMax,kMax+1)!shallow convective mass flux--m sub c
       concld      (1:iMax,1:kMax)      , &! REAL(r8), INTENT(out) :: state_concld(ibMax,kMax)  !convective cloud cover
       cld         (1:iMax,1:kMax)      , &! REAL(r8), INTENT(out) :: state_cld   (ibMax,kMax)  !cloud fraction
       sst         (1:iMax)             , &! REAL(r8), INTENT(in)  :: sst         (ibMax)       !sea surface temperature
       !zdu         (1:iMax,1:kMax)     , &! REAL(r8), INTENT(in)  :: state_zdu  (ibMax,kMax)  !detrainment rate from deep convection
       prec_str    (1:iMax)             , &! REAL(r8), INTENT(out)  :: prec_str   (ibMax)       ![Total] sfc flux of precip from stratiform (m/s) 
       snow_str    (1:iMax)             , &! REAL(r8), INTENT(out)  :: snow_str   (ibMax)       ![Total] sfc flux of snow from stratiform   (m/s)
       prec_sed    (1:iMax)             , &! REAL(r8), INTENT(out)  :: prec_sed   (ibMax)       !surface flux of total cloud water from sedimentation
       snow_sed    (1:iMax)             , &! REAL(r8), INTENT(out)  :: snow_sed   (ibMax)       !surface flux of cloud ice from sedimentation
       prec_pcw    (1:iMax)             , &! REAL(r8), INTENT(out)  :: prec_pcw   (ibMax)       !sfc flux of precip from microphysics(m/s)
       snow_pcw    (1:iMax)             )  ! REAL(r8), INTENT(out)  :: snow_pcw   (ibMax)       !sfc flux of snow from microphysics (m/s)
       Total_Rain=Total_Rain+MAX((prec_sed*0.5_r8*dt),0.0_r8)+MAX((prec_pcw*0.5_r8*dt),0.0_r8)
       Total_Snow=Total_Snow+MAX((snow_sed*0.5_r8*dt),0.0_r8)+MAX((snow_pcw*0.5_r8*dt),0.0_r8)
       CALL qnegat2 (q3, fdqn,(1.0_r8/dt), del, iMax, kMax)! time t+1

      ELSE IF(TRIM(ILCON).EQ.'HGFS') THEN
       ! grell mask
       snowh   =0.0_r8
       icefrac =0.0_r8
       landfrac=0.0_r8
       ocnfrac =0.0_r8
       DO i=1,iMax
          IF(mask(i) == 13 .or. mask(i) == 15 )snowh(i)=5.0_r8      
          IF(mask(i).GT.0_i8)THEN
             ! land
             icefrac (i)=0.0_r8
             landfrac(i)=1.0_r8
             ocnfrac (i)=0.0_r8
          ELSE
             ! water/ocean
             landfrac(i)=0.0_r8
             ocnfrac(i) =1.0_r8
             IF(ocnfrac(i).GT.0.01_r8.AND.tsfc(i).LT.260.0_r8) THEN
                icefrac(i)=1.0_r8
                ocnfrac(i) =0.0_r8
             ENDIF
          END IF
       END DO
       DO k=1,kMax
          DO i=1,iMax
             u2(i,k)=ub(i,k)
             v2(i,k)=vb(i,k)
             w2(i,k)=omgb(i,k)*1000.0_r8  ! (Pa/s)
          END DO
       END DO
!----------------------------
       DO k=1,kMax
          DO i=1,iMax
             press(i,k)=ps2(i)*1000_r8*sl(k)
             tv   (i,k)=t2 (i,k)*(1.0_r8+0.608_r8*q2(i,k))
          END DO
       END DO    
       !
       !  Calculate the distance between the surface and the first layer of the model
       !
       r1000=1000.0e0_r8 /gasr
       rbyg=gasr/grav*del(1)*0.5e0_r8
       DO i=1,iMax
          RHO     (i,1)=r1000*(ps2(i)*sl(1))/t2(i,1)
          delz    (i,1)=MAX((rbyg * tv(i,1)),0.5_r8)*0.75_r8
       END DO

       DO k=2,kMax
          DO i=1,iMax
             RHO (i,k)=r1000*(ps2(i)*sl(k))/t2(i,k)
             delz(i,k)=0.5_r8*gasr*(tv(i,k-1)+tv(i,k))* &
                       LOG(press(i,k-1)/press(i,k))/grav
          END DO
       END DO
  !-----------------------------------------------------------------------
       RAINNCV=0.0_r8
       SNOWNCV=0.0_r8
       CALL RunMicro_FERRIER (&
       iMax                            , &!INTEGER      , INTENT(IN)     :: nCols
       kMax                            , &!INTEGER      , INTENT(IN)     :: kMax
       DT                              , &!REAL(KIND=r8), INTENT(IN)     :: DT
       RHO        (1:iMax,1:kMax)      , &!REAL(KIND=r8), INTENT(IN),    :: DEL       (ims:ime, kms:kme, jms:jme)
       t3         (1:iMax,1:kMax)      , &!REAL(KIND=r8), INTENT(INOUT), :: gt (ims:ime, kms:kme, jms:jme)
       q3         (1:iMax,1:kMax)      , &!REAL(KIND=r8), INTENT(INOUT), :: qv (ims:ime, kms:kme, jms:jme)
       ql3        (1:iMax,1:kMax,latco), &!REAL(KIND=r8), INTENT(INOUT), :: qc (ims:ime, kms:kme, jms:jme)
       qi3        (1:iMax,1:kMax,latco), &!REAL(KIND=r8), INTENT(INOUT), :: qi (ims:ime, kms:kme, jms:jme)
       qr3        (1:iMax,1:kMax,latco), &!REAL(KIND=r8), INTENT(INOUT), :: qr (ims:ime, kms:kme, jms:jme)
       F_ICE_PHY  (1:iMax,1:kMax)      , &!REAL(KIND=r8), INTENT(INOUT), :: F_ICE_PHY  (ims:ime, kms:kme, jms:jme)            , &!REAL(KIND=r8), INTENT(INOUT), :: F_ICE_PHY  (ims:ime, kms:kme, jms:jme)
       F_RAIN_PHY (1:iMax,1:kMax)      , &!REAL(KIND=r8), INTENT(INOUT), :: F_RAIN_PHY (ims:ime, kms:kme, jms:jme)       , &!REAL(KIND=r8), INTENT(INOUT), :: F_RAIN_PHY (ims:ime, kms:kme, jms:jme)
       F_RIMEF_PHY(1:iMax,1:kMax)      , &!REAL(KIND=r8), INTENT(INOUT), :: F_RIMEF_PHY(ims:ime, kms:kme, jms:jme)
       RAINNCV    (1:iMax)             , &!REAL(KIND=r8), INTENT(INOUT), :: RAINNCV    (ims:ime,  jms:jme) 
       SNOWNCV    (1:iMax)             , &!REAL(KIND=r8), INTENT(INOUT), :: SNOWNCV    (ims:ime,  jms:jme) 
       ps2        (1:iMax)*1000_r8     , &
       colrad     (1:iMax)               )

!--------------------------------------------
        Total_Rain=Total_Rain + RAINNCV
        Total_Snow=Total_Snow + SNOWNCV
        snow_str=0.0_r8
        prec_sed=0.0_r8
        snow_sed=0.0_r8
        prec_pcw=0.0_r8
        snow_pcw=0.0_r8
       CALL qnegat2 (q3, fdqn,(1.0_r8/dt), del, iMax, kMax)! time t+1
      ELSE IF(TRIM(ILCON).EQ.'HWRF') THEN
       ! grell mask
       snowh   =0.0_r8
       icefrac =0.0_r8
       landfrac=0.0_r8
       ocnfrac =0.0_r8
       DO i=1,iMax
          IF(mask(i) == 13 .or. mask(i) == 15 )snowh(i)=5.0_r8      
          IF(mask(i).GT.0_i8)THEN
             ! land
             icefrac (i)=0.0_r8
             landfrac(i)=1.0_r8
             ocnfrac (i)=0.0_r8
          ELSE
             ! water/ocean
             landfrac(i)=0.0_r8
             ocnfrac(i) =1.0_r8
             IF(ocnfrac(i).GT.0.01_r8.AND.tsfc(i).LT.260.0_r8) THEN
                icefrac(i)=1.0_r8
                ocnfrac(i) =0.0_r8
             ENDIF
          END IF
       END DO
       DO k=1,kMax
          DO i=1,iMax
             u2(i,k)=ub(i,k)
             v2(i,k)=vb(i,k)
             w2(i,k)=omgb(i,k)*1000.0_r8  ! (Pa/s)
          END DO
       END DO
!----------------------------
       DO k=1,kMax
          DO i=1,iMax
             press(i,k)=ps2(i)*1000_r8*sl(k)
             tv   (i,k)=t2 (i,k)*(1.0_r8+0.608_r8*q2(i,k))
          END DO
       END DO    
       !
       !  Calculate the distance between the surface and the first layer of the model
       !
       r1000=1000.0e0_r8 /gasr
       rbyg=gasr/grav*del(1)*0.5e0_r8
       DO i=1,iMax
          RHO     (i,1)=r1000*(ps2(i)*sl(1))/t2(i,1)
          delz    (i,1)=MAX((rbyg * tv(i,1)),0.5_r8)*0.75_r8
       END DO

       DO k=2,kMax
          DO i=1,iMax
             RHO (i,k)=r1000*(ps2(i)*sl(k))/t2(i,k)
             delz(i,k)=0.5_r8*gasr*(tv(i,k-1)+tv(i,k))* &
                       LOG(press(i,k-1)/press(i,k))/grav
          END DO
       END DO
  !-----------------------------------------------------------------------
       RAINNCV=0.0_r8
       SNOWNCV=0.0_r8
       CALL RunMicro_HWRF (&
       iMax                             , &!INTEGER      , INTENT(IN)     :: nCols
       kMax                             , &!INTEGER      , INTENT(IN)     :: kMax
       DT                               , &!REAL(KIND=r8), INTENT(IN)     :: DT
       delz       (1:iMax,1:kMax)       , &!REAL(KIND=r8), INTENT(IN),    :: dz8w       (ims:ime, kms:kme, jms:jme)
       RHO        (1:iMax,1:kMax)       , &!REAL(KIND=r8), INTENT(IN),    :: rho_phy    (ims:ime, kms:kme, jms:jme)
       press      (1:iMax,1:kMax)       , &!REAL(KIND=r8), INTENT(IN),    :: p_phy      (ims:ime, kms:kme, jms:jme)
       LOWLYR     (1:iMax       )       , &!INTEGER      , INTENT(IN   )  :: LOWLYR     (ims:ime         , jms:jme)
       t3         (1:iMax,1:kMax)       , &!REAL(KIND=r8), INTENT(INOUT), :: gt         (ims:ime, kms:kme, jms:jme)
       q3         (1:iMax,1:kMax)       , &!REAL(KIND=r8), INTENT(INOUT), :: qv         (ims:ime, kms:kme, jms:jme)
       ql3        (1:iMax,1:kMax,latco) , &!REAL(KIND=r8), INTENT(INOUT), :: qc         (ims:ime, kms:kme, jms:jme)
       qi3        (1:iMax,1:kMax,latco) , &!REAL(KIND=r8), INTENT(INOUT), :: qi         (ims:ime, kms:kme, jms:jme)
       qr3        (1:iMax,1:kMax,latco) , &!REAL(KIND=r8), INTENT(INOUT), :: qr         (ims:ime, kms:kme, jms:jme)
       F_ICE_PHY  (1:iMax,1:kMax) , &!REAL(KIND=r8), INTENT(INOUT), :: F_ICE_PHY  (ims:ime, kms:kme, jms:jme)
       F_RAIN_PHY (1:iMax,1:kMax) , &!REAL(KIND=r8), INTENT(INOUT), :: F_RAIN_PHY (ims:ime, kms:kme, jms:jme)
       F_RIMEF_PHY(1:iMax,1:kMax) , &!REAL(KIND=r8), INTENT(INOUT), :: F_RIMEF_PHY(ims:ime, kms:kme, jms:jme)
       RAINNCV    (1:iMax)              , &!REAL(KIND=r8), INTENT(INOUT), :: RAINNCV    (ims:ime,          jms:jme) !GID
       SNOWNCV    (1:iMax)                )
!--------------------------------------------
        Total_Rain=Total_Rain + RAINNCV
        Total_Snow=Total_Snow + SNOWNCV
        snow_str=0.0_r8
        prec_sed=0.0_r8
        snow_sed=0.0_r8
        prec_pcw=0.0_r8
        snow_pcw=0.0_r8
       CALL qnegat2 (q3, fdqn,(1.0_r8/dt), del, iMax, kMax)! time t+1

      ELSE IF(TRIM(ILCON).EQ.'UKMO') THEN
       ! grell mask
       snowh   =0.0_r8
       icefrac =0.0_r8
       landfrac=0.0_r8
       ocnfrac =0.0_r8
       DO i=1,iMax
          IF(mask(i) == 13 .or. mask(i) == 15 )snowh(i)=5000.0_r8      
          IF(mask(i).GT.0_i8)THEN
             ! land
             icefrac (i)=0.0_r8
             landfrac(i)=1.0_r8
             ocnfrac (i)=0.0_r8
          ELSE
             ! water/ocean
             landfrac(i)=0.0_r8
             ocnfrac(i) =1.0_r8
             IF(ocnfrac(i).GT.0.01_r8.AND.tsfc(i).LT.260.0_r8) THEN
                icefrac(i)=1.0_r8
                ocnfrac(i) =0.0_r8
             ENDIF
          END IF
       END DO
       DO k=1,kMax
          DO i=1,iMax
             u2(i,k)=ub(i,k)
             v2(i,k)=vb(i,k)
             w2(i,k)=omgb(i,k)*1000.0_r8  ! (Pa/s)
          END DO
       END DO
!----------------------------
       DO k=1,kMax
          DO i=1,iMax
             press(i,k)=ps2(i)*1000_r8*sl(k)
             tv   (i,k)=t2 (i,k)*(1.0_r8+0.608_r8*q2(i,k))
          END DO
       END DO    
       !
       !  Calculate the distance between the surface and the first layer of the model
       !
       r1000=1000.0e0_r8 /gasr
       rbyg=gasr/grav*del(1)*0.5e0_r8
       DO i=1,iMax
          RHO     (i,1)=r1000*(ps2(i)*sl(1))/t2(i,1)
          delz    (i,1)=MAX((rbyg * tv(i,1)),0.5_r8)*0.75_r8
       END DO

       DO k=2,kMax
          DO i=1,iMax
             RHO (i,k)=r1000*(ps2(i)*sl(k))/t2(i,k)
             delz(i,k)=0.5_r8*gasr*(tv(i,k-1)+tv(i,k))* &
                       LOG(press(i,k-1)/press(i,k))/grav
          END DO
       END DO
  !-----------------------------------------------------------------------

       RAINNCV=0.0_r8
       SNOWNCV=0.0_r8
       CALL RunMicro_UKME(&
                        kMax                            , &
                        iMax                            , &
                        si         (1:kMax+1)           , &
                        sl         (1:kMax)             , &
                        dt                              , &
                        pblh       (1:iMax      )       , &
                        colrad     (1:iMax)             , &
                        kuo        (1:iMax)             , &
                        q3         (1:iMax,1:kMax)      , &!q              , &
                        qi3        (1:iMax,1:kMax,latco), &!QCF    , &
                        ql3        (1:iMax,1:kMax,latco), &!QCL    , &
                        qr3        (1:iMax,1:kMax,latco), &!qcf2   , &
                        t3         (1:iMax,1:kMax)      , &!T  , &
                        ps2        (1:iMax)*1000_r8     , &
                        F_RIMEF_PHY(1:iMax,1:kMax      ), &! CF     , &
                        F_RAIN_PHY (1:iMax,1:kMax      ), &! CFL    , &
                        F_ICE_PHY  (1:iMax,1:kMax      ), &! CFF    , &
                        snowh      (1:iMax)             , &
                        landfrac   (1:iMax)             , &
                        terr       (1:iMax)             , &
                        RAINNCV    (1:iMax)             , &
                        SNOWNCV    (1:iMax)               )


!--------------------------------------------
        Total_Rain=Total_Rain + RAINNCV
        Total_Snow=Total_Snow + SNOWNCV
        snow_str=0.0_r8
        prec_sed=0.0_r8
        snow_sed=0.0_r8
        prec_pcw=0.0_r8
        snow_pcw=0.0_r8
       CALL qnegat2 (q3, fdqn,(1.0_r8/dt), del, iMax, kMax)! time t+1

      ELSE IF(TRIM(ILCON).EQ.'MORR') THEN
       ! grell mask
       snowh   =0.0_r8
       icefrac =0.0_r8
       landfrac=0.0_r8
       ocnfrac =0.0_r8
       DO i=1,iMax
          IF(mask(i) == 13 .or. mask(i) == 15 )snowh(i)=5000.0_r8      
          IF(mask(i).GT.0_i8)THEN
             ! land
             icefrac (i)=0.0_r8
             landfrac(i)=1.0_r8
             ocnfrac (i)=0.0_r8
          ELSE
             ! water/ocean
             landfrac(i)=0.0_r8
             ocnfrac(i) =1.0_r8
             IF(ocnfrac(i).GT.0.01_r8.AND.tsfc(i).LT.260.0_r8) THEN
                icefrac(i)=1.0_r8
                ocnfrac(i) =0.0_r8
             ENDIF
          END IF
       END DO
       DO k=1,kMax
          DO i=1,iMax
             u2(i,k)=ub(i,k)
             v2(i,k)=vb(i,k)
             w2(i,k)=omgb(i,k)*1000.0_r8  ! (Pa/s)
          END DO
       END DO
  !-----------------------------------------------------------------------

       RAINNCV=0.0_r8
       SNOWNCV=0.0_r8
       CALL RunMicro_MORR( &
       iMax                                 , &!INTEGER      , INTENT(IN   ) :: nCols
       kMax                                 , &!INTEGER      , INTENT(IN   ) :: kMax 
       si          (1:kMax+1)               , &
       sl          (1:kMax)                 , &
       t3          (1:iMax,1:kMax)          , &!REAL(KIND=r8), INTENT(INOUT) :: Tc (1:nCols, 1:kMax)
       q3          (1:iMax,1:kMax)          , &!REAL(KIND=r8), INTENT(INOUT) :: qv (1:nCols, 1:kMax)
       ql3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qc (1:nCols, 1:kMax)
       qr3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qr (1:nCols, 1:kMax)
       qi3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qi (1:nCols, 1:kMax)
       QS3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qs (1:nCols, 1:kMax)
       QG3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qg (1:nCols, 1:kMax)
       NI3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: ni (1:nCols, 1:kMax)
       NS3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: ns (1:nCols, 1:kMax)
       NR3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: nr (1:nCols, 1:kMax)
       NG3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: NG (1:nCols, 1:kMax)   
       ps2         (1:iMax)*1000_r8         , & ! gps- AIR PRESSURE (PA)
       DT                                   , &!REAL(KIND=r8), INTENT(IN   ) :: dt_in
       w2         (1:iMax,1:kMax)           , &!REAL(KIND=r8), INTENT(IN   ) :: omega  ! omega (Pa/s)
       RAINNCV    (1:iMax)                  , &!REAL(KIND=r8), INTENT(OUT) :: LSRAIN(1:nCols)
       SNOWNCV    (1:iMax)                    )!REAL(KIND=r8), INTENT(OUT) :: LSSNOW(1:nCols)

!--------------------------------------------
        Total_Rain=Total_Rain + RAINNCV
        Total_Snow=Total_Snow + SNOWNCV
        snow_str=0.0_r8
        prec_sed=0.0_r8
        snow_sed=0.0_r8
        prec_pcw=0.0_r8
        snow_pcw=0.0_r8
       CALL qnegat2 (q3, fdqn,(1.0_r8/dt), del, iMax, kMax)! time t+1

      ELSE IF(TRIM(ILCON).EQ.'HUMO' .or. TRIM(ILCON).EQ.'HUMN') THEN
       ! grell mask
       snowh   =0.0_r8
       icefrac =0.0_r8
       landfrac=0.0_r8
       ocnfrac =0.0_r8
       DO i=1,iMax
          IF(mask(i) == 13 .or. mask(i) == 15 )snowh(i)=5000.0_r8      
          IF(mask(i).GT.0_i8)THEN
             ! land
             icefrac (i)=0.0_r8
             landfrac(i)=1.0_r8
             ocnfrac (i)=0.0_r8
          ELSE
             ! water/ocean
             landfrac(i)=0.0_r8
             ocnfrac(i) =1.0_r8
             IF(ocnfrac(i).GT.0.01_r8.AND.tsfc(i).LT.260.0_r8) THEN
                icefrac(i)=1.0_r8
                ocnfrac(i) =0.0_r8
             ENDIF
          END IF
       END DO
       DO k=1,kMax
          DO i=1,iMax
             u2(i,k)=ub(i,k)
             v2(i,k)=vb(i,k)
             w2(i,k)=omgb(i,k)*1000.0_r8  ! (Pa/s)
          END DO
       END DO
  !-----------------------------------------------------------------------

       print*, "==================> jdt = ", jdt, " timestep = ", jdt*dt, "TRIM(ILCON)", TRIM(ILCON)


       IF(TRIM(ILCON).EQ.'HUMN') THEN
         call hugMorr_spinup_reached(jdt, reached, passed)

         if (reached) then 
            print*, "filling reached, passed -------------", reached, passed
            call fill_vars_timesteps( &
            iMax                                 , &!INTEGER      , INTENT(IN   ) :: nCols
            kMax                                 , &!INTEGER      , INTENT(IN   ) :: kMax
            sl          (1:kMax)               , &
            t3          (1:iMax,1:kMax)          , &!REAL(KIND=r8), INTENT(INOUT) :: Tc (1:nCols, 1:kMax)
            q3          (1:iMax,1:kMax)          , &!REAL(KIND=r8), INTENT(INOUT) :: qv (1:nCols, 1:kMax)
            ql3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qc (1:nCols, 1:kMax)
            qr3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qr (1:nCols, 1:kMax)
            qi3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qi (1:nCols, 1:kMax)
            QS3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qs (1:nCols, 1:kMax)
            QG3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qg (1:nCols, 1:kMax)
            NI3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: ni (1:nCols, 1:kMax)
            NS3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: ns (1:nCols, 1:kMax)
            NR3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: nr (1:nCols, 1:kMax)
            NG3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: NG (1:nCols, 1:kMax)
            NC3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: NC (1:nCols, 1:kMax)
            TKE         (1:iMax,1:kMax)          , &!REAL(KIND=r8), INTENT(IN   ) :: KZH_TKE (1:nCols, 1:kMax) (M^2 S-1)
            w2          (1:iMax,1:kMax)          , &!REAL(KIND=r8), INTENT(IN   ) :: omega  ! omega (Pa/s)
            RAINNCV     (1:iMax)                 , &!REAL(KIND=r8), INTENT(OUT)   :: LSRAIN(1:nCols)
            SNOWNCV     (1:iMax)                   &!REAL(KIND=r8), INTENT(OUT)   :: LSSNOW(1:nCols)
            )
         endif

         RAINNCV=0.0_r8
         SNOWNCV=0.0_r8

         if (reached .and. passed) then

               ! ! ! Save output oringinal morrisson to compare with NN ===========================

               ! CALL RunMicro_HugMorr( &
               ! iMax                                 , &!INTEGER      , INTENT(IN   ) :: nCols
               ! kMax                                 , &!INTEGER      , INTENT(IN   ) :: kMax
               ! si          (1:kMax+1)               , &
               ! sl          (1:kMax)                 , &
               ! t3          (1:iMax,1:kMax)          , &!REAL(KIND=r8), INTENT(INOUT) :: Tc (1:nCols, 1:kMax)
               ! q3          (1:iMax,1:kMax)          , &!REAL(KIND=r8), INTENT(INOUT) :: qv (1:nCols, 1:kMax)
               ! ql3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qc (1:nCols, 1:kMax)
               ! qr3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qr (1:nCols, 1:kMax)
               ! qi3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qi (1:nCols, 1:kMax)
               ! QS3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qs (1:nCols, 1:kMax)
               ! QG3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qg (1:nCols, 1:kMax)
               ! NI3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: ni (1:nCols, 1:kMax)
               ! NS3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: ns (1:nCols, 1:kMax)
               ! NR3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: nr (1:nCols, 1:kMax)
               ! NG3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: NG (1:nCols, 1:kMax)
               ! NC3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: NC (1:nCols, 1:kMax)
               ! TKE         (1:iMax,1:kMax)          , &!REAL(KIND=r8), INTENT(IN   ) :: TKE (1:nCols, 1:kMax) (m^2 s-2)
               ! PBL_CoefKh  (1:iMax,1:kMax)          , &!REAL(KIND=r8), INTENT(IN   ) :: KZH (1:nCols, 1:kMax) (M^2 S-1)
               ! ps2         (1:iMax)*1000_r8         , &! gps- AIR PRESSURE (PA)
               ! DT                                   , &!REAL(KIND=r8), INTENT(IN   ) :: dt_in
               ! w2          (1:iMax,1:kMax)          , &!REAL(KIND=r8), INTENT(IN   ) :: omega  ! omega (Pa/s)
               ! EFFCS       (1:iMax,1:kMax)          , &!REAL(KIND=r8), INTENT(OUT  ) :: EFFCS (1:nCols, 1:kMax)   ! EFFCS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
               ! EFFIS       (1:iMax,1:kMax)          , &!REAL(KIND=r8), INTENT(OUT  ) :: EFFIS (1:nCols, 1:kMax)   ! EFFIS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
               ! RAINNCV    (1:iMax)                  , &!REAL(KIND=r8), INTENT(OUT) :: LSRAIN(1:nCols)
               ! SNOWNCV    (1:iMax)                  , &!REAL(KIND=r8), INTENT(OUT) :: LSSNOW(1:nCols)
               ! jdt*dt)

               ! do k=1, kMax
               !    call csv_write(111111, k, .false.)
               !    call csv_write(111111, t3(1,k), .false.)
               !    call csv_write(111111, q3(1,k), .false.)
               !    call csv_write(111111, ql3(1, k, latco), .false.)
               !    call csv_write(111111, qr3(1, k, latco), .false.)
               !    call csv_write(111111, qi3(1, k, latco), .false.)
               !    call csv_write(111111, QS3(1, k, latco), .false.)
               !    call csv_write(111111, QG3(1, k, latco), .false.)
               !    call csv_write(111111, NI3(1, k, latco), .false.)
               !    call csv_write(111111, NS3(1, k, latco), .false.)
               !    call csv_write(111111, NR3(1, k, latco), .false.)
               !    call csv_write(111111, NG3(1, k, latco), .false.)
               !    call csv_write(111111, NC3(1, k, latco), .false.)
               !    call csv_write(111111, EFFCS(1,k), .false.)
               !    call csv_write(111111, EFFIS(1,k), .false.)
               !    call csv_write(111111, RAINNCV(1:1), .false.)
               !    call csv_write(111111, SNOWNCV(1:1), .true.)
               ! end do
            
            ! OR NN ================

            print*, "calling NN -----------------------------------", reached, passed
            CALL RunMicro_HugMorr_NN( &
            iMax                                 , &!INTEGER      , INTENT(IN   ) :: nCols
            kMax                                 , &!INTEGER      , INTENT(IN   ) :: kMax
            sl          (1:kMax)               , &
            t3          (1:iMax,1:kMax)          , &!REAL(KIND=r8), INTENT(INOUT) :: Tc (1:nCols, 1:kMax)
            q3          (1:iMax,1:kMax)          , &!REAL(KIND=r8), INTENT(INOUT) :: qv (1:nCols, 1:kMax)
            ql3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qc (1:nCols, 1:kMax)
            qr3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qr (1:nCols, 1:kMax)
            qi3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qi (1:nCols, 1:kMax)
            QS3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qs (1:nCols, 1:kMax)
            QG3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qg (1:nCols, 1:kMax)
            NI3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: ni (1:nCols, 1:kMax)
            NS3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: ns (1:nCols, 1:kMax)
            NR3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: nr (1:nCols, 1:kMax)
            NG3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: NG (1:nCols, 1:kMax)
            NC3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: NC (1:nCols, 1:kMax)
            TKE         (1:iMax,1:kMax)          , &!REAL(KIND=r8), INTENT(IN   ) :: KZH (1:nCols, 1:kMax) (M^2 S-1)
            w2          (1:iMax,1:kMax)          , &!REAL(KIND=r8), INTENT(IN   ) :: omega  ! omega (Pa/s)
            RAINNCV    (1:iMax)                  , &!REAL(KIND=r8), INTENT(OUT) :: LSRAIN(1:nCols)
            SNOWNCV    (1:iMax)                    )!REAL(KIND=r8), INTENT(OUT) :: LSSNOW(1:nCols)

            print*, "Saving output to compare -----------------------------------", reached, passed
            do k=1, kMax
               call csv_write(222222, k, .false.)
               call csv_write(222222, t3(1,k), .false.)
               call csv_write(222222, q3(1,k), .false.)
               call csv_write(222222, ql3(1, k, latco), .false.)
               call csv_write(222222, qr3(1, k, latco), .false.)
               call csv_write(222222, qi3(1, k, latco), .false.)
               call csv_write(222222, QS3(1, k, latco), .false.)
               call csv_write(222222, QG3(1, k, latco), .false.)
               call csv_write(222222, NI3(1, k, latco), .false.)
               call csv_write(222222, NS3(1, k, latco), .false.)
               call csv_write(222222, NR3(1, k, latco), .false.)
               call csv_write(222222, NG3(1, k, latco), .false.)
               call csv_write(222222, NC3(1, k, latco), .false.)
               call csv_write(222222, EFFCS(1,k), .false.)
               call csv_write(222222, EFFIS(1,k), .false.)
               call csv_write(222222, RAINNCV(1:1), .false.)
               call csv_write(222222, SNOWNCV(1:1), .true.)
            end do
         end if
      else if(TRIM(ILCON).EQ.'HUMO') THEN
           CALL RunMicro_HugMorr( &
           iMax                                 , &!INTEGER      , INTENT(IN   ) :: nCols
           kMax                                 , &!INTEGER      , INTENT(IN   ) :: kMax
           si          (1:kMax+1)               , &
           sl          (1:kMax)                 , &
           t3          (1:iMax,1:kMax)          , &!REAL(KIND=r8), INTENT(INOUT) :: Tc (1:nCols, 1:kMax)
           q3          (1:iMax,1:kMax)          , &!REAL(KIND=r8), INTENT(INOUT) :: qv (1:nCols, 1:kMax)
           ql3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qc (1:nCols, 1:kMax)
           qr3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qr (1:nCols, 1:kMax)
           qi3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qi (1:nCols, 1:kMax)
           QS3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qs (1:nCols, 1:kMax)
           QG3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qg (1:nCols, 1:kMax)
           NI3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: ni (1:nCols, 1:kMax)
           NS3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: ns (1:nCols, 1:kMax)
           NR3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: nr (1:nCols, 1:kMax)
           NG3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: NG (1:nCols, 1:kMax)
           NC3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: NC (1:nCols, 1:kMax)
           TKE         (1:iMax,1:kMax)          , &!REAL(KIND=r8), INTENT(IN   ) :: TKE (1:nCols, 1:kMax) (m^2 s-2)
           PBL_CoefKh  (1:iMax,1:kMax)          , &!REAL(KIND=r8), INTENT(IN   ) :: KZH (1:nCols, 1:kMax) (M^2 S-1)
           ps2         (1:iMax)*1000_r8         , &! gps- AIR PRESSURE (PA)
           DT                                   , &!REAL(KIND=r8), INTENT(IN   ) :: dt_in
           w2          (1:iMax,1:kMax)          , &!REAL(KIND=r8), INTENT(IN   ) :: omega  ! omega (Pa/s)
           EFFCS       (1:iMax,1:kMax)          , &!REAL(KIND=r8), INTENT(OUT  ) :: EFFCS (1:nCols, 1:kMax)   ! EFFCS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
           EFFIS       (1:iMax,1:kMax)          , &!REAL(KIND=r8), INTENT(OUT  ) :: EFFIS (1:nCols, 1:kMax)   ! EFFIS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
           RAINNCV    (1:iMax)                  , &!REAL(KIND=r8), INTENT(OUT) :: LSRAIN(1:nCols)
           SNOWNCV    (1:iMax)                  , &!REAL(KIND=r8), INTENT(OUT) :: LSSNOW(1:nCols)
           jdt*dt)

      end if
!--------------------------------------------
        Total_Rain=Total_Rain + RAINNCV
        Total_Snow=Total_Snow + SNOWNCV
        snow_str=0.0_r8
        prec_sed=0.0_r8
        snow_sed=0.0_r8
        prec_pcw=0.0_r8
        snow_pcw=0.0_r8
       CALL qnegat2 (q3, fdqn,(1.0_r8/dt), del, iMax, kMax)! time t+1

    END IF

!    DO i=1,iMax
!       PRINT*,TRIM(ILCON),'   ',RAINNCV(i),SNOWNCV(i)
!    END DO

    ! Save humd/temp after large scale convection
    DO k=1,kMax
      DO i=1,iMax
         q3  (i,k)       = MAX(q3  (i,k)      ,1.0e-12_r8)
         IF (microphys) THEN
            ql3 (i,k,latco) = MAX(ql3 (i,k,latco),1.0e-12_r8)
            qi3 (i,k,latco) = MAX(qi3 (i,k,latco),1.0e-12_r8)

            ql2 (i,k,latco) = MAX(ql2 (i,k,latco),1.0e-12_r8)
            qi2 (i,k,latco) = MAX(qi2 (i,k,latco),1.0e-12_r8)
            IF((nClass+nAeros)>0 .and. PRESENT(gvarm))THEN
               qr3 (i,k,latco) = MAX(qr3 (i,k,latco),1.0e-12_r8)
               IF(TRIM(ILCON).EQ.'MORR'.or.TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HUMN')THEN
                  qs3 (i,k,latco) = MAX(qs3 (i,k,latco),1.0e-12_r8)
                  qg3 (i,k,latco) = MAX(qg3 (i,k,latco),1.0e-12_r8)
                  NI3 (i,k,latco) = MAX(NI3 (i,k,latco),0.0e-12_r8)
                  NS3 (i,k,latco) = MAX(NS3 (i,k,latco),0.0e-12_r8)
                  NR3 (i,k,latco) = MAX(NR3 (i,k,latco),0.0e-12_r8)
                  NG3 (i,k,latco) = MAX(NG3 (i,k,latco),0.0e-12_r8)
                  IF(TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HUMN') NC3 (i,k,latco) = MAX(NC3 (i,k,latco),0.0e-12_r8)

               END IF
               qr2 (i,k,latco) = MAX(qr2 (i,k,latco),1.0e-12_r8)
               IF(TRIM(ILCON).EQ.'MORR'.or.TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HUMN')THEN
                  qs2 (i,k,latco) = MAX(qs2 (i,k,latco),1.0e-12_r8)
                  qg2 (i,k,latco) = MAX(qg2 (i,k,latco),1.0e-12_r8)
                  NI2 (i,k,latco) = MAX(NI2 (i,k,latco),0.0e-12_r8)
                  NS2 (i,k,latco) = MAX(NS2 (i,k,latco),0.0e-12_r8)
                  NR2 (i,k,latco) = MAX(NR2 (i,k,latco),0.0e-12_r8)
                  NG2 (i,k,latco) = MAX(NG2 (i,k,latco),0.0e-12_r8)
                  IF(TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HUMN')NC2 (i,k,latco) = MAX(NC2 (i,k,latco),0.0e-12_r8)
               END IF
            END IF
         END IF
         tLrgs(i,k)=t3(i,k)
         qLrgs(i,k)=q3(i,k)
      END DO
    END DO

    IF (microphys) THEN
       DO k=1,kMax
          DO i=1,iMax
              gliqp(i,k) =ql3     (i,k,latco) 
              gicep(i,k) =qi3     (i,k,latco) 

              gliqm(i,k)   =ql2   (i,k,latco) 
              gicem(i,k)   =qi2   (i,k,latco) 
              IF((nClass+nAeros)>0 .and. PRESENT(gvarm))THEN
                 gvarp(i,k,1) =qr3   (i,k,latco) 
                 IF(TRIM(ILCON).EQ.'MORR'.or.TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HUMN')THEN
                    gvarp(i,k,2) =qs3   (i,k,latco) 
                    gvarp(i,k,3) =qg3   (i,k,latco) 
                    gvarp(i,k,4) =NI3   (i,k,latco) 
                    gvarp(i,k,5) =NS3   (i,k,latco) 
                    gvarp(i,k,6) =NR3   (i,k,latco) 
                    gvarp(i,k,7) =NG3   (i,k,latco) 
                    IF(TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HUMN')gvarp(i,k,8) =NC3   (i,k,latco) 

                 END IF
                 gvarm(i,k,1) =qr2   (i,k,latco) 
                 IF(TRIM(ILCON).EQ.'MORR'.or.TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HUMN')THEN
                    gvarm(i,k,2) =qs2   (i,k,latco) 
                    gvarm(i,k,3) =qg2   (i,k,latco) 
                    gvarm(i,k,4) =NI2   (i,k,latco) 
                    gvarm(i,k,5) =NS2   (i,k,latco) 
                    gvarm(i,k,6) =NR2   (i,k,latco) 
                    gvarm(i,k,7) =NG2   (i,k,latco) 
                    IF(TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HUMN')gvarm(i,k,8) =NC2   (i,k,latco) 

                 END IF 
              END IF
          END DO
       END DO   
    END IF


  END SUBROUTINE RunMicroPhysics

!**********************************************************************************
!**********************************************************************************

  SUBROUTINE  qnegat2 (fq, fdq, rdt, del, iMax, kMax)
    !
    ! input: fq  specific humidity (dimensionless mixing ratio)
    !        fp  surface pressure (cb)
    ! ouput: fq  adjusted specific humidity
    !        fp  unchanged
    !        fdq distribution of moisture modification
    !
    ! iMax......Number of grid points on a gaussian latitude circle   
    ! kMax......Number of sigma levels  
    ! imx.......=iMax+1 or iMax+2   :this dimension instead of iMax
    !              is used in order to avoid bank conflict of memory
    !              access in fft computation and make it efficient. the
    !              choice of 1 or 2 depends on the number of banks and
    !              the declared type of grid variable (real*4,real*8)
    !              to be fourier transformed.
    !              cyber machine has the symptom.
    !              cray machine has no bank conflict, but the argument
    !              'imx' in subr. fft991 cannot be replaced by iMax    
    ! del.......sigma spacing for each layer computed in routine "setsig".  
    ! dfact.....del(k+1)/del(k)
    !
    INTEGER, INTENT(in   ) :: iMax  
    INTEGER, INTENT(in   ) :: kMax
    REAL(KIND=r8)   , INTENT(in   ) :: rdt

    REAL(KIND=r8),    INTENT(inout) :: fq   (iMax,kMax)
    REAL(KIND=r8),    INTENT(inout) :: fdq  (iMax,kMax)  
    REAL(KIND=r8),    INTENT(in   ) :: del  (kMax)

    REAL(KIND=r8)   :: dfact(kMax)
    REAL(KIND=r8)   :: rdt2
    INTEGER :: klev
    INTEGER :: kblw
    INTEGER :: i
    INTEGER :: k  
    rdt2=rdt
    DO k=1,kMax-1
       dfact(k+1) = del(k+1)/del(k)
    END DO
    !     
    !     ecmwf vertical borrowing scheme
    !     fdq contains compensated borrowing above first level, uncompensated
    !     borrowing in first level
    !     
    DO k=1,kMax-1
       klev = kMax-k+1
       kblw = klev - 1
       DO i=1,iMax
          fdq(i,klev) = fq(i,klev)
          IF(fq(i,klev).LT.0.0e0_r8) fq(i,klev) = 1.0e-12_r8
          fdq(i,klev) = fq(i,klev) - fdq(i,klev)
          fq(i,kblw) = fq(i,kblw) - fdq(i,klev)*dfact(klev)
       END DO
    END DO

    DO i=1,iMax
       fdq(i,1) = fq(i,1)
       IF(fq(i,1).LT.0.0e0_r8) fq(i,1) = 1.0e-12_r8
       fdq(i,1) = fq(i,1) - fdq(i,1)
    END DO

  END SUBROUTINE qnegat2

!**********************************************************************************
!**********************************************************************************

  SUBROUTINE FinalizeMicroPhysics()
    IMPLICIT NONE
    DEALLOCATE(ql3)
    DEALLOCATE(qi3)
    DEALLOCATE(qr3)
    DEALLOCATE(qs3)
    DEALLOCATE(qg3)
    DEALLOCATE(NI3)
    DEALLOCATE(NS3)
    DEALLOCATE(NR3)
    DEALLOCATE(NG3)
    DEALLOCATE(NC3)
    
    DEALLOCATE(ql2)
    DEALLOCATE(qi2)
    DEALLOCATE(qr2)
    DEALLOCATE(qs2)
    DEALLOCATE(qg2)
    DEALLOCATE(NI2)
    DEALLOCATE(NS2)
    DEALLOCATE(NR2)
    DEALLOCATE(NG2)
    DEALLOCATE(NC2)

  END SUBROUTINE FinalizeMicroPhysics

END MODULE MicroPhysics
