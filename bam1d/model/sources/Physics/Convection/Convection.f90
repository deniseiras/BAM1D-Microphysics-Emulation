!
!  $Author: pkubota $
!  $Date: 2008/04/09 12:42:57 $
!  $Revision: 1.11 $
!
MODULE Convection

  !   InitConvection
  !
  !   cumulus_driver|--qnegat
  !                 |--Cu_Ara
  !                 |--Cu_Kuolcl
  !                 |--Cu_Grellens
  !                 |--Shall_Tied
  !                 |--Shall_Souza
  !                 |--lrgscl
  !                 


  USE Constants, ONLY :  &
       delq             ,&
       r8,i8,qmin,grav,gasr
       

  USE Diagnostics, ONLY:   &
       dodia             , &
       updia             , &
       StartStorDiag     , &
       nDiag_toprec      , & ! total precipiation
       nDiag_cvprec      , & ! convective precipitation
       nDiag_lsprec      , & ! large scale precipitation
       nDiag_snowfl      , & ! snowfall
       nDiag_clheat      , & ! convective latent heating
       nDiag_cmchan      , & ! convective moisture change
       nDiag_lslhea      , & ! large scale latent heating
       nDiag_lsmcha      , & ! large scale moisture change
       nDiag_sclhea      , & ! shallow convective latent heating
       nDiag_scmcha      , & ! shallow convective moisture change
       nDiag_nshcrm      , & ! negative specific humidity correction moisture source
       nDiag_qlicld      , & ! liquid water content in cloud after rainout
       nDiag_trcliq      , & ! Water Liquid Cloud kg/kg
       nDiag_trcice      , & ! Water Ice Cloud kg/kg
       nDiag_cape2d      , & ! CONVECTIVE AVAIL. POT.ENERGY M2/S2 
       nDiag_cine2d      , & ! CONVECTIVE INHIB. ENERGY M2/S2
       nDiag_sweath  !    , & ! SEVERE WEATHER THREAT
!       nDiag_qrmicr      , & ! = 230 ! QR - rain water mixing ratio  (kg/kg)
!       nDiag_qsmicr      , & ! = 231 ! QS - snow mixing ratio (kg/kg)
!       nDiag_qgmicr      , & ! = 232 ! QG - graupel mixing ratio (KG/KG)
!       nDiag_nimicr      , & ! = 233 ! NI - cloud ice number concentration (1/kg)
!       nDiag_nsmicr      , & ! = 234 ! NS - Snow Number concentration (1/kg)
!       nDiag_ncmicr      , & ! = 235 ! NC - Cloud droplet Number concentration (1/kg)
!       nDiag_nrmicr      , & ! = 236 ! NR - Rain Number concentration (1/kg)
!       nDiag_ngmicr      , & ! = 237 ! NG - Graupel number concentration (1/kg)
!       nDiag_co2aer          ! mixC02 - CO2 concentration (kg/kg)

  USE GridHistory, ONLY:   &
       IsGridHistoryOn   , &
       StoreGridHistory  , &
       dogrh             , &
       nGHis_cvprec     , &
       nGHis_clheat     , &
       nGHis_cvmosr     , &
       nGHis_sclhea     , &
       nGHis_shcvmo     , &
       nGHis_toprec     , &
       nGHis_snowfl     , &
       nGHis_sslaht     , &
       nGHis_spstms

  USE Options, ONLY :       &
       rccmbl            , &
       mlrg              , &
       iccon             , &
       ilcon             , &
       iscon             , &
       doprec            , &
       cflric            , &
       dt                , &
       kt                , & 
       ktp               , &  
       ktm               , & 
       jdt               , & 
       nfcnv0            , & 
       isimp             , & 
       nfctrl            , & 
       nfprt             , &
       microphys,nClass,nAeros

  USE wv_saturation, ONLY :       &
      gestbl,&
      findsp

 USE FieldsPhysics, ONLY:  &
       convc             , &
       convt             , &
       convb             , &
       prcp1             , &
       prcp2             , &
       prcp3             , &
       prcpt             , &
       toplv             , &
       botlv             , &
       taud              , &
       pblh              , &
       cu_hr             , &
       cu_kbot           , &
       cu_ktop           , &
       cu_Kuo            , &
       LOWLYR            , &
       F_ICE_PHY         , &
       F_RAIN_PHY        , &
       F_RIMEF_PHY       , &
       EFFCS             , &
       EFFIS             , &
       rVisDiff          , &
       rVisBeam          , &
       rNirDiff          , &
       rNirBeam          , &
       rVisDiffC         , &
       rVisBeamC         , &
       rNirDiffC         , &
       rNirBeamC         , &
       rSwToaDown        , &
       Dump

  USE Init, ONLY :       &
       nls

  USE Parallelism, ONLY: &
       MsgOne, FatalError

  USE DeepConvection, ONLY: &
      InitDeepConvection,RunDeepConvection,FinalizeDeepConvection

  USE MicroPhysics, ONLY: &
      InitMicroPhysics,RunMicroPhysics,FinalizeMicroPhysics
  
  
  USE ShallowConvection, ONLY: &
      InitShallowConvection,RunShallowConvection,FinalizeShallowConvection


  USE PhysicalFunctions, ONLY: calc_cape,SWEAT_index

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitConvection
  PUBLIC :: cumulus_driver
  PUBLIC :: FinalizeConvection

  PUBLIC :: InitCheckFileConvec
  PUBLIC :: ReStartConvec
  INTEGER, PARAMETER :: ppcnst=3
  REAL(KIND=r8), ALLOCATABLE  :: ql2(:,:,:)
  REAL(KIND=r8), ALLOCATABLE  :: qi2(:,:,:)
  REAL(KIND=r8), ALLOCATABLE  :: ql3(:,:,:)
  REAL(KIND=r8), ALLOCATABLE  :: qi3(:,:,:)

CONTAINS

  SUBROUTINE InitConvection(si,del        ,sl         ,cl         , &
                          kMax    ,iMax,jMax,ibMax,jbMax,&
                          trunc   ,ifdy       ,todcld     ,ids        , &
                          idc     ,ifday      ,tod           ,fNameMicro ,idate         )
      
    INTEGER, INTENT(IN) :: kMax
    INTEGER, INTENT(IN) :: iMax
    INTEGER, INTENT(IN) :: jMax
    INTEGER, INTENT(IN) :: ibMax
    INTEGER, INTENT(IN) :: jbMax
    INTEGER, INTENT(IN) :: trunc   
    REAL(KIND=r8),    INTENT(IN   ) :: si (kMax+1)
    REAL(KIND=r8),    INTENT(IN   ) :: del(kMax  )
    REAL(KIND=r8),    INTENT(IN   ) :: sl (kMax  )
    REAL(KIND=r8),    INTENT(IN   ) :: cl (kMax  )
    INTEGER         , INTENT(OUT  ) :: ifdy
    REAL(KIND=r8)   , INTENT(OUT  ) :: todcld
    INTEGER         , INTENT(OUT  ) :: ids(4)
    INTEGER         , INTENT(OUT  ) :: idc(4)
    INTEGER         , INTENT(IN   ) :: ifday
    REAL(KIND=r8)   , INTENT(IN   ) :: tod
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameMicro
    INTEGER         , INTENT(IN   ) :: idate(4)
    CHARACTER(LEN=*), PARAMETER     :: h='**(InitConvection)**'
    LOGICAL                         :: restart
    REAL(KIND=r8)    :: fhour

    IF(nfcnv0.NE.0) THEN
       CALL MsgOne(h,'Reading previous physics state for restart')
       restart=.TRUE.
    ELSE
       CALL MsgOne(h,'Initializing prec/cloud variables')
       restart=.FALSE.
    END IF

    CALL gestbl()
   fhour=ifdy*24.0_r8+tod/3600.0_r8

    CALL InitDeepConvection( &
                          dt      ,si   ,del  ,sl   , &
                          kMax    ,iMax ,jMax ,ibMax,jbMax, &
                          fhour   ,idate,iccon)

    CALL InitShallowConvection( &
                          si      ,del ,sl   ,cl   , &
                          kMax    ,jMax,ibMax,jbMax, &
                          trunc   ,ppcnst  ,ISCON        )

    CALL InitMicroPhysics( &
                          dt      ,si        ,del      ,sl       , &
                          kMax    ,iMax      ,jMax     ,ibMax    , &
                          jbMax   ,fNameMicro,ILCON   ,microphys, &
                          nClass  ,nAeros    ,EFFCS     ,EFFIS)
   
    ALLOCATE(ql3(ibMax,kMax,jbMax));ql3=0.00001e-12_r8
    ALLOCATE(qi3(ibMax,kMax,jbMax));qi3=0.00001e-12_r8
    
    ALLOCATE(ql2(ibMax,kMax,jbMax));ql2=0.00001e-12_r8
    ALLOCATE(qi2(ibMax,kMax,jbMax));qi2=0.00001e-12_r8

    IF(TRIM(isimp).NE.'YES') THEN
       CALL InitBoundCondConvec(&
           ifdy,todcld,ids,idc,ifday, &
           tod)
    END IF   

  END SUBROUTINE InitConvection

  SUBROUTINE InitBoundCondConvec(&
       ifdy,todcld,ids,idc,ifday, &
       tod)

    INTEGER         , INTENT(OUT  ) :: ifdy
    REAL(KIND=r8)   , INTENT(OUT  ) :: todcld
    INTEGER         , INTENT(OUT  ) :: ids(4)
    INTEGER         , INTENT(OUT  ) :: idc(4)
    INTEGER         , INTENT(IN   ) :: ifday
    REAL(KIND=r8)   , INTENT(IN   ) :: tod
    CHARACTER(LEN=*), PARAMETER :: h='**(InitBoundCondConvec)**'

    IF(nfcnv0.NE.0) THEN
       CALL MsgOne(h,'Reading previous physics state for restart')

       READ(UNIT=nfcnv0) ifdy,todcld,ids,idc
       READ(UNIT=nfcnv0) convc,convt,convb,prcp1,prcp2,prcp3, &
        prcpt,toplv,botlv,taud,pblh,cu_hr ,cu_kbot,cu_ktop,cu_Kuo 
       
       READ(UNIT=nfcnv0) F_ICE_PHY, F_RAIN_PHY, F_RIMEF_PHY,EFFCS,EFFIS
       
       IF(ifday.GT.0.OR.tod.GT.0.0_r8)READ(UNIT=nfcnv0)rVisDiff,rVisBeam,rNirDiff, &
        rNirBeam,rVisDiffC,rVisBeamC,rNirDiffC,rNirBeamC,rSwToaDown

       REWIND nfcnv0

       IF(nfctrl(4) .GE. 1)WRITE(UNIT=nfprt,FMT=555)ifdy,todcld,ids,idc
    ELSE
       CALL MsgOne(h,'Initializing prec/cloud variables')

       convc=0.0_r8
       convt=0.0_r8
       convb=0.0_r8
       prcp1=0.0_r8
       prcp2=0.0_r8
       prcp3=0.0_r8
       prcpt=0.0_r8
       toplv=0.0_r8
       botlv=0.0_r8
    END IF


555 FORMAT(' CLOUD PROGNOSTIC DATA READ IN. AT FORECAST DAY', &
         I8,' TOD ',F8.1/' STARTING',3I3,I5,' CURRENT',3I3,I5)
  END SUBROUTINE InitBoundCondConvec



  SUBROUTINE InitCheckFileConvec(ifdy  ,todcld,ids   ,idc    )
    INTEGER      , INTENT(OUT  ) :: ifdy
    REAL(KIND=r8), INTENT(OUT  ) :: todcld
    INTEGER      , INTENT(OUT  ) :: ids   (4)
    INTEGER      , INTENT(OUT  ) :: idc   (4)
    CHARACTER(LEN=*), PARAMETER :: h="**(InitCheckFileConvec)**"
    
    !
    !     read cloud dataset for cold start
    !
    IF(nfcnv0.NE.0) THEN
       CALL MsgOne(h,'Read prec/cloud variables')    
       READ(UNIT=nfcnv0) ifdy,todcld,ids,idc
       READ(UNIT=nfcnv0) convc,convt,convb,prcp1,prcp2,prcp3, &
            prcpt,toplv,botlv,taud,pblh,cu_hr ,cu_kbot,cu_ktop,cu_Kuo 
       READ(UNIT=nfcnv0) F_ICE_PHY, F_RAIN_PHY, F_RIMEF_PHY ,EFFCS,EFFIS

       REWIND nfcnv0

       IF(nfctrl(4) .GE. 1) WRITE(UNIT=nfprt,FMT=555)ifdy,todcld,ids,idc

    ELSE
       CALL MsgOne(h,'Initializing prec/cloud variables')
       convc=0.0_r8
       convt=0.0_r8
       convb=0.0_r8
       prcp1=0.0_r8
       prcp2=0.0_r8
       prcp3=0.0_r8
       prcpt=0.0_r8
       toplv=0.0_r8
       botlv=0.0_r8
    END IF
555 FORMAT(' CLOUD PROGNOSTIC DATA READ IN. AT FORECAST DAY', &
         I8,' TOD ',F8.1/' STARTING',3I3,I5,' CURRENT',3I3,I5)

  END SUBROUTINE InitCheckFileConvec
  
   SUBROUTINE ReStartConvec (ifday,tod,idate ,idatec,nfcnv1)

    INTEGER           ,INTENT(IN   ) :: ifday
    REAL(KIND=r8)     ,INTENT(IN   ) :: tod
    INTEGER           ,INTENT(IN   ) :: idate(4)
    INTEGER           ,INTENT(IN   ) :: idatec(4)
    INTEGER           ,INTENT(IN   ) :: nfcnv1

    IF(TRIM(isimp).NE.'YES') THEN
       CALL MsgOne('**(ReStartConvec)**','Saving physics state for restart')

       !$OMP SINGLE
       WRITE(UNIT=nfcnv1) ifday,tod,idate,idatec
       WRITE(UNIT=nfcnv1) convc,convt,convb,prcp1,prcp2,prcp3, &
            prcpt,toplv,botlv,taud,pblh,cu_hr  ,cu_kbot,cu_ktop,cu_Kuo 
       WRITE(UNIT=nfcnv1) F_ICE_PHY, F_RAIN_PHY, F_RIMEF_PHY,EFFCS,EFFIS 
       WRITE(UNIT=nfcnv1) rVisDiff,rVisBeam,rNirDiff,rNirBeam, &
            rVisDiffC,rVisBeamC,rNirDiffC,rNirBeamC,rSwToaDown
       !$OMP END SINGLE
    END IF

  END SUBROUTINE ReStartConvec
  
  
  SUBROUTINE cumulus_driver (&
      ! Run Flags
      ! Time info
            fac       , & !REAL(KIND=r8), INTENT(IN   ) :: fac
            fac2      , & !REAL(KIND=r8), INTENT(IN   ) :: fac2
            fac2x     , & !REAL(KIND=r8), INTENT(IN   ) :: fac2x  
      ! Model Geometry
            colrad    , & !REAL(KIND=r8), INTENT(IN   ) :: colrad(iMax)
            lonrad    , & !REAL(KIND=r8), INTENT(IN   ) :: lonrad(iMax)  
            del       , & !REAL(KIND=r8), INTENT(IN   ) :: del  (kMax)
            sl        , & !REAL(KIND=r8), INTENT(IN   ) :: sl   (kMax)
            si        , & !REAL(KIND=r8), INTENT(IN   ) :: si   (kMax+1)
      ! Model information
            iMax      , & !INTEGER      , INTENT(IN   ) :: iMax
            kMax      , & !INTEGER      , INTENT(IN   ) :: kMax
            latco     , & !INTEGER      , INTENT(IN   ) :: latco
            mask      , & !INTEGER(KIND=i8),INTENT(IN ) :: mask (iMax) 
            zs        , & !REAL(KIND=r8), INTENT(IN   ) :: zs   (iMax)
      ! CONVECTION: convective clouds
            convc     , & !REAL(KIND=r8), INTENT(INOUT) :: convc  (iMax)
            convt     , & !REAL(KIND=r8), INTENT(INOUT) :: convt  (iMax)
            convb     , & !REAL(KIND=r8), INTENT(INOUT) :: convb  (iMax)
            toplv     , & !REAL(KIND=r8), INTENT(INOUT) :: toplv  (iMax)
            botlv     , & !REAL(KIND=r8), INTENT(INOUT) :: botlv  (iMax)
            convts    , & !REAL(KIND=r8), INTENT(INOUT) :: convts (iMax)
            convcs    , & !REAL(KIND=r8), INTENT(INOUT) :: convcs (iMax)
            convbs    , & !REAL(KIND=r8), INTENT(INOUT) :: convbs (iMax)
            concld    , & !REAL(KIND=r8), INTENT(INOUT) :: concld (iMax,kMax)
            cld       , & !REAL(KIND=r8), INTENT(INOUT) :: cld    (iMax,kMax)
      ! SURFACE:  Fields
            tsfc      , & !REAL(KIND=r8), INTENT(IN   ) :: tsfc (iMax)
            tpert     , & !REAL(KIND=r8), INTENT(IN   ) :: tpert(iMax)
            qpert     , & !REAL(KIND=r8), INTENT(IN   ) :: qpert(iMax)
            sens      , & !REAL(KIND=r8), INTENT(IN   ) :: sens (iMax)
            evap      , & !REAL(KIND=r8), INTENT(IN   ) :: evap (iMax)
            ustar     , & !REAL(KIND=r8), INTENT(IN   ) :: ustar (iMax)
      ! Precipitation Field
            prcp1     , & !REAL(KIND=r8), INTENT(INOUT) :: prcp1  (iMax)
            prcp2     , & !REAL(KIND=r8), INTENT(INOUT) :: prcp2  (iMax)
            prcp3     , & !REAL(KIND=r8), INTENT(INOUT) :: prcp3  (iMax)
            prcpt     , & !REAL(KIND=r8), INTENT(INOUT) :: prcpt  (iMax)
            geshem    , & !REAL(KIND=r8), INTENT(INOUT) :: geshem (iMax)
            ppli      , & !REAL(KIND=r8), INTENT(INOUT) :: ppli   (iMax)
            ppci      , & !REAL(KIND=r8), INTENT(INOUT) :: ppci   (iMax)
            prct      , & !REAL(KIND=r8), INTENT(INOUT) :: prct   (iMax)
            prcc      , & !REAL(KIND=r8), INTENT(INOUT) :: prcc   (iMax)
            snowfl    , & !REAL(KIND=r8), INTENT(INOUT) :: snowfl (iMax)
      ! PBL:  Fields
            PBL_CoefKh, & !REAL(KIND=r8), INTENT(IN   ) :: PBL_CoefKh(iMax,kMax) 
            tke       , & !REAL(KIND=r8), INTENT(IN   ) :: tke  (iMax,kMax)
      ! Microphysics
            dudt      , &
            dvdt      , &
            qliq      , & !REAL(KIND=r8), INTENT(OUT  ) :: qliq(iMax,kMax)
            EFFCS     , & !REAL(KIND=r8), INTENT(INOUT) :: EFFCS   (iMax,kMax)  
            EFFIS     , & !REAL(KIND=r8), INTENT(INOUT) :: EFFIS   (iMax,kMax)   
      ! Atmospheric fields
            ta        , & !REAL(KIND=r8), INTENT(INOUT) :: ta (iMax,kMax)
            tb        , & !REAL(KIND=r8), INTENT(INOUT) :: tb (iMax,kMax)
            tc        , & !REAL(KIND=r8), INTENT(INOUT) :: tc (iMax,kMax)
            qa        , & !REAL(KIND=r8), INTENT(INOUT) :: qa (iMax,kMax)
            qb        , & !REAL(KIND=r8), INTENT(INOUT) :: qb (iMax,kMax)
            qc        , & !REAL(KIND=r8), INTENT(INOUT) :: qc (iMax,kMax)
            ub        , & !REAL(KIND=r8), INTENT(IN   ) :: ub   (iMax,kMax) ! (m/s) 
            vb        , & !REAL(KIND=r8), INTENT(IN   ) :: vb   (iMax,kMax) ! (m/s)
            omgb      , & !REAL(KIND=r8), INTENT(IN   ) :: omgb (iMax,kMax) ! (Pa/s)
            psb       , & !REAL(KIND=r8), INTENT(IN   ) :: psb  (iMax)
            psb2      , & !REAL(KIND=r8), INTENT(IN   ) :: psb2 (iMax)
            gicep     , & !REAL(KIND=r8), INTENT(INOUT) :: gicep  (iMax,kmax)
            gicem     , & !REAL(KIND=r8), INTENT(INOUT) :: gicem  (iMax,kmax)
            gliqp     , & !REAL(KIND=r8), INTENT(INOUT) :: gliqp  (iMax,kmax)
            gliqm     , & !REAL(KIND=r8), INTENT(INOUT) :: gliqm  (iMax,kmax)
            gvarp     , & !REAL(KIND=r8),OPTIONAL,   INTENT(INOUT) :: gvarp (iMax,kmax,nClass+nAeros)
            gvarm       ) !REAL(KIND=r8),OPTIONAL,   INTENT(INOUT) :: gvarm (iMax,kmax,nClass+nAeros)

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
    !    iccon                  ! cu schemes ex. KUO, ARA, GRE ..
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

    ! Dimensions
    INTEGER, INTENT(IN) :: iMax
    INTEGER, INTENT(IN) :: kMax

    ! Sizes
    REAL(KIND=r8), INTENT(IN) :: sl   (kMax)
    REAL(KIND=r8), INTENT(IN) :: del  (kMax)
    REAL(KIND=r8), INTENT(IN) :: si   (kMax+1)

    ! Fixed fields: latitudes, mask and topography
    INTEGER, INTENT(IN) :: latco
    INTEGER(KIND=i8), INTENT(IN) :: mask (iMax) 
    REAL(KIND=r8),    INTENT(IN) :: colrad(iMax)
    REAL(KIND=r8),    INTENT(IN) :: lonrad(iMax)  

    REAL(KIND=r8),    INTENT(IN) :: zs   (iMax)

    ! Temperature (K) and specific humidity (kg/kg) 
    ! at times (a) = T-1, (b) = T and (c) = T+1
    REAL(KIND=r8), INTENT(INOUT) :: ta (iMax,kMax)
    REAL(KIND=r8), INTENT(INOUT) :: tb (iMax,kMax)
    REAL(KIND=r8), INTENT(INOUT) :: tc (iMax,kMax)
    REAL(KIND=r8), INTENT(INOUT) :: qa (iMax,kMax)
    REAL(KIND=r8), INTENT(INOUT) :: qb (iMax,kMax)
    REAL(KIND=r8), INTENT(INOUT) :: qc (iMax,kMax)
    REAL(KIND=r8), INTENT(INOUT) :: dudt(iMax,kMax)
    REAL(KIND=r8), INTENT(INOUT) :: dvdt(iMax,kMax)

    REAL(KIND=r8), INTENT(OUT  ) :: qliq(iMax,kMax) !qrc liquid water content in cloud after rainout
    REAL(KIND=r8), INTENT(INOUT) :: EFFCS   (iMax,kMax)  
    REAL(KIND=r8), INTENT(INOUT) :: EFFIS   (iMax,kMax)   
    REAL(KIND=r8), INTENT(IN   ) :: PBL_CoefKh(iMax,kMax) 
    ! Wind at time T
    ! in half levels along with U,V,T,Q
    REAL(KIND=r8), INTENT(IN) :: ub   (iMax,kMax) ! (m/s) 
    REAL(KIND=r8), INTENT(IN) :: vb   (iMax,kMax) ! (m/s)
    REAL(KIND=r8), INTENT(IN) :: omgb (iMax,kMax) ! (Pa/s)

    ! Surface pressure (cb) at time T
    REAL(KIND=r8), INTENT(IN) :: psb  (iMax)
    REAL(KIND=r8), INTENT(IN) :: psb2 (iMax)
    REAL(KIND=r8), INTENT(IN) :: tsfc (iMax)
    REAL(KIND=r8), INTENT(IN) :: tpert(iMax)
    REAL(KIND=r8), INTENT(IN) :: qpert(iMax)
    REAL(KIND=r8), INTENT(IN) :: tke  (iMax,kMax)

    ! Heat/Water sfc fluxes
    REAL(KIND=r8), INTENT(IN) :: sens (iMax)
    REAL(KIND=r8), INTENT(IN) :: evap (iMax)
    REAL(KIND=r8), INTENT(IN) :: ustar(iMax)

    ! UNCLASSIFIED VARIABLES
    REAL(KIND=r8), INTENT(IN)      :: fac2x  
    REAL(KIND=r8), INTENT(IN)      :: fac
    REAL(KIND=r8), INTENT(IN)      :: fac2

    REAL(KIND=r8), INTENT(INOUT) :: convc  (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: convt  (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: convb  (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: prcp1  (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: prcp2  (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: prcp3  (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: prcpt  (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: toplv  (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: botlv  (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: convts (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: convcs (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: convbs (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: geshem (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: ppli   (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: ppci   (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: prcc   (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: snowfl (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: prct   (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: concld (iMax,kMax)
    REAL(KIND=r8), INTENT(INOUT) :: cld    (iMax,kMax)
    REAL(KIND=r8), INTENT(INOUT) :: gicem  (iMax,kmax)
    REAL(KIND=r8), INTENT(INOUT) :: gicep  (iMax,kmax)
    REAL(KIND=r8), INTENT(INOUT) :: gliqm  (iMax,kmax)
    REAL(KIND=r8), INTENT(INOUT) :: gliqp  (iMax,kmax)
    REAL(KIND=r8),OPTIONAL,   INTENT(INOUT) :: gvarm (iMax,kmax,nClass+nAeros)
    REAL(KIND=r8),OPTIONAL,   INTENT(INOUT) :: gvarp (iMax,kmax,nClass+nAeros)

    !
    !               L O C A L    V A R I A B L E S
    !              -------------------------------
    INTEGER       :: i
    INTEGER       :: k,n
    REAL(KIND=r8) :: ps1     (iMax)
    REAL(KIND=r8) :: ps2     (iMax)
    REAL(KIND=r8) :: PS_work (iMax)
    REAL(KIND=r8) :: terr    (iMax)  

    ! (T,q) before any convection
    REAL(KIND=r8) :: tBegin (iMax,kMax)
    REAL(KIND=r8) :: qBegin (iMax,kMax)

    ! (T,q) after deep convection
    REAL(KIND=r8) :: qDeep(iMax,kMax)
    REAL(KIND=r8) :: tDeep(iMax,kMax)
    ! (T,q) after shallow convection
    REAL(KIND=r8) :: tShal(iMax,kMax)
    REAL(KIND=r8) :: qShal(iMax,kMax)
    ! (T,q) after large scale adjustment
    REAL(KIND=r8) :: tLrgs(iMax,kMax)
    REAL(KIND=r8) :: qLrgs(iMax,kMax)

    ! Shallow heat/moist
    REAL(KIND=r8) :: sclhea(iMax,kMax)
    REAL(KIND=r8) :: scmcha(iMax,kMax)
    ! Convective heat/moist
    REAL(KIND=r8) :: clheat(iMax,kMax)
    REAL(KIND=r8) :: cmchan(iMax,kMax)
    ! Large scale heat/moist
    REAL(KIND=r8) :: lslhea(iMax,kMax)
    REAL(KIND=r8) :: lsmcha(iMax,kMax)

    ! Working copies of input fields (T,q) at times (a,b,c)
    REAL(KIND=r8) :: q1(iMax,kMax)
    REAL(KIND=r8) :: q2(iMax,kMax)
    REAL(KIND=r8) :: q3(iMax,kMax)


    REAL(KIND=r8) :: t1(iMax,kMax)
    REAL(KIND=r8) :: t2(iMax,kMax)
    REAL(KIND=r8) :: t3(iMax,kMax)

    ! Wind components for grell ensemble
    !REAL(KIND=r8) :: u2(iMax,kMax)
    !REAL(KIND=r8) :: v2(iMax,kMax)
    !REAL(KIND=r8) :: w2(iMax,kMax)

    REAL(KIND=r8) :: dlf          (iMax,kMax)    ! detrained water from ZM
    REAL(KIND=r8) :: rliq      (iMax)      ! vertical integral of liquid not yet in q(ixcldliq)
    real(KIND=r8) :: rliq2(iMax)                   ! vertical integral of liquid from shallow scheme
    REAL(KIND=r8) :: cmfmc  (iMax,kMax+1)   ! convective mass flux--m sub c
    REAL(KIND=r8) :: cmfmc2 (iMax,kMax+1)   ! shallow convective mass flux--m sub c
    REAL(KIND=r8) :: zdu          (iMax,kMax)   ! detrainment rate from deep convection

    !
    ! UNCLASSIFIED VARIABLES
    !
    REAL(KIND=r8) :: fdqn (iMax,kMax)
    REAL(KIND=r8) :: RAINCV     (iMax)
    REAL(KIND=r8) :: SNOWCV     (iMax)
    REAL(KIND=r8) :: Total_Rain (iMax)
    REAL(KIND=r8) :: Total_Snow (iMax)

    !*******************************************
    !               Ktopos nao usado fora
    !            kctop1  usado para ARA fora
    !*******************************************
    REAL(KIND=r8)    :: hrem  (iMax,kMax)
    REAL(KIND=r8)    :: qrem  (iMax,kMax)
    REAL(KIND=r8)    :: cldm  (iMax)
    INTEGER          :: kctop1(iMax)
    INTEGER          :: kcbot1(iMax)
    INTEGER          :: noshal(iMax)
    !**********
    !others
    !*********
    INTEGER          :: ktop (iMax)
    INTEGER          :: kuo  (iMax)
    INTEGER          :: ktops(iMax)
    REAL(KIND=r8)    :: plcl (iMax)
    INTEGER          :: kbot (iMax)
    INTEGER          :: kbots(iMax) 
    REAL(KIND=r8)    :: dq   (iMax,kMax)
    REAL(KIND=r8)    :: rdt
    LOGICAL          :: ghl_local
    REAL(KIND=r8)    :: snowflg(iMax)   
    REAL(KIND=r8)    :: prec_zmc(iMax)                ! total precipitation from ZM convection
    REAL(KIND=r8)    :: snow_zmc(iMax)                ! snow from ZM convection 
    real(KIND=r8)    :: prec_cmf(iMax)                ! total precipitation from Hack convection
    real(KIND=r8)    :: snow_cmf(iMax)                ! snow from Hack convection
    REAL(KIND=r8)    :: SCRa(iMax,kMax)
    REAL(KIND=r8)    :: SCRb(iMax,kMax)
    REAL(KIND=r8)    :: RAINNCV(iMax)
    REAL(KIND=r8)    :: SNOWNCV(iMax)
    REAL(KIND=r8)    :: cape(iMax) 
    REAL(KIND=r8)    :: cine (iMax) 
    REAL(KIND=r8)    :: LCL (iMax) 
    REAL(KIND=r8)    :: LFC (iMax) 
    REAL(KIND=r8)    :: SWEAT(iMax)
    INTEGER          :: i3dflag
!-----------
!nilo  IN-OUT FOR NEW_GRELL 
!-----------
    REAL(KIND=r8)  :: ddmu(iMax,kMax)
    REAL(KIND=r8)  :: ddql(iMax,kMax)

    REAL(KIND=r8) :: cape_old(iMax)


    !-----------------------------------------------------------------
    !-----------------------------------------------------------------

    rdt=1.0_r8/dt
    ghl_local = IsGridHistoryOn()

    ! Check for negative values of specific humidity
    ! Convert virtual temperature to thermodinamic temperature
    CALL qnegat (qa, fdqn, ta, (1.0_r8/dt), del, iMax, kMax)! time t-1
    CALL qnegat (qb, fdqn, tb, (1.0_r8/dt), del, iMax, kMax)! time t
    CALL qnegat (qc, fdqn, tc, (1.0_r8/dt), del, iMax, kMax)! time t+1

    ! Initialize cloud variables with unrealistic values
    DO i=1,iMax
       kbot (i) = 1
       ktop (i) = -1000
       ktops(i) = -1000
       kuo  (i) = -1000
       plcl (i) = -1.0e3_r8
       rliq (i) = 0.0_r8
       rliq2 (i) = 0.0_r8
    END DO   
    hrem= 0.0_r8
    qrem= 0.0_r8
    ddmu= 0.0_r8
    ddql= 0.0_r8
    dq  = 0.0_r8
    cmfmc  = 0.0_r8! convective mass flux--m sub c
    cmfmc2 = 0.0_r8! shallow convective mass flux--m sub c
    prec_cmf= 0.0_r8
    snow_cmf= 0.0_r8
    !ql3(1:iMax,1:kMax,latco)=0.0_r8
    !qi3(1:iMax,1:kMax,latco)=0.0_r8

    ! Initialize surface variables
    DO i=1,iMax
       !surface pressure
       !!T+1                ps2(i)=psb(i)
       !!T                  ps2(i)=psb2(i)
       ps1(i)    =psb2(i)!T+1 
       ps2(i)    =psb(i) !T   
       PS_work(i)=psb(i)

       terr(i)   =MAX(zs(i),0.0_r8)

       RAINCV(i)     = 0.0_r8
       SNOWCV(i)     = 0.0_r8
       Total_Rain(i) = 0.0_r8
       Total_Snow(i) = 0.0_r8
    END DO   

    ! Copy (T,q) at t-1, t and t+1 to work arrays
    DO i=1,iMax
       DO k=1,kMax
          dlf(i,k)=0.0_r8 
          zdu(i,k)=0.0_r8 
          T1(i,k)=ta(i,k)
          T2(i,k)=tb(i,k)
          T3(i,k)=tc(i,k)

          q1(i,k)=qa(i,k)
          q2(i,k)=qb(i,k)
          q3(i,k)=qc(i,k)
       END DO
    END DO
    IF (microphys) THEN
       DO k=1,kMax
          DO i=1,iMax
             ql3     (i,k,latco) = gliqp(i,k)
             qi3     (i,k,latco) = gicep(i,k)

             ql2     (i,k,latco) = gliqm(i,k) 
             qi2     (i,k,latco) = gicem(i,k)
          END DO
       END DO
    ELSE
       DO k=1,kMax
          DO i=1,iMax
             ql3     (i,k,latco) = 0.0_r8 
             qi3     (i,k,latco) = 0.0_r8 

             ql2     (i,k,latco) = 0.0_r8
             qi2     (i,k,latco) = 0.0_r8
          END DO
       END DO
    END IF
    DO i=1,iMax
       DO k=1,kMax
          tBegin(i,k)= tc(i,k)! time t+1
          qBegin(i,k)= qc(i,k)! time t+1
       END DO
    END DO
    !-----------------------------------------------------------------
    ! Calcule CAPE and CIN
    !-----------------------------------------------------------------
!CAPE em T-------------
    i3dflag=0
    SCRa=0.0_r8
    SCRb=0.0_r8
    cape_old=0.0_r8
    CALL calc_cape( &
       iMax                 , &!INTEGER      , INTENT(IN   ) :: nCols
       kMax                 , &!INTEGER      , INTENT(IN   ) :: kMax
       si                   , &!REAL(KIND=r8), INTENT(IN   ) :: si   (kMax+1)
       sl                   , &!REAL(KIND=r8), INTENT(IN   ) :: sl   (kMax)
       ps2*10.0_r8          , &!REAL(KIND=r8), INTENT(IN   ) :: PSFC (nCols)   ! psfc is pressure in hPa
       terr                 , &!REAL(KIND=r8), INTENT(IN   ) :: HGT  (nCols)   ! topography m
       t2                   , &!REAL(KIND=r8), INTENT(IN   ) :: TK    (nCols,kMax)     ! TK is temp in K, T is theta-300
       q2                   , &!REAL(KIND=r8), INTENT(IN   ) :: QV    (nCols,kMax)
       SCRa                 , &!REAL(KIND=r8), INTENT(OUT  ) :: SCRa  (nCols,kMax)
       SCRb                 , &!REAL(KIND=r8), INTENT(OUT  ) :: SCRb  (nCols,kMax)
       i3dflag                )!INTEGER      , INTENT(IN   ) :: i3dflag
       DO i=1,iMax
          cape_old(i) = MAX(SCRa(i,1),0.0_r8)
       END DO
!CAPE em T+1-------------
    i3dflag=0
    SCRa=0.0_r8
    SCRb=0.0_r8
    cape=0.0_r8
    cine=0.0_r8
    CALL calc_cape( &
       iMax                 , &!INTEGER      , INTENT(IN   ) :: nCols
       kMax                 , &!INTEGER      , INTENT(IN   ) :: kMax
       si                   , &!REAL(KIND=r8), INTENT(IN   ) :: si   (kMax+1)
       sl                   , &!REAL(KIND=r8), INTENT(IN   ) :: sl   (kMax)
       ps2*10.0_r8          , &!REAL(KIND=r8), INTENT(IN   ) :: PSFC (nCols)   ! psfc is pressure in hPa
       terr                 , &!REAL(KIND=r8), INTENT(IN   ) :: HGT  (nCols)   ! topography m
       t3                   , &!REAL(KIND=r8), INTENT(IN   ) :: TK    (nCols,kMax)     ! TK is temp in K, T is theta-300
       q3                   , &!REAL(KIND=r8), INTENT(IN   ) :: QV    (nCols,kMax)
       SCRa                 , &!REAL(KIND=r8), INTENT(OUT  ) :: SCRa  (nCols,kMax)
       SCRb                 , &!REAL(KIND=r8), INTENT(OUT  ) :: SCRb  (nCols,kMax)
       i3dflag                )!INTEGER      , INTENT(IN   ) :: i3dflag
       DO i=1,iMax
          cape(i) = MAX(SCRa(i,1),0.0_r8)
          cine(i) = MAX(SCRa(i,2),0.0_r8)
          LCL (i) = SCRa(i,3)
          LFC (i) = SCRa(i,4)
       END DO
       SWEAT=0.0_r8
       SWEAT=SWEAT_index(t3,ub,vb,iMax,kMax)

    !-----------------------------------------------------------------
    ! Deep Convection
    !-----------------------------------------------------------------
       IF((nClass+nAeros)>0 .and. PRESENT(gvarm))THEN

       CALL RunDeepConvection(&
      ! Run Flags
            iccon                          ,& !CHARACTER(LEN=*),    INTENT(IN   ) :: iccon
            cflric                         ,& !REAL(KIND=r8)   ,    INTENT(in   ) :: cflric
            microphys   , & !LOGICAL      , INTENT(in   ) :: microphys
      ! Time info
            dt                             ,& !REAL(KIND=r8)   ,    INTENT(in   ) :: dt
      ! Model Geometry
            sl        (1:kmax)             ,& !REAL(KIND=r8)   ,    INTENT(in   ) :: sl        (kmax)
            si        (1:kmax+1)           ,& !REAL(KIND=r8)   ,    INTENT(in   ) :: si        (kMax+1)
            del       (1:kmax)             ,& !REAL(KIND=r8)   ,    INTENT(in   ) :: del       (kmax)
      ! Model information
            jdt                             ,& !INTEGER         ,    INTENT(in   ) :: iMax  
            nClass                         ,& !INTEGER         ,    INTENT(in   ) :: iMax  
            nAeros                         ,&
            iMax                           ,& !INTEGER         ,    INTENT(in   ) :: iMax  
            kMax                           ,& !INTEGER         ,    INTENT(in   ) :: kMax
            nls                            ,& !INTEGER         ,    INTENT(in   ) :: nls 
            latco                          ,& !INTEGER         ,    INTENT(in   ) :: latco  
            mask      (1:iMax)             ,& !INTEGER(KIND=i8),    INTENT(IN   ) :: mask      (iMax)
            terr      (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: terr      (iMax)  
            colrad    (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: colrad    (iMax)
            lonrad    (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: lonrad    (iMax)  
      ! Surface field
            sens      (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: sens      (1:nCols)
            tpert     (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: tpert     (1:iMax)
            ustar     (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: ustar     (1:iMax)
      ! PBL field
            pblh      (1:iMax,latco)       ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: pblh      (1:iMax)
            tke       (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: tke       (iMax,kMax)
      ! CONVECTION: Cloud field
            ktop      (1:iMax)             ,& !INTEGER         ,    INTENT(inout) :: ktop      (iMax)
            kbot      (1:iMax)             ,& !INTEGER         ,    INTENT(inout) :: kbot      (iMax)
            ktops     (1:iMax)             ,& !INTEGER         ,    INTENT(INOUT) :: ktops     (iMax)
            kbots     (1:iMax)             ,& !INTEGER         ,    INTENT(INOUT) :: kbots     (iMax)
            cape_old  (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: cape_old  (1:nCols)
            cape      (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: cape      (1:nCols)
            cine      (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: cine      (1:nCols)
            hrem      (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(inout) :: hrem      (iMax,kmax)
            qrem      (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(inout) :: qrem      (iMax,kmax)
            kuo       (1:iMax)             ,& !INTEGER         ,    INTENT(inout) :: kuo       (iMax)
            cldm      (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(inout) :: cldm      (iMax)
            plcl      (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(INOUT) :: plcl      (iMax)
            ddql      (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(INOUT) :: ddql      (1:nCols,1:kMax)
            ddmu      (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(INOUT) :: ddmu      (1:nCols,1:kMax)
            fdqn      (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(INOUT) :: fdqn      (1:iMax,1:kMax)
            qliq      (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(OUT  ) :: qliq      (iMax,kMax) !qrc liquid water content in
            cld       (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(OUT  ) :: cld       (1:iMax,1:kMax)
            cmfmc     (1:iMax,1:kMax+1)    ,& !REAL(KIND=r8)   ,    INTENT(OUT  ) :: cmfmc     (1:iMax,1:kMax+1)
            dlf       (1:iMax,1:kMax)      ,& !REAL(KIND=r8)   ,    INTENT(OUT  ) :: dlf       (1:iMax,1:kMax)
            zdu       (1:iMax,1:kMax)      ,& !REAL(KIND=r8)   ,    INTENT(OUT  ) :: zdu       (1:iMax,1:kMax)
            rliq      (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(OUT  ) :: rliq      (1:iMax)
            RAINCV    (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(inout) :: RAINCV    (iMax)
            SNOWCV    (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(inout) :: SNOWCV    (iMax)
            snow_zmc  (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(OUT  ) :: snow_zmc  (1:iMax)
            prec_zmc  (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(OUT  ) :: prec_zmc  (1:iMax)
            Total_Rain(1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(INOUT) :: Total_Rain(iMax)
            Total_Snow(1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(INOUT) :: Total_Snow(iMax)
            dudt      (1:iMax,1:kMax)      ,& !
            dvdt      (1:iMax,1:kMax)      ,& !
      ! Atmospheric fields
            tc        (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(in   ) :: tc        (iMax,kmax)
            qc        (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(in   ) :: qc        (iMax,kmax)
            ps2       (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(in   ) :: ps2       (iMax)
            ub        (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: ub        (iMax,kMax) ! (m/s) 
            vb        (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: vb        (iMax,kMax) ! (m/s)
            omgb      (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: omgb      (iMax,kMax) ! (Pa/s)
            T2        (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(inout) :: T2        (iMax,kmax)
            T3        (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(inout) :: T3        (iMax,kmax)
            Q1        (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(inout) :: Q1        (iMax,kmax)
            Q2        (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(inout) :: Q2        (iMax,kmax)
            Q3        (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(inout) :: Q3        (iMax,kmax)
            dq        (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(inout) :: dq        (iMax,kmax)
            dump      (1:iMax,1:kmax,latco),& !REAL(KIND=r8)   , INTENT(inout) :: ql2     (1:iMax,1:kmax)
       ! Microphysics
            ql2       (1:iMax,1:kmax,latco),& !REAL(KIND=r8)   , INTENT(inout) :: ql2     (1:iMax,1:kmax)
            qi2       (1:iMax,1:kmax,latco),& !REAL(KIND=r8)   , INTENT(inout) :: qi2     (1:iMax,1:kmax)
            ql3       (1:iMax,1:kmax,latco),& !REAL(KIND=r8)   ,    INTENT(inout) :: ql3       (iMax,kmax)
            qi3       (1:iMax,1:kmax,latco),& !REAL(KIND=r8)   ,    INTENT(inout) :: qi3       (iMax,kmax)
            gvarm     (1:iMax,1:kmax,1:nClass+nAeros), &
            gvarp     (1:iMax,1:kmax,1:nClass+nAeros))  
    ! Save humd/temp after deep convection
    IF (microphys) THEN
       DO n=1,nClass+nAeros
          DO k=1,kMax
             DO i=1,iMax
                gvarp     (i,k,n)= MAX(gvarp     (i,k,n),0.0e-22_r8)
                gvarm     (i,k,n)= MAX(gvarm     (i,k,n),0.0e-22_r8)
             END DO
          END DO
      END DO
   END IF 

    ELSE
       CALL RunDeepConvection(&
      ! Run Flags
            iccon                          ,& !CHARACTER(LEN=*),    INTENT(IN   ) :: iccon
            cflric                         ,& !REAL(KIND=r8)   ,    INTENT(in   ) :: cflric
            microphys   , & !LOGICAL      , INTENT(in   ) :: microphys
      ! Time info
            dt                             ,& !REAL(KIND=r8)   ,    INTENT(in   ) :: dt
      ! Model Geometry
            sl        (1:kmax)             ,& !REAL(KIND=r8)   ,    INTENT(in   ) :: sl        (kmax)
            si        (1:kmax+1)           ,& !REAL(KIND=r8)   ,    INTENT(in   ) :: si        (kMax+1)
            del       (1:kmax)             ,& !REAL(KIND=r8)   ,    INTENT(in   ) :: del       (kmax)
      ! Model information
            jdt                             ,& !INTEGER         ,    INTENT(in   ) :: iMax  
            nClass                         ,& !INTEGER         ,    INTENT(in   ) :: iMax  
            nAeros                         ,&
            iMax                           ,& !INTEGER         ,    INTENT(in   ) :: iMax  
            kMax                           ,& !INTEGER         ,    INTENT(in   ) :: kMax
            nls                            ,& !INTEGER         ,    INTENT(in   ) :: nls 
            latco                          ,& !INTEGER         ,    INTENT(in   ) :: latco  
            mask      (1:iMax)             ,& !INTEGER(KIND=i8),    INTENT(IN   ) :: mask      (iMax)
            terr      (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: terr      (iMax)  
            colrad    (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: colrad    (iMax)
            lonrad    (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: lonrad    (iMax)  
      ! Surface field
            sens      (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: sens      (1:nCols)
            tpert     (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: tpert     (1:iMax)
            ustar     (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: ustar     (1:iMax)
      ! PBL field
            pblh      (1:iMax,latco)       ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: pblh      (1:iMax)
            tke       (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: tke       (iMax,kMax)
      ! CONVECTION: Cloud field
            ktop      (1:iMax)             ,& !INTEGER         ,    INTENT(inout) :: ktop      (iMax)
            kbot      (1:iMax)             ,& !INTEGER         ,    INTENT(inout) :: kbot      (iMax)
            ktops     (1:iMax)             ,& !INTEGER         ,    INTENT(INOUT) :: ktops     (iMax)
            kbots     (1:iMax)             ,& !INTEGER         ,    INTENT(INOUT) :: kbots     (iMax)
            cape_old  (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: cape_old  (1:nCols)
            cape      (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: cape      (1:nCols)
            cine      (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: cine      (1:nCols)
            hrem      (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(inout) :: hrem      (iMax,kmax)
            qrem      (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(inout) :: qrem      (iMax,kmax)
            kuo       (1:iMax)             ,& !INTEGER         ,    INTENT(inout) :: kuo       (iMax)
            cldm      (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(inout) :: cldm      (iMax)
            plcl      (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(INOUT) :: plcl      (iMax)
            ddql      (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(INOUT) :: ddql      (1:nCols,1:kMax)
            ddmu      (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(INOUT) :: ddmu      (1:nCols,1:kMax)
            fdqn      (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(INOUT) :: fdqn      (1:iMax,1:kMax)
            qliq      (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(OUT  ) :: qliq      (iMax,kMax) !qrc liquid water content in
            cld       (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(OUT  ) :: cld       (1:iMax,1:kMax)
            cmfmc     (1:iMax,1:kMax+1)    ,& !REAL(KIND=r8)   ,    INTENT(OUT  ) :: cmfmc     (1:iMax,1:kMax+1)
            dlf       (1:iMax,1:kMax)      ,& !REAL(KIND=r8)   ,    INTENT(OUT  ) :: dlf       (1:iMax,1:kMax)
            zdu       (1:iMax,1:kMax)      ,& !REAL(KIND=r8)   ,    INTENT(OUT  ) :: zdu       (1:iMax,1:kMax)
            rliq      (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(OUT  ) :: rliq      (1:iMax)
            RAINCV    (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(inout) :: RAINCV    (iMax)
            SNOWCV    (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(inout) :: SNOWCV    (iMax)
            snow_zmc  (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(OUT  ) :: snow_zmc  (1:iMax)
            prec_zmc  (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(OUT  ) :: prec_zmc  (1:iMax)
            Total_Rain(1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(INOUT) :: Total_Rain(iMax)
            Total_Snow(1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(INOUT) :: Total_Snow(iMax)
            dudt      (1:iMax,1:kMax)      ,& !
            dvdt      (1:iMax,1:kMax)      ,& !
      ! Atmospheric fields
            tc        (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(in   ) :: tc        (iMax,kmax)
            qc        (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(in   ) :: qc        (iMax,kmax)
            ps2       (1:iMax)             ,& !REAL(KIND=r8)   ,    INTENT(in   ) :: ps2       (iMax)
            ub        (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: ub        (iMax,kMax) ! (m/s) 
            vb        (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: vb        (iMax,kMax) ! (m/s)
            omgb      (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: omgb      (iMax,kMax) ! (Pa/s)
            T2        (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(inout) :: T2        (iMax,kmax)
            T3        (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(inout) :: T3        (iMax,kmax)
            Q1        (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(inout) :: Q1        (iMax,kmax)
            Q2        (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(inout) :: Q2        (iMax,kmax)
            Q3        (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(inout) :: Q3        (iMax,kmax)
            dq        (1:iMax,1:kmax)      ,& !REAL(KIND=r8)   ,    INTENT(inout) :: dq        (iMax,kmax)
            dump      (1:iMax,1:kmax,latco),& !REAL(KIND=r8)   , INTENT(inout) :: ql2     (1:iMax,1:kmax)
       ! Microphysics
            ql2       (1:iMax,1:kmax,latco), & !REAL(KIND=r8)   , INTENT(inout) :: ql2     (1:iMax,1:kmax)
            qi2       (1:iMax,1:kmax,latco),& !REAL(KIND=r8)   , INTENT(inout) :: qi2     (1:iMax,1:kmax)
            ql3       (1:iMax,1:kmax,latco),& !REAL(KIND=r8)   ,    INTENT(inout) :: ql3       (iMax,kmax)
            qi3       (1:iMax,1:kmax,latco) ) !REAL(KIND=r8)   ,    INTENT(inout) :: qi3       (iMax,kmax)

 
    END IF

    DO i=1,iMax
       cu_ktop(i,latco)= ktop(i)
       cu_kbot(i,latco)= kbot(i)
       cu_Kuo (i,latco)= kuo (i)
    END DO

    ! Save humd/temp after deep convection
    DO k=1,kMax
       DO i=1,iMax
          q3   (i,k)       = MAX(q3  (i,k)      ,1.0e-12_r8)
          IF (microphys) THEN
          ql3  (i,k,latco) = MAX(ql3 (i,k,latco),0.0e-12_r8)
          qi3  (i,k,latco) = MAX(qi3 (i,k,latco),0.0e-12_r8)
          ql2  (i,k,latco) = MAX(ql2 (i,k,latco),0.0e-12_r8)
          qi2  (i,k,latco) = MAX(qi2 (i,k,latco),0.0e-12_r8)
          END IF 
          qDeep(i,k)=q3(i,k)
          tDeep(i,k)=t3(i,k)
          cu_hr(i,k,latco) = fac*rdt*(tDeep(i,k)-tBegin(i,k))
       END DO
    END DO


    !-----------------------------------------------------------------
    ! Shallow Convection
    !-----------------------------------------------------------------

    CALL RunShallowConvection( &
      ! Run Flags
                      ISCON                           , & !CHARACTER(LEN=*), INTENT(IN   ) :: ISCON
                      iccon                           , & !CHARACTER(LEN=*), INTENT(IN   ) :: iccon
      ! Time info
                      jdt                             , & !INTEGER         , INTENT(in   ) :: jdt
                      dt                              , & !REAL(KIND=r8)   , INTENT(in   ) :: dt
      ! Model information
                      iMax                            , & !INTEGER         , INTENT(in   ) :: iMax
                      kmax                            , & !INTEGER         , INTENT(in   ) :: kmax
                      latco                           , & !INTEGER         , INTENT(in   ) :: latco
                      terr     (1:iMax)               , & !REAL(KIND=r8)   , INTENT(in   ) :: terr    (1:iMax)
      ! Model Geometry
                      si       (1:kMax+1)             , & !REAL(KIND=r8)   , INTENT(in   ) :: si      (1:kMax+1)
                      sl       (1:kmax)               , & !REAL(KIND=r8)   , INTENT(in   ) :: sl      (1:kmax)
                      del      (1:kmax)               , & !REAL(KIND=r8)   , INTENT(in   ) :: del     (1:kmax)
      ! Surface field
                      mask     (1:iMax)               , & !INTEGER(KIND=i8),INTENT(IN ) :: mask (iMax) 
                      sens     (1:iMax)               , & !
                      evap     (1:iMax)               , & !
                      qpert    (1:iMax)               , & !REAL(KIND=r8)   , INTENT(in   ) :: qpert   (1:iMax)
      ! PBL field
                      pblh     (1:iMax,latco)         , & !REAL(KIND=r8)   , INTENT(in   ) :: pblh    (1:iMax)
                      tke      (1:iMax,1:kMax)        , & !REAL(KIND=r8)   , INTENT(in   ) :: tke     (1:iMax,1:kMax)
      ! CONVECTION: Cloud field
                      dudt      (1:iMax,1:kMax)       ,& !
                      dvdt      (1:iMax,1:kMax)       ,& !
                      rliq     (1:iMax)               , & !REAL(KIND=r8)   , INTENT(in   ) :: rliq    (1:iMax)
                      ktop     (1:iMax)               , & !INTEGER         , INTENT(in   ) :: ktop    (1:iMax)
                      ktops    (1:iMax)               , & !INTEGER         , INTENT(in   ) :: ktops   (1:iMax)
                      kuo      (1:iMax)               , & !INTEGER         , INTENT(in   ) :: kuo     (1:iMax)
                      plcl     (1:iMax)               , & !REAL(KIND=r8)   , INTENT(inout) :: plcl    (1:iMax)
                      kcbot1   (1:iMax)               , & !INTEGER         , INTENT(inout) :: kcbot1  (1:iMax)
                      kctop1   (1:iMax)               , & !INTEGER         , INTENT(inout) :: kctop1  (1:iMax)
                      noshal   (1:iMax)               , & !INTEGER         , INTENT(inout) :: noshal  (1:iMax)
                      concld   (1:iMax,1:kMax)        , & !REAL(KIND=r8)   , INTENT(inout) :: concld  (1:iMax,1:kMax)
                      cld      (1:iMax,1:kMax)        , & !REAL(KIND=r8)   , INTENT(inout) :: cld     (1:iMax,1:kMax)
                      cmfmc    (1:iMax,1:kMax+1)      , & !REAL(KIND=r8)   , INTENT(inout) :: cmfmc   (1:iMax,1:kMax+1)
                      cmfmc2   (1:iMax,1:kMax+1)      , & !REAL(KIND=r8)   , INTENT(out  ) :: cmfmc2  (1:iMax,1:kMax+1)
                      dlf      (1:iMax,1:kMax)        , & !REAL(KIND=r8)   , INTENT(inout) :: dlf     (1:iMax,1:kMax)
                      fdqn     (1:iMax,1:kMax)        , & !REAL(KIND=r8)   , INTENT(inout) :: fdqn     (1:iMax,1:kMax)
                      rliq2    (1:iMax)               , & !REAL(KIND=r8)   , INTENT(out  ) :: rliq2   (1:iMax)
                      snow_cmf (1:iMax)               , & !REAL(KIND=r8)   , INTENT(out  ) :: snow_cmf(1:iMax)
                      prec_cmf (1:iMax)               , & !REAL(KIND=r8)   , INTENT(out  ) :: prec_cmf(1:iMax)
                      RAINCV   (1:iMax)               , & !REAL(KIND=r8)   , INTENT(out  ) :: RAINCV  (1:iMax)
                      SNOWCV   (1:iMax)               , & !REAL(KIND=r8)   , INTENT(out  ) :: SNOWCV  (1:iMax)
      ! Atmospheric fields
                      PS_work  (1:iMax)               , & !REAL(KIND=r8)   , INTENT(in   ) :: PS_work (1:iMax)
                      ps2      (1:iMax)               , & !REAL(KIND=r8)   , INTENT(in   ) :: ps2     (1:iMax)
                      ub       (1:iMax,1:kmax)        , & !REAL(KIND=r8)   , INTENT(in   ) :: ub      (1:iMax,1:kmax)
                      vb       (1:iMax,1:kmax)        , & !REAL(KIND=r8)   , INTENT(in   ) :: vb      (1:iMax,1:kmax)
                      omgb     (1:iMax,1:kmax)        , & !REAL(KIND=r8)   , INTENT(in   ) :: omgb    (1:iMax,1:kmax)
                      t3       (1:iMax,1:kmax)        , & !REAL(KIND=r8)   , INTENT(inout) :: t3      (1:iMax,1:kmax)
                      q3       (1:iMax,1:kmax)        , & !REAL(KIND=r8)   , INTENT(inout) :: q3      (1:iMax,1:kmax)
                      ql3      (1:iMax,1:kmax,latco)  , & !REAL(KIND=r8)   , INTENT(inout) :: ql3     (1:iMax,1:kmax)
                      qi3      (1:iMax,1:kmax,latco)    ) !REAL(KIND=r8)   , INTENT(inout) :: qi3     (1:iMax,1:kmax)

    !Total_Rain=RAINCV
    !Total_Snow=SNOWCV
    ! Save humd/temp after shallow convection
    DO k=1,kMax
      DO i=1,iMax
        q3  (i,k)       = MAX(q3  (i,k)      ,1.0e-12_r8)
        IF (microphys) THEN
        ql3 (i,k,latco) = MAX(ql3 (i,k,latco),1.0e-12_r8)
        qi3 (i,k,latco) = MAX(qi3 (i,k,latco),1.0e-12_r8)
        ql2 (i,k,latco) = MAX(ql2 (i,k,latco),1.0e-12_r8)
        qi2 (i,k,latco) = MAX(qi2 (i,k,latco),1.0e-12_r8)
         END IF
        tShal(i,k)=t3(i,k)
        qShal(i,k)=q3(i,k)
      END DO
    END DO



    IF (microphys) THEN
       DO k=1,kMax
          DO i=1,iMax
              gliqp(i,k)   =ql3     (i,k,latco) 
              gicep(i,k)   =qi3     (i,k,latco) 

              gliqm(i,k)   =ql2     (i,k,latco) 
              gicem(i,k)   =qi2     (i,k,latco) 
          END DO
       END DO   
    END IF
    !-----------------------------------------------------------------
    ! MicroPhysics Precipitation
    !-----------------------------------------------------------------
      IF((nClass+nAeros)>0 .and. PRESENT(gvarm))THEN
    
    CALL RunMicroPhysics(&
      ! Run Flags
               mlrg                                , &
               ILCON                               , &
               microphys                           , &
      ! Time info
               jdt                                 , &
               dt                                  , &
      ! Model Geometry
               sl          (1:kMax)                , &
               si          (1:kMax+1)              , &
               del         (1:kMax)                , &
               colrad      (1:iMax)                , &
               LOWLYR      (1:iMax,latco)          , &
               terr        (1:iMax)                , &
      ! Model information
               mask        (1:iMax)                , &
               nClass                              , &
               nAeros                              , &
               iMax                                , &
               kMax                                , &
               latco                               , &
      ! Surface field
               tsfc        (1:iMax)                , &
      ! PBL field
               PBL_CoefKh  (1:iMax,1:kMax)         , &
               tke         (1:iMax,1:kMax)         , &
               pblh        (1:iMax,latco)          , &
      ! CONVECTION: Cloud field
               kuo         (1:iMax)                , &
               cmfmc       (1:iMax,1:kMax+1)       , &
               cmfmc2      (1:iMax,1:kMax+1)       , &
               dlf         (1:iMax,1:kMax)         , &
               rliq        (1:iMax)                , &
               concld      (1:iMax,1:kMax)         , &
               cld         (1:iMax,1:kMax)         , &
               fdqn        (1:iMax,1:kMax)         , &
               EFFCS       (1:iMax,1:kMax)         , &
               EFFIS       (1:iMax,1:kMax)         , &
               Total_Rain  (1:iMax)                , &
               Total_Snow  (1:iMax)                , &
               RAINNCV     (1:iMax)                , &
               SNOWNCV     (1:iMax)                , &
               F_ICE_PHY   (1:iMax,1:kMax,latco)   , &
               F_RAIN_PHY  (1:iMax,1:kMax,latco)   , &
               F_RIMEF_PHY (1:iMax,1:kMax,latco)   , &
      ! Atmospheric fields
               ps2         (1:iMax)                , &
               t2          (1:iMax,1:kMax)         , &
               t3          (1:iMax,1:kMax)         , &
               q2          (1:iMax,1:kMax)         , &
               q3          (1:iMax,1:kMax)         , &
               ub          (1:iMax,1:kMax)         , &
               vb          (1:iMax,1:kMax)         , &
               omgb        (1:iMax,1:kMax)         , &
               dq          (1:iMax,1:kMax)         , &
               tLrgs       (1:iMax,1:kMax)         , & 
               qLrgs       (1:iMax,1:kMax)         , & 
      ! Microphysics
               gicem       (1:iMax,1:kmax)         , &
               gicep       (1:iMax,1:kmax)         , &
               gliqm       (1:iMax,1:kmax)         , &
               gliqp       (1:iMax,1:kmax)         , &
               gvarm       (1:iMax,1:kmax,1:nClass+nAeros), &
               gvarp       (1:iMax,1:kmax,1:nClass+nAeros)  )
    ELSE
    CALL RunMicroPhysics(&
      ! Run Flags
               mlrg                                 , &
               ILCON                                , &
               microphys                            , &
      ! Time info
               jdt                                  , &
               dt                                   , &
      ! Model Geometry
               sl           (1:kMax)                , &
               si           (1:kMax+1)              , &
               del          (1:kMax)                , &
               colrad       (1:iMax)                , &
               LOWLYR       (1:iMax,latco)          , &
               terr         (1:iMax)                , &
      ! Model information
               mask         (1:iMax)                , &
               nClass                               , &
               nAeros                               , &
               iMax                                 , &
               kMax                                 , &
               latco                                , &
      ! Surface field
               tsfc         (1:iMax)                , &
      ! PBL field
               PBL_CoefKh   (1:iMax,1:kMax)         , &
               tke          (1:iMax,1:kMax)         , &
               pblh         (1:iMax,latco)          , &
      ! CONVECTION: Cloud field
               kuo          (1:iMax)                , &
               cmfmc        (1:iMax,1:kMax+1)       , &
               cmfmc2       (1:iMax,1:kMax+1)       , &
               dlf          (1:iMax,1:kMax)         , &
               rliq         (1:iMax)                , &
               concld       (1:iMax,1:kMax)         , &
               cld          (1:iMax,1:kMax)         , &
               fdqn         (1:iMax,1:kMax)         , &
               EFFCS        (1:iMax,1:kMax)         , &
               EFFIS        (1:iMax,1:kMax)         , &
               Total_Rain   (1:iMax)                , &
               Total_Snow   (1:iMax)                , &
               RAINNCV      (1:iMax)                , &
               SNOWNCV      (1:iMax)                , &
               F_ICE_PHY    (1:iMax,1:kMax,latco)   , &
               F_RAIN_PHY   (1:iMax,1:kMax,latco)   , &
               F_RIMEF_PHY  (1:iMax,1:kMax,latco)   , &
      ! Atmospheric fields
               ps2          (1:iMax)                , &
               t2           (1:iMax,1:kMax)         , &
               t3           (1:iMax,1:kMax)         , &
               q2           (1:iMax,1:kMax)         , &
               q3           (1:iMax,1:kMax)         , &
               ub           (1:iMax,1:kMax)         , &
               vb           (1:iMax,1:kMax)         , &
               omgb         (1:iMax,1:kMax)         , &
               dq           (1:iMax,1:kMax)         , &
               tLrgs        (1:iMax,1:kMax)         , & 
               qLrgs        (1:iMax,1:kMax)         , & 
      ! Microphysics
               gicem        (1:iMax,1:kmax)         , &
               gicep        (1:iMax,1:kmax)         , &
               gliqm        (1:iMax,1:kmax)         , &
               gliqp        (1:iMax,1:kmax)           )
 
    END IF  
    !-----------------------------------------------------------------
    ! Convective Cloud Cover
    !-----------------------------------------------------------------
    
    IF(TRIM(iccon).EQ.'ARA'.OR.TRIM(iccon).EQ.'KUO'.OR.TRIM(iccon).EQ.'GRE'.OR.&
       TRIM(iccon).EQ.'RAS'.OR. TRIM(iccon).EQ.'ZMC'.OR.TRIM(iccon).EQ.'GEC'&
       .OR.TRIM(iccon).EQ.'OCL')THEN
      CALL CLOUD_COVER( &
             kt     ,ktp    ,iMax   ,kbot   ,ktop   ,noshal ,kctop1 , &
             kcbot1 ,RAINCV ,fac2x  ,rccmbl ,iccon  ,convc  ,convt  , &
             convb  ,prcp1  ,prcp2  ,prcp3  ,prcpt  ,toplv  ,botlv  , &
             convts ,convcs ,convbs )
    END IF

    !-----------------------------------------------------------------
    ! DIAGNOSTICS AND MODEL OUTPUT
    !-----------------------------------------------------------------

    !---------------------
    ! Calculate deep convection moistening and heating profiles
    !print*, 'Emarg: kmax=, max(tBegin),max(tDeep)=',kmax, max(tBegin(:,1:kmax/2),tBegin(:,kmax/2+1:kmax)),max(tDeep(:,1:kmax/2),tDeep(:,kmax/2+1:kmax))
    DO k=1,kMax
       DO i=1,iMax
          clheat(i,k)=fac*rdt*(tDeep(i,k)-tBegin(i,k))
          cmchan(i,k)=fac*rdt*(qDeep(i,k)-qBegin(i,k))
       END DO
    END DO

    !---------------------
    ! Calculate shallow convection moistening and heating profiles
    !print*, 'Emarg: kmax=,max(tShal),max(tDeep)=',kmax, max(tShal(:,1:kmax/2),tShal(:,kmax/2+1:kmax)),max(tDeep(:,1:kmax/2),tDeep(:,kmax/2+1:kmax))
    DO k=1,kMax
       DO i=1,iMax
          sclhea(i,k)=fac*rdt*(tShal(i,k)-tDeep(i,k))
          scmcha(i,k)=fac*rdt*(qShal(i,k)-qDeep(i,k))
       END DO
    END DO

    !---------------------
    ! Calculate large scale convection moistening and heating profiles
    DO k=1,kMax
       DO i=1,iMax
          lslhea(i,k)=fac*rdt*(tLrgs(i,k)-tShal(i,k))
          lsmcha(i,k)=fac*rdt*(qLrgs(i,k)-qShal(i,k))
       END DO
    END DO

    !***********************************
    ! move qDeep to qb and tDeep to tb
    ! move T3 to Tc and q3 to qc
    !***********************************
    DO i=1,iMax
       DO k=1,kMax
          ! Update T
          qb(i,k)=qDeep(i,k)!q  after deep convection
          tb(i,k)=tDeep(i,k)*(1.0_r8+delq*qDeep(i,k))!t  after deep convection
          ! Update T+1
          qc(i,k) = qLrgs(i,k)
          tc(i,k) = tLrgs(i,k)*(1.0_r8+delq*qLrgs(i,k))
       END DO
    END DO
    IF(TRIM(ISCON).EQ.'JHK' .OR. TRIM(ISCON).EQ.'UW' )THEN
       IF(TRIM(ILCON).EQ.'MIC'.or. TRIM(ILCON).EQ.'HWRF' .or. TRIM(ILCON).EQ.'HGFS'.or.&
          TRIM(ILCON).EQ.'UKMO' .or.TRIM(ILCON).EQ.'MORR' .or.TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HUMN' ) THEN
          ! Calculate precipiation in mm/s
          DO i=1,iMax
             ppci(i)=2.0e0_r8*1.0e3_r8*RAINCV(i) ! deep
             ppli(i)=2.0e0_r8*1.0e3_r8*(Total_Rain(i)+Total_Snow(i))-ppci(i) ! large
             !ppli(i)=2.0e0_r8*1.0e3_r8*(Total_Rain(i))-ppci(i) ! large
          END DO
       ELSE
          ! Calculate precipiation in mm/s
          DO i=1,iMax
             ppci(i)=2.0e0_r8*1.0e3_r8*RAINCV(i) ! deep
             ppli(i)=2.0e0_r8*1.0e3_r8*(Total_Rain(i)+Total_Snow(i))-ppci(i) ! large
             !ppli(i)=2.0e0_r8*1.0e3_r8*(Total_Rain(i))-ppci(i) ! large
          END DO
       END IF
    ELSE    
       IF(TRIM(ILCON).EQ.'MIC'.or. TRIM(ILCON).EQ.'HWRF'.or. TRIM(ILCON).EQ.'HGFS'.or.&
          TRIM(ILCON).EQ.'UKMO'.or. TRIM(ILCON).EQ.'MORR' .or. TRIM(ILCON).EQ.'HUMO' .or.TRIM(ILCON).EQ.'HUMN') THEN
       ! Calculate precipiation in mm/s
          DO i=1,iMax
             ppci(i)=2.0e0_r8*1.0e3_r8*RAINCV(i) ! deep
             ppli(i)=2.0e0_r8*1.0e3_r8*(Total_Rain(i)+Total_Snow(i))-ppci(i) ! large
             !ppli(i)=2.0e0_r8*1.0e3_r8*(Total_Rain(i))-ppci(i) ! large
          END DO
       ELSE
          ! Calculate precipiation in mm/s
          DO i=1,iMax
             ppci(i)=2.0e0_r8*1.0e3_r8*RAINCV(i) ! deep
             ppli(i)=2.0e0_r8*1.0e3_r8*(Total_Rain(i))-ppci(i) ! large
          END DO    
       END IF
    END IF

    ! Calculate precipiation in mm

    DO i=1,iMax
       geshem(i)=geshem(i)+fac2x*Total_Rain(i)
    END DO

    ! Time-step output of precipitation
    IF (doprec) THEN
       DO i=1,iMax
          prcc(i)=fac2*rdt*1.0e3_r8*RAINCV(i)
          prct(i)=fac2*rdt*1.0e3_r8*Total_Rain(i)
       END DO
    END IF
    
    ! Diagnose snow field
    IF(TRIM(ILCON).EQ.'MIC'.or. TRIM(ILCON).EQ.'HWRF'.or. TRIM(ILCON).EQ.'HGFS'.or. &
       TRIM(ILCON).EQ.'UKMO'.or.TRIM(ILCON).EQ.'MORR' .or.TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HUMN') THEN
       DO i=1,iMax          
             snowflg(i) = Total_Snow(i) !m 
             snowfl (i) = fac2*rdt*1.0e3_r8*Total_Snow(i) !mm/s
       END DO
    ELSE
       DO i=1,iMax          
          IF(0.35_r8*tLrgs(i,1)+0.65_r8*tLrgs(i,2).LE.273.2_r8)THEN
             snowflg(i) = Total_Rain(i)
             snowfl (i) = fac2*rdt*1.0e3_r8*Total_Rain(i)
          ELSE
             snowflg(i) = 0.0_r8
             snowfl (i) = 0.0_r8
          END IF
       END DO
    END IF

    !-----------------
    ! Storage Diagnostic Fields
    !------------------
    IF (StartStorDiag)THEN
      IF((nClass+nAeros)>0 .and. PRESENT(gvarp))THEN

       CALL ConvecDiagnStorage (&
            iMax     ,kMax     ,latco     ,rdt      ,fac2     , &
            fdqn     ,RAINCV   ,Total_Rain,snowfl   ,sclhea   , &
            scmcha   ,clheat   ,cmchan    ,lslhea   ,lsmcha   , &
            ql3 (1:iMax,1:kMax,latco),qi3 (1:iMax,1:kMax,latco),&
            cape,cine,LCL,LFC,SWEAT ,&
            gvarp  )
      ELSE
       CALL ConvecDiagnStorage (&
            iMax     ,kMax     ,latco     ,rdt      ,fac2     , &
            fdqn     ,RAINCV   ,Total_Rain,snowfl   ,sclhea   , &
            scmcha   ,clheat   ,cmchan    ,lslhea   ,lsmcha   , &
            ql3 (1:iMax,1:kMax,latco),qi3 (1:iMax,1:kMax,latco),&
            cape,cine,LCL,LFC,SWEAT   )
      END IF
    END IF
    !-----------------
    ! Storage GridHistory Fields
    !------------------
    IF(ghl_local)THEN 
       CALL ConvecGridHistStorage(&
            iMax     ,kMax      ,latco    ,rdt      ,fac2     , &
            RAINCV   ,Total_Rain,snowflg  ,sclhea   ,scmcha   , &
            clheat   ,cmchan    ,lslhea   ,lsmcha)
    END IF

  END SUBROUTINE cumulus_driver

  !-----------------------------------------------------------------
  ! FINALIZE CONVECTION
  !-----------------------------------------------------------------
  
  SUBROUTINE FinalizeConvection()

    DEALLOCATE(ql3)
    DEALLOCATE(qi3)
    
    DEALLOCATE(ql2)
    DEALLOCATE(qi2)
    CALL FinalizeDeepConvection()
    CALL FinalizeShallowConvection()
    CALL FinalizeMicroPhysics()
  END SUBROUTINE FinalizeConvection
  !-----------------------------------------------------------------
  !-----------------------------------------------------------------
  SUBROUTINE CLOUD_COVER( &
      kt     ,ktp    ,ncols  ,kbot   ,ktop   ,noshal1,kctop1 , &
      kcbot1 ,rrr    ,fac2x  ,rccmbl ,iccon  ,convc  ,convt  , &
      convb  ,prcp1  ,prcp2  ,prcp3  ,prcpt  ,toplv  ,botlv  , &
      convts ,convcs ,convbs )
   IMPLICIT NONE 
    INTEGER, INTENT(in   ) :: kt     
    INTEGER, INTENT(in   ) :: ktp    
    INTEGER, INTENT(in   ) :: ncols  
    INTEGER, INTENT(in   ) :: kbot   (ncols)
    INTEGER, INTENT(in   ) :: ktop   (ncols)
    INTEGER, INTENT(in   ) :: noshal1(ncols)
    INTEGER, INTENT(in   ) :: kctop1 (ncols)
    INTEGER, INTENT(in   ) :: kcbot1 (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: rrr    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: fac2x  
    REAL(KIND=r8),    INTENT(in   ) :: rccmbl ! radiative convective cloud minimum base layer index
    CHARACTER(LEN=*),INTENT(in) :: iccon

    REAL(KIND=r8),    INTENT(inout) :: convc  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: convt  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: convb  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: prcp1  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: prcp2  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: prcp3  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: prcpt  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: toplv  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: botlv  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: convts (ncols)
    REAL(KIND=r8),    INTENT(inout) :: convcs (ncols)
    REAL(KIND=r8),    INTENT(inout) :: convbs (ncols)
     
    REAL(KIND=r8)   , PARAMETER :: fp2457 = 0.2457_r8
    REAL(KIND=r8)   , PARAMETER :: fp1253 = 0.1253_r8
    REAL(KIND=r8)   , PARAMETER :: f0p8   = 0.8_r8
    REAL(KIND=r8)   , PARAMETER :: f8p0e3 = 8.0e3_r8
    INTEGER :: i      
    INTEGER :: is     
    INTEGER :: ijk    

    !--------------------------------------------------------
    !   CLOUD COVER, Cloud TOP-BOT FOR RADIATION (sub cldgen)
    !             DUE CONV PRECIPITATION
    !--------------------------------------------------------
    !   TO cldgen are necessary:
    !   a) cloud top and base   (convt, convb in cldgen)
    !   b) cloud amount is calculated convc (only rrr>0). It is calculate below.
    !     a+b used to defined high clouds due to strong convection
    !   prcpt=precipitation at each time step.(rrr)
    !   convt=ktop
    !   conbt=kbot (>=2) for radiation
    !*****************************************************************

    IF (kt .NE. ktp) THEN
       DO IJK=1,ncols
          convc(IJK) = 0.0_r8   ! call reset(convc(1),ncols)
          convt(IJK) = 0.0_r8   ! call reset(convt(1),ncols)
          convb(IJK) = 0.0_r8   ! call reset(convb(1),ncols)
       ENDDO

       DO i = 1, ncols
          prcpt(i) = prcpt(i) - prcp1(i) &
               + prcp3(i)
       END DO
       IF(TRIM(iccon).EQ.'ARA')THEN
         DO i = 1, ncols 
           IF (prcpt(i) .GT. 0.0e0_r8) THEN
             convc(i) = 0.2_r8+0.038_r8*prcpt(i)*23000.0_r8
             convc(i) = MAX(convc(i), 0.0e0_r8)
             convc(i) = MIN(convc(i), f0p8)
           END IF
         END DO
       ELSE
         DO i = 1, ncols 
           IF (prcpt(i) .GT. 0.0e0_r8) THEN
             convc(i) = fp2457 + fp1253 * LOG(prcpt(i) * f8p0e3)
             convc(i) = MAX(convc(i), 0.0e0_r8)
             convc(i) = MIN(convc(i), f0p8)
           END IF
         END DO        
       END IF  
       !--faltou   
       DO i = 1, ncols
          IF (prcp3(i) .GT. 0.0e0_r8) THEN
             convt(i)=toplv(i) / prcp3(i)
             convb(i)=botlv(i) / prcp3(i)
          END IF
       END DO
       DO i = 1, ncols
          convb(i) = MAX(convb(i),rccmbl)
          IF (convb(i) .GT. convt(i)) &
               convb(i) = convt(i)
       END DO
       !-----
       DO i = 1, ncols
          prcp1(i) = prcp2(i)
          prcp2(i) = prcp3(i)
       END DO
       DO IJK=1,ncols
          prcp3(IJK) = 0.0_r8   !call reset(prcp3(1),ncols)
          toplv(IJK) = 0.0_r8   !call reset(toplv(1),ncols)
          botlv(IJK) = 0.0_r8   !call reset(botlv(1),ncols)
       ENDDO
    END IF
    !*****************************************************
    IF(TRIM(iccon).EQ.'ARA') THEN
      DO i = 1, ncols
        IF (rrr(i) .GT. 0.0_r8) THEN
          prcp3(i) = prcp3(i) + fac2x * rrr(i)
          toplv(i) = toplv(i) + fac2x * rrr(i) * ktop(i)
          botlv(i) = botlv(i) + fac2x * rrr(i) * kbot(i)
        END IF
      END DO
    ELSE
      DO i = 1, ncols
        IF (rrr(i) .GT. 0.0_r8) THEN
          prcp3(i) = prcp3(i) + fac2x * rrr(i)
          toplv(i) = toplv(i) + fac2x * rrr(i) * ktop(i)
          botlv(i) = botlv(i) + fac2x * rrr(i) * kbot(i)
        END IF
      END DO
    END IF
    !
    IF(TRIM(iccon).EQ.'ARA') THEN
       DO IJK=1,ncols
          convts(IJK) = 0.0_r8 
          convbs(IJK) = 0.0_r8 
          convcs(IJK) = 0.0_r8 
       ENDDO

       DO is=1,ncols
          IF(noshal1(is).EQ.0) THEN
             convts(is)=kctop1(is)
             convbs(is)=kcbot1(is)
             !
             !     for mass flux   convcs(is)=0.5
             !
             convcs(is)= 0.3_r8
          ENDIF
       END DO
    ENDIF
  END SUBROUTINE CLOUD_COVER

  SUBROUTINE ConvecGridHistStorage(&
                   ncols    ,kMax     ,latco    ,rdt      ,fac2     , &
                   rrr      ,Total_Rain,snowflg  ,sclhea   ,scmcha   , &
                   clheat   ,cmchan   ,lslhea   ,lsmcha) 
                   
    IMPLICIT NONE
    INTEGER, INTENT(IN   )    :: ncols
    INTEGER, INTENT(IN   )    :: kMax
    INTEGER, INTENT(IN   )    :: latco    
    REAL(KIND=r8), INTENT(IN   ) :: rdt
    REAL(KIND=r8), INTENT(IN   ) :: fac2
    REAL(KIND=r8), INTENT(IN   ) :: rrr       (nCols)
    REAL(KIND=r8), INTENT(in   ) :: Total_Rain (nCols)
    REAL(KIND=r8), INTENT(in   ) :: snowflg   (nCols)
    REAL(KIND=r8), INTENT(in   ) :: sclhea    (nCols,kMax)
    REAL(KIND=r8), INTENT(in   ) :: scmcha    (nCols,kMax)
    REAL(KIND=r8), INTENT(in   ) :: clheat    (nCols,kMax)
    REAL(KIND=r8), INTENT(in   ) :: cmchan    (nCols,kMax)
    REAL(KIND=r8), INTENT(in   ) :: lslhea    (nCols,kMax)
    REAL(KIND=r8), INTENT(in   ) :: lsmcha    (nCols,kMax)

    IF(TRIM(iccon).EQ.'ARA'.OR.TRIM(iccon).EQ.'KUO'.OR.TRIM(iccon).EQ.'GRE'.OR.&
       TRIM(iccon).EQ.'RAS'.OR.TRIM(iccon).EQ.'ZMC'.OR.TRIM(iccon).EQ.'GEC')THEN
       IF(dogrh(nGHis_cvprec,latco)) &
            CALL StoreGridHistory(rrr(1:ncols),nGHis_cvprec,latco,fac2*rdt*1.0e3_r8)
       IF(dogrh(nGHis_clheat,latco)) &
            CALL StoreGridHistory(clheat(1:nCols,:),nGHis_clheat,latco)
       IF(dogrh(nGHis_cvmosr,latco)) &
            CALL StoreGridHistory(cmchan(1:nCols,:),nGHis_cvmosr,latco)
    END IF

    IF(TRIM(ISCON).EQ.'TIED' .OR. TRIM(ISCON).EQ.'MFLX'.or. TRIM(ISCON).EQ.'SOUZ' .or. TRIM(ISCON).EQ.'JHK'.or. TRIM(ISCON).EQ.'UW')THEN
       IF(dogrh(nGHis_sclhea,latco)) &
            CALL StoreGridHistory(sclhea,nGHis_sclhea,latco)
       IF(dogrh(nGHis_shcvmo,latco)) &
            CALL StoreGridHistory(scmcha,nGHis_shcvmo,latco)
    END IF
    !---------------
    !     gdivn,gtmpn,grotn,gun,gvn are temporary working space
    !     
    IF(TRIM(ILCON).EQ.'YES'.or.TRIM(ILCON).EQ.'LSC'.or.TRIM(ILCON).EQ.'MIC'.or. TRIM(ILCON).EQ.'HWRF'.or. &
       TRIM(ILCON).EQ.'HGFS'.or.TRIM(ILCON).EQ.'UKMO' .or.TRIM(ILCON).EQ.'MORR'.or.TRIM(ILCON).EQ.'HUMO' .or.TRIM(ILCON).EQ.'HUMN') THEN
       IF(dogrh(nGHis_toprec,latco)) &
            CALL StoreGridHistory(Total_Rain,nGHis_toprec,latco,fac2*rdt*1.0e3_r8)    
       IF(dogrh(nGHis_snowfl,latco)) &
            CALL StoreGridHistory(snowflg  ,nGHis_snowfl,latco,fac2*rdt*1.0e3_r8)
       IF(dogrh(nGHis_sslaht,latco)) &
            CALL StoreGridHistory(lslhea   ,nGHis_sslaht,latco)
       IF(dogrh(nGHis_spstms,latco)) &
            CALL StoreGridHistory(lsmcha   ,nGHis_spstms,latco)
    END IF
  END SUBROUTINE ConvecGridHistStorage

  SUBROUTINE ConvecDiagnStorage(&
                   ncols    ,kMax     ,latco    ,rdt      ,fac2     , &
                   fdqn     ,rrr      ,Total_Rain,snowfl   ,sclhea   , &
                   scmcha   ,clheat   ,cmchan   ,lslhea   ,lsmcha     , &
                   tracerLiq ,tracerIce,cape,cine,LCL,LFC,SWEAT       ,&
                   gvarp)    
                                      
    IMPLICIT NONE
    INTEGER, INTENT(IN   )    :: ncols
    INTEGER, INTENT(IN   )    :: kMax
    INTEGER, INTENT(IN   )    :: latco    
    REAL(KIND=r8), INTENT(IN   ) :: rdt
    REAL(KIND=r8), INTENT(IN   ) :: fac2
    REAL(KIND=r8), INTENT(IN   ) :: fdqn      (nCols,kMax)
    REAL(KIND=r8), INTENT(IN   ) :: rrr       (nCols)
    REAL(KIND=r8), INTENT(in   ) :: Total_Rain (nCols)
    REAL(KIND=r8), INTENT(in   ) :: snowfl    (nCols)
    REAL(KIND=r8), INTENT(in   ) :: sclhea    (nCols,kMax)
    REAL(KIND=r8), INTENT(in   ) :: scmcha    (nCols,kMax)
    REAL(KIND=r8), INTENT(in   ) :: clheat (nCols,kMax)
    REAL(KIND=r8), INTENT(in   ) :: cmchan (nCols,kMax)
    REAL(KIND=r8), INTENT(in   ) :: lslhea (nCols,kMax)
    REAL(KIND=r8), INTENT(in   ) :: lsmcha (nCols,kMax)
    REAL(KIND=r8), INTENT(in   ) :: tracerLiq (nCols,kMax)
    REAL(KIND=r8), INTENT(in   ) :: tracerIce(nCols,kMax)    
    REAL(KIND=r8), INTENT(in   ) :: cape(ncols) 
    REAL(KIND=r8), INTENT(in   ) :: cine (ncols) 
    REAL(KIND=r8), INTENT(in   ) :: LCL (ncols) 
    REAL(KIND=r8), INTENT(in   ) :: LFC (ncols) 
    REAL(KIND=r8), INTENT(in   ) :: SWEAT(ncols) 
    REAL(KIND=r8),OPTIONAL,   INTENT(IN) :: gvarp (nCols,kmax,nClass+nAeros)

    REAL(KIND=r8) :: qc (1:nCols, 1:kMax)
    REAL(KIND=r8) :: qr (1:nCols, 1:kMax)
    REAL(KIND=r8) :: qi (1:nCols, 1:kMax)
    REAL(KIND=r8) :: qs (1:nCols, 1:kMax)
    REAL(KIND=r8) :: qg (1:nCols, 1:kMax)
    REAL(KIND=r8) :: ni (1:nCols, 1:kMax)
    REAL(KIND=r8) :: ns (1:nCols, 1:kMax)
    REAL(KIND=r8) :: nr (1:nCols, 1:kMax)
    REAL(KIND=r8) :: NG (1:nCols, 1:kMax)  
    REAL(KIND=r8) :: nc (1:nCols, 1:kMax)
    REAL(KIND=r8) :: co2 (1:nCols, 1:kMax)

    REAL(KIND=r8)    :: bfr1   (nCols)
    REAL(KIND=r8)    :: bfr3   (nCols)
    INTEGER :: i,k,j,itr
    bfr1=LCL*0.0_r8 
    bfr3=LFC*0.0_r8 
    qr=0.0;qs=0.0;qg=0.0;ni=0.0;ns=0.0
    nr=0.0;NG=0.0;nc=0.0
    IF (microphys) THEN
       DO k=1,kMax
          DO i=1,nCols
             IF((nClass+nAeros)>0 .and. PRESENT(gvarp))THEN
                qr     (i,k) = gvarp(i,k,1)
                IF( TRIM(ILCON).EQ.'MORR'.or.TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HUMN') THEN
                   qs     (i,k) = gvarp(i,k,2)
                   qg     (i,k) = gvarp(i,k,3)
                   NI     (i,k) = gvarp(i,k,4)
                   NS     (i,k) = gvarp(i,k,5)
                   NR     (i,k) = gvarp(i,k,6)
                   NG     (i,k) = gvarp(i,k,7)
                   IF( TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HUMN') NC     (i,k) = gvarp(i,k,8)
                END IF
             END IF 
          END DO
       END DO
       
       IF((nClass+nAeros)>0 .and. PRESENT(gvarp))THEN
          DO itr=nClass+1,nClass+nAeros
             DO k=1,kMax
                DO i=1,nCols
                   co2     (i,k) = gvarp(i,k,itr)
                END DO
             END DO
          END DO
       END IF 

    END IF
!    !   nDiag_qrmicr      , & ! = 230 ! QR - rain water mixing ratio  (kg/kg)
!    IF(dodia(nDiag_qrmicr))CALL updia(qr,nDiag_qrmicr,latco)
!    !   nDiag_qsmicr      , & ! = 231 ! QS - snow mixing ratio (kg/kg)
!    IF(dodia(nDiag_qsmicr))CALL updia(qs,nDiag_qsmicr,latco)
!    !   nDiag_qgmicr      , & ! = 232 ! QG - graupel mixing ratio (KG/KG)
!    IF(dodia(nDiag_qgmicr))CALL updia(qg,nDiag_qgmicr,latco)
!    !   nDiag_nimicr      , & ! = 233 ! NI - cloud ice number concentration (1/kg)
!    IF(dodia(nDiag_nimicr))CALL updia(NI,nDiag_nimicr,latco)
!    !   nDiag_nsmicr      , & ! = 234 ! NS - Snow Number concentration (1/kg)
!    IF(dodia(nDiag_nsmicr))CALL updia(NS,nDiag_nsmicr,latco)
!    !   nDiag_ncmicr      , & ! = 235 ! NC - Cloud droplet Number concentration (1/kg)
!    IF(dodia(nDiag_ncmicr))CALL updia(NC,nDiag_ncmicr,latco)
!    !   nDiag_nrmicr      , & ! = 236 ! NR - Rain Number concentration (1/kg)
!    IF(dodia(nDiag_nrmicr))CALL updia(NR,nDiag_nrmicr,latco)
!    !   nDiag_ngmicr      , & ! = 237 ! NG - Graupel number concentration (1/kg)
!    IF(dodia(nDiag_ngmicr))CALL updia(NG,nDiag_ngmicr,latco)
!    !   nDiag_co2aer          ! mixC02 - CO2 concentration (kg/kg)
!    IF(dodia(nDiag_co2aer))CALL updia(co2,nDiag_co2aer,latco)

    ! "negative specific humidity" correction
    IF(dodia(nDiag_nshcrm))CALL updia(fdqn,nDiag_nshcrm,latco)
    IF(dodia(nDiag_trcliq))CALL updia(tracerLiq,nDiag_trcliq,latco)
    IF(dodia(nDiag_trcice))CALL updia(tracerIce,nDiag_trcice,latco)
    IF(dodia(nDiag_cape2d))CALL updia(cape,nDiag_cape2d,latco)
    IF(dodia(nDiag_cine2d))CALL updia(cine ,nDiag_cine2d,latco)
    IF(dodia(nDiag_sweath))CALL updia(SWEAT,nDiag_sweath,latco)

    IF(TRIM(iccon).EQ.'ARA'.OR.TRIM(iccon).EQ.'KUO'.OR.TRIM(iccon).EQ.'GRE'.OR.&
       TRIM(iccon).EQ.'RAS'.OR.TRIM(iccon).EQ.'ZMC'.OR.TRIM(iccon).EQ.'GEC')THEN
       IF(dodia(nDiag_cvprec)) THEN
          DO i=1,ncols
             bfr1(i)=fac2*rdt*1.0e3_r8*rrr(i)
          END DO
          CALL updia(bfr1,nDiag_cvprec,latco)
       END IF
       IF(dodia(nDiag_clheat))CALL updia(clheat,nDiag_clheat,latco)
       IF(dodia(nDiag_cmchan))CALL updia(cmchan,nDiag_cmchan,latco)
    END IF

    IF(TRIM(ISCON).EQ.'TIED' .OR. TRIM(ISCON).EQ.'MFLX'.or. TRIM(ISCON).EQ.'SOUZ' .or. TRIM(ISCON).EQ.'JHK'.or. TRIM(ISCON).EQ.'UW')THEN
       IF(dodia(nDiag_sclhea))CALL updia(sclhea,nDiag_sclhea,latco)
       IF(dodia(nDiag_scmcha))CALL updia(scmcha,nDiag_scmcha,latco)
    END IF
    !---------------
    !     gdivn,gtmpn,grotn,gun,gvn are temporary working space
    !     
    IF(TRIM(ILCON).EQ.'YES'.or.TRIM(ILCON).EQ.'LSC'.or.TRIM(ILCON).EQ.'MIC'.or. TRIM(ILCON).EQ.'HWRF'.or. &
       TRIM(ILCON).EQ.'HGFS'.or.TRIM(ILCON).EQ.'UKMO' .or.TRIM(ILCON).EQ.'MORR'.or.TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HUMN') THEN
       DO i=1,ncols
          IF(dodia(nDiag_toprec))bfr1(i)=fac2*rdt*1.0e3_r8*Total_Rain(i)
          IF(dodia(nDiag_lsprec))bfr3(i)=fac2*rdt*1.0e3_r8*(Total_Rain(i)-rrr(i))
       END DO

       IF(dodia(nDiag_toprec))CALL updia(bfr1,nDiag_toprec,latco)

       IF(dodia(nDiag_snowfl))CALL updia(snowfl,nDiag_snowfl,latco)

       IF(dodia(nDiag_lsprec))CALL updia(bfr3,nDiag_lsprec,latco)

       IF(dodia(nDiag_lslhea))CALL updia(lslhea,nDiag_lslhea,latco)

       IF(dodia(nDiag_lsmcha))CALL updia(lsmcha,nDiag_lsmcha,latco)
    END IF
  END SUBROUTINE ConvecDiagnStorage


  ! qnegat : routine for dealing with negative values of specific humidity
  !          for data on latitude circle.



  SUBROUTINE  qnegat (fq, fdq, fft, rdt, del, iMax, kMax)
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
    REAL(KIND=r8),    INTENT(inout) :: fft  (iMax,kMax)   
    REAL(KIND=r8),    INTENT(in   ) :: del  (kMax)

    REAL(KIND=r8)   :: dfact(kMax)

    INTEGER :: klev
    INTEGER :: kblw
    INTEGER :: i
    INTEGER :: k  

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

    DO k=1,kMax
       DO i=1,iMax
          fft(i,k)=fft(i,k)/(1.0_r8+delq*fq(i,k))
       END DO
    END DO

    IF(dodia(nDiag_nshcrm))THEN
       DO k=1,kMax
          DO i=1,iMax
             fdq(i,k)=fdq(i,k)*rdt
          END DO
       END DO
    END IF

  END SUBROUTINE qnegat

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
  SUBROUTINE qpart(kMax,iMax,t,ps,sl,q,ql,qi,opt)
    IMPLICIT NONE
    INTEGER      , INTENT(IN   ) :: kMax
    INTEGER      , INTENT(IN   ) :: iMax
    REAL(KIND=r8), INTENT(IN   ) :: t    (iMax,kMax)
    REAL(KIND=r8), INTENT(IN   ) :: ps   (iMax)
    REAL(KIND=r8), INTENT(IN   ) :: sl   (kMax)
    CHARACTER(LEN=*), INTENT(IN   ) :: opt

    REAL(KIND=r8), INTENT(INOUT) :: q    (iMax,kMax)
    REAL(KIND=r8), INTENT(INOUT) :: ql   (iMax,kMax)
    REAL(KIND=r8), INTENT(INOUT) :: qi   (iMax,kMax)
    REAL(KIND=r8), PARAMETER :: tmelt      = 273.16_r8             ! freezing T of fresh water ~ K
    REAL(KIND=r8), PARAMETER :: tmax_fice  = tmelt     - 10.0_r8   ! max temperature for cloud ice formation
    REAL(KIND=r8), PARAMETER :: tmin_fice  = tmax_fice - 30.0_r8   ! min temperature for cloud ice formation
    REAL(KIND=r8), PARAMETER :: tmax_fsnow = tmelt                 ! max temperature for transition to convective snow
    REAL(KIND=r8), PARAMETER :: tmin_fsnow = tmelt-5.0_r8           ! min temperature for transition to convective snow

    REAL(KIND=r8)                :: qes(iMax,kMax) 
    REAL(KIND=r8)                :: tsp(iMax,kMax) 
    REAL(KIND=r8)                :: p  (iMax,kMax)
    INTEGER ::i,k
    qi = 0.0_r8
    ql = 0.0_r8
    DO k=1,kMax
       DO i=1,iMax
              p  (i,k) = ps(i)*sl(k)*1000_r8              ! pressure in Pa
       END DO
    END DO   
    CALL findsp (iMax,kMax, q, t, p, tsp, qes)

    IF(TRIM(opt)=='subtraction')THEN
       DO k=1,kMax
          DO i=1,iMax
             !p  (i,k) = ps(i)*sl(k)*10.0_r8               ! pressure in mbar
             !
             ! sgb - IPH is for phase, dependent on TCRIT (water or ice)
             ! calculation of the pressure vapor
             !
             !esft=es5(t(i,k))
             !qes(i,k) = 0.622_r8*esft/(100.0_r8*p(i,k)-esft)
             IF(qes(i,k) <= 1.0e-12_r8  )qes(i,k)=1.0e-12_r8

             IF(q(i,k)   >  qes(i,k))THEN
               ql(i,k) = q (i,k) - qes(i,k)              
               q (i,k) = qes(i,k)
               IF (tsp(i,k) < tmin_fice)THEN
                  qi(i,k) = ql(i,k)
                  ql(i,k) = 0.0_r8
               ELSE
                  qi(i,k) = 0.0_r8
               END IF

             ELSE
               ql(i,k) = 0.0_r8
               qi(i,k) = 0.0_r8
             END IF
          END DO
       END DO
    ELSE IF (TRIM(opt)=='addition')THEN
       DO k=1,kMax
          DO i=1,iMax
               IF(qes(i,k) <= 1.0e-12_r8  )qes(i,k)=1.0e-12_r8
               !IF(q(i,k)   <  qes(i,k))THEN
                  q(i,k) = q (i,k) + ql(i,k) + qi(i,k)                 
                  !IF(q(i,k) > qes(i,k))THEN
                  !  q(i,k) = qes(i,k)
                  !END IF
                  !ql(i,k) = 0.0_r8
                  !qi(i,k) = 0.0_r8
               !END IF          
          END DO
       END DO
    END IF
  END SUBROUTINE qpart
  !---------------------------------
  REAL(KIND=r8) FUNCTION es5(t)
    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN) :: t
    REAL(KIND=r8)   , PARAMETER :: tcrit    =   273.15_r8
    REAL(KIND=r8)   , PARAMETER :: cp   =1004.0_r8
    REAL(KIND=r8)   , PARAMETER :: xl   =2.5e06_r8
    REAL(KIND=r8)   , PARAMETER :: rv   =461.9_r8
    REAL(KIND=r8)            :: ae  (2)
    REAL(KIND=r8)            :: be  (2)
    REAL(KIND=r8)            :: ht  (2)
    ht(1)=xl/cp
    ht(2)=2.834e6_r8/cp
    be(1)=0.622_r8*ht(1)/0.286_r8    
    ae(1)=be(1)/273.0_r8+LOG(610.71_r8)    
    be(2)=0.622_r8*ht(2)/0.286_r8  
    ae(2)=be(2)/273.0_r8+LOG(610.71_r8)

    IF (t <= tcrit) THEN
       es5 = EXP(ae(2)-be(2)/t)
    ELSE
       es5 = EXP(ae(1)-be(1)/t)
    END IF

  END FUNCTION es5
END MODULE Convection
