!
!  $Author: pkubota $
!  $Date: 2008/04/09 12:42:57 $
!  $Revision: 1.19 $
!
MODULE FieldsPhysics

  USE Sizes, ONLY: &
       sl

  USE Constants, ONLY: &
       rk,r8,r4,i8,i4

  IMPLICIT NONE

  ! Gaussian fields: 28 3D , 12 2D and 12 1D

!  PRIVATE

  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:)   :: sigki ! Fator de conversao de temperatura potencial

  !---------------------------------------------------------------------------------------------------------------
  ! SHORT WAVE RADIATION
  !---------------------------------------------------------------------------------------------------------------

  ! Coeficiente de transporte vertical

  ! Viscosity (turbulencia)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: PBL_CoefKm ! m2/s momentum
  ! Scalar diffusivity (water)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: PBL_CoefKh ! m2/s water and heat

  ! Radiation fields at next integer hour
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: yVisBeam ! Down Sfc SW flux visible beam    (all-sky)  (W/m^2)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: yVisDiff ! Down Sfc SW flux visible diffuse (all-sky)  (W/m^2)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: yNirBeam ! Down Sfc SW flux Near-IR beam    (all-sky)  (W/m^2)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: yNirDiff ! Down Sfc SW flux Near-IR diffuse (all-sky)  (W/m^2)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: yVisBeamC! Down Sfc SW flux visible beam    (clear)  (W/m^2) 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: yVisDiffC! Down Sfc SW flux visible diffuse (clear)  (W/m^2) 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: yNirBeamC! Down Sfc SW flux Near-IR beam    (clear)  (W/m^2) 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: yNirDiffC! Down Sfc SW flux Near-IR diffuse (clear)  (W/m^2) 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: ySwHeatRate ! Heating rate due to shortwave         (K/s)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: ySwHeatRatec! Heating rate due to shortwave (clear) (K/s)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: ySwToaDown! Incident SW at top (W/m^2)                        
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: ySwSfcNet ! Abs Sfc SW 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: ySwSfcNetC! Abs Sfc SW (clear) 

  ! Radiation fields at last integer hour
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: rVisBeam ! Down Sfc SW flux visible beam    (all-sky)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: rVisDiff ! Down Sfc SW flux visible diffuse (all-sky)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: rNirBeam ! Down Sfc SW flux Near-IR beam    (all-sky)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: rNirDiff ! Down Sfc SW flux Near-IR diffuse (all-sky)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: rVisBeamC! Down Sfc SW flux visible beam    (clear) 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: rVisDiffC! Down Sfc SW flux visible diffuse (clear) 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: rNirBeamC! Down Sfc SW flux Near-IR beam    (clear) 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: rNirDiffC! Down Sfc SW flux Near-IR diffuse (clear) 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: rSwToaDown! Incident SW at top (W/m^2)               
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: rSwSfcNet ! Abs Sfc SW 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: rSwSfcNetC! Abs Sfc SW (clear) 

  ! Surface albedo for SW 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: AlbVisBeam ! Visible beam surface albedo
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: AlbVisDiff ! Visible diffuse surface albedo
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: AlbNirBeam ! Near-ir beam surface albedo
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: AlbNirDiff ! Near-ir diffuse surface albedo

  !---------------------------------------------------------------------------------------------------------------
  ! LONG WAVE RADIATION
  !---------------------------------------------------------------------------------------------------------------

  ! LW Radiation fields at last integer hour
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: LwCoolRate ! Cooling rate due to longwave          (K/s)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: LwCoolRatec! Cooling rate due to longwave  (clear) (K/s)

  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: LwSfcDown! Down Sfc LW flux         (W/m2)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: LwSfcDownC! Down Sfc LW flux (clear) (W/m2)

  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: LwSfcNet    ! Net Sfc LW         (W/m2)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: LwSfcNetC ! Net Sfc LW (clear) (W/m2)

  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: LwToaUp! Longwave upward at top           (W/m2)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: LwToaUpC! Longwave upward at top (clear)   (W/m2)

  !---------------------------------------------------------------------------------------------------------------
  ! CLOUDS FOR RADIATION
  !---------------------------------------------------------------------------------------------------------------

  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   :: cldsav! Total Cloud cover
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: cldtot! Total cloud cover (at each layer)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: cldinv! Inversion clouds
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: cldsat! Saturation clouds
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: cldcon! Convection clouds
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: cldson! Shallow convective clouds

  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: clwd  ! Cloud liquid water path. 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: emisd ! emissivity
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: taud  ! Shortwave cloud optical depth
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: EFFCS ! CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: EFFIS ! CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
  !---------------------------------------------------------------------------------------------------------------
  ! CONVECTION
  !---------------------------------------------------------------------------------------------------------------
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) ::  qliq !liquid water content in cloud
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) ::  dudt !zonal wind tendency due deep convection
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) ::  dvdt !meridional wind tendency due deep convection

  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   ppli ! Precipitation rate ( large scale ) ( mm/s)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   ppci ! Precipitation rate ( cumulus ) ( mm/s )
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   prct ! cumulus and large scale scheme precipitation  (mm)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   prcc ! cumulus scheme precipitation (mm)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   prcp1! precipitation (cumulus) at each time step.(rrr) ( mm/s )
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   prcp2! precipitation (cumulus) at each time step.(rrr) ( mm/s )
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   prcp3! precipitation (cumulus) at each time step.(rrr) ( mm/s )
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   prcpt! precipitation(cumulus and large scale) at each time step(rrr)(mm/s)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   toplv! level of convective cloud top
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   botlv! level of convective cloud base
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   geshem! cumulus scheme precipitation (mm)

  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: convc ! Convective cloud cover in 3 hr. avrage
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: convt ! Convective cloud top  (sigma layer)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: convb ! Convective cloud base (sigma layer)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: convts! Level sigma of the convective cloud cover in 3 hr. avrage
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: convcs! Level sigma of the convective cloud top  (sigma layer)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: convbs! Level sigma of the convective cloud base (sigma layer)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: concld  ! convective cloud cover
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: cld     ! cloud fraction

  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: o3mix ! ozone mass mixing ratio (g/g)

  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: tracermix ! tracer mass mixing ratio (g/g)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: sheleg! Snow amount in (mm) (equivalent water depth)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: cu_hr
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::cu_kbot
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::cu_ktop
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::cu_Kuo

  !---------------------------------------------------------------------------------------------------------------
  ! PBL
  !---------------------------------------------------------------------------------------------------------------

  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: ustr  ! Surface zonal stress
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: vstr  ! Surface meridional stress
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: tkemyj! (IN) (*) turb_g(ngrid)%tkep(mzp,mxp,myp)

  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: sflux_t ! (IN) (*) turb_g(ngrid)%sflux_t(mxp,myp)
                                                                            !fluxo de calor sensivel [K m/s]
                                                                !          turbulent kinetic energy [m^2/s^2]
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: sflux_r! (IN) (*) turb_g(ngrid)%sflux_r(mxp,myp)
                                                               !          fluxo de vapor d'agua [kg/kg m/s]
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: sflux_u! (IN) (*) turb_g(ngrid)%sflux_u(mxp,myp)
                                                               !          fluxo de momento [m/s m/s]
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: sflux_v! (IN) (*) turb_g(ngrid)%sflux_v(mxp,myp)
                                                               !          fluxo de momento [m/s m/s]
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::  tstar ! (IN) (*) leaf_g(ngrid)%tstar(mxp,myp,npatch)
                                                               !          scale de temp. turbulenta [K]
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::  wstar ! (IN) (*) leaf_g(ngrid)%tstar(mxp,myp,npatch)
                                                               !          scale of vertical veloc. turbulency [m/s]
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: gl0   ! Maximum mixing length l0 in blackerdar's formula (m)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: Mmlen ! Maximum mixing length l0 in blackerdar's formula (m)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: tauresx ! Residual stress to be added in vdiff to correct
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: tauresy ! for turb stress mismatch between sfc and atm accumulated.

  !---------------------------------------------------------------------------------------------------------------
  ! SURFACE/SSIB
  !---------------------------------------------------------------------------------------------------------------
  TYPE LSM
    REAL(KIND=r8),PUBLIC, POINTER :: vcover(:,:)
    REAL(KIND=r8),PUBLIC, POINTER :: uve10m(:,:)
    REAL(KIND=r8),PUBLIC, POINTER :: vve10m(:,:)
  END TYPE LSM
  TYPE(LSM),PUBLIC :: sfc

  INTEGER(KIND=i8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: imask ! vegetation mask
  INTEGER(KIND=i8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: mlsi  ! vegetation mask !add solange 27-01-2012
  REAL   (KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: ssib  ! Fracao de umidade do solo
  REAL   (KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: bstar
  REAL   (KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: r_aer ! (IN) (*) leaf_g(ngrid)r_aer(mxp,myp,npatch)
  REAL   (KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   tcm! Temperatura da copa "dossel"(K)   modificada   t-1
  REAL   (KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   tgm! Temperatura da superficie do solo  (K)   modificada  t-1
  REAL   (KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   tdm! Temperatura do solo profundo (K)   modificada   t-1
  REAL   (KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: wm !w0(id)Grau de saturacao de umid do solo id=1 na camada superficial t+1
                                                             !w0(id)Grau de saturacao de umid do solo id=2 na camada de raizes t+1
                                                             !w0(id)Grau de saturacao de umid do solo id=3 na camada de drenagem t+1
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: capacm! capac0(iv). Agua interceptada iv=1 no dossel"water store
                                                                !             capacity t+1 of leaves"(m)  modificada t-1
                                                                ! capac0(iv). Agua interceptada iv=2 na cobertura do solo(m)modific.

  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: capacc! capac0(iv). Agua interceptada iv=1 no dossel"water store
                                                                !             capacity t+1 of leaves"(m)  modificada t-1
                                                                ! capac0(iv). Agua interceptada iv=2 na cobertura do solo(m)modific.
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   tc0! "dossel" canopy of temparature (K) t+1
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   tg0! surface soil temperature (K)    t+1
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   td0! dep soil temperature (K) t+1
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: w0 !w0(id)Grau de saturacao de umid do solo id=1 na camada superficial t+1
                                                             !w0(id)Grau de saturacao de umid do solo id=2 na camada de raizes t+1
                                                             !w0(id)Grau de saturacao de umid do solo id=3 na camada de drenagem t+1
  !add solange 27-01-2012
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: sm0 !sm0(id)Conteudo de umid do solo id=1 na camada superficial t+1  (m3/m3)
                                                              !sm0(id)Conteudo de umid do solo id=2 na camada de raizes t+1 ( m3/m3)
                                                              !sm0(id)Conteudo de umid do solo id=3 na camada de drenagem t+1 (m3/m3)
  !fim add solange

  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: capac0! capac0(iv). Agua interceptada iv=1 no dossel "water store capacity
                                                                !             of leaves"(m)  modificada t+1
                                                                ! capac0(iv). Agua interceptada iv=2 na cobertura do solo(m)modificada

  !OCEAN
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: tseam !  IF (tseam < 0)  sea surface temp. (K)
                                                              !  IF (tseam > 0)  ground temp. (K)

  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: gtsea ! IF (tseam < 0)  sea surface temp. (K)
                                                              ! IF (tseam > 0)  ground temp. (K)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: zorl  ! Aero. Roughness length  (m)

  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: HML    
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: HUML 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: HVML 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: TSK  
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: z0sea
  !---------------------------------------------------------------------------------------------------------------
  ! GRAVITY WAVE
  !---------------------------------------------------------------------------------------------------------------
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   var  !Surface height variance (m**2)

  !---------------------------------------------------------------------------------------------------------------
  ! OUTROS...
  !---------------------------------------------------------------------------------------------------------------

  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   tm0! prognostic surface of temparature (K) t+1
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   qm0! prognostic surface of specific humid (kg/kg) t+1
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   tmm! prognostic surface of temparature (K)  t-1
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   qmm! prognostic surface of specific humid (kg/kg)  t-1


  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   qsfc0! prognostic surface ocean temparature (K) t+1
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   tsfc0! prognostic surface ocean temparature (K) t+1
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   qsfcm! prognostic surface ocean temparature (K) t-1
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   tsfcm! prognostic surface ocean temparature (K) t-1


  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   tg1 ! deep soil temp (K)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   tg2 ! ground temperature (K)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   tg3 ! canopy temperature (K)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   soilm ! total soil water in (mm)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   sens  ! sensible heat flux (w/m^2)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   evap  ! latent heat flux(W/m^2)

  INTEGER, PUBLIC, PARAMETER :: nzg     =8                !- total number of soil layers
  INTEGER, PUBLIC, PARAMETER :: npatches=5                !- total number of veg patches
  INTEGER, PUBLIC, PARAMETER :: npatches_actual=2         !- actual number of veg patches
  ! (must be =< npatches and >= 2)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: wsib3d   ! initial soil wetness data
  ! at sib soil layers
  !
  !--common gl_sm - copie as mesmas linhas na gloop_slgm.f90 ----
  !
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: veg_type ! SIB veg type
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: frac_occ ! fractional area
  ! coverage

  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: xland! land mask (1 for land, 2 for water)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: seamask
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: xice
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: z0
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: ustar
  INTEGER      , PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: lowlyr
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: F_ICE_PHY  
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: F_RAIN_PHY 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: F_RIMEF_PHY

  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: snow
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: THZ0
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: QZ0
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: UZ0
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: VZ0
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: ZNT
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: PBLH
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: tpert
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: qpert  
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: AKHS
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: AKMS
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: CT
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: htdisp
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: temp2m
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: umes2m
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: gndvi
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: ndvim
  INTEGER(KIND=i8), PUBLIC, ALLOCATABLE, DIMENSION(:,:)    :: MskAnt


  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: Dump
  
  PUBLIC ::  InitFieldsPhyscs
CONTAINS






  SUBROUTINE InitFieldsPhyscs(ibMax, kMax, jbMax)
    INTEGER, INTENT(IN) :: ibMax
    INTEGER, INTENT(IN) :: kMax
    INTEGER, INTENT(IN) :: jbMax
    INTEGER :: k
    !
    !   map(ib,jb)
    !
    ALLOCATE(sigki (kMax))
    DO k=1,kmax
       sigki (k)=1.0e0_r8/EXP(rk*LOG(sl(k)))
    END DO


    ! Coeficiente de transporte vertical para a turbulencia
    ALLOCATE(PBL_CoefKm (ibMax,kMax+1,jbMax))
    PBL_CoefKm = 0.0_r8
    ALLOCATE(PBL_CoefKh (ibMax,kMax+1,jbMax))
    PBL_CoefKh = 0.0_r8


    ! Radiation fields at next integer hour
    ALLOCATE(ySwToaDown(ibMax,jbMax))
    ySwToaDown=0.0_r8
    ALLOCATE(yVisBeam (ibMax,jbMax))
    yVisBeam =0.0_r8
    ALLOCATE(yVisDiff (ibMax,jbMax))
    yVisDiff =0.0_r8
    ALLOCATE(yNirBeam (ibMax,jbMax))
    yNirBeam =0.0_r8
    ALLOCATE(yNirDiff (ibMax,jbMax))
    yNirDiff =0.0_r8
    ALLOCATE(yVisBeamC(ibMax,jbMax))
    yVisBeamC=0.0_r8
    ALLOCATE(yVisDiffC(ibMax,jbMax))
    yVisDiffC=0.0_r8
    ALLOCATE(yNirBeamC(ibMax,jbMax))
    yNirBeamC=0.0_r8
    ALLOCATE(yNirDiffC(ibMax,jbMax))
    yNirDiffC=0.0_r8
    ALLOCATE(ySwHeatRate   (ibMax,kMax,jbMax))
    ySwHeatRate   =0.0_r8
    ALLOCATE(ySwHeatRatec  (ibMax,kMax,jbMax))
    ySwHeatRatec  =0.0_r8
    ALLOCATE(ySwSfcNet (ibMax,jbMax))
    ySwSfcNet =0.0_r8
    ALLOCATE(ySwSfcNetC(ibMax,jbMax))
    ySwSfcNetC=0.0_r8

    ! Radiation fields at last integer hour
    ALLOCATE(rSwToaDown(ibMax,jbMax))
    rSwToaDown=0.0_r8
    ALLOCATE(rVisBeam (ibMax,jbMax))
    rVisBeam =0.0_r8
    ALLOCATE(rVisDiff (ibMax,jbMax))
    rVisDiff =0.0_r8
    ALLOCATE(rNirBeam (ibMax,jbMax))
    rNirBeam =0.0_r8
    ALLOCATE(rNirDiff (ibMax,jbMax))
    rNirDiff =0.0_r8
    ALLOCATE(rVisBeamC(ibMax,jbMax))
    rVisBeamC=0.0_r8
    ALLOCATE(rVisDiffC(ibMax,jbMax))
    rVisDiffC=0.0_r8
    ALLOCATE(rNirBeamC(ibMax,jbMax))
    rNirBeamC=0.0_r8
    ALLOCATE(rNirDiffC(ibMax,jbMax))
    rNirDiffC=0.0_r8
    ALLOCATE(rSwSfcNet (ibMax,jbMax))
    rSwSfcNet =0.0_r8
    ALLOCATE(rSwSfcNetC(ibMax,jbMax))
    rSwSfcNetC=0.0_r8

    ! Surface albedo for SW 
    ALLOCATE(AlbVisBeam (ibMax,jbMax))
    AlbVisBeam =0.0_r8
    ALLOCATE(AlbVisDiff (ibMax,jbMax))
    AlbVisDiff =0.0_r8
    ALLOCATE(AlbNirBeam (ibMax,jbMax))
    AlbNirBeam =0.0_r8
    ALLOCATE(AlbNirDiff (ibMax,jbMax))
    AlbNirDiff =0.0_r8

    ! LW Radiation fields at last integer hour
    ALLOCATE(LwCoolRate   (ibMax,kMax,jbMax))
    LwCoolRate   =0.0_r8
    ALLOCATE(LwCoolRatec  (ibMax,kMax,jbMax))
    LwCoolRatec  =0.0_r8
    ALLOCATE(LwSfcDown(ibMax,jbMax))
    LwSfcDown=0.0_r8
    ALLOCATE(LwSfcDownC(ibMax,jbMax))
    LwSfcDownC=0.0_r8
    ALLOCATE(LwSfcNet(ibMax,jbMax))
    LwSfcNet   =0.0_r8
    ALLOCATE(LwSfcNetC (ibMax,jbMax))
    LwSfcNetC =0.0_r8
    ALLOCATE(LwToaUp(ibMax,jbMax))
    LwToaUp=0.0_r8
    ALLOCATE(LwToaUpC(ibMax,jbMax))
    LwToaUpC=0.0_r8

    ! Clouds for Radiation
    ALLOCATE(cldsav(ibMax,jbMax))
    cldsav=0.0_r8
    ALLOCATE(cldtot(ibMax,kMax,jbMax))
    cldtot=0.0_r8
    ALLOCATE(cldinv(ibMax,kMax,jbMax))
    cldinv=0.0_r8
    ALLOCATE(cldsat(ibMax,kMax,jbMax))
    cldsat=0.0_r8
    ALLOCATE(cldcon(ibMax,kMax,jbMax))
    cldcon=0.0_r8
    ALLOCATE(cldson(ibMax,kMax,jbMax))
    cldson=0.0_r8
    ALLOCATE(clwd  (ibMax,kMax,jbMax))
    clwd  =0.0_r8
    ALLOCATE(emisd (ibMax,kMax,jbMax))
    emisd =0.0_r8
    ALLOCATE(taud  (ibMax,kMax,jbMax))
    taud =0.0_r8 
    ALLOCATE(EFFCS  (ibMax,kMax,jbMax))
    EFFCS =0.0_r8 
    ALLOCATE(EFFIS  (ibMax,kMax,jbMax))
    EFFIS =0.0_r8 

    ! Convection    
    ALLOCATE(qliq  (ibMax,kMax,jbMax))
    qliq  =0.0_r8
    ALLOCATE(dudt  (ibMax,kMax,jbMax))
    dudt  =0.0_r8
    ALLOCATE(dvdt  (ibMax,kMax,jbMax))
    dvdt  =0.0_r8
    ALLOCATE(ppli  (ibMax,jbMax))
    ppli  =0.0_r8
    ALLOCATE(ppci  (ibMax,jbMax))
    ppci  =0.0_r8
    ALLOCATE(prct  (ibMax,jbMax))
    prct   =0.0_r8
    ALLOCATE(prcc  (ibMax,jbMax))
    prcc    =0.0_r8
    ALLOCATE(prcp1 (ibMax,jbMax))
    prcp1  =0.0_r8
    ALLOCATE(prcp2 (ibMax,jbMax))
    prcp2  =0.0_r8
    ALLOCATE(prcp3 (ibMax,jbMax))
    prcp3  =0.0_r8
    ALLOCATE(prcpt (ibMax,jbMax))
    prcpt   =0.0_r8
    ALLOCATE(toplv (ibMax,jbMax))
    toplv   =0.0_r8
    ALLOCATE(botlv (ibMax,jbMax))
    botlv  =0.0_r8
    ALLOCATE(geshem(ibMax,jbMax))
    geshem =0.0_r8
    ALLOCATE(convc (ibMax,jbMax))
    convc=0.0_r8
    ALLOCATE(convt (ibMax,jbMax))
    convt=0.0_r8
    ALLOCATE(convb (ibMax,jbMax))
    convb=0.0_r8
    ALLOCATE(convts(ibMax,jbMax))
    convts=0.0_r8
    ALLOCATE(convcs(ibMax,jbMax))
    convcs=0.0_r8
    ALLOCATE(convbs(ibMax,jbMax))
    convbs=0.0_r8
    ALLOCATE(concld (ibMax,kMax,jbMax))
    concld =0.0_r8
    ALLOCATE(cld    (ibMax,kMax,jbMax))
    cld =0.0_r8

    ! PBL
    ALLOCATE(ustr  (ibMax,jbMax))
    ustr=0.0_r8
    ALLOCATE(vstr  (ibMax,jbMax))
    vstr=0.0_r8

    ! Surface/SSIB
    ALLOCATE(imask (ibMax,jbMax))
    imask =0
    ALLOCATE(ssib  (ibMax,jbMax))
    ssib=0.0_r8
    ALLOCATE(bstar (ibMax,jbMax))
    bstar=0.0_r8
    ! Outros (falta organizar)
    ALLOCATE(gl0   (ibMax,jbMax))
    gl0   =0.0_r8
    ALLOCATE(Mmlen   (ibMax,jbMax))
    Mmlen   =0.0_r8
    ALLOCATE(tauresx(ibMax,jbMax))
    tauresx =0.0_r8
    ALLOCATE(tauresy(ibMax,jbMax))
    tauresy=0.0_r8
    ALLOCATE(zorl  (ibMax,jbMax))
    zorl  =0.0_r8
    ALLOCATE(sheleg(ibMax,jbMax))
    sheleg=0.0_r8

    ALLOCATE(cu_hr   (ibMax,kMax,jbMax)) 
    cu_hr=0.0_r8
    ALLOCATE(cu_kbot (ibMax,jbMax))
    cu_kbot=1
    ALLOCATE(cu_ktop (ibMax,jbMax))
    cu_ktop=-100
    ALLOCATE(cu_Kuo  (ibMax,jbMax))
    cu_Kuo=0
    ALLOCATE(tseam (ibMax,jbMax))
    tseam=0.0_r8
! add solange 27-01-2012
    ALLOCATE(mlsi(ibMax,jbMax))
    mlsi=0_i8
    ALLOCATE(sm0   (ibMax,3,jbMax))
    sm0=0.0_r8
! fim add


    ALLOCATE(gtsea (ibMax,jbMax))
    gtsea =0.0_r8

    ALLOCATE(HML  (ibMax,jbMax)) 
    HML   =0.0_r8
    ALLOCATE(HUML(ibMax,jbMax)) 
    HUML=0.0_r8
    ALLOCATE(HVML(ibMax,jbMax)) 
    HVML=0.0_r8
    ALLOCATE(TSK (ibMax,jbMax)) 
    TSK =0.0_r8
    ALLOCATE(z0sea(ibMax,jbMax)) 
    z0sea=0.0_r8
    
    ALLOCATE(tm0   (ibMax,jbMax))
    tm0   =0.0_r8
    ALLOCATE(qm0   (ibMax,jbMax))
    qm0   =0.0_r8
    ALLOCATE(tc0   (ibMax,jbMax))
    tc0   =0.0_r8
    ALLOCATE(tg0   (ibMax,jbMax))
    tg0   =0.0_r8
    
    ALLOCATE(td0   (ibMax,jbMax))
    td0   =0.0_r8
    ALLOCATE(tdm   (ibMax,jbMax))
    tdm   =0.0_r8    
    
    ALLOCATE(w0    (ibMax,3,jbMax))
    w0    =0.0_r8
    ALLOCATE(capac0(ibMax,2,jbMax))
    capac0=0.0_r8
    ALLOCATE(tmm   (ibMax,jbMax))
    tmm   =0.0_r8
    ALLOCATE(qmm   (ibMax,jbMax))
    qmm   =0.0_r8
    ALLOCATE(tcm   (ibMax,jbMax))
    tcm   =0.0_r8
    ALLOCATE(tgm   (ibMax,jbMax))
    tgm   =0.0_r8
    ALLOCATE(wm    (ibMax,3,jbMax))
    wm    =0.0_r8
    ALLOCATE(capacm(ibMax,2,jbMax))
    capacm=0.0_r8
    ALLOCATE(capacc(ibMax,2,jbMax))
    capacc=0.0_r8
    ALLOCATE(var   (ibMax,jbMax))
    var   =0.0_r8
    ALLOCATE(tg1   (ibMax,jbMax))
    tg1    =0.0_r8
    ALLOCATE(tg2   (ibMax,jbMax))
    tg2    =0.0_r8
    ALLOCATE(tg3   (ibMax,jbMax))
    tg3    =0.0_r8
    ALLOCATE(soilm (ibMax,jbMax))
    soilm  =0.0_r8
    ALLOCATE(sens  (ibMax,jbMax))
    sens  =0.0_r8
    ALLOCATE(evap  (ibMax,jbMax))
    evap  =0.0_r8
    ALLOCATE(tracermix (ibMax,kMax,jbMax))
    tracermix =0.0_r8

    ALLOCATE(o3mix (ibMax,kMax,jbMax))
    o3mix =0.0_r8
    ALLOCATE(wsib3d   (ibMax,jbMax,8       ))
    wsib3d=0.0_r8
    ALLOCATE(veg_type (ibMax,jbMax,npatches))
    veg_type=0.0_r8
    ALLOCATE(frac_occ (ibMax,jbMax,npatches))
    frac_occ=0.0_r8
    ALLOCATE( XLAND(ibMax,jbMax))
    xland=0.0_r8
    ALLOCATE( SEAMASK(ibMax,jbMax))
    seamask=0.0_r8
    ALLOCATE( XICE(ibMax,jbMax))
    xice=0.0_r8
    ALLOCATE( Z0(ibMax,jbMax))
    z0=0.0_r8
    ALLOCATE( USTAR(ibMax,jbMax))
    ustar=0.0_r8
    ALLOCATE( LOWLYR(ibMax,jbMax))
    lowlyr=0

    ALLOCATE(F_ICE_PHY  (ibMax,kMax,jbMax))
    F_ICE_PHY=0.0_r8
    ALLOCATE(F_RAIN_PHY (ibMax,kMax,jbMax))
    F_RAIN_PHY=0.0_r8
    ALLOCATE(F_RIMEF_PHY(ibMax,kMax,jbMax))
    F_RIMEF_PHY=0.0_r8
    ALLOCATE(tkemyj(ibMax,kMax+1,jbMax))
    tkemyj=0.02_r8
    ALLOCATE(sflux_t(ibMax,jbMax))
    sflux_t=0.0_r8
    ALLOCATE(sflux_r(ibMax,jbMax))
    sflux_r=0.0_r8
    ALLOCATE(sflux_u(ibMax,jbMax))
    sflux_u=0.0_r8
    ALLOCATE(sflux_v(ibMax,jbMax))
    sflux_v=0.0_r8
    ALLOCATE(tstar(ibMax,jbMax))
    tstar=0.0_r8
    ALLOCATE(wstar(ibMax,jbMax))
    wstar=0.0_r8
    ALLOCATE(r_aer(ibMax,jbMax))
    r_aer=0.0_r8
    ALLOCATE( qsfc0 (ibMax,jbMax))
    qsfc0=0.0_r8
    ALLOCATE( tsfc0 (ibMax,jbMax))
    tsfc0=0.0_r8
    ALLOCATE( qsfcm (ibMax,jbMax))
    qsfcm=0.0_r8
    ALLOCATE( tsfcm (ibMax,jbMax))
    tsfcm=0.0_r8
    ALLOCATE( snow (ibMax,jbMax))
    snow=0.0_r8
    ALLOCATE( THZ0 (ibMax,jbMax))
    THZ0=0.0_r8
    ALLOCATE( QZ0  (ibMax,jbMax))
    QZ0=0.0_r8
    ALLOCATE( UZ0  (ibMax,jbMax))
    UZ0=0.0_r8
    ALLOCATE( VZ0  (ibMax,jbMax))
    VZ0=0.0_r8
    ALLOCATE( ZNT  (ibMax,jbMax))
    ZNT=0.0_r8
    ALLOCATE( PBLH (ibMax,jbMax))
    PBLH=0.0_r8
    ALLOCATE( tpert (ibMax,jbMax))
    tpert=0.0_r8
    ALLOCATE( qpert (ibMax,jbMax)) 
    qpert=0.0_r8
    ALLOCATE(AKHS (ibMax,jbMax))
    AKHS=0.0_r8
    ALLOCATE(AKMS (ibMax,jbMax))
    AKMS=0.0_r8
    ALLOCATE(CT(ibMax,jbMax))
    CT=0.0_r8
    ALLOCATE(htdisp(ibMax,jbMax))
    htdisp=0.0_r8
    ALLOCATE(temp2m(ibMax,jbMax))
    temp2m=0.0_r8
    ALLOCATE(umes2m(ibMax,jbMax))
    umes2m=0.0_r8
    ALLOCATE(gndvi(ibMax,jbMax))
    gndvi=0.0_r8
    ALLOCATE(ndvim(ibMax,jbMax))
    ndvim=0.0_r8    
    ALLOCATE(MskAnt(ibMax,jbMax))
    MskAnt=0_i8    
    ALLOCATE(Dump(ibMax,kMax,jbMax))
    Dump=0.0_r8
  END SUBROUTINE InitFieldsPhyscs
END MODULE FieldsPhysics
