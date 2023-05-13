MODULE Sfc_SeaFlux_Interface
  USE Constants, ONLY :     &
       stefan,cp,hl,gasr,grav,z0ice,r8,i8

  USE Options, ONLY : &
       oml_hml0,atmpbl,sfcpbl,OCFLUX,omlmodel,ICEMODEL

  USE Sfc_SeaFlux_COLA_Model, ONLY : InitSeaIceFlux_COLA_Model,SeaIceFlux_COLA_Model
  USE Sfc_SeaFlux_WGFS_Model, ONLY : Init_Sfc_SeaFlux_WGFS_Model,SeaFlux_WGFS_Model,OCEANML
  USE Sfc_SeaFlux_UKME_Model, ONLY : InitSfc_SeaFlux_UKME_Model,SF_EXCH
  USE Sfc_SeaIceFlux_WRF_Model, ONLY : InitSeaIce,GetFluxSeaIceModel

  USE SlabOceanModel        , ONLY : InitGetOceanAlb

  IMPLICIT NONE
  PRIVATE   

  REAL(KIND=r8)   :: rbyg
  PUBLIC :: Init_Sfc_SeaFlux_Interface
  PUBLIC :: seasfc
CONTAINS    


  SUBROUTINE Init_Sfc_SeaFlux_Interface(ibMax,jbMax,kMax,ibMaxPerJB,si,sl,del_in,RESTART,fNameSlabOcen,fgtmp)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ibMax
    INTEGER, INTENT(IN) :: jbMax
    INTEGER, INTENT(IN) :: kMax
    INTEGER         , INTENT(IN   ) :: ibMaxPerJB(jbMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: si(kMax+1)
    REAL(KIND=r8)   , INTENT(IN   ) :: sl(kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: del_in(kMax)
    LOGICAL         , INTENT(IN   ) ::  RESTART
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameSlabOcen
    REAL(KIND=r8)   , INTENT(IN   ) :: fgtmp(ibMax,kMax,jbMax)

    rbyg  =gasr/grav*del_in(1)*0.5_r8
    CALL InitSeaIceFlux_COLA_Model(kmax, sl, del_in, si)
    CALL Init_Sfc_SeaFlux_WGFS_Model()
    CALL InitSfc_SeaFlux_UKME_Model(si,sl,del_in,kMax)
    CALL InitGetOceanAlb(fNameSlabOcen)
    CALL InitSeaIce(ibMax,jbMax,kMax,ibMaxPerJB,fgtmp,RESTART)
  END SUBROUTINE Init_Sfc_SeaFlux_Interface


  SUBROUTINE seasfc( &
       tmtx  ,  &
       umtx  ,  &
       qmtx  ,  &
       kpbl  ,  &
       kMax ,  &
       slrad ,  &
       tsurf ,  &
       qsurf ,  &
                                !
       gu    ,  &
       gv    ,  &
       gt    ,  &
       gq    ,  &
       gps   ,  &
       tsea  ,  &
       dtc3x ,  &
       sinclt,  &
                                !
       sigki ,  &
       delsig,  &
       sens  ,  &
       evap  ,  &
       umom  ,  &
       vmom  ,  &
       rmi   ,  &
       rhi   ,  &
                                !
       cond  ,  &
       stor  ,  &
       zorl  ,  &
       ncols ,  &
       speedm,  &
       bstar ,  &
       Ustarm,  &
       z0sea ,  &
                                !
       rho   ,  &
       qsfc  ,  &
       tsfc  ,  &
       mskant,  &
                                !
       iMask ,  &
       zenith,  &
       ppli  ,  &
       ppci  ,  &
       LwSfcDown,&
       xvisb,&
       xvisd,&
       xnirb,&
       xnird,&
       HML   ,  &
       HUML  ,  &
       HVML  ,  &
       TSK   ,  &
       GSW   ,  &
       GLW   ,  &
                                !
       cldtot,  &
       ySwSfcNet , &
       month2     ,&
       LwSfcNet  , &
       pblh      , &
       QCF       , &
       QCL       , &
                                !        
       mlsi      , &
       latco     , &
       Mmlen )
    IMPLICIT NONE
    !
    !==========================================================================
    ! ncols......Number of grid points on a gaussian latitude circle
    ! kpbl.......Number of layers pbl process is included( for u v,t )
    ! kMax......Number of layers pbl process is included( for q     )
    ! tmtx.......Temperature related matrix
    !            gmt(i,k,1)*d(gt(i,k-1))/dt+gmt(i,k,2)*d(gt(i,k))/dt=gmt(i,k,3)
    !            gmt(i,1,1)=0.
    !            gmt(*,*,1)...dimensionless
    !            gmt(*,*,2)...dimensionless
    !            gmt(*,*,3)...deg/sec
    ! umtx.......Wind related matrix
    !            gmu(i,k,1)*d(gu(i,k-1))/dt+gmu(i,k,2)*d(gu(i,k))/dt=gmu(i,k,3)
    !            gmu(i,k,1)*d(gv(i,k-1))/dt+gmu(i,k,2)*d(gv(i,k))/dt=gmu(i,k,4)
    !            gmu(i,1,1)=0.
    !            gmu(*,*,1)...dimensionless
    !            gmu(*,*,2)...dimensionless
    !            gmu(*,*,3)...m/sec**2
    !            gmu(*,*,4)...m/sec**2
    ! qmtx.......specific humidity related matrix
    !            gmq(i,k,1)*d(gq(i,k-1))/dt+gmq(i,k,2)*d(gq(i,k))/dt=gmq(i,k,3)
    !            gmq(i,1,1)=0.
    !            gmq(*,*,1)...dimensionless
    !            gmq(*,*,2)...dimensionless
    !            gmq(*,*,3)...kg/kg/sec
    ! slrad......radiation interpolation
    ! tsurff.....earth's surface temperature used for radiation
    !            for the first time step when ground temperature is not yet
    !            computed (this is done by subr.tsinit ),
    ! qsurf......qsurf(i)=0.622e0*EXP(21.65605e0 -5418.0e0 /tsurf(i))/gps(i)
    ! gu.........(zonal      velocity)*sin(colat)
    ! gv.........(meridional velocity)*sin(colat)
    ! gt.........Temperature
    ! gq.........Specific humidity
    ! gps........Surface pressure in mb
    ! tsea.......effective surface radiative temperature ( tgeff )
    ! dtc3x......time increment dt
    ! sinclt.....sinclt=SIN(colrad(latitu))
    ! sigki......sigki (k)=1.0e0/EXP(akappa*LOG(sig(k))),  where "sig"
    !            sigma coordinate at middle of layer and akappa=gasr/cp
    ! delsig
    ! sens.......sensible heat flux
    ! evap.......latent heat flux  "evaporation"
    ! umom.......umom(i)=fmom*um(ncount),
    !            where .fmom  momentum flux      in n/m**2
    !            fmom= rhoair(ncount)*cu(ncount)*ustar(ncount)
    !            um  (ncount)=gu (i,1)/sinclt
    !            gu          = (zonal velocity)*sin(colat)
    ! vmom.......vmom(i)=rho(i)*gv(i)*rmi(i)
    !            rho  (i)=gps(i)/(gr100*gt(i))
    !            gr100 =gasr*0.01
    ! z0ice.......Roughness length of ice
    ! rmi.........rmi   (i)=cu(i)*ustar(i), where
    !             cu is friction  transfer coefficients
    !             ustar is surface friction velocity  (m/s)
    ! rhi.........rhi   (i)=ct(i)*ustar(i), where
    !             ct is heat transfer coefficients.
    !             ustar is surface friction velocity  (m/s)
    ! cond........cond(i)=gice*(tsurf(i)-tice) or
    !             cond(i)=(2.03/2.0)*(tsurf(i)-271.16)
    ! stor........stor(i)=hscap*c0(i)
    ! zorl........zorl (i)= 100.0 *zgrav*speedm(i)*rhi(i)
    !             zgrav =0.032 /grav
    ! cp..........specific heat of air           (j/kg/k)
    ! hl..........heat of evaporation of water     (j/kg)
    ! gasr........gas constant of dry air        (j/kg/k)
    ! grav........grav   gravity constant        (m/s**2)
    ! stefan......Stefan Boltzman constant
    !==========================================================================
    ! 
    INTEGER      ,INTENT(in   ) :: ncols
    INTEGER      ,INTENT(in   ) :: latco
    INTEGER      ,INTENT(IN   ) :: kpbl
    INTEGER      ,INTENT(IN   ) :: kMax
    REAL(KIND=r8),INTENT(INOUT) :: tmtx  (ncols,kpbl,3)
    REAL(KIND=r8),INTENT(INOUT) :: umtx  (ncols,kpbl,4)
    REAL(KIND=r8),INTENT(INOUT) :: qmtx  (ncols,kMax,3)
    REAL(KIND=r8),INTENT(IN   ) :: slrad (ncols)
    REAL(KIND=r8),INTENT(INOUT) :: tsurf (ncols)
    REAL(KIND=r8),INTENT(IN   ) :: qsurf (ncols)
    REAL(KIND=r8),INTENT(IN   ) :: gu   (ncols,kMax)
    REAL(KIND=r8),INTENT(IN   ) :: gv   (ncols,kMax)
    REAL(KIND=r8),INTENT(INOUT) :: gt   (ncols,kMax)
    REAL(KIND=r8),INTENT(INOUT) :: gq   (ncols,kMax)
    REAL(KIND=r8),INTENT(IN   ) :: gps   (ncols)
    REAL(KIND=r8),INTENT(INOUT) :: tsea  (ncols)
    REAL(KIND=r8),INTENT(IN   ) :: dtc3x
    REAL(KIND=r8),INTENT(IN   ) :: sinclt(ncols)
    REAL(KIND=r8),INTENT(IN   ) :: sigki (1)
    REAL(KIND=r8),INTENT(IN   ) :: delsig(1:kMax)
    REAL(KIND=r8),INTENT(INOUT) :: sens  (ncols)
    REAL(KIND=r8),INTENT(INOUT) :: evap  (ncols)
    REAL(KIND=r8),INTENT(INOUT) :: umom  (ncols)
    REAL(KIND=r8),INTENT(INOUT) :: vmom  (ncols)
    REAL(KIND=r8),INTENT(INOUT) :: rmi   (ncols)
    REAL(KIND=r8),INTENT(INOUT) :: rhi   (ncols)
    REAL(KIND=r8),INTENT(INOUT) :: cond  (ncols)
    REAL(KIND=r8),INTENT(INOUT) :: stor  (ncols)
    REAL(KIND=r8),INTENT(INOUT) :: zorl  (ncols)
    REAL(KIND=r8)               :: rnet  (ncols)
    REAL(KIND=r8),INTENT(INOUT) :: speedm(ncols)
    REAL(KIND=r8),INTENT(OUT  ) :: Ustarm(ncols)
    REAL(KIND=r8),INTENT(OUT  ) :: rho   (ncols)
    REAL(KIND=r8),INTENT(INOUT) :: z0sea (ncols)
    REAL(KIND=r8),INTENT(INOUT) :: qsfc  (ncols)
    REAL(KIND=r8),INTENT(INOUT) :: tsfc (ncols)
    INTEGER(KIND=i8), INTENT(IN   ) :: mskant(ncols)
    REAL(KIND=r8),INTENT(IN   ) :: zenith(ncols)
    REAL(KIND=r8), INTENT(IN   ) :: ppli   (nCols)! Precipitation rate ( large scale )       (mm)
    REAL(KIND=r8), INTENT(IN   ) :: ppci   (nCols)! Precipitation rate ( cumulus )           (mm)
    INTEGER(KIND=i8), INTENT(IN   ) :: imask (nCols)
    INTEGER(KIND=i8), INTENT(INOUT) :: mlsi  (ncols)
    REAL(KIND=r8),    INTENT(OUT  ) :: bstar (ncols)

    REAL(KIND=r8),    INTENT(INOUT) :: HML  (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: HUML (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: HVML (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: TSK  (ncols)
    REAL(KIND=r8)    :: emisd(ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: GSW(ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: GLW(ncols)

    REAL(KIND=r8),    INTENT(IN   ) :: cldtot (ncols,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: ySwSfcNet(ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: LwSfcNet(ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: pblh(ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: QCF(ncols,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: QCL(ncols,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: LwSfcDown(1:nCols)
    REAL(KIND=r8),    INTENT(IN   ) :: xvisb(1:nCols)
    REAL(KIND=r8),    INTENT(IN   ) :: xvisd(1:nCols)
    REAL(KIND=r8),    INTENT(IN   ) :: xnirb(1:nCols)
    REAL(KIND=r8),    INTENT(IN   ) :: xnird(1:nCols)
    INTEGER      ,    INTENT(IN   ) :: month2(1:nCols)
    REAL(KIND=r8),    INTENT(IN   ) :: Mmlen(1:nCols)

    REAL(KIND=r8) :: z0      (ncols)
    REAL(KIND=r8) :: z0_cola (ncols)

    REAL(KIND=r8)    :: ah    (ncols)
    REAL(KIND=r8)    :: al    (ncols)
    REAL(KIND=r8)    :: am    (ncols)
    REAL(KIND=r8)    :: cuni  (ncols)
    REAL(KIND=r8)    :: cui   (ncols)
    REAL(KIND=r8)    :: cu    (ncols)
    REAL(KIND=r8)    :: ct    (ncols)
    REAL(KIND=r8)    :: ctni  (ncols)
    REAL(KIND=r8)    :: cti   (ncols)
    REAL(KIND=r8)    :: um    (ncols)
    REAL(KIND=r8)    :: vm    (ncols)
    REAL(KIND=r8)    :: tha   (ncols)
    REAL(KIND=r8)    :: thm   (ncols)
    REAL(KIND=r8)    :: dzm   (ncols)
    REAL(KIND=r8)    :: thvgm (ncols)
    REAL(KIND=r8)    :: rib   (ncols)
    REAL(KIND=r8)    :: ustar (ncols)
    REAL(KIND=r8)    :: gtsav (ncols)
    REAL(KIND=r8)    :: gqsav (ncols)
    REAL(KIND=r8)    :: tmsav (ncols)
    REAL(KIND=r8)    :: qmsav (ncols)
    REAL(KIND=r8)    :: tssav (ncols)
    REAL(KIND=r8)    :: dqg0  (ncols)
    REAL(KIND=r8)    :: b00   (ncols)
    REAL(KIND=r8)    :: b03   (ncols)
    REAL(KIND=r8)    :: b04   (ncols)
    REAL(KIND=r8)    :: c0    (ncols)
    REAL(KIND=r8)    :: b30   (ncols)
    REAL(KIND=r8)    :: b33   (ncols)
    REAL(KIND=r8)    :: c3    (ncols)
    REAL(KIND=r8)    :: b40   (ncols)
    REAL(KIND=r8)    :: b44   (ncols)
    REAL(KIND=r8)    :: c4    (ncols)

    INTEGER :: i
    INTEGER :: ncount
    REAL(KIND=r8)    :: gbyhl
    REAL(KIND=r8)    :: gbycp
    REAL(KIND=r8)    :: gr100
    REAL(KIND=r8)    :: gb100
    REAL(KIND=r8)    :: zgrav
    REAL(KIND=r8)    :: gice
    REAL(KIND=r8)    :: hscap
    REAL(KIND=r8)    :: sl1kap
    REAL(KIND=r8)    :: st4
    REAL(KIND=r8)    :: dti
    REAL(KIND=r8)    :: dtm
    REAL(KIND=r8)    :: dtmdt
    REAL(KIND=r8)    :: dqm
    REAL(KIND=r8)    :: dqmdt   
    !real(r8)         :: excess(ncols)     ! Excess downward sfc latent heat flux
    !*JPB REAL(KIND=r8), PARAMETER :: dd=0.05_r8
    REAL(KIND=r8), PARAMETER :: dd=3.0_r8 ! Total depth of the ice slab (m), Using ECMWF value
    REAL(KIND=r8), PARAMETER :: tice=271.16_r8
    REAL(KIND=r8), PARAMETER :: dice=2.0_r8
    REAL(KIND=r8), PARAMETER :: hice=2.03_r8
    REAL(KIND=r8), PARAMETER :: rhoice=920.0_r8 ! Mean ice density (kg/m3)
    REAL(KIND=r8), PARAMETER :: cice=2093.0_r8  ! Heat Capacity of Ice (J/Kg)

    REAL(KIND=r8) :: ICE_ELATEN(1:ncols)  
    REAL(KIND=r8) :: ICE_HFLUX (1:ncols) 
    REAL(KIND=r8) :: ICE_rmi   (1:ncols)   
    REAL(KIND=r8) :: ICE_rhi   (1:ncols)   
    REAL(KIND=r8) :: IC3_tsurf (1:ncols)   
    REAL(KIND=r8) :: rmi_uk  (ncols) 
    REAL(KIND=r8) :: rhi_uk  (ncols) 
    REAL(KIND=r8) :: evap_uk (ncols) 
    REAL(KIND=r8) :: sens_uk (ncols)    
    REAL(KIND=r8) :: ustar_uk(ncols) 
    REAL(KIND=r8) :: ustar_wgfs(ncols) 

    REAL(KIND=r8) :: rmi_cola  (ncols) 
    REAL(KIND=r8) :: rhi_cola  (ncols) 
    REAL(KIND=r8) :: evap_cola (ncols) 
    REAL(KIND=r8) :: sens_cola (ncols)    
    REAL(KIND=r8) :: ustar_cola(ncols) 
    REAL(KIND=r8) :: cu_cola    (ncols)
    REAL(KIND=r8) :: ct_cola    (ncols)

    REAL(KIND=r8) :: RADNET  (ncols)
    REAL(KIND=r8) ::  U3D    (1:nCols)
    REAL(KIND=r8) ::  V3D    (1:nCols)
    REAL(KIND=r8) ::  T3D    (1:nCols)
    REAL(KIND=r8) ::  QV3D   (1:nCols)
    REAL(KIND=r8) ::  P3D    (1:nCols)
    REAL(KIND=r8) ::  PSFC   (1:nCols)
    REAL(KIND=r8) ::  CHS    (1:nCols)
    REAL(KIND=r8) ::  CHS2   (1:nCols)
    REAL(KIND=r8) ::  CQS2   (1:nCols)
    REAL(KIND=r8) ::  CPM    (1:nCols)
    !REAL(KIND=r8) ::  ZNT   (1:nCols)
    !REAL(KIND=r8) ::  UST   (1:nCols)
    REAL(KIND=r8) ::  PSIM   (1:nCols)
    REAL(KIND=r8) ::  PSIH   (1:nCols)
    REAL(KIND=r8) ::  XLAND  (1:nCols)
    REAL(KIND=r8) ::  HFX    (1:nCols)
    REAL(KIND=r8) ::  QFX    (1:nCols)
    REAL(KIND=r8) ::  LH     (1:nCols)
    !REAL(KIND=r8) ::  TSK    (1:nCols)
    REAL(KIND=r8) ::  FLHC   (1:nCols)
    REAL(KIND=r8) ::  FLQC   (1:nCols)
    REAL(KIND=r8) ::  QGH    (1:nCols)
    REAL(KIND=r8) ::  QSFC_1 (1:nCols)
    REAL(KIND=r8) :: U10    (1:nCols)
    REAL(KIND=r8) :: V10    (1:nCols)
    REAL(KIND=r8) :: GZ1OZ0 (1:nCols)
    REAL(KIND=r8) :: WSPD   (1:nCols)
    REAL(KIND=r8) :: BR     (1:nCols)
    REAL(KIND=r8) :: CHS_SEA (1:nCols)
    REAL(KIND=r8) :: CHS2_SEA(1:nCols)
    REAL(KIND=r8) :: CPM_SEA   (1:nCols)
    REAL(KIND=r8) :: CQS2_SEA  (1:nCols)
    REAL(KIND=r8) :: FLHC_SEA  (1:nCols)
    REAL(KIND=r8) :: FLQC_SEA  (1:nCols)
    REAL(KIND=r8) :: HFX_SEA   (1:nCols)
    REAL(KIND=r8) :: LH_SEA    (1:nCols)
    REAL(KIND=r8) :: QFX_SEA   (1:nCols)
    REAL(KIND=r8) :: QGH_SEA   (1:nCols)
    REAL(KIND=r8) :: QSFC_SEA  (1:nCols)
    REAL(KIND=r8) :: UST_SEA   (1:nCols)
    REAL(KIND=r8) :: ZNT_SEA   (1:nCols)
    REAL(KIND=r8) :: CM_SEA(1:nCols)
    REAL(KIND=r8) :: CH_SEA(1:nCols)
    REAL(KIND=r8) :: WSPD_SEA(1:nCols)
    REAL(KIND=r8) :: SST       (1:nCols)
    REAL(KIND=r8) :: XICE      (1:nCols)
    REAL(KIND=r8) :: H0ML      (1:nCols)
    REAL(KIND=r8), PARAMETER ::  xice_threshold = 0.5_r8

    sens =0.0_r8
    evap =0.0_r8    
    rmi_uk  =0.0_r8
    rhi_uk  =0.0_r8
    evap_uk=0.0_r8
    sens_uk=0.0_r8
    ustar_uk=0.0_r8
    ustar_wgfs=0.0_r8
    rmi_cola  =0.0_r8
    rhi_cola  =0.0_r8
    evap_cola =0.0_r8
    sens_cola =0.0_r8
    ustar_cola=0.0_r8
    cu_cola=0.0_r8
    ct_cola=0.0_r8
    ICE_ELATEN=0.0_r8
    ICE_HFLUX =0.0_r8
    ICE_rmi   =0.0_r8
    ICE_rhi   =0.0_r8
    IC3_tsurf =0.0_r8
    umom =0.0_r8
    vmom =0.0_r8
    rnet =0.0_r8
    z0   =0.0_r8
    z0_cola =0.0_r8
    qsfc =0.0_r8
    tsfc =0.0_r8
    bstar=0.0_r8

    speedm=0.0_r8
    ah    =0.0_r8
    al    =0.0_r8
    am    =0.0_r8
    cuni  =0.0_r8
    cui   =0.0_r8
    cu    =0.0_r8
    ctni  =0.0_r8
    cti   =0.0_r8
    ct    =0.0_r8
    um    =0.0_r8
    vm    =0.0_r8
    tha   =0.0_r8
    thm   =0.0_r8
    dzm   =0.0_r8
    thvgm =0.0_r8
    rib   =0.0_r8
    ustar =0.0_r8
    gtsav =0.0_r8
    gqsav =0.0_r8
    tmsav =0.0_r8
    qmsav =0.0_r8
    tssav =0.0_r8
    dqg0  =0.0_r8
    b00   =0.0_r8
    b03   =0.0_r8
    b04   =0.0_r8
    c0    =0.0_r8
    b30   =0.0_r8
    b33   =0.0_r8
    c3    =0.0_r8
    b40   =0.0_r8
    b44   =0.0_r8
    c4    =0.0_r8
    Ustarm=0.0_r8
    rho   =0.0_r8

    gbyhl =0.0_r8
    gbycp =0.0_r8
    gr100 =0.0_r8
    gb100 =0.0_r8
    zgrav =0.0_r8
    gice =0.0_r8
    hscap =0.0_r8
    sl1kap =0.0_r8
    st4 =0.0_r8
    dti =0.0_r8
    dtm =0.0_r8
    dtmdt =0.0_r8
    dqm =0.0_r8
    dqmdt =0.0_r8

    gr100 =gasr*0.01_r8
    gbycp =grav/(cp*delsig(1)*100.0_r8 *sigki(1))
    gbyhl =grav/(hl*delsig(1)*100.0_r8 )
    gb100 =grav/(   delsig(1)*100.0_r8 )
    zgrav =0.032_r8 /grav
    gice  =hice/dice ! 2.03_r8/2.0_r8
    hscap =rhoice*cice*dd/dtc3x
    sl1kap=sigki(1)! sigki ! Fator de conversao de temperatura potencial sigki (k)=1.0e0_r8 / EXP(rk*LOG(sl(k)))
    st4   =stefan*4.0_r8
    dti   =1.0_r8 /dtc3x

    DO i = 1, ncols
       RADNET(i)=ySwSfcNet(i)+LwSfcNet(i)
       IF(mskant(i) == 1_i8)THEN
          rnet (i)=-697.58_r8*slrad(i)
          rho  (i)=gps(i)/(gr100*gt(i,1))
          ah   (i)=gbycp/gps(i)
          al   (i)=gbyhl/gps(i)
          dqg0 (i)=0.622_r8 *EXP(30.25353_r8 -5418.0_r8 /tsurf(i)) &
               /(tsurf(i)*tsurf(i)*gps(i))
          gtsav(i)=gt   (i,1)
          gqsav(i)=gq   (i,1)
          tssav(i)=tsurf(i)
          tmsav(i)=tmtx (i,1,3)
          qmsav(i)=qmtx (i,1,3)
       END IF
    END DO

    c0  =0.0_r8
    cond=0.0_r8
    stor=0.0_r8
    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN
          z0(i)=0.001_r8
          IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) >= 271.17_r8) THEN
             ! 
             ! Solution of sea water
             !
             z0(i)=0.01_r8*zorl(i)
          ELSE IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) < 271.17_r8) THEN
             !
             ! Solution of sea ice
             !
             z0(i)=z0ice
          END IF
       END IF
    END DO

    ncount=0
8000 CONTINUE
    ncount=ncount+1

    IF(TRIM(OCFLUX) == 'COLA' .OR. TRIM(OCFLUX) == 'WGFS' .OR. TRIM(OCFLUX) == 'UKME')THEN

       CALL  SeaIceFlux_COLA_Model ( &
                                ! Run Flags
            atmpbl    , &
                                ! Model information
            ncols     ,kmax      ,sigki(1)     ,mskant(1:nCols)    ,mlsi(1:nCols)      , &
            delsig(1:kMax) , &
                                ! Model Geometry,&
            sinclt(1:nCols), &
                                ! Time info
            dtc3x          ,&
                                ! Atmospheric fields
            gu    (1:ncols,1:kMax),gv    (1:ncols,1:kMax),gt    (1:ncols,1:kMax),gq  (1:ncols,1:kMax)  ,&
            gps   (1:nCols)       ,qsurf (1:nCols)       ,tsurf (1:nCols)       ,tsea(1:nCols)         ,&
            speedm(1:nCols)       ,&
                                ! SSIB: Total radiation absorbed at ground
            rnet  (1:nCols)       , &
                                ! Turbulence fields
            Mmlen   (1:nCols)         ,cu_cola    (1:nCols)           ,ct_cola(1:nCols)           , &
            rmi_cola(1:nCols)         ,rhi_cola   (1:nCols)           ,z0_cola(1:nCols)           , &
            zorl    (1:nCols)         ,ustar_cola (1:nCols)           ,tmtx   (1:nCols,1:kMax,1:3), &
            qmtx    (1:nCols,1:kMax,1:3),umtx       (1:nCols,1:kMax,1:4), &
                                ! Heat and Vapor Flux
            sens_cola(1:nCols)        ,evap_cola  (1:nCols)          ,&
                                ! Solving a system of linear equations by Gauss elimination 
            c3      (1:nCols)         ,b30        (1:nCols)          ,c0     (1:nCols)            , &
            b33     (1:nCols)         ,c4         (1:nCols)          ,b40    (1:nCols)            , &
            b44     (1:nCols))
       DO i = 1, ncols
          IF(mskant(i) == 1_i8)THEN
             ustar(i)=ustar_cola(i)
             z0   (i)=z0_cola(i)
             cu   (i)=cu_cola(i)
             ct   (i)=ct_cola(i)
          END IF
       END DO

    END IF
    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN
          gt  (i,1)   = gtsav(i)
          gq  (i,1)   = gqsav(i)
          tsurf(i)    = tssav(i)
          tmtx(i,1,3) = tmsav(i)
          qmtx(i,1,3) = qmsav(i)
       END IF
    END DO

    IF(TRIM(OCFLUX) == 'WGFS')THEN

       DO i = 1, ncols
          emisd(i) =0.98_r8   
          IF(mskant(i) == 1)THEN
             U3D   (i) = gu (i,1)/sinclt(i)
             V3D   (i) = gv (i,1)/sinclt(i)        
             T3D   (i) = gt(i,1)
             QV3D  (i) = gq(i,1)
             P3D   (i) = gps(i)*100.0_r8 - gps(i)*delsig(1)*100.0_r8
             PSFC  (i) = gps(i)*100.0_r8
             IF(iMask(i) >= 1_i8  )THEN ! land mask (1 for land, 2 for water , 13 ice)
                XLAND (i) =  2
             ELSE
                XLAND (i) =  2
             END IF
             IF(omlmodel)THEN
                H0ML  (i) = oml_hml0 - 13.5_r8*LOG(MAX(ABS(tsea(i))-tice+0.01_r8,1.0_r8))
                IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) < tice+0.01_r8) THEN
                   TSK   (i) = ABS(tsurf(i))
                END IF
             ELSE
                H0ML  (i) = 0.0_r8
                TSK   (i) = ABS(tsurf(i))
             END IF
             SST   (i) = ABS(tsea(i))!   REAL, INTENT(in ) ::, SST    (1:nCols)
             IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) < tice+0.01_r8) THEN
                XICE  (i) = xice_threshold !   REAL, INTENT(in ) ::, XICE   (1:nCols)
             ELSE
                XICE  (i) = 0.0_r8
             END IF
             z0sea    (i) = MAX(z0   (i),0.1_r8*z0ice)
          END IF
       END DO

       CALL SeaFlux_WGFS_Model( &
            U3D       , &!   REAL, INTENT(IN   ) ::  U3D   (1:nCols) !-- U3D     3D u-velocity interpolated to theta points (m/s)
            V3D       , &!   REAL, INTENT(IN   ) ::  V3D   (1:nCols) !-- V3D     3D v-velocity interpolated to theta points (m/s)
            T3D       , &!   REAL, INTENT(IN   ) ::  T3D   (1:nCols) !-- T3D     temperature (K)
            QV3D      , &!   REAL, INTENT(IN   ) ::  QV3D  (1:nCols) !-- QV3D    3D water vapor mixing ratio (Kg/Kg)
            P3D       , &!   REAL, INTENT(IN   ) ::  P3D   (1:nCols) !-- P3D     3D pressure (Pa)
            PSFC      , &!   REAL, INTENT(IN   ) ::  PSFC  (1:nCols)                    surface pressure (Pa)
            CHS       , &!   REAL, INTENT(OUT  ) ::  CHS   (1:nCols)
            CHS2      , &!   REAL, INTENT(OUT  ) ::  CHS2  (1:nCols)
            CQS2      , &!   REAL, INTENT(OUT  ) ::  CQS2  (1:nCols)
            CPM       , &!   REAL, INTENT(OUT  ) ::  CPM   (1:nCols)
            z0sea     , &!   REAL, INTENT(INOUT) ::  ZNT   (1:nCols)
            ustar_wgfs, &!   REAL, INTENT(INOUT) ::  UST   (1:nCols)
            PSIM      , &!   REAL, INTENT(OUT  ) ::  PSIM  (1:nCols)
            PSIH      , &!   REAL, INTENT(OUT  ) ::  PSIH  (1:nCols)
            XLAND     , &!   REAL, INTENT(IN   ) ::  XLAND (1:nCols)
            HFX       , &!   REAL, INTENT(OUT  ) ::  HFX   (1:nCols)
            QFX       , &!   REAL, INTENT(OUT  ) ::  QFX   (1:nCols)
            LH        , &!   REAL, INTENT(OUT  ) ::  LH    (1:nCols)
            TSK       , &!   REAL, INTENT(IN   ) ::  TSK   (1:nCols)
            FLHC      , &!   REAL, INTENT(OUT  ) ::  FLHC  (1:nCols)
            FLQC      , &!   REAL, INTENT(OUT  ) ::  FLQC  (1:nCols)
            QGH       , &!   REAL, INTENT(OUT  ) ::  QGH    (1:nCols),
            QSFC_1    , &!   REAL, INTENT(OUT  ) ::  QSFC    (1:nCols),
            U10       , &!   REAL, INTENT(OUT  ) ::  U10     (1:nCols),
            V10       , &!   REAL, INTENT(OUT  ) ::  V10     (1:nCols),
            GZ1OZ0    , &!   REAL, INTENT(OUT  ) ::  GZ1OZ0  (1:nCols),
            WSPD      , &!   REAL, INTENT(OUT  ) ::  WSPD    (1:nCols),
            BR        , &!   REAL, INTENT(OUT  ) ::  BR     (1:nCols),
            CHS_SEA   , &!   REAL, INTENT(OUT) ::, CHS_SEA  (1:nCols)
            CHS2_SEA  , &!   REAL, INTENT(OUT) ::, CHS2_SEA (1:nCols)
            CPM_SEA   , &!   REAL, INTENT(OUT) ::, CPM_SEA  (1:nCols)
            CQS2_SEA  , &!   REAL, INTENT(OUT) ::, CQS2_SEA (1:nCols)
            FLHC_SEA  , &!   REAL, INTENT(OUT) ::, FLHC_SEA (1:nCols)
            FLQC_SEA  , &!   REAL, INTENT(OUT) ::, FLQC_SEA (1:nCols)
            HFX_SEA   , &!   REAL, INTENT(OUT) ::, HFX_SEA  (1:nCols)
            LH_SEA    , &!   REAL, INTENT(OUT) ::, LH_SEA   (1:nCols)
            QFX_SEA   , &!   REAL, INTENT(OUT) ::, QFX_SEA  (1:nCols)
            QGH_SEA   , &!   REAL, INTENT(OUT) ::, QGH_SEA  (1:nCols)
            QSFC_SEA  , &!   REAL, INTENT(OUT) ::, QSFC_SEA (1:nCols)
            UST_SEA   , &!   REAL, INTENT(OUT) ::, UST_SEA  (1:nCols)
            ZNT_SEA   , &!   REAL, INTENT(OUT) ::, ZNT_SEA  (1:nCols)
            CM_SEA    , &!   REAL, INTENT(OUT) ::, ZNT_SEA  (1:nCols)
            CH_SEA    , &!   REAL, INTENT(OUT) ::, ZNT_SEA  (1:nCols)
            WSPD_SEA  , &!   REAL, INTENT(OUT) ::, ZNT_SEA  (1:nCols)
            SST       , &!   REAL, INTENT(in ) ::, SST      (1:nCols)
            XICE      , &!   REAL, INTENT(in ) ::, XICE     (1:nCols)
            mskant    , & 
            delsig(1) , &
            nCols )

       IF(omlmodel)THEN
          CALL OCEANML( &
               ustar    , & !REAL,     INTENT(IN   ) :: UST  ( 1:nCols )
               U3D      , & !REAL,     INTENT(IN   ) :: U_PHY( 1:nCols )
               V3D      , & !REAL,     INTENT(IN   ) :: V_PHY( 1:nCols )
               mskant   , & !REAL,     INTENT(IN   ) :: XLAND( 1:nCols )
               HFX_SEA  , & !REAL,     INTENT(IN   ) :: HFX  ( 1:nCols )
               LH_SEA   , & !REAL,     INTENT(IN   ) :: LH   ( 1:nCols )
               ABS(tsea), & !REAL,     INTENT(IN   ) :: tsea ( 1:nCols )
               TSK      , & !REAL,     INTENT(INOUT) :: TSK  ( 1:nCols )
               HML      , & !REAL,     INTENT(INOUT) :: HML  ( 1:nCols )
               HUML     , & !REAL,     INTENT(INOUT) :: HUML ( 1:nCols )
               HVML     , & !REAL,     INTENT(INOUT) :: HVML ( 1:nCols )
               GSW      , & !REAL,     INTENT(IN   ) :: GSW  ( 1:nCols )
               GLW      , & !REAL,     INTENT(IN   ) :: GLW  ( 1:nCols )
               emisd    , & !REAL,     INTENT(IN   ) :: EMISS( 1:nCols )
               dtc3x    , & !REAL,     INTENT(IN   ) :: DT
               H0ML     , &
               nCols      ) !INTEGER,  INTENT(IN   ) :: nCols
       END IF

    ELSE IF (TRIM(OCFLUX) == 'UKME')THEN

       DO i = 1, ncols
          IF(mskant(i) == 1_i8)THEN
             IF(z0sea  (i)<= 0.0_r8) z0sea    (i) = MAX(z0   (i),0.1_r8*z0ice)
          END IF
       END DO

       CALL  SF_EXCH (&
            nCols     , &
            kMAx      , &
            cldtot    , &
            QCF       , &
            QCL       , &
            gq        , &
            gt        , &
            gu        , &
            gv        , &
            sinclt    , &
            gps       , &
            radnet    , &
            tsurf     , &
            pblh      , &
            mskant    , &
            iMask     , &
            ABS(tsurf), & !TSTAR_LAND , &
            ABS(tsea) , & !TSTAR_SSI , &
            ABS(tsea) , & !TSTAR_SEA ,&
            ABS(tsea) , & !TSTAR_SICE ,&
            rmi_uk    , & 
            rhi_uk    , &
            evap_uk   , &
            sens_uk   , &
            ustar_uk  , &
            z0sea       & !Z0MSEA &
            )
       !    ELSE IF (OCFLUX == 'COLA')THEN

       !    ELSE
       !       WRITE(*,*)'ERRO seasfc',OCFLUX
       !       STOP

    END IF

    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN
          gt   (i,1)   = gtsav(i)
          gq   (i,1)   = gqsav(i)
          tsurf(i)     = tssav(i)
          tmtx (i,1,3) = tmsav(i)
          qmtx (i,1,3) = qmsav(i)
       END IF
    END DO

    IF (TRIM(ICEMODEL)=='SSIB')THEN

       !
       !-----------------------------------------------------------------------
       !**********************************************
       CALL GetFluxSeaIceModel (&
                                ! Model information
            nCols,kMax,    latco,    delsig(1:1), &
                                ! Model Geometry
            zenith(1:ncols)      ,sinclt(1:ncols)      ,mskant(1:ncols),&
                                ! Time info
            dtc3x                ,month2(1:ncols)      , &
                                ! Atmospheric fields
            tsea (1:ncols)       ,gu   (1:ncols,1:kMax),gv   (1:ncols,1:kMax), &
            gt   (1:ncols,1:kMax),gq   (1:ncols,1:kMax),gps  (1:ncols)       , &
            ppli (1:ncols)       ,ppci (1:ncols)       , &
                                ! LW Radiation fields at last integer hour
            LwSfcDown(1:ncols)   , &
                                ! Radiation field (Interpolated) at time = tod
            xvisb(1:ncols)       ,xvisd(1:ncols)       ,xnirb(1:ncols)       , &
            xnird(1:ncols)       , &
                                ! Surface Flux of SeaIce Scheme
            ICE_ELATEN(1:ncols)  ,ICE_HFLUX (1:ncols)  ,ICE_rmi (1:ncols)    , &
            ICE_rhi   (1:ncols)  ,IC3_tsurf (1:ncols))

    END IF

    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN
          IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) < tice+0.01_r8) THEN
             !
             ! Solution of sea ice
             !
             IF(TRIM(OCFLUX) == 'COLA')THEN
                IF (TRIM(ICEMODEL)=='COLA')THEN
                   IF(omlmodel)THEN
                      rhi  (i) = rhi_cola(i)
                      rmi  (i) = rmi_cola(i)
                   ELSE
                      rhi  (i) = rhi_cola(i)
                      rmi  (i) = rmi_cola(i)
                   END IF
                ELSE IF (TRIM(ICEMODEL)=='SSIB')THEN
                   rhi  (i) = ICE_rhi(i)
                   rmi  (i) = ICE_rmi(i)
                END IF
                tmtx(i,1,3)=(c3(i)-b30(i)*c0(i))/(b33(i)*dtc3x)
                qmtx(i,1,3)=(c4(i)-b40(i)*c0(i))/(b44(i)*dtc3x)
             ELSE IF (TRIM(OCFLUX) == 'WGFS')THEN
                IF (TRIM(ICEMODEL)=='COLA')THEN                   
                   IF(omlmodel)THEN
                      !                     rhi  (i) = CHS_SEA(i)
                      !                     rmi  (i) = WSPD_SEA(i)*(CM_SEA (i))
                      rhi  (i) = rhi_cola(i)
                      rmi  (i) = rmi_cola(i)
                   ELSE
                      !                     rhi  (i) = CHS_SEA(i)
                      !                     rmi  (i) = WSPD_SEA(i)*(CM_SEA (i))
                      rhi  (i) = rhi_cola(i)
                      rmi  (i) = rmi_cola(i)
                   END IF
                ELSE IF (TRIM(ICEMODEL)=='SSIB')THEN
                   rhi  (i) = ICE_rhi(i)
                   rmi  (i) = ICE_rmi(i)
                END IF
                tmtx(i,1,3)=(c3(i)-b30(i)*c0(i))/(b33(i)*dtc3x)
                qmtx(i,1,3)=(c4(i)-b40(i)*c0(i))/(b44(i)*dtc3x)
             ELSE IF (TRIM(OCFLUX) == 'UKME')THEN
                IF (TRIM(ICEMODEL)=='COLA')THEN
                   rhi  (i) = rhi_uk(i)
                   rmi  (i) = rmi_uk(i)
                ELSE IF (TRIM(ICEMODEL)=='SSIB')THEN
                   rhi  (i) = ICE_rhi(i)
                   rmi  (i) = ICE_rmi(i)
                END IF
             ELSE
                STOP 'SEA AND ICE FLUX'
             END IF
             mlsi(i)=2_i8
             !              dTs
             ! rho*ci*Dzi* ---- = Rdown - rho*cp*ck* (Ts - Tl*Sigl) - rho*L*ck*(qs -ql) - cs*(Ts - Ti) - Sigma*Ts^4
             !              dt
             !
             !              ddTs     Rdown     rho*cp*ck* d(Ts - Tl*Sigl)   rho*L*ck*d(qs -ql)     cs*d(Ts - Ti)     Sigma*dTs^4
             ! rho*ci*Dzi* ------- = -----  - -------------------------- - --------------------- - --------------- - ----------
             !              dTs*dt   dTs                   dTs                       dTs                dTs             dTs
             !
             !              d  dTs   Rdown                  rho*L*ck*d(qs -ql)
             ! rho*ci*Dzi* ------- = -----  - rho*cp*ck  - -------------------- -  cs - 4Sigma*dTs^3
             !              dt dTs    dTs                        dTs          
             !
             !b00(i)=   hscap+cp*rho(i)*rhi(i)+hl*rho(i)*rhi(i)*dqg0(i) + gice+st4*tsurf(i)**3
             !b03(i)=        -cp*rho(i)*rhi(i)*sl1kap
             !b04(i)=-hl*rho(i)*rhi(i)
             ! b00 + b10 + b20 + b30 + b40    c0
             ! b01 +             b33          c3     
             ! b02                     b44    c4
             ! b03
             ! b04 
             ! Right side of eq.41 section III.A 
             ! COLA Physics Description Manual
             !c0 (i)=rnet(i) -cp*rho(i)*rhi(i)*(tsurf(i)-sl1kap*gt(i,1)) &
             !     -hl*rho(i)*rhi(i)*(qsurf(i)-       gq(i,1)) &
             !     -gice*(tsurf(i)-tice)-stefan*tsurf(i)**4
             !b30(i)=               -ah (i)*cp*rho(i)*rhi(i)
             !b33(i)=tmtx(i,1,2)*dti-b30(i)*          sl1kap
             !c3 (i)=tmtx(i,1,3)    -b30(i)*(tsurf(i)-sl1kap*gt(i,1))
             !b40(i)=               -al(i)*hl*rho(i)*rhi(i)* dqg0 (i)
             !b44(i)=qmtx(i,1,2)*dti+al(i)*hl*rho(i)*rhi(i)
             !c4 (i)=qmtx(i,1,3)    + &
             !     al(i)*hl*rho(i)*rhi(i)*(qsurf(i)-gq(i,1))
             !b00(i)=b00(i)-b30(i)*b03(i)/b33(i)-b40(i)*b04(i)/b44(i)
             !c0 (i)=c0 (i)-c3 (i)*b03(i)/b33(i)-c4 (i)*b04(i)/b44(i)
             !c0 (i)=c0 (i)/b00(i)
             !IF (TRIM(ICEMODEL)=='COLA')THEN
             !   IF(atmpbl == 1) THEN
             !      tsurf(i)=tsurf(i)+c0(i)
             !   END IF
             !ELSE IF (TRIM(ICEMODEL)=='SSIB')THEN
             IF (TRIM(ICEMODEL)=='SSIB')THEN
                tsurf(i)=IC3_tsurf (i)
             END IF

          ELSE IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) > tice+0.01_r8) THEN
             !
             ! Solution of sea water
             !
             IF(TRIM(OCFLUX) == 'WGFS')THEN
                rhi  (i) = CHS_SEA(i)
                rmi  (i) = WSPD_SEA(i)*(CM_SEA (i))
                zorl (i)= 100.0_r8 *zgrav*speedm(i)*rhi(i)
                sens (i) = HFX_SEA(i)
                evap (i) = LH_SEA(i)
             ELSE IF (TRIM(OCFLUX) == 'UKME')THEN
                rhi  (i) = rhi_uk(i)
                rmi  (i) = rmi_uk(i)
                zorl (i) = 100.0_r8    * zgrav*speedm(i)*rhi(i)
                sens (i) = sens_uk (i)
                evap (i) = evap_uk (i)    
             ELSE IF (TRIM(OCFLUX) == 'COLA')THEN
                rhi  (i) = rhi_cola(i)
                rmi  (i) = rmi_cola(i)
                zorl (i)= 100.0_r8 *zgrav*speedm(i)*rhi(i)
                !sens (i)= rho(i)*cp*(tsurf(i)-gt(i,1)*sigki(1))*rhi(i)
                !evap (i)= rho(i)*hl*(qsurf(i)-gq(i,1)         )*rhi(i)
                sens (i)=  sens_cola (i)
                evap (i)=  evap_cola (i)
             ELSE
                STOP 'SEA WATERFLUX'
             END IF
             mlsi(i)=0_i8
             IF(atmpbl /= 1)THEN
                tmtx(i,1,3)=(ah(i)*sens(i))/(dtc3x*ah(i)*rho(i)*cp*rhi(i))
                qmtx(i,1,3)=(al(i)*evap(i))/(dtc3x*al(i)*rho(i)*hl*rhi(i))
             ELSE
                tmtx(i,1,3)=(tmtx(i,1,3)+ah(i)*sens(i))/(tmtx(i,1,2)+dtc3x*ah(i)*rho(i)*cp*rhi(i))
                qmtx(i,1,3)=(qmtx(i,1,3)+al(i)*evap(i))/(qmtx(i,1,2)+dtc3x*al(i)*rho(i)*hl*rhi(i))
             END IF
          END IF
       END IF
    END DO

    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN
          gt(i,1)=gt(i,1)+tmtx(i,1,3)*dtc3x
          gq(i,1)=gq(i,1)+qmtx(i,1,3)*dtc3x
       END IF
    END DO

    IF (ncount == 1) go to 8000

    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN
          IF(sfcpbl == 2)THEN
             IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) < tice+0.01_r8) THEN
                !
                ! Solution of sea ice
                !
                IF(TRIM(OCFLUX) == 'COLA')THEN
                   IF (TRIM(ICEMODEL)=='COLA')THEN
                      IF(omlmodel)THEN
                         rhi  (i) = rhi_cola(i)
                         rmi  (i) = rmi_cola(i)
                         sens (i) = sens_cola (i)
                         evap (i) = evap_cola (i)
                      ELSE
                         rhi  (i) = rhi_cola(i)
                         rmi  (i) = rmi_cola(i)
                         sens (i) = sens_cola (i)
                         evap (i) = evap_cola (i)
                      END IF
                   ELSE IF (TRIM(ICEMODEL)=='SSIB')THEN
                      IF(omlmodel)THEN
                         rhi  (i) = ICE_rhi(i)
                         rmi  (i) = ICE_rmi(i)
                         sens (i) = ICE_HFLUX (i) 
                         evap (i) = ICE_ELATEN(i)
                      ELSE
                         rhi  (i) = ICE_rhi(i)
                         rmi  (i) = ICE_rmi(i)
                         sens (i) = ICE_HFLUX (i) 
                         evap (i) = ICE_ELATEN(i)
                      END IF
                   END IF
                ELSE IF(TRIM(OCFLUX) == 'WGFS')THEN
                   IF (TRIM(ICEMODEL)=='COLA')THEN
                      IF(omlmodel)THEN
                         rhi  (i) = rhi_cola(i)
                         rmi  (i) = rmi_cola(i)
                         sens (i) = sens_cola (i)
                         evap (i) = evap_cola (i)
                      ELSE
                         rhi  (i) = rhi_cola(i)
                         rmi  (i) = rmi_cola(i)
                         sens (i) = sens_cola (i)
                         evap (i) = evap_cola (i)
                      END IF
                   ELSE IF (TRIM(ICEMODEL)=='SSIB')THEN
                      IF(omlmodel)THEN
                         rhi  (i) = ICE_rhi(i)
                         rmi  (i) = ICE_rmi(i)
                         sens (i) = ICE_HFLUX (i) 
                         evap (i) = ICE_ELATEN(i)
                      ELSE
                         rhi  (i) = ICE_rhi(i)
                         rmi  (i) = ICE_rmi(i)
                         sens (i) = ICE_HFLUX (i) 
                         evap (i) = ICE_ELATEN(i)
                      END IF
                   END IF
                ELSE IF (TRIM(OCFLUX) == 'UKME')THEN
                   IF (TRIM(ICEMODEL)=='COLA')THEN
                      rhi  (i) = rhi_uk(i)
                      rmi  (i) = rmi_uk(i)
                      sens (i) = sens_uk(i)
                      evap (i) = evap_uk(i)
                   ELSE IF (TRIM(ICEMODEL)=='SSIB')THEN
                      rhi  (i) = ICE_rhi(i)
                      rmi  (i) = ICE_rmi(i)
                      sens (i) = ICE_HFLUX(i)
                      evap (i) = ICE_ELATEN(i)
                   END IF
                ELSE 
                   STOP 'SEA AND ICE FLUX'
                END IF
             ELSE IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) > tice+0.01_r8) THEN
                !
                ! Solution of sea water
                !
                IF (TRIM(OCFLUX) == 'COLA')THEN
                   rhi  (i) = rhi_cola  (i)
                   rmi  (i) = rmi_cola  (i)
                   sens (i) = sens_cola (i)
                   evap (i) = evap_cola (i)
                ELSE IF(TRIM(OCFLUX) == 'WGFS')THEN
                   sens (i) = HFX_SEA(i) 
                   evap (i) = LH_SEA(i)
                   rhi  (i) = CHS_SEA(i)
                   rmi  (i) = WSPD_SEA(i)*(CM_SEA(i) )
                ELSE IF (TRIM(OCFLUX) == 'UKME')THEN
                   rhi  (i) = rhi_uk(i)
                   rmi  (i) = rmi_uk(i)
                   sens (i) = sens_uk(i)
                   evap (i) = evap_uk(i)
                ELSE
                   STOP 'SEA WATERFLUX'
                END IF
             END IF
             IF(atmpbl /= 1)THEN
                dtmdt= (ah(i)*sens(i))/(dtc3x*ah(i)*rho(i)*cp*rhi(i))
                dqmdt= (al(i)*evap(i))/(dtc3x*al(i)*rho(i)*hl*rhi(i))
             ELSE
                dtmdt= (ah(i)*sens(i))/(tmtx(i,1,2)+dtc3x*ah(i)*rho(i)*cp*rhi(i))
                dqmdt= (al(i)*evap(i))/(qmtx(i,1,2)+dtc3x*al(i)*rho(i)*hl*rhi(i))
             END IF
             dtm=dtmdt*dtc3x
             dqm=dqmdt*dtc3x
             tsfc   (i)=gt(i,1)+dtm
             qsfc   (i)=gq(i,1)+dqm
          ELSE
             IF(atmpbl == 1)THEN

                IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) < tice+0.01_r8) THEN
                   !
                   ! Solution of sea ice
                   !
                   IF (TRIM(OCFLUX) == 'COLA')THEN
                      IF (TRIM(ICEMODEL)=='COLA')THEN
                         rhi  (i) = rhi_cola(i)
                         rmi  (i) = rmi_cola(i)
                         sens (i)=  sens_cola (i)
                         evap (i)=  evap_cola (i)
                      ELSE IF (TRIM(ICEMODEL)=='SSIB')THEN
                         rhi  (i) = ICE_rhi(i)
                         rmi  (i) = ICE_rmi(i)
                         sens (i) = ICE_HFLUX(i)
                         evap (i) = ICE_ELATEN(i)
                      END IF
                   ELSE IF(TRIM(OCFLUX) == 'WGFS')THEN
                      IF (TRIM(ICEMODEL)=='COLA')THEN
                         IF(omlmodel)THEN
                            rhi  (i) = rhi_cola(i)
                            rmi  (i) = rmi_cola(i)
                            sens (i)=  sens_cola (i)
                            evap (i)=  evap_cola (i)
                         ELSE
                            rhi  (i) = rhi_cola(i)
                            rmi  (i) = rmi_cola(i)
                            sens (i)=  sens_cola (i)
                            evap (i)=  evap_cola (i)
                         END IF
                      ELSE IF (TRIM(ICEMODEL)=='SSIB')THEN
                         rhi  (i) = ICE_rhi(i)
                         rmi  (i) = ICE_rmi(i)
                         sens (i) = ICE_HFLUX(i)
                         evap (i) = ICE_ELATEN(i)
                      END IF
                   ELSE IF (TRIM(OCFLUX) == 'UKME')THEN
                      IF (TRIM(ICEMODEL)=='COLA')THEN
                         rhi  (i) = rhi_uk(i)
                         rmi  (i) = rmi_uk(i)
                         sens (i) = sens_uk(i)
                         evap (i) = evap_uk(i)
                      ELSE IF (TRIM(ICEMODEL)=='SSIB')THEN
                         rhi  (i) = ICE_rhi(i)
                         rmi  (i) = ICE_rmi(i)
                         sens (i) = ICE_HFLUX(i)
                         evap (i) = ICE_ELATEN(i)
                      END IF
                   ELSE
                      STOP 'SEA AND ICE FLUX'
                   END IF
                ELSE IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) > tice+0.01_r8) THEN
                   !
                   ! Solution of sea water
                   !
                   IF (TRIM(OCFLUX) == 'COLA')THEN
                      rhi  (i) = rhi_cola(i)
                      rmi  (i) = rmi_cola(i)
                      sens (i)=  sens_cola (i)
                      evap (i)=  evap_cola (i)
                   ELSE IF(TRIM(OCFLUX) == 'WGFS')THEN
                      sens (i) = HFX_SEA(i) 
                      evap (i) = LH_SEA(i)
                      rhi  (i) = CHS_SEA(i)
                      rmi  (i) = WSPD_SEA(i)*(CM_SEA(i) )
                   ELSE IF (TRIM(OCFLUX) == 'UKME')THEN
                      rhi  (i) = rhi_uk(i)
                      rmi  (i) = rmi_uk(i)
                      sens (i) = sens_uk(i)
                      evap (i) = evap_uk(i)
                   ELSE
                      STOP 'SEA WATERFLUX'
                   END IF
                END IF
                IF(atmpbl /= 1)THEN
                   dtmdt= (ah(i)*sens(i))/(dtc3x*ah(i)*rho(i)*cp*rhi(i))
                   dqmdt= (al(i)*evap(i))/(dtc3x*al(i)*rho(i)*hl*rhi(i))
                ELSE
                   dtmdt= (ah(i)*sens(i))/(tmtx(i,1,2)+dtc3x*ah(i)*rho(i)*cp*rhi(i))
                   dqmdt= (al(i)*evap(i))/(qmtx(i,1,2)+dtc3x*al(i)*rho(i)*hl*rhi(i))
                END IF
                dtm=dtmdt*dtc3x
                dqm=dqmdt*dtc3x
                tsfc   (i)=gt(i,1)+dtm
                qsfc   (i)=gq(i,1)+dqm
             ELSE
                IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) < tice+0.01_r8) THEN
                   !
                   ! Solution of sea ice
                   !
                   IF (TRIM(OCFLUX) == 'COLA')THEN
                      IF (TRIM(ICEMODEL)=='COLA')THEN
                         rhi  (i) = rhi_cola(i)
                         rmi  (i) = rmi_cola(i)
                         sens (i)=  sens_cola (i)
                         evap (i)=  evap_cola (i)
                      ELSE IF (TRIM(ICEMODEL)=='SSIB')THEN
                         rhi  (i) = ICE_rhi(i)
                         rmi  (i) = ICE_rmi(i)
                         sens (i) = ICE_HFLUX(i)
                         evap (i) = ICE_ELATEN(i)
                      END IF
                   ELSE IF(TRIM(OCFLUX) == 'WGFS')THEN
                      IF (TRIM(ICEMODEL)=='COLA')THEN
                         IF(omlmodel)THEN
                            rhi  (i) = rhi_cola(i)
                            rmi  (i) = rmi_cola(i)
                            sens (i)=  sens_cola (i)
                            evap (i)=  evap_cola (i)
                         ELSE
                            rhi  (i) = rhi_cola(i)
                            rmi  (i) = rmi_cola(i)
                            sens (i)=  sens_cola (i)
                            evap (i)=  evap_cola (i)
                         END IF
                      ELSE IF (TRIM(ICEMODEL)=='SSIB')THEN
                         IF(omlmodel)THEN
                            rhi  (i) = ICE_rhi(i)
                            rmi  (i) = ICE_rmi(i)
                            sens (i) = ICE_HFLUX(i)
                            evap (i) = ICE_ELATEN(i)
                         ELSE
                            rhi  (i) = ICE_rhi(i)
                            rmi  (i) = ICE_rmi(i)
                            sens (i) = ICE_HFLUX(i)
                            evap (i) = ICE_ELATEN(i)
                         END IF
                      END IF
                   ELSE IF (TRIM(OCFLUX) == 'UKME')THEN
                      IF (TRIM(ICEMODEL)=='COLA')THEN
                         rhi  (i) = rhi_uk(i)
                         rmi  (i) = rmi_uk(i)
                         sens (i) = sens_uk(i)
                         evap (i) = evap_uk(i)
                      ELSE IF (TRIM(ICEMODEL)=='SSIB')THEN
                         rhi  (i) = ICE_rhi(i)
                         rmi  (i) = ICE_rmi(i)
                         sens (i) = ICE_HFLUX(i)
                         evap (i) = ICE_ELATEN(i)
                      END IF
                   ELSE
                      STOP 'SEA AND ICE FLUX'
                   END IF

                ELSE IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) > tice+0.01_r8) THEN
                   !
                   ! Solution of sea water
                   !
                   IF (TRIM(OCFLUX) == 'COLA')THEN
                      rhi  (i) = rhi_cola(i)
                      rmi  (i) = rmi_cola(i)
                      sens (i) = sens_cola (i)
                      evap (i) = evap_cola (i)
                   ELSE IF(TRIM(OCFLUX) == 'WGFS')THEN
                      sens (i) = HFX_SEA(i) 
                      evap (i) = LH_SEA(i)
                      rhi  (i) = CHS_SEA(i)
                      rmi  (i) = WSPD_SEA(i)*(CM_SEA(i) )
                   ELSE IF (TRIM(OCFLUX) == 'UKME')THEN
                      rhi  (i) = rhi_uk(i)
                      rmi  (i) = rmi_uk(i)
                      sens (i) = sens_uk(i)
                      evap (i) = evap_uk(i)
                   ELSE
                      STOP 'SEA WATERFLUX'
                   END IF
                END IF
                IF(atmpbl /= 1)THEN
                   dtmdt= (ah(i)*sens(i))/(dtc3x*ah(i)*rho(i)*cp*rhi(i))
                   dqmdt= (al(i)*evap(i))/(dtc3x*al(i)*rho(i)*hl*rhi(i))
                ELSE
                   dtmdt= (ah(i)*sens(i))/(tmtx(i,1,2)+dtc3x*ah(i)*rho(i)*cp*rhi(i))
                   dqmdt= (al(i)*evap(i))/(qmtx(i,1,2)+dtc3x*al(i)*rho(i)*hl*rhi(i))
                END IF
                dtm=dtmdt*dtc3x
                dqm=dqmdt*dtc3x
                tsfc   (i)=gt(i,1)+dtm
                qsfc   (i)=gq(i,1)+dqm
             END IF
             !dtmdt=(ah(i)*sens(i))/(tmtx(i,1,2)+dtc3x*ah(i)*rho(i)*cp*rhi(i))
             !dqmdt=(al(i)*evap(i))/(qmtx(i,1,2)+dtc3x*al(i)*rho(i)*hl*rhi(i))
             !dtm=dtmdt*dtc3x
             !dqm=dqmdt*dtc3x
             !tsfc   (i)=gt(i,1)+dtm
             !qsfc   (i)=gq(i,1)+dqm
             gt  (i,1)=gtsav(i)
             gq  (i,1)=gqsav(i)
             rnet(i)=rnet(i)-stefan*tsurf(i)**4
             IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) < tice+0.01_r8) THEN
                cond(i)=gice*(tsurf(i)-tice)
                stor(i)=hscap*c0(i)
                rnet(i)=rnet(i)-stefan*tsurf(i)**3*4.0_r8 *c0(i)
                tsurf(i)=MIN(tsurf(i),tice)
                tsea (i)=-   tsurf(i)
             END IF
          END IF
       END IF
    END DO

    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN
          umom(i)=rho(i)*gu(i,1)*rmi(i)
          vmom(i)=rho(i)*gv(i,1)*rmi(i)
          am  (i)=gb100/gps(i)
          IF(atmpbl /= 1)THEN
             umtx(i,1,3)=(-am(i)*umom(i))/(dtc3x*am(i)*rho(i)*rmi(i))
             umtx(i,1,4)=(-am(i)*vmom(i))/(dtc3x*am(i)*rho(i)*rmi(i))
          ELSE
             umtx(i,1,3)=(umtx(i,1,3)-am(i)*umom(i))/(umtx(i,1,2)+dtc3x*am(i)*rho(i)*rmi(i))
             umtx(i,1,4)=(umtx(i,1,4)-am(i)*vmom(i))/(umtx(i,1,2)+dtc3x*am(i)*rho(i)*rmi(i))
          END IF
          !
          !     set surface stress use of pseudo winds to true winds
          !     for output diagnostics
          !
          umom(i)=umom(i)/sinclt(i)
          vmom(i)=vmom(i)/sinclt(i)
          Ustarm(i) = SQRT(umom(i)**2 + vmom(i)**2)
          IF(Ustarm(i)==0.0_r8)Ustarm(i)=0.007_r8
          um  (i)=gu (i,1)/sinclt(i)
          vm  (i)=gv (i,1)/sinclt(i)
          speedm(i)=SQRT(um(i)**2 + vm(i)**2)
          speedm(i)=MAX(2.0_r8 , speedm(i))
          !
          dzm   (i)=rbyg*gt(i,1)

          !THAT  =effective_surface_skin_temperature [K]
          !BSTAR = (grav/(rhos*sqrt(CM*max(UU,1.e-30)/RHOS))) *  &
          !(CT*(TH-TA-(MAPL_GRAV/MAPL_CP)*DZ)/TA + MAPL_VIREPS*CQ*(QH-QA))
          !
          !                            grav
          !BSTAR = -------------------------------------------------------------- *
          !         (rho(i)*sqrt(CU(i)*max(speed(i),1.e-30_r8)/RHOS(i)))
          !
          !
          !
          !
          !         (CT(i)*(tsfc(i)-gt(i)-(GRAV/CP)*dzm(i))/gt(i) + MAPL_VIREPS*CT(i)*(qsfc-gq(i)))
          !
          !
          !
          !bstar(i) = (grav/(rho(i)*sqrt(cu(i)*max(speedm(i),1.e-30_r8)/rho(i)))) + &
          !           (ct(i)*(tsfc(i)-gt(i)-(grav/cp)*dzm(i))/gt(i) + mapl_vireps*ct(i)*(qsfc(i)-gq(i)))

          bstar(i)=cu(i)*grav*(ct(i)*(tsfc(i)-gt(i,1)-(grav/cp)*dzm(i))/gt(i,1))! + mapl_vireps*ct(i)*(qsfc(i)-gq(i)))

          !PRINT*,' seasfc(i)=' ,bstar(i),(grav/(rho(i)*sqrt(cu(i)*max(speedm(i),1.e-30_r8)/rho(i)))),&
          !                     ct(i)*(tsfc(i)-gt(i,1)-(grav/cp)*dzm(i))/gt(i),mapl_vireps*ct(i)*(qsfc(i)-gq(i))

       END IF
    END DO

  END SUBROUTINE seasfc

END MODULE Sfc_SeaFlux_Interface
