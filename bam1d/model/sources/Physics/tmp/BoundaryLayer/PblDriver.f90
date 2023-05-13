!
!  $Author: pkubota $
!  $Date: 2008/04/09 12:41:17 $
!  $Revision: 1.1 $
!
MODULE PblDriver

  USE Constants, ONLY :     &
       hl,            &
       cp,            &
       grav,          &
       gasr,          &
       r8

  USE Options, ONLY :  &
       nferr, nfprt,atmpbl,microphys

  USE Pbl_MellorYamada0, ONLY :     &
       InitPbl_MellorYamada0,&
       MellorYamada0

  USE Pbl_MellorYamada1, ONLY : &
       InitPbl_MellorYamada1,&
       MellorYamada1  

  USE CCM3PBL, ONLY : &
       vdinti,&
       vdintr

  USE Diagnostics, ONLY: &
        updia,dodia , &
        StartStorDiag,&
        nDiag_vdheat, & ! vertical diffusion heating
        nDiag_vdmois, & ! vertical diffusion moistening
        nDiag_vduzon, & ! vertical diffusion zonal momentum change
        nDiag_vdvmer, & ! vertical diffusion meridional momentum change
        nDiag_pblstr, & ! surface friction velocity
        nDiag_hghpbl, & ! planetary boundary layer height
        nDiag_khdpbl, & ! diffusion coefficient for heat
        nDiag_kmdpbl, & ! diffusion coefficient for momentum
        nDiag_ricpbl, & ! bulk Richardson no. from level to ref lev
        nDiag_ObuLen, & ! 203! Obukhov length                           (m)
        nDiag_InPhiM, & ! 204! inverse phi function for momentum     ()
        nDiag_tkemyj, & ! Turbulent Kinetic Energy
        nDiag_InPhiH    ! 205! inverse phi function for heat           ()


  USE GridHistory, ONLY:       &
       IsGridHistoryOn, StoreGridHistory, dogrh, &
       nGHis_vdheat, nGHis_vduzon, nGHis_vdvmer, nGHis_vdmois

  USE Parallelism, ONLY: &
       MsgOne, FatalError

  IMPLICIT NONE
  LOGICAL  :: TFlux=.FALSE.
  INTEGER, PARAMETER  :: pnats_lsc=0  ! number of non-advected trace species
  INTEGER, PARAMETER  :: pcnst_lsc=1 ! number of constituents (including water vapor)
  INTEGER, PARAMETER  :: pnats_mic=2  ! number of non-advected trace species
  INTEGER, PARAMETER  :: pcnst_mic=3 ! number of constituents (including water vapor)

  PRIVATE
  PUBLIC :: InitPBLDriver
  PUBLIC :: Pbl_Driver

CONTAINS
  SUBROUTINE InitPBLDriver(kmax, sl, delsig, si)
    IMPLICIT NONE
    INTEGER      , INTENT(IN) :: kmax
    REAL(KIND=r8), INTENT(IN) :: sl     (kmax)
    REAL(KIND=r8), INTENT(IN) :: delsig (kmax)
    REAL(KIND=r8), INTENT(IN) :: si     (kmax+1)
    INTEGER  :: pnats  ! number of non-advected trace species
    INTEGER  :: pcnst ! number of constituents (including water vapor)

    IF     (atmpbl==1)THEN
       CALL InitPbl_MellorYamada0(kmax, sl, delsig, si(1:kMax+1))
    ELSE IF(atmpbl==2)THEN
       CALL InitPbl_MellorYamada1()
    ELSE IF(atmpbl==3)THEN
       IF (microphys) THEN
          pnats=pnats_mic
          pcnst=pcnst_mic
          CALL vdinti(kmax,sl, delsig, si(1:kMax+1),pnats,pcnst)
       ELSE
          pnats=pnats_lsc
          pcnst=pcnst_lsc
          CALL vdinti(kmax,sl, delsig, si(1:kMax+1),pnats,pcnst)
       END IF
    ELSE   
        CALL MsgOne('**(InitPBLDriver)**','its not set pbl parametrization')

       WRITE(*,*)'its not set pbl parametrization'
       STOP
    END IF
  END SUBROUTINE InitPBLDriver

  !----------------------------------------------------------------------
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !----------------------------------------------------------------------
  SUBROUTINE pbl_driver ( &
       gu       ,gv      ,gt      ,gq      ,delsig  ,ncols  , &
       kmax     ,delt    ,colrad  ,tmtx    ,qmtx    ,umtx   , &
       tmsfc    ,qmsfc   ,umsfc   ,gl0     ,PBL     ,sigki  , &
       gps      ,sl      ,si      ,z0      ,xland   , &
       lowlyr   ,SICE    ,SNOW    ,sens    ,evap    ,thz0   , &
       qz0      ,uz0     ,vz0     ,tkemyj  ,USTAR   ,AKHS   , &
       AKMS     ,latco   ,QSFC    ,TSK     ,topog   ,ct     , &
       LwCoolRate,LwCoolRateC,cldtot,qliq, bstar,&
       htdisp   ,taux    ,tauy    ,kt      ,jdt     , &
       PBL_CoefKm, PBL_CoefKh ,pblh    ,tpert,qpert, tstar ,wstar, gice, gliq )
    IMPLICIT NONE
    INTEGER      ,    INTENT(in   ) :: ncols
    INTEGER      ,    INTENT(in   ) :: kmax
    INTEGER      ,    INTENT(in   ) :: PBL
    INTEGER      ,    INTENT(in   ) :: latco
    REAL(KIND=r8),    INTENT(in   ) :: gu     (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: gv     (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: gt     (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: gq     (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: delsig (kmax)
    REAL(KIND=r8),    INTENT(in   ) :: delt
    REAL(KIND=r8),    INTENT(in   ) :: colrad (ncols)
    REAL(KIND=r8),    INTENT(inout) :: tmtx   (ncols,kmax,3)
    REAL(KIND=r8),    INTENT(inout) :: qmtx   (ncols,kmax,5)
    REAL(KIND=r8),    INTENT(inout) :: umtx   (ncols,kmax,4)
    REAL(KIND=r8),    INTENT(in   ) :: tmsfc  (ncols,kmax,3)
    REAL(KIND=r8),    INTENT(in   ) :: qmsfc  (ncols,kmax,5)
    REAL(KIND=r8),    INTENT(in   ) :: umsfc  (ncols,kmax,4)
    REAL(KIND=r8),    INTENT(inout) :: gl0    (ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: sigki  (kmax)
    REAL(KIND=r8),    INTENT(IN   ) :: gps    (ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: sl     (kmax)
    REAL(KIND=r8),    INTENT(IN   ) :: si     (kmax+1)
    REAL(KIND=r8),    INTENT(inout) :: z0     (ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: XLAND  (ncols)
    INTEGER      ,    INTENT(in   ) :: LOWLYR (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: SICE   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: SNOW   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: sens   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: evap   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: THZ0   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: QZ0    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: UZ0    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: VZ0    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: tkemyj (ncols,kMax+1)
    REAL(KIND=r8),    INTENT(inout) :: USTAR  (ncols)  
    REAL(KIND=r8),    INTENT(inout) :: AKHS   (ncols)  
    REAL(KIND=r8),    INTENT(inout) :: AKMS   (ncols)  
    REAL(KIND=r8),    INTENT(inout) :: QSFC   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: TSK    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: topog  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: CT     (ncols) 
    REAL(KIND=r8),    INTENT(in   ) :: LwCoolRate  (ncols,kMax) 
    REAL(KIND=r8),    INTENT(in   ) :: LwCoolRateC (ncols,kMax) 
    REAL(KIND=r8),    INTENT(in   ) :: cldtot      (ncols,kMax) 
    REAL(KIND=r8),    INTENT(in   ) :: qliq        (ncols,kMax) 
    REAL(KIND=r8),    INTENT(in   ) :: bstar (ncols) 
    REAL(KIND=r8),    INTENT(in   ) :: htdisp(ncols) 
    REAL(KIND=r8),    INTENT(in   ) :: taux (1:nCols)
    REAL(KIND=r8),    INTENT(in   ) :: tauy (1:nCols)
    INTEGER      ,    INTENT(in   ) :: kt
    INTEGER      ,    INTENT(in   ) :: jdt
    REAL(KIND=r8),    INTENT(INOUT) :: PBL_CoefKm(ncols, kmax)
    REAL(KIND=r8),    INTENT(INOUT) :: PBL_CoefKh(ncols, kmax)

    REAL(KIND=r8), INTENT(INOUT) :: pblh   (1:nCols)  
    REAL(KIND=r8), INTENT(INOUT) :: tpert  (1:nCols)  
    REAL(KIND=r8), INTENT(INOUT) :: qpert  (1:nCols)  
    REAL(KIND=r8), INTENT(OUT  ) :: tstar  (1:nCols)
    REAL(KIND=r8), INTENT(OUT  ) :: wstar  (1:nCols)
    REAL(KIND=r8),    OPTIONAL,  INTENT(inout) :: gice (ncols,kMax)
    REAL(KIND=r8),    OPTIONAL,  INTENT(inout) :: gliq (ncols,kMax)

    REAL(KIND=r8) :: EXCH_H      (ncols,kMax)
    REAL(KIND=r8) :: ELFLX       (ncols) 
    REAL(KIND=r8) :: EL_MYJ      (nCols,kMax)
    REAL(KIND=r8) :: tkemyj_local(1:nCols,kMax+1) 
    REAL(KIND=r8) :: obklen      (nCols)           ! Obukhov length
    REAL(KIND=r8) :: phiminv     (nCols)        
    REAL(KIND=r8) :: phihinv     (nCols)        
    INTEGER       :: KPBL        (nCols)   
    INTEGER       :: pcnst                  ! number of constituents (including water vapor)      

    REAL(KIND=r8) :: RUBLTEN  (nCols,kMax)
    REAL(KIND=r8) :: RVBLTEN  (nCols,kMax)
    REAL(KIND=r8) :: RTHBLTEN (nCols,kMax)
    REAL(KIND=r8) :: RQVBLTEN (nCols,kMax,3)
    REAL(KIND=r8) :: RQCBLTEN (nCols,kMax)
    REAL(KIND=r8) :: DT
    REAL(KIND=r8) :: r100,rbyg
    REAL(KIND=r8) :: bps     (ncols,kmax)
    REAL(KIND=r8) :: PMID    (ncols,kmax)
    REAL(KIND=r8) :: PINT    (ncols,kmax)
    REAL(KIND=r8) :: T       (ncols,kmax)
    REAL(KIND=r8) :: TH      (ncols,kmax)
    REAL(KIND=r8) :: QV      (ncols,kmax)
    REAL(KIND=r8) :: QC      (ncols,kmax)
    REAL(KIND=r8) :: U       (ncols,kmax)
    REAL(KIND=r8) :: V       (ncols,kmax) 
    REAL(KIND=r8) :: RHO     (ncols,kmax)
    REAL(KIND=r8) :: EXNER   (nCols,kMax)
    REAL(KIND=r8) :: psur    (1:nCols)
    REAL(KIND=r8) :: delz    (nCols,kMax)
    REAL(KIND=r8) :: press   (nCols,kMax)
    REAL(KIND=r8) :: tv      (nCols,kMax)
    REAL(KIND=r8) :: ze      (nCols,kMax)
    REAL(KIND=r8) :: terr    (1:nCols)
    REAL(KIND=r8) :: hfx     (1:nCols)
    REAL(KIND=r8) :: qfx     (1:nCols)
    real(KIND=r8) :: rpdeli  (nCols,kMax)      ! 1./(pmid(k+1)-pmid(k))
    real(KIND=r8) :: rpdel   (nCols,kMax)      ! 1./(pint(k+1)-pint(k))
    real(KIND=r8) :: qm1     (nCols,kMax,3)
    REAL(KIND=r8) :: CHKLOWQ (nCols)!-- CHKLOWQ - is either 0 or 1 (so far set equal to 1).
    REAL(KIND=r8) :: rino    (ncols,kmax)
    REAL(KIND=r8) :: zhalf   (nCols,kMax + 1)        
    REAL(KIND=r8) :: delzhalf(nCols,kMax + 1)        
    REAL(KIND=r8) :: zehalf  (nCols,kMax + 1)        
    REAL(KIND=r8) :: kvh     (nCols,kMax + 1)         ! diffusion coefficient for heat
    REAL(KIND=r8) :: kvm     (nCols,kMax + 1)         ! diffusion coefficient for momentum

    !-- used only in MYJPBL. 
    REAL(KIND=r8) :: up1     (1:nCols,kMax)
    REAL(KIND=r8) :: vp1     (1:nCols,kMax)
    REAL(KIND=r8) :: tmv1    (1:nCols,kMax)
    REAL(KIND=r8) :: qpert_local(1:nCols,3)
    REAL(KIND=r8) :: cflx       (1:nCols,3)
    REAL(KIND=r8) :: pintm1     (1:nCols,kMax+1)
    REAL(KIND=r8) :: pstarln    (1:nCols)
    REAL(KIND=r8) :: zm         (1:nCols,kMax)
    REAL(KIND=r8) :: pmidm1     (1:nCols,kMax)
    REAL(KIND=r8) :: FRLAND     (1:nCols)!land_fraction (#)

    ! For RUC LSM CHKLOWQ needed for MYJPBL should 
    ! 1 because is actual specific humidity at the surface, and
    ! not the saturation value
    REAL    :: bldt 
    INTEGER ::i,k,l
    INTEGER :: STEPBL! bldt (max_dom)= 0,; minutes between boundary-layer physics calls
    !-- calculate pbl time step
    !   STEPBL = nint(BLDT*60./DT)
    !   STEPBL = max(STEPBL,1)
    kvh=0.0_r8
    kvm=0.0_r8
    rino=0.0_r8
    obklen  =0.0_r8
    phiminv =0.0_r8
    phihinv =0.0_r8
    KPBL    =1
    qpert_local=0.0_r8
    tkemyj_local=0.0_r8
    DO k=1,kMax+1
       DO i=1,nCols
          tkemyj_local(i,k)=tkemyj    (i,k)
       END DO
    END DO

    IF( atmpbl == 1) THEN
    IF (microphys) THEN
       pcnst=3
       CALL MellorYamada0 ( &
            gu    (1:ncols,1:kmax)    ,gv    (1:ncols,1:kmax   ) ,gt         (1:ncols,1:kmax)      ,gq           (1:ncols,1:kmax)      ,gps   (1:ncols)              ,&
            delsig(1:kmax)            ,ncols                     ,kmax                             ,delt                               ,colrad(1:ncols)              ,&
            tmtx  (1:ncols,1:kmax,1:3),qmtx  (1:ncols,1:kmax,1:5),umtx       (1:ncols,1:kmax,1:4)  ,tmsfc        (1:ncols,1:kmax,1:3)  ,qmsfc (1:ncols,1:kmax,1:5)   ,&
            umsfc (1:ncols,1:kmax,1:4),gl0   (1:ncols)           ,PBL_CoefKm (1:ncols,1:kmax)      ,PBL_CoefKh   (1:ncols,1:kmax)      ,taux  (1:nCols)              ,&
            tauy  (1:nCols)           ,sens  (1:nCols)           ,evap       (1:nCols)             ,pblh         (1:nCols)             ,QSFC  (1:nCols)              ,&
            TSK   (1:nCols)           ,tpert (1:nCols)           ,qpert_local(1:nCols,1:1)         ,tkemyj_local (1:nCols,1:kMax+1)    ,tstar (1:nCols)              ,&
            wstar (1:nCols)           ,gice  (1:ncols,1:kmax)    ,gliq       (1:ncols,1:kmax))
            qpert(1:nCols)=qpert_local(1:nCols,1)
    ELSE
       CALL MellorYamada0 ( &
            gu    (1:ncols,1:kmax)    ,gv    (1:ncols,1:kmax   ) ,gt         (1:ncols,1:kmax)      ,gq           (1:ncols,1:kmax)      ,gps   (1:ncols)              ,&
            delsig(1:kmax)            ,ncols                     ,kmax                             ,delt                               ,colrad(1:ncols)              ,&
            tmtx  (1:ncols,1:kmax,1:3),qmtx  (1:ncols,1:kmax,1:5),umtx       (1:ncols,1:kmax,1:4)  ,tmsfc        (1:ncols,1:kmax,1:3)  ,qmsfc (1:ncols,1:kmax,1:5)   ,&
            umsfc (1:ncols,1:kmax,1:4),gl0   (1:ncols)           ,PBL_CoefKm (1:ncols,1:kmax)      ,PBL_CoefKh   (1:ncols,1:kmax)      ,taux  (1:nCols)              ,&
            tauy  (1:nCols)           ,sens  (1:nCols)           ,evap       (1:nCols)             ,pblh         (1:nCols)             ,QSFC  (1:nCols)              ,&
            TSK   (1:nCols)           ,tpert (1:nCols)           ,qpert_local(1:nCols,1:1)         ,tkemyj_local (1:nCols,1:kMax+1)    ,tstar (1:nCols)              ,&
            wstar (1:nCols)           )
            qpert(1:nCols)=qpert_local(1:nCols,1)
    END IF
      
       DO k=1,kMax
          DO i=1,nCols
             tkemyj    (i,k)= tkemyj_local(i,k)
          END DO
       END DO

    ELSE IF( atmpbl == 2) THEN
       DT=2.0_r8*delt
       bldt=0.0_r8
       STEPBL = NINT(BLDT*60.0_r8/DT)
       STEPBL = 1!MAX(STEPBL,1)     
       r100=100.0e0_r8 /gasr
       DO k=1,kMAx
          DO i=1,nCols    
             !
             ! Factor conversion to potention temperature
             ! rk=gasr/cp
             !              1.0e0_r8
             !sigki (k)=---------------------
             !            sl(k) ** rk
             !
             bps  (i,k)=sigki(k) 
             PMID (i,k)=100.0_r8*gps(i)*sl(k) 
             PINT (i,k)=100.0_r8*gps(i)*si(k)
             T    (i,k)=gt (i,k)
             TH   (i,k)=T  (i,k)*bps(i,k)
             QV   (i,k)=gq (i,k)
             QC   (i,k)=0.0000_r8*gq (i,k)
             U    (i,k)=gu (i,k)/SIN( colrad(i))
             V    (i,k)=gv (i,k)/SIN( colrad(i))
             !                    (m*m)
             ! dzm   (i) = --------- = m
             !                      m
             rbyg=gasr/grav*delsig(k)*0.5e0_r8
             !DZ  (i,k)=rbyg * T(i,k)
             RHO (i,k)=r100*(gps(i)*sl(k))/T(i,k)
             ! pi     : psur function (PI) = cp(pres(z)/pres(super))**(R/cp) (x,y,z)
             !USTAR(i)=MIN(USTAR(i),10.9_r8)
             EXNER (i,k)=cp*((PINT(i,k)/PINT (i,1))**(gasr/cp))
          END DO
       END DO
    DO i=1,nCols
       psur(i)   =gps(i)*100.0_r8
       terr(i)   =MAX(topog(i),0.0_r8)   
    END DO
    DO k=1,kMax
       DO i=1,nCols
          press(i,k)=gps(i)*sl(k)
          press(i,k)=press(i,k)*100.0_r8
          tv(i,k)=gt(i,k)*(1.0_r8+0.608_r8*gq(i,k))
       END DO
    END DO    
    ze=0.0_r8
    DO i=1,nCols
      rbyg=gasr/grav*delsig(1)*0.5e0_r8
      IF(XLAND(i) >1.0_r8)THEN
         delz(i,1)=MAX((rbyg * tv(i,1) - htdisp(i)),0.5_r8)*0.75_r8
      ELSE
         delz(i,1)=MAX((rbyg * tv(i,1) - htdisp(i)),0.5_r8)      
      END IF
      ze  (i,1)= delz(i,1)!gt(i,1)*(gasr/grav)*(psur(i)-press(i,1))/psur(i)
      ze  (i,1)= terr(i  )+ze(i,1)
    END DO
    
    DO k=2,kMax
       DO i=1,nCols
          delz(i,k)=0.5_r8*gasr*(tv(i,k-1)+tv(i,k))* &
               LOG(press(i,k-1)/press(i,k))/grav
          ze(i,k)=ze(i,k-1)+ delz(i,k)
       END DO
    END DO

       DO i=1,nCols
          z0    (i)=MAX(z0(i),0.1e-7_r8)
          CHKLOWQ(i)=1
          ELFLX  (i) =  evap(i)
          hfx    (i) =  sens(i)
          qfx    (i) = (evap(i)/28.9_r8)/86400.0_r8

       ENDDO

       CALL MellorYamada1(nCols, & !INTENT(IN   )
            kMax    , & !INTENT(IN   ) 
            DT      , & !INTENT(IN   ) :: DT! time step (second)
            STEPBL  , & !INTENT(IN   ) :: STEPBL! bldt (max_dom)= 0,; minutes between boundary-layer physics calls
            ze(:,1) , & !INTENT(IN   ) :: HT   (1:nCols)!"HGT" "Terrain Height"   "m"
            delz    , & !INTENT(IN   ) :: DZ (1:nCols,1:kMax)!-- dz8w dz between full levels (m)
            PMID    , & !INTENT(IN   ) :: PMID(1:nCols,1:kMax)!-- p_phy pressure (Pa)
            PINT    , & !INTENT(IN   ) :: PINT  (1:nCols,1:kMax)!-- p8w pressure at full levels (Pa)
            TH      , & !INTENT(IN   ) :: TH(1:nCols,1:kMax)! potential temperature (K)
            T       , & !INTENT(IN   ) :: T(1:nCols,1:kMax) !t_phy         temperature (K)
            QV      , & !INTENT(IN   ) :: QV(1:nCols,1:kMax)! Qv         water vapor mixing ratio (kg/kg)
            QC      , & !INTENT(IN   ) :: CWM(1:nCols,1:kMax)! cloud water mixing ratio (kg/kg)
            U       , & !INTENT(IN   ) :: U(1:nCols,1:kMax)! u-velocity interpolated to theta points (m/s)
            V       , & !INTENT(IN   ) :: V(1:nCols,1:kMax)! v-velocity interpolated to theta points (m/s)
            TSK     , & !INTENT(IN   ) :: TSK    (1:nCols)!-- TSK surface temperature (K)
            CHKLOWQ , & !INTENT(IN   ) :: CHKLOWQ(1:nCols)!-- CHKLOWQ - is either 0 or 1 (so far set equal to 1).
            LOWLYR  , & !INTENT(IN   ) :: LOWLYR (1:nCols)!-- lowlyr index of lowest model layer above ground
            XLAND   , & !INTENT(IN   ) :: XLAND  (1:nCols)!-- XLAND         land mask (1 for land, 2 for water)
            SICE    , & !INTENT(IN   ) :: SICE   (1:nCols) !-- SICE liquid water-equivalent ice  depth (m)
            SNOW    , & !INTENT(IN   ) :: SNOW   (1:nCols) !-- SNOW liquid water-equivalent snow depth (m)
            ELFLX   , & !INTENT(IN   ) :: ELFLX  (1:nCols)!-- ELFLX--LH net upward latent heat flux at surface (W/m^2)
            bps     , & 
            colrad  , & 
            QSFC    , & !INTENT(INOUT) :: QSFC   (1:nCols)!-- qsfc specific humidity at lower boundary (kg/kg)
            THZ0    , & !INTENT(INOUT) :: THZ0   (1:nCols)!-- thz0 potential temperature at roughness length (K)
            QZ0     , & !INTENT(INOUT) :: QZ0    (1:nCols)!-- QZ0 specific humidity at roughness length (kg/kg)
            UZ0     , & !INTENT(INOUT) :: UZ0    (1:nCols)!-- uz0 u wind component at roughness length (m/s)
            VZ0     , & !INTENT(INOUT) :: VZ0    (1:nCols)!-- vz0 v wind component at roughness length (m/s)
            tkemyj_local(1:nCols,1:kMax) , & !INTENT(INOUT) :: !-- kinetic energy from Mellor-Yamada-Janjic (MYJ) (m^2/s^2)
            EXCH_H  , & !INTENT(INOUT) :: EXCH_H (1:nCols,1:kMax)
            USTAR   , & !INTENT(INOUT) :: USTAR  (1:nCols)!-- UST           u* in similarity theory (m/s)
            CT      , & !INTENT(INOUT) :: CT     (1:nCols)
            AKHS    , & !INTENT(INOUT) :: AKHS   (1:nCols)!-- akhs        sfc exchange coefficient of heat/moisture from MYJ
            AKMS    , & !INTENT(INOUT) :: AKMS   (1:nCols)!-- akms        sfc exchange coefficient of momentum from MYJ
            EL_MYJ  , & !INTENT(OUT  ) :: EL_MYJ   (1:nCols,1:kMax)! mixing length from Mellor-Yamada-Janjic (MYJ) (m)
            PBLH    , & !INTENT(OUT  ) :: PBLH     (1:nCols)       ! PBL height (m)
            KPBL    , & !INTENT(OUT  ) :: KPBL     (1:nCols)        !-- KPBL layer index of the PBL
            PBL_CoefKm,&
            PBL_CoefKh,&
            RUBLTEN , & !INTENT(OUT  ) :: RUBLTEN  (1:nCols,1:kMax)! U tendency due to PBL parameterization (m/s^2)
            RVBLTEN , & !INTENT(OUT  ) :: RVBLTEN  (1:nCols,1:kMax)! V tendency due to PBL parameterization (m/s^2)
            RTHBLTEN, & !INTENT(OUT  ) :: RTHBLTEN (1:nCols,1:kMax)! Theta tendency due to PBL parameterization (K/s)
            RQVBLTEN(1:nCols,1:kMax,1), & !INTENT(OUT  ) :: RQVBLTEN (1:nCols,1:kMax)! Qv tendency due to PBL parameterization (kg/kg/s)
            RQCBLTEN  ) !INTENT(OUT  ) :: RQCBLTEN (1:nCols,1:kMax)! Qc tendency due to PBL parameterization (kg/kg/s)
            qpert(1:nCols)=qpert_local(1:nCols,1)

       IF(TFlux)THEN
          DO i=1,nCols             
             qmtx(i,1,3) = qmsfc  (i,1,3)         
             tmtx(i,1,3) = tmsfc  (i,1,3)
             umtx(i,1,3) = umsfc  (i,1,3)
             umtx(i,1,4) = umsfc  (i,1,4)    
          END DO       
       ELSE
          DO i=1,nCols             

             qmtx(i,1,3) = RQVBLTEN(i,1,1)
             tmtx(i,1,3) = RTHBLTEN(i,1)/bps(i,1)
             umtx(i,1,3) = RUBLTEN (i,1)
             umtx(i,1,4) = RVBLTEN (i,1)
          END DO
       END IF
       DO k=1,kMax
          DO i=1,nCols
             tkemyj    (i,k) =tkemyj_local(i,k)
          END DO
       END DO       
       DO k=2,kMAx
          DO i=1,nCols 
             qmtx(i,k,3) = RQVBLTEN(i,k,1)
             tmtx(i,k,3) = RTHBLTEN(i,k)/bps(i,k)
             umtx(i,k,3) = RUBLTEN (i,k)
             umtx(i,k,4) = RVBLTEN (i,k)
          END DO
       END DO
       
       DO k=1,kMax
          DO i=1,nCols
             tkemyj    (i,k)= tkemyj_local(i,k)
          END DO
       END DO
       
    ELSE IF( atmpbl == 3) THEN
       !
       !identifies land points
       !
       DO i=1,nCols
          IF(XLAND(i) == 1.0_r8)THEN
             !
             ! LAND
             !
             FRLAND(i)= 1.0_r8
          ELSE
             !
             ! SEA
             !
             FRLAND(i)= 0.0_r8
          END IF
       END DO

       DT=2.0_r8 * delt
       r100=100.0e0_r8 /gasr
       l=kMax+1
       DO k=1,kMAx
          l=l-1
          DO i=1,nCols    
             !
             ! Factor conversion to potention temperature
             ! rk=gasr/cp
             !              1.0e0_r8
             !sigki (k)=---------------------
             !            sl(k) ** rk
             !
             bps  (i,l)=sigki(k)           ! Exner Function
             !
             ! Invert the layers of model
             !
             T    (i,l  )=gt (i,k)
             TH   (i,l  )=gt (i,k)*sigki(k) 
             qm1  (i,l,1)=gq (i,k)
             U    (i,l  )=gu (i,k)/SIN( colrad(i))
             V    (i,l  )=gv (i,k)/SIN( colrad(i))
          END DO
       END DO
       l=0
       DO k=kMAx,1,-1
          l=l+1
          DO i=1,nCols    
             pmidm1 (i,k) = 100.0_r8*gps(i)*sl(l) 
             tmv1   (i,k) = gt(i,l)*(1.0_r8 + 0.608_r8*gq(i,l))
          END DO
       END DO
       l=0
       DO k=kMAx+1,1,-1
          l=l+1
          DO i=1,nCols    
             pintm1 (i,k)=100.0_r8*gps(i)*si(l)
          END DO
       END DO
       DO i=1,nCols
          psur   (i) = gps(i)*100.0_r8
          pstarln(i) = log(psur(i))
          terr   (i) = MAX(topog(i),0.0_r8)   
       END DO
       DO k=1,kMax
          DO i=1,nCols
             press(i,k)=gps(i)*sl(k)
             press(i,k)=press(i,k)*100.0_r8
             tv(i,k)=gt(i,k)*(1.0_r8+0.608_r8*gq(i,k))
          END DO
       END DO    
       DO k=1,kMax
          DO i=1,nCols
             rpdel(i,k) = 1.0_r8/(pintm1(i,k+1) - pintm1(i,k))
          END DO
       END DO    

       do k=1,kMax-1
          do i=1,nCols
             rpdeli(i,k) = 1.0_r8/(pmidm1(i,k+1) - pmidm1(i,k))
          end do
       end do    
       !
       !  Convert j/kg to kg/mm3
       !

       DO i=1,nCols
         ELFLX  (i) =  evap(i)
         hfx    (i) =  sens(i)
         !
         !  J         N *m         Kg * m * m          m*m
         !------ = --------- =------------------- = ---------= hvap = 2.5104e+6! latent heat of vaporization of water (J kg-1)
         !  kg         kg          s*s Kg              s*s
         !
         ! (kg/m2/s)
         !   W        J         N *m       Kg * m * m       s*s         kg 
         ! ------ = ------- = --------- = ------------  * -------  = -------
         !  m*m      m*m*s      m*m*s       s*s*s*m*m       m*m       M*m*s
         !
         !           W         kg   
         !   =   ------ *  ------  =  LFlux / hvap
         !          m*m        J  
         qfx    (i)   = evap(i)/hl
         cflx   (i,1) = qfx(i)
       ENDDO
       rbyg=gasr/grav*delsig(1)*0.5e0_r8
       !
       !  Calculate the distance between the surface and the first layer of the model
       !
       ze=0.0_r8
       DO i=1,nCols
          terr(i)   =MAX(topog(i),0.0_r8)   
          IF(XLAND(i) >1.0_r8)THEN
             delz    (i,1)=MAX((rbyg * tv(i,1) - htdisp(i)),0.5_r8)*0.75_r8
             delzhalf(i,1)=MAX((rbyg * TSK(i)  - htdisp(i)),0.5_r8)*0.75_r8
          ELSE
             delz    (i,1)=MAX((rbyg * tv(i,1) - htdisp(i)),0.5_r8)      
             delzhalf(i,1)=MAX((rbyg * TSK(i)  - htdisp(i)),0.5_r8)     
          END IF
          ze    (i,1)= delz(i,1)
          ze    (i,1)=MAX(0.0_r8,ze(i,1))
          zehalf(i,1)=delzhalf(i,1)
          zehalf(i,1)=MAX(0.0_r8,zehalf(i,1))
       END DO
       DO k=2,kMax
         DO i=1,nCols
            delzhalf(i,k)=0.5_r8*gasr*(tv(i,k-1)+tv(i,k))* &
                 LOG(pintm1(i,kMax+3-k)/pintm1(i,kMax+2-k))/grav
            zehalf(i,k)=zehalf(i,k-1)+ delzhalf(i,k)
         END DO
       END DO
       DO i=1,nCols
          delzhalf(i,kMax+1)=gasr*(tv(i,kMax))* &
               LOG(pintm1(i,3)/pintm1(i,2))/grav
          zehalf(i,kMax+1)=zehalf(i,kMax)+ delzhalf(i,kMax+1)
       END DO
       DO k=1,kMax+1
          DO i=1,nCols
             zhalf(i,k) = zehalf(i,kMax+2-k)
          END DO
       END DO
       DO k=2,kMax
          DO i=1,nCols
             delz(i,k)=0.5_r8*gasr*(tv(i,k-1)+tv(i,k))* &
                  LOG(press(i,k-1)/press(i,k))/grav
             ze(i,k)=ze(i,k-1)+ delz(i,k)
          END DO
       END DO
       DO k=1,kMax
          DO i=1,nCols
              obklen (i) = 0.0_r8
              phiminv(i) = 0.0_r8
              phihinv(i) = 0.0_r8
              zm     (i,k) = ze(i,kMax+1-k)
             rino    (i,k) = 0.0_r8
          END DO
       END DO
       DO k=1,kMax+1
          DO i=1,nCols
             tkemyj_local(i,kMax+2-k)=tkemyj(i,k)
          END DO
       END DO
       IF (microphys) THEN
          pcnst=3
          DO i=1,nCols
                cflx(i,2)=0.0_r8 !surface constituent ice flux (kg/m2/s)
                cflx(i,3)=0.0_r8 !surface constituent liquid flux (kg/m2/s)
          END DO
          l=kMax+1
          DO k=1,kMax
             l=l-1
             DO i=1,nCols
                qm1(i,l,2)=  gice (i,k)!ice
                qm1(i,l,3)=  gliq (i,k)!liq
             END DO
          END DO 
       CALL vdintr(&
                 nCols                         , &! INTENT(IN   ) plon ! number of longitudes
                 nCols                         , &! INTENT(IN   ) plond ! slt extended domain longitude
                 kMax                          , &! INTENT(IN   ) plev ! number of vertical levels
                 pcnst                         , &! INTENT(IN   ) pcnst ! number of constituents (including water vapor)
                 DT                            , &! INTENT(IN   ) ztodt ! 2 delta-t
                 colrad      (1:nCols)         , &! INTENT(IN   ) Cosino the colatitude [radian]
                 gl0         (1:nCols)         , &! INTENT(IN   ) Maximum mixing length l0 in blackerdar's formula [m]
                 bstar       (1:nCols)         , &! INTENT(IN   ) surface_bouyancy_scale m s-2
                 FRLAND      (1:nCols)         , &! INTENT(IN   ) Fraction Land [%]
                 z0          (1:nCols)         , &! INTENT(IN   ) Rougosiness [m]
                 LwCoolRate  (1:nCols,1:kMax)  , &! INTENT(IN   ) air_temperature_tendency_due_to_longwave [K s-1]
                 LwCoolRateC (1:nCols,1:kMax)  , &! INTENT(IN   ) clear sky_air_temperature_tendency_lw [K s-1]
                 cldtot      (1:nCols,1:kMax)  , &! INTENT(IN   ) cloud fraction [%]
                 qliq        (1:nCols,1:kMax)  , &! INTENT(IN   ) cloud fraction [kg/kg]
                 pmidm1      (1:nCols,1:kMax)  , &! INTENT(IN   ) pmidm1(plond,plev) ! midpoint pressures
                 pintm1      (1:nCols,1:kMax+1), &! INTENT(IN   ) pintm1(plond,plev + 1)    ! interface pressures
                 bps         (1:nCols,1:kMax)  , &! INTENT(IN   ) psomc(plond,plev)         ! (psm1/pmidm1)**cappa
                 TH          (1:nCols,1:kMax)  , &! INTENT(IN   ) thm(plond,plev) ! potential temperature midpoints
                 zm          (1:nCols,1:kMax)  , &! INTENT(IN   ) zm(plond,plev) ! midpoint geopotential height above sfc
                 zhalf       (1:nCols,1:kMax+1), &! INTENT(IN   ) zhalf(plond,plev) ! interface pressures geopotential height
                 psur        (1:nCols)         , &! INTENT(IN   ) surface pressure [mb]
                 USTAR       (1:nCols)         , &! INTENT(IN   ) scale velocity turbulent [m/s]
                 TSK         (1:nCols)         , &! INTENT(IN   ) surface temperature
                 QSFC        (1:nCols)         , &! INTENT(IN   ) surface specific temperature
                 rpdel       (1:nCols,1:kMax)  , &! INTENT(IN   ) rpdel(plond,plev)! 1./pdel (thickness between interfaces)
                 rpdeli      (1:nCols,1:kMax)  , &! INTENT(IN   ) rpdeli(plond,plev)! 1./pdeli (thickness between midpoints)
                 U           (1:nCols,1:kMax)  , &! INTENT(IN   ) um1(plond,plev)         ! u-wind input
                 V           (1:nCols,1:kMax)  , &! INTENT(IN   ) vm1(plond,plev)         ! v-wind input
                 T           (1:nCols,1:kMax)  , &! INTENT(IN   ) tm1(plond,plev)         ! temperature input
                 taux        (1:nCols)         , &! INTENT(IN   ) taux(plond)                   ! x surface stress (n)
                 tauy        (1:nCols)         , &! INTENT(IN   ) tauy(plond)                   ! y surface stress (n)
                 HFX         (1:nCols)         , &! INTENT(IN   ) shflx(plond)! surface sensible heat flux (w/m2)
                 cflx        (1:nCols,1:pcnst) , &! INTENT(IN   ) cflx(plond,pcnst)! surface constituent flux (kg/m2/s)
                 qm1         (1:nCols,1:kMax,1:pcnst), &! INTENT(INOUT) qm1(plond,plev,pcnst)  ! initial/final constituent field
                 RTHBLTEN    (1:nCols,1:kMax)  , &! INTENT(OUT  ) dtv(plond,plev)     ! temperature tendency (heating)
                 RQVBLTEN    (1:nCols,1:kMax,1:pcnst)  , &! INTENT(OUT  ) dqv(plond,plev,pcnst)
                 RUBLTEN     (1:nCols,1:kMax)  , &! INTENT(OUT  ) duv(plond,plev)     ! u-wind tendency
                 RVBLTEN     (1:nCols,1:kMax)  , &! INTENT(OUT  ) dvv(plond,plev)     ! v-wind tendency
                 up1         (1:nCols,1:kMax)  , &! INTENT(OUT  ) up1(plond,plev)     ! u-wind after vertical diffusion
                 vp1         (1:nCols,1:kMax)  , &! INTENT(OUT  ) vp1(plond,plev)     ! v-wind after vertical diffusion
                 pblh        (1:nCols)         , &! INTENT(OUT  ) pblh(plond)! planetary boundary layer height
                 rino        (1:nCols,1:kMax)  , &! INTENT(INOUT) bulk Richardson no. from level to ref lev
                 tpert       (1:nCols)         , &! INTENT(OUT  ) tpert(plond)! convective temperature excess
                 qpert_local (1:nCols,1:pcnst) , &! INTENT(OUT  ) qpert(plond,pcnst)! convective humidity and constituent excess
                 tkemyj_local(1:nCols,1:kMax+1), &! INTENT(INOUT) Turbulent kinetic Energy [m/s]^2
                 kvh         (1:nCols,1:kMax+1), &! INTENT(OUT  ) Heat Coeficient Difusivity 
                 kvm         (1:nCols,1:kMax+1), &! INTENT(OUT  ) Momentun Coeficient Difusivity 
                 obklen      (1:nCols)         , &! INTENT(OUT  ) Heat Coeficient Difusivity 
                 phiminv     (1:nCols)         , &! INTENT(OUT  ) Momentum Stability Function 
                 phihinv     (1:nCols)         , &! INTENT(OUT  ) Heat Stability Function 
                 tstar       (1:nCols)         , &
                 wstar       (1:nCols)         )
        ELSE
         pcnst=1
            CALL vdintr(&
                 nCols                         , &! INTENT(IN   ) plon ! number of longitudes
                 nCols                         , &! INTENT(IN   ) plond ! slt extended domain longitude
                 kMax                          , &! INTENT(IN   ) plev ! number of vertical levels
                 pcnst                         , &! INTENT(IN   ) pcnst ! number of constituents (including water vapor)
                 DT                            , &! INTENT(IN   ) ztodt ! 2 delta-t
                 colrad      (1:nCols)         , &! INTENT(IN   ) Cosino the colatitude [radian]
                 gl0         (1:nCols)         , &! INTENT(IN   ) Maximum mixing length l0 in blackerdar's formula [m]
                 bstar       (1:nCols)         , &! INTENT(IN   ) surface_bouyancy_scale m s-2
                 FRLAND      (1:nCols)         , &! INTENT(IN   ) Fraction Land [%]
                 z0          (1:nCols)         , &! INTENT(IN   ) Rougosiness [m]
                 LwCoolRate  (1:nCols,1:kMax)  , &! INTENT(IN   ) air_temperature_tendency_due_to_longwave [K s-1]
                 LwCoolRateC (1:nCols,1:kMax)  , &! INTENT(IN   ) clear sky_air_temperature_tendency_lw [K s-1]
                 cldtot      (1:nCols,1:kMax)  , &! INTENT(IN   ) cloud fraction [%]
                 qliq        (1:nCols,1:kMax)  , &! INTENT(IN   ) cloud fraction [kg/kg]
                 pmidm1      (1:nCols,1:kMax)  , &! INTENT(IN   ) pmidm1(plond,plev) ! midpoint pressures
                 pintm1      (1:nCols,1:kMax+1), &! INTENT(IN   ) pintm1(plond,plev + 1)    ! interface pressures
                 bps         (1:nCols,1:kMax)  , &! INTENT(IN   ) psomc(plond,plev)         ! (psm1/pmidm1)**cappa
                 TH          (1:nCols,1:kMax)  , &! INTENT(IN   ) thm(plond,plev) ! potential temperature midpoints
                 zm          (1:nCols,1:kMax)  , &! INTENT(IN   ) zm(plond,plev) ! midpoint geopotential height above sfc
                 zhalf       (1:nCols,1:kMax+1), &! INTENT(IN   ) zhalf(plond,plev) ! interface pressures geopotential height
                 psur        (1:nCols)         , &! INTENT(IN   ) surface pressure [mb]
                 USTAR       (1:nCols)         , &! INTENT(IN   ) scale velocity turbulent [m/s]
                 TSK         (1:nCols)         , &! INTENT(IN   ) surface temperature
                 QSFC        (1:nCols)         , &! INTENT(IN   ) surface specific temperature
                 rpdel       (1:nCols,1:kMax)  , &! INTENT(IN   ) rpdel(plond,plev)! 1./pdel (thickness between interfaces)
                 rpdeli      (1:nCols,1:kMax)  , &! INTENT(IN   ) rpdeli(plond,plev)! 1./pdeli (thickness between midpoints)
                 U           (1:nCols,1:kMax)  , &! INTENT(IN   ) um1(plond,plev)         ! u-wind input
                 V           (1:nCols,1:kMax)  , &! INTENT(IN   ) vm1(plond,plev)         ! v-wind input
                 T           (1:nCols,1:kMax)  , &! INTENT(IN   ) tm1(plond,plev)         ! temperature input
                 taux        (1:nCols)         , &! INTENT(IN   ) taux(plond)                   ! x surface stress (n)
                 tauy        (1:nCols)         , &! INTENT(IN   ) tauy(plond)                   ! y surface stress (n)
                 HFX         (1:nCols)         , &! INTENT(IN   ) shflx(plond)! surface sensible heat flux (w/m2)
                 cflx        (1:nCols,1:pcnst) , &! INTENT(IN   ) cflx(plond,pcnst)! surface constituent flux (kg/m2/s)
                 qm1         (1:nCols,1:kMax,1:pcnst), &! INTENT(INOUT) qm1(plond,plev,pcnst)  ! initial/final constituent field
                 RTHBLTEN    (1:nCols,1:kMax)  , &! INTENT(OUT  ) dtv(plond,plev)     ! temperature tendency (heating)
                 RQVBLTEN    (1:nCols,1:kMax,1:pcnst) , &! INTENT(OUT  ) dqv(plond,plev,pcnst)
                 RUBLTEN     (1:nCols,1:kMax)  , &! INTENT(OUT  ) duv(plond,plev)     ! u-wind tendency
                 RVBLTEN     (1:nCols,1:kMax)  , &! INTENT(OUT  ) dvv(plond,plev)     ! v-wind tendency
                 up1         (1:nCols,1:kMax)  , &! INTENT(OUT  ) up1(plond,plev)     ! u-wind after vertical diffusion
                 vp1         (1:nCols,1:kMax)  , &! INTENT(OUT  ) vp1(plond,plev)     ! v-wind after vertical diffusion
                 pblh        (1:nCols)         , &! INTENT(OUT  ) pblh(plond)! planetary boundary layer height
                 rino        (1:nCols,1:kMax)  , &! INTENT(INOUT) bulk Richardson no. from level to ref lev
                 tpert       (1:nCols)         , &! INTENT(OUT  ) tpert(plond)! convective temperature excess
                 qpert_local (1:nCols,1:pcnst)       , &! INTENT(OUT  ) qpert(plond,pcnst)! convective humidity and constituent excess
                 tkemyj_local(1:nCols,1:kMax+1), &! INTENT(INOUT) Turbulent kinetic Energy [m/s]^2
                 kvh         (1:nCols,1:kMax+1), &! INTENT(OUT  ) Heat Coeficient Difusivity 
                 kvm         (1:nCols,1:kMax+1), &! INTENT(OUT  ) Momentun Coeficient Difusivity 
                 obklen      (1:nCols)         , &! INTENT(OUT  ) Heat Coeficient Difusivity 
                 phiminv     (1:nCols)         , &! INTENT(OUT  ) Momentum Stability Function 
                 phihinv     (1:nCols)         , &! INTENT(OUT  ) Heat Stability Function 
                 tstar       (1:nCols)         , &
                 wstar       (1:nCols)         )


        END IF

       IF (microphys) THEN!8121-8397
          l=kMax+1
          DO k=1,kMAx
             l=l-1
             DO i=1,nCols 
                qmtx(i,l,3) = RQVBLTEN(i,k,1)
                qmtx(i,l,4) = RQVBLTEN(i,k,2)  !ice
                qmtx(i,l,5) = RQVBLTEN(i,k,3)  !liq
                tmtx(i,l,3) = RTHBLTEN(i,k)/bps(i,k)
                umtx(i,l,3) = RUBLTEN (i,k)
                umtx(i,l,4) = RVBLTEN (i,k)
             END DO
          END DO
       ELSE
          l=kMax+1
          DO k=1,kMAx
             l=l-1
             DO i=1,nCols 
                qmtx(i,l,3) = RQVBLTEN(i,k,1)
                tmtx(i,l,3) = RTHBLTEN(i,k)/bps(i,k)
                umtx(i,l,3) = RUBLTEN (i,k)
                umtx(i,l,4) = RVBLTEN (i,k)
             END DO
          END DO
       END IF

       DO k=1,kMax
          DO i=1,nCols
             PBL_CoefKm(i,k)= kvm (i,kMax+1-k)
             PBL_CoefKh(i,k)= kvh (i,kMax+1-k)
          END DO
       END DO
       DO k=1,kMax+1
          DO i=1,nCols
             tkemyj    (i,k)=tkemyj_local(i,kMax+2-k)
          END DO
       END DO

       qpert(1:nCols)=qpert_local(1:nCols,1)

    ELSE
       WRITE(*,*)'its not set pbl parametrization'
       STOP

    END IF
    

    IF(PBL /= 1)THEN
       !-----------------
       ! Storage Diagnostic Fields
       !------------------
       IF( StartStorDiag)THEN
          CALL PblDiagnStorage(kt,jdt,latco,nCols,kMax,tmtx,qmtx,umtx,USTAR,pblh,&
            kvh,kvm  , rino ,obklen ,phiminv,phihinv ,tkemyj)
       END IF 
       !-----------------
       ! Storage GridHistory Fields
       !------------------
       IF(IsGridHistoryOn())THEN
          CALL PblGridHistoryStorage(kt,jdt,latco,nCols,kMax,tmtx,qmtx,umtx)
       END IF
    END IF
    
  END SUBROUTINE pbl_driver  

  SUBROUTINE PblGridHistoryStorage(kt,jdt,latco,nCols,kMax,tmtx,qmtx,umtx)
   IMPLICIT NONE
   INTEGER      , INTENT(IN   ) :: kt
   INTEGER      , INTENT(IN   ) :: jdt
   INTEGER      , INTENT(IN   ) :: latco
   INTEGER      , INTENT(IN   ) :: nCols
   INTEGER      , INTENT(IN   ) :: kMax
   REAL(KIND=r8), INTENT(in   ) :: tmtx (ncols,kmax,3)
   REAL(KIND=r8), INTENT(in   ) :: qmtx (ncols,kmax,5)
   REAL(KIND=r8), INTENT(in   ) :: umtx (ncols,kmax,4)
   IF( (kt.NE.0) .OR. (jdt.NE.1) ) THEN
      IF(dogrh(nGHis_vdheat,latco))CALL StoreGridHistory(tmtx(1:nCols,1:kMax,3),nGHis_vdheat,latco)
      IF(dogrh(nGHis_vduzon,latco))CALL StoreGridHistory(umtx(1:nCols,1:kMax,3),nGHis_vduzon,latco)
      IF(dogrh(nGHis_vdvmer,latco))CALL StoreGridHistory(umtx(1:nCols,1:kMax,4),nGHis_vdvmer,latco)
      IF(dogrh(nGHis_vdmois,latco))CALL StoreGridHistory(qmtx(1:nCols,1:kMax,3),nGHis_vdmois,latco)
   END IF
  
  END SUBROUTINE PblGridHistoryStorage
  
  
  SUBROUTINE PblDiagnStorage(kt,jdt,latco,nCols,kMax,tmtx,qmtx,umtx,USTAR,pblh,&
                             kvh,kvm  , rino,obklen,phiminv,phihinv,tkemyj)
   IMPLICIT NONE
   INTEGER      , INTENT(IN   ) :: kt
   INTEGER      , INTENT(IN   ) :: jdt
   INTEGER      , INTENT(IN   ) :: latco
   INTEGER      , INTENT(IN   ) :: nCols
   INTEGER      , INTENT(IN   ) :: kMax
   REAL(KIND=r8), INTENT(in   ) :: tmtx (ncols,kmax,3)
   REAL(KIND=r8), INTENT(in   ) :: qmtx (ncols,kmax,5)
   REAL(KIND=r8), INTENT(in   ) :: umtx (ncols,kmax,4)
   REAL(KIND=r8), INTENT(in   ) :: USTAR(ncols)
   REAL(KIND=r8), INTENT(in   ) :: pblh (ncols)
   REAL(KIND=r8), INTENT(in   ) :: kvh  (ncols,kmax+1)
   REAL(KIND=r8), INTENT(in   ) :: kvm  (ncols,kmax+1)
   REAL(KIND=r8), INTENT(in   ) :: rino (ncols,kmax)
   REAL(KIND=r8), INTENT(in   ) :: obklen (ncols)
   REAL(KIND=r8), INTENT(in   ) :: phiminv(ncols)
   REAL(KIND=r8), INTENT(in   ) :: phihinv(ncols)
   REAL(KIND=r8), INTENT(in   ) :: tkemyj(ncols,kMax+1)
   
   REAL(KIND=r8) :: brf(ncols,kmax)
   INTEGER        :: i,k
   IF( (kt.NE.0) .OR. (jdt.NE.1) ) THEN

      IF(dodia(nDiag_vdheat)) CALL updia(tmtx(1:nCols,1:kMax,3),nDiag_vdheat,latco)
   
      IF(dodia(nDiag_vduzon)) CALL updia(umtx(1:nCols,1:kMax,3),nDiag_vduzon,latco)
   
      IF(dodia(nDiag_vdvmer)) CALL updia(umtx(1:nCols,1:kMax,4),nDiag_vdvmer,latco)

      IF(dodia(nDiag_vdmois)) CALL updia(qmtx(1:nCols,1:kMax,3),nDiag_vdmois,latco)

      IF(dodia(nDiag_pblstr)) CALL updia(USTAR(1:nCols),nDiag_pblstr,latco)
     
      IF(dodia(nDiag_hghpbl)) CALL updia(pblh(1:nCols),nDiag_hghpbl,latco)
      
      IF(dodia(nDiag_ObuLen)) CALL updia(obklen (1:nCols),nDiag_ObuLen,latco)
      
      IF(dodia(nDiag_InPhiM)) CALL updia(phiminv(1:nCols),nDiag_InPhiM,latco)
      
      IF(dodia(nDiag_InPhiH)) CALL updia(phihinv(1:nCols),nDiag_InPhiH,latco)
      
      IF(dodia(nDiag_tkemyj)) CALL updia(tkemyj(1:nCols,1:kMax),nDiag_tkemyj,latco)

      IF(atmpbl == 3) THEN
         
         IF(dodia(nDiag_khdpbl)) THEN
            DO k=1,kMax
               DO i=1,nCols
                  brf(i,k) = kvh (i,kMax+1-k)
               END DO
            END DO
            CALL updia(brf (1:nCols,1:kMax),nDiag_khdpbl,latco)
         END IF 
         IF(dodia(nDiag_kmdpbl)) THEN
            DO k=1,kMax
               DO i=1,nCols
                  brf(i,k) = kvm (i,kMax+1-k)
               END DO
            END DO
            CALL updia(brf (1:nCols,1:kMax),nDiag_kmdpbl,latco)
         END IF
         IF(dodia(nDiag_ricpbl)) THEN
            DO k=1,kMax
               DO i=1,nCols
                  brf(i,k) = rino (i,kMax+1-k)
               END DO
            END DO
            CALL updia(brf(1:nCols,1:kMax),nDiag_ricpbl,latco)
         END IF   
      ELSE
         IF(dodia(nDiag_khdpbl))THEN
            CALL updia(kvh(1:nCols,1:kMax),nDiag_khdpbl,latco)
         END IF
         IF(dodia(nDiag_kmdpbl)) THEN
           CALL updia(kvm(1:nCols,1:kMax),nDiag_kmdpbl,latco)
         END IF
         IF(dodia(nDiag_ricpbl)) THEN
           CALL updia(rino(1:nCols,1:kMax),nDiag_ricpbl,latco)
         END IF 
      END IF
   
   END IF
  
  END SUBROUTINE PblDiagnStorage
  !----------------------------------------------------------------------
END MODULE PblDriver
