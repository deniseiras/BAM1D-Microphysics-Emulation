!
!  $Author: pkubota $
!  $Date: 2007/03/23 20:23:38 $
!  $Revision: 1.11 $
!
!
!MCGACPTEC : MEDIATION_LAYER:PHYSICS
!
MODULE SfcPBLDriver
  USE Constants, ONLY :     &
       cp,            &
       grav,          &
       gasr,          &
       r8

  USE Options, ONLY :     &
       sfcpbl,microphys,nClass

  USE Sfc_MellorYamada0, ONLY: &
       InitSfc_MellorYamada0, &
       SfcPbl_MY0

  USE Sfc_MellorYamada1, ONLY: &
       InitSfc_MellorYamada1, &
       SfcPbl_MYJ1

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: InitSfcPBL_Driver
  PUBLIC  :: SfcPBL_Driver

CONTAINS

  SUBROUTINE InitSfcPBL_Driver(kmax, sig, delsig, sigml,&
       RESTART,ibMax,jbMax,USTAR,LOWLYR )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: kmax
    REAL(KIND=r8),    INTENT(IN) :: sig   (kmax)
    REAL(KIND=r8),    INTENT(IN) :: delsig(kmax)
    REAL(KIND=r8),    INTENT(IN) :: sigml (kmax+1)
    LOGICAL        ,INTENT(IN   ) :: RESTART
    INTEGER        ,INTENT(IN   ) :: ibMax
    INTEGER        ,INTENT(IN   ) :: jbMax
    REAL(KIND=r8),INTENT(INOUT) :: USTAR  (1:ibMax,1:jbMax)
    INTEGER        ,INTENT(INOUT) :: LOWLYR (1:ibMax,1:jbMax)

    IF(sfcpbl == 1)THEN
       CALL InitSfc_MellorYamada0(kmax, sig, delsig, sigml)
    ELSE IF(sfcpbl == 2)THEN
       CALL InitSfc_MellorYamada1(kmax, sig, delsig, sigml,RESTART, &!(IN   )
            ibMax  , &!(IN   )
            jbMax  , &!(IN   )
            USTAR  , &!(INOUT)
            LOWLYR   )!(INOUT)
    ELSE IF(sfcpbl == 2)THEN
       WRITE(*,*)'CCM3 pbl parametrization'
    ELSE 
       WRITE(*,*)'its not set sfc pbl parametrization'
       STOP
    END IF

  END SUBROUTINE InitSfcPBL_Driver

  SUBROUTINE SfcPBL_Driver(gu     ,gv     ,gt   ,gq    ,&
       delsig ,ncols  ,kmax ,deltm ,&
       colrad ,tmsfc  ,qmsfc,umsfc ,&
       Mmlen  ,dt,&
       ITIMESTEP,topo   ,XLAND,LOWLYR,&
       gps    ,sl   ,si    ,&
       sigki  ,USTAR  ,tkemyj,qsfc ,&
       thz0   ,qz0    ,uz0  ,vz0   ,&
       znt    ,pblh   ,ELM,akhs ,akms  ,&
       latitu,TSK,CT,htdisp, &
       PBL_CoefKm, PBL_CoefKh, gice, gliq,gvar )
    IMPLICIT NONE
    INTEGER, INTENT(in        ) :: ncols
    INTEGER, INTENT(in        ) :: kmax
    REAL(KIND=r8),    INTENT(in   ) :: gu     (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: gv     (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: gt     (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: gq     (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: delsig (kmax)
    REAL(KIND=r8),    INTENT(in   ) :: deltm
    REAL(KIND=r8),    INTENT(in   ) :: colrad (ncols)
    REAL(KIND=r8),    INTENT(inout) :: tmsfc  (ncols,kmax,3)
    REAL(KIND=r8),    INTENT(inout) :: qmsfc  (ncols,kmax,5+nClass)
    REAL(KIND=r8),    INTENT(inout) :: umsfc  (ncols,kmax,4)
    REAL(KIND=r8),    INTENT(inout) :: Mmlen  (ncols)   
    REAL(KIND=r8),    INTENT(in   ) :: dt
    INTEGER      ,    INTENT(in         ) :: ITIMESTEP!-- itimestep         number of time steps
    REAL(KIND=r8),    INTENT(in   ) :: topo    (1:nCols)!"HGT" "Terrain Height"   "m"
    REAL(KIND=r8),    INTENT(in   ) :: XLAND (1:nCols)!land mask (1 for land, 2 for water)
    INTEGER        ,    INTENT(in   ) :: LOWLYR(1:nCols)!index of lowest model layer above ground
    REAL(KIND=r8),    INTENT(in   ) :: gps    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: sl     (kmax)
    REAL(KIND=r8),    INTENT(in   ) :: si     (kmax+1)
    REAL(KIND=r8),    INTENT(in   ) :: sigki  (kmax)
    REAL(KIND=r8),    INTENT(inout) :: USTAR  (nCols)! UST  u* in similarity theory (m/s)
    REAL(KIND=r8),    INTENT(inout) :: tkemyj (nCols,kMAx+1)! tkemyj           turbulence kinetic energy from Mellor-Yamada-Janjic (MYJ) (m^2/s^2)
    REAL(KIND=r8),    INTENT(inout) :: QSFC   (nCols)! qsfc specific humidity at lower boundary (kg/kg)
    REAL(KIND=r8),    INTENT(inout) :: THZ0   (nCols)! thz0 potential temperature at roughness length (K)
    REAL(KIND=r8),    INTENT(inout) :: QZ0    (nCols)
    REAL(KIND=r8),    INTENT(inout) :: UZ0    (nCols)! uz0  u wind component at roughness length (m/s)
    REAL(KIND=r8),    INTENT(inout) :: VZ0    (nCols)! vz0  v wind component at roughness length (m/s)
    REAL(KIND=r8),    INTENT(inout) :: ZNT    (nCols)! ZNT  time-varying roughness length (m)
    REAL(KIND=r8),    INTENT(inout) :: PBLH   (nCols)! PBLH PBL height (m)
    REAL(KIND=r8),    INTENT(inout) :: ELM    (nCols,kMAx)
    REAL(KIND=r8),    INTENT(inout) :: AKHS   (nCols)! akhs sfc exchange coefficient of heat/moisture from MYJ
    REAL(KIND=r8),    INTENT(inout) :: AKMS   (nCols)! akms sfc exchange coefficient of momentum from MYJ
    REAL(KIND=r8),    INTENT(inout) :: TSK  (1:nCols)!surface temperature (K)
    REAL(KIND=r8),    INTENT(inout) :: CT   (1:nCols)
    REAL(KIND=r8),    INTENT(in   ) :: htdisp (1:nCols)
    REAL(KIND=r8),    INTENT(INOUT) :: PBL_CoefKm(ncols, kmax)
    REAL(KIND=r8),    INTENT(INOUT) :: PBL_CoefKh(ncols, kmax)
    REAL(KIND=r8),    OPTIONAL,  INTENT(inout) :: gice (ncols,kMax)
    REAL(KIND=r8),    OPTIONAL,  INTENT(inout) :: gliq (ncols,kMax)
    REAL(KIND=r8),    OPTIONAL,  INTENT(inout) :: gvar (ncols,kMax,nClass)

    INTEGER, INTENT(in        ) :: latitu
    REAL(KIND=r8) :: MAVAIL(1:nCols)!surface moisture availability (between 0 and 1)
    REAL(KIND=r8) :: PMID  (1:nCols,1:kMAx)! p_phy         pressure (Pa)
    REAL(KIND=r8) :: PINT  (1:nCols,1:kMAx)! p8w           pressure at full levels (Pa)
    REAL(KIND=r8) :: QC    (1:nCols,1:kMAx)! "Cloud water mixing ratio"      " kg kg-1"
    REAL(KIND=r8) :: QV    (1:nCols,1:kMAx)! "Water vapor mixing ratio"      " kg kg-1"
    REAL(KIND=r8) :: T     (1:nCols,1:kMAx)!t_phy             temperature (K)
    REAL(KIND=r8) :: TH    (1:nCols,1:kMAx)!th_phy potential temperature (K)
    REAL(KIND=r8) :: U     (1:nCols,1:kMAx)!u_phy            u-velocity interpolated to theta points (m/s)
    REAL(KIND=r8) :: V     (1:nCols,1:kMAx)!v_phy            v-velocity interpolated to theta points (m/s)
    !
    REAL(KIND=r8) :: FLX_LH(1:nCols)!-- LH            net upward latent heat flux at surface (W/m^2)
    REAL(KIND=r8) :: HFX   (1:nCols)!-- HFX           net upward heat flux at the surface (W/m^2)
    REAL(KIND=r8) :: PSHLTR(1:nCols)!-- pshltr        diagnostic shelter (2m) pressure from MYJ (Pa)
    REAL(KIND=r8) :: QFX   (1:nCols)!-- QFX           net upward moisture flux at the surface (kg/m^2/s)
    REAL(KIND=r8) :: QSHLTR(1:nCols)!-- qshltr        diagnostic 2-m specific humidity from MYJ 
    REAL(KIND=r8) :: TSHLTR(1:nCols)!-- tshltr        diagnostic 2-m theta from MYJ
    REAL(KIND=r8) :: TH10  (1:nCols)!-- th10          diagnostic 10-m theta from MYJ
    REAL(KIND=r8) :: Q10   (1:nCols)!-- q10           diagnostic 10-m specific humidity from MYJ
    REAL(KIND=r8) :: U10   (1:nCols)! u10 diagnostic 10-m u component from surface layer
    REAL(KIND=r8) :: V10   (1:nCols)! v10 diagnostic 10-m v component from surface layer
    !
    !
    !
    REAL(KIND=r8) :: CHS  (1:nCols)! surface exchange coefficient for heat and moisture (m s-1)'
    REAL(KIND=r8) :: CHS2 (1:nCols)! 2m surface exchange coefficient for heat  (m s-1)'
    REAL(KIND=r8) :: CQS2 (1:nCols)! 2m surface exchange coefficient for moisture (m s-1)'
    REAL(KIND=r8) :: CPM  (1:nCols)
    REAL(KIND=r8) :: FLHC (1:nCols)
    REAL(KIND=r8) :: FLQC (1:nCols)
    REAL(KIND=r8) :: QGH  (1:nCols)
    REAL(KIND=r8) :: bps  (1:nCols,kMax)! Factor conversion to potention temperature
    REAL(KIND=r8) :: XMAVA,rbyg 
    REAL(KIND=r8) :: r100
    REAL(KIND=r8) :: con1  (kMax)
    REAL(KIND=r8) :: t0    (kMax)
    REAL(KIND=r8) :: t1    (kMax)
    REAL(KIND=r8) :: Pbl_ATemp(1:nCols,kMax)
    REAL(KIND=r8) :: Pbl_ITemp(1:nCols,kMax)
    REAL(KIND=r8) :: psur(1:nCols)
    REAL(KIND=r8) :: delz    (nCols,kMax)
    REAL(KIND=r8) :: press   (nCols,kMax)
    REAL(KIND=r8) :: tv      (nCols,kMax)
    REAL(KIND=r8) :: ze      (nCols,kMax)
    REAL(KIND=r8) :: terr(1:nCols)   
    REAL(KIND=r8) :: tkemyj_local  (1:nCols,kMax+1) 

    INTEGER :: i,k 
    REAL(KIND=r8) , PARAMETER :: rd=287.0_r8
    tkemyj_local=0.0_r8
    tkemyj=0.0_r8
    IF(sfcpbl == 1)THEN
    IF (microphys) THEN
       !PRINT*,' SfcPbl_MY0 begin'
      IF(nClass>0 .and. PRESENT(gvar))THEN
         CALL SfcPbl_MY0( &
            gu    (1:nCols,:) ,gv  (1:nCols,:)    ,gt(1:nCols,:)     ,gq(1:nCols,:)     ,&
            delsig(1:kmax   ) ,ncols              ,kmax              ,deltm                 ,&
            colrad(1:nCols  ) ,tmsfc(1:nCols,:,:) ,qmsfc(1:nCols,:,1:5+nClass),umsfc(1:nCols,:,:),&
            Mmlen (1:ncols  ) ,PBL_CoefKm(1:nCols,:),PBL_CoefKh(1:nCols,:),pblh(1:nCols  ),tkemyj_local,&
            gice  (1:nCols,:) ,gliq(1:nCols,:)     ,gvar(1:nCols,:,:))
       !PRINT*,' SfcPbl_MY0 end'
      ELSE
         CALL SfcPbl_MY0( &
            gu    (1:nCols,:) ,gv  (1:nCols,:)    ,gt(1:nCols,:)     ,gq(1:nCols,:)     ,&
            delsig(1:kmax   ) ,ncols              ,kmax              ,deltm                 ,&
            colrad(1:nCols  ) ,tmsfc(1:nCols,:,:) ,qmsfc(1:nCols,:,1:5+nClass),umsfc(1:nCols,:,:),&
            Mmlen (1:ncols  ) ,PBL_CoefKm(1:nCols,:),PBL_CoefKh(1:nCols,:),pblh(1:nCols  ),tkemyj_local,&
            gice  (1:nCols,:) ,gliq(1:nCols,:)    )

      END IF 
    ELSE
       CALL SfcPbl_MY0( &
            gu    (1:nCols,:) ,gv  (1:nCols,:)    ,gt(1:nCols,:)     ,gq(1:nCols,:)     ,&
            delsig(1:kmax   ) ,ncols                  ,kmax                     ,deltm                 ,&
            colrad(1:nCols  ) ,tmsfc(1:nCols,:,:) ,qmsfc(1:nCols,:,1:5+nClass),umsfc(1:nCols,:,:),&
            Mmlen (1:ncols  ) ,PBL_CoefKm(1:nCols,:),PBL_CoefKh(1:nCols,:),pblh(1:nCols  ),tkemyj_local)
    END IF

    DO k=1,kMax+1
       DO i=1,nCols
          tkemyj    (i,k)= tkemyj_local(i,k)
       END DO
    END DO

    ELSE IF(sfcpbl == 2)THEN

       r100=100.0e0_r8 /gasr
       DO k=1,kMAx-1
          con1  (        k)=grav*si(k+1)/(gasr*(sl(k)-sl(k+1)))
          con1  (        k)=con1(k)*con1(k)
          t0    (        k)=(sl(k+1)-si(k+1))/(sl(k+1)-sl(k))
          t1    (        k)=(si(k+1)-sl(k  ))/(sl(k+1)-sl(k))        
       END DO
       DO k = 1, kmax-1
          DO i = 1, ncols
             Pbl_ATemp(i,k)  = t0(k)*gt(i,k) + t1(k)*gt(i,k+1)
             Pbl_ITemp(i,k)  = 1.0_r8/(Pbl_ATemp(i,k)*Pbl_ATemp(i,k))
          END DO
       END DO
       !
       ! --  -- 2                 --                     -- 2
       !|        1  |         1        | g          si(k+1)      | 
       !| -----|  = ----- * |--- * ----------------| 
       !|  DZ  |         T*T        | R        sl(k) -sl(k+1) | 
       ! --  --                 --                     --  
       !
       ! --  --                     --                     -- 
       !|   1  |      1            | g          si(k+1)      | 
       !| -----|  = ----- * |--- * ------------| = SQRT((Pbl_ITemp(i,k) * con1(k)))
       !|  DZ  |      T            | R        sl(k) -sl(k+1) | 
       ! --  --                     --                     --  
       !
       !                  --                      -- 
       !           1         | g         si(k+1)        |                   1
       !DZ = ----- * |--- * ----------------| = ---------------------------------
       !           T         | R   sl(k) -sl(k+1)        |   SQRT((Pbl_ITemp(i,k) * con1(k)))
       !                  --                      --  
    DO i=1,nCols
       psur(i)   =gps(i)*100.0_r8
       terr(i)   =MAX(topo(i),0.0_r8)   
       !                --                         --  2
       !               |   grav        si(k+1)     |
       !   con1  (k) = |------- * ---------------- |             
       !               |  gasr      sl(k)-sl(k+1)) |
       !                --                         -- 
       !
       !               (m*m)
       ! dzm (i) = --------- = m
       !                 m
    END DO
    DO k=1,kMax
       DO i=1,nCols
          press(i,k)=gps(i)*sl(k)
          press(i,k)=press(i,k)*100.0_r8
          tv(i,k)=gt(i,k)*(1.0_r8+0.608_r8*gq(i,k))
       END DO
    END DO    
    DO i=1,nCols
      rbyg=gasr/grav*delsig(1)*0.5e0_r8
      IF(XLAND(i) >1.0_r8)THEN
         delz(i,1)=MAX((rbyg * tv(i,1) - htdisp(i)),0.5_r8)*0.75_r8
      ELSE
         delz(i,1)=MAX((rbyg * tv(i,1) - htdisp(i)),0.5_r8)      
      END IF
      ze  (i,1)= delz(i,1)!gt(i,1)*(gasr/grav)*(psur(i)-press(i,1))/psur(i)
      ze  (i,1)= terr(i)+ze(i,1)

    END DO
    DO k=2,kMax
       DO i=1,nCols
          delz(i,k)=0.5_r8*rd*(tv(i,k-1)+tv(i,k))* &
                    LOG(press(i,k-1)/press(i,k))/grav
          ze(i,k)=ze(i,k-1)+ delz(i,k)
       END DO
    END DO


       DO k=1,kMAx
          DO i=1,nCols
             !
             ! Factor conversion to potention temperature
             !
             bps (i,k)=sigki(k) 
             PMID(i,k)=100.0_r8*gps(i)*sl(k) 
             PINT(i,k)=100.0_r8*gps(i)*si(k)
             T   (i,k)=gt (i,k)
             TH  (i,k)=T  (i,k)*bps(i,k)
             QV  (i,k)=gq (i,k)
             QC  (i,k)=0.0_r8
             U   (i,k)=gu (i,k)/SIN( colrad(i))
             V   (i,k)=gv (i,k)/SIN( colrad(i))
          END DO
       END DO

       XMAVA=1.0_r8  
       DO i=1,nCols
          !..delsig     k=2  ****gu,gv,gt,gq,gyu,gyv,gtd,gqd,sl*** } delsig(2)
          !             k=3/2----si,ric,rf,km,kh,b,l -----------
          !             k=1  ****gu,gv,gt,gq,gyu,gyv,gtd,gqd,sl*** } delsig(1)
          !             k=1/2----si ----------------------------
          !
          ZNT (i)=MAX(ZNT(i),0.1e-7_r8)
          IF(XLAND(i) .LT. 1.5_r8)THEN
             MAVAIL(i)=XMAVA
          ELSE
             MAVAIL(i)=1.0_r8
          ENDIF
       ENDDO
       CALL SfcPbl_MYJ1(ITIMESTEP,&! INTENT(IN) :: !-- itimestep     number of time steps
            ze(:,1) ,&! INTENT(IN) :: HT     (iMax0:iMAx)!"HGT" "Terrain Height"   "m"   
            delz    ,&! INTENT(IN) :: DZ     (iMax0:iMAx,kMax0:kMAx)! dz between full levels (m)                               
            PMID    ,&! INTENT(IN) :: PMID   (iMax0:iMAx,kMax0:kMAx)! p_phy         pressure (Pa)
            PINT    ,&! INTENT(IN) :: PINT   (iMax0:iMAx,kMax0:kMAx)! p8w           pressure at full levels (Pa)
            TH      ,&! INTENT(IN) :: TH     (iMax0:iMAx,kMax0:kMAx)!th_phy potential temperature (K)
            T       ,&! INTENT(IN) :: T      (iMax0:iMAx,kMax0:kMAx)!t_phy            temperature (K)
            QV      ,&! INTENT(IN) :: QV     (iMax0:iMAx,kMax0:kMAx)! "Water vapor mixing ratio"      "kg kg-1"
            QC      ,&! INTENT(IN) :: QC     (iMax0:iMAx,kMax0:kMAx)! "Cloud water mixing ratio"      "kg kg-1"
            U         ,&! INTENT(IN) :: U         (iMax0:iMAx,kMax0:kMAx)!u_phy            u-velocity interpolated to theta points (m/s)
            V         ,&! INTENT(IN) :: V         (iMax0:iMAx,kMax0:kMAx)!v_phy            v-velocity interpolated to theta points (m/s)
            tkemyj_local(1:nCols,1:kMax)  ,&! INTENT(IN) :: tkemyj (iMax0:iMAx,kMax0:kMAx)! tke_myj turbulence kinetic energy from Mellor-Yamada-Janjic (MYJ) (m^2/s^2)               
            TSK     ,&! INTENT(IN) :: TSK    (iMax0:iMAx)!surface temperature (K)
            QSFC    ,&! INTENT(INOUT) :: QSFC (iMax0:iMAx)! qsfc specific humidity at lower boundary (kg/kg)
            THZ0    ,&! INTENT(INOUT) :: THZ0 (iMax0:iMAx)! thz0 potential temperature at roughness length (K)
            QZ0     ,&! INTENT(INOUT) :: QZ0 (iMax0:iMAx) ! QZ0 Water vapor mixing ratio at roughness length(Kg/kg)
            UZ0         ,&! INTENT(INOUT) :: UZ0  (iMax0:iMAx)! uz0  u wind component at roughness length (m/s)
            VZ0         ,&! INTENT(INOUT) :: VZ0  (iMax0:iMAx)! vz0  v wind component at roughness length (m/s)   
            LOWLYR  ,&! INTENT(IN) :: LOWLYR(iMax0:iMAx)!index of lowest model layer above ground
            XLAND         ,&! INTENT(IN) :: XLAND  (iMax0:iMAx)!land mask (1 for land, 2 for water)                               
            USTAR   ,&! INTENT(INOUT) :: USTAR(iMax0:iMAx)! UST  u* in similarity theory (m/s)
            ZNT     ,&! INTENT(INOUT) :: ZNT  (iMax0:iMAx)! ZNT  time-varying roughness length (m)
            PBLH    ,&! INTENT(INOUT) :: PBLH (iMax0:iMAx)! PBLH PBL height (m)
            ELM     ,&!
            MAVAIL  ,&! INTENT(IN) :: MAVAIL (iMax0:iMAx)!surface moisture availability (between 0 and 1)
            AKHS         ,&! INTENT(INOUT) :: AKHS (iMax0:iMAx)! akhs sfc exchange coefficient of heat/moisture from MYJ
            AKMS         ,&! INTENT(INOUT) :: AKMS (iMax0:iMAx)! akms sfc exchange coefficient of momentum from MYJ                             
            CHS         ,&! INTENT(OUT )  :: CHS  (iMax0:iMAx)
            CHS2         ,&! INTENT(OUT )  :: CHS2 (iMax0:iMAx)
            CQS2         ,&! INTENT(OUT )  :: CQS2 (iMax0:iMAx)
            HFX     ,&! INTENT(OUT) :: HFX   (iMax0:iMAx)!-- HFX net upward heat flux at the surface (W/m^2)
            QFX     ,&! INTENT(OUT) :: QFX   (iMax0:iMAx)!-- QFX net upward moisture flux at the surface (kg/m^2/s)
            FLX_LH  ,&! INTENT(OUT) :: FLX_LH(iMax0:iMAx)!-- LH  net upward latent heat flux at surface (W/m^2)
            FLHC         ,&! INTENT(OUT )  :: FLHC (iMax0:iMAx)
            FLQC    ,&! INTENT(OUT )  :: FLQC (iMax0:iMAx)
            QGH     ,&! INTENT(OUT )  :: QGH  (iMax0:iMAx)
            CPM     ,&! INTENT(OUT )  :: CPM  (iMax0:iMAx)
            CT         ,&! INTENT(OUT )  :: CT   (iMax0:iMAx)                               
            U10         ,&! INTENT(OUT) :: U10        (iMax0:iMAx)! u10 diagnostic 10-m u component from surface layer
            V10         ,&! INTENT(OUT) :: V10        (iMax0:iMAx)! v10 diagnostic 10-m v component from surface layer
            TSHLTR  ,&! INTENT(OUT) :: TSHLTR(iMax0:iMAx)!-- TH02        diagnostic 2-m theta from MYJ
            TH10    ,&! INTENT(OUT) :: TH10  (iMax0:iMAx)!-- th10          diagnostic 10-m theta from MYJ
            QSHLTR  ,&! INTENT(OUT) :: QSHLTR(iMax0:iMAx)!-- Q02        diagnostic 2-m specific humidity from MYJ 
            Q10     ,&! INTENT(OUT) :: Q10   (iMax0:iMAx)!-- q10           diagnostic 10-m specific humidity from MYJ
            PSHLTR         ,&! INTENT(OUT) :: PSHLTR(iMax0:iMAx)!-- pshltr        diagnostic shelter (2m) pressure from MYJ (Pa)    
            tmsfc   ,&
            qmsfc(1:nCols,:,1:5+nClass)   ,&
            umsfc   ,&
            gt      ,&
            gq      ,&
            gu      ,&
            gv      ,&
            bps     ,&
            PBL_CoefKm,&
            PBL_CoefKh,&
            latitu  ,&
            dt      ,&
            nCols   ,&
            kMax    ,&
            nClass)           
    DO k=1,kMax+1
       DO i=1,nCols
          tkemyj    (i,k)= tkemyj_local(i,k)
       END DO
    END DO

    ELSE
       WRITE(*,*)'its not set sfc pbl parametrization'
       STOP     
    END IF


  END SUBROUTINE SfcPBL_Driver
END MODULE SfcPBLDriver
