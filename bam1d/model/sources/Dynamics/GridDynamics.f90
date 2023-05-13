
!  $Author: pkubota $
!  $Date: 2010/04/20 20:18:04 $
!  $Revision: 1.18 $
!
!  modified by enver.ramirez
!  to read forcing file
!  Date: 2015/06/19
! 
!  modified by enver.ramirez
!  for compatibility with the GASS/DCP spec. sfc protocol
!  Date: Jan 2020
!
MODULE GridDynamics
       
  USE Constants,   ONLY : &
       i8,r8,i4,r4,&
       tov,               & ! intent(in)
       qmin                 ! intent(in)
   
  USE PhysicsDriver, ONLY : &
      DryPhysics

  USE GridHistory, ONLY:  &
       StoreGridHistory , &
       IsGridHistoryOn  , &
       nGHis_temper    , &
       nGHis_uzonal    , &
       nGHis_vmerid    , &
       nGHis_spchum    , &
       nGHis_tvirsf    , &
       nGHis_uzonsf    , &
       nGHis_vmersf    , &
       nGHis_sphusf    , &
       nGHis_snowdp    , &
       nGHis_rouglg    , &
       nGHis_tcanop    , &
       nGHis_tgfccv    , &
       nGHis_tgdeep    , &
       nGHis_swtsfz    , &
       nGHis_swtrtz    , &
       nGHis_swtrcz    , &
       nGHis_mostca    , &
       nGHis_mostgc    , &
       nGHis_vegtyp    , &
       dogrh

  USE Diagnostics, ONLY:  &
      pwater            , &
      updia             , &
      dodia             , &
      nDiag_sigdot            , &
      nDiag_omegav             , &
      nDiag_pwater            , &
      nDiag_divgxq              , &
      nDiag_vmoadv            , &
      nDiag_tmtdps              , &
      nDiag_tgfccv              , &
      nDiag_tcanop  , &
      nDiag_TEMPFORC              , &
      nDiag_MOISFORC


  USE FieldsPhysics, ONLY: &
       zorl              , &
       sheleg            , &
       imask             , &
       gtsea             , &
       tcm               , &
       tgm               , &
       tdm               , &
       wm                , &
       capacm

  USE Options, ONLY:       &
      microphys         , &
       nClass            , &
      isimp             , &
      istrt,  & !Enver
      intfor, ndg, specSfc, iffor, jdt, path_in, forcings_weight_d, forcings_weight_m, forcings_weight_t

!Enver additions 0 start
  USE Sizes, ONLY: si

  IMPLICIT NONE
  PRIVATE
   REAL(KIND=r8)   , ALLOCATABLE     :: pres (:)
   REAL(KIND=r8)   , ALLOCATABLE     :: omeg (:)
   REAL(KIND=r8)   , ALLOCATABLE     :: auvel (:)
   REAL(KIND=r8)   , ALLOCATABLE     :: avvel (:)
   REAL(KIND=r8)   , ALLOCATABLE     :: atemp (:)
   REAL(KIND=r8)   , ALLOCATABLE     :: aumid (:)
   REAL(KIND=r8)   , ALLOCATABLE     :: sigma(:)  !sigma = (P - Pt)/(Ps - Pt) 
   REAL(KIND=r8)   , ALLOCATABLE     :: sonda (:,:)
   REAL(KIND=r8)   , ALLOCATABLE     :: sond  (:,:)
   REAL(KIND=r8)   , ALLOCATABLE     :: sonda2 (:,:)
   REAL(KIND=r8)   , ALLOCATABLE     :: sond2  (:,:)

   REAL(KIND=r8)                     :: Ps1, Ps2

   REAL(KIND=r8)   , ALLOCATABLE     :: sigma_ndg(:)  !sigma = (P - Pt)/(Ps - Pt) 
   REAL(KIND=r8)   , ALLOCATABLE     :: sondg  (:,:)
   REAL(KIND=r8)   , ALLOCATABLE     :: sondg2  (:,:)
   REAL(KIND=r8)   , ALLOCATABLE     :: sondga (:,:)
   REAL(KIND=r8)   , ALLOCATABLE     :: sondga2 (:,:)


   REAL(KIND=r8)   , ALLOCATABLE     :: sigma_geo(:)  !sigma = (P - Pt)/(Ps - Pt) 
   REAL(KIND=r8)   , ALLOCATABLE     :: sogeo  (:,:)
   REAL(KIND=r8)   , ALLOCATABLE     :: sogeo2  (:,:)
   REAL(KIND=r8)   , ALLOCATABLE     :: sogeoa (:,:)
   REAL(KIND=r8)   , ALLOCATABLE     :: sogeoa2 (:,:)

   REAL(KIND=r8)   , ALLOCATABLE     :: gps  (:)  !gps  = sigma*(Ps - Pt) + Pt

  ! For specified Surface
   REAL(KIND=r8)                     :: specSfc_Tsfc
   REAL(KIND=r8)                     :: specSfc_H
   REAL(KIND=r8)                     :: specSfc_LE
   REAL(KIND=r8)                     :: specSfc_Tau
   REAL(KIND=r8)                     :: specSfc_Tsfc2
   REAL(KIND=r8)                     :: specSfc_H2
   REAL(KIND=r8)                     :: specSfc_LE2
   REAL(KIND=r8)                     :: specSfc_Tau2

   LOGICAL  :: OpenFile =.TRUE.
   LOGICAL  :: OpenFile_specSfc =.TRUE.
   LOGICAL  :: OpenFile_ndg =.TRUE.
   LOGICAL  :: OpenFile_geo =.TRUE.

  PUBLIC :: GrpComp, Forcing_In, Forcings_Ascii, Forcings_Ascii_lsf, Nudging_Ascii, Geostrophy_Ascii, specSfc_Ascii
!Enver additions 0 end
CONTAINS
  SUBROUTINE GrpComp(&
       gyu   , gyv   , gtd    , &
       gqd   , gtmp  , gq     , &
       gum   , gvm   , gtm    , &
       gqm   , omg   , ps     , &
       glnpm , colrad, rcl    , &
       gzs  ,&
       dt    , ifday , tod    , &
       idate , idatec, ibMax  , &
       kMax  , ibLim , slagr  , &
       jdt   , kt    , ktm    , &
       ktp   , ib    , jb     , &
       lonrad, cos2d , intcosz, &
       guo, gvo, gto, gqo,      &  ! for ndg
       gTsfc, gH, gLE, gTau,    &  ! for specSfc  
       gug, gvg,                &  ! for geo
       gicem  , gicet  , gliqm  , gliqt  , &
       gvarm  , gvart)
    !
    ! grpcomp: grid-point computations (all tendencies are computed) 
    !
    !
    ! slagr is the option for eulerian (slagr=.false.) or
    ! semi-Lagrangian integration (slagr=.true.)
    !
    !
    !
    CHARACTER(LEN=200)     :: ppath !Enver added ppath
    INTEGER, INTENT(IN   ) :: ibMax
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: ibLim
    LOGICAL(KIND=i8), INTENT(IN   ) :: slagr
    REAL(KIND=r8),    INTENT(IN   ) :: rcl    (ibMax)
    REAL(KIND=r8),    INTENT(IN   ) :: dt
    REAL(KIND=r8),    INTENT(IN   ) :: colrad (ibMax)
    REAL(KIND=r8),    INTENT(INOUT) :: gtmp   (ibMax,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: ps     (ibMax)
    REAL(KIND=r8),    INTENT(INOUT) :: gum    (ibMax,kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: gvm    (ibMax,kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: gtm    (ibMax,kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: gqm    (ibMax,kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: omg    (ibMax,kMax) 
    REAL(KIND=r8),    INTENT(INOUT) :: glnpm  (ibMax)
    REAL(KIND=r8),    INTENT(OUT  ) :: gyu    (ibMax,kMax)
    REAL(KIND=r8),    INTENT(OUT  ) :: gyv    (ibMax,kMax)
    REAL(KIND=r8),    INTENT(OUT  ) :: gtd    (ibMax,kMax)
    REAL(KIND=r8),    INTENT(OUT  ) :: gqd    (ibMax,kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: gq     (ibMax,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: gzs    (ibMax)
    INTEGER, INTENT(IN   ) :: ifday
    REAL(KIND=r8),    INTENT(IN   ) :: tod
    INTEGER, INTENT(IN   ) :: idate  (4)
    INTEGER, INTENT(IN   ) :: idatec (4)
    INTEGER, INTENT(IN   ) :: jdt
    INTEGER, INTENT(IN   ) :: kt
    INTEGER, INTENT(IN   ) :: ktm
    INTEGER, INTENT(IN   ) :: ktp
    INTEGER, INTENT(IN   ) :: jb
    INTEGER, INTENT(IN   ) :: ib
    REAL(KIND=r8)   , INTENT(IN   ) :: lonrad (ibMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: cos2d  (ibMax)    
    LOGICAL, INTENT(IN   ) :: intcosz
    REAL(KIND=r8),    OPTIONAL,   INTENT(INOUT) :: gicem  (ibMax, kMax)
    REAL(KIND=r8),    OPTIONAL,   INTENT(INOUT) :: gicet  (ibMax, kMax)
    REAL(KIND=r8),    OPTIONAL,   INTENT(INOUT) :: gliqm  (ibMax, kMax)
    REAL(KIND=r8),    OPTIONAL,   INTENT(INOUT) :: gliqt  (ibMax, kMax)
    REAL(KIND=r8),    OPTIONAL,   INTENT(INOUT) :: gvarm  (ibMax, kMax,nClass)
    REAL(KIND=r8),    OPTIONAL,   INTENT(INOUT) :: gvart  (ibMax, kMax,nClass)

    REAL(KIND=r8),    INTENT(INOUT) :: guo    (ibMax,kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: gvo    (ibMax,kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: gto    (ibMax,kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: gqo    (ibMax,kMax)
    !
    !   gTsfc, gH, gLE, gTau,    &  ! for specSfc  
    REAL(KIND=r8),    INTENT(INOUT) :: gTsfc (ibMax)
    REAL(KIND=r8),    INTENT(INOUT) :: gH    (ibMax)
    REAL(KIND=r8),    INTENT(INOUT) :: gLE   (ibMax)
    REAL(KIND=r8),    INTENT(INOUT) :: gTau  (ibMax)

    REAL(KIND=r8),    INTENT(INOUT) :: gug    (ibMax,kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: gvg    (ibMax,kMax)
    !
    !  local variables 
    !
    REAL(KIND=r8)   , DIMENSION(ibMax,kMax) :: gdt
    REAL(KIND=r8)   , DIMENSION(ibMax,kMax) :: psint
    REAL(KIND=r8)   , DIMENSION(ibMax,kMax) :: divint
    REAL(KIND=r8)   , DIMENSION(ibMax,kMax) :: dot   
    REAL(KIND=r8)                           :: gps(ibMax)
    REAL(KIND=r8)                           :: one
    INTEGER                        :: i
    INTEGER                        :: k
    INTEGER                        :: ncount
    INTEGER                        :: latco
    REAL(KIND=r8)                           :: xxday !Enver
    REAL(KIND=r8)                           :: xxday2 !Enver

!   for specSfc
    REAL(KIND=r8)                           :: specSfc_day  
    REAL(KIND=r8)                           :: specSfc_day2 
!--------------------------------------------------------------------------
!	 locations for available diagnostics in this subroutine
!-------------------------------------------------------------------------- 
! INTEGER, PARAMETER :: ndpsm  = 1  ! time mean surface pressure
! INTEGER, PARAMETER :: nDiag_sigdot = 8  ! sigma dot
! INTEGER, PARAMETER :: nDiag_omegav  = 7  ! omega
! INTEGER, PARAMETER :: nDiag_pwater = 14 ! precipitable water  			
! INTEGER, PARAMETER :: nDiag_divgxq   = 51 ! divergence * specific humidity
! INTEGER, PARAMETER :: nDiag_vmoadv = 52 ! vertical moisture advection
! INTEGER, PARAMETER :: nDiag_tmtdps   = 64 ! time mean deep soil temperature
! INTEGER, PARAMETER :: nDiag_tgfccv   = 65 ! ground/surface cover temperature
! INTEGER, PARAMETER :: nDiag_tcanop   = 66 ! canopy temperature
!--------------------------------------------------------------------------              
    one=1.0_r8
    psint=0.0_r8
    divint=0.0_r8
    dot=0.0_r8
    omg=0.0_r8
    gyu=0.0_r8
    gyv=0.0_r8
    gtd=0.0_r8
    gqd=0.0_r8
    glnpm=0.0_r8

    guo=0.0_r8
    gvo=0.0_r8
    gto=0.0_r8
    gqo=0.0_r8
   
    !for specSfc 
    !IF ( specSfc ) THEN
     gTsfc=0.0_r8
     gH=0.0_r8
     gLE=0.0_r8
     gTau=0.0_r8
    !ENDIF

    gug=0.0_r8
    gvg=0.0_r8

    IF(microphys)THEN
       IF(nClass>0)THEN
        gvart=0.0_r8;  gliqt=0.0_r8; gicet=0.0_r8
       ELSE
        gliqt=0.0_r8;gicet=0.0_r8
       END IF
    END IF
    latco=jb
    ppath = path_in
!
! Enver additions 1 start
  IF ( iffor == -1. ) THEN
    ! Forcing data set is fixed only one profile is read and
    ! use for the whole run
    CALL Forcing_In(one,psint,divint,dot,omg,gyu,gyv,gtd, &
                       gqd,ibmax,kmax,ppath)
  ELSEIF ( iffor == 0. ) THEN
    ! Forcing is linearly interpolated to current day and time forcing
    ! is assumed to be spaced each intfor seconds. Forcing file is
    ! in ASCII code
    ! print*, 'Inside GridDynamics, iffor == 0 '
    !CALL LocalJulian(idatec(4),idatec(2),idatec(3),idatec(1),tod,xxday)
    !print*, 'xxday=',xxday, 'tod=', tod, 'intfor=', intfor
    !  !open(500,file='times',form='formatted')
    !  ! write(500,'(i4,3i2,a4,4i4,f10.2,f10.2,f11.4,i3,i5)') idate(4), idate(2),
    !  ! idate(3), idate(1), "to",  idatec(4), idatec(2), idatec(3), idatec(1),
    !  ! tod, mod(tod,3600.0), xxday, iffor, intfor
    !  !close(500)
    CALL LocalJulian_Year(idate(4),idate(2),idate(3),idate(1), & 
                          idatec(4),idatec(2),idatec(3),idatec(1),tod,xxday)


    CALL Forcings_Ascii(one,psint,divint,dot,omg,gyu,gyv,gtd, &
                       gqd,glnpm,ibmax,kmax,ppath,tod,jdt,xxday,colrad) !removed xxday-1


    ! ndg = .TRUE.  !This line might be deleted when ndg option implemented 
    ndg = .TRUE.  !This line might be deleted when ndg option implemented - DENIS 
    IF ( ndg .eqv. .TRUE. ) THEN
    ! In addition to tendencies (gyu,gyv,gtd,gqd) we also read observed full
    ! values to perform nudging (Ndg) full values for those variables are stored in
    ! guo, gvo, gto, gqo
     CALL Nudging_Ascii(guo,gvo,gto,gqo,ibmax,kmax,ppath,tod,jdt,xxday,colrad) !removed xxday-1
    !
    ! CALL Geostrophy_Ascii(gug,gvg,ibmax,kmax,ppath,tod,jdt,xxday,colrad) !removed xxday-1
    ENDIF


    IF ( specSfc ) THEN
     CALL specSfc_Ascii(gTsfc,gH,gLE,gTau,ppath,ibmax,jdt,xxday) !removed xxday-1
    ENDIF

  ELSEIF ( iffor == 1. ) THEN
    ! Instead of u- and v- tendecies (gyu, gyv) we insert full u and v (gum, gvm)
    CALL LocalJulian_Year(idate(4),idate(2),idate(3),idate(1), & 
                          idatec(4),idatec(2),idatec(3),idatec(1),tod,xxday)
    !print*, 'xxday2=',xxday2, 'tod=', tod, 'intfor=', intfor

    CALL Forcings_Ascii_lsf(one,psint,divint,dot,omg,gum,gvm,gtd, &
                       gqd,ibmax,kmax,ppath,tod,jdt,xxday,colrad) !removed xxday-1
  
    ! Nudging not implemented here because u and v full values are being read
    ! another method must be used.
    ! Forcings_Ascii_lsf need a revision as: (gum, gvm) are time 'T-1', (gtd, gqd) are time 'T'

 ! for specSfc
    IF ( specSfc ) THEN
     CALL specSfc_Ascii(gTsfc,gH,gLE,gTau,ppath,ibmax,jdt,xxday) !removed xxday-1
    ENDIF

  ENDIF
! Enver additions 1 end
!
  !	
  !
  ! enforce humidity to be above a certain level (avoid negative values...)
  ! ----------------------------------------------------------------------
  !
  IF(isimp.ne.'YES ')THEN  
    gq=MAX(gq,qmin)
  ENDIF  
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  
  IF(dodia(nDiag_tmtdps).or.dodia(nDiag_tgfccv).or.dodia(nDiag_tcanop))THEN
      ncount=0
      DO i=1,ibLim
        gdt(i,1)=gtsea (i,jb)
        gdt(i,2)=gtsea (i,jb)
        gdt(i,3)=gtsea (i,jb)
        IF(imask(i,jb).ge.1)THEN
           ncount=ncount+1
           gdt(i,1)=tcm(i,jb)
           gdt(i,2)=tgm(i,jb)
           gdt(i,3)=tdm(i,jb)
        END IF
      END DO
  END IF
  
  IF(dodia(nDiag_tmtdps))CALL updia(gdt(:,3:3),nDiag_tmtdps,latco)
  IF(dodia(nDiag_tgfccv))CALL updia(gdt(:,2:2),nDiag_tgfccv,latco)
  IF(dodia(nDiag_tcanop))CALL updia(gdt(:,1:1),nDiag_tcanop,latco)
  !
  !     obtain grid history fields if requested
  !
  IF(IsGridHistoryOn())THEN
      DO k=1,kMax
        DO i=1,ibLim
           gtd(i,k)=tov(k)+gtmp(i,k)
        END DO
      END DO
      IF(dogrh(nGHis_temper,  latco)) CALL StoreGridHistory (gtd,                nGHis_temper,  latco)
!-1D      IF(dogrh(nGHis_uzonal,  latco)) CALL StoreGridHistory (gu,                 nGHis_uzonal,  latco, sqrt(rcl))
!-1D      IF(dogrh(nGHis_vmerid,  latco)) CALL StoreGridHistory (gv,                 nGHis_vmerid,  latco, sqrt(rcl))
      IF(dogrh(nGHis_spchum,  latco)) CALL StoreGridHistory (gq,                 nGHis_spchum,  latco)
      IF(dogrh(nGHis_tvirsf, latco)) CALL StoreGridHistory (gtd,                nGHis_tvirsf, latco)
!-1D       IF(dogrh(nGHis_uzonsf, latco)) CALL StoreGridHistory (gu,                 nGHis_uzonsf, latco, sqrt(rcl))
!-1D      IF(dogrh(nGHis_vmersf, latco)) CALL StoreGridHistory (gv,                 nGHis_vmersf, latco, sqrt(rcl))
      IF(dogrh(nGHis_sphusf, latco)) CALL StoreGridHistory (gq,                 nGHis_sphusf, latco)
      IF(dogrh(nGHis_snowdp, latco)) CALL StoreGridHistory (sheleg(:,latco), nGHis_snowdp, latco)
      IF(dogrh(nGHis_rouglg, latco)) CALL StoreGridHistory (zorl(:,latco),   nGHis_rouglg, latco)
      ncount=0
      DO i=1,ibLim
        gtd(i,1)=gtsea(i,latco)
        gtd(i,2)=gtsea(i,latco) 
        gtd(i,3)=gtsea(i,latco)
        gtd(i,4)=one
        gtd(i,5)=one
        gtd(i,6)=one
        gtd(i,7)=0.0e0_r8
        gtd(i,8)=0.0e0_r8
        gtd(i,9)=imask(i,latco)
        IF(imask(i,latco).ge.1)THEN
           ncount=ncount+1
           gtd(i,1)=tcm(i,latco)
           gtd(i,2)=tgm(i,latco)
           gtd(i,3)=tdm(i,latco)
           gtd(i,4)=wm (i,1,latco)
           gtd(i,5)=wm (i,2,latco)
           gtd(i,6)=wm (i,3,latco)
           gtd(i,7)=capacm(i,1,latco)
           gtd(i,8)=capacm(i,2,latco)
        END IF
      END DO
 
      IF(dogrh(nGHis_tcanop,    latco)) CALL StoreGridHistory (gtd(:,1), nGHis_tcanop,    latco)
      IF(dogrh(nGHis_tgfccv,    latco)) CALL StoreGridHistory (gtd(:,2), nGHis_tgfccv,    latco)
      IF(dogrh(nGHis_tgdeep,    latco)) CALL StoreGridHistory (gtd(:,3), nGHis_tgdeep,    latco)
      IF(dogrh(nGHis_swtsfz,    latco)) CALL StoreGridHistory (gtd(:,4), nGHis_swtsfz,    latco)
      IF(dogrh(nGHis_swtrtz,    latco)) CALL StoreGridHistory (gtd(:,5), nGHis_swtrtz,    latco)
      IF(dogrh(nGHis_swtrcz,    latco)) CALL StoreGridHistory (gtd(:,6), nGHis_swtrcz,    latco)
      IF(dogrh(nGHis_mostca, latco)) CALL StoreGridHistory (gtd(:,7), nGHis_mostca, latco, 1000.0_r8)
      IF(dogrh(nGHis_mostgc, latco)) CALL StoreGridHistory (gtd(:,8), nGHis_mostgc, latco, 1000.0_r8)
      IF(dogrh(nGHis_vegtyp,    latco)) CALL StoreGridHistory (gtd(:,9), nGHis_vegtyp,    latco)
  END IF
    IF(dodia(nDiag_divgxq))CALL updia  (psint ,nDiag_divgxq,latco)
    IF(dodia(nDiag_vmoadv))CALL updia(divint,nDiag_vmoadv,latco)
    IF(dodia(nDiag_omegav))CALL updia (omg   ,nDiag_omegav,latco)
    IF(dodia(nDiag_sigdot))CALL updia(dot,nDiag_sigdot,latco)
    
  IF(isimp.ne.'YES ')THEN  
  !
  !     gplam surface pressure in mb
  !
    DO k=1,kMax
      DO i=1,ibLim
       gtmp(i,k)= gtmp(i,k)+tov(k)
       gtm(i,k) = gtm(i,k)+tov(k)
      END DO
    END DO

       gps=10*ps
       IF (microphys) THEN
          IF(nClass>0 .and. PRESENT(gvarm))THEN 
           IF ( specSfc ) THEN
           CALL DryPhysics & 
             (gzs,  gtm, gqm, gum  ,gvm  ,gps,gyu    ,gyv    ,gtd   ,&
             gqd    ,colrad,ifday ,tod   ,gtmp   ,gq     ,omg    ,jb    ,&
             lonrad,glnpm  ,cos2d ,intcosz,kMax,ibMax, &
             gTsfc, gH, gLE, gTau, gicem, gicet, gliqm, gliqt , gvarm, gvart)
           ELSE
           CALL DryPhysics & 
             (gzs,  gtm, gqm, gum  ,gvm  ,gps,gyu    ,gyv    ,gtd   ,&
             gqd    ,colrad,ifday ,tod   ,gtmp   ,gq     ,omg    ,jb    ,&
             lonrad,glnpm  ,cos2d ,intcosz,kMax,ibMax, &
             gicem, gicet, gliqm, gliqt , gvarm, gvart)
           ENDIF
          ELSE
           IF ( specSfc ) THEN
           CALL DryPhysics & 
             (gzs,  gtm, gqm, gum  ,gvm  ,gps,gyu    ,gyv    ,gtd   ,&
             gqd    ,colrad,ifday ,tod   ,gtmp   ,gq     ,omg    ,jb    ,&
             lonrad,glnpm  ,cos2d ,intcosz,kMax,ibMax, &
             gTsfc, gH, gLE, gTau, gicem, gicet, gliqm, gliqt )
           ELSE
           CALL DryPhysics & 
             (gzs,  gtm, gqm, gum  ,gvm  ,gps,gyu    ,gyv    ,gtd   ,&
             gqd    ,colrad,ifday ,tod   ,gtmp   ,gq     ,omg    ,jb    ,&
             lonrad,glnpm  ,cos2d ,intcosz,kMax,ibMax, &
             gicem, gicet, gliqm, gliqt )
           ENDIF
          END IF
       ELSE
          IF ( specSfc ) THEN
          CALL DryPhysics & 
            (gzs,  gtm, gqm, gum  ,gvm  ,gps,gyu    ,gyv    ,gtd   ,&
            gqd    ,colrad,ifday ,tod   ,gtmp   ,gq     ,omg    ,jb    ,&
            lonrad,glnpm  ,cos2d ,intcosz,kMax,ibMax,gTsfc, gH, gLE, gTau)
          ELSE
          CALL DryPhysics & 
            (gzs,  gtm, gqm, gum  ,gvm  ,gps,gyu    ,gyv    ,gtd   ,&
            gqd    ,colrad,ifday ,tod   ,gtmp   ,gq     ,omg    ,jb    ,&
            lonrad,glnpm  ,cos2d ,intcosz,kMax,ibMax)
          ENDIF
       ENDIF

    DO k=1,kMax
      DO i=1,ibLim
       gtm(i,k)=gtm(i,k)-tov(k)
       gtmp(i,k)= gtmp(i,k)-tov(k)
      END DO
    END DO
      IF(dodia(nDiag_pwater))CALL updia(dot,nDiag_pwater,latco)
 END IF      	

  END SUBROUTINE GrpComp

  !CALL Forcing_In(one,psint,divint,dot,omg,gyu,gyv,gtd, &
   !                    gqd,ibmax,kmax,ppath)
 SUBROUTINE Forcing_In(fone,fpsint,fdivint,fdot,fomg,fgyu,fgyv,fgtd, &
                       fgqd,ibmax,kmax,path)
 !Enver Ramirez Gutierrez
 IMPLICIT NONE

 !integer, parameter :: ibmax=1 ! 
 INTEGER, INTENT(IN   ) :: ibMax !Enver is not the right ibMax
 !integer, parameter :: kmax=42 
 INTEGER, INTENT(IN   ) :: kMax
 !character(len=200)  :: path    
 CHARACTER(LEN=200 ), INTENT(IN   ) :: path 
 !real               :: si(kmax+1) !REAL   , INTENT(IN   )   :: si(:)
 !real               :: delsig(kmax) 
 !real               :: ci(kmax+1) 

! REAL   , ALLOCATABLE     :: sonda (:,:)
! REAL   , ALLOCATABLE     :: sond  (:,:)
 INTEGER(KIND=i8)                  :: nlev
 REAL(KIND=r8)                     :: aux (1000,6)
 REAL(KIND=r8)                     :: Ps
 REAL(KIND=r8)                     :: Pt
 INTEGER(KIND=i8)                  :: IL
 INTEGER(KIND=i8)                  :: k
 INTEGER(KIND=i8)                  :: i
 INTEGER(KIND=i8)                  :: j
 !-----------------------------------------------------------------
 REAL(KIND=r8)                           :: fone
 REAL(KIND=r8)   , DIMENSION(ibMax,kMax) :: fpsint
 REAL(KIND=r8)   , DIMENSION(ibMax,kMax) :: fdivint
 REAL(KIND=r8)   , DIMENSION(ibMax,kMax) :: fdot
 REAL(KIND=r8),    INTENT(INOUT) :: fomg    (ibMax,kMax)
 REAL(KIND=r8),    INTENT(OUT  ) :: fgyu    (ibMax,kMax)
 REAL(KIND=r8),    INTENT(OUT  ) :: fgyv    (ibMax,kMax)
 REAL(KIND=r8),    INTENT(OUT  ) :: fgtd    (ibMax,kMax)
 REAL(KIND=r8),    INTENT(OUT  ) :: fgqd    (ibMax,kMax)
 !-----------------------------------------------------------------
! REAL   , ALLOCATABLE     :: pres (:)
! REAL   , ALLOCATABLE     :: omeg (:)
! REAL   , ALLOCATABLE     :: auvel (:)
! REAL   , ALLOCATABLE     :: avvel (:)  
! REAL   , ALLOCATABLE     :: atemp (:)
! REAL   , ALLOCATABLE     :: aumid (:)
! REAL   , ALLOCATABLE     :: sigma(:)   ! sigma = (P - Pt)/(Ps - Pt) 
! REAL   , ALLOCATABLE     :: gps  (:)   !  gps  = sigma*(Ps - Pt) + Pt
 INTEGER(KIND=i8)                  :: ierr


!Enter forcing manually (start)
!   fone    = 1.0
!   fpsint  = 0.0 
!   fdivint = 0.0 
!   fdot    = 0.0 
!   fomg(ibmax,:)    = 0.0 
!   fgyu(ibmax,:)    = 0.0
!   fgyv(ibmax,:)    = 0.0
!   fgtd(ibmax,:)    = 0.0
!   fgqd(ibmax,:)    = 0.0
!Enter forcing manually (end)

 IF(.NOT.ALLOCATED(sonda)) ALLOCATE(sonda(kmax,6) )
 sonda=0.0
 aux=0.0
 !path = "./Data"
 OPEN(1234,file=TRIM(path)//'/FORCING_IN',form='formatted',&
      access='sequential',action='read',status='old' )

   nlev=0
   DO j=1,1000
      READ(1234,*,END=10)(aux (j,i),i=1,6)
      nlev = nlev + 1
   END DO
10 CONTINUE
 REWIND(1234)
 CLOSE(1234)

   IF(.NOT.ALLOCATED(pres)) ALLOCATE(pres(nlev))
   IF(.NOT.ALLOCATED(omeg)) ALLOCATE(omeg(nlev))
   IF(.NOT.ALLOCATED(auvel)) ALLOCATE(auvel(nlev))
   IF(.NOT.ALLOCATED(avvel)) ALLOCATE(avvel(nlev)) 
   IF(.NOT.ALLOCATED(atemp)) ALLOCATE(atemp(nlev))
   IF(.NOT.ALLOCATED(aumid)) ALLOCATE(aumid(nlev))

!   ALLOCATE(p_es(nlev))    
!   ALLOCATE(p_e (nlev)) 
!   ALLOCATE(u_ur(nlev)) 
!   ALLOCATE(u_ra(nlev)) 
   
   pres(1:nlev)=aux(1:nlev,1)         !pressao mb
   omeg(1:nlev)=aux(1:nlev,2)         !vertical velocity forcing
   !DELETE
   omeg(1:nlev)=aux(1:nlev,2)         !vertical velocity forcing
   !DELETE
   
   !Its difficult receive tendencies of U and V in polar coords
   !so we exclude that option for Forcing_In
   !IF(TRIM(OPTVELO)=='DIR') THEN    
   !   !Enver verify this wind rose i'm not sure is it is working well
   !   auvel(1:nlev)=-aux(1:nlev,3)*sin((aux(1:nlev,4)*3.14)/180.0)
   !   avvel(1:nlev)=-aux(1:nlev,3)*cos((aux(1:nlev,4)*3.14)/180.0)
   !ELSE   
      auvel(1:nlev)=aux(1:nlev,3)
      avvel(1:nlev)=aux(1:nlev,4)
   !END IF          
   
   atemp(1:nlev)=aux(1:nlev,5)         !temp. absoluta tendendy C/s
   !DELETE
   !atemp(1:nlev)=aux(1:nlev,5)*0.0         !temp. absoluta tendendy C/s
   !DELETE

   !IF(TRIM(OPTUMID)=='RELATIVA') THEN
   !   u_ur(1:nlev)=aux(1:nlev,3)
   !   p_es(1:nlev)=6.112*exp(17.67*temp(1:nlev)/(243.5+temp(1:nlev)))
   !   p_e (1:nlev)=(u_ur(1:nlev)/100.0)*p_es(1:nlev)
   !   u_ra(1:nlev)=(p_e(1:nlev)*622.0/(pres(1:nlev)-p_e(1:nlev)))
   !   umid(1:nlev)=u_ra(1:nlev)/1000.0
   !ELSE
      aumid(1:nlev)=aux(1:nlev,6)         !umidade especifica. kg/kg
   !DELETE
   !   aumid(1:nlev)=aux(1:nlev,6)*0.0         !umidade especifica
   !DELETE
   !END IF
   
   IF(.NOT.ALLOCATED(sond)) ALLOCATE(sond(nlev,6))
   
   sond(1:nlev,1)=pres(1:nlev)
   sond(1:nlev,2)=omeg(1:nlev)
   sond(1:nlev,3)=auvel(1:nlev)
   sond(1:nlev,4)=avvel(1:nlev)                
   sond(1:nlev,5)=atemp(1:nlev)   ! + 273.16)*(1.0+0.61*umid(1:nlev))
   sond(1:nlev,6)=aumid(1:nlev)


   print*, 'nlev=',nlev, path
   ALLOCATE(sigma(nlev),STAT=ierr) 
   if ( ierr /= 0 ) then
     print*, 'falha na alocacao'
   endif
   sigma=1.0; sigma(1) = 1.0
   IF(.NOT.ALLOCATED(gps)) ALLOCATE(gps  (nlev))

   print*, 'start print pres'
   print*, nlev, path
   print*, (pres(i),i=1,nlev)
   print*, (sigma(i),i=1,nlev)
   print*, 'end print pres'

   Ps=aux(1,1) !Test sond(1,1)
   Pt=aux(nlev,1) !Test sond(nlev,1)
   sigma(1)=1.0
   DO j=1,nlev-1
     sigma(j+1)= (sond(j,1) - Pt)/(Ps - Pt) !Test (sond(j,1) - Pt)/(Ps - Pt) 
     gps  (j)  = sigma(j)*(Ps - Pt) + Pt   !Test sigma(j)*(Ps - Pt) + Pt
   END DO   
    
   DO IL=1,6
      CALL INTER2(sond(1:nlev,IL),nlev,si,sonda(1:kMax,IL),kMax,sigma) 
   END DO
   print*, 'after INTER2'

   fone    = 1.0
   fpsint  = 0.0 
   fdivint = 0.0 
   fdot    = 0.0 
   fomg(ibmax,:)    = sonda(:,2)
   fgyu(ibmax,:)    = sonda(:,3)
   fgyv(ibmax,:)    = sonda(:,4)
   fgtd(ibmax,:)    = sonda(:,5)
   fgqd(ibmax,:)    = sonda(:,6)
   
   !To Cancel and very big number (Enver) DELETE THIS LINES
   !erg fomg(ibmax,kmax) = 0.0!; sonda(kmax,2) = 0.0
   !erg fgyu(ibmax,kmax) = 0.0!; sonda(kmax,3) = 0.0
   !erg fgyv(ibmax,kmax) = 0.0!; sonda(kmax,4) = 0.0
   !erg fgtd(ibmax,kmax) = 0.0!; sonda(kmax,5) = 0.0
   !erg fgqd(ibmax,kmax) = 0.0!; sonda(kmax,6) = 0.0
   !To Cancel and very big number (Enver) DELETE THIS LINES AFTER CORRECTION

   print*, 'omg'; print*, (sonda(i,2),i=1,kmax); print*, '----------------'
   !print*, 'gyu'; print*, (sonda(i,3),i=1,kmax); print*, '----------------'
   !print*, 'gyv'; print*, (sonda(i,4),i=1,kmax); print*, '----------------'
   print*, 'gtd'; print*, (sonda(i,5),i=1,kmax); print*, '----------------'
   print*, 'gqd'; print*, (sonda(i,6),i=1,kmax); print*, '----------------'
   



 END SUBROUTINE FORCING_IN

!CALL specSfc_Ascii(gTsfc,gH,gLE,gTau,ppath,jdt,xxday)
! SUBROUTINE specSfc_Ascii(fgTsfc,fgH,fgLE,fgTau,fpath,jdt,xM) 
! Evitamos fgTsfc,fgH,etc porque estan definidas
! de esa forma en FieldsDynamics (o FieldsPhysics) puede parecer
! que son las mismas, pero no son pues no se carga el espacio-nombrado con USE
! al inicio de este modulo => No estan accesibles
 SUBROUTINE specSfc_Ascii(sgTsfc,sgH,sgLE,sgTau,fpath,ibMax,jdt,xM)
   !Enver Ramirez Gutierrez
    IMPLICIT NONE
    CHARACTER(LEN=200 ), INTENT(IN   )      :: fpath 
    INTEGER,             INTENT(IN   )      :: ibMax !Enver is not the right ibMax
    REAL(KIND=r8),       INTENT(INOUT)      :: sgTsfc(ibMax)
    REAL(KIND=r8),       INTENT(INOUT)      :: sgH(ibMax)
    REAL(KIND=r8),       INTENT(INOUT)      :: sgLE(ibMax)
    REAL(KIND=r8),       INTENT(INOUT)      :: sgTau(ibMax)
    INTEGER    :: jdt
    REAL(KIND=r8),SAVE                      :: specSfc_Time
    REAL(KIND=r8),SAVE                      :: specSfc_Time2

    REAL(KIND=r8)       :: xM, p1, p2

   IF( OpenFile_specSfc )THEN    !OpenFile

       OPEN(12345,file=TRIM(fpath)//'/SFC_ASCII',form='formatted',&
           access='sequential',action='read',status='old' )

         !               day, K, W/m2 , W/m2, m2/s2
         READ(12345,*) specSfc_Time, specSfc_Tsfc, specSfc_H, specSfc_LE, specSfc_Tau


       OpenFile_specSfc = .FALSE.
   ENDIF


  ! IF( MOD(tod,intfor) == 0.0 .OR. jdt == 1 )THEN !Read second profile
  IF( jdt == 1 .OR. xM .ge. specSfc_Time2 )THEN !Read second sfc data

    ! Here we have to add code to ensure that
    !        x2 > xM                 (x2 be larger than xM)
    ! such that reading or parse the next forcing profile
    ! result in
    !        x1 < xM < x2
    !

      IF( jdt == 1 )THEN
      !     READ(12345,*) itmp, nlev, x2, Ps2, tmp
        READ(12345,*) specSfc_Time2, specSfc_Tsfc2, specSfc_H2, specSfc_LE2, specSfc_Tau2
      ELSE
       !aux = aux2
       !...
       specSfc_Time = specSfc_Time2
       specSfc_Tsfc = specSfc_Tsfc2
       specSfc_H    = specSfc_H2
       specSfc_LE   = specSfc_LE2
       specSfc_Tau  = specSfc_Tau2
     
      ! READ(1234,*) itmp, nlev, x2, Ps2, tmp
       READ(12345,*) specSfc_Time2, specSfc_Tsfc2, specSfc_H2, specSfc_LE2, specSfc_Tau2
     ENDIF

  ENDIF !Read second profile


  !Computing weights by linear interpolation
   p1 = (specSfc_Time2 - xM)/(specSfc_Time2 - specSfc_Time)
   p2 = 1.0 - p1

  !  specSfc_Ascii(fgTsfc,fgH,fgLE,fgTau,fpath,jdt,xM)
    sgTsfc(ibMax)  =  ( p1*specSfc_Tsfc +  p2*specSfc_Tsfc2 )
    sgH(ibMax)     =  ( p1*specSfc_H    +  p2*specSfc_H2 )
    sgLE(ibMax)    =  ( p1*specSfc_LE   +  p2*specSfc_LE2 ) 
    sgTau(ibMax)   =  ( p1*specSfc_Tau  +  p2*specSfc_Tau2 )
  
  !    IF(dodia(nDiag_TEMPFORC))CALL updia  (fgtd ,nDiag_TEMPFORC,1)
  !    IF(dodia(nDiag_MOISFORC))CALL updia  (fgqd ,nDiag_MOISFORC,1)

 END SUBROUTINE specSfc_Ascii



  !CALL Forcings_Ascii(one,psint,divint,dot,omg,gyu,gyv,gtd, &
   !                    gqd,glnpm,ibmax,kmax,ppath)
 SUBROUTINE Forcings_Ascii(fone,fpsint,fdivint,fdot,fomg,fgyu,fgyv,fgtd, &
                       fgqd,fglnpm,ibmax,kmax,path,tod,jdt,xM,colrad)
 !Enver Ramirez Gutierrez
 IMPLICIT NONE

 INTEGER, INTENT(IN   ) :: ibMax !Enver is not the right ibMax
 INTEGER, INTENT(IN   ) :: kMax
 CHARACTER(LEN=200 ), INTENT(IN   ) :: path 

! REAL   , ALLOCATABLE     :: sonda (:,:)
! REAL   , ALLOCATABLE     :: sond  (:,:)
! REAL   , ALLOCATABLE     :: sonda2 (:,:)
! REAL   , ALLOCATABLE     :: sond2  (:,:)
 INTEGER(KIND=i8)                  :: nlev
 REAL(KIND=r8)                     :: aux (1000,6)
 REAL(KIND=r8)                     :: aux2 (1000,6)
 !REAL(KIND=r8)                     :: Ps1, Ps2
 REAL(KIND=r8)                     :: Pt
 INTEGER(KIND=i8)                  :: IL
 INTEGER(KIND=i8)                  :: k
 INTEGER(KIND=i8)                  :: i
 INTEGER(KIND=i8)                  :: j
 !-----------------------------------------------------------------
 REAL(KIND=r8)                           :: fone
 REAL(KIND=r8)   , DIMENSION(ibMax,kMax) :: fpsint
 REAL(KIND=r8)   , DIMENSION(ibMax,kMax) :: fdivint
 REAL(KIND=r8)   , DIMENSION(ibMax,kMax) :: fdot
 REAL(KIND=r8),    INTENT(INOUT) :: fomg    (ibMax,kMax)
 REAL(KIND=r8),    INTENT(OUT  ) :: fgyu    (ibMax,kMax)
 REAL(KIND=r8),    INTENT(OUT  ) :: fgyv    (ibMax,kMax)
 REAL(KIND=r8),    INTENT(OUT  ) :: fgtd    (ibMax,kMax)
 REAL(KIND=r8),    INTENT(OUT  ) :: fgqd    (ibMax,kMax)
 REAL(KIND=r8),    INTENT(OUT  ) :: fglnpm    (ibMax,kMax)
 !-----------------------------------------------------------------
! REAL   , ALLOCATABLE     :: pres (:)
! REAL   , ALLOCATABLE     :: omeg (:) 
! REAL   , ALLOCATABLE     :: auvel (:)
! REAL   , ALLOCATABLE     :: avvel (:)
! REAL   , ALLOCATABLE     :: atemp (:)
! REAL   , ALLOCATABLE     :: aumid (:)
! REAL   , ALLOCATABLE     :: sigma(:)   ! sigma = (P - Pt)/(Ps - Pt) 
! REAL   , ALLOCATABLE     :: gps  (:)   !  gps  = sigma*(Ps - Pt) + Pt
 REAL(KIND=r8)   , INTENT(IN   ) :: tod
 INTEGER(KIND=i8)                  :: ierr
 REAL(KIND=r8),SAVE       :: x1, x2
 REAL(KIND=r8)       :: tmp, xM, p1, p2
 REAL(KIND=r8),    INTENT(IN   ) :: colrad (ibMax)
 INTEGER(KIND=i8)    :: itmp
 INTEGER    :: jdt
! LOGICAL, SAVE    :: OpenFile =.TRUE.
 LOGICAL    :: DynTend  =.FALSE.  ! .FALSE. !Associado somente com a tendencia
                                 !         !devido ao momento zonal e meridional
                                 !         !Se (DynTend=.FALSE.) => auvel=avvel=0.0 

 aux=0.0; aux2=0.0

 IF( OpenFile )THEN    !OpenFile

     OPEN(1234,file=TRIM(path)//'/FORCINGS_ASCII',form='formatted',&
         access='sequential',action='read',status='old' )

     READ(1234,*) itmp, nlev, x1, Ps1, tmp
    ! change days since to model days 
     x1 = x1 + 1
     DO j=1,nlev
       READ(1234,*)(aux (j,i),i=1,6)
     END DO

     ALLOCATE(pres(nlev))
     ALLOCATE(omeg(nlev)) 
     ALLOCATE(auvel(nlev))
     ALLOCATE(avvel(nlev))
     ALLOCATE(atemp(nlev))
     ALLOCATE(aumid(nlev))

     ALLOCATE(sond(nlev,6)); ALLOCATE(sond2(nlev,6))
     ALLOCATE(sigma(nlev)) 
     ALLOCATE(gps  (nlev))

     ALLOCATE(sonda(kmax,6)); ALLOCATE(sonda2(kmax,6))
     sonda=0.0; sonda2=0.0

      pres(1:nlev)=aux(1:nlev,1)         !pressao mb
      omeg(1:nlev)=aux(1:nlev,2)         !vertical velocity forcing
      IF(DynTend)THEN
       auvel(1:nlev)=aux(1:nlev,3)        ! u tendency
       avvel(1:nlev)=aux(1:nlev,4)        ! v tendency
      ELSE
       auvel(1:nlev)=aux(1:nlev,3)*0.0    !0.0 u tendency
       avvel(1:nlev)=aux(1:nlev,4)*0.0    !0.0 v tendency
      ENDIF
      atemp(1:nlev)=aux(1:nlev,5)        !temp. absoluta tendendy C/s
      aumid(1:nlev)=aux(1:nlev,6)         !umidade especifica. kg/kg

      sond(1:nlev,1)=pres(1:nlev)
      sond(1:nlev,2)=omeg(1:nlev)
      sond(1:nlev,3)=auvel(1:nlev)
      sond(1:nlev,4)=avvel(1:nlev)                
      sond(1:nlev,5)=atemp(1:nlev)   ! + 273.16)*(1.0+0.61*umid(1:nlev))
      sond(1:nlev,6)=aumid(1:nlev)
   
      sigma=1.0; sigma(1) = 1.0

      !Enver now Ps is read from the header of the forcing file (start)
      !   Ps=aux(1,1) !Test sond(1,1)
      Pt=aux(nlev,1) !Test sond(nlev,1)
      sigma(1)=1.0
      DO j=1,nlev-1
        sigma(j+1)= (sond(j,1) - Pt)/(Ps1 - Pt) !Test (sond(j,1) - Pt)/(Ps - Pt) 
        gps  (j)  = sigma(j)*(Ps1 - Pt) + Pt   !Test sigma(j)*(Ps - Pt) + Pt
      END DO   
    print*, 'niveis sigma'
    print*, sigma, si
    print*, 'niveis sigma'
      DO IL=1,6
         CALL INTER2(sond(1:nlev,IL),nlev,si,sonda(1:kMax,IL),kMax,sigma) 
      END DO
 
      print*, 'sonda',sonda(1:KMax,5),sonda(1:KMax,6) 
     OpenFile = .FALSE.
 ENDIF                 !OpenFile

! PRINT*,  'inside GridDynamics intfor=', intfor
!  PRINT*,  'MOD(tod,intfor)=',MOD(tod,intfor)
!  PRINT*,  'jdt=',jdt


! IF( MOD(tod,intfor) == 0.0 .OR. jdt == 1 )THEN !Read second profile
IF( jdt == 1 .OR. xM .ge. x2 )THEN !Read second profile !Enver 01Jul2013 to correct
                                                        !         a bug in time
                                                        !         interpolation.

  ! 28dec2014
  ! Here we have to add code to ensure that
  !        x2 > xM                 (x2 be larger than xM)
  ! such that reading or parse the next forcing profile
  ! result in
  !        x1 < xM < x2
  !
  !  28dec2014

   IF( jdt == 1 )THEN
     READ(1234,*) itmp, nlev, x2, Ps2, tmp
     ! change days since to model days 
     x2 = x2 + 1
     DO j=1,nlev
        READ(1234,*)(aux2 (j,i),i=1,6)
     END DO
   ELSE
     !aux = aux2
     !...
     sonda(:,2)    = sonda2(:,2)
     sonda(:,3)    = sonda2(:,3)
     sonda(:,4)    = sonda2(:,4)
     sonda(:,5)    = sonda2(:,5)
     sonda(:,6)    = sonda2(:,6)
     x1            = x2  !Enver 28Jun2013
     Ps1           = Ps2
     
     READ(1234,*) itmp, nlev, x2, Ps2, tmp
     ! change days since to model days 
     x2 = x2 + 1
     DO j=1,nlev
        READ(1234,*)(aux2 (j,i),i=1,6)
     END DO
   ENDIF
     print*, 'Ps1,Ps2', Ps1, Ps2
    !  print*, 'NLEV=',nlev, size(pres)
      pres(1:nlev)=aux2(1:nlev,1)         !pressao mb
      omeg(1:nlev)=aux2(1:nlev,2)         !vertical velocity forcing


      IF(DynTend)THEN
         auvel(1:nlev)=aux2(1:nlev,3)        ! u tendency
         avvel(1:nlev)=aux2(1:nlev,4)        ! v tendency
      ELSE
         auvel(1:nlev)=aux2(1:nlev,3)*0.0    !0.0 u tendency
         avvel(1:nlev)=aux2(1:nlev,4)*0.0    !0.0 v tendency
      ENDIF
      atemp(1:nlev)=aux2(1:nlev,5)        !temp. absoluta tendendy C/s
      aumid(1:nlev)=aux2(1:nlev,6)         !umidade especifica. kg/kg

      sond2(1:nlev,1)=pres(1:nlev)
      sond2(1:nlev,2)=omeg(1:nlev)
      sond2(1:nlev,3)=auvel(1:nlev)
      sond2(1:nlev,4)=avvel(1:nlev)                
      sond2(1:nlev,5)=atemp(1:nlev)   ! + 273.16)*(1.0+0.61*umid(1:nlev))
      sond2(1:nlev,6)=aumid(1:nlev)
   
      sigma=1.0; sigma(1) = 1.0

      !Enver now Ps is read from the header of the forcing file (start)
      !   Ps=aux2(1,1) !Test sond(1,1)
      Pt=aux2(nlev,1) !Test sond(nlev,1)
      sigma(1)=1.0
      DO j=1,nlev-1
        sigma(j+1)= (sond2(j,1) - Pt)/(Ps2 - Pt) !Test (sond(j,1) - Pt)/(Ps - Pt) 
        gps  (j)  = sigma(j)*(Ps2 - Pt) + Pt   !Test sigma(j)*(Ps - Pt) + Pt
      END DO   
      DO IL=1,6
         CALL INTER2(sond2(1:nlev,IL),nlev,si,sonda2(1:kMax,IL),kMax,sigma) 
      END DO
      print*, 'sonda2',sonda2(1:KMax,5),sonda2(1:KMax,6) 
 ENDIF !Read second profile


 !Computing weights by linear interpolation
  p1 = (x2 - xM)/(x2 - x1)
  p2 = 1.0 - p1

  print 123, 'x1 xM x2', x1, xM, x2, MOD(tod,intfor), p1, p2; 123 FORMAT(A,4F11.4,2F7.2)

  fone    = 1.0
  fpsint  = 0.0 
  fdivint = 0.0 
  fdot    = 0.0 
  fomg(ibmax,:)    =  ( p1*sonda(:,2) + p2*sonda2(:,2) ) / forcings_weight_d
  fgyu(ibmax,:)    =  ( p1*sonda(:,3) + p2*sonda2(:,3) ) * (SIN(colrad (ibmax))) / forcings_weight_d!sonda(:,3)
  fgyv(ibmax,:)    =  ( p1*sonda(:,4) + p2*sonda2(:,4) ) * (SIN(colrad (ibmax))) / forcings_weight_d!sonda(:,4)
  fgtd(ibmax,:)    =  ( p1*sonda(:,5) + p2*sonda2(:,5) ) / forcings_weight_t !sonda(:,5)
  fgqd(ibmax,:)    =  ( p1*sonda(:,6) + p2*sonda2(:,6) ) / forcings_weight_m !sonda(:,6)
  !teste
   !print*,  'teste=', p1*log(Ps1/10.0), p2*log(Ps2/10.0), ( p1*log(Ps1/10.0) + p2*log(Ps2/10.0) ), Ps1, Ps2
   !fglnpm(ibmax,:) = 4.6_r8
   print*,  'teste=', p1*log(Ps1/10.0), p2*log(Ps2/10.0), Ps1, Ps2, sonda(1,5), sonda2(1,5),shape(fglnpm)
  fglnpm(ibmax,1)  =  ( p1*log(Ps1/10.0) + p2*log(Ps2/10.0) )
  !teste
  
    IF(dodia(nDiag_TEMPFORC))CALL updia  (fgtd ,nDiag_TEMPFORC,1)
    IF(dodia(nDiag_MOISFORC))CALL updia  (fgqd ,nDiag_MOISFORC,1)


  !To Cancel and very big number (Enver) DELETE THIS LINES
  !erg fomg(ibmax,kmax) = 0.0!; sonda(kmax,2) = 0.0
  !erg fgyu(ibmax,kmax) = 0.0!; sonda(kmax,3) = 0.0
  !erg fgyv(ibmax,kmax) = 0.0!; sonda(kmax,4) = 0.0
  !erg fgtd(ibmax,kmax) = 0.0!; sonda(kmax,5) = 0.0
  !erg fgqd(ibmax,kmax) = 0.0!; sonda(kmax,6) = 0.0
  !To Cancel and very big number (Enver) DELETE THIS LINES AFTER CORRECTION

 END SUBROUTINE Forcings_Ascii

 SUBROUTINE Forcings_Ascii_lsf(fone,fpsint,fdivint,fdot,fomg,fgum,fgvm,fgtd, &
                       fgqd,ibmax,kmax,path,tod,jdt,xM,colrad)
 !Enver Ramirez Gutierrez
 IMPLICIT NONE

 INTEGER, INTENT(IN   ) :: ibMax !Enver is not the right ibMax
 INTEGER, INTENT(IN   ) :: kMax
 CHARACTER(LEN=200 ), INTENT(IN   ) :: path 

! REAL   , ALLOCATABLE     :: sonda (:,:)
! REAL   , ALLOCATABLE     :: sond  (:,:)
! REAL   , ALLOCATABLE     :: sonda2 (:,:)
! REAL   , ALLOCATABLE     :: sond2  (:,:)
 INTEGER(KIND=i8)                  :: nlev
 REAL(KIND=r8)                     :: aux (1000,6)
 REAL(KIND=r8)                     :: aux2 (1000,6)
 REAL(KIND=r8)                     :: Ps
 REAL(KIND=r8)                     :: Pt
 INTEGER(KIND=i8)                  :: IL
 INTEGER(KIND=i8)                  :: k
 INTEGER(KIND=i8)                  :: i
 INTEGER(KIND=i8)                  :: j
 !-----------------------------------------------------------------
 REAL(KIND=r8)                           :: fone
 REAL(KIND=r8)   , DIMENSION(ibMax,kMax) :: fpsint
 REAL(KIND=r8)   , DIMENSION(ibMax,kMax) :: fdivint
 REAL(KIND=r8)   , DIMENSION(ibMax,kMax) :: fdot
 REAL(KIND=r8),    INTENT(INOUT) :: fomg    (ibMax,kMax)
 REAL(KIND=r8),    INTENT(OUT  ) :: fgum    (ibMax,kMax)
 REAL(KIND=r8),    INTENT(OUT  ) :: fgvm    (ibMax,kMax)
 REAL(KIND=r8),    INTENT(OUT  ) :: fgtd    (ibMax,kMax)
 REAL(KIND=r8),    INTENT(OUT  ) :: fgqd    (ibMax,kMax)
 !-----------------------------------------------------------------
! REAL   , ALLOCATABLE     :: pres (:)
! REAL   , ALLOCATABLE     :: omeg (:) 
! REAL   , ALLOCATABLE     :: auvel (:)
! REAL   , ALLOCATABLE     :: avvel (:)
! REAL   , ALLOCATABLE     :: atemp (:)
! REAL   , ALLOCATABLE     :: aumid (:)
! REAL   , ALLOCATABLE     :: sigma(:)   ! sigma = (P - Pt)/(Ps - Pt) 
! REAL   , ALLOCATABLE     :: gps  (:)   !  gps  = sigma*(Ps - Pt) + Pt
 REAL(KIND=r8)   , INTENT(IN   ) :: tod
 INTEGER(KIND=i8)                  :: ierr
 REAL(KIND=r8),SAVE       :: x1, x2
 REAL(KIND=r8)       :: tmp, xM, p1, p2
 REAL(KIND=r8),    INTENT(IN   ) :: colrad (ibMax)
 INTEGER(KIND=i8)    :: itmp
 INTEGER    :: jdt
! LOGICAL, SAVE    :: OpenFile =.TRUE.
! LOGICAL    :: DynTend  =.TRUE.  ! .FALSE. !Associado somente com a tendencia
!                                 !         !devido ao momento zonal e meridional
!                                 !         !Se (DynTend=.FALSE.) => auvel=avvel=0.0 

 aux=0.0; aux2=0.0

 IF( OpenFile )THEN    !OpenFile

     OPEN(1234,file=TRIM(path)//'/FORCINGS_ASCII_lsf',form='formatted',&
         access='sequential',action='read',status='old' )

     READ(1234,*) itmp, nlev, x1, Ps, tmp
     DO j=1,nlev
       READ(1234,*)(aux (j,i),i=1,6)
     END DO

     ALLOCATE(pres(nlev))
     ALLOCATE(omeg(nlev)) 
     ALLOCATE(auvel(nlev))
     ALLOCATE(avvel(nlev))
     ALLOCATE(atemp(nlev))
     ALLOCATE(aumid(nlev))

     ALLOCATE(sond(nlev,6)); ALLOCATE(sond2(nlev,6))
     ALLOCATE(sigma(nlev)) 
     ALLOCATE(gps  (nlev))

     ALLOCATE(sonda(kmax,6)); ALLOCATE(sonda2(kmax,6))
     sonda=0.0; sonda2=0.0

      pres(1:nlev)=aux(1:nlev,1)         !pressao mb
      omeg(1:nlev)=aux(1:nlev,2)         !vertical velocity forcing
!      IF(DynTend)THEN
       auvel(1:nlev)=aux(1:nlev,3)        ! u tendency
       avvel(1:nlev)=aux(1:nlev,4)        ! v tendency
!      ELSE
!       auvel(1:nlev)=aux(1:nlev,3)*0.0    !0.0 u tendency
!       avvel(1:nlev)=aux(1:nlev,4)*0.0    !0.0 v tendency
!      ENDIF
      atemp(1:nlev)=aux(1:nlev,5)        !temp. absoluta tendendy C/s
      aumid(1:nlev)=aux(1:nlev,6)         !umidade especifica. kg/kg

      sond(1:nlev,1)=pres(1:nlev)
      sond(1:nlev,2)=omeg(1:nlev)
      sond(1:nlev,3)=auvel(1:nlev)
      sond(1:nlev,4)=avvel(1:nlev)                
      sond(1:nlev,5)=atemp(1:nlev)   ! + 273.16)*(1.0+0.61*umid(1:nlev))
      sond(1:nlev,6)=aumid(1:nlev)
   
      sigma=1.0; sigma(1) = 1.0

      !Enver now Ps is read from the header of the forcing file (start)
      !   Ps=aux(1,1) !Test sond(1,1)
      Pt=aux(nlev,1) !Test sond(nlev,1)
      sigma(1)=1.0
      DO j=1,nlev-1
        sigma(j+1)= (sond(j,1) - Pt)/(Ps - Pt) !Test (sond(j,1) - Pt)/(Ps - Pt) 
        gps  (j)  = sigma(j)*(Ps - Pt) + Pt   !Test sigma(j)*(Ps - Pt) + Pt
      END DO   
      print*, 'sond',sond(1:nlev,3),sond(1:nlev,4) 
     
      DO IL=1,6
         CALL INTER2(sond(1:nlev,IL),nlev,si,sonda(1:kMax,IL),kMax,sigma) 
      END DO
 
      !print*, 'sond',sond(1:nlev,3),sond(1:nlev,4) 
      print*, 'sonda',sonda(1:KMax,3),sonda(1:KMax,4) 
     OpenFile = .FALSE.
 ENDIF                 !OpenFile

! PRINT*,  'inside GridDynamics intfor=', intfor
!  PRINT*,  'MOD(tod,intfor)=',MOD(tod,intfor)
!  PRINT*,  'jdt=',jdt


! IF( MOD(tod,intfor) == 0.0 .OR. jdt == 1 )THEN !Read second profile
IF( jdt == 1 .OR. xM .ge. x2 )THEN !Read second profile !Enver 01Jul2013 to correct
                                                        !         a bug in time
                                                        !         interpolation.

  ! 28dec2014
  ! Here we have to add code to ensure that
  !        x2 > xM                 (x2 be larger than xM)
  ! such that reading or parse the next forcing profile
  ! result in
  !        x1 < xM < x2
  !
  !  28dec2014

   IF( jdt == 1 )THEN
     print*, 'forcings_weight_m', forcings_weight_m
     READ(1234,*) itmp, nlev, x2, Ps, tmp
      !print*, 'aux2'
     DO j=1,nlev
        READ(1234,*)(aux2 (j,i),i=1,6)
        !print*, (aux2 (j,i),i=1,6) 
     END DO
        !stop
   ELSE
     !aux = aux2
     !...
     sonda(:,2)    = sonda2(:,2)
     sonda(:,3)    = sonda2(:,3)
     sonda(:,4)    = sonda2(:,4)
     sonda(:,5)    = sonda2(:,5)
     sonda(:,6)    = sonda2(:,6)
     x1            = x2  !Enver 28Jun2013
     READ(1234,*) itmp, nlev, x2, Ps, tmp
     DO j=1,nlev
        READ(1234,*)(aux2 (j,i),i=1,6)
     END DO
   ENDIF
      print*, 'NLEV=',nlev, size(pres)
      pres(1:nlev)=aux2(1:nlev,1)         !pressao mb
      omeg(1:nlev)=aux2(1:nlev,2)         !vertical velocity forcing


!      IF(DynTend)THEN
         auvel(1:nlev)=aux2(1:nlev,3)        ! u tendency
         avvel(1:nlev)=aux2(1:nlev,4)        ! v tendency
!      ELSE
!         auvel(1:nlev)=aux(1:nlev,3)*0.0    !0.0 u tendency
!         avvel(1:nlev)=aux(1:nlev,4)*0.0    !0.0 v tendency
!      ENDIF
      atemp(1:nlev)=aux2(1:nlev,5)        !temp. absoluta tendendy C/s
      aumid(1:nlev)=aux2(1:nlev,6)         !umidade especifica. kg/kg

      sond2(1:nlev,1)=pres(1:nlev)
      sond2(1:nlev,2)=omeg(1:nlev)
      sond2(1:nlev,3)=auvel(1:nlev)
      sond2(1:nlev,4)=avvel(1:nlev)                
      sond2(1:nlev,5)=atemp(1:nlev)   ! + 273.16)*(1.0+0.61*umid(1:nlev))
      sond2(1:nlev,6)=aumid(1:nlev)
   
      sigma=1.0; sigma(1) = 1.0

      !Enver now Ps is read from the header of the forcing file (start)
      !   Ps=aux2(1,1) !Test sond(1,1)
      Pt=aux2(nlev,1) !Test sond(nlev,1)
      sigma(1)=1.0
      DO j=1,nlev-1
        sigma(j+1)= (sond2(j,1) - Pt)/(Ps - Pt) !Test (sond(j,1) - Pt)/(Ps - Pt) 
        gps  (j)  = sigma(j)*(Ps - Pt) + Pt   !Test sigma(j)*(Ps - Pt) + Pt
      END DO   
      print*, 'sond',sond2(1:nlev,3),sond2(1:nlev,4) 
      DO IL=1,6
         CALL INTER2(sond2(1:nlev,IL),nlev,si,sonda2(1:kMax,IL),kMax,sigma) 
      END DO
      !print*, 'sond',sond2(1:nlev,3),sond2(1:nlev,4) 
     ! print*, 'sonda2',sonda2(1:KMax,3),sonda2(1:KMax,4) 
 ENDIF !Read second profile


 !Computing weights by linear interpolation
  p1 = (x2 - xM)/(x2 - x1)
  p2 = 1.0 - p1

  print 123, 'x1 xM x2', x1, xM, x2, MOD(tod,intfor); 123 FORMAT(A,4F11.4)

  fone    = 1.0
  fpsint  = 0.0 
  fdivint = 0.0 
  fdot    = 0.0 
  fomg(ibmax,:)    =  ( p1*sonda(:,2) + p2*sonda2(:,2) ) / forcings_weight_d
  fgum(ibmax,:)    =  ( p1*sonda(:,3) + p2*sonda2(:,3) ) * (SIN(colrad (ibmax))) / forcings_weight_d!sonda(:,3)
  fgvm(ibmax,:)    =  ( p1*sonda(:,4) + p2*sonda2(:,4) ) * (SIN(colrad (ibmax))) / forcings_weight_d!sonda(:,4)
  fgtd(ibmax,:)    =  ( p1*sonda(:,5) + p2*sonda2(:,5) ) / forcings_weight_t !sonda(:,5)
  fgqd(ibmax,:)    =  ( p1*sonda(:,6) + p2*sonda2(:,6) ) / forcings_weight_m !sonda(:,6)
  
    IF(dodia(nDiag_TEMPFORC))CALL updia  (fgtd ,nDiag_TEMPFORC,1)
    IF(dodia(nDiag_MOISFORC))CALL updia  (fgqd ,nDiag_MOISFORC,1)


  !To Cancel and very big number (Enver) DELETE THIS LINES
  !erg fomg(ibmax,kmax) = 0.0!; sonda(kmax,2) = 0.0
  !erg fgyu(ibmax,kmax) = 0.0!; sonda(kmax,3) = 0.0
  !erg fgyv(ibmax,kmax) = 0.0!; sonda(kmax,4) = 0.0
  !erg fgtd(ibmax,kmax) = 0.0!; sonda(kmax,5) = 0.0
  !erg fgqd(ibmax,kmax) = 0.0!; sonda(kmax,6) = 0.0
  !To Cancel and very big number (Enver) DELETE THIS LINES AFTER CORRECTION

 END SUBROUTINE Forcings_Ascii_lsf

 !Nudging_Ascii(guo,gvo,gto,gqo,ibmax,kmax,ppath,tod,jdt,xxday,colrad)

 SUBROUTINE Nudging_Ascii(fguo,fgvo,fgto, &
                       fgqo,ibmax,kmax,path,tod,jdt,xM,colrad)
 !Enver Ramirez Gutierrez
 IMPLICIT NONE

 INTEGER, INTENT(IN   ) :: ibMax !Enver is not the right ibMax
 INTEGER, INTENT(IN   ) :: kMax
 CHARACTER(LEN=200 ), INTENT(IN   ) :: path 

 INTEGER(KIND=i8)                  :: nlev
 INTEGER(KIND=i8), PARAMETER       :: nvars=5
 REAL(KIND=r8)                     :: aux (1000,nvars)
 REAL(KIND=r8)                     :: aux2 (1000,nvars)
 REAL(KIND=r8)                     :: Ps
 REAL(KIND=r8)                     :: Pt
 INTEGER(KIND=i8)                  :: IL
 INTEGER(KIND=i8)                  :: k
 INTEGER(KIND=i8)                  :: i
 INTEGER(KIND=i8)                  :: j
 !-----------------------------------------------------------------
 REAL(KIND=r8),    INTENT(OUT  ) :: fguo    (ibMax,kMax)
 REAL(KIND=r8),    INTENT(OUT  ) :: fgvo    (ibMax,kMax)
 REAL(KIND=r8),    INTENT(OUT  ) :: fgto    (ibMax,kMax)
 REAL(KIND=r8),    INTENT(OUT  ) :: fgqo    (ibMax,kMax)
 !-----------------------------------------------------------------
 REAL(KIND=r8)   , INTENT(IN   ) :: tod
 INTEGER(KIND=i8)                  :: ierr
 REAL(KIND=r8),SAVE       :: nx1, nx2
 REAL(KIND=r8)       :: tmp, xM, p1, p2
 REAL(KIND=r8),    INTENT(IN   ) :: colrad (ibMax)
 INTEGER(KIND=i8)    :: itmp
 INTEGER    :: jdt
! LOGICAL, SAVE    :: OpenFile =.TRUE.
 LOGICAL    :: DynTend  =.TRUE.  ! .FALSE. !Associado somente com a tendencia
                                 !         !devido ao momento zonal e meridional
                                 !         !Se (DynTend=.FALSE.) => auvel=avvel=0.0 

 aux=0.0; aux2=0.0

 IF( OpenFile_ndg )THEN    !OpenFile

     OPEN(234,file=TRIM(path)//'/NUDGING_ASCII',form='formatted',&
         access='sequential',action='read',status='old' )

     OPEN(2704201,file=TRIM(path)//'../dataout/Nudging_Ascii_temp.out',form='formatted',&
         access='sequential')

     OPEN(2704202,file=TRIM(path)//'../dataout/Nudging_Ascii_q.out',form='formatted',&
         access='sequential')

     READ(234,*) itmp, nlev, nx1, Ps1, tmp
     DO j=1,nlev
       READ(234,*)(aux (j,i),i=1,nvars)
     END DO


     ALLOCATE(sondg(nlev,nvars)); ALLOCATE(sondg2(nlev,nvars))
     ALLOCATE(sigma_ndg(nlev)) 
     !ALLOCATE(gps  (nlev))

     ALLOCATE(sondga(kmax,nvars)); ALLOCATE(sondga2(kmax,nvars))
     sondga=0.0; sondga2=0.0

      sondg(1:nlev,1)=aux(1:nlev,1) !p
      sondg(1:nlev,2)=aux(1:nlev,2) !u
      sondg(1:nlev,3)=aux(1:nlev,3) !v
      sondg(1:nlev,4)=aux(1:nlev,4) !t
      sondg(1:nlev,5)=aux(1:nlev,5) !q
   
      sigma_ndg=1.0; sigma_ndg(1) = 1.0

      !Enver now Ps is read from the header of the forcing file (start)
      !   Ps=aux(1,1) !Test sond(1,1)
      Pt=aux(nlev,1) !Test sond(nlev,1)
      sigma_ndg(1)=1.0
      DO j=1,nlev-1
        sigma_ndg(j+1)= (sondg(j,1) - Pt)/(Ps1 - Pt) !Test (sond(j,1) - Pt)/(Ps - Pt) 
        !gps  (j)  = sigma_ndg(j)*(Ps - Pt) + Pt   !Test sigma(j)*(Ps - Pt) + Pt
      END DO   
!    print*, 'niveis sigma_ndg'
!    print*, sigma_ndg, si
!    print*, 'niveis sigma_ndg'
      DO IL=1,nvars
         CALL INTER2(sondg(1:nlev,IL),nlev,si,sondga(1:kMax,IL),kMax,sigma_ndg) 
      END DO
 
!      print*, 'sonda',sondga(1:KMax,5),sondga(1:KMax,6) 
     OpenFile_ndg = .FALSE.
 ENDIF                 !OpenFile


! IF( MOD(tod,intfor) == 0.0 .OR. jdt == 1 )THEN !Read second profile
IF( jdt == 1 .OR. xM .ge. nx2 )THEN !Read second profile !Enver 01Jul2013 to correct
                                                        !         a bug in time
                                                        !         interpolation.

  ! 28dec2014
  ! Here we have to add code to ensure that
  !        x2 > xM                 (x2 be larger than xM)
  ! such that reading or parse the next forcing profile
  ! result in
  !        x1 < xM < x2
  !
  !  28dec2014

   IF( jdt == 1 )THEN
     READ(234,*) itmp, nlev, nx2, Ps2, tmp
     DO j=1,nlev
        READ(234,*)(aux2 (j,i),i=1,nvars)
     END DO
   ELSE
     !aux = aux2
     !...
     sondga(:,2)    = sondga2(:,2)
     sondga(:,3)    = sondga2(:,3)
     sondga(:,4)    = sondga2(:,4)
     sondga(:,5)    = sondga2(:,5)
     nx1            = nx2  !Enver 28Jun2013
     Ps1           = Ps2
     READ(234,*) itmp, nlev, nx2, Ps2, tmp
     DO j=1,nlev
        READ(234,*)(aux2 (j,i),i=1,nvars)
     END DO
   ENDIF
!      print*, 'NLEV=',nlev, size(pres)
!      pres(1:nlev)=aux2(1:nlev,1)         !pressao mb
!      omeg(1:nlev)=aux2(1:nlev,2)         !vertical velocity forcing
!
!
!      IF(DynTend)THEN
!         auvel(1:nlev)=aux2(1:nlev,3)        ! u tendency
!         avvel(1:nlev)=aux2(1:nlev,4)        ! v tendency
!      ELSE
!         auvel(1:nlev)=aux2(1:nlev,3)*0.0    !0.0 u tendency
!         avvel(1:nlev)=aux2(1:nlev,4)*0.0    !0.0 v tendency
!      ENDIF
!      atemp(1:nlev)=aux2(1:nlev,5)        !temp. absoluta tendendy C/s
!      aumid(1:nlev)=aux2(1:nlev,6)         !umidade especifica. kg/kg

      sondg2(1:nlev,1)=aux2(1:nlev,1) !p
      sondg2(1:nlev,2)=aux2(1:nlev,2) !u
      sondg2(1:nlev,3)=aux2(1:nlev,3) !v
      sondg2(1:nlev,4)=aux2(1:nlev,4) !t                
      sondg2(1:nlev,5)=aux2(1:nlev,5) !q
   
      sigma_ndg=1.0; sigma_ndg(1) = 1.0

      !Enver now Ps is read from the header of the forcing file (start)
      !   Ps=aux2(1,1) !Test sond(1,1)
      Pt=aux2(nlev,1) !Test sond(nlev,1)
      sigma_ndg(1)=1.0
      DO j=1,nlev-1
        sigma_ndg(j+1)= (sondg2(j,1) - Pt)/(Ps2 - Pt) !Test (sond(j,1) - Pt)/(Ps - Pt) 
        !gps  (j)  = sigma_ndg(j)*(Ps - Pt) + Pt   !Test sigma(j)*(Ps - Pt) + Pt
      END DO   
      DO IL=1,nvars
         CALL INTER2(sondg2(1:nlev,IL),nlev,si,sondga2(1:kMax,IL),kMax,sigma_ndg) 
      END DO
 ENDIF !Read second profile


 !Computing weights by linear interpolation
  p1 = (nx2 - xM)/(nx2 - nx1)
  p2 = 1.0 - p1

  print 123, 'nx1 yM nx2', nx1, xM, nx2, MOD(tod,intfor); 123 FORMAT(A,4F11.4)
  ! print *, 'shape(tov)',shape(tov),'',shape(sondga),shape(fgvo),shape(fgto),shape(fgqo)

  fguo(ibmax,:)    =  ( p1*sondga(:,2) + p2*sondga2(:,2) ) * SIN(colrad (ibmax))
  fgvo(ibmax,:)    =  ( p1*sondga(:,3) + p2*sondga2(:,3) ) * SIN(colrad (ibmax))
  fgto(ibmax,:)    =  ( p1*sondga(:,4) + p2*sondga2(:,4) ) - tov
  fgqo(ibmax,:)    =  ( p1*sondga(:,5) + p2*sondga2(:,5) ) 

  !IF( jdt == 1 .OR. xM .ge. x2 )THEN   
!       write(2704201,*) fgto(ibmax,:)
!       write(2704202,*) fgqo(ibmax,:)
  !ENDIF

 END SUBROUTINE Nudging_Ascii

 SUBROUTINE Geostrophy_Ascii(fgug,fgvg, &
                       ibmax,kmax,path,tod,jdt,xM,colrad)
 !Enver Ramirez Gutierrez
 IMPLICIT NONE

 INTEGER, INTENT(IN   ) :: ibMax !Enver is not the right ibMax
 INTEGER, INTENT(IN   ) :: kMax
 CHARACTER(LEN=200 ), INTENT(IN   ) :: path 

 INTEGER(KIND=i8)                  :: nlev
 INTEGER(KIND=i8), PARAMETER       :: nvars=3
 REAL(KIND=r8)                     :: aux (1000,nvars)
 REAL(KIND=r8)                     :: aux2 (1000,nvars)
 REAL(KIND=r8)                     :: Ps
 REAL(KIND=r8)                     :: Pt
 INTEGER(KIND=i8)                  :: IL
 INTEGER(KIND=i8)                  :: k
 INTEGER(KIND=i8)                  :: i
 INTEGER(KIND=i8)                  :: j
 !-----------------------------------------------------------------
 REAL(KIND=r8),    INTENT(OUT  ) :: fgug    (ibMax,kMax)
 REAL(KIND=r8),    INTENT(OUT  ) :: fgvg    (ibMax,kMax)
 !-----------------------------------------------------------------
 REAL(KIND=r8)   , INTENT(IN   ) :: tod
 INTEGER(KIND=i8)                  :: ierr
 REAL(KIND=r8),SAVE       :: x1, x2
 REAL(KIND=r8)       :: tmp, xM, p1, p2
 REAL(KIND=r8),    INTENT(IN   ) :: colrad (ibMax)
 INTEGER(KIND=i8)    :: itmp
 INTEGER    :: jdt
! LOGICAL, SAVE    :: OpenFile =.TRUE.
 LOGICAL    :: DynTend  =.TRUE.  ! .FALSE. !Associado somente com a tendencia
                                 !         !devido ao momento zonal e meridional
                                 !         !Se (DynTend=.FALSE.) => auvel=avvel=0.0 

 aux=0.0; aux2=0.0

 IF( OpenFile_geo )THEN    !OpenFile

     OPEN(342,file=TRIM(path)//'/GEOSTROPHY_ASCII',form='formatted',&
         access='sequential',action='read',status='old' )

     READ(342,*) itmp, nlev, x1, Ps, tmp
     DO j=1,nlev
       READ(342,*)(aux (j,i),i=1,nvars)
     END DO


     ALLOCATE(sogeo(nlev,nvars)); ALLOCATE(sogeo2(nlev,nvars))
     ALLOCATE(sigma_geo(nlev)) 
     !ALLOCATE(gps  (nlev))

     ALLOCATE(sogeoa(kmax,nvars)); ALLOCATE(sogeoa2(kmax,nvars))
     sogeoa=0.0; sogeoa2=0.0

      sogeo(1:nlev,1)=aux(1:nlev,1) !p
      sogeo(1:nlev,2)=aux(1:nlev,2) !u
      sogeo(1:nlev,3)=aux(1:nlev,3) !v
   
      sigma_geo=1.0; sigma_geo(1) = 1.0

      !Enver now Ps is read from the header of the forcing file (start)
      !   Ps=aux(1,1) !Test sond(1,1)
      Pt=aux(nlev,1) !Test sond(nlev,1)
      sigma_geo(1)=1.0
      DO j=1,nlev-1
        sigma_geo(j+1)= (sogeo(j,1) - Pt)/(Ps - Pt) !Test (sond(j,1) - Pt)/(Ps - Pt) 
        !gps  (j)  = sigma_ndg(j)*(Ps - Pt) + Pt   !Test sigma(j)*(Ps - Pt) + Pt
      END DO   
!    print*, 'niveis sigma_ndg'
!    print*, sigma_ndg, si
!    print*, 'niveis sigma_ndg'
      DO IL=1,nvars
         CALL INTER2(sogeo(1:nlev,IL),nlev,si,sogeoa(1:kMax,IL),kMax,sigma_geo) 
      END DO
 
!      print*, 'sonda',sondga(1:KMax,5),sondga(1:KMax,6) 
     OpenFile_geo = .FALSE.
 ENDIF                 !OpenFile


! IF( MOD(tod,intfor) == 0.0 .OR. jdt == 1 )THEN !Read second profile
IF( jdt == 1 .OR. xM .ge. x2 )THEN !Read second profile !Enver 01Jul2013 to correct
                                                        !         a bug in time
                                                        !         interpolation.

  ! 28dec2014
  ! Here we have to add code to ensure that
  !        x2 > xM                 (x2 be larger than xM)
  ! such that reading or parse the next forcing profile
  ! result in
  !        x1 < xM < x2
  !
  !  28dec2014

   IF( jdt == 1 )THEN
     READ(342,*) itmp, nlev, x2, Ps, tmp
     DO j=1,nlev
        READ(342,*)(aux2 (j,i),i=1,nvars)
     END DO
   ELSE
     !aux = aux2
     !...
     sogeoa(:,2)    = sogeoa2(:,2)
     sogeoa(:,3)    = sogeoa2(:,3)
     x1            = x2  !Enver 28Jun2013
     READ(342,*) itmp, nlev, x2, Ps, tmp
     DO j=1,nlev
        READ(342,*)(aux2 (j,i),i=1,nvars)
     END DO
   ENDIF

      sogeo2(1:nlev,1)=aux2(1:nlev,1) !p
      sogeo2(1:nlev,2)=aux2(1:nlev,2) !u
      sogeo2(1:nlev,3)=aux2(1:nlev,3) !v
   
      sigma_geo=1.0; sigma_geo(1) = 1.0

      !Enver now Ps is read from the header of the forcing file (start)
      !   Ps=aux2(1,1) !Test sond(1,1)
      Pt=aux2(nlev,1) !Test sond(nlev,1)
      sigma_geo(1)=1.0
      DO j=1,nlev-1
        sigma_geo(j+1)= (sogeo2(j,1) - Pt)/(Ps - Pt) !Test (sond(j,1) - Pt)/(Ps - Pt) 
        !gps  (j)  = sigma_ndg(j)*(Ps - Pt) + Pt   !Test sigma(j)*(Ps - Pt) + Pt
      END DO   
      DO IL=1,nvars
         CALL INTER2(sogeo2(1:nlev,IL),nlev,si,sogeoa2(1:kMax,IL),kMax,sigma_geo) 
      END DO
 ENDIF !Read second profile


 !Computing weights by linear interpolation
  p1 = (x2 - xM)/(x2 - x1)
  p2 = 1.0 - p1

  fgug(ibmax,:)    =  ( p1*sogeoa(:,2) + p2*sogeoa2(:,2) ) * SIN(colrad (ibmax))
  fgvg(ibmax,:)    =  ( p1*sogeoa(:,3) + p2*sogeoa2(:,3) ) * SIN(colrad (ibmax))
  

 END SUBROUTINE Geostrophy_Ascii

 !CALL INTER(sond(:,IL),nlev,si,sonda(:,IL),kMax,sigma) 
 SUBROUTINE INTER2(Y,LEV,ALT,F,AN,XLEV)
  IMPLICIT NONE
  INTEGER(KIND=i8), INTENT(IN   ) :: LEV
  INTEGER, INTENT(IN   ) :: AN 
  REAL(KIND=r8)   , INTENT(IN   ) :: Y	 (LEV)
  REAL(KIND=r8)   , INTENT(OUT  ) :: F	 (AN )
  REAL(KIND=r8)   , INTENT(IN   ) :: ALT  (AN)
  REAL(KIND=r8)   , INTENT(IN   ) :: XLEV (LEV)

  INTEGER(KIND=i8)		 :: j
  INTEGER(KIND=i8)		 :: I
  REAL(KIND=r8)  		 :: undef=0.0 !undef=-9.99e33
  REAL(KIND=r8)  		 :: a(AN)
  REAL(KIND=r8)  		 :: b(AN)
  REAL(KIND=r8)  		 :: c(AN)

  a=0.0;b=0.0;c=0.0;F=0.0

  DO j=1,AN,1	       
    F(j)=undef     
  ENDDO
  DO j=1,AN,1
   DO 12 I=1,LEV-2,1			    
   IF(ALT(j).LE.XLEV(I).AND.ALT(j).GE.XLEV(I+2))THEN  
    IF((XLEV(I)-XLEV(I+1)).NE.0.0.AND.(XLEV(I)-XLEV(I+2)).NE. 0.0 )THEN	     

        a(j) = Y(I)/((XLEV(I)-XLEV(I+1))* &
  	       (XLEV(I)-XLEV(I+2)))
      
        b(j) = Y(I+1)/((XLEV (I+1)-XLEV(I))*&
  	       (XLEV(I+1)-XLEV(I+2)))
      
        c(j) = Y(I+2)/((XLEV(I+2)-XLEV(I))*&
  	       (XLEV(I+2)-XLEV(I+1)))
      
        F(j) =  a(j)*((alt(j)-XLEV(I+1))*(alt(j)-XLEV(I+2))) &
               +b(j)*((alt(j)-XLEV(I))*(alt(j)-XLEV(I+2)))	& 
               +c(j)*((alt(j)-XLEV(I))*(alt(j)-XLEV(I+1)))
  
       ENDIF		     
      ENDIF		    
12  ENDDO	
  ENDDO
      
 END SUBROUTINE INTER2

 SUBROUTINE LocalJulian(yi,mi,di,hi,tod,rxday)
    INTEGER, INTENT(IN   ) :: yi
    INTEGER, INTENT(IN   ) :: mi
    INTEGER, INTENT(IN   ) :: di
    INTEGER, INTENT(IN   ) :: hi
    REAL(KIND=r8)   , INTENT(OUT  ) :: rxday
    REAL(KIND=r8)   , INTENT(IN   ) :: tod
    REAL(KIND=r8)                   :: ss
    REAL(KIND=r8)                   :: yrl
    INTEGER(KIND=i8)                :: monl(12)
    INTEGER(KIND=i8)                :: monday(12)
    INTEGER(KIND=i8)                :: m
    REAL(KIND=r8)   , PARAMETER     :: f3600=3.6e3
   ! ss=0.0
    ss=tod
    yrl=365.25e0
    MONL = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    !
    !     id is now assumed to be the current date and hour
    !
    monday(1)=0
    DO m=2,12
       monday(m)=monday(m-1)+monl(m-1)
    END DO
    rxday=hi*f3600
    rxday=rxday+MOD(ss,f3600)
    rxday=monday(mi)+di+rxday/86400.0
    !rxday=rxday-MOD(yi+3,4)*0.25
    IF(MOD(yi,4).EQ.0.AND.mi.GT.2)rxday=rxday+1.0e0
    !rxday= MOD(rxday-1.0,yrl)
    !rxday= MOD(rxday,yrl)
 END SUBROUTINE LocalJulian

!LocalJulian_Year: Leads the model to deal with change of year during time integration
!28Dec2014

 SUBROUTINE LocalJulian_Year(yis,mis,dis,his,yi,mi,di,hi,tod,rxday)
    INTEGER, INTENT(IN   ) :: yis, yi
    INTEGER, INTENT(IN   ) :: mis, mi
    INTEGER, INTENT(IN   ) :: dis, di
    INTEGER, INTENT(IN   ) :: his, hi
    REAL(KIND=r8)   , INTENT(OUT  ) :: rxday
    REAL(KIND=r8)   , INTENT(IN   ) :: tod
    REAL(KIND=r8)                   :: ss
    REAL(KIND=r8)                   :: yrl
    INTEGER(KIND=i8)                :: monl(12)
    INTEGER(KIND=i8)                :: monday(12)
    INTEGER(KIND=i8)                :: m
    REAL(KIND=r8)   , PARAMETER     :: f3600=3.6e3
    INTEGER               :: yidx
    LOGICAL               :: isLeap

   ! ss=0.0
    ss=tod
    yrl=365.25e0
    MONL = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    !
    !     id is now assumed to be the current date and hour
    !
    monday(1)=0
    DO m=2,12
       monday(m)=monday(m-1)+monl(m-1)
    END DO
    rxday=hi*f3600
    rxday=rxday+MOD(ss,f3600)
    rxday=monday(mi)+di+rxday/86400.0
    !rxday=rxday-MOD(yi+3,4)*0.25
    !IF(MOD(yi,4).EQ.0.AND.mi.GT.2)rxday=rxday+1.0e0
       isLeap=(MOD(yi,4) .EQ. 0 .and. .not. MOD(yi,100) .EQ. 0).OR.(MOD(yi,4).EQ.400)
    IF(isLeap.AND.mi.GT.2)rxday=rxday+1.0e0
    !rxday= MOD(rxday-1.0,yrl)
    !rxday= MOD(rxday,yrl)

    !Adiciona el #ro de dias entre los anos yis e yi, tiene una verificacion
    !simple de los anos bisiestos

    yidx = yis
    do while ( yidx .lt. yi )
     rxday = rxday + yrl      !* ( yi - yidx )
      !Pope Gregory XIII criterium for leap year
      !https://alternativefocus.wordpress.com/2011/01/27/why-2100-will-not-be-a-leap-year/
       isLeap=(MOD(yidx,4) .EQ. 0 .and. .not. MOD(yidx,100) .EQ. 0).OR.(MOD(yidx,4).EQ.400)
      !IF( MOD(yidx,4) .EQ. 0 ) rxday = rxday + 1.0e0
       IF(isLeap) rxday = rxday + 1.0e0
     yidx = yidx + 1 
    end do    
 END SUBROUTINE LocalJulian_Year


END MODULE GridDynamics
