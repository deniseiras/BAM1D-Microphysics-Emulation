PROGRAM model
  
  USE ModRadiationDriver, ONLY:        &
       InitRadiationDriver         

  USE Surface, ONLY:          &
      InitSurface

  USE FieldsPhysics, ONLY:     &
       InitFieldsPhyscs,  &
       xland,capac0,w0 , prct, prcc,geshem ,gtsea,td0,zorl,&
       iMask,xice,lowlyr,ustar,z0,temp2m,umes2m,seamask,tracermix
             
  USE Convection, ONLY:       &
      InitConvection
              
  USE Sizes, ONLY:            &
      ibMax                 , &
      jbMax                 , &
      mnMax                 , &
      mMax                  , &
      nMax                  , &
      mnMap                 , &
      kMax                  , &
      imax                  , &
      jmax                  , &
      ijMax                 , &
      ijMaxGauQua           , &
      ibPerIJ               , &
      jbPerIJ               , &
      iPerIJB               , &
      jPerIJB               , &
      ibMaxPerJB,&
      si                    , & 
      sl                    , &
      cl                    , &
      ci                    , &
      del                   , &
      rpi  

  USE GridHistory, ONLY:      &
      InitGridHistory       , &
      WriteGridHistory      , &
      TurnOnGridHistory     , &
      IsGridHistoryOn

  USE InputOutput, ONLY:      &
      InitInputOutput       , &
      fsbc

             
  USE Diagnostics, ONLY:      &
       StartStorDiag, &
      InitDiagnostics       , &
      rsdiag                , &
      wridia                , &
      wrprog                , &
      reqdg                 , &
      combf                 , &
      dodia                 , &
      itavl                 , &
      iavrq                 , &
      nucf                  , &
      lvrq                  , &
      nurq                  , &
      lvcf                  , &
      itcf                  , &
      mxavl                 , &
      icf  

!  USE PlanBoundLayer, ONLY:   &
!      InitPlanBoundLayer
  
  USE Constants, ONLY:        &
      i8,r8,i4,r4,&
      ndavl,numx,ngrp,nunits,&
      ncf   ,&
      ncf2  ,&
      ngrmx ,&
      ndrq ,&
      ncdg ,&
      jxavl,&
      jxcdg,&
      tbar                  , &
      tov                   , &
      pai                   , &
      gasr                  , &
      cp

  
  USE FieldsDynamics, ONLY:        &
        fgtmpp,fgtmp,fgq,fgu,fgv,fglnps,fgps,fgrot,fgdiv,fgqp,fgpsp
  
  USE Init, ONLY :            &
      InitAll
  
  USE ModTimeStep, ONLY:      &
      TimeStep
      

  USE IOLowLevel, ONLY:    &
  InitReadWriteSpec         , &
  ReadProgHead              

  USE Utils, ONLY: &
       CreateGridValues, &
       InitTimeStamp, &
       TimeStamp, &
       tmstmp2
  
  USE  Options, ONLY:          &
       ReadNameList         , &
       SetOutPut,             &
       DumpOptions          , &
       SOND_IN              , &
       istrt                , &
       filta                , &
       filtb                , &
       delt                 , &
       nfin0                , &
       nfin1                , &
       initlz               , &
       nfcnv0               , &
       dogwd                , &
       isimp                , &
       grhflg               , &
       grhfl                , &
       igwd                 , &
       ifsst                , &
       sstlag               , &
       intsst               , &
       igrfu                , &
       iptu                 , &
       ighdr                , &
       monl                 , &
       ighou                , &
       enhdif               , &
       dt                   , &
       idate                , &
       idatec               , &
       idatef               , &
       kt                   , &
       ktm                  , &
       ktp                  , &
       jdt                  , &
       ddelt                , &
       rhdifd               , &
       rhdift               , &
       nfsibo               , &
       nfsibi               , &
       nf3d                 , &
       nstep                , &
       maxstp               , &
       isteps               , &
       nfcnv1               , &
       allghf               , &
       dodyn                , &
       maxtim               , &
       maxtfm               , &    
       nfvar                , &
       nferr                , &
       nfprt                , &
       ifprt                , &
       jdt2                 , &
       ifilt                , &
       trunc                , &
       vert                 , &
       nfdrct               , &
       nfdiag               , &
       WriteFaked           , &
       fNameInput0          , &
       fNameInput1          , &
       fNameSSTAOI          , &
       fNameSnow            , &
       fNameSoilms          , &
       fNameAlbedo          , &
       fNameCnftBl          , &
       fNameCnf2Tb          , &
       fNameLookTb          , &
       fNameSibVeg          , &
       fNameSibAlb          , &
       fNameDTable          , &
       fNameRDT             , &
       fNameOrgvar          , &
       fNameSibmsk          , &
       fNameTg3zrl          , &
       fNameIBISDeltaTemp   , &
       fNameIBISMask        , &
       fNameSandMask        , &
       fNameClayMask        , &
       fNameClimaTemp       , &
       fNameSlabOcen        , &
       fnamemicro           ,&
       start                , &
       deltOut              , &
       slagr                , &
       reducedGrid          , &
       linearGrid           , &
       nlnminit             , &
       diabatic             , &
       eigeninit            , &
       rsettov              , &
       ct_in                , &
       cq_in                , &
       del_in               , &
       nls                  , &
       nlcs                 , &
       iwrkdm               , &
       maxtid               , &
       cthl                 , &
       PREFY                , & 
       PREFX                , &  
       EXTF                 , & 
       EXDF                 , & 
       EXTH                 , & 
       EXDH                 , & 
       EXTW                 , & 
       EXDW                 , & 
       EXTS                 , &
       yrl                  , &
       nfsibd               , &
       nfsibt               , &
       sonda                , &
       path_in              , &
       dirfNameOutput
   USE PhysicalFunctions, ONLY: &
      InitPhysicalFunctions     

  USE PblDriver, ONLY:   &
       InitPBLDriver
  USE GwddDriver, ONLY:   &
       InitGWDDDriver
  USE SfcPBLDriver, ONLY: &
       InitSfcPBL_Driver
      
  IMPLICIT NONE
  INTEGER, PARAMETER   :: mgaus       = 100 
  INTEGER, PARAMETER   :: ngaus       = 64
  INTEGER, PARAMETER   :: mspec       = 4
  INTEGER, PARAMETER   :: nspec       = 0
  INTEGER, PARAMETER   :: nfsf        = 62
  INTEGER, PARAMETER   :: nfkm        = 17
  INTEGER, PARAMETER   :: ityp        = 13
  INTEGER, PARAMETER   :: imon        = 12
  INTEGER, PARAMETER   :: icg         = 2
  INTEGER, PARAMETER   :: iwv         = 3
  INTEGER, PARAMETER   :: ild         = 2
  INTEGER, PARAMETER   :: idp         = 3
  INTEGER, PARAMETER   :: ibd         = 2
  CHARACTER(LEN=200)   :: roperm
  CHARACTER(LEN=  7)   :: namef
  CHARACTER(LEN= 10)   :: labeli
  CHARACTER(LEN= 10)   :: labelc
  CHARACTER(LEN= 10)   :: labelf
  CHARACTER(LEN=  4)   :: TRC
  CHARACTER(LEN=  4)   :: LV
  LOGICAL              :: restart=.FALSE.
  LOGICAL(KIND=i8)     :: inirestart=.TRUE.
  REAL(KIND=r8)   , ALLOCATABLE :: lsmk(:)
  INTEGER              :: ifday
  INTEGER              :: it
  REAL(KIND=r8)                 :: tod
  REAL(KIND=r8)                 :: fa 
  REAL(KIND=r8)                 :: fb
  REAL(KIND=r8)                 :: fb1
  INTEGER              :: istart
  INTEGER              :: istep
  INTEGER              :: iend
  INTEGER              :: i
  INTEGER              :: ij
  INTEGER              :: j
  INTEGER              :: l
  INTEGER              :: ncount
  INTEGER              :: PlotStep 
  CHARACTER (LEN=10)   :: DateInit_s 
  CHARACTER (LEN=10)   :: DateNow_s
  CHARACTER (LEN=10)   :: cc1
  CHARACTER (LEN=10)   :: ModelName
  INTEGER              :: jhr
  INTEGER              :: jmon
  INTEGER              :: jday
  INTEGER              :: jyr
  REAL(KIND=r8)                 :: ahour
  INTEGER              :: iafday
  REAL(KIND=r8)                 :: itod
  INTEGER              :: jahr
  INTEGER              :: jaday
  INTEGER              :: jamon
  INTEGER              :: jayr
  INTEGER              :: kdt
  INTEGER              :: irstl
  LOGICAL(KIND=i8)              :: masbd
  INTEGER              :: mascd
  INTEGER              :: ifdy
  REAL(KIND=r8)                 :: todcld
  INTEGER              :: ids(4)
  INTEGER              :: idc(4)
  REAL(KIND=r8)                 :: totm 
  REAL(KIND=r8)                 :: todsib
  INTEGER              :: limlow
  INTEGER              :: ierr
  REAL(KIND=r8)                 :: fdh
  REAL(KIND=r8)                 :: dth
  REAL(KIND=r8)                 :: delth
  REAL(KIND=r8)                 :: fdayh
  REAL(KIND=r8)                 :: zero=0.0
  LOGICAL(KIND=i8)              :: idiaten
  LOGICAL(KIND=i8)              :: enhdifl
  LOGICAL(KIND=i8)              :: enhdifl2
  CHARACTER(LEN=4 )    :: nexp
  CHARACTER(LEN=40)    :: jttl
  CHARACTER(LEN=4)     :: grid='1D  '
  INTEGER              :: maxt0
  REAL(KIND=r8)   , ALLOCATABLE :: radwrk(:,:)  !(iwrkdm,4)
  REAL(KIND=r8)   , ALLOCATABLE :: rlsm(:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: GaussColat(:)   
  REAL(KIND=r8)   , ALLOCATABLE :: AuxGaussColat(:)   
  REAL(KIND=r8)   , ALLOCATABLE :: colradD(:)   
  REAL(KIND=r8)   , ALLOCATABLE :: rclD   (:)
  REAL(KIND=r8)   , ALLOCATABLE :: gzs_orig(:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: lsmk1    (:,:)
  REAL(KIND=r8)   , SAVE        :: long  
  REAL(KIND=r8)   , SAVE        :: longitude
  REAL(KIND=r8)   , SAVE        :: latitude 
  INTEGER              :: k 
  INTEGER              :: nffcst   
  NAMELIST /MODEL_1D/longitude,latitude 

  !*************************** Denis *****************************
  ! data for training
  ! Ferrier
  !open(unit=9876,file='/home/denis/_COM_BACKUP/NN_BAM1d/bam1d_data/ferrier_inputs.csv',status='unknown')
  !open(unit=8456,file='/home/denis/_COM_BACKUP/NN_BAM1d/bam1d_data/ferrier_outputs.csv',status='unknown')
  !write(9876, '(A)') 'k,rho,gt,qv,qc,qi,qr,f_ice_phy,f_rain_phy,f_rimef_phy,pgr,colrad'
  !write(8456, '(A)') 'k,rho,gt,qv,qc,qi,qr,f_ice_phy,f_rain_phy,f_rimef_phy,rainncv,snowncv'
  
  ! HUMO
  open(unit=9876,file='/home/denis/_COM_BACKUP/NN_BAM1d/bam1d_data/hug_morr_inputs.csv',status='unknown')
  open(unit=8456,file='/home/denis/_COM_BACKUP/NN_BAM1d/bam1d_data/hug_morr_outputs.csv',status='unknown')
  write(9876, '(A)') 'k,sl,Tc,qv,qc,qr,qi,qs,qg,ni,ns,nr,NG,NC,tke,omega'
  write(8456, '(A)') 'k,Tc,qv,qc,qr,qi,qs,qg,ni,ns,nr,NG,NC,LSRAIN,LSSNOW'
  
  ! open(unit=84562,file='/home/denis/_COM_BACKUP/NN_BAM1d/bam1d_data/hug_morr_rainsnow_puts.csv',status='unknown')
  ! write(84562, '(A)') 'LSRAIN,LSSNOW'

  ! Save output oringinal morrisson to compare with NN
  ! open(unit=111111,file='./hug_morr_outputs_ORIGINAL_PHYSICS.csv',status='unknown')
  ! write(111111, '(A)') '"k","Tc","qv","qc","qr","qi","qs","qg","ni","ns","nr","NG","NC","EFFCS","EFFIS","LSRAIN","LSSNOW"'

  ! Save output NN morrisson 
  open(unit=222222,file='./hug_morr_outputs_NN_PHYSICS.csv',status='unknown')
  write(222222, '(A)') '"k","Tc","qv","qc","qr","qi","qs","qg","ni","ns","nr","NG","NC","EFFCS","EFFIS","LSRAIN","LSSNOW"'

  NEXP    ='0003'
  grhfl   =.false.
  fsbc    =.true.
  ifday=0
  !
  ! input case
  !
  CALL ReadNameList() 
  READ (111,MODEL_1D)

  IF(trunc < 100)WRITE(TRC,'(a1,i3.3)')'T',trunc
  IF(trunc >= 100 .and. trunc < 1000)WRITE(TRC,'(a1,i3)')'T',trunc
  IF(trunc >= 1000)WRITE(TRC,'(a1,i4.4)')'T',trunc   
  IF(vert < 100)WRITE(LV,'(a1,i2.2)')'L',vert
  IF(vert >= 100 .and. vert < 1000)WRITE(LV,'(a1,i3.3)')'L',vert
     
  ALLOCATE(radwrk(iwrkdm,4))     

  PlotStep = NINT(deltOut/delt)
  
  !
  ! echo problem size
  !
  WRITE(*,"('** timestep output is ',i4)")PlotStep
  WRITE(ModelName,"('T',i3.3,'L',i2.2)") trunc, vert
  JTTL    =  'CPTEC AGCM R1.2 2001  '//TRIM(ModelName)//'  '//TRIM(start)
  WRITE(*,"('**')")
  WRITE(cc1,"(i10)") maxtim
  WRITE(*,"('** model ',a,' runs for ',a,' timesteps')") ModelName,TRIM(ADJUSTL(cc1))
  WRITE(cc1,"(i10)") INT(delt)
  WRITE(*,"('** timestep is ',a,' (s)')") TRIM(ADJUSTL(cc1))
  WRITE(cc1,"(i10)") PlotStep
  WRITE(*,"('** outputs are spaced by ',a,' timesteps')") TRIM(ADJUSTL(cc1))
  WRITE(*,"('** output file directory ',a)") TRIM(dirFNameOutput)
  WRITE(*,"('** model configuration: ')",ADVANCE='NO') 
  !
  ! initialize modules
  !
print*, 'before InitAll Enver'
  CALL InitAll(trunc, vert, reducedGrid, linearGrid, & 
               del_in, rhdifd, rhdift)
print*, 'after InitAll Enver'

  IF(grid == '2D  ') THEN
    !OPEN (197, FILE='/gfs/dk05/pkubota/PHYSCS-1.0.0/datain/GaussColat-1D.T062L28', &
    !        ACTION='read')

    ALLOCATE (GaussColat    (jMax/2))
    ALLOCATE (AuxGaussColat (jMax/2))
    ALLOCATE (colradD    (jMax)  )
    ALLOCATE (rclD       (jMax)  )
    ALLOCATE (gzs_orig      (iMax,jMax))
    ALLOCATE (lsmk1    (iMax,jMax))
    READ(197)GaussColat

    AuxGaussColat = 1.0_r8/(SIN(GaussColat)*SIN(GaussColat))

    DO j=1,jMax/2
      colradD(  j)     =   GaussColat(j)
      colradD(jMax+1-j)  =   pai-GaussColat(j)
      rclD   (j) =   AuxGaussColat(j)
      rclD   (jMax+1-j)  =   AuxGaussColat(j)
    END DO
  ELSE
    ALLOCATE (GaussColat    (1))
    ALLOCATE (AuxGaussColat (1))
    ALLOCATE (colradD    (1))
    ALLOCATE (rclD       (1))
    ALLOCATE (gzs_orig      (1,1))
    ALLOCATE (lsmk1    (1,1))
    long               = longitude
    GaussColat         = latitude
    GaussColat         = ((GaussColat + 90.0_r8)*pai)/180.0_r8
    AuxGaussColat      = 1.0_r8/(SIN(GaussColat)*SIN(GaussColat))
    colradD(1)         = GaussColat(1)
    rclD   (1)       = AuxGaussColat(1)
  END IF
  CALL CreateGridValues(colradD)

  CALL InitPhysicalFunctions(si,sl,kMax)

  CALL InitInputOutput(nfprt, nferr, ifprt, ngrmx,ncf   ,ncf2  , mMax  , &
                       nMax , mnMax, mnMap,kmax,path_in ,fNameSnow,fNameSSTAOI , &
                       fNameSSTAOI,fNameSoilms,fNameAlbedo,fNameCnftBl , &
                       fNameCnf2Tb,fNameLookTb)
StartStorDiag=.TRUE.
  IF(grid == '2D  ') THEN
      
  CALL InitDiagnostics(mgaus  ,ngaus  ,mspec  ,nspec  ,dodyn  , colradD , &
                       mMax   ,nMax   ,mnMax  ,mnMap  ,iMax   ,jMax    , &
                       kMax   ,ibMax  ,jbMax  ,ibMaxPerJB,grid,fNameDTable,&
                       fNameRDT)
  ELSE

  CALL InitDiagnostics(mgaus  ,ngaus  ,mspec  ,nspec  ,dodyn  , colradD , &
                       mMax   ,nMax   ,mnMax  ,mnMap  ,1   ,1    , &
                       kMax   ,1  ,1  ,ibMaxPerJB,grid,fNameDTable,&
                       fNameRDT)
  
  END IF

  CALL InitReadWriteSpec(&
             ndrq    ,ncdg    ,ndavl   ,mxavl   ,icf     ,mMax    , &
             mnMax   ,kMax    ,ijMax   ,iMax    ,jMax    ,ibMax   , &
             jbMax   ,nfprt   ,nferr   ,ifprt   , &
             reqdg   ,combf   ,dodia   ,itavl   ,iavrq   ,&
             nucf    ,lvrq    ,nurq    ,lvcf    ,itcf    )  

  IF(grid == '2D  ') THEN
    !CALL InitFieldsPhyscs(ibMax, kMax, jbMax,iMax,jMax)
    CALL InitFieldsPhyscs(ibMax, kMax, jbMax)
  ELSE

    !CALL InitFieldsPhyscs(1    , kMax,     1,   1,   1) 
    CALL InitFieldsPhyscs(1, kMax, 1)
  END IF
  DO jdt=0,maxtim
      !cehl(jdt)=.FALSE.
      !cdhl(jdt)=.FALSE.
     cthl(jdt)=.FALSE.
  END DO
  !     
  !     set isteps=1 due to common syntax
  !     
  dogwd=1
  !IF(grid == '2D  ') THEN
    !CALL InitVariancia(igwd  ,nfvar ,fNameOrgvar)
  !ELSE
    !CALL InitVariancia1D(igwd  ,nfvar ,fNameOrgvar)  
  !END IF    
  delt     = ddelt
  IF (nstep.eq.1) nstep=7

  IF(enhdif.eq.'YES ') THEN
      enhdifl = .TRUE.
  ELSE
      enhdifl = .FALSE.
  ENDIF
  !
  ! Initialize modules 
  !
  idatec(1)=idate(1)
  idatec(2)=idate(3)
  idatec(3)=idate(2)
  idatec(4)=idate(4)
  idate=idatec
  CALL InitTimeStamp (DateInit_s, idate)   
  labeli  = DateInit_s
  CALL TimeStamp     (DateNow_s, idatec, maxtim, delt)
  labelf  = DateNow_s
  labelc  = labelf 
  idatec  = idate
  REWIND (nfin1)
  !
  ! prepare output files
  !
  roperm  = dirfNameOutput
  namef   = "GFCT"//TRIM(PREFY)     
  !
  !     read cloud dataset - logic assumes that initialization not
  !     performed for warm start
  !
  jdt=0  
!  IF(grid == '2D  ') THEN
!    CALL InitBoundCond(ibMax,jbMax ,iwrkdm,ityp  ,ifdy  ,todcld,ids ,idc   ,ifday ,&
!                       tod  ,totm  ,todsib,radwrk,idate ,idatec,jdt ,si    ,sl    ,&
!                       fNameSibmsk,fNameTg3zrl ,ibMaxPerJB)
!  ELSE
print*, 'before SOND_IN Enver'
    CALL SOND_IN(path_in,gasr,cp,kmax,si,si,ci,del, sl ,cl ,rpi)
print*, 'after SOND_IN Enver'

!    CALL InitBoundCond1D(1,1 ,iwrkdm,ityp  ,ifdy  ,todcld,ids ,idc   ,ifday ,&
!                         tod  ,totm  ,todsib,radwrk,idate ,idatec,jdt ,si    ,sl    ,&
!                         fNameSibmsk,fNameTg3zrl ,ibMaxPerJB)  
!  END IF		       

  !
  !     write diabatic heating rate for nstep into file nfdbh
  !     
  tod=0.0_r8
  !
  !     start reading initial values (t=t   )
  !
  IF (nfin0 .eq. nfin1) THEN      
    !
    ! cold inicialization
    !  
    IF(grid == '2D  ') THEN
     ! READ(199)fgtmp  , fgq     ,fgu   , fgv    , fglnps
     !  fgps = EXP(fglnps)
    ELSE      
      fgtmpp(1,1:kmax,1) = sonda(1:kmax,2)-tov(1:kmax)
      fgtmp (1,1:kmax,1) = sonda(1:kmax,2)-tov(1:kmax)
      fgqp  (1,1:kmax,1) = sonda(1:kmax,3)  
      fgq   (1,1:kmax,1) = sonda(1:kmax,3)  
      fgu   (1,1:kmax,1) = sonda(1:kmax,4)
      fgv   (1,1:kmax,1) = sonda(1:kmax,5)
      fglnps(1,1)   = log(sonda(1,1)/10.0_r8)
      fgps  (1,1)   = EXP(fglnps(1,1))
      fgpsp (1,1)   = EXP(fglnps(1,1))
    END IF
  END IF
print*, 'in Main fglnps=', fglnps, fgps, fgpsp
  !
  !

  IF(ifsst.gt.3.and.sstlag.lt.zero)THEN
      WRITE(nfprt,336)ifsst,sstlag
      STOP 336
  END IF

  CALL InitPBLDriver(ibMax,jbMax,kmax, sl, del_in, si)

  CALL InitGWDDDriver(ibMax,jbMax,iMax,jMax,kmax,  si,sl,del,ibMaxPerJB,grid)
  CALL InitConvection(si,del_in        ,sl         ,cl         , &
                          kMax    ,iMax,jMax,ibMax,jbMax,&
                          trunc, ifdy       ,todcld     ,ids        , &
                          idc     ,ifday      ,tod           ,fNameMicro, idate          )
  CALL InitSurface(ibMax             ,jbMax             ,iMax              ,jMax          , &
                   kMax              ,del_in            ,path_in           ,fNameSibVeg   , &
                   fNameSibAlb       ,idate             ,idatec            ,dt            , &
                   nfsibd            ,nfprt             ,nfsibt            ,fNameSibmsk   , &
                   ifday             ,ibMaxPerJB        ,tod               ,ids           , &
                   idc               ,ifdy              ,todsib            ,fNameIBISMask , &
                   fNameIBISDeltaTemp,fNameSandMask     ,fNameClayMask     ,fNameClimaTemp, &
                   si                ,sl                ,RESTART           ,fgtmp         , &
                   fgq               ,fNameSlabOcen )


  CALL InitRadiationDriver(monl,yrl,kmax,sl,si,dt,nls)

  !
  !     write diagnostics/prognostics directory
  !
  CALL InitGridHistory (del_in ,nexp  ,jttl  ,idate ,allghf,grhflg,igrfu  , &
                        iptu   ,ighdr ,iMax  ,jMax  ,ibMax ,jbMax ,ibPerIJ, &
                        jbPerIJ,kMax)
  !
  CALL InitSfcPBL_Driver(kmax, sl, del_in, si,&
                        RESTART,ibMax  ,jbMax  ,USTAR  , &!(INOUT)
                        LOWLYR)
  !
  ! Write problem options to stdio
  !
  CALL DumpOptions()

  !
  fdh=24.0_r8
  dth=86400.0_r8/fdh
  delth=delt/dth
  fdayh=float(ifday)*fdh+tod/dth
  !
  maxt0=nint((float(ifday)*fdh+tod/dth)/delth) 
  CALL SetOutPut (tod,idatec)

  !     
  !     this is to remove accumulations from nlnmi
  !
  geshem=0.0_r8
  !
  !     clear all diagnostic accumulators
  !
  CALL rsdiag 
  !
  !     check files
  !     if nfin0=nfin1   then  cold start
  !
  limlow=1
  IF(nfin0.eq.nfin1)THEN
      !
      !     read cloud dataset for cold start
      !     
      ! CALL InitCheckfile(&
      !            ibMax  ,jbMax  ,ifdy  ,todcld,ids   ,idc   ,ifday , &
      !            tod   ,idate ,idatec,jdt   ,todsib  ,ibMaxPerJB)
      !
      !     cold start (at first delt/4 ,then delt/2 )
      !
      limlow =2
      dt= delt /4.0_r8
      !
      ! filter arguments for first time step
      !
      fa  = 0.0_r8
      fb  = 1.0_r8
      fb1 = 1.0_r8

  DO jdt=1,2
 
         IF (jdt.eq.2) THEN
            CALL rsdiag()
            CALL TurnOnGridHistory()
         END IF

         IF(jdt.eq.2.and.grhflg) grhfl=.true.

         istrt=jdt      
         kdt=jdt
         IF(ifprt(7).ge.1)WRITE(nfprt,104) kdt

         !enver   WRITE(*,"(' cold start step ',i2)") jdt
         !     
         !     calculate matrices for semi-implicit integration
         !     
         ! perform time step 
         !
         ifilt=0
         WRITE(*,*)ifday,idatec,idate
         idiaten=.FALSE.
         enhdifl2=.FALSE.

         CALL TimeStep(fb1,fa,fb,slagr,nlnminit,inirestart,idiaten,enhdifl2,dt,kt,ktm,ktp,jdt, & 
                       ifday,tod,idate, idatec,colradD,rclD,grid,long)

         ! prepare next time step, including filter arguments

      dt  =dt*2.0_r8
      jdt2=0
      ktm =kt    
      fsbc=.false.
      fb1 = 0.0_r8
  END DO
  tod=dt !PYK
  END IF
  !
  !     smooth start
  !
  IF(igwd.eq.'YES ')dogwd=0
  ahour=ifday*24.0_r8+tod/3600.0_r8
  istrt=0
  !
  ! filter arguments for all remaining time steps
  !
  fa     = filta
  fb     = filtb
  fb1    = 0.0_r8
  isteps = 1
  istep  = 1
  !
  ! time step loop
  !
  it=maxtfm   !Enver maxtim
  DO jdt=limlow,maxtim !Enver maxtim
      WRITE(*,"(' integration time step ',i5)") jdt
      !
      !     step loop starts
      !     
      IF(grhflg)grhfl=.true.
      
      IF(isimp.ne.'YES ') THEN  
        IF(grid == '2D  ') THEN
        !  CALL InitGetsbc  (ibMax,jbMax,ifday ,tod   ,idate ,idatec,jdt   )
        ELSE
	 ! CALL InitGetsbc1D(ibMax,jbMax,ifday ,tod   ,idate ,idatec,jdt   )
        END IF    
      END IF      
      !     
      !
      !        
      kdt=jdt      
      tod=tod+dt
      IF(abs( mod(tod+0.03125_r8,86400.0_r8)-0.03125_r8).lt.0.0625_r8)THEN
        tod=0.0_r8
        ifday=ifday+1
      END IF
      
      CALL tmstmp2(idate,ifday,tod,jhr,jday,jmon,jyr)
      idatec(1)=jhr
      idatec(2)=jmon
      idatec(3)=jday
      idatec(4)=jyr
      ahour=(ifday*24.0e0_r8)+(tod/3.6e3_r8)      
      kt   =int(ahour-(1.0e-2_r8))
      ktp  =int(ahour+(dt/3.6e3_r8)-(1.0e-2_r8))
      IF(jdt.eq.maxtim) THEN
         ktm=kt
      END IF
      !      
      ! perform time step 
      !
      ifilt=1
      jdt2=1
      idiaten=.FALSE.
      CALL TimeStep(fb1,fa,fb,slagr,nlnminit,inirestart,idiaten,enhdifl,dt,kt,ktm,ktp,jdt, &
                    ifday,tod,idate,idatec,colradD,rclD,grid,long)

      IF(IsGridHistoryOn())THEN
         CALL WriteGridHistory (ighou, ifday, tod, idate)
      END IF
      fsbc=.false.
      ktm=kt
      !
      ! output, if desired
      !	
      IF(MOD(jdt,plotstep)==0) THEN

        maxstp=nint((float(ifday)*fdh+tod/dth)/delth)-maxt0
       
        WRITE(*,"('Write file at timestep ',i5)") jdt
        
        IF (WriteFaked) THEN 
          !
          !     reset precip. every maxstp time steps
          !     
          IF(grid == '2D  ') THEN
             CALL wrprog (nfdrct ,nfdiag,ifday     ,tod   ,idate ,idatec,fgrot , &
                         fgdiv   ,fgq   ,fglnps    ,fgtmp ,zorl  ,gtsea ,td0   , &
                         capac0  ,w0    ,imask     ,nexp  ,jttl  ,nf3d  ,del_in, &
                         gzs_orig,lsmk1 ,iMax*jMax ,kmax  ,imax  ,jmax  , &
                         roperm  ,namef ,labeli    ,labelf,extw  ,exdw  ,TRIM(TRC), &
                         TRIM(LV),longitude,latitude,si,it,fgu   ,fgv   )  
             CALL wridia(nfdiag, maxstp)
          ELSE
             CALL wrprog (nfdrct ,nfdiag,ifday  ,tod   ,idate ,idatec,fgrot(1,:,1) , &
                         fgdiv(1,:,1)   ,fgq(1,:,1)   ,fglnps(1,1) ,fgtmp(1,:,1) ,zorl  ,gtsea ,td0   , &
                         capac0  ,w0    ,imask  ,nexp  ,jttl  ,nf3d  ,del_in, &
                         gzs_orig,lsmk1  ,1*1    ,kmax  ,1     ,1     , &
                         roperm  ,namef ,labeli ,labelf,extw  ,exdw  ,TRIM(TRC), &
                         TRIM(LV),longitude,latitude,si,it,fgu(1,:,1),fgv(1,:,1))
             CALL wridia(nfdiag, maxstp)
         END IF

          !
          !     zero reset diagnostic fields
          !     
          CALL rsdiag
           geshem=0.0_r8
           limlow=1     
        END IF

        IF(jdt.ne.maxtim)THEN
             maxt0=nint((float(ifday)*fdh+tod/dth)/delth)
        END IF
      END IF
      !
      fb1 = fb
  ENDDO

  ! **************** Denis
  close(9876)
  close(8456)
  ! close(84562)
  
  ! close(111111)
  close(222222)

  STOP
7 FORMAT( ' SIMPLIFIED PHYSICS OPTION IN EFFECT'/ &
       '  ALL OTHER OPTIONS OVERRIDDEN'/ &
       '  INITIAL CONDITIONS ISOTHERMAL RESTING ATMOSPHERE'/ &
       '  FLAT LOWER BOUNDARY'/ &
       ' PHYSICS:'/ &
       '  NEWTONIAN COOLING'/ &
       '  RALEIGH DAMPING'/ &
       '  NO WATER VAPOR OR CONDENSATION EFFECTS')
102 FORMAT(' DT=',G13.6,' IFDAY=',I10,' TOD=',G13.6, &
         ' SEC: WEPROG AT=',G13.6,'  STEP=',I10,' MAXSTP=',I5)
104 FORMAT(' ITERATION COUNT FOR THE COLD START=',I2)
336 FORMAT(' FOR IFSST=',I5,' SSTLAG MUST BE SET NONNEGATIVE.  NOT ',G12.5)         
END PROGRAM model
