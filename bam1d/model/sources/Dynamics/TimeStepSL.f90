!
!  $Author: pkubota $
!  $Date: 2010/03/22 14:30:49 $
!  $Revision: 1.16 $
!  Enver Ramirez 
!     modified to variate surface pressure along integration
!  Enver Ramirez
!     Jan 2020 modified to allow specified surface as in GASS/DCP protocol
!
MODULE ModTimeStep
  USE Constants, ONLY:     &
      r8,i8,&
      coriol             , &
      root2              , &
      tov          

  USE Sizes, ONLY:         &
      ibMax              , &
      jbMax              , &
      ibMaxPerJB         , &
      jPerIJB            , &
      iPerIJB            , & 
      mnMax              , &
      jMax               , &
      kMax               , &
      iMax               , &
      sl
    USE ModRadiationDriver,  ONLY:    &
         coszmed


  USE FieldsDynamics , ONLY :    &
               fguo, fgvo, fgto, fgqo,          & !variables for Nudging
               fgTsfc, fgH, fgLE, fgTau,        & !for specSfc
               fgug, fgvg,                      & !variables for Geostrophy
               fgyu   , fgyv        , fgtd    , &
               fgqd   , fgtmp        , fgq          , &
               fgum   , fgvm        , fgtmpm  , &
               fgqm   , omg        , fgps    , &
               fglnpm , fgzs    , fgqmm   , &
               fgumm  , fgvmm   , fgicem   ,fgtmpmm, &
               fgicet  , fgliqm   , fgliqt,fgtmpp , &
               fgqp   , &
               fgu    , &
               fgv    , &
               fgicep , &
               fgliqp ,fgpsp, &
         qice, fgice, fgicep, fgicem, fgicet, &
         qvar, fgvar, fgvarp, fgvarm, fgvart, &
         qliq, fgliq, fgliqp, fgliqm, fgliqt, &
         fgpass_scalars, adr_scalars
       

    USE FieldsPhysics, ONLY:  &
         ustr               , &
         vstr               , &
         imask              , &
         gtsea              , &
         PBL_CoefKm         , &
         PBL_CoefKh !hmjb            


  USE ModTimeFilter, ONLY: &
      TimeFilterStep1    , &
      TimeFilterStep2

  USE ModAddTend, ONLY:    &
      AddTend

  USE GridDynamics,  ONLY:   &
      GrpComp

  USE PhysicsDriver, ONLY :  &
      HumidPhysics 

  USE Options, ONLY :       &
      nClass             , &
      microphys          , &
      isimp              , &
      cdhl               , &
      first              , &
      nfdhn              , &
      vcrit              , &
      alpha              , &
      nfprt              , &
      nferr              , &
      yrl                , &
      monl               , &
      intcosz            , &
      initlz
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: TimeStep
CONTAINS

  !     timestep: performs a time-step of the model. Through the values
  !               of fb1, fa and fb we control if this is a initial
  !               time step,  part of a cold start or if it is a normal
  !               time step. 

  SUBROUTINE TimeStep (fb1, fa, fb, slagr, nlnminit, inirestart, idiaten, &
                       enhdif,dt,kt,ktm,ktp, jdt,ifday, tod, idate, idatec,&
                       colrad,rcs2 ,grid,longitude)
    REAL(KIND=r8)            , INTENT(IN)               :: fb1
    REAL(KIND=r8)            , INTENT(IN)               :: fa
    REAL(KIND=r8)            , INTENT(IN)               :: fb
    LOGICAL(KIND=i8)         , INTENT(IN)               :: slagr
    LOGICAL(KIND=i8)         , INTENT(IN)               :: nlnminit    
    LOGICAL(KIND=i8)         , INTENT(IN)               :: inirestart
    LOGICAL(KIND=i8)         , INTENT(IN)               :: idiaten   
    LOGICAL(KIND=i8)         , INTENT(IN)               :: enhdif
    REAL(KIND=r8)            , INTENT(IN)               :: dt
    INTEGER         , INTENT(IN)               :: jdt
    INTEGER         , INTENT(IN)               :: ifday
    REAL(KIND=r8)            , INTENT(IN)               :: tod
    INTEGER         , INTENT(IN)               :: idate(4)
    INTEGER         , INTENT(IN)               :: idatec (4)
    INTEGER         , INTENT(IN)               :: kt
    INTEGER         , INTENT(IN)               :: ktm
    INTEGER         , INTENT(IN)               :: ktp
    REAL(KIND=r8)            , INTENT(IN)               :: colrad(:)
    REAL(KIND=r8)            , INTENT(IN)               :: rcs2  (:)
    CHARACTER(len=*)         , INTENT(IN)               :: grid     
    REAL(KIND=r8)            , INTENT(IN)               :: longitude
    INTEGER         , SAVE                     :: ifp = 1    
    CHARACTER(len=4)         , SAVE                     :: siph = 'NO  '
    REAL(KIND=r8)            , ALLOCATABLE, SAVE        :: colrad2D(:,:) 
    REAL(KIND=r8)            , ALLOCATABLE, SAVE        :: rcl     (:,:) 
    REAL(KIND=r8)            , ALLOCATABLE, SAVE        :: lonrad  (:,:)
    REAL(KIND=r8)            , ALLOCATABLE, SAVE        :: lati    (:,:)       
    REAL(KIND=r8)            , ALLOCATABLE, SAVE        :: long    (:,:)   
    REAL(KIND=r8)                  :: cosz   
    REAL(KIND=r8)            , ALLOCATABLE, SAVE        :: cos2d   (:,:)   
    REAL(KIND=r8)            , ALLOCATABLE, SAVE        :: guqn(:,:)
    INTEGER                                    :: jb 
    INTEGER                                    :: ib 
    INTEGER                                    :: j 
    INTEGER                                    :: jhalf
    INTEGER                                    :: i,k
    INTEGER                                    :: ij
    REAL(KIND=r8)                                         :: vmax(kMax)
    REAL(KIND=r8)                                       :: g
    !
    !
    IF(ifp.EQ.1) THEN   
      IF(grid == '2D  ') THEN
       ALLOCATE (colrad2D(ibMax,jbMax))
       ALLOCATE (guqn    (ibMax,kMax))
       ALLOCATE (rcl     (ibMax,jbMax))
       ALLOCATE (lonrad  (ibMax,jbMax))
       ALLOCATE (long    (iMax,jMax))
       ALLOCATE (lati    (iMax,jMax))
       ALLOCATE (cos2d   (ibMax,jbMax))
       guqn=0.0_r8
       long  =0.0_r8
        DO j=1,jMax
         g=0.0_r8
          DO i=1,iMax
            long(i,j)=g
            g=g+(360.0_r8/REAL(iMax,kind=r8))
          END DO
        END DO
       
       colrad2D=0.0_r8
       lonrad=0.0_r8
       rcl=0.0_r8
       DO jb = 1, jbMax
          DO ib = 1, ibMaxPerJB(jb)
             j = jPerIJB(ib,jb)
             i = iPerIJB(ib,jb)
             jhalf = MIN(j, jMax-j+1)
             colrad2D(ib,jb)=colrad(j)
             rcl     (ib,jb)=rcs2(jhalf)
             lonrad  (ib,jb)=long(i,j)
             lati      (i,j)=colrad(j)
          END DO
       END DO 
      ELSE
       ALLOCATE (colrad2D(1,1))
       ALLOCATE (guqn    (1,kMax))
       ALLOCATE (rcl     (1,1))
       ALLOCATE (lonrad  (1,1))
       ALLOCATE (long    (1,1))
       ALLOCATE (lati    (1,1))
       ALLOCATE (cos2d   (1,1))
       guqn         =0.0_r8
       long         =longitude
       colrad2D(1,1)=colrad(1)
       rcl     (1,1)=rcs2  (1)
       lonrad  (1,1)=long  (1,1)
       lati    (1,1)=colrad(1)
      END IF
       ifp=0
    END IF

    IF(grid == '2D  ') THEN
      DO jb = 1, jmax
        CALL COSZMED(idatec,tod,yrl,colrad(jb),lonrad(:,jb),cosz,iMax)
      END DO
      DO jb = 1, jbMax
         DO ib = 1, ibMax
                cos2d  (ib,jb)=cosz
         END DO
      END DO
    ELSE
      CALL COSZMED(idatec,tod,yrl,colrad(1),lonrad(1,1),cosz,1)
      cos2d  (1,1)=cosz
    END IF
    
    IF (cdhl(jdt)) THEN
      ustr = 0.0_r8
      vstr = 0.0_r8
    END IF
    vmax = 0.0_r8
    !
    !  Complete filtering of previous time-step variables
    !  --------------------------------------------------
    !

    print*,'-1fgzs=', fgzs(1,1) 
    CALL TimeFilterStep2(fb1,grid)
    !IF (inirestart) RETURN    

    !
    ! Grid-point computations over latitudes
    ! --------------------------------------
    !

    print*,'fgzs=', fgzs(1,1)
    print*, 'fgu,v,t,q TimeStepSL', fgtmpm(1,1,1), fgtmpm(1,KMax-4,1), fglnpm, fgps, fgzs

    IF(grid == '2D  ') THEN
      DO jb = 1, jbMax
        DO ib = 1,ibMax
          CALL GrpComp ( &
               fgyu   (ib,:,jb), fgyv    (ib,:,jb), fgtd    (ib,:,jb), &
               fgqd   (ib,:,jb), fgtmp   (ib,:,jb), fgq     (ib,:,jb), &
               fgum   (ib,:,jb), fgvm    (ib,:,jb), fgtmpm  (ib,:,jb), &
               fgqm   (ib,:,jb), omg     (ib,:,jb), fgps    (ib  ,jb), &
               fglnpm (ib  ,jb), colrad2D(ib,  jb), rcl     (ib  ,jb), &
               fgzs   (ib  ,jb), & 
               dt              , ifday            , tod              , &
               idate           , idatec           , 1                , &
               kMax            , 1                , slagr            , &
               jdt             , kt               , ktm              , &
               ktp             , ib               , jb               , &
               lonrad(ib,jb)   , cos2d(ib,jb)     , intcosz          , &
               fguo(ib,:,jb), fgvo(ib,:,jb), fgto(ib,:,jb), fgqo(ib,:,jb)    , & !variables for Nudging
               fgTsfc, fgH, fgLE, fgTau, & !fgTsfc(ib,jb), fgH(ib,jb), fgLE(ib,jb), fgTau(ib,jb),                                       &  ! for specSfc
               fgug(ib,:,jb), fgvg(ib,:,jb)                                  , & !variables for Geostrophy
               fgicem(ib,:,jb) , fgicet(ib,:,jb)  , fgliqm (ib,:,jb)  , &
               fgliqt (ib,:,jb) )
        END DO
      END DO
    ELSE
       IF (.NOT. microphys) THEN
          CALL GrpComp ( &
               fgyu   (1,:,1), fgyv    (1,:,1), fgtd    (1,:,1), &
               fgqd   (1,:,1), fgtmp   (1,:,1), fgq     (1,:,1), &
               fgum   (1,:,1), fgvm    (1,:,1), fgtmpm  (1,:,1), &
               fgqm   (1,:,1), omg     (1,:,1), fgps    (1  ,1), &
               fglnpm (1  ,1), colrad2D(1,  1), rcl     (1  ,1), &
               fgzs   (1  ,1), & 
               dt            , ifday          , tod            , &
               idate         , idatec         , 1           , &
               kMax          , 1              , slagr          , &
               jdt           , kt             , ktm            , &
               ktp           , 1              , 1           , &
               lonrad(1,1)   , cos2d(1,1)     , intcosz        , &  
               fguo(1,:,1), fgvo(1,:,1), fgto(1,:,1), fgqo(1,:,1)  ,& !variables for Nudging
               fgTsfc, fgH, fgLE, fgTau, & !fgTsfc(1,1), fgH(1,1), fgLE(1,1), fgTau(1,1),                                       &  ! for specSfc
               fgug(1,:,1), fgvg(1,:,1)                         ) !variables for Geostrophy
       ELSE
          IF(nClass>0 )THEN
          CALL GrpComp ( &
               fgyu   (1,:,1), fgyv    (1,:,1), fgtd    (1,:,1), &
               fgqd   (1,:,1), fgtmp   (1,:,1), fgq     (1,:,1), &
               fgum   (1,:,1), fgvm    (1,:,1), fgtmpm  (1,:,1), &
               fgqm   (1,:,1), omg     (1,:,1), fgps    (1  ,1), &
               fglnpm (1  ,1), colrad2D(1,  1), rcl     (1  ,1), &
               fgzs   (1  ,1), & 
               dt            , ifday          , tod            , &
               idate         , idatec         , 1           , &
               kMax          , 1              , slagr          , &
               jdt           , kt             , ktm            , &
               ktp           , 1              , 1           , &
               lonrad(1,1)   , cos2d(1,1)     , intcosz        , &
               fguo(1,:,1), fgvo(1,:,1), fgto(1,:,1), fgqo(1,:,1)    , & !variables for Nudging
               fgTsfc, fgH, fgLE, fgTau, & !fgTsfc(1,1), fgH(1,1), fgLE(1,1), fgTau(1,1),                                       &  ! for specSfc
               fgug(1,:,1), fgvg(1,:,1)                          , & !variables for Geostrophy
               fgicem(1,:,1)  , fgicet(1,:,1)   , fgliqm(1,:,1)    ,&
               fgliqt(1,:,1)  ,&
               fgvarm(1,:,1,:)  , fgvart(1,:,1,:)  )
           ELSE
          CALL GrpComp ( &
               fgyu   (1,:,1), fgyv    (1,:,1), fgtd    (1,:,1), &
               fgqd   (1,:,1), fgtmp   (1,:,1), fgq     (1,:,1), &
               fgum   (1,:,1), fgvm    (1,:,1), fgtmpm  (1,:,1), &
               fgqm   (1,:,1), omg     (1,:,1), fgps    (1  ,1), &
               fglnpm (1  ,1), colrad2D(1,  1), rcl     (1  ,1), &
               fgzs   (1  ,1), & 
               dt            , ifday          , tod            , &
               idate         , idatec         , 1           , &
               kMax          , 1              , slagr          , &
               jdt           , kt             , ktm            , &
               ktp           , 1              , 1           , &
               lonrad(1,1)   , cos2d(1,1)     , intcosz        , &
               fguo(1,:,1), fgvo(1,:,1), fgto(1,:,1), fgqo(1,:,1)    , & !variables for Nudging
               fgTsfc, fgH, fgLE, fgTau, &!fgTsfc(1,1), fgH(1,1), fgLE(1,1), fgTau(1,1),                                       &  ! for specSfc
               fgug(1,:,1), fgvg(1,:,1)                          , & !variables for Geostrophy
               fgicem(1,:,1)  , fgicet(1,:,1)   , fgliqm(1,:,1)    ,&
               fgliqt(1,:,1)    )

           END IF
    END IF

    END IF
    !
    
    first = .FALSE.
    !
    !  Perform semi-Lagrangian computations
    !  ------------------------------------
    ! IF (slagr)  CALL semilagr ...
    !
    !   Em semilagr devem ser calculados os pontos de partida e pontos medios
    !   da trajetoria ( guv e gum contem os valores do vento).
    !   Os valores das tendencias (gyu, gyv, gtd) devem ser interpolados
    !   para as posicoes dos pontos medios da trajetoria. As tendencias 
    !   (gyum, gyvm, gtdm, gvdlnpm)  do tempo anterior devem ser
    !     interpoladas para os pontos de partida da trajetoria.
    !   Os respectivos valores interpolados devem ser armazenados nas mesmas variaveis
    !   (na rotina seguinte addtend,  as tendencias serao finalizadas.
    !   (Uma formulacao alternativa e calcular a tendencia no tempo intermediario como
    !   a media entre seus valores nos pontos de grade (ja armazenados em gyu, gyv e gtd)
    !   e estes mesmos valores interpolados para os pontos de partida 
    !
    !  Finish tendencies
    !  -----------------
    ! 
    IF (slagr)  THEN 
!-1D       CALL SemiLagr (dt)  ! computes departure points and add tendencies.
    ELSE
       CALL AddTend  (dt, nlnminit,grid)
    END IF
    !
    !  Begin filtering of previous time-step variables
    !  -----------------------------------------------
    !
    CALL TimeFilterStep1(fa, fb,grid)
    !
    ! Grid-point computations for water
    ! ---------------------------------
    !
    !     
    !     perform moist ,large scale & dry convection
    !     
    IF(TRIM(isimp) .ne.'YES ') THEN
      IF(grid == '2D  ') THEN

       DO jb = 1, jbMax   
        DO  k = 1, kMax 
         DO ib = 1, ibMax
          fgtmp(ib,k,jb)=  fgtmp(ib,k,jb)+tov(k)
         END DO 
        END DO
       END DO
       DO jb = 1, jbMax 
        DO ib = 1, ibMax
         !CALL HumidPhysics(jb, ib,guqn(ib,:),fgqmm(ib,:,jb), fgtmp(ib,:,jb),&
         !                  fgq(ib,:,jb), fgps(ib,  jb),fgumm(ib,:,jb)  , &
         !                  fgvmm(ib,:,jb),omg(ib,:,jb),1,kMax)
        END DO
       END DO
       DO jb = 1, jbMax   
        DO  k = 1, kMax 
         DO ib = 1, ibMax
          fgtmp(ib,k,jb)=  fgtmp(ib,k,jb)-tov(k)
         END DO 
        END DO
       END DO
      ELSE
         fgtmp(1,1:kMax,1)=  fgtmp(1,1:kMax,1)+tov(1:kMax)
         jb=1
         ib=1
           IF (.NOT. microphys) THEN

         CALL HumidPhysics(1 , 1  ,kMax  , fgqmm(ib,1:kMax,jb), fgtmp(ib,1:kMax,jb) , &
                           fgq    (ib,1:kMax,jb), fgpsp   (ib,       jb), fgum   (ib,1:kMax,jb) , &
                           fgvm   (ib,1:kMax,jb), omg     (ib,1:kMax,jb), fgtmpmm(ib,1:kMax,jb) , &
                           fgtmpm (ib,1:kMax,jb), fgqm     (ib,1:kMax,jb), fgps  (ib       ,jb) , &
                           fgzs  (ib       ,jb), colrad2D(ib,       jb),lonrad(1,jb) )

        else
         IF(nClass>0 )THEN

         CALL HumidPhysics(1 , 1  ,kMax  , fgqmm(ib,1:kMax,jb), fgtmp(ib,1:kMax,jb) , &
                           fgq    (ib,1:kMax,jb), fgpsp   (ib,       jb), fgum   (ib,1:kMax,jb) , &
                           fgvm   (ib,1:kMax,jb), omg     (ib,1:kMax,jb), fgtmpmm(ib,1:kMax,jb) , &
                           fgtmpm (ib,1:kMax,jb), fgqm     (ib,1:kMax,jb), fgps  (ib       ,jb) , &
                           fgzs  (ib       ,jb), colrad2D(ib,       jb) ,lonrad(1,jb), fgicem(ib,1:kMax,jb) , &
                           fgicep(ib,1:kMax,jb), fgliqm  (ib,1:kMax,jb), fgliqp(ib,1:kMax,jb)  ,&
                           fgvarm(1:ibMax,1:kMax,jb,1:nClass),fgvarp(1:ibMax,1:kMax,jb,1:nClass)   )
         ELSE
         CALL HumidPhysics(1 , 1  ,kMax  , fgqmm(ib,1:kMax,jb), fgtmp(ib,1:kMax,jb) , &
                           fgq    (ib,1:kMax,jb), fgpsp   (ib,       jb), fgum   (ib,1:kMax,jb) , &
                           fgvm   (ib,1:kMax,jb), omg     (ib,1:kMax,jb), fgtmpmm(ib,1:kMax,jb) , &
                           fgtmpm (ib,1:kMax,jb), fgqm     (ib,1:kMax,jb), fgps  (ib       ,jb) , &
                           fgzs  (ib       ,jb), colrad2D(ib,       jb) ,lonrad(1,jb), fgicem(ib,1:kMax,jb) , &
                           fgicep(ib,1:kMax,jb), fgliqm  (ib,1:kMax,jb), fgliqp(ib,1:kMax,jb)  )

         END IF
         endif
        fgtmp(1,1:kMax,1)=  fgtmp(1,1:kMax,1)-tov(1:kMax)  
      END IF
    END IF
  END SUBROUTINE TimeStep
END MODULE ModTimeStep
