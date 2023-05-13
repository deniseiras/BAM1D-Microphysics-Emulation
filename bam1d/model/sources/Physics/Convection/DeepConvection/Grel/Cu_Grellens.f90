!
!  $Author: pkubota $
!  $Date: 2011/05/27 23:15:12 $
!  $Revision: 1.6 $
!
MODULE Cu_Grellens
  !
  ! grellens--|---grellens2-|
  !                         |cup_env
  !                         |  
  !                         |cup_env
  !                         |  
  !                         |cup_env_clev
  !                         |  
  !                         |cup_env_clev
  !                         |  
  !                         |cup_maximi
  !                         |  
  !                         |cup_kbcon
  !                         |  
  !                         |cup_minimi
  !                         |  
  !                         |cup_up_he
  !                         |  
  !                         |cup_up_he
  !                         |  
  !                         |cup_ktop
  !                         |  
  !                         |cup_minimi
  !                         |  
  !                         |cup_up_nms
  !                         |  
  !                         |cup_up_nms
  !                         |  
  !                         |cup_dd_nms
  !                         |  
  !                         |cup_dd_nms
  !                         |  
  !                         |cup_dd_he
  !                         |  
  !                         |cup_dd_he
  !                         |  
  !                         |cup_dd_moisture
  !                         |  
  !                         |cup_dd_moisture
  !                         |  
  !                         |cup_up_moisture
  !                         |  
  !                         |cup_up_moisture
  !                         |  
  !                         |cup_up_aa0
  !                         |  
  !                         |cup_up_aa0
  !                         |  
  !                         |cup_dd_edt
  !                         |  
  !                         |cup_dellabot
  !                         |  
  !                         |cup_dellabot
  !                         |  
  !                         |cup_dellas
  !                         |  
  !                         |cup_dellas
  !                         |  
  !                         |cup_env
  !                         |  
  !                         |cup_env_clev
  !                         |  
  !                         |cup_up_he
  !                         |  
  !                         |cup_up_nms
  !                         |  
  !                         |cup_dd_nms
  !                         |  
  !                         |cup_dd_he
  !                         |  
  !                         |cup_dd_moisture
  !                         |  
  !                         |cup_up_moisture
  !                         |  
  !                         |cup_up_aa0
  !                         |  
  !                         |cup_maximi
  !                         |  
  !                         |cup_kbcon
  !                         |  
  !                         |cup_forcing_ens_16
  !                         |  
  !                         |cup_output_ens
  ! 
  USE Constants, ONLY :  &
       grav          , &
       undef           , &
       r8  
  USE Options, ONLY :       &
       grepar1           , &! integer: 0 ensemble 1 GRE   4 OMG   7 KUO  10 Chappel 13 ARA   24 ensemble2     
       grepar2           , &
       grepar3           , &
       grepar4 
 USE Parallelism, ONLY: &
      myId,&
      MsgOne, &
      FatalError

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: grellens,InitGrellens
  REAL(KIND=r8)   , PARAMETER :: cp   =1004.0_r8
  REAL(KIND=r8)   , PARAMETER :: xl   =2.5e06_r8
  REAL(KIND=r8)   , PARAMETER :: rv   =461.9_r8
  !----------------------------------------------------------------------------
  REAL(KIND=r8)   , PARAMETER :: detra = 10000.0_r8! detrainment value, original 12000. (10000-12000)
  REAL(KIND=r8)   , PARAMETER :: xfmax = 0.10_r8! xff max value ocean. 0.0---0.2  (best 0.05-0.10) 
  !----------------------------------------------------------------------------
  INTEGER, PARAMETER :: iens     =  1 
  INTEGER, PARAMETER :: iens_tmp =  1
  INTEGER, PARAMETER :: mjx      =  1 
  INTEGER, PARAMETER :: maxens   =  3 ! ensemble one on mbdt from PARAME(mm5)
  INTEGER, PARAMETER :: maxens2  =  3 ! ensemble two on precip efficiency
  INTEGER, PARAMETER :: maxens3  = 20 ! ensemble three done in cup_forcing  !!maxens3  = 16
  INTEGER, PARAMETER :: ensdim   = 1*maxens*maxens2*maxens3 !144
  REAL(KIND=r8)   , PARAMETER :: tcrit    =   273.15_r8
  !
  ! workfunctions for downdraft
  !
  REAL(KIND=r8)   , PARAMETER :: beta1=0.0_r8
  REAL(KIND=r8)   , PARAMETER :: beta2=0.25_r8
  REAL(KIND=r8)   , PARAMETER :: beta3=0.50_r8
  REAL(KIND=r8)   , PARAMETER :: beta( 1:maxens3)=(/&
       1.00_r8,1.00_r8,1.00_r8,1.00_r8,&
       1.00_r8,1.00_r8,0.00_r8,0.25_r8,&
       0.50_r8,1.00_r8,1.00_r8,1.00_r8,&
       1.00_r8,1.00_r8,1.00_r8,1.00_r8,&
       1.00_r8,1.00_r8,1.00_r8,1.00_r8 /)
  INTEGER, PARAMETER :: mkxcrt=15
  REAL(KIND=r8)   , PARAMETER :: pcrit( 1:mkxcrt)=(/&
       850.0_r8,800.0_r8,750.0_r8,700.0_r8,650.0_r8,600.0_r8,550.0_r8,500.0_r8,450.0_r8,400.0_r8, &
        35.0_r8,300.0_r8,250.0_r8,200.0_r8,150.0_r8/)

  REAL(KIND=r8)   , PARAMETER :: acrit( 1:mkxcrt)=(/&
       0.0633_r8,0.0445_r8,0.0553_r8,0.0664_r8,0.0750_r8,0.1082_r8,0.1521_r8,0.2216_r8, &
       0.3151_r8,0.3677_r8,0.4100_r8,0.5255_r8,0.7663_r8,1.1686_r8,1.6851_r8/)

  ! GDAS derived acrit

  REAL(KIND=r8)   , PARAMETER :: acritt (1:mkxcrt)=(/&
       0.203_r8,0.515_r8,0.521_r8,0.566_r8,0.625_r8,0.665_r8,0.659_r8,0.688_r8, &
       0.743_r8,0.813_r8,0.886_r8,0.947_r8,1.138_r8,1.377_r8,1.896_r8/)

  REAL(KIND=r8)   , PARAMETER :: dec_fudge = -0.20_r8
  REAL(KIND=r8)            :: ae  (2)
  REAL(KIND=r8)            :: be  (2)
  REAL(KIND=r8)            :: ht  (2)

  REAL(KIND=r8)               :: radius
  REAL(KIND=r8)               :: entr_rate
  REAL(KIND=r8)               :: mentrd_rate
  REAL(KIND=r8)               :: mentr_rate
  REAL(KIND=r8)               :: edtmin
  REAL(KIND=r8)               :: edtmax
  REAL(KIND=r8)               :: edtmax1
  REAL(KIND=r8)               :: effmax
  REAL(KIND=r8)               :: depth_min
  REAL(KIND=r8)               :: cap_maxs 
  REAL(KIND=r8)               :: cap_maxs_land 
  REAL(KIND=r8)               :: cap_max_increment
  REAL(KIND=r8)               :: zkbmax
  REAL(KIND=r8)               :: zcutdown
  REAL(KIND=r8)               :: z_detr


CONTAINS
  SUBROUTINE InitGrellens()
    !
    ! specify entrainmentrate and detrainmentrate
    ! Larger radius will give less mass fluix and make cloud grow taller
    ! and shift heating. Recomend 10 km.
    ! snf        radius=12000.
    !
    radius=detra
    !
    !  gross entrainment rate
    !
    entr_rate=0.2_r8/radius
    !
    ! entrainment of mass
    !
    mentrd_rate=entr_rate
    mentr_rate =entr_rate
    !
    ! initial detrainmentrates
    !
    !
    ! max/min allowed value for epsilon (ratio downdraft base 
    ! mass flux/updraft base mass flux
    !
    edtmin=0.2_r8
    !
    ! snf   edtmax=0.60_r8
    ! snf   edtmax=0.75_r8
    !
    edtmax=0.99_r8    ! ok over land snf may 2004 
    edtmax1=0.99_r8   ! Ok testado no oceano... com 0.8_r8 fica mal na NINA8889 
    !---------------
    effmax=0.99_r8
    edtmax=effmax
    edtmax1=effmax 
    !
    !  minimum depth (m), clouds must have
    !
    depth_min=500.0_r8
    !
    ! maximum depth (mb) of capping
    ! larger cap = no convection    !!!!!!!not true
    !     if(iens == 3)cap_max2=50.0_r8
    !     if(iens == 2)cap_max2=75.0_r8 
    !     if(iens == 1)cap_max2=100.0_r8
    ! original    cap_maxs=125.0_r8 cap_max_increment=50.0_r8      !new
    !
    cap_maxs=grepar3
    cap_maxs_land=grepar3
    !--------------------------------
    !!    cap_max_increment=50.0_r8
    !!snf modified   
    cap_max_increment=grepar4

    !
    ! max height(m) above ground where updraft air can originate
    !
    zkbmax=4000.0_r8
    !
    ! height(m) above which no downdrafts are allowed to originate
    !
    zcutdown=3000.0_r8
    !
    ! depth(m) over which downdraft detrains all its mass
    !
    z_detr=1250.0_r8
    !
    ht(1)=xl/cp
    
    ht(2)=2.834e6_r8/cp
    
    be(1)=0.622_r8*ht(1)/0.286_r8
    
    ae(1)=be(1)/273.0_r8+LOG(610.71_r8)
    
    be(2)=0.622_r8*ht(2)/0.286_r8
    
    ae(2)=be(2)/273.0_r8+LOG(610.71_r8)

  END SUBROUTINE InitGrellens

  SUBROUTINE grellens(&
       ps     ,sl     ,delsig ,ua     ,va     ,omg    ,t2     ,tn1    ,q2     ,&
       qn1    ,zz     ,xland  ,dtime  ,RAINCV ,kuo    ,ktop   ,kbot   ,plcl   ,&
       qliq   ,nCols  ,kMax   , tke   ,latco  ,lon    ,lat    ,ddql   ,ddmu   ,&
       cape   ,cine_LS,sens   ,cape_old)

    !--------------------------------------------------------------------
    !  This convection  f90 program subroutine has been modified by figueroa 
    !  in dec2003-jan2003 from mm5-Chemestry model and it has been adapted
    !  to CPTEC-GCM. 
    !- Some modification are:
    !- adapting the original parameters for T062L28 and T126L28 GCM.
    !- modified Betas in Kuo-clousure after testing 
    !- addition  some new changes from RAMS version (after testing)
    !- In April 2004 was added caps ensemble from Grell new code
    !- Omega closure was removed from ensemble closure. However it can be 
    !  use as single closure.
    !- This convection scheme was tested for different closures and 
    !  different
    !  parameters.  After several experiments and long integrations we 
    !  have suggested some
    !  critical values for CPTEC GCM. You can find them in the namelist.
    !
    !  NOTE: It is not official version yet!! It can be used only for tests
    !  more informations: nilo@cptec.inpe.br
    !-------------------------------------------------------------------
    !--------------------------------------------------------------------
    ! nilo.figueroa 2012 sep6-2012
    !-------------------------------------------------------------------

    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: nCols
    INTEGER, INTENT(IN   ) :: kMax
    REAL(KIND=r8)   , INTENT(IN   ) :: dtime
    REAL(KIND=r8)   , INTENT(IN   ) :: ua    (1:nCols,1:kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: va    (1:nCols,1:kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: omg   (1:nCols,1:kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: t2    (1:nCols,1:kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: q2    (1:nCols,1:kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: tn1   (1:nCols,1:kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: qn1   (1:nCols,1:kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: ps    (1:nCols       )
    REAL(KIND=r8)   , INTENT(IN   ) :: zz    (1:nCols       )
    INTEGER         , INTENT(IN   ) :: xland (1:nCols       )
    REAL(KIND=r8)   , INTENT(IN   ) :: sl    (1:kMax        )
    REAL(KIND=r8)   , INTENT(IN   ) :: delsig(1:kMax        )
    INTEGER         , INTENT(INOUT) :: kuo   (1:nCols       )
    INTEGER         , INTENT(INOUT) :: ktop  (1:nCols       )
    INTEGER         , INTENT(INOUT) :: kbot  (1:nCols       )
    REAL(KIND=r8)   , INTENT(OUT  ) :: qliq (1:nCols,1:kMax)
    !
    ! output variables after cumulus parameterization
    !
    REAL(KIND=r8)   , INTENT(INOUT) :: RAINCV(1:nCols       )
    REAL(KIND=r8)   , INTENT(INOUT) :: plcl  (1:nCols       ) 
!---------nilo
    INTEGER         , INTENT(IN   ) :: latco
    REAL(KIND=r8)   , INTENT(IN   ) :: lon(nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: lat(nCols)
    REAL(KIND=r8)                   :: lat2(nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: tke   (1:nCols,1:kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: ddql (1:nCols,1:kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: ddmu (1:nCols,1:kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: cape_old(1:nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: cape(1:nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: cine_LS(1:nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: sens(1:nCols)
!------------------
    !
    ! local variables
    !
    REAL(KIND=r8)                   :: t        (nCols,kMax)
    REAL(KIND=r8)                   :: q        (nCols,kMax)
    REAL(KIND=r8)                   :: tn       (nCols,kMax)
    REAL(KIND=r8)                   :: qo       (nCols,kMax)
    REAL(KIND=r8)                   :: p        (nCols,kMax)
    REAL(KIND=r8)                   :: po       (nCols,kMax)
    REAL(KIND=r8)                   :: us       (nCols,kMax)
    REAL(KIND=r8)                   :: vs       (nCols,kMax)
    REAL(KIND=r8)                   :: omeg     (nCols,kMax)
!    REAL(KIND=r8)                   :: ter11    (nCols     )
    REAL(KIND=r8)                   :: psur     (nCols     )
    REAL(KIND=r8)                   :: outt     (nCols,kMax) 
    REAL(KIND=r8)                   :: outq     (nCols,kMax) 
    REAL(KIND=r8)                   :: outqc    (nCols,kMax)
    REAL(KIND=r8)                   :: pre1     (nCols     )
!    REAL(KIND=r8)                   :: prec     (nCols     )
!    REAL(KIND=r8)                   :: dt
    INTEGER                         :: i
    INTEGER                         :: k
    INTEGER                         :: kk
    !
    REAL(KIND=r8)                   :: cupclw   (nCols,kMax)
    REAL(KIND=r8)                   :: bncy     (nCols,kMax)
    REAL(KIND=r8)                   :: massfln  (nCols,ensdim)
!    REAL(KIND=r8)                   :: down_massflx(nCols)
    INTEGER                         :: ierr     (nCols)

!new 2011-dec
    REAL(KIND=r8)                   :: xmb      (nCols     )
    REAL(KIND=r8)                   :: out1_mass_flux_up(nCols,kMax)
    REAL(KIND=r8)                   :: out2_mass_flux_d(nCols,kMax)
    REAL(KIND=r8)                   :: out3_detr(nCols,kMax)
    REAL(KIND=r8)                   :: out4_liq_water(nCols,kMax)
    REAL(KIND=r8)                   :: out5_qc(nCols,kMax)
    REAL(KIND=r8)                   :: out6_wf(nCols)
    REAL(KIND=r8)                   :: out7_entr(nCols,kMax)
    REAL(KIND=r8)                   :: qlmedia(nCols)
    REAL(KIND=r8)                   :: ddmd(nCols,kMax)  
    REAL(KIND=r8)                   :: dden(nCols,kMax)
    REAL(KIND=r8)                   :: ddde(nCols,kMax)
    REAL(KIND=r8)                   :: ddca(nCols)
    REAL(KIND=r8)                   :: dp1, dp2,aa,bb,cc,dd,ee
    REAL(KIND=r8)                   :: dcape(nCols)
    REAL(KIND=r8)                   :: cine (nCols)
    
!new 2012
    REAL(KIND=r8)                   :: pw(nCols)


    !
    !    prepare input, erase output
    !
    DO i=1,nCols
       kuo (i)=0
       kbot(i)=1
       ktop(i)=1
       ierr(i)=0
       plcl(i) = 1.0e0_r8
       lat2(i)=90.0_r8-lat(i)*180.0_r8/3.1415926_r8
       cine(i)=0.0_r8
    END DO
    !
    !move  variables from GCM to local variables
    !
    !    dt=dtime
    !    ter11(1:nCols) = zz(1:nCols)
    psur (1:nCols) = ps(1:nCols)*10.0_r8
    DO k = 1, kMax
       DO i = 1,nCols
          po  (i,k) = ps(i)*sl(k)*10.0_r8               ! pressure in mbar
          p   (i,k) = po(i,k)
          us  (i,k) = ua(i,k) 
          vs  (i,k) = va(i,k) 
          omeg(i,k) = omg(i,k)   !test 
          !
          !
          q   (i,k) = q2(i,k)
          qo  (i,k) = qn1(i,k)
          !----------------------------------------------------------q to r
          ! q  (i,k)=q2(i,k)/(1.0_r8-q2(i,k))
          ! qo(i,k)=qn1(i,k)/(1.0_r8-qn1(i,k))
          !-------------------------------------------------------------
          t   (i,k) = t2(i,k)
          tn  (i,k) = tn1(i,k)
          IF(TN(I,K) < 200.0_r8)    TN(I,K) = T(I,K)
          IF(QO(I,K) < 1.E-08_r8)   QO(I,K) = 1.E-08_r8
       END DO
    END DO
    !
    ! call cumulus parameterization
    !
    pre1   = 0.0_r8
    outt   = 0.0_r8
    outq   = 0.0_r8
    outqc = 0.0_r8
    dcape=0.0_r8
    !new-2007
    cupclw =0.0_r8  !new  cloud liquid water
    xmb=0.0_r8               !total base mass flux
    massfln=0.0_r8  !downdraft mass flux  
    !-----
    out1_mass_flux_up=0.0_r8
    out2_mass_flux_d=0.0_r8
    out3_detr=0.0_r8
    out4_liq_water=0.0_r8
    out5_qc=0.0_r8
    out6_wf=0.0_r8
    out7_entr=0.0_r8
    ddql=0.0_r8
    ddmu=0.0_r8
    ddmd=0.0_r8
    dden=0.0_r8
    ddde=0.0_r8
    ddca=0.0_r8
!---------------------------------------------
    pw=0.0_r8    ! precipitable water
    DO k=1,kmax
       DO i = 1,ncols
          pw(i) = pw(i) + delsig(k)*qo(i,k)
       END DO
    END DO
    DO i = 1,ncols
       pw(i)=100.0_r8*pw(i)*psur(i)/grav    !9.81_r8
       dd=cape(i)-cape_old(i)
       IF(cape(i)>0.0_r8.AND.dd>0.0_r8)THEN
       dcape(i)=dd
       ELSE
       dcape(i)=0.0_r8
       END IF
    END DO
!---------------------------
    CALL grellens2(&
         t      ,q      ,tn     ,qo      ,po     ,p      ,psur   ,us     ,vs    , &
         omeg   ,zz     ,dtime  ,nCols   ,kMax   ,xland  ,qliq   ,pw     ,lat2  , &
         lon    ,tke    ,dcape  ,cine    ,sens   ,ierr   ,outt   ,outq   ,outqc , &
         pre1   ,ktop   ,kbot   ,xmb     ,cupclw ,massfln,out1_mass_flux_up,    out2_mass_flux_d, &
         out3_detr,  out4_liq_water     , &
         out5_qc,    out6_wf ,out7_entr                                          )
    !
    ! after cumulus parameterization
    ! out  tn1, qn1, prec, kuo,ktop, kbot
    !
    DO i = 1,nCols
        IF(ierr(i) == 0)THEN
           RAINCV(i)=dtime*pre1(i)             !in mm/sec(ditme),by 0.5_r8(if leap-frog or 2dt)
           RAINCV(i)=RAINCV(i)/1000.0_r8       !in m for gcm
        END IF   
        IF(RAINCV(i) > 0.0_r8)kuo(i)=1
        kk=kbot(i)
        plcl(i)=p(i,kk)/10.0_r8  ! from mb to cb for Shallow convection
    END DO
!------------------------------------------
        ! out1_mass_flux_up(i,k)=zuo(i,k)  ! flux de massa normalizado
        ! out1_mass_flux_up(i,k)=xmb(i)*out1_mass_flux_up(i,k) Fluxo de massa  (onde xmb=fluxo de massa na basae).
        ! out2_mass_flux_d(i,k)=zdo(i,k)*edt(i) fluxo de massa descendente.normalizado..
        ! out2_mass_flux_d(i,k)=xmb(i)*out2_mass_flux_d(i,k)  fluxo de massa descendente (xmb= na base)
        ! out3_detr(i,k)=cd(i,k). detrainment
        ! qrc = liquid water content in cloud after rainout
        ! out4_liq_water(i,k)=qrc(i,k)      
        ! qc = cloud q (including liquid water) after entrainment
        ! out5_qc(i,k)=qc(i,k)
        ! out6=workFuntion
        ! out7-entrainment
!-----------------------------------------------
    DO k=1,kMax
       DO i=1,nCols   
          IF(ierr(i) == 0)THEN
             IF(RAINCV(i) > 0.0_r8)THEN
                 tn1(i,k) = tn1(i,k)+ 2.0_r8*outt(i,k)*dtime
                 !outq(i,k)=(outq(i,k)+outqc(i,k))/(1.0_r8+outq(i,k)+outqc(i,k))
                 qn1(i,k) = qn1(i,k)+ 2.0_r8*outq(i,k)*dtime
                 IF (cupclw(i,k).lt.1.0e-8_r8) cupclw(i,k)=0.0_r8  !new
             END IF
          END IF
       END DO
    END DO
!-------------------------------------------------up
       DO i=1,nCols
             IF(RAINCV(i) > 0.0_r8)THEN
                 ddca(i)=out6_wf(i)
                DO k=kbot(i),ktop(i)
                   ddmu(i,k)=xmb(i)*out1_mass_flux_up(i,k)
                   IF(ddmu(i,k).lt.1.0e-8_r8)ddmu(i,k)=0.0_r8

                   !!!!!ddql(i,k)=out5_qc(i,k)
                   ddql(i,k)=out4_liq_water(i,k)
                   IF (ddql(i,k).lt.1.0e-8_r8)ddql(i,k)=0.0_r8
                   dden(i,k)=MIN(MAX(0.0_r8,out7_entr(i,k)),10.0_r8)
                   ddde(i,k)=MIN(MAX(0.0_r8,out3_detr(i,k)),10.0_r8)
                END DO
              
                DO k=1,ktop(i)
                  ddmd(i,k)=xmb(i)*out2_mass_flux_d(i,k)
                  IF(ddmd(i,k).lt.1.0e-8_r8)ddmd(i,k)=0.0_r8
               END DO

             END IF
       END DO

      !---------------------------------------------------
      !   SALVAR   SO QUANDO CHOVE  (RAINCV(i)>=0.0_r8)
      !-------------------------------------------------
      ! CAPE E CINE vem de FORA
      !xmb     total base mass flux
      !------
      DO i=1,nCols
         IF(RAINCV(i)<=0.0_r8)THEN
            !cine(i)=undef
            ddca(i)=undef
            xmb(i)=undef
            ddmu(i,:)=undef
            ddmd(i,:)=undef
            ddde(i,:)=undef
            dden(i,:)=undef
            ddql(i,:)=undef
          END IF
      END DO
!      if (dodia(nDiag_deepcape))  call updia(ddca                 ,nDiag_deepcape   ,latco,jdt)
!      !if (dodia(nDiag_deepxx1d)) call updia(xmb                   ,nDiag_deepxx1d ,latco,jdt)
!      if (dodia(nDiag_deepxx1d)) call updia(cine                   ,nDiag_deepxx1d ,latco,jdt)
!
!      if (dodia(nDiag_deepupm))  call updia(ddmu                  ,nDiag_deepupm  ,latco,jdt) !!!!!!IMPOETANT
!      if (dodia(nDiag_deepdwm))  call updia(ddmd                  ,nDiag_deepdwm  ,latco,jdt)
!      if (dodia(nDiag_deepdetr)) call updia(ddde                  ,nDiag_deepdetr ,latco,jdt)
!      if (dodia(nDiag_deepentr))  call updia(dden                 ,nDiag_deepentr   ,latco,jdt)
!      !if (dodia(nDiag_deepql ))  call updia(ddql                  ,nDiag_deepql   ,latco,jdt)
!      if (dodia(nDiag_deepql ))  call updia(omeg                  ,nDiag_deepql   ,latco,jdt)

    RETURN
  END SUBROUTINE grellens


  !*-------------------------

  SUBROUTINE grellens2(& 
       t      ,q      ,tn     ,qo     ,po     ,p      ,psur   ,us     ,vs       , &
       omeg   ,z1     ,dtime  ,nCols  ,kMax   ,mask   ,qliq   ,pwt    ,lat     , &
       lon    ,tke    ,dcape  ,cine   ,sens   ,ierr   ,outt   ,outq   ,outqc   , &
       pre    ,ktop   ,kbcon  ,xmb    ,cupclw ,massfln,out1_mass_flux_up,   out2_mass_flux_d, &
       out3_detr,  out4_liq_water     ,&
       out5_qc, out6_wf,out7_entr)

    IMPLICIT NONE
    !
    ! input variables
    ! IN
    INTEGER         , INTENT(IN   )       :: nCols
    INTEGER         , INTENT(IN   )       :: kMax
    REAL(KIND=r8)   , INTENT(IN   )       :: dtime
    REAL(KIND=r8)   , INTENT(IN   )       :: pwt  (nCols)
    REAL(KIND=r8)   , INTENT(IN   )       :: lat  (nCols)
    REAL(KIND=r8)   , INTENT(IN   )       :: lon  (nCols)
    REAL(KIND=r8)   , INTENT(IN   )       :: t    (nCols,kMax)
    REAL(KIND=r8)   , INTENT(INOUT)       :: q    (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )       :: tn   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(INOUT)       :: qo   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )       :: p    (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )       :: po   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )       :: us   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )       :: vs   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )       :: omeg (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )       :: z1   (nCols)
    REAL(KIND=r8)   , INTENT(IN   )       :: tke  (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )       :: dcape(nCols)
    REAL(KIND=r8)   , INTENT(INOUT)       :: cine (nCols)
    REAL(KIND=r8)   , INTENT(IN   )       :: sens (nCols)
    REAL(KIND=r8)   , INTENT(IN   )       :: psur (nCols)
    INTEGER         , INTENT(IN   )       :: mask (nCols)
    REAL(KIND=r8)   , INTENT(OUT  )       :: qliq(nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  )       :: xmb  (nCols)
    INTEGER         , INTENT(INOUT)       :: ierr      (nCols)
!----
    REAL(KIND=r8)   , INTENT(INOUT)       :: massfln(nCols,ensdim) !new
    REAL(KIND=r8)   , INTENT(INOUT)       :: cupclw (nCols,kMax)
    REAL(KIND=r8)                         :: buoy(nCols,kMax)
    
    !
    ! output variables
    ! OUT


    REAL(KIND=r8)   , INTENT(INOUT)       :: outt (nCols,kMax)
    REAL(KIND=r8)   , INTENT(INOUT)       :: outq (nCols,kMax)
    REAL(KIND=r8)   , INTENT(INOUT)       :: outqc(nCols,kMax)
    REAL(KIND=r8)   , INTENT(INOUT)       :: pre  (nCols)
    INTEGER         , INTENT(INOUT)       :: ktop (nCols)
    INTEGER         , INTENT(INOUT)       :: kbcon(nCols)

  !new-dec2011
!----------------------------------------------------------------------------------
    REAL(KIND=r8)   , INTENT(OUT  )      :: out1_mass_flux_up(nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  )      :: out2_mass_flux_d(nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  )      :: out3_detr(nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  )      :: out4_liq_water(nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  )      :: out5_qc(nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  )      :: out6_wf(nCols)
    REAL(KIND=r8)   , INTENT(OUT  )      :: out7_entr(nCols, kMax)

    !
    ! LOCAL VARIABLES
    !
    INTEGER            :: i
    INTEGER            :: k
    INTEGER            :: j
!not used    INTEGER            :: iedt
    INTEGER            :: istart 
    INTEGER            :: iend

    INTEGER            :: kdet1     (nCols)
    INTEGER            :: kdet      (nCols)
    REAL(KIND=r8)      :: mconv     (nCols)
    INTEGER            :: kzdown    (nCols)
    INTEGER            :: kbmax     (nCols)
    INTEGER            :: k22       (nCols)
    INTEGER            :: jmin      (nCols)
    INTEGER            :: kstabi    (nCols)
    INTEGER            :: kstabm    (nCols)
    INTEGER            :: KZI       (nCols)

    REAL(KIND=r8)               :: aaeq      (nCols)
    REAL(KIND=r8)               :: edt       (nCols)
    REAL(KIND=r8)               :: aa1       (nCols)
    REAL(KIND=r8)               :: aa0       (nCols)
    REAL(KIND=r8)               :: hkb       (nCols)
    REAL(KIND=r8)               :: hkbo      (nCols)
    !! snf REAL(KIND=r8)               :: xmb       (nCols)
    REAL(KIND=r8)               :: pwav      (nCols)
    REAL(KIND=r8)               :: pwev      (nCols)
    REAL(KIND=r8)               :: pwavo     (nCols)
    REAL(KIND=r8)               :: pwevo     (nCols)
    REAL(KIND=r8)               :: bu        (nCols)
    REAL(KIND=r8)               :: cap_max   (nCols)
    REAL(KIND=r8)               :: vshear    (nCols)
    REAL(KIND=r8)               :: sdp       (nCols)
    REAL(KIND=r8)               :: vws       (nCols)
    REAL(KIND=r8)               :: he        (nCols,kMax)
    REAL(KIND=r8)               :: hes       (nCols,kMax)
    REAL(KIND=r8)               :: qes       (nCols,kMax)
    REAL(KIND=r8)               :: z         (nCols,kMax)
    REAL(KIND=r8)               :: dby       (nCols,kMax)
    REAL(KIND=r8)               :: qc        (nCols,kMax)
    REAL(KIND=r8)               :: qrcd      (nCols,kMax)
    REAL(KIND=r8)               :: pwd       (nCols,kMax)
    REAL(KIND=r8)               :: pw        (nCols,kMax)
    REAL(KIND=r8)               :: heo       (nCols,kMax)
    REAL(KIND=r8)               :: heso      (nCols,kMax)
    REAL(KIND=r8)               :: qeso      (nCols,kMax)
    REAL(KIND=r8)               :: zo        (nCols,kMax)
    REAL(KIND=r8)               :: dbyo      (nCols,kMax)
    REAL(KIND=r8)               :: qco       (nCols,kMax)
    REAL(KIND=r8)               :: qrcdo     (nCols,kMax)
    REAL(KIND=r8)               :: pwdo      (nCols,kMax)
    REAL(KIND=r8)               :: pwo       (nCols,kMax)
    REAL(KIND=r8)               :: hcd       (nCols,kMax)
    REAL(KIND=r8)               :: hcdo      (nCols,kMax)
    REAL(KIND=r8)               :: qcd       (nCols,kMax)
    REAL(KIND=r8)               :: qcdo      (nCols,kMax)
    REAL(KIND=r8)               :: dbyd      (nCols,kMax)
    REAL(KIND=r8)               :: dbydo     (nCols,kMax)
    REAL(KIND=r8)               :: hc        (nCols,kMax)
    REAL(KIND=r8)               :: hco       (nCols,kMax)
    REAL(KIND=r8)               :: qrc       (nCols,kMax)
    REAL(KIND=r8)               :: qrco      (nCols,kMax)
    REAL(KIND=r8)               :: zu        (nCols,kMax)
    REAL(KIND=r8)               :: zuo       (nCols,kMax)
    REAL(KIND=r8)               :: zd        (nCols,kMax)
    REAL(KIND=r8)               :: zdo       (nCols,kMax)
    REAL(KIND=r8)               :: qes_cup   (nCols,kMax)
    REAL(KIND=r8)               :: q_cup     (nCols,kMax)
    REAL(KIND=r8)               :: he_cup    (nCols,kMax)
    REAL(KIND=r8)               :: hes_cup   (nCols,kMax)
    REAL(KIND=r8)               :: z_cup     (nCols,kMax)
    REAL(KIND=r8)               :: p_cup     (nCols,kMax)
    REAL(KIND=r8)               :: gamma_cup (nCols,kMax)
    REAL(KIND=r8)               :: t_cup     (nCols,kMax)
    REAL(KIND=r8)               :: qeso_cup  (nCols,kMax)
    REAL(KIND=r8)               :: qo_cup    (nCols,kMax)
    REAL(KIND=r8)               :: heo_cup   (nCols,kMax)
    REAL(KIND=r8)               :: heso_cup  (nCols,kMax)
    REAL(KIND=r8)               :: zo_cup    (nCols,kMax)
    REAL(KIND=r8)               :: po_cup    (nCols,kMax)
    REAL(KIND=r8)               :: gammao_cup(nCols,kMax)
    REAL(KIND=r8)               :: tn_cup    (nCols,kMax)
    REAL(KIND=r8)               :: cd        (nCols,kMax)
    REAL(KIND=r8)               :: cdd       (nCols,kMax)

    REAL(KIND=r8)               :: dellat_ens (nCols,kMax, maxens2)
    REAL(KIND=r8)               :: dellaq_ens (nCols,kMax, maxens2)
    REAL(KIND=r8)               :: dellaqc_ens(nCols,kMax, maxens2)
    REAL(KIND=r8)               :: pwo_ens    (nCols,kMax, maxens2)
    REAL(KIND=r8)               :: xf_ens     (nCols,ensdim)
    REAL(KIND=r8)               :: outt_ens   (nCols,ensdim)
    REAL(KIND=r8)               :: pr_ens     (nCols,ensdim)
    !!!snf REAL(KIND=r8)               :: massfln    (nCols,ensdim)

    REAL(KIND=r8)               :: edtc     (nCols,maxens2)

    REAL(KIND=r8)               :: dq 
    REAL(KIND=r8)               :: mbdt
    REAL(KIND=r8)               :: zktop
    !-------
    !new
    !---------
    REAL(KIND=r8)               :: dh2         (nCols) 
    REAL(KIND=r8)               :: xfac1       (nCols) 
    REAL(KIND=r8)               :: xfac_for_dn (nCols)
    INTEGER                     :: left        (nCols)
    INTEGER                     :: nLeft,ib   !not used ,nNewLeft
    INTEGER                     :: maxens22
!---------------------------------
!2011
    INTEGER                     :: kk1     (nCols),kk
    REAL(KIND=r8)               :: tke1D   (kMax)
    REAL(KIND=r8)               :: tkeMax  (nCols)
    REAL(KIND=r8)               :: tkeMedia(nCols)
    REAL(KIND=r8)               :: den0    (nCols,kMax)
    REAL(KIND=r8)               :: den1    (nCols,kMax)
    REAL(KIND=r8)               :: da,dz
    REAL(KIND=r8)               :: hh      (nCols)

!----------------new 2012
    REAL(KIND=r8)               :: HR    (nCols,kMax)
    REAL(KIND=r8)               :: bb,cc,dd,ee,ff
    REAL(KIND=r8)               :: entr2D(nCols,kMax)
    REAL(KIND=r8)               :: detr2D(nCols,kMax)
    REAL(KIND=r8)               :: entr2


    !
    ! Compress Local Variable 
    !
    INTEGER :: nCols_gz
    REAL(KIND=r8)    :: edtc_gz       (nCols,maxens2)
    INTEGER :: ierr_gz       (nCols)
    REAL(KIND=r8)    :: dellat_ens_gz (nCols,kMax, maxens2)
    REAL(KIND=r8)    :: dellaq_ens_gz (nCols,kMax, maxens2)
    REAL(KIND=r8)    :: dellaqc_ens_gz(nCols,kMax, maxens2)
    REAL(KIND=r8)    :: pwo_ens_gz    (nCols,kMax, maxens2)
    REAL(KIND=r8)    :: heo_cup_gz    (nCols,kMax)
    REAL(KIND=r8)    :: zo_cup_gz     (nCols,kMax)
    REAL(KIND=r8)    :: po_cup_gz     (nCols,kMax)
    REAL(KIND=r8)    :: hcdo_gz       (nCols,kMax)
    REAL(KIND=r8)    :: zdo_gz        (nCols,kMax)
    REAL(KIND=r8)    :: cdd_gz        (nCols,kMax)
    REAL(KIND=r8)    :: heo_gz        (nCols,kMax)
    REAL(KIND=r8)    :: qo_cup_gz     (nCols,kMax)
    REAL(KIND=r8)    :: qrcdo_gz      (nCols,kMax)
    REAL(KIND=r8)    :: qo_gz         (nCols,kMax)
    REAL(KIND=r8)    :: zuo_gz        (nCols,kMax)
    REAL(KIND=r8)    :: cd_gz         (nCols,kMax)
    REAL(KIND=r8)    :: hco_gz        (nCols,kMax)
    INTEGER :: ktop_gz       (nCols)
    INTEGER :: k22_gz        (nCols)
    INTEGER :: kbcon_gz      (nCols)
    INTEGER :: jmin_gz       (nCols)
    INTEGER :: kdet_gz       (nCols)
    REAL(KIND=r8)    :: qco_gz        (nCols,kMax)
    REAL(KIND=r8)    :: qrco_gz       (nCols,kMax)
    REAL(KIND=r8)    :: tn_gz         (nCols,kMax)
    REAL(KIND=r8)    :: po_gz         (nCols,kMax)
    REAL(KIND=r8)    :: z1_gz         (nCols)
    REAL(KIND=r8)    :: psur_gz       (nCols)
    REAL(KIND=r8)    :: gamma_cup_gz  (nCols,kMax)
    REAL(KIND=r8)    :: pr_ens_gz     (nCols,ensdim)
    REAL(KIND=r8)    :: pwo_gz        (nCols,kMax)
    REAL(KIND=r8)    :: pwdo_gz       (nCols,kMax)
    REAL(KIND=r8)    :: outt_ens_gz   (nCols,ensdim)    
    REAL(KIND=r8)    :: he_cup_gz     (nCols,kMax)
    INTEGER :: kbmax_gz      (nCols)
    REAL(KIND=r8)    :: heso_cup_gz   (nCols,kMax)
    REAL(KIND=r8)    :: cap_max_gz    (nCols)
    REAL(KIND=r8)    :: aa0_gz        (nCols)
    REAL(KIND=r8)    :: aa1_gz        (nCols)
    REAL(KIND=r8)    :: xmb_gz        (nCols)
    REAL(KIND=r8)    :: xf_ens_gz     (nCols,ensdim)
    INTEGER :: mask_gz       (nCols)
    REAL(KIND=r8)    :: mconv_gz      (nCols)
    REAL(KIND=r8)    :: omeg_gz       (nCols,kMax)
    REAL(KIND=r8)    :: massfln_gz    (nCols,ensdim)
    REAL(KIND=r8)    :: p_cup_gz      (nCols,kMax)
    INTEGER :: listim  (nCols)
    LOGICAL :: bitx    (nCols)
    INTEGER :: litx    (nCols)
    integer , parameter :: i_use_stat_prop=0
!new 2012
    REAL(KIND=r8)    :: entr2D_gz   (nCols,kMax)
    REAL(KIND=r8)    :: detr2D_gz   (nCols,kMax)
    REAL(KIND=r8)    :: dcape_gz      (nCols)
    REAL(KIND=r8)    :: cine_gz      (nCols)
    REAL(KIND=r8)    :: sens_gz      (nCols)
    REAL(KIND=r8)    :: tkeMax_gz     (nCols)
    REAL(KIND=r8)    :: tkeMedia_gz     (nCols)
    REAL(KIND=r8)    :: den1_gz       (nCols,kMax)

!
!***************** the following are your basic environmental
!                  variables. They carry a "_cup" if they are
!                  on model cloud levels (staggered). They carry
!                  an "o"-ending (z becomes zo), if they are the forced
!                  variables. They are preceded by x (z becomes xz)
!                  to indicate modification by some typ of cloud
!
  ! z           = heights of model levels
  ! q           = environmental mixing ratio
  ! qes         = environmental saturation mixing ratio
  ! t           = environmental temp
  ! p           = environmental pressure
  ! he          = environmental moist static energy
  ! hes         = environmental saturation moist static energy
  ! z_cup       = heights of model cloud levels
  ! q_cup       = environmental q on model cloud levels
  ! qes_cup     = saturation q on model cloud levels
  ! t_cup       = temperature (Kelvin) on model cloud levels
  ! p_cup       = environmental pressure
  ! he_cup = moist static energy on model cloud levels
  ! hes_cup = saturation moist static energy on model cloud levels
  ! gamma_cup = gamma on model cloud levels
!
!
  ! hcd = moist static energy in downdraft
  ! zd normalized downdraft mass flux
  ! dby = buoancy term
  ! entr = entrainment rate
  ! zd   = downdraft normalized mass flux
  ! entr= entrainment rate
  ! hcd = h in model cloud
  ! bu = buoancy term
  ! zd = normalized downdraft mass flux
  ! gamma_cup = gamma on model cloud levels
  ! mentr_rate = entrainment rate
  ! qcd = cloud q (including liquid water) after entrainment
  ! qrch = saturation q in cloud
  ! pwd = evaporate at that level
  ! pwev = total normalized integrated evaoprate (I2)
  ! entr= entrainment rate
  ! z1 = terrain elevation
  ! entr = downdraft entrainment rate
  ! jmin = downdraft originating level
  ! kdet = level above ground where downdraft start detraining
  ! psur        = surface pressure
  ! z1          = terrain elevation
  ! pr_ens = precipitation ensemble
  ! xf_ens = mass flux ensembles
  ! massfln = downdraft mass flux ensembles used in next timestep
  ! omeg = omega from large scale model
  ! mconv = moisture convergence from large scale model
  ! zd      = downdraft normalized mass flux
  ! zu      = updraft normalized mass flux
  ! dir     = "storm motion"
  ! mbdt    = arbitrary numerical parameter
  ! dtime   = dt over which forcing is applied
  ! iact_gr_old = flag to tell where convection was active
  ! kbcon       = LFC of parcel from k22
  ! k22         = updraft originating level
  ! icoic       = flag if only want one closure (usually set to zero!)
  ! dby = buoancy term
  ! ktop = cloud top (output)
  ! xmb    = total base mass flux
  ! hc = cloud moist static energy
  ! hkb = moist static energy at originating level
  ! mentr_rate = entrainment rate

    !
    !snf parameter from namelist 
    !
    ! begin executable
    !
    qliq    =0.0_r8
    outqc   =0.0_r8
    kzdown  =0
    kbmax   =0
    k22     =0
    jmin    =0
    kstabi  =0
    kstabm  =0
    KZI     =0
    left    =0
    aaeq    =0.0_r8
    edt     =0.0_r8
    aa1     =0.0_r8
    aa0     =0.0_r8
    hkb     =0.0_r8
    hkbo    =0.0_r8
    xmb     =0.0_r8
    pwav    =0.0_r8
    pwev    =0.0_r8
    pwavo   =0.0_r8
    pwevo   =0.0_r8
    bu      =0.0_r8
    cap_max =0.0_r8
    vshear  =0.0_r8
    sdp     =0.0_r8
    vws     =0.0_r8
    he      =0.0_r8
    hes     =0.0_r8
    qes     =0.0_r8
    z       =0.0_r8
    dby     =0.0_r8
    qc      =0.0_r8
    qrcd    =0.0_r8
    pwd     =0.0_r8
    pw      =0.0_r8
    heo     =0.0_r8
    heso    =0.0_r8
    qeso    =0.0_r8
    zo      =0.0_r8
    dbyo    =0.0_r8
    qco     =0.0_r8
    qrcdo  =0.0_r8
    pwdo  =0.0_r8
    pwo  =0.0_r8
    hcd  =0.0_r8
    hcdo  =0.0_r8
    qcd  =0.0_r8
    qcdo  =0.0_r8
    dbyd  =0.0_r8
    dbydo  =0.0_r8
    hc  =0.0_r8
    hco  =0.0_r8
    qrc  =0.0_r8
    qrco  =0.0_r8
    zu  =0.0_r8
    zuo  =0.0_r8
    zd  =0.0_r8
    zdo  =0.0_r8
    qes_cup   =0.0_r8
    q_cup  =0.0_r8
    he_cup    =0.0_r8
    hes_cup   =0.0_r8
    z_cup  =0.0_r8
    p_cup  =0.0_r8
    gamma_cup =0.0_r8
    t_cup  =0.0_r8
    qeso_cup  =0.0_r8
    qo_cup    =0.0_r8
    heo_cup   =0.0_r8
    heso_cup  =0.0_r8
    zo_cup    =0.0_r8
    po_cup    =0.0_r8
    gammao_cup=0.0_r8
    tn_cup    =0.0_r8
    cd  =0.0_r8
    cdd  =0.0_r8
    dellat_ens  =0.0_r8
    dellaq_ens  =0.0_r8
    dellaqc_ens =0.0_r8
    pwo_ens  =0.0_r8   
    xf_ens   =0.0_r8   
    outt_ens  =0.0_r8  
    pr_ens   =0.0_r8   
    edtc =0.0_r8
    dq =0.0_r8
    mbdt=0.0_r8
    zktop=0.0_r8



    dh2   =0.0_r8
    xfac1   =0.0_r8
    xfac_for_dn=0.0_r8
    dellat_ens =0.0_r8
    dellaq_ens =0.0_r8
    dellaqc_ens=0.0_r8
    pwo_ens    =0.0_r8
    hkbo=0
    istart=1
    iend=nCols
    maxens22=grepar2
    !
    !snf  is it necessary to save massfln for next step?.
    ! no it is not necessary
    !
    !!!! snf massfln=0.0_r8
    mconv=0.0_r8
    qrco=0.0_r8
    qrco_gz=0.0_r8
    DO i=istart,iend
       !
       ! prepare input, erase output
       !
       kdet  (i) =2
       kdet1 (i) =0
       pre   (I) =0.0_r8
    END DO
    !
    ! calculate moisture convergence mconv
    !
    DO k=2,kMax-1
       DO i = istart,iend
          dq      = 0.5_r8*(q(i,k+1)-q(i,k-1))
          mconv(i) = mconv(i) + omeg(i,k)*dq/9.81_r8
       END DO
    END DO

    DO I = istart,iend
       IF(mconv(I) < 0.0_r8)  mconv(I) = 0.0_r8
    END DO
    !
    ! initial detrainmentrates
    !
    DO k=1,kMax
       DO i=istart,iend
          !
          ! snf with 0.5_r8 and 0.1_r8 does not difference..why?? 
          !
          cd (i,k) = 0.1_r8*entr_rate             !!!new2
          cdd(i,k) = 0.0_r8
       END DO
    END DO

    DO i=istart,iend
       aa0   (i)=0.0_r8                                   !snf1
       aa1   (i)=0.0_r8
       hh(i)=0.0_r8
       cine (i)=0.0_r8
       kstabm(i)=kMax-2
       aaeq  (i)=0.0_r8                            !added snf  
       IF(aaeq(i) <  0.0_r8)THEN
          ierr(i)=20
       ELSE
          ierr (i)=0
       END IF
       !nilo 2012
       tkeMedia(i)=0.0_r8
    END DO
    !
    !snf3
    !--- initialize cap_max
    !
    DO i=istart,iend
       cap_max(i)=cap_maxs
       IF(mask(i).NE.1)cap_max(i)=cap_maxs_land
    END DO
    mbdt=(REAL(1,kind=r8)-3.0_r8)*dtime*1.e-3_r8 + dtime*5.e-03_r8  !new
    !
    ! environmental conditions, FIRST HEIGHTS
    !
    DO k=1,maxens*maxens22*maxens3
       DO i=istart,iend
          IF(ierr(i).NE.20)THEN
             xf_ens  (i,(iens-1)*maxens*maxens22*maxens3+k)= 0.0_r8
             pr_ens  (i,(iens-1)*maxens*maxens22*maxens3+k)= 0.0_r8
             outt_ens(i,(iens-1)*maxens*maxens22*maxens3+k)= 0.0_r8
          END IF
       END DO
    END DO
    !
    ! calculate moist static energy, heights, qes
    !
    CALL cup_env(z      , & ! z      (out)
         qes    , & ! qes    (out)
         he     , & ! he     (inout)
         hes    , & ! hes    (out)
         t      , & ! t      (in)
         q      , & ! q      (inout)
         p      , & ! p      (in)
         z1     , & ! z1     (in)
         nCols  , & ! nCols  (in)
         kMax   , & ! kMax   (in)
         istart , & ! istart (in)
         iend   , & ! iend   (in)
         psur   , & ! psur   (in)
         ierr   , & ! ierr   (in)
         0      ,&  !tcrit  (in)
         den0     ) ! den    (out)

    CALL cup_env(zo     , & ! zo     (out)
         qeso   , & ! qeso   (out)
         heo    , & ! heo    (inout)
         heso   , & ! heso   (out)
         tn     , & ! tn     (in)
         qo     , & ! qo     (inout)
         po     , & ! po     (in)
         z1     , & ! z1     (in)
         nCols  , & ! nCols  (in)
         kMax   , & ! kMax   (in)
         istart , & ! istart (in)
         iend   , & ! iend   (in)
         psur   , & ! psur   (in)
         ierr   , & ! ierr   (in)
         0      , & ! tcrit  (in)
         den1      ) ! den    (out)
    !
    ! environmental values on cloud levels
    !
    CALL cup_env_clev(t        , & ! t         (in)
         qes      , & ! qes       (in)
         q        , & ! q         (in)
         he       , & ! he        (in)
         hes      , & ! hes       (in)
         z        , & ! z         (in)
         p        , & ! p         (in)
         qes_cup  , & ! qes_cup   (out)
         q_cup    , & ! q_cup     (out)
         he_cup   , & ! he_cup    (out)
         hes_cup  , & ! hes_cup   (out)
         z_cup    , & ! z_cup     (out)
         p_cup    , & ! p_cup     (out)
         gamma_cup, & ! gamma_cup (out)
         t_cup    , & ! t_cup     (out)
         psur     , & ! psur      (in)
         nCols    , & ! nCols     (in)
         kMax     , & ! kMax      (in)
         istart   , & ! istart    (in)
         iend     , & ! iend      (in)
         ierr     , & ! ierr      (in)
         z1         ) ! z1        (in)

    CALL cup_env_clev(tn        , &! tn        (in)
         qeso      , &! qeso      (in)
         qo        , &! qo        (in)
         heo       , &! heo       (in)
         heso      , &! heso      (in)
         zo        , &! zo        (in)
         po        , &! po        (in)
         qeso_cup  , &! qeso_cup  (out)
         qo_cup    , &! qo_cup    (out)
         heo_cup   , &! heo_cup   (out)
         heso_cup  , &! heso_cup  (out)
         zo_cup    , &! zo_cup    (out)
         po_cup    , &! po_cup    (out)
         gammao_cup, &! gammao_cup(out)
         tn_cup    , &! tn_cup    (out)
         psur      , &! psur      (in)
         nCols     , &! nCols     (in)
         kMax      , &! kMax      (in)
         istart    , &! istart    (in)
         iend      , &! iend      (in)
         ierr      , &! ierr      (in)
         z1          )! z1        (in)
    !
    !
    !
    kbmax=0
    DO k=1,kMax
       DO i=istart,iend
          IF(ierr(i) == 0 .AND. zo_cup(i,k) >  zkbmax+z1(i) .AND. kbmax(i) ==0)THEN
             kbmax(i)=k
          END IF
       END DO
    END DO
    !
    ! level where detrainment for downdraft starts
    !
    kdet1=0
    DO k=1,kMax
       DO i=istart,iend
          IF(ierr(i) == 0 .AND. zo_cup(i,k) >  z_detr+z1(i) .AND.kdet1(i) ==0)THEN
             kdet (i)=k
             kdet1(i)=k
          END IF
       END DO
    END DO
    !
    !-------------------------------
    ! USE TKE  (it is still no used)
    !---------------------------------
    ! Determine PBL top using TKE (TKEG) and 
    ! liquid water mixing ratio (RCPG)
    !------------------------------------
    !!       call get_zi(nCols,kMax,istart,iend,j,ierr,kzi,TKEG &
    !!                   ,RCPG,zo,z1,tkmin)
    !
    !       DO  I=ISTART,IEND   
    !         IF(ierr(I) == 0)THEN
    !          tkemax(i) = 0.0_r8
    !          do k=1,kzi(i)
    !           tkemax(i) = max(tkemax(i),tkeg(i,k))
    !          enddo
    !          if(tkemax(i) < 1.5_r8 )  cap_max(i) = 75.0_r8  
    !          if(tkemax(i) < 0.5_r8 )  cap_max(i) = 25.0_r8
    !         endif
    !       enddo
    !------------------------------------
    !
    ! determine level with highest moist static energy content - k22
    ! kstart = 3 
    !
    CALL cup_maximi(heo_cup  , &  ! heo_cup (in)
         nCols    , &  ! nCols   (in)
         kMax     , &  ! kMax    (in) 
         3        , &  ! ks      (in) !era 3  
         kbmax    , &  ! kbmax   (in)
         k22      , &  ! k22     (out)
         istart   , &  ! istart  (in)
         iend     , &  ! iend    (in)
         ierr       )  ! ierr    (in)
    !                
    DO i=istart,iend
       IF(ierr(i) == 0)THEN
          kzi(i) = 1
          IF(k22(i) >= kbmax(i))ierr(i)=2
       END IF
    END DO
    !
    ! determine the level of convective cloud base  - kbcon
    !
    !snf8 cup_KBCOn is called first timE
    !  Cap_mas for 
    !-------------------------------------
    !snf  call cup_kbcon for cap_max=cap_max-(1-1)*cap_max_increment
    !---------------
    !
    CALL cup_kbcon(&
         cap_max_increment, & ! cap_max_increment (in)
         1                , & ! iloop             (in)
         k22              , & ! k22               (inout)
         kbcon            , & ! kbcon             (out)
         heo_cup          , & ! heo_cup           (in)
         heso_cup         , & ! heso_cup          (in)
         nCols            , & ! nCols             (in)
         kMax             , & ! kMax              (in)
         istart           , & ! istart            (in)
         iend             , & ! iend              (in)
         ierr             , & ! ierr              (inout)
         kbmax            , & ! kbmax             (in)
         po_cup           , & ! po_cup            (in)
         cap_max            ) ! cap_max           (in)
    DO I=ISTART,IEND
       IF(ierr(I) == 0)THEN
          hkb(i)=hkbo(i)
       END IF
    END DO
!----------------------------------new
         tkeMax=0.0_r8
    DO i=istart,iend
        tke1D=0.0_r8
       IF(ierr(I) == 0)THEN
          DO k=2,kbcon(i)
          tke1D(k)=tke(i,k)
          END DO
          tkeMax(i)=MAXVAL(tke1D)
        END IF
    END DO
!-------------------------------------------
! calculate PBL TOP level
!-----------------------------
    DO i=istart,iend
       kk1(i)=kbcon(i)
       DO k=2,kbcon(i)
          IF(tke(i,k)<=0.17_r8)THEN   !!
            kk1(i)=k
            EXIT
          END IF
       END DO
    END DO
!---------------------------------------

    tkeMedia=0.0_r8
    DO I=ISTART,IEND
       IF(ierr(I) == 0)THEN
          DO k=1,kbcon(i)
           tkeMedia(i)=tkeMedia(i)+tke(i,k)
          ENDDO
           if(kbcon(i)>=1)tkeMedia(i)=tkeMedia(i)/float(kbcon(i))
       END IF
    END DO
!---------------------------------------------------------------
    !
    ! increase detrainment in stable layers
    !
    CALL cup_minimi( &
         heso_cup , &  ! heso_cup (in)
         nCols    , &  ! nCols    (in)
         kMax     , &  ! kMax     (in)
         kbcon    , &  ! kbcon    (in)
         kstabm   , &  ! kstabm   (in)
         kstabi   , &  ! kstabi   (out)
         istart   , &  ! istart   (in)
         iend     , &  ! iend     (in)
         ierr       )  ! ierr     (in)
!--------------------------------------
!--------------------------
! define entr and detr rate    nilo1
!---------------------------
    entr2D=0.0_r8
    detr2D=0.0_r8
    entr2=0.0_r8
    entr2D=entr_rate
    detr2D=mentr_rate
    entr2=entr_rate
!--------------------------------------
    DO k=1,kMax
       DO i=istart,iend
          IF( ierr(i) == 0 .AND. kstabm(i)-1 >= kstabi(i) .AND. &
               k >= kstabi(i) .AND. k<= kstabm(i)-1 )THEN
                 entr2=entr2D(i,k)
             cd(i,k)=cd(i,k-1)+1.5_r8*entr2  
             IF(iens >  4)THEN
                cd(i,k)=cd(i,k-1)+REAL(iens-4,kind=r8)*entr2&
                     /REAL(kstabm(i)-kstabi(i),kind=r8)
             ELSE
                cd(i,k)=cd(i,k)
             END IF
             IF(cd(i,k) >  10.0_r8*entr2) cd(i,k)=10.0_r8*entr2 !new
          END IF
       END DO
    END DO
    !
    ! calculate incloud moist static energy
    !
        CALL cup_up_he(  &
         k22       , & ! k22        (in)
         hkb       , & ! hkb        (out)
         z_cup     , & ! z_cup      (in)
         cd        , & ! cd         (in)
         mentr_rate, & ! mentr_rate (in)
         he_cup    , & ! he_cup     (in)
         hc        , & ! hc         (out)
         nCols     , & ! nCols      (in)
         kMax      , & ! kMax       (in)
         kbcon     , & ! kbcon      (in)
         ierr      , & ! ierr       (in)
         istart    , & ! istart     (in)
         iend      , & ! iend       (in)
         dby       , & ! dby        (out)
         he        , & ! he         (in)
         hes_cup   , & ! ) ! hes_cup    (in)
         entr2D)
    !
    CALL cup_up_he(  &
         k22       , & ! k22        (in)
         hkbo      , & ! hkbo       (out)
         zo_cup    , & ! zo_cup     (in)
         cd        , & ! cd         (in)
         mentr_rate, & ! mentr_rate (in)
         heo_cup   , & ! heo_cup    (in)
         hco       , & ! hco        (out)
         nCols     , & ! nCols      (in)
         kMax      , & ! kMax       (in)
         kbcon     , & ! kbcon      (in)
         ierr      , & ! ierr       (in)
         istart    , & ! istart     (in)
         iend      , & ! iend       (in)
         dbyo      , & ! dbyo       (out)
         heo       , & ! heo        (in)
         heso_cup  , & !
         entr2D)
    !
    ! determine cloud top - ktop
    !
    CALL cup_ktop(  &
         1        , & ! ilo    (in)
         dbyo     , & ! dbyo   (inout)
         kbcon    , & ! kbcon  (in)
         ktop     , & ! ktop   (out)
         nCols    , & ! nCols  (in)
         kMax     , & ! kMax   (in)
         istart   , & ! istart (in)
         iend     , & ! iend   (in)
         ierr       ) ! ierr   (inout)

    kzdown(istart:iend)=0
    DO k=1,kMax
       DO i=istart,iend
          IF(ierr(i) == 0)THEN
             zktop=(zo_cup(i,ktop(i))-z1(i))*0.6_r8
             zktop=MIN(zktop+z1(i),zcutdown+z1(i))
             IF(zo_cup(i,k) >  zktop .AND. kzdown(i) == 0 )THEN
                kzdown(i) = k
             END IF
          END IF
       END DO
    END DO
    !
    ! downdraft originating level - jmin
    ! jmin output from cup_minimi
    !
    CALL cup_minimi( &
         heso_cup  , &! heso_cup (in)
         nCols     , &! nCols    (in)
         kMax      , &! kMax     (in)
         k22       , &! k22      (in)
         kzdown    , &! kzdown   (in)
         jmin      , &! jmin     (out)
         istart    , &! istart   (in)
         iend      , &! iend     (in)
         ierr        )! ierr     (in)
    !
    ! check whether it would have buoyancy, if there where
    ! no entrainment/detrainment
    !
    dh2(istart:iend)=0.0_r8    
    nLeft = 0
    DO i=istart,iend
       IF(ierr(i) == 0)THEN
          nLeft = nLeft + 1
          left(nLeft) = i
          IF (jmin(i)-1 <  kdet(i)) kdet(i)=jmin(i)-1
          IF (jmin(i) >= ktop(i)-1) jmin(i)=ktop(i)-2
       END IF
    END DO

    DO ib=1,nLeft
       i=left(ib)
101    CONTINUE
       DO k=jmin(i)-1,1,-1
          dh2(i)    = dh2(i) + (zo_cup  (i,k+1) - zo_cup  (i,k)) &
               * (heso_cup(i,jmin(i)) - heso_cup(i,k))
          IF(dh2(i) >  0.0_r8)THEN
             jmin(i)=jmin(i)-1
             IF(jmin(i) > 3)THEN
                IF (jmin(i)-1 <  kdet(i)  ) kdet(i)=jmin(i)-1
                IF (jmin(i)   >= ktop(i)-1) jmin(i)=ktop(i)-2
                dh2(i)=0.0_r8
                go to 101
             ELSE IF(jmin(i) <= 3 .AND. ierr(i) /= 9 )THEN
                ierr(i)=9
             END IF
          END IF
       END DO
    END DO

    DO i=istart,iend
       IF(ierr(i) == 0)THEN
          IF(jmin(i) <= 3 .AND. ierr(i) /= 4 .AND. dh2(i) <= 0.0_r8)THEN
             ierr(i)=4
          END IF
       END IF
    END DO
    !
    ! Must have at least depth_min m between cloud convective base
    ! and cloud top
    !
    DO i=istart,iend
       IF(ierr(i) == 0)THEN
          IF(-zo_cup(i,kbcon(i))+zo_cup(i,ktop(i)) <  depth_min)THEN
             ierr(i)=6
          END IF
       END IF
    END DO
    !
    ! normalized updraft mass flux profile
    !
    CALL cup_up_nms( &
         zu        , & ! zu         (out)
         z_cup     , & ! z_cup      (in)
         mentr_rate, & ! mentr_rate (in)
         cd        , & ! cd         (in)
         kbcon     , & ! kbcon      (in)
         ktop      , & ! ktop       (in)
         nCols     , & ! nCols      (in)
         kMax      , & ! kMax       (in)
         istart    , & ! istart     (in)
         iend      , & ! iend       (in)
         ierr      , & ! ierr       (in)
         k22       , & ! ) ! k22        (in)
         entr2D)

    CALL cup_up_nms( &
         zuo       , & ! zuo        (out)
         zo_cup    , & ! zo_cup     (in)
         mentr_rate, & ! mentr_rate (in)
         cd        , & ! cd         (in)
         kbcon     , & ! kbcon      (in)
         ktop      , & ! ktop       (in)
         nCols     , & ! nCols      (in)
         kMax      , & ! kMax       (in)
         istart    , & ! istart     (in)
         iend      , & ! iend       (in)
         ierr      , & ! ierr       (in)
         k22       , & !  ) ! k22        (in)
         entr2D)
    !
    ! normalized downdraft mass flux profile,also work on bottom
    ! detrainment in this routin
    !
    CALL cup_dd_nms(  &
         zd         , & ! zd          (out)
         z_cup      , & ! z_cup       (in)
         cdd        , & ! cdd         (out)
         mentrd_rate, & ! mentrd_rate (in)
         jmin       , & ! jmin        (in)
         ierr       , & ! ierr        (in)
         nCols      , & ! nCols       (in)
         kMax       , & ! kMax        (in)
         istart     , & ! istart      (in)
         iend       , & ! iend        (in)
         0          , & ! itest       (in)
         kdet       , & ! kdet        (in)
         z1         , & !  ) ! z1          (in)
         detr2D)

    CALL cup_dd_nms(  &
         zdo        , & ! zdo         (out)
         zo_cup     , & ! zo_cup      (in)
         cdd        , & ! cdd         (out)
         mentrd_rate, & ! mentrd_rate (in)
         jmin       , & ! jmin        (in)
         ierr       , & ! ierr        (in)
         nCols      , & ! nCols       (in)
         kMax       , & ! kMax        (in)
         istart     , & ! istart      (in)
         iend       , & ! iend        (in)
         1          , & ! itest       (in)
         kdet       , & ! kdet        (in)
         z1         , & !  ) ! z1          (in)
         detr2D)
    !
    !  downdraft moist static energy
    !
    CALL cup_dd_he (  &
         hes_cup    , &! hes_cup     (in)
         hcd        , &! hcd         (out)
         z_cup      , &! z_cup       (in)
         cdd        , &! cdd         (in)
         mentrd_rate, &! mentrd_rate (in)
         jmin       , &! jmin        (in)
         ierr       , &! ierr        (in)
         nCols      , &! nCols       (in)
         kMax       , &! kMax        (in)
         istart     , &! istart      (in)
         iend       , &! iend        (in)
         he         , &! he          (in)
         dbyd       , & !  )! dbyd        (out)
         detr2D)

    CALL cup_dd_he (  &
         heso_cup   , &! heso_cup    (in)
         hcdo       , &! hcdo        (out)
         zo_cup     , &! zo_cup      (in)
         cdd        , &! cdd         (in)
         mentrd_rate, &! mentrd_rate (in)
         jmin       , &! jmin        (in)
         ierr       , &! ierr        (in)
         nCols      , &! nCols       (in)
         kMax       , &! kMax        (in)
         istart     , &! istart      (in)
         iend       , &! iend        (in)
         heo        , &! heo         (in)
         dbydo       , &! )! dbydo       (out)
         detr2D)

    !
    !  calculate moisture properties of downdraft
    !
    !
    !snf out  qcd = cloud q (including liquid water) after entrainment
    CALL cup_dd_moisture( &
         zd         , & ! zd          (in)
         hcd        , & ! hcd         (in)
         hes_cup    , & ! hes_cup     (in)
         qcd        , & ! qcd         (out)
         qes_cup    , & ! qes_cup     (in)
         pwd        , & ! pwd         (out)
         q_cup      , & ! q_cup       (in)
         z_cup      , & ! z_cup       (in)
         cdd        , & ! cdd         (in)
         mentrd_rate, & ! mentrd_rate (in)
         jmin       , & ! jmin        (in)
         ierr       , & ! ierr        (inout)
         gamma_cup  , & ! gamma_cup   (in)
         pwev       , & ! pwev        (out)
         nCols      , & ! nCols       (in)
         kMax       , & ! kMax        (in)
         istart     , & ! istart      (in)
         iend       , & ! iend        (in)
         bu         , & ! bu          (out)
         qrcd       , & ! qrcd        (out)
         q          , & ! q           (in)
         2          , & !  ) ! iloop       (in)
         detr2D)
    CALL cup_dd_moisture( &
         zdo        , & ! zdo         (in)
         hcdo       , & ! hcdo        (in)
         heso_cup   , & ! heso_cup    (in)
         qcdo       , & ! qcdo        (out)
         qeso_cup   , & ! qeso_cup    (in)
         pwdo       , & ! pwdo        (out)
         qo_cup     , & ! qo_cup      (in)
         zo_cup     , & ! zo_cup      (in)
         cdd        , & ! cdd         (in)
         mentrd_rate, & ! mentrd_rate (in)
         jmin       , & ! jmin        (in)
         ierr       , & ! ierr        (inout)
         gammao_cup , & ! gammao_cup  (in)
         pwevo      , & ! pwevo       (out)
         nCols      , & ! nCols       (in)
         kMax       , & ! kMax        (in)
         istart     , & ! istart      (in)
         iend       , & ! iend        (in)
         bu         , & ! bu          (out)
         qrcdo      , & ! qrcdo       (out)
         qo         , & ! qo          (in)
         1          , & !  ) ! iloop       (in)
         detr2D)
    !
    ! calculate moisture properties of updraft
    !
    !snf 
    !OUT
  ! qc = cloud q (including liquid water) after entrainment
  ! qrc = liquid water content in cloud after rainout
  ! pw = condensate that will fall out at that level

   CALL cup_up_moisture( &
         ierr       , & ! ierr       (in)
         z_cup      , & ! z_cup      (in)
         qc         , & ! qc         (out)
         qrc        , & ! qrc        (out)
         pw         , & ! pw         (out)
         pwav       , & ! pwav       (out)
         kbcon      , & ! kbcon      (in)
         ktop       , & ! ktop       (in)
         nCols      , & ! nCols      (in)
         kMax       , & ! kMax       (in)
         istart     , & ! istart     (in)
         iend       , & ! iend       (in)
         cd         , & ! cd         (in)
         dby        , & ! dby        (inout)
         mentr_rate , & ! mentr_rate (in)
         q          , & ! q          (in)
         gamma_cup  , & ! gamma_cup  (in)
         zu         , & ! zu         (in)
         qes_cup    , & ! qes_cup    (in)
         k22        , & ! k22        (in)
         q_cup      , & !  ) ! q_cup      (in)
         entr2D)
   !snf new
      DO k=1,kMax
         DO i=istart,iend
            cupclw(i,k)=qrc(i,k)
         END DO
      END DO

    CALL cup_up_moisture( &
         ierr       , & ! ierr       (in)
         zo_cup     , & ! zo_cup     (in)
         qco        , & ! qco        (out)
         qrco       , & ! qrco       (out)
         pwo        , & ! pwo        (out)
         pwavo      , & ! pwavo      (out)
         kbcon      , & ! kbcon      (in)
         ktop       , & ! ktop       (in)
         nCols      , & ! nCols      (in)
         kMax       , & ! kMax       (in)
         istart     , & ! istart     (in)
         iend       , & ! iend       (in)
         cd         , & ! cd         (in)
         dbyo       , & ! dbyo       (inout)
         mentr_rate , & ! mentr_rate (in)
         q          , & ! q          (in)
         gammao_cup , & ! gammao_cup (in)
         zuo        , & ! zuo        (in)
         qeso_cup   , & ! qeso_cup   (in)
         k22        , & ! k22        (in)
         qo_cup     , & !  ) ! qo_cup     (in)
         entr2D)

    DO k=1,kMax
       DO i=1,nCols
          IF(ierr(i) == 0 .AND. k<= ktop (i).and.k >= kbcon(i))THEN
             IF(qrco(i,k) < 0.0_r8) THEN
                qliq(i,k) = 0.0_r8
             ELSE IF(qrco(i,k) >= 1.0_r8) THEN
                qliq(i,k) = 0.0_r8
             ELSE
                qliq(i,k) = qrco(i,k)
             END IF
          END IF
       END DO
    END DO

    !
    ! calculate workfunctions for updrafts
    !
    CALL cup_up_aa0(  &
         aa0        , & ! aa0       (inout)
         z          , & ! z         (in)
         zu         , & ! zu        (in)
         dby        , & ! dby       (in)
         gamma_cup  , & ! gamma_cup (in)
         t_cup      , & ! t_cup     (in)
         kbcon      , & ! kbcon     (in)
         ktop       , & ! ktop      (in)
         kMax       , & ! kMax      (in)
         nCols      , & ! nCols     (in)
         istart     , & ! istart    (in)
         iend       , & ! iend      (in)
         ierr         ) ! ierr      (inout)

    CALL cup_up_aa0(  &
         aa1        , & ! aa1       (inout)
         zo         , & ! z0        (in)
         zuo        , & ! zu0       (in)
         dbyo       , & ! dbyo      (in)
         gammao_cup , & ! gammao_cup(in)
         tn_cup     , & ! tn_cup    (in)
         kbcon      , & ! kbcon     (in)
         ktop       , & ! ktop      (in)
         kMax       , & ! kMax      (in)
         nCols      , & ! nCols     (in)
         istart     , & ! istart    (in)
         iend       , & ! iend      (in)
         ierr         ) ! ierr      (inout)

    CALL cup_up_CIN(  &
         hh        , & ! hh       (inout
         cine      , & ! aa1       (inout)
         zo         , & ! z0        (in)
         zuo        , & ! zu0       (in)
         dbyo       , & ! dbyo      (in)
         gammao_cup , & ! gammao_cup(in)
         tn_cup     , & ! tn_cup    (in)
         kbcon      , & ! kbcon     (in)
         ktop       , & ! ktop      (in)
         kMax       , & ! kMax      (in)
         nCols      , & ! nCols     (in)
         istart     , & ! istart    (in)
         iend       , & ! iend      (in)
         ierr       , & ! ierr      (inout)
         kk1        , & ! pbl k     (inout)
         z1         , &
         buoy        )


!------------
    DO i=istart,iend
       IF(ierr(i) == 0)THEN
          IF(aa1(i) == 0.0_r8)THEN
             ierr(i)=17
          END IF
       END IF
    END DO
    !
    ! determine downdraft strength in terms of windshear
    !
    CALL cup_dd_edt( &
         ierr      , &! ierr    (in)
         us        , &! us      (in)
         vs        , &! vs      (in)
         zo        , &! zo      (in)
         ktop      , &! ktop    (in)
         kbcon     , &! kbcon   (in)
         edt       , &! edt     (out)
         po        , &! po      (in)
         pwavo     , &! pwavo   (in)
         pwevo     , &! pwevo   (in)
         nCols     , &! nCols   (in)
         kMax      , &! kMax    (in)
         istart    , &! istart  (in)
         iend      , &! iend    (in)
         edtmax    , &! edtmax  (in)
         edtmin    , &! edtmin  (in)
         maxens2   , &! maxens2 (in)
         edtc      , &! edtc    (out)
         vshear    , &! vshear  (out)
         sdp       , &! sdp     (out)
         vws       , &! vws     (out)
         mask      , &! mask    (in)
         edtmax1   , &! edtmax1 (in)
         maxens22    )! maxens22(in)

!nilo new
!------------------------------------------
    out1_mass_flux_up=0.0_r8
    out2_mass_flux_d=0.0_r8
    out3_detr=0.0_r8
    out4_liq_water=0.0_r8
    out5_qc=0.0_r8
    out6_wf=0.0_r8
    out7_entr=0.0_r8
          
    DO i=istart,iend
       IF(ierr(i) ==0)THEN
          out6_wf(i)=aa1(i)
          DO k=1,kMax
             out1_mass_flux_up(i,k)=zuo(i,k)
             out2_mass_flux_d(i,k)=zdo(i,k)*edt(i)
             out3_detr(i,k)=cd(i,k)
             out4_liq_water(i,k)=qrc(i,k) 
             out5_qc(i,k)=qc(i,k)
             out7_entr(i,k)=entr2D(i,k)
          END DO
       END IF
    END DO
!-------------------------------------------

    DO i = 1, ncols
       listim(i)=i
    END DO
    nCols_gz=0
    DO i=1,ncols
       IF ((ierr(i) == 0 ) .OR. (ierr(i) > 995)) THEN
          nCols_gz=nCols_gz+1
          bitx(i)=.TRUE.
       END IF
    END DO

    IF (nCols_gz > 0 ) THEN

        nCols_gz=0
        DO  i=1,ncols
           IF(bitx(i))THEN
             nCols_gz=nCols_gz+1
            litx(nCols_gz) = listim(i)
           END IF  
        END DO

        DO k=1,maxens2
           DO i=1,nCols_gz
              edtc_gz (i,k) = edtc  (litx(i),k)
           END DO
        END DO

        DO k=1,ensdim
              DO i=1,nCols_gz
                 pr_ens_gz   (i,k) =pr_ens   (litx(i),k) 
                 outt_ens_gz (i,k) =outt_ens (litx(i),k)
                 xf_ens_gz   (i,k) =xf_ens   (litx(i),k) 
                 massfln_gz  (i,k) =massfln  (litx(i),k) 
              END DO
        END DO


        DO j=1,kMax
              DO i=1,nCols_gz
                     heo_cup_gz  (i,j)=heo_cup    (litx(i),j)
                     zo_cup_gz   (i,j)=zo_cup     (litx(i),j)
                     po_cup_gz   (i,j)=po_cup     (litx(i),j)
                     hcdo_gz     (i,j)=hcdo       (litx(i),j)
                     zdo_gz      (i,j)=zdo        (litx(i),j)
                     cdd_gz      (i,j)=cdd        (litx(i),j)
                     heo_gz      (i,j)=heo        (litx(i),j)
                     qo_cup_gz   (i,j)=qo_cup     (litx(i),j)
                     qrcdo_gz    (i,j)=qrcdo      (litx(i),j)
                     qo_gz       (i,j)=qo         (litx(i),j)
                     zuo_gz      (i,j)=zuo        (litx(i),j)
                     cd_gz       (i,j)=cd         (litx(i),j)
                     hco_gz      (i,j)=hco        (litx(i),j)
                     qco_gz      (i,j)= qco       (litx(i),j)
                     qrco_gz     (i,j)= qrco      (litx(i),j)
                     tn_gz       (i,j)= tn        (litx(i),j)
                     po_gz       (i,j)= po        (litx(i),j)
                     gamma_cup_gz(i,j)= gamma_cup (litx(i),j)
                     pwo_gz      (i,j)= pwo       (litx(i),j)
                     pwdo_gz     (i,j)= pwdo      (litx(i),j)
                     he_cup_gz   (i,j)= he_cup    (litx(i),j)
                     heso_cup_gz (i,j)= heso_cup  (litx(i),j)
                     omeg_gz     (i,j)= omeg      (litx(i),j)
                     p_cup_gz    (i,j)= p_cup     (litx(i),j)
                     entr2D_gz   (i,j)= entr2D    (litx(i),j)
                     detr2D_gz   (i,j)= detr2D    (litx(i),j)
                     den1_gz     (i,j)=den1       (litx(i),j)
              END DO
        END DO


       DO i=1,nCols_gz
               ierr_gz    (i)=ierr    (litx(i))
               ktop_gz    (i)=ktop    (litx(i))
               k22_gz     (i)=k22     (litx(i))
               kbcon_gz   (i)=kbcon   (litx(i))
               jmin_gz    (i)=jmin    (litx(i))
               kdet_gz    (i)=kdet    (litx(i))
               z1_gz      (i)=z1      (litx(i))
               psur_gz    (i)=psur    (litx(i))
               kbmax_gz   (i)=kbmax   (litx(i))
               cap_max_gz (i)=cap_max (litx(i))
               aa0_gz     (i)=aa0     (litx(i))
               aa1_gz     (i)=aa1     (litx(i))
               xmb_gz     (i)=xmb     (litx(i))
               mask_gz    (i)=mask    (litx(i))
               mconv_gz   (i)=mconv   (litx(i))
               dcape_gz   (i)=dcape   (litx(i))
               cine_gz    (i)=cine    (litx(i))
               sens_gz    (i)=cine    (litx(i))
               tkeMax_gz  (i)=tkeMax  (litx(i))
               tkeMedia_gz(i)=tkeMedia(litx(i))
       END DO

       CALL Ensemble(&
            istart                        , &     !INTEGER,(IN   )
            nCols_gz                      , &     !INTEGER,(IN   )
            nCols                         , &     !INTEGER,(IN   )
            kMax                          , &     !INTEGER,(IN   )
            maxens                        , &     !INTEGER,(IN   )
            maxens2                       , &     !INTEGER,(IN   )
            maxens22                      , &     !INTEGER,(IN   )
            maxens3                       , &     !INTEGER,(IN   )
            ensdim                        , &     !INTEGER,(IN   )
            mbdt                          , &     !REAL   ,(IN   )
            dtime                         , &     !REAL   ,(IN   )
            edtc_gz       , &          !REAL   ,(IN   )(nCols,     maxens2)        edtc            (1:nCols   ,1:maxens2        ), &         !REAL   ,(IN        )(nCols,     maxens2)
            ierr_gz       , &          !INTEGER,(INOUT)(nCols             )        ierr            (1:nCols                        ), &         !INTEGER,(INOUT)(nCols             )
            dellat_ens_gz , &          !REAL   ,(OUT  )(nCols,kMax,maxens2)     dellat_ens         (1:nCols   ,1:kMax,1:maxens2), &     !REAL   ,(OUT  )(nCols,kMax,maxens2)
            dellaq_ens_gz , &          !REAL   ,(OUT  )(nCols,kMax,maxens2)             dellaq_ens         (1:nCols   ,1:kMax,1:maxens2), &     !REAL   ,(OUT  )(nCols,kMax,maxens2)
            dellaqc_ens_gz, &          !REAL   ,(OUT  )(nCols,kMax,maxens2)             dellaqc_ens   (1:nCols   ,1:kMax,1:maxens2), &     !REAL   ,(OUT  )(nCols,kMax,maxens2)
            pwo_ens_gz    , &          !REAL   ,(OUT  )(nCols,kMax,maxens2)             pwo_ens         (1:nCols   ,1:kMax,1:maxens2), &     !REAL   ,(OUT  )(nCols,kMax,maxens2)
            heo_cup_gz    , &          !REAL   ,(IN   )(nCols,kMax             )             heo_cup         (1:nCols   ,1:kMax  ), &     !REAL   ,(IN   )(nCols,kMax         )
            zo_cup_gz     , &          !REAL   ,(IN   )(nCols,kMax             )             zo_cup         (1:nCols   ,1:kMax  ), &     !REAL   ,(IN   )(nCols,kMax         )
            po_cup_gz     , &          !REAL   ,(INOUT)(nCols,kMax             )             po_cup         (1:nCols   ,1:kMax  ), &     !REAL   ,(INOUT)(nCols,kMax         )
            hcdo_gz       , &          !REAL   ,(IN   )(nCols,kMax             )             hcdo          (1:nCols   ,1:kMax  ), &     !REAL   ,(IN   )(nCols,kMax         )
            zdo_gz        , &          !REAL   ,(IN   )(nCols,kMax             )             zdo           (1:nCols   ,1:kMax  ), &     !REAL   ,(IN   )(nCols,kMax         )
            cdd_gz        , &          !REAL   ,(INOUT)(nCols,kMax             )             cdd           (1:nCols   ,1:kMax  ), &     !REAL   ,(INOUT)(nCols,kMax         )
            heo_gz        , &          !REAL   ,(IN   )(nCols,kMax             )             heo           (1:nCols   ,1:kMax  ), &     !REAL   ,(IN   )(nCols,kMax         )
            qo_cup_gz     , &          !REAL   ,(IN   )(nCols,kMax             )             qo_cup         (1:nCols   ,1:kMax  ), &     !REAL   ,(IN   )(nCols,kMax         )
            qrcdo_gz      , &          !REAL   ,(IN   )(nCols,kMax             )             qrcdo         (1:nCols   ,1:kMax  ), &     !REAL   ,(IN   )(nCols,kMax         )
            qo_gz         , &          !REAL   ,(IN   )(nCols,kMax             )             qo                 (1:nCols   ,1:kMax  ), &     !REAL   ,(IN   )(nCols,kMax         )
            zuo_gz        , &          !REAL   ,(IN   )(nCols,kMax             )             zuo           (1:nCols   ,1:kMax  ), &     !REAL   ,(IN   )(nCols,kMax         )
            cd_gz         , &          !REAL   ,(IN   )(nCols,kMax             )             cd                 (1:nCols   ,1:kMax  ), &     !REAL   ,(IN   )(nCols,kMax         )
            hco_gz        , &          !REAL   ,(IN   )(nCols,kMax             )             hco           (1:nCols   ,1:kMax  ), &     !REAL   ,(IN   )(nCols,kMax         )
            ktop_gz       , &          !INTEGER,(IN   )(nCols             )        ktop            (1:nCols                ), &         !INTEGER,(IN        )(nCols             )
            k22_gz        , &     !INTEGER,(IN   )(nCols              )       k22            (1:nCols                ), &         !INTEGER,(IN        )(nCols             )
            kbcon_gz      , &          !INTEGER,(IN   )(nCols             )        kbcon            (1:nCols                ), &         !INTEGER,(IN        )(nCols             )
            jmin_gz       , &          !INTEGER,(IN   )(nCols             )        jmin            (1:nCols                ), &         !INTEGER,(IN        )(nCols             )
            kdet_gz       , &          !INTEGER,(IN   )(nCols             )        kdet            (1:nCols                ), &         !INTEGER,(IN        )(nCols             )
            qco_gz        , &          !REAL   ,(IN   )(nCols,kMax             )             qco           (1:nCols   ,1:kMax  ), &     !REAL   ,(IN   )(nCols,kMax         )
            qrco_gz       , &          !REAL   ,(IN   )(nCols,kMax             )             qrco          (1:nCols   ,1:kMax  ), &     !REAL   ,(IN   )(nCols,kMax         )
            tn_gz         , &          !REAL   ,(IN   )(1:nCols_gz,kMax   )             tn                 (1:nCols   ,1:kMax  ), &     !REAL   ,(IN   )(nCols,kMax         )
            po_gz         , &          !REAL   ,(IN   )(nCols,kMax             )             po                 (1:nCols   ,1:kMax  ), &     !REAL   ,(IN   )(nCols,kMax         )
            z1_gz         , &          !REAL   ,(IN   )(nCols             )        z1            (1:nCols                ), &         !REAL   ,(IN        )(nCols             )
            psur_gz       , &          !REAL   ,(IN   )(nCols             )        psur            (1:nCols                ), &         !REAL   ,(IN        )(nCols             )
            gamma_cup_gz  , &            !REAL   ,(INOUT)(nCols,kMax        )             gamma_cup     (1:nCols   ,1:kMax  ), &        !REAL        ,(INOUT)(nCols,kMax           )
            pr_ens_gz     , &            !REAL   ,(INOUT)(nCols,ensdim      )             pr_ens           (1:nCols   ,1:ensdim), &        !REAL        ,(INOUT)(nCols,ensdim           )
            pwo_gz        , &            !REAL   ,(IN   )(nCols,kMax        )             pwo           (1:nCols   ,1:kMax  ), &        !REAL        ,(IN   )(nCols,kMax           )
            pwdo_gz       , &            !REAL   ,(IN   )(nCols,kMax        )             pwdo           (1:nCols   ,1:kMax  ), &        !REAL        ,(IN   )(nCols,kMax           )
            outt_ens_gz   , &            !REAL   ,(OUT  )(nCols,ensdim      )             outt_ens           (1:nCols   ,1:ensdim), &        !REAL        ,(OUT  )(nCols,ensdim           )        
            he_cup_gz     , &            !REAL   ,(IN   )(nCols,kMax        )             he_cup           (1:nCols   ,1:kMax  ), &        !REAL        ,(IN   )(nCols,kMax           )
            kbmax_gz      , &          !INTEGER,(IN   )(nCols             )        kbmax            (1:nCols                ), &         !INTEGER,(IN        )(nCols             )
            heso_cup_gz   , &          !REAL   ,(IN   )(nCols,kMax             )        heso_cup      (1:nCols   ,1:kMax  ), &         !REAL   ,(IN        )(nCols,kMax            )
            cap_max_gz    , &          !REAL   ,(IN   )(nCols             )        cap_max            (1:nCols                ), &         !REAL   ,(IN        )(nCols             )
            aa0_gz        , &          !REAL   ,(INOUT)(nCols             )        aa0            (1:nCols                ), &         !REAL   ,(INOUT)(nCols             )
            aa1_gz        , &          !REAL   ,(IN   )(nCols             )        aa1            (1:nCols                ), &         !REAL   ,(IN        )(nCols             )
            xmb_gz        , &          !REAL   ,(OUT  )(nCols             )        xmb            (1:nCols                ), &         !REAL   ,(OUT  )(nCols             )
            xf_ens_gz     , &          !REAL   ,(OUT  )(nCols,ensdim      )        xf_ens            (1:nCols   ,1:ensdim), &         !REAL   ,(OUT  )(nCols,ensdim      )
            mask_gz       , &          !INTEGER,(IN   )(nCols             )        mask            (1:nCols                ), &         !INTEGER,(IN        )(nCols             )
            mconv_gz      , &          !REAL   ,(IN   )(nCols             )        mconv            (1:nCols                ), &         !REAL   ,(IN        )(nCols             )
            omeg_gz       , &            !REAL   ,(IN   )(nCols,kMax        )             omeg           (1:nCols   ,1:kMax  ), &        !REAL        ,(IN   )(nCols,kMax           )
            massfln_gz    , &            !REAL   ,(OUT  )(nCols,ensdim      )             massfln           (1:nCols   ,1:ensdim), &        !REAL        ,(OUT  )(nCols,ensdim           )
            p_cup_gz      , & !!!nilo )  !REAL   ,(IN   )(nCols,kMax        )             p_cup           (1:nCols   ,1:kMax)  )        !REAL        ,(IN   )(nCols,kMax           )
            entr2D_gz     , &    !
            detr2D_gz     , &  !!)
            den1_gz       , &
            dcape_gz      , &
            cine_gz       , &
            sens_gz       , &
            tkeMax_gz     , &
            tkeMedia_gz   )


       DO i=1,nCols_gz
          ierr(litx(i)) = ierr_gz(i)
       END DO

       DO i=1,nCols_gz
              aa0  (litx(i)) = aa0_gz (i)
              xmb  (litx(i)) = xmb_gz( i)
       END DO
       
       DO k=1,kMax
          DO i=1,nCols_gz
                 po_cup   (litx(i),k) = po_cup_gz   (i,k)
                 cdd      (litx(i),k) = cdd_gz      (i,k) 
                 gamma_cup(litx(i),k) = gamma_cup_gz(i,k)
                 qrco     (litx(i),k) = qrco_gz     (i,k) 
          END DO
       END DO
       
       DO j=1,ensdim
          DO i=1,nCols_gz
                 pr_ens   (litx(i),j) = pr_ens_gz   (i,j) 
                 outt_ens (litx(i),j) = outt_ens_gz (i,j) 
                 xf_ens   (litx(i),j) = xf_ens_gz   (i,j) 
                 massfln  (litx(i),j) = massfln_gz  (i,j) 
          END DO
       END DO
       DO j=1,maxens2
          DO k=1,kMax
             DO i=1,nCols_gz
                   dellat_ens  (litx(i),k,j)=dellat_ens_gz (i,k,j)
                   dellaq_ens  (litx(i),k,j)=dellaq_ens_gz (i,k,j)
                   dellaqc_ens (litx(i),k,j)=dellaqc_ens_gz(i,k,j)
                   pwo_ens     (litx(i),k,j)=pwo_ens_gz    (i,k,j)
             END DO
          END DO
       END DO

    DO k=1,kMax
       DO i=1,nCols
          IF( ierr(i) == 0 .and.k<= ktop (i).and.k >= kbcon(i))THEN
             IF(qrco(i,k) < 0.0_r8) THEN
                qliq(i,k) = 0.0_r8
             ELSE IF(qrco(i,k) >= 1.0_r8) THEN
                qliq(i,k) = 0.0_r8
             ELSE
                qliq(i,k) = qrco(i,k)
             END IF
          END IF
       END DO
    END DO

    !
    !--- FEEDBACK
    !
    IF(i_use_stat_prop == 0 ) THEN

    CALL cup_output_ens( &
         xf_ens     , & ! xf_ens      (in)
         ierr       , & ! ierr        (inout)
         dellat_ens , & ! dellat_ens  (in)
         dellaq_ens , & ! dellaq_ens  (in)
         dellaqc_ens, & ! dellaqc_ens (in)
         outt       , & ! outt        (inout) hmjb
         outq       , & ! outq        (inout) hmjb
         outqc      , & ! outqc       (out)
         pre        , & ! pre         (out)
         pwo_ens    , & ! pwo_ens     (in)
         xmb        , & ! xmb         (out)
         ktop       , & ! ktop        (in)
         nCols      , & ! nCols       (in)
         kMax       , & ! kMax        (in)
         istart     , & ! istart      (in)
         iend       , & ! iend        (in)
         maxens2    , & ! maxens2     (in)
         maxens     , & ! maxens      (in)
         iens       , & ! iens        (in)
         pr_ens     , & ! pr_ens      (inout)
         outt_ens   , & ! outt_ens    (inout)
         maxens3    , & ! maxens3     (in)
         ensdim     , & ! ensdim      (in)
         massfln    , & ! massfln     (inout)
         xfac1      , & ! xfac1       (out)
         xfac_for_dn, & ! xfac_for_dn (out) 
         maxens22    ) ! maxens22    (in)
    ELSE
    CALL cup_output_ens_catt(              &
        xf_ens         , &! intent(in   )
        ierr           , &! intent(inout)
        dellat_ens     , &! intent(in        )
        dellaq_ens     , &! intent(in        )
        dellaqc_ens    , &! intent(in        )
        outt           , &! intent(out  )
        outq           , &! intent(out  )
        outqc          , &! intent(out  )
        pre            , &! intent(out  )
        pwo_ens        , &! intent(in   )
        xmb            , &! intent(out  )
        ktop           , &! intent(in   )
        nCols          , &! intent(in   )
        kMax           , &! intent(in   ) 
        istart         , & ! istart      (in)
        iend           , & ! iend        (in)
        maxens2        , &! intent(in   )
        maxens         , &! intent(in   )
        iens           , &! intent(in   )
        pr_ens         , &! intent(inout)
        outt_ens       , &! intent(inout)
        maxens3        , &! intent(in   )
        ensdim         , &! intent(in   )
        massfln        , &! intent(inout)
        xfac1          , &! intent(out  )
        xfac_for_dn       )! intent(out  )

    END IF               
    DO i=istart,iend
       pre(i)=MAX(pre(i),0.0_r8)
       !snf
    END DO
    ELSE
        outt =0.0_r8
        outq  =0.0_r8
        outqc =0.0_r8
        pre=0.0_r8
        xmb=0.0_r8
        pr_ens=0.0_r8
        outt_ens=0.0_r8
        massfln=0.0_r8
        xfac1=0.0_r8
        xfac_for_dn=0.0_r8
        qliq=0.0_r8
    END IF

    RETURN
  END SUBROUTINE grellens2
  !
  !END CUP
  !

  SUBROUTINE Ensemble(&
       istart     , &
       iend       , &
       nCols      , &
       kMax       , &
       maxens     , & 
       maxens2    , &
       maxens22   , &
       maxens3    , &
       ensdim     , &
       mbdt       , &
       dtime      , & 
       edtc       , &
       ierr       , &
       dellat_ens , &
       dellaq_ens , &
       dellaqc_ens, &
       pwo_ens    , &
       heo_cup    , &
       zo_cup     , &
       po_cup     , & 
       hcdo       , &
       zdo        , &
       cdd        , &
       heo        , &
       qo_cup     , &
       qrcdo      , &
       qo         , &
       zuo        , &
       cd         , &
       hco        , &
       ktop       , &
       k22        , &
       kbcon      , &
       jmin       , &
       kdet       , &
       qco        , &
       qrco       , &
       tn         , &
       po         , & 
       z1         , & 
       psur       , & 
       gamma_cup  , & 
       pr_ens     , &
       pwo        , &
       pwdo       , &
       outt_ens   , &
       he_cup     , &
       kbmax      , &
       heso_cup   , &    
       cap_max    , &
       aa0        , &
       aa1        , &
       xmb        , &
       xf_ens     , &
       mask       , &
       mconv      , &
       omeg       , &
       massfln    , &
       p_cup      , &  !!)
       entr2D     , &
       detr2D     , & !! )
       den1       , &
       dcape       , &
       cine       , &
       sens       , &
       tkeMax     , &
       tkeMedia   )


    INTEGER :: iedt
    INTEGER, INTENT(IN   ) :: istart
    INTEGER, INTENT(IN   ) :: iend
    INTEGER, INTENT(IN   ) :: nCols
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: maxens
    INTEGER, INTENT(IN   ) :: maxens2,maxens22
    INTEGER, INTENT(IN   ) :: maxens3
    INTEGER, INTENT(IN   ) :: ensdim 
    REAL(KIND=r8)   , INTENT(IN   ) :: mbdt
    REAL(KIND=r8)   , INTENT(IN   ) :: dtime
    REAL(KIND=r8)   , INTENT(IN   ) :: edtc       (nCols,maxens2)
    INTEGER, INTENT(INOUT) :: ierr       (nCols)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dellat_ens (nCols,kMax, maxens2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dellaq_ens (nCols,kMax, maxens2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dellaqc_ens(nCols,kMax, maxens2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: pwo_ens    (nCols,kMax, maxens2)
    REAL(KIND=r8)   , INTENT(IN   ) :: heo_cup    (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: zo_cup     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: po_cup     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: hcdo       (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: zdo        (nCols,kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: cdd        (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: heo        (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: qo_cup     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: qrcdo      (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: qo         (nCols,kMax)! water vapor mixing ratio (kg/kg) at time t+1
    REAL(KIND=r8)   , INTENT(IN   ) :: zuo        (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: cd         (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: hco        (nCols,kMax)
    INTEGER, INTENT(IN   ) :: ktop       (nCols)
    INTEGER, INTENT(IN   ) :: k22        (nCols)
    INTEGER, INTENT(IN   ) :: kbcon      (nCols)! level of convective cloud base       
    INTEGER, INTENT(IN   ) :: jmin       (nCols)
    INTEGER, INTENT(IN   ) :: kdet       (nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: qco        (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: qrco       (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: tn         (nCols,kMax)! temperature (K) at time t+1
    REAL(KIND=r8)   , INTENT(IN   ) :: po         (nCols,kMax)! pressao de superficie no tempo t mb
    REAL(KIND=r8)   , INTENT(IN   ) :: z1         (nCols)! topography (m)
    REAL(KIND=r8)   , INTENT(IN   ) :: psur       (nCols)! pressao de superficie no tempo t mb
    REAL(KIND=r8)   , INTENT(INOUT) :: gamma_cup  (nCols,kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: pr_ens     (nCols,ensdim)
    REAL(KIND=r8)   , INTENT(IN   ) :: pwo        (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: pwdo       (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: outt_ens   (nCols,ensdim)    
    REAL(KIND=r8)   , INTENT(IN   ) :: he_cup     (nCols,kMax)
    INTEGER, INTENT(IN   ) :: kbmax      (nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: heso_cup   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: cap_max    (nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: aa0        (nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: aa1        (nCols)
    REAL(KIND=r8)   , INTENT(OUT  ) :: xmb        (nCols)
    REAL(KIND=r8)   , INTENT(OUT  ) :: xf_ens     (nCols,ensdim)
    INTEGER, INTENT(IN   ) :: mask       (nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: mconv      (nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: omeg       (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: massfln    (nCols,ensdim)
    REAL(KIND=r8)   , INTENT(IN   ) :: p_cup      (nCols,kMax)
    
    !
    ! LOCAL VARIABLE
    !
    REAL(KIND=r8)    :: dellat    (nCols,kMax)
    REAL(KIND=r8)    :: dellaq    (nCols,kMax)
    REAL(KIND=r8)    :: dellah    (nCols,kMax)
    REAL(KIND=r8)    :: dellaqc   (nCols,kMax)
    REAL(KIND=r8)    :: xhe       (nCols,kMax)
    REAL(KIND=r8)    :: edt       (nCols)
    REAL(KIND=r8)    :: bu        (nCols)    
    INTEGER          :: ierr2     (nCols)
    INTEGER          :: ierr3     (nCols)    
    REAL(KIND=r8)    :: dbyd      (nCols,kMax)
    REAL(KIND=r8)    :: xq        (nCols,kMax)
    REAL(KIND=r8)    :: xt        (nCols,kMax)
    REAL(KIND=r8)    :: xqes      (nCols,kMax)
    REAL(KIND=r8)    :: edto      (nCols)
    REAL(KIND=r8)    :: xhes      (nCols,kMax)
    REAL(KIND=r8)    :: xff_ens3  (nCols,maxens3)
    REAL(KIND=r8)    :: xk        (nCols,maxens) 
    REAL(KIND=r8)    :: xaa0_ens  (nCols,maxens)
    INTEGER          :: k22x      (nCols)    
    INTEGER          :: kbconx    (nCols)
    INTEGER          :: nallp     
    REAL(KIND=r8)    :: xaa0      (nCols)
    REAL(KIND=r8)    :: denx      (nCols,kMax)
    REAL(KIND=r8)    :: xt_cup    (nCols,kMax)
    REAL(KIND=r8)    :: xdby      (nCols,kMax)
    REAL(KIND=r8)    :: xzu       (nCols,kMax)
    REAL(KIND=r8)    :: xz        (nCols,kMax)
    REAL(KIND=r8)    :: xq_cup    (nCols,kMax)
    REAL(KIND=r8)    :: xqes_cup  (nCols,kMax)
    REAL(KIND=r8)    :: xpwav     (nCols)
    REAL(KIND=r8)    :: xpw       (nCols,kMax)
    REAL(KIND=r8)    :: xqrc      (nCols,kMax)
    REAL(KIND=r8)    :: xqc       (nCols,kMax)
    REAL(KIND=r8)    :: xz_cup    (nCols,kMax)
    REAL(KIND=r8)    :: xqrcd     (nCols,kMax)
    REAL(KIND=r8)    :: xpwev     (nCols)
    REAL(KIND=r8)    :: xpwd      (nCols,kMax)
    REAL(KIND=r8)    :: xqcd      (nCols,kMax)
    REAL(KIND=r8)    :: xhes_cup  (nCols,kMax)
    REAL(KIND=r8)    :: xhcd      (nCols,kMax)
    REAL(KIND=r8)    :: xzd       (nCols,kMax)
    REAL(KIND=r8)    :: xhc       (nCols,kMax)
    REAL(KIND=r8)    :: xhe_cup   (nCols,kMax)
    REAL(KIND=r8)    :: xhkb      (nCols)
    REAL(KIND=r8)    :: scr1      (nCols,kMax)
    REAL(KIND=r8)    :: dz
    REAL(KIND=r8)    :: massfld    
    INTEGER :: i
    INTEGER :: k
    INTEGER :: nens
    INTEGER :: nens3
    REAL(KIND=r8)   , INTENT(IN   ) :: entr2D   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: detr2D   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: den1     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: dcape (nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: cine (nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: sens (nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: tkeMax(nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: tkeMedia(nCols)

    dellat    =0.0_r8;dellaq    =0.0_r8
    dellah    =0.0_r8;dellaqc   =0.0_r8
    xhe       =0.0_r8;edt       =0.0_r8
    bu        =0.0_r8;ierr2     =0
    ierr3     =0;dbyd           =0.0_r8
    xq        =0.0_r8;xt        =0.0_r8
    xqes      =0.0_r8;edto      =0.0_r8
    xhes      =0.0_r8;xff_ens3  =0.0_r8
    xk        =0.0_r8;xaa0_ens  =0.0_r8
    k22x      =0;kbconx         =0
    nallp     =0;xaa0           =0.0_r8
    xt_cup    =0.0_r8;xdby      =0.0_r8
    xzu       =0.0_r8;xz        =0.0_r8
    xq_cup    =0.0_r8;xqes_cup  =0.0_r8
    xpwav     =0.0_r8;xpw       =0.0_r8
    xqrc      =0.0_r8;xqc       =0.0_r8
    xz_cup    =0.0_r8;xqrcd     =0.0_r8
    xpwev     =0.0_r8;xpwd      =0.0_r8
    xqcd      =0.0_r8;xhes_cup  =0.0_r8
    xhcd      =0.0_r8;xzd       =0.0_r8
    xhc       =0.0_r8;xhe_cup   =0.0_r8
    xhkb      =0.0_r8;scr1      =0.0_r8
    dz        =0.0_r8;massfld   =0.0_r8
    !
    ! LOOP FOR ENSEMBLE MAXENS2
    !
    DO 250 iedt=1,maxens22

       DO i=istart,iend
          IF(ierr(i) == 0)THEN
             edt (i)=edtc(i,iedt)
             edto(i)=edtc(i,iedt)
          END IF
       END DO

       DO k=1,kMax
          DO i=istart,iend
             dellat_ens (i,k,iedt)=0.0_r8
             dellaq_ens (i,k,iedt)=0.0_r8
             dellaqc_ens(i,k,iedt)=0.0_r8
             pwo_ens    (i,k,iedt)=0.0_r8
          END DO
       END DO
       !
       !--- downdraft workfunctions
       !
       !
       !--- change per unit mass that a model cloud would modify the environment
       !
       !--- 1.0_r8 in bottom layer
       !
       CALL cup_dellabot( &
            heo_cup    , &  ! heo_cup     (in)
            ierr       , &  ! ierr        (in)
            zo_cup     , &  ! zo_cup      (in)
            po_cup     , &  ! po_cup      (in)
            hcdo       , &  ! hcdo        (in)
            edto       , &  ! edto        (in)
            zdo        , &  ! zdo         (in)
            cdd        , &  ! cdd         (in)
            heo        , &  ! heo         (in)
            nCols      , &  ! nCols       (in)
            kMax       , &  ! kMax        (in)
            istart     , &  ! istart      (in)
            iend       , &  ! iend        (in)
            dellah     , &  ! dellah      (out)
            mentrd_rate, &  !  mentrd_rate (in)
            detr2D)

       CALL cup_dellabot(&
            qo_cup     , &  ! qo_cup      (in)
            ierr       , &  ! ierr        (in)
            zo_cup     , &  ! zo_cup      (in)
            po_cup     , &  ! po_cup      (in)
            qrcdo      , &  ! qrcdo       (in)
            edto       , &  ! edto        (in)
            zdo        , &  ! zdo         (in)
            cdd        , &  ! cdd         (in)
            qo         , &  ! qo          (in)
            nCols      , &  ! nCols       (in)
            kMax       , &  ! kMax        (in)
            istart     , &  ! istart      (in)
            iend       , &  ! iend        (in)
            dellaq     , &  ! dellaq      (out)
            mentrd_rate, &  !  mentrd_rate (in)
            detr2D)
       !
       !--- 2. everywhere else
       !
       CALL cup_dellas(&
            ierr       , &  ! ierr        (in)
            zo_cup     , &  ! zo_cup      (in)
            po_cup     , &  ! po_cup      (in)
            hcdo       , &  ! hcdo        (in)
            edto       , &  ! edto        (in)
            zdo        , &  ! zdo         (in)
            cdd        , &  ! cdd         (in)
            heo        , &  ! heo         (in)
            nCols      , &  ! nCols       (in)
            kMax       , &  ! kMax        (in)
            istart     , &  ! istart      (in)
            iend       , &  ! iend        (in)
            dellah     , &  ! dellah      (out)
            mentrd_rate, &  ! mentrd_rate (in)
            zuo        , &  ! zuo         (in)
            cd         , &  ! cd          (in)
            hco        , &  ! hco         (in)
            ktop       , &  ! ktop        (in)
            k22        , &  ! k22         (in)
            kbcon      , &  ! kbcon       (in)
            mentr_rate , &  ! mentr_rate  (in)
            jmin       , &  ! jmin        (in)
            heo_cup    , &  ! heo_cup     (in)
            kdet       , &  ! kdet        (in)
            k22        , &  !)    ! k22         (in)
            entr2D     , &  !  )  ! mentrd_rate (in)
            detr2D      )

       !
       !-- take out cloud liquid water for detrainment
       !
       DO k=1,kMax
          DO i=istart,iend
             scr1   (i,k)=0.0_r8
             dellaqc(i,k)=0.0_r8
             IF(ierr(i) == 0)THEN
                scr1(i,k)=qco(i,k)-qrco(i,k)
                IF(k == ktop(i)-0)dellaqc(i,k)=                 &
                     0.01_r8*zuo(i,ktop(i))*qrco(i,ktop(i))*        &
                     9.81_r8/(po_cup(i,k  )-po_cup(i,k+1))

                IF(k <  ktop(i)  .AND.k >  kbcon(i))THEN
                   dz=zo_cup(i,k+1)-zo_cup(i,k)
                   dellaqc(i,k)=0.01_r8*9.81_r8*cd(i,k)*dz*zuo(i,k)    &
                        *0.5_r8*(qrco(i,k)+qrco(i,k+1))/            &
                        (po_cup(i,k  )-po_cup(i,k+1))
                END IF
             END IF
          END DO
       END DO
       !
       CALL cup_dellas( &
            ierr       , &  ! ierr        (in)
            zo_cup     , &  ! zo_cup      (in)
            po_cup     , &  ! po_cup      (in)
            qrcdo      , &  ! qrcdo       (in)
            edto       , &  ! edto        (in)
            zdo        , &  ! zdo         (in)
            cdd        , &  ! cdd         (in)
            qo         , &  ! qo          (in)
            nCols      , &  ! nCols       (in)
            kMax       , &  ! kMax        (in)
            istart     , &  ! istart      (in)
            iend       , &  ! iend        (in)
            dellaq     , &  ! dellaq      (out)
            mentrd_rate, &  ! mentrd_rate (in)
            zuo        , &  ! zuo         (in)
            cd         , &  ! cd          (in)
            scr1       , &  ! scr1        (in)
            ktop       , &  ! ktop        (in)
            k22        , &  ! k22         (in)
            kbcon      , &  ! kbcon       (in)
            mentr_rate , &  ! mentr_rate  (in)
            jmin       , &  ! jmin        (in)
            qo_cup     , &  ! qo_cup      (in)
            kdet       , &  ! kdet        (in)
            k22        , &  !  )  ! k22         (in)
            entr2D     , &  !  )  ! mentrd_rate (in)
            detr2D)
       !
       !--- using dellas, calculate changed environmental profiles
       !
       !................second loop..........................start 200
!!!old       do 200 nens=1,maxens
!!!old           mbdt=mbdt_ens(nens)
       !
       DO k=1,maxens
          DO i=istart,iend
             xaa0_ens(i,k)=0.0_r8
          END DO
       END DO
       !
       !-----------------------------
       !
       DO k=1,kMax-1
          DO i=istart,iend
             dellat(i,k)=0.0_r8
             IF(ierr(i) == 0)THEN
                xhe   (i,k)=dellah(i,k)*mbdt+heo(i,k)
                xq    (i,k)=dellaq(i,k)*mbdt+qo (i,k)
                dellat(i,k)=(1.0_r8/1004.0_r8)*(dellah(i,k)-2.5e06_r8*dellaq(i,k))
                xt    (i,k)= dellat(i,k)*mbdt+tn(i,k)
                IF(xq(i,k) <= 0.0_r8)xq(i,k)=1.e-08_r8
             END IF
          END DO
       END DO
       !
       !
       DO i=istart,iend
          IF(ierr(i) == 0)THEN
             xhe(i,kMax)=heo(i,kMax)
             xq (i,kMax)=qo (i,kMax)
             xt (i,kMax)=tn (i,kMax)
             IF(xq(i,kMax) <= 0.0_r8)xq(i,kMax)=1.e-08_r8
          END IF
       END DO
       !
       ! calculate moist static energy, heights, qes
       !
       CALL cup_env(&
            xz        , &  ! xz     (out)
            xqes      , &  ! xqes   (out)
            xhe       , &  ! xhe    (inout)
            xhes      , &  ! xhes   (out)
            xt        , &  ! xt     (in)
            xq        , &  ! xq     (inout)
            po        , &  ! po     (in)
            z1        , &  ! z1     (in)
            nCols     , &  ! nCols  (in)
            kMax      , &  ! kMax   (in)
            istart    , &  ! istart (in)
            iend      , &  ! iend   (in)
            psur      , &  ! psur   (in)
            ierr      , &  ! ierr   (in)
             2         , &  !  2     (in)
            denx        )  ! den    (out)
       !
       ! environmental values on cloud levels
       !
       CALL cup_env_clev( &
            xt        , &  ! xt        (in)
            xqes      , &  ! xqes      (in)
            xq        , &  ! xq        (in)
            xhe       , &  ! xhe       (in)
            xhes      , &  ! xhes      (in)
            xz        , &  ! xz        (in)
            po        , &  ! po        (in)
            xqes_cup  , &  ! xqes_cup  (out)
            xq_cup    , &  ! xq_cup    (out)
            xhe_cup   , &  ! xhe_cup   (out)
            xhes_cup  , &  ! xhes_cup  (out)
            xz_cup    , &  ! xz_cup    (out)
            po_cup    , &  ! po_cup    (out)
            gamma_cup , &  ! gamma_cup (out)
            xt_cup    , &  ! xt_cup    (out)
            psur      , &  ! psur      (in)
            nCols     , &  ! nCols     (in)
            kMax      , &  ! kMax      (in)
            istart    , &  ! istart    (in)
            iend      , &  ! iend      (in)
            ierr      , &  ! ierr      (in)
            z1          )  ! z1        (in)
       !
       !STATIC CONTROL
       !
       ! moist static energy inside cloud
       !
       DO i=istart,iend
          IF(ierr(i) == 0)THEN
             xhkb(i)=xhe(i,k22(i))
          END IF
       END DO

          CALL cup_up_he(&
            k22       , & ! k22        (in)
            xhkb      , & ! xhkb       (out)
            xz_cup    , & ! xz_cup     (in)
            cd        , & ! cd         (in)
            mentr_rate, & ! mentr_rate (in)
            xhe_cup   , & ! xhe_cup    (in)
            xhc       , & ! xhc        (out)
            nCols     , & ! nCols      (in)
            kMax      , & ! kMax       (in)
            kbcon     , & ! kbcon      (in)
            ierr      , & ! ierr       (in)
            istart    , & ! istart     (in)
            iend      , & ! iend       (in)
            xdby      , & ! xdby       (out)
            xhe       , & ! xhe        (in)
            xhes_cup  , & !  ) ! xhes_cup   (in)
            entr2D)
       !
       ! normalized mass flux profile
       !
       CALL cup_up_nms(&
            xzu       , & ! xzu        (out)
            xz_cup    , & ! xz_cup     (in)
            mentr_rate, & ! mentr_rate (in)
            cd        , & ! cd         (in)
            kbcon     , & ! kbcon      (in)
            ktop      , & ! ktop       (in)
            nCols     , & ! nCols      (in)
            kMax      , & ! kMax       (in)
            istart    , & ! istart     (in)
            iend      , & ! iend       (in)
            ierr      , & ! ierr       (in)
            k22       , & ! ) ! k22        (in)
            entr2D)

       CALL cup_dd_nms(xzd        , &! xzd        (out)
            xz_cup     , &! xz_cup     (in)
            cdd        , &! cdd        (out)
            mentrd_rate, &! mentrd_rate(in)
            jmin       , &! jmin       (in)
            ierr       , &! ierr       (in)
            nCols      , &! nCols      (in)
            kMax       , &! kMax       (in)
            istart     , &! istart     (in)
            iend       , &! iend       (in)
            1          , &! 1
            kdet       , &! kdet       (in)
            z1         , & !  )! z1         (in)
            detr2D)

       !
       ! moisture downdraft
       !
       CALL cup_dd_he(xhes_cup   , &! xhes_cup    (in)
            xhcd       , &! xhcd        (out)
            xz_cup     , &! xz_cup      (in)
            cdd        , &! cdd         (in)
            mentrd_rate, &! mentrd_rate (in)
            jmin       , &! jmin        (in)
            ierr       , &! ierr        (in)
            nCols      , &! nCols       (in)
            kMax       , &! kMax        (in)
            istart     , &! istart      (in)
            iend       , &! iend        (in)
            xhe        , &! xhe         (in)
            dbyd       , & !  )! dbyd        (out)
            detr2D)

       CALL cup_dd_moisture(xzd        , &  ! xzd         (in)
            xhcd       , &  ! xhcd        (in)
            xhes_cup   , &  ! xhes_cup    (in)
            xqcd       , &  ! xqcd        (out)
            xqes_cup   , &  ! xqes_cup    (in)
            xpwd       , &  ! xpwd        (out)
            xq_cup     , &  ! xq_cup      (in)
            xz_cup     , &  ! xz_cup      (in)
            cdd        , &  ! cdd         (in)
            mentrd_rate, &  ! mentrd_rate (in)
            jmin       , &  ! jmin        (in)
            ierr       , &  ! ierr        (inout)
            gamma_cup  , &  ! gamma_cup   (in)
            xpwev      , &  ! xpwev       (out)
            nCols      , &  ! nCols       (in)
            kMax       , &  ! kMax        (in)
            istart     , &  ! istart      (in)
            iend       , &  ! iend        (in)
            bu         , &  ! bu          (out)
            xqrcd      , &  ! xqrcd       (out)
            xq         , &  ! xq          (in)
            3          , & !  )  ! 3
            detr2D)

       !
       ! moisture updraft
       !
       
         CALL cup_up_moisture(ierr       , &  ! ierr       (in)
            xz_cup     , &  ! xz_cup     (in)
            xqc        , &  ! xqc        (out)
            xqrc       , &  ! xqrc       (out)
            xpw        , &  ! xpw        (out)
            xpwav      , &  ! xpwav      (out)
            kbcon      , &  ! kbcon      (in)
            ktop       , &  ! ktop       (in)
            nCols      , &  ! nCols      (in)
            kMax       , &  ! kMax       (in)
            istart     , &  ! istart     (in)
            iend       , &  ! iend       (in)
            cd         , &  ! cd         (in)
            xdby       , &  ! xdby       (inout)
            mentr_rate , &  ! mentr_rate (in)
            xq         , &  ! xq         (in)
            gamma_cup  , &  ! gamma_cup  (in)
            xzu        , &  ! xzu        (in)
            xqes_cup   , &  ! xqes_cup   (in)
            k22        , &  ! k22        (in)
            xq_cup     , &  !  )  ! xq_cup     (in)
            detr2D)

       !
       ! workfunctions for updraft
       !
       CALL cup_up_aa0(xaa0       , & ! xaa0      (inout)
            xz         , & ! xz        (in)
            xzu        , & ! xzu       (in)
            xdby       , & ! xdby      (in)
            gamma_cup  , & ! gamma_cup (in)
            xt_cup     , & ! xt_cup    (in)
            kbcon      , & ! kbcon     (in)
            ktop       , & ! ktop      (in)
            kMax       , & ! kMax      (in)
            nCols      , & ! nCols     (in)
            istart     , & ! istart    (in)
            iend       , & ! iend      (in)
            ierr         ) ! ierr      (in)

       !
       ! workfunctions for downdraft
       !---------0--------------
       ! 
       DO 200 nens=1,maxens
          DO i=istart,iend 
             IF(ierr(i) == 0)THEN
                xaa0_ens(i,nens)=xaa0(i)
             END IF
          END DO
          nallp=(iens-1)*maxens3*maxens*maxens22 &
               +(iedt-1)*maxens*maxens3 &
               +(nens-1)*maxens3
          DO nens3=1,maxens3
             DO k=1,kMax
                DO i=istart,iend
                   IF( k <= ktop(i) .AND. ierr(i) == 0 )THEN                
                      pr_ens(i,nallp+nens3)=pr_ens(i,nallp+nens3)+&
                           pwo(i,k)+beta(nens3)*edto(i)*pwdo(i,k)
                   END IF
                END DO
             END DO
             DO i=istart,iend 
                IF(ierr(i) == 0)THEN
                   outt_ens (i,nallp+nens3)=dellat(i,1)
                   IF(pr_ens(i,nallp+nens3) < 0.0_r8)THEN
                      pr_ens(i,nallp+nens3)= 0.0_r8
                   END IF
                END IF
             END DO
          END DO
200    END DO
       !...............end 200
       !
       ! LARGE SCALE FORCING
       !
       CALL cup_maximi(he_cup    , & ! he_cup (in)
            nCols     , & ! nCols  (in)
            kMax      , & ! kMax   (in) 
            3         , & ! 3      (in)
            kbmax     , & ! kbmax  (in)
            k22x      , & ! k22x   (out)
            istart    , & ! istart (in)
            iend      , & ! iend   (in)
            ierr        ) ! ierr   (in)

       DO i=istart,iend
          IF(ierr(i) == 0)THEN
             k22x (i)=k22(i)
          END IF
          ierr2(i)=ierr(i)
          ierr3(i)=ierr(i)
       END DO
       !
       ! --- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
       ! snf  call cup_kbcon for cap_max=cap_max-(2-1)*cap_max_increment
       !
       CALL cup_kbcon(&
            cap_max_increment, &  ! cap_max_increment (in)
            2         , &         ! 2
            k22x      , &         ! k22x              (inout)
            kbconx    , &         ! kbconx            (out)
            heo_cup   , &         ! heo_cup           (in)
            heso_cup  , &         ! heso_cup          (in)
            nCols     , &         ! nCols             (in)
            kMax      , &         ! kMax              (in)
            istart    , &         ! istart            (in)
            iend      , &         ! iend              (in)
            ierr2     , &         ! ierr2             (inout)
            kbmax     , &         ! kbmax             (in)
            po_cup    , &         ! po_cup            (in)
            cap_max     )         ! cap_max           (in)
       !
       ! snf  call cup_kbcon for cap_max=cap_max-(3-1)*cap_max_increment
       !
       CALL cup_kbcon(&
            cap_max_increment, &  ! cap_max_increment (in)
            3                , &  ! 3                 (in)
            k22x             , &  ! k22x              (inout)
            kbconx           , &  ! kbconx            (out)
            heo_cup          , &  ! heo_cup           (in)
            heso_cup         , &  ! heso_cup          (in)
            nCols            , &  ! nCols             (in)
            kMax             , &  ! kMax              (in)
            istart           , &  ! istart            (in)
            iend             , &  ! iend              (in)
            ierr3            , &  ! ierr3             (inout)
            kbmax            , &  ! kbmax             (in)
            po_cup           , &  ! po_cup            (in)
            cap_max            )  ! cap_max           (in)

          CALL cup_forcing_ens( &
               aa0       , & ! aa0      (inout)
               aa1       , & ! aa1      (in)
               xaa0_ens  , & ! xaa0_ens (in)
               mbdt      , & ! mbdt     (in)
               dtime     , & ! dtime    (in)
               xmb       , & ! xmb      (out)
               ierr      , & ! ierr     (inout)
               nCols     , & ! nCols    (in)
               kMax      , & ! kMax     (in)
               istart    , & ! istart   (in)
               iend      , & ! iend     (in)
               xf_ens    , & ! xf_ens   (out)
               'deeps'   , & ! 'deeps'  (in)
               mask      , & ! mask     (in)
               maxens    , & ! maxens   (in)
               iens      , & ! iens     (in)
               iedt      , & ! iedt     (in)
               maxens3   , & ! maxens3  (in)
               mconv     , & ! mconv    (in)
               omeg      , & ! omeg     (in)
               k22       , & ! k22      (in)
               pr_ens    , & ! pr_ens   (in)
               edto      , & ! edto     (in)
               kbcon     , & ! kbcon    (in)
               ensdim    , & ! ensdim   (in)
               massfln   , & ! massfln  (out)
               massfld   , & ! massfld  (inout)
               xff_ens3  , & ! xff_ens3 (out)
               xk        , & ! xk       (out)
               p_cup     , & ! p_cup    (in)
               ktop      , & ! ktop     (in)
               ierr2     , & ! ierr2    (in)
               ierr3     , & ! ierr3    (in)
               grepar1   , & ! grepar1  (in)
               xfmax     , & ! xfmax    (in)
               maxens22  , & !  ) !maxens22  (in)
               den1      , &
               dcape      , &
               cine      , &
               sens      , &
               tkeMax    , &
               tkeMedia  )

       DO k=1,kMax
          DO i=istart,iend
             IF(ierr(i) == 0)THEN
                dellat_ens (i,k,iedt)=dellat (i,k)
                dellaq_ens (i,k,iedt)=dellaq (i,k)
                dellaqc_ens(i,k,iedt)=dellaqc(i,k)
                pwo_ens    (i,k,iedt)=pwo    (i,k)+edt(i)*pwdo(i,k)
             ELSE 
                dellat_ens (i,k,iedt)=0.0_r8
                dellaq_ens (i,k,iedt)=0.0_r8
                dellaqc_ens(i,k,iedt)=0.0_r8
                pwo_ens    (i,k,iedt)=0.0_r8
             END IF
          END DO
       END DO
250 END DO

  END SUBROUTINE Ensemble
  !
  !END CUP
  !-----------------------------------------------------------------------subroutines
  !*-----------
  SUBROUTINE cup_env( &
       z      ,qes    ,he     ,hes    ,t      ,q      , &
       p      ,z1     ,nCols  ,kMax   ,istart ,iend   ,psur   , &
       ierr   ,itest,den                                    )

    IMPLICIT NONE
    INTEGER, INTENT(IN   )    :: nCols
    INTEGER, INTENT(IN   )    :: kMax
    REAL(KIND=r8)   , INTENT(OUT  )    :: z   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  )    :: qes (nCols,kMax)! pressure vapor
    REAL(KIND=r8)   , INTENT(OUT  )    :: den (nCols,kMax)!density
    REAL(KIND=r8)   , INTENT(INOUT)    :: he  (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  )    :: hes (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )    :: t   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(INOUT)    :: q   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )    :: p   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )    :: z1  (nCols)
    INTEGER, INTENT(IN   )    :: istart
    INTEGER, INTENT(IN   )    :: iend
    REAL(KIND=r8)   , INTENT(IN   )    :: psur(nCols)
    INTEGER, INTENT(IN   )    :: ierr(nCols)
    INTEGER, INTENT(IN   )    :: itest
    REAL(KIND=r8)             :: esft
    !
    ! local variables
    !
    INTEGER                   :: i
    INTEGER                   :: k
    REAL(KIND=r8)                      :: tv  (nCols,kMax) ! virtual temperature
    den=0.0_r8
    DO k=1,kMax
       DO i=istart,iend
          IF(ierr(i) == 0)THEN
             !
             ! sgb - IPH is for phase, dependent on TCRIT (water or ice)
             ! calculation of the pressure vapor
             !
             esft=es5(t(i,k))
             !qes(i,k) = 0.622_r8*esft/(100.0_r8*p(i,k)-esft)
             qes(i,k) = 0.622_r8*esft/MAX((100.0_r8*p(i,k) - esft),1.0e-12_r8)
             IF(qes(i,k) <= 1.0e-08_r8  )      qes(i,k)=1.0e-08_r8
             IF(q(i,k)   >  qes(i,k))        q(i,k)=qes(i,k)
             !
             ! calculation of virtual temperature
             !
             tv(i,k) = t(i,k)+0.608_r8*q(i,k)*t(i,k)
             den(i,k)=100.0_r8*p(i,k)/(287.0_r8*tv(i,k))  !nilo
          END IF
       END DO
    END DO
    !
    ! z's are calculated with changed h's and q's and t's
    ! if itest=2
    !
    !
    ! calculate heights geopotential
    !
    IF(itest.NE.2)THEN
       DO i=istart,iend
          IF(ierr(i) == 0)THEN
             z(i,1) = MAX(0.0_r8,z1(i))-(LOG(p(i,1))-LOG(psur(i)) )*287.0_r8 &
             * tv(i,1)/9.81_r8
             he  (i,1)=9.81_r8*z(i,1)+1004.0_r8*t(i,1)+2.5e06_r8*q  (i,1)
             hes (i,1)=9.81_r8*z(i,1)+1004.0_r8*t(i,1)+2.5e06_r8*qes(i,1)
             IF(he(i,1) >= hes(i,1))he(i,1)=hes(i,1)
          END IF
       END DO
       DO k=2,kMax
          DO i=istart,iend
             IF(ierr(i) == 0)THEN
                z(i,k) = z(i,k-1)     -(LOG(p(i,k))-LOG(p(i,k-1)))*287.0_r8 &
                        * (0.5_r8*tv(i,k)+0.5_r8*tv(i,k-1) )/9.81_r8
                he  (i,k)=9.81_r8*z(i,k)+1004.0_r8*t(i,k)+2.5e06_r8*q  (i,k)
                hes (i,k)=9.81_r8*z(i,k)+1004.0_r8*t(i,k)+2.5e06_r8*qes(i,k)
                IF(he(i,k) >= hes(i,k))he(i,k)=hes(i,k)
             END IF
          END DO
       END DO
    ELSE
       DO k=1,kMax
          DO i=istart,iend
             IF(ierr(i) == 0)THEN
                z(i,k)=(he(i,k)-1004.0_r8*t(i,k)-2.5e6_r8*q(i,k))/9.81_r8
                z(i,k)=MAX(1.0e-3_r8,z(i,k))
                hes(i,k)=9.81_r8*z(i,k)+1004.0_r8*t(i,k)+2.5e06_r8*qes(i,k)
                IF(he(i,k) >= hes(i,k))he(i,k)=hes(i,k)
             END IF
          END DO
       END DO
    END IF
    RETURN
  END SUBROUTINE cup_env

  !*--------
  SUBROUTINE cup_env_clev( &
       t        ,qes      ,q      ,he       ,hes     ,z      , &
       p        ,qes_cup  ,q_cup  ,he_cup   ,hes_cup ,z_cup  ,p_cup  , &
       gamma_cup,t_cup    ,psur   ,nCols    ,kMax    ,istart ,iend   , &
       ierr     ,z1                                                 )

    IMPLICIT NONE
    INTEGER, INTENT(IN   )                  :: nCols
    INTEGER, INTENT(IN   )                  :: kMax
    INTEGER, INTENT(IN   )                  :: istart
    INTEGER, INTENT(IN   )                  :: iend
    INTEGER, INTENT(IN   )                  :: ierr(nCols)
    REAL(KIND=r8)   , INTENT(IN   )                  :: t        (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )                  :: qes      (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )                  :: q        (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )                  :: he       (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )                  :: hes      (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )                  :: z        (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )                  :: p        (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  )                  :: qes_cup  (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  )                  :: q_cup    (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  )                  :: he_cup   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  )                  :: hes_cup  (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  )                  :: z_cup    (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  )                  :: p_cup    (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  )                  :: gamma_cup(nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  )                  :: t_cup    (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )                  :: psur     (nCols)    
    REAL(KIND=r8)   , INTENT(IN   )                  :: z1       (nCols)    

    INTEGER                  :: i
    INTEGER                  :: k

  ! ierr error value, maybe modified in this routine
  ! q           = environmental mixing ratio
  ! q_cup       = environmental mixing ratio on cloud levels
  ! qes         = environmental saturation mixing ratio
  ! qes_cup     = environmental saturation mixing ratio on cloud levels
  ! t           = environmental temp
  ! t_cup       = environmental temp on cloud levels
  ! p           = environmental pressure
  ! p_cup       = environmental pressure on cloud levels
  ! z           = environmental heights
  ! z_cup       = environmental heights on cloud levels
  ! he          = environmental moist static energy
  ! he_cup      = environmental moist static energy on cloud levels
  ! hes         = environmental saturation moist static energy
  ! hes_cup     = environmental saturation moist static energy on cloud levels
  ! gamma_cup   = gamma on cloud levels
  ! psur        = surface pressure
  ! z1          = terrain elevation


    DO k=2,kMax
       DO i=istart,iend
          IF(ierr(i) == 0)THEN
             qes_cup(i,k) = 0.5_r8*(qes(i,k-1) + qes(i,k))
             q_cup  (i,k) = 0.5_r8*(  q(i,k-1) +   q(i,k))
             hes_cup(i,k) = 0.5_r8*(hes(i,k-1) + hes(i,k))
             he_cup (i,k) = 0.5_r8*( he(i,k-1) +  he(i,k))

             IF(he_cup(i,k)  >   hes_cup(i,k)) he_cup(i,k) = hes_cup(i,k)
             z_cup    (i,k) = 0.5_r8*(z(i,k-1) + z(i,k))
             p_cup    (i,k) = 0.5_r8*(p(i,k-1) + p(i,k))
             t_cup    (i,k) = 0.5_r8*(t(i,k-1) + t(i,k))
             gamma_cup(i,k) =(xl/cp)*(xl/(rv*t_cup(i,k)    &
                  *t_cup(i,k)))*qes_cup(i,k)
          END IF
       END DO
    END DO
    !
    DO i=istart,iend
       IF(ierr(i) == 0)THEN
          qes_cup  (i,1) =  qes(i,1)
          q_cup    (i,1) =    q(i,1)
          hes_cup  (i,1) =  hes(i,1)
          he_cup   (i,1) =   he(i,1)

          z_cup    (i,1) = 0.5_r8*( z(i,1) +   z1(i))
          p_cup    (i,1) = 0.5_r8*( p(i,1) + psur(i))
          t_cup    (i,1) =      t(i,1)
          gamma_cup(i,1) = xl/cp*(xl/(rv*t_cup(i,1)               &
               *t_cup(i,1)))*qes_cup(i,1)
       END IF
    END DO
  END SUBROUTINE cup_env_clev

  !*--------

  SUBROUTINE cup_maximi( &
       array    ,nCols    ,kMax      ,ks       ,ke       , &
       maxx     ,istart   ,iend     ,ierr)

    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: nCols 
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: ks
    INTEGER, INTENT(IN   ) :: istart
    INTEGER, INTENT(IN   ) :: iend
    REAL(KIND=r8)   , INTENT(IN   ) :: array(nCols, kMax)
    INTEGER, INTENT(IN   ) :: ierr (nCols) 
    INTEGER, INTENT(OUT  ) :: maxx (nCols) 
    INTEGER, INTENT(IN   ) :: ke   (nCols)

    REAL(KIND=r8)                   :: x    (nCols)
    INTEGER                :: i
    INTEGER                :: k

    DO i=istart,iend
       maxx(i)=ks
       IF(ierr(i) == 0)THEN
          x(i)=array(i,ks)
       END IF
    END DO
    DO k=1,kMax
       DO i=istart,iend
          IF(ierr(i) == 0 .AND. k >= ks .AND.  k <= ke(i) ) THEN
             IF(array(i,k) >= x(i)) THEN
                x(i)=array(i,k)
                maxx(i)=k
             END IF
          END IF
       END DO
    END DO
    RETURN
  END SUBROUTINE cup_maximi

  !*------

  SUBROUTINE cup_minimi(array     ,nCols    ,kMax       ,ks        ,kend      , &
       kt        ,istart   ,iend      ,ierr)

    IMPLICIT NONE
    INTEGER, INTENT(IN   )   :: nCols 
    INTEGER, INTENT(IN   )   :: kMax
    INTEGER, INTENT(IN   )   :: istart
    INTEGER, INTENT(IN   )   :: iend
    REAL(KIND=r8)   , INTENT(IN   )   :: array(nCols, kMax)
    INTEGER, INTENT(OUT  )   :: kt   (nCols)
    INTEGER, INTENT(IN   )   :: ks   (nCols)
    INTEGER, INTENT(IN   )   :: kend (nCols)
    INTEGER, INTENT(IN   )   :: ierr (nCols)

    REAL(KIND=r8)                     :: x    (nCols     )
    INTEGER                  :: kstop(nCols     )
    INTEGER                  :: i
    INTEGER                  :: k 

    DO i=istart,iend
       kt(i)=ks(i)
       IF(ierr(i) == 0)THEN
          x    (i)=array(i,ks(i))
          kstop(i)=MAX(ks(i)+1,kend(i))
       END IF
    END DO

    DO k=1,kMax
       DO i=istart,iend
          IF(ierr(i) == 0)THEN
             IF (k >= ks(i)+1 .AND. k <= kstop(i)) THEN
                IF(array(i,k) <  x(i)) THEN
                   x(i)=array(i,k)
                   kt(i)=k
                END IF
             END IF
          END IF
       END DO
    END DO
    RETURN
  END SUBROUTINE cup_minimi

  SUBROUTINE cup_kbcon(cap_inc   ,&
       iloop     ,k22       ,kbcon     ,he_cup    ,hes_cup   , &
       nCols     ,kMax      ,istart    ,iend      ,ierr      , &
       kbmax     ,p_cup     ,cap_max)
    IMPLICIT NONE
    INTEGER, INTENT(IN   )        :: nCols
    INTEGER, INTENT(IN   )        :: kMax
    INTEGER, INTENT(IN   )        :: istart
    INTEGER, INTENT(IN   )        :: iend
    INTEGER, INTENT(IN   )        :: iloop
    INTEGER, INTENT(OUT  )        :: kbcon   (nCols)
    INTEGER, INTENT(INOUT)        :: k22     (nCols)
    INTEGER, INTENT(INOUT)        :: ierr    (nCols)
    INTEGER, INTENT(IN   )        :: kbmax   (nCols)
    REAL(KIND=r8)   , INTENT(IN   )        :: cap_max (nCols) 
    REAL(KIND=r8)   , INTENT(IN   )        :: he_cup  (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   )        :: hes_cup (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   )        :: p_cup   (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   )        :: cap_inc
    !
    ! new
    !
    REAL(KIND=r8)                          :: plus
    REAL(KIND=r8)                          :: pbcdif
    INTEGER                       :: i
    INTEGER :: left(iend-istart+1)
    INTEGER :: nLeft
    INTEGER :: nNewLeft
    INTEGER :: toContinue(iend-istart+1)
    INTEGER :: nToContinue
    INTEGER :: cnt
    !
    ! determine the level of convective cloud base  - kbcon
    !
    kbcon=1
    nLeft = 0
    DO i=istart,iend
       IF(ierr(i) == 0 ) THEN
          kbcon(i)=k22(i)
          nLeft = nLeft + 1
          left(nLeft) = i
       ELSE
          kbcon(i)=1
       END IF
    END DO
    DO
       IF (nLeft == 0) THEN
          EXIT
       ELSE

          nNewLeft = 0
          nToContinue = 0
          !CDIR NODEP
          DO cnt = 1, nLeft
             i = left(cnt)
             IF(he_cup(i,k22(i)) <  hes_cup(i,kbcon(i))) THEN
                kbcon(i)=kbcon(i)+1
                IF(kbcon(i) >  kbmax(i)+2)THEN
                   IF (iloop <  4) THEN
                      ierr(i) =   3
                   ELSE IF (iloop == 4) THEN
                      ierr(i) = 997
                   END IF
                ELSE
                   nNewLeft = nNewLeft + 1
                   left(nNewLeft) = i
                END IF
             ELSE
                nToContinue = nToContinue + 1
                toContinue(nToContinue) = i
             END IF
          END DO

          !CDIR NODEP
          DO cnt = 1, nToContinue
             i = toContinue(cnt)
             IF(kbcon(i)-k22(i) /= 1) THEN
                !
                ! cloud base pressure and max moist static energy pressure
                !
                ! i.e., the depth (in mb) of the layer of negative buoyancy                  
                !
                pbcdif=-p_cup(i,kbcon(i))+p_cup(i,k22(i))
                plus  =MAX(25.0_r8, cap_max(i)-REAL(iloop-1,kind=r8)*cap_inc)   !new
                IF(pbcdif > plus)THEN
                   k22  (i)=k22(i)+1
                   kbcon(i)=k22(i)
                   nNewLeft = nNewLeft + 1
                   left(nNewLeft) = i
                END IF
             END IF
          END DO

          nLeft = nNewLeft

       END IF
    END DO
  END SUBROUTINE cup_kbcon


  SUBROUTINE cup_up_he(k22       ,hkb       ,z_cup     ,cd        ,entr      , &
       he_cup    ,hc        ,nCols     ,kMax      ,kbcon     , &
       ierr      ,istart    ,iend      ,dby       ,he        , &
       hes_cup,entr2D)

    IMPLICIT NONE
    INTEGER, INTENT(IN   )                       :: nCols
    INTEGER, INTENT(IN   )                       :: kMax
    INTEGER, INTENT(IN   )                       :: istart
    INTEGER, INTENT(IN   )                       :: iend
    REAL(KIND=r8)   , INTENT(IN   )                       :: entr
    INTEGER, INTENT(IN   )                       :: kbcon    (nCols)     
    INTEGER, INTENT(IN   )                       :: ierr     (nCols)     
    INTEGER, INTENT(IN   )                       :: k22      (nCols)     
    REAL(KIND=r8)   , INTENT(OUT  )                       :: hkb      (nCols)     
    REAL(KIND=r8)   , INTENT(IN   )                       :: he_cup   (nCols, kMax)
    REAL(KIND=r8)   , INTENT(OUT  )                       :: hc       (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   )                       :: z_cup    (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   )                       :: cd       (nCols, kMax)
    REAL(KIND=r8)   , INTENT(OUT  )                       :: dby      (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   )                       :: he       (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   )                       :: hes_cup  (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   )                       :: entr2D (nCols, kMax)
    REAL(KIND=r8):: entr2


    INTEGER                       :: i
    INTEGER                       :: k
    REAL(KIND=r8)                          :: dz
!
  ! hc = cloud moist static energy
  ! hkb = moist static energy at originating level
  ! he = moist static energy on model levels
  ! he_cup = moist static energy on model cloud levels
  ! hes_cup = saturation moist static energy on model cloud levels
  ! dby = buoancy term
  ! cd= detrainment function
  ! z_cup = heights of model cloud levels
  ! entr = entrainment rate

    !
    ! moist static energy inside cloud
    !
    DO i=istart,iend
       IF(ierr(i) == 0)THEN
          hkb(i)=he_cup(i,k22(i))
       END IF
    END DO

    DO k=1,kMax
       DO i=istart,iend
          IF( ierr(i) == 0 .AND. k < k22(i)  .AND. k<= kbcon(i)-1)THEN
             hc(i,k)=he_cup(i,k)
             !!!dby(i,k)=0.0_r8
             dby(i,k)=hc(i,k)-hes_cup(i,k)
          END IF
          IF(ierr(i) == 0 .AND. k >= k22(i)  .AND.k <= kbcon(i)-1)THEN
             hc(i,k)=hkb(i)
             !!!!dby(i,k)=0.0_r8
            dby(i,k)=hc(i,k)-hes_cup(i,k)
          END IF
       END DO
    END DO

    DO i=istart,iend
       IF(ierr(i) == 0)THEN
          hc (i,kbcon(i))= hkb(i)
          dby(i,kbcon(i))= hkb(i)-hes_cup(i,kbcon(i))
       END IF
    END DO
    DO k=1,kMax-1
       DO i=istart,iend
          IF(k >= 2 .AND. k >= kbcon(i).AND.ierr(i) == 0)THEN
             dz=z_cup(i,k)-z_cup(i,k-1)
             entr2=entr2D(i,k)
             hc(i,k)=(hc(i,k-1)*(1.0_r8-0.5_r8*cd(i,k)*dz)+entr2*          &
                  dz*he(i,k-1))/(1.0_r8+entr2*dz-0.5_r8*cd(i,k)*dz)

             dby(i,k)=hc(i,k)-hes_cup(i,k)
          END IF
       END DO
    END DO
    RETURN
  END SUBROUTINE cup_up_he


  SUBROUTINE cup_ktop( &
       ilo      ,dby      ,kbcon    ,ktop     ,nCols     ,kMax     , &
       istart   ,iend     ,ierr)

    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: nCols
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: istart
    INTEGER, INTENT(IN   ) :: iend
    INTEGER, INTENT(IN   ) :: ilo
    INTEGER, INTENT(INOUT) :: ierr   (nCols)
    INTEGER, INTENT(IN   ) :: kbcon  (nCols)
    INTEGER, INTENT(OUT  ) :: ktop   (nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: dby    (nCols,kMax)

    INTEGER :: i
    INTEGER :: k

    ktop (istart:iend)=1
    DO k=1,kMax-2  
       DO i=istart,iend
          IF(ierr(i) == 0 .AND. k >= kbcon(i)+1 )THEN
             IF(dby(i,k) <= 0.0_r8 .AND. ktop(i) == 1) THEN
                ktop(i)=k-1
             END IF
          END IF
       END DO
    END DO

    DO i=istart,iend      
       IF(ierr(i) == 0 .AND. ktop(i) == 1)THEN    
          IF (ilo == 1) ierr(i)=5
          IF (ilo == 2) ierr(i)=998
       END IF
    END DO

    DO  k=1,kMax
       DO i=istart,iend
          IF(ierr(i) == 0 .AND.k >= ktop(i)+1)THEN
             dby(i,k)=0.0_r8
          END IF
       END DO
    END DO
    RETURN
  END SUBROUTINE cup_ktop

  !*------

  SUBROUTINE cup_up_nms( &
       zu        ,z_cup     ,entr      ,cd        ,kbcon     , &
       ktop      ,nCols     ,kMax      ,istart    ,iend      , &
       ierr      ,k22,entr2D                                          )

    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: nCols
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: istart
    INTEGER, INTENT(IN   ) :: iend
    REAL(KIND=r8)   , INTENT(IN   ) :: entr
    REAL(KIND=r8)   , INTENT(OUT  ) :: zu   (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: z_cup(nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: cd   (nCols, kMax)
    INTEGER, INTENT(IN   ) :: kbcon(nCols  )
    INTEGER, INTENT(IN   ) :: ktop (nCols  )
    INTEGER, INTENT(IN   ) :: k22  (nCols  )
    INTEGER, INTENT(IN   ) :: ierr (nCols  )
    REAL(KIND=r8)   , INTENT(IN   ) :: entr2D (nCols, kMax)
    REAL(KIND=r8)::entr2

    INTEGER                :: i
    INTEGER                :: k
    REAL(KIND=r8)                   :: dz

    DO k=1,kMax
       DO i=istart,iend
          IF(ierr(i) == 0) THEN
             zu(i,k)=0.0_r8
          END IF
       END DO
    END DO
    DO k=1,kMax 
       DO i=istart,iend
          IF(ierr(i) == 0 .AND. k >= k22(i) .AND. k <= kbcon(i))THEN
             zu(i,k)=1.0_r8
          END IF
       END DO
    END DO
    DO k=1,kMax
       DO i=istart,iend
          IF(ierr(i) == 0 .AND. k >= kbcon(i)+1 .AND. k <= ktop(i))THEN
             dz=z_cup(i,k)-z_cup(i,k-1)
             entr2=entr2D(i,k)
             zu(i,k)=zu(i,k-1)*(1.0_r8+(entr2-cd(i,k))*dz)
          END IF
       END DO
    END DO

    RETURN
  END SUBROUTINE cup_up_nms

  !*--------

  SUBROUTINE cup_dd_nms(zd       ,z_cup     ,cdd       ,entr       ,jmin     , &
       ierr     ,nCols     ,kMax       ,istart     ,iend     , &
       itest    ,kdet      ,z1,detr2D)

    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: nCols
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: istart
    INTEGER, INTENT(IN   ) :: iend
    INTEGER, INTENT(IN   ) :: itest
    INTEGER, INTENT(IN   ) :: jmin   (nCols)
    INTEGER, INTENT(IN   ) :: ierr   (nCols)
    INTEGER, INTENT(IN   ) :: kdet   (nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: entr
    REAL(KIND=r8)   , INTENT(OUT  ) :: zd     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: z_cup  (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: cdd    (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: z1     (nCols    )
     REAL(KIND=r8)   , INTENT(IN   ) :: detr2D(nCols,kMax)
    REAL(KIND=r8)  ::entr2
    INTEGER                :: i
    INTEGER                :: k
    INTEGER                :: ki
    REAL(KIND=r8)                   :: dz
    REAL(KIND=r8)                   :: a
    REAL(KIND=r8)                   :: perc
    !
  ! z_cup = height of cloud model level
  ! z1 = terrain elevation
  ! entr = downdraft entrainment rate
  ! jmin = downdraft originating level
  ! kdet = level above ground where downdraft start detraining
  ! itest = flag to whether to calculate cdd

    !
    ! perc is the percentage of mass left when hitting the ground
    !
    !perc=0.2_r8
    perc=0.03_r8       !it is ok - new.
    DO k=1,kMax
       DO i=istart,iend
          IF(ierr(i) == 0) THEN
             zd(i,k)=0.0_r8
             IF(itest == 0)cdd(i,k)=0.0_r8
          END IF
       END DO
    END DO

    a=1.0_r8-perc

    DO i=istart,iend
       IF(ierr(i) == 0)THEN
          zd(i,jmin(i))=1.0_r8
       END IF
    END DO
    DO ki=kMax-1,1,-1
       DO i=istart,iend
          IF(ierr(i) == 0 .AND. ki <= jmin(i)-1 .AND. ki >= 1)THEN
             !
             ! integrate downward, specify detrainment(cdd)!
             !
             dz=z_cup(i,ki+1)-z_cup(i,ki)
             entr2=detr2D(i,ki)
             IF(ki <= kdet(i).AND.itest == 0)THEN
                cdd(i,ki)=entr2+(1.0_r8- (a*(z_cup(i,ki)-z1(i))          &
                     +perc*(z_cup(i,kdet(i))-z1(i)) )           &
                     /(a*(z_cup(i,ki+1)-z1(i))                  &
                     +perc*(z_cup(i,kdet(i))-z1(i))))/dz
             END IF
             zd(i,ki)=zd(i,ki+1)*(1.0_r8+(entr2-cdd(i,ki))*dz)
             !
             !----------------------
             !
          END IF
       END DO
    END DO
    RETURN
  END SUBROUTINE cup_dd_nms

  !*------

  SUBROUTINE cup_dd_he(hes_cup    ,hcd       ,z_cup     ,cdd        , &
       entr       ,jmin      ,ierr      ,nCols      ,kMax    , &
       istart     ,iend      ,he        ,dby, detr2D     )

    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: nCols
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: istart
    INTEGER, INTENT(IN   ) :: iend
    REAL(KIND=r8)   , INTENT(IN   ) :: z_cup   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: cdd     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: he      (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dby     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: hcd     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: hes_cup (nCols,kMax)
    INTEGER, INTENT(IN   ) :: jmin    (nCols)
    INTEGER, INTENT(IN   ) :: ierr    (nCols)
    INTEGER                :: i
    INTEGER                :: k
    INTEGER                :: ki
    REAL(KIND=r8)                   :: dz
    REAL(KIND=r8)   , INTENT(IN   ) :: entr
    REAL(KIND=r8)   , INTENT(IN   ) :: detr2D(nCols,kMax)
    REAL(KIND=r8):: entr2
  !-------
  ! hcd = downdraft moist static energy
  ! he = moist static energy on model levels
  ! he_cup = moist static energy on model cloud levels
  ! hes_cup = saturation moist static energy on model cloud levels
  ! dby = buoancy term
  ! cdd= detrainment function
  ! z_cup = heights of model cloud levels
  ! entr = entrainment rate
  ! zd   = downdraft normalized mass flux



    DO k=2,kMax
       DO i=istart,iend
          dby(i,k)=0.0_r8
          IF(ierr(I) == 0)THEN
             hcd(i,k)=hes_cup(i,k)
          END IF
       END DO
    END DO

    DO i=istart,iend
       IF(ierr(i) == 0)THEN
          hcd(i,jmin(i)) = hes_cup(i,jmin(i))
          dby(i,jmin(i)) = hcd    (i,jmin(i)) - hes_cup(i,jmin(i))
       END IF
    END DO

    DO ki=kMax-1,1,-1
       DO i=istart,iend
          IF(ierr(i) == 0 .AND. ki <= jmin(i)-1 )THEN
             dz        = z_cup(i,ki+1)-z_cup(i,ki)
             entr2=detr2D(i,ki)
             hcd(i,ki) = (hcd(i,ki+1)*(1.0_r8-0.5_r8*cdd(i,ki)*dz)     &
                  + entr2*dz*he(i,ki)  )                   &
                  / (1.0_r8+entr2*dz-0.5_r8*cdd(i,ki)*dz)
             dby(i,ki) = hcd(i,ki)-hes_cup(i,ki)
          END IF
       END DO
    END DO
    RETURN
  END SUBROUTINE cup_dd_he

  !*------

  SUBROUTINE cup_dd_moisture( &
       zd        ,hcd       ,hes_cup   ,qcd       , &
       qes_cup    ,pwd       ,q_cup     ,z_cup     ,cdd       , &
       entr       ,jmin      ,ierr      ,gamma_cup ,pwev      , &
       nCols      ,kMax      ,istart    ,iend      ,bu        , &
       qrcd       ,q         , &
       iloop, detr2D                                                    )

    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: nCols
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: istart
    INTEGER, INTENT(IN   ) :: iend
    INTEGER, INTENT(IN   ) :: iloop
    REAL(KIND=r8)   , INTENT(IN   ) :: entr
    INTEGER, INTENT(IN   ) :: jmin      (nCols     )
    INTEGER, INTENT(INOUT) :: ierr      (nCols     )
    REAL(KIND=r8)   , INTENT(OUT  ) :: bu        (nCols     )
    REAL(KIND=r8)   , INTENT(OUT  ) :: pwev      (nCols     )
    REAL(KIND=r8)   , INTENT(IN   ) :: zd        (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: qcd       (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: pwd       (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: qrcd      (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: hes_cup   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: hcd       (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: qes_cup   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: q_cup     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: z_cup     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: cdd       (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: gamma_cup (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: q         (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: detr2D(nCols,kMax)
    REAL(KIND=r8):: entr2

    INTEGER                :: i
    INTEGER                :: k
    INTEGER                :: ki
    REAL(KIND=r8)                   :: dz
    REAL(KIND=r8)                   :: dqeva
    REAL(KIND=r8)                   :: dh
  !------
  ! cdd= detrainment function
  ! q = environmental q on model levels
  ! q_cup = environmental q on model cloud levels
  ! qes_cup = saturation q on model cloud levels
  ! hes_cup = saturation h on model cloud levels
  ! hcd = h in model cloud
  ! bu = buoancy term
  ! zd = normalized downdraft mass flux
  ! gamma_cup = gamma on model cloud levels
  ! mentr_rate = entrainment rate
  ! qcd = cloud q (including liquid water) after entrainment
  ! qrch = saturation q in cloud
  ! pwd = evaporate at that level
  ! pwev = total normalized integrated evaoprate (I2)
  ! entr= entrainment rate


    DO i=istart,iend
       IF(ierr(i) == 0) THEN
          bu  (i) = 0.0_r8
          pwev(i) = 0.0_r8
       END IF
    END DO

    DO k=1,kMax
       DO i=istart,iend
          IF(ierr(i) == 0) THEN
             qcd (i,k) = 0.0_r8
             qrcd(i,k) = 0.0_r8
             pwd (i,k) = 0.0_r8
          END IF
       END DO
    END DO

    DO i=istart,iend
       IF(ierr(i) == 0)THEN
          k              = jmin(i)
          dz             = z_cup(i,k+1)-z_cup(i,k)
          qcd (i,k)      = q_cup(i,k)
          qrcd(i,k)      = qes_cup(i,k)
          pwd (i,jmin(i))= MIN(0.0_r8,qcd(i,k)-qrcd(i,k))
          pwev(i)        = pwev(i)+pwd(i,jmin(i))
          qcd (i,k)      = qes_cup(i,k)
          dh             = hcd(i,k)-hes_cup(i,k)
          bu(i)          = dz*dh
       END IF
    END DO

    DO ki=kMax-1,1,-1
       DO i=istart,iend
          IF(ierr(i) == 0 .AND. ki <= jmin(i)-1)THEN
             dz        = z_cup(i,ki+1)-z_cup(i,ki)
             entr2=detr2D(i,ki)
             qcd(i,ki) = (qcd(i,ki+1)*(1.0_r8-0.5_r8*cdd(i,ki)*dz)     &
                       + entr2*dz*q(i,ki)   )                   &
                       / (1.0_r8+entr2*dz-0.5_r8*cdd(i,ki)*dz)
             !
             ! to be negatively buoyant, hcd should be smaller than hes!
             !
             dh         = hcd(i,ki)-hes_cup(i,ki)
             bu  (i)    = bu(i)+dz*dh
             qrcd(i,ki) = qes_cup(i,ki)+(1.0_r8/xl)*(gamma_cup(i,ki)           &
                        / (1.0_r8+gamma_cup(i,ki)))*dh
             dqeva      = qcd(i,ki)-qrcd(i,ki)

             IF(dqeva >  0.0_r8) dqeva=0.0_r8

             pwd (i,ki) = zd  (i,ki)*dqeva
             qcd (i,ki) = qrcd(i,ki)
             pwev(i)    = pwev(i)+pwd(i,ki)
          END IF
       END DO
    END DO

    DO i=istart,iend
       IF(ierr(i) == 0)THEN
          IF(bu(i) >= 0.0_r8 .AND. iloop == 1)THEN
             ierr(i)=7
          END IF
       END IF
    END DO
    RETURN
  END SUBROUTINE cup_dd_moisture

  !*--------

  SUBROUTINE cup_up_moisture( &
       ierr       ,z_cup     ,qc        ,qrc        ,pw        , &
       pwav       ,kbcon     ,ktop      ,nCols      ,kMax      , &
       istart     ,iend      ,cd        ,dby        ,mentr_rate, &
       q          ,gamma_cup ,zu        ,qes_cup    ,k22       , &
       qe_cup, entr2D                                                    )

    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: istart
    INTEGER, INTENT(IN   ) :: iend
    INTEGER, INTENT(IN   ) :: nCols
    INTEGER, INTENT(IN   ) :: kMax
    REAL(KIND=r8)   , INTENT(IN   ) :: mentr_rate
    INTEGER, INTENT(IN   ) :: kbcon      (nCols)
    INTEGER, INTENT(IN   ) :: ktop       (nCols)
    INTEGER, INTENT(IN   ) :: ierr       (nCols)
    INTEGER, INTENT(IN   ) :: k22        (nCols)
    REAL(KIND=r8)   , INTENT(OUT  ) :: pwav       (nCols)     
    REAL(KIND=r8)   , INTENT(IN   ) :: q          (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: zu         (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: gamma_cup  (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: qe_cup     (nCols, kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: dby        (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: cd         (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: z_cup      (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: qes_cup    (nCols, kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: qc         (nCols, kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: qrc        (nCols, kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: pw         (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: entr2D(nCols, kMax)
    REAL(KIND=r8)  :: entr2


    INTEGER                    :: i
    INTEGER                    :: k
    INTEGER                    :: IALL
    REAL(KIND=r8)                       :: radius2
    REAL(KIND=r8)                       :: dz
    REAL(KIND=r8)                       :: qrch
    REAL(KIND=r8)                       :: c0
!---------
  ! cd= detrainment function
  ! q = environmental q on model levels
  ! qe_cup = environmental q on model cloud levels
  ! qes_cup = saturation q on model cloud levels
  ! dby = buoancy term
  ! cd= detrainment function
  ! zu = normalized updraft mass flux
  ! gamma_cup = gamma on model cloud levels
  ! mentr_rate = entrainment rate
  !
  ! qc = cloud q (including liquid water) after entrainment
  ! qrch = saturation q in cloud
  ! qrc = liquid water content in cloud after rainout
  ! pw = condensate that will fall out at that level
  ! pwav = totan normalized integrated condensate (I1)
  ! c0 = conversion rate (cloud to rain)



    IALL=0
    c0=0.002_r8
    !
    ! no precip for small clouds
    !
    IF(mentr_rate >  0.0_r8)THEN
       radius2=0.2_r8/mentr_rate
       IF(radius2 <  900.0_r8)c0=0.0_r8
    ENDIF


    DO i=istart,iend
       IF(ierr(i) == 0) THEN
          pwav(i)=0.0_r8
       END IF
    END DO

    DO k=1,kMax
       DO i=istart,iend
          IF(ierr(i) == 0) THEN
             pw(i,k) =0.0_r8
             !
             !snf        qc(i,k) =qes_cup(i,k)
             !
             qc(i,k) =qe_cup(i,k)   !new
             qrc(i,k)=0.0_r8
          END IF
       END DO
    END DO
    DO k=1,kMax
       DO i=istart,iend
          IF(ierr(i) == 0 .AND. k >= k22(i) .AND. k<= kbcon(i)-1)THEN
             qc(i,k)=qe_cup(i,k22(i))
          END IF
       END DO
    END DO

    !

    DO k=1,kMax
       DO  i=istart,iend
          IF(ierr(i) == 0 .AND. k >= kbcon(i) .AND. k <= ktop(i) ) THEN
             dz=z_cup(i,k)-z_cup(i,k-1)
             entr2=entr2D(i,k)
             !
             ! 1. steady state plume equation, for what could
             !    be in cloud without condensation
             !
             qc(i,k)=(qc(i,k-1)*(1.0_r8-0.5_r8*cd(i,k)*dz)+entr2*       &
                  dz*q(i,k-1))/(1.0_r8+entr2*dz-0.5_r8*cd(i,k)*dz)
             !
             !2. saturation  in cloud, this is what is allowed to be in it
             !
             qrch=qes_cup(i,k)+(1.0_r8/xl)*(gamma_cup(i,k)       &
                  /(1.0_r8+gamma_cup(i,k)))*dby(i,k)
             !
             ! liquid water content in cloud after rainout
             !
             qrc(i,k)=(qc(i,k)-qrch)/(1.0_r8+c0*dz)
             IF(qrc(i,k) <  0.0_r8)THEN
                qrc(i,k)=0.0_r8
             END IF
             !
             ! 3.Condensation
             !
             pw(i,k)=c0*dz*qrc(i,k)*zu(i,k)
             IF(IALL == 1)THEN
                qrc(i,k)=0.0_r8
                pw(i,k)=(qc(i,k)-qrch)*zu(i,k)
                IF(pw(i,k) <  0.0_r8)pw(i,k)=0.0_r8
             END IF
             !
             ! set next level
             !
             qc(i,k)=qrc(i,k)+qrch
             !
             ! integrated normalized ondensate
             !
             pwav(i)=pwav(i)+pw(i,k)
          END IF
       END DO
    END DO

    RETURN
  END SUBROUTINE cup_up_moisture

  !*---------- 

  SUBROUTINE cup_up_aa0( &
       aa0        ,z         ,zu        ,dby        ,gamma_cup , &
       t_cup      ,kbcon     ,ktop      ,kMax       ,nCols     , &
       istart     ,iend      ,ierr                              )

    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: nCols
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: istart
    INTEGER, INTENT(IN   ) :: iend
    INTEGER, INTENT(IN   ) :: kbcon     (nCols)     
    INTEGER, INTENT(IN   ) :: ktop      (nCols)     
    INTEGER, INTENT(IN   ) :: ierr      (nCols)     
    REAL(KIND=r8)   , INTENT(INOUT) :: aa0       (nCols)     
    REAL(KIND=r8)   , INTENT(IN   ) :: z         (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: zu        (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: gamma_cup (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: t_cup     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: dby       (nCols,kMax)

    REAL(KIND=r8)                         :: dz
    REAL(KIND=r8)                         :: da
    INTEGER                      :: i 
    INTEGER                      :: k

!
  ! aa0 cloud work function
  ! gamma_cup = gamma on model cloud levels
  ! t_cup = temperature (Kelvin) on model cloud levels
  ! dby = buoancy term
  ! zu= normalized updraft mass flux
  ! z = heights of model levels
  ! ierr error value, maybe modified in this routine
  !


    DO i=istart,iend
       IF(ierr(i) == 0) THEN
          aa0(i)=0.0_r8
       END IF
    END DO

    DO  k=2,kMax-1
       DO  i=istart,iend
          IF(ierr(i) == 0 .AND. k > kbcon(i) .AND. k <= ktop(i) ) THEN
             dz=z(i,k)-z(i,k-1)
             da=zu(i,k)*dz*(9.81_r8/(1004.0_r8*(   &
                  (t_cup(i,k)))))*dby(i,k-1)/ &
                  (1.0_r8+gamma_cup(i,k))
             IF (k == ktop(i) .AND. da <= 0.0_r8) THEN
                CYCLE
             ELSE
                aa0(i)=aa0(i)+da
                IF(aa0(i) <  0.0_r8)aa0(i)=0.0_r8
             END IF
          END IF
       END DO
    END DO
    RETURN
  END SUBROUTINE cup_up_aa0

  !*----------

  SUBROUTINE cup_up_CIN( &
       hh, cine        ,z         ,zu        ,dby        ,gamma_cup , &
       t_cup      ,kbcon     ,ktop      ,kMax       ,nCols     , &
       istart     ,iend      ,ierr      ,kk1,z1, buoy)

    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: nCols
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: istart
    INTEGER, INTENT(IN   ) :: iend
    INTEGER, INTENT(IN   ) :: kbcon     (nCols)
    INTEGER, INTENT(IN   ) :: ktop      (nCols)
    INTEGER, INTENT(IN   ) :: ierr      (nCols)
    INTEGER, INTENT(INOUT) :: kk1       (nCols)
    REAL(KIND=r8)   , INTENT(IN) ::z1 (nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: hh       (nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: buoy     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: cine       (nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: z         (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: zu        (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: gamma_cup (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: t_cup     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: dby       (nCols,kMax)

    REAL(KIND=r8)                         :: dz
    REAL(KIND=r8)                         :: da
    INTEGER                      :: i
    INTEGER                      :: k

!
  ! aa0 cloud work function
  ! gamma_cup = gamma on model cloud levels
  ! t_cup = temperature (Kelvin) on model cloud levels
  ! dby = buoancy term
  ! zu= normalized updraft mass flux
  ! z = heights of model levels
  ! ierr error value, maybe modified in this routine
  !

   da=0.0_r8
   cine=0.0_r8
   hh=0.0_r8
    DO i=istart,iend
          if(kk1(i)>kbcon(i))kk1(i)=kbcon(i)
    END DO

    DO  i=istart,iend
           hh(i)=z(i,kk1(i))-MAX(0.0_r8,z1(i))
           hh(i)=MIN(3500.0_r8,hh(i))
           hh(i)=MAX(1.0e-8_r8, hh(i))
    ENDDO

     DO  i=istart,iend
         DO k=kk1(i), ktop(i)-2
             dz=z(i,k)-z(i,k-1)
             da=dz*(9.81_r8/(1004.0_r8*(   &
                  (t_cup(i,k)))))*dby(i,k-1)/ &
                  (1.0_r8+gamma_cup(i,k))
             IF (da< 0.0_r8)cine(i)=cine(i)+da
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!buoy(i,k)=dby(i,k)
       END DO
             cine(i)=-cine(i)
             cine(i)=MIN(500.0_r8,cine(i))
             cine(i)=MAX(1.0e-8_r8, cine(i))
    END DO

    RETURN
  END SUBROUTINE cup_up_CIN

  !*----------

  SUBROUTINE cup_dd_edt(ierr      ,us        ,vs        ,z         ,ktop      , &
       kbcon     ,edt       ,p         ,pwav      ,pwev      , &
       nCols     ,kMax      ,istart    ,iend      ,edtmax    , &
       edtmin    ,maxens2   ,edtc      ,vshear    ,sdp       , &
       vws       ,mask      ,edtmax1,maxens22                        )


    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: nCols
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: istart
    INTEGER, INTENT(IN   ) :: iend
    INTEGER, INTENT(IN   ) :: maxens2,maxens22
    REAL(KIND=r8)   , INTENT(IN   ) :: edtmax
    REAL(KIND=r8)   , INTENT(IN   ) :: edtmin
    INTEGER, INTENT(IN   ) :: ktop   (nCols)
    INTEGER, INTENT(IN   ) :: kbcon  (nCols)
    INTEGER, INTENT(IN   ) :: ierr   (nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: us     (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: vs     (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: z      (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: p      (nCols, kMax)     
    REAL(KIND=r8)   , INTENT(OUT  ) :: edt    (nCols     )
    REAL(KIND=r8)   , INTENT(IN   ) :: pwav   (nCols     )
    REAL(KIND=r8)   , INTENT(IN   ) :: pwev   (nCols     )
    REAL(KIND=r8)   , INTENT(OUT  ) :: vshear (nCols     )
    REAL(KIND=r8)   , INTENT(OUT  ) :: sdp    (nCols     )
    REAL(KIND=r8)   , INTENT(OUT  ) :: vws    (nCols     )
    REAL(KIND=r8)   , INTENT(OUT  ) :: edtc   (nCols, maxens2)                !new 
    REAL(KIND=r8)   , INTENT(IN   ) :: edtmax1
    INTEGER, INTENT(IN   ) :: mask   (nCols)    

    REAL(KIND=r8)                         :: pefb
    REAL(KIND=r8)                         :: prezk
    REAL(KIND=r8)                         :: zkbc
    REAL(KIND=r8)                         :: pef
    REAL(KIND=r8)                         :: einc
    REAL(KIND=r8)                         :: aa1
    REAL(KIND=r8)                         :: aa2
    INTEGER                      :: i
    INTEGER                      :: kk
    INTEGER                      :: k
    !
    ! determine downdraft strength in terms of windshear
    ! calculate an average wind shear over the depth of the cloud
    !
    DO i=istart,iend
       IF(ierr(i) == 0) THEN
          edt   (i)=0.0_r8
          vws   (i)=0.0_r8
          sdp   (i)=0.0_r8
          vshear(i)=0.0_r8
       END IF
    END DO

    DO kk = 1,kMax-1
       DO i=istart,iend
          IF(ierr(i) == 0) THEN
             IF(kk  <=  min(ktop(i),kMax-1) .AND. kk  >=  kbcon(i)) THEN
                aa1=ABS((us(i,kk+1)-us(i,kk))/(z(i,kk+1)-z(i,kk)))
                aa2=ABS((vs(i,kk+1)-vs(i,kk))/(z(i,kk+1)-z(i,kk)))
                vws(i) = vws(i)+(aa1+aa2)*(p(i,kk) - p(i,kk+1))
                sdp(i) = sdp(i) + p(i,kk) - p(i,kk+1)
             END IF
             IF (kk  ==  kMax-1)vshear(i) = 1.0e3_r8 * vws(i) / sdp(i)
          END IF
       END DO
    END DO

    DO i=istart,iend
       IF(ierr(i) == 0)THEN
          pef=(1.591_r8-0.639_r8*vshear(i)+0.0953_r8*(vshear(i)**2)        &
               -0.00496_r8*(vshear(i)**3))
          !-------------------------------------------
          !snf
          IF(mask(i) == 1)THEN
             IF(pef >  edtmax1)pef=edtmax1          
          ELSE
             IF(pef >  edtmax )pef=edtmax
          END IF
          !
          !------------------
          !
          IF(pef <  edtmin)pef=edtmin
          !
          ! cloud base precip efficiency
          !
          zkbc=z(i,kbcon(i))*3.281e-3_r8
          prezk=0.02_r8
          IF(zkbc >  3.0_r8)THEN
             prezk= 0.96729352_r8+zkbc*(-0.70034167_r8+zkbc*(0.162179896_r8+zkbc*(-  &
                  1.2569798e-2_r8+zkbc*(4.2772e-4_r8-zkbc*5.44e-6_r8))))
          END IF

          IF(zkbc >  25.0_r8)THEN
             prezk=2.4_r8
          END IF
          pefb=1.0_r8/(1.0_r8+prezk)
          !
          !          if(pefb >  edtmax)pefb=edtmax
          !-------------------------------------------
          !snf
          IF(mask(i) == 1)THEN
             IF(pefb >  edtmax1)pefb=edtmax1
          ELSE
             IF(pefb >  edtmax)pefb=edtmax
          END IF
          !------------------

          IF(pefb <  edtmin)pefb=edtmin
          edt(i)=1.0_r8-0.5_r8*(pefb+pef)
          !
          ! edt here is 1-precipeff
          !
       END IF
    END DO

    DO k=1,maxens22
       DO i=istart,iend
          IF(ierr(i) == 0)THEN
             einc      = edt(i) / REAL(maxens22+1  ,kind=r8  )
             edtc(i,k) = edt(i) - REAL(k,kind=r8)*einc
             !edtc(i,1) = edt(i)*0.75_r8
             !edtc(i,2) = edt(i)*0.50_r8
!snf-new     !edtc(i,3) =edt(i)*0.25_r8
! forcando usar 0.25_r8 quando apenas 1 ensamble 
!        if(maxens22.eq.1)edtc(i,k) = edt(i)*0.25_r8
          END IF
       END DO
    END DO

    DO k=1,maxens22
       DO i=istart,iend
          IF(ierr(i) == 0)THEN
             edtc(i,k)=-edtc(i,k)*pwav(i)/pwev(i)
             !
             !             if(edtc(i,k) >  edtmax)edtc(i,k)=edtmax
             !-------------------------------------------
             !snf
             IF(mask(i) == 1)THEN
                IF(edtc(i,k) >  edtmax1)edtc(i,k)=edtmax1
             ELSE
                IF(edtc(i,k) >  edtmax )edtc(i,k)=edtmax
             END IF
             !------------------
             IF(edtc(i,k) <  edtmin)edtc(i,k)=edtmin
          END IF
       END DO
    END DO
    RETURN
  END SUBROUTINE cup_dd_edt

  !*--------

  SUBROUTINE cup_dd_aa0( &
       edt        ,ierr      ,aa0       ,jmin      ,gamma_cup , &
       t_cup      ,hcd       ,hes_cup   ,z         ,nCols     , &
       kMax       ,istart    ,iend      ,zd                     )

    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: nCols
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: istart
    INTEGER, INTENT(IN   ) :: iend
    INTEGER, INTENT(IN   ) :: jmin      (nCols)     
    INTEGER, INTENT(IN   ) :: ierr      (nCols)     
    REAL(KIND=r8)   , INTENT(IN   ) :: gamma_cup (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: t_cup     (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: z         (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: hes_cup   (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: zd        (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: hcd       (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: edt       (nCols     )
    REAL(KIND=r8)   , INTENT(INOUT) :: aa0       (nCols     )
    REAL(KIND=r8)    :: dz
    INTEGER :: i
    INTEGER :: k
    INTEGER :: kk


    DO k=1,kMax-1
       DO i=istart,iend
          IF(ierr(i) == 0 .AND. k <  jmin(i))THEN
             kk=jmin(i)-k
             !
             ! original
             !
             dz=(z(i,kk)-z(i,kk+1))
             aa0(i)=aa0(i)+zd(i,kk)*edt(i)*dz*(9.81_r8/(1004.0_r8*t_cup(i,kk)))  &
                  *((hcd(i,kk)-hes_cup(i,kk))/(1.0_r8+gamma_cup(i,kk)))
          END IF
       END DO
    END DO
    RETURN
  END SUBROUTINE cup_dd_aa0

  !*--------

  SUBROUTINE cup_dellabot( &
       he_cup    ,ierr       ,z_cup     ,p_cup     ,hcd         , &
       edt       ,zd         ,cdd       ,he        ,nCols       , &
       kMax       ,istart     ,iend      ,della    ,mentrd_rate,detr2D)

    IMPLICIT NONE
    INTEGER, INTENT(IN   )                   :: nCols
    INTEGER, INTENT(IN   )                   :: kMax
    INTEGER, INTENT(IN   )                   :: istart
    INTEGER, INTENT(IN   )                   :: iend
    REAL(KIND=r8)   , INTENT(IN   )                   :: z_cup  (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )                   :: p_cup  (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )                   :: hcd    (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )                   :: zd     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )                   :: cdd    (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )                   :: he     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  )                   :: della  (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )                   :: he_cup (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )                   :: edt    (nCols)
    INTEGER, INTENT(IN   )                   :: ierr   (nCols)
    REAL(KIND=r8)   , INTENT(IN   )                   :: mentrd_rate
    REAL(KIND=r8)   , INTENT(IN   )                   :: detr2D(nCols,kMax)
    REAL(KIND=r8):: entr2


    REAL(KIND=r8)                      :: detdo1
    REAL(KIND=r8)                      :: detdo2
    REAL(KIND=r8)                      :: entdo
    REAL(KIND=r8)                      :: g
    REAL(KIND=r8)                      :: dp
    REAL(KIND=r8)                      :: dz
    REAL(KIND=r8)                      :: subin
    REAL(KIND=r8)                      :: detdo
    INTEGER                   :: i
    g=9.81_r8
    DO i=istart,iend
       della(i,1)=0.0_r8
       IF(ierr(i) == 0) THEN
          dz        =       z_cup(i,2)-z_cup(i,1)
          entr2     = detr2D(i,1)
          dp        = 100.0_r8*(p_cup(i,1)-p_cup(i,2))
          detdo1    = edt(i)*zd(i,2)*cdd(i,1)*dz
          detdo2    = edt(i)*zd(i,1)
          entdo     = edt(i)*zd(i,2)*entr2*dz
          !snf
          subin     =-edt(i)*zd(i,2)
          detdo     = detdo1+detdo2-entdo+subin
          della(i,1)= (  detdo1*0.5_r8*(hcd(i,1)+hcd(i,2))      &
               + detdo2*    hcd(i,1)                   &
               + subin *    he_cup(i,2)                &
               - entdo *    he(i,1)    )*g/dp
       END IF
    END DO
    RETURN
  END SUBROUTINE cup_dellabot

  !*--------

  SUBROUTINE cup_dellas( &
       ierr       ,z_cup      ,p_cup     ,hcd       ,edt       , &
       zd         ,cdd        ,he        ,nCols     ,kMax      , &
       istart     ,iend       ,della     , &
       mentrd_rate,zu         ,cd        ,hc        ,ktop      , &
       k22        ,kbcon      ,mentr_rate,jmin      ,he_cup    , &
       kdet       ,kpbl       ,entr2D,       detr2D              )

    IMPLICIT NONE
    INTEGER, INTENT(IN   )           :: nCols
    INTEGER, INTENT(IN   )           :: kMax
    INTEGER, INTENT(IN   )           :: istart
    INTEGER, INTENT(IN   )           :: iend
    REAL(KIND=r8)   , INTENT(IN   )           :: z_cup (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )           :: p_cup (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )           :: hcd   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )           :: zd    (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )           :: cdd   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )           :: he    (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  )           :: della (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )           :: hc    (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )           :: cd    (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )           :: zu    (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )           :: he_cup(nCols,kMax)
    INTEGER, INTENT(IN   )           :: kbcon (nCols )
    INTEGER, INTENT(IN   )           :: ktop  (nCols )
    INTEGER, INTENT(IN   )           :: k22   (nCols )
    INTEGER, INTENT(IN   )           :: jmin  (nCols )
    INTEGER, INTENT(IN   )           :: ierr  (nCols )
    INTEGER, INTENT(IN   )           :: kdet  (nCols )
    INTEGER, INTENT(IN   )           :: kpbl  (nCols )
    REAL(KIND=r8)   , INTENT(IN   )           :: edt   (nCols )
    REAL(KIND=r8)   , INTENT(IN   )           :: mentrd_rate
    REAL(KIND=r8)   , INTENT(IN   )           :: mentr_rate
    REAL(KIND=r8)   , INTENT(IN   )           :: entr2D (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )           :: detr2D (nCols,kMax)


    REAL(KIND=r8)                        :: entdo
    REAL(KIND=r8)                        :: g
    REAL(KIND=r8)                        :: dp
    REAL(KIND=r8)                        :: dz
    REAL(KIND=r8)                        :: subin
    REAL(KIND=r8)                        :: detdo
    REAL(KIND=r8)                        :: entup
    REAL(KIND=r8)                        :: detup
    REAL(KIND=r8)                        :: subdown
    REAL(KIND=r8)                        :: entdoj
    REAL(KIND=r8)                        :: entupk
    REAL(KIND=r8)                        :: detupk
    REAL(KIND=r8)                        :: totmas
    INTEGER                     :: ier  (nCols)
    INTEGER                     :: i
    INTEGER                     :: k

    g=9.81_r8
    ier=0
    DO k=2,kMax
       DO i=istart,iend
          della(i,k)=0.0_r8
       END DO
    END DO

    DO k=2,kMax
       DO i=istart,iend
          IF(ierr(i) == 0 .AND. k <= ktop(i) ) THEN 
             !
             ! specify detrainment of downdraft, has to be consistent
             ! with zd calculations in soundd
             !
             dz    = z_cup(i,k+1)-z_cup(i,k)
             detdo = edt(i)*cdd(i,k)   *dz*zd(i,k+1)
             !! entdo = edt(i)*mentrd_rate*dz*zd(i,k+1)
             entdo = edt(i)*detr2D(i,k)*dz*zd(i,k+1)
             subin = zu(i,k+1)-zd(i,k+1)*edt(i)
             entup = 0.0_r8
             detup = 0.0_r8
             subdown = (zu(i,k)-zd(i,k)*edt(i))                 !new
             entdoj  = 0.0_r8
             entupk  = 0.0_r8
             detupk  = 0.0_r8
             !
             !         if(k >= kbcon(i))entup=mentr_rate*dz*zu(i,k)     !old
             !         subdown=(zu(i,k)-zd(i,k)*edt(i))                 !old
             !
             IF(k >= kbcon(i) .AND. k <  ktop(i))THEN
                !!entup = mentr_rate*dz*zu(i,k)
                 entup = entr2D(i,k)*dz*zu(i,k)
                detup = cd(i,k+1) *dz*zu(i,k)
             END IF

             IF(k == jmin(i))THEN
                entdoj  =edt(i)*zd(i,k)
             END IF
             !
             !         if(k == kpbl(i)-1)then                           !old
             !
             IF(k == k22(i)-1)THEN                            !new
                entupk  = zu(i,kpbl(i))
             END IF

             IF(k >  kdet(i))THEN
                detdo   = 0.0_r8
             END IF
             !
             !         if(k == ktop(i)-1)then                          !old
             !
             IF(k == ktop(i)-0)THEN                          !new
                detupk  = zu(i,ktop(i))
                subin   = 0.0_r8
             END IF
             IF(k <  kbcon(i))THEN
                detup   = 0.0_r8
             END IF
             !
             !changed due to subsidence and entrainment
             !
             totmas=subin-subdown+detup-entup-entdo+                       &
                  detdo-entupk-entdoj+detupk

             IF(ABS(totmas) >  1.e-6_r8)THEN                  !test new
                !!! nilo3 ier(i)=1
             END IF
             !--
             dp =  100.0_r8*( p_cup(i,k)-p_cup(i,k+1) )
             della(i,k)=(                                                  &
                  subin  * he_cup(i,k+1)                                 &
                  - subdown* he_cup(i,k  )                                 &
                  + detup  * 0.5_r8*( hc(i,k+1)+ hc(i,k))                     &
                  + detdo  * 0.5_r8*(hcd(i,k+1)+hcd(i,k))                     &
                  - entup  * he    (i,k)                                       &
                  - entdo  * he    (i,k)                                       &
                  - entupk * he_cup(i,k22(i))                              &
                  - entdoj * he_cup(i,jmin(i))                             &
                  + detupk * hc(i,ktop(i))                                 &
                  )*g/dp
          END IF
       END DO
    END DO
    IF (ANY(ier /=0)) THEN
       WRITE (0, '( " some ier /= 0; will stop")')
       STOP "** ERROR AT cup_dellas **"
    END IF
    RETURN
  END SUBROUTINE cup_dellas

  !*------

  SUBROUTINE cup_forcing_ens( &
       aa0       ,aa1       ,xaa0      ,mbdt      ,dtime     , &
       xmb       ,ierr      ,nCols     ,kMax      ,istart    , &
       iend      ,xf        ,name      ,mask      ,maxens    , &
       iens      ,iedt      ,maxens3   ,mconv     , &
       omeg      ,k22       ,pr_ens    ,edt       ,kbcon     , &
       ensdim    ,massfln   ,massfld   ,xff_ens3  ,xk        , &
       p_cup     ,ktop      ,ierr2     ,ierr3     ,grepar1   , &
       xfmax,maxens22,den1,dcape,cine,sens,tkeMax,tkeMedia )

    IMPLICIT NONE
    CHARACTER (LEN=*), INTENT(IN) ::  name
    INTEGER, INTENT(IN)    :: istart
    INTEGER, INTENT(IN)    :: iend
    INTEGER, INTENT(IN)    :: nCols
    INTEGER, INTENT(IN)    :: kMax
    INTEGER, INTENT(IN)    :: maxens
    INTEGER, INTENT(IN)    :: maxens3
    INTEGER, INTENT(IN)    :: ensdim
    INTEGER, INTENT(IN)    :: iens
    INTEGER, INTENT(IN)    :: iedt
    INTEGER, INTENT(IN)    :: maxens22
    INTEGER, INTENT(IN)    :: grepar1
    INTEGER, INTENT(IN)    :: mask      (nCols)
    INTEGER, INTENT(IN)    :: k22       (nCols)
    INTEGER, INTENT(IN)    :: kbcon     (nCols)
    INTEGER, INTENT(INOUT) :: ierr      (nCols)
    INTEGER, INTENT(IN)    :: ierr2     (nCols)
    INTEGER, INTENT(IN)    :: ierr3     (nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: aa0       (nCols)       
    REAL(KIND=r8)   , INTENT(IN)    :: aa1       (nCols)       
    REAL(KIND=r8)   , INTENT(OUT)   :: xmb       (nCols)       
    REAL(KIND=r8)   , INTENT(IN)    :: edt       (nCols)       
    REAL(KIND=r8)   , INTENT(IN)    :: mconv     (nCols)       
    REAL(KIND=r8)   , INTENT(IN)    :: mbdt        
    REAL(KIND=r8)   , INTENT(OUT)   :: xk        (nCols,maxens )   
    REAL(KIND=r8)   , INTENT(OUT)   :: xff_ens3  (nCols,maxens3)   
    REAL(KIND=r8)   , INTENT(IN)    :: omeg      (nCols,kMax)   
    REAL(KIND=r8)   , INTENT(IN)    :: xaa0      (nCols,maxens)
    REAL(KIND=r8)   , INTENT(OUT)   :: xf        (nCols,ensdim)
    REAL(KIND=r8)   , INTENT(IN)    :: pr_ens    (nCols,ensdim)
    REAL(KIND=r8)   , INTENT(OUT)   :: massfln   (nCols,ensdim)
    INTEGER, INTENT(IN)    :: ktop      (nCols)    
    REAL(KIND=r8)   , INTENT(IN)    :: p_cup     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN)    :: xfmax
    REAL(KIND=r8)   , INTENT(INOUT) :: massfld
    REAL(KIND=r8)   , INTENT(IN)    :: dtime
!-new
    REAL(KIND=r8)   , INTENT(IN)    :: den1 (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN)    :: dcape (nCols)
    REAL(KIND=r8)   , INTENT(IN)    :: cine (nCols)
    REAL(KIND=r8)   , INTENT(IN)    :: sens (nCols)
    REAL(KIND=r8)   , INTENT(IN)    :: tkeMax (nCols)
    REAL(KIND=r8)   , INTENT(IN)    :: tkeMedia (nCols)

    !
    ! new
    !

    !
    !--
    !
    REAL(KIND=r8)                     :: xomg(nCols) 
    REAL(KIND=r8)                     :: xff0(nCols) 
    INTEGER                  :: nens
    INTEGER                  :: n
    INTEGER                  :: nens3
    INTEGER                  :: iresult
    INTEGER                  :: iresultd
    INTEGER                  :: iresulte
    INTEGER                  :: i
    INTEGER                  :: k
    INTEGER                  :: nall(nCols,maxens)
    REAL(KIND=r8)                     :: xff_max
    !
    !------ ensemble 3 dimension = 16
    !
    REAL(KIND=r8)    :: a1,aa,bb,cc,ee,ff   
    INTEGER :: kclim(nCols) 

    REAL(KIND=r8)    :: aclim1
    REAL(KIND=r8)    :: aclim2
    REAL(KIND=r8)    :: aclim3
    REAL(KIND=r8)    :: aclim4
    LOGICAL :: teste2(nCols,maxens3)
    LOGICAL :: teste3(nCols,maxens3)
    INTEGER                  ::k_omeg

!   
!
  ! ierr error value, maybe modified in this routine
  ! pr_ens = precipitation ensemble
  ! xf_ens = mass flux ensembles
  ! massfln = downdraft mass flux ensembles used in next timestep
  ! omeg = omega from large scale model
  ! mconv = moisture convergence from large scale model
  ! zd      = downdraft normalized mass flux
  ! zu      = updraft normalized mass flux
  ! aa0     = cloud work function without forcing effects
  ! aa1     = cloud work function with forcing effects
  ! xaa0    = cloud work function with cloud effects (ensemble dependent)
  ! edt     = epsilon
  ! dir     = "storm motion"
  ! mbdt    = arbitrary numerical parameter
  ! dtime   = dt over which forcing is applied
  ! iact_gr_old = flag to tell where convection was active
  ! kbcon       = LFC of parcel from k22
  ! k22         = updraft originating level
  ! icoic       = flag if only want one closure (usually set to zero!)
  ! name        = deep or shallow convection flag


    IF(name == 'deeps' ) THEN 
       !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       !----------
       ! below new
       !-----------
       ! dded for ensemble 3 with dimension = 16
       !
       kclim(istart:iend)=0
       DO k=mkxcrt,1,-1
          DO i=istart,iend
             IF(ierr(i) == 0)THEN
                IF(p_cup(i,ktop(i)) < pcrit(k) .AND. kclim(i)==0)THEN
                   kclim(i)=k
                END IF
                IF(p_cup(i,ktop(i)) > pcrit(1)  .AND. kclim(i)==0)kclim(i)=1
             END IF
          END DO
       END DO
       !
       !--- large scale forcing
       !
       DO i=istart,iend
          xmb(i)=0.0_r8
          xomg(i) = -omeg(i,1)/9.81_r8
          IF(ierr(i) > 995)THEN
             aa0 (i)=0.0_r8
             ierr(i)=0
          END IF
       END DO

       DO k=2,kMax
          DO i=istart,iend
             IF(ierr(i) == 0 .AND. k <= kbcon(i))THEN
                IF(-omeg(i,k)/9.81_r8 > -omeg(i,1)/9.81_r8) THEN
                   xomg(i) = -omeg(i,k)/9.81_r8
                END IF
             END IF
          END DO
       END DO


       DO i=istart,iend
          IF(ierr(i) ==  0 )THEN
             !
             ! treatment different for this closure
             !
             !
             !
             !          grell scheme closure
             !---------
             !  GRE alone is better than others in GCM
             !----------
             xff0    (i)  =       (aa1(i)-aa0(i))/dtime
             xff_ens3(i,1)=       (aa1(i)-aa0(i))/dtime
             xff_ens3(i,2)=0.9_r8 * ((aa1(i)-aa0(i))/dtime) 
             xff_ens3(i,3)=1.1_r8 * ((aa1(i)-aa0(i))/dtime)
             !
             !          omeg is in bar/s, mconv done with omeg in Pa/s
             !
             !           more like brown (1979), or frank-cohen (199?)
             !---------
             ! all omegas are very bad in GCM 
             !  F4 is too bad 
             !--------- 
             !
             xff_ens3(i,4)=     -omeg(i,k22(i))/9.81_r8
             xff_ens3(i,5)=     -omeg(i,kbcon(i))/9.81_r8
             xff_ens3(i,6)=      xomg(i)               !new
             !
             ! more like Krishnamurti et al.
             ! for different betas 
             !  beta=1 is found better than others in GCM
             !  all betas are beta=1.0_r8
             !--------
             !
             xff_ens3(i,7)=mconv(i)
             xff_ens3(i,8)=mconv(i)
             xff_ens3(i,9)=mconv(i)
             !
             ! more like Fritsch Chappel or Kain Fritsch (plus triggers)
             ! here it was tested different DTs
             !--------
             !
             xff_ens3(i,10)=aa1(i)/(60.0_r8*20.0_r8)
             xff_ens3(i,11)=aa1(i)/(60.0_r8*40.0_r8)
             !xff_ens3(i,12)=aa1(i)/(60.0_r8*90.0_r8) !good
             xff_ens3(i,12)=aa1(i)/(60.0_r8*60.0_r8)
             !
             !-------------
             k= MAX(kclim(i)-1,1)
             aclim1=- dec_fudge*acrit (kclim(i))*1.0e3_r8
             aclim2=- dec_fudge*acrit (k)       *1.0e3_r8
             aclim3=- dec_fudge*acritt(kclim(i))*1.0e3_r8
             aclim4=- dec_fudge*acritt(k)       *1.0e3_r8

             xff_ens3(i,13)=MAX(0.0_r8,(AA1(I)-aclim1)/dtime)
             xff_ens3(i,14)=MAX(0.0_r8,(AA1(I)-aclim2)/dtime)
             xff_ens3(i,15)=MAX(0.0_r8,(AA1(I)-aclim3)/dtime)
             xff_ens3(i,16)=MAX(0.0_r8,(AA1(I)-aclim4)/dtime)

!
!new closures, may20,20313
!----------------------------------------------------------------
! USE DCAPE/CINE/TKEMAX
!xff_ens3=den1(i,1)*Const1*tkeMax(i)**0.5_r8*EXP(-cine/tkemax)
!xff_ens3=den1(i,1)*Const2*(2*cape(i)**0.5_r8*EXP(-cine/tkemax) 
!xff_ens3=dcape/dtime 
!xff_ens3=dcape/dtime*EXP(-cine/tkemax) 
! Const1, Const2, cloud fraction, 0.05, 0.10...
!--------------------------------------------------------------
             aa=MIN(MAX(0.01_r8,tkeMax(i)),5.0_r8)
             bb=MIN(MAX(0.0_r8,cine(i)),200.0_r8)
             cc=MAX( 0.0_r8, exp(-bb/aa))
             cc=MIN(cc,1.0_r8)
             ee=(MAX(0.0001_r8,tkeMax(i)))**0.5_r8

             xff_ens3(i,17)=MAX(den1(i,1)*0.25_r8*ee*cc,0.0_r8)  !Turbulence+cine impact
             !xff_ens3(i,17)=MAX(den1(i,1)*0.10_r8*ee*cc,0.0_r8)
             !xff_ens3(i,18)=MAX(den1(i,1)*0.10_r8*ee,0.0_r8)     !Turbulence+cine impact
!nilo             xff_ens3(i,17)=0.9_r8*MAX(dcape(i)/dtime,0.0_r8)
             xff_ens3(i,18)=1.1_r8*MAX(dcape(i)/dtime,0.0_r8)
             xff_ens3(i,19)=MAX(dcape(i)/dtime,0.0_r8)
             !xff_ens3(i,20)=MAX(dcape(i)/dtime*cc ,0.0_r8)        !
             xff_ens3(i,20)=MAX(dcape(i)/dtime,0.0_r8)        !

!----------------------------------------------------------end 
     END IF
   END DO

       DO nens=1,maxens
          DO   i=istart,iend
             IF(ierr(i) == 0)THEN
                xk(i,nens)=(xaa0(i,nens)-aa1(i))/mbdt
                IF(xk(i,nens) <= 0.0_r8 .AND. xk(i,nens) > -1.0e-6_r8)xk(i,nens)=-1.0e-6_r8
                IF(xk(i,nens) >  0.0_r8 .AND. xk(i,nens) <  1.0e-6_r8)xk(i,nens)= 1.0e-6_r8
                nall(i,nens)=(iens-1)*maxens3*maxens*maxens22       &
                            +(iedt-1)*maxens3*maxens               &
                            +(nens-1)*maxens3
             END IF
          END DO
       END DO
       !
       !  add up all ensembles
       !
       !*********************************************************************35
       !-----------------------------------------------
       ! observe the mass flux calculation:
       !-----------------------------------------------!
       ! ne   |     ierr     | mass flux               !
       ! 1    |     ierr =0  |  mf1 = xff_ens3/xk (ne) !
       ! 1    |     ierr >0  |  mf1 =  0               !
       ! 2    |     ierr2=0  |  mf2 = mf1              !
       ! 2    |     ierr2>0  |  mf2 =  0               !
       ! 3    |     ierr3=0  |  mf3 = mf1              !
       ! 3    |     ierr3>0  |  mf3 =  0               !
       ! 
       !
       ! xk(ne) is the same for any 'ne'.    
       !
       ! if ierr2 > 0 (convection was not permited for that cap_max)
       ! then equal to zero the mass flux for the second member of the ensemble (maxens)
       !
       teste2=.TRUE.
       DO nens3=1,maxens3
          DO   i=istart,iend
             IF(ierr(i) == 0 .AND. ierr2(i) > 0)THEN
                xf     (i,nall(i,2)+nens3)=0.0_r8
                massfln(i,nall(i,2)+nens3)=0.0_r8
                teste2(i,nens3)=.FALSE.
             END IF
          END DO
       END DO

       teste3=.TRUE.
       DO nens3=1,maxens3
          DO   i=istart,iend
             IF(ierr(i) == 0 .AND. ierr3(i) > 0)THEN
                xf(i,nall(i,3)+nens3)=0.0_r8
                massfln(i,nall(i,3)+nens3)=0.0_r8
                teste3(i,nens3)=.FALSE.
             END IF
          END DO
       END DO

       DO nens=1,maxens
          DO i=istart,iend
             IF(ierr(i) == 0 .AND. teste2(i,nens) .AND. teste3(i,nens) ) THEN
                !
                !---------------------------------------------------
                !! for every xk, we have maxens3 xffs,
                !! iens is from outermost ensemble (most expensive!
                !! iedt (maxens2 belongs to it)
                !! is from second, next outermost, not so expensive
                !! so, for every outermost loop, we have maxens*maxens2*3
                !! ensembles!!! nall would be 0, if everything is on first
                !! loop index, then nens would start counting, then iedt, then iensi...
                !------------------------------------------------
                !
                iresultd=0
                iresulte=0
                !
                ! check for upwind convection
                !
                iresult=0
                massfld=0.0_r8
               
                IF(xk(i,nens) <  0.0_r8 .AND. xff0(i) >  0.0_r8)iresultd=1
                iresulte=MAX(iresult,iresultd)
                IF(iresulte == 1)THEN
                   !
                   !  snf--------
                   !  xff_max=0.0_r8, 0.01_r8 0.05_r8 0.075_r8 0.10_r8
                   !  xff_max=0.05_r8     !!land good
                   !  xff_max=0.10_r8     !!ocean good
                   !
                   xff_max=0.05_r8
                   IF(mask(i) == 1)xff_max=xfmax    ! if ocean, copy value passed to subroutine (0.10_r8)   
                   !
                   !------------
                   ! xff_max=0.0_r8     original from Grell old1..increase preci over ITZ and Amazon
                   ! xff_max=0.01_r8  Pacific ZIT is strong, decrease Prec SPZ...does not increase in SACZ
                   ! xff_max=0.05_r8  SPZ is strong and eastward, weaker ZIT... 
                   ! xff_max=0.10_r8-0.20_r8  Increase precipitation on Amazon, decrease prec SPZ, SACZ And Indian
                   ! We suggest between 0.025_r8-0.10_r8
                   !
                   IF(xff0(i) >  xff_max)THEN 
                      xf(i,nall(i,nens)+1)=MAX(0.0_r8,-xff_ens3(i,1)/xk(i,nens))   &
                           +massfld
                      xf(i,nall(i,nens)+2)=MAX(0.0_r8,-xff_ens3(i,2)/xk(i,nens))   &
                           +massfld
                      xf(i,nall(i,nens)+3)=MAX(0.0_r8,-xff_ens3(i,3)/xk(i,nens))   &
                           +massfld
                      !
                      !below is new
                      !-----------
                      !
                      xf(i,nall(i,nens)+13)=MAX(0.0_r8,-xff_ens3(i,13)/xk(i,nens)) &
                           +massfld
                      xf(i,nall(i,nens)+14)=MAX(0.0_r8,-xff_ens3(i,14)/xk(i,nens)) &
                           +massfld
                      xf(i,nall(i,nens)+15)=MAX(0.0_r8,-xff_ens3(i,15)/xk(i,nens)) &
                           +massfld
                      xf(i,nall(i,nens)+16)=MAX(0.0_r8,-xff_ens3(i,16)/xk(i,nens)) &
                           +massfld
                      xf(i,nall(i,nens)+17)=MAX(0.0_r8,-xff_ens3(i,17)/xk(i,nens)) &
                           +massfld
                      xf(i,nall(i,nens)+18)=MAX(0.0_r8,-xff_ens3(i,18)/xk(i,nens)) &
                           +massfld
                      xf(i,nall(i,nens)+19)=MAX(0.0_r8,-xff_ens3(i,19)/xk(i,nens)) &
                           +massfld
                      xf(i,nall(i,nens)+20)=MAX(0.0_r8,-xff_ens3(i,20)/xk(i,nens)) &
                           +massfld


                   ELSE
                      xf(i,nall(i,nens)+1)=massfld
                      xf(i,nall(i,nens)+2)=massfld
                      xf(i,nall(i,nens)+3)=massfld
                      !
                      !below is new
                      !--------------
                      !
                      xf(i,nall(i,nens)+13)=massfld
                      xf(i,nall(i,nens)+14)=massfld
                      xf(i,nall(i,nens)+15)=massfld
                      xf(i,nall(i,nens)+16)=massfld
                      xf(i,nall(i,nens)+17)=massfld
                      xf(i,nall(i,nens)+18)=massfld
                      xf(i,nall(i,nens)+19)=massfld
                      xf(i,nall(i,nens)+20)=massfld
                   END IF

                   !
                   xf(i,nall(i,nens)+4)=MAX(0.0_r8,xff_ens3(i,4)+massfld)
                   xf(i,nall(i,nens)+5)=MAX(0.0_r8,xff_ens3(i,5)+massfld)
                   xf(i,nall(i,nens)+6)=MAX(0.0_r8,xff_ens3(i,6)+massfld)

                   !
                   a1 = MAX(1.e-9_r8,pr_ens(i,nall(i,nens)+7))
                   xf(i,nall(i,nens)+7)=MAX(0.0_r8,xff_ens3(i,7))/a1
                   a1 = MAX(1.e-9_r8,pr_ens(i,nall(i,nens)+8))
                   xf(i,nall(i,nens)+8)=MAX(0.0_r8,xff_ens3(i,8))/a1
                   a1 = MAX(1.e-9_r8,pr_ens(i,nall(i,nens)+9))
                   xf(i,nall(i,nens)+9)=MAX(0.0_r8,xff_ens3(i,9))/a1
                   !
                   IF(XK(i,nens) <  0.0_r8 .AND. xff0(i) >  xff_max)THEN

                      xf(i,nall(i,nens)+10)=MAX(0.0_r8,-xff_ens3(i,10)/xk(i,nens))+massfld
                      xf(i,nall(i,nens)+11)=MAX(0.0_r8,-xff_ens3(i,11)/xk(i,nens))+massfld
                      xf(i,nall(i,nens)+12)=MAX(0.0_r8,-xff_ens3(i,12)/xk(i,nens))+massfld
                   ELSE
                      xf(i,nall(i,nens)+10)=massfld
                      xf(i,nall(i,nens)+11)=massfld
                      xf(i,nall(i,nens)+12)=massfld
                   END IF
!new nilo2
               IF(grepar1==0)then
                      !xf(i,nall(i,nens)+1)=xf(i,nall(i,nens)+17)
                      !xf(i,nall(i,nens)+2)=xf(i,nall(i,nens)+18)
                      !xf(i,nall(i,nens)+3)=xf(i,nall(i,nens)+19)
                      !xf(i,nall(i,nens)+4)=xf(i,nall(i,nens)+20)
                      !xf(i,nall(i,nens)+5)=xf(i,nall(i,nens)+17)
                      !xf(i,nall(i,nens)+6)=xf(i,nall(i,nens)+18)
                      !xf(i,nall(i,nens)+7)=xf(i,nall(i,nens)+19)
                      !xf(i,nall(i,nens)+8)=xf(i,nall(i,nens)+20)
                      !xf(i,nall(i,nens)+9)=xf(i,nall(i,nens)+17)
                      !xf(i,nall(i,nens)+10)=xf(i,nall(i,nens)+18)
                      !xf(i,nall(i,nens)+11)=xf(i,nall(i,nens)+19)
                      !xf(i,nall(i,nens)+12)=xf(i,nall(i,nens)+20)
                      !xf(i,nall(i,nens)+13)=xf(i,nall(i,nens)+17)
                      !xf(i,nall(i,nens)+14)=xf(i,nall(i,nens)+18)
                      !xf(i,nall(i,nens)+15)=xf(i,nall(i,nens)+19)
                      !xf(i,nall(i,nens)+16)=xf(i,nall(i,nens)+20)
                      !xf(i,nall(i,nens)+17)=xf(i,nall(i,nens)+17)
                      !xf(i,nall(i,nens)+18)=xf(i,nall(i,nens)+17)
                      !xf(i,nall(i,nens)+19)=xf(i,nall(i,nens)+17)
                      !xf(i,nall(i,nens)+20)=xf(i,nall(i,nens)+17)
                 ENDIF

                   !
                END IF
             END IF
          END DO
       END DO



       DO nens=1,maxens      
          DO nens3=1,maxens3
             DO i=istart,iend
                IF(ierr(i) == 0)THEN         
                   IF(teste2(i,nens) .AND. teste3(i,nens) ) THEN
                      iresultd=0
                      iresulte=0
                      !
                      ! check for upwind convection
                      !
                      iresult=0
                      massfld=0.0_r8
                      IF(xk(i,nens) <  0.0_r8 .AND. xff0(i) >  0.0_r8)iresultd=1
                      iresulte=MAX(iresult,iresultd)
                      IF(iresulte == 1)THEN
                         !
                         !****************************************************************
                         !----- 1d closure ensemble -------------
                         !
                         IF(grepar1 >= 1 .AND. grepar1 <= maxens3)THEN
                            xf(i,nall(i,nens)+nens3)=xf(i,nall(i,nens)+grepar1)
                         END IF
                         !
                         !-------------------------
                         ! store new for next time step
                         !-------
                         !
                         massfln(i,nall(i,nens)+nens3)=edt(i)*xf(i,nall(i,nens)+nens3)
                         massfln(i,nall(i,nens)+nens3)=MAX(0.0_r8,massfln(i,nall(i,nens)+nens3))
                      END IF
                   END IF
                END IF
             END DO
          END DO
       END DO

    END IF

    IF( name /=  'deeps')THEN 
       DO n=1,ensdim
          DO i=istart,iend
             IF(ierr(i) .NE. 20 .AND. ierr(i).NE.0)THEN
                xf     (i,n)=0.0_r8
                massfln(i,n)=0.0_r8
             END IF
          END DO
       END DO
    END IF
    RETURN
  END SUBROUTINE cup_forcing_ens

  !------------------------------------------------------------------------

  SUBROUTINE cup_output_ens( &
       xf_ens    ,ierr      ,dellat    ,dellaq    ,dellaqc    , &
       outt      ,outq      ,outqc     ,pre       ,pw         , &
       xmb       ,ktop      ,nCols     ,kMax       , &
       istart    ,iend      ,maxens2    , &
       maxens    ,iens      ,pr_ens    ,outt_ens  ,maxens3    , &
       ensdim    ,massfln   ,xfac1     ,xfac_for_dn,maxens22   )

    IMPLICIT NONE
    INTEGER, INTENT(IN   )         :: nCols
    INTEGER, INTENT(IN   )         :: kMax
    INTEGER, INTENT(IN   )         :: istart
    INTEGER, INTENT(IN   )         :: iend
    INTEGER, INTENT(IN   )         :: ensdim
    INTEGER, INTENT(IN   )         :: maxens2,maxens22
    INTEGER, INTENT(IN   )         :: maxens
    INTEGER, INTENT(IN   )         :: iens
    INTEGER, INTENT(IN   )         :: maxens3
    REAL(KIND=r8)   , INTENT(OUT  )         :: pre     (nCols)            
    REAL(KIND=r8)   , INTENT(OUT  )         :: xmb     (nCols)            
    REAL(KIND=r8)   , INTENT(OUT  )         :: xfac1   (nCols)            
    REAL(KIND=r8)   , INTENT(INOUT)         :: outt    (nCols,kMax)        
    REAL(KIND=r8)   , INTENT(INOUT)         :: outq    (nCols,kMax)   
!hmjb    REAL(KIND=r8)   , INTENT(OUT  )         :: outt    (nCols,kMax)        
!hmjb    REAL(KIND=r8)   , INTENT(OUT  )         :: outq    (nCols,kMax)                 
    REAL(KIND=r8)   , INTENT(OUT  )         :: outqc   (nCols,kMax)        
    REAL(KIND=r8)   , INTENT(IN   )         :: dellat  (nCols,kMax,maxens2)
    REAL(KIND=r8)   , INTENT(IN   )         :: dellaq  (nCols,kMax,maxens2)
    REAL(KIND=r8)   , INTENT(IN   )         :: pw      (nCols,kMax,maxens2)
    REAL(KIND=r8)   , INTENT(IN   )         :: xf_ens  (nCols,ensdim)     
    REAL(KIND=r8)   , INTENT(INOUT)         :: pr_ens  (nCols,ensdim)     
    REAL(KIND=r8)   , INTENT(INOUT)         :: massfln (nCols,ensdim)     
    REAL(KIND=r8)   , INTENT(INOUT)         :: outt_ens(nCols,ensdim)
    INTEGER, INTENT(IN   )         :: ktop    (nCols)
    INTEGER, INTENT(INOUT)         :: ierr    (nCols)
    REAL(KIND=r8)   , INTENT(OUT  )         :: xfac_for_dn(nCols)
    INTEGER          :: ncount  (nCols)
    !
    ! new
    !
    INTEGER                     :: i
    INTEGER                     :: k
    INTEGER                     :: n
    REAL(KIND=r8)                        :: outtes
    REAL(KIND=r8)                        :: ddtes
    REAL(KIND=r8)                        :: dtt     (nCols,kMax)
    REAL(KIND=r8)                        :: dtq     (nCols,kMax)
    REAL(KIND=r8)                        :: dtqc    (nCols,kMax)
    REAL(KIND=r8)                        :: dtpw    (nCols,kMax)
    REAL(KIND=r8)                        :: dellaqc (nCols,kMax,maxens2)
!
! xf_ens = ensemble mass fluxes
  ! pr_ens = precipitation ensembles
  ! dellat = change of temperature per unit mass flux of cloud ensemble
  ! dellaq = change of q per unit mass flux of cloud ensemble
  ! dellaqc = change of qc per unit mass flux of cloud ensemble
  ! outtem = output temp tendency (per s)
  ! outq   = output q tendency (per s)
  ! outqc  = output qc tendency (per s)
  ! pre    = output precip
  ! xmb    = total base mass flux
  ! xfac1  = correction factor
  ! pw = pw -epsilon*pd (ensemble dependent)
  ! ierr error value, maybe modified in this routine


    DO k=1,kMax
       DO i=istart,iend
          IF(ierr(i) == 0)THEN
             outt (i,k) = 0.0_r8
             outq (i,k) = 0.0_r8
             outqc(i,k) = 0.0_r8
             dtt  (i,k) = 0.0_r8
             dtq  (i,k) = 0.0_r8
             dtqc (i,k) = 0.0_r8
             dtpw (i,k) = 0.0_r8
          END IF
       END DO
    END DO

    DO i=istart,iend
       IF(ierr(i) == 0)THEN
          pre  (i)      = 0.0_r8
          xmb  (i)      = 0.0_r8
          xfac1(i)      = 1.0_r8
          xfac_for_dn(i)= 1.0_r8
          ncount(i)     = 0
       END IF
    END DO
    !
    ! calculate mass fluxes
    !
    ! Simple average  (OLD)
    ! --------------
    !
    DO n=(iens-1)*maxens22*maxens*maxens3+1, iens*maxens22*maxens*maxens3 
       DO i=istart,iend
          IF(ierr(i) == 0)THEN
             pr_ens  (i,n) = pr_ens  (i,n)*xf_ens (i,n)
             outt_ens(i,n) = outt_ens(i,n)*xf_ens (i,n)
             IF(xf_ens(i,n) >= 0.0_r8)THEN
                xmb   (i) = xmb   (i) + xf_ens(i,n)
                ncount(i) = ncount(i) + 1
             END IF
          END IF
       END DO
    END DO
    DO i=istart,iend
       IF(ierr(i) == 0)THEN          
          IF(ncount(i) >  0)THEN
             xmb (i)=xmb(i)/REAL(ncount(i),kind=r8)
          ELSE
             xmb (i)=0.0_r8
             ierr(i)=13
          END IF
          xfac1(i)=xmb(i)!new1
       END IF
    END DO
    !
    !--------------------
    !! now do feedback 
    !----------------------
    !
    ddtes=250.0_r8   !new        For shall ddtes=500.0_r8
    !
    DO n=1,maxens22
       DO k=1,kMax
          DO i=istart,iend
             IF(ierr(i) == 0 .AND. k <= ktop(i))THEN
                dtt (i,k)  = dtt (i,k) + dellat (i,k,n)
                dtq (i,k)  = dtq (i,k) + dellaq (i,k,n)
                dtqc(i,k)  = dtqc(i,k) + dellaqc(i,k,n)
                dtpw(i,k)  = dtpw(i,k) + pw     (i,k,n)
             END IF
          END DO
       END DO
    END DO
    DO k=1,kMax
       DO i=istart,iend
          IF(ierr(i) == 0 .AND. k <= ktop(i))THEN
             outtes = dtt(i,k)*xmb(i)*86400.0_r8/REAL(maxens22,kind=r8)
             IF(outtes  >   2.0_r8*ddtes .AND. k >  2)THEN
                xmb(i) = 2.0_r8*ddtes/outtes * xmb(i)
                outtes = 1.0_r8*ddtes
             END IF


             IF(outtes  <   -ddtes)THEN
                xmb(i) = -ddtes/outtes * xmb(i)
                outtes = -ddtes
             END IF

             IF(outtes  >   0.5_r8*ddtes .AND. k <= 2)THEN
                xmb(i) =    ddtes/outtes * xmb(i)
                outtes = 0.5_r8*ddtes
             END IF
             outt (i,k) = outt (i,k) + xmb(i)*dtt (i,k)/REAL(maxens22,kind=r8)
             outq (i,k) = outq (i,k) + xmb(i)*dtq (i,k)/REAL(maxens22,kind=r8)
             outqc(i,k) = outqc(i,k) + xmb(i)*dtqc(i,k)/REAL(maxens22,kind=r8)
             pre  (i)   = pre  (i)   + xmb(i)*dtpw(i,k)/REAL(maxens22,kind=r8)
          END IF
       END DO
    END DO
    !
    ! below is new  it is only for statistics?
    !
    DO k=(iens-1)*maxens22*maxens*maxens3+1,iens*maxens22*maxens*maxens3
       DO i=istart,iend
          IF(ierr(i) == 0)THEN
             xfac1   (i) = xmb(i)        / (xfac1(i)+1.e-16_r8)
             massfln (i,k) = massfln (i,k) * xfac1(i)*xfac_for_dn(i)
             pr_ens  (i,k) = pr_ens  (i,k) * xfac1(i)
             outt_ens(i,k) = outt_ens(i,k) * xfac1(i)
          END IF
       END DO
    END DO

    RETURN
  END SUBROUTINE cup_output_ens


subroutine cup_output_ens_catt(              &
                                xf_ens     , &! intent(in   )
                                ierr       , &! intent(inout)
                                dellat     , &! intent(in   )
                                dellaq     , &! intent(in   )
                                dellaqc    , &! intent(in   )
                                outt       , &! intent(out  )
                                outq       , &! intent(out  )
                                outqc      , &! intent(out  )
                                pre        , &! intent(out  )
                                pwo_ens    , &! intent(in   )
                                xmb        , &! intent(out  )
                                ktop       , &! intent(in   )
                                iMax       , &! intent(in   )
                                kMax       , &! intent(in   ) 
                                istart     , & ! istart      (in)
                                iend       , & ! iend        (in)
                                maxens2    , &! intent(in   )
                                maxens     , &! intent(in   )
                                iens       , &! intent(in   )
                                pr_ens     , &! intent(inout)
                                outt_ens   , &! intent(inout)
                                maxens3    , &! intent(in   )
                                ensdim     , &! intent(in   )
                                massfln    , &! intent(inout)
                                xfac1      , &! intent(out  )
                                xfac_for_dn  )! intent(out  )


  implicit none
  !
  ! xf = ensemble mass fluxes
  ! pr_ens = surface precipitation ensembles
  ! dellat = change of temperature per unit mass flux of cloud ensemble
  ! dellaq = change of q per unit mass flux of cloud ensemble
  ! dellaqc = change of qc per unit mass flux of cloud ensemble
  ! outt = output temp tendency (per s)
  ! outq   = output q tendency (per s)
  ! outqc  = output qc tendency (per s)
  ! pre    = output precip
  ! xmb    = total base mass flux
  ! xfac1  = correction factor
  ! pwo_ens = pw -epsilon*pd (ensemble dependent)
  ! ierr error value, maybe modified in this routine
  !
  integer,intent(in   ) :: kMax
  integer,intent(in   ) :: iMax
  INTEGER, INTENT(IN   )         :: istart
  INTEGER, INTENT(IN   )         :: iend
  integer,intent(in   ) :: ensdim
  integer,intent(in   ) :: maxens2
  integer,intent(in   ) :: maxens
  integer,intent(in   ) :: iens
  integer,intent(in   ) :: maxens3  
  integer,intent(in   ) :: ktop    (iMax)
  integer,intent(inout) :: ierr    (iMax)
  real(KIND=r8)   ,intent(in   ) :: xf_ens  (iMax,ensdim) 
  real(KIND=r8)   ,intent(inout) :: pr_ens  (iMax,ensdim)
  real(KIND=r8)   ,intent(inout) :: outt_ens(iMax,ensdim)
  real(KIND=r8)   ,intent(inout) :: massfln (iMax,ensdim)
  real(KIND=r8)   ,intent(out  ) :: outt    (iMax,kMax)
  real(KIND=r8)   ,intent(out  ) :: outq    (iMax,kMax)
  real(KIND=r8)   ,intent(out  ) :: outqc   (iMax,kMax)
  real(KIND=r8)   ,intent(in   ) :: dellat  (iMax,kMax,maxens2)
  real(KIND=r8)   ,intent(in   ) :: dellaq  (iMax,kMax,maxens2)
  real(KIND=r8)   ,intent(in   ) :: pwo_ens (iMax,kMax,maxens2)
  real(KIND=r8)   ,intent(in   ) :: dellaqc (iMax,kMax,maxens2)
  real(KIND=r8)   ,intent(out  ) :: pre     (iMax)
  real(KIND=r8)   ,intent(out  ) :: xmb     (iMax)
  real(KIND=r8)   ,intent(out  ) :: xfac1   (iMax)
  real(KIND=r8)   ,intent(out  ) :: xfac_for_dn(iMax) 
  real(KIND=r8)   :: dti
  integer :: ncount  (iMax)
  integer :: i
  integer :: k
  integer :: n
  real(KIND=r8)    :: outtes
  real(KIND=r8)    :: ddtes
  real(KIND=r8)    :: dtt     (iMax,kMax)
  real(KIND=r8)    :: dtq     (iMax,kMax)
  real(KIND=r8)    :: dtqc    (iMax,kMax)
  real(KIND=r8)    :: dtpw    (iMax,kMax)


  real(KIND=r8)         :: xmb_ave         (iMax)
  real(KIND=r8)         :: pr_ave         (iMax)
  
  real(KIND=r8)         :: xmb_std         (iMax)
   real(KIND=r8)         :: pr_std         (iMax)
 
  real(KIND=r8)         :: xmb_ske         (iMax)
  real(KIND=r8)         :: pr_ske         (iMax)

  real(KIND=r8)         :: xmb_cur         (iMax)
   real(KIND=r8)         :: pr_cur         (iMax)
 
  real(KIND=r8)         :: x_ave        (iMax,maxens3)
  real(KIND=r8)         :: x_std        (iMax,maxens3)
  real(KIND=r8)         :: x_ske        (iMax,maxens3)
  real(KIND=r8)         :: x_cur        (iMax,maxens3)
  real(KIND=r8)         :: x_ave_cap        (iMax,maxens)
  real(KIND=r8)         :: x_ave_cap1(iMax,maxens)
  real(KIND=r8)         :: x_ave_cap2(iMax,maxens)
  real(KIND=r8)         :: x_ave_cap3(iMax,maxens)
  real(KIND=r8)         :: x_ave_cap4(iMax,maxens)
  real(KIND=r8)         :: x_ave_cap5(iMax,maxens)


  !
  !srf -dec/2003
  ! tunning parameter = use it to calibrate the precipitation field.
  !
  real(KIND=r8)    , parameter :: tunning=0.0_r8 !data tunning /0.7/
  !
  ! statistics properties
  !
  integer , parameter :: i_use_stat_prop=1

        !
        !-time average (hardwire to franl=10800 confrq=600)
        !         dti = = franl/confrq
  dti=0.05555_r8! (for franl=10800 confrq=600)

  do k=1,kMax
     do i=istart,iend
        outt  (i,k) = 0.0_r8
        outq  (i,k) = 0.0_r8
        outqc (i,k) = 0.0_r8
        dtt   (i,k) = 0.0_r8
        dtq   (i,k) = 0.0_r8
        dtqc  (i,k) = 0.0_r8
        dtpw  (i,k) = 0.0_r8
        ncount(i) = 0
     enddo
  enddo

  do i=istart,iend
     pre        (i)  =0.0_r8
     xmb        (i)  =0.0_r8
     xfac1      (i)=1.0_r8
     xfac_for_dn(i) = 1.0_r8
  enddo
  !
  !--- calculate ensemble average mass fluxes
  !
  !lufla begin -------------
  if(i_use_stat_prop == 1 ) then

      !if(.not. cup_output_vars_alloc) call alloc_cup_output_vars(iMax,maxens,maxens3)
      !
      !1st time: updraft mass fluxes statistics properties
      !          sgrell1_3d and sgrell2_3d arrays will store them for output
      !
     call massflx_stats_catt( &
          xf_ens    , & !01  intent(in   ) !01 intent(in   ) real   , intent(in   ) :: xf_ens           (iMax,ensdim)
          iMax     , & !03  intent(in   ) !03 intent(in   ) integer, intent(in   ) :: iMax
          istart    , & !  INTEGER, INTENT(IN   )    :: istart
          iend      , & !  INTEGER, INTENT(IN   )    :: iend
          ensdim    , & !06  intent(in   ) !06 intent(in   ) integer, intent(in   ) :: ensdim
          maxens    , & !07  intent(in   ) !07 intent(in   ) integer, intent(in   ) :: maxens
          maxens2   , & !08  intent(in   ) !08 intent(in   ) integer, intent(in   ) :: maxens2
          maxens3   , & !09  intent(in   ) !09 intent(in   ) integer, intent(in   ) :: maxens3
          xmb_ave   , & !10  intent(out  ) !10 intent(out  ) real   , intent(out  ) :: xt_ave           (iMax)
          xmb_std   , & !11  intent(out  ) !11 intent(out  ) real   , intent(out  ) :: xt_std           (iMax)
          xmb_cur   , & !12  intent(out  ) !12 intent(out  ) real   , intent(out  ) :: xt_cur           (iMax)
          xmb_ske   , & !13  intent(out  ) !13 intent(out  ) real   , intent(out  ) :: xt_ske           (iMax)
          x_ave     , & !14  intent(out  ) !14 intent(out  ) real   , intent(out  ) :: x_ave           (iMax,maxens3)
          x_std     , & !15  intent(out  ) !15 intent(out  ) real   , intent(out  ) :: x_std           (iMax,maxens3)
          x_ske     , & !16  intent(out  ) !16 intent(out  ) real   , intent(out  ) :: x_ske           (iMax,maxens3)
          x_cur     , & !17  intent(out  ) !17 intent(out  ) real   , intent(out  ) :: x_cur           (iMax,maxens3)
          x_ave_cap , & !18  intent(out  ) !18 intent(out  ) real   , intent(out  ) :: x_ave_cap   (iMax,maxens)
          ierr      , &        !22  intent(in   ) !22 intent(in   )  integer, intent(in   ) :: ierr            (iMax)
          !kMax        , & !25  intent(in   ) !25 intent(in   )  integer, intent(in   ) :: kMax
          1         , & !28  intent(in   ) !28 intent(in   )  integer, intent(in   ) :: itest
          x_ave_cap1, & !29  intent(out  ) !29 intent(out  ) real   , intent(out  ) :: x_ave_cap1(iMax,maxens)
          x_ave_cap2, & !30  intent(out  ) !30 intent(out  ) real   , intent(out  ) :: x_ave_cap2(iMax,maxens)
          x_ave_cap3, & !31  intent(out  ) !31 intent(out  ) real   , intent(out  ) :: x_ave_cap3(iMax,maxens)
          x_ave_cap4, & !32  intent(out  ) !32 intent(out  ) real   , intent(out  ) :: x_ave_cap4(iMax,maxens)
          x_ave_cap5, & !33  intent(out  ) !33 intent(out  ) real   , intent(out  ) :: x_ave_cap5(iMax,maxens)
          dti       )        !34  intent(in   ) !34 intent(in   )real   , intent(in   ) :: dti
     !
     !2nd time: precipitation statistics properties
     !sgrell1_3d array will store them for output
     !
     call massflx_stats_catt( &           
          pr_ens    , & !01  intent(in   ) !01 intent(in   ) real   , intent(in   ) :: xf_ens           (iMax,ensdim)
          iMax     , &        !03  intent(in   ) !03 intent(in   ) integer, intent(in   ) :: iMax
          istart    , & !  INTEGER, INTENT(IN   )    :: istart
          iend      , & !  INTEGER, INTENT(IN   )    :: iend
          ensdim    , &        !06  intent(in   ) !06 intent(in   ) integer, intent(in   ) :: ensdim
          maxens    , &        !07  intent(in   ) !07 intent(in   ) integer, intent(in   ) :: maxens        
          maxens2   , & !08  intent(in   ) !08 intent(in   ) integer, intent(in   ) :: maxens2
          maxens3   , &        !09  intent(in   ) !09 intent(in   ) integer, intent(in   ) :: maxens3
          pr_ave    , &        !10  intent(out  ) !10 intent(out  ) real   , intent(out  ) :: xt_ave           (iMax)
          pr_std    , &        !11  intent(out  ) !11 intent(out  ) real   , intent(out  ) :: xt_std           (iMax)
          pr_cur    , &        !12  intent(out  ) !12 intent(out  ) real   , intent(out  ) :: xt_cur           (iMax)
          pr_ske    , &        !13  intent(out  ) !13 intent(out  ) real   , intent(out  ) :: xt_ske           (iMax)
          x_ave     , &        !14  intent(out  ) !14 intent(out  ) real   , intent(out  ) :: x_ave           (iMax,maxens3)
          x_std     , &        !15  intent(out  ) !15 intent(out  ) real   , intent(out  ) :: x_std           (iMax,maxens3)
          x_ske     , &        !16  intent(out  ) !16 intent(out  ) real   , intent(out  ) :: x_ske           (iMax,maxens3)
          x_cur     , &        !17  intent(out  ) !17 intent(out  ) real   , intent(out  ) :: x_cur           (iMax,maxens3)
          x_ave_cap , &        !18  intent(out  ) !18 intent(out  ) real   , intent(out  ) :: x_ave_cap   (iMax,maxens)
          ierr      , &        !22  intent(in   ) !22 intent(in   )  integer, intent(in   ) :: ierr            (iMax)
          !kMax        , & !25  intent(in   ) !25 intent(in   )  integer, intent(in   ) :: kMax
          2         , & !28  intent(in   ) !28 intent(in   )  integer, intent(in   ) :: itest
          x_ave_cap1, & !29  intent(out  ) !29 intent(out  ) real   , intent(out  ) :: x_ave_cap1(iMax,maxens)
          x_ave_cap2, & !30  intent(out  ) !30 intent(out  ) real   , intent(out  ) :: x_ave_cap2(iMax,maxens)
          x_ave_cap3, & !31  intent(out  ) !31 intent(out  ) real   , intent(out  ) :: x_ave_cap3(iMax,maxens)
          x_ave_cap4, & !32  intent(out  ) !32 intent(out  ) real   , intent(out  ) :: x_ave_cap4(iMax,maxens)
          x_ave_cap5, & !33  intent(out  ) !33 intent(out  ) real   , intent(out  ) :: x_ave_cap5(iMax,maxens)
          dti       )   !34  intent(in   ) !34 intent(in   )real   , intent(in   ) :: dti


     do n=(iens-1)*maxens2*maxens*maxens3+1,iens*maxens2*maxens*maxens3
        do i=istart,iend
           if(ierr(i) == 0)then
              pr_ens  (i,n) = pr_ens  (i,n)*xf_ens(i,n)
              outt_ens(i,n) = outt_ens(i,n)*xf_ens(i,n)
              xmb     (i    ) = xmb     (i) + xf_ens(i,n)
              ncount  (i    ) = ncount  (i) + 1
           endif
        enddo
     enddo

     do i=istart,iend
        if(ierr(i) == 0)then     
           xfac_for_dn(i) = xmb(i)/real(ncount(i),kind=r8)
           !tunning process
           xmb(i) = xmb_ave(i) - tunning * xmb_std(i)
           if(xmb(i).lt.0.0_r8) xmb(i)= 0.1_r8 * xmb_ave(i)
           xfac1(i) = xmb(i)
           !srf - inconsistency at downdraft mass flux due tunning process
           !      on xmb - correction factor:
           xfac_for_dn(i) = xmb(i)/(xfac_for_dn(i) + 1.0e-16_r8)
           if(xmb_ave(i) .lt. 1.0e-10_r8) ierr(i) = 11
        endif
     enddo

  else
     !
     !----  simple average       
     !
     do n=(iens-1)*maxens2*maxens*maxens3+1,iens*maxens2*maxens*maxens3
        do i=istart,iend
           if(ierr(i) == 0)then
              pr_ens  (i,n) = pr_ens  (i,n)*xf_ens(i,n)
              outt_ens(i,n) = outt_ens(i,n)*xf_ens(i,n)
              !
              !srf- esta restricao gera incompatibilidade com o calculo de 
              !     do fluxo de massa do  downdraft calculado fora
              !
              if(xf_ens(i,n) >= 0.0_r8)then
                 xmb   (i) = xmb   (i) + xf_ens(i,n)
                 ncount(i) = ncount(i) + 1 
              endif
           endif
        enddo
     enddo

     do i=istart,iend
        if(ierr(i) == 0)then
           if(ncount(i) >= 0)then
              xmb(i)=xmb(i)/real(ncount(i),kind=r8)
           else
              xmb(i)=0.0_r8
              ierr(i)=13
           endif
        endif
        xfac1(i)=xmb(i) 
     enddo
  end if
  !
  !-- now do feedback
  !
  ddtes=250.0_r8    

  do n=1,maxens2
     do k=1,maxval(ktop)
        do i=istart,iend
           if(ierr(i) == 0.and.k <= ktop(i))then
              dtt (i,k) = dtt (i,k) + dellat (i,k,n)
              dtq (i,k) = dtq (i,k) + dellaq (i,k,n)
              dtqc(i,k) = dtqc(i,k) + dellaqc(i,k,n)
              dtpw(i,k) = dtpw(i,k) + pwo_ens(i,k,n)
           end if
        end do
     end do
  end do

  do k=1,kMax  
     do i=istart,iend
        if(ierr(i) == 0 .and.k <= ktop(i))then
           outtes = dtt(i,k)*xmb(i)*86400.0_r8/real(maxens2,kind=r8)
           if (outtes  >   2.0_r8*ddtes .and. k >  2)then
              xmb(i) = 2.0_r8*ddtes/outtes * xmb(i)
              outtes = 1.0_r8*ddtes
           endif

           if(outtes  <   -ddtes)then
              xmb(i) = -ddtes/outtes * xmb(i)
              outtes = -ddtes
           endif

           if(outtes  >   0.5_r8*ddtes .and. k <= 2)then
              xmb(i) =    ddtes/outtes * xmb(i)
              outtes = 0.5_r8*ddtes
           endif

           outt (i,k) = outt (i,k) + xmb(i)*dtt(i,k) /real(maxens2,kind=r8)
           outq (i,k) = outq (i,k) + xmb(i)*dtq(i,k) /real(maxens2,kind=r8)
           outqc(i,k) = outqc(i,k) + xmb(i)*dtqc(i,k)/real(maxens2,kind=r8)
           pre  (i)   = pre  (i)   + xmb(i)*dtpw(i,k)/real(maxens2,kind=r8) ! unit : kg[liq water]/(m^2 s)
        endif
     enddo
  enddo
  !srf 20fev2003
  !   xfac1(i) = xmb(i)/xfac1(i)  = sometimes getting nan
  ! 
  do i=istart,iend
     if(ierr(i) == 0)then
        xfac1(i)=xmb(i)/(xfac1(i)+1.0e-16_r8)
     end if
  end do

  !
  ! below is new  it is only for statistics?
  !
  do k=(iens-1)*maxens2*maxens*maxens3+1,iens*maxens2*maxens*maxens3
     do i=istart,iend
        if(ierr(i) == 0)then
           massfln (i,k) = massfln (i,k)*xfac1(i)*xfac_for_dn(i)
           pr_ens  (i,k) = pr_ens  (i,k)*xfac1(i)
           outt_ens(i,k) = outt_ens(i,k)*xfac1(i)
        end if
     end do
  end do
end subroutine cup_output_ens_catt

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

subroutine massflx_stats_catt( &
     xf_ens    , & !01 intent(in   ) real   , intent(in   ) :: xf_ens      (iMax,ensdim)
     iMax      , & !03 intent(in   ) integer, intent(in   ) :: iMax
     istart    , & !  INTEGER, INTENT(IN   )    :: istart
     iend      , & !  INTEGER, INTENT(IN   )    :: iend
     ensdim    , & !06 intent(in   ) integer, intent(in   ) :: ensdim
     maxens    , & !07 intent(in   ) integer, intent(in   ) :: maxens  
     maxens2   , & !08 intent(in   ) integer, intent(in   ) :: maxens2
     maxens3   , & !09 intent(in   ) integer, intent(in   ) :: maxens3
     xt_ave    , & !10 intent(out  ) real   , intent(out  ) :: xt_ave      (iMax)
     xt_std    , & !11 intent(out  ) real   , intent(out  ) :: xt_std      (iMax)
     xt_cur    , & !12 intent(out  ) real   , intent(out  ) :: xt_cur      (iMax)
     xt_ske    , & !13 intent(out  ) real   , intent(out  ) :: xt_ske      (iMax)
     x_ave     , & !14 intent(out  ) real   , intent(out  ) :: x_ave       (iMax,maxens3)
     x_std     , & !15 intent(out  ) real   , intent(out  ) :: x_std       (iMax,maxens3)
     x_ske     , & !16 intent(out  ) real   , intent(out  ) :: x_ske       (iMax,maxens3)
     x_cur     , & !17 intent(out  ) real   , intent(out  ) :: x_cur       (iMax,maxens3)
     x_ave_cap , & !18 intent(out  ) real   , intent(out  ) :: x_ave_cap   (iMax,maxens)
     ierr      , & !22 intent(in   )  integer, intent(in   ) :: ierr        (iMax)
!     kMax        , & !25 intent(in   )  integer, intent(in   ) :: kMax
     itest     , & !28 intent(in   )  integer, intent(in   ) :: itest
     x_ave_cap1, & !29 intent(out  ) real   , intent(out  ) :: x_ave_cap1(iMax,maxens)
     x_ave_cap2, & !30 intent(out  ) real   , intent(out  ) :: x_ave_cap2(iMax,maxens)
     x_ave_cap3, & !31 intent(out  ) real   , intent(out  ) :: x_ave_cap3(iMax,maxens)
     x_ave_cap4, & !32 intent(out  ) real   , intent(out  ) :: x_ave_cap4(iMax,maxens)
     x_ave_cap5, & !33 intent(out  ) real   , intent(out  ) :: x_ave_cap5(iMax,maxens)
     dti          )!34 intent(in   )real   , intent(in   ) :: dti

  implicit none
  integer, intent(in   ) :: iMax
  INTEGER, INTENT(IN   )    :: istart
  INTEGER, INTENT(IN   )    :: iend
  integer, intent(in   ) :: ensdim
  integer, intent(in   ) :: maxens  
  integer, intent(in   ) :: maxens2
  integer, intent(in   ) :: maxens3
!  integer, intent(in   ) :: kMax
  integer, intent(in   ) :: itest
  integer, intent(in   ) :: ierr        (iMax)
  real(KIND=r8)   , intent(in   ) :: xf_ens      (iMax,ensdim)
  real(KIND=r8)   , intent(out  ) :: xt_ave      (iMax)
  real(KIND=r8)   , intent(out  ) :: xt_std      (iMax)
  real(KIND=r8)   , intent(out  ) :: xt_ske      (iMax)
  real(KIND=r8)   , intent(out  ) :: xt_cur      (iMax)
  real(KIND=r8)   , intent(out  ) :: x_ave       (iMax,maxens3)
  real(KIND=r8)   , intent(out  ) :: x_std       (iMax,maxens3)
  real(KIND=r8)   , intent(out  ) :: x_ske       (iMax,maxens3)
  real(KIND=r8)   , intent(out  ) :: x_cur       (iMax,maxens3)
  real(KIND=r8)   , intent(out  ) :: x_ave_cap   (iMax,maxens)


  !08-03-2004
  !   extra arrays for outputting precip from cap ensembles in dependence of closures
  !
  real(KIND=r8)   , intent(out  ) :: x_ave_cap1(iMax,maxens)
  real(KIND=r8)   , intent(out  ) :: x_ave_cap2(iMax,maxens)
  real(KIND=r8)   , intent(out  ) :: x_ave_cap3(iMax,maxens)
  real(KIND=r8)   , intent(out  ) :: x_ave_cap4(iMax,maxens)
  real(KIND=r8)   , intent(out  ) :: x_ave_cap5(iMax,maxens)
  real(KIND=r8)   , intent(in   ) :: dti

  !- for output at analysis only
  real(KIND=r8) :: sgrell1_3d (100,iMax)
  real(KIND=r8) :: sgrell2_3d (100,iMax)


  real(KIND=r8)    :: small_number
  real(KIND=r8)    :: rxxx
  integer :: i
  integer :: k
  integer :: num
  integer :: kk
  integer :: num2
  integer :: num3
  integer :: iedt
  sgrell1_3d=0.0_r8
  sgrell2_3d=0.0_r8
  num =maxens *maxens2 ! old way: ensdim/maxens3
  num2=maxens2*maxens3 ! old way: ensdim/maxens
  num3=maxens2         ! old way: num2/maxens3
  small_number=1.e-6_r8


  do k=1,maxens
     do i=istart,iend
        x_ave_cap (i,k)=0.0_r8
        x_ave_cap1(i,k)=0.0_r8
        x_ave_cap2(i,k)=0.0_r8
        x_ave_cap3(i,k)=0.0_r8
        x_ave_cap4(i,k)=0.0_r8
        x_ave_cap5(i,k)=0.0_r8
     enddo
  enddo

  do k=1,maxens3
     do i=istart,iend
        x_ave(i,k)=0.0_r8
        x_std(i,k)=0.0_r8
        x_ske(i,k)=0.0_r8
        x_cur(i,k)=0.0_r8
     enddo
  enddo

  do i=istart,iend
     xt_ave(i)=0._r8
     xt_std(i)=0._r8
     xt_ske(i)=0._r8
     xt_cur(i)=0._r8
  enddo
  !
  !srf-                the mean for each cap_max (1,2 and 3) :
  !     average in maxen2 and maxens3 ensemble elements for each
  !     maxens (cap_max) ensemble element 
  do iedt=1,maxens2
     do   k=1,maxens
        do kk=1,maxens3
           do i=istart,iend
              if(ierr(i).eq.0) &
                   x_ave_cap(i,k) = x_ave_cap(i,k) + &
                   xf_ens(i,maxens3*(k-1)+(iedt-1)*maxens*maxens3+kk)
           enddo
        enddo
     enddo
  enddo
  !
  !srf-  the mean for each cap_max mean on the closure group (grell, as, kf, etc)
  !                                     and the precip efficiency
  !      average in maxen2 and maxens3 ensemble elements for each
  !      maxens (cap_max) ensemble element 
  ! x_ave_cap1(k) = fluxo de massa media sobre os tres elementos
  !                 de fechamento grell e sobre os tres elementos 
  !                 da eficiencia de precipitacao para cada 
  !                 cap_max (k=1,2,3)
  !
  do iedt=1,maxens2
     do   k=1,maxens
        do i=istart,iend
           if(ierr(i).eq.0) then
              x_ave_cap1(i,k) = x_ave_cap1(i,k) +                        & 
                   (xf_ens(i,maxens3*(k-1)+(iedt-1)*maxens*maxens3+1)+     & 
                   xf_ens(i,maxens3*(k-1)+(iedt-1)*maxens*maxens3+2)+     & 
                   xf_ens(i,maxens3*(k-1)+(iedt-1)*maxens*maxens3+3))*0.33333_r8

              x_ave_cap2(i,k) = x_ave_cap2(i,k) +                        &
                   (xf_ens(i,maxens3*(k-1)+(iedt-1)*maxens*maxens3+4)+    &
                   xf_ens(i,maxens3*(k-1)+(iedt-1)*maxens*maxens3+5)+    &
                   xf_ens(i,maxens3*(k-1)+(iedt-1)*maxens*maxens3+6))*0.33333_r8

              x_ave_cap3(i,k) = x_ave_cap3(i,k) +                       &
                   (xf_ens(i,maxens3*(k-1)+(iedt-1)*maxens*maxens3+7)+   &
                   xf_ens(i,maxens3*(k-1)+(iedt-1)*maxens*maxens3+8)+   &
                   xf_ens(i,maxens3*(k-1)+(iedt-1)*maxens*maxens3+9))*0.33333_r8

              x_ave_cap4(i,k) = x_ave_cap4(i,k) +                       &
                   (xf_ens(i,maxens3*(k-1)+(iedt-1)*maxens*maxens3+10)+  &
                   xf_ens(i,maxens3*(k-1)+(iedt-1)*maxens*maxens3+11)+  &
                   xf_ens(i,maxens3*(k-1)+(iedt-1)*maxens*maxens3+12))*0.33333_r8

              x_ave_cap5(i,k) = x_ave_cap5(i,k) +                        &
                   (xf_ens(i,maxens3*(k-1)+(iedt-1)*maxens*maxens3+13)+    &
              xf_ens(i,maxens3*(k-1)+(iedt-1)*maxens*maxens3+15)+    &
              xf_ens(i,maxens3*(k-1)+(iedt-1)*maxens*maxens3+16))*0.33333_r8
           endif
        enddo
     enddo
  enddo

  do k=1,maxens
     do i=istart,iend
        if(ierr(i).eq.0) then
           x_ave_cap (i,k)=x_ave_cap (i,k)/real(num2,kind=r8)
           x_ave_cap1(i,k)=x_ave_cap1(i,k)/real(num3,kind=r8)
           x_ave_cap2(i,k)=x_ave_cap2(i,k)/real(num3,kind=r8)
           x_ave_cap3(i,k)=x_ave_cap3(i,k)/real(num3,kind=r8)
           x_ave_cap4(i,k)=x_ave_cap4(i,k)/real(num3,kind=r8)
           x_ave_cap5(i,k)=x_ave_cap5(i,k)/real(num3,kind=r8)
        endif
     enddo
  enddo
  !
  !srf- mass flux average for each closure:
  !     average in maxens and maxens2 ensemble elements for each
  !     maxens3 (closure) ensemble element
  !     
  do kk=1,num ! num = maxens * maxens2
     do k=1,maxens3
        do i=istart,iend
           if(ierr(i).eq.0)then
              x_ave(i,k)=x_ave(i,k) + xf_ens(i,maxens3*(kk-1)+k)
           endif
        enddo
     enddo
  enddo

  do k=1,maxens3
     do i=istart,iend
        if(ierr(i).eq.0) x_ave(i,k) = x_ave(i,k)/real(num,kind=r8)
     enddo
  enddo
  !
  !srf- total average in maxens, maxen2 and maxens3 ensemble elements
  !new way
  !
  do k=1,maxens3
     do i=istart,iend
        if(ierr(i).eq.0) xt_ave(i) = xt_ave(i) + x_ave(i,k)/real(maxens3,kind=r8)
     enddo
  enddo
  !
  !old way
  !  do k=1,maxens3
  !     do i=istart,iend
  !        if(ierr(i).eq.0) xt_ave(i)=xt_ave(i)+x_ave(i,k)
  !     enddo
  !  enddo
  !  do i=istart,iend
  !     if(ierr(i).eq.0) xt_ave(i)=xt_ave(i)/real(maxens3)
  !  enddo
  !
  !--- now do std, skewness,curtosis
  !srf- stats properties in maxens and maxens2 ensemble elements for each
  !     maxens3 (closure) ensemble element (for each closure)
  !
  do kk=1,num
     do k=1,maxens3
        do i=istart,iend
           if(ierr(i).eq.0.and.x_ave(i,k).gt.0.0_r8)then
              !-----------------------    bug- underflow
              rxxx       = max( small_number, xf_ens(i,maxens3*(kk-1)+k)-x_ave(i,k))
              x_std(i,k) = x_std(i,k)+(rxxx)**2
              x_ske(i,k) = x_ske(i,k)+(rxxx)**3
              x_cur(i,k) = x_cur(i,k)+(rxxx)**4
              xt_std(i)  = xt_std(i) +(xf_ens(i,maxens3*(kk-1)+k)-xt_ave(i) )**2
           endif
        enddo
     enddo
  enddo



  do k=1,maxens3
     do i=istart,iend
        if(ierr(i).eq.0.and.xt_ave(i).gt.0.0_r8)then
           xt_ske(i) = xt_ske(i)+(x_ave(i,k)-xt_ave(i))**3
           !bug- underflow
           if(x_ave(i,k)-xt_ave(i) > small_number) &
                xt_cur(i) = xt_cur(i)+(x_ave(i,k)-xt_ave(i))**4
        endif
     enddo
  enddo
  do k=1,maxens3
     do i=istart,iend
        !srf-21/02/2005
        !       if(ierr(i).eq.0.and.x_std(i,k).gt.0.)then
        if(ierr(i).eq.0.and.x_std(i,k).gt.small_number)then
           x_std(i,k) = x_std(i,k)/real(num,kind=r8)
           x_std(i,k) = sqrt(x_std(i,k))
           if(xt_std(i) .gt. 0.0_r8 ) then
              x_ske(i,k) = x_ske(i,k)/real(num,kind=r8)/x_std(i,k)**3
              x_cur(i,k) = x_cur(i,k)/real(num,kind=r8)/x_std(i,k)**4

           endif
        endif
     enddo
  enddo

  do i=istart,iend
     if(ierr(i).eq.0.and.xt_ave(i).gt.0.0_r8)then
        xt_std(i)=xt_std(i)/real(num*maxens3,kind=r8)
        xt_std(i)=sqrt(xt_std(i))
        if(xt_std(i) .gt. small_number ) then
           xt_ske(i)=xt_ske(i)/real(maxens3,kind=r8)/xt_std(i)**3
           xt_cur(i)=xt_cur(i)/real(maxens3,kind=r8)/xt_std(i)**4
        endif
     endif
  end do
  !srf -- store at analysis files (for output)
  !-------------------for deep cumulus -----------------------

  if(itest.eq.1)then
     if(maxens3.ne.16 ) then
      !if(maxens3.ne.16 .or. iMax .lt. 30) then
        print*,'**** maxens3 for deep different of 16 *****'
        print*,'****               or                 *****'
        print*,'****       nzp   less than 30         *****'
        stop
     else
        do i=istart,iend
           !
           !this the total mean mass flux (average over all ensemble)
           !
           sgrell1_3d(2,i) = xt_ave(i)
           sgrell1_3d(3,i) = xt_std(i)
           sgrell1_3d(4,i) = xt_ske(i)
           sgrell1_3d(5,i) = xt_cur(i)
           !
           !this the mass flux for the closure type 1  - grell
           !
           sgrell1_3d(6,i) = 0.333_r8*(x_ave(i,1)+x_ave(i,2)+x_ave(i,3))
           !
           !this the mass flux for the closure type 4  - low level omega
           !
           sgrell1_3d(7,i) = 0.333_r8*(x_ave(i,4)+x_ave(i,5)+x_ave(i,6))
           !
           !this the mass flux for the closure type 7  - moisture convergence
           !
           sgrell1_3d(8,i) = 0.333_r8*(x_ave(i,7)+x_ave(i,8)+x_ave(i,9))
           !
           !this the mass flux for the closure type 10  - stability closure (fc-kf)
           !
           sgrell1_3d(9,i) = 0.333_r8*(x_ave(i,10)+x_ave(i,11)+x_ave(i,12))
           !
           !this the mass flux for the closure type 13  - arakawa-schubert
           !
           sgrell1_3d(10,i) = 0.25_r8*(x_ave(i,13)+x_ave(i,14)+x_ave(i,15)+x_ave(i,16))
           !
           !        !this the average mass flux for each cap_max ensemble element
           !        do k=1,maxens
           !           sgrell1_3d(16+(k-1),i) = x_ave_cap(i,k)
           !        enddo
           !
           sgrell2_3d(1,i)  =x_ave_cap1(i,1)
           sgrell2_3d(2,i)  =x_ave_cap1(i,2)
           sgrell2_3d(3,i)  =x_ave_cap1(i,3)
           sgrell2_3d(4,i)  =x_ave_cap2(i,1)
           sgrell2_3d(5,i)  =x_ave_cap2(i,2)
           sgrell2_3d(6,i)  =x_ave_cap2(i,3)
           sgrell2_3d(7,i)  =x_ave_cap3(i,1)
           sgrell2_3d(8,i)  =x_ave_cap3(i,2)
           sgrell2_3d(9,i)  =x_ave_cap3(i,3)
           sgrell2_3d(10,i) =x_ave_cap4(i,1)
           sgrell2_3d(11,i) =x_ave_cap4(i,2)
           sgrell2_3d(12,i) =x_ave_cap4(i,3)
           sgrell2_3d(13,i) =x_ave_cap5(i,1)
           sgrell2_3d(14,i) =x_ave_cap5(i,2)
           sgrell2_3d(15,i) =x_ave_cap5(i,3)
        end do
        !
        !this the average mass flux for each cap_max ensemble element  
        !
        do k=1,maxens
           do i=istart,iend
              if(ierr(i).eq.0) then
                 sgrell1_3d(16+(k-1),i) = x_ave_cap(i,k)
              endif
           end do
        enddo
     end if
  elseif(itest.eq.2)then
     do i=istart,iend

        sgrell1_3d(11,i) = 0.333_r8*(x_ave(i,1)+x_ave(i,2)+x_ave(i,3)) &
             *sgrell1_3d(6,i)*3600.0_r8
        sgrell1_3d(12,i) = 0.333_r8*(x_ave(i,4)+x_ave(i,5)+x_ave(i,6)) &
             *sgrell1_3d(7,i)*3600.0_r8
        sgrell1_3d(13,i) = 0.333_r8*(x_ave(i,7)+x_ave(i,8)+x_ave(i,9)) &
             *sgrell1_3d(8,i)*3600.0_r8
        sgrell1_3d(14,i) = 0.333_r8*(x_ave(i,10)+x_ave(i,11)+x_ave(i,12)) &
             *sgrell1_3d(9,i)*3600.0_r8
        sgrell1_3d(15,i) = 0.250_r8*(x_ave(i,13)+x_ave(i,14)+x_ave(i,15)+x_ave(i,16)) &
             *sgrell1_3d(10,i)*3600.0_r8 
        !
        !    convert from kg/m2/s to kg/m2/hour
        !
        sgrell2_3d(1,i)  =x_ave_cap1(i,1)*sgrell2_3d(1 ,i)*3600.0_r8
        sgrell2_3d(2,i)  =x_ave_cap1(i,2)*sgrell2_3d(2 ,i)*3600.0_r8
        sgrell2_3d(3,i)  =x_ave_cap1(i,3)*sgrell2_3d(3 ,i)*3600.0_r8
        sgrell2_3d(4,i)  =x_ave_cap2(i,1)*sgrell2_3d(4 ,i)*3600.0_r8
        sgrell2_3d(5,i)  =x_ave_cap2(i,2)*sgrell2_3d(5 ,i)*3600.0_r8
        sgrell2_3d(6,i)  =x_ave_cap2(i,3)*sgrell2_3d(6 ,i)*3600.0_r8
        sgrell2_3d(7,i)  =x_ave_cap3(i,1)*sgrell2_3d(7 ,i)*3600.0_r8
        sgrell2_3d(8,i)  =x_ave_cap3(i,2)*sgrell2_3d(8 ,i)*3600.0_r8
        sgrell2_3d(9,i)  =x_ave_cap3(i,3)*sgrell2_3d(9 ,i)*3600.0_r8
        sgrell2_3d(10,i) =x_ave_cap4(i,1)*sgrell2_3d(10,i)*3600.0_r8
        sgrell2_3d(11,i) =x_ave_cap4(i,2)*sgrell2_3d(11,i)*3600.0_r8
        sgrell2_3d(12,i) =x_ave_cap4(i,3)*sgrell2_3d(12,i)*3600.0_r8
        sgrell2_3d(13,i) =x_ave_cap5(i,1)*sgrell2_3d(13,i)*3600.0_r8
        sgrell2_3d(14,i) =x_ave_cap5(i,2)*sgrell2_3d(14,i)*3600.0_r8
        sgrell2_3d(15,i) =x_ave_cap5(i,3)*sgrell2_3d(15,i)*3600.0_r8
        !
        !-time average (hardwire to franl=10800 confrq=600)
        !         dti = = franl/confrq
        !     dti=0.05555 (for franl=10800 confrq=600)
        !
        !- precipitation closure 1 (grell, mm/h) 
        sgrell2_3d(16,i)  =sgrell2_3d(16,i)+ sgrell2_3d(1 ,i)*dti ! large cap
        sgrell2_3d(17,i)  =sgrell2_3d(17,i)+ sgrell2_3d(2 ,i)*dti ! medium cap
        sgrell2_3d(18,i)  =sgrell2_3d(18,i)+ sgrell2_3d(3 ,i)*dti ! low cap'
        !
        !- precipitation closure 2 (omega, mm/h) 
        !
        sgrell2_3d(19,i)  =sgrell2_3d(19,i)+ sgrell2_3d(4 ,i)*dti ! large cap
        sgrell2_3d(20,i)  =sgrell2_3d(20,i)+ sgrell2_3d(5 ,i)*dti ! medium cap
        sgrell2_3d(21,i)  =sgrell2_3d(21,i)+ sgrell2_3d(6 ,i)*dti ! low cap'
        !
        !- precipitation closure 3 (kuo, mm/h) 
        !
        sgrell2_3d(22,i)  =sgrell2_3d(22,i)+ sgrell2_3d(7 ,i)*dti ! large cap
        sgrell2_3d(23,i)  =sgrell2_3d(23,i)+ sgrell2_3d(8 ,i)*dti ! medium cap
        sgrell2_3d(24,i)  =sgrell2_3d(24,i)+ sgrell2_3d(9 ,i)*dti ! low cap'
        !
        !- precipitation closure 4 (k&f, mm/h) 
        !
        sgrell2_3d(25,i)  =sgrell2_3d(25,i)+ sgrell2_3d(10,i)*dti ! large cap
        sgrell2_3d(26,i)  =sgrell2_3d(26,i)+ sgrell2_3d(11,i)*dti ! medium cap
        sgrell2_3d(27,i)  =sgrell2_3d(27,i)+ sgrell2_3d(12,i)*dti ! low cap'
        !
        !- precipitation closure 4 (a&s, mm/h) 
        !
        sgrell2_3d(28,i)  =sgrell2_3d(28,i)+ sgrell2_3d(13,i)*dti ! large cap
        sgrell2_3d(29,i)  =sgrell2_3d(29,i)+ sgrell2_3d(14,i)*dti ! medium cap
        sgrell2_3d(30,i)  =sgrell2_3d(30,i)+ sgrell2_3d(15,i)*dti ! low cap'
        !---------------for shallow cumulus ---------------------------------
     end do
  elseif(itest.eq.3)then
     if(maxens3.ne.10) then
        print*,'**** maxens3 for shallow different of 10 *****'
        stop
     else
        !
        !this the mass flux for the closure type 1  - grell
        !
        do i=istart,iend
           sgrell1_3d(21,i) = 0.333_r8*(x_ave(i,1)+x_ave(i,2)+x_ave(i,3))
           !this the mass flux for the closure type 4 -  arakawa-schubert
           sgrell1_3d(22,i) = 0.025_r8*(x_ave(i,4)+x_ave(i,5)+x_ave(i,6)+x_ave(i,7) )
           !this the mass flux for the closure type 7 - stability closure (fc-kf)
           sgrell1_3d(23,i) = 0.333_r8*(x_ave(i,8)+x_ave(i,9)+x_ave(i,10))
           !this the total mean mass flux (average ovell all ensemble)
           !       sgrell1_3d(24,i) = xt_ave(i)
           !       sgrell1_3d(25,i) = xt_std(i)
           !       sgrell1_3d(26,i) = xt_ske(i)
           !       sgrell1_3d(27,i) = xt_cur(i)

        enddo
     endif
  endif
end subroutine massflx_stats_catt










  !---------------------------------
  REAL(KIND=r8) FUNCTION es5(t)
    REAL(KIND=r8), INTENT(IN) :: t

    IF (t <= tcrit) THEN
       es5 = EXP(ae(2)-be(2)/t)
    ELSE
       es5 = EXP(ae(1)-be(1)/t)
    END IF
  END FUNCTION es5
END MODULE Cu_Grellens
