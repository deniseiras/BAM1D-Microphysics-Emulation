MODULE Cu_RAS3PHASE
  USE Utils, ONLY: &
       IJtoIBJB ,&
       NearestIJtoIBJB,&
       linearijtoibjb
  
  USE Options, ONLY: &
       reducedGrid

  USE Mod_GET_PRS, ONLY: GET_PRS
  IMPLICIT NONE
  PRIVATE
  INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(13,60) ! the '60' maps to 64-bit real
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)   ! Kind for 32-bits Integer Numbers

  !  --- ...  Geophysics/Astronomy constants

  REAL(kind=r8),PARAMETER:: con_g      =9.80665e+0_r8     ! gravity           (m/s2)
  REAL(kind=r8),PARAMETER:: con_pi     =3.1415926535897931 ! pi
  REAL(kind=r8),PARAMETER:: con_rerth  =6.3712e+6      ! radius of earth   (m)
 
  !  --- ...  Thermodynamics constants


  REAL(kind=r8),PARAMETER:: con_cp     =1.0046e+3_r8      ! spec heat air @p    (J/kg/K)
  REAL(kind=r8),PARAMETER:: con_rv     =4.6150e+2_r8      ! gas constant H2O    (J/kg/K)
  REAL(kind=r8),PARAMETER:: con_hvap   =2.5000e+6_r8      ! lat heat H2O cond   (J/kg)
  REAL(kind=r8),PARAMETER:: con_hfus   =3.3358e+5_r8      ! lat heat H2O fusion (J/kg)
  REAL(kind=r8),PARAMETER:: con_rd     =2.8705e+2_r8      ! gas constant air    (J/kg/K)
  REAL(kind=r8),PARAMETER:: con_cvap   =1.8460e+3_r8      ! spec heat H2O gas   (J/kg/K)
  REAL(kind=r8),PARAMETER:: con_cliq   =4.1855e+3_r8      ! spec heat H2O liq   (J/kg/K)
  REAL(kind=r8),PARAMETER:: con_csol   =2.1060e+3_r8      ! spec heat H2O ice   (J/kg/K)
  REAL(kind=r8),PARAMETER:: con_ttp    =2.7316e+2_r8      ! temp at H2O 3pt     (K)
  REAL(kind=r8),PARAMETER:: con_psat   =6.1078e+2_r8      ! pres at H2O 3pt     (Pa)  


  !  Secondary constants

  REAL(kind=r8),PARAMETER:: con_rocp   =con_rd/con_cp
  REAL(kind=r8),PARAMETER:: con_fvirt  =con_rv/con_rd-1.0_r8
  REAL(kind=r8),PARAMETER:: con_eps    =con_rd/con_rv
  REAL(kind=r8),PARAMETER:: con_epsm1  =con_rd/con_rv-1.0_r8

  !     module module_ras

  !     real, parameter :: con_FVirt=0.0
  !
  INTEGER, PARAMETER :: nrcmax=15 ! Maximum # of random clouds per 1200s
  !INTEGER, PARAMETER :: nrcmax=12 ! Maximum # of random clouds per 1200s
  !     integer, parameter :: nrcmax=15 ! Maximum # of random clouds per 1200s
  !     integer, parameter :: nrcmax=20
  INTEGER     ::  nrcm!     nrcm     - integer, number of random clouds                  1    !
  REAL (kind=r8), ALLOCATABLE :: xkt2(:,:,:)

  REAL (kind=r8), PARAMETER :: delt_c=1800.0_r8
  LOGICAL       , PARAMETER :: fix_ncld_hr=.TRUE.
  !
  REAL(kind=r8), PARAMETER :: ZERO=0.0_r8
  REAL(kind=r8), PARAMETER :: HALF=0.5_r8
  REAL(kind=r8), PARAMETER :: ONE=1.0_r8
  REAL(kind=r8), PARAMETER :: TWO=2.0_r8
  REAL(kind=r8), PARAMETER :: FOUR_P2=4.E2_r8
  REAL(kind=r8), PARAMETER :: FOUR=4.0_r8
  REAL(kind=r8), PARAMETER :: ONE_M1=1.E-1_r8
  REAL(kind=r8), PARAMETER :: ONE_M2=1.E-2_r8
  REAL(kind=r8), PARAMETER :: ONE_M5=1.E-5_r8
  REAL(kind=r8), PARAMETER :: ONE_M6=1.E-6_r8
  REAL(kind=r8), PARAMETER :: ONE_M10=1.E-10_r8
  !
  REAL(kind=r8), PARAMETER :: cmb2pa = 100.0_r8  ! Conversion from MB to PA
  REAL(kind=r8), PARAMETER :: onebg       = ONE / con_g
  REAL(kind=r8), PARAMETER :: gravcon     = cmb2pa * ONEBG
  REAL(kind=r8), PARAMETER :: gravfac     = con_g / CMB2PA
  REAL(kind=r8), PARAMETER :: elocp       = con_hvap  / con_cp
  REAL(kind=r8), PARAMETER :: elfocp      = (con_hvap +con_hfus) / con_cp
  REAL(kind=r8), PARAMETER :: rkapi       = ONE / con_rocp
  REAL(kind=r8), PARAMETER :: rkpp1i      = ONE / (ONE+con_rocp)
  REAL(kind=r8), PARAMETER :: zfac        = 0.28888889E-4_r8 * ONEBG
  !
  !
  !     logical, parameter :: advcld=.false. advups=.true.
  !     logical, parameter :: advcld=.true., advups=.true., advtvd=.false.
  LOGICAL, PARAMETER :: advcld=.TRUE.
  LOGICAL, PARAMETER :: advups=.FALSE.
  LOGICAL, PARAMETER :: advtvd=.TRUE.
  LOGICAL, PARAMETER :: flipv = .TRUE.!     flipv    - logical, flag for vertical direction              1    !

  !     logical, parameter :: advcld=.false., advups=.false.
  !
  REAL(kind=r8), ALLOCATABLE  ::  RASAL(:)
  REAL(kind=r8) :: AFC
  REAL(kind=r8) :: facdt

  !     PARAMETER (DD_DP=1000.0, RKNOB=1.0, EKNOB=1.0)   ! No downdraft!
  !     PARAMETER (DD_DP=100.0,  RKNOB=1.0, EKNOB=1.0)
  !     PARAMETER (DD_DP=200.0,  RKNOB=1.0, EKNOB=1.0)
  !     PARAMETER (DD_DP=250.0,  RKNOB=1.0, EKNOB=1.0)
  !     PARAMETER (DD_DP=300.0,  RKNOB=1.0, EKNOB=1.0)
  !     PARAMETER (DD_DP=450.0,  RKNOB=1.0, EKNOB=1.0)
  !     PARAMETER (DD_DP=500.0,  RKNOB=0.5, EKNOB=1.0)
  !     PARAMETER (DD_DP=500.0,  RKNOB=0.70, EKNOB=1.0)
  !     PARAMETER (DD_DP=500.0,  RKNOB=0.75, EKNOB=1.0)
  REAL(kind=r8), PARAMETER :: DD_DP=500.0_r8,  RKNOB=1.0_r8, EKNOB=1.0_r8
  !     PARAMETER (DD_DP=500.0,  RKNOB=1.5, EKNOB=1.0)
!!!!! PARAMETER (DD_DP=450.0,  RKNOB=1.5, EKNOB=1.0)
  !     PARAMETER (DD_DP=450.0,  RKNOB=2.0, EKNOB=1.0)
  !     PARAMETER (DD_DP=450.0,  RKNOB=0.5, EKNOB=1.0)
  !     PARAMETER (DD_DP=350.0,  RKNOB=0.5, EKNOB=1.0)
  !     PARAMETER (DD_DP=350.0,  RKNOB=1.0, EKNOB=1.0)
  !     PARAMETER (DD_DP=350.0,  RKNOB=2.0, EKNOB=1.0)
  !     PARAMETER (DD_DP=350.0,  RKNOB=3.0, EKNOB=1.0)
  !
  REAL(kind=r8), PARAMETER :: RHMAX=1.0_r8     !  MAX RELATIVE HUMIDITY
  REAL(kind=r8), PARAMETER :: QUAD_LAM=1.0_r8  !  MASK FOR QUADRATIC LAMBDA
  !     PARAMETER (RHRAM=0.15)    !  PBL RELATIVE HUMIDITY RAMP
  REAL(kind=r8), PARAMETER :: RHRAM=0.05_r8    !  PBL RELATIVE HUMIDITY RAMP
  !     PARAMETER (RHRAM=0.10)    !  PBL RELATIVE HUMIDITY RAMP
  REAL(kind=r8), PARAMETER :: HCRIT=4000.0_r8  !  Critical Moist Static Energy
  REAL(kind=r8), PARAMETER :: qudfac=QUAD_LAM*half
  !     parameter (qudfac=QUAD_LAM*0.25)    ! Yogesh's
  REAL(kind=r8), PARAMETER :: testmb=0.1_r8
  REAL(kind=r8), PARAMETER :: tstmbi=one/testmb
  !
  !
  !     PARAMETER (ALMIN1=0.00E-6, ALMIN2=0.00E-5, ALMAX=1.0E-1)
  !     PARAMETER (ALMIN1=0.00E-6, ALMIN2=0.00E-5, ALMAX=2.0E-2)
  !     PARAMETER (ALMIN1=0.00E-6, ALMIN2=1.00E-6, ALMAX=2.0E-2)
  !     PARAMETER (ALMIN1=5.00E-6, ALMIN2=2.50E-5, ALMAX=2.0E-2)
  !     PARAMETER (ALMIN1=0.00E-6, ALMIN2=2.50E-5, ALMAX=2.0E-2)
!!!   PARAMETER (ALMIN1=0.00E-6, ALMIN2=2.50E-5, ALMAX=1.0E-2)
  !!    PARAMETER (ALMIN1=0.00E-6, ALMIN2=2.50E-5, ALMAX=1.0E-3)
  !     PARAMETER (ALMIN1=0.00E-6, ALMIN2=1.00E-5, ALMAX=1.0E-2)
  !     PARAMETER (ALMIN1=0.00E-6, ALMIN2=2.00E-5, ALMAX=1.0E-2)
  !     PARAMETER (ALMIN1=0.00E-6, ALMIN2=2.50E-5, ALMAX=1.0E-2)
  REAL(kind=r8), PARAMETER :: ALMIN1=0.00E-6_r8, ALMIN2=4.00E-5_r8, ALMAX=1.0E-2_r8
  !cnt  PARAMETER (ALMIN1=0.00E-6, ALMIN2=2.50E-5, ALMAX=5.0E-3)
  !LL   PARAMETER (ALMIN1=0.00E-6, ALMIN2=2.50E-5, ALMAX=4.0E-3)
  !     PARAMETER (ALMIN1=0.00E-6, ALMIN2=1.00E-5, ALMAX=2.0E-2)
  !     PARAMETER (ALMIN1=0.00E-6, ALMIN2=5.00E-4, ALMAX=2.0E-2)
  !     PARAMETER (ALMIN1=0.10E-4, ALMIN2=0.15E-4, ALMAX=1.0E-1)
  !     PARAMETER (ALMIN1=0.00E-4, ALMIN2=0.40E-4, ALMAX=2.0E-2)
  !     PARAMETER (ALMIN1=0.20E-4, ALMIN2=0.40E-4, ALMAX=2.0E-2)
  !     PARAMETER (ALMIN1=0.25E-4, ALMIN2=0.50E-4, ALMAX=2.0E-2)
  !     PARAMETER (ALMIN1=0.40E-4, ALMIN2=0.50E-4, ALMAX=2.0E-2)
  !
  REAL(kind=r8), PARAMETER :: BLDMAX = 200.0_r8 !mb
  !
  INTEGER KBLMX

  !     PARAMETER (QI0=1.0E-4, QW0=1.0E-4)
  !     PARAMETER (QI0=0.0E-5, QW0=0.0E-0)
  REAL(kind=r8), PARAMETER :: QI0=1.0E-5_r8, QW0=1.0E-5_r8
  !     PARAMETER (QI0=1.0E-4, QW0=1.0E-5) ! 20050509
  !     PARAMETER (QI0=1.0E-5, QW0=1.0E-6)
  !     PARAMETER (QI0=0.0E-5, QW0=0.0E-5)
!!!   PARAMETER (QI0=5.0E-4, QW0=1.0E-5)
  !     PARAMETER (QI0=5.0E-4, QW0=5.0E-4)
  !     PARAMETER (QI0=2.0E-4, QW0=2.0E-5)
  !     PARAMETER (QI0=2.0E-5, QW0=2.0E-5)
  !     PARAMETER (QI0=2.0E-4, QW0=1.0E-4)
  !     PARAMETER (QI0=2.0E-4, QW0=1.0E-5)
  !     PARAMETER (QI0=1.0E-3, QW0=2.0E-5)
  !     PARAMETER (QI0=1.0E-3, QW0=7.0E-4)

  !     PARAMETER (C0I=5.0E-4)
  !     PARAMETER (C0I=4.0E-4)
  REAL(kind=r8), PARAMETER :: C0I=1.0E-3_r8

  !     parameter (c0=1.0e-3)
  !     parameter (c0=1.5e-3)
  REAL(kind=r8), PARAMETER ::  c0=2.0e-3_r8
  !     parameter (c0=1.0e-3, KBLMX=10, ERRMIN=0.0001, ERRMI2=0.1*ERRMIN)
  !     parameter (c0=2.0e-3, KBLMX=10, ERRMIN=0.0001, ERRMI2=0.1*ERRMIN)
  !
  !     parameter (TF=130.16, TCR=160.16, TCRF=1.0/(TCR-TF),TCL=2.0)
  !     parameter (TF=230.16, TCR=260.16, TCRF=1.0/(TCR-TF))
  REAL(kind=r8), PARAMETER ::  TF=233.16_r8, TCR=263.16_r8, TCRF=1.0_r8/(TCR-TF),TCL=2.0_r8
  !
  !     For Tilting Angle Specification
  !
  REAL(kind=r8) :: TLBPL(7) 
  REAL(kind=r8) :: drdp(5)
  REAL(kind=r8) :: VTP
  !
  REAL(KIND=r8), PARAMETER :: PLAC(1:8)=(/100.0_r8, 200.0_r8, 300.0_r8, 400.0_r8, 500.0_r8, 600.0_r8, 700.0_r8, 800.0_r8/)
  !     DATA TLAC/ 37.0,  25.0,  17.0,  12.0,  10.0,  8.0,  6.0,  5.0/
  !     DATA TLAC/ 35.0,  24.0,  17.0,  12.0,  10.0,  8.0,  6.0,  5.0/
  !     DATA TLAC/ 35.0,  25.0,  20.0,  17.5,  15.0,  12.5,  10.0,  5.0/
  REAL(KIND=r8), PARAMETER :: TLAC(1:8)=(/ 35.0_r8,  25.0_r8,  20.0_r8,  17.5_r8,  15.0_r8,  12.5_r8,  10.0_r8,  7.5_r8/)
  !     DATA TLAC/ 37.0,  26.0,  18.0,  14.0,  10.0,  8.0,  6.0,  5.0/
  !     DATA TLAC/ 25.0,  22.5,  20.0,  17.5,  15.0,  12.5,  10.0,  10.0/
  REAL(KIND=r8), PARAMETER :: REFP(1:6)=(/500.0_r8, 300.0_r8, 250.0_r8, 200.0_r8, 150.0_r8, 100.0_r8/)
  !     DATA REFR/ 0.25,   0.5,  0.75,   1.0,   1.5,   2.0/
  !     DATA REFR/ 0.5,   1.0,  1.5,   2.0,   3.0,   4.0/
  REAL(KIND=r8), PARAMETER :: REFR(1:6)=(/ 1.0_r8,   2.0_r8,  3.0_r8,   4.0_r8,   6.0_r8,   8.0_r8/)
  !
  REAL(kind=r8) :: AC(16)
  REAL(kind=r8) :: AD(16)
  !
  INTEGER, PARAMETER :: NQRP=500001

  REAL(kind=r8) C1XQRP, C2XQRP, TBQRP(NQRP), TBQRA(NQRP)     &
       &,                    TBQRB(NQRP)
  !
  REAL(kind=r8) :: rasalf

  INTEGER, PARAMETER :: NVTP=10001
  REAL(kind=r8)      :: C1XVTP, C2XVTP, TBVTP(NVTP)
  ! 
  !  END    module module_ras
  !
  !
  !      module module_rascnv
  !
  REAL(KIND=r8), PARAMETER :: frac=0.5_r8, crtmsf=0.0_r8
  !     PARAMETER (MAX_NEG_BOUY=0.25, REVAP=.TRUE., CUMFRC=.true.)
  !     PARAMETER (MAX_NEG_BOUY=0.20, REVAP=.TRUE., CUMFRC=.true.)
  !     PARAMETER (MAX_NEG_BOUY=0.15, REVAP=.true., CUMFRC=.true.)
  !LL3  PARAMETER (MAX_NEG_BOUY=0.10, REVAP=.TRUE., CUMFRC=.true.)
  !     PARAMETER (MAX_NEG_BOUY=0.15, REVAP=.true., CUMFRC=.false.)
  REAL(KIND=r8), PARAMETER :: MAX_NEG_BOUY=0.15_r8
  LOGICAL      , PARAMETER :: REVAP=.TRUE., CUMFRC=.TRUE.
  !     PARAMETER (MAX_NEG_BOUY=0.15, REVAP=.true., CUMFRC=.false.)
  !     PARAMETER (MAX_NEG_BOUY=0.05, REVAP=.true., CUMFRC=.true.)


  LOGICAL, PARAMETER :: WRKFUN = .FALSE.,  UPDRET = .FALSE.
  LOGICAL, PARAMETER :: CRTFUN = .TRUE.,   CALKBL = .FALSE., BOTOP=.TRUE.
  ! LOGICAL, PARAMETER :: CRTFUN = .TRUE.,    BOTOP=.TRUE.
  !
  !     parameter (rhfacs=0.70, rhfacl=0.70)
  !     parameter (rhfacs=0.75, rhfacl=0.75)
  REAL(KIND=r8), PARAMETER ::  rhfacs=0.80_r8, rhfacl=0.80_r8
  !     parameter (rhfacs=0.80, rhfacl=0.85)
  REAL(KIND=r8), PARAMETER :: FACE=5.0_r8, DELX=10000.0_r8, DDFAC=FACE*DELX*0.001_r8
  !
  !     real (kind=r8), parameter :: pgftop=0.7, pgfbot=0.3        &
  !     real (kind=r8), parameter :: pgftop=0.75, pgfbot=0.35      &
  REAL (kind=r8), PARAMETER :: pgftop=0.80_r8, pgfbot=0.30_r8      
  REAL (kind=r8), PARAMETER :: pgfgrad=(pgfbot-pgftop)*0.001_r8
  !
  !
  INTEGER       , ALLOCATABLE :: nlons(:)
  REAL (kind=r8), ALLOCATABLE ::si(:)
  REAL (kind=r8), ALLOCATABLE ::sl(:)
  INTEGER      :: latg

  !      end module module_rascnv
  !
  !   module funcphys
  INTEGER,PARAMETER:: nxpvs=7501

  REAL(r8) c1xpvs,c2xpvs,tbpvs(nxpvs)

  !   END module funcphys
!  Parameters
        INTEGER,PARAMETER:: n=624
        INTEGER,PARAMETER:: m=397
        INTEGER(KIND=i8),PARAMETER:: mata=-1727483681_i8 ! constant vector a
        INTEGER(KIND=i8),PARAMETER:: umask=-2147483648_i8 ! most significant w-r bits
        INTEGER(KIND=i8),PARAMETER:: lmask =2147483647_i8 ! least significant r bits
        INTEGER(KIND=i8),PARAMETER:: tmaskb=-1658038656_i8 ! tempering parameter
        INTEGER(KIND=i8),PARAMETER:: tmaskc=-272236544_i8 ! tempering parameter
        INTEGER(KIND=i8),PARAMETER:: mag01(0:1)=(/0_i8,mata/)
        INTEGER(KIND=i8),PARAMETER:: iseed=4357_i8

!  Defined types
        TYPE random_stat
          PRIVATE
          INTEGER:: mti=n+1
          INTEGER:: mt(0:n-1)
          INTEGER:: iset
          REAL   :: gset
        END TYPE
!  Saved data
        TYPE(random_stat),SAVE:: sstat
!  Overloaded interfaces
        INTERFACE random_setseed
          MODULE PROCEDURE random_setseed_s
          MODULE PROCEDURE random_setseed_t
        END INTERFACE
        INTERFACE random_number
          MODULE PROCEDURE random_number_i
          MODULE PROCEDURE random_number_s
          MODULE PROCEDURE random_number_t
        END INTERFACE

  PUBLIC :: Init_Cu_RAS3PHASE
  PUBLIC :: Run_Cu_RAS3PHASE
  PUBLIC :: Finalize_Cu_RAS3PHASE
  
  
CONTAINS


!  Subprogram random_setseed_s
!  Sets seed in saved mode.
        SUBROUTINE random_setseed_s(inseed)
          IMPLICIT NONE
          INTEGER,INTENT(in):: inseed
          CALL random_setseed_t(inseed,sstat)
        END SUBROUTINE random_setseed_s
!  Subprogram random_setseed_t
!  Sets seed in thread-safe mode.
        SUBROUTINE random_setseed_t(inseed,stat)
          IMPLICIT NONE
          INTEGER,INTENT(in):: inseed
          TYPE(random_stat),INTENT(out):: stat
          INTEGER ii,mti
          ii=inseed
          IF(ii.EQ.0) ii=iseed
          stat%mti=n
          stat%mt(0)=IAND(ii,-1)
          DO mti=1,n-1
            stat%mt(mti)=IAND(69069*stat%mt(mti-1),-1)
          ENDDO
          stat%iset=0
          stat%gset=0.
        END SUBROUTINE random_setseed_t
!  Subprogram random_number_i
!  Generates random numbers in interactive mode.
        SUBROUTINE random_number_i(harvest,inseed)
          IMPLICIT NONE
          REAL,INTENT(out):: harvest(:)
          INTEGER,INTENT(in):: inseed
          TYPE(random_stat) stat
          CALL random_setseed_t(inseed,stat)
          CALL random_number_t(harvest,stat)
        END SUBROUTINE random_number_i
!  Generates random numbers in saved mode; overloads Fortran 90 standard.
        SUBROUTINE random_number_s(harvest)
          IMPLICIT NONE
          REAL,INTENT(out):: harvest(:)
          IF(sstat%mti.EQ.n+1) CALL random_setseed_t(int(iseed,kind=i4),sstat)
          CALL random_number_t(harvest,sstat)
        END SUBROUTINE random_number_s
!  Subprogram random_number_t
!  Generates random numbers in thread-safe mode.
        SUBROUTINE random_number_t(harvest,stat)
          IMPLICIT NONE
          REAL, PARAMETER  :: twop32=2.0**32
          REAL, PARAMETER  :: twop32m1i=1.0/(twop32-1.0)
          REAL,INTENT(out):: harvest(:)
          TYPE(random_stat),INTENT(inout):: stat
          INTEGER j,kk,y
          INTEGER tshftu,tshfts,tshftt,tshftl
          tshftu(y)=ISHFT(y,-11)
          tshfts(y)=ISHFT(y,7)
          tshftt(y)=ISHFT(y,15)
          tshftl(y)=ISHFT(y,-18)
          DO j=1,SIZE(harvest)
             IF(stat%mti.GE.n) THEN
                DO kk=0,n-m-1
                   y=INT(IOR(IAND(stat%mt(kk),umask),IAND(stat%mt(kk+1),lmask)),KIND=i4)
                   stat%mt(kk)=INT(IEOR(IEOR(stat%mt(kk+m),ISHFT(y,-1)), &
                        &                           mag01(IAND(y,1))),kind=i4)
                ENDDO
                DO kk=n-m,n-2
                   y=INT(IOR(IAND(stat%mt(kk),umask),IAND(stat%mt(kk+1),lmask)),KIND=i4)
                   stat%mt(kk)=INT(IEOR(IEOR(stat%mt(kk+(m-n)),ISHFT(y,-1)), &
                        &                           mag01(IAND(y,1))),KIND=i4)
                ENDDO
                y=INT(IOR(IAND(stat%mt(n-1),umask),IAND(stat%mt(0),lmask)),KIND=i4)
                stat%mt(n-1)=INT(IEOR(IEOR(stat%mt(m-1),ISHFT(y,-1)), &
                     &                      mag01(IAND(y,1))),KIND=i4)
                stat%mti=0
             ENDIF
             y=stat%mt(stat%mti)
             y=IEOR(y,tshftu(y))
             y=INT(IEOR(y,IAND(tshfts(y),tmaskb)),KIND=i4)
             y=INT(IEOR(y,IAND(tshftt(y),tmaskc)),KIND=i4)
             y=IEOR(y,tshftl(y))
             IF(y.LT.0) THEN
                harvest(j)=(REAL(y)+twop32)*twop32m1i
             ELSE
                harvest(j)=REAL(y)*twop32m1i
             ENDIF
             stat%mti=stat%mti+1
          ENDDO
        END SUBROUTINE random_number_t

  !-----------------------------------------------------------------------------------------
  SUBROUTINE Init_Cu_RAS3PHASE(kMax,dt,fhour,idate,iMax,jMax,ibMax,jbMax,si_in,sl_in)
    IMPLICIT NONE
    INTEGER      , INTENT(IN   ) :: kMax
    REAL(kind=r8), INTENT(IN   ) :: dt
    REAL(kind=r8), INTENT(IN   ) ::  fhour
    INTEGER      , INTENT(IN   ) :: idate(4)!=(/00,01,29,2015/)
    INTEGER      , INTENT(IN   ) :: iMax,jMax,ibMax,jbMax
    REAL(kind=r8), INTENT(IN   ) ::  si_in(kMax+1)
    REAL(kind=r8), INTENT(IN   ) ::  sl_in(kMax)
    !real (kind=r8) :: rannum(lonr*latr*nrcm)

    INTEGER       :: nrc,i,j,ij,ib,jb
    INTEGER       :: iseed
    INTEGER       :: seed0
    REAL(kind=r8),ALLOCATABLE ::  brf(:,:,:)
    REAL(kind=r8) ::  wrk(1)
    REAL(kind=r8), PARAMETER :: cons_0=0.0_r8,   cons_24=24.0_r8
    REAL(kind=r8), PARAMETER :: cons_99=99.0_r8, cons_1p0d9=1.0E9_r8
    REAL (kind=r8), ALLOCATABLE :: rannum(:)
    PRINT*,'Init_Cu_RAS3PHASE'
    ALLOCATE (si(kMax+1));si=si_in
    ALLOCATE (sl(kMax));sl=sl_in

    nrcm = MAX(nrcmax, NINT((nrcmax*dt)/600.0_r8))
    ALLOCATE(rannum(iMax*jMax*nrcm));rannum=0.0_r8
    ALLOCATE(brf(iMax,jMax,nrcm));brf=0.0_r8

    wrk=0.0_r8

    ALLOCATE(xkt2  (ibMax,jbMax,nrcm));xkt2=0.0_r8
    ALLOCATE(nlons(ibMax));nlons=iMax
    !
    latg=jMax
    seed0 = idate(1) + idate(2) + idate(3) + idate(4)
    CALL random_setseed(seed0)
    CALL RANDOM_NUMBER(wrk)
    seed0 = seed0 + NINT(wrk(1)*1000.0_r8)

    !
    iseed = INT( MOD(100.0_r8*SQRT(fhour*3600.0_r8),cons_1p0d9) + 1 + seed0 ,KIND=i4)

    CALL random_setseed(iseed)
    CALL RANDOM_NUMBER(rannum)
    ij=0
    do nrc=1,nrcm
       do j=1,jmax
          do i=1,iMax
             ij=ij+1
             brf(i,j,nrc) = rannum(ij)
          enddo
       enddo
    enddo

    DO nrc=1,nrcm
       IF (reducedGrid) THEN
          CALL LinearIJtoIBJB(brf(:,:,nrc),xkt2(:,:,nrc))
       ELSE
          CALL IJtoIBJB(brf(:,:,nrc) ,xkt2(:,:,nrc) )
       END IF
    END DO
    
    DEALLOCATE(rannum)
    DEALLOCATE(brf)
    CALL gfuncphys()
    CALL ras_init(kMax)

  END SUBROUTINE Init_Cu_RAS3PHASE
  !
  SUBROUTINE set_ras_afc(dt)
    IMPLICIT NONE
    REAL(kind=r8) :: DT
    !     AFC = -(1.04E-4*DT)*(3600./DT)**0.578
    AFC = -(1.01097E-4_r8*DT)*(3600.0_r8/DT)**0.57777778_r8
  END SUBROUTINE set_ras_afc

  SUBROUTINE ras_init(levs)
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN   ) :: levs

    !
    REAL(kind=r8) :: tem
    REAL(kind=r8) :: actop
    REAL(kind=r8) :: tem1
    REAL(kind=r8) :: tem2
    REAL(KIND=r8) :: A(15)
    INTEGER       :: i, l
    !     PARAMETER (ACTP=1.7,   FACM=1.20)
    REAL(kind=r8), PARAMETER :: ACTP=1.7_r8
    REAL(kind=r8), PARAMETER :: FACM=1.00_r8
    !     PARAMETER (ACTP=1.7,   FACM=0.90)
    !     PARAMETER (ACTP=1.7,   FACM=0.75)
    !     PARAMETER (ACTP=1.7,   FACM=0.60)
    !     PARAMETER (ACTP=1.7,   FACM=0.5)   ! cnt
    !     PARAMETER (ACTP=1.7,   FACM=0.4)
    !     PARAMETER (ACTP=1.7,   FACM=0.0)
    !
    !
    REAL(kind=r8), PARAMETER :: PH(1:15)=(/150.0_r8, 200.0_r8, 250.0_r8, 300.0_r8, 350.0_r8, 400.0_r8, 450.0_r8, 500.0_r8,    &
         550.0_r8, 600.0_r8, 650.0_r8, 700.0_r8, 750.0_r8, 800.0_r8, 850.0_r8/)
    !
    REAL(kind=r8), PARAMETER ::  AUX(1:15)=(/ 1.6851_r8, 1.1686_r8, 0.7663_r8, 0.5255_r8, 0.4100_r8, 0.3677_r8,           &
         &       0.3151_r8, 0.2216_r8, 0.1521_r8, 0.1082_r8, 0.0750_r8, 0.0664_r8,            &
         &       0.0553_r8, 0.0445_r8, 0.0633_r8/)
    !
    !LOGICAL :: first=.TRUE.
    !
    !IF (first) THEN
    !
    ALLOCATE (rasal(levs))
    !                                   set critical workfunction arrays
    ACTOP = ACTP*FACM
    DO L=1,15
       A(L) = AUX(L)*FACM
    ENDDO
    DO L=2,15
       TEM   = 1.0_r8 / (PH(L) - PH(L-1))
       AC(L) = (PH(L)*A(L-1) - PH(L-1)*A(L)) * TEM
       AD(L) = (A(L) - A(L-1)) * TEM
    ENDDO
    AC(1)  = ACTOP
    AC(16) = A(15)
    AD(1)  = 0.0_r8
    AD(16) = 0.0_r8
    !
    !       CALL SETES
    CALL SETQRP()
    CALL SETVTP()
    !
    !       kblmx = levs / 2
    !
    !       RASALF  = 0.10
    !       RASALF  = 0.20
    RASALF  = 0.30_r8
    !       RASALF  = 0.35
    !
    DO L=1,LEVS
       RASAL(L) = RASALF
    ENDDO
    !
    !
    DO i=1,7
       tlbpl(i) = (tlac(i)-tlac(i+1)) / (plac(i)-plac(i+1))
    ENDDO
    DO i=1,5
       drdp(i)  = (REFR(i+1)-REFR(i)) / (REFP(i+1)-REFP(i))
    ENDDO
    !
    VTP    = 36.34_r8*SQRT(1.2_r8)* (0.001_r8)**0.1364_r8
    !
    !IF (me .EQ. 0) PRINT *,' NO DOWNDRAFT FOR CLOUD TYPES'          &
    !     &,        ' DETRAINING WITHIN THE BOTTOM ',DD_DP,' hPa LAYERS'
    !
    !first = .FALSE.
    !ENDIF
    !
  END SUBROUTINE ras_init
  !-------------------------------------------------------------------------------
  SUBROUTINE gfuncphys()
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: gfuncphys    Compute all physics function tables
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Compute all physics function tables.  Lookup tables are
    !   set up for computing saturation vapor pressure, dewpoint temperature,
    !   equivalent potential temperature, moist adiabatic temperature and humidity,
    !   pressure to the kappa, and lifting condensation level temperature.
    !
    ! Program History Log:
    ! 1999-03-01  Iredell             f90 module
    !
    ! Usage:  call gfuncphys
    !
    ! Subprograms called:
    !   gpvs        compute saturation vapor pressure table
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    CALL gpvs()
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE gfuncphys
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------
  SUBROUTINE gpvs()
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: gpvs         Compute saturation vapor pressure table
    !   Author: N Phillips            W/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Computes saturation vapor pressure table as a function of
    !   temperature for the table lookup function fpvs.
    !   Exact saturation vapor pressures are calculated in subprogram fpvsx.
    !   The current implementation computes a table with a length
    !   of 7501 for temperatures ranging from 180. to 330. Kelvin.
    !
    ! Program History Log:
    !   91-05-07  Iredell
    !   94-12-30  Iredell             expand table
    ! 1999-03-01  Iredell             f90 module
    ! 2001-02-26  Iredell             ice phase
    !
    ! Usage:  call gpvs
    !
    ! Subprograms called:
    !   (fpvsx)    inlinable function to compute saturation vapor pressure
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    INTEGER jx
    REAL(r8) :: xmin,xmax,xinc,x,t
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xmin=180.0_r8
    xmax=330.0_r8
    xinc=(xmax-xmin)/(nxpvs-1)
    !   c1xpvs=1.-xmin/xinc
    c2xpvs=1.0_r8/xinc
    c1xpvs=1.0_r8-xmin*c2xpvs
    DO jx=1,nxpvs
       x=xmin+(jx-1)*xinc
       t=x
       tbpvs(jx)=fpvsx(t)
    ENDDO
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE gpvs
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------
  ELEMENTAL FUNCTION fpvs(t)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: fpvs         Compute saturation vapor pressure
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Compute saturation vapor pressure from the temperature.
    !   A linear interpolation is done between values in a lookup table
    !   computed in gpvs. See documentation for fpvsx for details.
    !   Input values outside table range are reset to table extrema.
    !   The interpolation accuracy is almost 6 decimal places.
    !   On the Cray, fpvs is about 4 times faster than exact calculation.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             expand table
    ! 1999-03-01  Iredell             f90 module
    ! 2001-02-26  Iredell             ice phase
    !
    ! Usage:   pvs=fpvs(t)
    !
    !   Input argument list:
    !     t          Real(r8) temperature in Kelvin
    !
    !   Output argument list:
    !     fpvs       Real(r8) saturation vapor pressure in Pascals
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8) :: fpvs
    REAL(r8),INTENT(in):: t
    INTEGER jx
    REAL(r8) xj
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xj=MIN(MAX(c1xpvs+c2xpvs*t,1.0_r8),REAL(nxpvs,r8))
    jx=INT(MIN(xj,nxpvs-1.0_r8),KIND=i4)
    fpvs=tbpvs(jx)+(xj-jx)*(tbpvs(jx+1)-tbpvs(jx))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION fpvs
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------
  ELEMENTAL FUNCTION fpvsx(t)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: fpvsx        Compute saturation vapor pressure
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Exactly compute saturation vapor pressure from temperature.
    !   The saturation vapor pressure over either liquid and ice is computed
    !   over liquid for temperatures above the triple point,
    !   over ice for temperatures 20 degress below the triple point,
    !   and a linear combination of the two for temperatures in between.
    !   The water model assumes a perfect gas, constant specific heats
    !   for gas, liquid and ice, and neglects the volume of the condensate.
    !   The model does account for the variation of the latent heat
    !   of condensation and sublimation with temperature.
    !   The Clausius-Clapeyron equation is integrated from the triple point
    !   to get the formula
    !       pvsl=con_psat*(tr**xa)*exp(xb*(1.-tr))
    !   where tr is ttp/t and other values are physical constants.
    !   The reference for this computation is Emanuel(1994), pages 116-117.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             exact computation
    ! 1999-03-01  Iredell             f90 module
    ! 2001-02-26  Iredell             ice phase
    !
    ! Usage:   pvs=fpvsx(t)
    !
    !   Input argument list:
    !     t          Real(r8) temperature in Kelvin
    !
    !   Output argument list:
    !     fpvsx      Real(r8) saturation vapor pressure in Pascals
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8) :: fpvsx
    REAL(r8),INTENT(in) :: t
    REAL(r8),PARAMETER :: tliq=con_ttp
    REAL(r8),PARAMETER :: tice=con_ttp-20.0
    REAL(r8),PARAMETER :: dldtl=con_cvap-con_cliq
    REAL(r8),PARAMETER :: heatl=con_hvap
    REAL(r8),PARAMETER :: xponal=-dldtl/con_rv
    REAL(r8),PARAMETER :: xponbl=-dldtl/con_rv+heatl/(con_rv*con_ttp)
    REAL(r8),PARAMETER :: dldti=con_cvap-con_csol
    REAL(r8),PARAMETER :: heati=con_hvap+con_hfus
    REAL(r8),PARAMETER :: xponai=-dldti/con_rv
    REAL(r8),PARAMETER :: xponbi=-dldti/con_rv+heati/(con_rv*con_ttp)
    REAL(r8) tr,w,pvl,pvi
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tr=con_ttp/t
    IF(t.GE.tliq) THEN
       fpvsx=con_psat*(tr**xponal)*EXP(xponbl*(1.0_r8-tr))
    ELSEIF(t.LT.tice) THEN
       fpvsx=con_psat*(tr**xponai)*EXP(xponbi*(1.0_r8-tr))
    ELSE
       w=(t-tice)/(tliq-tice)
       pvl=con_psat*(tr**xponal)*EXP(xponbl*(1.0_r8-tr))
       pvi=con_psat*(tr**xponai)*EXP(xponbi*(1.0_r8-tr))
       fpvsx=w*pvl+(1.0_r8-w)*pvi
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION  fpvsx

  !-------------------------------------------------------------------------------

      SUBROUTINE sig2press(nCols    ,& !
                           kMax     ,&!
                           pgr      ,&!
                           sl       ,&!
                           si       ,&!
                           prsi     ,&!
                           prsl     ,&!
                           prsik    ,&!
                           prslk      )
 
      IMPLICIT NONE
 
      INTEGER      , INTENT(IN   ) ::  nCols
      INTEGER      , INTENT(IN   ) ::  kMax
      REAL(kind=r8), INTENT(IN   ) ::  pgr(nCols)    !cb
      REAL(kind=r8), INTENT(IN   ) ::  sl(kMax)
      REAL(kind=r8), INTENT(IN   ) ::  si(kMax+1)
      REAL(kind=r8), INTENT(OUT  ) ::  prsi(nCols,kMax+1)
      REAL(kind=r8), INTENT(OUT  ) ::  prsl(nCols,kMax)
      REAL(kind=r8), INTENT(OUT  ) ::  prsik(nCols,kMax+1)
      REAL(kind=r8), INTENT(OUT  ) ::  prslk(nCols,kMax)
      
      REAL(kind=r8) :: pgrk(nCols)
      REAL(kind=r8) :: tem, rkapi, rkapp1
      REAL(kind=r8) :: slk(kMax) 
      REAL(kind=r8) :: sik(kMax+1)
      REAL(kind=r8) :: sikp1(kMax+1)
      INTEGER       :: i,k

      REAL(kind=r8), PARAMETER :: PT01=0.01
      REAL(kind=r8), PARAMETER:: con_rd     =2.8705e+2      ! gas constant air    (J/kg/K)
      REAL(kind=r8), PARAMETER:: con_cp     =1.0046e+3      ! spec heat air @p    (J/kg/K)
      REAL(kind=r8), PARAMETER:: con_rocp   =con_rd/con_cp
      REAL(kind=r8), PARAMETER:: rkap = con_rocp
      REAL(kind=r8), PARAMETER:: rk = con_rocp
      REAL(kind=r8), PARAMETER :: cb2mb   = 10.0

      RKAPI  = 1.0 / RKAP
      RKAPP1 = 1.0 + RKAP
      DO k=1,kMax+1
        sik(k)   = si(k) ** rkap
        sikp1(k) = si(k) ** rkapp1
      END DO
      DO k=1,kMax
        tem      = rkapp1 * (si(k) - si(k+1))
        slk(k)   = (sikp1(k)-sikp1(k+1))/tem
        !sl(k)    = slk(k) ** rkapi
      END DO

      DO i=1,nCols
         prsi(i,kMax+1)  = si(kMax+1)*pgr(i)      ! prsi are now pressures
         pgrk(i)         = (pgr(i)*pt01) ** rk
         prsik(i,kMax+1) = sik(kMax+1) * pgrk(i)
      END DO
      DO k=1,kMax
        DO i=1,nCols
          prsi(i,k)  = si(k)*pgr(i)               ! prsi are now pressures
          prsl(i,k)  = sl(k)*pgr(i)
          prsik(i,k) = sik(k) * pgrk(i)
          prslk(i,k) = slk(k) * pgrk(i)
        END DO
      END DO
      RETURN
      END SUBROUTINE sig2press

  !-----------------------------------------------------------------------------------------
  SUBROUTINE Run_Cu_RAS3PHASE(microphys,nClass,nCols,latco,ntrac, kMax,dt,&
                              t2,t3,q2,q3,ql2,ql3,qi2,qi3,u2,v2,ps2,colrad,ustar,pblh,tke,mask2,&!,KPBL,ustar,&
                              kbot ,ktop ,kuo  ,raincv,gvarm,gvarp)
    IMPLICIT NONE
    LOGICAL      , INTENT(IN   ) :: microphys
    INTEGER      , INTENT(IN   ) :: nClass
    INTEGER      , INTENT(IN   ) :: nCols
    INTEGER      , INTENT(IN   ) :: latco
    INTEGER      , INTENT(IN   ) :: ntrac
    INTEGER      , INTENT(IN   ) :: kMax
    REAL(KIND=r8), INTENT(IN   ) :: dt
    REAL(KIND=r8), INTENT(IN   ) :: t2 (nCols,kMax)!     tgrs     - real, layer mean temperature ( k )      ix,levs !
    REAL(KIND=r8), INTENT(INOUT) :: t3 (nCols,kMax)!     tgrs     - real, layer mean temperature ( k )      ix,levs !
    REAL(KIND=r8), INTENT(IN   ) :: q2 (nCols,kMax)!     tgrs     - real, layer mean temperature ( k )      ix,levs !
    REAL(KIND=r8), INTENT(INOUT) :: q3 (nCols,kMax)!     tgrs     - real, layer mean temperature ( k )      ix,levs !
    REAL(KIND=r8), INTENT(IN   ) :: ql2 (nCols,kMax)!     tgrs     - real, layer mean temperature ( k )      ix,levs !
    REAL(KIND=r8), INTENT(INOUT) :: ql3 (nCols,kMax)!     tgrs     - real, layer mean temperature ( k )      ix,levs !
    REAL(KIND=r8), INTENT(IN   ) :: qi2 (nCols,kMax)!     tgrs     - real, layer mean temperature ( k )      ix,levs !
    REAL(KIND=r8), INTENT(INOUT) :: qi3 (nCols,kMax)!     tgrs     - real, layer mean temperature ( k )      ix,levs !
    REAL(KIND=r8), INTENT(IN   ) :: u2 (nCols,kMax)!     tgrs     - real, layer mean temperature ( k )      ix,levs !
    REAL(KIND=r8), INTENT(IN   ) :: v2 (nCols,kMax)!     tgrs     - real, layer mean temperature ( k )      ix,levs !
    REAL(KIND=r8), INTENT(IN   ) :: ps2(nCols)        !surface pressure cb
    REAL(KIND=r8), INTENT(IN   ) :: colrad (nCols)
    REAL(KIND=r8), INTENT(IN   ) :: ustar(nCols)
    REAL(KIND=r8), INTENT(IN   ) :: pblh(nCols)
    REAL(KIND=r8), INTENT(IN   ) :: tke   (1:nCols,1:kMax)
    INTEGER      , INTENT(IN   ) :: mask2(nCols)
    REAL(KIND=r8)   , OPTIONAL,   INTENT(INOUT) :: gvarm (nCols,kmax,nClass)
    REAL(KIND=r8)   , OPTIONAL,   INTENT(INOUT) :: gvarp (nCols,kmax,nClass)

    
!    INTEGER      , INTENT(IN   ) ::  KPBL(nCols)
!    REAL(KIND=r8), INTENT(IN   ) :: ustar(nCols)
!    INTEGER      , INTENT(OUT  ) ::  kbot (nCols)
!    INTEGER      , INTENT(OUT  ) ::  ktop (nCols)
!    INTEGER      , INTENT(OUT  ) ::  kuo  (nCols)
    REAL(kind=r8) , INTENT(OUT  ) :: raincv(nCols)
    INTEGER       :: KPBL(nCols)
    INTEGER       :: kbot (nCols)
    INTEGER       :: ktop (nCols)
    INTEGER       :: kuo  (nCols)
    REAL(kind=r8) :: rain1(nCols)

    REAL(KIND=r8) :: tgrs (nCols,kMax)!     tgrs     - real, layer mean temperature ( k )ix,levs !
    REAL(KIND=r8) :: qgrs (nCols,kMax,3)! qgrs     - real, layer mean tracer concentration     ix,levs,ntrac!
    REAL(KIND=r8) :: qliq (nCols,kMax)!     qgrs     - real, layer mean tracer concentrationix,levs,ntrac!
    REAL(KIND=r8) :: qice (nCols,kMax)!     qgrs     - real, layer mean tracer concentrationix,levs,ntrac!

    REAL(KIND=r8) :: ugrs (nCols,kMax)!     ugrs,vgrs- real, u/v component of layer wind ix,levs !
    REAL(KIND=r8) :: vgrs (nCols,kMax)
    REAL(KIND=r8) :: ps(nCols)        !cb

    REAL (kind=r8), PARAMETER ::                               &
!    & dxmax=-8.8818363,   dxmin=-5.2574954,   dxinv=1.0/(dxmax-dxmin)
!    & dxmax=-9.5468126,   dxmin=-5.2574954,   dxinv=1.0/(dxmax-dxmin)
!    & dxmax=-18.40047804, dxmin=-9.800790154, dxinv=1.0/(dxmax-dxmin)
!    & dxmax=-15.8949521,  dxmin=-9.800790154, dxinv=1.0/(dxmax-dxmin)
!    & dxmax=-15.31958795, dxmin=-9.800790154, dxinv=1.0/(dxmax-dxmin)
!    & dxmax=-14.95494484, dxmin=-9.800790154, dxinv=1.0/(dxmax-dxmin)
!    & dxmax=-14.50865774, dxmin=-9.800790154, dxinv=1.0/(dxmax-dxmin)
     & dxmax=-16.118095651,dxmin=-9.800790154, dxinv=1.0/(dxmax-dxmin)
      REAL (kind=r8), PARAMETER :: rhc_max = 0.9999      ! 20060512
!     real (kind=r8), parameter :: rhc_max = 0.999       ! for pry
      REAL (kind=r8), PARAMETER :: wg = 0.5_r8 !wg - gustiness factor m/s
      REAL (kind=r8), PARAMETER :: ccwf     = 1.0  !         print *,' RAS Convection scheme used with ccwf=',ccwf
      REAL(kind=r8), PARAMETER :: cb2mb   = 10.0


    INTEGER, PARAMETER ::   num_p3d=3
    !     lonf,latg- integer, number of lon/lat points                 1    !
    REAL(KIND=r8) :: prsi(nCols,kMax+1)!     prsi     - real, pressure at layer interfaces             ix,levs+1
    REAL(KIND=r8) :: prsl(nCols,kMax)!     prsl     - real, mean layer presure                       ix,levs !
    REAL(KIND=r8) :: prsik(nCols,kMax+1)
    REAL(KIND=r8) :: prslk(nCols,kMax)
    REAL(KIND=r8) :: phii(nCols,kMax+1)
    REAL(KIND=r8) :: phil(nCols,kMax)
    REAL(KIND=r8) :: del(nCols,kMax)
    !LOCAL 
    REAL(KIND=r8) :: coslat(nCols)!     coslat   - real, cos of latitude                             nCols   !
    REAL(KIND=r8) :: work1(nCols)
    REAL(KIND=r8) :: work2(nCols)
    REAL(KIND=r8) :: rhc(nCols,kMax)
    REAL(KIND=r8) :: clw(nCols,kMax,2+nClass)
    LOGICAL       :: CALKBLMsk(nCols)
    REAL(kind=r8) :: rannum(nCols,nrcm)
    INTEGER       :: i,k,nn
    INTEGER       :: itc
    LOGICAL       :: trans_trac       = .TRUE.    ! This is effective only for RAS
    INTEGER       :: tottracer=0
    INTEGER       :: ntcw=3
    INTEGER       :: ntoz=0
    INTEGER       :: ncld=1
    INTEGER       :: itrc
    INTEGER       :: trc_shft=0
    INTEGER       :: tracers
    REAL(kind=r8) :: crtrh(3)
    REAL(kind=r8) :: rhbbot 
    REAL(kind=r8) :: rhpbl  
    REAL(kind=r8) :: theta
    REAL(kind=r8) :: tem1
    REAL(kind=r8) :: tem2
    REAL(kind=r8) :: rhbtop 
    REAL(kind=r8) :: DDVEL(nCols)
    REAL(kind=r8) :: garea(nCols)
    REAL(kind=r8) :: ua(nCols)
    REAL(kind=r8) :: cd(nCols)
    REAL(kind=r8) :: SS
    INTEGER       :: lmh(nCols) 
    REAL(kind=r8) :: ccwfac(nCols)
    REAL(kind=r8) :: ud_mf(nCols,kMax)
    REAL(kind=r8) :: dd_mf(nCols,kMax)
    REAL(kind=r8) :: det_mf(nCols,kMax)
    LOGICAL       :: lprnt
    !IF(nClass>0 .and. PRESENT(gvarm))THEN
       tottracer=nClass
    !ELSE
    !   tottracer=ntrac
    !END IF
    !tottracer=ntrac
    IF ( ccwf >= 0.0 ) THEN
      ccwfac = ccwf
    ELSE
      ccwfac = -999.0
    ENDIF
    DO i=1, nCols
       ps(i)=ps2(i) !cb
       lmh(i) = kMax
       DDVEL(i)=0.0_r8
       rain1(i)=0.0_r8

    END DO
    DO k=1,kMax
       DO i=1, nCols
           tgrs (i,k)= t3(i,k) ! layer mean temperature ( k )	   K
           ugrs (i,k)= u2 (i,k)
           vgrs (i,k)= v2 (i,k)
           qliq (i,k)=ql3(i,k)
           qice (i,k)=qi3(i,k)
       END DO      
    END DO
    DO itrc=1,3
       DO k=1,kMax
          DO i=1, nCols
             IF(itrc ==1)THEN
                qgrs (i,k,itrc)= q3(i,k)!qgrs     - real, layer mean tracer concentration     ix,levs,ntrac!
             END IF 

             IF(itrc ==2)THEN
                qgrs (i,k,itrc)= qi3(i,k)!qgrs     - real, layer mean tracer concentration     ix,levs,ntrac!
             END IF 
             IF(itrc ==3)THEN
                qgrs (i,k,itrc)= ql3(i,k)!qgrs     - real, layer mean tracer concentration     ix,levs,ntrac!
             END IF 
          END DO      
       END DO
    END DO

!     allocate ( clw(nCols,kMax,tottracer+2) )
      CALL sig2press(nCols    ,& !
                           kMax     ,&!
                           ps       ,&! cb
                           sl       ,&!
                           si       ,&!
                           prsi     ,&!!cb
                           prsl     ,&!!cb
                           prsik    ,&!!cb
                           prslk      )!cb
      phii= 0.0_r8
      phil= 0.0_r8
      CALL GET_PRS(nCols     ,&!
                     nCols   ,&!
                     kMax    ,&!
                     1       ,&!ntrac
                     tgrs  (1:nCols,1:kMax)     ,&!
                     qgrs  (1:nCols,1:kMax,1),&!
                     prsi  (1:nCols,1:kMax+1)  ,&!
                     prsik (1:nCols,1:kMax+1)  ,&!
                     prsl  (1:nCols,1:kMax)    ,&!
                     prslk (1:nCols,1:kMax)    ,&!
                     phii  (1:nCols,1:kMax+1)  ,&!
                     phil  (1:nCols,1:kMax)    ,&!
                     del   (1:nCols,1:kMax)      )

!       PRINT*,MAXVAL(phii),MINVAL(phii),MAXVAL(phil),MINVAL(phil)
!DO k = 1, kMax
!      DO i = 1, nCols
!         IF(phii(i,k)/con_g < pblh(i))THEN
!            kpbl(i) =k
!         END IF 
!      END DO
!    END DO

!-------------------------------------------
! calculate PBL TOP level
!-----------------------------
       ! grell mask
       !      mask2(i)=0 ! land
       !      mask2(i)=1 ! water/ocean

    DO i=1,nCols
       IF(mask2(i) == 0)THEN
          !      mask2(i)=0 ! land
          CALKBLMsk(i)=.TRUE.
       ELSE
          CALKBLMsk(i)=.FALSE.
       END IF
    END DO


    kpbl=1
    DO k=1,kMax/2
       DO i = 1, nCols
         ! IF(tke(i,k) >= 0.17_r8 .and. phii(i,k)/con_g <  pblh(i))THEN   !!
         IF(phii(i,k)/con_g <  pblh(i))THEN   !!

            kpbl(i)=k
          END IF
       END DO
    END DO
!---------------------------------------

    lprnt=.FALSE.
    crtrh(:)         = 0.95_r8
    DO itc=1,nrcm
       DO i=1, nCols
          rannum(i,itc) =xkt2  (i,latco,itc)
       END DO
    END DO
      DO i = 1, nCols
         ua(i)=SQRT((ugrs(i,1)**2)+(vgrs (i,1)**2))
         SS=(ua(i) + wg)           !velocity incl. gustiness param.

         CD(i)=(ustar(i)/SS)**2

       ! colrad.....colatitude  colrad=0-3.14 from np to sp in radians
       ! colrad.....colatitude  colrad=0-pi   from np to sp in radians

       !theta = 90.0_r8   -(180.0_r8/pai)*colrad(i) ! colatitude -> latitude
       theta = (con_pi/2)-colrad(i) ! colatitude -> latitude
       coslat(i)   = cos(((colrad(i)))-(3.1415926e0_r8/2.0_r8))

       !coslat(i) = COS(theta)
       ! the 180 degrees are divided into 37 bands with 5deg each
       ! except for the first and last, which have 2.5 deg
       ! The centers of the bands are located at:
       !   90, 85, 80, ..., 5, 0, -5, ..., -85, -90 (37 latitudes)
 
        tem1       = con_rerth * (con_pi + con_pi)*coslat(i) / nlons(i)
        tem2       = con_rerth * con_pi/latg
        garea(i)   = tem1 * tem2
        !dlength(i) = sqrt( tem1*tem1 + tem2*tem2 )
      ENDDO
!
!  --- ...  figure out how many extra tracers are there
!
!      IF ( trans_trac ) THEN
!        IF ( ntcw > 0 ) THEN
!          IF ( ntoz < ntcw ) THEN
!            trc_shft = ntcw + ncld
!          ELSE
!            trc_shft = ntoz
!          ENDIF
!        ELSEIF ( ntoz > 0 ) THEN
!          trc_shft = ntoz
!        ELSE
!          trc_shft = 1
!        ENDIF
!
!        tracers   = ntrac - trc_shft
!        IF ( ntoz > 0 ) tottracer = tracers + 1  ! ozone is added separately
!      ELSE
!        tottracer = 0                            ! no convective transport of tracers
!      ENDIF

!     allocate ( clw(nCols,kMax,tottracer+2) )
      DO k = 1, kMax
        DO i = 1, nCols
          clw(i,k,1) = 0.0
          clw(i,k,2) = -999.9
        ENDDO
      ENDDO
!  --- ...  for convective tracer transport (while using ras)

!      if ( ras ) then
!        IF ( tottracer > 0 ) THEN
!
!          IF ( ntoz > 0 ) THEN
!            clw(:,:,3) = qgrs(:,:,ntoz)
!
!            IF ( tracers > 0 ) THEN
!              DO nn = 1, tracers
!                clw(:,:,3+nn) = qgrs(:,:,nn+trc_shft)
!              ENDDO
!            ENDIF
!          ELSE
!            DO nn = 1, tracers
!              clw(:,:,2+nn) = qgrs(:,:,nn+trc_shft)
!            ENDDO
!          ENDIF
!
!        ENDIF
!      endif   ! end if_ras
!  --- ...  calling precipitation processes

      DO i = 1, nCols
        work1(i) = (LOG(coslat(i) / (nlons(i)*latg)) - dxmin) * dxinv
        work1(i) = MAX( 0.0_r8, MIN( 1.0_r8, work1(i) ) )
        work2(i) = 1.0_r8 - work1(i)
      ENDDO

!  --- ...  calling convective parameterization

      IF ( ntcw > 0 ) THEN
      rhbbot = crtrh(1)
      rhpbl  = crtrh(2)
      rhbtop = crtrh(3)

        DO k = 1, kMax
          DO i = 1, nCols
            rhc(i,k) = rhbbot - (rhbbot - rhbtop) * (1.0_r8 - prslk(i,k))
            rhc(i,k) = rhc_max*work1(i) + rhc(i,k)*work2(i)
            rhc(i,k) = MAX( 0.0_r8, MIN( 1.0_r8, rhc(i,k) ) )
          ENDDO
        ENDDO

        IF ( num_p3d == 3 ) THEN    ! call brad ferrier's microphysics
!  --- ...  algorithm to separate different hydrometeor species

          DO k = 1, kMax
            DO i = 1, nCols
              clw(i,k,1)  = qice(i,k)
              clw(i,k,2)  = qliq(i,k)
!  --- ...  array to track fraction of "cloud" in the form of ice
            ENDDO
          ENDDO
          IF(nClass>0 .and. PRESENT(gvarm))THEN
             DO itrc=1,nClass
                DO k = 1, kMax
                   DO i = 1, nCols
                      clw(i,k,2+itrc)  = gvarp (i,k,itrc)
                   END DO
                END DO
            END DO
          END IF
        ELSE   ! if_num_p3d
        

          DO k = 1, kMax
            DO i = 1, nCols
              clw(i,k,1) = qgrs(i,k,1)
            ENDDO
          ENDDO

        ENDIF  ! end if_num_p3d

      ELSE    ! if_ntcw

        rhc(:,:) = 1.0

      ENDIF   ! end if_ntcw
    !    prepare input, erase output
    !
    DO i=1,nCols
       kuo (i)=0
       kbot(i)=1
       ktop(i)=1
    END DO


        call rascnv(   nCols,    kMax,   2*dt, 2*dt, rannum,         &
                 tgrs ,    qgrs(:,:,1:1),   ugrs,    vgrs, qice,qliq,clw,tottracer ,&
                 prsi ,   prsl,   prsik,  prslk, phil,  phii,&
                 kpbl ,   cd,     rain1,  kbot,  ktop,  kuo,&
                 DDVEL, flipv, cb2mb,garea, lmh, ccwfac, &
                 nrcm, rhc,CALKBLMsk,  ud_mf, dd_mf, det_mf, lprnt)


    DO i = 1, nCols
        raincv(i) =rain1(i)*0.5_r8
       IF(RAINCV(i) > 0.0_r8)kuo(i)=1
    ENDDO

    DO k=1,kMax
       DO i=1, nCols
            IF(rain1(i) > 0.0_r8)THEN
              t3(i,k) = tgrs (i,k  ) !t3 (i,k) + (tgrs (i,k  ) - t2  (i,k))/(2*dt)! layer mean temperature ( k )K
              q3 (i,k)= qgrs (i,k,1)!q3 (i,k) + (qgrs (i,k,1) - q2  (i,k))/(2*dt)
              ql3(i,k)= qliq (i,k)!ql3(i,k) + (clw  (i,k,2) - qliq(i,k))/(2*dt)
              qi3(i,k)= qice (i,k)!qi3(i,k) + (clw  (i,k,1) - qice(i,k))/(2*dt)
            END IF   
       END DO      
    END DO
          IF(nClass>0 .and. PRESENT(gvarm))THEN
             DO itrc=1,nClass
                DO k = 1, kMax
                   DO i = 1, nCols
                       !gvarp (i,k,itrc) =clw(i,k,2+itrc)
                   END DO
                END DO
            END DO
          END IF

!      deallocate ( clw )

  END SUBROUTINE Run_Cu_RAS3PHASE

  !
  !  !-----------------------------------------------------------------------------------------

  SUBROUTINE rascnv(   nCols,     k,      dt,    dtf,  rannum      &
       &,                 tin,   qin,    uin,    vin,  qice,qliq,ccin,  trac &
       &,                 prsi,  prsl,   prsik,  prslk, phil,  phii       &
       &,                 KPBL,  CDRAG,  RAINC,  kbot,  ktop,  kuo        &
       &,                 DDVEL, FLIPV,  facmb,    garea, lmh, ccwfac     &
       &,                 nrcm,  rhc,  CALKBLMsk,  ud_mf, dd_mf,  det_mf,lprnt)
    !
    !*********************************************************************
    !*********************************************************************
    !************         Relaxed Arakawa-Schubert      ******************
    !************             Parameterization          ******************
    !************          Plug Compatible Driver       ******************
    !************               23 May 2002             ******************
    !************                                       ******************
    !************               Developed By            ******************
    !************                                       ******************
    !************             Shrinivas Moorthi         ******************
    !************                                       ******************
    !************                  EMC/NCEP             ******************
    !*********************************************************************
    !*********************************************************************
    !
    !
    !      use module_ras, DPD => DD_DP
    !      use module_rascnv
    IMPLICIT NONE
    !
    !
    !      input
    !
    INTEGER      , INTENT(IN   ) ::  nrcm!     nrcm     - integer, number of random clouds                  1    !
    INTEGER      , INTENT(IN   ) ::  trac!     ntrac    - integer, number of tracers                        1    !
    INTEGER      , INTENT(IN   ) :: nCols!     nCols, IX   - integer, horiz dimention and num of used pts      1    !
    INTEGER      , INTENT(IN   ) :: k !     levs     - integer, vertical layer dimension                 1    !
    REAL(kind=r8), INTENT(IN   ) :: DT !     dtp,dtf  - real, time interval (second)                      1    !
    REAL(kind=r8), INTENT(IN   ) :: dtf!     dtp,dtf  - real, time interval (second)                      1    !
    REAL(kind=r8), INTENT(IN   ) :: rannum(nCols,nrcm)
    REAL(kind=r8), INTENT(INOUT) :: tin(nCols,k)
    REAL(kind=r8), INTENT(INOUT) :: qin(nCols,k)
    REAL(kind=r8), INTENT(INOUT) :: uin(nCols,k)
    REAL(kind=r8), INTENT(INOUT) :: vin(nCols,k)
    REAL(kind=r8), INTENT(INOUT) :: qice(nCols,k)
    REAL(kind=r8), INTENT(INOUT) :: qliq(nCols,k)
    REAL(kind=r8), INTENT(INOUT) :: ccin(nCols,k,trac+2)
    REAL(kind=r8), INTENT(IN   ) :: prsi(nCols,k+1) 
    REAL(kind=r8), INTENT(IN   ) :: prsik(nCols,k+1)
    REAL(kind=r8), INTENT(IN   ) :: prsl(nCols,k)
    REAL(kind=r8), INTENT(IN   ) :: prslk(nCols,k+1)
    REAL(kind=r8), INTENT(IN   ) :: phil(nCols,k)
    REAL(kind=r8), INTENT(IN   ) :: phii(nCols,k+1)
    INTEGER      , INTENT(IN   ) :: KPBL(nCols)
    REAL(kind=r8), INTENT(IN   ) :: CDRAG(nCols)

    REAL(kind=r8), INTENT(OUT  ) :: RAINC(nCols) 
    INTEGER      , INTENT(OUT  ) ::  kbot(nCols)
    INTEGER      , INTENT(OUT  ) ::  ktop(nCols)
    INTEGER      , INTENT(OUT  ) ::  kuo(nCols)
    REAL(kind=r8), INTENT(INOUT) :: DDVEL(nCols)
    LOGICAL      , INTENT(IN   ) :: FLIPV
    REAL(kind=r8), INTENT(IN   ) :: facmb
    REAL(kind=r8), INTENT(IN   ) :: garea(nCols)
    INTEGER      , INTENT(IN   ) :: lmh(nCols)
    REAL(kind=r8), INTENT(IN   ) :: ccwfac(nCols)  
    REAL(kind=r8), INTENT(IN   ) :: rhc(nCols,k)
    REAL(kind=r8), INTENT(OUT  ) :: ud_mf(nCols,k)
    REAL(kind=r8), INTENT(OUT  ) :: dd_mf(nCols,k)
    REAL(kind=r8), INTENT(OUT  ) :: det_mf(nCols,k)
    LOGICAL      , INTENT(IN   ) :: lprnt
    LOGICAL      , INTENT(IN   ) :: CALKBLMsk(nCols)






    !
    !     locals
    !
    INTEGER ::  ncrnd

    REAL(kind=r8) :: RAIN
    REAL(kind=r8) :: toi(k)
    REAL(kind=r8) :: qoi(k)
    REAL(kind=r8) :: uvi(k,trac+2)   
    REAL(kind=r8) :: TCU(k)
    REAL(kind=r8) :: QCU(k)
    REAL(kind=r8) :: PCU(k)
    REAL(kind=r8) :: clw(k)
    REAL(kind=r8) :: cli(k)  
    REAL(kind=r8) :: QII(k)
    REAL(kind=r8) :: QLI(k)
    REAL(kind=r8) :: PRS(k+1)
    REAL(kind=r8) :: PSJ(k+1)     
    REAL(kind=r8) :: phi_l(k)
    REAL(kind=r8) :: phi_h(k+1)              
    REAL(kind=r8) :: RCU(k,trac+2)
    REAL(kind=r8) :: wfnc
    REAL(kind=r8) :: flx(k+1)
    REAL(kind=r8) :: FLXD(K+1)

    REAL(kind=r8) :: tla
    !REAL(kind=r8) :: pl
    INTEGER       :: irnd,ib

    INTEGER      , PARAMETER :: ICM=100
    !REAL(KIND=r8), PARAMETER :: DAYLEN=86400.0_r8
    !REAL(KIND=r8), PARAMETER :: PFAC=1.0_r8/450.0_r8
    REAL(KIND=r8), PARAMETER :: clwmin=1.0e-10_r8
    INTEGER       :: IC(ICM)
    !
    REAL(kind=r8), ALLOCATABLE ::  ALFINT(:,:)
    !REAL(kind=r8) :: ALFINQ(K)
    REAL(kind=r8) :: PRSM(K)
    REAL(kind=r8) :: PSJM(K)  
    REAL(kind=r8) :: trcfac(trac+2,k)
    REAL(kind=r8) :: alfind(K)
    REAL(kind=r8) :: rhc_l(k)
    REAL(kind=r8) :: dtvd(2,4)
    !REAL(kind=r8) :: CFAC
    REAL(kind=r8) :: TEM
    !REAL(kind=r8) :: dpi
    REAL(kind=r8) :: sgc
    REAL(kind=r8) :: ccwf
    REAL(kind=r8) :: tem1
    REAL(kind=r8) :: tem2
    !
    INTEGER       :: KCR
    INTEGER       :: KFX
    INTEGER       :: NCMX
    INTEGER       :: NC
    INTEGER       :: KTEM
    INTEGER       :: I
    INTEGER       :: L
    INTEGER       :: lm1
    INTEGER       :: ntrc
    !INTEGER       :: ia
    INTEGER       :: ll
    INTEGER       :: km1
    INTEGER       :: kp1
    INTEGER       :: ipt
    !INTEGER       :: lv
    INTEGER       :: KBL
    INTEGER       :: n 
    INTEGER       :: lmhij
    INTEGER       :: KRMIN
    INTEGER       :: KRMAX
    INTEGER       :: KFMAX
    !
    LOGICAL       :: DNDRFT, lprint

    !
    !     locals
    !
    ncrnd=0

    RAIN=0.0_r8
    toi=0.0_r8
    qoi=0.0_r8
    uvi=0.0_r8
    TCU=0.0_r8
    QCU=0.0_r8
    PCU=0.0_r8
    clw=0.0_r8
    cli=0.0_r8
    QII=0.0_r8
    QLI=0.0_r8
    PRS=0.0_r8
    PSJ=0.0_r8  
    phi_l=0.0_r8
    phi_h=0.0_r8
    RCU=0.0_r8
    wfnc=0.0_r8
    flx=0.0_r8
    FLXD=0.0_r8

    tla=0.0_r8
    IC=0
    !
    PRSM=0.0_r8
    PSJM=0.0_r8
    trcfac=0.0_r8
    alfind=0.0_r8
    rhc_l=0.0_r8
    dtvd=0.0_r8
    TEM=0.0_r8
    sgc=0.0_r8
    ccwf=0.0_r8
    tem1=0.0_r8
    tem2=0.0_r8
    !
    KCR=0
    KFX=0
    NCMX=0
    NC=0
    KTEM=0
    I=0
    L=0
    lm1=0
    ntrc=0
    ll=0
    km1=0
    kp1=0
    ipt=0
    KBL=0
    n =0
    lmhij=0
    KRMIN=0
    KRMAX=0
    KFMAX=0
    !
    km1    = k - 1
    kp1    = k + 1
    !
    ntrc = trac
    trcfac(:,:) = 1.0_r8             !  For other tracers
    IF (CUMFRC) THEN
       ntrc = ntrc + 2
       !       trcfac(trac+1) = 0.45_r8       !  For press grad correction c=0.55_r8
       !       trcfac(trac+2) = 0.45_r8       !  in momentum mixing calculations
    ENDIF
    !
    IF (.NOT. ALLOCATED(alfint))THEN
       ALLOCATE(alfint(k,ntrc+4))
       alfint=0.0_r8
    END IF
    !
    CALL set_ras_afc(dt)
    !
    ccwf = 0.5_r8
    DO IPT=1,nCols
       !CALKBL=CALKBLMsk(ipt)
       !
       ! Resolution dependent press grad correction momentum mixing
       !
       IF (CUMFRC) THEN
!!!!      tem = max(0.0, min(1.0, sqrt(sqrt(garea(ipt)*0.25E-10))))
          !         tem = max(0.0, min(1.0, sqrt(garea(ipt)*0.25E-10)))
          !         tem = max(0.1, min(1.0, sqrt(garea(ipt)*0.25E-10)))
          !         tem = max(0.25, min(1.0, sqrt(garea(ipt)*0.25E-10)))
          !         tem = max(0.50, min(1.0, sqrt(garea(ipt)*0.25E-10)))
          !         tem = max(0.45, min(1.0, sqrt(garea(ipt)*0.25E-10))) ! for r2 and rf exp
          !         tem = 1.0                  ! for r1 exp
          !         tem = 0.45                 ! for r6 exp

          !         trcfac(trac+1,l) = tem       !  For press grad correction c=0.55
          !         trcfac(trac+2,l) = tem       !  in momentum mixing calculations
          !
          IF (ccwfac(ipt) >= 0.0_r8) ccwf = ccwfac(ipt)
       ENDIF
       DO l=1,k
          ud_mf(ipt,l)  = 0.0_r8
          dd_mf(ipt,l)  = 0.0_r8
          det_mf(ipt,l) = 0.0_r8
       ENDDO
       !
       !     Compute NCRND  : here LMH is the number of layers above the
       !                      bottom surface.  For sigma coordinate LMH=K.
       !
       LMHIJ = LMH(ipt)
       IF (flipv) THEN
          ll  = kp1 - LMH(ipt)
          tem = 1.0_r8 / prsi(ipt,ll)
       ELSE
          ll  = LMH(ipt)
          tem = 1.0_r8 / prsi(ipt,ll+1)
       ENDIF
       KRMIN = 1
       KRMAX = km1
       KFMAX = KRMAX
       DO L=1,LMHIJ-1
          ll = l
          IF (flipv) ll = kp1 -l ! Input variables are bottom to top!
          SGC = prsl(ipt,ll) * tem
          IF (SGC .LE. 0.050_r8) KRMIN = L
          !         IF (SGC .LE. 0.600_r8) KRMAX = L
          IF (SGC .LE. 0.700_r8) KRMAX = L
          !         IF (SGC .LE. 0.800_r8) KRMAX = L
          !!        IF (SGC .LE. 0.760_r8) KRMAX = L
          !         IF (SGC .LE. 0.930_r8) KFMAX = L
          IF (SGC .LE. 0.970_r8) KFMAX = L    ! Commented on 20060202
          !LL2      IF (SGC .LE. 0.950_r8) KFMAX = L
       ENDDO
       !     if (lprnt .and. ipt .eq. ipr) print *,' krmin=',krmin,' krmax=',
       !    &krmax,' kfmax=',kfmax,' lmhij=',lmhij,' tem=',tem
       !
       IF (fix_ncld_hr) THEN
          NCRND = INT(MIN(nrcmax, (KRMAX-KRMIN+1)) * (DTF/1200) + 0.50001_r8,KIND=i4)
          !         NCRND = min(nrcmax, (KRMAX-KRMIN+1)) * (DTF/1800) + 0.50001_r8
          facdt = delt_c / dt
       ELSE
          NCRND = INT(MIN(nrcmax, (KRMAX-KRMIN+1)),KIND=i4)
          facdt = 1.0_r8 
       ENDIF
       IF (DT .GT. DTF) NCRND = (5*NCRND) / 4
       NCRND   = MAX(NCRND, 1)
       !
       KCR    = MIN(LMHIJ,KRMAX)
       KTEM   = MIN(LMHIJ,KFMAX)
       KFX    = KTEM - KCR
       !     if(lprnt)print*,' enter RASCNV k=',k,' ktem=',ktem,' LMHIJ='
       !    &,                 LMHIJ
       !    &,               ' krmax=',krmax,' kfmax=',kfmax
       !    &,               ' kcr=',kcr, ' cdrag=',cdrag(ipr)

       IF (KFX .GT. 0) THEN
          IF (BOTOP) THEN
             DO NC=1,KFX
                IC(NC) = KTEM + 1 - NC
             ENDDO
          ELSE
             DO NC=KFX,1,-1
                IC(NC) = KTEM + 1 - NC
             ENDDO
          ENDIF
       ENDIF
       !
       NCMX  = KFX + NCRND
       IF (NCRND .GT. 0) THEN
          DO I=1,NCRND
             IRND = INT((RANNUM(ipt,I)-0.0005_r8)*(KCR-KRMIN+1),KIND=i4)
             IC(KFX+I) = IRND + KRMIN
          ENDDO
       ENDIF
       !
       !     ia = 1
       !
       !     print *,' in rascnv: k=',k,'lat=',lat,' lprnt=',lprnt
       !     if (lprnt) then
       !        if (me .eq. 0) then
       !        print *,' tin',(tin(ia,l),l=k,1,-1)
       !        print *,' qin',(qin(ia,l),l=k,1,-1)
       !     endif
       !
       !
       !lprint = lprnt .AND. ipt .EQ. ipr
       lprint = lprnt
       !       kuo(ipt)  = 0
       DO l=1,k
          ll = l
          IF (flipv) ll = kp1 -l ! Input variables are bottom to top!
          CLW(l)     = 0.0_r8       ! Assumes initial value of Cloud water
          CLI(l)     = 0.0_r8       ! Assumes initial value of Cloud ice
          ! to be zero i.e. no environmental condensate!!!
          !         CLT(ipt,l) = 0.0_r8
          QII(l)     = 0.0_r8
          QLI(l)     = 0.0_r8
          !                          Initialize heating, drying, cloudiness etc.
          tcu(l)     = 0.0_r8
          qcu(l)     = 0.0_r8
          pcu(l)     = 0.0_r8
          flx(l)     = 0.0_r8
          flxd(l)    = 0.0_r8
          rcu(l,1)   = 0.0_r8
          rcu(l,2)   = 0.0_r8
          !                          Transfer input prognostic data into local variable
          toi(l)     = tin(ipt,ll)
          qoi(l)     = qin(ipt,ll)
          uvi(l,trac+1) = uin(ipt,ll)
          uvi(l,trac+2) = vin(ipt,ll)
          !
          DO n=1,trac
             uvi(l,n) = ccin(ipt,ll,n+2)
          ENDDO
          !
       ENDDO
       flx(k+1)  = 0.0_r8
       flxd(k+1) = 0.0_r8
       !
       IF (ccin(ipt,1,2) .LE. -999.0_r8) THEN
          DO l=1,k
             ll = l
             IF (flipv) ll = kp1 -l ! Input variables are bottom to top!
             !PK tem = ccin(ipt,ll,1)                                      &
             !PK     &            * MAX(ZERO, MIN(ONE, (TCR-toi(L))*TCRF))
             !PK ccin(ipt,ll,2) = ccin(ipt,ll,1) - tem
             !PK ccin(ipt,ll,1) = tem
             tem = qice(ipt,ll)                                      &
                  &            * MAX(ZERO, MIN(ONE, (TCR-toi(L))*TCRF))
             qliq(ipt,ll) = qice(ipt,ll) - tem
             qice(ipt,ll) = tem

          ENDDO
       ENDIF

       IF (advcld) THEN
          DO l=1,k
             ll = l
             IF (flipv) ll = kp1 -l ! Input variables are bottom to top!
             !PK QII(L) = ccin(ipt,ll,1)
             !PK QLI(L) = ccin(ipt,ll,2)
             QII(L) = qice(ipt,ll)
             QLI(L) = qliq(ipt,ll)

          ENDDO
       ENDIF
       !
       KBL  = KPBL(ipt)
       IF (flipv) KBL  = MAX(MIN(k, kp1-KPBL(ipt)), k/2)
       rain = 0.0_r8
       !
       DO L=1,kp1
          ll = l
          IF (flipv) ll = kp1 + 1 - l      ! Input variables are bottom to top!
          PRS(LL)   = prsi(ipt, L) * facmb ! facmb is for conversion to MB
          PSJ(LL)   = prsik(ipt,L)
          phi_h(LL) = phii(ipt,L)
       ENDDO
       !
       DO L=1,k
          ll = l
          IF (flipv) ll = kp1 - l          ! Input variables are bottom to top!
          PRSM(LL)  = prsl(ipt, L) * facmb ! facmb is for conversion to MB
          PSJM(LL)  = prslk(ipt,L)
          phi_l(LL) = phil(ipt,L)
          rhc_l(LL) = rhc(ipt,L)
          !
          !         rhc_l(ll) = 1.0_r8
       ENDDO
       !
       !     if(lprint) print *,' PRS=',PRS
       !     if(lprint) print *,' PRSM=',PRSM
       !     if (lprint) then
       !        print *,' qns=',qns(ia),' qoi=',qn0(ia,k),'qin=',qin(ia,1)
       !        if (me .eq. 0) then
       !        print *,' toi',(tn0(ia,l),l=1,k)
       !        print *,' qoi',(qn0(ia,l),l=1,k),' kbl=',kbl
       !     endif
       !
       !
       !!      PSJP(KP1) = PSJ(KP1) * PRS(KP1) * RKPP1I
       !       do l=k,kctop(1),-1
       !         DPI     = RKPP1I / (PRS(L+1) - PRS(L))
       !         PSJM(L) = (PSJ(L+1)*PRS(L+1) - PSJ(L)*PRS(L)) * DPI
       !!        PSJP(L) = PSJ(L) * PRS(L) * RKPP1I
       !!        DPI(L)  = 1.0 / (PRS(L+1) - PRS(L))
       !!        PSJM(L) = (PSJP(L+1) - PSJP(L)) * DPI(L)
       !         PRSM(L) = 1000.0_r8 * PSJM(L) ** (1.0_r8/rkap)
       !!        PRSM(L) = 1000.0_r8 * PSJM(L) ** rkapi
       !!        PRSM(L) = 0.5_r8 * (prs(L+1)+prs(L))
       !       enddo
       !
       !
       !     
       !     print *,' ipt=',ipt
       alfint(:,:) = 0.5_r8              ! For second order scheme
       alfind(:)   = 0.5_r8
       IF (advups) THEN               ! For first order upstream for updraft
          alfint(:,:) = 1.0_r8
       ELSEIF (advtvd) THEN           ! TVD flux limiter scheme for updraft
          alfint(:,:) = 1.0_r8
          l   = krmin
          lm1 = l - 1
          dtvd(1,1) = con_cp*(toi(l)-toi(lm1)) + phi_l(l)-phi_l(lm1)        &
               &              + con_hvap *(qoi(l)-qoi(lm1))
          dtvd(1,2) = qoi(l) - qoi(lm1)
          dtvd(1,3) = qli(l) - qli(lm1)
          dtvd(1,4) = qii(l) - qii(lm1)
          DO l=krmin+1,k
             lm1 = l - 1
             !     print *,' toi=',toi(l),toi(lm1),' phi_l=',phi_l(l),phi_l(lm1)
             !    &,' qoi=',qoi(l),qoi(lm1),' con_cp=',con_cp,' con_hvap =',con_hvap 
             dtvd(2,1)   = con_cp*(toi(l)-toi(lm1)) + phi_l(l)-phi_l(lm1)    &
                  &                  + con_hvap *(qoi(l)-qoi(lm1))
             !     print *,' l=',l,' dtvd=',dtvd(:,1)
             IF (ABS(dtvd(2,1)) > 1.0e-10_r8) THEN
                tem1        = dtvd(1,1) / dtvd(2,1)
                tem2        = ABS(tem1)
                alfint(l,1) = 1.0_r8 - 0.5_r8*(tem1 + tem2)/(1.0_r8 + tem2)   ! for h
             ENDIF
             !     print *,' alfint=',alfint(l,1),' l=',l,' ipt=',ipt
             dtvd(1,1)   = dtvd(2,1)
             !
             dtvd(2,2)   = qoi(l) - qoi(lm1)
             !     print *,' l=',l,' dtvd2=',dtvd(:,2)
             IF (ABS(dtvd(2,2)) > 1.0e-10_r8) THEN
                tem1        = dtvd(1,2) / dtvd(2,2)
                tem2        = ABS(tem1)
                alfint(l,2) = 1.0_r8 - 0.5_r8*(tem1 + tem2)/(1.0_r8 + tem2)   ! for q
             ENDIF
             dtvd(1,2)   = dtvd(2,2)
             !
             dtvd(2,3)   = qli(l) - qli(lm1)
             !     print *,' l=',l,' dtvd3=',dtvd(:,3)
             IF (ABS(dtvd(2,3)) > 1.0e-10_r8) THEN
                tem1        = dtvd(1,3) / dtvd(2,3)
                tem2        = ABS(tem1)
                alfint(l,3) = 1.0_r8 - 0.5_r8*(tem1 + tem2)/(1.0_r8 + tem2)   ! for ql
             ENDIF
             dtvd(1,3)   = dtvd(2,3)
             !
             dtvd(2,4)   = qii(l) - qii(lm1)
             !     print *,' l=',l,' dtvd4=',dtvd(:,4)
             IF (ABS(dtvd(2,4)) > 1.0e-10_r8) THEN
                tem1        = dtvd(1,4) / dtvd(2,4)
                tem2        = ABS(tem1)
                alfint(l,4) = 1.0_r8 - 0.5_r8*(tem1 + tem2)/(1.0_r8 + tem2)   ! for qi
             ENDIF
             dtvd(1,4)   = dtvd(2,4)
          ENDDO
          !
          IF (ntrc > 0) THEN
             DO n=1,ntrc
                l = krmin
                dtvd(1,1)   = uvi(l,n) - uvi(l-1,n)
                DO l=krmin+1,k
                   dtvd(2,1)     = uvi(l,n) - uvi(l-1,n)
                   !     print *,' l=',l,' dtvdn=',dtvd(:,1),' n=',n,' l=',l
                   IF (ABS(dtvd(2,1)) > 1.0e-10_r8) THEN
                      tem1          = dtvd(1,1) / dtvd(2,1)
                      tem2          = ABS(tem1)
                      alfint(l,n+4) = 1.0_r8 - 0.5_r8*(tem1 + tem2)/(1.0_r8 + tem2) ! for tracers
                   ENDIF
                   dtvd(1,1)     = dtvd(2,1)
                ENDDO
             ENDDO
          ENDIF
       ENDIF
       !
       !     print *,' after alfint for ipt=',ipt
       IF (CUMFRC) THEN

          DO l=krmin,k
             tem = 1.0_r8 - MAX(pgfbot, MIN(pgftop, pgftop+pgfgrad*prsm(l)))
             trcfac(trac+1,l) = tem
             trcfac(trac+2,l) = tem
          ENDDO
       ENDIF
       !
       lprint = lprnt !.AND. ipt .EQ. ipr
       !     if (lprint) then
       !       print *,' trcfac=',trcfac(1+trac,krmin:k)
       !       print *,' alfint=',alfint(krmin:k,1)
       !       print *,' alfinq=',alfint(krmin:k,2)
       !       print *,' alfini=',alfint(krmin:k,4)
       !       print *,' alfinu=',alfint(krmin:k,5)
       !     endif
       !
       IF (calkbl) kbl = k
       DO NC=1,NCMX
          !
          IB = IC(NC)
          IF (ib .GT. kbl) CYCLE
          !       lprint = lprnt .and. ipt .eq. ipr
          !       lprint = lprnt .and. ipt .eq. ipr .and. ib .eq. 41
          !
          DNDRFT = DD_DP .GT. 0.0_r8
          !
          !     if (lprint) print *,' calling cloud type ib=',ib,' kbl=',kbl
          !    *,' kpbl=',kpbl,' alfint=',alfint,' frac=',frac
          !    *,' ntrc=',ntrc,' ipt=',ipt
          !
          !       if (advtvd) then           ! TVD flux limiter scheme for updraft
          !         l   = ib
          !         lm1 = l - 1
          !         dtvd(1,1) = con_cp*(toi(l)-toi(lm1)) + phi_l(l)-phi_l(lm1)
          !    &              + con_hvap *(qoi(l)-qoi(lm1))
          !         dtvd(1,2) = qoi(l) - qoi(lm1)
          !         dtvd(1,3) = qli(l) - qli(lm1)
          !         dtvd(1,4) = qii(l) - qii(lm1)
          !         do l=ib+1,k
          !           lm1 = l - 1
          !           dtvd(2,1)   = con_cp*(toi(l)-toi(lm1)) + phi_l(l)-phi_l(lm1)
          !    &                  + con_hvap *(qoi(l)-qoi(lm1))
          !           if (abs(dtvd(2,1)) > 1.0e-10_r8) then
          !             tem1        = dtvd(1,1) / dtvd(2,1)
          !             tem2        = abs(tem1)
          !             alfint(l,1) = 1.0_r8 - 0.5_r8*(tem1 + tem2)/(1.0_r8 + tem2)   ! for h
          !           endif
          !           dtvd(1,1)   = dtvd(2,1)
          !
          !           dtvd(2,2)   = qoi(l) - qoi(lm1)
          !           if (abs(dtvd(2,2)) > 1.0e-10_r8) then
          !             tem1        = dtvd(1,2) / dtvd(2,2)
          !             tem2        = abs(tem1)
          !             alfint(l,2) = 1.0_r8 - 0.5_r8*(tem1 + tem2)/(1.0_r8 + tem2)   ! for q
          !           endif
          !           dtvd(1,2)   = dtvd(2,2)
          !
          !           dtvd(2,3)   = qli(l) - qli(lm1)
          !           if (abs(dtvd(2,3)) > 1.0e-10_r8) then
          !             tem1        = dtvd(1,3) / dtvd(2,3)
          !             tem2        = abs(tem1)
          !             alfint(l,3) = 1.0_r8 - 0.5_r8*(tem1 + tem2)/(1.0_r8 + tem2)   ! for ql
          !           endif
          !           dtvd(1,3)   = dtvd(2,3)
          !
          !           dtvd(2,4)   = qii(l) - qii(lm1)
          !           if (abs(dtvd(2,4)) > 1.0e-10_r8) then
          !             tem1        = dtvd(1,4) / dtvd(2,4)
          !             tem2        = abs(tem1)
          !             alfint(l,4) = 1.0_r8 - 0.5_r8*(tem1 + tem2)/(1.0_r8 + tem2)   ! for qi
          !           endif
          !           dtvd(1,4)   = dtvd(2,4)
          !         enddo
          !
          !         if (ntrc > 0) then
          !           do n=1,ntrc
          !             l = ib
          !             dtvd(1,1)   = uvi(l,n) - uvi(l-1,n)
          !             do l=ib+1,k
          !               dtvd(2,1)     = uvi(l,n) - uvi(l-1,n)
          !               if (abs(dtvd(2,1)) > 1.0e-10_r8) then
          !                 tem1        = dtvd(1,1) / dtvd(2,1)
          !                 tem2          = abs(tem1)
          !                 alfint(l,n+4) = 1.0_r8 - 0.5_r8*(tem1 + tem2)/(1.0_r8 + tem2) ! for tracers
          !               endif
          !               dtvd(1,1)     = dtvd(2,1)
          !             enddo
          !           enddo
          !         endif
          !       endif
          !
          !
          !     if (lprint) then
          !     ia = ipt
          !     print *,' toi=',(toi(ia,l),l=1,K)
          !     print *,' qoi=',(qoi(ia,l),l=1,K),' kbl=',kbl
          !     print *,' toi=',(toi(l),l=1,K)
          !     print *,' qoi=',(qoi(l),l=1,K),' kbl=',kbl
          !     print *,' prs=',(prs(l),l=1,K)
          !     endif
          !
          WFNC = 0.0_r8
          DO L=IB,K+1
             FLX(L)    = 0.0_r8
             FLXD(L)   = 0.0_r8
          ENDDO
          !
          !
          !     if (me .eq. 0) then
          !     if(lprint)then
          !     print *,' CALLING CLOUD TYPE IB= ', IB,' DT=',DT,' K=',K
          !    &, 'ipt=',ipt
          !     print *,' TOI=',(TOI(L),L=IB,K)
          !     print *,' QOI=',(QOI(L),L=IB,K)
          !     endif
          !     print *,' alft=',alfint
          !
          TLA = -10.0_r8
          !
          !     if (lprint) print *,' qliin=',qli
          !     if (lprint) print *,' qiiin=',qii
          CALL CLOUD(lmhij, IB, ntrc                                    &
               &,               FRAC,  MAX_NEG_BOUY                     &
!               &,              RASAL(IB), FRAC,  MAX_NEG_BOUY                     &

               &,              ALFINT, rhfacs, garea(ipt)         &
               !    &,              ALFINT, ALFINQ, rhfacl, rhfacs, garea(ipt)         &
               &,              alfind, rhc_l                                      &
               !
               &,              TOI, QOI, UVI, PRS, PRSM, phi_l, phi_h             &
               !    &,              TOI, QOI, UVI, PRS, PRSM, PSJ, PSJM
               !    &,              TOI, QOI, UVI, PRS, PRSM, PSJ, PSJM, DPI
               !    &,              TOI, QOI, UVI, PRS, PSJ
               &,              QLI, QII, KBL, DDVEL(ipt)                          &
               &,              CDRAG(ipt),lprint, trcfac, ccwf                    &
               !    &,              IDIAG, lprnt
               &,              TCU, QCU, RCU, PCU, FLX, FLXD                      &
               &,              RAIN, REVAP, DT                                    &
               &,              WFNC, WRKFUN, CALKBL, CRTFUN, TLA, DNDRFT, DD_DP)     
          !    &,              WFNC, WRKFUN, CALKBL, CRTFUN, TLA, DNDRFT, UPDRET)
          !     if (lprint) print *,' rain=',rain,' ipt=',ipt
          !     if (me .eq. 0) then
          IF (lprint) THEN
             PRINT *,' after calling CLOUD TYPE IB= ', IB                      &
                  &,' rain=',rain,' prskd=',prs(ib),' qli=',qli(ib),' qii=',qii(ib)
             !     print *,' TOI=',(TOI(L),L=1,K),' me=',me,' ib=',ib
             !     print *,' QOI=',(QOI(L),L=1,K)
          ENDIF
          !     if (lprint) print *,' qliou=',qli
          !     if (lprint) print *,' qiiou=',qii
          !
          DO L=IB,K
             ll = l
             IF (flipv) ll  = kp1 -l    ! Input variables are bottom to top!
             ud_mf(ipt,ll)  = ud_mf(ipt,ll)  + flx(l+1)
             dd_mf(ipt,ll)  = dd_mf(ipt,ll)  + flxd(l+1)
          ENDDO
          ll = ib
          IF (flipv) ll  = kp1 - ib
          det_mf(ipt,ll) = det_mf(ipt,ll) + flx(ib)
          ! 
          !     Compute cloud amounts for the Goddard radiation
          !
          !         IF (FLX(KBL) .GT. 0.0_r8) THEN
          !           PL   = 0.5_r8 * (PRS(IB) + PRS(IB+1))
          !           CFAC = MIN(1.0_r8, MAX(0.0_r8, (850.0_r8-PL)*PFAC))
          !         ELSE
          !           CFAC = 0.0_r8
          !         ENDIF
          !
          !   Warining!!!!
          !   ------------
          !   By doing the following, CLOUD does not contain environmental
          !   condensate!
          !
          IF (.NOT. advcld) THEN
             DO l=1,K
                !             clw(l ) = clw(l) + QLI(L) + QII(L)
                clw(l ) = clw(l) + QLI(L)
                cli(l ) = cli(l) + QII(L)
                QLI(L)  = 0.0_r8
                QII(L)  = 0.0_r8
             ENDDO
          ENDIF
          !
       ENDDO                      ! End of the NC loop!
       !
       RAINC(ipt) = rain * 0.001_r8    ! Output rain is in meters
       !     if(lprint)print*,' convective precip=',rain*86400/dt,' mm/day'
       !    1,               ' ipt=',ipt
       !
       !     if (lprint) then
       !        print *,' toi',(tn0(imax,l),l=1,k)
       !        print *,' qoi',(qn0(imax,l),l=1,k)
       !     endif
       !
       DO l=1,k
          ll = l
          IF (flipv) ll  = kp1 - l
          tin(ipt,ll)    = toi(l)                   ! Temperature
          qin(ipt,ll)    = qoi(l)                   ! Specific humidity
          uin(ipt,ll)    = uvi(l,trac+1)            ! U momentum
          vin(ipt,ll)    = uvi(l,trac+2)            ! V momentum
          !         clw(l)         = clw(l) + qli(l) + qii(l) ! Cloud condensate
          !         ccin(ipt,ll,1) = ccin(ipt,ll,1) + clw(l)
          DO n=1,trac
             ccin(ipt,ll,n+2) = uvi(l,n)             ! Tracers
          ENDDO
       ENDDO
       IF (advcld) THEN
          DO l=1,k
             ll = l
             IF (flipv) ll  = kp1 - l
             !           ccin(ipt,ll,1) = qli(l) + qii(l) ! Cloud condensate
             !PK ccin(ipt,ll,1) = qii(l)          ! Cloud ice
             !PK ccin(ipt,ll,2) = qli(l)          ! Cloud water
             qice(ipt,ll) = qii(l)          ! Cloud ice
             qliq(ipt,ll) = qli(l)          ! Cloud water

          ENDDO
       ELSE
          DO l=1,k
             ll = l
             IF (flipv) ll  = kp1 - l
             !           ccin(ipt,ll,1) = ccin(ipt,ll,1) + clw(l)
             !PK ccin(ipt,ll,1) = ccin(ipt,ll,1) + cli(l)
             !PK ccin(ipt,ll,2) = ccin(ipt,ll,2) + clw(l)
             qice(ipt,ll) = qice(ipt,ll) + cli(l)
             qliq(ipt,ll) = qliq(ipt,ll) + clw(l)

          ENDDO
       ENDIF
       !
       !       kuo(ipt)  = 0
       !
       ktop(ipt) = kp1
       kbot(ipt) = 0

       DO l=lmhij-1,1,-1
          IF (prs(lmhij+1)-prs(l) .GT. 250.0_r8 .AND. tcu(l) .NE. 0.0_r8) THEN ! for r1 &rf
             !         if (prs(lmhij+1)-prs(l) .gt. 100.0_r8 .and. tcu(l) .ne. 0.0_r8) then ! for r1 &rf
             !1        if (prs(lmhij+1)-prs(l) .gt. 300.0_r8 .and. tcu(l) .ne. 0.0_r8) then
             !         if (prs(lmhij+1)-prs(l) .gt. 500.0_r8 .and. tcu(l) .ne. 0.0_r8) then
             !         if (prs(lmhij+1)-prs(l) .gt. 400.0_r8 .and. tcu(l) .ne. 0.0_r8) then
             !         if (prs(kp1)-prs(l) .gt. 500.0_r8 .and. tcu(l) .ne. 0.0_r8) then ! for r2 exp
             !         if (prs(kp1)-prs(l) .gt. 400.0_r8 .and. tcu(l) .ne. 0.0_r8) then
             !         if (prs(lmhij+1)-prs(l) .gt. 200.0_r8 .and. tcu(l) .ne. 0.0_r8) then
             !         if (prsm(l) .lt. 900.0_r8 .and. tcu(l) .ne. 0.0_r8) then
             !         if (phi_l(l) .gt. 10000.0_r8 .and. tcu(l) .ne. 0.0_r8) then
             kuo(ipt) = 1
          ENDIF
          !  New test for convective clouds ! added in 08/21/96
          IF (clw(l)+cli(l) .GT. 0.0_r8 .OR.                               &
               &        qli(l)+qii(l) .GT. clwmin) ktop(ipt) = l
       ENDDO
       DO l=1,km1
          IF (clw(l)+cli(l) .GT. 0.0_r8 .OR.                               &
               &        qli(l)+qii(l) .GT. clwmin) kbot(ipt) = l
       ENDDO
       IF (flipv) THEN
          ktop(ipt) = kp1 - ktop(ipt)
          kbot(ipt) = kp1 - kbot(ipt)
       ENDIF
       !
       !     if (lprint) then
       !        print *,' tin',(tin(ia,l),l=k,1,-1)
       !        print *,' qin',(qin(ia,l),l=k,1,-1)
       !     endif
       !
       !     Velocity scale from the downdraft!
       !
       DDVEL(ipt) = DDVEL(ipt) * DDFAC * con_g / (prs(K+1)-prs(k))
       !
    ENDDO                            ! End of the IPT Loop!
    DEALLOCATE (alfint)
    !
    RETURN
  END  SUBROUTINE rascnv

  !-----------------------------------------------------------------------------------------

  SUBROUTINE CRTWRK(PL, CCWF, ACR)
    !use module_ras , only : ac, ad
    IMPLICIT NONE
    !
    REAL(kind=r8), INTENT(IN   ) :: PL, CCWF
    REAL(kind=r8), INTENT(OUT  ) :: ACR
    INTEGER :: IWK
    !
    ACR=0.0_r8
    IWK = INT(PL * 0.02_r8 - 0.999999999_r8,KIND=i4)
    IWK = MAX(1, MIN(IWK,16))
    ACR = (AC(IWK) + PL * AD(IWK)) * CCWF
    !
    RETURN
  END SUBROUTINE CRTWRK


  !-----------------------------------------------------------------------------------------

  SUBROUTINE CLOUD(                                                 &
       &                  K, KD, M                                        &
       &,                  FRACBL, MAX_NEG_BOUY                    &
!       &,                 RASALF, FRACBL, MAX_NEG_BOUY                    &

       &,                 ALFINT,       RHFACS, garea           &
                                !    &,                 ALFINT, ALFINQ, RHFACL, RHFACS, garea           &
       &,                 alfind, rhc_ls                                  &
       
       &,                 TOI, QOI, ROI, PRS, PRSM, phil, phih            &
                                !    &,                 TOI, QOI, ROI, PRS, PRSM, PRJ, PRJM, DPI        &
                                !    &,                 TOI, QOI, ROI, PRS, PRJ                         &
       &,                 QLI, QII, KPBL, DSFC                            &
       &,                 CD,lprnt, trcfac,ccwf                           &
                                !    &,                 IDIAG, lprnt                                    &
       
       &,                 TCU, QCU, RCU, PCU, FLX, FLXD                   &
                                !    &,                 TCD, QCD                                        &
       &,                 CUP, REVAP, DT                                  &
       &,                 WFNC, WRKFUN, CALKBL, CRTFUN, TLA, DNDRFT, DPD)  

    !
    !***********************************************************************
    !******************** Relaxed  Arakawa-Schubert ************************
    !****************** Plug Compatible Scalar Version *********************
    !************************ SUBROUTINE CLOUD  ****************************
    !************************  October 2004     ****************************
    !********************  VERSION 2.0  (modified) *************************
    !************* Shrinivas.Moorthi@noaa.gov (301) 763 8000(X7233) ********
    !***********************************************************************
    !*Reference:
    !-----------
    !     NOAA Technical Report NWS/NCEP 99-01:
    !     Documentation of Version 2 of Relaxed-Arakawa-Schubert
    !     Cumulus Parameterization with Convective Downdrafts, June 1999.
    !     by S. Moorthi and M. J. Suarez.
    !
    !***********************************************************************
    !
    !===>    UPDATES CLOUD TENDENCIES DUE TO A SINGLE CLOUD
    !===>    DETRAINING AT LEVEL KD.
    !
    !***********************************************************************
    !
    !===>  TOI(K)     INOUT   TEMPERATURE             KELVIN
    !===>  QOI(K)     INOUT   SPECIFIC HUMIDITY       NON-DIMENSIONAL
    !===>  ROI(K,M)   INOUT   TRACER                  ARBITRARY
    !===>  QLI(K)     INOUT   LIQUID WATER            NON-DIMENSIONAL
    !===>  QII(K)     INOUT   ICE                     NON-DIMENSIONAL

    !===>  PRS(K+1)   INPUT   PRESSURE @ EDGES        MB
    !===>  PRSM(K)    INPUT   PRESSURE @ LAYERS       MB
    !===>  PHIH(K+1)  INPUT   GEOPOTENTIAL @ EDGES  IN MKS units
    !===>  PHIL(K)    INPUT   GEOPOTENTIAL @ LAYERS IN MKS units
    !===>  PRJ(K+1)   INPUT   (P/P0)^KAPPA  @ EDGES   NON-DIMENSIONAL
    !===>  PRJM(K)    INPUT   (P/P0)^KAPPA  @ LAYERS  NON-DIMENSIONAL

    !===>  K      INPUT   THE RISE & THE INDEX OF THE SUBCLOUD LAYER
    !===>  KD     INPUT   DETRAINMENT LEVEL ( 1<= KD < K )          
    !===>  M      INPUT   NUMBER OF TRACERS. MAY BE ZERO.
    !===>  DNDRFT INPUT   LOGICAL .TRUE. OR .FALSE.
    !===>  DPD    INPUT   Minumum Cloud Depth for DOWNDRFAT Computation hPa
    !
    !===>  TCU(K  )   UPDATE  TEMPERATURE TENDENCY       DEG
    !===>  QCU(K  )   UPDATE  WATER VAPOR TENDENCY       (G/G)
    !===>  RCU(K,M)   UPDATE  TRACER TENDENCIES          ND
    !===>  PCU(K-1)   UPDATE  PRECIP @ BASE OF LAYER     KG/M^2
    !===>  FLX(K  )   UPDATE  MASS FLUX @ TOP OF LAYER   KG/M^2
    !===>  CUP        UPDATE  PRECIPITATION AT THE SURFACE KG/M^2
    !
    !      use module_ras
    IMPLICIT NONE
    !
    !  INPUT ARGUMENTS
    INTEGER      , INTENT(IN   ) :: K
    INTEGER      , INTENT(IN   ) :: KD
    INTEGER      , INTENT(IN   ) :: M   !trac +2
    !REAL(kind=r8), INTENT(IN   ) :: RASALF   ! not used
    REAL(kind=r8), INTENT(IN   ) :: FRACBL
    REAL(kind=r8), INTENT(IN   ) :: MAX_NEG_BOUY
    REAL(kind=r8), INTENT(IN   ) :: ALFINT(K,M+4)
    !REAL(kind=r8), INTENT(IN   ) :: RHFACL ! not used 
    REAL(kind=r8), INTENT(IN   ) :: RHFACS 
    REAL(kind=r8), INTENT(IN   ) :: garea
    REAL(kind=r8), INTENT(IN   ) :: alfind(k)
    REAL(kind=r8), INTENT(IN   ) :: rhc_ls(k)
    REAL(kind=r8), INTENT(INOUT) :: TOI(K) 
    REAL(kind=r8), INTENT(INOUT) :: QOI(K )
    REAL(kind=r8), INTENT(INOUT) :: ROI(K,M)
    REAL(kind=r8), INTENT(IN   ) :: PRS(K+1)
    REAL(kind=r8), INTENT(IN   ) :: PRSM(K) 
    REAL(kind=r8), INTENT(IN   ) :: PHIL(K)
    REAL(kind=r8), INTENT(IN   ) :: PHIH(K+1) 
    REAL(kind=r8), INTENT(INOUT) :: QLI(K)
    REAL(kind=r8), INTENT(INOUT) :: QII(K)
    INTEGER      , INTENT(INOUT) :: KPBL
    REAL(kind=r8), INTENT(INOUT) :: DSFC
    REAL(kind=r8), INTENT(IN   ) :: CD
    LOGICAL      , INTENT(IN   ) :: lprnt
    REAL(kind=r8), INTENT(IN   ) :: trcfac(M,k)
    REAL(kind=r8), INTENT(IN   ) :: ccwf

    REAL(kind=r8), INTENT(INOUT) :: TCU(K)
    REAL(kind=r8), INTENT(INOUT) :: QCU(K)
    REAL(kind=r8), INTENT(INOUT) :: RCU(K,M)
    REAL(kind=r8), INTENT(INOUT) :: PCU(K) 
    REAL(kind=r8), INTENT(INOUT) :: FLX(K+1)
    REAL(kind=r8), INTENT(INOUT) :: FLXD(K+1)
    REAL(kind=r8), INTENT(INOUT) :: CUP

    LOGICAL      , INTENT(IN   ) :: REVAP
    REAL(kind=r8), INTENT(IN   ) :: DT
    REAL(kind=r8), INTENT(INOUT) :: WFNC
    LOGICAL      , INTENT(IN   ) :: WRKFUN
    LOGICAL      , INTENT(IN   ) :: CALKBL
    LOGICAL      , INTENT(IN   ) :: CRTFUN
    REAL(kind=r8), INTENT(INOUT) :: TLA
    LOGICAL      , INTENT(IN   ) ::  DNDRFT
    REAL(kind=r8), INTENT(IN   ) :: DPD


    LOGICAL       :: CALCUP


    !REAL(kind=r8) :: UFN
    INTEGER       :: KBL
    INTEGER       :: KB1

    !     real(kind=r8) RASALF, FRACBL, MAX_NEG_BOUY, ALFINT(K),     &
    !     real(kind=r8) ALFINQ(K), DPD, alfind(k), rhc_ls(k)

    !  UPDATE ARGUMENTS

    REAL(kind=r8)  ::    TCD(K),   QCD(K)

    !  TEMPORARY WORK SPACE

    REAL(kind=r8) :: HOL(KD:K)
    REAL(kind=r8) :: QOL(KD:K)
    REAL(kind=r8) :: GAF(KD:K+1)
    REAL(kind=r8) :: HST(KD:K)
    REAL(kind=r8) :: QST(KD:K)
    REAL(kind=r8) :: TOL(KD:K)
    REAL(kind=r8) :: GMH(KD:K)
    REAL(kind=r8) :: GMS(KD:K+1)
    REAL(kind=r8) :: GAM(KD:K+1)
    REAL(kind=r8) :: AKT(KD:K)
    REAL(kind=r8) :: AKC(KD:K)
    REAL(kind=r8) :: BKC(KD:K)
    REAL(kind=r8) :: LTL(KD:K)
    REAL(kind=r8) :: RNN(KD:K)
    REAL(kind=r8) :: FCO(KD:K)
    REAL(kind=r8) :: PRI(KD:K)
    !REAL(kind=r8) :: PRH(KD:K)
    REAL(kind=r8) :: QIL(KD:K)
    REAL(kind=r8) :: QLL(KD:K)
    REAL(kind=r8) :: ZET(KD:K)
    REAL(kind=r8) :: XI(KD:K)
    REAL(kind=r8) :: RNS(KD:K)
    REAL(kind=r8) :: Q0U(KD:K)
    REAL(kind=r8) :: Q0D(KD:K)
    REAL(kind=r8) :: vtf(KD:K)
    REAL(kind=r8) :: DLB(KD:K+1)
    REAL(kind=r8) :: DLT(KD:K+1)
    REAL(kind=r8) :: ETA(KD:K+1)
    REAL(kind=r8) :: PRL(KD:K+1)
    REAL(kind=r8) :: CIL(KD:K)
    REAL(kind=r8) :: CLL(KD:K)
    REAL(kind=r8) :: ETAI(KD:K)

    REAL(kind=r8) :: ALM
    REAL(kind=r8) :: DET
    REAL(kind=r8) :: HCC
    REAL(kind=r8) :: CLP
    REAL(kind=r8) :: HSU
    REAL(kind=r8) :: HSD
    REAL(kind=r8) :: QTL
    REAL(kind=r8) :: QTV
    REAL(kind=r8) :: AKM
    REAL(kind=r8) :: WFN
    REAL(kind=r8) :: HOS
    REAL(kind=r8) :: QOS
    REAL(kind=r8) :: AMB
    REAL(kind=r8) :: TX1
    REAL(kind=r8) :: TX2
    REAL(kind=r8) :: TX3
    REAL(kind=r8) :: TX4
    REAL(kind=r8) :: TX5
    REAL(kind=r8) :: QIS
    REAL(kind=r8) :: QLS
    REAL(kind=r8) :: HBL
    REAL(kind=r8) :: QBL
    REAL(kind=r8) :: RBL(M)
    REAL(kind=r8) :: QLB
    REAL(kind=r8) :: QIB
    REAL(kind=r8) :: PRIS
    !REAL(kind=r8) :: TX6
    REAL(kind=r8) :: ACR
    !REAL(kind=r8) :: TX7
    !EAL(kind=r8) :: TX8
    !REAL(kind=r8) :: TX9
    REAL(kind=r8) :: RHC
    REAL(kind=r8) :: hstkd
    REAL(kind=r8) :: qstkd
    REAL(kind=r8) :: ltlkd
    REAL(kind=r8) :: q0ukd
    REAL(kind=r8) :: q0dkd
    REAL(kind=r8) :: dlbkd
    REAL(kind=r8) :: qtp
    REAL(kind=r8) :: qw00
    REAL(kind=r8) :: qi00
    REAL(kind=r8) :: qrbkd
    REAL(kind=r8) :: hstold
    REAL(kind=r8) :: rel_fac

    !     INTEGER IA,  I1,  I2, ID1, ID2
    !     INTEGER IB,  I3

    LOGICAL :: UNSAT
    LOGICAL ::  ep_wfn

    LOGICAL :: LOWEST
    !LOGICAL ::  SKPDD

    REAL(kind=r8) :: TL
    REAL(kind=r8) :: PL
    REAL(kind=r8) :: QL
    REAL(kind=r8) :: QS
    REAL(kind=r8) :: DQS
    REAL(kind=r8) :: ST1
    !REAL(kind=r8) :: SGN
    REAL(kind=r8) :: TAU
    REAL(kind=r8) :: QTVP
    REAL(kind=r8) :: HB
    REAL(kind=r8) :: QB
    !REAL(kind=r8) :: TB
    !REAL(kind=r8) :: QQQ
    REAL(kind=r8) :: HCCP
    REAL(kind=r8) :: DS
    REAL(kind=r8) :: DH
    REAL(kind=r8) :: AMBMAX
    REAL(kind=r8) :: X00
    REAL(kind=r8) :: EPP
    REAL(kind=r8) :: QTLP
    REAL(kind=r8) :: DPI
    REAL(kind=r8) :: DPHIB
    REAL(kind=r8) :: DPHIT
    REAL(kind=r8) :: DEL_ETA
    REAL(kind=r8) :: DETP
    REAL(kind=r8) :: TEM
    REAL(kind=r8) :: TEM1
    REAL(kind=r8) :: TEM2
    REAL(kind=r8) :: TEM3
    REAL(kind=r8) :: TEM4
    REAL(kind=r8) :: ST2
    REAL(kind=r8) :: ST3
    REAL(kind=r8) :: ST4
    REAL(kind=r8) :: ST5
    !REAL(kind=r8) :: ERRH
    !REAL(kind=r8) :: ERRW
    !REAL(kind=r8) :: ERRE
    REAL(kind=r8) :: TEM5
    REAL(kind=r8) :: TEM6
    REAL(kind=r8) :: HBD
    REAL(kind=r8) :: QBD
    REAL(kind=r8) :: st1s
    !REAL(kind=r8), PARAMETER :: ERRMIN=0.0001_r8!, ERRMI2=0.1_r8*ERRMIN

    !     parameter (c0=1.0e-3_r8, KBLMX=20, ERRMIN=0.0001_r8, ERRMI2=0.1_r8*ERRMIN)

    !INTEGER    :: I
    INTEGER    :: L
    INTEGER    :: N
    INTEGER    :: KD1
    INTEGER    :: II 
    INTEGER    :: KP1
    !INTEGER    :: IT
    INTEGER    :: KM1
    INTEGER    :: KTEM
    !INTEGER    :: KK
    !INTEGER    :: KK1
    !INTEGER    :: LM1
    !INTEGER    :: LL
    !INTEGER    :: LP1
    INTEGER    :: kbls
    INTEGER    :: kmxh

    REAL(kind=r8) ::  avt
    REAL(kind=r8) ::  avq
    REAL(kind=r8) ::  avr
    REAL(kind=r8) ::  avh
    !
    !     REEVAPORATION
    !
    !     real(kind=r8), parameter ::
    !    &                   clfa = -0.452550814376093547E-03_r8
    !    &,                  clfb =  0.161398573159240791E-01_r8
    !    &,                  clfc = -0.163676268676807096_r8
    !    &,                  clfd =  0.447988962175259131_r8
    !    &,                  point3 = 0.3, point01=0.01_r8

    !     real(kind=r8), parameter :: rainmin=1.0e-9_r8
    REAL(kind=r8), PARAMETER :: rainmin=1.0e-8_r8
    !REAL(kind=r8), PARAMETER :: oneopt9=1.0_r8/0.09_r8
    !REAL(kind=r8), PARAMETER :: oneopt4=1.0_r8/0.04_r8

    REAL(kind=r8) :: CLFRAC

    REAL(kind=r8) :: ACTEVAP
    !REAL(kind=r8) :: AREARAT
    REAL(kind=r8) :: DELTAQ
    !REAL(kind=r8) :: MASS
    !REAL(kind=r8) :: MASSINV
    REAL(kind=r8) :: POTEVAP  
    REAL(kind=r8) :: TEQ
    REAL(kind=r8) :: QSTEQ
    REAL(kind=r8) :: DQDT
    REAL(kind=r8) :: QEQ
    !
    !     Temporary workspace and parameters needed for downdraft
    !
    !REAL(kind=r8) :: GMF
    !
    REAL(kind=r8) :: BUY(KD:K+1)
    REAL(kind=r8) :: QRB(KD:K)
    REAL(kind=r8) :: QRT(KD:K) 
    REAL(kind=r8) :: ETD(KD:K+1)
    REAL(kind=r8) :: HOD(KD:K+1)
    REAL(kind=r8) :: QOD(KD:K+1) 
    REAL(kind=r8) :: GHD(KD:K)
    REAL(kind=r8) :: GSD(KD:K)
    REAL(kind=r8) :: EVP(KD:K)    
    REAL(kind=r8) :: ETZ(KD:K)
    REAL(kind=r8) :: CLDFR(KD:K)
    REAL(kind=r8) :: TRAIN
    REAL(kind=r8) :: DOF
    REAL(kind=r8) :: CLDFRD
    !REAL(kind=r8) :: FAC
    !REAL(kind=r8) :: RSUM1
    !REAL(kind=r8) :: RSUM2
    !REAL(kind=r8) :: RSUM3
    REAL(kind=r8) :: dpneg
    INTEGER       :: IDH
    LOGICAL       :: DDFT
    !     real(kind=r8) eps, epsm1, rvi, facw, faci, hsub, tmix, DEN
    !     real(kind=r8) eps, epsm1, rv, rd, depth
    !     real(kind=r8) eps, epsm1, rv, rd, fpvs, depth
    !
    !
    !***********************************************************************
    !
    !  UPDATE ARGUMENTS

    TCD=0.0_r8;QCD=0.0_r8

    !  TEMPORARY WORK SPACE

    HOL=0.0_r8
    QOL=0.0_r8
    GAF=0.0_r8
    HST=0.0_r8
    QST=0.0_r8
    TOL=0.0_r8
    GMH=0.0_r8
    GMS=0.0_r8
    GAM=0.0_r8
    AKT=0.0_r8
    AKC=0.0_r8
    BKC=0.0_r8
    LTL=0.0_r8
    RNN=0.0_r8
    FCO=0.0_r8
    PRI=0.0_r8
    !REAL(kind=r8) :: PRH(KD:K)
    QIL=0.0_r8
    QLL=0.0_r8
    ZET=0.0_r8
    XI =0.0_r8
    RNS=0.0_r8
    Q0U=0.0_r8
    Q0D=0.0_r8
    vtf=0.0_r8
    DLB=0.0_r8
    DLT=0.0_r8
    ETA=0.0_r8
    PRL=0.0_r8
    CIL=0.0_r8
    CLL=0.0_r8
    ETAI=0.0_r8

    ALM=0.0_r8
    DET=0.0_r8
    HCC=0.0_r8
    CLP=0.0_r8
    HSU=0.0_r8
    HSD=0.0_r8
    QTL=0.0_r8
    QTV=0.0_r8
    AKM=0.0_r8
    WFN=0.0_r8
    HOS=0.0_r8
    QOS=0.0_r8
    AMB=0.0_r8
    TX1=0.0_r8
    TX2=0.0_r8
    TX3=0.0_r8
    TX4=0.0_r8
    TX5=0.0_r8
    QIS=0.0_r8
    QLS=0.0_r8
    HBL=0.0_r8
    QBL=0.0_r8
    RBL=0.0_r8
    QLB=0.0_r8
    QIB=0.0_r8
    PRIS=0.0_r8
    ACR=0.0_r8
    RHC=0.0_r8
    hstkd=0.0_r8
    qstkd=0.0_r8
    ltlkd=0.0_r8
    q0ukd=0.0_r8
    q0dkd=0.0_r8
    dlbkd=0.0_r8
    qtp=0.0_r8
    qw00=0.0_r8
    qi00=0.0_r8
    qrbkd=0.0_r8
    hstold=0.0_r8
    rel_fac=0.0_r8

    !     INTEGER IA,  I1,  I2, ID1, ID2
    !     INTEGER IB,  I3


    !LOGICAL ::  SKPDD

    TL=0.0_r8
    PL=0.0_r8
    QL=0.0_r8
    QS=0.0_r8
    DQS=0.0_r8
    ST1=0.0_r8
    TAU=0.0_r8
    QTVP=0.0_r8
    HB=0.0_r8
    QB=0.0_r8
    HCCP=0.0_r8
    DS=0.0_r8
    DH=0.0_r8
    AMBMAX=0.0_r8
    X00=0.0_r8
    EPP=0.0_r8
    QTLP=0.0_r8
    DPI=0.0_r8
    DPHIB=0.0_r8
    DPHIT=0.0_r8
    DEL_ETA=0.0_r8
    DETP=0.0_r8
    TEM=0.0_r8
    TEM1=0.0_r8
    TEM2=0.0_r8
    TEM3=0.0_r8
    TEM4=0.0_r8
    ST2=0.0_r8
    ST3=0.0_r8
    ST4=0.0_r8
    ST5=0.0_r8
    TEM5=0.0_r8
    TEM6=0.0_r8
    HBD=0.0_r8
    QBD=0.0_r8
    st1s=0.0_r8

    L=0
    N=0
    KD1=0
    II =0
    KP1=0
    KM1=0
    KTEM=0
    kbls=0
    kmxh=0

    avt=0.0_r8
    avq=0.0_r8
    avr=0.0_r8
    avh=0.0_r8
    CLFRAC=0.0_r8
    ACTEVAP=0.0_r8
    DELTAQ=0.0_r8
    POTEVAP =0.0_r8 
    TEQ=0.0_r8
    QSTEQ=0.0_r8
    DQDT=0.0_r8
    QEQ=0.0_r8
    !
    !     Temporary workspace and parameters needed for downdraft
    !
    !
    BUY=0.0_r8
    QRB=0.0_r8
    QRT =0.0_r8
    ETD=0.0_r8
    HOD=0.0_r8
    QOD =0.0_r8
    GHD=0.0_r8
    GSD=0.0_r8
    EVP =0.0_r8
    ETZ=0.0_r8
    CLDFR=0.0_r8
    TRAIN=0.0_r8
    DOF=0.0_r8
    CLDFRD=0.0_r8
    dpneg=0.0_r8
    IDH=0

    !
    DO l=1,K
       tcd(L) = 0.0_r8
       qcd(L) = 0.0_r8
    ENDDO
    !
    KP1     = K  + 1
    KM1     = K  - 1
    KD1     = KD + 1
    kblmx   = k / 2
    !
    !     if (lprnt) print *,' IN CLOUD for KD=',kd
    !     if (lprnt) print *,' prs=',prs(Kd:K+1)
    !     if (lprnt) print *,' phil=',phil(KD:K)
    !     if (lprnt) print *,' phih=',phih(KD:K+1)
    !     if (lprnt) print *,' toi=',toi
    !     if (lprnt) print *,' qoi=',qoi
    !
    !     do l=kd1,k
    !       alfint(l) = (prjm(l)-prj(l)) / (prjm(l)-prjm(l-1))
    !       alfinq(l) = alfint(l)
    !     enddo
    !
    CLDFRD   = 0.0_r8
    DOF      = 0.0_r8
    PRL(KP1) = PRS(KP1)
    !
    DO L=KD,K
       RNN(L) = 0.0_r8
       ZET(L) = 0.0_r8
       XI(L)  = 0.0_r8
       !
       TOL(L) = TOI(L)
       QOL(L) = QOI(L)
       PRL(L) = PRS(L)
       BUY(L) = 0.0_r8
       CLL(L) = QLI(L)
       CIL(L) = QII(L)
    ENDDO
    !
    DO L=KD, K
       DPI    = ONE / (PRL(L+1) - PRL(L))
       PRI(L) = GRAVFAC * DPI
       !
       PL     = PRSM(L)
       TL     = TOL(L)

       !     if (lprnt) print *,' l=',l,' prl=',prl(l+1),prl(l),' pl=',pl,
       !    &' dpi=',dpi,' prsm=',prsm(l)

       AKT(L) = (PRL(L+1) - PL) * DPI
       !
       !     if (lprnt) print *,' l=',l,' prl=',prl(l+1),prl(l),' pl=',pl,
       !    &' dpi=',dpi,' prsm=',prsm(l),' akt=',akt(l)
       !
       CALL QSATCN(TL, PL, QS, DQS)
       !
       !     if(lprnt)print*,' qs=',qs,' tl=',tl,' pl=',pl
       !    1,               ' dqs=',dqs,' qol=',qol(l)
       !
       !
       QST(L) = QS
       GAM(L) = DQS * ELOCP
       ST1    = ONE + GAM(L)
       GAF(L) = (ONE/con_hvap ) * (GAM(L)/(ONE + GAM(L)))

       QL     = MAX(MIN(QS*RHMAX,QOL(L)), ONE_M10)
       QOL(L) = QL

       TEM    = con_cp * TL
       LTL(L) = TEM * ST1 / (ONE+con_FVirt*(QST(L)+TL*DQS))
       vtf(L) = 1.0_r8 + con_FVirt * QL
       ETA(L) = ONE / (LTL(L) * VTF(L))

       HOL(L) = TEM + QL * con_hvap 
       HST(L) = TEM + QS * con_hvap 
       !
       !     if(lprnt)print*,' l=',l,' hst=',hst(l),' tem=',tem
       !    1,               ' qs=',qs,' con_hvap =',con_hvap 
       !     if (lprnt) print *,' L=',L,' tem=',tem,' ql=',ql,' con_hvap =',con_hvap 
       !    &,' con_hfus=',con_hfus,' qii=',qii(l),' con_cp=',con_cp,' tl=',tl
       !    &,' qs=',qs,' qol=',qol(l),' rhmax=',rhmax,' hol=',hol(l)
       !    &,' pl=',pl

    ENDDO
    !
    ETA(K+1) = ZERO
    GMS(K)   = ZERO
    !
    AKT(KD)  = HALF
    GMS(KD)  = ZERO
    !
    CLP      = ZERO
    !
    GAM(K+1) = GAM(K)
    GAF(K+1) = GAF(K)
    !
    DO L=K,KD1,-1
       !       TEM1   = con_cp * TOL(L) * VTF(L) / PRH(L)

       DPHIB  = PHIL(L) - PHIH(L+1)
       DPHIT  = PHIH(L) - PHIL(L)
       !
       DLB(L) = DPHIB * ETA(L)
       DLT(L) = DPHIT * ETA(L)
       !
       QRB(L) = DPHIB
       QRT(L) = DPHIT
       !
       ETA(L) = ETA(L+1) + DPHIB

       !     if (lprnt) print *,' L=',L,' dphib=',dphib,' dphit=',dphit
       !    &,' eta=',eta(l),' hol_new=',hol(l)+eta(l)
       !    &,' con_cp=',con_cp,' tol=',tol(l),' vtf=',vtf(l)
       !
       HOL(L) = HOL(L) + ETA(L)
       hstold = hst(l)
       HST(L) = HST(L) + ETA(L)
       !
       !     if(lprnt)print*,' l=',l,' hst=',hst(l),' eta=',eta(l)
       !    1,               ' hstold=',hstold

       ETA(L) = ETA(L) + DPHIT
    ENDDO
    !
    !     For the cloud top layer
    !
    L = KD

    DPHIB  = PHIL(L) - PHIH(L+1)
    !
    DLB(L) = DPHIB * ETA(L)
    !
    QRB(L) = DPHIB
    QRT(L) = DPHIB
    !
    ETA(L) = ETA(L+1) + DPHIB

    HOL(L) = HOL(L) + ETA(L)
    HST(L) = HST(L) + ETA(L)
    !
    !     if (kd .eq. 12) then
    !     if (lprnt) print *,' IN CLOUD for KD=',KD,' K=',K
    !     if (lprnt) print *,' l=',l,' hol=',hol(l),' hst=',hst(l)
    !     if (lprnt) print *,' TOL=',tol
    !     if (lprnt) print *,' qol=',qol
    !     if (lprnt) print *,' hol=',hol
    !     if (lprnt) print *,' hst=',hst
    !     endif
    !
    !     To determine KBL internally -- If KBL is defined externally
    !     the following two loop should be skipped
    !
    !     if (lprnt) print *,' calkbl=',calkbl

    IF (CALKBL) THEN
       KTEM = MAX(KD, K-KBLMX-2)
       kmxh = k

       !        DO L=KM1,KTEM,-1
       !     if(lprnt) print *,' l=',l,' kmxh=',kmxh,' prl=',prl(l)
       !    &, prl(k),' hol=',hol(l),hol(kmxh)
       !          if (prl(k) - prl(l) .gt. 100.0_r8) exit
       !          if (hol(l) .gt. hol(kmxh)) kmxh = l
       !     if(lprnt) print *,' l=',l,' kmxh=',kmxh,' prl=',prl(l)
       !        ENDDO

       DO L=kmxh,KTEM+1,-1
          kbls = l
          IF (hst(l-1) .GT. hst(l)) EXIT
       ENDDO
       KBL   = Kmxh
       TX1   = ZERO
       UNSAT = .FALSE.
       DO L=kmxh-1,KTEM,-1
          TEM = HOL(K) - HOL(L)
          TX3 = (HOL(L) - HOL(L+1)) / (PRL(L+2) - PRL(L))

          !     if (lprnt) print *,' l=',l,' kbl=',kbl,' tx3=',tx3,' tx1=',tx1
          IF (TX3 .LT. TX1 .AND. TEM .LT. HCRIT) THEN
             TX1   = TX3
             KBL   = L
             !            KBL   = L+1
             UNSAT = .TRUE.
          ELSEIF (UNSAT .AND.                                          &
               &           ( ((KBL .LT. K-1) .AND. TX3 .GT. 0.5_r8*TX1)              &
               &              .OR. TEM .GT. HCRIT) ) THEN
             TX1 = -1.0E20_r8
          ENDIF
       ENDDO
       !     if(lprnt) print *,' kbl=',kbl,' kbls=',kbls,' kmxh=',kmxh
       !
       !        ii = min(kbl,kbls)
       ii = kbl
       DO l=ktem,kmxh-1
          !          if (hol(kmxh) .gt. hst(l)) kbl = l
          IF (hol(kmxh) .GT. hst(l)) kbl = l+1    ! Commented on 09/20/04
       ENDDO
       !        if(lprnt) print *,' kblhst=',kbl,' ii=',ii

       IF (prl(K+1) - prl(ii) .GT. 50.0_r8 .AND. ii .GT. kbl) kbl = ii
       !        if(lprnt) print *,' kbl2=',kbl,' ii=',ii
       IF (kbl .NE. ii) THEN
          !          kbl = min(K, max(kbl+1, kd-1))
!!!        kbl = min(K, max(kbl, kd-1))
          IF (PRL(K+1)-PRL(KBL) .GT. bldmax) kbl = MAX(kbl,ii)
       ENDIF
       !        if (ii .gt. kbl) then
       !           if (hol(Kmxh)-hol(kbl) .gt. hcrit) kbl = ii
       !        endif
       !
       !        ii = kbl
       !        do l=ii,k
       !          if (hol(k) .gt. hst(l)) kbl = l
       !        enddo
!!!      kbl = min(K, max(kbl, kd-1))
       !
       KBL  = MIN(k, MAX(KBL,K-KBLMX))
!!!!!    kbl = K - 2
!!!
       !        tem1 = max(10.0_r8, min(50.0_r8,(prl(k+1) - prl(kd))*0.05_r8))
       !LL2     tem1 = max(10.0_r8, min(50.0_r8,(prl(k+1) - prl(kd))*0.066_r8))
       !        do l=k,k-kblmx,-1
       !          tem = prl(k+1) - prl(l)
       !LL        if (tem .gt. 20.0_r8) then
       !LL        if (tem .gt. 40.0_r8) then
       !!r1       if (tem .gt. 50.0_r8) then
       !          if (tem .gt. tem1) then
       !            kbl = min(kbl,l)
       !            exit
       !          endif
       !        enddo
       ! 
       tem1 = MAX(prl(k+1)-prl(k),                                    &
            &                     MIN((prl(kbl) - prl(kd))*0.05_r8, 20.0_r8))
       !    &                     min((prl(kbl) - prl(kd))*0.05_r8, 30.0_r8))
       IF (prl(k+1)-prl(kbl) .LT. tem1) THEN
          KTEM = MAX(KD+1, K-KBLMX)
          DO l=k,KTEM,-1
             tem = prl(k+1) - prl(l)
             IF (tem .GT. tem1) THEN
                kbl = MIN(kbl,l)
                EXIT
             ENDIF
          ENDDO
       ENDIF
!!!

       KPBL = KBL
       !     if(lprnt)print*,' 1st kbl=',kbl,' kblmx=',kblmx,' kd=',kd
       !     if(lprnt)print*,' tx3=',tx3,' tx1=',tx1,' tem=',tem
       !    1,               ' hcrit=',hcrit
    ELSE
       KBL  = KPBL
       !     if(lprnt)print*,' 2nd kbl=',kbl
    ENDIF
    !     if(lprnt)print*,' after CALKBL l=',l,' hol=',hol(l)
    !    1,               ' hst=',hst(l)
    !
    KBL      = MAX(KBL,KD)
    KB1      = KBL - 1
    !!
    !!
    !     do kbl=k,kd1,-1
    !       st1  = 1.0_r8 / (PRL(K+1) - PRL(KBL))
    !       tem1 = (PRL(K+1)-PRL(K)) * st1
    !       HBL = HOL(K) * tem1
    !
    !       DO L=KM1,KBL,-1
    !         tem2 = (PRL(K+1)-PRL(L)) * st1
    !         TEM  = tem2 - tem1
    !         HBL  = HBL + HOL(L) * TEM
    !         tem1 = tem2
    !       enddo
    !     if(lprnt) print *,' HBL=',HBL,' KBL=',KBL
    !       KB1 = KBL - 1
    !       st2 = 0.5_r8 * (hst(kbl)+hst(kb1))
    !       if (st2 .le. hbl) exit
    !     ENDDO
    !     if (hst(kbl) .le. hbl) kbl = kbl + 1
    !     kbl = min(kbl+1, K)
    !     KB1 = KBL - 1
    !!    if (lprnt) print *,' HBL=',HBL,' HST=',st2,' KBL=',KBL,' kb1=',kb1
    !     if(lprnt)print *,' HBL=',HBL,' HST=',hst(l),' KBL=',KBL
    !    1,                ' kb1=',kb1
    !     if (PRL(K+1)-PRL(KBL) .gt. bldmax) return
    !!
    !!
    !
    !     if (lprnt) print *,' kbl=',kbl,' prlkbl=',prl(kbl),prl(k+1)
    IF(kb1 .LE. kd)THEN
       !       if(lprnt)print*,' kb1=',kb1,' kd=',kd,' EXIT CLOUD'
       RETURN
    ENDIF
    IF(PRL(K+1)-PRL(KBL) .GT. bldmax)THEN
       !       if(lprnt)print*,' prl(k+1)=',prl(k+1),' prl(kbl)=',prl(kbl)
       !    1,               ' bldmax=',bldmax,' k+1=',k+1,' kbl=',kbl
       !    2,               ' EXIT CLOUD'
       RETURN
    ENDIF
    !
    !     if (lprnt) print *,' kbl=',kbl
    !
    PRIS     = ONE / (PRL(K+1)-PRL(KBL))
    TX1      = ETA(KBL)
    !
    GMS(KBL) = 0.0_r8
    XI(KBL)  = 0.0_r8
    ZET(KBL) = 0.0_r8
    !     DEPTH    = ETA(KD) - ETA(KBL)
    !
    DO L=K,KD,-1
       IF (L .GE. KBL) THEN
          ETA(L) = (PRL(K+1)-PRL(L)) * PRIS
       ELSE
          ZET(L) = (ETA(L) - TX1) * ONEBG
          XI(L)  =  ZET(L) * ZET(L) * QUDFAC
          ETA(L) =  ZET(L) - ZET(L+1)
          GMS(L) =  XI(L)  - XI(L+1)
       ENDIF
    ENDDO
    !
    HBL = HOL(K) * ETA(K)
    QBL = QOL(K) * ETA(K)
    QLB = CLL(K) * ETA(K)
    QIB = CIL(K) * ETA(K)
    !     TX1 = QOL(K) / QST(K) * ETA(K)
    TX1 = QST(K) * ETA(K)
    !
    DO L=KM1,KBL,-1
       TEM = ETA(L) - ETA(L+1)
       HBL = HBL + HOL(L) * TEM
       !     if(lprnt)print*,' l=',l,' qbl=',qbl,' qol=',qol(l)
       !    1,               ' tem=',tem
       QBL = QBL + QOL(L) * TEM
       QLB = QLB + CLL(L) * TEM
       QIB = QIB + CIL(L) * TEM
       !        TX1 = TX1 + QOL(L) / QST(L) * TEM
       TX1 = TX1 + QST(L) * TEM
    ENDDO
    !     if (lprnt) print *,' hbl=',hbl,' qbl=',qbl
    !                                   Find Min value of HOL in TX2
    TX2 = HOL(KD)
    IDH = KD1
    DO L=KD1,KB1
       IF (HOL(L) .LT. TX2) THEN
          TX2 = HOL(L)
          IDH = L             ! Level of minimum moist static energy!
       ENDIF
    ENDDO
    IDH = 1
    IDH = MAX(KD1, IDH)
    !
    TEM1 = HBL - HOL(KD)
    TEM  = HBL - HST(KD1)                                             &
         &             - LTL(KD1) *( con_FVirt *(QOL(KD1)-QST(KD1)))
    LOWEST = KD .EQ. KB1
    !

    !     TX1   = QBL / TX1
    TX1   = RHFACS - QBL / TX1       !     Average RH
    !     TX1   = RHFACS - TX1             !     Average of each layer RH
    UNSAT = (TEM .GT. ZERO .OR. (LOWEST .AND. TEM1 .GE. ZERO))        &
         &         .AND. (TX1 .LT. RHRAM)                                   &
         &         .AND. (KBL .GT. KD)

    !     if(lprnt) print *,' unsat=',unsat,' tem=',tem,' tem1=',tem1
    !    &,' tx1=',tx1,' rhram=',rhram,' kbl=',kbl,' kd=',kd,' lowest='
    !    &,lowest,' rhfacs=',rhfacs,' ltl=',ltl(kd1),' qol=',qol(kd1)
    !    &,' qst=',qst(kd1),' hst=',hst(kd1),' con_FVirt=',con_FVirt

    !
    !===>  IF NO SOUNDING MEETS FIRST CONDITION, RETURN
    !     if(lprnt .and. (.not. unsat)) print *,' tx1=',tx1,' rhfacs='
    !    &,rhfacs, ' tem=',tem,' hst=',hst(kd1)

    IF (.NOT. UNSAT) RETURN
    !
    !     TEM1   = TX1 - RHFACS
    !     RHC    = MAX(ZERO, MIN(ONE, EXP(20.0_r8*TEM1) ))
    RHC    = MAX(ZERO, MIN(ONE, EXP(-20.0_r8*TX1) ))
    !
    DO N=1,M
       RBL(N) = ROI(K,N) * ETA(K)
    ENDDO
    DO N=1,M
       DO L=KM1,KBL,-1
          RBL(N) = RBL(N) + ROI(L,N)*(ETA(L)-ETA(L+1))
       ENDDO
    ENDDO
    !
    TX4    = 0.0_r8
    TX5    = 0.0_r8
    !
    TX3      = QST(KBL) - GAF(KBL) * HST(KBL)
    QIL(KBL) = MAX(ZERO, MIN(ONE, (TCR-TCL-TOL(KBL))*TCRF))
    !
    DO L=KB1,KD1,-1
       TEM      = QST(L) - GAF(L) * HST(L)
       TEM1     = (TX3 + TEM) * 0.5_r8
       ST2      = (GAF(L)+GAF(L+1)) * 0.5_r8
       !
       FCO(L+1) =            TEM1 + ST2 * HBL

       !     if(lprnt) print *,' fco=',fco(l+1),' tem1=',tem1,' st2=',st2
       !    &,' hbl=',hbl,' tx3=',tx3,' tem=',tem,' gaf=',gaf(l),' l=',l

       RNN(L+1) = ZET(L+1) * TEM1 + ST2 * TX4
       GMH(L+1) = XI(L+1)  * TEM1 + ST2 * TX5
       !
       TX3      = TEM
       TX4      = TX4 + ETA(L) * HOL(L)
       TX5      = TX5 + GMS(L) * HOL(L)
       !
       QIL(L)   = MAX(ZERO, MIN(ONE, (TCR-TCL-TOL(L))*TCRF))
       QLL(L+1) = (0.5_r8*con_hfus) * ST2 * (QIL(L)+QIL(L+1)) + ONE
    ENDDO
    !
    !     FOR THE CLOUD TOP -- L=KD
    !
    L = KD
    !
    TEM      = QST(L) - GAF(L) * HST(L)
    TEM1     = (TX3 + TEM) * 0.5_r8
    ST2      = (GAF(L)+GAF(L+1)) * 0.5_r8
    !
    FCO(L+1) =            TEM1 + ST2 * HBL
    RNN(L+1) = ZET(L+1) * TEM1 + ST2 * TX4
    GMH(L+1) = XI(L+1)  * TEM1 + ST2 * TX5
    !
    FCO(L)   = TEM + GAF(L) * HBL
    RNN(L)   = TEM * ZET(L) + (TX4 + ETA(L)*HOL(L)) * GAF(L)
    GMH(L)   = TEM * XI(L)  + (TX5 + GMS(L)*HOL(L)) * GAF(L)
    !
    !   Replace FCO for the Bottom
    !
    FCO(KBL) = QBL
    RNN(KBL) = 0.0_r8
    GMH(KBL) = 0.0_r8
    !
    QIL(KD)  =  MAX(ZERO, MIN(ONE, (TCR-TCL-TOL(KD))*TCRF))
    QLL(KD1) = (0.5_r8*con_hfus) * ST2 * (QIL(KD) + QIL(KD1)) + ONE
    QLL(KD ) = con_hfus * GAF(KD) * QIL(KD) + ONE
    !
    !     if (lprnt) print *,' fco=',fco(kd:kbl)
    !     if (lprnt) print *,' qil=',qil(kd:kbl)
    !     if (lprnt) print *,' qll=',qll(kd:kbl)
    !
    st1  = qil(kd)
    st2  = c0i * st1
    tem  = c0  * (1.0_r8-st1)
    tem2 = st2*qi0 + tem*qw0
    !
    DO L=KD,KB1
       tx2    = akt(l) * eta(l)
       tx1    = tx2 * tem2
       q0u(l) = tx1
       FCO(L) = FCO(L+1) - FCO(L) + tx1
       RNN(L) = RNN(L+1) - RNN(L)                                     &
            &          + ETA(L)*(QOL(L)+CLL(L)+CIL(L)) + tx1*zet(l)
       GMH(L) = GMH(L+1) - GMH(L)                                     &
            &          + GMS(L)*(QOL(L)+CLL(L)+CIL(L)) + tx1*xi(l)
       !
       tem1   = (1.0_r8-akt(l)) * eta(l)

       !     if(lprnt) print *,' qll=',qll(l),' st2=',st2,' tem=',tem
       !    &,' tx2=',tx2,' akt=',akt(l),' eta=',eta(l)

       AKT(L) = QLL(L)   + (st2 + tem) * tx2

       !     if(lprnt) print *,' akt==',akt(l),' l==',l

       AKC(L) = 1.0_r8 / AKT(L)
       !
       st1    = 0.5_r8 * (qil(l)+qil(l+1))
       st2    = c0i * st1
       tem    = c0  * (1.0_r8-st1)
       tem2   = st2*qi0 + tem*qw0
       !
       BKC(L) = QLL(L+1) - (st2 + tem) * tem1
       !
       tx1    = tem1*tem2
       q0d(l) = tx1
       FCO(L) = FCO(L) + tx1
       RNN(L) = RNN(L) + tx1*zet(l+1)
       GMH(L) = GMH(L) + tx1*xi(l+1)
    ENDDO

    !     if(lprnt) print *,' akt=',akt(kd:kb1)
    !     if(lprnt) print *,' akc=',akc(kd:kb1)

    qw00 = qw0
    qi00 = qi0
    ii = 0
777 CONTINUE
    !
    !     if (lprnt) print *,' after 777 ii=',ii,' ep_wfn=',ep_wfn
    !
    ep_wfn = .FALSE.
    RNN(KBL) = 0.0_r8
    TX3      = bkc(kb1) * (QIB + QLB)
    TX4      = 0.0_r8
    TX5      = 0.0_r8
    DO L=KB1,KD1,-1
       TEM    = BKC(L-1)       * AKC(L)
       !     if (lprnt) print *,' tx3=',tx3,' fco=',fco(l),' akc=',akc(l)
       !    &,' bkc=',bkc(l-1), ' l=',l
       TX3    = (TX3 + FCO(L)) * TEM
       TX4    = (TX4 + RNN(L)) * TEM
       TX5    = (TX5 + GMH(L)) * TEM
    ENDDO
    IF (KD .LT. KB1) THEN
       HSD   = HST(KD1)                                               &
            &         + LTL(KD1) *  con_FVirt *(QOL(KD1)-QST(KD1))
    ELSE
       HSD   = HBL
    ENDIF
    !
    !     if (lprnt) print *,' tx3=',tx3,' fco=',fco(kd),' akc=',akc(kd)
    TX3 = (TX3 + FCO(KD)) * AKC(KD)
    TX4 = (TX4 + RNN(KD)) * AKC(KD)
    TX5 = (TX5 + GMH(KD)) * AKC(KD)
    ALM = con_hfus*QIL(KD) - LTL(KD) * VTF(KD)
    !
    HSU = HST(KD) + LTL(KD) * con_FVirt * (QOL(KD)-QST(KD))

    !     if (lprnt) print *,' hsu=',hsu,' hst=',hst(kd),
    !    &' ltl=',ltl(kd),' qol=',qol(kd),' qst=',qst(kd)
    !
    !===> VERTICAL INTEGRALS NEEDED TO COMPUTE THE ENTRAINMENT PARAMETER
    !
    TX1 = ALM * TX4
    TX2 = ALM * TX5

    DO L=KD,KB1
       TAU = HOL(L) - HSU
       TX1 = TX1 + TAU * ETA(L)
       TX2 = TX2 + TAU * GMS(L)
    ENDDO
    !
    !     MODIFY HSU TO INCLUDE CLOUD LIQUID WATER AND ICE TERMS
    !
    !     if (lprnt) print *,' hsu=',hsu,' alm=',alm,' tx3=',tx3

    HSU   = HSU - ALM * TX3
    !
    CLP   = ZERO
    ALM   = -100.0_r8
    HOS   = HOL(KD)
    QOS   = QOL(KD)
    QIS   = CIL(KD)
    QLS   = CLL(KD)
    UNSAT = HBL .GT. HSU .AND. ABS(tx1) .GT. 1.0e-4_r8

    !     if (lprnt) print *,' ii=',ii,' unsat=',unsat,' hsu=',hsu
    !    &,' hbl=',hbl,' tx1=',tx1,' hsd=',hsd


    !***********************************************************************


    ST1  = HALF*(HSU + HSD)
    IF (UNSAT) THEN
       !
       !  STANDARD CASE:
       !   CLOUD CAN BE NEUTRALLY BOUYANT AT MIDDLE OF LEVEL KD W/ +VE LAMBDA.
       !   EPP < .25 IS REQUIRED TO HAVE REAL ROOTS.
       !
       clp = 1.0_r8
       st2 = hbl - hsu

       !     if(lprnt) print *,' tx2=',tx2,' tx1=',tx1,' st2=',st2
       !
       IF (tx2 .EQ. 0.0_r8) THEN
          alm = - st2 / tx1
          IF (alm .GT. almax) alm = -100.0_r8
       ELSE
          x00 = tx2 + tx2
          epp = tx1 * tx1 - (x00+x00)*st2
          IF (epp .GT. 0.0_r8) THEN
             x00  = 1.0_r8 / x00
             tem  = SQRT(epp)
             tem1 = (-tx1-tem)*x00
             tem2 = (-tx1+tem)*x00
             IF (tem1 .GT. almax) tem1 = -100.0_r8
             IF (tem2 .GT. almax) tem2 = -100.0_r8
             alm  = MAX(tem1,tem2)

             !     if (lprnt) print *,' tem1=',tem1,' tem2=',tem2,' alm=',alm
             !    &,' tx1=',tx1,' tem=',tem,' epp=',epp,' x00=',x00,' st2=',st2

          ENDIF
       ENDIF

       !     if (lprnt) print *,' almF=',alm,' ii=',ii,' qw00=',qw00
       !    &,' qi00=',qi00
       !
       !  CLIP CASE:
       !   NON-ENTRAINIG CLOUD DETRAINS IN LOWER HALF OF TOP LAYER.
       !   NO CLOUDS ARE ALLOWED TO DETRAIN BELOW THE TOP LAYER.
       !
    ELSEIF ( (HBL .LE. HSU) .AND.                                    &
         &          (HBL .GT. ST1   )     ) THEN
       ALM = ZERO
       CLP = (HBL-ST1) / (HSU-ST1)
    ENDIF
    !
    UNSAT = .TRUE.
    IF (ALMIN1 .GT. 0.0_r8) THEN
       IF (ALM .GE. ALMIN1) UNSAT = .FALSE.
    ELSE
       LOWEST   = KD .EQ. KB1
       IF ( (ALM .GT. ZERO) .OR.                                       &
            &      (.NOT. LOWEST .AND. ALM .EQ. ZERO) ) UNSAT = .FALSE.
    ENDIF
    !
    !     if (alm*depth/con_g .ge. 1.0) UNSAT = .TRUE.
    !
    !===>  IF NO SOUNDING MEETS SECOND CONDITION, RETURN
    !
    IF (UNSAT) THEN
       IF (ii .GT. 0 .OR. (qw00 .EQ. 0.0_r8 .AND. qi00 .EQ. 0.0_r8)) RETURN
       CLP = 1.0_r8
       ep_wfn = .TRUE.
       GO TO 888
    ENDIF
    !
    !     if (lprnt) print *,' hstkd=',hst(kd),' qstkd=',qst(kd)
    !    &,' ii=',ii,' clp=',clp

    st1s = ONE
    IF(CLP.GT.ZERO .AND. CLP.LT.ONE) THEN
       ST1     = HALF*(ONE+CLP)
       ST2     = ONE - ST1
       st1s    = st1
       hstkd   = hst(kd)
       qstkd   = qst(kd)
       ltlkd   = ltl(kd)
       q0ukd   = q0u(kd)
       q0dkd   = q0d(kd)
       dlbkd   = dlb(kd)
       qrbkd   = qrb(kd)
       !
       HST(KD) = HST(KD)*ST1 + HST(KD1)*ST2
       HOS     = HOL(KD)*ST1 + HOL(KD1)*ST2
       QST(KD) = QST(KD)*ST1 + QST(KD1)*ST2
       QOS     = QOL(KD)*ST1 + QOL(KD1)*ST2
       QLS     = CLL(KD)*ST1 + CLL(KD1)*ST2
       QIS     = CIL(KD)*ST1 + CIL(KD1)*ST2
       LTL(KD) = LTL(KD)*ST1 + LTL(KD1)*ST2
       !
       DLB(KD) = DLB(KD)*CLP
       qrb(KD) = qrb(KD)*CLP
       ETA(KD) = ETA(KD)*CLP
       GMS(KD) = GMS(KD)*CLP
       Q0U(KD) = Q0U(KD)*CLP
       Q0D(KD) = Q0D(KD)*CLP
    ENDIF
    !
    !
    !***********************************************************************
    !
    !    Critical workfunction is included in this version
    !
    ACR = 0.0_r8
    TEM = PRL(KD1) - (PRL(KD1)-PRL(KD)) * CLP * HALF
    tx1 = PRL(KBL) - TEM
!!!   tx2 = min(700.0_r8,max(tx1,100.0_r8))
    !!    tx2 = min(900.0_r8,max(tx1,100.0_r8))
    !     tx2 = min(800.0_r8,max(tx1,200.0_r8))
    !     rel_fac =  dt * 600.0_r8 / (3600.0_r8*((800.0_r8-tx2)*0.5_r8+(tx2-200.0_r8)*3.0_r8))
    !     rel_fac =  dt         / (6.0_r8*((800.0_r8-tx2)*0.5_r8+(tx2-200.0_r8)*3.0_r8))
    !     rel_fac =  dt * facdt / (4.5_r8*((900.0_r8-tx2)*0.5_r8+(tx2-100.0_r8)*3.0_r8))
    !!    rel_fac =  dt * facdt / (4.5_r8*((900.0_r8-tx2)*0.5_r8+(tx2-100.0_r8)*6.0_r8))
!!!   rel_fac =  dt * facdt / (6.0_r8*((700.0_r8-tx2)*1.0_r8+(tx2-100.0_r8)*3.0_r8))
!!!!  rel_fac =  dt * facdt / (6.0_r8*((700.0_r8-tx2)*0.5_r8+(tx2-100.0_r8)*2.0_r8))
    !
    !
    !     tx2 = min(800.0_r8,max(tx1,100.0_r8))
    !     tem1    = log(tx2*0.01_r8) / log(8.0_r8)
    tx2 = MIN(900.0_r8,MAX(tx1,100.0_r8))
    tem1    = LOG(tx2*0.01_r8) / LOG(10.0_r8)
    !     rel_fac = (dt * facdt)  / (3600.0_r8 * (tem1*4.0_r8 + (1-tem1)*1.0_r8))
    rel_fac = (dt * facdt)  / (3600.0_r8 * (tem1*3.0_r8 + (1-tem1)*1.0_r8))
    !     rel_fac = (dt * facdt)  / (3600.0_r8 * (tem1*2.0_r8 + (1-tem1)*1.0_r8))
    !cnt  rel_fac = (dt * facdt)  / (3600.0_r8 * 1.5_r8)
    !     rel_fac = 0.3_r8 
    !
    rel_fac = MAX(zero, MIN(one,rel_fac))

    IF (CRTFUN) THEN
       CALL CRTWRK(TEM, CCWF, ST1)
       !       ACR = (PRL(K) - TEM) * ST1
       ACR = TX1 * ST1
    ENDIF
    !
    !===>  NORMALIZED MASSFLUX
    !
    !  ETA IS THE THICKNESS COMING IN AND THE MASS FLUX GOING OUT.
    !  GMS IS THE THICKNESS OF THE SQUARE; IT IS LATER REUSED FOR GAMMA_S
    !
    !     ETA(K) = ONE

    DO L=KB1,KD,-1
       ETA(L)  = ETA(L+1) + ALM * (ETA(L) + ALM * GMS(L))
    ENDDO
    DO L=KD,KBL
       ETAI(L) = 1.0_r8 / ETA(L)
    ENDDO

    !     if (lprnt) print *,' eta=',eta,' ii=',ii,' alm=',alm
    !
    !===>  CLOUD WORKFUNCTION
    !
    WFN   = ZERO
    AKM   = ZERO
    DET   = ZERO
    HCC   = HBL
    UNSAT = .FALSE.
    QTL   = QST(KB1) - GAF(KB1)*HST(KB1)
    TX1   = HBL
    !
    !     tem   = qst(kbl) - gaf(kbl)*hst(kbl)
    !     qtv   = 0.5_r8 * ((tem+qtl) + (gaf(kbl)+gaf(kb1))*hbl)
    !     det   = max(ZERO, qbl-qtv)
    !     qtv   = qbl - det
    !     det   = det + qlb + qib
    !!
    qtv   = qbl
    det   = qlb + qib
    !
    tx2   = 0.0_r8
    dpneg = 0.0_r8
    !
    DO L=KB1,KD1,-1
       DEL_ETA = ETA(L) - ETA(L+1)
       HCCP = HCC + DEL_ETA*HOL(L)
       !
       QTLP = QST(L-1) - GAF(L-1)*HST(L-1)
       QTVP = 0.5_r8 * ((QTLP+QTL)*ETA(L)                                &
            &              + (GAF(L)+GAF(L-1))*HCCP)
       ST1  = ETA(L)*Q0U(L) + ETA(L+1)*Q0D(L)
       DETP = (BKC(L)*DET - (QTVP-QTV)                                &
            &        + DEL_ETA*(QOL(L)+CLL(L)+CIL(L)) + ST1)  * AKC(L)

       !     if(lprnt) print *,' detp=',detp,' bkc=',bkc(l),' det=',det
       !     if (lprnt .and. kd .eq. 15) 
       !    &          print *,' detp=',detp,' bkc=',bkc(l),' det=',det
       !    &,' qtvp=',qtvp,' qtv=',qtv,' del_eta=',del_eta,' qol='
       !    &,qol(l),' st1=',st1,' akc=',akc(l)
       !
       TEM1   = AKT(L)   - QLL(L)
       TEM2   = QLL(L+1) - BKC(L)
       RNS(L) = TEM1*DETP  + TEM2*DET - ST1

       qtp    = 0.5_r8 * (qil(L)+qil(L-1))
       tem2   = MIN(qtp*(detp-eta(l)*qw00),                           &
            &               (1.0_r8-qtp)*(detp-eta(l)*qi00))
       st1    = MIN(tx2,tem2)
       tx2    = tem2
       !
       IF (rns(l) .LT. zero .OR. st1 .LT. zero) ep_wfn = .TRUE.
       IF (DETP .LE. ZERO) UNSAT = .TRUE.
       !        IF (DETP .LE. ZERO .or. rns(l) .lt. zero) UNSAT = .TRUE.

       ST1  = HST(L) - LTL(L)*con_FVirt*(QST(L)-QOL(L))


       TEM2 = HCCP   + DETP   * QTP * con_hfus
       !
       !     if(lprnt) print *,' hst=',hst(l),' ltl=',ltl(l),' con_FVirt=',con_FVirt
       !     if (lprnt .and. kd .eq. 15) 
       !    &          print *,' hst=',hst(l),' ltl=',ltl(l),' con_FVirt=',con_FVirt
       !    &,' qst=',qst(l),' qol=',qol(l),' hccp=',hccp,' detp=',detp
       !    *,' qtp=',qtp,' con_hfus=',con_hfus,' vtf=',vtf(l)

       ST2  = LTL(L) * VTF(L)
       TEM5 = CLL(L) + CIL(L)
       TEM3 = (TX1  - ETA(L+1)*ST1 - ST2*(DET-TEM5*eta(l+1))) * DLB(L)
       TEM4 = (TEM2 - ETA(L  )*ST1 - ST2*(DETP-TEM5*eta(l)))  * DLT(L)
       !
       !     if (lprnt) then
       !     if (lprnt .and. kd .eq. 12) then 
       !       print *,' tem3=',tem3,' tx1=',tx1,' st1=',st1,' eta1=',eta(l+1)
       !    &, ' st2=',st2,' det=',det,' tem5=',tem5,' dlb=',dlb(l)
       !       print *,' tem4=',tem4,' tem2=',tem2,' detp=',detp
       !    &, ' eta=',eta(l),' dlt=',dlt(l),' rns=',rns(l),' l=',l
       !       print *,' bt1=',tem3/(eta(l+1)*qrb(l))
       !    &,         ' bt2=',tem4/(eta(l)*qrt(l))
       !      endif

       ST1  = TEM3 + TEM4

       !     if (lprnt) print *,' wfn=',wfn,' st1=',st1,' l=',l,' ep_wfn=',
       !    &ep_wfn,' akm=',akm

       WFN = WFN + ST1       
       AKM = AKM - MIN(ST1,ZERO)

       !     if (lprnt) print *,' wfn=',wfn,' akm=',akm
       !     if (lprnt .and. kd .eq. 12) print *,' wfn=',wfn,' akm=',akm

       IF (st1 .LT. zero .AND. wfn .LT. zero) THEN
          dpneg = dpneg + prl(l+1) - prl(l)
       ENDIF

       !        BUY(L) = 0.5_r8 * (ETA(L+1) + ETA(L)) * ST1
       !        BUY(L) = ETA(L+1)*tem3 + ETA(L)*tem4
       !!       BUY(L) = tem3*ETAI(L+1) + tem4*ETAI(L)
       BUY(L) = 0.5_r8 * (tem3/(eta(l+1)*qrb(l)) + tem4/(eta(l)*qrt(l)))
       !        BUY(L) = 0.5_r8 * st1 / ((eta(l)+eta(l+1))*(qrb(l)+qrt(l)))
       !
       HCC = HCCP
       DET = DETP
       QTL = QTLP
       QTV = QTVP
       TX1 = TEM2

    ENDDO

    DEL_ETA = ETA(KD) - ETA(KD1)
    HCCP    = HCC + DEL_ETA*HOS
    !
    QTLP = QST(KD) - GAF(KD)*HST(KD)
    QTVP = QTLP*ETA(KD) + GAF(KD)*HCCP
    ST1  = ETA(KD)*Q0U(KD) + ETA(KD1)*Q0D(KD)
    DETP = (BKC(KD)*DET - (QTVP-QTV)                                  &
         &     + DEL_ETA*(QOS+QLS+QIS) + ST1) * AKC(KD)
    !
    TEM1    = AKT(KD)  - QLL(KD)
    TEM2    = QLL(KD1) - BKC(KD)
    RNS(KD) = TEM1*DETP  + TEM2*DET - ST1
    !
    IF (rns(kd) .LT. zero) ep_wfn = .TRUE.
    IF (DETP.LE.ZERO) UNSAT = .TRUE.
    !
888 CONTINUE

    !     if (lprnt) print *,' ep_wfn=',ep_wfn,' ii=',ii,' rns=',rns(kd)
    !    &,' clp=',clp,' hst(kd)=',hst(kd)

    IF (ep_wfn) THEN
       IF ((qw00 .EQ. 0.0_r8 .AND. qi00 .EQ. 0.0_r8)) RETURN
       IF (ii .EQ. 0) THEN
          ii  = 1
          IF (clp .GT. 0.0_r8 .AND. clp .LT. 1.0_r8) THEN
             hst(kd) = hstkd
             qst(kd) = qstkd
             ltl(kd) = ltlkd
             q0u(kd) = q0ukd
             q0d(kd) = q0dkd
             dlb(kd) = dlbkd
             qrb(kd) = qrbkd
          ENDIF
          DO l=kd,kb1
             FCO(L) = FCO(L) - q0u(l) - q0d(l)
             RNN(L) = RNN(L) - q0u(l)*zet(l) - q0d(l)*zet(l+1)
             GMH(L) = GMH(L) - q0u(l)*xi(l)  - q0d(l)*zet(l+1)
             ETA(L) = ZET(L) - ZET(L+1)
             GMS(L) = XI(L)  - XI(L+1)
             Q0U(L) = 0.0_r8
             Q0D(L) = 0.0_r8
          ENDDO
          qw00 = 0.0_r8
          qi00 = 0.0_r8

          !     if (lprnt) print *,' returning to 777 : ii=',ii,' qw00=',qw00,qi00
          !    &,' clp=',clp,' hst(kd)=',hst(kd)

          go to 777
       ELSE
          unsat = .TRUE.
       ENDIF
    ENDIF
    !
    !
    !     ST1 = 0.5 * (HST(KD)  - LTL(KD)*con_FVirt*(QST(KD)-QOS)
    !    &          +  HST(KD1) - LTL(KD1)*con_FVirt*(QST(KD1)-QOL(KD1)))
    !
    ST1 = HST(KD)  - LTL(KD)*con_FVirt*(QST(KD)-QOS)
    ST2 = LTL(KD)  * VTF(KD)
    TEM5 = (QLS + QIS) * eta(kd1)
    ST1  = HALF * (TX1-ETA(KD1)*ST1-ST2*(DET-TEM5))*DLB(KD)
    !
    !     if (lprnt) print *,' st1=',st1,' st2=',st2,' ltl=',ltl(kd)
    !    *,ltl(kd1),' qos=',qos,qol(kd1)

    WFN = WFN + ST1
    AKM = AKM - MIN(ST1,ZERO)   ! Commented on 08/26/02 - does not include top
    !

    !     BUY(KD) = 0.5_r8 * (ETA(KD1) + ETA(KD)) * ST1
    !     BUY(KD) = ETA(KD1) * ST1
    !!    BUY(KD) = ST1 * ETAI(KD1)
    BUY(KD) = ST1 / (ETA(KD1)*qrb(kd))
    !     BUY(KD) = 0.5_r8 * ST1 / (qrb(kd) * (eta(kd)+eta(kd1)))
    !
    !     if (lprnt) print *,' wfn=',wfn,' akm=',akm,' st1=',st1
    !    &,' dpneg=',dpneg

    DET = DETP
    HCC = HCCP
    AKM = AKM / WFN


    !***********************************************************************
    !
    !     If only to calculate workfunction save it and return
    !
    IF (WRKFUN) THEN
       IF (WFN .GE. 0.0_r8) WFNC = WFN
       RETURN
    ELSEIF (.NOT. CRTFUN) THEN
       ACR = WFNC
    ENDIF
    !
    !===>  THIRD CHECK BASED ON CLOUD WORKFUNCTION
    !
    CALCUP = .FALSE.

    !     TEM  =  MIN(CD*100.0_r8, MAX_NEG_BOUY)
    TEM  =  MIN(CD*200.0_r8, MAX_NEG_BOUY)
    !LL2  tem  = max_neg_bouy
    !     tem1 = dpneg / (prl(kbl)-prsm(kd))
    IF (WFN .GT. ACR .AND.  (.NOT. UNSAT)                             &
                                !    & .and. tem1 .le. 0.1_r8  .AND. AKM .LE. TEM) THEN
                                !    & .and. dpneg .lt. 100.0_r8  .AND. AKM .LE. TEM) THEN
         & .AND. dpneg .LT. 150.0_r8  .AND. AKM .LE. TEM) THEN
       !    & .and. dpneg .lt. 200.0_r8  .AND. AKM .LE. TEM) THEN
       !
       CALCUP = .TRUE.
    ENDIF

    !     if (lprnt) print *,' calcup=',calcup,' akm=',akm,' tem=',tem
    !    *,' unsat=',unsat,' clp=',clp,' rhc=',rhc,' cd=',cd,' acr=',acr
    !
    !===>  IF NO SOUNDING MEETS THIRD CONDITION, RETURN
    !
    !     if (lprnt .and. kd .eq. 15) stop
    IF (.NOT. CALCUP) RETURN
    !
    ! This is for not LL - 20050601
    IF (ALMIN2 .NE. 0.0_r8) THEN
       !       ST1 = 0.0_r8
       IF (ALMIN1 .NE. ALMIN2) ST1 = 1.0_r8 / MAX(ONE_M10,(ALMIN2-ALMIN1))
       IF (ALM .LT. ALMIN2) THEN
          !!         CLP = CLP * (ALM - ALMIN1) * ST1
          !!         CLP = CLP * (0.1_r8 + 0.9_r8*(ALM - ALMIN1) * ST1)
          CLP = CLP * MAX(0.0_r8, MIN(1.0_r8,(0.3_r8 + 0.7_r8*(ALM-ALMIN1)*ST1)))
          !          CLP = CLP * max(0.0_r8, min(1.0_r8,(0.2_r8 + 0.8_r8*(ALM-ALMIN1)*ST1)))
          !          CLP = CLP * max(0.0_r8, min(1.0_r8,(0.1_r8 + 0.9_r8*(ALM-ALMIN1)*ST1)))
       ENDIF
    ENDIF
    !
    !     if (lprnt) print *,' clp=',clp
    !
    CLP = CLP * RHC
    DO l=kd,kb1
       rnn(l) = rns(l)
    ENDDO
    DO L=KBL,K 
       RNN(L) = 0.0_r8 
    ENDDO
    !     if (lprnt) print *,' rnn=',rnn
    !
    !     If downdraft is to be invoked, do preliminary check to see
    !     if enough rain is available and then call DDRFT.
    !
    DDFT = .FALSE.
    IF (DNDRFT) THEN
       !
       TRAIN = 0.0_r8
       IF (CLP .GT. 0.0_r8) THEN
          DO L=KD,KB1
             TRAIN = TRAIN + RNN(L)
          ENDDO
       ENDIF

       PL = (PRL(KD1) + PRL(KD))*HALF
       TEM = PRL(K+1)*(1.0_r8-DPD*0.001_r8)
       !cnt    TEM = MIN(PRL(K+1)-DPD, PRL(KBL)-50.0_r8)
       !       TEM = MIN(PRL(K+1)-DPD, PRL(KBL)-300.0_r8)
       IF (TRAIN .GT. 1.0E-4_r8 .AND. PL .LE. TEM) DDFT  = .TRUE.
       !
       !       if (ddft) then
       !!        DO L=KBL-3,KD,-1
       !!        DO L=KBL-3,KD1,-1
       !         DO L=KBL-3,KD1+1,-1
       !           IF (BUY(L) .LT. 0.1_r8) THEN
       !             DDFT = .FALSE.
       !             EXIT
       !           ENDIF
       !         ENDDO
       !       endif
    ENDIF
    !
    !     if (lprnt) print *,' BEFORE CALLING DDRFT KD=',kd,' DDFT=',DDFT
    !    &,                  ' PL=',PL,' TRAIN=',TRAIN
    !     if (lprnt) print *,' buy=',(buy(l),l=kd,kb1)

    IF (DDFT) THEN
       !
       !     Call Downdraft scheme based on (Cheng and Arakawa, 1997)
       !
       CALL DDRFT(                                                     &
            &              K, KD                                               &
            &,             TLA, ALFIND                                         &
!            &,             TOL, QOL, HOL,   PRL, QST, HST, GAM, GAF, HBL, QBL  &
            &,             TOL, QOL, HOL,   PRL, QST, HST, GAM, GAF  &
            &,             QRB, QRT, BUY,   KBL, IDH, ETA, RNN, ETAI           &
            &,             ALM, TRAIN, DDFT                               &
 !           &,             ALM, WFN, TRAIN, DDFT                               &
            &,             ETD, HOD, QOD,   EVP, DOF, CLDFR, ETZ               &
                                !    &,             ETD, HOD, QOD,   EVP, DOF, CLDFRD, ETZ              &
            &,             GMS, GSD, GHD)               
       !    &,             TX1, TX2, TX3, TX4, TX5, TX6, TX7, TX8, TX9)

    ENDIF
    !
    !  No Downdraft case (including case with no downdraft soln)
    !  ---------------------------------------------------------
    !
    IF (.NOT. DDFT) THEN
       DO L=KD,K+1
          ETD(L) = 0.0_r8
          HOD(L) = 0.0_r8
          QOD(L) = 0.0_r8
       ENDDO
       DO L=KD,K
          EVP(L) = 0.0_r8
          ETZ(L) = 0.0_r8
       ENDDO

    ENDIF
    !     if (lprnt) print *,' hod=',hod
    !     if (lprnt) print *,' etd=',etd
    !
    !
    !===> CALCULATE GAMMAS  i.e. TENDENCIES PER UNIT CLOUD BASE MASSFLUX
    !           Includes downdraft terms!

    avh = 0.0_r8

    !
    !     Fraction of detrained condensate evaporated
    !
    !     tem1 = max(ZERO, min(HALF, (prl(kd)-FOUR_P2)*ONE_M2))
    !     tem1 = max(ZERO, min(HALF, (prl(kd)-300.0_r8)*0.005_r8))
    tem1 = 0.0_r8
    !     tem1 = 1.0_r8
    !     if (kd1 .eq. kbl) tem1 = 0.0_r8
    !
    tem2    = 1.0_r8 - tem1
    TEM = DET * QIL(KD)

    !     st1 = (HCC-ETA(KD)*HST(KD)) / (1.0_r8+gam(KD))

    st1 = (HCC+con_hfus*TEM-ETA(KD)*HST(KD)) / (1.0_r8+gam(KD))
    DS  = ETA(KD1) * (HOS- HOL(KD)) - con_hvap *(QOS - QOL(KD))
    DH  = ETA(KD1) * (HOS- HOL(KD))

    !     GMS(KD) = (DS + st1 - tem1*det*con_hvap ) * PRI(KD)
    !     GMS(KD) = (DS + st1 - tem1*(det*con_hvap +tem*con_hfus)) * PRI(KD)

    GMS(KD) = (DS + st1 - tem1*det*con_hvap -tem*con_hfus) * PRI(KD)
    !     GMS(KD) = (DS + st1 - det*con_hvap  - tem*con_hfus) * PRI(KD)
    GMH(KD) = PRI(KD) * (HCC-ETA(KD)*HOS + DH)

    !     GMH(KD) = PRI(KD) * (HCC-ETA(KD)*HOS + con_hfus*tem + DH)
    !     GMH(KD) = PRI(KD) * (HCC-ETA(KD)*HOS + con_hfus*tem*tem1 + DH)

    !     if (lprnt) print *,' gmhkd=',gmh(kd),' gmskd=',gms(kd)
    !    &,' det=',det,' tem=',tem,' tem1=',tem1,' tem2=',tem2
    !
    !      TENDENCY FOR SUSPENDED ENVIRONMENTAL ICE AND/OR LIQUID WATER
    !
    !     tem2    = 1.0_r8 - tem1
    !     tem3    = tem2 * (1.0_r8+con_hfus/con_hvap )
    !     QIL(KD) =     (tem3*TEM + ETA(KD1)*(QIS-CIL(KD))                  &

    QIL(KD) =     (tem2*TEM + ETA(KD1)*(QIS-CIL(KD))                  &
         &                        - ETA(KD)*QIS ) * PRI(KD)
    QLL(KD) = (tem2*(DET-TEM) + ETA(KD1)*(QLS-CLL(KD))                &
         &                          - ETA(KD)*QLS ) * PRI(KD)
    !
    GHD(KD) = 0.0_r8
    GSD(KD) = 0.0_r8
    !
    DO L=KD1,K
       ST1 = ONE - ALFINT(L,1)
       ST2 = ONE - ALFINT(L,2)
       ST3 = ONE - ALFINT(L,3)
       ST4 = ONE - ALFINT(L,4)
       ST5 = ONE - ALFIND(L)
       !      IF (L .LT. KBL) THEN
       !      IF (L .LT. K) THEN
       HB       = ALFINT(L,1)*HOL(L-1) + ST1*HOL(L)
       QB       = ALFINT(L,2)*QOL(L-1) + ST2*QOL(L)

       TEM      = ALFINT(L,4)*CIL(L-1) + ST4*CIL(L)
       TEM2     = ALFINT(L,3)*CLL(L-1) + ST3*CLL(L)

       TEM1     = ETA(L) * (TEM - CIL(L))
       TEM3     = ETA(L) * (TEM2 - CLL(L))

       HBD      = ALFIND(L)*HOL(L-1) + ST5*HOL(L)
       QBD      = ALFIND(L)*QOL(L-1) + ST5*QOL(L)

       TEM5     = ETD(L) * (HOD(L) - HBD)
       TEM6     = ETD(L) * (QOD(L) - QBD)
       !
       DH       = ETA(L) * (HB - HOL(L)) + TEM5
       DS       = DH - con_hvap  * (ETA(L) * (QB - QOL(L)) + TEM6)

       GMH(L)   = DH * PRI(L)
       GMS(L)   = DS * PRI(L)

       !     if (lprnt) print *,' gmh=',gmh(l),' gms=',gms(l)
       !    &,' dh=',dh,' ds=',ds,' qb=',qb,' qol=',qol(l),' eta=',eta(l)
       !    &,' hb=',hb,' hol=',hol(l),' l=',l,' hod=',hod(l)
       !    &,' etd=',etd(l),' qod=',qod(l),' tem5=',tem5,' tem6=',tem6
       !
       GHD(L)   = TEM5 * PRI(L)
       GSD(L)   = (TEM5 - con_hvap  * TEM6) * PRI(L)
       !
       QIL(L)   = TEM1 * PRI(L)
       QLL(L)   = TEM3 * PRI(L)

       TEM1     = ETA(L) * (CIL(L-1) - TEM)
       TEM3     = ETA(L) * (CLL(L-1) - TEM2)

       DH       = ETA(L) * (HOL(L-1) - HB) - TEM5
       DS       = DH - con_hvap  * ETA(L) * (QOL(L-1) - QB)                &
            &                 + con_hvap  * (TEM6 - EVP(L-1))

       GMH(L-1) = GMH(L-1) + DH * PRI(L-1)
       GMS(L-1) = GMS(L-1) + DS * PRI(L-1)
       !
       !     if (lprnt) print *,' gmh1=',gmh(l-1),' gms1=',gms(l-1)
       !    &,' dh=',dh,' ds=',ds,' qb=',qb,' qol=',qol(l-1)
       !    &,' hb=',hb,' hol=',hol(l-1),' evp=',evp(l-1)
       !
       GHD(L-1) = GHD(L-1) - TEM5 * PRI(L-1)
       GSD(L-1) = GSD(L-1) - (TEM5-con_hvap *(TEM6-EVP(L-1))) * PRI(L-1)

       QIL(L-1) = QIL(L-1) + TEM1 * PRI(L-1)
       QLL(L-1) = QLL(L-1) + TEM3 * PRI(L-1)
       !      ELSEIF (L .EQ. KBL) THEN
       !!       HB       = ALFINT(L)*HOL(L-1) + ST1*HBL
       !!       QB       = ALFINT(L)*QOL(L-1) + ST1*QBL
       !        HB       = ALFINT(L)*HOL(L-1) + ST1*HOL(L)
       !        QB       = ALFINT(L)*QOL(L-1) + ST1*QOL(L)

       !!       HB       = HBL
       !!       QB       = QBL
       !        HBD      = ALFINT(L)*HOL(L-1) + ST1*HOL(L)
       !        QBD      = ALFINT(L)*QOL(L-1) + ST1*QOL(L)

       !!       TEM      = ALFINQ(L)*CIL(L-1) + ST2*QIB
       !!       TEM2     = ALFINQ(L)*CLL(L-1) + ST2*QLB
       !        TEM      = ALFINQ(L)*CIL(L-1) + ST2*CIL(L)
       !        TEM2     = ALFINQ(L)*CLL(L-1) + ST2*CLL(L)

       !        TEM1     = ETA(L) * (TEM - QIB)
       !        TEM3     = ETA(L) * (TEM2 - QLB)

       !        TEM5     =  ETD(L) * (HOD(L) - HBD)
       !        TEM6     =  ETD(L) * (QOD(L) - QBD)

       !        tem4     = GRAVFAC * pris
       !        TX1      = ETA(L) * (HB - HBL) * TEM4
       !        TX2      = TX1 - con_hvap  * ETA(L) * (QB - QBL) * TEM4
       !        DH       = TEM5

       !        DS       =  DH - con_hvap  * (TEM6 + EVP(L))


       !        GMH(L)   = TX1 + DH * PRI(L)
       !        GMS(L)   = TX2 + DS * PRI(L)
       !
       !        GHD(L)   = TEM5 * PRI(L)
       !        GSD(L)   = (TEM5 - con_hvap  * (TEM6+EVP(L))) * PRI(L)
       !
       !        QIL(L)   = TEM1 * tem4
       !        QLL(L)   = TEM3 * tem4

       !        TEM1     = ETA(L) * (CIL(L-1) - TEM)
       !        TEM3     = ETA(L) * (CLL(L-1) - TEM2)

       !        DH       = ETA(L) * (HOL(L-1) - HB) - TEM5
       !        DS       = DH - con_hvap  * ETA(L) * (QOL(L-1) - QB)
       !    *                 + con_hvap  * (TEM6 - EVP(L-1))

       !        GMH(L-1) = GMH(L-1) + DH * PRI(L-1)
       !        GMS(L-1) = GMS(L-1) + DS * PRI(L-1)
       !
       !        GHD(L-1) = GHD(L-1) - TEM5 * PRI(L-1)
       !        GSD(L-1) = GSD(L-1) - (TEM5-con_hvap *(TEM6-EVP(L-1)))
       !    *                                  * PRI(L-1)

       !        QIL(L-1) = QIL(L-1) + TEM1 * PRI(L-1)
       !        QLL(L-1) = QLL(L-1) + TEM3 * PRI(L-1)
       !      ELSE
       !!       HB       = ALFINT(L)*HOL(L-1) + ST1*HOL(L)
       !!       QB       = ALFINT(L)*QOL(L-1) + ST1*QOL(L)
       !!                                                                                                                
       !!       TEM      = ALFINQ(L)*CIL(L-1) + ST2*CIL(L)
       !!       TEM2     = ALFINQ(L)*CLL(L-1) + ST2*CLL(L)
       !!                                                                                                                
       !!       TEM1     = ETA(L) * (TEM - CIL(L))
       !!       TEM3     = ETA(L) * (TEM2 - CLL(L))
       !!                                                                                                                
       !!       TEM5     = ETD(L) * (HOD(L) - HB)
       !!       TEM6     = ETD(L) * (QOD(L) - QB)
       !
       !!       DH       = ETA(L) * (HB - HOL(L)) + TEM5
       !!       DS       = DH - con_hvap  * (ETA(L) * (QB - QOL(L)) + TEM6)
       !                                                                                                                 
       !!       GMH(L)   = DH * PRI(L)
       !1       GMS(L)   = DS * PRI(L)

       !     if (lprnt) print *,' gmh=',gmh(l),' gms=',gms(l)
       !    &,' dh=',dh,' ds=',ds,' qb=',qb,' qol=',qol(l),' eta=',eta(l)
       !    &,' hb=',hb,' hol=',hol(l),' l=',l
       !
       !!       GHD(L)   = TEM5 * PRI(L)
       !!       GSD(L)   = (TEM5 - con_hvap  * TEM6) * PRI(L)
       !
       !!       QIL(L)   = TEM1 * PRI(L)
       !!       QLL(L)   = TEM3 * PRI(L)
       !
       !
       !
       !        HBD      = ALFINT(L)*HOL(L-1) + ST1*HOL(L)
       !        QBD      = ALFINT(L)*QOL(L-1) + ST1*QOL(L)
       !        TEM5     =  ETD(L) * (HOD(L) - HBD)
       !        TEM6     =  ETD(L) * (QOD(L) - QBD)
       !        DH       =  TEM5
       !        DS       =  DH - con_hvap  * (TEM6 + EVP(L))
       !
       !        GMH(L)   = TX1 + DH * PRI(L)
       !        GMS(L)   = TX2 + DS * PRI(L)
       !        GHD(L)   = DH * PRI(L)
       !        GSD(L)   = DS * PRI(L)
       !
       !        DH       = - TEM5
       !        DS       = DH  + con_hvap  * TEM6
       !        GMH(L-1) = GMH(L-1) + DH * PRI(L-1)
       !        GMS(L-1) = GMS(L-1) + DS * PRI(L-1)
       !
       !        GHD(L-1) = GHD(L-1) + DH * PRI(L-1)
       !        GSD(L-1) = GSD(L-1) + DS * PRI(L-1)
       !
       !        QIL(L)   = QIL(L-1)
       !        QLL(L)   = QLL(L-1)
       !      ENDIF

       avh = avh + gmh(l-1)*(prs(l)-prs(l-1))

    ENDDO
    !
    !!    TEM2   =  - con_hvap  * EVP(K) * PRI(K)
    !!    GMS(K) = GMS(K) + TEM2
    !!    GSD(K) = GSD(K) + TEM2
    !

    HBD  = HOL(K)
    QBD  = QOL(K)
    TEM5 =  ETD(K+1) * (HOD(K+1) - HBD)
    TEM6 =  ETD(K+1) * (QOD(K+1) - QBD)
    DH   = - TEM5
    DS   = DH  + con_hvap  * TEM6
    TEM1 = DH * PRI(K)
    TEM2 = (DS - con_hvap  * EVP(K)) * PRI(K)
    !!    TEM2 = - con_hvap  * EVP(K) * PRI(K)
    GMH(K) = GMH(K) + TEM1
    GMS(K) = GMS(K) + TEM2
    GHD(K) = GHD(K) + TEM1
    GSD(K) = GSD(K) + TEM2

    !     if (lprnt) print *,' gmhk=',gmh(k),' gmsk=',gms(k)
    !    &,' tem1=',tem1,' tem2=',tem2,' dh=',dh,' ds=',ds
    !
    avh = avh + gmh(K)*(prs(KP1)-prs(K))
    !
    tem4   = - GRAVFAC * pris
    TX1    = DH * tem4
    TX2    = DS * tem4
    !
    DO L=KBL,K
       GMH(L) = GMH(L) + TX1
       GMS(L) = GMS(L) + TX2
       GHD(L) = GHD(L) + TX1
       GSD(L) = GSD(L) + TX2
       !
       avh = avh + tx1*(prs(l+1)-prs(l))
    ENDDO

    !     DO L=KBL,K
    !       tem = (eta(l+1) - eta(l)) * pri(l)
    !       tx1 = dh * tem
    !       tx2 = ds * tem
    !       GMH(L) = GMH(L) + TX1
    !       GMS(L) = GMS(L) + TX2
    !       GHD(L) = GHD(L) + TX1
    !       GSD(L) = GSD(L) + TX2

    !       avh = avh + tx1*(prs(l+1)-prs(l))
    !     ENDDO
    !
    !     if (lprnt) then
    !        print *,' gmh=',gmh
    !        print *,' gms=',gms(KD:K)
    !     endif
    !
    !***********************************************************************
    !***********************************************************************

    !===>  KERNEL (AKM) CALCULATION BEGINS

    !===>  MODIFY SOUNDING WITH UNIT MASS FLUX
    !
    !     TESTMB = 0.01_r8

    DO L=KD,K

       TEM1   = GMH(L)
       TEM2   = GMS(L)
       HOL(L) = HOL(L) +  TEM1*TESTMB
       QOL(L) = QOL(L) + (TEM1-TEM2)  * (TESTMB/con_hvap )
       !    &                   +  con_hfus*QIL(L))  * (TESTMB/con_hvap )
       HST(L) = HST(L) +  TEM2*(ONE+GAM(L))*TESTMB
       QST(L) = QST(L) +  TEM2*GAM(L)*(TESTMB/con_hvap )
       CLL(L) = CLL(L) + QLL(L) * TESTMB
       CIL(L) = CIL(L) + QIL(L) * TESTMB
    ENDDO
    !

    IF (alm .GT. 0.0_r8) THEN
       HOS = HOS + GMH(KD)  * TESTMB
       QOS = QOS + (GMH(KD)-GMS(KD)) * (TESTMB/con_hvap )
       !    &          + con_hfus*QIL(KD)) * (TESTMB/con_hvap )

       QLS     = QLS + QLL(KD) * TESTMB
       QIS     = QIS + QIL(KD) * TESTMB
    ELSE
       st2 = 1.0_r8 - st1s
       HOS = HOS + (st1s*GMH(KD)+st2*GMH(KD1))  * TESTMB
       QOS = QOS + (st1s * (GMH(KD)-GMS(KD))                           &
            &            +  st2  * (GMH(KD1)-GMS(KD1))) * (TESTMB/con_hvap )
       !       QOS = QOS + (st1s * (GMH(KD)-GMS(KD)+con_hfus*QIL(KD))              &
       !    &            +  st2  * (GMH(KD1)-GMS(KD1)+con_hfus*QIL(KD1)))          &
       !    &            * (TESTMB/con_hvap )
       HST(kd) = HST(kd) + (st1s*GMS(kd)*(ONE+GAM(kd))                 &
            &                    +  st2*gms(kd1)*(ONE+GAM(kd1))) * TESTMB
       QST(kd) = QST(kd) + (st1s*GMS(kd)*GAM(kd)                       &
            &                    +  st2*gms(kd1)*gam(kd1)) * (TESTMB/con_hvap )

       QLS     = QLS + (st1s*QLL(KD)+st2*QLL(KD1)) * TESTMB
       QIS     = QIS + (st1s*QIL(KD)+st2*QIL(KD1)) * TESTMB
    ENDIF

    !
    TEM = PRL(K+1) - PRL(K)
    HBL = HOL(K) * TEM
    QBL = QOL(K) * TEM
    QLB = CLL(K) * TEM
    QIB = CIL(K) * TEM
    DO L=KM1,KBL,-1
       TEM = PRL(L+1) - PRL(L)
       HBL = HBL + HOL(L) * TEM
       QBL = QBL + QOL(L) * TEM
       QLB = QLB + CLL(L) * TEM
       QIB = QIB + CIL(L) * TEM
    ENDDO
    HBL = HBL * PRIS
    QBL = QBL * PRIS
    QLB = QLB * PRIS
    QIB = QIB * PRIS

    !     if (lprnt) print *,' hbla=',hbl,' qbla=',qbl

    !***********************************************************************

    !===>  CLOUD WORKFUNCTION FOR MODIFIED SOUNDING, THEN KERNEL (AKM)
    !
    AKM = ZERO
    TX1 = ZERO
    QTL = QST(KB1) - GAF(KB1)*HST(KB1)
    QTV = QBL
    HCC = HBL
    TX2 = HCC
    TX4 = (con_hfus*0.5_r8)*MAX(ZERO,MIN(ONE,(TCR-TCL-TOL(KB1))*TCRF))
    !     TX4 = (con_hfus*0.5_r8)*MAX(0.0_r8,MIN(1.0_r8,(TCR-TOL(KB1))*TCRF))
    !
    !     tem   = qst(kbl) - gaf(kbl)*hst(kbl)
    !     qtv   = 0.5_r8 * ((tem+qtl) + (gaf(kbl)+gaf(kb1))*hbl)
    !     tx1   = max(ZERO, qbl-qtv)
    !     qtv   = qbl - tx1
    !     tx1   = tx1 + qib + qlb
    !
    qtv   = qbl
    tx1   = qib + qlb
    !

    DO L=KB1,KD1,-1
       DEL_ETA = ETA(L) - ETA(L+1)
       HCCP = HCC + DEL_ETA*HOL(L)
       !
       QTLP = QST(L-1) - GAF(L-1)*HST(L-1)
       QTVP = 0.5_r8 * ((QTLP+QTL)*ETA(L)                                &
            &                +(GAF(L)+GAF(L-1))*HCCP)

       DETP = (BKC(L)*TX1 - (QTVP-QTV)                                &
            &        +  DEL_ETA*(QOL(L)+CLL(L)+CIL(L))                         &
            &        +  ETA(L)*Q0U(L) + ETA(L+1)*Q0D(L)) * AKC(L)
       IF (DETP .LE. ZERO) UNSAT = .TRUE.

       ST1  = HST(L) - LTL(L)*con_FVirt*(QST(L)-QOL(L))

       TEM2 = (con_hfus*0.5_r8)*MAX(ZERO,MIN(ONE,(TCR-TCL-TOL(L-1))*TCRF))
       !        TEM2 = (con_hfus*0.5_r8)*MAX(0.0_r8,MIN(1.0_r8,(TCR-TOL(L-1))*TCRF))
       TEM1 = HCCP + DETP * (TEM2+TX4)

       ST2  = LTL(L) * VTF(L)
       TEM5 = CLL(L) + CIL(L)
       AKM  = AKM +                                                   &
            &     (  (TX2  -ETA(L+1)*ST1-ST2*(TX1-TEM5*eta(l+1))) * DLB(L)     &
            &      + (TEM1 -ETA(L  )*ST1-ST2*(DETP-TEM5*eta(l)))  * DLT(L) )
       !
       HCC  = HCCP
       TX1  = DETP
       TX2  = TEM1
       QTL  = QTLP
       QTV  = QTVP
       TX4  = TEM2
    ENDDO
    !
    IF (unsat) RETURN
    !
    !  Eventhough we ignore the change in lambda, we still assume
    !  that the cLoud-top contribution is zero; as though we still
    !  had non-bouyancy there.
    !
    !
    ST1 = HST(KD)  - LTL(KD)*con_FVirt*(QST(KD)-QOS)
    ST2 = LTL(KD)  * VTF(KD)
    TEM5 = (QLS + QIS) * eta(kd1)
    AKM  = AKM + HALF * (TX2-ETA(KD1)*ST1-ST2*(TX1-TEM5)) * DLB(KD)
    !
    AKM = (AKM - WFN) * (ONE/TESTMB)


    !***********************************************************************

    !===>   MASS FLUX

    !     if (acr .gt. 0.0_r8) then
    !       tem = max(0.01_r8, min(0.05_r8, (wfn-acr)/acr))
    !     else
    !       tem = max(0.01_r8, 0.05_r8*min(1.0_r8, wfn))
    !     endif
!!!   tem2 = (rasalf*(tem-0.01_r8) + 0.05_r8 - tem) * oneopt4
    !     tem2 = (rel_fac*(tem-0.01_r8) + 0.05_r8 - tem) * oneopt4
    tem2 = rel_fac
    !
    AMB = - (WFN-ACR) / AKM
    !
    IF(lprnt) PRINT *,' wfn=',wfn,' acr=',acr,' akm=',akm             &
         &,' amb=',amb,' KD=',kd,' cldfrd=',cldfrd,' tem2=',tem2            &
         &,' rel_fac=',rel_fac,' prskd=',prs(kd)

    !===>   RELAXATION AND CLIPPING FACTORS
    !
    AMB = AMB * CLP * tem2

!!!   if (DDFT) AMB = MIN(AMB, ONE/CLDFRD)

    !===>   SUB-CLOUD LAYER DEPTH LIMIT ON MASS FLUX

    AMBMAX = (PRL(KP1)-PRL(KBL))*(FRACBL*GRAVCON)
    AMB    = MAX(MIN(AMB, AMBMAX),ZERO)


    !     if(lprnt) print *,' AMB=',amb,' clp=',clp,' ambmax=',ambmax
    !***********************************************************************
    !*************************RESULTS***************************************
    !***********************************************************************

    !===>  PRECIPITATION AND CLW DETRAINMENT
    !
    avt = 0.0_r8
    avq = 0.0_r8
    avr = dof

    !
    DSFC = DSFC + AMB * ETD(K) * (1.0_r8/DT)
    !
    !     DO L=KBL,KD,-1
    DO L=K,KD,-1
       PCU(L) = PCU(L) + AMB*RNN(L)      !  (A40)
       avr = avr + rnn(l)
       !     if(lprnt) print *,' avr=',avr,' rnn=',rnn(l),' l=',l
    ENDDO
    !
    !===> TEMPARATURE AND Q CHANGE AND CLOUD MASS FLUX DUE TO CLOUD TYPE KD
    !
    TX1  = AMB * (ONE/con_cp)
    TX2  = AMB * (ONE/con_hvap )
    DO L=KD,K
       ST1    = GMS(L)*TX1
       TOI(L) = TOI(L) + ST1
       TCU(L) = TCU(L) + ST1
       TCD(L) = TCD(L) + GSD(L) * TX1
       !
       !       st1 = st1 - (con_hfus/con_cp) * QIL(L) * AMB
       st1 = st1 - (con_hvap /con_cp) * (QIL(L) + QLL(L)) * AMB

       avt = avt + st1 * (prs(l+1)-prs(l))

       FLX(L)  = FLX(L)  + ETA(L)*AMB
       FLXD(L) = FLXD(L) + ETD(L)*AMB
       !
       QII(L)  = QII(L) + QIL(L) * AMB
       TEM     = 0.0_r8

       QLI(L)  = QLI(L) + QLL(L) * AMB + TEM

       ST1     = (GMH(L)-GMS(L)) * TX2
       !       ST1     = (GMH(L)-GMS(L)+con_hfus*QIL(L)) * TX2

       QOI(L)  = QOI(L) + ST1
       QCU(L)  = QCU(L) + ST1
       QCD(L)  = QCD(L) + (GHD(L)-GSD(L)) * TX2
       !
       avq = avq + (st1+(QLL(L)+QIL(L))*amb) * (prs(l+1)-prs(l))
       !       avq = avq + st1 * (prs(l+1)-prs(l))
       !       avr = avr + (QLL(L) + QIL(L)*(1+con_hfus/con_hvap ))
       !       avr = avr + (QLL(L) + QIL(L))
       !    *                  * (prs(l+1)-prs(l)) * gravcon

       !     if(lprnt) print *,' avr=',avr,' qll=',qll(l),' l=',l
       !    &,' qil=',qil(l)

    ENDDO
    avr = avr * amb
    !
    !      Correction for negative condensate!
    !     if (advcld) then
    !       do l=kd,k
    !         if (qli(l) .lt. 0.0_r8) then
    !           qoi(l) = qoi(l) + qli(l)
    !           toi(l) = toi(l) - (con_hvap /con_cp) * qli(l)
    !           qli(l) = 0.0_r8
    !         endif
    !         if (qii(l) .lt. 0.0_r8) then
    !           qoi(l) = qoi(l) + qii(l)
    !           toi(l) = toi(l) - ((con_hvap +con_hfus)/con_cp) * qii(l)
    !           qii(l) = 0.0_r8
    !         endif
    !       enddo
    !     endif

    !
    !
    !     if (lprnt) then
    !       print *,' For KD=',KD
    !       avt = avt * con_cp * 100.0_r8*86400.0_r8 / (con_hvap *DT*con_g)
    !       avq = avq *  100.0_r8*86400.0_r8 / (DT*con_g)
    !       avr = avr * 86400.0_r8 / DT
    !       print *,' avt=',avt,' avq=',avq,' avr=',avr,' avh='
    !    *   ,avh,' alm=',alm,' DDFT=',DDFT,' KD=',KD
    !    &,' TOIK-',toi(k),' TOIK-1=',toi(k-1),' TOIK-2=',toi(k-2)
    !        if (kd .eq. 12 .and. .not. ddft) stop
    !       if (avh .gt. 0.1_r8 .or. abs(avt+avq) .gt. 1.0e-5_r8 .or.
    !    &      abs(avt-avr) .gt. 1.0e-5_r8 .or. abs(avr+avq) .gt. 1.0e-5_r8) stop
    !
    !     if (lprnt) then
    !       print *,' For KD=',KD
    !       print *,' TCU=',(tcu(l),l=kd,k)
    !       print *,' QCU=',(Qcu(l),l=kd,k)
    !     endif
    !
    TX1 = 0.0_r8
    TX2 = 0.0_r8
    !
    !     REEVAPORATION OF FALLING CONVECTIVE RAIN
    !
    IF (REVAP) THEN
       !      AFC     = -(1.04E-4_r8*DT)*(3600._r8/DT)**0.578_r8
       !      rknob   = 5.0_r8
       !      rknob   = 3.0_r8
       !      rknob   = 0.0_r8
       !      rknob   = 1.9_r8
       !      rknob   = 2.0_r8
       !
       tem = 0.0_r8
       DO l=kd,kbl
          !        tem = tem + pcu(l)
          IF (L .LT. IDH .OR. (.NOT. DDFT)) THEN
             tem = tem + amb * rnn(l)
          ENDIF
       ENDDO
       tem = tem + amb * dof
       tem = tem * (3600.0_r8/dt)
       !      tem1 = 4.0E10_r8/max(garea,one) * sqrt((prl(kbl)-prl(kd))/prl(K+1))
       !LLLL  tem1 = rknob * sqrt(sqrt(4.0E10_r8/max(garea,one)))
       !!     tem2 = sqrt(4.0E10_r8/max(garea,one))
       !!     tem1 = sqrt(tem2)
       !!     tem1 = rknob * max(one, tem1*tem2)
       tem1 = MAX(1.0_r8, MIN(100.0_r8,SQRT((5.0E10_r8/MAX(garea,one)))))
       !Cntnewtem1 = max(1.0_r8, min(100.0_r8,(5.0E10_r8/max(garea,one))))
       !Cnt   tem1 = max(1.0_r8, min(100.0_r8,(4.0E10_r8/max(garea,one))))
       !      tem1 = rknob * sqrt(4.0E10_r8/max(garea,one))
       !      clfrac = ((clfa*tem + clfb)*tem + clfc)*tem + clfd
       !      clfrac = min(point3, max(point01,clfrac))

       !      if (lprnt) print *,' clfr0=',clf(tem),' tem=',tem,' tem1=',tem1

       clfrac = MAX(ZERO, MIN(ONE, rknob*clf(tem)*tem1))

       !      if (lprnt) print *,' cldfrd=',cldfrd,' amb=',amb
       !    &,' clfrac=',clfrac

       !      TX3    = AMB*ETA(KD)*PRI(KD)
       !      CLDFRD = MIN(AMB*CLDFRD, ONE)
       !      if(lprnt) print *,' cldfrd=',cldfrd,' amb=',amb
       !      CLDFRD = MIN(AMB*CLDFRD, 1.0_r8)
       !
       !     if(lprnt) print *,' tx3=',tx3,' etakd=',eta(kd),' pri=',pri(kd)
       !     if(lprnt) print *,' RNN=',RNN(kd:k)
       !
       !cnt   DO L=KD,K
       DO L=KD,KBL         ! Testing on 20070926
          !        clvfr = (prl(l)+prl(l+1)) / (prl(k)+prl(k+1))
!!!      clvfr = 0.5_r8 * (prl(l)+prl(l+1)) / prl(k+1)
          !        clvfr = min(1.0_r8, clvfr * clvfr)
          !                                                 for L=KD,K
          IF (L .GE. IDH .AND. DDFT) THEN
             TX2    = TX2 + AMB * RNN(L)
             CLDFRD = MIN(AMB*CLDFR(L), clfrac)
!!!        CLDFRD = MIN(AMB*CLDFR(L), ONE)
             !          if (l .eq. kbl) tx2 = tx2 + amb * dof
             !          if (l .eq. kbl) tx1 = tx1 + amb * dof
          ELSE
             TX1 = TX1 + AMB * RNN(L)
          ENDIF
          tx4 = zfac * phil(l)
          tx4 = (one - tx4 * (one - half*tx4)) * afc
          !
          !        CLFRAC = MIN(TX3*rknob*1.1_r8, ONE)
          !        CLFRAC = ONE
          !        CLFRAC = MIN(TX3*rknob, ONE)
          !        CLFRAC = MIN(TX3*rknob, 1.0_r8)

          IF (TX1 .GT. 0.0_r8 .OR. TX2 .GT. 0.0_r8) THEN
             TEQ     = TOI(L)
             QEQ     = QOI(L)
             PL      = 0.5_r8 * (PRL(L+1)+PRL(L))

             ST1     = MAX(ZERO, MIN(ONE, (TCR-TEQ)*TCRF))
             !         ST1     = MAX(0.0_r8, MIN(1.0_r8, (TCR-TEQ)*TCRF))
             ST2     = ST1*ELFOCP + (1.0_r8-ST1)*ELOCP

             CALL QSATCN ( TEQ,PL,QSTEQ,DQDT)
             !
             !         tx8  = 10.0_r8 * fpvs(teq)                  ! fpvs is in centibars!
             !         tx9  = 1.0_r8 / max(pl + epsm1 * tx8, 1.0e-10_r8)
             !         qsteq  = MIN(eps*tx8*tx9, 1.0_r8)
             !         dqdt = pl * qsteq * con_hvap  *  tx9 / (teq*teq*rv)
             !
             DELTAQ = 0.5_r8 * (QSTEQ*rhc_ls(l)-QEQ) / (1.0_r8+ST2*DQDT)
             !
             QEQ    = QEQ + DELTAQ
             TEQ    = TEQ - DELTAQ*ST2
             !
             TEM1   = MAX(ZERO, MIN(ONE, (TCR-TEQ)*TCRF))
             !         TEM1   = MAX(0.0_r8, MIN(1.0_r8, (TCR-TEQ)*TCRF))
             TEM2   = TEM1*ELFOCP + (1.0_r8-TEM1)*ELOCP

             CALL QSATCN ( TEQ,PL,QSTEQ,DQDT)
             !
             !         tx8  = 10.0_r8 * fpvs(teq)                  ! fpvs is in centibars!
             !         tx9  = 1.0_r8 / max(pl + epsm1 * tx8, 1.0e-10_r8)
             !         qsteq  = MIN(eps*tx8*tx9, 1.0_r8)
             !         dqdt = pl * qsteq * con_hvap  *  tx9 / (teq*teq*rv)
             !
             DELTAQ = (QSTEQ*rhc_ls(l)-QEQ) / (1.0_r8+TEM2*DQDT)
             !
             QEQ    = QEQ + DELTAQ
             TEQ    = TEQ - DELTAQ*TEM2

             IF (QEQ .GT. QOI(L)) THEN
                POTEVAP = (QEQ-QOI(L))*(PRL(L+1)-PRL(L))*GRAVCON
                !           POTEVAP = (QEQ-QOI(L))*(PRL(L+1)-PRL(L))*GRAVCON * 0.85_r8

                !           TEM3    = SQRT(PL*0.001_r8)
                tem4    = 0.0_r8
                IF (tx1 .GT. 0.0_r8)                                           &
                     &      TEM4    = POTEVAP * (1.0_r8 - EXP( tx4*TX1**0.57777778_r8 ) )
                !    &      TEM4    = POTEVAP * (1.0_r8 - EXP( AFC*tx4*SQRT(TX1) ) )
                !    &      TEM4    = POTEVAP * (1.0_r8 - EXP( AFC*SQRT(TX1*TEM3) ) )
                !    &      TEM4    = POTEVAP * (1.0_r8 - EXP(-0.32_r8*SQRT(DT*TX1*0.001_r8) ) )
                ACTEVAP = MIN(TX1, TEM4*CLFRAC)
                !           ACTEVAP = MIN(TX1, TEM4*CLFRAC*clvfr)
                !     if(lprnt) print *,' L=',L,' actevap=',actevap,' tem4=',tem4,
                !    &' clfrac='
                !    &,clfrac,' potevap=',potevap,'efac=',AFC*SQRT(TX1*TEM3)
                !    &,' tx1=',tx1
                IF (tx1 .LT. rainmin*dt) actevap = MIN(tx1, potevap)
                !
                tem4    = 0.0_r8
                IF (tx2 .GT. 0.0_r8)                                           &
                     &      TEM4    = POTEVAP * (1.0_r8 - EXP( tx4*TX2**0.57777778_r8 ) )
                !    &      TEM4    = POTEVAP * (1.0_r8 - EXP( AFC*tx4*SQRT(TX2) ) )
                !    &      TEM4    = POTEVAP * (1.0_r8 - EXP( AFC*SQRT(TX2*TEM3) ) )
                !    &      TEM4    = POTEVAP * (1.0_r8 - EXP(-0.32_r8*SQRT(DT*TX2*0.001_r8) ) )
                TEM4    = MIN(MIN(TX2, TEM4*CLDFRD), potevap-actevap)
                IF (tx2 .LT. rainmin*dt) tem4 = MIN(tx2, potevap-actevap)
                !
                TX1     = TX1 - ACTEVAP
                TX2     = TX2 - TEM4
                ST1     = (ACTEVAP+TEM4) * PRI(L)
                QOI(L)  = QOI(L) + ST1
                QCU(L)  = QCU(L) + ST1
                !

                ST1     = ST1 * ELOCP
                TOI(L)  = TOI(L) - ST1 
                TCU(L)  = TCU(L) - ST1
             ENDIF
          ENDIF
       ENDDO
       !
!!!    CUP = CUP + TX1 + TX2
       CUP = CUP + TX1 + TX2 + DOF * AMB
    ELSE
       DO L=KD,K
          TX1 = TX1 + AMB * RNN(L)
       ENDDO
       CUP = CUP + TX1 + DOF * AMB
    ENDIF

    !     CUP = CUP + TX1 + TX2 + DOF * AMB
    !     CUP = CUP + TX1 + TX2
    !     if (lprnt) print *,' tx1=',tx1,' tx2=',tx2,' dof=',dof
    !    &,' cup=',cup*86400/dt,' amb=',amb
    !    &,' amb=',amb,' cup=',cup,' clfrac=',clfrac,' cldfrd=',cldfrd
    !    &,' ddft=',ddft,' kd=',kd,' kbl=',kbl,' k=',k
    !
    !    MIXING OF PASSIVE TRACERS
    !
    DO N=1,M

       DO L=KD,K
          HOL(L) = ROI(L,N)
       ENDDO
       !
       HCC     = RBL(N)
       HOD(KD) = HOL(KD)
       !      Compute downdraft properties for the tracer
       DO L=KD1,K
          ST1 = ONE - ALFIND(L)
          HB  = ALFIND(L)  * HOL(L-1) + ST1 * HOL(L)
          IF (ETZ(L-1) .NE. 0.0_r8) THEN
             DEL_ETA = ETD(L) - ETD(L-1)
             TEM     = 1.0_r8 / ETZ(L-1)
             IF (DEL_ETA .GT. 0.0_r8) THEN
                HOD(L) = (ETD(L-1)*(HOD(L-1)-HOL(L-1))                     &
                     &                +  ETD(L)  *(HOL(L-1)-HB)                         &
                     &                +  ETZ(L-1)*HB) * TEM
             ELSE
                HOD(L) = (ETD(L-1)*(HOD(L-1)-HB) + ETZ(L-1)*HB) * TEM
             ENDIF
          ELSE
             HOD(L) = HB
          ENDIF
       ENDDO

       DO L=KB1,KD,-1
          HCC = HCC + (ETA(L)-ETA(L+1))*HOL(L)
       ENDDO
       !
       GMH(KD) = PRI(KD) * (HCC-ETA(KD)*HOL(KD))
       DO L=KD1,K
          ST1 = ONE - ALFINT(L,N+4)
          ST2 = ONE - ALFIND(L)
          !        IF (L .LT. KBL) THEN
          HB       = ALFINT(L,N+4) * HOL(L-1) + ST1 * HOL(L)
          HBD      = ALFIND(L) * HOL(L-1) + ST2 * HOL(L)
          TEM5     = ETD(L)    * (HOD(L) - HBD)
          !!         DH       = ETA(L)    * (HB - HOL(L)) * trcfac(n) + TEM5
          DH       = ETA(L)    * (HB - HOL(L)) + TEM5
          !!         GMH(L  ) = DH * PRI(L)
          GMH(L  ) = DH * PRI(L) * trcfac(n,l)
          !!         DH       = ETA(L)    * (HOL(L-1) - HB) * trcfac(n) - TEM5
          !!         GMH(L-1) = GMH(L-1)  + DH * PRI(L-1)
          DH       = ETA(L)    * (HOL(L-1) - HB) - TEM5
          GMH(L-1) = GMH(L-1)  + DH * PRI(L-1) * trcfac(n,l)
          !        ELSEIF (L .EQ. KBL) THEN
          !          HB       = ALFINT(L) * HOL(L-1) + ST1 * RBL(N)
          !          HBD      = ALFINT(L) * HOL(L-1) + ST1 * HOL(L)
          !          DH       = ETD(L)    * (HOD(L) - HBD)
          !          tem4     = GRAVFAC * pris
          !          TX1      = ETA(L)    * (HB - RBL(N)) * TEM4
          !          GMH(L)   = (TX1      + DH * PRI(L)) * trcfac(n)
          !          DH       = ETA(L)    * (HOL(L-1) - HB) - DH
          !          GMH(L-1) = GMH(L-1)  + DH * PRI(L-1) * trcfac(n)
          !        ELSE
          !          HBD      = ALFINT(L) * HOL(L-1) + ST1 * HOL(L)
          !          DH       = ETD(L)    * (HOD(L) - HBD)
          !          GMH(L)   = (TX1      + DH * PRI(L)) * trcfac(n)
          !          GMH(L-1) = GMH(L-1)  - DH * PRI(L-1) * trcfac(n)
          !        ENDIF
       ENDDO
       !
       DO L=KD,K
          ST1      = GMH(L)*AMB
          ROI(L,N) = HOL(L)   + ST1
          RCU(L,N) = RCU(L,N) + ST1
       ENDDO
    ENDDO                             ! Tracer loop M

    !     if (lprnt) print *,' toio=',toi
    !     if (lprnt) print *,' qoio=',qoi
    !     if (lprnt .and. kd .eq. 41) stop
    !     if (toi(K)-toi(k-1) .lt. 20.0_r8) stop
    !***********************************************************************
    !***********************************************************************
    !***********************************************************************

    RETURN
  END SUBROUTINE CLOUD


  !-----------------------------------------------------------------------------------------

  SUBROUTINE DDRFT(                                                 &
       &                  K, KD                                           &
       &,                 TLA, ALFIND                                     &
       &,                 TOL, QOL, HOL, PRL, QST, HST, GAM, GAF&
!       &,                 TOL, QOL, HOL, PRL, QST, HST, GAM, GAF, HBL, QBL&
       &,                 QRB, QRT, BUY, KBL, IDH, ETA, RNN, ETAI         &
       &,                 ALM, TRAIN, DDFT                           &
!       &,                 ALM, WFN, TRAIN, DDFT                           &
       &,                 ETD, HOD, QOD, EVP, DOF, CLDFRD, WCB            &
       &,                 GMS, GSD, GHD)                   




    !    &,                 GMS, GSD, GHD,lprnt)
    !    &,                 TX1, TX2, TX3, TX4, TX5, TX6, TX7, TX8, TX9)

    !
    !***********************************************************************
    !******************** Cumulus Downdraft Subroutine *********************
    !****************** Based on Cheng and Arakawa (1997)  ****** **********
    !************************ SUBROUTINE DDRFT  ****************************
    !*************************  October 2004  ******************************
    !***********************************************************************
    !***********************************************************************
    !************* Shrinivas.Moorthi@noaa.gov (301) 763 8000(X7233) ********
    !***********************************************************************
    !***********************************************************************
    !23456789012345678901234567890123456789012345678901234567890123456789012
    !
    !===>  TOL(K)     INPUT   TEMPERATURE            KELVIN
    !===>  QOL(K)     INPUT   SPECIFIC HUMIDITY      NON-DIMENSIONAL

    !===>  PRL(K+1)   INPUT   PRESSURE @ EDGES       MB

    !===>  K     INPUT   THE RISE & THE INDEX OF THE SUBCLOUD LAYER
    !===>  KD    INPUT   DETRAINMENT LEVEL ( 1<= KD < K )          
    !     
    !      use module_ras
    IMPLICIT NONE
    !
    !  INPUT ARGUMENTS
    !
    INTEGER      , INTENT(IN   ) :: K
    INTEGER      , INTENT(IN   ) :: KD
    REAL(kind=r8), INTENT(INOUT) :: TLA
    REAL(kind=r8), INTENT(IN   ) :: ALFIND(K)
    REAL(kind=r8), INTENT(IN   ) :: TOL(KD:K)
    REAL(kind=r8), INTENT(IN   ) :: QOL(KD:K)
    REAL(kind=r8), INTENT(IN   ) :: HOL(KD:K)
    REAL(kind=r8), INTENT(IN   ) :: PRL(KD:K+1)
    REAL(kind=r8), INTENT(IN   ) :: QST(KD:K)  
    REAL(kind=r8), INTENT(IN   ) :: HST(KD:K) 
    REAL(kind=r8), INTENT(IN   ) :: GAM(KD:K+1)
    REAL(kind=r8), INTENT(IN   ) :: GAF(KD:K+1)
!    REAL(kind=r8), INTENT(IN   ) :: HBL  ! not used
!    REAL(kind=r8), INTENT(IN   ) :: QBL ! not used
    REAL(kind=r8), INTENT(IN   ) :: QRB(KD:K)
    REAL(kind=r8), INTENT(IN   ) :: QRT(KD:K) 
    REAL(kind=r8), INTENT(INOUT) :: BUY(KD:K+1)
    INTEGER,       INTENT(IN   ) :: KBL
    INTEGER,       INTENT(IN   ) :: IDH
    REAL(kind=r8), INTENT(IN   ) :: ETA(KD:K+1)
    REAL(kind=r8), INTENT(INOUT) :: RNN(KD:K)
    REAL(kind=r8), INTENT(IN   ) :: ETAI(KD:K)
    REAL(kind=r8), INTENT(IN   ) :: ALM
    REAL(kind=r8), INTENT(IN   ) :: TRAIN
    LOGICAL,       INTENT(INOUT) :: DDFT
    REAL(kind=r8), INTENT(INOUT) :: ETD(KD:K+1)
    REAL(kind=r8), INTENT(INOUT) :: HOD(KD:K+1)
    REAL(kind=r8), INTENT(INOUT) :: QOD(KD:K+1)
    REAL(kind=r8), INTENT(INOUT) :: EVP(KD:K)
    REAL(kind=r8), INTENT(OUT  ) :: DOF
    REAL(kind=r8), INTENT(OUT  ) :: CLDFRD(KD:K)  
    REAL(kind=r8), INTENT(OUT  ) :: WCB(KD:K)
    REAL(kind=r8), INTENT(OUT  ) :: GMS(KD:K+1)
    REAL(kind=r8), INTENT(OUT  ) :: GSD(KD:K)
    REAL(kind=r8), INTENT(OUT  ) :: GHD(KD:K)


    LOGICAL :: SKPDD, SKPUP
    INTEGER :: KB1

    REAL(kind=r8) :: RNS(KD:K) 
    !
    !REAL(kind=r8) :: PRIS 
    !
    !     TEMPORARY WORK SPACE
    !
    REAL(kind=r8) :: TX1
    REAL(kind=r8) :: TX2
    REAL(kind=r8) :: TX3
    REAL(kind=r8) :: TX4
    REAL(kind=r8) :: TX5
    REAL(kind=r8) :: TX6
    !REAL(kind=r8) :: TX7
    REAL(kind=r8) :: TX8
    REAL(kind=r8) :: TX9
    LOGICAL       :: UNSAT

    !REAL(kind=r8) :: TL
    !REAL(kind=r8) :: PL
    !REAL(kind=r8) :: QL
    !REAL(kind=r8) :: QS
    !REAL(kind=r8) :: DQS
    REAL(kind=r8) :: ST1
    !REAL(kind=r8) :: HB
    !REAL(kind=r8) :: QB
    !REAL(kind=r8) :: TB  
    REAL(kind=r8) :: QQQ
    REAL(kind=r8) :: PICON
    REAL(kind=r8) :: DEL_ETA     
    REAL(kind=r8) :: TEM
    REAL(kind=r8) :: TEM1
    REAL(kind=r8) :: TEM2
    REAL(kind=r8) :: TEM3
    REAL(kind=r8) :: TEM4
    REAL(kind=r8) :: ST2   
    REAL(kind=r8) :: ERRH
    REAL(kind=r8) :: ERRW
    REAL(kind=r8) :: ERRE
    !REAL(kind=r8) :: TEM5  
    REAL(kind=r8) :: TEM6
    !REAL(kind=r8) :: HBD
    !REAL(kind=r8) :: QBD

    INTEGER  L,  N,  KD1, II                                     &
         &,       KP1,  KM1, KTEM, KK, KK1, LM1, LL, LP1                 &
         &,        ntla

    !
    INTEGER, PARAMETER :: NUMTLA=2
    !     integer, parameter :: NUMTLA=4
    REAL(KINd=r8),PARAMETER :: ERRMIN=0.0001_r8, ERRMI2=0.1_r8*ERRMIN
    !     parameter (ERRMIN=0.00001_r8, ERRMI2=0.1_r8*ERRMIN)
    !
    REAL(kind=r8) :: STLA
    REAL(kind=r8) :: CTL2
    REAL(kind=r8) :: CTL3
    !REAL(kind=r8) :: CTLA
    !REAL(kind=r8) :: VTRM
    REAL(kind=r8) :: VTPEXP 
    REAL(kind=r8) :: WCMIN
    REAL(kind=r8) :: WCBASE
    REAL(kind=r8) :: F2
    REAL(kind=r8) :: QRAF
    REAL(kind=r8) :: QRBF
    REAL(kind=r8) :: CMPOR 
    REAL(kind=r8) :: del_tla
    !REAL(kind=r8) :: sialf
    !
    REAL(kind=r8),PARAMETER :: ONPG=1.0_r8+0.5_r8, GMF=1.0_r8/ONPG, RPART=0.0_r8
    !     parameter (ONPG=1.0_r8+0.5_r8, GMF=1.0_r8/ONPG, RPART=1.0_r8)
    !     parameter (ONPG=1.0_r8+0.5_r8, GMF=1.0_r8/ONPG, RPART=0.5_r8)
    !     PARAMETER (AA1=1.0_r8, BB1=1.5_r8, CC1=1.1_r8, DD1=0.85_r8, F3=CC1, F5=2.5_r8)
    !     PARAMETER (AA1=2.0_r8, BB1=1.5_r8, CC1=1.1_r8, DD1=0.85_r8, F3=CC1, F5=2.5_r8)
    REAL(kind=r8),PARAMETER ::  AA1=1.0_r8, BB1=1.0_r8, CC1=1.0_r8, DD1=1.0_r8, F3=CC1,  F5=1.0_r8
    REAL(kind=r8),PARAMETER ::  QRMIN=1.0E-6_r8, WC2MIN=0.01_r8, GMF1=GMF/AA1, GMF5=GMF/F5
    !     parameter (QRMIN=1.0E-6_r8, WC2MIN=1.00_r8, GMF1=GMF/AA1, GMF5=GMF/F5)
    !     parameter (sialf=0.5_r8)
    !
    REAL(kind=r8),PARAMETER :: PI=3.1415926535897931_r8, PIINV=1.0_r8/PI
    INTEGER :: ITR
    !     PARAMETER (ITRMU=25, ITRMD=25, ITRMIN=7)
    INTEGER,PARAMETER :: ITRMU=25, ITRMD=25, ITRMIN=12, ITRMND=12
    !     PARAMETER (ITRMU=25, ITRMD=25, ITRMIN=12)
    !     PARAMETER (ITRMU=14, ITRMD=18, ITRMIN=7)
    !     PARAMETER (ITRMU=10, ITRMD=10, ITRMIN=5)
    REAL(kind=r8) :: QRP(KD:K+1)
    REAL(kind=r8) :: WVL(KD:K+1)
    REAL(kind=r8) :: AL2
    REAL(kind=r8) :: WVLO(KD:K+1)
    !
    REAL(kind=r8) :: RNF(KD:K)
    REAL(kind=r8) :: ROR(KD:K+1)
    REAL(kind=r8) :: STLT(KD:K)
    REAL(kind=r8) :: RNT
    REAL(kind=r8) :: RNB
    REAL(kind=r8) :: ERRQ
    REAL(kind=r8) :: RNTP
    INTEGER       :: IDW 
    INTEGER       :: IDN(K)
    INTEGER       :: idnm
    !REAL(kind=r8) :: ELM(K)
    !     real(kind=r8) EM(K*K), ELM(K)
    REAL(kind=r8) :: EDZ
    REAL(kind=r8) :: DDZ
    REAL(kind=r8) :: CE
    REAL(kind=r8) :: QHS
    REAL(kind=r8) :: FAC
    REAL(kind=r8) :: FACG
   ! REAL(kind=r8) :: ASIN
    REAL(kind=r8) :: RSUM1
    REAL(kind=r8) :: RSUM2
!    REAL(kind=r8) :: RSUM3
!    REAL(kind=r8) :: CEE

    LOGICAL       :: DDLGK
    !
    REAL(kind=r8) :: AA(KD:K,KD:K+1)
    REAL(kind=r8) :: QW(KD:K,KD:K)          
    REAL(kind=r8) :: BUD(KD:K)
    REAL(kind=r8) :: VT(2)
    REAL(kind=r8) :: VRW(2)
    REAL(kind=r8) :: TRW(2)    
    REAL(kind=r8) :: GQW(KD:K)         
    REAL(kind=r8) :: QA(3)
    REAL(kind=r8) :: WA(3)
    REAL(kind=r8) :: DOFW         
    REAL(kind=r8) :: QRPI(KD:K)
    !REAL(kind=r8) :: QRPS(KD:K)


    !    &,                    GQW(KD:K), WCB(KD:K)

    !***********************************************************************

   ! REAL(kind=r8) :: QRPF
    DOF=0.0_r8
    CLDFRD(KD:K)  =0.0_r8
    WCB(KD:K)=0.0_r8
    GMS(KD:K+1)=0.0_r8
    GSD(KD:K)=0.0_r8
    GHD(KD:K)=0.0_r8


    KB1=0

    RNS(KD:K) =0.0_r8
    TX1=0.0_r8
    TX2=0.0_r8
    TX3=0.0_r8
    TX4=0.0_r8
    TX5=0.0_r8
    TX6=0.0_r8
    TX8=0.0_r8
    TX9=0.0_r8

    ST1=0.0_r8
    QQQ=0.0_r8
    PICON=0.0_r8
    DEL_ETA=0.0_r8
    TEM=0.0_r8
    TEM1=0.0_r8
    TEM2=0.0_r8
    TEM3=0.0_r8
    TEM4=0.0_r8
    ST2=0.0_r8
    ERRH=0.0_r8
    ERRW=0.0_r8
    ERRE=0.0_r8
    TEM6=0.0_r8

    L=0  
    N=0  
    KD1=0 
    II=0      
    KP1=0  
    KM1=0 
    KTEM=0 
    KK=0 
    KK1=0 
    LM1=0 
    LL=0 
    LP1=0
    ntla=0
    STLA=0.0_r8
    CTL2=0.0_r8
    CTL3=0.0_r8
    VTPEXP =0.0_r8
    WCMIN=0.0_r8
    WCBASE=0.0_r8
    F2=0.0_r8
    QRAF=0.0_r8
    QRBF=0.0_r8
    CMPOR =0.0_r8
    del_tla=0.0_r8
    !
    ITR=0
    QRP=0.0_r8
    WVL=0.0_r8
    AL2=0.0_r8
    WVLO=0.0_r8
    !
    RNF=0.0_r8
    ROR=0.0_r8
    STLT=0.0_r8
    RNT=0.0_r8
    RNB=0.0_r8
    ERRQ=0.0_r8
    RNTP=0.0_r8
    IDW =0
    IDN(K)=0
    idnm=0
    EDZ=0.0_r8
    DDZ=0.0_r8
    CE=0.0_r8
    QHS=0.0_r8
    FAC=0.0_r8
    FACG=0.0_r8
    RSUM1=0.0_r8
    RSUM2=0.0_r8
    AA=0.0_r8
    QW=0.0_r8
    BUD=0.0_r8
    VT=0.0_r8
    VRW=0.0_r8
    TRW=0.0_r8
    GQW=0.0_r8
    QA=0.0_r8
    WA=0.0_r8
    DOFW=0.0_r8
    QRPI=0.0_r8
    !

    !     if(lprnt) print *,' K=',K,' KD=',KD,' In Downdrft'

    KD1    = KD + 1
    KP1    = K  + 1
    KM1    = K  - 1
    KB1    = KBL - 1
    !
    CMPOR  = CMB2PA / con_rd
    !
    !     VTP    = 36.34_r8*SQRT(1.2_r8)* (0.001_r8)**0.1364_r8
    VTPEXP = -0.3636_r8
    !     PIINV  = 1.0_r8 / PI
    PICON  = PI * ONEBG * 0.5_r8
    !
    !
    !     Compute Rain Water Budget of the Updraft (Cheng and Arakawa, 1997)
    !
    CLDFRD = 0.0_r8
    RNTP   = 0.0_r8
    DOF    = 0.0_r8
    ERRQ   = 10.0_r8
    RNB    = 0.0_r8
    RNT    = 0.0_r8
    TX2    = PRL(KBL)
    !
    TX1      = (PRL(KD) + PRL(KD1)) * 0.5_r8
    ROR(KD)  = CMPOR*TX1 / (TOL(KD)*(1.0_r8+con_FVirt*QOL(KD)))
    !     GMS(KD)  = VTP * ROR(KD) ** VTPEXP
    GMS(KD)  = VTP * VTPF(ROR(KD))
    !
    QRP(KD)  = QRMIN
    !
    TEM      = TOL(K) * (1.0_r8 + con_FVirt * QOL(K))
    ROR(K+1) = 0.5_r8 * CMPOR * (PRL(K+1)+PRL(K)) / TEM
    GMS(K+1) = VTP * VTPF(ROR(K+1))
    QRP(K+1) = QRMIN
    !!    BUY(KD)  = MAX(BUY(KD),ONE_M1)
    !     BUY(KD)  = MAX(BUY(KD), 0.1_r8)
    !     BUY(KD)  = MAX(BUY(KD), 0.0_r8)
    !
    kk = kbl
    DO L=KD1,K
       TEM = 0.5_r8 * (TOL(L)+TOL(L-1))                                   &
            &      * (1.0_r8 + (0.5_r8*con_FVirt) * (QOL(L)+QOL(L-1)))
       ROR(L) = CMPOR * PRL(L) / TEM
       !       GMS(L) = VTP * ROR(L) ** VTPEXP
       GMS(L) = VTP * VTPF(ROR(L))
       QRP(L) = QRMIN
       !!      BUY(L) = MAX(BUY(L),ONE_M1)
       !       BUY(L) = MAX(BUY(L), 0.1_r8)
       !       BUY(L) = MAX(BUY(L), 1.0E-5_r8)
       IF (buy(l) .LE. 0.0_r8 .AND. kk .EQ. KBL) THEN
          kk = l
          !       if (buy(l) .le. 0.0_r8) then
          !         if (buy(l-1) .gt. 0.0_r8 .and. buy(l+1) .gt. 0.0_r8) then
          !           buy(l) = 0.5_r8 * (buy(l+1) + buy(l-1))
          !         elseif (buy(l-1) .gt. 0.0_r8) then
          !           buy(l) = 0.5_r8*buy(l-1)
          !           buy(l) = 0.25_r8 * buy(l-1)
          !         else
          !            BUY(L) = 1.0E-4_r8
          !            BUY(L) = 5.0E-4_r8
          !            BUY(L) = 1.0E-5_r8
          !         endif
       ENDIF
       !       BUY(L) = MAX(BUY(L), 1.0E-4_r8)
       !       BUY(L) = MAX(BUY(L), 1.0E-5_r8)
       !       BUY(L) = MAX(BUY(L), 5.0E-4_r8)
    ENDDO
    IF (kk .NE. kbl) THEN
       DO l=kk,kbl
          buy(l) = 0.9_r8 * buy(l-1)
       ENDDO
    ENDIF
    !
    DO l=kd,k
       qrpi(l) = buy(l)
    ENDDO
    DO l=kd1,kb1
       buy(l) = 0.25_r8 * (qrpi(l-1)+qrpi(l)+qrpi(l)+qrpi(l+1))
       !       tem = 0.5_r8 * (eta(l)+eta(l+1))
       !       buy(l) = buy(l) * tem * tem
    ENDDO
    !     tem = 0.5_r8 * (eta(KD)+eta(kd1))
    !     buy(kd) = buy(kd) * tem * tem

    !
    !     CALL ANGRAD(TX1, ALM, STLA, CTL2, AL2, PI, TLA, TX2, WFN, TX3)
    tx1 = 1000.0_r8 + tx1 - prl(k+1)
    CALL ANGRAD(TX1, ALM,  AL2, TLA)
    !
    !    Following Ucla approach for rain profile
    !
    F2      = 2.0_r8*BB1*ONEBG/(PI*0.2_r8)
    WCMIN   = SQRT(WC2MIN)
    WCBASE  = WCMIN
    !
    !     del_tla = TLA * 0.2_r8
    !     del_tla = TLA * 0.25_r8
    del_tla = TLA * 0.3_r8
    TLA     = TLA - DEL_TLA
    !
    !     do ntla=1,numtla
    !
    !     if (errq .lt. 1.0_r8 .or. tla .gt. 45.0_r8) cycle
    !
    !     tla = tla + del_tla
    !     STLA = SIN(TLA*PI/180.0_r8)
    !     CTL2 = 1.0_r8 - STLA * STLA
    !
    !     if (lprnt) print *,' tla=',tla,' al2=',al2,' ptop='
    !    &,0.5_r8*(prl(kd)+prl(kd1)),' ntla=',ntla
    !     if (lprnt) print *,' buy=',(buy(l),l=kd,kbl)
    !
    !     STLA = F2     * STLA * AL2
    !     CTL2 = DD1    * CTL2
    !     CTL3 = 0.1364_r8 * CTL2
    !
    DO L=KD,K
       RNF(L)   = 0.0_r8
       RNS(L)   = 0.0_r8
       WVL(L)   = 0.0_r8
       STLT(L)  = 0.0_r8
       GQW(L)   = 0.0_r8
       QRP(L)   = QRMIN
       DO N=KD,K
          QW(N,L) = 0.0_r8
       ENDDO
    ENDDO
    !
    !-----QW(N,L) = D(W(N)*W(N))/DQR(L)
    !
    KK = KBL
    !     WVL(KK)    = WCBASE
    QW(KD,KD)  = -QRB(KD)  * GMF1
    GHD(KD)    = ETA(KD)   * ETA(KD)
    GQW(KD)    = QW(KD,KD) * GHD(KD)
    GSD(KD)    = ETAI(KD)  * ETAI(KD)
    !     GSD(KD)    = 1.0_r8 / GHD(KD)
    !
    GQW(KK)    = -  QRB(KK-1) * (GMF1+GMF1)
    !
    WCB(KK)    = WCBASE * WCBASE
    !     WVL(KK)    = WCBASE
    !     STLT(KBL)  = 1.0_r8 / WCBASE

    TX1        = WCB(KK)
    GSD(KK)    = 1.0_r8
    GHD(KK)    = 1.0_r8
    !
    TEM        = GMF1 + GMF1
    !!    TX1        = WCB(KK) + buy(kb1)*tem*qrb(kb1)
    DO L=KB1,KD1,-1
       GHD(L)  = ETA(L)  * ETA(L)
       GSD(L)  = ETAI(L) * ETAI(L)
       !        GSD(L)  = 1.0_r8 / GHD(L)
       GQW(L)  = - GHD(L) * (QRB(L-1)+QRT(L)) * TEM
       QW(L,L) = - QRT(L) * TEM
       !
       !        TX1     = TX1 + BUY(L) * TEM
       !!       TX1     = TX1 + BUY(L) * TEM * (qrb(l-1)+qrt(l)) * ghd(l)
       st1     = 0.5_r8 * (eta(l) + eta(l+1))
       TX1     = TX1 + BUY(L) * TEM * (qrb(l)+qrt(l)) * st1 * st1
       WCB(L)  = TX1 * GSD(L)
    ENDDO
    !
    TEM1        = (QRB(KD) + QRT(KD1) + QRT(KD1)) * GMF1
    GQW(KD1)    = - GHD(KD1) * TEM1
    !     QW(L,KD1)   = - QRT(KD1) * TEM
    QW(KD1,KD1) = - QRT(KD1) * TEM
    !     WCB(KD)     = (TX1 + BUY(KD)*TEM) * GSD(KD)
    st1     = 0.5_r8 * (eta(kd) + eta(kd1))
    WCB(KD)     = (TX1 + BUY(KD)*TEM*qrb(kd)*st1*st1) * GSD(KD)
    !
    DO L=KD1,KBL
       DO N=KD,L-1
          QW(N,L) = GQW(L) * GSD(N)
       ENDDO
    ENDDO
    QW(KBL,KBL) = 0.0_r8
    !
    !     WVL(KBL)    = WCBASE
    !     STLT(KBL)   = 1.0_r8 / WCBASE
    !
    !
    DO ntla=1,numtla
       !
       !     if (errq .lt. 1.0_r8 .or. tla .gt. 45.0_r8) cycle
       IF (errq .LT. 0.1_r8 .OR. tla .GT. 45.0_r8) CYCLE
       !
       tla = tla + del_tla
       STLA = SIN(TLA*PI/180.0_r8)
       CTL2 = 1.0_r8 - STLA * STLA
       !
       !     if (lprnt) print *,' tla=',tla,' al2=',al2,' ptop='
       !    &,0.5_r8*(prl(kd)+prl(kd1)),' ntla=',ntla,' f2=',f2,' stla=',stla
       !     if (lprnt) print *,' buy=',(buy(l),l=kd,kbl)
       !
       STLA = F2     * STLA * AL2
       CTL2 = DD1    * CTL2
       CTL3 = 0.1364_r8 * CTL2
       !
       DO L=KD,K
          RNF(L)   = 0.0_r8
          WVL(L)   = 0.0_r8
          STLT(L)  = 0.0_r8
          QRP(L)   = QRMIN
       ENDDO
       WVL(KBL)    = WCBASE
       STLT(KBL)   = 1.0_r8 / WCBASE
       !
       DO L=KD,K+1
          DO N=KD,K
             AA(N,L) = 0.0_r8
          ENDDO
       ENDDO
       !
       SKPUP = .FALSE.
       !
       DO ITR=1,ITRMU               ! Rain Profile Iteration starts!
          IF (.NOT. SKPUP) THEN
             wvlo = wvl
             !
             !-----CALCULATING THE VERTICAL VELOCITY
             !
             TX1      = 0.0_r8
             QRPI(KBL) = 1.0_r8 / QRP(KBL)
             DO L=KB1,KD,-1
                TX1     = TX1    + QRP(L+1) * GQW(L+1)
                ST1     = WCB(L) + QW(L,L)  * QRP(L)                        &
                     &                       + TX1      * GSD(L)
                !           if (st1 .gt. 0.0_r8) then
                IF (st1 .GT. wc2min) THEN
                   !             WVL(L)  = SQRT(ST1)
                   WVL(L)  = 0.5_r8 * (SQRT(ST1) + WVL(L))
                   !             if (itr .eq. 1) wvl(l) = wvl(l) * 0.25_r8
                ELSE
                   !     if (lprnt)  print *,' l=',l,' st1=',st1,' wcb=',wcb(l),' qw='
                   !    &,qw(l,l),' qrp=',qrp(l),' tx1=',tx1,' gsd=',gsd(l),' ite=',itr
                   !             wvl(l) = 0.5_r8*(wcmin+wvl(l))
                   wvl(l) = 0.5_r8 * (wvl(l) + wvl(l+1))
                   qrp(l) = 0.5_r8 * ((wvl(l)*wvl(l)-wcb(l)-tx1*gsd(l))/qw(l,l) &
                        &                     + qrp(l))
                   !!            wvl(l) = 0.5_r8 * (wvl(l) + wvl(l+1))
                ENDIF
                !           wvl(l)  = 0.5_r8 * (wvl(l) + wvlo(l))
                !           WVL(L)  = SQRT(MAX(ST1,WC2MIN))
                wvl(l)  = MAX(wvl(l), wcbase)
                STLT(L) = 1.0_r8 / WVL(L)
                QRPI(L) = 1.0_r8 / QRP(L)
             ENDDO
             !         qrps = qrp
             !         do l=kd1,kb1
             !           qrp(l) = 0.25_r8 * (qrps(l-1)+qrps(l)+qrps(l)+qrps(l+1))
             !           qrpi(l) = 1.0_r8 / qrp(l)
             !         enddo
             !         qrpi(kd) = 1.0_r8 / qrp(kd)
             !
             !     if (lprnt) then
             !     print *,' ITR=',ITR,' ITRMU=',ITRMU
             !     print *,' WVL=',(WVL(L),L=KD,KBL)
             !     print *,' qrp=',(qrp(L),L=KD,KBL)
             !     print *,' qrpi=',(qrpi(L),L=KD,KBL)
             !     print *,' rnf=',(rnf(L),L=KD,KBL)
             !     endif
             !
             !-----CALCULATING TRW, VRW AND OF
             !
             !         VT(1)   = GMS(KD) * QRP(KD)**0.1364_r8
             VT(1)   = GMS(KD) * QRPF(QRP(KD))
             TRW(1)  = ETA(KD) * QRP(KD) * STLT(KD)
             TX6     = TRW(1) * VT(1)
             VRW(1)  = F3*WVL(KD) - CTL2*VT(1)
             BUD(KD) = STLA * TX6 * QRB(KD) * 0.5_r8
             RNF(KD) = BUD(KD)
             DOF     = 1.1364_r8 * BUD(KD) * QRPI(KD)
             DOFW    = -BUD(KD) * STLT(KD)
             !
             RNT     = TRW(1) * VRW(1)
             TX2     = 0.0_r8
             TX4     = 0.0_r8
             RNB     = RNT
             TX1     = 0.5_r8
             TX8     = 0.0_r8
             !
             IF (RNT .GE. 0.0_r8) THEN
                TX3 = (RNT-CTL3*TX6) * QRPI(KD)
                TX5 = CTL2 * TX6 * STLT(KD)
             ELSE
                TX3 = 0.0_r8
                TX5 = 0.0_r8
                RNT = 0.0_r8
                RNB = 0.0_r8
             ENDIF
             !
             DO L=KD1,KB1
                KTEM    = MAX(L-2, KD)
                LL      = L - 1
                !
                !           VT(2)   = GMS(L) * QRP(L)**0.1364_r8
                VT(2)   = GMS(L) * QRPF(QRP(L))
                TRW(2)  = ETA(L) * QRP(L) * STLT(L)
                VRW(2)  = F3*WVL(L) - CTL2*VT(2)
                QQQ     = STLA * TRW(2) * VT(2)
                ST1     = TX1  * QRB(LL)
                BUD(L)  = QQQ * (ST1 + QRT(L))
                !
                QA(2)   = DOF
                WA(2)   = DOFW
                DOF     = 1.1364_r8 * BUD(L) * QRPI(L)
                DOFW    = -BUD(L) * STLT(L)
                !
                RNF(LL) = RNF(LL) + QQQ * ST1
                RNF(L)  =           QQQ * QRT(L)
                !
                TEM3    = VRW(1) + VRW(2)
                TEM4    = TRW(1) + TRW(2)
                !
                TX6     = .25_r8 * TEM3 * TEM4
                TEM4    = TEM4 * CTL3
                !
                !-----BY QR ABOVE
                !
                !           TEM1    = .25_r8*(TRW(1)*TEM3 - TEM4*VT(1))*TX7
                TEM1    = .25_r8*(TRW(1)*TEM3 - TEM4*VT(1))*QRPI(LL)
                ST1     = .25_r8*(TRW(1)*(CTL2*VT(1)-VRW(2))                   &
                     &                  * STLT(LL) + F3*TRW(2))
                !-----BY QR BELOW
                TEM2    = .25_r8*(TRW(2)*TEM3 - TEM4*VT(2))*QRPI(L)
                ST2     = .25_r8*(TRW(2)*(CTL2*VT(2)-VRW(1))                   &
                     &                 * STLT(L)  + F3*TRW(1))
                !
                !      From top to  the KBL-2 layer
                !
                QA(1)   = TX2
                QA(2)   = QA(2) + TX3 - TEM1
                QA(3)   = -TEM2
                !
                WA(1)   = TX4
                WA(2)   = WA(2) + TX5 - ST1
                WA(3)   = -ST2
                !
                TX2     = TEM1
                TX3     = TEM2
                TX4     = ST1
                TX5     = ST2
                !
                VT(1)   = VT(2)
                TRW(1)  = TRW(2)
                VRW(1)  = VRW(2)
                !
                IF (WVL(KTEM) .EQ. WCMIN) WA(1) = 0.0_r8
                IF (WVL(LL)   .EQ. WCMIN) WA(2) = 0.0_r8
                IF (WVL(L)    .EQ. WCMIN) WA(3) = 0.0_r8
                DO N=KTEM,KBL
                   AA(LL,N) = (WA(1)*QW(KTEM,N) * STLT(KTEM)                 &
                        &                 +  WA(2)*QW(LL,N)   * STLT(LL)                   &
                        &                 +  WA(3)*QW(L,N)    * STLT(L) ) * 0.5_r8
                ENDDO
                AA(LL,KTEM) = AA(LL,KTEM) + QA(1)
                AA(LL,LL)   = AA(LL,LL)   + QA(2)
                AA(LL,L)    = AA(LL,L)    + QA(3)
                BUD(LL)     = (TX8 + RNN(LL)) * 0.5_r8                         &
                     &                    - RNB + TX6 - BUD(LL)
                AA(LL,KBL+1) = BUD(LL)
                RNB = TX6
                TX1 = 1.0_r8
                TX8 = RNN(LL)
             ENDDO
             L  = KBL
             LL = L - 1
             !         VT(2)   = GMS(L) * QRP(L)**0.1364_r8
             VT(2)   = GMS(L) * QRPF(QRP(L))
             TRW(2)  = ETA(L) * QRP(L) * STLT(L)
             VRW(2)  = F3*WVL(L) - CTL2*VT(2)
             ST1     = STLA * TRW(2) * VT(2) * QRB(LL)
             BUD(L)  = ST1

             QA(2)   = DOF
             WA(2)   = DOFW
             DOF     = 1.1364_r8 * BUD(L) * QRPI(L)
             DOFW    = -BUD(L) * STLT(L)
             !
             RNF(LL) = RNF(LL) + ST1
             !
             TEM3    = VRW(1) + VRW(2)
             TEM4    = TRW(1) + TRW(2)
             !
             TX6     = .25_r8 * TEM3 * TEM4
             TEM4    = TEM4 * CTL3
             !
             !-----BY QR ABOVE
             !
             TEM1    = .25_r8*(TRW(1)*TEM3 - TEM4*VT(1))*QRPI(LL)
             ST1     = .25_r8*(TRW(1)*(CTL2*VT(1)-VRW(2))                     &
                  &                * STLT(LL) + F3*TRW(2))
             !-----BY QR BELOW
             TEM2    = .25_r8*(TRW(2)*TEM3 - TEM4*VT(2))*QRPI(L)
             ST2     = .25_r8*(TRW(2)*(CTL2*VT(2)-VRW(1))                     &
                  &                 * STLT(L)  + F3*TRW(1))
             !
             !      For the layer next to the top of the boundary layer
             !
             QA(1)   = TX2
             QA(2)   = QA(2) + TX3 - TEM1
             QA(3)   = -TEM2
             !
             WA(1)   = TX4
             WA(2)   = WA(2) + TX5 - ST1
             WA(3)   = -ST2
             !
             TX2     = TEM1
             TX3     = TEM2
             TX4     = ST1
             TX5     = ST2
             !
             IDW     = MAX(L-2, KD)
             !
             IF (WVL(IDW) .EQ. WCMIN) WA(1) = 0.0_r8
             IF (WVL(LL)  .EQ. WCMIN) WA(2) = 0.0_r8
             IF (WVL(L)   .EQ. WCMIN) WA(3) = 0.0_r8
             !
             KK = IDW
             DO N=KK,L
                AA(LL,N) = (WA(1)*QW(KK,N) * STLT(KK)                       &
                     &               +  WA(2)*QW(LL,N) * STLT(LL)                       &
                     &               +  WA(3)*QW(L,N)  * STLT(L) ) * 0.5_r8

             ENDDO
             !
             AA(LL,IDW) = AA(LL,IDW) + QA(1)
             AA(LL,LL)  = AA(LL,LL)  + QA(2)
             AA(LL,L)   = AA(LL,L)   + QA(3)
             BUD(LL)    = (TX8+RNN(LL)) * 0.5_r8 - RNB + TX6 - BUD(LL)
             !
             AA(LL,L+1) = BUD(LL)
             !
             RNB        = TRW(2) * VRW(2)
             !
             !      For the top of the boundary layer
             !
             IF (RNB .LT. 0.0_r8) THEN
                KK    = KBL
                TEM   = VT(2) * TRW(2)
                QA(2) = (RNB - CTL3*TEM) * QRPI(KK)
                WA(2) = CTL2 * TEM * STLT(KK)
             ELSE
                RNB   = 0.0_r8
                QA(2) = 0.0_r8
                WA(2) = 0.0_r8
             ENDIF
             !
             QA(1) = TX2
             QA(2) = DOF + TX3 - QA(2)
             QA(3) = 0.0_r8
             !
             WA(1) = TX4
             WA(2) = DOFW + TX5 - WA(2)
             WA(3) = 0.0_r8
             !
             KK = KBL
             IF (WVL(KK-1) .EQ. WCMIN) WA(1) = 0.0_r8
             IF (WVL(KK)   .EQ. WCMIN) WA(2) = 0.0_r8
             !
             DO II=1,2
                N = KK + II - 2
                AA(KK,N) = (WA(1)*QW(KK-1,N) * STLT(KK-1)                  &
                     &                +  WA(2)*QW(KK,N)   * STLT(KK)) * 0.5_r8
             ENDDO
             FAC = 0.5_r8
             LL  = KBL
             L   = LL + 1
             LM1 = LL - 1
             AA(LL,LM1)  = AA(LL,LM1) + QA(1)
             AA(LL,LL)   = AA(LL,LL)  + QA(2)
             BUD(LL)     = 0.5_r8*RNN(LM1) - TX6 + RNB - BUD(LL)
             AA(LL,LL+1) = BUD(LL)
             !
             !-----SOLVING THE BUDGET EQUATIONS FOR DQR
             !
             DO L=KD1,KBL
                LM1  = L - 1
                UNSAT = ABS(AA(LM1,LM1)) .LT. ABS(AA(L,LM1))
                DO  N=LM1,KBL+1
                   IF (UNSAT) THEN
                      TX1       = AA(LM1,N)
                      AA(LM1,N) = AA(L,N)
                      AA(L,N)   = TX1
                   ENDIF
                ENDDO
                TX1 = AA(L,LM1) / AA(LM1,LM1)
                DO  N=L,KBL+1
                   AA(L,N) = AA(L,N) - TX1 * AA(LM1,N)
                ENDDO
             ENDDO
             !
             !-----BACK SUBSTITUTION AND CHECK IF THE SOLUTION CONVERGES
             !
             KK = KBL
             KK1 = KK + 1
             AA(KK,KK1) = AA(KK,KK1) / AA(KK,KK)      !   Qr correction !
             TX2        = ABS(AA(KK,KK1)) * QRPI(KK)  !   Error Measure !
             !     if (lprnt) print *,' tx2a=',tx2,' aa1=',aa(kk,kk1)
             !    &,' qrpi=',qrpi(kk)
             !
             KK = KBL + 1
             DO L=KB1,KD,-1
                LP1   = L + 1
                TX1  = 0.0_r8
                DO N=LP1,KBL
                   TX1  = TX1 + AA(L,N) * AA(N,KK)
                ENDDO
                AA(L,KK) = (AA(L,KK) - TX1) / AA(L,L)       ! Qr correction !
                TX2      = MAX(TX2, ABS(AA(L,KK))*QRPI(L))  ! Error Measure !
                !     if (lprnt) print *,' tx2b=',tx2,' aa1=',aa(l,kk)
                !    &,' qrpi=',qrpi(l),' L=',L
             ENDDO
             !
             !         tem = 0.5_r8
             IF (tx2 .GT. 1.0_r8 .AND. ABS(errq-tx2) .GT. 0.1_r8) THEN
                tem = 0.5_r8
                !!        elseif (tx2 .lt. 0.1_r8) then
                !!          tem = 1.2_r8
             ELSE
                tem = 1.0_r8
             ENDIF
             !
             DO L=KD,KBL
                !            QRP(L) = MAX(QRP(L)+AA(L,KBL+1), QRMIN)
                QRP(L) = MAX(QRP(L)+AA(L,KBL+1)*tem, QRMIN)
             ENDDO
             !
             !     if (lprnt) print *,' itr=',itr,' tx2=',tx2
             IF (ITR .LT. ITRMIN) THEN
                TEM = ABS(ERRQ-TX2) 
                IF (TEM .GE. ERRMI2 .AND. TX2 .GE. ERRMIN) THEN 
                   ERRQ  = TX2                              ! Further iteration !
                ELSE 
                   SKPUP = .TRUE.                           ! Converges      !
                   ERRQ  = 0.0_r8                              ! Rain profile exists!
                   !     print *,' here1',' tem=',tem,' tx2=',tx2,' errmi2=',
                   !    *errmi2,' errmin=',errmin
                ENDIF
             ELSE
                TEM = ERRQ - TX2
                !            IF (TEM .LT. ZERO .AND. ERRQ .GT. 0.1_r8) THEN
                IF (TEM .LT. ZERO .AND. ERRQ .GT. 0.5_r8) THEN
                   !            IF (TEM .LT. ZERO .and.                                    &
                   !    &          (ntla .lt. numtla .or. ERRQ .gt. 0.5_r8)) THEN
                   !     if (lprnt) print *,' tx2=',tx2,' errq=',errq,' tem=',tem
                   SKPUP = .TRUE.                           ! No convergence !
                   ERRQ = 10.0_r8                              ! No rain profile!
!!!!         ELSEIF (ABS(TEM).LT.ERRMI2 .OR. TX2.LT.ERRMIN) THEN
                ELSEIF (TX2.LT.ERRMIN) THEN
                   SKPUP = .TRUE.                           ! Converges      !
                   ERRQ = 0.0_r8                               ! Rain profile exists!
                   !     print *,' here2'
                ELSEIF (tem .LT. zero .AND. errq .LT. 0.1_r8) THEN
                   skpup = .TRUE.
                   !              if (ntla .eq. numtla .or. tem .gt. -0.003) then
                   errq  = 0.0_r8
                   !              else
                   !                errq = 10.0
                   !              endif
                ELSE
                   ERRQ = TX2                               ! Further iteration !
                   !     if (lprnt) print *,' itr=',itr,' errq=',errq
                   !              if (itr .eq. itrmu .and. ERRQ .GT. ERRMIN*10             &
                   !    &            .and. ntla .eq. 1) ERRQ = 10.0 
                ENDIF
             ENDIF
             !
             !         if (lprnt) print *,' ERRQ=',ERRQ

          ENDIF                                           ! SKPUP  ENDIF!
          !
       ENDDO                                          ! End of the ITR Loop!!
       !     enddo                                          ! End of ntla loop
       !
       !     if(lprnt) then
       !       print *,' QRP=',(QRP(L),L=KD,KBL)
       !       print *,'RNF=',(RNF(L),L=KD,KBL),' RNT=',RNT,' RNB=',RNB
       !    &,' errq=',errq
       !     endif
       !
       IF (ERRQ .LT. 0.1_r8) THEN
          DDFT = .TRUE.
          RNB  = - RNB
          !    do l=kd1,kb1-1
          !      if (wvl(l)-wcbase .lt. 1.0E-9) ddft = .false.
          !    enddo
       ELSE
          DDFT = .FALSE.
       ENDIF
       !
       !     Caution !! Below is an adjustment to rain flux to maintain
       !                conservation of precip!
       !
       IF (DDFT) THEN
          TX1 = 0.0_r8
          DO L=KD,KB1
             TX1 = TX1 + RNF(L)
          ENDDO
          !     if (lprnt) print *,' tx1+rnt+rnb=',tx1+rnt+rnb, ' train=',train
          TX1 = TRAIN / (TX1+RNT+RNB)
          IF (ABS(TX1-1.0_r8) .LT. 0.2_r8) THEN
             RNT = MAX(RNT*TX1,ZERO)
             RNB = RNB * TX1
          ELSE
             DDFT = .FALSE.
             ERRQ = 10.0_r8
          ENDIF
       ENDIF
    ENDDO                                          ! End of ntla loop
    !
    DOF = 0.0_r8
    IF (.NOT. DDFT) RETURN     ! Rain profile did not converge!
    !

    DO L=KD,KB1
       RNF(L) = RNF(L) * TX1

    ENDDO
    !     if (lprnt) print *,' TRAIN=',TRAIN
    !     if (lprnt) print *,' RNF=',RNF
    !
    !     Adjustment is over
    !
    !     Downdraft
    !
    DO L=KD,K
       WCB(L) = 0.0_r8
    ENDDO
    !
    SKPDD = .NOT. DDFT
    !
    ERRQ  = 10.0_r8
    IF (.NOT. SKPDD) THEN
       !
       !     Calculate Downdraft Properties
       !

       KK = MAX(KB1,KD1)
       DO L=KK,K
          STLT(L) = STLT(L-1)
       ENDDO
       TEM1 = 1.0_r8 / BB1
       !
       DO L=KD,K
          IF (L .LE. KBL) THEN
             TEM     = STLA * TEM1
             STLT(L) = ETA(L) * STLT(L) * TEM / ROR(L)
          ELSE
             STLT(L) = 0.0_r8
          ENDIF
       ENDDO
       !       if (lprnt) print *,' STLT=',stlt

       rsum1 = 0.0_r8
       rsum2 = 0.0_r8

       !
       IDN      = 99
       DO L=KD,K+1
          ETD(L)  = 0.0_r8
          WVL(L)  = 0.0_r8
          !         QRP(L)  = 0.0_r8
       ENDDO
       DO L=KD,K
          EVP(L)   = 0.0_r8
          BUY(L)   = 0.0_r8
          QRP(L+1) = 0.0_r8
       ENDDO
       HOD(KD)  = HOL(KD)
       QOD(KD)  = QOL(KD)
       TX1      = 0.0_r8                               ! sigma at the top
!!!     TX1      = STLT(KD)*QRB(KD)*ONE              ! sigma at the top
       !       TX1      = MIN(STLT(KD)*QRB(KD)*ONE, ONE)    ! sigma at the top
       !       TX1      = MIN(STLT(KD)*QRB(KD)*0.5_r8, ONE)    ! sigma at the top
       RNTP     = 0.0_r8
       TX5      = TX1
       QA(1)    = 0.0_r8
       !     if(lprnt) print *,' stlt=',stlt(kd),' qrb=',qrb(kd)
       !    *,' tx1=',tx1,' ror=',ror(kd),' gms=',gms(kd),' rpart=',rpart
       !    *,' rnt=',rnt
       !
       !       Here we assume RPART of detrained rain RNT goes to Pd
       !
       IF (RNT .GT. 0.0_r8) THEN
          IF (TX1 .GT. 0.0_r8) THEN
             QRP(KD) = (RPART*RNT / (ROR(KD)*TX1*GMS(KD)))               &
                  &                                          ** (1.0_r8/1.1364_r8)
          ELSE
             tx1 = RPART*RNT / (ROR(KD)*GMS(KD)*QRP(KD)**1.1364_r8)
          ENDIF
          RNTP    = (1.0_r8 - RPART) * RNT
          BUY(KD) = - ROR(KD) * TX1 * QRP(KD)
       ELSE
          QRP(KD) = 0.0_r8
       ENDIF
       !
       !     L-loop for the downdraft iteration from KD1 to K+1 (bottom surface)
       !
       !     BUD(KD) = ROR(KD)
       idnm = 1
       DO L=KD1,K+1

          QA(1) = 0.0_r8
          ddlgk = idn(idnm) .EQ. 99
          IF (.NOT. ddlgk) CYCLE
          IF (L .LE. K) THEN
             ST1   = 1.0_r8 - ALFIND(L)
             WA(1) = ALFIND(L)*HOL(L-1) + ST1*HOL(L)
             WA(2) = ALFIND(L)*QOL(L-1) + ST1*QOL(L)
             WA(3) = ALFIND(L)*TOL(L-1) + ST1*TOL(L)
             QA(2) = ALFIND(L)*HST(L-1) + ST1*HST(L)
             QA(3) = ALFIND(L)*QST(L-1) + ST1*QST(L)
          ELSE
             WA(1) = HOL(K)
             WA(2) = QOL(K)
             WA(3) = TOL(K)
             QA(2) = HST(K)
             QA(3) = QST(K)
          ENDIF
          !
          FAC = 2.0_r8
          IF (L .EQ. KD1) FAC = 1.0_r8

          FACG    = FAC * 0.5_r8 * GMF5     !  12/17/97
          !
          !         DDLGK   =  IDN(idnm) .EQ. 99
          BUD(KD) = ROR(L)

          !         IF (DDLGK) THEN
          TX1    = TX5
          WVL(L) = MAX(WVL(L-1),ONE_M1)

          QRP(L) = MAX(QRP(L-1),QRP(L))
          !
          !           VT(1)  = GMS(L-1) * QRP(L-1) ** 0.1364_r8
          VT(1)  = GMS(L-1) * QRPF(QRP(L-1))
          RNT    = ROR(L-1) * (WVL(L-1)+VT(1))*QRP(L-1)
          !     if(lprnt) print *,' l=',l,' qa=',qa(1), ' tx1RNT=',RNT*tx1,
          !    *' wvl=',wvl(l-1)
          !    *,' qrp=',qrp(l-1),' tx5=',tx5,' tx1=',tx1,' rnt=',rnt

          !

          !           TEM    = MAX(ALM, 2.5E-4) * MAX(ETA(L), 1.0)
          TEM    = MAX(ALM,ONE_M6) * MAX(ETA(L), ONE)
          !           TEM    = MAX(ALM, 1.0E-5) * MAX(ETA(L), 1.0)
          TRW(1) = PICON*TEM*(QRB(L-1)+QRT(L-1))
          TRW(2) = 1.0_r8 / TRW(1)
          !
          VRW(1) = 0.5_r8 * (GAM(L-1) + GAM(L))
          VRW(2) = 1.0_r8 / (VRW(1) + VRW(1))
          !
          TX4    =  (QRT(L-1)+QRB(L-1))*(ONEBG*FAC*500.00_r8*EKNOB)
          !
          DOFW   = 1.0_r8 / (WA(3) * (1.0_r8 + con_FVirt*WA(2)))      !  1.0_r8 / TVbar!
          !
          ETD(L) = ETD(L-1)
          HOD(L) = HOD(L-1)
          QOD(L) = QOD(L-1)
          !
          ERRQ   = 10.0_r8

          !
          IF (L .LE. KBL) THEN
             TX3 = STLT(L-1) * QRT(L-1) * (0.5_r8*FAC)
             TX8 = STLT(L)   * QRB(L-1) * (0.5_r8*FAC)
             TX9 = TX8 + TX3
          ELSE
             TX3 = 0.0_r8
             TX8 = 0.0_r8
             TX9 = 0.0_r8
          ENDIF
          !
          TEM  = WVL(L-1) + VT(1)
          IF (TEM .GT. 0.0_r8) THEN
             TEM1 = 1.0_r8 / (TEM*ROR(L-1))
             TX3 = VT(1) * TEM1 * ROR(L-1) * TX3
             TX6 = TX1 * TEM1
          ELSE
             TX6 = 1.0_r8
          ENDIF
          !         ENDIF
          !
          IF (L .EQ. KD1) THEN
             IF (RNT .GT. 0.0_r8) THEN
                TEM    = MAX(QRP(L-1),QRP(L))
                WVL(L) = TX1 * TEM * QRB(L-1)*(FACG*5.0_r8)
             ENDIF
             WVL(L) = MAX(ONE_M2, WVL(L))
             TRW(1) = TRW(1) * 0.5_r8
             TRW(2) = TRW(2) + TRW(2)
          ELSE
             IF (DDLGK) EVP(L-1) = EVP(L-2)
          ENDIF
          !
          !       No downdraft above level IDH
          !

          IF (L .LT. IDH) THEN

             ETD(L)   = 0.0_r8
             HOD(L)   = WA(1)
             QOD(L)   = WA(2)
             EVP(L-1) = 0.0_r8
             WVL(L)   = 0.0_r8
             QRP(L)   = 0.0_r8
             BUY(L)   = 0.0_r8
             TX5      = TX9
             ERRQ     = 0.0_r8
             RNTP     = RNTP + RNT * TX1
             RNT      = 0.0_r8
             WCB(L-1) = 0.0_r8
          ENDIF
          !         BUD(KD) = ROR(L)
          !
          !       Iteration loop for a given level L begins
          !
          !         if (lprnt) print *,' tx8=',tx8,' tx9=',tx9,' tx5=',tx5
          !    &,                      ' tx1=',tx1
          DO ITR=1,ITRMD
             !
             !           UNSAT =  DDLGK .AND. (ERRQ .GT. ERRMIN)
             UNSAT =  ERRQ .GT. ERRMIN
             IF (UNSAT) THEN
                !
                !             VT(1)  = GMS(L) * QRP(L) ** 0.1364
                VT(1)  = GMS(L) * QRPF(QRP(L))
                TEM    =  WVL(L) + VT(1)
                !
                IF (TEM .GT. 0.0_r8) THEN
                   ST1    = ROR(L) * TEM * QRP(L) + RNT
                   IF (ST1 .NE. 0.0_r8) ST1 = 2.0_r8 * EVP(L-1) / ST1
                   TEM1   = 1.0_r8 / (TEM*ROR(L))
                   TEM2   = VT(1) * TEM1 * ROR(L) * TX8
                ELSE
                   TEM1   = 0.0_r8
                   TEM2   = TX8
                   ST1    = 0.0_r8
                ENDIF
                !     if (lprnt) print *,' st1=',st1,' tem=',tem,' ror=',ror(l)
                !    &,' qrp=',qrp(l),' rnt=',rnt,' ror1=',ror(l-1),' wvl=',wvl(l)
                !    &,' wvl1=',wvl(l-1),' tem2=',tem2,' vt=',vt(1),' tx3=',tx3
                !
                st2 = tx5
                TEM = ROR(L)*WVL(L) - ROR(L-1)*WVL(L-1)
                IF (tem .GT. 0.0_r8) THEN
                   TX5 = (TX1 - ST1 + TEM2 + TX3)/(1.0_r8+tem*tem1)
                ELSE
                   TX5 = TX1 - tem*tx6 - ST1 + TEM2 + TX3
                ENDIF
                TX5   = MAX(TX5,ZERO)
                tx5 = 0.5_r8 * (tx5 + st2)
                !
                !             qqq = 1.0_r8 + tem * tem1 * (1.0_r8 - sialf)
                !
                !             if (qqq .gt. 0.0_r8) then
                !               TX5   = (TX1 - sialf*tem*tx6 - ST1 + TEM2 + TX3) / qqq
                !             else
                !               TX5   = (TX1 - tem*tx6 - ST1 + TEM2 + TX3)
                !             endif
                !
                !     if(lprnt) print *,' tx51=',tx5,' tx1=',tx1,' st1=',st1,' tem2='
                !     if(tx5 .le. 0.0_r8 .and. l .gt. kd+2)
                !    * print *,' tx51=',tx5,' tx1=',tx1,' st1=',st1,' tem2='
                !    *,tem2,' tx3=',tx3,' tem=',tem,' tem1=',tem1,' wvl=',wvl(l-1),
                !    &wvl(l),' l=',l,' itr=',itr,' evp=',evp(l-1),' vt=',vt(1)
                !    *,' qrp=',qrp(l),' rnt=',rnt,' kd=',kd
                !     if (lprnt) print *,' etd=',etd(l),' wvl=',wvl(l)
                !    &,' trw=',trw(1),trw(2),' ror=',ror(l),' wa=',wa


                !
                TEM1   = ETD(L)
                ETD(L) = ROR(L) * TX5 * MAX(WVL(L),ZERO)
                !
                IF (etd(l) .GT. 0.0_r8) etd(l) = 0.5_r8 * (etd(l) + tem1)
                !

                DEL_ETA = ETD(L) - ETD(L-1)

                !               TEM       = DEL_ETA * TRW(2)
                !               TEM2      = MAX(MIN(TEM, 1.0_r8), -1.0_r8)
                !               IF (ABS(TEM) .GT. 1.0_r8 .AND. ETD(L) .GT. 0.0_r8 ) THEN
                !                 DEL_ETA = TEM2 * TRW(1)
                !                 ETD(L)  = ETD(L-1) + DEL_ETA
                !               ENDIF
                !               IF (WVL(L) .GT. 0.0_r8) TX5 = ETD(L) / (ROR(L)*WVL(L))
                !
                ERRE  = ETD(L) - TEM1
                !
                tem  = MAX(ABS(del_eta), trw(1))
                tem2 = del_eta / tem
                TEM1 = SQRT(MAX((tem+DEL_ETA)*(tem-DEL_ETA),ZERO))
                !               TEM1 = SQRT(MAX((TRW(1)+DEL_ETA)*(TRW(1)-DEL_ETA),0.0_r8))

                EDZ  = (0.5_r8 + ASIN(TEM2)*PIINV)*DEL_ETA + TEM1*PIINV

                DDZ   = EDZ - DEL_ETA
                WCB(L-1) = ETD(L) + DDZ
                !
                TEM1  = HOD(L)
                IF (DEL_ETA .GT. 0.0_r8) THEN
                   QQQ    = 1.0_r8 / (ETD(L) + DDZ)
                   HOD(L) = (ETD(L-1)*HOD(L-1) + DEL_ETA*HOL(L-1)          &
                        &                                            + DDZ*WA(1)) * QQQ
                   QOD(L) = (ETD(L-1)*QOD(L-1) + DEL_ETA*QOL(L-1)          &
                        &                                            + DDZ*WA(2)) * QQQ
                ELSEIF((ETD(L-1) + EDZ) .GT. 0.0_r8) THEN
                   QQQ    = 1.0_r8 / (ETD(L-1) + EDZ)
                   HOD(L) = (ETD(L-1)*HOD(L-1) + EDZ*WA(1)) * QQQ
                   QOD(L) = (ETD(L-1)*QOD(L-1) + EDZ*WA(2)) * QQQ
                ENDIF
                ERRH  = HOD(L) - TEM1
                ERRQ  = ABS(ERRH/HOD(L))  + ABS(ERRE/MAX(ETD(L),ONE_M5))
                !     if (lprnt) print *,' ERRQP=',errq,' errh=',errh,' hod=',hod(l)
                !    &,' erre=',erre,' etd=',etd(l),' del_eta=',del_eta
                DOF   = DDZ
                VT(2) = QQQ

                !
                DDZ  = DOF
                TEM4 = QOD(L)
                TEM1 = VRW(1)
                !
                QHS  = QA(3) + 0.5_r8 * (GAF(L-1)+GAF(L))                    &
                     &                           * (HOD(L)-QA(2))
                !
                !                                           First iteration       !
                !
                ST2  = PRL(L) * (QHS + TEM1 * (QHS-QOD(L)))
                TEM2 = ROR(L) * QRP(L)
                CALL QRABF(TEM2,QRAF,QRBF)
                TEM6 = TX5 * (1.6_r8 + 124.9_r8 * QRAF) * QRBF * TX4
                !
                CE   = TEM6 * ST2 / ((5.4E5_r8*ST2 + 2.55E6_r8)*(ETD(L)+DDZ))
                !
                TEM2   = - ((1.0_r8+TEM1)*(QHS+CE) + TEM1*QOD(L))
                TEM3   = (1.0_r8 + TEM1) * QHS * (QOD(L)+CE)
                TEM    = MAX(TEM2*TEM2 - 4.0_r8*TEM1*TEM3,ZERO)
                QOD(L) = MAX(TEM4, (- TEM2 - SQRT(TEM)) * VRW(2))
                !

                !
                !                                            second iteration   !
                !
                ST2  = PRL(L) * (QHS + TEM1 * (QHS-QOD(L)))
                CE   = TEM6 * ST2 / ((5.4E5_r8*ST2 + 2.55E6_r8)*(ETD(L)+DDZ))
                !             CEE  = CE * (ETD(L)+DDZ)
                !


                TEM2   = - ((1.0_r8+TEM1)*(QHS+CE) + TEM1*tem4)
                TEM3   = (1.0_r8 + TEM1) * QHS * (tem4+CE)
                TEM    = MAX(TEM2*TEM2 - 4.0_r8*TEM1*TEM3,ZERO)
                QOD(L) = MAX(TEM4, (- TEM2 - SQRT(TEM)) * VRW(2))
                !                                              Evaporation in Layer L-1
                !

                EVP(L-1) = (QOD(L)-TEM4) * (ETD(L)+DDZ)
                !                                              Calculate Pd (L+1/2)
                QA(1)    = TX1*RNT + RNF(L-1) - EVP(L-1)
                !
                !     if(lprnt) print *,' etd=',etd(l),' tx5=',tx5,' rnt=',rnt
                !    *,' rnf=',rnf(l-1),' evp=',evp(l-1),' itr=',itr,' L=',L

                !
                IF (qa(1) .GT. 0.0_r8) THEN
                   IF (ETD(L) .GT. 0.0_r8) THEN
                      TEM    = QA(1) / (ETD(L)+ROR(L)*TX5*VT(1))
                      QRP(L) = MAX(TEM,ZERO)
                   ELSEIF (TX5 .GT. 0.0_r8) THEN
                      QRP(L) = (MAX(ZERO,QA(1)/(ROR(L)*TX5*GMS(L))))           &
                           &                                          ** (1.0_r8/1.1364_r8)
                   ELSE
                      QRP(L) = 0.0_r8
                   ENDIF
                ELSE
                   qrp(l) = 0.5_r8 * qrp(l)
                ENDIF
                !                                              Compute Buoyancy
                TEM1   = WA(3)+(HOD(L)-WA(1)-con_hvap*(QOD(L)-WA(2)))         &
                     &                                                  * (1.0_r8/con_cp)
                !             if (lprnt) print *,' tem1=',tem1,' wa3=',wa(3),' hod='
                !    &,hod(l),' wa1=',wa(1),' qod=',qod(l),' wa2=',wa(2),' con_hvap=',con_hvap
                !    &,' cmpor=',cmpor,' dofw=',dofw,' prl=',prl(l),' qrp=',qrp(l)
                TEM1   = TEM1 * (1.0_r8 + con_FVirt*QOD(L))
                ROR(L) = CMPOR * PRL(L) / TEM1
                TEM1   = TEM1 * DOFW
!!!           TEM1   = TEM1 * (1.0_r8 + con_FVirt*QOD(L)) * DOFW

                BUY(L) = (TEM1 - 1.0_r8 - QRP(L)) * ROR(L) * TX5
                !                                              Compute W (L+1/2)

                TEM1   = WVL(L)
                !             IF (ETD(L) .GT. 0.0_r8) THEN
                WVL(L) = VT(2) * (ETD(L-1)*WVL(L-1) - FACG                &
                     &                 * (BUY(L-1)*QRT(L-1)+BUY(L)*QRB(L-1)))
                !
                !             if (lprnt) print *,' wvl=',wvl(l),'vt2=',vt(2),' buy1='
                !    &,buy(l-1),' buy=',buy(l),' qrt1=',qrt(l-1),' qrb1=',qrb(l-1)
                !    &,' etd1=',etd(l-1),' wvl1=',wvl(l-1)
                !             ENDIF
                !
                IF (wvl(l) .LT. 0.0_r8) THEN
                   !               WVL(L) = max(wvl(l), 0.1_r8*tem1)
                   !               WVL(L) = 0.5_r8*tem1
                   !               WVL(L) = 0.1_r8*tem1
                   !               WVL(L) = 0.0_r8
                   WVL(L) = 1.0e-10_r8
                ELSE
                   WVL(L) = 0.5_r8*(WVL(L)+TEM1)
                ENDIF

                !
                !             WVL(L) = max(0.5_r8*(WVL(L)+TEM1), 0.0_r8)

                ERRW   = WVL(L) - TEM1
                !
                ERRQ   = ERRQ + ABS(ERRW/MAX(WVL(L),ONE_M5))

                !     if (lprnt) print *,' errw=',errw,' wvl=',wvl(l)
                !     if(lprnt .or. tx5 .eq. 0.0_r8) then
                !     if(tx5 .eq. 0.0_r8 .and. l .gt. kbl) then
                !        print *,' errq=',errq,' itr=',itr,' l=',l,' wvl=',wvl(l)
                !    &,' tx5=',tx5,' idnm=',idnm,' etd1=',etd(l-1),' etd=',etd(l)
                !    &,' kbl=',kbl
                !     endif
                !
                !     if(lprnt) print *,' itr=',itr,' itrmnd=',itrmnd,' itrmd=',itrmd
                !             IF (ITR .GE. MIN(ITRMIN,ITRMD/2)) THEN
                IF (ITR .GE. MIN(ITRMND,ITRMD/2)) THEN
                   !     if(lprnt) print *,' itr=',itr,' etd1=',etd(l-1),' errq=',errq
                   IF (ETD(L-1) .EQ. 0.0_r8 .AND. ERRQ .GT. 0.2_r8) THEN
                      !     if(lprnt) print *,' bud=',bud(kd),' wa=',wa(1),wa(2)
                      ROR(L)   = BUD(KD)
                      ETD(L)   = 0.0_r8
                      WVL(L)   = 0.0_r8
                      ERRQ     = 0.0_r8
                      HOD(L)   = WA(1)
                      QOD(L)   = WA(2)
                      !                 TX5      = TX1 + TX9
                      IF (L .LE. KBL) THEN
                         TX5      = TX9
                      ELSE
                         TX5 = (STLT(KB1) * QRT(KB1)                         &
                              &                  +  STLT(KBL) * QRB(KB1)) * (0.5_r8*FAC)
                      ENDIF

                      !     if(lprnt) print *,' tx1=',tx1,' rnt=',rnt,' rnf=',rnf(l-1)
                      !    *,' evp=',evp(l-1),' l=',l
                      EVP(L-1) = 0.0_r8
                      TEM      = MAX(TX1*RNT+RNF(L-1),ZERO)
                      QA(1)    = TEM - EVP(L-1)
                      !                 IF (QA(1) .GT. 0.0_r8) THEN
                      !     if(lprnt) print *,' ror=',ror(l),' tx5=',tx5,' tx1=',tx1
                      !    *,' tx9=',tx9,' gms=',gms(l),' qa=',qa(1)
                      !     if(lprnt) call mpi_quit(13)
                      !     if (tx5 .eq. 0.0_r8 .or. gms(l) .eq. 0.0_r8)
                      !     if (lprnt) 
                      !    *  print *,' Atx5=',tx5,' gms=',gms(l),' ror=',ror(l)
                      !    *,' L=',L,' QA=',QA(1),' tx1=',tx1,' tx9=',tx9
                      !    *,' kbl=',kbl,' etd1=',etd(l-1),' idnm=',idnm,' idn=',idn(idnm)
                      !    *,' errq=',errq
                      QRP(L)   = (QA(1) / (ROR(L)*TX5*GMS(L)))              &
                           &                                            ** (1.0_r8/1.1364_r8)
                      !                 endif
                      BUY(L)   = - ROR(L) * TX5 * QRP(L)
                      WCB(L-1) = 0.0_r8
                   ENDIF
                   !
                   DEL_ETA = ETD(L) - ETD(L-1)
                   IF(DEL_ETA .LT. 0.0_r8 .AND. ERRQ .GT. 0.1_r8) THEN
                      ROR(L)   = BUD(KD)
                      ETD(L)   = 0.0_r8
                      WVL(L)   = 0.0_r8
!!!!!             TX5      = TX1 + TX9
                      CLDFRD(L-1) = TX5
                      !
                      DEL_ETA  = - ETD(L-1)
                      EDZ      = 0.0_r8
                      DDZ      = -DEL_ETA
                      WCB(L-1) = DDZ

                      !
                      HOD(L)   = HOD(L-1)
                      QOD(L)   = QOD(L-1)

                      !
                      TEM4     = QOD(L)
                      TEM1     = VRW(1)
                      !
                      QHS      = QA(3) + 0.5_r8 * (GAF(L-1)+GAF(L))            &
                           &                                   * (HOD(L)-QA(2))

                      !
                      !                                           First iteration       !
                      !
                      ST2  = PRL(L) * (QHS + TEM1 * (QHS-QOD(L)))
                      TEM2 = ROR(L) * QRP(L-1)
                      CALL QRABF(TEM2,QRAF,QRBF)
                      TEM6 = TX5 * (1.6_r8 + 124.9_r8 * QRAF) * QRBF * TX4
                      !
                      CE   = TEM6*ST2/((5.4E5_r8*ST2 + 2.55E6_r8)*(ETD(L)+DDZ))
                      !

                      TEM2   = - ((1.0_r8+TEM1)*(QHS+CE) + TEM1*QOD(L))
                      TEM3   = (1.0_r8 + TEM1) * QHS * (QOD(L)+CE)
                      TEM    = MAX(TEM2*TEM2 -FOUR*TEM1*TEM3,ZERO)
                      QOD(L) = MAX(TEM4, (- TEM2 - SQRT(TEM)) * VRW(2))
                      !
                      !                                            second iteration   !
                      !
                      ST2  = PRL(L) * (QHS + TEM1 * (QHS-QOD(L)))
                      CE   = TEM6*ST2/((5.4E5_r8*ST2 + 2.55E6_r8)*(ETD(L)+DDZ))
                      !                 CEE  = CE * (ETD(L)+DDZ)
                      !


                      TEM2   = - ((1.0_r8+TEM1)*(QHS+CE) + TEM1*tem4)
                      TEM3   = (1.0_r8 + TEM1) * QHS * (tem4+CE)
                      TEM    = MAX(TEM2*TEM2 -FOUR*TEM1*TEM3,ZERO)
                      QOD(L) = MAX(TEM4, (- TEM2 - SQRT(TEM)) * VRW(2))

                      !                                              Evaporation in Layer L-1
                      !
                      EVP(L-1) = (QOD(L)-TEM4) * (ETD(L)+DDZ)

                      !                                               Calculate Pd (L+1/2)
                      !                 RNN(L-1) = TX1*RNT + RNF(L-1) - EVP(L-1)
                      QA(1)    = TX1*RNT + RNF(L-1)
                      EVP(L-1) = MIN(EVP(L-1), QA(1))
                      QA(1)    = QA(1) - EVP(L-1)
                      qrp(l)   = 0.0_r8
                      !
                      !     if (tx5 .eq. 0.0_r8 .or. gms(l) .eq. 0.0_r8)
                      !     if (lprnt)
                      !    *  print *,' Btx5=',tx5,' gms=',gms(l),' ror=',ror(l)
                      !    *,' L=',L,' QA=',QA(1),' tx1=',tx1,' tx9=',tx9
                      !    *,' kbl=',kbl,' etd1=',etd(l-1),' DEL_ETA=',DEL_ETA
                      !    &,' evp=',evp(l-1)
                      !
                      !                 IF (QA(1) .GT. 0.0_r8) THEN
                      !!                  RNS(L-1) = QA(1)
!!!                 tx5      = tx9
                      !                   QRP(L) = (QA(1) / (ROR(L)*TX5*GMS(L)))              &
                      !    &                                         ** (1.0_r8/1.1364_r8)
                      !                 endif
                      !                 ERRQ   = 0.0_r8
                      !                                              Compute Buoyancy
                      !                 TEM1   = WA(3)+(HOD(L)-WA(1)-con_hvap*(QOD(L)-WA(2)))     &
                      !    &                                                  * (1.0_r8/con_cp)
                      !                 TEM1   = TEM1 * (1.0_r8 + con_FVirt*QOD(L)) * DOFW

                      !                 BUY(L) = (TEM1 - 1.0_r8 - QRP(L)) * ROR(L) * TX5
                      !
                      !                 IF (QA(1) .GT. 0.0_r8) RNS(L) = QA(1)
                      IF (L .LE. K) THEN
                         RNS(L) = QA(1)
                         QA(1)  = 0.0_r8
                      ENDIF
                      tx5      = tx9
                      ERRQ     = 0.0_r8
                      QRP(L)   = 0.0_r8
                      BUY(L)   = 0.0_r8

                      !
                   ENDIF
                ENDIF
             ENDIF
             !

          ENDDO                ! End of the iteration loop  for a given L!
          !       if (kd .eq. 13 .and. .not. ddft) stop
          IF (L .LE. K) THEN
             IF (ETD(L-1) .EQ. 0.0_r8                                       &
                  &         .AND. ERRQ .GT. 0.1_r8 .AND. l .LE. kbl) THEN
!!!  &         .AND. ERRQ .GT. ERRMIN*10.0 .and. l .le. kbl) THEN
                !    &         .AND. ERRQ .GT. ERRMIN*10.0) THEN
                ROR(L)   = BUD(KD)
                HOD(L)   = WA(1)
                QOD(L)   = WA(2)
                TX5      =       TX9     ! Does not make too much difference!
                !              TX5      = TX1 + TX9
                EVP(L-1) = 0.0_r8
                !              EVP(L-1) = CEE * (1.0_r8 - qod(l)/qa(3))
                QA(1)    = TX1*RNT + RNF(L-1)
                EVP(L-1) = MIN(EVP(L-1), QA(1))
                QA(1)    = QA(1) - EVP(L-1)
                !              QRP(L)   = 0.0_r8
                IF (tx5 .EQ. 0.0_r8 .OR. gms(l) .EQ. 0.0_r8) THEN
                   PRINT *,' Ctx5=',tx5,' gms=',gms(l),' ror=',ror(l)              &
                        &,' L=',L,' QA=',QA(1),' tx1=',tx1,' tx9=',tx9                     &
                        &,' kbl=',kbl,' etd1=',etd(l-1),' DEL_ETA=',DEL_ETA
                ENDIF
                !              IF (QA(1) .GT. 0.0_r8) THEN
                QRP(L) = (QA(1) / (ROR(L)*TX5*GMS(L)))                 &
                     &                                         ** (1.0_r8/1.1364_r8)
                !              ENDIF
                ETD(L)   = 0.0_r8
                WVL(L)   = 0.0_r8
                ST1      = 1.0_r8 - ALFIND(L)

                ERRQ     = 0.0_r8
                BUY(L)   = - ROR(L) * TX5 * QRP(L)
                WCB(L-1) = 0.0_r8
             ENDIF
          ENDIF

          !
          !
          LL = MIN(IDN(idnm), K+1)
          IF (ERRQ .LT. 1.0_r8 .AND. L .LE. LL) THEN
             IF (ETD(L-1) .GT. 0.0_r8 .AND. ETD(L) .EQ. 0.0_r8) THEN
                IDN(idnm) = L
                wvl(l)    = 0.0_r8
                IF (L .LT. KBL .OR. tx5 .GT. 0.0_r8) idnm  = idnm + 1
                errq      = 0.0_r8
             ENDIF
             IF (etd(l) .EQ. 0.0_r8 .AND. l .GT. kbl) THEN
                idn(idnm) = l
                IF (tx5 .GT. 0.0_r8) idnm  = idnm + 1
             ENDIF
          ENDIF

          !       if (lprnt) then
          !       print *,' ERRQ=',ERRQ,' IDN=',IDN(idnm),' idnm=',idnm
          !       print *,' L=',L,' QRP=',QRP(L),' ETD=',ETD(L),' QA=',QA(1)
          !    *,' evp=',evp(l-1),' rnf=',rnf(l-1)
          !       endif

          ! 
          !     If downdraft properties are not obtainable, (i.e.solution does
          !      not converge) , no downdraft is assumed
          !
          !          IF (ERRQ .GT. ERRMIN*100.0_r8 .AND. IDN(idnm) .EQ. 99)          &
          IF (ERRQ .GT. 0.1_r8 .AND. IDN(idnm) .EQ. 99)                   &
               &                          DDFT = .FALSE.
          !
          !
          DOF = 0.0_r8
          IF (.NOT. DDFT) RETURN
          !
          !     if (ddlgk .or. l .le. idn(idnm)) then
          !     rsum2 = rsum2 + evp(l-1)
          !     print *,' rsum1=',rsum1,' rsum2=',rsum2,' L=',L,' qa=',qa(1)
          !    *,' evp=',evp(l-1)
          !     else
          !     rsum1 = rsum1 + rnf(l-1)
          !     print *,' rsum1=',rsum1,' rsum2=',rsum2,' L=',L,' rnf=',rnf(l-1)
          !     endif

       ENDDO                      ! End of the L Loop of downdraft !

       TX1 = 0.0_r8

       DOF = QA(1)
       !
       !     print *,' dof=',dof,' rntp=',rntp,' rnb=',rnb
       !     print *,' total=',(rsum1+dof+rntp+rnb)

    ENDIF                       ! SKPDD endif
    !

    RNN(KD) = RNTP
    TX1     = EVP(KD)
    TX2     = RNTP + RNB + DOF

    !     if (lprnt) print *,' tx2=',tx2
    II = IDH
    IF (II .GE. KD1+1) THEN
       RNN(KD)   = RNN(KD) + RNF(KD)
       TX2       = TX2 + RNF(KD)
       RNN(II-1) = 0.0_r8
       TX1       = EVP(II-1)
    ENDIF
    !     if (lprnt) print *,' tx2=',tx2,' idnm=',idnm,' idn=',idn(idnm)
    DO L=KD,K
       II = IDH

       IF (L .GT. KD1 .AND. L .LT. II) THEN
          RNN(L-1) = RNF(L-1)
          TX2      = TX2 + RNN(L-1)

       ELSEIF (L .GE. II .AND. L .LT. IDN(idnm)) THEN
!!!       ELSEIF (L .GE. II .AND. L .LE. IDN(idnm)) THEN

          !           do jj=2,idnm
          !             if (l .ge. idn(jj-1) .and. l .lt. idn(jj)) then
!!!             RNN(L)   = 0.0
!!!             TX1      = TX1 + EVP(L)
          !             endif
          !           enddo
          !
          rnn(l) = rns(l)
          tx2    = tx2 + rnn(l)
          TX1    = TX1 + EVP(L)

       ELSEIF (L .GE. IDN(idnm)) THEN
          ETD(L+1) = 0.0_r8
          HOD(L+1) = 0.0_r8
          QOD(L+1) = 0.0_r8
          EVP(L)   = 0.0_r8
          RNN(L)   = RNF(L) + RNS(L)
          TX2      = TX2    + RNN(L)
       ENDIF
       !     if (lprnt) print *,' tx2=',tx2,' L=',L,' rnn=',rnn(l)
    ENDDO
    !       IF (K+1 .GT. IDN(idnm)) THEN
    !         ETD(K+1) = 0.0_r8
    !         HOD(K+1) = 0.0_r8
    !         QOD(K+1) = 0.0_r8
    !         EVP(K)   = 0.0_r8
    !         RNN(K)   = RNF(K)
    !         TX2      = TX2 + RNN(K)
    !       ENDIF
    !
    !      For Downdraft case the rain is that falls thru the bottom

    L = KBL

    RNN(L)    = RNN(L) + RNB
    CLDFRD(L) = TX5

    !
    !     Caution !! Below is an adjustment to rain flux to maintain
    !                conservation of precip!

    !
    !     if (lprnt) print *,' train=',train,' tx2=',tx2,' tx1=',tx1

    IF (TX1 .GT. 0.0_r8) THEN
       TX1 = (TRAIN - TX2) / TX1
    ELSE
       TX1 = 0.0_r8
    ENDIF

    !       TX5      = EVP(KBL)

    !!      EVP(KBL) = EVP(KBL) * TX1

    !       TX3      = RNN(KBL) + EVP(KBL) + DOF
    !       TX2      = RNN(KBL)
    !       TX4      = EVP(KBL)

    !       DO L=KD,KB1
    DO L=KD,K

       !         TX5    = TX5 + EVP(L)
       EVP(L) = EVP(L) * TX1
       !         TX3    = TX3 + EVP(L) + RNN(L)
       !         TX2    = TX2 + RNN(L)
       !         TX4    = TX4 + EVP(L)
    ENDDO
    !
    !     if (lprnt .and. kd .eq. 52) stop
    !***********************************************************************
    !***********************************************************************

    RETURN
  END   SUBROUTINE DDRFT

  !-----------------------------------------------------------------------------------------

  SUBROUTINE QSATCN(TT,P,Q,DQDT)

    !      USE FUNCPHYS , ONLY : fpvs
    IMPLICIT NONE
    !     include 'constant.h'
    !
    REAL(kind=r8) TT, P, Q, DQDT
    !
    REAL(kind=r8) rvi, facw, faci, hsub, tmix, DEN
    REAL(kind=r8) ZERO,ONE,ONE_M10
    PARAMETER (RVI=1.0_r8/con_rv)
    PARAMETER (FACW=con_cvap-con_cliq, FACI=con_cvap-con_CSOL)
    PARAMETER (HSUB=con_hvap+con_hfus, tmix=con_ttp-20.0_r8, DEN=1.0_r8/(con_ttp-TMIX))
    PARAMETER (ZERO=0.0_r8,ONE=1.0_r8,ONE_M10=1.E-10_r8)
    !
    !CFPP$ NOCONCUR R
    REAL(kind=r8) es, d, hlorv, W
    !
    !     es    = 10.0_r8 * fpvs(tt)                ! fpvs is in centibars!
    es    = 0.01_r8 * fpvs(tt)                ! fpvs is in Pascals!
    D     = 1.0_r8 / MAX(p+con_epsm1*es,ONE_M10)
    !
    q     = MIN(con_eps*es*D, ONE)
    !
    W     = MAX(ZERO, MIN(ONE, (TT - TMIX)*DEN))
    hlorv = ( W      * (con_hvap + FACW * (tt-con_ttp))                       &
         &       + (1.0_r8-W) * (HSUB + FACI * (tt-con_ttp)) ) * RVI
    dqdt  = p * q * hlorv *  D / (tt*tt)
    !
    RETURN
  END SUBROUTINE QSATCN
  !-----------------------------------------------------------------------------------------

  SUBROUTINE ANGRAD( PRES, ALM,  AL2, TLA)
    !     SUBROUTINE ANGRAD( PRES, ALM, STLA, CTL2, AL2                     &
    !    &,                  PI, TLA, PRB, WFN, UFN)
    !      use module_ras , only : refp, refr, tlac, plac, tlbpl, drdp, almax
    IMPLICIT NONE

    !     real(kind=r8) PRES, STLA, CTL2, pi,  pifac                 &
    REAL(kind=r8) ::PRES                                         &
         &,                    ALM,  AL2,  TLA,  TEM 
    !
    !
    !     pifac = pi / 180.0_r8
    !     print *,' pres=',pres
    IF (TLA .LT. 0.0_r8) THEN
       IF (PRES .LE. PLAC(1)) THEN
          TLA = TLAC(1)
       ELSEIF (PRES .LE. PLAC(2)) THEN
          TLA = TLAC(2) + (PRES-PLAC(2))*tlbpl(1)
       ELSEIF (PRES .LE. PLAC(3)) THEN
          TLA = TLAC(3) + (PRES-PLAC(3))*tlbpl(2)
       ELSEIF (PRES .LE. PLAC(4)) THEN
          TLA = TLAC(4) + (PRES-PLAC(4))*tlbpl(3)
       ELSEIF (PRES .LE. PLAC(5)) THEN
          TLA = TLAC(5) + (PRES-PLAC(5))*tlbpl(4)
       ELSEIF (PRES .LE. PLAC(6)) THEN
          TLA = TLAC(6) + (PRES-PLAC(6))*tlbpl(5)
       ELSEIF (PRES .LE. PLAC(7)) THEN
          TLA = TLAC(7) + (PRES-PLAC(7))*tlbpl(6)
       ELSEIF (PRES .LE. PLAC(8)) THEN
          TLA = TLAC(8) + (PRES-PLAC(8))*tlbpl(7)
       ELSE
          TLA = TLAC(8)
       ENDIF
       !         tla = tla * 1.5

       !         STLA = SIN(TLA*PIFAC)
       !         TEM1 = COS(TLA*PIFAC)
       !         CTL2 = TEM1 * TEM1

    ELSE
       !         STLA = SIN(TLA*PIFAC)
       !         TEM1 = COS(TLA*PIFAC)
       !         CTL2 = TEM1 * TEM1

    ENDIF
    IF (PRES .GE. REFP(1)) THEN
       TEM = REFR(1)
    ELSEIF (PRES .GE. REFP(2)) THEN
       TEM = REFR(1) + (PRES-REFP(1)) * drdp(1)
    ELSEIF (PRES .GE. REFP(3)) THEN
       TEM = REFR(2) + (PRES-REFP(2)) * drdp(2)
    ELSEIF (PRES .GE. REFP(4)) THEN
       TEM = REFR(3) + (PRES-REFP(3)) * drdp(3)
    ELSEIF (PRES .GE. REFP(5)) THEN
       TEM = REFR(4) + (PRES-REFP(4)) * drdp(4)
    ELSEIF (PRES .GE. REFP(6)) THEN
       TEM = REFR(5) + (PRES-REFP(5)) * drdp(5)
    ELSE
       TEM = REFR(6)
    ENDIF
    !!      AL2 = min(ALMAX, MAX(ALM, 2.0E-4/TEM))
    !       AL2 = min(2.0E-3, MAX(ALM, 2.0E-4/TEM))
    !
    tem = 2.0E-4_r8 / tem
    al2 = MIN(4.0_r8*tem, MAX(alm, tem))
    !
    RETURN
  END SUBROUTINE ANGRAD
  !-----------------------------------------------------------------------------------------

  SUBROUTINE SETQRP()
    ! use module_ras , only : NQRP,C1XQRP,C2XQRP,TBQRP,TBQRA,TBQRB
    IMPLICIT NONE

    REAL(kind=r8) :: tem2,tem1,x,xinc,xmax,xmin
    INTEGER :: jx
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !CFPP$ NOCONCUR R
    !     XMIN=1.0E-6
    XMIN=0.0_r8
    XMAX=5.0_r8
    XINC=(XMAX-XMIN)/(NQRP-1)
    C1XQRP=1.0_r8-XMIN/XINC
    C2XQRP=1.0_r8/XINC
    TEM1 = 0.001_r8 ** 0.2046_r8
    TEM2 = 0.001_r8 ** 0.525_r8
    DO JX=1,NQRP
       X         = XMIN + (JX-1)*XINC
       TBQRP(JX) =        X ** 0.1364_r8
       TBQRA(JX) = TEM1 * X ** 0.2046_r8
       TBQRB(JX) = TEM2 * X ** 0.525_r8
    ENDDO
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    RETURN
  END  SUBROUTINE SETQRP

  !-----------------------------------------------------------------------------------------

  REAL(kind=r8) FUNCTION QRPF(QRP)
    !
    !      use module_ras , only : NQRP,C1XQRP,C2XQRP,TBQRP,TBQRA,TBQRB
    IMPLICIT NONE

    REAL(kind=r8) QRP, XJ, REAL_NQRP
    REAL(kind=r8), PARAMETER :: ONE=1.0_r8
    INTEGER JX
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    REAL_NQRP=REAL(NQRP)
    XJ   = MIN(MAX(C1XQRP+C2XQRP*QRP,ONE),REAL_NQRP)
    !     XJ   = MIN(MAX(C1XQRP+C2XQRP*QRP,ONE),FLOAT(NQRP))
    JX   = INT(MIN(XJ,NQRP-ONE),KIND=i4)
    QRPF = TBQRP(JX)  + (XJ-JX) * (TBQRP(JX+1)-TBQRP(JX))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    RETURN
  END FUNCTION QRPF
  !-----------------------------------------------------------------------------------------

  SUBROUTINE QRABF(QRP,QRAF,QRBF)
    !      use module_ras , only : NQRP,C1XQRP,C2XQRP,TBQRP,TBQRA,TBQRB
    IMPLICIT NONE
    !
    REAL(kind=r8) QRP, QRAF, QRBF, XJ, REAL_NQRP
    REAL(kind=r8), PARAMETER :: ONE=1.0_r8
    INTEGER JX
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    REAL_NQRP=REAL(NQRP)
    XJ   = MIN(MAX(C1XQRP+C2XQRP*QRP,ONE),REAL_NQRP)
    JX   = INT(MIN(XJ,NQRP-ONE),KIND=i4)
    XJ   = XJ - JX
    QRAF = TBQRA(JX)  + XJ * (TBQRA(JX+1)-TBQRA(JX))
    QRBF = TBQRB(JX)  + XJ * (TBQRB(JX+1)-TBQRB(JX))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    RETURN
  END SUBROUTINE QRABF

  !-----------------------------------------------------------------------------------------

  SUBROUTINE SETVTP()
    !use module_ras , only : NVTP,C1XVTP,C2XVTP,TBVTP
    IMPLICIT NONE

    REAL(kind=r8) :: xinc,x,xmax,xmin
    INTEGER :: jx
    REAL(kind=r8), PARAMETER :: VTPEXP=-0.3636_r8
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !CFPP$ NOCONCUR R
    XMIN=0.05_r8
    XMAX=1.5_r8
    XINC=(XMAX-XMIN)/(NVTP-1)
    C1XVTP=1.0_r8-XMIN/XINC
    C2XVTP=1.0_r8/XINC
    DO JX=1,NVTP
       X         = XMIN + (JX-1)*XINC
       TBVTP(JX) =        X ** VTPEXP
    ENDDO
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    RETURN
  END  SUBROUTINE SETVTP


  !-----------------------------------------------------------------------------------------

 REAL(kind=r8)  FUNCTION VTPF(ROR)
    !
    !use module_ras , only : NVTP,C1XVTP,C2XVTP,TBVTP
    IMPLICIT NONE
    REAL(kind=r8) :: ROR, XJ, REAL_NVTP
    REAL(kind=r8), PARAMETER :: ONE=1.0_r8
    INTEGER :: JX
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    REAL_NVTP=REAL(NVTP)
    XJ   = MIN(MAX(C1XVTP+C2XVTP*ROR,ONE),REAL_NVTP)
    JX   = INT(MIN(XJ,NVTP-ONE),KIND=i4)
    VTPF = TBVTP(JX)  + (XJ-JX) * (TBVTP(JX+1)-TBVTP(JX))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    RETURN
  END  FUNCTION VTPF
  !-----------------------------------------------------------------------------------------


    REAL(kind=r8) FUNCTION CLF(PRATE)
    !
    IMPLICIT NONE
    REAL(kind=r8) PRATE
    !
    REAL (kind=r8), PARAMETER :: ccf1=0.30_r8, ccf2=0.09_r8                 &
         &,                                   ccf3=0.04_r8, ccf4=0.01_r8          &
         &,                                   pr1=1.0_r8,   pr2=5.0_r8            &
         &,                                   pr3=20.0_r8
    !
    IF (prate .LT. pr1) THEN
       clf = ccf1
    ELSEIF (prate .LT. pr2) THEN
       clf = ccf2
    ELSEIF (prate .LT. pr3) THEN
       clf = ccf3
    ELSE
       clf = ccf4
    ENDIF
    !
    RETURN
  END FUNCTION CLF



  !-----------------------------------------------------------------------------------------
  SUBROUTINE Finalize_Cu_RAS3PHASE()
    IMPLICIT NONE
    DEALLOCATE (rasal)

  END SUBROUTINE Finalize_Cu_RAS3PHASE
  !-----------------------------------------------------------------------------------------

END MODULE Cu_RAS3PHASE


!PROGRAM Main
! USE Cu_RAS3PHASE, ONLY: Init_Cu_RAS3PHASE
! IMPLICIT NONE
! INTEGER :: kMax=28
! REAL(KIND=8) :: dt=1200.0
! REAL(KIND=8) :: fhour=0.0
! INTEGER :: idate(1:4)=(/00,01,29,2015/)
! INTEGER :: iMax=2
! INTEGER :: jMax=1
!! INTEGER :: ibMax=2
! INTEGER :: jbMax=1
! REAL(kind=8) ::  si_in(kMax+1)
! REAL(kind=8) ::  sl_in(kMax)
!
! CALL Init_Cu_RAS3PHASE(kMax,dt,fhour,idate,iMax,jMax,ibMax,jbMax,si_in,sl_in)
! PRINT*,' --  '
! CALL Init_Cu_RAS3PHASE(kMax,dt,fhour,idate,iMax,jMax,ibMax,jbMax,si_in,sl_in)
!
!END PROGRAM Main
