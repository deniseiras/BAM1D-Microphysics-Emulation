!
!  $Author: pkubota $
!  $Date: 2011/04/07 16:00:31 $
!  $Revision: 1.18 $
!
MODULE Utils

  ! CreateAssocLegFunc
  ! DestroyAssocLegFunc
  !
  ! DumpAssocLegFunc   ------------------| DumpMatrix(interface)
  !
  ! CreateGaussQuad    ------------------| CreateLegPol
  !
  ! CreateGridValues
  !
  ! DestroyGaussQuad   ------------------| DestroyLegPol
  !
  ! DumpGaussQuad      ------------------| DumpMatrix(interface)
  !
  ! iminv
  !
  ! Rg                 ------------------| Balanc
  !                                  |
  !                                  | Orthes
  !                                  |
  !                                  | Ortran
  !                                  |
  !                                  | Hqr2    ------|Hqr3
  !                                  |
  !                                  | Balbak
  !                                  |
  !                                  | Znorma
  !
  ! Tql2
  ! Tred2
  ! tmstmp2
  ! InitTimeStamp
  !
  ! TimeStamp  --------------------------| caldat
  !
  ! IBJBtoIJ_R (Interface)
  ! IJtoIBJB_R (Interface)
  ! IBJBtoIJ_I (Interface)
  ! IJtoIBJB_I (Interface)
  !
  ! SplineIJtoIBJB_R2D (Interface) ------| CyclicCubicSpline
  !
  ! SplineIBJBtoIJ_R2D (Interface) ------| CyclicCubicSpline
  !
  ! LinearIJtoIBJB_R2D (Interface) ------| CyclicLinear
  !
  ! LinearIBJBtoIJ_R2D (Interface) ------| CyclicLinear
  !
  ! NearestIJtoIBJB_I2D (Interface)------| CyclicNearest_i
  !
  ! NearestIJtoIBJB_R2D (Interface)------| CyclicNearest_r
  !
  ! NearestIBJBtoIJ_I2D (Interface)------| CyclicNearest_i
  !
  ! NearestIBJBtoIJ_R2D (Interface)------| CyclicNearest_r
  !
  ! FreqBoxIJtoIBJB_I2D (Interface)------| CyclicFreqBox_i
  !
  ! FreqBoxIJtoIBJB_R2D (Interface)------| CyclicFreqBox_r
  !
  ! SeaMaskIJtoIBJB_R2D (Interface)------| CyclicSeaMask_r
  !
  ! SeaMaskIBJBtoIJ_R2D (Interface)------| CyclicSeaMask_r
  !
  ! AveBoxIJtoIBJB_R2D (Interface) ------| CyclicAveBox_r
  !
  ! AveBoxIBJBtoIJ_R2D (Interface) ------| CyclicAveBox_r
  !
  ! vfirec     --------------------------| vfinit
  !
  !
  !-------------------------------------------------------------------
  !  ASSOCIATED LEGENDRE FUNCTIONS
  !  Module computes and stores Associated Legendre
  !  Functions and Epslon.
  !
  !  Module exports three routines:
  !     CreateAssocLegFunc  initializes module and compute functions;
  !     DestroyAssocLegFunc destroys module;
  !     DumpAssocLegFunc    Dumps module info
  !
  !  Module usage:
  !  CreateAssocLegFunc should be invoked once, prior to any other
  !  module routine. It computes and hides the function values.
  !  DestroyAssocLegFunc destroys all internal info.
  !
  !  Module use values exported by Sizes and procedures from Auxiliary


  USE Sizes, ONLY:  &
       iMax,        &
       jMax,        &
       kMax,        &
       jbMax,       &
       ibMax


  USE Constants, ONLY : &
       i4,         &
       i8,         &
       r8,         &
       pai,        &
       twomg,      &
       er,         &
       r16


  USE Parallelism, ONLY: &
       MsgOne, &
       FatalError

 IMPLICIT NONE

  PRIVATE

  !
  !  LEGANDRE POLINOMIAL AND ITS ROOTS
  !
  !  Module exports four routines:
  !     CreateLegPol  initializes module;
  !     DestroyLegPol destroys module;
  !     LegPol        computes polinomial
  !     LegPolRoots   computes roots of even degree Legandre Pol
  !
  !  Module does not export (or require) any data value.
  !
  PUBLIC :: InitTimeStamp
  PUBLIC :: TimeStamp
  PUBLIC :: IJtoIBJB
  PUBLIC :: NearestIJtoIBJB
  PUBLIC :: LinearIJtoIBJB
  PUBLIC :: FreqBoxIJtoIBJB
  PUBLIC :: AveBoxIJtoIBJB

! Enver start to include Coriolis and pressure gradient effects in the momentum equation
!  PUBLIC :: fcor
!  PUBLIC :: colrad
! Enver end 

  INTEGER,  PRIVATE :: JulianDayInitIntegration  
  
  INTERFACE IJtoIBJB
     MODULE PROCEDURE &
          IJtoIBJB_R, IJtoIBJB_I, &
          IJtoIBJB3_R, IJtoIBJB3_I
  END INTERFACE
  INTERFACE NearestIJtoIBJB
     MODULE PROCEDURE &
          NearestIJtoIBJB_I2D, NearestIJtoIBJB_R2D, &
          NearestIJtoIBJB_I3D, NearestIJtoIBJB_R3D
  END INTERFACE
  INTERFACE LinearIJtoIBJB
     MODULE PROCEDURE LinearIJtoIBJB_R2D
  END INTERFACE
  INTERFACE FreqBoxIJtoIBJB
     MODULE PROCEDURE FreqBoxIJtoIBJB_I2D, FreqBoxIJtoIBJB_R2D
  END INTERFACE
  INTERFACE AveBoxIJtoIBJB
     MODULE PROCEDURE AveBoxIJtoIBJB_R2D
  END INTERFACE

  PUBLIC :: lati
  PUBLIC :: CreateGridValues
  PUBLIC :: vfirec
  PUBLIC :: tmstmp2
  REAL(KIND=r8), ALLOCATABLE :: lati(:)
!Enver start
  REAL(KIND=r8), PUBLIC  :: fcor         !single real in 1D for Coriolis parameter   , ALLOCATABLE :: fcor
!  REAL(KIND=r8)  :: colrad(jMax) !single real in 1D   ,  ALLOCATABLE :: colrad
!Enver end


CONTAINS
  SUBROUTINE CreateGridValues(colrad)
    REAL(KIND=r8) :: colrad(jMax)
    !REAL(KIND=r8) :: fcor
    INTEGER :: j
    ALLOCATE (lati    (jMax))

    DO j=1,jMax
       lati(j)=colrad(j)
    END DO
    
     fcor = twomg * COS(colrad(1))

  END SUBROUTINE CreateGridValues

  SUBROUTINE FreqBoxIJtoIBJB_I2D(FieldIn,FieldOut)
    INTEGER(KIND=i8), INTENT(IN  ) :: FieldIn (iMax,jMax)
    INTEGER(KIND=i8), INTENT(OUT ) :: FieldOut(ibMax,jbMax)

    INTEGER(KIND=i8)   :: FOut(iMax)
    INTEGER            :: i
    INTEGER            :: iFirst
    INTEGER            :: ilast
    INTEGER            :: ib
    INTEGER            :: j
    INTEGER            :: jb

    CHARACTER(LEN=*), PARAMETER :: h="**FreqBoxIJtoIBJB**"

    !IF (reducedGrid) THEN
    !   !$OMP PARALLEL DO PRIVATE(iFirst,ilast,fout,i,ib,jb)
    !   DO j = myfirstlat,mylastlat
    !      iFirst = myfirstlon(j)
    !      ilast  = mylastlon(j)
    !      CALL CyclicFreqBox_i(iMax, iMaxPerJ(j), &
    !      FieldIn(1:iMax,j), FOut(1:iMax), ifirst, ilast)
    !      DO i = ifirst,ilast
    !         ib = ibperij(i,j)
    !         jb = jbperij(i,j)
    !         FieldOut(ib,jb) = Fout(i)
    !      END DO
    !   END DO
    !   !$OMP END PARALLEL DO
    !ELSE
    !   !$OMP PARALLEL DO PRIVATE(ib,i,j)
       DO jb = 1, jbMax
          DO ib = 1, ibMax!PerJB(jb)
             !i = iPerIJB(ib,jb)
             !j = jPerIJB(ib,jb)
             FieldOut(ib,jb)=FieldIn(ib,jb)
          END DO
       END DO
    !   !$OMP END PARALLEL DO
    !END IF
  END SUBROUTINE FreqBoxIJtoIBJB_I2D

  SUBROUTINE FreqBoxIJtoIBJB_R2D(FieldIn,FieldOut)
    REAL(KIND=r8)   , INTENT(IN  ) :: FieldIn (iMax,jMax)
    REAL(KIND=r8)   , INTENT(OUT ) :: FieldOut(ibMax,jbMax)

    REAL(KIND=r8)      :: FOut(iMax)
    INTEGER            :: i
    INTEGER            :: iFirst
    INTEGER            :: ilast
    INTEGER            :: ib
    INTEGER            :: j
    INTEGER            :: jb

    CHARACTER(LEN=*), PARAMETER :: h="**FreqBoxIJtoIBJB**"

    !IF (reducedGrid) THEN
    !   !$OMP PARALLEL DO PRIVATE(iFirst,ilast,fout,i,ib,jb)
    !   DO j = myfirstlat,mylastlat
    !      iFirst = myfirstlon(j)
    !      ilast  = mylastlon(j)
    !      CALL CyclicFreqBox_r(iMax, iMaxPerJ(j), &
    !      FieldIn(1:,j), FOut, ifirst, ilast)
    !      DO i = ifirst,ilast
    !         ib = ibperij(i,j)
    !         jb = jbperij(i,j)
    !         FieldOut(ib,jb) = Fout(i)
    !      END DO
    !   END DO
    !   !$OMP END PARALLEL DO
    !ELSE
    !   !$OMP PARALLEL DO PRIVATE(ib,i,j)
       DO jb = 1, jbMax
          DO ib = 1, ibMax!PerJB(jb)
             !i = iPerIJB(ib,jb)
             !j = jPerIJB(ib,jb)
             FieldOut(ib,jb)=FieldIn(ib,jb)
          END DO
       END DO
    !   !$OMP END PARALLEL DO
    !END IF
  END SUBROUTINE FreqBoxIJtoIBJB_R2D


  SUBROUTINE LinearIJtoIBJB_R2D(FieldIn,FieldOut)
    REAL(KIND=r8), INTENT(IN  ) :: FieldIn (iMax,jMax)
    REAL(KIND=r8), INTENT(OUT ) :: FieldOut(ibMax,jbMax)

    REAL(KIND=r8) :: Fout(imax)
    INTEGER            :: i
    INTEGER            :: iFirst
    INTEGER            :: ilast
    INTEGER            :: ib
    INTEGER            :: j
    INTEGER            :: jb

    CHARACTER(LEN=*), PARAMETER :: h="**LinearIJtoIBJB**"

    !IF (reducedGrid) THEN
    !   !$OMP PARALLEL DO PRIVATE(iFirst,ilast,i,ib,fout,jb)
    !   DO j = myfirstlat,mylastlat
    !      ifirst = myfirstlon(j)
    !      ilast  = mylastlon(j)
    !      CALL CyclicLinear(iMax, iMaxPerJ(j), &
    !           FieldIn(1,j), FOut, ifirst,ilast)
    !      DO i = ifirst,ilast
    !         ib = ibperij(i,j)
    !         jb = jbperij(i,j)
    !         FieldOut(ib,jb) = Fout(i)
    !      END DO
    !   END DO
    !   !$OMP END PARALLEL DO
    !ELSE
    !   !$OMP PARALLEL DO PRIVATE(ib,i,j)
       DO jb = 1, jbMax
          DO ib = 1, ibMax!PerJB(jb)
             !i = iPerIJB(ib,jb)
             !j = jPerIJB(ib,jb)
             FieldOut(ib,jb)=FieldIn(ib,jb)
          END DO
       END DO
    !  !$OMP END PARALLEL DO
    !END IF
  END SUBROUTINE LinearIJtoIBJB_R2D



  SUBROUTINE NearestIJtoIBJB_I2D(FieldIn,FieldOut)
    INTEGER, INTENT(IN  ) :: FieldIn (iMax,jMax)
    INTEGER, INTENT(OUT ) :: FieldOut(ibMax,jbMax)

    INTEGER            :: Fout(iMax)
    INTEGER            :: i
    INTEGER            :: iFirst
    INTEGER            :: ilast
    INTEGER            :: ib
    INTEGER            :: j
    INTEGER            :: jb

    CHARACTER(LEN=*), PARAMETER :: h="**NearestIJtoIBJB**"

    !IF (reducedGrid) THEN
    !   !$OMP PARALLEL DO PRIVATE(iFirst,ilast,fout,i,ib,jb)
    !   DO j = myfirstlat,mylastlat
    !      iFirst = myfirstlon(j)
    !      ilast  = mylastlon(j)
    !      CALL CyclicNearest_i(iMax, iMaxPerJ(j), &
    !           FieldIn(1,j), FOut, ifirst, ilast)
    !      DO i = ifirst,ilast
    !         ib = ibperij(i,j)
    !         jb = jbperij(i,j)
    !         FieldOut(ib,jb) = Fout(i)
    !      END DO
    !   END DO
    !   !$OMP END PARALLEL DO
    !ELSE
    !   !$OMP PARALLEL DO PRIVATE(ib,i,j)
       DO jb = 1, jbMax
          DO ib = 1, ibMax!PerJB(jb)
             !i = iPerIJB(ib,jb)
             !j = jPerIJB(ib,jb)
             FieldOut(ib,jb)=FieldIn(ib,jb)
          END DO
       END DO
    !   !$OMP END PARALLEL DO
    !END IF
  END SUBROUTINE NearestIJtoIBJB_I2D



  SUBROUTINE NearestIJtoIBJB_R2D(FieldIn,FieldOut)
    REAL(KIND=r8), INTENT(IN  ) :: FieldIn (iMax,jMax)
    REAL(KIND=r8), INTENT(OUT ) :: FieldOut(ibMax,jbMax)

    REAL(KIND=r8)      :: Fout(iMax)
    INTEGER            :: i
    INTEGER            :: iFirst
    INTEGER            :: ilast 
    INTEGER            :: ib
    INTEGER            :: j
    INTEGER            :: jb

    CHARACTER(LEN=*), PARAMETER :: h="**NearestIJtoIBJB**"

    !IF (reducedGrid) THEN
    !   !$OMP PARALLEL DO PRIVATE(iFirst,ilast,fout,i,ib,jb)
    !   DO j = myfirstlat,mylastlat
    !      iFirst = myfirstlon(j)
    !      ilast  = mylastlon(j)
    !      CALL CyclicNearest_r(iMax, iMaxPerJ(j), &
    !           FieldIn(1,j), FOut, ifirst, ilast)
    !      DO i = ifirst,ilast
    !         ib = ibperij(i,j)
    !         jb = jbperij(i,j)
    !         FieldOut(ib,jb) = Fout(i)
    !      END DO
    !   END DO
    !   !$OMP END PARALLEL DO
    !ELSE
    !   !$OMP PARALLEL DO PRIVATE(ib,i,j)
       DO jb = 1, jbMax
          DO ib = 1, ibMax!PerJB(jb)
             !i = iPerIJB(ib,jb)
             !j = jPerIJB(ib,jb)
             FieldOut(ib,jb)=FieldIn(ib,jb)
          END DO
       END DO
    !   !$OMP END PARALLEL DO
    !END IF
  END SUBROUTINE NearestIJtoIBJB_R2D

  ! 3D version by hmjb
  SUBROUTINE NearestIJtoIBJB_I3D(FieldIn,FieldOut)
    INTEGER, INTENT(IN  ) :: FieldIn (iMax,kMax,jMax)
    INTEGER, INTENT(OUT ) :: FieldOut(ibMax,kMax,jbMax)

    INTEGER            :: FOut(iMax)
    INTEGER            :: i,k
    INTEGER            :: ifirst
    INTEGER            :: ilast 
    INTEGER            :: ib
    INTEGER            :: j
    INTEGER            :: jb

    CHARACTER(LEN=*), PARAMETER :: h="**NearestIJtoIBJB**"

    !IF (reducedGrid) THEN
    !   !$OMP PARALLEL DO PRIVATE(iFirst,ilast,fout,i,ib,jb,k)
    !   DO j = myfirstlat,mylastlat
    !      iFirst = myfirstlon(j)
    !      ilast  = mylastlon(j)
    !      DO k=1,kMax 
    !         CALL CyclicNearest_i(iMax, iMaxPerJ(j), &
    !              FieldIn(1,k,j), FOut, ifirst, ilast)
    !         DO i = ifirst,ilast
    !            ib = ibperij(i,j)
    !            jb = jbperij(i,j)
    !            FieldOut(ib,k,jb) = Fout(i)
    !         END DO
    !      ENDDO
    !   END DO
    !   !$OMP END PARALLEL DO
    !ELSE
    !   !$OMP PARALLEL DO PRIVATE(ib,i,j)
       DO jb = 1, jbMax
          DO ib = 1, ibMax!PerJB(jb)
            ! i = iPerIJB(ib,jb)
            ! j = jPerIJB(ib,jb)
             FieldOut(ib,:,jb)=FieldIn(ib,:,jb)
          END DO
       END DO
    !   !$OMP END PARALLEL DO
    !END IF
  END SUBROUTINE NearestIJtoIBJB_I3D


  ! 3D version by hmjb
  SUBROUTINE NearestIJtoIBJB_R3D(FieldIn,FieldOut)
    REAL(KIND=r8), INTENT(IN  ) :: FieldIn (iMax,kMax,jMax)
    REAL(KIND=r8), INTENT(OUT ) :: FieldOut(ibMax,kMax,jbMax)

    REAL(KIND=r8)      :: FOut(iMax)
    INTEGER            :: i,k
    INTEGER            :: iFirst
    INTEGER            :: ilast
    INTEGER            :: ib
    INTEGER            :: j
    INTEGER            :: jb

    CHARACTER(LEN=*), PARAMETER :: h="**NearestIJtoIBJB**"

    !IF (reducedGrid) THEN
    !   !$OMP PARALLEL DO PRIVATE(iFirst,ilast,fout,i,ib,jb,k)
    !   DO j = myfirstlat,mylastlat
    !      iFirst = myfirstlon(j)
    !      ilast  = mylastlon(j)
    !      DO k=1,kMax 
    !         CALL CyclicNearest_r(iMax, iMaxPerJ(j), &
    !              FieldIn(1,k,j), FOut, ifirst, ilast)
    !         DO i = ifirst,ilast
    !            ib = ibperij(i,j)
    !            jb = jbperij(i,j)
    !            FieldOut(ib,k,jb) = Fout(i)
    !         END DO
    !      ENDDO
    !   END DO
    !   !$OMP END PARALLEL DO
    !ELSE
    !   !$OMP PARALLEL DO PRIVATE(ib,i,j)
       DO jb = 1, jbMax
          DO ib = 1, ibMax!PerJB(jb)
             !i = iPerIJB(ib,jb)
             !j = jPerIJB(ib,jb)
             FieldOut(ib,:,jb)=FieldIn(ib,:,jb)
          END DO
       END DO
    !   !$OMP END PARALLEL DO
    !END IF
  END SUBROUTINE NearestIJtoIBJB_R3D

 SUBROUTINE AveBoxIJtoIBJB_R2D(FieldIn,FieldOut)
    REAL(KIND=r8)   , INTENT(IN  ) :: FieldIn (iMax,jMax)
    REAL(KIND=r8)   , INTENT(OUT ) :: FieldOut(ibMax,jbMax)

    REAL(KIND=r8)      :: FOut(iMax)
    INTEGER            :: i
    INTEGER            :: iFirst
    INTEGER            :: ilast
    INTEGER            :: ib
    INTEGER            :: j
    INTEGER            :: jb

    CHARACTER(LEN=*), PARAMETER :: h="**AveBoxIJtoIBJB_R2D**"

    !IF (reducedGrid) THEN
    !   !$OMP PARALLEL DO PRIVATE(iFirst,ilast,fout,i,ib,jb)
    !   DO j = myfirstlat,mylastlat
    !      iFirst = myfirstlon(j)
    !      ilast  = mylastlon(j)
    !      CALL CyclicAveBox_r(iMax, iMaxPerJ(j), &
    !           FieldIn(1,j), FOut, ifirst, ilast)
    !      DO i = ifirst,ilast
    !         ib = ibperij(i,j)
    !         jb = jbperij(i,j)
    !         FieldOut(ib,jb) = Fout(i)
    !      END DO
    !   END DO
    !   !$OMP END PARALLEL DO
    !ELSE
    !   !$OMP PARALLEL DO PRIVATE(ib,i,j)
       DO jb = 1, jbMax
          DO ib = 1, ibMax!PerJB(jb)
             !i = iPerIJB(ib,jb)
             !j = jPerIJB(ib,jb)
             FieldOut(ib,jb)=FieldIn(ib,jb)
          END DO
       END DO
    !   !$OMP END PARALLEL DO
    !END IF
 END SUBROUTINE AveBoxIJtoIBJB_R2D

 !
 ! maps (i,j) into (ib,jb)
 !
 SUBROUTINE IJtoIBJB_R(var_in,var_out)
    REAL(KIND=r8)   , INTENT(IN  ) :: var_in (:,:)
    REAL(KIND=r8)   , INTENT(OUT ) :: var_out(:,:)
    INTEGER               :: i
    INTEGER               :: j
    INTEGER               :: ib
    INTEGER               :: jb
   !$OMP PARALLEL DO PRIVATE(i,j,ib)
    DO jb = 1, jbMax
      DO ib = 1,ibMax
         var_out(ib,jb)=var_in(ib,jb)
       END DO
    END DO
   !$OMP END PARALLEL DO
 END SUBROUTINE IJtoIBJB_R

 ! 3D version by hmjb
 SUBROUTINE IJtoIBJB3_R(var_in,var_out)
    REAL(KIND=r8)   , INTENT(IN  ) :: var_in (iMax,kMax,jMax)
    REAL(KIND=r8)   , INTENT(OUT ) :: var_out(ibMax,kMax,jbMax)
    INTEGER               :: i
    INTEGER               :: j
    INTEGER               :: ib
    INTEGER               :: jb
   !$OMP PARALLEL DO PRIVATE(i,j,ib)
    DO jb = 1, jbMax
      DO ib = 1,ibMax
         var_out(ib,:,jb)=var_in(ib,:,jb)
       END DO
    END DO
   !$OMP END PARALLEL DO
 END SUBROUTINE IJtoIBJB3_R

 !
 ! maps (i,j) into (ib,jb)
 !
 SUBROUTINE IJtoIBJB_I(var_in,var_out)
    INTEGER(KIND=i8), INTENT(IN  ) :: var_in (:,:)
    INTEGER(KIND=i8), INTENT(OUT ) :: var_out(:,:)
    INTEGER               :: i
    INTEGER               :: j
    INTEGER               :: ib
    INTEGER               :: jb
   !$OMP PARALLEL DO PRIVATE(i,j,ib)
    DO jb = 1, jbMax
      DO ib = 1, ibMax
         var_out(ib,jb)=var_in(ib,jb)
       END DO
    END DO
   !$OMP END PARALLEL DO
 END SUBROUTINE IJtoIBJB_I
 ! 3D version by hmjb
 SUBROUTINE IJtoIBJB3_I(var_in,var_out)
    INTEGER(KIND=i8), INTENT(IN  ) :: var_in (iMax,kMax,jMax)
    INTEGER(KIND=i8), INTENT(OUT ) :: var_out(ibMax,kMax,jbMax)
    INTEGER               :: i
    INTEGER               :: j
    INTEGER               :: ib
    INTEGER               :: jb
   !$OMP PARALLEL DO PRIVATE(i,j,ib)
    DO jb = 1, jbMax
      DO ib = 1, ibMax
         var_out(ib,:,jb)=var_in(ib,:,jb)
       END DO
    END DO
   !$OMP END PARALLEL DO
 END SUBROUTINE IJtoIBJB3_I



SUBROUTINE tmstmp2(id    ,ifday ,tod   ,ihr   ,iday  ,mon   ,iyr)
  !
  ! $Author: cptec $
  ! $Date: 2000/01/10 22:27:56 $
  ! $Revision: 1.1 $
  !
  !==========================================================================
  !    id(4).......date of current data
  !                id(1)....hour(00/12)
  !                id(2)....month
  !                id(3)....day of month
  !                id(4)....year
  !    ifday.......model forecast day
  !    tod.........todx=tod+swint*f3600, model forecast time of 
  !                day in seconds
  !                swint....sw subr. call interval in hours
  !                swint has to be less than or equal to trint
  !                              and mod(trint,swint)=0 
  !                f3600=3.6e3
  !    ihr.........hour(00/12)
  !    iday........day of month
  !    mon.........month
  !    iyr.........year
  !    yrl.........length of year in days      
  !    monl(12)....length of each month in days 
  !==========================================================================
  !
  IMPLICIT NONE
  INTEGER, INTENT(in ) :: id(4)
  INTEGER, INTENT(in ) :: ifday
  REAL(KIND=r8),    INTENT(in ) :: tod
  INTEGER, INTENT(out) :: ihr
  INTEGER, INTENT(out) :: iday
  INTEGER, INTENT(out) :: mon
  INTEGER, INTENT(out) :: iyr

  INTEGER :: kday 
  INTEGER :: idaymn
  REAL(KIND=r8)    :: ctim 
  REAL(KIND=r8)    :: hrmodl
  INTEGER :: monl(12)           ! from common comtim

  REAL(KIND=r8), PARAMETER :: yrl =   365.2500_r8
  REAL(KIND=r8), PARAMETER ::  ep = 0.015625_r8
  DATA MONL/31,28,31,30,31,30,&
            31,31,30,31,30,31/

  ctim=tod+id(1)*3600.0_r8
  
  IF (ctim >= 86400.e0_r8) THEN
      kday=1
      ctim=ctim-86400.e0_r8
  ELSE
      kday=0
  END IF
  !
  !     adjust time to reduce round off error in divsion
  !
  iday = id(3) + ifday + kday
  hrmodl = (ctim+ep)/3600.0_r8
  ihr = hrmodl
  mon = id(2)
  iyr = id(4)
  DO
      idaymn = monl(mon)
      IF (yrl == 365.25e0_r8 .AND. MOD(iyr,4) == 0 .AND. mon == 2) &
          idaymn=29
      IF (iday <= idaymn) RETURN
      iday = iday - idaymn
      mon = mon + 1
      IF (mon < 13) CYCLE
      mon = 1
      iyr = iyr + 1
  END DO
END SUBROUTINE tmstmp2


       Subroutine InitTimeStamp(dateInit_s,idate) 
         Character(len = 10), intent(out)    ::  DateInit_s
         Integer, intent(in) :: idate(4) 

!local variables
         Integer hhi, mmi, ddi, yyyyi 

         yyyyi = idate(4) 
         mmi   = idate(2) 
         ddi   = idate(3) 
         hhi   = idate(1) 

         WRITE(DateInit_s,'(I4.4,3I2.2)') yyyyi, mmi, ddi, hhi

! computes the julian day of this calendar date 

         JulianDayInitIntegration = julday(mmi, ddi, yyyyi) 

       end Subroutine InitTimeStamp

       Subroutine TimeStamp(DateNow_s, idatec, jdt, dt)

         character(len =10) , intent(out) :: DateNow_s
         Integer,             intent(out) ::idatec(4)  

         Integer, intent(in)  :: jdt 
         real(KIND=r8),    intent(in)  :: dt

!local variables          
         Integer          :: hhc, mmc, ddc, yyyyc, juliandaynow    

         JulianDayNow = JulianDayInitIntegration + (int(dt)*jdt)/(24*3600)
         call caldat(JulianDayNow, mmc, ddc, yyyyc) 


!         write(*,*) juliandaynow
!         write(*,*) juliandayinitintegration 
!         write(*,*) mmc
!         write(*,*) ddc
!         write(*,*) yyyyc

         hhc = mod(int(dt)*jdt/3600,24) 

         write(DateNow_s,'(I4.4, 3I2.2)' ) yyyyc, mmc, ddc, hhc


         idatec = (/hhc,mmc,ddc,yyyyc /) 

       end Subroutine TimeStamp

       SUBROUTINE CALDAT(JULIAN,MM,ID,IYYY)
         IMPLICIT NONE
	 Integer,             intent(in) :: JULIAN
	 Integer,             intent(out) :: mm
	 Integer,             intent(out) :: ID
	 Integer,             intent(out) :: IYYY
         Integer,PARAMETER :: IGREG=2299161
         Integer  :: jalpha
         Integer  :: ja
         Integer  :: jb
         Integer  :: jc
         Integer  :: jd
         Integer  :: je

         ! input:  julian day 
         ! output: mm = mes ; id = dia, iyyy = ano 
         !PARAMETER (IGREG=2299161)
         IF(JULIAN.GE.IGREG)THEN
            JALPHA=INT(((JULIAN-1867216)-0.25)/36524.25)
            JA=JULIAN+1+JALPHA-INT(0.25*JALPHA)
         ELSE
            JA=JULIAN
         ENDIF
         JB=JA+1524
         JC=INT(6680.+((JB-2439870)-122.1)/365.25)
         JD=365*JC+INT(0.25*JC)
         JE=INT((JB-JD)/30.6001)
         ID=JB-JD-INT(30.6001*JE)
         MM=JE-1
         IF(MM.GT.12)MM=MM-12
         IYYY=JC-4715
         IF(MM.GT.2)IYYY=IYYY-1
         IF(IYYY.LE.0)IYYY=IYYY-1
       END SUBROUTINE CALDAT

 !________________________________________________________________________
       FUNCTION JULDAY(MM,ID,IYYY)
         IMPLICIT NONE
	 Integer,             intent(in) :: MM
	 Integer,             intent(in) :: ID
	 Integer,             intent(inout) :: IYYY
	 Integer,PARAMETER :: IGREG=15+31*(10+12*1582)
         Integer  :: julday
         Integer  :: jy
         Integer  :: jm
         Integer  :: ja

	 ! input:  mm = mes;  id = dia;  iyyy = ano 
         ! output: julian day dessa data. 
         
         !PARAMETER (IGREG=15+31*(10+12*1582))
         IF (IYYY.EQ.0) STOP 'There is no Year Zero.'
         IF (IYYY.LT.0) IYYY=IYYY+1
         IF (MM.GT.2) THEN
            JY=IYYY
            JM=MM+1
         ELSE
            JY=IYYY-1
            JM=MM+13
         ENDIF
         JULDAY=INT(365.25*JY)+INT(30.6001*JM)+ID+1720995
         IF (ID+31*(MM+12*IYYY).GE.IGREG) THEN
            JA=INT(0.01*JY)
            JULDAY=JULDAY+2-JA+INT(0.25*JA)
         ENDIF
       END FUNCTION JULDAY
!________________________________________________________________________




!
!------------------------------- VFORMAT ----------------------------------
!
 SUBROUTINE vfirec(iunit,a,n,type)

  INTEGER, INTENT(IN)  :: iunit  !#TO deve ser kind default
  INTEGER, INTENT(IN)  :: n
  REAL(KIND=r8), INTENT(OUT)    :: a(n)
  CHARACTER(len=* ), INTENT(IN) :: type
  !
  ! local
  !
  CHARACTER(len=1 ) :: vc(0:63)
  CHARACTER(len=80) :: line
  CHARACTER(len=1 ) :: cs
  INTEGER           :: ich0
  INTEGER           :: ich9
  INTEGER           :: ichcz
  INTEGER           :: ichca
  INTEGER           :: ichla
  INTEGER           :: ichlz
  INTEGER           :: i
  INTEGER           :: nvalline
  INTEGER           :: nchs
  INTEGER           :: ic
  INTEGER           :: ii
  INTEGER           :: isval
  INTEGER           :: iii
  INTEGER           :: ics
  INTEGER           :: nn
  INTEGER           :: nbits
  INTEGER           :: nc
  REAL(KIND=r8)              :: bias
  REAL(KIND=r8)              :: fact
  REAL(KIND=r8)              :: facti
  REAL(KIND=r8)              :: scfct

  vc='0'
  IF (vc(0).ne.'0') CALL vfinit(vc)

  ich0 =ichar('0')
  ich9 =ichar('9')
  ichcz=ichar('Z')
  ichlz=ichar('z')
  ichca=ichar('A')
  ichla=ichar('a')

  READ (iunit,'(2i8,2e20.10)')nn,nbits,bias,fact

  IF (nn.ne.n) THEN
    PRINT*,' Word count mismatch on vfirec record '
    PRINT*,' Words on record - ',nn
    PRINT*,' Words expected  - ',n
    STOP 'vfirec'
  END IF

  nvalline=(78*6)/nbits
  nchs=nbits/6

  DO i=1,n,nvalline
    READ(iunit,'(a78)') line
    ic=0
    DO ii=i,i+nvalline-1
      isval=0
      IF(ii.gt.n) EXIT
      DO iii=1,nchs
         ic=ic+1
         cs=line(ic:ic)
         ics=ichar(cs)
         IF (ics.le.ich9) THEN
            nc=ics-ich0
         ELSE IF (ics.le.ichcz) THEN
            nc=ics-ichca+10
         ELSE
            nc=ics-ichla+36
         END IF
         isval=ior(ishft(nc,6*(nchs-iii)),isval)
      END DO ! loop iii
        a(ii)=isval
    END DO ! loop ii

  END DO ! loop i

  facti=1.0_r8/fact

  IF (type.eq.'LIN') THEN
    DO i=1,n

      a(i)=a(i)*facti-bias

      !print*,'VFM=',i,a(i)
    END DO
  ELSE IF (type.eq.'LOG') THEN
    scfct=2.0_r8**(nbits-1)
    DO i=1,n
        a(i)=sign(1.0_r8,a(i)-scfct)  &
           *(10.0_r8**(abs(20.0_r8*(a(i)/scfct-1.0_r8))-10.0_r8))
    END DO
  END IF
 END SUBROUTINE vfirec
!--------------------------------------------------------
 SUBROUTINE vfinit(vc)
   CHARACTER(len=1), INTENT(OUT  ) :: vc   (*)
   CHARACTER(len=1)                :: vcscr(0:63)
   INTEGER                         :: n

   DATA vcscr/'0','1','2','3','4','5','6','7','8','9'   &
              ,'A','B','C','D','E','F','G','H','I','J'  &
              ,'K','L','M','N','O','P','Q','R','S','T'  &
              ,'U','V','W','X','Y','Z','a','b','c','d'  &
              ,'e','f','g','h','i','j','k','l','m','n'  &
              ,'o','p','q','r','s','t','u','v','w','x'  &
              ,'y','z','{','|'/

  DO n=0,63
      vc(n)=vcscr(n)
  END DO
 END SUBROUTINE vfinit


END MODULE Utils

