MODULE TimesTmp
 USE Constants, ONLY: rk,i8,r8,i4,r4

 IMPLICIT NONE


CONTAINS
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
  INTEGER(KIND=i8), INTENT(in ) :: id(4)
  INTEGER(KIND=i8), INTENT(in ) :: ifday
  REAL(KIND=r8),    INTENT(in ) :: tod
  INTEGER(KIND=i8), INTENT(out) :: ihr
  INTEGER(KIND=i8), INTENT(out) :: iday
  INTEGER(KIND=i8), INTENT(out) :: mon
  INTEGER(KIND=i8), INTENT(out) :: iyr

  INTEGER(KIND=i8) :: kday 
  INTEGER(KIND=i8) :: idaymn
  REAL(KIND=r8)    :: ctim 
  REAL(KIND=r8)    :: hrmodl
  INTEGER(KIND=i8) :: monl(12)           ! from common comtim

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
END MODULE TimesTmp
