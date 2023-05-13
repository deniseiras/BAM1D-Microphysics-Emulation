     Module ModTimeStamp
      USE Constants, ONLY: rk,i8,r8,i4,r4
      IMPLICIT NONE
! date are always in the form yyyymmddhh ( year, month, day, hour). hour in 
! is in 0-24 form 

       Integer, save, private :: JulianDayInitIntegration  
contains 

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



     END Module ModTimeStamp
