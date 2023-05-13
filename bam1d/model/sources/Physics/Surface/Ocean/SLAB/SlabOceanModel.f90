MODULE SlabOceanModel
 IMPLICIT NONE  
  PRIVATE 
  ! Selecting Kinds
  INTEGER, PARAMETER :: r4 = SELECTED_REAL_KIND(6)  ! Kind for 32-bits Real Numbers
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)   ! Kind for 32-bits Integer Numbers
  INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers
  INTEGER, PARAMETER :: r16 = SELECTED_REAL_KIND(15)! Kind for 128-bits Real Numbers


  !  This is an example to use the ocean albedo look-up-table (ocnalbtab.bin) 
  !   to obtain the albedo at your specified spectral band, optical depth,
  !   cosine of solar zenith, wind and chlorophyll concentration.
  !  
  INTEGER, PARAMETER :: nb=24
  INTEGER, PARAMETER :: nt=16
  INTEGER, PARAMETER :: ns=15
  INTEGER, PARAMETER :: nw=7
  INTEGER, PARAMETER :: nc=5
  INTEGER, PARAMETER :: nline=8640
  REAL(KIND=r8), ALLOCATABLE  :: alb (:,:)
  REAL(KIND=r8), ALLOCATABLE  :: rflx(:,:)

  REAL(KIND=r8), PARAMETER  :: taunode(16)=(/0.0_r8, 0.05_r8, 0.1_r8, 0.16_r8, 0.24_r8,  0.35_r8, 0.5_r8, 0.7_r8, 0.99_r8, &
       1.3_r8, 1.80_r8, 2.5_r8, 5.00_r8, 9.00_r8, 15.0_r8, 25.0_r8 /)
  REAL(KIND=r8), PARAMETER  :: szanode(15)=(/0.05_r8, 0.09_r8, 0.15_r8, 0.21_r8, 0.27_r8, 0.33_r8, 0.39_r8, 0.45_r8, &
       0.52_r8, 0.60_r8, 0.68_r8, 0.76_r8, 0.84_r8, 0.92_r8, 1.0_r8 /)
  REAL(KIND=r8), PARAMETER  :: windnode(7)=(/0.0_r8, 3.0_r8, 6.0_r8, 9.0_r8, 12.0_r8, 15.0_r8, 18.0_r8 /)
  REAL(KIND=r8), PARAMETER  :: chlnode(5) =(/0.0_r8, 0.1_r8, 0.5_r8, 2.0_r8, 12.0_r8/)
  REAL(KIND=r8), PARAMETER  :: wlnode(25) =(/0.25_r8, 0.30_r8, 0.33_r8, 0.36_r8, 0.40_r8, 0.44_r8, 0.48_r8, 0.52_r8, &
       0.57_r8, 0.64_r8, 0.69_r8, 0.752_r8, 0.780_r8, 0.867_r8, 1.0_r8, 1.096_r8,&
       1.19_r8, 1.276_r8, 1.534_r8, 1.645_r8, 2.128_r8, 2.381_r8, 2.907_r8, &
       3.425_r8, 4.0_r8/)

 PUBLIC :: InitGetOceanAlb 
 PUBLIC :: GetOceanAlb
CONTAINS
  SUBROUTINE InitGetOceanAlb(fNameSlabOcen)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameSlabOcen

    !    !!!Note: if the record unit in your sestem is WORDS instead of BYTE,
    !    !!!      you MUST change the "recl" (record length) in the following
    !    !!!      OPEN statement from 24*4 to 24 (recl=24).        
    INTEGER :: lrec
    INTEGER :: irec
    REAL(KIND=r4)  :: alb2 (nb)
    ALLOCATE(alb (nline,nb));alb=0.0_r8
    ALLOCATE(rflx(nline,nb));rflx=0.0_r8
    alb2=0.0_r4
    INQUIRE(IOLENGTH=lrec)alb2
      !print*, 'lrec=',lrec, nb, fNameSlabOcen !delete this line Enver
    OPEN(1,file=TRIM(fNameSlabOcen),FORM='UNFORMATTED',ACCESS='DIRECT',recl=lrec)
      !print*, 'Enver after open'
    DO irec=1,nline
       READ(1,rec=irec)alb2
       !print*, 'Enver after ',irec,nline,'read fNameSlabOcen'
       rflx(irec,1:nb)=alb2(1:nb)
       alb (irec,1:nb)=alb2(1:nb)
    END DO
    CLOSE(1)
    RETURN

  END SUBROUTINE InitGetOceanAlb


  FUNCTION GetOceanAlb(tau_in,sza_in,wind_in,chl_in,wls_in,wle_in)

    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN   ) :: tau_in 
    REAL(KIND=r8), INTENT(IN   ) :: sza_in 
    REAL(KIND=r8), INTENT(IN   ) :: wind_in
    REAL(KIND=r8), INTENT(IN   ) :: chl_in 
    REAL(KIND=r8), INTENT(IN   ) :: wls_in 
    REAL(KIND=r8), INTENT(IN   ) :: wle_in 
    REAL(KIND=r8) :: GetOceanAlb
    REAL(KIND=r8) :: albb  (nb)
    REAL(KIND=r8) :: wtflx (nb)
    REAL(KIND=r8) :: rt(2)
    REAL(KIND=r8) :: rs(2)
    REAL(KIND=r8) :: rw(2)
    REAL(KIND=r8) :: rc(2)
    REAL(KIND=r8) :: tau
    REAL(KIND=r8) :: sza
    REAL(KIND=r8) :: wind
    REAL(KIND=r8) :: chl
    REAL(KIND=r8) :: wls 
    REAL(KIND=r8) :: wle 
    REAL(KIND=r8) :: twtb
    REAL(KIND=r8) :: wtb
    REAL(KIND=r8) :: wtf
    REAL(KIND=r8) :: wtt
    REAL(KIND=r8) :: albedo
    INTEGER :: ib
    INTEGER :: ib1
    INTEGER :: ib2
    INTEGER :: icc
    INTEGER :: ic
    INTEGER :: irec
    INTEGER :: irr
    INTEGER :: iss
    INTEGER :: is
    INTEGER :: itt
    INTEGER :: it
    INTEGER :: iww
    INTEGER :: iw
    !     --------------------------- INPUT ---------------------------------------
    !   
    !    specify the parameters for albedo here:

    !        tau = 0.40              !aerosol/cloud optical depth
    !        sza = 0.50              !cosine of solar zenith angle
    !        wind = 10.00            !wind speed in m/s
    !        chl = 0.10              !chlorophyll concentration in mg/m3
    !        wls = 0.69              !start wavelength (um) of your band
    !        wle = 1.19              !end wavelength (um) of your band
    !
    !   Output albedo for these input is 0.068. If not, your system may use
    !   different record unit and you need change the record length.
    !     -------------------------------------------------------------------------

    !    now find the albedo corresponding to the 4 parameters above:
    !-----------------------------------------------------------------------
    ! Computes surface albedos over ocean for Slab Ocean Model (SOM)
    !
    ! Two spectral surface albedos for direct (dir) and diffuse (dif)
    ! incident radiation are calculated. The spectral intervals are:
    !   s (shortwave)  = 0.2-0.7 micro-meters
    !   l (longwave)   = 0.7-5.0 micro-meters
    !

    !real asdir(plond)     ! Srf alb for direct  rad   0.2-0.7 micro-ms (0.20-0.69 )
    !real asdif(plond)     ! Srf alb for diffuse rad   0.2-0.7 micro-ms (1.19-2.38 )

    !real aldir(plond)     ! Srf alb for direct rad   0.7-5.0 micro-ms (0.69-1.19 )
    !real aldif(plond)     ! Srf alb for diffuse rad  0.7-5.0 micro-ms (2.38-4.0um)
    GetOceanAlb=0.8e0_r8
    albb  (1:nb)=0.0_r8
    wtflx (1:nb)=0.0_r8
    rt(1:2)=0.0_r8
    rs(1:2)=0.0_r8
    rw(1:2)=0.0_r8
    rc(1:2)=0.0_r8
    tau=0.0_r8
    sza=0.0_r8
    wind=0.0_r8
    chl=0.0_r8
    wls =0.0_r8
    wle =0.0_r8
    twtb=0.0_r8
    wtb=0.0_r8
    wtf=0.0_r8
    wtt=0.0_r8
    albedo=0.0_r8
    ib=0
    ib1=0
    ib2=0
    icc=0
    ic=0
    irec=0
    irr=0
    iss=0
    is=0
    itt=0
    it=0
    iww=0
    iw=0
    
    tau = tau_in 
    sza = sza_in 
    wind= wind_in
    chl = chl_in 
    wls = wls_in 
    wle = wle_in 
    PRINT*,tau.LT.0.0_r8 , sza.LT.0.0_r8 , sza.GT.1.0_r8, wind.LT.0.0_r8,chl.LT.0.0_r8 
    PRINT*,tau,sza,wind,chl
    IF(tau.LT.0.0_r8 .OR. (sza.LT.0.0_r8 .OR. sza.GT.1.0_r8) .OR. wind.LT.0.0_r8 &
         .OR. chl.LT.0.0_r8)STOP 'Err: input parameters wrong!'
    IF(tau  .GT. taunode(nt))tau=taunode(nt)
    IF(wind .GT. windnode(nw))wind=windnode(nw)
    IF(chl  .GT. 45.0_r8)chl=45.0_r8
    CALL locate(taunode,nt,tau,it)
    CALL locate(szanode,ns,sza,is)
    CALL locate(windnode,nw,wind,iw)
    CALL locate(chlnode,nc,chl,ic)
    rt(2) = (tau-taunode(it))/(taunode(it+1)-taunode(it))
    rt(1) = 1.0_r8-rt(2)
    rs(2) = (sza-szanode(is))/(szanode(is+1)-szanode(is))
    rs(1) = 1.0_r8-rs(2)
    rw(2) = (wind-windnode(iw))/(windnode(iw+1)-windnode(iw))
    rw(1) = 1.0_r8-rw(2)
    rc(2) = (chl-chlnode(ic))/(chlnode(ic+1)-chlnode(ic))
    rc(1) = 1.0_r8-rc(2)

    IF(wls .GT. wle)STOP 'Err: Start wavelength should be smaller.'
    IF(wls .LT. wlnode(1))wls=wlnode(1)
    IF(wle .GT. wlnode(25))wle=wlnode(25)
    CALL locate(wlnode,25,wls,ib1)
    CALL locate(wlnode,25,wle,ib2)

    !             ** get alb(ib1:ib2) by 4 dimensional linear interpolaton **
    DO ib=ib1,ib2
       albb(ib) = 0.0_r8
       wtflx(ib) = 0.0_r8
    ENDDO
    DO  itt=it,it+1
       DO  iss=is,is+1
          DO  iww=iw,iw+1
             DO  icc=ic,ic+1
                irec = (itt-1)*15*7*5 + (iss-1)*7*5 + (iww-1)*5 + icc
                !read(1,rec=irec) alb
                wtt = rt(itt-it+1)*rs(iss-is+1)*rw(iww-iw+1)*rc(icc-ic+1)
                DO ib=ib1,ib2
                   albb(ib) = albb(ib) + wtt*alb(irec,ib)
                ENDDO
             END DO
          END DO

          !                  *** get 24 band down flux weights ***
          irr = 8400 + (itt-1)*15 + iss
          !read(1,rec=irr) rflx

          wtf = rt(itt-it+1)*rs(iss-is+1)
          !         do ib=ib1,ib2
          DO ib=1,24
             wtflx(ib) = wtflx(ib) + wtf*rflx(irr,ib)
          ENDDO
       END DO
    END DO

    !             ** get albedo in the specified band by weighted sum **
    twtb = 0.0_r8
    albedo = 0.0_r8
    DO ib=ib1,ib2
       IF(ib .EQ. ib1)THEN
          wtb = (wlnode(ib1+1)-wls)/(wlnode(ib1+1)-wlnode(ib1))
       ELSE IF (ib .EQ. ib2)THEN
          wtb = (wle-wlnode(ib2))/(wlnode(ib2+1)-wlnode(ib2))
       ELSE
          wtb = 1.0_r8
       ENDIF
       albedo = albedo + wtb*wtflx(ib)*albb(ib)
       twtb = twtb + wtb*wtflx(ib)
    ENDDO
    IF (twtb.EQ.0.0_r8 .AND. ib1.EQ.ib2)THEN
       albedo = albb(ib1)
    ELSE
       albedo = albedo/twtb
    ENDIF
    GetOceanAlb=albedo
    !WRITE(*,'(6a6,/,6f6.2,/,a11,f6.3,/)')'Tau','COSUN','Wind', &
    !     'Chl','WL1','WL2',tau,sza,wind,chl,wls,wle,&
    !     'Albedo =',albedo

  END FUNCTION GetOceanAlb

  !=======================================================================
  SUBROUTINE locate(xx,n,x,j)
    IMPLICIT NONE
    !
    ! purpose:  given an array xx of length n, and given a value X, returns
    !           a value J such that X is between xx(j) and xx(j+1). xx must
    !           be monotonic, either increasing of decreasing. this function
    !           returns j=1 or j=n-1 if x is out of range.
    !
    ! input:
    !   xx      monitonic table
    !   n       size of xx
    !   x       single floating point value perhaps within the range of xx
    !
    INTEGER, INTENT(IN   ) :: n
    REAL(KIND=r8), INTENT(IN   ) :: x,xx(n)
    INTEGER, INTENT(OUT  ) :: j

    INTEGER :: jl,jm,ju
    j=0;jl=0;jm=0;ju=0
    
    IF(x.EQ.xx(1)) THEN
       j=1
       RETURN
    ENDIF
    IF(x.EQ.xx(n)) THEN
       j=n-1
       RETURN
    ENDIF
    jl=1
    ju=n
10  IF(ju-jl.GT.1) THEN
       jm=(ju+jl)/2
       IF((xx(n).GT.xx(1)).EQV.(x.GT.xx(jm)))THEN
          jl=jm
       ELSE
          ju=jm
       ENDIF
       GOTO 10
    ENDIF
    j=jl
    RETURN
  END SUBROUTINE locate
  !=======================================================================

END MODULE SlabOceanModel
