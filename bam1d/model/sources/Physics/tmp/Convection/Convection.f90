MODULE wv_saturation
 IMPLICIT NONE
    INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
    REAL(KIND=r8),PARAMETER :: SHR_CONST_MWWV   = 18.016_r8       ! molecular weight water vapor
    REAL(KIND=r8),PARAMETER :: SHR_CONST_MWDAIR = 28.966_r8       ! molecular weight dry air ~ kg/kmole
    REAL(KIND=r8),PARAMETER :: epsilo = shr_const_mwwv/shr_const_mwdair ! ratio of h2o to dry air molecular weights 
    REAL(KIND=r8),PARAMETER :: SHR_CONST_LATVAP = 2.501e6_r8      ! latent heat of evaporation ~ J/kg
    REAL(KIND=r8),PARAMETER :: SHR_CONST_TKFRZ  = 273.16_r8       ! freezing T of fresh water ~ K (intentionally made == to TKTRIP)
    REAL(KIND=r8),PARAMETER :: SHR_CONST_AVOGAD = 6.02214e26_r8   ! Avogadro's number ~ molecules/kmole
    REAL(KIND=r8),PARAMETER :: SHR_CONST_BOLTZ  = 1.38065e-23_r8  ! Boltzmann's constant ~ J/K/molecule
    REAL(KIND=r8),PARAMETER :: SHR_CONST_LATICE = 3.337e5_r8      ! latent heat of fusion ~ J/kg
    REAL(KIND=r8), PARAMETER :: latice = shr_const_latice ! Latent heat of fusion

    REAL(KIND=r8),PARAMETER :: trice  =  20.00_r8         ! Trans range from es over h2o to es over ice
    REAL(KIND=r8),PARAMETER :: ttrice=trice
    REAL(KIND=r8),PARAMETER :: SHR_CONST_CPDAIR = 1.00464e3_r8    ! specific heat of dry air ~ J/kg/K
    REAL(KIND=r8),PARAMETER :: cpair = shr_const_cpdair  ! specific heat of dry air (J/K/kg)
    REAL(KIND=r8),PARAMETER :: cp    =cpair
    REAL(KIND=r8),PARAMETER :: latvap = shr_const_latvap ! Latent heat of vaporization
    REAL(KIND=r8),PARAMETER :: SHR_CONST_RGAS   = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ ! Universal gas constant ~ J/K/kmole
    REAL(KIND=r8),PARAMETER :: SHR_CONST_RWV    = SHR_CONST_RGAS/SHR_CONST_MWWV    ! Water vapor gas constant ~ J/K/kg
    INTEGER, PARAMETER :: plenest = 250! length of saturation vapor pressure table
    REAL(KIND=r8), PUBLIC, PARAMETER :: tmelt = shr_const_tkfrz   ! Freezing point of water

    REAL(KIND=r8)           :: estbl(plenest)      ! table values of saturation vapor pressure
    REAL(KIND=r8),PARAMETER :: tmn  = 173.16_r8          ! Minimum temperature entry in table
    REAL(KIND=r8),PARAMETER :: tmx  = 375.16_r8          ! Maximum temperature entry in table
    REAL(KIND=r8),PARAMETER :: tmin=tmn       ! min temperature (K) for table
    REAL(KIND=r8),PARAMETER :: tmax= tmx      ! max temperature (K) for table
    LOGICAL ,PARAMETER :: icephs=.TRUE.  ! false => saturation vapor press over water only
    INTEGER ,PARAMETER ::  iterp =2             ! #iterations for precipitation calculation
    INTEGER  ::  k1mb  =1  ! index of the eta level near 1 mb
    REAL(KIND=r8)           :: pcf(6)     ! polynomial coeffs -> es transition water to ice
    REAL(KIND=r8), PRIVATE :: hlatf  = latice
    REAL(KIND=r8), PRIVATE :: hlatv  = latvap
    REAL(KIND=r8), PRIVATE :: rgasv  = SHR_CONST_RWV    ! Gas constant for water vapor
    REAL(KIND=r8), PRIVATE :: t0 = tmelt                ! approximate freezing temp

CONTAINS

  SUBROUTINE findsp (nCols,pver, q, t, p, tsp, qsp)
    IMPLICIT NONE
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    !     find the wet bulb temperature for a given t and q
    !     in a longitude height section
    !     wet bulb temp is the temperature and spec humidity that is 
    !     just saturated and has the same enthalpy
    !     if q > qs(t) then tsp > t and qsp = qs(tsp) < q
    !     if q < qs(t) then tsp < t and qsp = qs(tsp) > q
    !
    ! Method: 
    ! a Newton method is used
    ! first guess uses an algorithm provided by John Petch from the UKMO
    ! we exclude points where the physical situation is unrealistic
    ! e.g. where the temperature is outside the range of validity for the
    !      saturation vapor pressure, or where the water vapor pressure
    !      exceeds the ambient pressure, or the saturation specific humidity is 
    !      unrealistic
    ! 
    ! Author: P. Rasch
    ! 
    !-----------------------------------------------------------------------
    !
    !     input arguments
    !
    INTEGER, INTENT(in) :: nCols                 ! number of columns (max)
    INTEGER, INTENT(in) :: pver                  ! number of vertical levels

    REAL(KIND=r8), INTENT(in) :: q(nCols,pver)        ! water vapor (kg/kg)
    REAL(KIND=r8), INTENT(in) :: t(nCols,pver)        ! temperature (K)
    REAL(KIND=r8), INTENT(in) :: p(nCols,pver)        ! pressure    (Pa)
    !
    ! output arguments
    !
    REAL(KIND=r8), INTENT(out) :: tsp(nCols,pver)      ! saturation temp (K)
    REAL(KIND=r8), INTENT(out) :: qsp(nCols,pver)      ! saturation mixing ratio (kg/kg)
    !
    ! local variables
    !
    INTEGER i                 ! work variable
    INTEGER k                 ! work variable
    LOGICAL lflg              ! work variable
    INTEGER iter              ! work variable
    INTEGER l                 ! work variable
    LOGICAL :: error_found

    REAL(KIND=r8) omeps                ! 1 minus epsilon
    REAL(KIND=r8) trinv                ! work variable
    REAL(KIND=r8) es                   ! sat. vapor pressure
    REAL(KIND=r8) desdt                ! change in sat vap pressure wrt temperature
    !     real(KIND=r8) desdp                ! change in sat vap pressure wrt pressure
    REAL(KIND=r8) dqsdt                ! change in sat spec. hum. wrt temperature
    REAL(KIND=r8) dgdt                 ! work variable
    REAL(KIND=r8) g                    ! work variable
    REAL(KIND=r8) weight(nCols)        ! work variable
    REAL(KIND=r8) hlatsb               ! (sublimation)
    REAL(KIND=r8) hlatvp               ! (vaporization)
    REAL(KIND=r8) hltalt(nCols,pver)   ! lat. heat. of vap.
    REAL(KIND=r8) tterm                ! work var.
    REAL(KIND=r8) qs                   ! spec. hum. of water vapor
    REAL(KIND=r8) tc                   ! crit temp of transition to ice

    ! work variables
    REAL(KIND=r8) t1, q1, dt, dq
    REAL(KIND=r8) dtm, dqm
    REAL(KIND=r8) qvd, a1, tmp
    REAL(KIND=r8) rair
    REAL(KIND=r8) r1b, c1, c2, c3
    REAL(KIND=r8) denom
    REAL(KIND=r8) dttol
    REAL(KIND=r8) dqtol
    INTEGER doit(nCols) 
    REAL(KIND=r8) enin(nCols), enout(nCols)
    REAL(KIND=r8) tlim(nCols)
    REAL(KIND=r8) epsqs
    epsqs = epsilo
    k1mb=1
    omeps = 1.0_r8 - epsqs
    trinv = 1.0_r8/ttrice
    a1 = 7.5_r8*LOG(10.0_r8)
    rair =  287.04_r8
    c3 = rair*a1/cp
    dtm = 0.0_r8    ! needed for iter=0 blowup with f90 -ei
    dqm = 0.0_r8    ! needed for iter=0 blowup with f90 -ei
    dttol = 1.e-4_r8 ! the relative temp error tolerance required to quit the iteration
    dqtol = 1.e-4_r8 ! the relative moisture error tolerance required to quit the iteration
    !  tmin = 173.16_r8 ! the coldest temperature we can deal with
    !
    ! max number of times to iterate the calculation
    iter = 50
    !
    DO k = k1mb,pver

       !
       ! first guess on the wet bulb temperature
       !
       DO i = 1,ncols

          ! limit the temperature range to that relevant to the sat vap pres tables

          tlim(i) = MIN(MAX(t(i,k),173.0_r8),373.0_r8)

          es = estblf(tlim(i))

          denom = p(i,k) - omeps*es
          qs = epsqs*es/denom
          doit(i) = 0
          enout(i) = 1.0_r8
          ! make sure a meaningful calculation is possible
          IF (p(i,k) > 5.0_r8*es .AND. qs > 0.0_r8 .AND. qs < 0.5_r8) THEN
             !
             ! Saturation specific humidity
             !
             qs = MIN(epsqs*es/denom,1.0_r8)
             !
             ! "generalized" analytic expression for t derivative of es
             !  accurate to within 1 percent for 173.16 < t < 373.16
             !
             ! Weighting of hlat accounts for transition from water to ice
             ! polynomial expression approximates difference between es over
             ! water and es over ice from 0 to -ttrice (C) (min of ttrice is
             ! -40): required for accurate estimate of es derivative in transition
             ! range from ice to water also accounting for change of hlatv with t
             ! above freezing where const slope is given by -2369 j/(kg c) = cpv - cw
             !
             tc     = tlim(i) - t0
             lflg   = (tc >= -ttrice .AND. tc < 0.0_r8)
             weight(i) = MIN(-tc*trinv,1.0_r8)
             hlatsb = hlatv + weight(i)*hlatf
             hlatvp = hlatv - 2369.0_r8*tc
             IF (tlim(i) < t0) THEN
                hltalt(i,k) = hlatsb
             ELSE
                hltalt(i,k) = hlatvp
             END IF
             enin(i) = cp*tlim(i) + hltalt(i,k)*q(i,k)

             ! make a guess at the wet bulb temp using a UKMO algorithm (from J. Petch)
             tmp =  q(i,k) - qs
             c1 = hltalt(i,k)*c3
             c2 = (tlim(i) + 36.0_r8)**2
             r1b    = c2/(c2 + c1*qs)
             qvd   = r1b*tmp
             tsp(i,k) = tlim(i) + ((hltalt(i,k)/cp)*qvd)
             !#ifdef DEBUG
             !             if ( (lchnk == lchnklook(nlook) ) .and. (i == icollook(nlook) ) ) then
             !                write (6,*) ' relative humidity ', q(i,k)/qs
             !                write (6,*) ' first guess ', tsp(i,k)
             !             endif
             !#endif
             es = estblf(tsp(i,k))
             qsp(i,k) = MIN(epsqs*es/(p(i,k) - omeps*es),1.0_r8)
          ELSE
             doit(i) = 1
             tsp(i,k) = tlim(i)
             qsp(i,k) = q(i,k)
             enin(i) = 1.0_r8
          ENDIF
       END DO   ! end do i
       !
       ! now iterate on first guess
       !
       DO l = 1, iter
          dtm = 0.0_r8
          dqm = 0.0_r8
          DO i = 1,ncols
             IF (doit(i) == 0) THEN
                es = estblf(tsp(i,k))
                !
                ! Saturation specific humidity
                !
                qs = MIN(epsqs*es/(p(i,k) - omeps*es),1.0_r8)
                !
                ! "generalized" analytic expression for t derivative of es
                ! accurate to within 1 percent for 173.16 < t < 373.16
                !
                ! Weighting of hlat accounts for transition from water to ice
                ! polynomial expression approximates difference between es over
                ! water and es over ice from 0 to -ttrice (C) (min of ttrice is
                ! -40): required for accurate estimate of es derivative in transition
                ! range from ice to water also accounting for change of hlatv with t
                ! above freezing where const slope is given by -2369 j/(kg c) = cpv - cw
                !
                tc     = tsp(i,k) - t0
                lflg   = (tc >= -ttrice .AND. tc < 0.0_r8)
                weight(i) = MIN(-tc*trinv,1.0_r8)
                hlatsb = hlatv + weight(i)*hlatf
                hlatvp = hlatv - 2369.0_r8*tc
                IF (tsp(i,k) < t0) THEN
                   hltalt(i,k) = hlatsb
                ELSE
                   hltalt(i,k) = hlatvp
                END IF
                IF (lflg) THEN
                   tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3)+tc*(pcf(4) + tc*pcf(5))))
                ELSE
                   tterm = 0.0_r8
                END IF
                desdt = hltalt(i,k)*es/(rgasv*tsp(i,k)*tsp(i,k)) + tterm*trinv
                dqsdt = (epsqs + omeps*qs)/(p(i,k) - omeps*es)*desdt
                !              g = cp*(tlim(i)-tsp(i,k)) + hltalt(i,k)*q(i,k)- hltalt(i,k)*qsp(i,k)
                g = enin(i) - (cp*tsp(i,k) + hltalt(i,k)*qsp(i,k))
                dgdt = -(cp + hltalt(i,k)*dqsdt)
                t1 = tsp(i,k) - g/dgdt
                dt = ABS(t1 - tsp(i,k))/t1
                tsp(i,k) = MAX(t1,tmin)
                es = estblf(tsp(i,k))
                q1 = MIN(epsqs*es/(p(i,k) - omeps*es),1.0_r8)
                q1=MAX(q1,1.e-12_r8)
                dq = ABS(q1 - qsp(i,k))/MAX(q1,1.e-12_r8)
                qsp(i,k) = q1
                !#ifdef DEBUG
                !               if ( (lchnk == lchnklook(nlook) ) .and. (i == icollook(nlook) ) ) then
                !                  write (6,*) ' rel chg lev, iter, t, q ', k, l, dt, dq, g
                !               endif
                !#endif
                dtm = MAX(dtm,dt)
                dqm = MAX(dqm,dq)
                ! if converged at this point, exclude it from more iterations
                IF (dt < dttol .AND. dq < dqtol) THEN
                   doit(i) = 2
                ENDIF
                enout(i) = cp*tsp(i,k) + hltalt(i,k)*qsp(i,k)
                ! bail out if we are too near the end of temp range
                IF (tsp(i,k) < 174.16_r8) THEN
                   doit(i) = 4
                ENDIF
             ELSE
             ENDIF
          END DO              ! do i = 1,ncols

          IF (dtm < dttol .AND. dqm < dqtol) THEN
             go to 10
          ENDIF

       END DO                 ! do l = 1,iter
10     CONTINUE

       error_found = .FALSE.
       IF (dtm > dttol .OR. dqm > dqtol) THEN
          DO i = 1,ncols
             IF (doit(i) == 0) error_found = .TRUE.
          END DO
          IF (error_found) THEN
             DO i = 1,ncols
                IF (doit(i) == 0) THEN
                   WRITE (6,*) ' findsp not converging at point i, k ', i, k
                   WRITE (6,*) ' t, q, p, enin ', t(i,k), q(i,k), p(i,k), enin(i)
                   WRITE (6,*) ' tsp, qsp, enout ', tsp(i,k), qsp(i,k), enout(i)
                   STOP 'FINDSP'
                ENDIF
             END DO
          ENDIF
       ENDIF
       DO i = 1,ncols
          IF (doit(i) == 2 .AND. ABS((enin(i)-enout(i))/(enin(i)+enout(i))) > 1.e-4_r8) THEN
             error_found = .TRUE.
          ENDIF
       END DO
       IF (error_found) THEN
          DO i = 1,ncols
             IF (doit(i) == 2 .AND. ABS((enin(i)-enout(i))/(enin(i)+enout(i))) > 1.e-4_r8) THEN
                WRITE (6,*) ' the enthalpy is not conserved for point ', &
                     i, k, enin(i), enout(i)
                WRITE (6,*) ' t, q, p, enin ', t(i,k), q(i,k), p(i,k), enin(i)
                WRITE (6,*) ' tsp, qsp, enout ', tsp(i,k), qsp(i,k), enout(i)
                STOP 'FINDSP'
             ENDIF
          END DO
       ENDIF

    END DO                    ! level loop (k=1,pver)

    RETURN
  END SUBROUTINE findsp



  SUBROUTINE gestbl(       )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Builds saturation vapor pressure table for later lookup procedure.
    ! 
    ! Method: 
    ! Uses Goff & Gratch (1946) relationships to generate the table
    ! according to a set of free parameters defined below.  Auxiliary
    ! routines are also included for making rapid estimates (well with 1%)
    ! of both es and d(es)/dt for the particular table configuration.
    ! 
    ! Author: J. Hack
    ! 
    !-----------------------------------------------------------------------
    !   use pmgrid, only: masterproc
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    !
    !---------------------------Local variables-----------------------------
    !
    REAL(KIND=r8)  epsqs
    REAL(KIND=r8) t             ! Temperature
    INTEGER n          ! Increment counter
    INTEGER lentbl     ! Calculated length of lookup table
    INTEGER itype      ! Ice phase: 0 -> no ice phase
    !            1 -> ice phase, no transitiong
    !           -x -> ice phase, x degree transition
    !
    !-----------------------------------------------------------------------
    !
    ! Set es table parameters
    !
    !   tmin   = tmn       ! Minimum temperature entry in table
    !   tmax   = tmx       ! Maximum temperature entry in table
    !   ttrice = trice     ! Trans. range from es over h2o to es over ice
    !   icephs = ip        ! Ice phase (true or false)
    !
    ! Set physical constants required for es calculation
    !
    epsqs  = epsilo
    !
    lentbl = INT(tmax-tmin+2.000001_r8)
    IF (lentbl .GT. plenest) THEN
       WRITE(6,9000) tmax, tmin, plenest
       STOP 'GESTBL'    ! Abnormal termination
    END IF
    !
    ! Begin building es table.
    ! Check whether ice phase requested.
    ! If so, set appropriate transition range for temperature
    !
    IF (icephs) THEN
       IF (ttrice /= 0.0_r8) THEN
          itype = -ttrice
       ELSE
          itype = 1
       END IF
    ELSE
       itype = 0
    END IF
    !
    t = tmin - 1.0_r8
    DO n=1,lentbl
       t = t + 1.0_r8
       CALL gffgch(t,estbl(n),itype)
    END DO
    !
    DO n=lentbl+1,plenest
       estbl(n) = -99999.0_r8
    END DO
    !
    ! Table complete -- Set coefficients for polynomial approximation of
    ! difference between saturation vapor press over water and saturation
    ! pressure over ice for -ttrice < t < 0 (degrees C). NOTE: polynomial
    ! is valid in the range -40 < t < 0 (degrees C).
    !
    !                  --- Degree 5 approximation ---
    !
    pcf(1) =  5.04469588506e-01_r8
    pcf(2) = -5.47288442819e+00_r8
    pcf(3) = -3.67471858735e-01_r8
    pcf(4) = -8.95963532403e-03_r8
    pcf(5) = -7.78053686625e-05_r8
    !
    !                  --- Degree 6 approximation ---
    !
    !-----pcf(1) =  7.63285250063e-02
    !-----pcf(2) = -5.86048427932e+00
    !-----pcf(3) = -4.38660831780e-01
    !-----pcf(4) = -1.37898276415e-02
    !-----pcf(5) = -2.14444472424e-04
    !-----pcf(6) = -1.36639103771e-06
    !
    !   if (masterproc) then
    !      write(6,*)' *** SATURATION VAPOR PRESSURE TABLE COMPLETED ***'
    !   end if

    RETURN
    !
9000 FORMAT('GESTBL: FATAL ERROR *********************************',/, &
         ' TMAX AND TMIN REQUIRE A LARGER DIMENSION ON THE LENGTH', &
         ' OF THE SATURATION VAPOR PRESSURE TABLE ESTBL(PLENEST)',/, &
         ' TMAX, TMIN, AND PLENEST => ', 2f7.2, i3)
    !
  END SUBROUTINE gestbl
  
  
  REAL(KIND=r8) FUNCTION estblf( td )
    !
    ! Saturation vapor pressure table lookup
    !
    REAL(KIND=r8), INTENT(in) :: td         ! Temperature for saturation lookup
    !
    REAL(KIND=r8) :: e          ! intermediate variable for es look-up
    REAL(KIND=r8) :: tmin       ! min temperature (K) for table
    REAL(KIND=r8) :: tmax       ! max temperature (K) for table

    REAL(KIND=r8) :: ai
    INTEGER  :: i

    tmin=tmn  
    tmax=tmx  

    !
    e = MAX(MIN(td,tmax),tmin)   ! partial pressure
    i = INT(e-tmin)+1
    ai = AINT(e-tmin)
    estblf = (tmin+ai-e+1.0_r8)* &
         estbl(i)-(tmin+ai-e)* &
         estbl(i+1)
  END FUNCTION estblf
  
  
  SUBROUTINE aqsat(t       ,p       ,es      ,qs        ,ii      , &
       ILEN    ,kk      ,kstart  ,kend      )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Utility procedure to look up and return saturation vapor pressure from
    ! precomputed table, calculate and return saturation specific humidity
    ! (g/g),for input arrays of temperature and pressure (dimensioned ii,kk)
    ! This routine is useful for evaluating only a selected region in the
    ! vertical.
    ! 
    ! Method: 
    ! <Describe the algorithm(s) used in the routine.> 
    ! <Also include any applicable external references.> 
    ! 
    ! Author: J. Hack
    ! 
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    INTEGER  , INTENT(in) :: ii             ! I dimension of arrays t, p, es, qs
    INTEGER  , INTENT(in) :: kk             ! K dimension of arrays t, p, es, qs
    INTEGER  , INTENT(in) :: ILEN           ! Length of vectors in I direction which
    INTEGER  , INTENT(in) :: kstart         ! Starting location in K direction
    INTEGER  , INTENT(in) :: kend           ! Ending location in K direction
    REAL(KIND=r8), INTENT(in) :: t(ii,kk)          ! Temperature
    REAL(KIND=r8), INTENT(in) :: p(ii,kk)          ! Pressure
    !
    ! Output arguments
    !
    REAL(KIND=r8), INTENT(out) :: es(ii,kk)         ! Saturation vapor pressure
    REAL(KIND=r8), INTENT(out) :: qs(ii,kk)         ! Saturation specific humidity
    !
    !---------------------------Local workspace-----------------------------
    !
    REAL(KIND=r8) epsqs      ! Ratio of h2o to dry air molecular weights 
    REAL(KIND=r8) omeps             ! 1 - 0.622
    INTEGER i, k           ! Indices
    !
    !-----------------------------------------------------------------------
    !
    epsqs = epsilo
    omeps = 1.0_r8 - epsqs
    DO k=kstart,kend
       DO i=1,ILEN
          es(i,k) = estblf(t(i,k))
          !
          ! Saturation specific humidity
          !
          qs(i,k) = epsqs*es(i,k)/(p(i,k) - omeps*es(i,k))
          !
          ! The following check is to avoid the generation of negative values
          ! that can occur in the upper stratosphere and mesosphere
          !
          qs(i,k) = MIN(1.0_r8,qs(i,k))
          !
          IF (qs(i,k) < 0.0_r8) THEN
             qs(i,k) = 1.0_r8
             es(i,k) = p(i,k)
          END IF
       END DO
    END DO
    !
    RETURN
  END SUBROUTINE aqsat
  SUBROUTINE vqsatd(t       ,p       ,es      ,qs      ,gam      , &
       len     )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Utility procedure to look up and return saturation vapor pressure from
    ! precomputed table, calculate and return saturation specific humidity
    ! (g/g), and calculate and return gamma (l/cp)*(d(qsat)/dT).  The same
    ! function as qsatd, but operates on vectors of temperature and pressure
    ! 
    ! Method: 
    ! 
    ! Author: J. Hack
    ! 
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    INTEGER, INTENT(in) :: len       ! vector length
    REAL(KIND=r8), INTENT(in) :: t(len)       ! temperature
    REAL(KIND=r8), INTENT(in) :: p(len)       ! pressure
    !
    ! Output arguments
    !
    REAL(KIND=r8), INTENT(out) :: es(len)   ! saturation vapor pressure
    REAL(KIND=r8), INTENT(out) :: qs(len)   ! saturation specific humidity
    REAL(KIND=r8), INTENT(out) :: gam(len)  ! (l/cp)*(d(qs)/dt)
    !
    !--------------------------Local Variables------------------------------
    !
    LOGICAL lflg        ! true if in temperature transition region
    !
    INTEGER i           ! index for vector calculations
    !
    REAL(KIND=r8) omeps     ! 1. - 0.622
    REAL(KIND=r8) trinv     ! reciprocal of ttrice (transition range)
    REAL(KIND=r8) tc        ! temperature (in degrees C)
    REAL(KIND=r8) weight    ! weight for es transition from water to ice
    REAL(KIND=r8) hltalt    ! appropriately modified hlat for T derivatives
    !
    REAL(KIND=r8) hlatsb    ! hlat weighted in transition region
    REAL(KIND=r8) hlatvp    ! hlat modified for t changes above freezing
    REAL(KIND=r8) tterm     ! account for d(es)/dT in transition region
    REAL(KIND=r8) desdt     ! d(es)/dT
    REAL(KIND=r8) epsqs
    !
    !-----------------------------------------------------------------------
    !
    epsqs = epsilo

    omeps = 1.0_r8 - epsqs
    DO i=1,len
       es(i) = estblf(t(i))
       !
       ! Saturation specific humidity
       !
       qs(i) = epsqs*es(i)/(p(i) - omeps*es(i))
       !
       ! The following check is to avoid the generation of negative
       ! values that can occur in the upper stratosphere and mesosphere
       !
       qs(i) = MIN(1.0_r8,qs(i))
       !
       IF (qs(i) < 0.0_r8) THEN
          qs(i) = 1.0_r8
          es(i) = p(i)
       END IF
    END DO
    !
    ! "generalized" analytic expression for t derivative of es
    ! accurate to within 1 percent for 173.16 < t < 373.16
    !
    trinv = 0.0_r8
    IF ((.NOT. icephs) .OR. (ttrice.EQ.0.0_r8)) go to 10
    trinv = 1.0_r8/ttrice
    DO i=1,len
       !
       ! Weighting of hlat accounts for transition from water to ice
       ! polynomial expression approximates difference between es over
       ! water and es over ice from 0 to -ttrice (C) (min of ttrice is
       ! -40): required for accurate estimate of es derivative in transition
       ! range from ice to water also accounting for change of hlatv with t
       ! above freezing where const slope is given by -2369 j/(kg c) = cpv - cw
       !
       tc     = t(i) - tmelt
       lflg   = (tc >= -ttrice .AND. tc < 0.0_r8)
       weight = MIN(-tc*trinv,1.0_r8)
       hlatsb = hlatv + weight*hlatf
       hlatvp = hlatv - 2369.0_r8*tc
       IF (t(i) < tmelt) THEN
          hltalt = hlatsb
       ELSE
          hltalt = hlatvp
       END IF
       IF (lflg) THEN
          tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3) + tc*(pcf(4) + tc*pcf(5))))
       ELSE
          tterm = 0.0_r8
       END IF
       desdt  = hltalt*es(i)/(rgasv*t(i)*t(i)) + tterm*trinv
       gam(i) = hltalt*qs(i)*p(i)*desdt/(cp*es(i)*(p(i) - omeps*es(i)))
       IF (qs(i) == 1.0_r8) gam(i) = 0.0_r8
    END DO
    RETURN
    !
    ! No icephs or water to ice transition
    !
10  DO i=1,len
       !
       ! Account for change of hlatv with t above freezing where
       ! constant slope is given by -2369 j/(kg c) = cpv - cw
       !
       hlatvp = hlatv - 2369.0_r8*(t(i)-tmelt)
       IF (icephs) THEN
          hlatsb = hlatv + hlatf
       ELSE
          hlatsb = hlatv
       END IF
       IF (t(i) < tmelt) THEN
          hltalt = hlatsb
       ELSE
          hltalt = hlatvp
       END IF
       desdt  = hltalt*es(i)/(rgasv*t(i)*t(i))
       gam(i) = hltalt*qs(i)*p(i)*desdt/(cp*es(i)*(p(i) - omeps*es(i)))
       IF (qs(i) == 1.0_r8) gam(i) = 0.0_r8
    END DO
    !
    RETURN
    !
  END SUBROUTINE vqsatd

 SUBROUTINE gffgch(t       ,es      ,itype   )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Computes saturation vapor pressure over water and/or over ice using
    ! Goff & Gratch (1946) relationships. 
    ! <Say what the routine does> 
    ! 
    ! Method: 
    ! T (temperature), and itype are input parameters, while es (saturation
    ! vapor pressure) is an output parameter.  The input parameter itype
    ! serves two purposes: a value of zero indicates that saturation vapor
    ! pressures over water are to be returned (regardless of temperature),
    ! while a value of one indicates that saturation vapor pressures over
    ! ice should be returned when t is less than freezing degrees.  If itype
    ! is negative, its absolute value is interpreted to define a temperature
    ! transition region below freezing in which the returned
    ! saturation vapor pressure is a weighted average of the respective ice
    ! and water value.  That is, in the temperature range 0 => -itype
    ! degrees c, the saturation vapor pressures are assumed to be a weighted
    ! average of the vapor pressure over supercooled water and ice (all
    ! water at 0 c; all ice at -itype c).  Maximum transition range => 40 c
    ! 
    ! Author: J. Hack
    ! 
    !-----------------------------------------------------------------------
    !   use shr_kind_mod, only: r8 => shr_kind_r8
    !   use physconst, only: tmelt
    !   use abortutils, only: endrun

    IMPLICIT NONE
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    REAL(KIND=r8), INTENT(in) :: t          ! Temperature
    !
    ! Output arguments
    !
    INTEGER, INTENT(inout) :: itype   ! Flag for ice phase and associated transition

    REAL(KIND=r8), INTENT(out) :: es         ! Saturation vapor pressure
    !
    !---------------------------Local variables-----------------------------
    !
    REAL(KIND=r8) e1         ! Intermediate scratch variable for es over water
    REAL(KIND=r8) e2         ! Intermediate scratch variable for es over water
    REAL(KIND=r8) eswtr      ! Saturation vapor pressure over water
    REAL(KIND=r8) f          ! Intermediate scratch variable for es over water
    REAL(KIND=r8) f1         ! Intermediate scratch variable for es over water
    REAL(KIND=r8) f2         ! Intermediate scratch variable for es over water
    REAL(KIND=r8) f3         ! Intermediate scratch variable for es over water
    REAL(KIND=r8) f4         ! Intermediate scratch variable for es over water
    REAL(KIND=r8) f5         ! Intermediate scratch variable for es over water
    REAL(KIND=r8) ps         ! Reference pressure (mb)
    REAL(KIND=r8) t0         ! Reference temperature (freezing point of water)
    REAL(KIND=r8) term1      ! Intermediate scratch variable for es over ice
    REAL(KIND=r8) term2      ! Intermediate scratch variable for es over ice
    REAL(KIND=r8) term3      ! Intermediate scratch variable for es over ice
    REAL(KIND=r8) tr         ! Transition range for es over water to es over ice
    REAL(KIND=r8) ts         ! Reference temperature (boiling point of water)
    REAL(KIND=r8) weight     ! Intermediate scratch variable for es transition
    INTEGER itypo   ! Intermediate scratch variable for holding itype
    !
    !-----------------------------------------------------------------------
    !
    ! Check on whether there is to be a transition region for es
    !
    IF (itype < 0) THEN
       tr    = ABS(real(itype,kind=r8))
       itypo = itype
       itype = 1
    ELSE
       tr    = 0.0_r8
       itypo = itype
    END IF
    IF (tr > 40.0_r8) THEN
       WRITE(6,900) tr
       STOP 'GFFGCH'                ! Abnormal termination
    END IF
    !
    IF(t < (tmelt - tr) .AND. itype == 1) go to 10
    !
    ! Water
    !
    ps = 1013.246_r8
    ts = 373.16_r8
    e1 = 11.344_r8*(1.0_r8 - t/ts)
    e2 = -3.49149_r8*(ts/t - 1.0_r8)
    f1 = -7.90298_r8*(ts/t - 1.0_r8)
    f2 = 5.02808_r8*LOG10(ts/t)
    f3 = -1.3816_r8*(10.0_r8**e1 - 1.0_r8)/10000000.0_r8
    f4 = 8.1328_r8*(10.0_r8**e2 - 1.0_r8)/1000.0_r8
    f5 = LOG10(ps)
    f  = f1 + f2 + f3 + f4 + f5
    es = (10.0_r8**f)*100.0_r8
    eswtr = es
    !
    IF(t >= tmelt .OR. itype == 0) go to 20
    !
    ! Ice
    !
10  CONTINUE
    t0    = tmelt
    term1 = 2.01889049_r8/(t0/t)
    term2 = 3.56654_r8*LOG(t0/t)
    term3 = 20.947031_r8*(t0/t)
    es    = 575.185606e10_r8*EXP(-(term1 + term2 + term3))
    !
    IF (t < (tmelt - tr)) go to 20
    !
    ! Weighted transition between water and ice
    !
    weight = MIN((tmelt - t)/tr,1.0_r8)
    es = weight*es + (1.0_r8 - weight)*eswtr
    !
20  CONTINUE
    itype = itypo
    RETURN
    !
900 FORMAT('GFFGCH: FATAL ERROR ******************************',/, &
         'TRANSITION RANGE FOR WATER TO ICE SATURATION VAPOR', &
         ' PRESSURE, TR, EXCEEDS MAXIMUM ALLOWABLE VALUE OF', &
         ' 40.0 DEGREES C',/, ' TR = ',f7.2)
    !
  END SUBROUTINE gffgch
END MODULE wv_saturation

!
!  $Author: pkubota $
!  $Date: 2008/04/09 12:42:57 $
!  $Revision: 1.11 $
!
MODULE Convection

  !   InitConvection
  !
  !   cumulus_driver|--qnegat
  !                 |--Cu_Ara
  !                 |--Cu_Kuolcl
  !                 |--Cu_Grellens
  !                 |--Shall_Tied
  !                 |--Shall_Souza
  !                 |--lrgscl
  !                 


  USE Constants, ONLY :  &
       delq             ,&
       r8,i8,qmin,grav,undef
       

  USE Diagnostics, ONLY:   &
       dodia             , &
       updia             , &
       StartStorDiag     , &
       nDiag_toprec      , & ! total precipiation
       nDiag_cvprec      , & ! convective precipitation
       nDiag_lsprec      , & ! large scale precipitation
       nDiag_snowfl      , & ! snowfall
       nDiag_clheat      , & ! convective latent heating
       nDiag_cmchan      , & ! convective moisture change
       nDiag_lslhea      , & ! large scale latent heating
       nDiag_lsmcha      , & ! large scale moisture change
       nDiag_sclhea      , & ! shallow convective latent heating
       nDiag_scmcha      , & ! shallow convective moisture change
       nDiag_nshcrm      , & ! negative specific humidity correction moisture source
       nDiag_qlicld      , & ! liquid water content in cloud after rainout
       nDiag_trcliq      , & ! Water Liquid Cloud kg/kg
       nDiag_trcice      , & ! Water Ice Cloud kg/kg
       nDiag_cape2d      , & ! CONVECTIVE AVAIL. POT.ENERGY M2/S2 
       nDiag_cine2d      , & ! CONVECTIVE INHIB. ENERGY M2/S2
       nDiag_sweath      , & ! SEVERE WEATHER THREAT
       nDiag_covtcl      

  USE GridHistory, ONLY:   &
       IsGridHistoryOn   , &
       StoreGridHistory  , &
       dogrh             , &
       nGHis_cvprec     , &
       nGHis_clheat     , &
       nGHis_cvmosr     , &
       nGHis_sclhea     , &
       nGHis_shcvmo     , &
       nGHis_toprec     , &
       nGHis_snowfl     , &
       nGHis_sslaht     , &
       nGHis_spstms

  USE Options, ONLY :       &
       rccmbl            , &
       mlrg              , &
       iccon             , &
       ilcon             , &
       iscon             , &
       doprec            , &
       cflric            , &
       dt                , &
       kt                , & 
       ktp               , &  
       ktm               , & 
       jdt               , & 
       nfcnv0            , & 
       isimp             , & 
       nfctrl            , & 
       nfprt             , &
       microphys

  USE wv_saturation, ONLY :       &
      gestbl,&
      findsp

 USE FieldsPhysics, ONLY:  &
       convc             , &
       convt             , &
       convb             , &
       prcp1             , &
       prcp2             , &
       prcp3             , &
       prcpt             , &
       toplv             , &
       botlv             , &
       rVisDiff          , &
       rVisBeam          , &
       rNirDiff          , &
       rNirBeam          , &
       rVisDiffC         , &
       rVisBeamC         , &
       rNirDiffC         , &
       rNirBeamC         , &
       rSwToaDown  

  USE Init, ONLY :       &
       nls

  USE Cu_Kuolcl, ONLY :       &
      InitCu_Kuolcl         , &
      kuolcl

  USE Shall_Tied, ONLY:       &
      InitShall_Tied        , &
      shalv2

  USE MicroPhysics, ONLY:       &
      InitMicroPhysics ,stratiform_tend   
      
  USE Cu_ZhangMcFarlane, ONLY:       &
      Init_Cu_ZhangMcFarlane ,convect_deep_tend   
      
   USE Shall_JHack, ONLY:       &
      convect_shallow_init ,convect_shallow_tend   

  USE Parallelism, ONLY: &
       MsgOne, FatalError

       
  USE Cu_Grellens, ONLY: grellens,InitGrellens
  USE Shall_Souza, ONLY: shallsouza
  USE ModConRas, ONLY: arprep,shllcl,InitArakawa
  USE LrgSclPrec, ONLY: InitLrgScl,lrgscl 
  USE PhysicalFunctions, ONLY: calc_cape,SWEAT_index

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitConvection
  PUBLIC :: cumulus_driver
  PUBLIC :: InitCheckFileConvec
  PUBLIC :: ReStartConvec
  INTEGER, PARAMETER :: ppcnst=3
  REAL(KIND=r8), ALLOCATABLE  :: ql2(:,:,:)
  REAL(KIND=r8), ALLOCATABLE  :: qi2(:,:,:)
  REAL(KIND=r8), ALLOCATABLE  :: ql3(:,:,:)
  REAL(KIND=r8), ALLOCATABLE  :: qi3(:,:,:)
  
CONTAINS

  SUBROUTINE InitConvection(si,del        ,sl         ,cl         , &
                          kMax    ,iMax,jMax,ibMax,jbMax,&
                          ifdy       ,todcld     ,ids        , &
                          idc     ,ifday      ,tod                     )
      
    INTEGER, INTENT(IN) :: kMax
    INTEGER, INTENT(IN) :: iMax
    INTEGER, INTENT(IN) :: jMax
    INTEGER, INTENT(IN) :: ibMax
    INTEGER, INTENT(IN) :: jbMax
    REAL(KIND=r8),    INTENT(IN) :: si (kMax+1)
    REAL(KIND=r8),    INTENT(IN) :: del(kMax  )
    REAL(KIND=r8),    INTENT(IN) :: sl (kMax  )
    REAL(KIND=r8),    INTENT(IN) :: cl (kMax  )     
    INTEGER         , INTENT(OUT  ) :: ifdy
    REAL(KIND=r8)   , INTENT(OUT  ) :: todcld
    INTEGER         , INTENT(OUT  ) :: ids(4)
    INTEGER         , INTENT(OUT  ) :: idc(4)
    INTEGER         , INTENT(IN   ) :: ifday
    REAL(KIND=r8)   , INTENT(IN   ) :: tod
    CALL gestbl()
    IF(TRIM(iccon).EQ.'ARA')  CALL InitArakawa   ()
    IF(TRIM(iccon).EQ.'KUO')  CALL InitCu_Kuolcl () 
    IF(TRIM(iccon).EQ.'GRE')  CALL InitGrellens()
    IF(TRIM(iccon).EQ.'ZMC')  CALL Init_Cu_ZhangMcFarlane(dt,si,sl,del,ibMax,kMax,jbMax,jMax)
    IF(TRIM(ISCON).EQ.'TIED') CALL InitShall_Tied(si, del, sl, cl, kMax)
    IF(TRIM(ISCON).EQ.'JHK'.or. TRIM(ISCON).EQ.'UW')  CALL convect_shallow_init(ISCON,kMax,jMax,ibMax,jbMax,ppcnst,si,sl,del)
    IF(TRIM(ILCON).EQ.'YES'.or.TRIM(ILCON).EQ.'LSC') CALL InitLrgScl()    
    IF(TRIM(ILCON).EQ.'YES'.or.TRIM(ILCON).EQ.'LSC') CALL InitLrgScl()    
    IF(TRIM(ILCON).EQ.'MIC') CALL InitMicroPhysics(kMax,jMax,ibMax,jbMax,ppcnst,si,sl,del)

    ALLOCATE(ql3(ibMax,kMax,jbMax));ql3=0.00001e-12_r8
    ALLOCATE(qi3(ibMax,kMax,jbMax));qi3=0.00001e-12_r8
    ALLOCATE(ql2(ibMax,kMax,jbMax));ql2=0.00001e-12_r8
    ALLOCATE(qi2(ibMax,kMax,jbMax));qi2=0.00001e-12_r8

    IF(TRIM(isimp).NE.'YES') THEN
       CALL InitBoundCondConvec(&
           ifdy,todcld,ids,idc,ifday, &
           tod)
    END IF   
  END SUBROUTINE InitConvection

  SUBROUTINE InitBoundCondConvec(&
       ifdy,todcld,ids,idc,ifday, &
       tod)

    INTEGER         , INTENT(OUT  ) :: ifdy
    REAL(KIND=r8)   , INTENT(OUT  ) :: todcld
    INTEGER         , INTENT(OUT  ) :: ids(4)
    INTEGER         , INTENT(OUT  ) :: idc(4)
    INTEGER         , INTENT(IN   ) :: ifday
    REAL(KIND=r8)   , INTENT(IN   ) :: tod
    CHARACTER(LEN=*), PARAMETER :: h='**(InitBoundCondConvec)**'

    IF(nfcnv0.NE.0) THEN
       CALL MsgOne(h,'Reading previous physics state for restart')

       READ(UNIT=nfcnv0) ifdy,todcld,ids,idc
       READ(UNIT=nfcnv0) convc,convt,convb,prcp1,prcp2,prcp3, &
        prcpt,toplv,botlv
       IF(ifday.GT.0.OR.tod.GT.0.0_r8)READ(UNIT=nfcnv0)rVisDiff,rVisBeam,rNirDiff, &
        rNirBeam,rVisDiffC,rVisBeamC,rNirDiffC,rNirBeamC,rSwToaDown

       REWIND nfcnv0

       IF(nfctrl(4) .GE. 1)WRITE(UNIT=nfprt,FMT=555)ifdy,todcld,ids,idc
    ELSE
       CALL MsgOne(h,'Initializing prec/cloud variables')

       convc=0.0_r8
       convt=0.0_r8
       convb=0.0_r8
       prcp1=0.0_r8
       prcp2=0.0_r8
       prcp3=0.0_r8
       prcpt=0.0_r8
       toplv=0.0_r8
       botlv=0.0_r8
    END IF


555 FORMAT(' CLOUD PROGNOSTIC DATA READ IN. AT FORECAST DAY', &
         I8,' TOD ',F8.1/' STARTING',3I3,I5,' CURRENT',3I3,I5)
  END SUBROUTINE InitBoundCondConvec



  SUBROUTINE InitCheckFileConvec(ifdy  ,todcld,ids   ,idc    )
    INTEGER      , INTENT(OUT  ) :: ifdy
    REAL(KIND=r8), INTENT(OUT  ) :: todcld
    INTEGER      , INTENT(OUT  ) :: ids   (4)
    INTEGER      , INTENT(OUT  ) :: idc   (4)
    CHARACTER(LEN=*), PARAMETER :: h="**(InitCheckFileConvec)**"
    
    !
    !     read cloud dataset for cold start
    !
    IF(nfcnv0.NE.0) THEN
       CALL MsgOne(h,'Read prec/cloud variables')    
       READ(UNIT=nfcnv0) ifdy,todcld,ids,idc
       READ(UNIT=nfcnv0) convc,convt,convb,prcp1,prcp2,prcp3, &
            prcpt,toplv,botlv


       REWIND nfcnv0

       IF(nfctrl(4) .GE. 1) WRITE(UNIT=nfprt,FMT=555)ifdy,todcld,ids,idc

    ELSE
       CALL MsgOne(h,'Initializing prec/cloud variables')
       convc=0.0_r8
       convt=0.0_r8
       convb=0.0_r8
       prcp1=0.0_r8
       prcp2=0.0_r8
       prcp3=0.0_r8
       prcpt=0.0_r8
       toplv=0.0_r8
       botlv=0.0_r8
    END IF
555 FORMAT(' CLOUD PROGNOSTIC DATA READ IN. AT FORECAST DAY', &
         I8,' TOD ',F8.1/' STARTING',3I3,I5,' CURRENT',3I3,I5)

  END SUBROUTINE InitCheckFileConvec
  
   SUBROUTINE ReStartConvec (ifday,tod,idate ,idatec,nfcnv1)

    INTEGER           ,INTENT(IN   ) :: ifday
    REAL(KIND=r8)     ,INTENT(IN   ) :: tod
    INTEGER           ,INTENT(IN   ) :: idate(4)
    INTEGER           ,INTENT(IN   ) :: idatec(4)
    INTEGER           ,INTENT(IN   ) :: nfcnv1

    IF(TRIM(isimp).NE.'YES') THEN
       CALL MsgOne('**(ReStartConvec)**','Saving physics state for restart')

       !$OMP SINGLE
       WRITE(UNIT=nfcnv1) ifday,tod,idate,idatec
       WRITE(UNIT=nfcnv1) convc,convt,convb,prcp1,prcp2,prcp3, &
            prcpt,toplv,botlv
       WRITE(UNIT=nfcnv1) rVisDiff,rVisBeam,rNirDiff,rNirBeam, &
            rVisDiffC,rVisBeamC,rNirDiffC,rNirBeamC,rSwToaDown
       !$OMP END SINGLE
    END IF

  END SUBROUTINE ReStartConvec 
  
  
  SUBROUTINE cumulus_driver (&
            iMax     ,kMax     ,ta       ,tb       ,tc       ,qa       , &
            qb       ,qc       ,ub       ,vb       ,omgb     ,psb      , &
            psb2     ,del      ,sl       ,si       ,zs       ,sens     , &
            evap     ,mask     ,latco    ,fac2x    ,convc    ,convt    , &
            convb    ,prcp1    ,prcp2    ,prcp3    ,prcpt    ,toplv    , &
            botlv    ,convts   ,convcs   ,convbs   ,fac2     ,fac      , &
            geshem   ,ppli     ,ppci     ,prct     ,prcc     ,snowfl   , &
            qliq     ,tsfc     ,pblh     ,tpert    ,qpert    ,tke      , &
            concld   ,cld      ,gliqp    ,gicep    ,gliqm    , &
            gicem)
    !************************************************************************
    !   The cumulus_driver subroutine calls deep and shallow cumulus
    !   parameterization schemes.
    !   more information nilo@cptec.inpe.br
    !   NOTE: This version is not official. You can use only for test.
    !************************************************************************
    !
    !  Definition/
    !---------------
    !             I N P U T  O U T P U T  F O R   G C M
    !             -------------------------------------
    ! INPUT
    !
    !** integer
    !    iMax                   ! end index for longitude domain
    !    kMax                   ! end index for u,v,t,p sigma levels
    !    jdt                   ! number of time step
    !    iccon                  ! cu schemes ex. KUO, ARA, GRE ..
    !   kuo                      ! convection yes(1) or not(0) for shallow convection
    !
    !** real
    !    dt                  ! time step (s)
    !    ta                     ! temperature (K) at time t-1
    !    tb                     ! temperature (K) at time t
    !    tc                     ! temperature (K) at time t+1
    !    qa                     ! water vapor mixing ratio (kg/kg) at time t-1
    !    qb                     ! water vapor mixing ratio (kg/kg) at time t
    !    qc                     ! water vapor mixing ratio (kg/kg) at time t+1
    !    psb                    ! surface pressure (cb)     at time t
    !    ub                     ! u-velocity (m/s) at time t
    !    vb                     ! v-velocity (m/s) at time t
    !    omgb                   ! vertical omega velocity (Pa/s) at time t
    !                           ! it is in half levels along with U,V,T,Q
    !    sl                     ! half sigma layers
    !    si                     ! full sigma layers
    !    del                    ! difference between full sigma layers
    !    xland                  ! land-sea mask (1 for land; 0 for water)
    !    zs                     ! topography (m)
    !    DX                     ! horizontal space interval (m)
    !    qrem,cldm              ! local variables for  RAS-Scheme
    !
    !    hrem,qrem              ! these arrays are needed for the heating 
    !                           ! and mostening from ras  scheme
    !
    !
    !    ktops, kbots           ! these arrays are needed for the new 
    !                           ! shallow convection scheme
    !    cldm                   ! needed for cloud fraction based on mass 
    !                           ! flux
    !    noshal1, kctop1, kcbot1! needed for cloud fraction based on mass 
    !                           ! flux new arrays needed for shallow 
    !                           ! convective clouds
    !     
    !
    !
    !   OUTPUT
    !**  integer
    !    kuo                    ! indix for shalow convection KUO,RAS,KUOG, GRELL
    !    ktop                   ! level of convective cloud top
    !    kbot                   ! level of convective cloud base
    !    plcl                   ! pressure at convection levl for shallow convection
    !                           ! in Kuo 
    !
    !** real
    !   RAINCV                  ! cumulus scheme precipitation (mm)
    !   tc                      ! new temperature (K) at time t+1  after CU precip
    !   qc                      ! new  water vapor mixing ratio (kg/kg) at time t+1.
    !
    !
    !*********************************************************************************
    IMPLICIT NONE
    !              I N P U T     O U T P U T    V A R I A B L E S
    !              ----------------------------------------------
    !              Xa at t-1   Xb at t    Xc at t+1

    ! Dimensions
    INTEGER, INTENT(IN) :: iMax
    INTEGER, INTENT(IN) :: kMax

    ! Sizes
    REAL(KIND=r8), INTENT(IN) :: sl   (kMax)
    REAL(KIND=r8), INTENT(IN) :: del  (kMax)
    REAL(KIND=r8), INTENT(IN) :: si   (kMax+1)

    ! Fixed fields: latitudes, mask and topography
    INTEGER, INTENT(IN) :: latco
    INTEGER(KIND=i8), INTENT(IN) :: mask (iMax) 
    REAL(KIND=r8),    INTENT(IN) :: zs   (iMax)

    ! Temperature (K) and specific humidity (kg/kg) 
    ! at times (a) = T-1, (b) = T and (c) = T+1
    REAL(KIND=r8), INTENT(INOUT) :: ta (iMax,kMax)
    REAL(KIND=r8), INTENT(INOUT) :: tb (iMax,kMax)
    REAL(KIND=r8), INTENT(INOUT) :: tc (iMax,kMax)
    REAL(KIND=r8), INTENT(INOUT) :: qa (iMax,kMax)
    REAL(KIND=r8), INTENT(INOUT) :: qb (iMax,kMax)
    REAL(KIND=r8), INTENT(INOUT) :: qc (iMax,kMax)
    REAL(KIND=r8), INTENT(OUT  ) :: qliq(iMax,kMax) !qrc liquid water content in cloud after rainout

    ! Wind at time T
    ! in half levels along with U,V,T,Q
    REAL(KIND=r8), INTENT(IN) :: ub   (iMax,kMax) ! (m/s) 
    REAL(KIND=r8), INTENT(IN) :: vb   (iMax,kMax) ! (m/s)
    REAL(KIND=r8), INTENT(IN) :: omgb (iMax,kMax) ! (Pa/s)

    ! Surface pressure (cb) at time T
    REAL(KIND=r8), INTENT(IN) :: psb  (iMax)
    REAL(KIND=r8), INTENT(IN) :: psb2 (iMax)
    REAL(KIND=r8), INTENT(IN) :: tsfc (iMax)
    REAL(KIND=r8), INTENT(IN) :: pblh (iMax)
    REAL(KIND=r8), INTENT(IN) :: tpert(iMax)
    REAL(KIND=r8), INTENT(IN) :: qpert(iMax)
    REAL(KIND=r8), INTENT(IN) :: tke  (iMax,kMax)

    ! Heat/Water sfc fluxes
    REAL(KIND=r8), INTENT(IN) :: sens (iMax)
    REAL(KIND=r8), INTENT(IN) :: evap (iMax)

    ! UNCLASSIFIED VARIABLES
    REAL(KIND=r8), INTENT(IN)      :: fac2x  
    REAL(KIND=r8), INTENT(IN)      :: fac
    REAL(KIND=r8), INTENT(IN)      :: fac2

    REAL(KIND=r8), INTENT(INOUT) :: convc  (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: convt  (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: convb  (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: prcp1  (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: prcp2  (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: prcp3  (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: prcpt  (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: toplv  (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: botlv  (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: convts (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: convcs (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: convbs (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: geshem (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: ppli   (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: ppci   (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: prcc   (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: snowfl (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: prct   (iMax)
    REAL(KIND=r8), INTENT(INOUT) :: concld    (iMax,kMax)
    REAL(KIND=r8), INTENT(INOUT) :: cld       (iMax,kMax)
!    REAL(KIND=r8), INTENT(INOUT) :: tracerLiq (iMax,kMax)
!    REAL(KIND=r8), INTENT(INOUT) :: tracerIce (iMax,kMax)
!    REAL(KIND=r8), INTENT(INOUT) :: tracerLiqm (iMax,kMax)
!    REAL(KIND=r8), INTENT(INOUT) :: tracerIcem (iMax,kMax)
    REAL(KIND=r8),    OPTIONAL,   INTENT(INOUT) :: gicem (iMax,kmax)
    REAL(KIND=r8),    OPTIONAL,   INTENT(INOUT) :: gicep (iMax,kmax)
    REAL(KIND=r8),    OPTIONAL,   INTENT(INOUT) :: gliqm (iMax,kmax)
    REAL(KIND=r8),    OPTIONAL,   INTENT(INOUT) :: gliqp (iMax,kmax)

    !
    !               L O C A L    V A R I A B L E S
    !              -------------------------------
    INTEGER       :: i
    INTEGER       :: k
    INTEGER       :: mask2   (iMax)
    REAL(KIND=r8) :: ps1     (iMax)
    REAL(KIND=r8) :: ps2     (iMax)
    REAL(KIND=r8) :: PS_work (iMax)
    REAL(KIND=r8) :: terr    (iMax)  

    ! (T,q) before any convection
    REAL(KIND=r8) :: tBegin (iMax,kMax)
    REAL(KIND=r8) :: qBegin (iMax,kMax)

    ! (T,q) after deep convection
    REAL(KIND=r8) :: qDeep(iMax,kMax)
    REAL(KIND=r8) :: tDeep(iMax,kMax)
    ! (T,q) after shallow convection
    REAL(KIND=r8) :: tShal(iMax,kMax)
    REAL(KIND=r8) :: qShal(iMax,kMax)
    ! (T,q) after large scale adjustment
    REAL(KIND=r8) :: tLrgs(iMax,kMax)
    REAL(KIND=r8) :: qLrgs(iMax,kMax)

    ! Shallow heat/moist
    REAL(KIND=r8) :: sclhea(iMax,kMax)
    REAL(KIND=r8) :: scmcha(iMax,kMax)
    ! Convective heat/moist
    REAL(KIND=r8) :: clheat(iMax,kMax)
    REAL(KIND=r8) :: cmchan(iMax,kMax)
    ! Large scale heat/moist
    REAL(KIND=r8) :: lslhea(iMax,kMax)
    REAL(KIND=r8) :: lsmcha(iMax,kMax)

    ! Working copies of input fields (T,q) at times (a,b,c)
    REAL(KIND=r8) :: q1(iMax,kMax)
    REAL(KIND=r8) :: q2(iMax,kMax)
    REAL(KIND=r8) :: q3(iMax,kMax)


    REAL(KIND=r8) :: t1(iMax,kMax)
    REAL(KIND=r8) :: t2(iMax,kMax)
    REAL(KIND=r8) :: t3(iMax,kMax)

    ! Wind components for grell ensemble
    REAL(KIND=r8) :: u2(iMax,kMax)
    REAL(KIND=r8) :: v2(iMax,kMax)
    REAL(KIND=r8) :: w2(iMax,kMax)

    REAL(KIND=r8) :: icefrac  (iMax)
    REAL(KIND=r8) :: landfrac (iMax)
    REAL(KIND=r8) :: ocnfrac  (iMax)
    REAL(KIND=r8) :: landm      (iMax)     ! land fraction ramped over water
    REAL(KIND=r8) :: snowh      (iMax)     ! Snow depth over land, water equivalent (m)
    REAL(KIND=r8) :: dlf          (iMax,kMax)    ! detrained water from ZM
    REAL(KIND=r8) :: rliq      (iMax)      ! vertical integral of liquid not yet in q(ixcldliq)
    real(KIND=r8) :: rliq2(iMax)                   ! vertical integral of liquid from shallow scheme
    REAL(KIND=r8) :: cmfmc  (iMax,kMax+1)   ! convective mass flux--m sub c
    REAL(KIND=r8) :: cmfmc2 (iMax,kMax+1)   ! shallow convective mass flux--m sub c
    REAL(KIND=r8) :: ts      (iMax)      ! surface temperature
    REAL(KIND=r8) :: sst      (iMax)      ! sea surface temperature
    REAL(KIND=r8) :: zdu          (iMax,kMax)   ! detrainment rate from deep convection
    REAL(KIND=r8) :: prec_str     (iMax)     ! [Total] sfc flux of precip from stratiform (m/s) 
    REAL(KIND=r8) :: snow_str     (iMax)  ! [Total] sfc flux of snow from stratiform   (m/s)
    REAL(KIND=r8) :: prec_sed     (iMax)  ! surface flux of total cloud water from sedimentation
    REAL(KIND=r8) :: snow_sed     (iMax)  ! surface flux of cloud ice from sedimentation
    REAL(KIND=r8) :: prec_pcw     (iMax)  ! sfc flux of precip from microphysics(m/s)
    REAL(KIND=r8) :: snow_pcw     (iMax)  ! sfc flux of snow from microphysics (m/s)

    !
    ! UNCLASSIFIED VARIABLES
    !
    REAL(KIND=r8) :: fdqn (iMax,kMax)
    REAL(KIND=r8) :: RAINCV     (iMax)
    REAL(KIND=r8) :: SNOWCV     (iMax)
    REAL(KIND=r8) :: Total_Rain (iMax)
    REAL(KIND=r8) :: Total_Snow (iMax)

    !*******************************************
    !               Ktopos nao usado fora
    !            kctop1  usado para ARA fora
    !*******************************************
    REAL(KIND=r8)    :: hrem  (iMax,kMax)
    REAL(KIND=r8)    :: qrem  (iMax,kMax)
    REAL(KIND=r8)    :: cldm  (iMax)
    INTEGER          :: kctop1(iMax)
    INTEGER          :: kcbot1(iMax)
    INTEGER          :: noshal(iMax)
    !**********
    !others
    !*********
    INTEGER          :: ktop (iMax)
    INTEGER          :: kuo  (iMax)
    INTEGER          :: ktops(iMax)
    REAL(KIND=r8)    :: plcl (iMax)
    INTEGER          :: kbot (iMax)
    INTEGER          :: kbots(iMax) 
    REAL(KIND=r8)    :: dq   (iMax,kMax)
    REAL(KIND=r8)    :: rdt
    LOGICAL          :: newr
    LOGICAL          :: ghl_local
    REAL(KIND=r8)    :: snowflg(iMax)   
    REAL(KIND=r8)    :: prec_zmc(iMax)                ! total precipitation from ZM convection
    REAL(KIND=r8)    :: snow_zmc(iMax)                ! snow from ZM convection 
    real(KIND=r8)    :: prec_cmf(iMax)                ! total precipitation from Hack convection
    real(KIND=r8)    :: snow_cmf(iMax)                ! snow from Hack convection
    REAL(KIND=r8)    :: SCRa(iMax,kMax)
    REAL(KIND=r8)    :: SCRb(iMax,kMax)
    REAL(KIND=r8)    :: cape(iMax) 
    REAL(KIND=r8)    :: cin (iMax) 
    REAL(KIND=r8)    :: LCL (iMax) 
    REAL(KIND=r8)    :: LFC (iMax) 
    REAL(KIND=r8)    :: SWEAT(iMax)
    INTEGER          :: i3dflag
    !-----------------------------------------------------------------
    !-----------------------------------------------------------------

    rdt=1.0_r8/dt
    ghl_local = IsGridHistoryOn()

    ! Check for negative values of specific humidity
    ! Convert virtual temperature to thermodinamic temperature
    CALL qnegat (qa, fdqn, ta, (1.0_r8/dt), del, iMax, kMax)! time t-1
    CALL qnegat (qb, fdqn, tb, (1.0_r8/dt), del, iMax, kMax)! time t
    CALL qnegat (qc, fdqn, tc, (1.0_r8/dt), del, iMax, kMax)! time t+1
    ! Initialize cloud variables with unrealistic values
    DO i=1,iMax
       ktop (i) = -1000
       ktops(i) = -1000
       kuo  (i) = -1000
       plcl (i) = -1.0e3_r8
       rliq (i) = 0.0_r8
       rliq2 (i) = 0.0_r8
    END DO   
    cmfmc  = 0.0_r8! convective mass flux--m sub c
    cmfmc2 = 0.0_r8! shallow convective mass flux--m sub c
    prec_cmf= 0.0_r8
    snow_cmf= 0.0_r8
    !ql3(1:iMax,1:kMax,latco)=0.0_r8
    !qi3(1:iMax,1:kMax,latco)=0.0_r8

    ! Initialize surface variables
    DO i=1,iMax
       !surface pressure
       !!T+1                ps2(i)=psb(i)
       !!T                  ps2(i)=psb2(i)
       ps1(i)    =psb2(i)!T+1 
       ps2(i)    =psb(i) !T   
       PS_work(i)=psb(i)

       terr(i)   =MAX(zs(i),0.0_r8)

       RAINCV(i)     = 0.0_r8
       SNOWCV(i)     = 0.0_r8
       Total_Rain(i) = 0.0_r8
       Total_Snow(i) = 0.0_r8
    END DO   

    ! Copy (T,q) at t-1, t and t+1 to work arrays
    DO i=1,iMax
       DO k=1,kMax
          dlf(i,k)=0.0_r8 
          zdu(i,k)=0.0_r8 
          T1(i,k)=ta(i,k)
          T2(i,k)=tb(i,k)
          T3(i,k)=tc(i,k)

          q1(i,k)=qa(i,k)
          q2(i,k)=qb(i,k)
          q3(i,k)=qc(i,k)
       END DO
    END DO
    IF (microphys) THEN
       DO k=1,kMax
          DO i=1,iMax
             ql3     (i,k,latco) = gliqp(i,k)
             qi3     (i,k,latco) = gicep(i,k)
             ql2     (i,k,latco) = gliqm(i,k)
             qi2     (i,k,latco) = gicem(i,k)
          END DO
       END DO
    ELSE
       DO k=1,kMax
          DO i=1,iMax
             ql3     (i,k,latco) = 0.0_r8 
             qi3     (i,k,latco) = 0.0_r8 
             ql2     (i,k,latco) = 0.0_r8 
             qi2     (i,k,latco) = 0.0_r8 
          END DO
       END DO
    END IF
    DO i=1,iMax
       DO k=1,kMax
          tBegin(i,k)= tc(i,k)! time t+1
          qBegin(i,k)= qc(i,k)! time t+1
       END DO
    END DO
    !-----------------------------------------------------------------
    ! Calcule CAPE and CIN
    !-----------------------------------------------------------------
    i3dflag=0
    SCRa=0.0_r8
    SCRb=0.0_r8
    CALL calc_cape( &
       iMax                 , &!INTEGER      , INTENT(IN   ) :: nCols
       kMax                 , &!INTEGER      , INTENT(IN   ) :: kMax
       si                   , &!REAL(KIND=r8), INTENT(IN   ) :: si   (kMax+1)
       sl                   , &!REAL(KIND=r8), INTENT(IN   ) :: sl   (kMax)
       ps2*10.0_r8          , &!REAL(KIND=r8), INTENT(IN   ) :: PSFC (nCols)   ! psfc is pressure in hPa
       terr                 , &!REAL(KIND=r8), INTENT(IN   ) :: HGT  (nCols)   ! topography m
       t3                   , &!REAL(KIND=r8), INTENT(IN   ) :: TK    (nCols,kMax)     ! TK is temp in K, T is theta-300
       q3                   , &!REAL(KIND=r8), INTENT(IN   ) :: QV    (nCols,kMax)
       SCRa                 , &!REAL(KIND=r8), INTENT(OUT  ) :: SCRa  (nCols,kMax)
       SCRb                 , &!REAL(KIND=r8), INTENT(OUT  ) :: SCRb  (nCols,kMax)
       i3dflag                )!INTEGER      , INTENT(IN   ) :: i3dflag
       DO i=1,iMax
          cape(i) = SCRa(i,1)
          cin (i) = SCRa(i,2)
          LCL (i) = SCRa(i,3)
          LFC (i) = SCRa(i,4)
       END DO
       SWEAT=0.0_r8
       SWEAT=SWEAT_index(t3,ub,vb,iMax,kMax)
 
    !-----------------------------------------------------------------
    ! Deep Convection
    !-----------------------------------------------------------------
    IF(TRIM(iccon).EQ.'ARA') THEN
       CALL arprep( &
            tc    ,qc    ,sl    ,si   ,ps2  ,ktop  ,kbot  ,RAINCV, &
            hrem  ,qrem  ,dt    ,T3   ,q3   ,del   ,kuo   ,cldm  , &
            cflric,kMax+1,kMax-1,iMax ,kMax ,nls)

       CALL shllcl(dt    ,ps2   ,sl   ,si   ,q3    ,T3    ,kuo   , &
            plcl  ,ktops ,kbots ,iMax ,kMax  )
       CALL qnegat2 (q3, fdqn,(1.0_r8/dt), del, iMax, kMax)! time t+1
    END IF
    
    IF(TRIM(iccon).EQ.'KUO')THEN
       CALL  kuolcl(dt, ps2, del, sl, si, q1, q3, &
            T3, dq, RAINCV, kuo, plcl, ktop, &
            kbot, iMax, kMax)
       CALL qnegat2 (q3, fdqn,(1.0_r8/dt), del, iMax, kMax)! time t+1
    END IF

    IF(TRIM(iccon).EQ.'GRE')THEN

       ! grell mask
       DO i=1,iMax      
          IF(mask(i).GT.0_i8)THEN
             mask2(i)=0 ! land
          ELSE
             mask2(i)=1 ! water/ocean
          END IF
       END DO

       ! grell wind
       DO i=1,iMax
          DO k=1,kMax
             u2(i,k)=ub(i,k)
             v2(i,k)=vb(i,k)
             w2(i,k)=omgb(i,k)*1000.0_r8  ! (Pa/s)
          END DO
       END DO

       CALL grellens(ps2, sl,u2,v2,w2,T2,T3, q2,q3,         &
            terr,mask2, dt, RAINCV,kuo,ktop,kbot,plcl,qliq,      &
            iMax,kMax)
       CALL qnegat2 (q3, fdqn,(1.0_r8/dt), del, iMax, kMax)! time t+1
    END IF

    IF(TRIM(iccon).EQ.'ZMC')THEN
       DO k=1,kMax
          DO i=1,iMax
             u2(i,k)=ub(i,k)
             v2(i,k)=vb(i,k)
             w2(i,k)=omgb(i,k)*1000.0_r8  ! (Pa/s)
          END DO
       END DO
       CALL convect_deep_tend( &
       iMax                         , &! INTEGER, INTENT(in) :: pcols   ! number of columns (max)
       kMax+1                       , &! INTEGER, INTENT(in) :: pverp   ! number of vertical levels + 1
       kMax                         , &! INTEGER, INTENT(in) :: pver    ! number of vertical levels
       latco                        , &! INTEGER, INTENT(in) :: latco    ! number of latitudes
       1                            , &! INTEGER, INTENT(in) :: pcnst=1 ! number of advected constituents (including water vapor)
       2                            , &! INTEGER, INTENT(in) :: pnats=2 ! number of non-advected constituents
       ps2     (1:iMax)*1000_r8     , &! REAL(r8), INTENT(in)  :: state_ps  (pcols)    !(pcols) ! surface pressure(Pa)
       terr    (1:iMax)*grav        , &! REAL(r8), INTENT(in)  :: state_phis   (pcols)    !(pcols)  ! surface geopotential
       t3      (1:iMax,1:kMax)      , &! REAL(r8), INTENT(in)  :: state_t   (pcols,pver)!(pcols,pver)! temperature (K)
       q3      (1:iMax,1:kMax)      , &! REAL(r8), INTENT(in)  :: state_qv  (pcols,pver)!(pcols,pver,ppcnst)! vapor  mixing ratio (kg/kg moist or dry air depending on type)
       ql3     (1:iMax,1:kMax,latco), &! REAL(r8), INTENT(in)  :: state_ql (pcols,pver)!(pcols,pver,ppcnst)! liquid mixing ratio (kg/kg moist or dry air depending on type)
       qi3     (1:iMax,1:kMax,latco), &! REAL(r8), INTENT(in)  :: state_qi (pcols,pver)!(pcols,pver,ppcnst)! ice  mixing ratio (kg/kg moist or dry air depending on type)
       w2      (1:iMax,1:kMax)      , &! REAL(r8), INTENT(in)  :: state_omega  (pcols,pver)!(pcols,pver)! vertical pressure velocity (Pa/s) 
       cld     (1:iMax,1:kMax)      , &! REAL(r8), INTENT(out) :: state_cld(ibMax,kMax)          !  cloud fraction
       prec_zmc(1:iMax)             , &! real(r8), intent(out) :: prec(pcols)   ! total precipitation(m/s)
       pblh    (1:iMax)             , &! real(r8), intent(in)  :: pblh(pcols)
       cmfmc   (1:iMax,1:kMax+1)    , &! real(r8), intent(out) :: cmfmc(pcols,pverp)
       tpert   (1:iMax)             , &! real(r8), intent(in)  :: tpert(pcols)
       dlf     (1:iMax,1:kMax)      , &! real(r8), intent(out) :: dlf(pcols,pver)! scattrd version of the detraining cld h2o tend
       zdu     (1:iMax,1:kMax)      , &! real(r8), intent(out) :: zdu(pcols,pver)
       rliq    (1:iMax)             , &! real(r8), intent(out) :: rliq(pcols) ! reserved liquid (not yet in cldliq) for energy integrals
       2*dt                         , &! real(r8), intent(in)  :: ztodt       ! 2 delta t (model time increment)
       snow_zmc(1:iMax)             , &! real(r8), intent(out) :: snow(pcols)   ! snow from ZM convection 
       ktop    (1:iMax)             , &
       kbot    (1:iMax)             , &
       kuo     (1:iMax)               )
       RAINCV=MAX(prec_zmc*dt,0.0_r8)
       SNOWCV=MAX(snow_zmc*dt,0.0_r8)
       CALL qnegat2 (q3, fdqn,(1.0_r8/dt), del, iMax, kMax)! time t+1
    END IF


    ! Save humd/temp after deep convection
    DO k=1,kMax
       DO i=1,iMax
          q3   (i,k)       = MAX(q3  (i,k)      ,1.0e-12_r8)
          IF (microphys) THEN
          ql3  (i,k,latco) = MAX(ql3 (i,k,latco),0.0e-12_r8)
          qi3  (i,k,latco) = MAX(qi3 (i,k,latco),0.0e-12_r8)
          ql2  (i,k,latco) = MAX(ql2 (i,k,latco),0.0e-12_r8)
          qi2  (i,k,latco) = MAX(qi2 (i,k,latco),0.0e-12_r8)
          END IF 
          qDeep(i,k)=q3(i,k)
          tDeep(i,k)=t3(i,k)
       END DO
    END DO

    ! Copy deep convection rain to total rain
    Total_Rain=RAINCV
    Total_Snow=SNOWCV    
    !-----------------------------------------------------------------
    ! Shallow Convection
    !-----------------------------------------------------------------
    IF(TRIM(ISCON).EQ.'TIED')THEN
       IF(TRIM(iccon).EQ.'KUO'.OR.TRIM(iccon).EQ.'GRE' .OR.TRIM(iccon).EQ.'ZMC')THEN
          newr=.FALSE.
          CALL shalv2(si, sl, t3, q3, PS_work, dt, &
               ktop, plcl, kuo, kMax+1, kctop1, kcbot1, noshal, &
               newr, iMax, kMax)
       END IF

       IF(TRIM(iccon).EQ.'ARA') THEN
          newr=.TRUE.
          CALL shalv2(si, sl, t3, q3, PS_work, dt, &
               ktops, plcl, kuo, kMax+1, kctop1, kcbot1, noshal, &
               newr, iMax, kMax)
       END IF
       CALL qnegat2 (q3, fdqn,(1.0_r8/dt), del, iMax, kMax)! time t+1
    END IF

     IF(TRIM(ISCON).EQ.'JHK' .OR. TRIM(ISCON).EQ.'UW' )THEN
       !IF (nscalars <= 0) THEN
       !    dlf   =0.0_r8
       !    rliq  =0.0_r8
       !    cmfmc =0.0_r8
       !END IF
       DO i=1,iMax
          DO k=1,kMax
             u2(i,k)=ub(i,k)
             v2(i,k)=vb(i,k)
             w2(i,k)=omgb(i,k)*1000.0_r8  ! (Pa/s)
          END DO
       END DO
       call convect_shallow_tend ( &
        latco      , &!INTEGER, INTENT(IN   )  :: latco
        iMax       , &!INTEGER, INTENT(IN   )  :: pcols
        kMax       , &!INTEGER, INTENT(IN   )  :: pver
        kMax+1     , &!INTEGER, INTENT(IN   )  :: pverp
        jdt        , &!INTEGER, INTENT(IN   )  :: nstep
        1          , &! INTEGER, INTENT(in):: pcnst=1 ! number of advected constituents (including water vapor)
        2          , &! INTEGER, INTENT(in):: pnats=2 ! number of non-advected constituents
        2*dt       , &!real(r8), intent(in) :: ztodt        ! 2 delta-t (seconds)
        ps2     (1:iMax)*1000_r8     , &! REAL(r8), INTENT(in)  :: state_ps   (pcols)    !(pcols)   ! surface pressure(Pa)      
        terr    (1:iMax)*grav        , &! REAL(r8), INTENT(in)  :: state_phis   (pcols)  !(pcols)   ! surface geopotential
        t3      (1:iMax,1:kMax)      , &! REAL(r8), INTENT(in)  :: state_t   (pcols,pver)!(pcols,pver)! temperature (K)
        q3      (1:iMax,1:kMax)      , &! REAL(r8), INTENT(in)  :: state_qv   (pcols,pver)!(pcols,pver,ppcnst)! vapor  mixing ratio (kg/kg moist or dry air depending on type)
        ql3     (1:iMax,1:kMax,latco), &! REAL(r8), INTENT(in)  :: state_ql    (pcols,pver)!(pcols,pver,ppcnst)! liquid mixing ratio (kg/kg moist or dry air depending on type)
        qi3     (1:iMax,1:kMax,latco), &! REAL(r8), INTENT(in)  :: state_qi    (pcols,pver)!(pcols,pver,ppcnst)! ice    mixing ratio (kg/kg moist or dry air depending on type)
        u2      (1:iMax,1:kMax)      , &
        v2      (1:iMax,1:kMax)      , &
        w2      (1:iMax,1:kMax)      , &!REAL(r8), INTENT(in   )  :: state_omega  (pcols,pver)!(pcols,pver)! vertical pressure velocity (Pa/s) 
        concld  (1:iMax,1:kMax)      , &!REAL(r8), INTENT(in   )  ::
        cld     (1:iMax,1:kMax)      , &!REAL(r8), INTENT(in   )  ::
        qpert   (1:iMax)             , &!real(r8), intent(in   ) :: qpert(pcols)  ! PBL perturbation specific humidity
        pblh    (1:iMax)             , &!real(r8), intent(in   ) :: pblht(pcols)    ! PBL height (provided by PBL routine)
        cmfmc   (1:iMax,1:kMax+1)    , &!real(r8), intent(inout) :: cmfmc(pcols,pverp)  ! moist deep + shallow convection cloud mass flux
        cmfmc2  (1:iMax,1:kMax+1)    , &!real(r8), intent(out  ) :: cmfmc2(pcols,pverp)  ! moist shallow convection cloud mass flux
        prec_cmf(1:iMax)             , &!real(r8), intent(out  ) :: precc(pcols)     ! convective precipitation rate
        dlf     (1:iMax,1:kMax)      , &!real(r8), intent(inout) :: qc(pcols,pver)      ! dq/dt due to export of cloud water  ! detrained water 
        rliq    (1:iMax)             , &!real(r8), intent(inout) :: rliq(pcols) ! vertical integral of liquid not yet in q(ixcldliq)
        rliq2   (1:iMax)             , &!real(r8), intent(out  ) :: rliq2(pcols) 
        snow_cmf(1:iMax)             , &!real(r8), intent(out  ) :: snow(pcols)  ! snow from this convection     
        tke     (1:iMax,1:kMax)      , &
        ktop    (1:iMax)             , &
        KUO     (1:iMax)             , &
        kctop1  (1:iMax)             , &
        kcbot1  (1:iMax)             , &
        noshal  (1:iMax)               )
       prec_cmf=(prec_cmf*dt)
       snow_cmf=(snow_cmf*dt)
       CALL qnegat2 (q3, fdqn,(1.0_r8/dt), del, iMax, kMax)! time t+1
       
        RAINCV = RAINCV + MAX(prec_cmf,0.0_r8)
        SNOWCV = SNOWCV + MAX(snow_cmf,0.0_r8)
    END IF
    
    IF(TRIM(ISCON).EQ.'SOUZ')THEN
       WRITE(*,*)"it is not available yet. Problem in the water balance'"
       STOP "ERROR"
       !CALL Shallsouza(t3,q3,PS_work,sl,sens,evap,dt,iMax,kMax,kuo, &
       !     noshal, kctop1, kcbot1, 560.0_r8, 1.6_r8)
    END IF

    Total_Rain=RAINCV
    Total_Snow=SNOWCV
    ! Save humd/temp after shallow convection
    DO k=1,kMax
      DO i=1,iMax
        q3  (i,k)       = MAX(q3  (i,k)      ,1.0e-12_r8)
        IF (microphys) THEN
        ql3 (i,k,latco) = MAX(ql3 (i,k,latco),1.0e-12_r8)
        qi3 (i,k,latco) = MAX(qi3 (i,k,latco),1.0e-12_r8)
        ql2 (i,k,latco) = MAX(ql2 (i,k,latco),1.0e-12_r8)
        qi2 (i,k,latco) = MAX(qi2 (i,k,latco),1.0e-12_r8)
         END IF
        tShal(i,k)=t3(i,k)
        qShal(i,k)=q3(i,k)
      END DO
    END DO

    !-----------------------------------------------------------------
    ! Large Scale Precipitation
    !-----------------------------------------------------------------
    
    IF(TRIM(ILCON).EQ.'LSC' .OR. TRIM(ILCON).EQ.'YES' ) THEN
      CALL lrgscl(Total_Rain, t3, dq, q3, ps2, del, sl, dt, &
                mlrg, latco, iMax, kMax)
      CALL qnegat2 (q3, fdqn,(1.0_r8/dt), del, iMax, kMax)! time t+1
        snow_str=0.0_r8
        prec_sed=0.0_r8
        snow_sed=0.0_r8
        prec_pcw=0.0_r8
        snow_pcw=0.0_r8
    END IF
    IF(TRIM(ILCON).EQ.'MIC') THEN
       ! grell mask
       snowh          =0.0_r8
       icefrac=0.0_r8
       landfrac=0.0_r8
       ocnfrac=0.0_r8
       DO i=1,iMax
          IF(mask(i) == 13 .or. mask(i) == 15 )snowh(i)=5.0_r8      
          IF(mask(i).GT.0_i8)THEN
             ! land
             icefrac(i)=0.0_r8
             landfrac(i)=1.0_r8
             ocnfrac(i)=0.0_r8
          ELSE
             ! water/ocean
             landfrac(i)=0.0_r8
             ocnfrac(i) =1.0_r8
             IF(ocnfrac(i).GT.0.01_r8.AND.tsfc(i).LT.260.0_r8) THEN
                icefrac(i)=1.0_r8
                ocnfrac(i) =0.0_r8
             ENDIF
          END IF
       END DO
       DO i=1,iMax
          DO k=1,kMax
             u2(i,k)=ub(i,k)
             v2(i,k)=vb(i,k)
             w2(i,k)=omgb(i,k)*1000.0_r8  ! (Pa/s)
          END DO
       END DO
       !CALL qpart(kMax,iMax,t1(1:iMax,1:kMax),ps2(1:iMax),sl(1:kMax),q1(1:iMax,1:kMax),ql3(1:iMax,1:kMax,latco),qi3(1:iMax,1:kMax,latco),'subtraction')              
       !CALL qpart(kMax,iMax,t2(1:iMax,1:kMax),ps2(1:iMax),sl(1:kMax),q2(1:iMax,1:kMax),ql3(1:iMax,1:kMax,latco),qi3(1:iMax,1:kMax,latco),'subtraction')              
       !CALL qpart(kMax,iMax,t3(1:iMax,1:kMax),ps2(1:iMax),sl(1:kMax),q3(1:iMax,1:kMax),ql3(1:iMax,1:kMax,latco),qi3(1:iMax,1:kMax,latco),'subtraction')       
       landm          =0.0_r8  
       ts             =tsfc
       sst            =tsfc
       CALL stratiform_tend(&
       jdt                              , &! INTEGER , INTENT(in) :: ibMax               ! number of columns (max)
       iMax                             , &! INTEGER , INTENT(in) :: ibMax               ! number of columns (max)
       kMax                             , &! INTEGER , INTENT(in) :: kMax               ! number of vertical levels
       kMax+1                           , &! INTEGER , INTENT(in) :: kMax+1               ! number of vertical levels + 1
       ppcnst                           , &! INTEGER , INTENT(in) :: ppcnst               ! number of constituent
       latco                            , &! INTEGER , INTENT(in) :: latco                      ! latitude
       dt                               , &! REAL(r8), INTENT(in)  :: dtime                        ! timestep
       ps2         (1:iMax)*1000_r8     , &! REAL(r8), INTENT(in)  :: state_ps    (ibMax)        !(ibMax)     ! surface pressure(Pa)
       terr        (1:iMax)*grav        , &! REAL(r8), INTENT(in)  :: state_phis  (ibMax)        !(ibMax)     ! surface geopotential
       t2          (1:iMax,1:kMax)      , &! REAL(r8), INTENT(in)  :: state_t     (ibMax,kMax)!(ibMax,kMax)! temperature (K)
       t3          (1:iMax,1:kMax)      , &! REAL(r8), INTENT(in)  :: state_t     (ibMax,kMax)!(ibMax,kMax)! temperature (K)
       q2          (1:iMax,1:kMax)      , &! REAL(r8), INTENT(in)  :: state_qv    (ibMax,kMax)!(ibMax,kMax,ppcnst)! vapor  mixing ratio (kg/kg moist or dry air depending on type)
       q3          (1:iMax,1:kMax)      , &! REAL(r8), INTENT(in)  :: state_qv    (ibMax,kMax)!(ibMax,kMax,ppcnst)! vapor  mixing ratio (kg/kg moist or dry air depending on type)
       ql2         (1:iMax,1:kMax,latco), &! REAL(r8), INTENT(in)  :: state_ql    (pcols,pver)!(pcols,pver,ppcnst)! liquid mixing ratio (kg/kg moist or dry air depending on type)
       ql3         (1:iMax,1:kMax,latco), &! REAL(r8), INTENT(in)  :: state_ql    (pcols,pver)!(pcols,pver,ppcnst)! liquid mixing ratio (kg/kg moist or dry air depending on type)
       qi2         (1:iMax,1:kMax,latco), &! REAL(r8), INTENT(in)  :: state_qi    (pcols,pver)!(pcols,pver,ppcnst)! ice    mixing ratio (kg/kg moist or dry air depending on type)
       qi3         (1:iMax,1:kMax,latco), &! REAL(r8), INTENT(in)  :: state_qi    (pcols,pver)!(pcols,pver,ppcnst)! ice    mixing ratio (kg/kg moist or dry air depending on type)
       w2          (1:iMax,1:kMax)      , &! REAL(r8), INTENT(in)  :: state_omega (ibMax,kMax)!(ibMax,kMax)! vertical pressure velocity (Pa/s) 
       icefrac     (1:iMax)             , &! REAL(r8), INTENT(in)  :: icefrac     (ibMax)            ! sea ice fraction (fraction)
       landfrac    (1:iMax)             , &! REAL(r8), INTENT(in)  :: landfrac    (ibMax)            ! land fraction (fraction)
       ocnfrac     (1:iMax)             , &! REAL(r8), INTENT(in)  :: ocnfrac     (ibMax)            ! ocean fraction (fraction)
       landm       (1:iMax)             , &! REAL(r8), INTENT(in)  :: landm       (ibMax)            ! land fraction ramped over water
       snowh       (1:iMax)             , &! REAL(r8), INTENT(in)  :: snowh       (ibMax)            ! Snow depth over land, water equivalent (m)
       dlf         (1:iMax,1:kMax)      , &! REAL(r8), INTENT(in)  :: state_dlf   (ibMax,kMax)    ! detrained water from ZM
       rliq        (1:iMax)             , &! REAL(r8), INTENT(in)  :: rliq        (ibMax)             ! vertical integral of liquid not yet in q(ixcldliq)
       cmfmc       (1:iMax,1:kMax+1)    , &! REAL(r8), INTENT(in)  :: state_cmfmc (ibMax,kMax+1)   ! convective mass flux--m sub c
       cmfmc2      (1:iMax,1:kMax+1)    , &! REAL(r8), INTENT(in)  :: state_cmfmc2(ibMax,kMax+1)   ! shallow convective mass flux--m sub c
       concld      (1:iMax,1:kMax)      , &! REAL(r8), INTENT(out) :: state_concld(ibMax,kMax)    ! convective cloud cover
       cld         (1:iMax,1:kMax)      , &! REAL(r8), INTENT(out) :: state_cld   (ibMax,kMax)          !  cloud fraction
       sst         (1:iMax)             , &! REAL(r8), INTENT(in)  :: sst         (ibMax)             ! sea surface temperature
       !zdu         (1:iMax,1:kMax)      , &! REAL(r8), INTENT(in)  :: state_zdu  (ibMax,kMax)          ! detrainment rate from deep convection
       prec_str    (1:iMax)             , &! REAL(r8), INTENT(out)  :: prec_str   (ibMax)  ! [Total] sfc flux of precip from stratiform (m/s) 
       snow_str    (1:iMax)             , &! REAL(r8), INTENT(out)  :: snow_str   (ibMax)  ! [Total] sfc flux of snow from stratiform   (m/s)
       prec_sed    (1:iMax)             , &! REAL(r8), INTENT(out)  :: prec_sed   (ibMax)  ! surface flux of total cloud water from sedimentation
       snow_sed    (1:iMax)             , &! REAL(r8), INTENT(out)  :: snow_sed   (ibMax)  ! surface flux of cloud ice from sedimentation
       prec_pcw    (1:iMax)             , &! REAL(r8), INTENT(out)  :: prec_pcw   (ibMax)  ! sfc flux of precip from microphysics(m/s)
       snow_pcw    (1:iMax)             )  ! REAL(r8), INTENT(out)  :: snow_pcw   (ibMax)  ! sfc flux of snow from microphysics (m/s)
       Total_Rain=Total_Rain+MAX((prec_sed*0.5_r8*dt),0.0_r8)+MAX((prec_pcw*0.5_r8*dt),0.0_r8)
       Total_Snow=Total_Snow+MAX((snow_sed*0.5_r8*dt),0.0_r8)+MAX((snow_pcw*0.5_r8*dt),0.0_r8)
       !CALL qpart(kMax,iMax,t3(1:iMax,1:kMax),ps2(1:iMax),sl(1:kMax),q3(1:iMax,1:kMax),ql3(1:iMax,1:kMax,latco),qi3(1:iMax,1:kMax,latco),'addition')
      CALL qnegat2 (q3, fdqn,(1.0_r8/dt), del, iMax, kMax)! time t+1
    END IF

    ! Save humd/temp after large scale convection
    DO k=1,kMax
      DO i=1,iMax
         q3  (i,k)       = MAX(q3  (i,k)      ,1.0e-12_r8)
         IF (microphys) THEN
         ql3 (i,k,latco) = MAX(ql3 (i,k,latco),1.0e-12_r8)
         qi3 (i,k,latco) = MAX(qi3 (i,k,latco),1.0e-12_r8)
         ql2 (i,k,latco) = MAX(ql2 (i,k,latco),1.0e-12_r8)
         qi2 (i,k,latco) = MAX(qi2 (i,k,latco),1.0e-12_r8)
         END IF
        tLrgs(i,k)=t3(i,k)
        qLrgs(i,k)=q3(i,k)
      END DO
    END DO

    IF (microphys) THEN
       DO k=1,kMax
          DO i=1,iMax
              gliqp(i,k) =ql3     (i,k,latco) 
              gicep(i,k) =qi3     (i,k,latco) 
              gliqm(i,k) =ql2     (i,k,latco) 
              gicem(i,k) =qi2     (i,k,latco) 
          END DO
       END DO   
    END IF

    !-----------------------------------------------------------------
    ! Convective Cloud Cover
    !-----------------------------------------------------------------
    
    IF(TRIM(iccon).EQ.'ARA'.OR.TRIM(iccon).EQ.'KUO'.OR.TRIM(iccon).EQ.'GRE'.OR.TRIM(iccon).EQ.'ZMC')THEN
      CALL CLOUD_COVER( &
             kt     ,ktp    ,iMax   ,kbot   ,ktop   ,noshal ,kctop1 , &
             kcbot1 ,RAINCV ,fac2x  ,rccmbl ,iccon  ,convc  ,convt  , &
             convb  ,prcp1  ,prcp2  ,prcp3  ,prcpt  ,toplv  ,botlv  , &
             convts ,convcs ,convbs )
    END IF

    !-----------------------------------------------------------------
    ! DIAGNOSTICS AND MODEL OUTPUT
    !-----------------------------------------------------------------

    !---------------------
    ! Calculate deep convection moistening and heating profiles
    DO k=1,kMax
       DO i=1,iMax
          clheat(i,k)=fac*rdt*(tDeep(i,k)-tBegin(i,k))
          cmchan(i,k)=fac*rdt*(qDeep(i,k)-qBegin(i,k))
       END DO
    END DO

    !---------------------
    ! Calculate shallow convection moistening and heating profiles
    DO k=1,kMax
       DO i=1,iMax
          sclhea(i,k)=fac*rdt*(tShal(i,k)-tDeep(i,k))
          scmcha(i,k)=fac*rdt*(qShal(i,k)-qDeep(i,k))
       END DO
    END DO

    !---------------------
    ! Calculate large scale convection moistening and heating profiles
    DO k=1,kMax
       DO i=1,iMax
          lslhea(i,k)=fac*rdt*(tLrgs(i,k)-tShal(i,k))
          lsmcha(i,k)=fac*rdt*(qLrgs(i,k)-qShal(i,k))
       END DO
    END DO

    !***********************************
    ! move qDeep to qb and tDeep to tb
    ! move T3 to Tc and q3 to qc
    !***********************************
    DO i=1,iMax
       DO k=1,kMax
          ! Update T
          qb(i,k)=qDeep(i,k)!q  after deep convection
          tb(i,k)=tDeep(i,k)*(1.0_r8+delq*qDeep(i,k))!t  after deep convection
          ! Update T+1
          qc(i,k) = qLrgs(i,k)
          tc(i,k) = tLrgs(i,k)*(1.0_r8+delq*qLrgs(i,k))
       END DO
    END DO
    IF(TRIM(ISCON).EQ.'JHK' .OR. TRIM(ISCON).EQ.'UW' )THEN
       IF(TRIM(ILCON).EQ.'MIC') THEN
          ! Calculate precipiation in mm/s
          DO i=1,iMax
             ppci(i)=2.0e0_r8*1.0e3_r8*RAINCV(i) ! deep
             ppli(i)=2.0e0_r8*1.0e3_r8*(Total_Rain(i)+Total_Snow(i))-ppci(i) ! large
             !ppli(i)=2.0e0_r8*1.0e3_r8*(Total_Rain(i))-ppci(i) ! large
          END DO
       ELSE
          ! Calculate precipiation in mm/s
          DO i=1,iMax
             ppci(i)=2.0e0_r8*1.0e3_r8*RAINCV(i) ! deep
             ppli(i)=2.0e0_r8*1.0e3_r8*(Total_Rain(i)+Total_Snow(i))-ppci(i) ! large
             !ppli(i)=2.0e0_r8*1.0e3_r8*(Total_Rain(i))-ppci(i) ! large
          END DO
       END IF
    ELSE    
       IF(TRIM(ILCON).EQ.'MIC') THEN
       ! Calculate precipiation in mm/s
          DO i=1,iMax
             ppci(i)=2.0e0_r8*1.0e3_r8*RAINCV(i) ! deep
             ppli(i)=2.0e0_r8*1.0e3_r8*(Total_Rain(i)+Total_Snow(i))-ppci(i) ! large
             !ppli(i)=2.0e0_r8*1.0e3_r8*(Total_Rain(i))-ppci(i) ! large
          END DO
       ELSE
          ! Calculate precipiation in mm/s
          DO i=1,iMax
             ppci(i)=2.0e0_r8*1.0e3_r8*RAINCV(i) ! deep
             ppli(i)=2.0e0_r8*1.0e3_r8*(Total_Rain(i))-ppci(i) ! large
          END DO    
       END IF
    END IF

    ! Calculate precipiation in mm

    DO i=1,iMax
       geshem(i)=geshem(i)+fac2x*Total_Rain(i)
    END DO

    ! Time-step output of precipitation
    IF (doprec) THEN
       DO i=1,iMax
          prcc(i)=fac2*rdt*1.0e3_r8*RAINCV(i)
          prct(i)=fac2*rdt*1.0e3_r8*Total_Rain(i)
       END DO
    END IF
    
    ! Diagnose snow field
    IF(TRIM(ILCON).EQ.'MIC') THEN
       DO i=1,iMax          
             snowflg(i) = snow_str(i) !m 
             snowfl (i) = fac2*rdt*1.0e3_r8*Total_Snow(i) !mm/s
       END DO
    ELSE
       DO i=1,iMax          
          IF(0.35_r8*tLrgs(i,1)+0.65_r8*tLrgs(i,2).LE.273.2_r8)THEN
             snowflg(i) = Total_Rain(i)
             snowfl (i) = fac2*rdt*1.0e3_r8*Total_Rain(i)
          ELSE
             snowflg(i) = 0.0_r8
             snowfl (i) = 0.0_r8
          END IF
       END DO
    END IF

    !-----------------
    ! Storage Diagnostic Fields
    !------------------
    IF (StartStorDiag)THEN
       CALL ConvecDiagnStorage (&
            iMax     ,kMax     ,latco     ,rdt      ,fac2     , &
            fdqn     ,RAINCV   ,Total_Rain,snowfl   ,sclhea   , &
            scmcha   ,clheat   ,cmchan    ,lslhea   ,lsmcha   , &
            ql3 (1:iMax,1:kMax,latco),qi3 (1:iMax,1:kMax,latco),&
            cape,cin,LCL,LFC,SWEAT,jdt   )
    END IF
    !-----------------
    ! Storage GridHistory Fields
    !------------------
    IF(ghl_local)THEN 
       CALL ConvecGridHistStorage(&
            iMax     ,kMax      ,latco    ,rdt      ,fac2     , &
            RAINCV   ,Total_Rain,snowflg  ,sclhea   ,scmcha   , &
            clheat   ,cmchan    ,lslhea   ,lsmcha)
    END IF

  END SUBROUTINE cumulus_driver
  


  SUBROUTINE CLOUD_COVER( &
      kt     ,ktp    ,ncols  ,kbot   ,ktop   ,noshal1,kctop1 , &
      kcbot1 ,rrr    ,fac2x  ,rccmbl ,iccon  ,convc  ,convt  , &
      convb  ,prcp1  ,prcp2  ,prcp3  ,prcpt  ,toplv  ,botlv  , &
      convts ,convcs ,convbs )
   IMPLICIT NONE 
    INTEGER, INTENT(in   ) :: kt     
    INTEGER, INTENT(in   ) :: ktp    
    INTEGER, INTENT(in   ) :: ncols  
    INTEGER, INTENT(in   ) :: kbot   (ncols)
    INTEGER, INTENT(in   ) :: ktop   (ncols)
    INTEGER, INTENT(in   ) :: noshal1(ncols)
    INTEGER, INTENT(in   ) :: kctop1 (ncols)
    INTEGER, INTENT(in   ) :: kcbot1 (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: rrr    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: fac2x  
    REAL(KIND=r8),    INTENT(in   ) :: rccmbl ! radiative convective cloud minimum base layer index
    CHARACTER(LEN=*),INTENT(in) :: iccon

    REAL(KIND=r8),    INTENT(inout) :: convc  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: convt  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: convb  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: prcp1  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: prcp2  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: prcp3  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: prcpt  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: toplv  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: botlv  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: convts (ncols)
    REAL(KIND=r8),    INTENT(inout) :: convcs (ncols)
    REAL(KIND=r8),    INTENT(inout) :: convbs (ncols)
     
    REAL(KIND=r8)   , PARAMETER :: fp2457 = 0.2457_r8
    REAL(KIND=r8)   , PARAMETER :: fp1253 = 0.1253_r8
    REAL(KIND=r8)   , PARAMETER :: f0p8   = 0.8_r8
    REAL(KIND=r8)   , PARAMETER :: f8p0e3 = 8.0e3_r8
    INTEGER :: i      
    INTEGER :: is     
    INTEGER :: ijk    

    !--------------------------------------------------------
    !   CLOUD COVER, Cloud TOP-BOT FOR RADIATION (sub cldgen)
    !             DUE CONV PRECIPITATION
    !--------------------------------------------------------
    !   TO cldgen are necessary:
    !   a) cloud top and base   (convt, convb in cldgen)
    !   b) cloud amount is calculated convc (only rrr>0). It is calculate below.
    !     a+b used to defined high clouds due to strong convection
    !   prcpt=precipitation at each time step.(rrr)
    !   convt=ktop
    !   conbt=kbot (>=2) for radiation
    !*****************************************************************

    IF (kt .NE. ktp) THEN
       DO IJK=1,ncols
          convc(IJK) = 0.0_r8   ! call reset(convc(1),ncols)
          convt(IJK) = 0.0_r8   ! call reset(convt(1),ncols)
          convb(IJK) = 0.0_r8   ! call reset(convb(1),ncols)
       ENDDO

       DO i = 1, ncols
          prcpt(i) = prcpt(i) - prcp1(i) &
               + prcp3(i)
       END DO
       IF(TRIM(iccon).EQ.'ARA')THEN
         DO i = 1, ncols 
           IF (prcpt(i) .GT. 0.0e0_r8) THEN
             convc(i) = 0.2_r8+0.038_r8*prcpt(i)*23000.0_r8
             convc(i) = MAX(convc(i), 0.0e0_r8)
             convc(i) = MIN(convc(i), f0p8)
           END IF
         END DO
       ELSE
         DO i = 1, ncols 
           IF (prcpt(i) .GT. 0.0e0_r8) THEN
             convc(i) = fp2457 + fp1253 * LOG(prcpt(i) * f8p0e3)
             convc(i) = MAX(convc(i), 0.0e0_r8)
             convc(i) = MIN(convc(i), f0p8)
           END IF
         END DO        
       END IF  
       !--faltou   
       DO i = 1, ncols
          IF (prcp3(i) .GT. 0.0e0_r8) THEN
             convt(i)=toplv(i) / prcp3(i)
             convb(i)=botlv(i) / prcp3(i)
          END IF
       END DO
       DO i = 1, ncols
          convb(i) = MAX(convb(i),rccmbl)
          IF (convb(i) .GT. convt(i)) &
               convb(i) = convt(i)
       END DO
       !-----
       DO i = 1, ncols
          prcp1(i) = prcp2(i)
          prcp2(i) = prcp3(i)
       END DO
       DO IJK=1,ncols
          prcp3(IJK) = 0.0_r8   !call reset(prcp3(1),ncols)
          toplv(IJK) = 0.0_r8   !call reset(toplv(1),ncols)
          botlv(IJK) = 0.0_r8   !call reset(botlv(1),ncols)
       ENDDO
    END IF
    !*****************************************************
    IF(TRIM(iccon).EQ.'ARA') THEN
      DO i = 1, ncols
        IF (rrr(i) .GT. 0.0_r8) THEN
          prcp3(i) = prcp3(i) + fac2x * rrr(i)
          toplv(i) = toplv(i) + fac2x * rrr(i) * ktop(i)
          botlv(i) = botlv(i) + fac2x * rrr(i) * kbot(i)
        END IF
      END DO
    ELSE
      DO i = 1, ncols
        IF (rrr(i) .GT. 0.0_r8) THEN
          prcp3(i) = prcp3(i) + fac2x * rrr(i)
          toplv(i) = toplv(i) + fac2x * rrr(i) * ktop(i)
          botlv(i) = botlv(i) + fac2x * rrr(i) * kbot(i)
        END IF
      END DO
    END IF
    !
    IF(TRIM(iccon).EQ.'ARA') THEN
       DO IJK=1,ncols
          convts(IJK) = 0.0_r8 
          convbs(IJK) = 0.0_r8 
          convcs(IJK) = 0.0_r8 
       ENDDO

       DO is=1,ncols
          IF(noshal1(is).EQ.0) THEN
             convts(is)=kctop1(is)
             convbs(is)=kcbot1(is)
             !
             !     for mass flux   convcs(is)=0.5
             !
             convcs(is)= 0.3_r8
          ENDIF
       END DO
    ENDIF
  END SUBROUTINE CLOUD_COVER

  SUBROUTINE ConvecGridHistStorage(&
                   ncols    ,kMax     ,latco    ,rdt      ,fac2     , &
                   rrr      ,Total_Rain,snowflg  ,sclhea   ,scmcha   , &
                   clheat   ,cmchan   ,lslhea   ,lsmcha) 
                   
    IMPLICIT NONE
    INTEGER, INTENT(IN   )    :: ncols
    INTEGER, INTENT(IN   )    :: kMax
    INTEGER, INTENT(IN   )    :: latco    
    REAL(KIND=r8), INTENT(IN   ) :: rdt
    REAL(KIND=r8), INTENT(IN   ) :: fac2
    REAL(KIND=r8), INTENT(IN   ) :: rrr       (nCols)
    REAL(KIND=r8), INTENT(in   ) :: Total_Rain (nCols)
    REAL(KIND=r8), INTENT(in   ) :: snowflg   (nCols)
    REAL(KIND=r8), INTENT(in   ) :: sclhea    (nCols,kMax)
    REAL(KIND=r8), INTENT(in   ) :: scmcha    (nCols,kMax)
    REAL(KIND=r8), INTENT(in   ) :: clheat    (nCols,kMax)
    REAL(KIND=r8), INTENT(in   ) :: cmchan    (nCols,kMax)
    REAL(KIND=r8), INTENT(in   ) :: lslhea    (nCols,kMax)
    REAL(KIND=r8), INTENT(in   ) :: lsmcha    (nCols,kMax)

    IF(TRIM(iccon).EQ.'ARA'.OR.TRIM(iccon).EQ.'KUO'.OR.TRIM(iccon).EQ.'GRE'.OR.TRIM(iccon).EQ.'ZMC')THEN
       IF(dogrh(nGHis_cvprec,latco)) &
            CALL StoreGridHistory(rrr(1:ncols),nGHis_cvprec,latco,fac2*rdt*1.0e3_r8)
       IF(dogrh(nGHis_clheat,latco)) &
            CALL StoreGridHistory(clheat(1:nCols,:),nGHis_clheat,latco)
       IF(dogrh(nGHis_cvmosr,latco)) &
            CALL StoreGridHistory(cmchan(1:nCols,:),nGHis_cvmosr,latco)
    END IF

    IF(TRIM(ISCON).EQ.'TIED' .or. TRIM(ISCON).EQ.'SOUZ' .or. TRIM(ISCON).EQ.'JHK'.or. TRIM(ISCON).EQ.'UW')THEN
       IF(dogrh(nGHis_sclhea,latco)) &
            CALL StoreGridHistory(sclhea,nGHis_sclhea,latco)
       IF(dogrh(nGHis_shcvmo,latco)) &
            CALL StoreGridHistory(scmcha,nGHis_shcvmo,latco)
    END IF
    !---------------
    !     gdivn,gtmpn,grotn,gun,gvn are temporary working space
    !     
    IF(TRIM(ILCON).EQ.'YES'.or.TRIM(ILCON).EQ.'LSC'.or.TRIM(ILCON).EQ.'MIC') THEN
       IF(dogrh(nGHis_toprec,latco)) &
            CALL StoreGridHistory(Total_Rain,nGHis_toprec,latco,fac2*rdt*1.0e3_r8)    
       IF(dogrh(nGHis_snowfl,latco)) &
            CALL StoreGridHistory(snowflg  ,nGHis_snowfl,latco,fac2*rdt*1.0e3_r8)
       IF(dogrh(nGHis_sslaht,latco)) &
            CALL StoreGridHistory(lslhea   ,nGHis_sslaht,latco)
       IF(dogrh(nGHis_spstms,latco)) &
            CALL StoreGridHistory(lsmcha   ,nGHis_spstms,latco)
    END IF
  END SUBROUTINE ConvecGridHistStorage

  SUBROUTINE ConvecDiagnStorage(&
                   ncols    ,kMax     ,latco    ,rdt      ,fac2     , &
                   fdqn     ,rrr      ,Total_Rain,snowfl   ,sclhea   , &
                   scmcha   ,clheat   ,cmchan   ,lslhea   ,lsmcha     , &
                   tracerLiq ,tracerIce,cape,cin,LCL,LFC,SWEAT,jdt)    
                                      
    IMPLICIT NONE
    INTEGER, INTENT(IN   )    :: ncols
    INTEGER, INTENT(IN   )    :: kMax
    INTEGER, INTENT(IN   )    :: latco    
    REAL(KIND=r8), INTENT(IN   ) :: rdt
    REAL(KIND=r8), INTENT(IN   ) :: fac2
    REAL(KIND=r8), INTENT(IN   ) :: fdqn      (nCols,kMax)
    REAL(KIND=r8), INTENT(IN   ) :: rrr       (nCols)
    REAL(KIND=r8), INTENT(in   ) :: Total_Rain (nCols)
    REAL(KIND=r8), INTENT(in   ) :: snowfl    (nCols)
    REAL(KIND=r8), INTENT(in   ) :: sclhea    (nCols,kMax)
    REAL(KIND=r8), INTENT(in   ) :: scmcha    (nCols,kMax)
    REAL(KIND=r8), INTENT(in   ) :: clheat (nCols,kMax)
    REAL(KIND=r8), INTENT(in   ) :: cmchan (nCols,kMax)
    REAL(KIND=r8), INTENT(in   ) :: lslhea (nCols,kMax)
    REAL(KIND=r8), INTENT(in   ) :: lsmcha (nCols,kMax)
    REAL(KIND=r8), INTENT(in   ) :: tracerLiq (nCols,kMax)
    REAL(KIND=r8), INTENT(in   ) :: tracerIce(nCols,kMax)    
    REAL(KIND=r8), INTENT(in   ) :: cape(ncols) 
    REAL(KIND=r8), INTENT(in   ) :: cin (ncols) 
    REAL(KIND=r8), INTENT(in   ) :: LCL (ncols) 
    REAL(KIND=r8), INTENT(in   ) :: LFC (ncols) 
    REAL(KIND=r8), INTENT(in   ) :: SWEAT(ncols) 
    INTEGER      , INTENT(IN   ) :: jdt
    
    REAL(KIND=r8)    :: bfr1   (nCols)
    REAL(KIND=r8)    :: bfr3   (nCols)
    INTEGER :: i
    ! "negative specific humidity" correction
    IF(dodia(nDiag_nshcrm))CALL updia(fdqn,nDiag_nshcrm,latco)
    IF(dodia(nDiag_trcliq))CALL updia(tracerLiq,nDiag_trcliq,latco)
    IF(dodia(nDiag_trcice))CALL updia(tracerIce,nDiag_trcice,latco)
    IF(dodia(nDiag_cape2d))CALL updia(cape,nDiag_cape2d,latco)
    IF(dodia(nDiag_cine2d))CALL updia(cin ,nDiag_cine2d,latco)
    IF(dodia(nDiag_sweath))CALL updia(SWEAT,nDiag_sweath,latco)
    DO i=1,ncols
       IF(rrr(i) <= 0.0_r8)THEN
          bfr1(i)=undef
       ELSE
          bfr1(i)=1.0
       END IF
    END DO
    IF(dodia(nDiag_covtcl))CALL updia(bfr1,nDiag_covtcl,latco,jdt)
                                     !(field, loca, lat,jdt)
    IF(TRIM(iccon).EQ.'ARA'.OR.TRIM(iccon).EQ.'KUO'.OR.TRIM(iccon).EQ.'GRE'.OR.TRIM(iccon).EQ.'ZMC')THEN
       IF(dodia(nDiag_cvprec)) THEN
          DO i=1,ncols
             bfr1(i)=fac2*rdt*1.0e3_r8*rrr(i)
          END DO
          CALL updia(bfr1,nDiag_cvprec,latco)
       END IF
       IF(dodia(nDiag_clheat))CALL updia(clheat,nDiag_clheat,latco)
       IF(dodia(nDiag_cmchan))CALL updia(cmchan,nDiag_cmchan,latco)
    END IF

    IF(TRIM(ISCON).EQ.'TIED' .or. TRIM(ISCON).EQ.'SOUZ' .or. TRIM(ISCON).EQ.'JHK'.or. TRIM(ISCON).EQ.'UW')THEN
       IF(dodia(nDiag_sclhea))CALL updia(sclhea,nDiag_sclhea,latco)
       IF(dodia(nDiag_scmcha))CALL updia(scmcha,nDiag_scmcha,latco)
    END IF
    !---------------
    !     gdivn,gtmpn,grotn,gun,gvn are temporary working space
    !     
    IF(TRIM(ILCON).EQ.'YES'.or.TRIM(ILCON).EQ.'LSC'.or.TRIM(ILCON).EQ.'MIC') THEN
       DO i=1,ncols
          IF(dodia(nDiag_toprec))bfr1(i)=fac2*rdt*1.0e3_r8*Total_Rain(i)

          IF(dodia(nDiag_lsprec))bfr3(i)=fac2*rdt*1.0e3_r8*(Total_Rain(i)-rrr(i))
       END DO

       IF(dodia(nDiag_toprec))CALL updia(bfr1,nDiag_toprec,latco)

       IF(dodia(nDiag_snowfl))CALL updia(snowfl,nDiag_snowfl,latco)

       IF(dodia(nDiag_lsprec))CALL updia(bfr3,nDiag_lsprec,latco)

       IF(dodia(nDiag_lslhea))CALL updia(lslhea,nDiag_lslhea,latco)

       IF(dodia(nDiag_lsmcha))CALL updia(lsmcha,nDiag_lsmcha,latco)
    END IF
  END SUBROUTINE ConvecDiagnStorage


  ! qnegat : routine for dealing with negative values of specific humidity
  !          for data on latitude circle.



  SUBROUTINE  qnegat (fq, fdq, fft, rdt, del, iMax, kMax)
    !
    ! input: fq  specific humidity (dimensionless mixing ratio)
    !        fp  surface pressure (cb)
    ! ouput: fq  adjusted specific humidity
    !        fp  unchanged
    !        fdq distribution of moisture modification
    !
    ! iMax......Number of grid points on a gaussian latitude circle   
    ! kMax......Number of sigma levels  
    ! imx.......=iMax+1 or iMax+2   :this dimension instead of iMax
    !              is used in order to avoid bank conflict of memory
    !              access in fft computation and make it efficient. the
    !              choice of 1 or 2 depends on the number of banks and
    !              the declared type of grid variable (real*4,real*8)
    !              to be fourier transformed.
    !              cyber machine has the symptom.
    !              cray machine has no bank conflict, but the argument
    !              'imx' in subr. fft991 cannot be replaced by iMax    
    ! del.......sigma spacing for each layer computed in routine "setsig".  
    ! dfact.....del(k+1)/del(k)
    !
    INTEGER, INTENT(in   ) :: iMax  
    INTEGER, INTENT(in   ) :: kMax
    REAL(KIND=r8)   , INTENT(in   ) :: rdt

    REAL(KIND=r8),    INTENT(inout) :: fq   (iMax,kMax)
    REAL(KIND=r8),    INTENT(inout) :: fdq  (iMax,kMax)  
    REAL(KIND=r8),    INTENT(inout) :: fft  (iMax,kMax)   
    REAL(KIND=r8),    INTENT(in   ) :: del  (kMax)

    REAL(KIND=r8)   :: dfact(kMax)

    INTEGER :: klev
    INTEGER :: kblw
    INTEGER :: i
    INTEGER :: k  

    DO k=1,kMax-1
       dfact(k+1) = del(k+1)/del(k)
    END DO
    !     
    !     ecmwf vertical borrowing scheme
    !     fdq contains compensated borrowing above first level, uncompensated
    !     borrowing in first level
    !     
    DO k=1,kMax-1
       klev = kMax-k+1
       kblw = klev - 1
       DO i=1,iMax
          fdq(i,klev) = fq(i,klev)
          IF(fq(i,klev).LT.0.0e0_r8) fq(i,klev) = 1.0e-12_r8
          fdq(i,klev) = fq(i,klev) - fdq(i,klev)
          fq(i,kblw) = fq(i,kblw) - fdq(i,klev)*dfact(klev)
       END DO
    END DO

    DO i=1,iMax
       fdq(i,1) = fq(i,1)
       IF(fq(i,1).LT.0.0e0_r8) fq(i,1) = 1.0e-12_r8
       fdq(i,1) = fq(i,1) - fdq(i,1)
    END DO

    DO k=1,kMax
       DO i=1,iMax
          fft(i,k)=fft(i,k)/(1.0_r8+delq*fq(i,k))
       END DO
    END DO

    IF(dodia(nDiag_nshcrm))THEN
       DO k=1,kMax
          DO i=1,iMax
             fdq(i,k)=fdq(i,k)*rdt
          END DO
       END DO
    END IF

  END SUBROUTINE qnegat

  SUBROUTINE  qnegat2 (fq, fdq, rdt, del, iMax, kMax)
    !
    ! input: fq  specific humidity (dimensionless mixing ratio)
    !        fp  surface pressure (cb)
    ! ouput: fq  adjusted specific humidity
    !        fp  unchanged
    !        fdq distribution of moisture modification
    !
    ! iMax......Number of grid points on a gaussian latitude circle   
    ! kMax......Number of sigma levels  
    ! imx.......=iMax+1 or iMax+2   :this dimension instead of iMax
    !              is used in order to avoid bank conflict of memory
    !              access in fft computation and make it efficient. the
    !              choice of 1 or 2 depends on the number of banks and
    !              the declared type of grid variable (real*4,real*8)
    !              to be fourier transformed.
    !              cyber machine has the symptom.
    !              cray machine has no bank conflict, but the argument
    !              'imx' in subr. fft991 cannot be replaced by iMax    
    ! del.......sigma spacing for each layer computed in routine "setsig".  
    ! dfact.....del(k+1)/del(k)
    !
    INTEGER, INTENT(in   ) :: iMax  
    INTEGER, INTENT(in   ) :: kMax
    REAL(KIND=r8)   , INTENT(in   ) :: rdt

    REAL(KIND=r8),    INTENT(inout) :: fq   (iMax,kMax)
    REAL(KIND=r8),    INTENT(inout) :: fdq  (iMax,kMax)  
    REAL(KIND=r8),    INTENT(in   ) :: del  (kMax)

    REAL(KIND=r8)   :: dfact(kMax)

    INTEGER :: klev
    INTEGER :: kblw
    INTEGER :: i
    INTEGER :: k  

    DO k=1,kMax-1
       dfact(k+1) = del(k+1)/del(k)
    END DO
    !     
    !     ecmwf vertical borrowing scheme
    !     fdq contains compensated borrowing above first level, uncompensated
    !     borrowing in first level
    !     
    DO k=1,kMax-1
       klev = kMax-k+1
       kblw = klev - 1
       DO i=1,iMax
          fdq(i,klev) = fq(i,klev)
          IF(fq(i,klev).LT.0.0e0_r8) fq(i,klev) = 1.0e-12_r8
          fdq(i,klev) = fq(i,klev) - fdq(i,klev)
          fq(i,kblw) = fq(i,kblw) - fdq(i,klev)*dfact(klev)
       END DO
    END DO

    DO i=1,iMax
       fdq(i,1) = fq(i,1)
       IF(fq(i,1).LT.0.0e0_r8) fq(i,1) = 1.0e-12_r8
       fdq(i,1) = fq(i,1) - fdq(i,1)
    END DO

  END SUBROUTINE qnegat2
  SUBROUTINE qpart(kMax,iMax,t,ps,sl,q,ql,qi,opt)
    IMPLICIT NONE
    INTEGER      , INTENT(IN   ) :: kMax
    INTEGER      , INTENT(IN   ) :: iMax
    REAL(KIND=r8), INTENT(IN   ) :: t    (iMax,kMax)
    REAL(KIND=r8), INTENT(IN   ) :: ps   (iMax)
    REAL(KIND=r8), INTENT(IN   ) :: sl   (kMax)
    CHARACTER(LEN=*), INTENT(IN   ) :: opt

    REAL(KIND=r8), INTENT(INOUT) :: q    (iMax,kMax)
    REAL(KIND=r8), INTENT(INOUT) :: ql   (iMax,kMax)
    REAL(KIND=r8), INTENT(INOUT) :: qi   (iMax,kMax)
    REAL(KIND=r8), PARAMETER :: tmelt      = 273.16_r8             ! freezing T of fresh water ~ K
    REAL(KIND=r8), PARAMETER :: tmax_fice  = tmelt     - 10.0_r8   ! max temperature for cloud ice formation
    REAL(KIND=r8), PARAMETER :: tmin_fice  = tmax_fice - 30.0_r8   ! min temperature for cloud ice formation
    REAL(KIND=r8), PARAMETER :: tmax_fsnow = tmelt                 ! max temperature for transition to convective snow
    REAL(KIND=r8), PARAMETER :: tmin_fsnow = tmelt-5.0_r8           ! min temperature for transition to convective snow

    REAL(KIND=r8)                :: qes(iMax,kMax) 
    REAL(KIND=r8)                :: tsp(iMax,kMax) 
    REAL(KIND=r8)                :: p  (iMax,kMax)
    REAL(KIND=r8)                :: esft
    INTEGER ::i,k
    qi = 0.0_r8
    ql = 0.0_r8
    DO k=1,kMax
       DO i=1,iMax
              p  (i,k) = ps(i)*sl(k)*1000_r8              ! pressure in Pa
       END DO
    END DO   
    CALL findsp (iMax,kMax, q, t, p, tsp, qes)

    IF(TRIM(opt)=='subtraction')THEN
       DO k=1,kMax
          DO i=1,iMax
             !p  (i,k) = ps(i)*sl(k)*10.0_r8               ! pressure in mbar
             !
             ! sgb - IPH is for phase, dependent on TCRIT (water or ice)
             ! calculation of the pressure vapor
             !
             !esft=es5(t(i,k))
             !qes(i,k) = 0.622_r8*esft/(100.0_r8*p(i,k)-esft)
             IF(qes(i,k) <= 1.0e-12_r8  )qes(i,k)=1.0e-12_r8

             IF(q(i,k)   >  qes(i,k))THEN
               ql(i,k) = q (i,k) - qes(i,k)              
               q (i,k) = qes(i,k)
               IF (tsp(i,k) < tmin_fice)THEN
                  qi(i,k) = ql(i,k)
                  ql(i,k) = 0.0_r8
               ELSE
                  qi(i,k) = 0.0_r8
               END IF

             ELSE
               ql(i,k) = 0.0_r8
               qi(i,k) = 0.0_r8
             END IF
          END DO
       END DO
    ELSE IF (TRIM(opt)=='addition')THEN
       DO k=1,kMax
          DO i=1,iMax
               IF(qes(i,k) <= 1.0e-12_r8  )qes(i,k)=1.0e-12_r8
               !IF(q(i,k)   <  qes(i,k))THEN
                  q(i,k) = q (i,k) + ql(i,k) + qi(i,k)                 
                  !IF(q(i,k) > qes(i,k))THEN
                  !  q(i,k) = qes(i,k)
                  !END IF
                  !ql(i,k) = 0.0_r8
                  !qi(i,k) = 0.0_r8
               !END IF          
          END DO
       END DO
    END IF
  END SUBROUTINE qpart
  !---------------------------------
  REAL(KIND=r8) FUNCTION es5(t)
    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN) :: t
    REAL(KIND=r8)   , PARAMETER :: tcrit    =   273.15_r8
    REAL(KIND=r8)   , PARAMETER :: cp   =1004.0_r8
    REAL(KIND=r8)   , PARAMETER :: xl   =2.5e06_r8
    REAL(KIND=r8)   , PARAMETER :: rv   =461.9_r8
    REAL(KIND=r8)            :: ae  (2)
    REAL(KIND=r8)            :: be  (2)
    REAL(KIND=r8)            :: ht  (2)
    ht(1)=xl/cp
    ht(2)=2.834e6_r8/cp
    be(1)=0.622_r8*ht(1)/0.286_r8    
    ae(1)=be(1)/273.0_r8+LOG(610.71_r8)    
    be(2)=0.622_r8*ht(2)/0.286_r8  
    ae(2)=be(2)/273.0_r8+LOG(610.71_r8)

    IF (t <= tcrit) THEN
       es5 = EXP(ae(2)-be(2)/t)
    ELSE
       es5 = EXP(ae(1)-be(1)/t)
    END IF

  END FUNCTION es5
END MODULE Convection
