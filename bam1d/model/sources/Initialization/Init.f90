MODULE Init   ! Version 0 of Nov 25th, 2001



  !  INITIALIZATION PROCEDURE

  !  Initializes all modules required by the computation, at
  !  a precise order. Exports a single procedure, InitAll.
  !  That procedure ought to be invoked before any other
  !  procedure in the computation, since it initializes all
  !  other modules.
  !  See procedure header for arguments description.



  USE Sizes, ONLY : ibMax, jbMax,kMax, &
       jMax,jbMax_ext, jMaxHalf, jMinPerM, iMaxPerJ, &
       mMax, nExtMax, mnMax, mnExtMax, mnExtMap,si, &
       RegisterBasicSizes, RegisterOtherSizes, DumpSizes,InitVerSizes
  USE Transform, ONLY : NextSizeFFT
  USE FieldsDynamics, ONLY: InitFields
  USE Constants, ONLY: InitConstants,i8,r8
  USE Options, ONLY: &
       nscalars

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: InitAll

  INTEGER           , PUBLIC           :: nls
  INTEGER           , PUBLIC           :: nlcs

  LOGICAL(KIND=i8), PARAMETER :: dumpLocal=.FALSE.
  REAL(KIND=r8), PARAMETER :: Sigma_Estrat=0.099999_r8
  INTEGER :: k


CONTAINS



  !InitAll: Initializes all modules used by the computation.
  !         There are four mandatory arguments and four
  !         optional arguments. The four mandatory arguments are:
  !
  !         trunc       model truncation
  !         vert        number of model vertical layers
  !         reducedGrid TRUE iff model uses reduced Grid
  !         linearGrid  TRUE iff model uses linear Grid
  !
  !         If no other argument is given, InitAll computes the
  !         maximum number of longitudes and the number of
  !         latitudes required by the
  !         combination of (reduced or gaussian) with
  !         (linear or quadratic) grids. It proceeds by 
  !         computing Gaussian Points and Weights, placing
  !         latitudes as Gaussian Points and computing Associated
  !         Legendre Functions.
  !         In the case of reduced grid, it proceeds by computing 
  !         the number of longitudes per latitude according to
  !         Courtier-Naughton criteria.
  !         In any case, it follows by accomodating the number of
  !         longitudes per latitude to the acceptable FFT sizes and
  !         ends by initializing MultiFFT and  
  !         MultiLegTrans modules.
  !
  !         Parts of the computation can be bypassed if any of
  !         the optional arguments were given. The optional arguments
  !         and their effect are:
  !
  !Obs. do Jairo: os opcionais serao revistos e alterados, pois
  !permitem inconsistencias e contem informacoes duplicadas.
  !
  !         lonPerLat: an integer vector, indexed by latitudes, 
  !         containing the number of longitudes at each latitude.
  !         Its presence bypass the Courtier-Naughton criteria
  !         computation, but not the accomodation to FFT sizes.
  !           
  !         wavesPerLat: an integer vector, indexed by latitudes,
  !         containing the number of Fourier waves at each latitude.
  !         Its presence bypass the accomodation of Courtier-Naughton
  !         results to the possible sizes of FFTs.
  !
  !         gausPt: a real vector containing all Gaussian Points. Its
  !         size is the number of model latitudes. The presence of this
  !         argument bypasses the computation of Gaussian Points.
  !
  !         gausWg: a real vector containing all Gaussian Weights. Its
  !         size is the number of model latitudes. The presence of this
  !         argument bypasses the computation of Gaussian Weights.

  !Obs do Jairo: esta rotina sera aumentada para conter todas as 
  !inicializacoes do modelo. Esta longe de sua forma final.


  SUBROUTINE InitAll (trunc, vert, reducedGrid, linearGrid, &
         del_in, rhdifd_in, rhdift_in, &
       lonPerLat, wavesPerLat, gausPt, gausWg)
    INTEGER, INTENT(IN)           :: trunc
    INTEGER, INTENT(IN)           :: vert
    LOGICAL, INTENT(IN)           :: reducedGrid 
    LOGICAL(KIND=i8), INTENT(IN)           :: linearGrid 
    REAL(KIND=r8),    INTENT(IN)           :: del_in(:)
    REAL(KIND=r8),    INTENT(IN)           :: rhdifd_in
    REAL(KIND=r8),    INTENT(IN)           :: rhdift_in
    INTEGER, INTENT(IN), OPTIONAL :: lonPerLat(:)
    INTEGER, INTENT(IN), OPTIONAL :: wavesPerLat(:)
    REAL(KIND=r8),    INTENT(IN), OPTIONAL :: gausPt(:)
    REAL(KIND=r8),    INTENT(IN), OPTIONAL :: gausWg(:)

    CHARACTER(LEN=*), PARAMETER   :: h="**(InitAll)**"
    CHARACTER(LEN=10) :: c1, c2
    INTEGER :: wavesFactor
    INTEGER :: nLon
    INTEGER :: nLat, left,mymnMax, mymnExtMax, &
                    MNMax_si
    INTEGER :: j
    INTEGER, ALLOCATABLE :: mPerLat(:)
    INTEGER, ALLOCATABLE :: iMaxPerLat(:)
    mymnMax=1
    mymnExtMax=1 
    MNMax_si=1
    IF (dumpLocal) THEN
       WRITE(c1,"(i10)") trunc
       WRITE(c2,"(i10)") vert
       WRITE(*,"(a,' trunc=',a,'; vert=',a)",ADVANCE='NO') &
            h, TRIM(ADJUSTL(c1)), TRIM(ADJUSTL(c2))
       IF (reducedGrid) THEN
          WRITE(*,"('; Reduced and')",ADVANCE='NO')
       ELSE
          WRITE(*,"('; Gaussian and')",ADVANCE='NO')
       END IF
       IF (linearGrid) THEN
          WRITE(*,"(' Linear Grid')")
       ELSE
          WRITE(*,"(' Quadratic Grid')")
       END IF
      ! CALL flush(*)
    END IF

    ! linear or quadratic grid

    IF (linearGrid) THEN
       wavesFactor = 2
    ELSE
       wavesFactor = 3
    END IF

    ! maximum number of longitudes
    ! for nLat=nLon/2 with even nLat, it suffices nLon multiple of 4

    IF (PRESENT(lonPerLat)) THEN
       nLon = NextSizeFFT(MAXVAL(lonPerLat))
    ELSE 
       nLon = NextSizeFFT(wavesFactor*trunc + 1)
       DO
          left = MOD(nLon,4)
          IF (left == 0) EXIT
          nLon = NextSizeFFT(nLon + 4 - left)
       END DO
    END IF

    ! number of latitudes

    IF (PRESENT(lonPerLat)) THEN
       nLat = SIZE(lonPerLat)
    ELSE
       nLat = nLon/2
    END IF

    ! verify consistency and register basic sizes

    CALL RegisterBasicSizes(trunc, nLat, nLon, vert)
    IF (dumpLocal) THEN
       WRITE(*,"(a,' will DumpSizes just after RegisterBasicSizes')") h
       CALL DumpSizes()
       !CALL flush(6)
    END IF

    CALL InitVerSizes(del_in)

    nls=0
    DO k=1,kMax
       IF (si(k) <= Sigma_Estrat) nls=nls+1
    END DO
    IF (nls == 0) nls=1
    nlcs=kMax+2

    ! computes mPerLat

    ALLOCATE (mPerLat(jMax))
    IF (PRESENT(wavesPerLat)) THEN
       mPerLat = wavesPerLat
    ELSE
       mPerLat = mMax
    END IF

    ! compute iMaxPerLat

    ALLOCATE(iMaxPerLat(jMax))
    IF (PRESENT(lonPerLat)) THEN
       iMaxPerLat = lonPerLat
    ELSE 
       DO j = 1, jMaxHalf
          iMaxPerLat(j) = 1!NextSizeFFT(wavesFactor*(mPerLat(j)-1) + 1)
       END DO
          iMaxPerLat(jMaxHalf+1:jMax)=1!iMaxPerLat(jMaxHalf:1:-1)
    END IF
    
    ! finish building sizes
    
    CALL RegisterOtherSizes(iMaxPerLat, mPerLat)
    
    IF (dumpLocal) THEN
       WRITE(*,"(a,' will DumpSizes just after RegisterOtherSizes')") h
       CALL DumpSizes()
    END IF

    ! for now on, all global constants defined at sizes can be used
    
    ! initialize Fields, mnMax, mnExtMax

    CALL InitFields(ibMax, kMax, jbMax, mymnMax, mymnExtMax, kmax, &
                    MNMax_si,jbMax_ext, nscalars)

    ! initialize Constants
    
    CALL InitConstants(kMax)
  END SUBROUTINE InitAll
END MODULE Init
