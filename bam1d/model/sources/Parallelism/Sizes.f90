MODULE Sizes   ! Version 0 of Nov 25th, 2001

 USE Constants, ONLY: rk,i8,r8,i4,r4


  ! SINGLE REPOSITORY OF SIZES AND MAPPINGS OF GRID, FOURIER AND SPECTRAL
  ! FIELD REPRESENTATIONS. SHOULD BE USED BY ALL MODEL MODULES THAT REQUIRE
  ! SUCH INFORMATION.
  !
  ! Module exports three procedures:
  ! RegisterBasicSizes: set all repository values for spectral representation.
  !                     It also sets latitudes and maximum number of longitudes
  !                     for grid representation.
  ! RegisterOtherSizes: set remaining repository values.
  ! DumpSizes:          detailed dumping of repository values
  !
  ! Module exports a set of values and arrays, described bellow.
  !
  ! Notation used throughout the code:
  ! m  Legendre Order (also Fourier wave number). Goes from 1 to mMax.
  ! n  Legendre Degree. Goes from m to nMax on regular spectral fields
  !    and from m to nExtMax on extended spectral fields.
  ! mn index to store the spectral triangle (m,n) into a single dimension.
  ! k  vertical index, from 1 to kMax
  ! j  latitude index, from 1 to jMax (full field)
  !    or from 1 to jMaxHalf (hemisphere)
  ! i  longitude index, from 1 to iMax (gaussian grid) or
  !    from 1 to iMaxPerJ(j) (at latitude j of the reduced grid)
  ! ib index of a longitude on a block of longitudes packed in one
  !    dimension of grids;
  ! jb index for which block of longitudes to take, packed in another
  !    dimension of grids

  ! SPECTRAL REPRESENTATION:
  ! A regular spectral field should be dimensioned (2*mnMax,kMax), where
  ! the first dimension accomodates pairs of (real,imaginary) spectral
  ! coefficients for a fixed vertical and the second dimension varies
  ! verticals.
  ! Extended spectral fields should use (2*mnExtMax,kMax).
  ! This module provides mapping functions to map (m,n) into mn:
  ! For a regular spectral field with complex entries, Legendre Order m
  ! and Degree n has real coefficient at position 2*mnMap(m,n)-1 and
  ! imaginary coefficient at position 2*mnMap(m,n). Inverse mappings 
  ! (from mn to m and n) are also provided. Maps for the extended spectral 
  ! field are also provided.

  ! GRID REPRESENTATION:
  ! Since the number of longitudes per latitude may vary with latitude
  ! (on reduced grids), sets of latitudes (with all longitudes) 
  ! are packed together in the first dimension of Grids. 
  ! Near the pole, many latitudes are packed; near the equator, just a few. 
  ! Second dimension is vertical.
  ! Third dimension is the number of packed latitudes required to represent
  ! a full field.
  ! A grid field should be dimensioned (ibMax, kMax, jbMax).
  ! This module provides mapping functions to map latitude j and longitude
  ! i into (ib,jb).
  ! Map ibPerIJ(i,j) gives index of first dimension that stores longitude
  ! i of latitude j. Map jbPerIJ(i,j) gives index of third dimension that
  ! stores longitude i of latitude j. Consequently, the point of longitude
  ! i, latitude j and vertical k is stored at (ibPerIJ(i,j), k, jbPerIJ(i,j)).
  ! Inverse mappings (from (ib,jb) to (i,j)) are also provided.

  ! FOURIER REPRESENTATION:
  ! For the moment, Fourier fields are represented externally to the
  ! transform. That should not happen in the future.
  ! First dimension contains pairs of (real,imaginary) fourier 
  ! coefficients. Second dimension is latitude. Third dimension is
  ! vertical. A full fourier field should be dimensioned 
  ! (iMax+1, jMax, kMax)


  IMPLICIT NONE


  
  ! SPECTRAL REPRESENTATION:
  ! mMax is the maximum Legendre Order (also maximum Fourier wave number).
  !      it is set to model truncation + 1 (Legendre Order 0 is m=1)
  ! nMax is the maximum Legendre Degree for regular spectral fields.
  !      it is set to model truncation + 1 (Legendre Degree 1 is n=1)
  ! mnMax is the amount of spectral coeficients per latitude, for 
  !       regular spectral fields. It is the number of points in the
  !       regular triangular spectral plane.
  ! mMap indexed by mn in [1:mnMax]; returns m at this mn for regular
  !      spectral fields;
  ! nMap indexed by mn in [1:mnMax]; returns n at this mn for regular
  !      spectral fields;
  ! mnMap indexed by (m,n); returns mn that stores spectral coefficient
  !       (m,n), when spectral coefficients are real (e.g. Associated 
  !       Legendre Functions). For complex spectral coefficients, 
  !       mn=2*mnMap(m,n) for the imaginary part and mn=2*mnMap(m,n)-1
  !       for the real part.
  ! Remaining variables with Ext extension have the same meaning, but
  ! for the extended spectral space (where nExtMax=trunc+2), that is,
  ! where each Legendre Order has one more Legendre Degree than usual.



  INTEGER              :: mMax=-1
  INTEGER              :: nMax=-1
  INTEGER              :: mnMax=-1
  INTEGER, ALLOCATABLE :: mMap(:)
  INTEGER, ALLOCATABLE :: nMap(:)
  INTEGER, ALLOCATABLE :: mnMap(:,:)
  INTEGER, ALLOCATABLE :: snnp1(:)
  INTEGER              :: nExtMax=-1
  INTEGER              :: mnExtMax=-1
  INTEGER, ALLOCATABLE :: mExtMap(:)
  INTEGER, ALLOCATABLE :: nExtMap(:)
  INTEGER, ALLOCATABLE :: mnExtMap(:,:)



  ! LATITUDES REPRESENTATION:
  ! jMax is the number of latitudes for full Fourier and Grid representations.
  ! jMaxHalf is the number of latitudes at each hemisphere.
  ! mMaxPerJ is an array indexed by latitudes that stores the maximum value
  !          of m at each latitude. Consequently, the contribution of 
  !          latitude j for the Legendre Transform should be taken only
  !          up to Legendre Order mMaxPerJ(j). By the same token, 
  !          FFTs at latitude j should be computed only up to 
  !          Fourier wave number mMaxPerJ(j).
  !          For regular grids, mMaxPerJ(j) = mMax for all j.
  ! jMinPerM is the inverse mapping - an array indexed by Legendre Orders (m)
  !          containing the smallest latitude (value of j) that has that
  !          order. Latitudes that contain Legendre Order (and Fourier
  !          wave number) m are jMinPerM(m) <= j <= jMax-jMinPerM(m)+1
  !          For regular grids, jMinPerM(m) = 1 for all m.



  INTEGER              :: jMax=-1
  INTEGER              :: jMaxHalf=-1
  INTEGER, ALLOCATABLE :: mMaxPerJ(:)
  INTEGER, ALLOCATABLE :: jMinPerM(:)



  ! LONGITUDES REPRESENTATION:
  ! iMax is the maximum number of longitudes per latitude;
  !      it is the number of longitudes per latitude for all latitudes
  !      on regular grids;
  !      it is only the maximum number of longitudes per latitude on 
  !      reduced grids; it is the actual number of longitudes per latitude
  !      close to the equator on reduced grids;
  ! ijMax is the number of horizontal grid points at regular or reduced
  !       grids.
  ! iMaxPerJ is the actual number of longitudes per latitude on regular
  !          and reduced grids; latitude j has iMaxPerJ(j) longitudes.



  INTEGER              :: iMax=-1
  INTEGER              :: ijMax=-1
  INTEGER              :: ijMaxGauQua=-1
  INTEGER, ALLOCATABLE :: iMaxPerJ(:)



  ! GRID REPRESENTATION:
  ! All longitudes of a set of latitudes are packed together in the
  ! first dimension of grids. That decreases the waste of memory when
  ! comparing to store one latitude only, for the case of reduced
  ! grid. 
  ! ibMax is the maximum number of longitudes packed into the first dimension
  !       of grids. The actual number of longitudes vary with the third
  !       dimension of the grid representation.
  ! jbMax is the number of sets of longitudes required to store an
  !       entire field.
  ! ibPerIJ maps longitude i and latitude j into the first dimension
  !         of grid representation. It is indexed ibPerIJ(i,j).
  ! jbPerIJ maps longitude i and latitude j into the third dimension
  !         of grid representation. It is indexed jbPerIJ(i,j).
  ! iPerIJB gives which longitude is stored at first dimension index
  !         i and third dimension index j of Grid representations.
  !         It is indexed iPerIJB(i,j)
  ! jPerIJB gives which latitude is stored at first dimension index
  !         i and third dimension index j of Grid representations.
  !         It is indexed jPerIJB(i,j)
  ! ibMaxPerJB gives how many latitudes are actually stored at third
  !            dimension jb. Since the number of longitudes vary with
  !            latitudes, the amount of space actually used in the first
  !            dimension of grid representations vary with the third
  !            dimension. Array ibMaxPerJB, indexed by jb, accounts for
  !            such variation.


  INTEGER              :: ibMax=-1
  INTEGER              :: jbMax=-1
  INTEGER              :: jbMax_ext
  INTEGER, ALLOCATABLE :: ibPerIJ(:,:)
  INTEGER, ALLOCATABLE :: jbPerIJ(:,:)
  INTEGER, ALLOCATABLE :: iPerIJB(:,:)
  INTEGER, ALLOCATABLE :: jPerIJB(:,:)
  INTEGER, ALLOCATABLE :: ibMaxPerJB(:)

  INTEGER              :: kMax=-1

  REAL(KIND=r8), ALLOCATABLE :: ci(:)      ! 1 - sigma each level (level 1 at surface)
  REAL(KIND=r8), ALLOCATABLE :: si(:)      ! sigma
  REAL(KIND=r8), ALLOCATABLE :: del(:)     ! layer thickness (in sigma)
!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
  REAL(KIND=r8), ALLOCATABLE :: delcl(:)     ! layer thickness (in cl)
!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

  REAL(KIND=r8), ALLOCATABLE :: rdel2(:)
  REAL(KIND=r8), ALLOCATABLE :: sl(:)      ! sigma at layer midpoint
  REAL(KIND=r8), ALLOCATABLE :: cl(:)      ! 1.0 - sl
  REAL(KIND=r8), ALLOCATABLE :: rpi(:)     ! 'pi' ratios at adjacent layers


  LOGICAL(KIND=i8), PARAMETER, PRIVATE :: dumpLocal=.TRUE.
CONTAINS

  SUBROUTINE InitVerSizes(del_in)
    REAL(KIND=r8), INTENT(IN) :: del_in(:)
    INTEGER :: kMax
    INTEGER :: i, j, k
    REAL(KIND=r8)    :: rk1, sirk, sirk1, dif
    kMax = SIZE(del_in)
    ALLOCATE(del(kMax))
!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
    ALLOCATE(delcl(kMax-1))
!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
    ALLOCATE(rdel2(kMax))
    ALLOCATE(ci(kMax+1))
    ALLOCATE(si(kMax+1))
    ALLOCATE(sl(kMax))
    ALLOCATE(cl(kMax))
    ALLOCATE(rpi(kMax-1))

    del = del_in
    rdel2 = 0.5_r8/del

    rk1 = rk + 1.0_r8

    !cdir novector
    ci(1) = 0.0_r8
    DO k=1, kmax-1
       ci(k+1)=ci(k)+del(k)
    END DO
    ci(kmax+1)=1.0_r8

    DO k=1, kmax+1
       si(k) = 1.0_r8 - ci(k)
    END DO

    DO k=1, kmax
       sirk =EXP(rk1*LOG(si(k)))
       IF(k.LE.kmax-1) THEN
          sirk1=EXP(rk1*LOG(si(k+1)))
       ELSE
          sirk1=0.0_r8
       END IF
       dif = sirk-sirk1
       dif = dif / (rk1*(si(k)-si(k+1)))
       sl(k) = EXP(LOG(dif)/rk)
       cl(k) = 1.0_r8 - sl(k)
    END DO
!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
    DO k = 1, kmax-1
       delcl(k) = cl(k+1) - cl(k)
    END DO
!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

    !     
    !     compute pi ratios for temp. matrix.
    !
    DO k=1, kmax-1
       rpi(k) = EXP(rk*LOG(sl(k+1)/sl(k)))
    END DO
  END SUBROUTINE InitVerSizes

  SUBROUTINE RegisterBasicSizes(trunc, nLat, nLon, vert)
    INTEGER, INTENT(IN) :: trunc
    INTEGER, INTENT(IN) :: nLat
    INTEGER, INTENT(IN) :: nLon
    INTEGER, INTENT(IN) :: vert
    INTEGER :: m, n, mn, diag, ele
    CHARACTER(LEN=*), PARAMETER :: h="**(RegisterBasicSizes)**"

    jMax     = 1!nLat
    jMaxHalf = 1!nLat/2

    iMax     = 1!nLon

    kMax     = vert

    mMax     = trunc + 1
    nMax     = mMax
    nExtMax  = mMax + 1
    mnExtMax = (nExtMax+2)*(nExtMax-1)/2
    ALLOCATE (mExtMap(mnExtMax))
    ALLOCATE (nExtMap(mnExtMax))
    ALLOCATE (mnExtMap(mMax,nExtMax))
    nExtMap = -1  ! flag mapping error
    mExtMap = -1  ! flag mapping error
    mnExtMap = -1 ! flag mapping error
    mn = 0
    DO diag = 1, nExtMax ! diagonal
       DO ele = 1, MIN(mMax-diag+2, mMax)
          m = ele
          n = ele+diag-1
          mn = mn + 1
          mnExtMap(m,n) = mn
          mExtMap(mn)   = m
          nExtMap(mn)   = n
       END DO
    END DO
    mnMax = (mMax * (nMax+1))/2
    ALLOCATE (mnMap(mMax,nMax))
    ALLOCATE (mMap(mnMax))
    ALLOCATE (nMap(mnMax))
    mnMap = -1  ! flag mapping error
    mMap = -1   ! flag mapping error
    nMap = -1   ! flag mapping error
    mn = 0
    DO diag = 1, nMax 
       DO ele = 1, mMax-diag+1
          m = ele
          n = m+diag-1
          mn = mn + 1
          mnMap(m,n) = mn
          mMap(mn)   = m
          nMap(mn)   = n
       END DO
    END DO
    ALLOCATE (snnp1(2*mnMax))
    DO m = 1, mMax
       DO n = m, nMax
          mn = mnMap(m,n)
          snnp1(2*mn-1) = (n-1)*n
          snnp1(2*mn  ) = (n-1)*n
       END DO
    END DO
    IF (dumpLocal) THEN
       WRITE(*,"(a,' Dump at the end ')") h
       CALL DumpSizes()
    END IF
    IF (trunc == 021 .OR. trunc == 030 .OR. trunc == 042 .OR. &
        trunc == 047 .OR. trunc == 062 .OR. trunc == 079 .OR. &
        trunc == 085 .OR. trunc == 094 .OR. trunc == 106 .OR. &
        trunc == 126 .OR. trunc == 159 .OR. trunc == 170 .OR. &
        trunc == 213 .OR. trunc == 254 .OR. trunc == 319) THEN
        ijMaxGauQua = iMax*jMax
    ELSE
       WRITE(*,"(a,' nao ha ijMaxGauQua para esta resolucao do modelo')") h
       STOP h
    END IF
!    IF (trunc == 62 .AND. vert == 28) THEN
!       ijMaxGauQua = 18432
!    ELSE IF (trunc == 126 .AND. vert == 28) THEN
!       ijMaxGauQua = 73728
!    ELSE IF (trunc == 170 .AND. vert == 42) THEN
!       ijMaxGauQua = 131072
!    ELSE
!       WRITE(*,"(a,' nao ha ijMaxGauQua para esta resolucao do modelo')") h
!       STOP h
!    END IF

  END SUBROUTINE RegisterBasicSizes
  SUBROUTINE RegisterOtherSizes(iMaxPerLat, mPerLat)
    INTEGER, INTENT(IN) :: iMaxPerLat(jMax)
    INTEGER, INTENT(IN) :: mPerLat(jMax)
    CHARACTER(LEN=*), PARAMETER :: h="**(RegisterOtherSizes)**"
    INTEGER, PARAMETER :: MinLatPerBlk=1
    INTEGER :: i, j, m, ib, jb, cnt
    ALLOCATE (mMaxPerJ(jMax))
    mMaxPerJ = mPerLat
    ALLOCATE (iMaxPerJ(jMax))
    iMaxPerJ = iMaxPerLat
    !IF (iMax /= MAXVAL(iMaxPerJ)) THEN
    !   STOP ' imax and imaxperj disagree'6923
    !END IF
    ijMax = SUM(iMaxPerJ)
    ALLOCATE (jMinPerM(mMax))
    jMinPerM = jMaxHalf
    DO j = 1, jMaxHalf
       m = mMaxPerJ(j)
       jMinPerM(1:m) = MIN(j, jMinPerM(1:m))
    END DO

    ! # longitudes per block

    ibMax = MinLatPerBlk*iMax
    
    ! # blocks

    jbMax = 1
    cnt = 0
    DO j = 1, jMax
       cnt = cnt + iMaxPerJ(j)
       IF (cnt > ibMax) THEN
          jbMax = jbMax + 1
          cnt = iMaxPerJ(j)
       ELSE IF (cnt == ibMax) THEN
          jbMax = jbMax + 1
          cnt = 0
       END IF
    END DO
    IF (cnt == 0) THEN
       jbMax = jbMax - 1
    END IF

    ! maps (i,j) into (ib,jb) and vice-versa

    ALLOCATE (ibPerIJ(iMax ,jMax ))
    ibPerIJ=-1
    ALLOCATE (jbPerIJ(iMax ,jMax ))
    jbPerIJ=-1
    ALLOCATE (iPerIJB(ibMax,jbMax))
    iPerIJB=-1
    ALLOCATE (jPerIJB(ibMax,jbMax))
    jPerIJB=-1
    ALLOCATE (ibMaxPerJB(jbMax))
    ibMaxPerJB=-1

    jb = 1
    ib = 0
    DO j = 1, jMax
       IF (ib + iMaxPerJ(j) > ibMax) THEN
          jb = jb + 1
          ib = 0
       END IF
       ibPerIJ(   1:   iMaxPerJ(j),  j) =1! (/ (i, i=ib+1,ib+iMaxPerJ(j)) /)
       jbPerIJ(   1:   iMaxPerJ(j),  j) =1! jb
       iPerIJB(ib+1:ib+iMaxPerJ(j), jb) =1! (/ (i, i=1,iMaxPerJ(j)) /)
       jPerIJB(ib+1:ib+iMaxPerJ(j), jb) =1! j
       ib = ib + iMaxPerJ(j)
    END DO

    DO j = 1, jMax
       ibMaxPerJB(jbPerIJ(iMaxPerJ(j),j)) = 1!ibPerIJ(iMaxPerJ(j),j)
    END DO
    jbmax_ext = jb

    IF (dumpLocal) THEN
       WRITE(*,"(a,' Dump at the end ')") h
       CALL DumpSizes()
    END IF
  END SUBROUTINE RegisterOtherSizes
  SUBROUTINE DumpSizes()
    CHARACTER(LEN=*), PARAMETER :: h="**(DumpSizes)**"
    CHARACTER(LEN=10) :: c1, c2, c3, c4, c5, c6
    LOGICAL(KIND=i8) :: Mask(jMaxHalf)
    INTEGER :: first(1)
    INTEGER :: jMaxPerM(mMax)
    INTEGER :: i, j, jj, lastjj, jb, m
    INTEGER :: firstLat, lastLat, firstM, lastM, firstJ, lastJ, lastIB
    IF (mMax == -1) THEN
       WRITE(*,"(a,' Sizes not created')") h
    ELSE
       WRITE(c1,"(i10)") mMax
       WRITE(c2,"(i10)") nMax
       WRITE(c3,"(i10)") nExtMax
       WRITE(c4,"(i10)") mnExtMax
       WRITE(*,"(a,' mMax=',a,'; nMax=',a,'; nExtMax=',a,'; mnExtMax=',a)") &
            h, TRIM(ADJUSTL(c1)), TRIM(ADJUSTL(c2)), TRIM(ADJUSTL(c3)), &
            TRIM(ADJUSTL(c4))
       WRITE(c1,"(i10)") jMax
       WRITE(c2,"(i10)") jMaxHalf
       WRITE(*,"(a,' jMax=',a,'; jMaxHalf=',a)") &
            h, TRIM(ADJUSTL(c1)), TRIM(ADJUSTL(c2))
       WRITE(c1,"(i10)") iMax
       WRITE(c2,"(i10)") kMax
       WRITE(c3,"(i10)") ijMax
       IF (ijMax == -1) THEN
          WRITE(*,"(a,' iMax=',a,'; kMax=',a)") h, &
               TRIM(ADJUSTL(c1)), TRIM(ADJUSTL(c2))
          WRITE(*,"(a,' Sizes not fully created')") h
       ELSE
          WRITE(*,"(a,' iMax=',a,'; kMax=',a,'; ijMax=',a)") h, &
               TRIM(ADJUSTL(c1)), TRIM(ADJUSTL(c2)), TRIM(ADJUSTL(c3))
          WRITE(*,"(' latitudes   longitudes   m')")
          DO i = MINVAL(iMaxPerJ), MAXVAL(iMaxPerJ)-1
             Mask = iMaxPerJ(1:jMaxHalf) == i
             IF (ANY(Mask)) THEN
                first = MINLOC(iMaxPerJ(1:jMaxHalf), Mask)
                firstLat = first(1)
                lastLat = firstLat + COUNT(Mask) - 1
                WRITE(*,"(i5,':',i3,5x,i5,6x,'0:',i3)") &
                     firstLat, lastLat, i, mMaxPerJ(lastLat)-1
             END IF
          END DO
          i = MAXVAL(iMaxPerJ)
          first = MINLOC(iMaxPerJ, iMaxPerJ==i)
          firstLat = first(1)
          lastLat = jMax-firstLat+1
          WRITE(*,"(i5,':',i3,5x,i5,6x,'0:',i3)") &
               firstLat, lastLat, i, mMaxPerJ(lastLat)-1
          DO i = MAXVAL(iMaxPerJ)-1, MINVAL(iMaxPerJ), -1
             Mask = iMaxPerJ(jMaxHalf+1:jMax) == i
             IF (ANY(Mask)) THEN
                first = MINLOC(iMaxPerJ(jMaxHalf+1:jMax), Mask)
                firstLat = first(1) + jMaxHalf
                lastLat = firstLat + COUNT(Mask) - 1
                WRITE(*,"(i5,':',i3,5x,i5,6x,'0:',i3)") &
                     firstLat, lastLat, i, mMaxPerJ(lastLat)-1
             END IF
          END DO


          WRITE(*,"('     m       latitudes')")
          jMaxPerM = jMax - jMinPerM + 1

          lastM = 0
          DO
             IF (lastM == mMax) THEN
                EXIT
             ELSE
                firstM=lastM+1
                lastM=firstM
                firstJ=jMinPerM(firstM)
                lastJ=jMaxPerM(firstM)
                DO m = firstM+1, mMax
                   IF ((firstJ == jMinPerM(m)) .AND. &
                        (lastJ == jMaxPerM(m))) THEN
                      lastM = m
                   ELSE
                      EXIT
                   END IF
                END DO
                WRITE(*,"(i5,':',i3,5x,i3,':',i3)") &
                     firstM-1, lastM-1, firstJ, lastJ
             END IF
          END DO

          WRITE(c1,"(i10)") ibMax
          WRITE(c2,"(i10)") jbMax
          WRITE(c3,"(i10)") jbMax*ibMax
          WRITE(*,"(a,' ibMax=',a,'; jbMax=',a,'; ijbMax=',a)") h, &
               TRIM(ADJUSTL(c1)), TRIM(ADJUSTL(c2)), TRIM(ADJUSTL(c3))
          DO j = 1, jMax
             IF (jbPerIJ(1,j) /= jbPerIJ(iMaxPerJ(j),j)) THEN
                STOP "error in mapping jbPerIJ"
             END IF
             WRITE(c1,"(i10)") 1
             WRITE(c2,"(i10)") iMaxPerJ(j)
             WRITE(c3,"(i10)") j
             WRITE(c4,"(i10)") ibPerIJ(1,j)
             WRITE(c5,"(i10)") ibPerIJ(iMaxPerJ(j),j)
             WRITE(c6,"(i10)") jbPerIJ(1,j)
             WRITE(*,"(a,'(Long,Lat)=(',a,':',a,',',a,') &
                  &mapped into(ib,jb)=(',&
                  &a,':',a,',',a,')')") h, &
                  TRIM(ADJUSTL(c1)), TRIM(ADJUSTL(c2)), &
                  TRIM(ADJUSTL(c3)), TRIM(ADJUSTL(c4)), &
                  TRIM(ADJUSTL(c5)), TRIM(ADJUSTL(c6))
          END DO


          DO jb = 1, jbMax
             lastIb = ibMaxPerJB(jb)
             WRITE(c1,"(i10)") 1
             WRITE(c2,"(i10)") lastIb
             WRITE(c3,"(i10)") jb
             WRITE(*,"(a,'(ib,jb)=(',a,':',a,',',a,') maps (Long,Lat)=')") &
                  h, TRIM(ADJUSTL(c1)), TRIM(ADJUSTL(c2)), TRIM(ADJUSTL(c3))
             firstJ = jPerIJB(1,jb)
             lastJ  = jPerIJB(lastIb,jb)
             DO jj = firstJ, lastJ, 5
                lastjj = MIN(jj+4,lastJ)
                DO j = jj, lastjj
                   WRITE(c4,"(i10)") 1
                   WRITE(c5,"(i10)") iMaxPerJ(j)
                   WRITE(c6,"(i10)") j
                   WRITE(*,"('(',a,':',a,',',a,')')",ADVANCE='no') &
                        TRIM(ADJUSTL(c4)), TRIM(ADJUSTL(c5)), TRIM(ADJUSTL(c6))
                   IF (j /= lastjj) THEN
                      WRITE(*,"(',')", ADVANCE='no')
                   ELSE IF (j /= lastj) THEN
                      WRITE(*,"(',')")
                   ELSE
                      WRITE(*,"('')")
                   END IF
                END DO
             END DO
          END DO
       END IF
    END IF
  END SUBROUTINE DumpSizes
END MODULE Sizes
