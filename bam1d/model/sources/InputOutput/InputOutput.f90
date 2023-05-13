MODULE InputOutput
  USE Constants, ONLY: i8,r8,i4,r4

  USE IOLowLevel, ONLY: &
      ReadHead     , &  
      GReadHead    , & 
      ReadField    , &
      GReadField   , &
      WriteHead    , &
      GWriteHead   , &
      WriteField   , &
      GWriteField  , &
      ReadGetALB   , &
      ReadGetSST   , &
      ReadGetSLM   , &
      ReadGetSLM3D , &
      ReadGetSNW   , &
      ReadGetSST2  , &
      ReadOzone    , &
      ReadTracer


  USE Utils, ONLY: &
       IJtoIBJB      ,&
       AveBoxIJtoIBJB,&
       NearestIJtoIBJB

   USE Options, ONLY: &
       nfprt, &
       nfctrl, &
       nfsst, &
       nfndvi, &
       nfozone, &
       nfSoilMostSib2,&
       reducedGrid,&
       labelsi,&
       labelsj,&
       fNameSnow, &
       fNameSibAlb,&
       fNameAlbedo, &
       fNameSoilms,&
       fNameSSTAOI, &
       fNameNDVIAOI,&
       fNameSoilMoistSib2,&
       fNameOzone , &
       fNametracer, &
       ifco2, ifozone, & !hmjb
       nfco2, nfozone, & !hmjb
       nftrc,iftracer, &
       schemes    , &
       co2val, &             !hmjb for new co2 values
       monl     

!-1D  USE ModTransDiagCol, ONLY: &
!-1D      TransDiagCol

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: InitInputOutput
  PUBLIC :: cnvray
  PUBLIC :: wrfld
  PUBLIC :: sclout
  PUBLIC :: transp
  PUBLIC :: transw
  PUBLIC :: getsbc
  PUBLIC :: getsbc1D
  
  PUBLIC :: nferr
  PUBLIC :: ifprt
  PUBLIC :: fsbc
  INTEGER              :: nferr
  INTEGER              :: ifprt(100)
  INTEGER              :: mMax
  INTEGER              :: nMax
  INTEGER              :: mnMax
  INTEGER, ALLOCATABLE :: mnMap(:,:)
  INTEGER              :: kMax

  LOGICAL(KIND=i8)              :: fsbc       
  INTEGER              :: ngrmx
  INTEGER              :: ncf
  INTEGER              :: ncf2
  CHARACTER(LEN=100)   :: path
!  CHARACTER(LEN=200)   :: fNameSnow
   CHARACTER(LEN=200)   :: fNameSSTWKL
!  CHARACTER(LEN=200)   :: fNameSSTAOI
!  CHARACTER(LEN=200)   :: fNameSoilms
!  CHARACTER(LEN=200)   :: fNameSibAlb
  INTEGER, ALLOCATABLE :: looku (:,:,:)
  REAL(KIND=r8),    ALLOCATABLE :: cnfac (:)
  REAL(KIND=r8),    ALLOCATABLE :: cnfac2(:)
  
  REAL(KIND=r8),    PARAMETER   :: undef =1.0e53_r8
  REAL(KIND=r8),    PARAMETER   :: eeemin=1.0e-35_r8
  REAL(KIND=r8),    PARAMETER   :: eeemax=1.0e35_r8
!
! unit number of input output file 
!
  INTEGER, PARAMETER   ::  nfin0   =18 ! input  file at time level t-dt
  INTEGER, PARAMETER   ::  nfin1   =18 ! input  file at time level t
  INTEGER, PARAMETER   ::  nfout0  =20 ! output file at time level t-dt
  INTEGER, PARAMETER   ::  nfout1  =21 ! output file at time level t
  INTEGER, PARAMETER   ::  nfout2  =21 ! output file at t=0 ( normal-mode initialized )
  INTEGER, PARAMETER   ::  nfclm0  =10 ! sst,soil moisture etc.  input
  INTEGER, PARAMETER   ::  nfclm1  =11 ! sst,soil moisture etc.  output
  INTEGER, PARAMETER   ::  nftgz0  =61 ! ground temperature and roughness length input
  INTEGER, PARAMETER   ::  nftgz1  =61 ! ground temperature and roughness length output
  INTEGER, PARAMETER   ::  nfsibt  =99 ! sib surface vegetation type
  INTEGER, PARAMETER   ::  nfsibd  =88 ! sib vegetation parameter
  INTEGER, PARAMETER   ::  nfsibi  =77 ! sib prognostic variable input  file
  INTEGER, PARAMETER   ::  nfsibo  =66 ! sib prognostic variable output file
  INTEGER, PARAMETER   ::  nfcnv0  = 0 ! initial information on convective clouds for int. radiation
  INTEGER, PARAMETER   ::  nfcnv1  =32 ! output information on convective clouds for int. radiation
  INTEGER, PARAMETER   ::  nfvar   =33 ! surface height variance
  INTEGER, PARAMETER   ::  nfsnw   =51 ! snow        file
  INTEGER, PARAMETER   ::  nfalb   =52 ! albedo file
  INTEGER, PARAMETER   ::  nfslm   =53 ! soil moisture file
  INTEGER, PARAMETER   ::  nfcldr  =74 ! cloud radiation fields
  INTEGER, PARAMETER   ::  nfnmi   =80 ! normal modes
  INTEGER, PARAMETER   ::  igrfu   =45
  INTEGER, PARAMETER   ::  iptu    =42
  INTEGER, PARAMETER   ::  nfdbh   =75 ! heating rate used for diabatic nlnmi


  INTEGER, ALLOCATABLE :: TMap(:)

CONTAINS



  ! InitInputOutput: Initializes module



  SUBROUTINE InitInputOutput(nfprt_in, nferr_in, ifprt_in, ngrmx_in, ncf_in , &
                             ncf2_in , mMax_in , nMax_in , mnMax_in,mnMap_in, &
                             kmax_in ,path_in  ,fNameSnow_in,fNameSSTWKL_in , &
                             fNameSSTAOI_in,fNameSoilms_in,fNameSibAlb_in   , &
                             fNameCnftBl   ,fNameCnf2Tb   ,fNameLookTb)
    INTEGER,          INTENT(IN) :: mMax_in
    INTEGER,          INTENT(IN) :: nMax_in
    INTEGER,          INTENT(IN) :: mnMax_in
    INTEGER,          INTENT(IN) :: mnMap_in(:,:)
    INTEGER,          INTENT(IN) :: kmax_in
    INTEGER,          INTENT(IN) :: nfprt_in
    INTEGER,          INTENT(IN) :: nferr_in
    INTEGER,          INTENT(IN) :: ifprt_in(100)
    INTEGER,          INTENT(IN) :: ngrmx_in
    INTEGER,          INTENT(IN) :: ncf_in
    INTEGER,          INTENT(IN) :: ncf2_in
    CHARACTER(LEN=*), INTENT(IN) :: path_in
    CHARACTER(LEN=*), INTENT(IN) :: fNameSnow_in
    CHARACTER(LEN=*), INTENT(IN) :: fNameSSTWKL_in
    CHARACTER(LEN=*), INTENT(IN) :: fNameSSTAOI_in
    CHARACTER(LEN=*), INTENT(IN) :: fNameSoilms_in
    CHARACTER(LEN=*), INTENT(IN) :: fNameSibAlb_in
    CHARACTER(LEN=*), INTENT(IN) :: fNameCnftBl
    CHARACTER(LEN=*), INTENT(IN) :: fNameCnf2Tb
    CHARACTER(LEN=*), INTENT(IN) :: fNameLookTb
    INTEGER :: mm
    INTEGER :: nn
    INTEGER :: l
    ALLOCATE (mnMap(mMax_in,nMax_in))
    
    OPEN(37, file=TRIM(fNameCnftBl),FORM="FORMATTED")
    OPEN(38, file=TRIM(fNameCnf2Tb),FORM="FORMATTED")
    OPEN(39, file=TRIM(fNameLookTb),FORM="FORMATTED")
    
    path  = path_in       
    nferr = nferr_in
    ifprt = ifprt_in
    mMax  = mMax_in 
    nMax  = nMax_in 
    mnMax = mnMax_in
    mnMap = mnMap_in
    kmax  = kmax_in
    !fNameSnow  =fNameSnow_in
    fNameSSTWKL=fNameSSTWKL_in
    !fNameSSTAOI=fNameSSTAOI_in
    !fNameSoilms=fNameSoilms_in
    !fNameSibAlb=fNameSibAlb_in
    ncf   = ncf_in
    ALLOCATE(cnfac(ncf))
    REWIND 37
    READ(37,"(5e16.8)") cnfac
    REWIND 37
    CLOSE(37,status='KEEP')


    ncf2  = ncf2_in
    ALLOCATE(cnfac2(ncf2))  
    REWIND 38  
    READ(38,"(5e16.8)") cnfac2
    REWIND 38  
    CLOSE(38,status='KEEP')


    ngrmx = ngrmx_in
    ALLOCATE(looku(0:9,0:9,0:ngrmx))
    REWIND 39
    READ(39,"(20i4)") looku
    REWIND 39
    CLOSE(39,status='KEEP')

    ALLOCATE(TMap(2*mnMax))
    l=0
    DO mm=1,mMax
       DO nn=mm,nMax
          l=l+1
          TMap(2*l-1)=2*mnMap(mm,nn)-1
          TMap(2*l  )=2*mnMap(mm,nn)
       END DO
    END DO
  END SUBROUTINE InitInputOutput



  ! cnvray: convert array



  SUBROUTINE cnvray (array, idim, ifr, ito)
    INTEGER, INTENT(IN   ) :: idim
    REAL(KIND=r8),    INTENT(INOUT) :: array(idim)
    INTEGER, INTENT(IN   ) :: ifr
    INTEGER, INTENT(IN   ) :: ito

    CHARACTER(LEN=20) :: c0
    CHARACTER(LEN=20) :: c1
    INTEGER           :: i
    INTEGER           :: icf
    INTEGER           :: igpf
    INTEGER           :: iuf
    INTEGER           :: igpt
    INTEGER           :: iut
    REAL(KIND=r8)              :: cf
    REAL(KIND=r8)              :: cf2

    CHARACTER(LEN=*), PARAMETER :: h="**(cnvray)**"

    ! consistency

    IF (ifr <= -1) THEN
       WRITE(c0,"(i20)") ifr
       WRITE(*,"(a)") h//" ERROR: ifr ("//TRIM(ADJUSTL(c0))//") <= -1 "
       STOP h
    ELSE IF (ito <= -1) THEN
       WRITE(c0,"(i20)") ito
       WRITE(*,"(a)") h//" ERROR: ito ("//TRIM(ADJUSTL(c0))//") <= -1 "
       STOP h
    ELSE IF (idim <= 0) THEN
       WRITE(c0,"(i20)") idim
       WRITE(*,"(a)") h//" ERROR: idim ("//TRIM(ADJUSTL(c0))//") <= 0 "
       STOP h
    ELSE IF (ito /= ifr) THEN
       igpf=ifr/10
       igpt=ito/10
       IF (igpf /= igpt) THEN
          WRITE(c0,"(i20)") igpf
          WRITE(c1,"(i20)") igpt
          WRITE(*,"(a)") h//" ERROR: igpf ("//TRIM(ADJUSTL(c0))//&
               &") /= igpt ("//TRIM(ADJUSTL(c1))//")"
          STOP h
       ELSE IF (igpf > ngrmx) THEN
          WRITE(c0,"(i20)") igpf
          WRITE(c1,"(i20)") ngrmx
          WRITE(*,"(a)") h//" ERROR: igpf ("//TRIM(ADJUSTL(c0))//&
               &") > ngrmx ("//TRIM(ADJUSTL(c1))//")"
          STOP h
       ELSE

          ! table look-up

          iuf=MOD(ifr,10)
          iut=MOD(ito,10)
          icf=looku(iuf,iut,igpf)

          ! consistency, again

          IF (icf < 1 .OR. icf > ncf) THEN
             WRITE(c0,"(i20)") icf
             WRITE(c1,"(i20)") ncf
             WRITE(*,"(a)") h//" ERROR: icf ("//TRIM(ADJUSTL(c0))//&
                  &") < 1 or > ncf ("//TRIM(ADJUSTL(c1))//")"
             STOP h
          END IF

          ! get coeficients

          cf=cnfac(icf)
          IF (icf <= ncf2) THEN
             cf2=cnfac2(icf)
          ELSE
             cf2=0.0_r8
          END IF

          ! convert array

          DO i = 1, idim
             IF (array(i) /= undef) THEN
                array(i)=cf*array(i)+cf2
             END IF
          END DO
       END IF
    END IF
  END SUBROUTINE cnvray


  ! wrfld: writes an unformatted 2-d field on file unit iudiag.


  SUBROUTINE wrfld(field,nwv,unit)
    !
    ! wrfld  :writes an unformatted 2-d field on file unit unit.
    !
    INTEGER     , INTENT(IN   ) :: nwv
    INTEGER     , INTENT(IN   ) :: unit
    REAL(KIND=r8)        , INTENT(IN   ) :: field(nwv)
    WRITE(unit)REAL(field,r4)
  END SUBROUTINE wrfld


  ! transp : after input, transposes arrays of spectral coefficients
  !          by swapping the order of the subscripts representing the
  !          degree and order of the associated legendre functions.



  SUBROUTINE transp(a,kmax,isign)
    INTEGER, INTENT(IN   ) :: kmax
    INTEGER, INTENT(IN   ) :: isign
    REAL(KIND=r8)   , INTENT(INOUT) :: a(2*mnMax,kmax)

    INTEGER                :: k         
    INTEGER                :: mn
    REAL(KIND=r8)                   :: w(2*mnMax)

    IF(isign.EQ.1) THEN
       DO k=1,kmax
          DO mn = 1, 2*mnMax
             w(mn) = a(TMap(mn),k)
          END DO
          DO mn = 1, 2*mnMax
             a(mn,k) = w(mn)
          END DO
       END DO
    ELSE
       DO k=1,kmax
          DO mn = 1, 2*mnMax
             w(TMap(mn)) = a(mn,k)
          END DO
          DO mn = 1, 2*mnMax
             a(mn,k) = w(mn)
          END DO
       END DO
    END IF
  END SUBROUTINE transp



  ! transw : writes out an array containing spectral coefficients
  !          representing one level of a global field after
  !          transposing the array by swapping the order of
  !          the subscripts for the degree and the order of
  !          the associated legendre functions.



  SUBROUTINE transw(a, unit)
    REAL(KIND=r8),    INTENT(IN) :: a(2*mnMax)
    INTEGER, INTENT(IN) :: unit

    INTEGER :: mn
    REAL(KIND=r8)    :: w(2*mnMax)

    DO mn = 1, 2*mnMax
       w(mn) = a(TMap(mn))
    END DO
    WRITE(unit) w
  END SUBROUTINE transw
  !
  ! scale, convert to 32 bits and output field
  !
  SUBROUTINE sclout(unit, field, nharm, levs, fact1, ivar, nufr, nuto)
    INTEGER, INTENT(IN) :: unit
    INTEGER, INTENT(IN) :: nharm
    INTEGER, INTENT(IN) :: levs
    REAL(KIND=r8),    INTENT(IN) :: field(nharm,levs)
    REAL(KIND=r8),    INTENT(IN) :: fact1
    CHARACTER(LEN=4), INTENT(IN) :: ivar
    INTEGER, INTENT(IN) :: nufr
    INTEGER, INTENT(IN) :: nuto
    REAL(KIND=r8) :: fldaux(nharm,1)

    INTEGER :: k, i

    DO k=1,levs
       DO i=1,nharm
          fldaux(i,1) = fact1 * field(i,k)
       END DO

       CALL cnvray(fldaux,nharm,nufr,nuto)

       IF(ivar.EQ.'SPEC') THEN
          CALL transp(fldaux, 1, 1)
       END IF
       !     
       !     check for values outside iee range ( +- 10**38, +- 10**-40 )
       !     
       DO i=1,nharm
          IF(ABS(fldaux(i,1)).LT.eeemin) &
               fldaux(i,1)=SIGN(eeemin,fldaux(i,1))
          IF(ABS(fldaux(i,1)).GT.eeemax)THEN
             fldaux(i,1)=SIGN(eeemax,fldaux(i,1))
          ENDIF
       END DO
       CALL WriteField(unit, fldaux(:,1))
    END DO
  END SUBROUTINE sclout

  !
  ! getsbc :read surface/atmosphere boundary conditions.
  !
  SUBROUTINE getsbc (imax,jmax,kmax,galb ,gsst ,gndvi,gslm,gsnw,gozo,tracermix,wsib3d,&
       ifday,tod,idate,idatec,&
       ifalb,ifsst,ifndvi,ifslm,ifslmSib2,ifsnw,ifozone,iftracer,&
       sstlag,intsst,intndvi,intsoilm,fint,tice,&
       yrl ,monl,ibMax,jbMax,ibMaxPerJB)
    IMPLICIT NONE
    !
    ! INPUT/OUTPUT VARIABLES
    !
    ! Real size of the grid
    INTEGER, INTENT(in   ) :: imax
    INTEGER, INTENT(in   ) :: jmax
    INTEGER, INTENT(in   ) :: kmax
    ! Size of block divided grid
    INTEGER, INTENT(in   ) :: ibMax
    INTEGER, INTENT(in   ) :: jbMax
    INTEGER, INTENT(in   ) :: ibMaxPerJB(jbMax)

    ! Boundary fields output
    REAL(KIND=r8), INTENT(out  ) :: galb(ibMax,jbMax) ! albedo
    REAL(KIND=r8), INTENT(out  ) :: gndvi(ibMax,jbMax) ! ndvi    
    REAL(KIND=r8), INTENT(out  ) :: gsst(ibMax,jbMax) ! sst
    REAL(KIND=r8), INTENT(out  ) :: gslm(ibMax,jbMax) ! soil moisture
    REAL(KIND=r8), INTENT(out  ) :: gsnw(ibMax,jbMax) ! snow
    REAL(KIND=r8), INTENT(out  ) :: wsib3d(ibMax,jbMax,3) ! moisture
    !hmjb o ozonio nao pode ser apenas 'out' pois, no caso de usar a antiga
    !  getoz(), ele sairia daqui com valores indefinidos... Com inout,
    !  ele entra e,  se nao for alterado, sai como entrou
    REAL(KIND=r8), INTENT(inout) :: gozo(ibMax,kMax,jbMax) ! ozone
    REAL(KIND=r8), INTENT(inout) :: tracermix(ibMax,kMax,jbMax) ! tracer

    ! Options for reading boundary fields
    INTEGER, INTENT(inout) :: ifalb
    INTEGER, INTENT(inout) :: ifsst
    INTEGER, INTENT(inout) :: ifndvi
    INTEGER, INTENT(inout) :: ifslm
    INTEGER, INTENT(inout) :: ifsnw
    INTEGER, INTENT(inout) :: ifslmSib2
    INTEGER, INTENT(inout) :: ifozone
    INTEGER, INTENT(inout) :: iftracer
    ! Time
    INTEGER, INTENT(in   ) :: ifday
    REAL(KIND=r8), INTENT(in   ) :: tod
    INTEGER, INTENT(in   ) :: idate(4)
    INTEGER, INTENT(in   ) :: idatec(4)
    REAL(KIND=r8), INTENT(in   ) :: sstlag
    INTEGER, INTENT(in   ) :: intsst
    INTEGER, INTENT(in   ) :: intndvi
    INTEGER, INTENT(in   ) :: intsoilm
    REAL(KIND=r8), INTENT(in   ) :: fint
    REAL(KIND=r8), INTENT(in   ) :: tice
    REAL(KIND=r8), INTENT(in   ) :: yrl
    INTEGER, INTENT(in   ) :: monl(12)
    !
    ! LOCAL VARIABLES
    !
    REAL(KIND=r8)                :: xndvi   (ibMax,jbMax)
    REAL(KIND=r8)                :: xsst    (ibMax,jbMax)
    REAL(KIND=r8)                :: bfr_in  (imax,jmax)
    REAL(KIND=r8)                :: bfrw_in  (imax,jmax,3)
    REAL(KIND=r8)                :: bfrw_out  (ibmax,jbmax,3)
    REAL(KIND=r4)                :: rbrfw3d    (iMax,jMax,3)

    REAL(KIND=r8)                :: bfr_in3 (imax,kmax,jmax)
    REAL(KIND=r4)                :: bfr_3r4 (imax,kmax,jmax)
    REAL(KIND=r8)                :: bfr_out (ibMax,jbMax)
    REAL(KIND=r8)                :: bfr_out3(ibMax,kmax,jbMax)
    REAL(KIND=r4)                :: rbrf    (iMax,jMax)
    REAL(KIND=r4)                :: rbrf3   (iMax,kmax,jMax)

    !
    !
    INTEGER                :: lrecl,LRecIn
    REAL(KIND=r8)          :: fhr
    INTEGER                :: mf
    INTEGER                :: mn
    INTEGER                :: mf_ndvi
    INTEGER                :: mn_ndvi

    INTEGER                :: month
    INTEGER                :: mm
    INTEGER                :: i
    INTEGER                :: j
    INTEGER                :: k
    INTEGER                :: irec
    INTEGER                :: irec_ndvi

    REAL(KIND=r8)                :: f1
    REAL(KIND=r8)                :: f2
    REAL(KIND=r8)                :: f1_ndvi
    REAL(KIND=r8)                :: f2_ndvi
    REAL(KIND=r8)                :: gmax
    REAL(KIND=r8)                :: gmin
    REAL(KIND=r8)                :: fsst
    REAL(KIND=r8)                :: fndvi
    REAL(KIND=r8)                :: fisst
    REAL(KIND=r8)                :: findvi
    REAL(KIND=r8)                :: xx1
    REAL(KIND=r8)                :: xx2
    REAL(KIND=r8)                :: xday
    INTEGER :: ierr
    !
    !   ifxxx=0    xxx is not processed
    !   ifxxx=1    xxx is set to month=idate(2) in the first call,
    !              but not processed from the subsequent calls.
    !              ifxxx is set to zero after interpolation
    !   ifxxx=2    xxx is interpolated to current day and time every fint
    !              hours synchronized to 00z regardless of initial time.
    !              interpolation is continuous (every time step) if fint<0.
    !   ifxxx=3    xxx is interpolated to current day and time when ifday=0
    !              and tod=0.0 but not processed otherwise
    !              ( appropriate only when xxx is predicted )
    !
    !              the following are for sst only (fint applies as in
    !              ifxxx=2):
    !   ifsst=4    sst is linearly interpolated from continuous direct
    !              access data set to current day and time.  data set
    !              is assumed to be spaced every intsst days or every
    !              calendar month is intsst < 0.
    !   ifsst=5    sst is expanded from piecewise cubic coefficients in
    !              direct access data set to current day and time.  data set
    !              is assumed to be spaced every intsst days.
    !   note:      for ifsst=4 or 5 sstlag must be set.  sstlag is the
    !              number of days plus any fraction prior to the initial
    !              condition date and time the data set begins if intsst > 0
    !              sstlag is the number of whole months prior to the initial
    !              condition date the data set begins if intsst < 0.
    !
    !     ifsst=-1 for numerical weather forecasting using mean weekly sst:
    !              sst is read in the first call as the second record of
    !              the archieve but not processed from the subsequent calls.
    !
    rbrf=0.0_r4
    rbrfw3d=0.0_r4
    rbrf3=0.0_r4
    IF (fint > 0.0_r8) THEN
       fhr=REAL(idate(1),r8)+tod/3600.0_r8+1.0e-3_r8
       IF (.NOT. fsbc .AND. ABS( MOD(fhr,fint)) > 1.0e-2_r8) THEN
          RETURN
       END IF
    END IF
    IF (ifsst == 4 .AND. intsst <= 0) THEN
       CALL GetRecWgtMonthlySST &
            (idate, idatec, tod, labelsi, labelsj, &
            irec, f1, f2, mf, mn,monl)

!!$       WRITE (UNIT=nfprt, FMT='(A)') ' GetRecWgtMonthlySST'
!!$       WRITE (UNIT=nfprt, FMT='(/,4(A,I5),/)') &
!!$            ' reci = ', irec, ' recf = ', irec+1, &
!!$            ' mra = ', mf, ' mrb = ', mf+1
!!$       WRITE (UNIT=nfprt, FMT=*) ' fa  (*mra) = ', f1, ' fb  (*mrb) = ', f2
    ELSE
       CALL GetWeightsOld(yrl,monl,idatec, tod, f1, f2,mf)
!!$       WRITE (UNIT=nfprt, FMT=*) ' fa  (*mra) = ', f1, ' fb  (*mrb) = ', f2
    END IF


    IF (ifndvi == 4 .AND. intndvi <= 0) THEN
       CALL GetRecWgtMonthlySST &
            (idate, idatec, tod, labelsi, labelsj, &
            irec_ndvi, f1_ndvi, f2_ndvi, mf_ndvi, mn_ndvi,monl)

!!$       WRITE (UNIT=nfprt, FMT='(A)') ' GetRecWgtMonthlySST'
!!$       WRITE (UNIT=nfprt, FMT='(/,4(A,I5),/)') &
!!$            ' reci = ', irec, ' recf = ', irec+1, &
!!$            ' mra = ', mf, ' mrb = ', mf+1
!!$       WRITE (UNIT=nfprt, FMT=*) ' fa  (*mra) = ', f1, ' fb  (*mrb) = ', f2
    ELSE
       CALL GetWeightsOld(yrl,monl,idatec, tod, f1_ndvi, f2_ndvi,mf_ndvi)
!!$       WRITE (UNIT=nfprt, FMT=*) ' fa  (*mra) = ', f1, ' fb  (*mrb) = ', f2
    END IF

    !
    ! process albedo file
    !
    IF (ifalb /= 0) THEN
       IF (ifalb == 1) THEN
          month=idate(2)
          INQUIRE (IOLENGTH=LRecIn) rbrf
          OPEN (UNIT=nfalb,FILE=TRIM(fNameAlbedo),FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIn, &
               ACTION='read', STATUS='old', IOSTAT=ierr)
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameAlbedo), ierr
             STOP "**(ERROR)**"
          END IF
          irec=month
          CALL ReadGetALB(nfalb,irec,bfr_in)

          IF (reducedGrid) THEN
             CALL AveBoxIJtoIBJB(bfr_in,galb)
          ELSE
             CALL IJtoIBJB(bfr_in ,galb)
          END IF
          CLOSE(UNIT=nfalb)
          ifalb=0
       ELSE IF (&
            (ifalb == 2) .OR. &
            (ifalb == 3 .AND. tod == 0.0_r8 .AND. ifday == 0)) THEN
          rbrf=0.0_r4
          INQUIRE (IOLENGTH=LRecIn) rbrf
          OPEN (UNIT=nfalb,FILE=TRIM(fNameAlbedo),FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIn, &
               ACTION='read', STATUS='old', IOSTAT=ierr)
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameAlbedo), ierr
             STOP "**(ERROR)**"
          END IF
          irec=mf
          CALL ReadGetALB(nfalb,irec,bfr_in)
          IF (reducedGrid) THEN
             CALL AveBoxIJtoIBJB(bfr_in,galb)
          ELSE
             CALL IJtoIBJB(bfr_in ,galb)
          END IF
          IF (irec == 12) THEN
             irec=1
          ELSE   
             irec=irec+1
          END IF
          CALL ReadGetALB(nfalb,irec,bfr_in)
          IF (reducedGrid) THEN
             CALL AveBoxIJtoIBJB(bfr_in,bfr_out)
          ELSE
             CALL IJtoIBJB(bfr_in ,bfr_out)
          END IF
          CLOSE(UNIT=nfalb)
          gmax=-1.0e10_r8
          gmin=+1.0e10_r8

          DO j=1,jbMax
             DO i=1,ibMaxPerJB(j)
                galb(i,j)=f2*galb(i,j)+f1*bfr_out(i,j)
                gmax=MAX(gmax,galb(i,j))
                gmin=MIN(gmin,galb(i,j))
             END DO
          END DO

          IF (ifalb == 3 .AND. tod == 0.0_r8 .AND. ifday == 0) THEN
             ifalb=0
          END IF

          IF (nfctrl(23) >= 1) THEN
             WRITE(UNIT=nfprt,FMT=888) mf,f1,f2,gmax,gmin
          END IF

       ELSE
          WRITE(UNIT=nfprt,FMT=999)
          STOP
       END IF
    END IF

    !
    ! process sst file
    !

    IF (ifsst /= 0) THEN
       IF (ifsst == -1) THEN
          INQUIRE (IOLENGTH=LRecIn) rbrf
          OPEN (UNIT=nfsst, FILE=TRIM(fNameSSTAOI),FORM='FORMATTED', ACCESS='SEQUENTIAL', &
               ACTION='READ',STATUS='OLD', IOSTAT=ierr)
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameSSTAOI), ierr
             STOP "**(ERROR)**"
          END IF
          irec=1
          CALL ReadGetSST(nfsst,irec,bfr_in,'txt')
          irec=2
          CALL ReadGetSST(nfsst,irec,bfr_in,'txt')
          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(bfr_in ,gsst)
          ELSE
             CALL IJtoIBJB(bfr_in ,gsst)
          END IF
          CLOSE(UNIT=nfsst)
          PRINT*,gsst
          gmax=-1.0e10_r8
          gmin=+1.0e10_r8

          DO j=1,jbMax
             DO i=1,ibMaxPerJB(j)
                IF (gsst(i,j) > 10.0_r8) THEN
                   gsst(i,j)=-gsst(i,j)
                ELSE IF (gsst(i,j) < 0.0_r8) THEN
                   gsst(i,j)=290.0_r8
                ELSE
                   PRINT *, " OPTION ifsst=-1 INCORRECT VALUE OF SST "
                   STOP "**(ERROR)**"
                END IF
                gmax=MAX(gmax,gsst(i,j))
                gmin=MIN(gmin,gsst(i,j))
             END DO
          END DO
          WRITE(UNIT=nfprt,FMT=667) ifsst,gmax,gmin
          ifsst=0
       ELSE IF (ifsst == 1) THEN
          OPEN(UNIT=nfsst, FILE=TRIM(fNameSSTAOI), FORM='formatted', ACCESS='sequential',&
               ACTION='read', STATUS='old', IOSTAT=ierr)
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameSSTAOI), ierr
             STOP "**(ERROR)**"
          END IF
          READ(UNIT=nfsst)
          month=idate(2)
          DO mm=1,month
             CALL ReadGetSST(nfsst,irec,bfr_in,'txt')
          END DO
          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(bfr_in ,gsst)
          ELSE
             CALL IJtoIBJB(bfr_in ,gsst)
          END IF
          CLOSE(UNIT=nfsst)
          DO j=1,jbMax
             DO i=1,ibMaxPerJB(j)
                IF (gsst(i,j) > 10.0_r8) THEN
                   gsst(i,j)=gsst(i,j)
                ELSE IF (gsst(i,j) < 0.0_r8) THEN
                   gsst(i,j)=290.0_r8
                ELSE
                   PRINT *, " OPTION ifsst=-1 INCORRECT VALUE OF SST "
                   STOP "**(ERROR)**"
                END IF
             END DO
          END DO
          ifsst=0
       ELSE IF (ifsst == 2.OR. &
            (ifsst == 3.AND.tod == 0.0_r8.AND.ifday == 0)) THEN
          INQUIRE (IOLENGTH=LRecIn) rbrf
          OPEN (UNIT=nfsst, FILE=TRIM(fNameSSTAOI),FORM='FORMATTED', ACCESS='SEQUENTIAL', &
             ACTION='READ',STATUS='OLD', IOSTAT=ierr)
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameSSTAOI), ierr
             STOP "**(ERROR)**"
          END IF

          irec = mf+1
          CALL ReadGetSST(nfsst,irec,bfr_in,'txt')

          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(bfr_in ,xsst)
          ELSE
             CALL IJtoIBJB(bfr_in ,xsst)
          END IF
          
          IF (irec == 13) THEN
             irec=2
          ELSE
             irec=irec+1
          END IF
          CALL ReadGetSST(nfsst,irec,bfr_in,'txt')
          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(bfr_in ,bfr_out)
          ELSE
             CALL IJtoIBJB(bfr_in ,bfr_out)
          END IF

          CLOSE(UNIT=nfsst)
          gmax=-1.0e10_r8
          gmin=+1.0e10_r8
          DO j=1,jbMax
             DO i=1,ibMaxPerJB(j)

                fsst=f2*xsst(i,j)+f1*bfr_out(i,j)
                IF (fsst > gmax) THEN
                   gmax=fsst
                END IF
                IF (fsst < gmin) THEN
                   gmin=fsst
                END IF
                IF (fsst > 10.0_r8) THEN
                   xsst(i,j)=-fsst
                ELSE IF (fsst < 0.0_r8) THEN
                   xsst(i,j)=290.0_r8
                ELSE
                   PRINT *, " OPTION ifsst=-1 INCORRECT VALUE OF SST "
                   STOP "**(ERROR)**"
                END IF
             END DO
          END DO
 
          DO j=1,jbMax
             DO i=1,ibMaxPerJB(j)
                IF (tod == 0.0_r8.AND.ifday == 0) THEN
                   gsst(i,j)=xsst(i,j)
                ELSE IF (xsst(i,j) < 0.0_r8.AND.ABS(xsst(i,j)) >= tice) THEN
                   gsst(i,j)=xsst(i,j)
                ELSE IF (xsst(i,j) < 0.0_r8.AND.ABS(gsst(i,j)) >= tice) THEN
                   gsst(i,j)=-tice+1.0e-2_r8
                END IF
             END DO
          END DO

          IF (ifsst == 3.AND.tod == 0.0_r8.AND.ifday == 0) THEN
             ifsst=0
          END IF
          IF (nfctrl(23) >= 1) THEN
             WRITE(UNIT=nfprt,FMT=666) mf,f1,f2,gmax,gmin
          END IF
          PRINT*,gsst
       ELSE IF (ifsst == 4) THEN
          IF (intsst > 0) THEN
             fisst=REAL(intsst,r8)
             xday=ifday+tod/86400.0_r8+sstlag
             irec=INT(xday/fisst+1.0e-3_r8+1.0_r8)
             xx1= MOD(xday,fisst)/fisst
             xx2=1.0_r8-xx1
          ELSE
             xx1=f1
             xx2=f2
          END IF
          INQUIRE (IOLENGTH=lrecl) rbrf
          !lrecl=lrecl/2
          OPEN(UNIT=nfsst,FILE=TRIM(fNameSSTAOI),FORM='unformatted',ACCESS='direct',&
               RECL=lrecl,ACTION='read', STATUS='old', IOSTAT=ierr)
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameSSTAOI), ierr
             STOP "**(ERROR)**"
          END IF
          CALL ReadGetSST2(nfsst,bfr_in,irec)
          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(bfr_in ,xsst)
          ELSE
             CALL IJtoIBJB(bfr_in ,xsst)
          END IF
          CALL ReadGetSST2(nfsst,bfr_in,irec+1)
          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(bfr_in ,bfr_out)
          ELSE
             CALL IJtoIBJB(bfr_in ,bfr_out)
          END IF
          CLOSE(UNIT=nfsst)
          gmax=-1.0e10_r8
          gmin=+1.0e10_r8
          DO j=1,jbMax
             DO i=1,ibMaxPerJB(j)
                fsst=xx2*xsst(i,j)+xx1*bfr_out(i,j)
                IF (fsst > gmax) THEN
                   gmax=fsst
                END IF
                IF (fsst < gmin) THEN
                   gmin=fsst
                END IF
                IF (fsst > 10.0_r8) THEN
                   xsst(i,j)=-fsst
                ELSE IF (fsst < 0.0_r8) THEN
                   xsst(i,j)=290.0_r8
                ELSE
                   PRINT *, " OPTION ifsst=-1 INCORRECT VALUE OF SST "
                   STOP "**(ERROR)**"
                END IF
             END DO
          END DO

          IF (nfctrl(23) >= 1) THEN
             WRITE(UNIT=nfprt,FMT=666) irec,xx1,xx2,gmax,gmin
          END IF

          DO j=1,jbMax
             DO i=1,ibMaxPerJB(j)
                IF (tod == 0.0_r8.AND.ifday == 0) THEN
                   gsst(i,j)=xsst(i,j)
                ELSE IF (xsst(i,j) < 0.0_r8.AND.ABS(xsst(i,j)) >= tice) THEN
                   gsst(i,j)=xsst(i,j)
                ELSE IF (xsst(i,j) < 0.0_r8.AND.ABS(gsst(i,j)) >= tice) THEN
                   gsst(i,j)=-tice+1.0e-2_r8
                END IF
             END DO
          END DO

       ELSE IF (ifsst == 5.AND.intsst > 0) THEN

          !*(JP)* Eliminei este caso pelas obs do Bonatti e minhas

          PRINT *, " OPTION ifsst=5 NOT CORRECTLY IMPLEMENTED "
          STOP "**(ERROR)**"

       ELSE
          WRITE(UNIT=nfprt,FMT=1999)
          STOP
       END IF
    END IF


    IF (schemes ==2) THEN

       !
       ! process ndvi file
       !

       IF (ifndvi /= 0) THEN
          IF (ifndvi == -1) THEN
             INQUIRE (IOLENGTH=LRecIn) rbrf
             OPEN (UNIT=nfndvi, FILE=TRIM(fNameNDVIAOI),FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIn,&
             ACTION='READ',STATUS='OLD', IOSTAT=ierr)
             IF (ierr /= 0) THEN
                WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                     TRIM(fNameNDVIAOI), ierr
                STOP "**(ERROR)**"
             END IF
             irec_ndvi=1
             CALL ReadGetSST(nfndvi,irec_ndvi,bfr_in)
             IF (reducedGrid) THEN
                CALL NearestIJtoIBJB(bfr_in ,gndvi)
             ELSE
                CALL IJtoIBJB(bfr_in ,gndvi)
             END IF
             CLOSE(UNIT=nfndvi)
             gmax=0.0e0_r8
             gmin=+1.0e0_r8

             DO j=1,jbMax
                DO i=1,ibMaxPerJB(j)
                   IF (gndvi(i,j) > 0.0_r8) THEN
                      gndvi(i,j)=gndvi(i,j)
                   ELSE IF (gndvi(i,j) <=0.0_r8) THEN
                      gndvi(i,j)=0.0_r8
                   ELSE
                      PRINT *, " OPTION ifndvi=-1 INCORRECT VALUE OF SST "
                      STOP "**(ERROR)**"
                   END IF
                   gmax=MAX(gmax,gndvi(i,j))
                   gmin=MIN(gmin,gndvi(i,j))
                END DO
             END DO

             WRITE(UNIT=nfprt,FMT=667) ifndvi,gmax,gmin
             ifndvi=0
          ELSE IF (ifndvi == 1) THEN
             OPEN(UNIT=nfndvi, FILE=TRIM(fNameNDVIAOI), FORM='unformatted', ACCESS='sequential',&
                  ACTION='read', STATUS='old', IOSTAT=ierr)
             IF (ierr /= 0) THEN
                WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                     TRIM(fNameNDVIAOI), ierr
                STOP "**(ERROR)**"
             END IF
             month=idate(2)
             DO mm=1,month
                CALL ReadGetSST(nfndvi,irec_ndvi,bfr_in)
             END DO
             IF (reducedGrid) THEN
                CALL NearestIJtoIBJB(bfr_in ,gndvi)
             ELSE
                CALL IJtoIBJB(bfr_in ,gndvi)
             END IF
             CLOSE(UNIT=nfndvi)
             DO j=1,jbMax
                DO i=1,ibMaxPerJB(j)
                   IF (gndvi(i,j) > 0.0_r8) THEN
                         gndvi(i,j)=gndvi(i,j)
                   ELSE IF (gndvi(i,j) <=0.0_r8) THEN
                      gndvi(i,j)=0.0_r8
                   ELSE
                      PRINT *, " OPTION ifndvi=-1 INCORRECT VALUE OF SST "
                      STOP "**(ERROR)**"
                   END IF
                END DO
             END DO
             ifndvi=0
          ELSE IF (ifndvi == 2.OR. &
               (ifndvi == 3.AND.tod == 0.0_r8.AND.ifday == 0)) THEN
             INQUIRE (IOLENGTH=LRecIn) rbrf
             OPEN (UNIT=nfndvi, FILE=TRIM(fNameNDVIAOI),FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIn,&
                   ACTION='READ',STATUS='OLD', IOSTAT=ierr)
             IF (ierr /= 0) THEN
                WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                     TRIM(fNameNDVIAOI), ierr
                STOP "**(ERROR)**"
             END IF

             irec_ndvi = mf
             CALL ReadGetSST(nfndvi,irec_ndvi,bfr_in)

             IF (reducedGrid) THEN
                CALL NearestIJtoIBJB(bfr_in ,xndvi)
             ELSE
                CALL IJtoIBJB(bfr_in ,xndvi)
             END IF
          
             IF (irec_ndvi == 12) THEN
                irec_ndvi=1
             ELSE
                irec_ndvi=irec_ndvi+1
             END IF
             CALL ReadGetSST(nfndvi,irec_ndvi,bfr_in)
             IF (reducedGrid) THEN
                CALL NearestIJtoIBJB(bfr_in ,bfr_out)
             ELSE
                CALL IJtoIBJB(bfr_in ,bfr_out)
             END IF

             CLOSE(UNIT=nfndvi)
             gmax=0.0e0_r8
             gmin=+1.0e0_r8
             DO j=1,jbMax
                DO i=1,ibMaxPerJB(j)

                   fndvi=f2_ndvi*xndvi(i,j)+f1_ndvi*bfr_out(i,j)
                   IF (fndvi > gmax) THEN
                      gmax=fndvi
                   END IF
                   IF (fndvi < gmin) THEN
                      gmin=fndvi
                   END IF
                   IF (fndvi > 0.0_r8) THEN
                      xndvi(i,j)=fndvi
                   ELSE IF (fndvi <= 0.0_r8) THEN
                      xndvi(i,j)=0.0_r8
                   ELSE
                      PRINT *, " OPTION ifndvi=-1 INCORRECT VALUE OF SST "
                      STOP "**(ERROR)**"
                   END IF
                END DO
             END DO
             DO j=1,jbMax
                DO i=1,ibMaxPerJB(j)
                   IF (tod == 0.0_r8.AND.ifday == 0) THEN
                      gndvi(i,j)=xndvi(i,j)
                   ELSE
                      gndvi(i,j)=xndvi(i,j)
                   END IF
                END DO
             END DO

             IF (ifndvi == 3.AND.tod == 0.0_r8.AND.ifday == 0) THEN
                ifndvi=0
             END IF
             IF (nfctrl(23) >= 1) THEN
                WRITE(UNIT=nfprt,FMT=666) mf,f1_ndvi,f2_ndvi,gmax,gmin
             END IF
          ELSE IF (ifndvi == 4) THEN
             IF (intndvi > 0) THEN
                findvi=REAL(intndvi,r8)
                xday=ifday+tod/86400.0_r8+sstlag
                irec_ndvi=INT(xday/findvi+1.0e-3_r8+1.0_r8)
                xx1= MOD(xday,findvi)/findvi
                xx2=1.0_r8-xx1
             ELSE
                xx1=f1_ndvi
                xx2=f2_ndvi
             END IF
             INQUIRE (IOLENGTH=lrecl) bfr_in
             lrecl=lrecl/2
             OPEN(UNIT=nfndvi,FILE=TRIM(fNameNDVIAOI),FORM='unformatted',ACCESS='direct',&
                  RECL=lrecl,ACTION='read', STATUS='old', IOSTAT=ierr)
             IF (ierr /= 0) THEN
                WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                     TRIM(fNameNDVIAOI), ierr
                STOP "**(ERROR)**"
             END IF
             CALL ReadGetSST2(nfndvi,bfr_in,irec_ndvi)
             IF (reducedGrid) THEN
                CALL NearestIJtoIBJB(bfr_in ,xndvi)
             ELSE
                CALL IJtoIBJB(bfr_in ,xndvi)
             END IF
             CALL ReadGetSST2(nfndvi,bfr_in,irec_ndvi+1)
             IF (reducedGrid) THEN
                CALL NearestIJtoIBJB(bfr_in ,bfr_out)
             ELSE
                CALL IJtoIBJB(bfr_in ,bfr_out)
             END IF
             CLOSE(UNIT=nfndvi)
             gmax=0.0e0_r8
             gmin=+1.0e0_r8
             DO j=1,jbMax
                DO i=1,ibMaxPerJB(j)
                   fndvi=xx2*xndvi(i,j)+xx1*bfr_out(i,j)
                   IF (fndvi > gmax) THEN
                      gmax=fndvi
                   END IF
                   IF (fndvi < gmin) THEN
                         gmin=fndvi
                   END IF
                   IF (fndvi > 0.0_r8) THEN
                      xndvi(i,j)=fndvi
                   ELSE IF (fndvi <= 0.0_r8) THEN
                      xndvi(i,j)=0.0_r8
                   ELSE
                      PRINT *, " OPTION ifndvi=-1 INCORRECT VALUE OF SST "
                      STOP "**(ERROR)**"
                   END IF
                END DO
             END DO

             IF (nfctrl(23) >= 1) THEN
                WRITE(UNIT=nfprt,FMT=666) irec_ndvi,xx1,xx2,gmax,gmin
             END IF

             DO j=1,jbMax
                DO i=1,ibMaxPerJB(j)
                   IF (tod == 0.0_r8.AND.ifday == 0) THEN
                      gndvi(i,j)=xndvi(i,j)
                   ELSE
                      gndvi(i,j)=xndvi(i,j)
                   END IF
                END DO
             END DO

          ELSE IF (ifndvi == 5.AND.intndvi > 0) THEN

             !*(JP)* Eliminei este caso pelas obs do Bonatti e minhas

             PRINT *, " OPTION ifndvi=5 NOT CORRECTLY IMPLEMENTED "
             STOP "**(ERROR)**"

          ELSE
             WRITE(UNIT=nfprt,FMT=1999)
             STOP
          END IF
       END IF
    END IF
    !
    ! process snow file
    !

    IF (ifsnw /= 0) THEN
       IF (ifsnw == 1) THEN
          INQUIRE (IOLENGTH=LRecIn) rbrf
          OPEN (UNIT=nfsnw,FILE=TRIM(fNameSnow), FORM='FORMATTED', ACCESS='SEQUENTIAL', &
                ACTION='READ',STATUS='OLD', IOSTAT=ierr)
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameSnow), ierr
             STOP "**(ERROR)**"
          END IF
          irec=1
          CALL ReadGetSNW(nfsnw,irec,bfr_in,'txt')
          IF (reducedGrid) THEN
             CALL AveBoxIJtoIBJB(bfr_in,gsnw)
          ELSE
             CALL IJtoIBJB(bfr_in,gsnw)
          END IF
          CLOSE(UNIT=nfsnw)
          ifsnw=0
       ELSE IF (ifsnw == 2.OR. &
            (ifsnw == 3.AND.tod == 0.0_r8.AND.ifday == 0)) THEN
          INQUIRE (IOLENGTH=LRecIn) rbrf
          OPEN (UNIT=nfsnw,FILE=TRIM(fNameSnow), FORM='FORMATTED', ACCESS='SEQUENTIAL', &
                ACTION='READ',STATUS='OLD', IOSTAT=ierr)
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameSnow), ierr
             STOP "**(ERROR)**"
          END IF
          irec=1
          CALL ReadGetSNW(nfsnw,irec,bfr_in,'txt')
          IF (reducedGrid) THEN
             CALL AveBoxIJtoIBJB(bfr_in,gsnw)
          ELSE
             CALL IJtoIBJB(bfr_in,gsnw)
          END IF
          CLOSE(UNIT=nfsnw)
          gmax=-1.0e10_r8
          gmin=+1.0e10_r8
          DO j=1,jbMax
             DO i=1,ibMaxPerJB(j)
                gmax=MAX(gmax,gsnw(i,j))
                gmin=MIN(gmin,gsnw(i,j))
             END DO
          END DO

          IF (ifsnw == 3.AND.tod == 0.0_r8.AND.ifday == 0) THEN
             ifsnw=0
          END IF
          IF (nfctrl(23) >= 1) THEN
             WRITE(UNIT=nfprt,FMT=444) gmax,gmin
          END IF
       ELSE
          WRITE(UNIT=nfprt,FMT=555)
          STOP
       END IF
    END IF

    IF (schemes ==5) THEN
       !
       ! process soil moisture file
       !

       IF (ifslmSib2 /= 0) THEN
          IF (ifslmSib2 == 1) THEN
             rbrfw3d=0.0_r4
             INQUIRE (IOLENGTH=LRecIn) rbrfw3d
             OPEN (UNIT=nfSoilMostSib2,FILE=TRIM(fNameSoilMoistSib2),FORM='UNFORMATTED', ACCESS='DIRECT', &
                  ACTION='read', RECL=LRecIn, STATUS='OLD', IOSTAT=ierr) 
             IF (ierr /= 0) THEN
                WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                     TRIM(fNameSoilMoistSib2), ierr
                STOP "**(ERROR)**"
             END IF
          
             irec=idate(2)
             CALL ReadGetSLM3D(nfSoilMostSib2,irec,bfrw_in)
     
             DO k=1,3
                IF (reducedGrid) THEN
                   CALL AveBoxIJtoIBJB(bfrw_in(:,:,k),wsib3d(:,:,k))
                ELSE
                   CALL IJtoIBJB(bfrw_in(:,:,k),wsib3d(:,:,k))
                END IF
             END DO
     
             CLOSE(UNIT=nfSoilMostSib2)
             ifslmSib2=0
          ELSE IF (ifslmSib2 == 2.OR. &
               (ifslmSib2 == 3.AND.tod == 0.0_r8.AND.ifday == 0)) THEN
             rbrfw3d=0.0_r4
             INQUIRE (IOLENGTH=LRecIn) rbrfw3d
             OPEN (UNIT=nfSoilMostSib2,FILE=TRIM(fNameSoilMoistSib2),FORM='UNFORMATTED', ACCESS='DIRECT', &
                  ACTION='read', RECL=LRecIn, STATUS='OLD', IOSTAT=ierr) 
             IF (ierr /= 0) THEN
                WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                     TRIM(fNameSoilMoistSib2), ierr
                STOP "**(ERROR)**"
             END IF
             irec=mf
             CALL ReadGetSLM3D(nfSoilMostSib2,irec,bfrw_in)
 
             DO k=1,3
                IF (reducedGrid) THEN
                   CALL AveBoxIJtoIBJB(bfrw_in(:,:,k),wsib3d(:,:,k))
                ELSE
                   CALL IJtoIBJB(bfrw_in(:,:,k),wsib3d(:,:,k))
                END IF
             END DO

             IF (irec == 12) THEN
                irec=1
             ELSE
                irec=irec+1    
             END IF
  
             CALL ReadGetSLM3D(nfSoilMostSib2,irec,bfrw_in)

             DO k=1,3
                IF (reducedGrid) THEN
                   CALL AveBoxIJtoIBJB(bfrw_in(:,:,k),bfrw_out(:,:,k))
                ELSE
                   CALL IJtoIBJB(bfrw_in(:,:,k),bfrw_out(:,:,k))
                END IF
             END DO
             CLOSE(UNIT=nfSoilMostSib2)
             gmax=-1.0e0_r8
             gmin=+1.0e0_r8
             DO k=1,3
                DO j=1,jbMax
                   DO i=1,ibMaxPerJB(j)
                      wsib3d(i,j,k)=f2*wsib3d(i,j,k)+f1*bfrw_out(i,j,k)
                      gmax=MAX(gmax,wsib3d(i,j,k))
                      gmin=MIN(gmin,wsib3d(i,j,k))
                   END DO
                END DO
             END DO
             IF (ifslmSib2 == 3.AND.tod == 0.0_r8.AND.ifday == 0) THEN
                ifslmSib2=0
             END IF
             IF (nfctrl(23) >= 1) THEN
                WRITE(UNIT=nfprt,FMT=222) mf,f1,f2,gmax,gmin
             END IF
          ELSE
             WRITE(UNIT=nfprt,FMT=333)
             STOP
          END IF
       END IF
       
    ELSE
       !
       ! process soil moisture file
       !
       IF (ifslm /= 0) THEN
          IF (ifslm == 1) THEN
             INQUIRE (IOLENGTH=LRecIn) rbrf
             OPEN (UNIT=nfslm,FILE=TRIM(fNameSoilms),FORM='FORMATTED', ACCESS='SEQUENTIAL', &
                  ACTION='read', STATUS='OLD', IOSTAT=ierr) 
             IF (ierr /= 0) THEN
                   WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameSoilms), ierr
                STOP "**(ERROR)**"
             END IF
             
             irec=idate(2)
             CALL ReadGetSLM(nfslm,irec,bfr_in,'txt')
             IF (reducedGrid) THEN
                CALL AveBoxIJtoIBJB(bfr_in,gslm)
             ELSE
                CALL IJtoIBJB(bfr_in,gslm)
             END IF
             CLOSE(UNIT=nfslm)
             ifslm=0
          ELSE IF (ifslm == 2.OR. &
               (ifslm == 3.AND.tod == 0.0_r8.AND.ifday == 0)) THEN
             INQUIRE (IOLENGTH=LRecIn) rbrf
             OPEN (UNIT=nfslm,FILE=TRIM(fNameSoilms),FORM='FORMATTED', ACCESS='SEQUENTIAL', &
                  ACTION='read', STATUS='OLD', IOSTAT=ierr) 
             IF (ierr /= 0) THEN
                WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                     TRIM(fNameSoilms), ierr
                STOP "**(ERROR)**"
             END IF
             irec=mf
             CALL ReadGetSLM(nfslm,irec,bfr_in,'txt')
 
             IF (reducedGrid) THEN
                CALL AveBoxIJtoIBJB(bfr_in,gslm)
             ELSE
                CALL IJtoIBJB(bfr_in,gslm)
             END IF
             IF (irec == 12) THEN
                irec=1
             ELSE
                irec=irec+1    
             END IF
  
             CALL ReadGetSLM(nfslm,irec,bfr_in,'txt')
          
             IF (reducedGrid) THEN
                CALL AveBoxIJtoIBJB(bfr_in,bfr_out)
             ELSE
                CALL IJtoIBJB(bfr_in,bfr_out)
             END IF
             CLOSE(UNIT=nfslm)
             gmax=-1.0e10_r8
             gmin=+1.0e10_r8
             DO j=1,jbMax
                DO i=1,ibMaxPerJB(j)
                   gslm(i,j)=f2*gslm(i,j)+f1*bfr_out(i,j)
                   gmax=MAX(gmax,gslm(i,j))
                   gmin=MIN(gmin,gslm(i,j))
                END DO
             END DO
             IF (ifslm == 3.AND.tod == 0.0_r8.AND.ifday == 0) THEN
                ifslm=0
             END IF
             IF (nfctrl(23) >= 1) THEN
                WRITE(UNIT=nfprt,FMT=222) mf,f1,f2,gmax,gmin
             END IF
          ELSE
             WRITE(UNIT=nfprt,FMT=333)
             STOP
          END IF
       END IF
    END IF

    !
    ! Process CO2 file/field/value
    !

    IF(ifco2.EQ.-1) THEN
       CALL getco2(idatec,co2val)
    ELSEIF(ifco2.EQ.1) THEN
       !CALL READ_MONTH_CO2
    ELSEIF(ifco2.EQ.2) THEN
    ELSEIF(ifco2.EQ.3) THEN
    ELSEIF(ifco2.EQ.4) THEN
    ENDIF

    !
    ! Process ozone file
    !

    IF (ifozone /= 0) THEN
       !   =1    read field from single month file (first call only)
       IF (ifozone == 1) THEN
          INQUIRE (IOLENGTH=LRecIn) rbrf
          OPEN (UNIT=nfozone, FILE=TRIM(fNameOzone), FORM='UNFORMATTED', &
          ACCESS='DIRECT', RECL=LRecIn*kMax, ACTION='READ', STATUS='OLD', IOSTAT=ierr)
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameOzone), ierr
             STOP "**(ERROR)**"
          END IF
          CALL ReadOzone(nfozone,bfr_in3,1)
          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(bfr_in3 ,gozo)
          ELSE
             CALL IJtoIBJB(bfr_in3 ,gozo)
          END IF
          CLOSE(UNIT=nfozone)
          ifozone=-1
          !   =2    interpolated to current day and time from 12 month clim
          !   =3    interpolated to current day and time from 12 month predicted field
       ELSE IF (ifozone == 2.OR. &
            (ifozone == 3.AND.tod == 0.0_r8.AND.ifday == 0)) THEN
          INQUIRE (IOLENGTH=lrecl) bfr_in3
          lrecl=lrecl/2
          OPEN(UNIT=nfozone,file=TRIM(fNameOzone),ACCESS='direct',&
               FORM='unformatted',RECL=lrecl,STATUS='old')
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameOzone), ierr
             STOP "**(ERROR)**"
          END IF

          CALL ReadOzone(nfozone,bfr_in3,mf)
          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(bfr_in3 ,gozo)
          ELSE
             CALL IJtoIBJB(bfr_in3 ,gozo)
          END IF

          mf=mf+1
          IF (mf == 13) mf=1
          CALL ReadOzone(nfozone,bfr_in3,mf)
          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(bfr_in3 ,bfr_out3)
          ELSE
             CALL IJtoIBJB(bfr_in3 ,bfr_out3)
          END IF

          CLOSE (UNIT=nfozone)
          gmax=-1.0e10_r8
          gmin=+1.0e10_r8
          DO j=1,jbMax
             DO i=1,ibMaxPerJB(j)
                DO k=1,kMax
                   gozo(i,k,j)=f2*gozo(i,k,j)+f1*bfr_out3(i,k,j)
                   gmax=MAX(gmax,gozo(i,k,j))
                   gmin=MIN(gmin,gozo(i,k,j))
                END DO
             END DO
          END DO
          IF (ifozone == 3) THEN
             ifozone=-3
          END IF
          IF (nfctrl(23) >= 1) THEN
             WRITE(UNIT=nfprt,FMT=223) mf,f1,f2,gmax,gmin
          END IF
          !   =4    interpolated from continuous direct access data set to current day and time
       ELSE IF (ifozone == 4) THEN
          WRITE(UNIT=nfprt,FMT=*) 'ERROR: DIRECT ACCESS OZONE FILE NOT IMPLEMENTED! ABORTING...'
          STOP
       END IF
    END IF


    !
    ! Process tracer file
    !

    IF (iftracer /= 0) THEN
       !   =1    read field from single month file (first call only)
       IF (iftracer == 1) THEN
          !INQUIRE (IOLENGTH=LRecIn) rbrf
          !OPEN (UNIT=nftrc, FILE=TRIM(fNametracer), FORM='UNFORMATTED', &
          !ACCESS='DIRECT', RECL=LRecIn, ACTION='READ', STATUS='OLD', IOSTAT=ierr)
          OPEN (UNIT=nftrc, FILE=TRIM(fNametracer), FORM='UNFORMATTED', &
          ACCESS='SEQUENTIAL', ACTION='READ', STATUS='OLD', IOSTAT=ierr)
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNametracer), ierr
             STOP "**(ERROR)**"
          END IF
          CALL ReadTracer(nftrc,bfr_in3)
          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(bfr_in3 ,tracermix)
          ELSE
             CALL IJtoIBJB(bfr_in3 ,tracermix)
          END IF
          CLOSE(UNIT=nftrc)
          iftracer=-1
          !   =2    interpolated to current day and time from 12 month clim
          !   =3    interpolated to current day and time from 12 month predicted field
       ELSE IF (iftracer == 2.OR. &
            (iftracer == 3.AND.tod == 0.0_r8.AND.ifday == 0)) THEN
          INQUIRE (IOLENGTH=lrecl) bfr_in3
          lrecl=lrecl/2
          OPEN(UNIT=nftrc,file=TRIM(fNametracer),ACCESS='direct',&
               FORM='unformatted',RECL=lrecl,STATUS='old')
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNametracer), ierr
             STOP "**(ERROR)**"
          END IF

          CALL ReadTracer(nftrc,bfr_in3,mf)
          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(bfr_in3 ,tracermix)
          ELSE
             CALL IJtoIBJB(bfr_in3 ,tracermix)
          END IF

          mf=mf+1
          IF (mf == 13) mf=1
          CALL ReadTracer(nftrc,bfr_in3,mf)
          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(bfr_in3 ,bfr_out3)
          ELSE
             CALL IJtoIBJB(bfr_in3 ,bfr_out3)
          END IF

          CLOSE (UNIT=nftrc)
          gmax=-1.0e10_r8
          gmin=+1.0e10_r8
          DO j=1,jbMax
             DO i=1,ibMaxPerJB(j)
                DO k=1,kMax
                   tracermix(i,k,j)=f2*tracermix(i,k,j)+f1*bfr_out3(i,k,j)
                   gmax=MAX(gmax,tracermix(i,k,j))
                   gmin=MIN(gmin,tracermix(i,k,j))
                END DO
             END DO
          END DO
          IF (iftracer == 3) THEN
             iftracer=-3
          END IF
          IF (nfctrl(23) >= 1) THEN
             WRITE(UNIT=nfprt,FMT=223) mf,f1,f2,gmax,gmin
          END IF
          !   =4    interpolated from continuous direct access data set to current day and time
       ELSE IF (iftracer == 4) THEN
          WRITE(UNIT=nfprt,FMT=*) 'ERROR: DIRECT ACCESS OZONE FILE NOT IMPLEMENTED! ABORTING...'
          STOP
       END IF
    END IF


222 FORMAT(' SOILM   START MONTH=',i2,'  F1,F2=',2f6.3,'  MAX,MIN=',2e12.5)
223 FORMAT(' OZONE   START MONTH=',i2,'  F1,F2=',2f6.3,'  MAX,MIN=',2e12.5)
333 FORMAT(' ABNORMAL END IN SUBR.GETSBC AT SOILM  INTERPOLATION')
444 FORMAT(' SNOW HAS ONLY ONE FILE','  MAX,MIN=',2E12.5)
555 FORMAT(' ABNORMAL END IN SUBR.GETSBC AT SNOW   INTERPOLATION')
666 FORMAT(' SST START REC (MONTH+2) =',I5, &
         '  F1,F2=',2G13.6,'  MAX,MIN=',2G12.5)
667 FORMAT(' SST:  IFSST=',I2,'  MAX,MIN=',2G12.5)
888 FORMAT(' ALBEDO  START MONTH=',I2, &
         '  F1,F2=',2F6.3,'  MAX,MIN=',2E12.5)
999 FORMAT(' ABNORMAL END IN SUBR.GETSBC AT ALBEDO INTERPOLATION')
1999 FORMAT('ABNORMAL END IN SUBR.GETSBC AT SST   INTERPOLATION')
  END SUBROUTINE getsbc



  !
  ! getsbc :read surface boundary conditions.
  !
  SUBROUTINE getsbc_old (ncols,jmax,galb,gsst,gslm,gsnw,ifday,tod,&
                     idate,idatec,nfprt,ifprt,jdt,ifalb,ifsst,& 
                     ifslm,ifsnw,sstlag,intsst,fint,tice,yrl ,& 
                     monl) 
    INTEGER, INTENT(in   ) :: ncols  
    INTEGER, INTENT(in   ) :: jmax  
    REAL   (KIND=r8), INTENT(out  ) :: galb(ncols*jmax)
    REAL   (KIND=r8), INTENT(out  ) :: gsst(ncols*jmax)
    REAL   (KIND=r8), INTENT(out  ) :: gslm(ncols*jmax)
    REAL   (KIND=r8), INTENT(out  ) :: gsnw(ncols*jmax)
    INTEGER, INTENT(in   ) :: ifday
    REAL   (KIND=r8), INTENT(in   ) :: tod
    INTEGER, INTENT(in   ) :: idate(4) 
    INTEGER, INTENT(in   ) :: idatec(4)
    INTEGER, INTENT(in   ) :: nfprt       ! intent(in)
    INTEGER, INTENT(in   ) :: ifprt(100)  ! not used
    INTEGER, INTENT(in   ) :: jdt         ! not used
    INTEGER, INTENT(inout) :: ifalb 
    INTEGER, INTENT(inout) :: ifsst 
    INTEGER, INTENT(inout) :: ifslm 
    INTEGER, INTENT(inout) :: ifsnw 
    REAL   (KIND=r8), INTENT(in   ) :: sstlag 
    INTEGER, INTENT(in   ) :: intsst
    REAL   (KIND=r8), INTENT(in   ) :: fint  
    REAL   (KIND=r8), INTENT(in   ) :: tice
    REAL   (KIND=r8), INTENT(in   ) :: yrl 
    INTEGER, INTENT(in   ) :: monl(12)
    !
    !  common comtim
    !
    REAL   (KIND=r8)                :: xsst  (ncols*jmax)
    REAL   (KIND=r8)                :: bfr1  (ncols*jmax)
    REAL   (KIND=r8)                :: bfr2  (ncols*jmax)
    REAL   (KIND=r4)                :: bfra  (ncols*jmax)
    REAL   (KIND=r4)                :: bfrb  (ncols*jmax)
    REAL   (KIND=r4)                :: bfrc  (ncols*jmax)
    REAL   (KIND=r4)                :: bfrd  (ncols*jmax)
    INTEGER(KIND=i4)                :: lrecl
    CHARACTER(len=5)                :: type2
    CHARACTER(len=14)               :: fmtadj
    LOGICAL(KIND=i8)                         :: ly
    REAL   (KIND=r8)                :: fhr
    INTEGER                :: mon
    INTEGER                :: mf
    INTEGER                :: mnl
    INTEGER                :: mn
    INTEGER                :: mnlf
    INTEGER                :: month
    INTEGER                :: mm
    INTEGER                :: ij
    INTEGER                :: irec
    INTEGER                :: mfx
    INTEGER                :: mnln
    REAL   (KIND=r8)                :: yday
    REAL   (KIND=r8)                :: f1
    REAL   (KIND=r8)                :: f2
    REAL   (KIND=r8)                :: add
    REAL   (KIND=r8)                :: gmax
    REAL   (KIND=r8)                :: gmin
    REAL   (KIND=r8)                :: fsst
    REAL   (KIND=r8)                :: fisst
    REAL   (KIND=r8)                :: xx1
    REAL   (KIND=r8)                :: xx2
    REAL   (KIND=r8)                :: xday
    REAL   (KIND=r8)                :: zday
    !
    !   ifxxx=0    xxx is not processed
    !   ifxxx=1    xxx is set to month=idate(2) in the first call,
    !              but not processed from the subsequent calls.
    !              ifxxx is set to zero after interpolation
    !   ifxxx=2    xxx is interpolated to current day and time every fint
    !              hours synchronized to 00z regardless of initial time.
    !              interpolation is continuous (every time step) if fint<0.
    !   ifxxx=3    xxx is interpolated to current day and time when ifday=0
    !              and tod=0.0 but not processed otherwise
    !              ( appropriate only when xxx is predicted )
    !
    !              the following are for sst only (fint applies as in
    !              ifxxx=2):
    !   ifsst=4    sst is linearly interpolated from continuous direct
    !              access data set to current day and time.  data set
    !              is assumed to be spaced every intsst days or every 
    !              calendar month is intsst < 0.
    !   ifsst=5    sst is expanded from piecewise cubic coefficients in
    !              direct access data set to current day and time.  data set
    !              is assumed to be spaced every intsst days.
    !   note:      for ifsst=4 or 5 sstlag must be set.  sstlag is the
    !              number of days plus any fraction prior to the initial 
    !              condition date and time the data set begins if intsst > 0
    !              sstlag is the number of whole months prior to the initial 
    !              condition date the data set begins if intsst < 0.
    !
    !     ifsst=-1 for numerical weather forecasting using mean weekly sst:
    !              sst is read in the first call as the second record of 
    !              the archieve but not processed from the subsequent calls. 
    !

    REWIND(nfsnw) 
    REWIND(nfalb) 
    REWIND(nfslm) 
      
    IF (ifsst >= 4) THEN
      IF(fsbc) THEN
       INQUIRE (IOLENGTH=lrecl) bfra 
       lrecl=lrecl/2
       OPEN(nfsst,file=TRIM(fNameSSTAOI),ACCESS='DIRECT',&
       FORM='UNFORMATTED',RECL=lrecl,STATUS='unknown')
      END IF
    ELSE
    REWIND(nfsst) 
    OPEN(nfsst, file=TRIM(fNameSSTAOI), &
         ACTION="read",FORM="UNFORMATTED")
    END IF
    
    OPEN(nfsnw, file=TRIM(fNameSnow), &
        ACTION="read",FORM="UNFORMATTED")  
        
    OPEN(nfalb, file=TRIM(fNameSibAlb), &
        ACTION="read",FORM="UNFORMATTED")   
       
    OPEN(nfslm, file=TRIM(fNameSoilms), &
        ACTION="read",FORM="UNFORMATTED")         


    IF (fint > 0.0_r8) THEN
       fhr=float(idate(1))+tod/3600.0_r8+1.0e-3_r8
       IF (.NOT. fsbc .AND. ABS( MOD(fhr,fint)) > 1.0e-2_r8) THEN
          RETURN
       END IF
    END IF

    mon=idatec(2)
    yday=FLOAT(idatec(3))+float(idatec(1))/24.0_r8+ MOD(tod,3600.0_r8)/86400.0_r8
    mf=mon-1
    ly=yrl == 365.25_r8 .AND. MOD(idatec(4),4) == 0
    mnl=monl(mon)

    IF (ly .AND. mon == 2) THEN
       mnl=29
    END IF
    IF (yday > 1.0_r8+float(mnl)/2.0_r8) THEN
       mf=mon
    END IF

    mn=mf+1

    IF (mf <  1) THEN
       mf=12
    END IF
    IF (mn > 12) THEN
       mn=1
    END IF

    mnlf=monl(mf)

    IF (ly .AND. mf == 2) THEN
       mnlf=29
    END IF

    add=float(mnlf)/2.0_r8-1.0_r8

    IF (mf == mon) THEN
       add=-add-2.0_r8
    END IF

    mnln=monl(mn)

    IF (ly .AND. mn == 2) THEN
       mnln=29
    END IF

    f1=2.0_r8*(yday+add)/float(mnlf+mnln)
    f2=1.0_r8-f1

    IF (ifalb /= 0) THEN
       IF (ifalb == 1) THEN
          REWIND nfalb
          month=idate(2)
          DO mm=1,month
             CALL ReadGetALB(nfalb,bfr1)
          END DO
          DO ij=1,ncols*jmax
             galb(ij)=bfr1(ij)
          END DO
          ifalb=0
       ELSE IF (&
            (ifalb == 2) .OR. &
            (ifalb == 3 .AND. tod == 0.0_r8 .AND. ifday == 0)) THEN
          REWIND nfalb
          DO mm=1,mf
              CALL ReadGetALB(nfalb,bfr1)
          END DO
          IF (mf == 12) THEN
             REWIND nfalb
          END IF
          CALL ReadGetALB(nfalb,bfr2)
          gmax=-1.0e10_r8
          gmin=+1.0e10_r8
          DO ij=1,ncols*jmax
             galb(ij)=f2*bfr1(ij)+f1*bfr2(ij)
             gmax=MAX(gmax,galb(ij))
             gmin=MIN(gmin,galb(ij))
          END DO
          IF (ifalb == 3 .AND. tod == 0.0_r8 .AND. ifday == 0) THEN
             ifalb=0
          END IF
          REWIND nfalb
          IF (ifprt(23) >= 1) THEN
             WRITE(nfprt,888) mf,f1,f2,gmax,gmin
          END IF
       ELSE
          WRITE(nfprt,999)
       END IF
    END IF


    IF (ifsst /= 0) THEN
       IF (ifsst == -1) THEN
          REWIND nfsst
          CALL ReadGetSST(nfsst,bfrb)
          CALL ReadGetSST(nfsst,bfra)
          gmax=-1.0e10_r8
          gmin=+1.0e10_r8
          DO ij=1,ncols*jmax
             IF (bfra(ij) > 10.0_r8) THEN
                gsst(ij)=-bfra(ij)
             ELSE IF (bfra(ij) < 0.0_r8) THEN
                gsst(ij)=290.0_r8
             END IF
             gmax=MAX(gmax,gsst(ij))
             gmin=MIN(gmin,gsst(ij))
          END DO
          WRITE(6,667) ifsst,gmax,gmin
          ifsst=0
       ELSE IF (ifsst == 1) THEN
          REWIND nfsst
          READ(nfsst)         
          month=idate(2)
          DO mm=1,month
             CALL ReadGetSST(nfsst,bfra)
          END DO
          DO ij=1,ncols*jmax
             IF (bfra(ij) > 10.0_r8) THEN
                gsst(ij)=-bfra(ij)
             ELSE IF (bfra(ij) < 0.0_r8) THEN
                gsst(ij)=290.0_r8
             END IF
          END DO
          ifsst=0
       ELSE IF (ifsst == 2.OR. &
            (ifsst == 3.AND.tod == 0.0.AND.ifday == 0)) THEN
          REWIND nfsst
          READ(nfsst)
          DO  mm=1,mf
             CALL ReadGetSST(nfsst,bfra)
          END DO
          IF (mf == 12) THEN
             REWIND nfsst
             READ(nfsst)
          END IF
          CALL ReadGetSST(nfsst,bfrb)
          gmax=-1.0e10_r8
          gmin=+1.0e10_r8
          DO ij=1,ncols*jmax
             fsst=f2*bfra(ij)+f1*bfrb(ij)
             IF (fsst > gmax) THEN
                gmax=fsst
             END IF
             IF (fsst < gmin) THEN
                gmin=fsst
             END IF
             IF (fsst > 10.0_r8) THEN
                xsst(ij)=-fsst
             ELSE IF (fsst < 0.0_r8) THEN
                xsst(ij)=290.0_r8
             END IF
          END DO
          DO ij=1,ncols*jmax
             IF (tod == 0.0_r8.AND.ifday == 0) THEN
                gsst(ij)=xsst(ij)
             ELSE IF (xsst(ij) < 0.0_r8.AND.ABS(xsst(ij)) >= tice) THEN
                gsst(ij)=xsst(ij)
             ELSE IF (xsst(ij) < 0.0_r8.AND.ABS(gsst(ij)) >= tice) THEN
                gsst(ij)=-tice+1.0e-2_r8
             END IF
          END DO
          IF (ifsst == 3.AND.tod == 0.0.AND.ifday == 0) THEN
             ifsst=0
          END IF
          REWIND nfsst
          IF (ifprt(23) >= 1) THEN
             WRITE(nfprt,666) mf,f1,f2,gmax,gmin
          END IF
       ELSE IF (ifsst == 4) THEN
          IF (intsst > 0) THEN
             fisst=float(intsst)
             xday=ifday+tod/86400.0_r8+sstlag
             irec=xday/fisst+1.0e-3_r8+1
             xx1= MOD(xday,fisst)/fisst
             xx2=1.0_r8-xx1
          ELSE
             mon=idatec(2)
             xday=idatec(3)+float(idatec(1))/24.0_r8+MOD(tod,3600.0_r8)/86400.0_r8
             zday=1.0_r8+float(monl(mon))/2.0_r8           
             mfx=mon-1
             IF (xday > zday) THEN
                mfx=mon           
             END IF
             irec=12*(idatec(4)-idate(4))+mfx-idate(2)+INT(sstlag)+1
             xx1=f1
             xx2=f2
          END IF
          irec=irec+2
          CALL ReadGetSST2(nfsst,bfra,irec)
          CALL ReadGetSST2(nfsst,bfrb,irec+1)
          gmax=-1.0e10_r8
          gmin=+1.0e10_r8
          DO ij=1,ncols*jmax
             fsst=xx2*bfra(ij)+xx1*bfrb(ij)
             IF (fsst > gmax) THEN
                gmax=fsst
             END IF
             IF (fsst < gmin) THEN
                gmin=fsst
             END IF
             IF (fsst > 10.0_r8) THEN
                xsst(ij)=-fsst
             ELSE IF (fsst < 0.0_r8) THEN
                xsst(ij)=290.0_r8
             END IF
          END DO
          IF (ifprt(23) >= 1) THEN
             WRITE(nfprt,666) irec,xx1,xx2,gmax,gmin
          END IF
          DO ij=1,ncols*jmax
             IF (tod == 0.0_r8.AND.ifday == 0) THEN
                gsst(ij)=xsst(ij)
             ELSE IF (xsst(ij) < 0.0_r8.AND.ABS(xsst(ij)) >= tice) THEN
                gsst(ij)=xsst(ij)
             ELSE IF (xsst(ij) < 0.0_r8.AND.ABS(gsst(ij)) >= tice) THEN
                gsst(ij)=-tice+1.0e-2_r8
             END IF
          END DO
       ELSE IF (ifsst == 5.AND.intsst > 0) THEN
          fisst=float(intsst)
          zday=ifday+tod/86400.0_r8
          xday=zday+sstlag
          irec=4*INT(xday/fisst+1.0e-3_r8)+1
          !*jpb  +2 due to land sea mask
          irec=irec+2
          CALL ReadGetSST2(nfsst,bfra,irec)
          CALL ReadGetSST2(nfsst,bfra,irec+1)
          CALL ReadGetSST2(nfsst,bfra,irec+2)
          CALL ReadGetSST2(nfsst,bfra,irec+3)
          gmax=-1.0e10_r8
          gmin=+1.0e10_r8
          !*jpb   i do not understand this case (ifsst=5):
          !*jpb   xday is not a piecewise cubic coefficient ...
          !*jpb   it seems that this approach does not work
          DO ij=1,ncols*jmax
             fsst=bfra(ij)+ &
                  xday*(bfrb(ij)+xday*(bfrc(ij)+xday*bfrd(ij)))
             IF (fsst > gmax) THEN
                gmax=fsst
             END IF
             IF (fsst < gmin) THEN
                gmin=fsst
             END IF
             IF (fsst > 10.0_r8) THEN
                xsst(ij)=-fsst
             ELSE IF (fsst < 0.0_r8) THEN
                xsst(ij)=290.0_r8
             END IF
          END DO
          xx1=zday
          xx2=xday
          IF (ifprt(23) >= 1) THEN
             WRITE(nfprt,666) irec,xx1,xx2,gmax,gmin
          END IF
          DO ij=1,ncols*jmax
             IF (tod == 0.0_r8.AND.ifday == 0) THEN
                gsst(ij)=xsst(ij)
             ELSE IF (xsst(ij) < 0.0_r8.AND.ABS(xsst(ij)) >= tice) THEN
                gsst(ij)=xsst(ij)
             ELSE IF (xsst(ij) < 0.0_r8.AND.ABS(gsst(ij)) >= tice) THEN
                gsst(ij)=-tice+1.0e-2_r8
             END IF
          END DO
       ELSE
          WRITE(nfprt,1999)
       END IF
    END IF

    IF (ifsnw /= 0) THEN
       IF (ifsnw == 1) THEN
          REWIND nfsnw
          CALL ReadGetSNW(nfsnw,bfra)
          DO ij=1,ncols*jmax
             gsnw(ij)=bfra(ij)
          END DO
          ifsnw=0
       ELSE IF (ifsnw == 2.OR. &
            (ifsnw == 3.AND.tod == 0.0_r8.AND.ifday == 0)) THEN
          REWIND nfsnw
          CALL ReadGetSNW(nfsnw,bfra)
          gmax=-1.0e10_r8
          gmin=+1.0e10_r8
          DO ij=1,ncols*jmax
             gsnw(ij)=bfra(ij)
             gmax=MAX(gmax,gsnw(ij))
             gmin=MIN(gmin,gsnw(ij))
          END DO
          IF (ifsnw == 3.AND.tod == 0.0_r8.AND.ifday == 0) THEN
             ifsnw=0
          END IF
          REWIND nfsnw
          IF (ifprt(23) >= 1) THEN
             WRITE(nfprt,444) gmax,gmin
          END IF
       ELSE
          WRITE(nfprt,555)
       END IF
    END IF


    IF (ifslm /= 0) THEN
       IF (ifslm == 1) THEN
          REWIND nfslm
          month=idate(2)
          DO mm=1,month
             CALL ReadGetSLM(nfslm,bfr1)
          END DO
          DO ij=1,ncols*jmax
             gslm(ij)=bfr1(ij)
          END DO
          ifslm=0
       ELSE IF (ifslm == 2.OR. &
            (ifslm == 3.AND.tod == 0.0_r8.AND.ifday == 0)) THEN
          REWIND nfslm
          DO mm=1,mf
             CALL ReadGetSLM(nfslm,bfr1)
          END DO
          IF (mf == 12) THEN
             REWIND nfslm
          END IF
           CALL ReadGetSLM(nfslm,bfr2)
          gmax=-1.0e10_r8
          gmin=+1.0e10_r8
          DO ij=1,ncols*jmax
             gslm(ij)=f2*bfr1(ij)+f1*bfr2(ij)
             gmax=MAX(gmax,gslm(ij))
             gmin=MIN(gmin,gslm(ij))
          END DO
          IF (ifslm == 3.AND.tod == 0.0.AND.ifday == 0) THEN
             ifslm=0
          END IF
          REWIND nfslm
          IF (ifprt(23) >= 1) THEN
             WRITE(nfprt,222) mf,f1,f2,gmax,gmin
          END IF
       ELSE
          WRITE(nfprt,333)
       END IF
    END IF

222 FORMAT(' SOILM   START MONTH=',i2,'  F1,F2=',2f6.3,'  MAX,MIN=',2e12.5)
333 FORMAT(' ABNORMAL END IN SUBR.GETSBC AT SOILM  INTERPOLATION')
444 FORMAT(' SNOW HAS ONLY ONE FILE','  MAX,MIN=',2E12.5)
555 FORMAT(' ABNORMAL END IN SUBR.GETSBC AT SNOW   INTERPOLATION')
666 FORMAT(' SST     START MONTH=',I2, &
         '  F1,F2=',2G13.6,'  MAX,MIN=',2G12.5)
667 FORMAT(' SST:  IFSST=',I2,'  MAX,MIN=',2G12.5)
888 FORMAT(' ALBEDO  START MONTH=',I2, &
         '  F1,F2=',2F6.3,'  MAX,MIN=',2E12.5)
999 FORMAT(' ABNORMAL END IN SUBR.GETSBC AT ALBEDO INTERPOLATION')
1999 FORMAT('ABNORMAL END IN SUBR.GETSBC AT SST   INTERPOLATION')
  END SUBROUTINE getsbc_old


  SUBROUTINE getsbc1D (ncols,jmax,galb,gsst,gslm,gsnw,ifday,tod,&
                       idate,idatec,nfprt,ifprt,jdt,ifalb,ifsst,& 
                       ifslm,ifsnw,sstlag,intsst,fint,tice,yrl ,& 
                       monl) 
    INTEGER, INTENT(in   ) :: ncols  
    INTEGER, INTENT(in   ) :: jmax  
    REAL   (KIND=r8), INTENT(out  ) :: galb(ncols*jmax)
    REAL   (KIND=r8), INTENT(out  ) :: gsst(ncols*jmax)
    REAL   (KIND=r8), INTENT(out  ) :: gslm(ncols*jmax)
    REAL   (KIND=r8), INTENT(out  ) :: gsnw(ncols*jmax)
    INTEGER, INTENT(in   ) :: ifday
    REAL   (KIND=r8), INTENT(in   ) :: tod
    INTEGER, INTENT(in   ) :: idate(4) 
    INTEGER, INTENT(in   ) :: idatec(4)
    INTEGER, INTENT(in   ) :: nfprt       ! intent(in)
    INTEGER, INTENT(in   ) :: ifprt(100)  ! not used
    INTEGER, INTENT(in   ) :: jdt         ! not used
    INTEGER, INTENT(inout) :: ifalb 
    INTEGER, INTENT(inout) :: ifsst 
    INTEGER, INTENT(inout) :: ifslm 
    INTEGER, INTENT(inout) :: ifsnw 
    REAL   (KIND=r8), INTENT(in   ) :: sstlag 
    INTEGER, INTENT(in   ) :: intsst
    REAL   (KIND=r8), INTENT(in   ) :: fint  
    REAL   (KIND=r8), INTENT(in   ) :: tice
    REAL   (KIND=r8), INTENT(in   ) :: yrl 
    INTEGER, INTENT(in   ) :: monl(12)
    !
    !  common comtim
    !
    REAL   (KIND=r8)                :: xsst  (ncols*jmax)
    REAL   (KIND=r8)                :: bfr1  
    REAL   (KIND=r8)                :: bfr2  
    REAL   (KIND=r4)                :: bfra
    REAL   (KIND=r4)                :: bfrb
    REAL   (KIND=r4)                :: bfrc 
    REAL   (KIND=r4)                :: bfrd 
    INTEGER                :: lrecl
    CHARACTER(len=5)                :: type2
    CHARACTER(len=14)               :: fmtadj
    LOGICAL(KIND=i8)                         :: ly
    REAL   (KIND=r8)                :: fhr
    INTEGER                :: mon
    INTEGER                :: mf
    INTEGER                :: mnl
    INTEGER                :: mn
    INTEGER                :: mnlf
    INTEGER                :: month
    INTEGER                :: mm
    INTEGER                :: ij
    INTEGER                :: irec
    INTEGER                :: mfx
    INTEGER                :: mnln
    REAL   (KIND=r8)                :: yday
    REAL   (KIND=r8)                :: f1
    REAL   (KIND=r8)                :: f2
    REAL   (KIND=r8)                :: add
    REAL   (KIND=r8)                :: gmax
    REAL   (KIND=r8)                :: gmin
    REAL   (KIND=r8)                :: fsst
    REAL   (KIND=r8)                :: fisst
    REAL   (KIND=r8)                :: xx1
    REAL   (KIND=r8)                :: xx2
    REAL   (KIND=r8)                :: xday
    REAL   (KIND=r8)                :: zday
    REAL   (KIND=r8)                :: alb(12)
    REAL   (KIND=r8)                :: snw(12)
    REAL   (KIND=r8)                :: slm(12)
    REAL   (KIND=r8)                :: sst(12)
    NAMELIST /MODEL_SBC/alb,snw,slm,sst

    !
    !   ifxxx=0    xxx is not processed
    !   ifxxx=1    xxx is set to month=idate(2) in the first call,
    !              but not processed from the subsequent calls.
    !              ifxxx is set to zero after interpolation
    !   ifxxx=2    xxx is interpolated to current day and time every fint
    !              hours synchronized to 00z regardless of initial time.
    !              interpolation is continuous (every time step) if fint<0.
    !   ifxxx=3    xxx is interpolated to current day and time when ifday=0
    !              and tod=0.0 but not processed otherwise
    !              ( appropriate only when xxx is predicted )
    !
    !              the following are for sst only (fint applies as in
    !              ifxxx=2):
    !   ifsst=4    sst is linearly interpolated from continuous direct
    !              access data set to current day and time.  data set
    !              is assumed to be spaced every intsst days or every 
    !              calendar month is intsst < 0.
    !   ifsst=5    sst is expanded from piecewise cubic coefficients in
    !              direct access data set to current day and time.  data set
    !              is assumed to be spaced every intsst days.
    !   note:      for ifsst=4 or 5 sstlag must be set.  sstlag is the
    !              number of days plus any fraction prior to the initial 
    !              condition date and time the data set begins if intsst > 0
    !              sstlag is the number of whole months prior to the initial 
    !              condition date the data set begins if intsst < 0.
    !
    !     ifsst=-1 for numerical weather forecasting using mean weekly sst:
    !              sst is read in the first call as the second record of 
    !              the archieve but not processed from the subsequent calls. 
    !
    REWIND 111
    READ (111,MODEL_SBC)
    IF (fint > 0.0_r8) THEN
       fhr=float(idate(1))+tod/3600.0_r8+1.0e-3_r8
       IF (.NOT. fsbc .AND. ABS( MOD(fhr,fint)) > 1.0e-2_r8) THEN
          RETURN
       END IF
    END IF

    mon=idatec(2)
    yday=idatec(3)+float(idatec(1))/24.0_r8+ MOD(tod,3600.0_r8)/86400.0_r8
    mf=mon-1
    ly=yrl == 365.25_r8 .AND. MOD(idatec(4),4) == 0
    mnl=monl(mon)

    IF (ly .AND. mon == 2) THEN
       mnl=29
    END IF
    IF (yday > 1.0_r8+float(mnl)/2.0_r8) THEN
       mf=mon
    END IF

    mn=mf+1

    IF (mf <  1) THEN
       mf=12
    END IF
    IF (mn > 12) THEN
       mn=1
    END IF

    mnlf=monl(mf)

    IF (ly .AND. mf == 2) THEN
       mnlf=29
    END IF

    add=float(mnlf)/2.0_r8-1.0_r8

    IF (mf == mon) THEN
       add=-add-2.0_r8
    END IF

    mnln=monl(mn)

    IF (ly .AND. mn == 2) THEN
       mnln=29
    END IF

    f1=2.0_r8*(yday+add)/float(mnlf+mnln)
    f2=1.0_r8-f1

    IF (ifalb /= 0) THEN
       IF (ifalb == 1) THEN
          month=idate(2)
          DO mm=1,month
             bfr1=alb(mm)
          END DO
          DO ij=1,ncols*jmax
             galb(ij)=bfr1
          END DO
          ifalb=0
       ELSE IF (&
            (ifalb == 2) .OR. &
            (ifalb == 3 .AND. tod == 0.0_r8 .AND. ifday == 0)) THEN
          DO mm=1,mf
              bfr1=alb(mm)
          END DO
          IF (mf == 12) THEN
             bfr2=alb(1)
          END IF
          gmax=-1.0e10_r8
          gmin=+1.0e10_r8
          DO ij=1,ncols*jmax
             galb(ij)=f2*bfr1+f1*bfr2
             gmax=MAX(gmax,galb(ij))
             gmin=MIN(gmin,galb(ij))
          END DO
          IF (ifalb == 3 .AND. tod == 0.0_r8 .AND. ifday == 0) THEN
             ifalb=0
          END IF
          REWIND nfalb
          IF (ifprt(23) >= 1) THEN
             WRITE(nfprt,888) mf,f1,f2,gmax,gmin
          END IF
       ELSE
          WRITE(nfprt,999)
       END IF
    END IF


    IF (ifsst /= 0) THEN
       IF (ifsst == -1) THEN
          bfrb=REAL(sst(2),kind=r4) 
          bfra=REAL(sst(1),kind=r4) 
          gmax=-1.0e10_r8
          gmin=+1.0e10_r8
          DO ij=1,ncols*jmax
             IF (bfra > 10.0_r8) THEN
                gsst(ij)=-bfra
             ELSE IF (bfra < 0.0_r8) THEN
                gsst(ij)=290.0_r8
             END IF
             gmax=MAX(gmax,gsst(ij))
             gmin=MIN(gmin,gsst(ij))
          END DO
          WRITE(6,667) ifsst,gmax,gmin
          ifsst=0
       ELSE IF (ifsst == 1) THEN
          month=idate(2)
          DO mm=1,month 
             bfra=REAL(sst(mm),kind=r4)
          END DO
          DO ij=1,ncols*jmax
             IF (bfra > 10.0_r8) THEN
                gsst(ij)=-bfra
             ELSE IF (bfra < 0.0_r8) THEN
                gsst(ij)=290.0_r8
             END IF
          END DO
          ifsst=0
       ELSE IF (ifsst == 2.OR. &
            (ifsst == 3.AND.tod == 0.0_r8.AND.ifday == 0)) THEN
          DO  mm=1,mf
             bfra=REAL(sst(mm) ,kind=r4)
          END DO
          IF (mf == 12) THEN
             bfrb=REAL(sst(1) ,kind=r4)
          END IF
          gmax=-1.0e10_r8
          gmin=+1.0e10_r8
          DO ij=1,ncols*jmax
             fsst=f2*bfra+f1*bfrb
             IF (fsst > gmax) THEN
                gmax=fsst
             END IF
             IF (fsst < gmin) THEN
                gmin=fsst
             END IF
             IF (fsst > 10.0_r8) THEN
                xsst(ij)=-fsst
             ELSE IF (fsst < 0.0_r8) THEN
                xsst(ij)=290.0_r8
             END IF
          END DO
          DO ij=1,ncols*jmax
             IF (tod == 0.0_r8.AND.ifday == 0) THEN
                gsst(ij)=xsst(ij)
             ELSE IF (xsst(ij) < 0.0_r8.AND.ABS(xsst(ij)) >= tice) THEN
                gsst(ij)=xsst(ij)
             ELSE IF (xsst(ij) < 0.0_r8.AND.ABS(gsst(ij)) >= tice) THEN
                gsst(ij)=-tice+1.0e-2_r8
             END IF
          END DO
          IF (ifsst == 3.AND.tod == 0.0_r8.AND.ifday == 0) THEN
             ifsst=0
          END IF
          IF (ifprt(23) >= 1) THEN
             WRITE(nfprt,666) mf,f1,f2,gmax,gmin
          END IF
       ELSE IF (ifsst == 4) THEN
          IF (intsst > 0) THEN
             fisst=float(intsst)
             xday=ifday+tod/86400.0_r8+sstlag
             irec=xday/fisst+1.0e-3_r8+1
             xx1= MOD(xday,fisst)/fisst
             xx2=1.0_r8-xx1
          ELSE
             mon=idatec(2)
             xday=idatec(3)+float(idatec(1))/24.0_r8+MOD(tod,3600.0_r8)/86400.0_r8
             zday=1.0_r8+float(monl(mon))/2.0_r8           
             mfx=mon-1
             IF (xday > zday) THEN
                mfx=mon           
             END IF
             irec=12*(idatec(4)-idate(4))+mfx-idate(2)+INT(sstlag)+1
             xx1=f1
             xx2=f2
          END IF
          irec=irec+1
          bfra=REAL(sst(irec  ),kind=r4)
          bfrb=REAL(sst(irec+1),kind=r4)
          gmax=-1.0e10_r8
          gmin=+1.0e10_r8
          DO ij=1,ncols*jmax
             fsst=xx2*bfra+xx1*bfrb
             IF (fsst > gmax) THEN
                gmax=fsst
             END IF
             IF (fsst < gmin) THEN
                gmin=fsst
             END IF
             IF (fsst > 10.0_r8) THEN
                xsst(ij)=-fsst
             ELSE IF (fsst < 0.0_r8) THEN
                xsst(ij)=290.0_r8
             END IF
          END DO
          IF (ifprt(23) >= 1) THEN
             WRITE(nfprt,666) irec,xx1,xx2,gmax,gmin
          END IF
          DO ij=1,ncols*jmax
             IF (tod == 0.0_r8.AND.ifday == 0) THEN
                gsst(ij)=xsst(ij)
             ELSE IF (xsst(ij) < 0.0_r8.AND.ABS(xsst(ij)) >= tice) THEN
                gsst(ij)=xsst(ij)
             ELSE IF (xsst(ij) < 0.0_r8.AND.ABS(gsst(ij)) >= tice) THEN
                gsst(ij)=-tice+1.0e-2_r8
             END IF
          END DO
       ELSE IF (ifsst == 5.AND.intsst > 0) THEN
          fisst=float(intsst)
          zday=ifday+tod/86400.0_r8
          xday=zday+sstlag
          irec=4*INT(xday/fisst+1.0e-3_r8)+1
          !*jpb  +2 due to land sea mask
          irec=irec+1
          bfra=REAL(sst(irec  ),KIND=r4)
          bfrb=REAL(sst(irec+1),KIND=r4)
          bfrc=REAL(sst(irec+2),KIND=r4)
          bfrd=REAL(sst(irec+3),KIND=r4)
          gmax=-1.0e10_r8
          gmin=+1.0e10_r8
          !*jpb   i do not understand this case (ifsst=5):
          !*jpb   xday is not a piecewise cubic coefficient ...
          !*jpb   it seems that this approach does not work
          DO ij=1,ncols*jmax
             fsst=bfra+ &
                  xday*(bfrb+xday*(bfrc+xday*bfrd))
             IF (fsst > gmax) THEN
                gmax=fsst
             END IF
             IF (fsst < gmin) THEN
                gmin=fsst
             END IF
             IF (fsst > 10.0_r8) THEN
                xsst(ij)=-fsst
             ELSE IF (fsst < 0.0_r8) THEN
                xsst(ij)=290.0_r8
             END IF
          END DO
          xx1=zday
          xx2=xday
          IF (ifprt(23) >= 1) THEN
             WRITE(nfprt,666) irec,xx1,xx2,gmax,gmin
          END IF
          DO ij=1,ncols*jmax
             IF (tod == 0.0_r8.AND.ifday == 0) THEN
                gsst(ij)=xsst(ij)
             ELSE IF (xsst(ij) < 0.0_r8.AND.ABS(xsst(ij)) >= tice) THEN
                gsst(ij)=xsst(ij)
             ELSE IF (xsst(ij) < 0.0_r8.AND.ABS(gsst(ij)) >= tice) THEN
                gsst(ij)=-tice+1.0e-2_r8
             END IF
          END DO
       ELSE
          WRITE(nfprt,1999)
       END IF
    END IF

    IF (ifsnw /= 0) THEN
       IF (ifsnw == 1) THEN 
          month=idate(2)          
          DO mm=1,month
            bfra=REAL(snw(mm),KIND=r4)
          END DO
          DO ij=1,ncols*jmax
             gsnw(ij)=bfra
          END DO
          ifsnw=0
       ELSE IF (ifsnw == 2.OR. &
            (ifsnw == 3.AND.tod == 0.0_r8.AND.ifday == 0)) THEN
          DO mm=1,mf
             bfra=REAL(snw(mm),KIND=r4)
          END DO
          IF (mf == 12) THEN
            bfr2=REAL(snw(1),KIND=r4)
          END IF
          gmax=-1.0e10_r8
          gmin=+1.0e10_r8
          DO ij=1,ncols*jmax
             gsnw(ij)=f2*bfra+f1*bfr2
             gmax=MAX(gmax,gsnw(ij))
             gmin=MIN(gmin,gsnw(ij))
          END DO
          IF (ifsnw == 3.AND.tod == 0.0_r8.AND.ifday == 0) THEN
             ifsnw=0
          END IF
          REWIND nfsnw
          IF (ifprt(23) >= 1) THEN
             WRITE(nfprt,444) gmax,gmin
          END IF
       ELSE
          WRITE(nfprt,555)
       END IF
    END IF


    IF (ifslm /= 0) THEN
       IF (ifslm == 1) THEN
          month=idate(2)
          DO mm=1,month
            bfr1=slm(mm)
          END DO
          DO ij=1,ncols*jmax
             gslm(ij)=bfr1
          END DO
          ifslm=0
       ELSE IF (ifslm == 2.OR. &
            (ifslm == 3.AND.tod == 0.0_r8.AND.ifday == 0)) THEN
          DO mm=1,mf
             bfr1=slm(mm)
          END DO
          IF (mf == 12) THEN
            bfr2=slm(1)
          END IF
          gmax=-1.0e10_r8
          gmin=+1.0e10_r8
          DO ij=1,ncols*jmax
             gslm(ij)=f2*bfr1+f1*bfr2
             gmax=MAX(gmax,gslm(ij))
             gmin=MIN(gmin,gslm(ij))
          END DO
          IF (ifslm == 3.AND.tod == 0.0_r8.AND.ifday == 0) THEN
             ifslm=0
          END IF
          REWIND nfslm
          IF (ifprt(23) >= 1) THEN
             WRITE(nfprt,222) mf,f1,f2,gmax,gmin
          END IF
       ELSE
          WRITE(nfprt,333)
       END IF
    END IF
222 FORMAT(' SOILM   START MONTH=',i2,'  F1,F2=',2f6.3,'  MAX,MIN=',2e12.5)
333 FORMAT(' ABNORMAL END IN SUBR.GETSBC AT SOILM  INTERPOLATION')
444 FORMAT(' SNOW HAS ONLY ONE FILE','  MAX,MIN=',2E12.5)
555 FORMAT(' ABNORMAL END IN SUBR.GETSBC AT SNOW   INTERPOLATION')
666 FORMAT(' SST     START MONTH=',I2, &
         '  F1,F2=',2G13.6,'  MAX,MIN=',2G12.5)
667 FORMAT(' SST:  IFSST=',I2,'  MAX,MIN=',2G12.5)
888 FORMAT(' ALBEDO  START MONTH=',I2, &
         '  F1,F2=',2F6.3,'  MAX,MIN=',2E12.5)
999 FORMAT(' ABNORMAL END IN SUBR.GETSBC AT ALBEDO INTERPOLATION')
1999 FORMAT('ABNORMAL END IN SUBR.GETSBC AT SST   INTERPOLATION')
  END SUBROUTINE getsbc1D


  SUBROUTINE GetRecWgtMonthlySST &
       (idate, idatec, tod, labelsi, labelsj, &
       irec, f1, f2, mra, mrb,monl)

    IMPLICIT NONE

    ! Computes the Corresponding Records to do Linear
    ! Time Interpolation and the Respectives Weights.

    INTEGER, INTENT (IN) :: idate(4), idatec(4),monl(12)
    REAL (KIND=r8), INTENT (IN) :: tod
    CHARACTER (LEN=10), INTENT (IN) :: labelsi, labelsj

    INTEGER, INTENT (OUT) :: irec, mra, mrb
    REAL (KIND=r8), INTENT (OUT) :: f1, f2

    ! Local Constants
    INTEGER :: ysi, msi, dsi, ysj, msj, dsj, ndij, nd, &
         tmca, tmcb, tmcf
    REAL (KIND=r8) :: xday, zdayf, zdaya, zdayb, tc

    ! Get Year, Month and Day of the Initial and Second Medium Date
    ! for SST Direct Access File Data

    READ (labelsi(1:4), '(I4)') ysi
    READ (labelsi(5:6), '(I2)') msi
    READ (labelsi(7:8), '(I2)') dsi
    READ (labelsj(1:4), '(I4)') ysj
    READ (labelsj(5:6), '(I2)') msj
    READ (labelsj(7:8), '(I2)') dsj

    ! Lag of Days for SST Data:
    ! Just for Checking if the Scale is a Month
    ndij=0
    IF (msi+1 <= msj-1) THEN
       DO nd=msi+1,msj-1
          ndij=ndij+monl(nd)
       END DO
    ELSE
       DO nd=msi+1,12
          ndij=ndij+monl(nd)
       END DO
       DO nd=1,msj-1
          ndij=ndij+monl(nd)
       END DO
    END IF
    ndij=ndij+monl(msi)-dsi+dsj+365*(ysj-ysi-1)

    ! Check for Monthly Scale SST Data
    IF (ABS(ndij) <= 27 .OR. ABS(ndij) >= 32) THEN
       WRITE (UNIT=0, FMT='(/,A)') ' *** Error: The SST Data Is Not On Monthly Scale   ***'
       WRITE (UNIT=0, FMT='(/,A,I8,12X,A,/)') ' *** Lag Of Days For SST Data: ', ndij, '***'
       WRITE (UNIT=0, FMT='(A,/)') ' *** Program STOP: SUBROUTINE GetRecWgtMonthlySST  ***'
       STOP
    END IF

    ! Length in Days of the Date of Forecasting
    tmcf=monl(idatec(2))
    IF (idatec(2) == 2 .AND. MOD(idatec(4),4) == 0) tmcf=29
    ! Medium Day of the Month of Forecasting
    zdayf=0.5_r8*REAL(tmcf,r8)+1.0_r8
    ! Fractional Day of Forecasting
    tc=REAL(idate(1),r8)/24.0_r8+tod/86400.0_r8
    ! Correcting Factor if Necessary (tc is in Days)
    IF (tc >= 1.0_r8) tc=tc-1.0_r8
    xday=REAL(idatec(3),r8)+tc
    ! Getting the Corresponding Record in SST Data
    irec=12-msi+idatec(2)+12*(idatec(4)-ysi-1)+2
    IF (xday >= zdayf) irec=irec+1

    ! Months for the Linear Time Interpolation Related to the Records
    mra=MOD(irec-3+msi,12)
    IF (mra == 0) mra=12
    mrb=mra+1
    IF (mrb > 12) mrb=1

    ! Length in Days for the First Month of Interpolation
    tmca=monl(mra)
    IF (mra == 2 .AND. MOD(ysi,4) == 0) tmca=29
    ! Medium Fracitonal Day for the First Month of Interpolation
    zdaya=0.5_r8*REAL(tmca,r8)+1.0_r8-REAL(tmca,r8)
    ! Length in Days for the Second Month of Interpolation
    tmcb=monl(mrb)
    IF (mrb == 2 .AND. MOD(ysj,4) == 0) tmcb=29
    ! Medium Fracitonal Day for the Second Month of Interpolation
    zdayb=0.5_r8*REAL(tmcb,r8)+1.0_r8
    ! Scaling Fractional Day of Forecasting, if Necessary
    IF (xday >= zdayf) xday=xday-REAL(tmca,r8)
    ! Interpolation Factors
    f1=(xday-zdaya)/(zdayb-zdaya)
    f2=1.0_r8-f1

  END SUBROUTINE GetRecWgtMonthlySST

  SUBROUTINE GetWeightsOld (yrl,monl,idatec, tod, f1, f2,mf)

    IMPLICIT NONE

    ! Computes Weights as in getsbc:

    INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15)
    INTEGER, INTENT (IN) :: idatec(4)
    INTEGER, INTENT (IN) :: monl(12)
    REAL (KIND=r8), INTENT (IN) :: tod
    REAL (KIND=r8), INTENT (IN) :: yrl
    REAL (KIND=r8), INTENT (OUT) :: f1, f2
    INTEGER,  INTENT (OUT):: mf
    INTEGER :: mon, mnl, mn, mnlf, mnln
    REAL (KIND=r8) :: yday, add
    LOGICAL :: ly

    mon=idatec(2)
    yday=REAL(idatec(3),r8)+REAL(idatec(1),r8)/24.0_r8+MOD(tod,3600.0_r8)/86400.0_r8
    mf=mon-1
    ly= yrl == 365.25_r8 .AND. MOD(idatec(4),4) == 0
    mnl=monl(mon)
    IF (ly .AND. mon == 2) mnl=29
    ! Em getsbc seria apenas >
    ! As consideracoes de interpolacao leva a >=
    IF (yday >= 1.0_r8+0.5_r8*REAL(mnl,r8)) mf=mon
    mn=mf+1
    IF (mf < 1) mf=12
    IF (mn > 12) mn=1
    mnlf=monl(mf)
    IF (ly .AND. mf == 2) mnlf=29
    add=0.5_r8*REAL(mnlf,r8)-1.0_r8
    IF (mf == mon) add=-add-2.0_r8
    mnln=monl(mn)
    IF (ly .AND. mn == 2) mnln=29
    f1=2.0_r8*(yday+add)/REAL(mnlf+mnln,r8)
    f2=1.0_r8-f1

  END SUBROUTINE GetWeightsOld


  !hmjb
  SUBROUTINE getco2(time,co2val)
    !==========================================================================
    ! getco2: Interpolates Mauna Loa data for a given time
    !
    ! *** Atmospheric CO2 concentrations (ppmv) derived from in situ  ***
    ! *** air samples collected at Mauna Loa Observatory, Hawaii      ***
    !
    ! Data:
    !
    !   http://cdiac.ornl.gov/trends/co2/contents.htm
    !   http://cdiac.ornl.gov/ftp/trends/co2/maunaloa.co2
    !
    ! Parabolic fitting by hbarbosa@cptec.inpe.br, 17 Jan 2007:
    !
    !   co2val = a*(time-2000)^2 + b*(time-2000) + c
    !
    !       a  = 0.0116696   +/- 0.0005706    (4.89%)
    !       b  = 1.79984     +/- 0.022        (1.222%)
    !       c  = 369         +/- 0.1794       (0.04863%)
    !
    !==========================================================================
    !     time.......date of current data
    !     time(1)....hour(00/12)
    !     time(2)....month
    !     time(3)....day of month
    !     time(4)....year
    !
    !    co2val....co2val is wgne standard value in ppm "co2val = /345.0/
    !==========================================================================

    IMPLICIT NONE
    REAL(KIND=r8), PARAMETER :: A = 0.0116696
    REAL(KIND=r8), PARAMETER :: B = 1.79984
    REAL(KIND=r8), PARAMETER :: C = 369.0

    INTEGER,       INTENT(IN ) :: time(4)
    REAL(KIND=r8), INTENT(OUT) :: co2val

    REAL(KIND=r8) :: TDIF

    tdif=time(4) + (time(2)-1.)/12. + (time(3)-1.+ time(1)/24.)/365. - 2000.

    co2val = A*tdif**2 + B*tdif + C

    !    WRITE(*,123) time,tdif+2000.,co2val
    !123 format('hmjb co2val date=',3(I2,1x),I4,' fyear=',F10.5,' val=',F7.3)

    RETURN
  END SUBROUTINE getco2
  !hmjb

END MODULE InputOutput
