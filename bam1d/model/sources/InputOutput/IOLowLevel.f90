MODULE IOLowLevel      
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ReadHead  
  PUBLIC :: GReadHead
  PUBLIC :: ReadField
  PUBLIC :: GReadField
  PUBLIC :: WriteHead  
  PUBLIC :: GWriteHead
  PUBLIC :: WriteField
  PUBLIC :: GWriteField
  PUBLIC :: WriteDir
  PUBLIC :: WriteDire
  PUBLIC :: ReadProgHead 
  PUBLIC :: GReadProgHead
  PUBLIC :: WriteProgHead
  PUBLIC :: ReadLandSeaMask
  PUBLIC :: ReadLandSeaMask2
  PUBLIC :: GetUnit
  PUBLIC :: ReadVar
  PUBLIC :: ReadGetALB
  PUBLIC :: ReadGetSST
  PUBLIC :: ReadGetSST2
  PUBLIC :: ReadGetSLM  
  PUBLIC :: ReadGetSNW
  PUBLIC :: ReadGetNFTGZ
  PUBLIC :: InitReadWriteSpec
  PUBLIC :: WriteDHNHead
  PUBLIC :: LandSeaMask
  PUBLIC :: ReadOzone !hmjb
  PUBLIC :: ReadTracer !hmjb
  PUBLIC :: ReadGetSLM3D

  INTERFACE WriteDHNHead
     MODULE PROCEDURE WriteDHNHead4, WriteDHNHead8
  END INTERFACE     
  INTERFACE ReadOzone
     MODULE PROCEDURE ReadOzone8
  END INTERFACE
  INTERFACE ReadTracer
     MODULE PROCEDURE ReadTracer8,ReadTracer8s
  END INTERFACE
  INTERFACE ReadGetNFTGZ
     MODULE PROCEDURE ReadNFTGZ4, ReadNFTGZ8,ReadNFTGZ4_rec, ReadNFTGZ8_rec, &
                      ReadNFTGZ4_txt, ReadNFTGZ8_txt
  END INTERFACE 
  INTERFACE ReadGetSNW 
     MODULE PROCEDURE ReadSNW4, ReadSNW8, ReadSNW4_rec, ReadSNW8_rec, &
                      ReadSNW4_txt, ReadSNW8_txt
  END INTERFACE   
  INTERFACE ReadGetSLM 
     MODULE PROCEDURE ReadSLM4, ReadSLM8, ReadSLM4_rec, ReadSLM8_rec, &
                      ReadSLM4_txt, ReadSLM8_txt
  END INTERFACE 
  INTERFACE ReadGetSLM3D
     MODULE PROCEDURE ReadSLM43D, ReadSLM83D
  END INTERFACE
  INTERFACE ReadGetALB
     MODULE PROCEDURE ReadAlb4, ReadAlb8, ReadAlb4_rec, ReadAlb8_rec
  END INTERFACE
  INTERFACE ReadGetSST 
     MODULE PROCEDURE ReadSST4, ReadSST8,  ReadSST4_rec, ReadSST8_rec, &
                      ReadSST4_txt, ReadSST8_txt
  END INTERFACE
  INTERFACE ReadGetSST2 
     MODULE PROCEDURE ReadSST4Rec, ReadSST8Rec, ReadSST4Rec2, ReadSST8Rec2
  END INTERFACE  
  INTERFACE ReadVar 
     MODULE PROCEDURE ReadVar4, ReadVar8, ReadVar4_txt, ReadVar8_txt
  END INTERFACE  
  INTERFACE GReadHead
     MODULE PROCEDURE GReadHead4, GReadHead8
  END INTERFACE
  INTERFACE ReadHead
     MODULE PROCEDURE ReadHead4, ReadHead8
  END INTERFACE
  INTERFACE WriteHead
     MODULE PROCEDURE WriteHead4, WriteHead8
  END INTERFACE
  INTERFACE GWriteHead
     MODULE PROCEDURE GWriteHead4, GWriteHead8
  END INTERFACE
  INTERFACE ReadProgHead
     MODULE PROCEDURE ReadProgHead4, ReadProgHead8
  END INTERFACE  
  INTERFACE GReadProgHead
     MODULE PROCEDURE GReadProgHead4, GReadProgHead8
  END INTERFACE
  INTERFACE WriteProgHead
     MODULE PROCEDURE WriteProgHead4, WriteProgHead8
  END INTERFACE
  INTERFACE WriteDir
     MODULE PROCEDURE WriteDir4, WriteDir8
  END INTERFACE  
  INTERFACE WriteDire
     MODULE PROCEDURE WriteDire4, WriteDire8 
  END INTERFACE
  INTERFACE GReadField
     MODULE PROCEDURE GReadField41D, GReadField42D, GReadField81D, GReadField82D
  END INTERFACE
  INTERFACE ReadField
     MODULE PROCEDURE ReadField41D, ReadField42D, ReadField81D, ReadField82D
  END INTERFACE
  INTERFACE WriteField
     MODULE PROCEDURE WriteField41D, WriteField42D, WriteField81D, WriteField82D
  END INTERFACE
  INTERFACE GWriteField
     MODULE PROCEDURE GWriteField41D, GWriteField42D, GWriteField81D, GWriteField82D
  END INTERFACE
  INTERFACE ReadLandSeaMask
     MODULE PROCEDURE ReadLandSeaMask4, ReadLandSeaMask8
  END INTERFACE
  INTERFACE ReadLandSeaMask2
     MODULE PROCEDURE ReadLandSeaMask2_4, ReadLandSeaMask2_8
  END INTERFACE  
  INTERFACE LandSeaMask
     MODULE PROCEDURE LandSeaMask4, LandSeaMask8
  END INTERFACE
  INTERFACE GetUnit
     MODULE PROCEDURE GetUnit4, GetUnit8
  END INTERFACE
  INTEGER, PARAMETER :: r4=SELECTED_REAL_KIND(6)
  INTEGER, PARAMETER :: i4=SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: r8=SELECTED_REAL_KIND(15)
  INTEGER, PARAMETER :: i8=SELECTED_INT_KIND(14)
  !INTEGER, PARAMETER :: r8=SELECTED_REAL_KIND(14)
  !INTEGER, PARAMETER :: i8=SELECTED_INT_KIND(8)
  REAL, PARAMETER :: MaxIEEEExpo=36  ! maximum IEEE Expo (close)
  REAL, PARAMETER :: MinIEEEExpo=-36 ! minimum IEEE Expo (close)
  CHARACTER(len=40), ALLOCATABLE :: reqdg(:)
  CHARACTER(len=40), ALLOCATABLE :: combf(:)
  LOGICAL   , ALLOCATABLE :: dodia(:)
  INTEGER   , ALLOCATABLE :: itavl(:)
  INTEGER   , ALLOCATABLE :: iavrq(:)
  INTEGER   , ALLOCATABLE :: nucf (:)
  INTEGER   , ALLOCATABLE :: lvrq (:)
  INTEGER   , ALLOCATABLE :: nurq (:)
  INTEGER   , ALLOCATABLE :: lvcf (:)
  INTEGER   , ALLOCATABLE :: itcf (:)
  INTEGER                         :: ndrq  
  INTEGER                         :: ncdg  
  INTEGER                         :: ndavl 
  INTEGER                         :: mxavl   
  INTEGER                         :: icf  
  INTEGER                         :: mMax 
  INTEGER                         :: mnMax
  INTEGER                         :: kMax 
  INTEGER                         :: ijMax
  INTEGER                         :: iMax 
  INTEGER                         :: jMax 
  INTEGER                         :: ibMax
  INTEGER                         :: jbMax
  INTEGER                         :: nfprt
  INTEGER                         :: nferr
  INTEGER                         :: ifprt(100)
  
CONTAINS
  SUBROUTINE InitReadWriteSpec(&
             ndrq_in ,ncdg_in ,ndavl_in,mxavl_in,icf_in  ,mMax_in , &
             mnMax_in,kMax_in ,ijMax_in,iMax_in ,jMax_in ,ibMax_in, &
             jbMax_in,nfprt_in,nferr_in,ifprt_in, &
             reqdg_in,combf_in,dodia_in,itavl_in,iavrq_in,&
             nucf_in ,lvrq_in ,nurq_in ,lvcf_in ,itcf_in )

    INTEGER , INTENT(IN   ) :: ndrq_in
    INTEGER , INTENT(IN   ) :: ncdg_in
    INTEGER , INTENT(IN   ) :: ndavl_in
    CHARACTER(len=40), INTENT(IN   ) :: reqdg_in(:)
    CHARACTER(len=40), INTENT(IN   ) :: combf_in(:)
    LOGICAL(KIND=i8), INTENT(IN   ) :: dodia_in(:)
    INTEGER , INTENT(IN   ) :: itavl_in(:)
    INTEGER , INTENT(IN   ) :: iavrq_in(:)
    INTEGER , INTENT(IN   ) :: nucf_in (:)
    INTEGER , INTENT(IN   ) :: lvrq_in (:)
    INTEGER , INTENT(IN   ) :: nurq_in (:)
    INTEGER , INTENT(IN   ) :: lvcf_in (:)
    INTEGER , INTENT(IN   ) :: itcf_in (:)
    INTEGER , INTENT(IN   ) :: mxavl_in   
    INTEGER , INTENT(IN   ) :: icf_in  
    INTEGER , INTENT(IN   ) :: mMax_in 
    INTEGER , INTENT(IN   ) :: mnMax_in
    INTEGER , INTENT(IN   ) :: kMax_in 
    INTEGER , INTENT(IN   ) :: ijMax_in
    INTEGER , INTENT(IN   ) :: iMax_in 
    INTEGER , INTENT(IN   ) :: jMax_in 
    INTEGER , INTENT(IN   ) :: ibMax_in
    INTEGER , INTENT(IN   ) :: jbMax_in
    INTEGER , INTENT(IN   ) :: nfprt_in
    INTEGER , INTENT(IN   ) :: nferr_in
    INTEGER , INTENT(IN   ) :: ifprt_in(:)

    ALLOCATE(reqdg(ndrq_in))
    ALLOCATE(combf(ncdg_in))  
    ALLOCATE(dodia(ndavl_in))
    ALLOCATE(itavl(ndavl_in))
    ALLOCATE(iavrq(ndavl_in))
    ALLOCATE(nucf (ncdg_in ))
    ALLOCATE(lvrq (ndrq_in ))
    ALLOCATE(nurq (ndrq_in ))
    ALLOCATE(lvcf (ncdg_in ))
    ALLOCATE(itcf (ncdg_in )) 
    mxavl =   mxavl_in
    icf   =     icf_in  
    mMax  =    mMax_in 
    mnMax =   mnMax_in
    kMax  =    kMax_in 
    ijMax =   ijMax_in
    iMax  =    iMax_in 
    jMax  =    jMax_in 
    ibMax =   ibMax_in
    jbMax =   jbMax_in
    nfprt =   nfprt_in
    nferr =   nferr_in
    ifprt =   ifprt_in
    ndrq  =  ndrq_in
    ncdg  =  ncdg_in
    ndavl =  ndavl_in
    
    reqdg =  reqdg_in
    combf =  combf_in
    dodia =  dodia_in
    itavl =  itavl_in
    iavrq =  iavrq_in
    nucf  =  nucf_in 
    lvrq  =  lvrq_in 
    nurq  =  nurq_in 
    lvcf  =  lvcf_in 
    itcf  =  itcf_in 
    
  END SUBROUTINE InitReadWriteSpec

  SUBROUTINE ReadHead4(n, ifday, tod, idate, idatec, si, sl)
    INTEGER(KIND=i8), INTENT(IN)  :: n
    INTEGER(KIND=i4), INTENT(OUT) :: ifday
    REAL   (KIND=r4), INTENT(OUT) :: tod
    INTEGER(KIND=i4), INTENT(OUT) :: idate(4)
    INTEGER(KIND=i4), INTENT(OUT) :: idatec(4)
    REAL   (KIND=r4), INTENT(OUT) :: si(:)
    REAL   (KIND=r4), INTENT(OUT) :: sl(:)
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadHead4)**"
    READ(n)ifday, tod, idate, idatec, si, sl
  END SUBROUTINE ReadHead4
  SUBROUTINE ReadHead8(n, ifday, tod, idate, idatec, si, sl)
    INTEGER(KIND=i8), INTENT(IN)  :: n
    INTEGER(KIND=i8), INTENT(OUT) :: ifday
    REAL   (KIND=r8), INTENT(OUT) :: tod
    INTEGER(KIND=i8), INTENT(OUT) :: idate(4)
    INTEGER(KIND=i8), INTENT(OUT) :: idatec(4)
    REAL   (KIND=r8), INTENT(OUT) :: si(:)
    REAL   (KIND=r8), INTENT(OUT) :: sl(:)
    INTEGER(KIND=i4) :: iaux(10)
    REAL   (KIND=r4) :: raux1(kmax), raux2(kmax+1)
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadHead8)**"
    READ(n)iaux, raux2, raux1
    ifday  = INT(iaux(  1 ), i8)
    tod    = INT(iaux(  2 ), r8)
    idate  = INT(iaux(3:6 ), i8)
    idatec = INT(iaux(7:10), i8)
    si     = REAL(raux2 , r8)
    sl     = REAL(raux1 , r8)
  END SUBROUTINE ReadHead8
  SUBROUTINE GReadHead4(n, ifday, tod, idate, idatec, si, sl)
    INTEGER(KIND=i8), INTENT(IN)  :: n
    INTEGER(KIND=i4), INTENT(OUT) :: ifday
    REAL   (KIND=r4), INTENT(OUT) :: tod
    INTEGER(KIND=i4), INTENT(OUT) :: idate(4)
    INTEGER(KIND=i4), INTENT(OUT) :: idatec(4)
    REAL   (KIND=r4), INTENT(OUT) :: si(:)
    REAL   (KIND=r4), INTENT(OUT) :: sl(:)
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadHead4)**"
    READ(n)ifday, tod, idate, idatec, si, sl
  END SUBROUTINE GReadHead4  
  SUBROUTINE GReadHead8(n, ifday, tod, idate, idatec, si, sl)
    INTEGER(KIND=i8), INTENT(IN)  :: n
    INTEGER(KIND=i8), INTENT(OUT) :: ifday
    REAL   (KIND=r8), INTENT(OUT) :: tod
    INTEGER(KIND=i8), INTENT(OUT) :: idate(4)
    INTEGER(KIND=i8), INTENT(OUT) :: idatec(4)
    REAL   (KIND=r8), INTENT(OUT) :: si(:)
    REAL   (KIND=r8), INTENT(OUT) :: sl(:)
    INTEGER(KIND=i8) :: iaux(10)
    REAL   (KIND=r8) :: raux1(kmax), raux2(kmax+1)
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadHead8)**"
    READ(n)iaux, raux2, raux1
    ifday  = INT(iaux(  1 ), i8)
    tod    = INT(iaux(  2 ), r8)
    idate  = INT(iaux(3:6 ), i8)
    idatec = INT(iaux(7:10), i8)
    si     = REAL(raux2 , r8)
    sl     = REAL(raux1 , r8)
  END SUBROUTINE GReadHead8
  
  SUBROUTINE ReadField42D(n, field)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    REAL   (KIND=r4), INTENT(OUT) :: field(:,:)
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadField42D)**"
    INTEGER :: k
    INTEGER :: d2
    d2 = SIZE(field,2)
    DO k = 1, d2
       READ(n)field(:,k)
    END DO
  END SUBROUTINE ReadField42D
  SUBROUTINE ReadField82D(n, field)
    INTEGER(KIND=i8), INTENT(IN)  :: n
    REAL   (KIND=r8), INTENT(OUT) :: field(:,:)
    REAL   (KIND=r4) :: raux3(SIZE(field,1))
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadField82D)**"
    INTEGER :: k
    INTEGER :: d2
    d2 = SIZE(field,2)
    DO k=1, d2
       READ(n)raux3
       field(:,k) = REAL(raux3, r8)
    END DO
  END SUBROUTINE ReadField82D
  SUBROUTINE ReadField41D(n, field)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    REAL   (KIND=r4), INTENT(OUT) :: field(:)
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadField41D)**"
    READ(n)field
  END SUBROUTINE ReadField41D
  SUBROUTINE ReadField81D(n, field)
    INTEGER(KIND=i8), INTENT(IN)  :: n
    REAL   (KIND=r8), INTENT(OUT) :: field(:)
    REAL   (KIND=r4) :: raux3(SIZE(field))
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadField81D)**"
    READ(n)raux3
    field = REAL(raux3, r8)
  END SUBROUTINE ReadField81D
  
  
  SUBROUTINE GReadField42D(n, field)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    REAL   (KIND=r4), INTENT(OUT) :: field(:,:)
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadField42D)**"
    INTEGER :: k
    INTEGER :: d2
    d2 = SIZE(field,2)
    DO k = 1, d2
       READ(n)field(:,k)
    END DO
  END SUBROUTINE GReadField42D
  SUBROUTINE GReadField82D(n, field)
    INTEGER(KIND=i8), INTENT(IN)  :: n
    REAL   (KIND=r8), INTENT(OUT) :: field(:,:)
    REAL   (KIND=r8) :: raux3(SIZE(field,1))
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadField82D)**"
    INTEGER :: k
    INTEGER :: d2
    d2 = SIZE(field,2)
    DO k=1, d2
       READ(n)raux3
       field(:,k) =raux3
    END DO
  END SUBROUTINE GReadField82D
  SUBROUTINE GReadField41D(n, field)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    REAL   (KIND=r4), INTENT(OUT) :: field(:)
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadField41D)**"
    READ(n)field
  END SUBROUTINE GReadField41D
  SUBROUTINE GReadField81D(n, field)
    INTEGER(KIND=i8), INTENT(IN)  :: n
    REAL   (KIND=r8), INTENT(OUT) :: field(:)
    REAL   (KIND=r8) :: raux3(SIZE(field))
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadField81D)**"
    READ(n)field
  END SUBROUTINE GReadField81D
  
  SUBROUTINE WriteHead4(n, ifday, tod, idate, idatec, si, sl)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    INTEGER(KIND=i4), INTENT(IN)  :: ifday
    REAL   (KIND=r4), INTENT(IN)  :: tod
    INTEGER(KIND=i4), INTENT(IN)  :: idate(4)
    INTEGER(KIND=i4), INTENT(IN)  :: idatec(4)
    REAL   (KIND=r4), INTENT(IN)  :: si(:)
    REAL   (KIND=r4), INTENT(IN)  :: sl(:)
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteHead4)**"
    WRITE(n)ifday, tod, idate, idatec, si, sl
  END SUBROUTINE WriteHead4
  SUBROUTINE WriteHead8(n, ifday, tod, idate, idatec, si, sl)
    INTEGER(KIND=i8), INTENT(IN)  :: n
    INTEGER(KIND=i8), INTENT(IN)  :: ifday
    REAL   (KIND=r8), INTENT(IN)  :: tod
    INTEGER(KIND=i8), INTENT(IN)  :: idate(4)
    INTEGER(KIND=i8), INTENT(IN)  :: idatec(4)
    REAL   (KIND=r8), INTENT(IN)  :: si(:)
    REAL   (KIND=r8), INTENT(IN)  :: sl(:)
    INTEGER(KIND=i4) :: iaux(10)
    REAL   (KIND=r4) :: raux1(kmax), raux2(kmax+1)
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteHead8)**"
    iaux(  1 ) = INT (ifday , i4)
    iaux(  2 ) = INT (tod   , i4)
    iaux(3:6 ) = INT (idate , i4)
    iaux(7:10) = INT (idatec, i4)
    raux2      = REAL(si    , r4)
    raux1      = REAL(sl    , r4)
    WRITE(n)iaux, raux2, raux1
  END SUBROUTINE WriteHead8

  SUBROUTINE GWriteHead4(n, ifday, tod, idate, idatec, si, sl)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    INTEGER(KIND=i4), INTENT(IN)  :: ifday
    REAL   (KIND=r4), INTENT(IN)  :: tod
    INTEGER(KIND=i4), INTENT(IN)  :: idate(4)
    INTEGER(KIND=i4), INTENT(IN)  :: idatec(4)
    REAL   (KIND=r4), INTENT(IN)  :: si(:)
    REAL   (KIND=r4), INTENT(IN)  :: sl(:)
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteHead4)**"
    WRITE(n)ifday, tod, idate, idatec, si, sl
  END SUBROUTINE GWriteHead4
  SUBROUTINE GWriteHead8(n, ifday, tod, idate, idatec, si, sl)
    INTEGER(KIND=i8), INTENT(IN)  :: n
    INTEGER(KIND=i8), INTENT(IN)  :: ifday
    REAL   (KIND=r8), INTENT(IN)  :: tod
    INTEGER(KIND=i8), INTENT(IN)  :: idate(4)
    INTEGER(KIND=i8), INTENT(IN)  :: idatec(4)
    REAL   (KIND=r8), INTENT(IN)  :: si(:)
    REAL   (KIND=r8), INTENT(IN)  :: sl(:)
    INTEGER(KIND=i8) :: iaux(10)
    REAL   (KIND=r8) :: raux1(kmax), raux2(kmax+1)
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteHead8)**"
    iaux(  1 ) = INT (ifday , i8)
    iaux(  2 ) = INT (tod   , r8)
    iaux(3:6 ) = INT (idate , i8)
    iaux(7:10) = INT (idatec, i8)
    raux2      = REAL(si    , r8)
    raux1      = REAL(sl    , r8)
    WRITE(n)iaux, raux2, raux1
  END SUBROUTINE GWriteHead8
  
  SUBROUTINE WriteField42D(n, field)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    REAL   (KIND=r4), INTENT(IN)  :: field(:,:)
    REAL   (KIND=r4) :: raux3(SIZE(field,1))
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteField42D)**"
    REAL :: MaxIEEE, MinIEEE
    INTEGER :: k, l
    INTEGER :: d1, d2
    MaxIEEE = 10.0**MaxIEEEExpo
    MinIEEE = 10.0**MinIEEEExpo
    d2=SIZE(field,2);d1=SIZE(field,1)
    DO k = 1, d2
       DO l = 1, d1
          IF (field(l,k) >= 0.0) THEN
             raux3(l) = MAX(MIN(field(l,k),MaxIEEE),MinIEEE)
          ELSE
             raux3(l) = MIN(MAX(field(l,k),-MaxIEEE),-MinIEEE)
          END IF
       END DO
       WRITE(n)raux3(:)
    END DO
  END SUBROUTINE WriteField42D
  SUBROUTINE WriteField82D(n, field)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    REAL   (KIND=r8), INTENT(IN)  :: field(:,:)
    REAL   (KIND=r4) :: raux3(SIZE(field,1))
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteField82D)**"
    REAL :: MaxIEEE, MinIEEE
    INTEGER :: k, l
    INTEGER :: d1, d2
    MaxIEEE = 10.0**MaxIEEEExpo
    MinIEEE = 10.0**MinIEEEExpo
    d2=SIZE(field,2);d1=SIZE(field,1)
    DO k = 1, d2
       DO l = 1, d1
          IF (field(l,k) >= 0.0) THEN
             raux3(l) = MAX(MIN(field(l,k),MaxIEEE),MinIEEE)
          ELSE
             raux3(l) = MIN(MAX(field(l,k),-MaxIEEE),-MinIEEE)
          END IF
       END DO
       WRITE(n)raux3(:)
    END DO
  END SUBROUTINE WriteField82D


  SUBROUTINE WriteField41D(n, field)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    REAL   (KIND=r4), INTENT(IN)  :: field(:)
    REAL   (KIND=r4) :: raux3(SIZE(field,1))
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteField41D)**"
    REAL :: MaxIEEE, MinIEEE
    INTEGER :: l
    INTEGER :: d1
    MaxIEEE = 10.0**MaxIEEEExpo
    MinIEEE = 10.0**MinIEEEExpo
    d1=SIZE(field,1)
    DO l = 1, d1
       IF (field(l) >= 0.0) THEN
          raux3(l) = MAX(MIN(field(l),MaxIEEE),MinIEEE)
       ELSE
          raux3(l) = MIN(MAX(field(l),-MaxIEEE),-MinIEEE)
       END IF
    END DO
    WRITE(n)raux3(:)
  END SUBROUTINE WriteField41D
  SUBROUTINE WriteField81D(n, field)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    REAL   (KIND=r8), INTENT(IN)  :: field(:)
    REAL   (KIND=r4) :: raux3(SIZE(field,1))
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteField81D)**"
    REAL :: MaxIEEE, MinIEEE
    INTEGER :: l
    INTEGER :: d1
    MaxIEEE = 10.0**MaxIEEEExpo
    MinIEEE = 10.0**MinIEEEExpo
    d1=SIZE(field,1)
    DO l = 1, d1
       IF (field(l) >= 0.0) THEN
          raux3(l) = MAX(MIN(field(l),MaxIEEE),MinIEEE)
       ELSE
          raux3(l) = MIN(MAX(field(l),-MaxIEEE),-MinIEEE)
       END IF
    END DO
    WRITE(n)raux3(:)
  END SUBROUTINE WriteField81D


  SUBROUTINE GWriteField42D(n, field)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    REAL   (KIND=r4), INTENT(IN)  :: field(:,:)
    REAL   (KIND=r4) :: raux3(SIZE(field,1))
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteField42D)**"
    REAL :: MaxIEEE, MinIEEE
    INTEGER :: k, l
    INTEGER :: d1, d2
    MaxIEEE = 10.0**MaxIEEEExpo
    MinIEEE = 10.0**MinIEEEExpo
    d2=SIZE(field,2);d1=SIZE(field,1)
    DO k = 1, d2
       DO l = 1, d1
          IF (field(l,k) >= 0.0) THEN
             raux3(l) = MAX(MIN(field(l,k),MaxIEEE),MinIEEE)
          ELSE
             raux3(l) = MIN(MAX(field(l,k),-MaxIEEE),-MinIEEE)
          END IF
       END DO
       WRITE(n)raux3(:)
    END DO
  END SUBROUTINE GWriteField42D
  SUBROUTINE GWriteField82D(n, field)
    INTEGER(KIND=i8), INTENT(IN)  :: n
    REAL   (KIND=r8), INTENT(IN)  :: field(:,:)
    REAL   (KIND=r8) :: raux3(SIZE(field,1))
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteField82D)**"
    REAL :: MaxIEEE, MinIEEE
    INTEGER :: k, l
    INTEGER :: d1, d2
    MaxIEEE = 10.0**MaxIEEEExpo
    MinIEEE = 10.0**MinIEEEExpo
    d2=SIZE(field,2);d1=SIZE(field,1)
    DO k = 1, d2
       DO l = 1, d1
          IF (field(l,k) >= 0.0) THEN
             raux3(l) = MAX(MIN(field(l,k),MaxIEEE),MinIEEE)
          ELSE
             raux3(l) = MIN(MAX(field(l,k),-MaxIEEE),-MinIEEE)
          END IF
       END DO
       WRITE(n)raux3(:)
    END DO
  END SUBROUTINE GWriteField82D


  SUBROUTINE GWriteField41D(n, field)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    REAL   (KIND=r4), INTENT(IN)  :: field(:)
    REAL   (KIND=r4) :: raux3(SIZE(field,1))
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteField41D)**"
    REAL :: MaxIEEE, MinIEEE
    INTEGER :: l
    INTEGER :: d1
    MaxIEEE = 10.0**MaxIEEEExpo
    MinIEEE = 10.0**MinIEEEExpo
    d1=SIZE(field,1)
    DO l = 1, d1
       IF (field(l) >= 0.0) THEN
          raux3(l) = MAX(MIN(field(l),MaxIEEE),MinIEEE)
       ELSE
          raux3(l) = MIN(MAX(field(l),-MaxIEEE),-MinIEEE)
       END IF
    END DO
    WRITE(n)raux3(:)
  END SUBROUTINE GWriteField41D
  SUBROUTINE GWriteField81D(n, field)
    INTEGER(KIND=i8), INTENT(IN)  :: n
    REAL   (KIND=r8), INTENT(IN)  :: field(:)
    REAL   (KIND=r8) :: raux3(SIZE(field,1))
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteField81D)**"
    REAL :: MaxIEEE, MinIEEE
    INTEGER :: l
    INTEGER :: d1
    MaxIEEE = 10.0**MaxIEEEExpo
    MinIEEE = 10.0**MinIEEEExpo
    d1=SIZE(field,1)
    DO l = 1, d1
       IF (field(l) >= 0.0) THEN
          raux3(l) = MAX(MIN(field(l),MaxIEEE),MinIEEE)
       ELSE
          raux3(l) = MIN(MAX(field(l),-MaxIEEE),-MinIEEE)
       END IF
    END DO
    WRITE(n)raux3(:)
  END SUBROUTINE GWriteField81D
  
  SUBROUTINE WriteDir4(n, idate, ihr, iday, mon, iyr, del,tod)
    INTEGER (KIND=i4), INTENT(IN ) :: n
    INTEGER (KIND=i4), INTENT(IN ) :: idate(4)
    INTEGER (KIND=i4), INTENT(IN ) :: ihr
    INTEGER (KIND=i4), INTENT(IN ) :: iday
    INTEGER (KIND=i4), INTENT(IN ) :: mon
    INTEGER (KIND=i4), INTENT(IN ) :: iyr
    REAL    (KIND=r4), INTENT(IN ) :: del(kMax)
    REAL    (KIND=r4), INTENT(IN ) :: tod
    INTEGER (KIND=i4)              :: isg(2)
    
    CHARACTER (LEN= 4) :: imdl
    CHARACTER (LEN=40) :: jttl
    CHARACTER (LEN=20), PARAMETER :: ittl='CPTEC SIGMA VERS 2.0'
    CHARACTER (LEN= 4), PARAMETER :: nexp='0001'
    CHARACTER (LEN=40), PARAMETER :: orogra = 'TOPOGRAPHY' 
    CHARACTER (LEN=40), PARAMETER :: lseamk = 'LAND SEA MASK'
    CHARACTER (LEN=40), PARAMETER :: lnsurf = 'LN SURFACE PRESSURE'
    CHARACTER (LEN=40), PARAMETER :: divrgn = 'DIVERGENCE'
    CHARACTER (LEN=40), PARAMETER :: vortic = 'VORTICITY'
    CHARACTER (LEN=40), PARAMETER :: spechu = 'SPECIFIC HUMIDITY'
    CHARACTER (LEN=40), PARAMETER :: tempvi = 'VIRTUAL TEMPERATURE'
    CHARACTER (LEN=40), PARAMETER :: surfte = 'SURFACE TEMPERATURE'
    CHARACTER (LEN=40), PARAMETER :: srfrou = 'ROUGHNESS LENGTH'
    CHARACTER (LEN=40), PARAMETER :: deepte = 'DEEP SOIL TEMPERATURE'
    CHARACTER (LEN=40), PARAMETER :: stcnpy = 'STORAGE ON CANOPY'
    CHARACTER (LEN=40), PARAMETER :: stgrnd = 'STORAGE ON GROUND'
    CHARACTER (LEN=40), PARAMETER :: wt1soi = 'SOIL WETNESS OF SURFACE'
    CHARACTER (LEN=40), PARAMETER :: wt2soi = 'SOIL WETNESS OF ROOT ZONE'
    CHARACTER (LEN=40), PARAMETER :: wt3soi = 'SOIL WETNESS OF DRAINAGE ZONE' 
    CHARACTER (LEN=29), PARAMETER :: fmt1='(A40,2X,A4,2X,I8,3X,I4,4X,I3)'
    CHARACTER (LEN=4 ), PARAMETER :: diag='DIAG'
    INTEGER :: m
    INTEGER :: nn
    INTEGER :: ix

    isg(1)=iMax*jMax
    isg(2)=2*mnMax

    jttl='CPTEC AGCM REVIS 1.0 2000  T   L    COLD'
    WRITE (jttl(29:31), '(i3.3)') mMax-1
    WRITE (jttl(33:34), '(i2.2)') kMax
    WRITE (imdl, '(A1,I3.3)') 'T', mMax-1

    WRITE (n, '(A20)')   ittl
    WRITE (n, '(A4,1X,A4,1X,A4,1X,11I5,1X,A4)') &
         nexp, 'SEQU', imdl, mMax, kmax, kmax, &
         ihr, iday, mon, iyr, idate, 'TAPE'
    WRITE (n, '(A40)')   jttl
    WRITE (n, '(5E16.8)')   del
    WRITE (n, fmt1) orogra, 'FIXD', 2*mnMax, 1, 10
    WRITE (n, fmt1) lseamk, 'FIXD', ijmax, 1, 0
    WRITE (n, fmt1) lnsurf, 'PROG', 2*mnMax, 1, 142
    WRITE (n, fmt1) divrgn, 'PROG', 2*mnMax, kmax, 50
    WRITE (n, fmt1) vortic, 'PROG', 2*mnMax, kmax, 50
    WRITE (n, fmt1) spechu, 'PROG', 2*mnMax, kmax, 0
    WRITE (n, fmt1) tempvi, 'PROG', 2*mnMax, kmax, 40
    WRITE (n, fmt1) srfrou, 'PROG', ijmax, 1, 10
    WRITE (n, fmt1) surfte, 'PROG', ijmax, 1, 40
    WRITE (n, fmt1) deepte, 'PROG', ijmax, 1, 40
    WRITE (n, fmt1) stcnpy, 'PROG', ijmax, 1, 10
    WRITE (n, fmt1) stgrnd, 'PROG', ijmax, 1, 10
    WRITE (n, fmt1) wt1soi, 'PROG', ijmax, 1, 0
    WRITE (n, fmt1) wt2soi, 'PROG', ijmax, 1, 0
    WRITE (n, fmt1) wt3soi, 'PROG', ijmax, 1, 0

    IF (iday == 0 .and. tod == 0.0) RETURN
  
    DO m=1,mxavl
      IF (dodia(m) .and. (iavrq(m) > 0)) THEN
        nn=iavrq(m)
        WRITE(n,160)reqdg(nn),diag,isg(itavl(m)),lvrq(nn),nurq(nn)
        IF(ifprt(91) >= 1)WRITE(nfprt,161) reqdg(nn),diag, &
             isg(itavl(m)),lvrq(nn),nurq(nn)
      END IF
    END DO

    IF(icf.ne.0)THEN
      DO ix=1,icf
        WRITE(n,160)combf(ix),diag,isg(itcf(ix)),lvcf(ix),nucf(ix)
        IF(ifprt(91) >= 1)WRITE(nfprt,161) combf(ix),diag, &
                         isg(itcf(ix)),lvcf(ix),nucf(ix)
      END DO
    END IF
160 FORMAT(A40,2X,A4,2X,I8,3X,I4,4X,I3)
161 FORMAT(' ',A40,2X,A4,2X,I8,3X,I4,4X,I3)    
  END SUBROUTINE WriteDir4
  SUBROUTINE WriteDir8(n, idate, ihr, iday, mon, iyr, del,tod)
    INTEGER (KIND=i4), INTENT(IN ) :: n
    INTEGER (KIND=i4), INTENT(IN ) :: idate(4)
    INTEGER (KIND=i4), INTENT(IN ) :: ihr
    INTEGER (KIND=i4), INTENT(IN ) :: iday
    INTEGER (KIND=i4), INTENT(IN ) :: mon
    INTEGER (KIND=i4), INTENT(IN ) :: iyr
    REAL    (KIND=r8), INTENT(IN ) :: del(kMax) 
    REAL    (KIND=r8), INTENT(IN ) :: tod
    INTEGER (KIND=i4)              :: isg(2)

    CHARACTER (LEN= 4) :: imdl
    CHARACTER (LEN=40) :: jttl
    CHARACTER (LEN=20), PARAMETER :: ittl='CPTEC SIGMA VERS 2.0'
    CHARACTER (LEN= 4), PARAMETER :: nexp='0001'
    CHARACTER (LEN=40), PARAMETER :: orogra = 'TOPOGRAPHY' 
    CHARACTER (LEN=40), PARAMETER :: lseamk = 'LAND SEA MASK'
    CHARACTER (LEN=40), PARAMETER :: lnsurf = 'LN SURFACE PRESSURE'
    CHARACTER (LEN=40), PARAMETER :: divrgn = 'DIVERGENCE'
    CHARACTER (LEN=40), PARAMETER :: vortic = 'VORTICITY'
    CHARACTER (LEN=40), PARAMETER :: spechu = 'SPECIFIC HUMIDITY'
    CHARACTER (LEN=40), PARAMETER :: tempvi = 'VIRTUAL TEMPERATURE'
    CHARACTER (LEN=40), PARAMETER :: surfte = 'SURFACE TEMPERATURE'
    CHARACTER (LEN=40), PARAMETER :: srfrou = 'ROUGHNESS LENGTH'
    CHARACTER (LEN=40), PARAMETER :: deepte = 'DEEP SOIL TEMPERATURE'
    CHARACTER (LEN=40), PARAMETER :: stcnpy = 'STORAGE ON CANOPY'
    CHARACTER (LEN=40), PARAMETER :: stgrnd = 'STORAGE ON GROUND'
    CHARACTER (LEN=40), PARAMETER :: wt1soi = 'SOIL WETNESS OF SURFACE'
    CHARACTER (LEN=40), PARAMETER :: wt2soi = 'SOIL WETNESS OF ROOT ZONE'
    CHARACTER (LEN=40), PARAMETER :: wt3soi = 'SOIL WETNESS OF DRAINAGE ZONE' 
    CHARACTER (LEN=29), PARAMETER :: fmt1='(A40,2X,A4,2X,I8,3X,I4,4X,I3)'
    CHARACTER (LEN=4 ), PARAMETER :: diag='DIAG'
    INTEGER :: m
    INTEGER :: nn
    INTEGER :: ix

    isg(1)=iMax*jMax
    isg(2)=2*mnMax
    
    jttl='CPTEC AGCM REVIS 1.0 2000  T   L    COLD'
    WRITE (jttl(29:31), '(i3.3)') mMax-1
    WRITE (jttl(33:34), '(i2.2)') kMax
    WRITE (imdl, '(A1,I3.3)') 'T', mMax-1

    WRITE (n, '(A20)')   ittl
    WRITE (n, '(A4,1X,A4,1X,A4,1X,11I5,1X,A4)') &
         nexp, 'SEQU', imdl, mMax, kmax, kmax, &
         ihr, iday, mon, iyr, idate, 'TAPE'
    WRITE (n, '(A40)')   jttl
    WRITE (n, '(5E16.8)')   del
    WRITE (n, fmt1) orogra, 'FIXD', 2*mnMax, 1, 10
    WRITE (n, fmt1) lseamk, 'FIXD', ijmax, 1, 0
    WRITE (n, fmt1) lnsurf, 'PROG', 2*mnMax, 1, 142
    WRITE (n, fmt1) divrgn, 'PROG', 2*mnMax, kmax, 50
    WRITE (n, fmt1) vortic, 'PROG', 2*mnMax, kmax, 50
    WRITE (n, fmt1) spechu, 'PROG', 2*mnMax, kmax, 0
    WRITE (n, fmt1) tempvi, 'PROG', 2*mnMax, kmax, 40
    WRITE (n, fmt1) srfrou, 'PROG', ijmax, 1, 10
    WRITE (n, fmt1) surfte, 'PROG', ijmax, 1, 40
    WRITE (n, fmt1) deepte, 'PROG', ijmax, 1, 40
    WRITE (n, fmt1) stcnpy, 'PROG', ijmax, 1, 10
    WRITE (n, fmt1) stgrnd, 'PROG', ijmax, 1, 10
    WRITE (n, fmt1) wt1soi, 'PROG', ijmax, 1, 0
    WRITE (n, fmt1) wt2soi, 'PROG', ijmax, 1, 0
    WRITE (n, fmt1) wt3soi, 'PROG', ijmax, 1, 0

    IF (iday == 0 .and. tod == 0.0) RETURN
  
    DO m=1,mxavl
      IF (dodia(m) .and. (iavrq(m) > 0)) THEN
        nn=iavrq(m)
        WRITE(n,160)reqdg(nn),diag,isg(itavl(m)),lvrq(nn),nurq(nn)
        IF(ifprt(91) >= 1)WRITE(nfprt,161) reqdg(nn),diag, &
             isg(itavl(m)),lvrq(nn),nurq(nn)
      END IF
    END DO

    IF(icf.ne.0)THEN
      DO ix=1,icf
        WRITE(n,160)combf(ix),diag,isg(itcf(ix)),lvcf(ix),nucf(ix)
        IF(ifprt(91) >= 1)WRITE(nfprt,161) combf(ix),diag, &
                         isg(itcf(ix)),lvcf(ix),nucf(ix)
      END DO
    END IF
  
160 FORMAT(A40,2X,A4,2X,I8,3X,I4,4X,I3)
161 FORMAT(' ',A40,2X,A4,2X,I8,3X,I4,4X,I3)
  END SUBROUTINE WriteDir8
  SUBROUTINE WriteDire4(n, idate, ihr, iday, mon, iyr, del,tod)
    INTEGER (KIND=i4), INTENT(IN ) :: n
    INTEGER (KIND=i4), INTENT(IN ) :: idate(4)
    INTEGER (KIND=i4), INTENT(IN ) :: ihr
    INTEGER (KIND=i4), INTENT(IN ) :: iday
    INTEGER (KIND=i4), INTENT(IN ) :: mon
    INTEGER (KIND=i4), INTENT(IN ) :: iyr
    REAL    (KIND=r4), INTENT(IN ) :: del(kMax)
    REAL    (KIND=r4), INTENT(IN ) :: tod
    
    CHARACTER (LEN= 4) :: imdl
    CHARACTER (LEN=40) :: jttl
    CHARACTER (LEN=20), PARAMETER :: ittl='CPTEC SIGMA VERS 2.0'
    CHARACTER (LEN= 4), PARAMETER :: nexp='0001'
    CHARACTER (LEN=40), PARAMETER :: orogra = 'TOPOGRAPHY' 
    CHARACTER (LEN=40), PARAMETER :: lseamk = 'LAND SEA MASK'
    CHARACTER (LEN=40), PARAMETER :: lnsurf = 'LN SURFACE PRESSURE'
    CHARACTER (LEN=40), PARAMETER :: divrgn = 'DIVERGENCE'
    CHARACTER (LEN=40), PARAMETER :: vortic = 'VORTICITY'
    CHARACTER (LEN=40), PARAMETER :: spechu = 'SPECIFIC HUMIDITY'
    CHARACTER (LEN=40), PARAMETER :: tempvi = 'VIRTUAL TEMPERATURE'
    CHARACTER (LEN=40), PARAMETER :: surfte = 'SURFACE TEMPERATURE'
    CHARACTER (LEN=40), PARAMETER :: srfrou = 'ROUGHNESS LENGTH'
    CHARACTER (LEN=40), PARAMETER :: deepte = 'DEEP SOIL TEMPERATURE'
    CHARACTER (LEN=40), PARAMETER :: stcnpy = 'STORAGE ON CANOPY'
    CHARACTER (LEN=40), PARAMETER :: stgrnd = 'STORAGE ON GROUND'
    CHARACTER (LEN=40), PARAMETER :: wt1soi = 'SOIL WETNESS OF SURFACE'
    CHARACTER (LEN=40), PARAMETER :: wt2soi = 'SOIL WETNESS OF ROOT ZONE'
    CHARACTER (LEN=40), PARAMETER :: wt3soi = 'SOIL WETNESS OF DRAINAGE ZONE' 
    CHARACTER (LEN=29), PARAMETER :: fmt1='(A40,2X,A4,2X,I8,3X,I4,4X,I3)'

    jttl='CPTEC AGCM REVIS 1.0 2000  T   L    COLD'
    WRITE (jttl(29:31), '(i3.3)') mMax-1
    WRITE (jttl(33:34), '(i2.2)') kMax
    WRITE (imdl, '(A1,I3.3)') 'T', mMax-1

    WRITE (n, '(A20)')   ittl
    WRITE (n, '(A4,1X,A4,1X,A4,1X,11I5,1X,A4)') &
         nexp, 'SEQU', imdl, mMax, kmax, kmax, &
         ihr, iday, mon, iyr, idate, 'TAPE'
    WRITE (n, '(A40)')   jttl
    WRITE (n, '(5E16.8)')   del
    WRITE (n, fmt1) orogra, 'FIXD', 2*mnMax, 1, 10
    WRITE (n, fmt1) lseamk, 'FIXD', ijmax, 1, 0
    WRITE (n, fmt1) lnsurf, 'PROG', 2*mnMax, 1, 142
    WRITE (n, fmt1) divrgn, 'PROG', 2*mnMax, kmax, 50
    WRITE (n, fmt1) vortic, 'PROG', 2*mnMax, kmax, 50
    WRITE (n, fmt1) spechu, 'PROG', 2*mnMax, kmax, 0
    WRITE (n, fmt1) tempvi, 'PROG', 2*mnMax, kmax, 40
    WRITE (n, fmt1) srfrou, 'PROG', ijmax, 1, 10
    WRITE (n, fmt1) surfte, 'PROG', ijmax, 1, 40
    WRITE (n, fmt1) deepte, 'PROG', ijmax, 1, 40
    WRITE (n, fmt1) stcnpy, 'PROG', ijmax, 1, 10
    WRITE (n, fmt1) stgrnd, 'PROG', ijmax, 1, 10
    WRITE (n, fmt1) wt1soi, 'PROG', ijmax, 1, 0
    WRITE (n, fmt1) wt2soi, 'PROG', ijmax, 1, 0
    WRITE (n, fmt1) wt3soi, 'PROG', ijmax, 1, 0
160 FORMAT(A40,2X,A4,2X,I8,3X,I4,4X,I3)
161 FORMAT(' ',A40,2X,A4,2X,I8,3X,I4,4X,I3)    
  END SUBROUTINE WriteDire4
  SUBROUTINE WriteDire8(n, idate, ihr, iday, mon, iyr, del,tod)
    INTEGER (KIND=i4), INTENT(IN ) :: n
    INTEGER (KIND=i8), INTENT(IN ) :: idate(4)
    INTEGER (KIND=i8), INTENT(IN ) :: ihr
    INTEGER (KIND=i8), INTENT(IN ) :: iday
    INTEGER (KIND=i8), INTENT(IN ) :: mon
    INTEGER (KIND=i8), INTENT(IN ) :: iyr
    REAL    (KIND=r8), INTENT(IN ) :: del(kMax) 
    REAL    (KIND=r8), INTENT(IN ) :: tod

    CHARACTER (LEN= 4) :: imdl
    CHARACTER (LEN=40) :: jttl
    CHARACTER (LEN=20), PARAMETER :: ittl='CPTEC SIGMA VERS 2.0'
    CHARACTER (LEN= 4), PARAMETER :: nexp='0001'
    CHARACTER (LEN=40), PARAMETER :: orogra = 'TOPOGRAPHY' 
    CHARACTER (LEN=40), PARAMETER :: lseamk = 'LAND SEA MASK'
    CHARACTER (LEN=40), PARAMETER :: lnsurf = 'LN SURFACE PRESSURE'
    CHARACTER (LEN=40), PARAMETER :: divrgn = 'DIVERGENCE'
    CHARACTER (LEN=40), PARAMETER :: vortic = 'VORTICITY'
    CHARACTER (LEN=40), PARAMETER :: spechu = 'SPECIFIC HUMIDITY'
    CHARACTER (LEN=40), PARAMETER :: tempvi = 'VIRTUAL TEMPERATURE'
    CHARACTER (LEN=40), PARAMETER :: surfte = 'SURFACE TEMPERATURE'
    CHARACTER (LEN=40), PARAMETER :: srfrou = 'ROUGHNESS LENGTH'
    CHARACTER (LEN=40), PARAMETER :: deepte = 'DEEP SOIL TEMPERATURE'
    CHARACTER (LEN=40), PARAMETER :: stcnpy = 'STORAGE ON CANOPY'
    CHARACTER (LEN=40), PARAMETER :: stgrnd = 'STORAGE ON GROUND'
    CHARACTER (LEN=40), PARAMETER :: wt1soi = 'SOIL WETNESS OF SURFACE'
    CHARACTER (LEN=40), PARAMETER :: wt2soi = 'SOIL WETNESS OF ROOT ZONE'
    CHARACTER (LEN=40), PARAMETER :: wt3soi = 'SOIL WETNESS OF DRAINAGE ZONE' 
    CHARACTER (LEN=29), PARAMETER :: fmt1='(A40,2X,A4,2X,I8,3X,I4,4X,I3)'


    
    jttl='CPTEC AGCM REVIS 1.0 2000  T   L    COLD'
    WRITE (jttl(29:31), '(i3.3)') mMax-1
    WRITE (jttl(33:34), '(i2.2)') kMax
    WRITE (imdl, '(A1,I3.3)') 'T', mMax-1

    WRITE (n, '(A20)')   ittl
    WRITE (n, '(A4,1X,A4,1X,A4,1X,11I5,1X,A4)') &
         nexp, 'SEQU', imdl, mMax, kmax, kmax, &
         ihr, iday, mon, iyr, idate, 'TAPE'
    WRITE (n, '(A40)')   jttl
    WRITE (n, '(5E16.8)')   del
    WRITE (n, fmt1) orogra, 'FIXD', 2*mnMax, 1, 10
    WRITE (n, fmt1) lseamk, 'FIXD', ijmax, 1, 0
    WRITE (n, fmt1) lnsurf, 'PROG', 2*mnMax, 1, 142
    WRITE (n, fmt1) divrgn, 'PROG', 2*mnMax, kmax, 50
    WRITE (n, fmt1) vortic, 'PROG', 2*mnMax, kmax, 50
    WRITE (n, fmt1) spechu, 'PROG', 2*mnMax, kmax, 0
    WRITE (n, fmt1) tempvi, 'PROG', 2*mnMax, kmax, 40
    WRITE (n, fmt1) srfrou, 'PROG', ijmax, 1, 10
    WRITE (n, fmt1) surfte, 'PROG', ijmax, 1, 40
    WRITE (n, fmt1) deepte, 'PROG', ijmax, 1, 40
    WRITE (n, fmt1) stcnpy, 'PROG', ijmax, 1, 10
    WRITE (n, fmt1) stgrnd, 'PROG', ijmax, 1, 10
    WRITE (n, fmt1) wt1soi, 'PROG', ijmax, 1, 0
    WRITE (n, fmt1) wt2soi, 'PROG', ijmax, 1, 0
    WRITE (n, fmt1) wt3soi, 'PROG', ijmax, 1, 0
160 FORMAT(A40,2X,A4,2X,I8,3X,I4,4X,I3)
161 FORMAT(' ',A40,2X,A4,2X,I8,3X,I4,4X,I3)
  END SUBROUTINE WriteDire8
  SUBROUTINE ReadProgHead4(n, ifday, tod, idate, idatec)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    INTEGER(KIND=i4), INTENT(OUT) :: ifday
    REAL   (KIND=r4), INTENT(OUT) :: tod
    INTEGER(KIND=i4), INTENT(OUT) :: idate(4)
    INTEGER(KIND=i4), INTENT(OUT) :: idatec(4)
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadProgHead4)**"
    READ(n)ifday, tod, idate, idatec
  END SUBROUTINE ReadProgHead4
  SUBROUTINE ReadProgHead8(n, ifday, tod, idate, idatec)
    INTEGER(KIND=i8), INTENT(IN)  :: n
    INTEGER(KIND=i8), INTENT(OUT) :: ifday
    REAL   (KIND=r8), INTENT(OUT) :: tod
    INTEGER(KIND=i8), INTENT(OUT) :: idate(4)
    INTEGER(KIND=i8), INTENT(OUT) :: idatec(4)
    INTEGER(KIND=i4) :: iaux(10)
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadProgHead8)**"
    READ(n)iaux
    ifday  = INT(iaux(  1 ), i8)
    tod    = INT(iaux(  2 ), i8)
    idate  = INT(iaux(3:6 ), i8)
    idatec = INT(iaux(7:10), i8)
  END SUBROUTINE ReadProgHead8
  SUBROUTINE GReadProgHead4(n, ifday, tod, idate, idatec)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    INTEGER(KIND=i4), INTENT(OUT) :: ifday
    REAL   (KIND=r4), INTENT(OUT) :: tod
    INTEGER(KIND=i4), INTENT(OUT) :: idate(4)
    INTEGER(KIND=i4), INTENT(OUT) :: idatec(4)
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadProgHead4)**"
    READ(n)ifday, tod, idate, idatec
  END SUBROUTINE GReadProgHead4
  SUBROUTINE GReadProgHead8(n, ifday, tod, idate, idatec)
    INTEGER(KIND=i8), INTENT(IN)  :: n
    INTEGER(KIND=i8), INTENT(OUT) :: ifday
    REAL   (KIND=r8), INTENT(OUT) :: tod
    INTEGER(KIND=i8), INTENT(OUT) :: idate(4)
    INTEGER(KIND=i8), INTENT(OUT) :: idatec(4)
    INTEGER(KIND=i8) :: iaux(10)
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadProgHead8)**"
    READ(n)iaux
    ifday  = INT(iaux(  1 ), i8)
    tod    = INT(iaux(  2 ), i8)
    idate  = INT(iaux(3:6 ), i8)
    idatec = INT(iaux(7:10), i8)
  END SUBROUTINE GReadProgHead8
  SUBROUTINE WriteProgHead4(n, ifday, tod, idate, idatec)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    INTEGER(KIND=i4), INTENT(IN)  :: ifday
    REAL   (KIND=r4), INTENT(IN)  :: tod
    INTEGER(KIND=i4), INTENT(IN)  :: idate(4)
    INTEGER(KIND=i4), INTENT(IN)  :: idatec(4)
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteProgHead4)**"
    WRITE(n)ifday, tod, idate, idatec
  END SUBROUTINE WriteProgHead4
  SUBROUTINE WriteProgHead8(n, ifday, tod, idate, idatec)
    INTEGER(KIND=i8), INTENT(IN)  :: n
    INTEGER(KIND=i8), INTENT(IN)  :: ifday
    REAL   (KIND=r8), INTENT(IN)  :: tod
    INTEGER(KIND=i8), INTENT(IN)  :: idate(4)
    INTEGER(KIND=i8), INTENT(IN)  :: idatec(4)
    INTEGER(KIND=i4) :: iaux(10)
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteProgHead8)**"
    iaux(  1 ) = INT (ifday , i4)
    iaux(  2 ) = INT (tod   , i4)
    iaux(3:6 ) = INT (idate , i4)
    iaux(7:10) = INT (idatec, i4)
    WRITE(n)iaux
  END SUBROUTINE WriteProgHead8
  SUBROUTINE ReadLandSeaMask4(n, iMax, lsmk)
    INTEGER(KIND=i4), INTENT(IN ) :: n
    INTEGER(KIND=i4), INTENT(IN ) :: iMax
    REAL(KIND=r4),    INTENT(OUT) :: lsmk(:)
    INTEGER(KIND=i4)              :: int_lsmk(SIZE(lsmk))
    CHARACTER (LEN=7) :: fmti
    INTEGER :: l
    WRITE (fmti, '(''('',i3,''i1)'')') iMax
    READ (n, fmti) int_lsmk
    lsmk(:)=REAL(1-2*int_lsmk(:),r4)
  END SUBROUTINE ReadLandSeaMask4
  SUBROUTINE ReadLandSeaMask8(n, iMax, lsmk)
    INTEGER(KIND=i8), INTENT(IN ) :: n
    INTEGER(KIND=i8), INTENT(IN ) :: iMax
    REAL(KIND=r8),    INTENT(OUT) :: lsmk(:)
    INTEGER(KIND=i4)              :: int_lsmk(SIZE(lsmk))
    CHARACTER (LEN=7) :: fmti
    INTEGER :: l
    WRITE (fmti, '(''('',i3,''i1)'')') iMax
    READ (n, fmti) int_lsmk
    lsmk(:)=REAL(1-2*int_lsmk(:),r8)
  END SUBROUTINE ReadLandSeaMask8

  SUBROUTINE ReadLandSeaMask2_4(fname, lsmk)
    CHARACTER(LEN=*), INTENT(IN ) :: fname
    REAL(KIND=r4),    INTENT(OUT) :: lsmk(imax*jmax)
    INTEGER(KIND=i4) :: unit
    INTEGER(KIND=i4) :: ierr
    CALL GetUnit(unit)
    OPEN (unit, FILE=TRIM(fname), FORM="UNFORMATTED", ACCESS="SEQUENTIAL",&
         ACTION="READ", STATUS="OLD", IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(*,"('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(fname), ierr
       STOP "**(ERROR)**"
    END IF
    READ (unit, IOSTAT=ierr) lsmk
    IF (ierr /= 0) THEN
       WRITE(*,"('**(ERROR)** Read file ',a,' returned iostat=',i4)") &
            TRIM(fname), ierr
       STOP "**(ERROR)**"
    END IF
    CLOSE (unit)
  END SUBROUTINE ReadLandSeaMask2_4

  SUBROUTINE ReadLandSeaMask2_8(fname, lsmk)
    CHARACTER(LEN=*  ), INTENT(IN ) :: fname
    REAL     (KIND=r8),    INTENT(OUT) :: lsmk(imax*jmax)
    REAL     (KIND=r4) :: aux(imax*jmax)
    INTEGER  (KIND=i8) :: unit
    INTEGER  (KIND=i4) :: ierr
    CALL GetUnit(unit)
    OPEN (unit, FILE=TRIM(fname), FORM="UNFORMATTED", ACCESS="SEQUENTIAL",&
         ACTION="READ", STATUS="OLD", IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(*,"('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(fname), ierr
       STOP "**(ERROR)**"
    END IF
    READ (unit, IOSTAT=ierr) aux
    IF (ierr /= 0) THEN
       WRITE(*,"('**(ERROR)** Read file ',a,' returned iostat=',i4)") &
            TRIM(fname), ierr
       STOP "**(ERROR)**"
    END IF
    lsmk = REAL(aux,r8)
    CLOSE(unit)
  END SUBROUTINE ReadLandSeaMask2_8
SUBROUTINE  LandSeaMask4(ifsst,nfsst,labeli,intsst,sstlag,fNameSSTAOI,rlsm)
  IMPLICIT NONE
  INTEGER  (KIND=i8), INTENT(IN   )  :: ifsst
  INTEGER  (KIND=i8), INTENT(IN   )  :: nfsst
  INTEGER  (KIND=i8), INTENT(INOUT)  :: intsst
  REAL     (KIND=r8), INTENT(OUT  )  :: sstlag
  REAL     (KIND=r4), INTENT(OUT  )  :: rlsm(:,:)
  CHARACTER(LEN=*  ), INTENT(IN   )  :: labeli
  CHARACTER(LEN=*  ), INTENT(IN   )  :: fNameSSTAOI
  
  CHARACTER(LEN=10 )                 :: labelsi
  CHARACTER(LEN=10 )                 :: labelsj
  INTEGER  (KIND=i4)                 :: lrecl
  INTEGER  (KIND=i4)                 :: nsst
  REAL     (KIND=r8)                 :: dlag
  REAL     (KIND=r8),ALLOCATABLE     :: var4(:,:)  
  INTEGER  (KIND=i8)                 :: j
  INTEGER  (KIND=i8)                 :: jmax
  INTEGER  (KIND=i8)                 :: i  
  INTEGER  (KIND=i8)                 :: imax
  jmax = SIZE(rlsm,2)   
  imax = SIZE(rlsm,1)
  ALLOCATE(var4(imax,jmax))


    ! Use open statement to open direct access file when ifsst .ge. 4:

    INQUIRE (IOLENGTH=lrecl) rlsm
    lrecl=lrecl/2
    
    OPEN (nfsst,file=TRIM(fNameSSTAOI),ACCESS='DIRECT',FORM='UNFORMATTED',&
          RECL=lrecl,STATUS='UNKNOWN')
    READ (nfsst,REC=1)nsst,labelsi,labelsj
    WRITE(*,*)nsst,labelsi,labelsj
    CALL daylag (labelsi,labelsj,dlag,intsst)
         intsst=NINT(dlag)
    IF (intsst > 10) intsst=-intsst
    CALL daylag (labelsi,labeli,dlag,intsst)
         sstlag=dlag
    WRITE(*,'(/,a)')' Direct Access SST File:'
    WRITE(*,'(1x,a4,3(1x,a2))') &
               labelsi(1:4),labelsi(5:6),labelsi(7:8),labelsi(9:10)
    WRITE(*,'(1x,a4,3(1x,a2))') &
               labelsj(1:4),labelsj(5:6),labelsj(7:8),labelsj(9:10)
    WRITE(*,'(i7,a)')intsst,'  Days'
    WRITE(*,'(1x,a4,3(1x,a2))') &
               labelsi(1:4),labelsi(5:6),labelsi(7:8),labelsi(9:10)
    WRITE(*,'(1x,a4,3(1x,a2))') &
               labeli(1:4),labeli(5:6),labeli(7:8),labeli(9:10)
    IF (intsst > 0) THEN
       WRITE(*,'(f7.1,a,/)')sstlag,'  Days'
    ELSE
       WRITE(*,'(f7.1,a,/)')sstlag,'  Months'
    END IF
    IF (sstlag < 0) THEN
        WRITE (*, 336) ifsst, sstlag
        STOP 336
    END IF
    READ (nfsst,REC=2) var4
    CLOSE (nfsst)

  DO j=1,jmax
      DO i=1,imax
         rlsm(i,j)=var4(i,j)
      END DO
  END DO
  DEALLOCATE(var4)
336 FORMAT(' FOR IFSST=',I5,' SSTLAG MUST BE SET NONNEGATIVE.  NOT ',G12.5)  
END SUBROUTINE  LandSeaMask4
SUBROUTINE  LandSeaMask8(ifsst,nfsst,labeli,intsst,sstlag,fNameSSTAOI,rlsm)
  IMPLICIT NONE
  INTEGER  (KIND=i8), INTENT(IN   )  :: ifsst
  INTEGER  (KIND=i8), INTENT(IN   )  :: nfsst
  INTEGER  (KIND=i8), INTENT(INOUT)  :: intsst
  REAL     (KIND=r8), INTENT(OUT  )  :: sstlag
  REAL     (KIND=r8), INTENT(OUT  )  :: rlsm(:,:)
  CHARACTER(LEN=*  ), INTENT(IN   )  :: labeli
  CHARACTER(LEN=*  ), INTENT(IN   )  :: fNameSSTAOI
  
  CHARACTER(LEN=10 )                 :: labelsi
  CHARACTER(LEN=10 )                 :: labelsj
  INTEGER  (KIND=i4)                 :: lrecl
  INTEGER  (KIND=i4)                 :: nsst
  REAL     (KIND=r8)                 :: dlag
  REAL     (KIND=r4),ALLOCATABLE     :: var4(:,:)  
  INTEGER  (KIND=i8)                 :: j
  INTEGER  (KIND=i8)                 :: jmax
  INTEGER  (KIND=i8)                 :: i  
  INTEGER  (KIND=i8)                 :: imax
  jmax = SIZE(rlsm,2)   
  imax = SIZE(rlsm,1)
  ALLOCATE(var4(imax,jmax))


    ! Use open statement to open direct access file when ifsst .ge. 4:

    INQUIRE (IOLENGTH=lrecl) rlsm
    lrecl=lrecl/2
    
    OPEN (nfsst,file=TRIM(fNameSSTAOI),ACCESS='DIRECT',FORM='UNFORMATTED',&
          RECL=lrecl,STATUS='UNKNOWN')
    READ (nfsst,REC=1)nsst,labelsi,labelsj
    WRITE(*,*)nsst,labelsi,labelsj
    CALL daylag (labelsi,labelsj,dlag,intsst)
         intsst=NINT(dlag)
    IF (intsst > 10) intsst=-intsst
    CALL daylag (labelsi,labeli,dlag,intsst)
         sstlag=dlag
    WRITE(*,'(/,a)')' Direct Access SST File:'
    WRITE(*,'(1x,a4,3(1x,a2))') &
               labelsi(1:4),labelsi(5:6),labelsi(7:8),labelsi(9:10)
    WRITE(*,'(1x,a4,3(1x,a2))') &
               labelsj(1:4),labelsj(5:6),labelsj(7:8),labelsj(9:10)
    WRITE(*,'(i7,a)')intsst,'  Days'
    WRITE(*,'(1x,a4,3(1x,a2))') &
               labelsi(1:4),labelsi(5:6),labelsi(7:8),labelsi(9:10)
    WRITE(*,'(1x,a4,3(1x,a2))') &
               labeli(1:4),labeli(5:6),labeli(7:8),labeli(9:10)
    IF (intsst > 0) THEN
       WRITE(*,'(f7.1,a,/)')sstlag,'  Days'
    ELSE
       WRITE(*,'(f7.1,a,/)')sstlag,'  Months'
    END IF
    IF (sstlag < 0) THEN
        WRITE (*, 336) ifsst, sstlag
        STOP 336
    END IF
    READ (nfsst,REC=2) var4
    CLOSE (nfsst)

  DO j=1,jmax
      DO i=1,imax
         rlsm(i,j)=var4(i,j)
      END DO
  END DO
  DEALLOCATE(var4)
336 FORMAT(' FOR IFSST=',I5,' SSTLAG MUST BE SET NONNEGATIVE.  NOT ',G12.5)  
END SUBROUTINE  LandSeaMask8



SUBROUTINE daylag (labeli,labelf,dlag,intsst)
  IMPLICIT NONE
  CHARACTER (LEN=*  ) :: labeli
  CHARACTER (LEN=*  ) :: labelf
  REAL      (KIND=r8) :: dlag
  REAL                :: xday
  REAL                :: yday
  INTEGER   (KIND=i8) :: intsst
  INTEGER :: yi,mi,di,hi,yf,mf,df,hf,ndy,y,n,ndi,ndf
  INTEGER, DIMENSION (12) :: ndm = &
           (/31,28,31,30,31,30,31,31,30,31,30,31/)
  INTEGER, DIMENSION (12) :: ndmi = &
           (/31,28,31,30,31,30,31,31,30,31,30,31/)
  INTEGER, DIMENSION (12) :: ndmf = &
           (/31,28,31,30,31,30,31,31,30,31,30,31/)

  READ (labeli(1:4), '(i4)') yi
  READ (labeli(5:6), '(i2)') mi
  READ (labeli(7:8), '(i2)') di
  READ (labeli(9:10),'(i2)') hi
  READ (labelf(1:4), '(i4)') yf
  READ (labelf(5:6), '(i2)') mf
  READ (labelf(7:8), '(i2)') df
  READ (labelf(9:10),'(i2)') hf
  IF (MOD(yi,4) .EQ. 0) ndmi(2)=29
  IF (MOD(yf,4) .EQ. 0) ndmf(2)=29

  IF (intsst > 0) THEN

    ndy=0
    DO y=yi+1,yf-1
       DO n=1,12
          ndy=ndy+ndm(n)
       END DO
       IF (MOD(y,4) .EQ. 0) ndy=ndy+1
    END DO

    ndi=di
    DO n=1,mi-1
       ndi=ndi+ndmi(n)
    END DO
    ndf=df
    DO n=1,mf-1
       ndf=ndf+ndmf(n)
    END DO

    IF (yf .EQ. yi) THEN
       dlag=REAL(ndf-ndi)+REAL(hf-hi)/24.0
    ELSE IF(yf .GT. yi) THEN
       ndi=365-ndi
       IF (ndmi(2) .EQ. 29) ndi=ndi+1
       dlag=REAL(ndf+ndi)+REAL(hf-hi)/24.0+ndy
    ELSE
       dlag=-1.0
    END IF

  ELSE

    IF (mf >= mi) THEN
      dlag=mf-mi+12*(yf-yi)
    ELSE
      dlag=12+mf-mi+12*(yf-yi-1)
    END IF
    IF (MOD(yf,4) .EQ. 0) ndmf(2)=29
    xday=REAL(df)+REAL(hf)/24.0
    yday=1.0+REAL(ndmf(mf))/2.0
    IF (xday <= yday) dlag=dlag-1.0

  END IF

END SUBROUTINE daylag
  SUBROUTINE GetUnit4(unit)
    INTEGER(KIND=i4), INTENT(OUT) :: unit
    LOGICAL :: op
    INTEGER(KIND=i4), PARAMETER :: LBUnit=1   ! lower bound unit number
    INTEGER(KIND=i4), PARAMETER :: UBUnit=99  ! upper bound unit number
    DO unit = LBUnit, UBUnit
       INQUIRE(unit, OPENED=op)
       IF (.NOT. op) THEN
          EXIT
       END IF
    END DO
    IF (unit > UBUnit) THEN
       WRITE(*,"('**(ERROR)** All file units are opened')")
       STOP "**(ERROR)**"
    END IF
  END SUBROUTINE GetUnit4

  SUBROUTINE GetUnit8(unit)
    INTEGER(KIND=i8), INTENT(OUT) :: unit
    LOGICAL :: op
    INTEGER(KIND=i8), PARAMETER :: LBUnit=1   ! lower bound unit number
    INTEGER(KIND=i8), PARAMETER :: UBUnit=99  ! upper bound unit number
    DO unit = LBUnit, UBUnit
       INQUIRE(unit, OPENED=op)
       IF (.NOT. op) THEN
          EXIT
       END IF
    END DO
    IF (unit > UBUnit) THEN
       WRITE(*,"('**(ERROR)** All file units are opened')")
       STOP "**(ERROR)**"
    END IF
  END SUBROUTINE GetUnit8

  SUBROUTINE ReadVar4(nfvar,irec,var,fclose)
    INTEGER(KIND=i4) , INTENT(in   ) :: nfvar
    INTEGER(KIND=i4) , INTENT(IN   ) :: irec
    REAL   (KIND=r4) , INTENT(out  ) :: var (:,:)
    INTEGER           , INTENT(IN   ) :: fclose
    REAL   (KIND=r4)  :: var4(iMax,jMax)
    INTEGER(KIND=i4)  :: i
    INTEGER(KIND=i4)  :: j
     
    READ(UNIT=nfvar,rec=irec) var4
    DO j=1,jMax
       DO i=1,iMax
          var(i,j)=var4(i,j)
       END DO
    END DO
    IF(fclose == 0)CLOSE(UNIT=nfvar,STATUS='KEEP')
  END SUBROUTINE ReadVar4
  SUBROUTINE ReadVar8(nfvar,irec,var,fclose)
    INTEGER           , INTENT(in   ) :: nfvar
    INTEGER           , INTENT(IN   ) :: irec
    REAL     (KIND=r8), INTENT(out  ) :: var (:,:)  
    INTEGER           , INTENT(IN   ) :: fclose  
    REAL     (KIND=r4) :: var8(iMax,jMax)
    INTEGER            :: i
    INTEGER            :: j 
    READ(UNIT=nfvar,rec=irec) var8
    DO j=1,jMax
       DO i=1,iMax
          var(i,j)=REAL(var8(i,j),r8)     
       END DO
    END DO
    IF(fclose == 0) CLOSE(UNIT=nfvar,STATUS='KEEP') 
  END SUBROUTINE ReadVar8

  SUBROUTINE ReadVar4_txt(nfvar,irec,var,fclose,type2)
    INTEGER(KIND=i4) , INTENT(in   ) :: nfvar
    INTEGER(KIND=i4) , INTENT(IN   ) :: irec
    CHARACTER(LEN=*) , INTENT(IN   ) :: type2
    REAL   (KIND=r4) , INTENT(out  ) :: var (:,:)
    INTEGER           , INTENT(IN   ) :: fclose
    REAL   (KIND=r4)  :: var4(iMax,jMax)
    INTEGER(KIND=i4)  :: i
    INTEGER(KIND=i4)  :: j
    PRINT*,type2     
    DO i=1,irec
       READ(nfvar,*) var4
    END DO
    REWIND(nfvar)
    DO j=1,jMax
       DO i=1,iMax
          var(i,j)=var4(i,j)
       END DO
    END DO
    IF(fclose == 0)CLOSE(UNIT=nfvar,STATUS='KEEP')
  END SUBROUTINE ReadVar4_txt
  SUBROUTINE ReadVar8_txt(nfvar,irec,var,fclose,type2)
    INTEGER           , INTENT(in   ) :: nfvar
    INTEGER           , INTENT(IN   ) :: irec
    CHARACTER(LEN=*)  , INTENT(IN   ) :: type2
    REAL     (KIND=r8), INTENT(out  ) :: var (:,:)  
    INTEGER           , INTENT(IN   ) :: fclose  
    REAL     (KIND=r4) :: var8(iMax,jMax)
    INTEGER            :: i
    INTEGER            :: j 
    PRINT*,type2
    DO i=1,irec
       READ(nfvar,*) var8
    END DO
    REWIND(nfvar)
    DO j=1,jMax
       DO i=1,iMax
          var(i,j)=REAL(var8(i,j),r8)     
       END DO
    END DO
    IF(fclose == 0) CLOSE(UNIT=nfvar,STATUS='KEEP') 
  END SUBROUTINE ReadVar8_txt


 SUBROUTINE ReadAlb4(n,field)
   INTEGER  (kind=i4)         , INTENT(in   ) :: n
   REAL     (KIND=r4)         , INTENT(out  ) :: field  (:)
   REAL     (KIND=r8)                         :: raux3(SIZE(field))
   CHARACTER(LEN=*), PARAMETER :: h="**(ReadAlb4)**"
   READ(n) raux3
   field = REAL(raux3, r4)
 END SUBROUTINE ReadAlb4
 SUBROUTINE ReadAlb8(n,field)
   INTEGER  (kind=i4)         , INTENT(in   ) :: n
   REAL     (KIND=r8)         , INTENT(out  ) :: field  (:)
   REAL     (KIND=r8)                         :: raux3(SIZE(field))
   CHARACTER(LEN=*), PARAMETER :: h="**(ReadAlb8)**"   
   READ(n) raux3
   field = raux3
 END SUBROUTINE ReadAlb8 
 SUBROUTINE ReadAlb4_rec(n,irec,field)
    INTEGER(KIND=i4) , INTENT(in   ) :: n
    INTEGER(KIND=i4) , INTENT(IN   ) :: irec
    REAL   (KIND=r4) , INTENT(out  ) :: field  (:,:)
    REAL   (KIND=r8) :: raux3(SIZE(field,1),SIZE(field,2))
    INTEGER(KIND=i4) :: i,j
    INTEGER(KIND=i4) :: idim,jdim
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadAlb4)**"
    idim=SIZE(field,1) ;jdim=SIZE(field,2)
    READ(UNIT=n,rec=irec) raux3
    DO j=1,jdim
      DO i=1,idim
        field(i,j) = REAL(raux3(i,j), r4)
      END DO
    END DO  
  END SUBROUTINE ReadAlb4_rec
  SUBROUTINE ReadAlb8_rec(n,irec,field)
    INTEGER , INTENT(in   ) :: n
    INTEGER , INTENT(IN   ) :: irec
    REAL     (KIND=r8) , INTENT(out  ) :: field  (:,:)
    REAL     (KIND=r8) :: raux3(SIZE(field,1),SIZE(field,2))
    INTEGER(KIND=i4) :: i,j
    INTEGER(KIND=i4) :: idim,jdim
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadAlb8)**"   
    idim=SIZE(field,1) ;jdim=SIZE(field,2)
    READ(UNIT=n,rec=irec) raux3
    DO j=1,jdim
      DO i=1,idim
        field(i,j) = REAL(raux3(i,j), r8)
      END DO
    END DO  
  END SUBROUTINE ReadAlb8_rec


 SUBROUTINE ReadSST4(n,field)
   INTEGER  (kind=i4)         , INTENT(in   ) :: n
   REAL     (KIND=r4)         , INTENT(out  ) :: field  (:)
   REAL     (KIND=r8)                         :: raux3(SIZE(field))
   CHARACTER(LEN=*), PARAMETER :: h="**(ReadSST4)**"   
   READ(n) raux3
   field = REAL(raux3, r4)
 END SUBROUTINE ReadSST4
 SUBROUTINE ReadSST8(n,field)
    INTEGER(KIND=i4)          , INTENT(IN   ) :: n
    REAL   (KIND=r8)          , INTENT(OUT  ) :: field(:)
    REAL   (KIND=r4)                          :: raux3(SIZE(field))
   CHARACTER(LEN=*), PARAMETER :: h="**(ReadSST8)**"   
    READ(n)raux3
    field = REAL(raux3, r8)
 END SUBROUTINE ReadSST8 
  SUBROUTINE ReadSST4_rec(n,irec,field)
    INTEGER , INTENT(in   ) :: n
    INTEGER , INTENT(IN   ) :: irec
    REAL     (KIND=r4)         , INTENT(out  ) :: field  (:,:)
    REAL     (KIND=r4)                         :: raux3(SIZE(field,1),SIZE(field,2))
    INTEGER(KIND=i4)                         :: i,j
    INTEGER(KIND=i4)                         :: idim,jdim
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadSST4)**" 
    idim=SIZE(field,1) ;jdim=SIZE(field,2)
    READ (UNIT=n, REC=irec) raux3
    DO j=1,jdim
      DO i=1,idim
        field(i,j) = REAL(raux3(i,j), r4)
      END DO
    END DO  
  END SUBROUTINE ReadSST4_rec
  SUBROUTINE ReadSST8_rec(n,irec,field)
    INTEGER , INTENT(IN   ) :: n
    INTEGER , INTENT(IN   ) :: irec
    REAL   (KIND=r8)          , INTENT(OUT  ) :: field(:,:)
    REAL   (KIND=r4)                          :: raux3(SIZE(field,1),SIZE(field,2))
    INTEGER(KIND=i4)                          :: i,j
    INTEGER(KIND=i4)                          :: idim,jdim
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadSST8)**"   
    idim=SIZE(field,1) ;jdim=SIZE(field,2)
    READ (UNIT=n, REC=irec) raux3
    DO j=1,jdim
      DO i=1,idim
        field(i,j) = REAL(raux3(i,j), r8)
      END DO
    END DO  
  END SUBROUTINE ReadSST8_rec

  SUBROUTINE ReadSST4_txt(n,irec,field,type2)
    INTEGER , INTENT(in   ) :: n
    INTEGER , INTENT(IN   ) :: irec
    CHARACTER(LEN=*), INTENT(IN   ) ::   type2
    REAL     (KIND=r4)         , INTENT(out  ) :: field  (:,:)
    REAL     (KIND=r4)                         :: raux3(SIZE(field,1),SIZE(field,2))
    INTEGER(KIND=i4)                         :: i,j
    INTEGER(KIND=i4)                         :: idim,jdim
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadSST4)**" 
    idim=SIZE(field,1) ;jdim=SIZE(field,2)
    PRINT*,type2
    DO i=1,irec
       READ (n,*) raux3
    END DO   
    REWIND(n)
    DO j=1,jdim
      DO i=1,idim
        field(i,j) = REAL(raux3(i,j), r4)
      END DO
    END DO  
  END SUBROUTINE ReadSST4_txt
  SUBROUTINE ReadSST8_txt(n,irec,field,type2)
    INTEGER , INTENT(IN   ) :: n
    INTEGER , INTENT(IN   ) :: irec
    CHARACTER(LEN=*), INTENT(IN   ) ::   type2
    REAL   (KIND=r8)          , INTENT(OUT  ) :: field(:,:)
    REAL   (KIND=r4)                          :: raux3(SIZE(field,1),SIZE(field,2))
    INTEGER(KIND=i4)                          :: i,j
    INTEGER(KIND=i4)                          :: idim,jdim
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadSST8)**"   
    idim=SIZE(field,1) ;jdim=SIZE(field,2)
    PRINT*,type2
    DO i=1,irec
       READ (n,*) raux3
    END DO
    REWIND(n)   
    DO j=1,jdim
      DO i=1,idim
        field(i,j) = REAL(raux3(i,j), r8)
      END DO
    END DO  
  END SUBROUTINE ReadSST8_txt
 
 SUBROUTINE ReadSST4Rec(n,field,irec)
   INTEGER  (kind=i4)         , INTENT(in   ) :: n
   INTEGER  (kind=i4)         , INTENT(in   ) :: irec
   REAL     (KIND=r4)         , INTENT(out  ) :: field  (:)
   REAL     (KIND=r8)                         :: raux3(SIZE(field))
   READ(n,rec=irec) raux3
   field = REAL(raux3, r4)
 END SUBROUTINE ReadSST4Rec
 SUBROUTINE ReadSST8Rec(n,field,irec)
    INTEGER(KIND=i4)          , INTENT(IN   ) :: n
    INTEGER(kind=i4)          , INTENT(in   ) :: irec
    REAL   (KIND=r8)          , INTENT(OUT  ) :: field(:)
    REAL   (KIND=r4)                          :: raux3(SIZE(field))
    READ(n,rec=irec)raux3
    field = REAL(raux3, r8)
 END SUBROUTINE ReadSST8Rec

  SUBROUTINE ReadSST4Rec2(n,field,irec)
    INTEGER           , INTENT(in   ) :: n
    INTEGER           , INTENT(in   ) :: irec
    REAL     (KIND=r4)         , INTENT(out  ) :: field  (:,:)
    REAL     (KIND=r4)                         :: raux3(SIZE(field,1)*SIZE(field,2))
    INTEGER(KIND=i4)                         :: i,j,ij
    INTEGER(KIND=i4)                         :: idim,jdim
    idim=SIZE(field,1) ;jdim=SIZE(field,2)
    READ(UNIT=n,rec=irec) raux3
    ij=0
    DO j=1,jdim
      DO i=1,idim
        ij=ij+1
        field(i,j) = REAL(raux3(ij), r4)
      END DO
    END DO
  END SUBROUTINE ReadSST4Rec2
  SUBROUTINE ReadSST8Rec2(n,field,irec)
    INTEGER          , INTENT(IN   ) :: n
    INTEGER          , INTENT(in   ) :: irec
    REAL   (KIND=r8)          , INTENT(OUT  ) :: field(:,:)
    REAL   (KIND=r4)                          :: raux3(SIZE(field,1)*SIZE(field,2))
    INTEGER(KIND=i4)                          :: i,j,ij
    INTEGER(KIND=i4)                          :: idim,jdim
    idim=SIZE(field,1) ;jdim=SIZE(field,2)
    READ(UNIT=n,rec=irec)raux3
    ij=0
    DO j=1,jdim
      DO i=1,idim
        ij=ij+1    
        field(i,j) = REAL(raux3(ij), r8)
      END DO
    END DO  
  END SUBROUTINE ReadSST8Rec2


 SUBROUTINE ReadSLM4(n,field)
   INTEGER  (kind=i4)         , INTENT(in   ) :: n
   REAL     (KIND=r4)         , INTENT(out  ) :: field  (:) 
   REAL     (KIND=r4)                         :: raux3(SIZE(field))
   READ(n) raux3
   field = raux3
 END SUBROUTINE ReadSLM4
 SUBROUTINE ReadSLM8(n,field)
   INTEGER  (kind=i4)         , INTENT(in   ) :: n
   REAL     (KIND=r8)         , INTENT(out  ) :: field  (:)
   REAL     (KIND=r8)                         :: raux3(SIZE(field))
   READ(n) raux3
   field = raux3
 END SUBROUTINE ReadSLM8 

  SUBROUTINE ReadSLM4_rec(n,irec,field)
    INTEGER(KIND=i4)  , INTENT(in   ) :: n
    INTEGER           , INTENT(in   ) :: irec
    REAL     (KIND=r4), INTENT(out  ) :: field  (:,:) 
    REAL     (KIND=r4)                :: raux3(SIZE(field,1),SIZE(field,2))
    INTEGER(KIND=i4)                  :: i,j
    INTEGER(KIND=i4)                  :: idim,jdim
    idim=SIZE(field,1) ;jdim=SIZE(field,2)
    READ(UNIT=n,rec=irec) raux3
    DO j=1,jdim
      DO i=1,idim
        field(i,j) = raux3(i,j)
      END DO
    END DO  
  END SUBROUTINE ReadSLM4_rec
  SUBROUTINE ReadSLM8_rec(n,irec,field)
    INTEGER           , INTENT(in   ) :: n
    INTEGER           , INTENT(in   ) :: irec
    REAL     (KIND=r8), INTENT(out  ) :: field  (:,:)
    REAL     (KIND=r4)                :: raux3(SIZE(field,1),SIZE(field,2))
    INTEGER(KIND=i4)                  :: i,j
    INTEGER(KIND=i4)                  :: idim,jdim
    idim=SIZE(field,1) ;jdim=SIZE(field,2)
    READ(UNIT=n,rec=irec) raux3
    DO j=1,jdim
      DO i=1,idim
        field(i,j) = raux3(i,j)
      END DO
    END DO
  END SUBROUTINE ReadSLM8_rec

  SUBROUTINE ReadSLM4_txt(n,irec,field,type2)
    INTEGER(KIND=i4)  , INTENT(in   ) :: n
    INTEGER           , INTENT(in   ) :: irec
    CHARACTER(LEN=*)  , INTENT(IN   ) :: type2
    REAL     (KIND=r4), INTENT(out  ) :: field  (:,:) 
    REAL     (KIND=r4)                :: raux3(SIZE(field,1),SIZE(field,2))
    INTEGER(KIND=i4)                  :: i,j
    INTEGER(KIND=i4)                  :: idim,jdim
    idim=SIZE(field,1) ;jdim=SIZE(field,2)
    PRINT*, type2
    DO i=1,irec
       READ(n,*) raux3
    END DO
    REWIND(n)
    DO j=1,jdim
      DO i=1,idim
        field(i,j) = raux3(i,j)
      END DO
    END DO  
  END SUBROUTINE ReadSLM4_txt
  SUBROUTINE ReadSLM8_txt(n,irec,field,type2)
    INTEGER           , INTENT(in   ) :: n
    INTEGER           , INTENT(in   ) :: irec
    CHARACTER(LEN=*)  , INTENT(IN   ) :: type2
    REAL     (KIND=r8), INTENT(out  ) :: field  (:,:)
    REAL     (KIND=r4)                :: raux3(SIZE(field,1),SIZE(field,2))
    INTEGER(KIND=i4)                  :: i,j
    INTEGER(KIND=i4)                  :: idim,jdim
    idim=SIZE(field,1) ;jdim=SIZE(field,2)
    PRINT*, type2
    DO i=1,irec
       READ(n,*) raux3
    END DO
    REWIND(n)
    DO j=1,jdim
      DO i=1,idim
        field(i,j) = raux3(i,j)
      END DO
    END DO
  END SUBROUTINE ReadSLM8_txt

 SUBROUTINE ReadSLM43D(n,irec,field)
    INTEGER(KIND=i4)  , INTENT(in   ) :: n
    INTEGER           , INTENT(in   ) :: irec
    REAL     (KIND=r4), INTENT(out  ) :: field  (:,:,:) 
    REAL     (KIND=r4)                :: raux3(SIZE(field,1),SIZE(field,2),SIZE(field,3))
    INTEGER(KIND=i4)                  :: i,j,k
    INTEGER(KIND=i4)                  :: idim,jdim,kdim
    idim=SIZE(field,1) ;jdim=SIZE(field,2);kdim=SIZE(field,3)
    READ(UNIT=n,rec=irec) raux3
    DO k=1,kdim
       DO j=1,jdim
          DO i=1,idim
             field(i,j,k) = raux3(i,j,k)
          END DO
      END DO
    END DO  
  END SUBROUTINE ReadSLM43D
  SUBROUTINE ReadSLM83D(n,irec,field)
    INTEGER           , INTENT(in   ) :: n
    INTEGER           , INTENT(in   ) :: irec
    REAL     (KIND=r8), INTENT(out  ) :: field  (:,:,:)
    REAL     (KIND=r4)                :: raux3(SIZE(field,1),SIZE(field,2),SIZE(field,3))
    INTEGER(KIND=i4)                  :: i,j,k
    INTEGER(KIND=i4)                  :: idim,jdim,kdim
    idim=SIZE(field,1) ;jdim=SIZE(field,2);kdim=SIZE(field,3)
    READ(UNIT=n,rec=irec) raux3
    DO k=1,kdim
       DO j=1,jdim
          DO i=1,idim
             field(i,j,k) = raux3(i,j,k)
          END DO
      END DO
    END DO  
  END SUBROUTINE ReadSLM83D


 SUBROUTINE ReadSNW4(n,field)
   INTEGER  (kind=i4)         , INTENT(in   ) :: n
   REAL     (KIND=r4)         , INTENT(out  ) :: field  (:)
    REAL    (KIND=r8)                         :: raux3(SIZE(field))
    READ(n) raux3
    field = REAL(raux3, r4)
 END SUBROUTINE ReadSNW4
 SUBROUTINE ReadSNW8(n,field)
    INTEGER(KIND=i4),             INTENT(IN)  :: n
    REAL   (KIND=r8),             INTENT(OUT) :: field(:)
    REAL   (KIND=r4)                          :: raux3(SIZE(field))
    READ(n)raux3
    field = REAL(raux3, r8)
 END SUBROUTINE ReadSNW8 



  SUBROUTINE ReadSNW4_rec(n,irec,field)
    INTEGER           , INTENT(IN   ) :: n
    INTEGER           , INTENT(IN   ) :: irec
    REAL     (KIND=r4)         , INTENT(out  ) :: field  (:,:)
    REAL     (KIND=r4)                         :: raux3(SIZE(field,1),SIZE(field,2))
    INTEGER(KIND=i4)                         :: i,j
    INTEGER(KIND=i4)                         :: idim,jdim
    idim=SIZE(field,1) ;jdim=SIZE(field,2)
    READ(UNIT=n,rec=irec) raux3
    DO j=1,jdim
      DO i=1,idim
        field(i,j) = REAL(raux3(i,j), r4)
      END DO
    END DO  
  END SUBROUTINE ReadSNW4_rec
  SUBROUTINE ReadSNW8_rec(n,irec,field)
    INTEGER         , INTENT(IN   ) :: n
    INTEGER         , INTENT(IN   ) :: irec
    REAL   (KIND=r8), INTENT(OUT  ) :: field(:,:)
    REAL   (KIND=r4)                          :: raux3(SIZE(field,1),SIZE(field,2))
    INTEGER(KIND=i4)                          :: i,j
    INTEGER(KIND=i4)                          :: idim,jdim
    idim=SIZE(field,1) ;jdim=SIZE(field,2)
    READ(UNIT=n,rec=irec)raux3
    DO j=1,jdim
      DO i=1,idim
         field(i,j) = REAL(raux3(i,j), r8)
      END DO
    END DO
  END SUBROUTINE ReadSNW8_rec

  SUBROUTINE ReadSNW4_txt(n,irec,field,type2)
    INTEGER           , INTENT(IN   ) :: n
    INTEGER           , INTENT(IN   ) :: irec
    CHARACTER(LEN=*)  , INTENT(IN   ) :: type2
    REAL     (KIND=r4)         , INTENT(out  ) :: field  (:,:)
    REAL     (KIND=r4)                         :: raux3(SIZE(field,1),SIZE(field,2))
    INTEGER(KIND=i4)                         :: i,j
    INTEGER(KIND=i4)                         :: idim,jdim
    idim=SIZE(field,1) ;jdim=SIZE(field,2)
    PRINT*,type2
    DO i=1,irec
       READ(n,*) raux3
    END DO
    REWIND(n)
    DO j=1,jdim
      DO i=1,idim
        field(i,j) = REAL(raux3(i,j), r4)
      END DO
    END DO  
  END SUBROUTINE ReadSNW4_txt
  SUBROUTINE ReadSNW8_txt(n,irec,field,type2)
    INTEGER         , INTENT(IN   ) :: n
    INTEGER         , INTENT(IN   ) :: irec
    CHARACTER(LEN=*), INTENT(IN   ) :: type2
    REAL   (KIND=r8), INTENT(OUT  ) :: field(:,:)
    REAL   (KIND=r4)                          :: raux3(SIZE(field,1),SIZE(field,2))
    INTEGER(KIND=i4)                          :: i,j
    INTEGER(KIND=i4)                          :: idim,jdim
    idim=SIZE(field,1) ;jdim=SIZE(field,2)
    PRINT*,type2
    DO i=1,irec
       READ(n,*)raux3
    END DO
    REWIND(n)
    DO j=1,jdim
      DO i=1,idim
         field(i,j) = REAL(raux3(i,j), r8)
      END DO
    END DO
  END SUBROUTINE ReadSNW8_txt

  ! The Ozone files are written in a way that they can be read from GrADS.
  ! Therefore,  the order of the varibles inside the binary had to be x, y, z
  !   what is different from the global model x, z, y order.
  !   Moreover, the file is open in GrADS with a yrev option. So y=1 
  !   means the north pole as exepected by the global model.
  SUBROUTINE ReadOzone8(n,field,irec)
     INTEGER          , INTENT(IN   ) :: n,irec
     REAL   (KIND=r8) , INTENT(OUT  ) :: field(:,:,:)
     REAL   (KIND=r4)                 :: raux3(SIZE(field,1),SIZE(field,3),SIZE(field,2))
     CHARACTER(LEN=*), PARAMETER :: h="**(ReadOzone8)**"   
     INTEGER :: i,j,k,im,jm,km

     READ(UNIT=n,REC=irec) raux3
     ! input field is (i,k,j)
     im=size(field,1)
     jm=size(field,3)
     km=size(field,2)
     do j=1,jm
        do k=1,km
           do i=1,im
              field(i,k,j) = REAL(raux3(i,j,k), r8)
           enddo
        enddo
     enddo

  END SUBROUTINE ReadOzone8
  ! The tracer files are written in a way that they can be read from GrADS.
  ! Therefore,  the order of the varibles inside the binary had to be x, y, z
  !   what is different from the global model x, z, y order.
  !   Moreover, the file is open in GrADS with a yrev option. So y=1 
  !   means the north pole as exepected by the global model.
  SUBROUTINE ReadTracer8(n,field,irec)
     INTEGER          , INTENT(IN   ) :: n,irec
     REAL   (KIND=r8) , INTENT(OUT  ) :: field(:,:,:)
     REAL   (KIND=r4)                 :: raux3(SIZE(field,1),SIZE(field,3))
     CHARACTER(LEN=*), PARAMETER :: h="**(ReadOzone8)**"   
     INTEGER :: i,j,k,im,jm,km

     !READ(UNIT=n,REC=irec) raux3
     ! input field is (i,k,j)
     im=size(field,1)
     jm=size(field,3)
     km=size(field,2)
     do k=1,km
        READ(UNIT=n,REC=k) raux3
        PRINT*,'pkubota',k
        PRINT*,raux3
        do j=1,jm
           do i=1,im
              field(i,k,j) = REAL(raux3(i,j),kind= r8)
              IF(field(i,k,j).LT.0.0e0_r8)field(i,k,j) = 1.0e-12_r8
           enddo
        enddo
     enddo

  END SUBROUTINE ReadTracer8
  SUBROUTINE ReadTracer8s(n,field)
     INTEGER          , INTENT(IN   ) :: n
     REAL   (KIND=r8) , INTENT(OUT  ) :: field(:,:,:)
     REAL   (KIND=r4)                 :: raux3(SIZE(field,1),size(field,3))
     CHARACTER(LEN=*), PARAMETER :: h="**(ReadOzone8)**"
     INTEGER :: i,j,k,im,jm,km

     !READ(UNIT=n,REC=irec) raux3
     ! input field is (i,k,j)
     im=size(field,1)
     jm=size(field,3)
     km=size(field,2)
     do k=1,km
        READ(UNIT=n) raux3
        do j=1,jm
           do i=1,im
              field(i,k,j) = REAL(raux3(i,j),kind= r8)
              IF(field(i,k,j).LT.0.0e0_r8)field(i,k,j) = 1.0e-12_r8
           enddo
        enddo
     enddo

  END SUBROUTINE ReadTracer8s

 SUBROUTINE ReadNFTGZ4(n,field1,field2,field3,field4)
  INTEGER(kind=i8), INTENT(in   ) :: n 
  REAL   (kind=r8), INTENT(out  ) :: field1 (:,:)
  REAL   (kind=r8), INTENT(out  ) :: field2 (:,:)
  REAL   (kind=r8), INTENT(out  ) :: field3 (:,:)
  REAL   (kind=r8), INTENT(out  ) :: field4 (:,:)
  REAL   (KIND=r8)                :: raux1(SIZE(field1,1),SIZE(field1,2))
  REAL   (KIND=r8)                :: raux2(SIZE(field2,1),SIZE(field2,2))
  REAL   (KIND=r8)                :: raux3(SIZE(field3,1),SIZE(field3,2))
  REAL   (KIND=r8)                :: raux4(SIZE(field4,1),SIZE(field4,2))
  
  READ  (n) raux1,raux2,raux3,raux4
  field1 = REAL(raux1, r8)
  field2 = REAL(raux2, r8)
  field3 = REAL(raux3, r8)
  field4 = REAL(raux4, r8)

  REWIND n
 END SUBROUTINE ReadNFTGZ4
 SUBROUTINE ReadNFTGZ8(n,field1,field2,field3,field4)
  INTEGER(kind=i4), INTENT(in   ) :: n 
  REAL   (kind=r4), INTENT(out  ) :: field1 (:,:)
  REAL   (kind=r4), INTENT(out  ) :: field2 (:,:)
  REAL   (kind=r4), INTENT(out  ) :: field3 (:,:)
  REAL   (kind=r4), INTENT(out  ) :: field4 (:,:)
  REAL   (KIND=r8)                :: raux1(SIZE(field1,1),SIZE(field1,2))
  REAL   (KIND=r8)                :: raux2(SIZE(field2,1),SIZE(field2,2))
  REAL   (KIND=r8)                :: raux3(SIZE(field3,1),SIZE(field3,2))
  REAL   (KIND=r8)                :: raux4(SIZE(field4,1),SIZE(field4,2))
  
  READ  (n) raux1,raux2,raux3,raux4
  field1 = REAL(raux1, r4)
  field2 = REAL(raux2, r4)
  field3 = REAL(raux3, r4)
  field4 = REAL(raux4, r4)

  REWIND n
 END SUBROUTINE ReadNFTGZ8 
  SUBROUTINE ReadNFTGZ4_rec(n,irec,field1,field2,field3)
    INTEGER, INTENT(in   ) :: n 
    INTEGER, INTENT(in   ) :: irec
    REAL   (KIND=r8), INTENT(out  ) :: field1 (:,:)
    REAL   (KIND=r8), INTENT(out  ) :: field2 (:,:)
    REAL   (KIND=r8), INTENT(out  ) :: field3 (:,:)
    REAL   (KIND=r4)                :: raux1(SIZE(field1,1),SIZE(field1,2))
    REAL   (KIND=r4)                :: raux2(SIZE(field2,1),SIZE(field2,2))
    REAL   (KIND=r4)                :: raux3(SIZE(field3,1),SIZE(field3,2))

    READ (UNIT=n, Rec=irec)   raux1
    READ (UNIT=n, Rec=irec+1) raux2
    READ (UNIT=n, Rec=irec+2) raux3

    field1 = REAL(raux1, r8)
    field2 = REAL(raux2, r8)
    field3 = REAL(raux3, r8)

    !REWIND(UNIT=n)
  END SUBROUTINE ReadNFTGZ4_rec
  SUBROUTINE ReadNFTGZ8_rec(n,irec,field1,field2,field3)
    INTEGER(KIND=i4), INTENT(in   ) :: n 
    INTEGER         , INTENT(in   ) :: irec
    REAL   (KIND=r4), INTENT(out  ) :: field1 (:,:)
    REAL   (KIND=r4), INTENT(out  ) :: field2 (:,:)
    REAL   (KIND=r4), INTENT(out  ) :: field3 (:,:)
    REAL   (KIND=r4)                :: raux1(SIZE(field1,1),SIZE(field1,2))
    REAL   (KIND=r4)                :: raux2(SIZE(field2,1),SIZE(field2,2))
    REAL   (KIND=r4)                :: raux3(SIZE(field3,1),SIZE(field3,2))

    READ (UNIT=n, Rec=irec ) raux1
    READ (UNIT=n, Rec=irec+1) raux2
    READ (UNIT=n, Rec=irec+2) raux3

    field1 = REAL(raux1, r4)
    field2 = REAL(raux2, r4)
    field3 = REAL(raux3, r4)

    REWIND(UNIT=n)
  END SUBROUTINE ReadNFTGZ8_rec


  SUBROUTINE ReadNFTGZ4_txt(n,irec,field1,field2,field3,type2)
    INTEGER, INTENT(in   ) :: n 
    INTEGER, INTENT(in   ) :: irec
    REAL   (KIND=r8), INTENT(out  ) :: field1 (:,:)
    REAL   (KIND=r8), INTENT(out  ) :: field2 (:,:)
    REAL   (KIND=r8), INTENT(out  ) :: field3 (:,:)
    CHARACTER(LEN=3), INTENT(in   ) :: type2
    REAL   (KIND=r4)                :: raux1(SIZE(field1,1),SIZE(field1,2))
    REAL   (KIND=r4)                :: raux2(SIZE(field2,1),SIZE(field2,2))
    REAL   (KIND=r4)                :: raux3(SIZE(field3,1),SIZE(field3,2))
    PRINT*,type2
    READ (n, *) raux1
    READ (n, *) raux2
    READ (n, *) raux3
    REWIND(n)
    field1 = REAL(raux1, r8)
    field2 = REAL(raux2, r8)
    field3 = REAL(raux3, r8)

    REWIND(UNIT=n)
  END SUBROUTINE ReadNFTGZ4_txt
  SUBROUTINE ReadNFTGZ8_txt(n,irec,field1,field2,field3,type2)
    INTEGER(KIND=i4), INTENT(in   ) :: n 
    INTEGER         , INTENT(in   ) :: irec
    REAL   (KIND=r4), INTENT(out  ) :: field1 (:,:)
    REAL   (KIND=r4), INTENT(out  ) :: field2 (:,:)
    REAL   (KIND=r4), INTENT(out  ) :: field3 (:,:)
    CHARACTER(LEN=3), INTENT(in   ) :: type2
    REAL   (KIND=r4)                :: raux1(SIZE(field1,1),SIZE(field1,2))
    REAL   (KIND=r4)                :: raux2(SIZE(field2,1),SIZE(field2,2))
    REAL   (KIND=r4)                :: raux3(SIZE(field3,1),SIZE(field3,2))
    PRINT*,type2
    READ (n, *) raux1
    READ (n, *) raux2
    READ (n, *) raux3
    REWIND(n)
    field1 = REAL(raux1, r4)
    field2 = REAL(raux2, r4)
    field3 = REAL(raux3, r4)

    REWIND(UNIT=n)
  END SUBROUTINE ReadNFTGZ8_txt

  SUBROUTINE WriteDHNHead4(n, ifday, tod)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    INTEGER(KIND=i4), INTENT(IN)  :: ifday
    REAL   (KIND=r4), INTENT(IN)  :: tod
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteProgHead4)**"
    WRITE(n)ifday, tod
  END SUBROUTINE WriteDHNHead4
  SUBROUTINE WriteDHNHead8(n, ifday, tod)
    INTEGER(KIND=i8), INTENT(IN)  :: n
    INTEGER(KIND=i8), INTENT(IN)  :: ifday
    REAL   (KIND=r8), INTENT(IN)  :: tod
    INTEGER(KIND=i4) :: iaux(2)
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteProgHead8)**"
    iaux(  1 ) = INT (ifday , i4)
    iaux(  2 ) = INT (tod   , i4)
    WRITE(n)iaux
  END SUBROUTINE WriteDHNHead8
END MODULE IOLowLevel
