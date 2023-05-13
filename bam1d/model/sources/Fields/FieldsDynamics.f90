!
!  $Author: pkubota $
!  $Date: 2011/04/20 14:34:46 $
!  $Revision: 1.10 $
!
MODULE FieldsDynamics

  USE Constants, Only: r8
  USE Options,   Only: &
      slagr,           &
      slhum,           &
      microphys,       &
      nClass,          &
      SL_twotime_scheme

  IMPLICIT NONE


  ! Spectral fields: 13 2D and 7 1D


  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qtmpp
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qtmpt
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qtmpt_si
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qrotp
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qrott
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qrott_si
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qdivp
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qdivt
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qdivt_si
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qqp
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qice
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: qvar
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qliq
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qup
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qvp
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qtphi
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qqphi
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qdiaten
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ) :: qlnpp
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ) :: qlnpl
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ) :: qlnpt
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ) :: qlnpt_si
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ) :: qgzslap
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ) :: qgzs
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ) :: qplam
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ) :: qpphi
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ) :: qgzsphi
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qozon  ! add solange 27-01-2012

  ! Gaussian fields: 28 + 2 * nscalars 3D and 12 2D

  INTEGER       :: adr_scalars
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:,:,:) :: fgpass_scalars

  ! (du/dt) tendency of zonal wind * cos(latitude)
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgyu !T  
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgyum!T-1
  ! -(dv/dt) negative of tendency of v*cos(latitude)
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgyv !T  
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgyvm!T-1
  ! specific humidity tendency (kg/kg/s)
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgqd !T  
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgqdm!T-1 
  ! specific humidity (kg/kg)
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgq  !T
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgqmm!T-1 (depois da convec)
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgqm!T-1 (antes da convec)
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgqp!T+1

  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgicemm !T-1
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgicem !T-1
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgice  !T  
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgicep ! T+1
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgicet ! tendency 

  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:,:) :: fgvarmm !T-1
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:,:) :: fgvarm !T-1
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:,:) :: fgvar  !T  
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:,:) :: fgvarp ! T+1
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:,:) :: fgvart !tendency  

  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgliqmm  !T-1
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgliqm  !T-1
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgliq   !T  
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgliqp  ! T+1
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgliqt  !tendency 
  ! virtual temperature tendency (K/s)
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgtd !T  
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgtdm!T-1
  ! virtual temperature (K)
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgtmp !T
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgtmpm!T-1
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgtmpp!T+1
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgtmpmm
  ! zonal wind (?)
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgu !T  
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgum!T-1
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgumm!T-1
  ! meridional wind (?)
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgv !T  
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgvm!T-1
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgvmm!T-1

  ! vertical wind (?)
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgw !T  
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgwm!T-1 
  ! 
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgdiv !T  
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgdivm!T-1
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgdivmm!T-1

  !
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgtphi !T  
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgtphim!T-1
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgtphimm!T-1
  !
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fgpphi
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fgpphim
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fgpphimm
  !
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgtlam
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgtlamm
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgtlammm
  !
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fgplam
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fgplamm
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fgplammm
  !
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgqlam
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgulam
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgvlam

  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgrot
  ! vertical omega(?) velocity (cb/s)
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: omg!T

  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgqphi
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fguphi
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgvphi
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fgumean
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fgvmean
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fgvdlnp
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fglnps
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fglnpm
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fglnpmm
  ! surface pressure (cb) 
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fgps !T
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fgpsp!T+1
  ! geopotential at the surface, fgzs/grav == topography
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fgzs
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fgzslam
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fgzsphi
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fgvdlnpm

  ! variables for Nudging, though dimensions differ from those inside Nudging_Ascii
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:, :  ,:) :: fguo
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:, :  ,:) :: fgvo
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:, :  ,:) :: fgto
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:, :  ,:) :: fgqo

  ! variables for specified surface (specSfc)
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION (:,:) :: fgTsfc
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION (:,:) :: fgH
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION (:,:) :: fgLE
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION (:,:) :: fgTau

  ! variables for Geostrophy, though dimensions differ from those inside Geostrophy_Ascii
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:, :  ,:) :: fgug
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:, :  ,:) :: fgvg

CONTAINS
  SUBROUTINE InitFields(ibMax, kMax, jbMax, mnMax, mnExtMax, &
                        kmaxloc, MNMax_si, jbMax_ext, nscalars)
    INTEGER, INTENT(IN) :: ibMax
    INTEGER, INTENT(IN) :: kMax
    INTEGER, INTENT(IN) :: kMaxloc
    INTEGER, INTENT(IN) :: jbMax
    INTEGER, INTENT(IN) :: jbMax_ext
    INTEGER, INTENT(IN) :: mnMax
    INTEGER, INTENT(IN) :: MNMax_si
    INTEGER, INTENT(IN) :: mnExtMax
    INTEGER, INTENT(IN) :: nscalars

    ALLOCATE (qtmpp(2*mnMax,kMaxloc))
    qtmpp=0.0_r8
    ALLOCATE (qtmpt(2*mnMax,kMaxloc))
    qtmpt=0.0_r8
    ALLOCATE (qtmpt_si(2*mnMax_si,kMax))
    qtmpt_si=0.0_r8
    ALLOCATE (qrotp(2*mnMax,kMaxloc))
    qrotp=0.0_r8
    ALLOCATE (qrott(2*mnMax,kMaxloc))
    qrott=0.0_r8
    ALLOCATE (qrott_si(2*mnMax_si,kMax))
    qrott_si=0.0_r8
    ALLOCATE (qdivp(2*mnMax,kMaxloc))
    qdivp=0.0_r8
    ALLOCATE (qdivt(2*mnMax,kMaxloc))
    qdivt=0.0_r8
    ALLOCATE (qdivt_si(2*mnMax_si,kMax))
    qdivt_si=0.0_r8
    ALLOCATE (qqp(2*mnMax,kMaxloc))
    qqp=0.0_r8
    IF (microphys) THEN
       IF(nClass>0)THEN
          ALLOCATE (qvar(2*mnMax,kMaxloc,nClass))
          qvar=0.0_r8
       END IF
       ALLOCATE (qice(2*mnMax,kMaxloc))
       qice=0.0_r8
       ALLOCATE (qliq(2*mnMax,kMaxloc))
       qliq=0.0_r8
    ENDIF
    ALLOCATE (qdiaten(2*mnMax,kMaxloc))
    qdiaten=0.0_r8
    ALLOCATE (qup(2*mnExtMax,kMaxloc))
    qup=0.0_r8
    ALLOCATE (qvp(2*mnExtMax,kMaxloc))
    qvp=0.0_r8
    ALLOCATE (qtphi(2*mnExtMax,kMaxloc))
    qtphi=0.0_r8
    ALLOCATE (qqphi(2*mnExtMax,kMaxloc))
    qqphi=0.0_r8
    ALLOCATE (qlnpp(2*mnMax))
    qlnpp=0.0_r8
    ALLOCATE (qlnpl(2*mnMax))
    qlnpl=0.0_r8
    ALLOCATE (qlnpt(2*mnMax))
    qlnpt=0.0_r8
    ALLOCATE (qlnpt_si(2*mnMax_si))
    qlnpt_si=0.0_r8
    ALLOCATE (qgzslap(2*mnMax))
    qgzslap=0.0_r8
    ALLOCATE (qgzs(2*mnMax))
    qgzs=0.0_r8
    ALLOCATE (qplam(2*mnMax))
    qplam=0.0_r8
    ALLOCATE (qpphi(2*mnExtMax))
    qpphi=0.0_r8
    ALLOCATE (qgzsphi(2*mnExtMax))
    qgzsphi=0.0_r8
    ALLOCATE (qozon(2*mnMax,kMaxloc))  ! solange add
    qozon=0.0_r8

    IF (nscalars .gt. 0) THEN
       ALLOCATE (fgpass_scalars(ibMax,kMax,jbmax_ext,nscalars,2))
       fgpass_scalars=0.0_r8
      ELSE
       ALLOCATE (fgpass_scalars(1,1,1,1,1))
       fgpass_scalars=0.0_r8
    ENDIF
    adr_scalars = 1
    
    ALLOCATE (fguo(ibMax,kMax,jbMax))
    fguo=0.0_r8
    ALLOCATE (fgvo(ibMax,kMax,jbMax))
    fgvo=0.0_r8
    ALLOCATE (fgto(ibMax,kMax,jbMax))
    fgto=0.0_r8
    ALLOCATE (fgqo(ibMax,kMax,jbMax))
    fgqo=0.0_r8

    ALLOCATE (fgTsfc(ibMax,jbMax))
    fgTsfc=0.0_r8
    ALLOCATE (fgH(ibMax,jbMax))
    fgH=0.0_r8
    ALLOCATE (fgLE(ibMax,jbMax))
    fgLE=0.0_r8
    ALLOCATE (fgTau(ibMax,jbMax))
    fgTau=0.0_r8

    ALLOCATE (fgug(ibMax,kMax,jbMax))
    fgug=0.0_r8
    ALLOCATE (fgvg(ibMax,kMax,jbMax))
    fgvg=0.0_r8

    ALLOCATE (fgyu(ibMax,kMax,jbMax))
    fgyu=0.0_r8
    ALLOCATE (fgyv(ibMax,kMax,jbMax))
    fgyv=0.0_r8
    ALLOCATE (fgtd(ibMax,kMax,jbMax))
    fgtd=0.0_r8
    ALLOCATE (fgqd(ibMax,kMax,jbMax_ext))
    fgqd=0.0_r8
    ALLOCATE (fgdiv(ibMax,kMax,jbMax))
    fgdiv=0.0_r8
    ALLOCATE (fgrot(ibMax,kMax,jbMax))
    fgrot=0.0_r8
    ALLOCATE (fgq(ibMax,kMax,jbMax))
    fgq=0.0_r8
    ALLOCATE (fgtmp(ibMax,kMax,jbMax))
    fgtmp=0.0_r8
    IF (SL_twotime_scheme) THEN
       ALLOCATE (fgu(ibMax,kMax,jbMax))
       fgu=0.0_r8
       ALLOCATE (fgv(ibMax,kMax,jbMax))
       fgv=0.0_r8    
       ALLOCATE (fgw(ibMax,kMax,jbMax))
       fgw=0.0_r8
       ALLOCATE (fgum(ibMax,kMax,jbMax_ext))
       fgum=0.0_r8
       ALLOCATE (fgumm(ibMax,kMax,jbMax_ext))
       fgumm=0.0_r8
       ALLOCATE (fgvm(ibMax,kMax,jbMax_ext))
       fgvm=0.0_r8
       ALLOCATE (fgvmm(ibMax,kMax,jbMax_ext))
       fgvmm=0.0_r8
       ALLOCATE (fgwm(ibMax,kMax,jbMax_ext))
       fgwm=0.0_r8
       ALLOCATE (fgtmpp(ibMax,kMax,jbMax_ext))
       fgtmpp=0.0_r8
       ALLOCATE (fgqp(ibMax,kMax,jbMax_ext))
       fgqp=0.0_r8
       ALLOCATE (fgyum(ibMax,kMax,jbMax))
       fgyum=0.0_r8
       ALLOCATE (fgyvm(ibMax,kMax,jbMax))
       fgyvm=0.0_r8
       ALLOCATE (fgtdm(ibMax,kMax,jbMax))
       fgtdm=0.0_r8
       ALLOCATE (fgqdm(ibMax,kMax,jbMax))
       fgqdm=0.0_r8
      ELSE
       ALLOCATE (fgu(ibMax,kMax,jbMax_ext))
       fgu=0.0_r8
       ALLOCATE (fgv(ibMax,kMax,jbMax_ext))
       fgv=0.0_r8    
       ALLOCATE (fgw(ibMax,kMax,jbMax_ext))
       fgw=0.0_r8
       ALLOCATE (fgum(ibMax,kMax,jbMax))
       fgum=0.0_r8
       ALLOCATE (fgumm(ibMax,kMax,jbMax))
       fgumm=0.0_r8
       ALLOCATE (fgvm(ibMax,kMax,jbMax))
       fgvm=0.0_r8    
       ALLOCATE (fgvmm(ibMax,kMax,jbMax))
       fgvmm=0.0_r8    
       ALLOCATE (fgtmpp(ibMax,kMax,jbMax))
       fgtmpp=0.0_r8
       ALLOCATE (fgqp(ibMax,kMax,jbMax))
       fgqp=0.0_r8
       ALLOCATE (fgyum(ibMax,kMax,jbMax_ext))
       fgyum=0.0_r8
       ALLOCATE (fgyvm(ibMax,kMax,jbMax_ext))
       fgyvm=0.0_r8
       ALLOCATE (fgtdm(ibMax,kMax,jbMax_ext))
       fgtdm=0.0_r8
       ALLOCATE (fgqdm(ibMax,kMax,jbMax))
       fgqdm=0.0_r8
    ENDIF 
    IF (microphys) THEN 
       ALLOCATE (fgice(ibMax,kMax,jbMax))
       fgice=0.0_r8
       ALLOCATE (fgicep(ibMax,kMax,jbMax))
       fgicep=0.0_r8
       ALLOCATE (fgicemm(ibMax,kMax,jbMax))
       fgicemm=0.0_r8
       ALLOCATE (fgicem(ibMax,kMax,jbMax))
       fgicem=0.0_r8
       ALLOCATE (fgicet(ibMax,kMax,jbMax_ext))
       fgicet=0.0_r8
       IF(nClass>0)THEN
          ALLOCATE (fgvar(ibMax,kMax,jbMax,nClass))
          fgvar=0.0_r8
          ALLOCATE (fgvarp(ibMax,kMax,jbMax,nClass))
          fgvarp=0.0_r8
          ALLOCATE (fgvarm(ibMax,kMax,jbMax,nClass))
          fgvarm=0.0_r8
          ALLOCATE (fgvarmm(ibMax,kMax,jbMax,nClass))
          fgvarmm=0.0_r8

          ALLOCATE (fgvart(ibMax,kMax,jbMax_ext,nClass))
          fgvart=0.0_r8
       END IF 
       ALLOCATE (fgliq(ibMax,kMax,jbMax))
       fgliq=0.0_r8
       ALLOCATE (fgliqp(ibMax,kMax,jbMax))
       fgliqp=0.0_r8
       ALLOCATE (fgliqmm(ibMax,kMax,jbMax))
       fgliqmm=0.0_r8
       ALLOCATE (fgliqm(ibMax,kMax,jbMax))
       fgliqm=0.0_r8
       ALLOCATE (fgliqt(ibMax,kMax,jbMax_ext))
       fgliqt=0.0_r8
    ENDIF
    ALLOCATE (fgqmm(ibMax,kMax,jbMax))
    fgqmm=0.0_r8
    ALLOCATE (fgqm(ibMax,kMax,jbMax))
    fgqm=0.0_r8 
    ALLOCATE (omg(ibMax,kMax,jbMax))
    omg=0.0_r8
    ALLOCATE (fgtmpm(ibMax,kMax,jbMax))
    fgtmpm=0.0_r8
    ALLOCATE (fgtmpmm(ibMax,kMax,jbMax))
    fgtmpmm=0.0_r8    
    ALLOCATE (fgdivm(ibMax,kMax,jbMax))
    fgdivm=0.0_r8
    ALLOCATE (fgdivmm(ibMax,kMax,jbMax))
    fgdivmm=0.0_r8
    ALLOCATE (fgtphi(ibMax,kMax,jbMax))
    fgtphi=0.0_r8
    ALLOCATE (fgqphi(ibMax,kMax,jbMax))
    fgqphi=0.0_r8
    ALLOCATE (fguphi(ibMax,kMax,jbMax))
    fguphi=0.0_r8
    ALLOCATE (fgvphi(ibMax,kMax,jbMax))
    fgvphi=0.0_r8
    ALLOCATE (fgtphim(ibMax,kMax,jbMax))
    fgtphim=0.0_r8
    ALLOCATE (fgtphimm(ibMax,kMax,jbMax))
    fgtphimm=0.0_r8
    ALLOCATE (fgtlam(ibMax,kMax,jbMax))
    fgtlam=0.0_r8
    ALLOCATE (fgqlam(ibMax,kMax,jbMax))
    fgqlam=0.0_r8
    ALLOCATE (fgulam(ibMax,kMax,jbMax))
    fgulam=0.0_r8
    ALLOCATE (fgvlam(ibMax,kMax,jbMax))
    fgvlam=0.0_r8
    ALLOCATE (fgtlamm(ibMax,kMax,jbMax))
    fgtlamm=0.0_r8
    ALLOCATE (fgtlammm(ibMax,kMax,jbMax))
    fgtlammm=0.0_r8
    IF (SL_twotime_scheme) THEN
       ALLOCATE (fgumean(ibMax,     jbMax))
       fgumean=0.0_r8
       ALLOCATE (fgvmean(ibMax,     jbMax))
       fgvmean=0.0_r8
       ALLOCATE (fgpsp(ibMax,     jbMax_ext))
       fgpsp=0.0_r8
       ALLOCATE (fgvdlnpm(ibMax,     jbMax))
       fgvdlnpm=0.0_r8
      ELSE
       ALLOCATE (fgumean(ibMax,     jbMax_ext))
       fgumean=0.0_r8
       ALLOCATE (fgvmean(ibMax,     jbMax_ext))
       fgvmean=0.0_r8
       ALLOCATE (fgpsp(ibMax,     jbMax))
       fgpsp=0.0_r8
       ALLOCATE (fgvdlnpm(ibMax,     jbMax_ext))
       fgvdlnpm=0.0_r8
    ENDIF 
    ALLOCATE (fgvdlnp(ibMax,     jbMax))
    fgvdlnp=0.0_r8
    ALLOCATE (fglnps(ibMax,     jbMax))
    fglnps=0.0_r8
    ALLOCATE (fgps(ibMax,     jbMax))
    fgps=0.0_r8
    ALLOCATE (fglnpm(ibMax,     jbMax))
    fglnpm=0.0_r8
    ALLOCATE (fglnpmm(ibMax,     jbMax))
    fglnpmm=0.0_r8
    ALLOCATE (fgpphi(ibMax,     jbMax))
    fgpphi=0.0_r8
    ALLOCATE (fgplam(ibMax,     jbMax))
    fgplam=0.0_r8
    ALLOCATE (fgplamm(ibMax,     jbMax))
    fgplamm=0.0_r8
    ALLOCATE (fgplammm(ibMax,     jbMax))
    fgplammm=0.0_r8
    ALLOCATE (fgpphim(ibMax,     jbMax))
    fgpphim=0.0_r8
    ALLOCATE (fgpphimm(ibMax,     jbMax))
    fgpphimm=0.0_r8
    ALLOCATE (fgzs(ibMax,     jbMax))
    fgzs=0.0_r8
    ALLOCATE (fgzslam(ibMax,     jbMax))
    fgzslam=0.0_r8
    ALLOCATE (fgzsphi(ibMax,     jbMax))
    fgzsphi=0.0_r8
  END SUBROUTINE InitFields
END MODULE FieldsDynamics
