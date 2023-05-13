MODULE ModTimeFilter 
  
  USE Constants, ONLY: i8,r8,i4,r4

  USE Sizes, ONLY : &
       ibMax,       & ! intent(in)
       jbMax,       & ! intent(in)
       kMax           ! intent(in)
  USE FieldsDynamics, ONLY : &
       fgtmpm ,  & ! intent(inout)
       fgtmpmm,  &
       fgtmp  ,  & ! intent(in)
       fgdivm ,  & ! intent(inout)
       fgdivmm,  &
       fgdiv  ,  & ! intent(in)
       fgum   ,  & ! intent(inout)
       fgumm  ,  & ! intent(inout)     
       fgu    ,  & ! intent(in)
       fgvm   ,  & ! intent(inout)
       fgvmm  ,  & ! intent(inout)       
       fgv    ,  & ! intent(in)
       fgqm   ,  & ! intent(inout)
       fgqmm  ,  & ! intent(out)
       fgq    ,  & ! intent(in)
       fgtlamm,  & ! intent(inout)
       fgtlam ,  & ! intent(in)
       fgtphim,  & ! intent(inout)
       fgtphimm, &
       fgtphi ,  & ! intent(in)
       fglnpm ,  & ! intent(inout)
       fglnpmm,  &
       fglnps ,  & ! intent(in)
       fgplamm,  & ! intent(inout)
       fgtlammm, &
       fgplam ,  & ! intent(in)
       fgplammm, &
       fgpphim,  & ! intent(inout)
       fgpphimm,  & ! intent(inout)
       fgpphi  ,  &
       fgicemm , &
       fgicem , &
       fgice , &
       fgicep , &
       fgliqmm , &
       fgliqm , &
       fgliq , &
       fgliqp , &
       fgvarmm , &
       fgvarm , &
       fgvar , &
       fgvarp,&
       fgqp
  USE Options, ONLY :       &
      nClass  ,microphys           , slhum 

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: TimeFilterStep1
  PUBLIC :: TimeFilterStep2
CONTAINS

  ! TimeFilterStep1: First part of the asselin/robert time-filter
  ! (computes a partially filtered value of fold)

  SUBROUTINE TimeFilterStep1(fa, fb,grid)
    REAL(KIND=r8)            , INTENT(IN) :: fa
    REAL(KIND=r8)            , INTENT(IN) :: fb
    CHARACTER(len=*), INTENT(IN) :: grid  
    INTEGER :: ib, jb, k      ,kk
    IF(grid == '1D  ') THEN
      jbMax=1
      ibMax=1
    END IF
    fgqmm   =fgqm
    fgumm   =fgum
    fgvmm   =fgvm
    fgtmpmm =fgtmpm
    fgdivmm =fgdivm
    fglnpmm =fglnpm
    fgtlammm=fgtlamm
    fgtphimm=fgtphim
    fgplammm=fgplamm 
    fgpphimm=fgpphim
    IF(microphys)THEN
      fgicemm=fgicem
      fgliqmm=fgliqm
      IF(nClass > 0)THEN
         fgvarmm=fgvarm
      END IF
    END IF
    DO jb = 1, jbMax
       DO k = 1, kMax
          DO ib = 1, ibMax
             fgtmpm (ib,k,jb) = fa*fgtmp (ib,k,jb) + fb*fgtmpm (ib,k,jb)
             fgdivm (ib,k,jb) = fa*fgdiv (ib,k,jb) + fb*fgdivm (ib,k,jb)
             fgum   (ib,k,jb) = fa*fgu   (ib,k,jb) + fb*fgum   (ib,k,jb)
             fgvm   (ib,k,jb) = fa*fgv   (ib,k,jb) + fb*fgvm   (ib,k,jb)
!             fgqm   (ib,k,jb) = fa*fgq   (ib,k,jb) + fb*fgqm   (ib,k,jb)
             fgtlamm(ib,k,jb) = fa*fgtlam(ib,k,jb) + fb*fgtlamm(ib,k,jb)
             fgtphim(ib,k,jb) = fa*fgtphi(ib,k,jb) + fb*fgtphim(ib,k,jb)
          END DO
       END DO
    END DO
    IF(.not.slhum) THEN
       DO jb = 1, jbMax
          DO k = 1, kMax
             DO ib = 1, ibMax
                fgqm(ib,k,jb) = fa*fgq(ib,k,jb) + fb*fgqm(ib,k,jb)
             END DO
          END DO
      END DO
    ELSE 
       IF (fa.ne.0.0_r8) THEN 
          DO jb = 1, jbMax
             DO k = 1, kMax
                DO ib = 1, ibMax
                   fgqm(ib,k,jb) = fa*fgq(ib,k,jb) + &
                                   fb*(fgqm(ib,k,jb)+fgqp(ib,k,jb))
                END DO
             END DO
             IF (microphys) THEN
                DO k = 1, kMax
                   DO ib = 1, ibMax
                      fgicem(ib,k,jb) = fa*fgice(ib,k,jb) + &
                                        fb*(fgicem(ib,k,jb)+fgicep(ib,k,jb))
                      fgliqm(ib,k,jb) = fa*fgliq(ib,k,jb) + &
                                        fb*(fgliqm(ib,k,jb)+fgliqp(ib,k,jb))
                   END DO
                END DO
                DO kk=1,nClass
                   DO k = 1, kMax
                      DO ib = 1, ibMax
                         fgvarm(ib,k,jb,kk) = fa*fgvar(ib,k,jb,kk) + &
                                       fb*(fgvarm(ib,k,jb,kk)+fgvarp(ib,k,jb,kk))
                      END DO
                   END DO
                END DO
             ENDIF
          END DO
       ENDIF
       DO jb = 1, jbMax
          DO k = 1, kMax
             DO ib = 1, ibMax
                fgq(ib,k,jb) = fgqp(ib,k,jb)
             END DO
          END DO
       END DO
       IF (microphys) THEN
          DO jb = 1, jbMax
             DO k = 1, kMax
                DO ib = 1, ibMax
                   fgice(ib,k,jb) = fgicep(ib,k,jb)
                   fgliq(ib,k,jb) = fgliqp(ib,k,jb)
                END DO
             END DO
             DO kk=1,nClass
                DO k = 1, kMax
                   DO ib = 1, ibMax
                      fgvar(ib,k,jb,kk) = fgvarp(ib,k,jb,kk)
                   END DO
                END DO
             END DO
          END DO
       ENDIF
    ENDIF


    DO jb = 1, jbMax
       DO ib = 1, ibMax
          fglnpm (ib,jb) = fa*fglnps(ib,jb) + fb*fglnpm (ib,jb)
          fgplamm(ib,jb) = fa*fgplam(ib,jb) + fb*fgplamm(ib,jb)
          fgpphim(ib,jb) = fa*fgpphi(ib,jb) + fb*fgpphim(ib,jb)
       END DO
    END DO
  END SUBROUTINE TimeFilterStep1


  ! TimeFilterStep2: Second part of the asselin/robert time-filter
  ! (the partially filtered value of fold is filtered completely)

  SUBROUTINE TimeFilterStep2(fb1,grid)
    REAL(KIND=r8)            , INTENT(IN) :: fb1
    CHARACTER(len=*), INTENT(IN) :: grid  
    INTEGER :: ib, jb, k,kk
    IF(grid == '1D  ') THEN
      jbMax=1
      ibMax=1
    END IF
    DO jb = 1, jbMax
       DO k = 1, kMax
          DO ib = 1, ibMax
             fgtmpm (ib,k,jb) = fgtmpm (ib,k,jb) + fb1*fgtmp (ib,k,jb)
             fgdivm (ib,k,jb) = fgdivm (ib,k,jb) + fb1*fgdiv (ib,k,jb)
             fgum   (ib,k,jb) = fgum   (ib,k,jb) + fb1*fgu   (ib,k,jb)
             fgvm   (ib,k,jb) = fgvm   (ib,k,jb) + fb1*fgv   (ib,k,jb)
             fgqm   (ib,k,jb) = fgqm   (ib,k,jb) + fb1*fgq   (ib,k,jb)
             fgtlamm(ib,k,jb) = fgtlamm(ib,k,jb) + fb1*fgtlam(ib,k,jb)
             fgtphim(ib,k,jb) = fgtphim(ib,k,jb) + fb1*fgtphi(ib,k,jb)
          END DO
       END DO
    END DO
    IF(.not.slhum) THEN
       IF (microphys) THEN
          DO jb = 1, jbMax
             DO k = 1, kMax
                DO ib = 1, ibMax
                   fgicem(ib,k,jb) = fgicem(ib,k,jb) + fb1*(fgice(ib,k,jb))
                   fgliqm(ib,k,jb) = fgliqm(ib,k,jb) + fb1*(fgliq(ib,k,jb))
                END DO
             END DO
             DO kk=1,nClass
                DO k = 1, kMax
                   DO ib = 1, ibMax
                      fgvarm  (ib,k,jb,kk)   = fgvarm(ib,k,jb,kk)   + fb1 * ( fgvar(ib,k,jb,kk)) 
                   END DO
                END DO
             END DO
          END DO
       ENDIF
    ENDIF
    DO jb = 1, jbMax
       DO ib = 1, ibMax
          fglnpm (ib,jb) = fglnpm (ib,jb) + fb1*fglnps(ib,jb)
          fgplamm(ib,jb) = fgplamm(ib,jb) + fb1*fgplam(ib,jb)
          fgpphim(ib,jb) = fgpphim(ib,jb) + fb1*fgpphi(ib,jb)
       END DO
    END DO
  END SUBROUTINE TimeFilterStep2
END MODULE ModTimeFilter
