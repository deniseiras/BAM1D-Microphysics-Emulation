MODULE ModAddTend

  USE Constants, ONLY :	  &
       r8,i8,r4,i4

  USE Sizes,  ONLY : &
       ibMax,   & ! intent(in)
       jbMax,        & ! intent(in)
       kMax            ! intent(in)
  USE FieldsDynamics, ONLY : &
       fgyum,        & ! intent(in)
       fgyvm,        & ! intent(in)
       fgtdm,        & ! intent(in)
       fgvdlnpm,     & ! intent(in)
       fgum,         & ! intent(in)
       fgvm,         & ! intent(in)
       fgtmpm,       & ! intent(in)
       fgqm,         & ! intent(in)
       fglnpm,       & ! intent(in)
       fgyu,         & ! intent(inout)
       fgyv,         & ! intent(inout)
       fgtd,         & ! intent(inout)
       fgqd,         & ! intent(inout)
       fgvdlnp,      &         ! intent(inout)
       fgu,          &
       fgv,          &
       fgtmp,        &
       fgq  ,        &
       fglnps ,      &
       fgice ,       &
       fgicem ,      &
       fgicet,       &
       fgliq,        &
       fgliqm,       &
       fgliqt,       &
       fgvar,        &
       fgvarm,       &
       fgvart

  USE Options, ONLY :       &
      nClass  ,microphys           , slhum 

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: AddTend_RK4   !4th order Runge Kutta
CONTAINS

  ! addtend: finish tendency computations, adding contributions from
  !          old and current time step.

  SUBROUTINE AddTend_RK4(dt, nlnminit,grid)
    REAL(KIND=r8)            , INTENT(IN) :: dt
    LOGICAL(KIND=i8)         , INTENT(IN) :: nlnminit    
    CHARACTER(len=*), INTENT(IN) :: grid     
    INTEGER :: ib, jb, k,kk
    REAL(KIND=r8) :: dt2
    IF(grid=='1D  ') THEN
      jbMax=1
      ibMax=1
    END IF
    dt2 = dt + dt
       DO jb = 1, jbMax
          DO k = 1, kMax
             DO ib = 1, ibMax
                fgu  (ib,k,jb)   = fgum(ib,k,jb)   + dt2 * ( fgyu(ib,k,jb)) 
                fgv  (ib,k,jb)   = fgvm(ib,k,jb)   + dt2 * ( fgyv(ib,k,jb)) 
                fgtmp(ib,k,jb)   = fgtmpm(ib,k,jb) + dt2 * ( fgtd(ib,k,jb)) 
                fgq  (ib,k,jb)   = fgqm(ib,k,jb)   + dt2 * ( fgqd(ib,k,jb)) 
             END DO
          END DO
       END DO
       IF (microphys) THEN
          DO jb = 1, jbMax
             DO k = 1, kMax
                DO ib = 1, ibMax
                   fgice  (ib,k,jb)   = fgicem(ib,k,jb)   + dt2 * ( fgicet(ib,k,jb)) 
                   fgliq  (ib,k,jb)   = fgliqm(ib,k,jb)   + dt2 * ( fgliqt(ib,k,jb)) 
                END DO
             END DO
             DO kk=1,nClass
                DO k = 1, kMax
                   DO ib = 1, ibMax
                      fgvar  (ib,k,jb,kk)   = fgvarm(ib,k,jb,kk)   + dt2 * ( fgvart(ib,k,jb,kk)) 
                   END DO
                END DO
             END DO
          END DO
       ENDIF
       
       DO jb = 1, jbMax
          DO ib = 1, ibMax
             fglnps(ib,jb) = fglnpm(ib,jb) !+ dt2 * ( fgvdlnp(ib,jb) + fgvdlnpm(ib,jb) )
          END DO
       END DO
  END SUBROUTINE AddTend_RK4
END MODULE ModAddTend
