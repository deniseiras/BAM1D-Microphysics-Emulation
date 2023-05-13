MODULE Parallelism

  ! Deals with parallelism; exports
  !   MPI procedures for
  !       initialization 
  !       normal finalization 
  !       error finalization
  !   MPI process data:
  !       how many processes in MPI computation
  !       which processes is this
  ! Creates one file per process for printouts,
  ! unscrambling outputs. Provides procedures for:
  !       writing messages to output and to the file
  !       writing messages only to the file


  IMPLICIT NONE

  PRIVATE

  ! public data

  INTEGER,          PUBLIC :: maxNodes      ! # MPI processes in the computation
  INTEGER,          PUBLIC :: myId          ! MPI process rank
  INTEGER,          PUBLIC :: maxNodes_four ! # MPI processes in fourier group  
  INTEGER,          PUBLIC :: myId_four     ! MPI process rank in fourier group
  INTEGER,          PUBLIC :: mygroup_four  ! fourier group
  INTEGER,          PUBLIC :: unitDump      ! this process dumping file unit
  INTEGER,          PUBLIC :: COMM_FOUR     ! Communicator of Fourier Group

  ! public procedures

  PUBLIC :: MsgOne
  PUBLIC :: FatalError
  ! private data

  CHARACTER(LEN=4) :: cNProc      ! maxNodes in characters, for printing
  CHARACTER(LEN=4) :: cThisProc   ! myId in characters, for printing
  INTEGER, PARAMETER :: stdout=6


CONTAINS


  !*** Dump message at stdout and dump file ***



  SUBROUTINE MsgOne(h, message)
    CHARACTER(LEN=*), INTENT(IN) :: h, message
       WRITE(stdout,*)  TRIM(h)//' '//TRIM(message)
  END SUBROUTINE MsgOne



  !*** Dump error message everywhere and destroy parallelism ***



  SUBROUTINE FatalError(message)
    CHARACTER(LEN=*), INTENT(IN) :: message
    INTEGER :: ierror=-1
    INTEGER :: ierr
    CHARACTER(LEN=11) :: h="**(ERROR)**"
    WRITE(stdout,*)  TRIM(h)//' '//TRIM(message)
    STOP
  END SUBROUTINE FatalError

END MODULE Parallelism
