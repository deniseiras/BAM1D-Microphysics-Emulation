MODULE Micro_HugMorr_NN

  USE csv_file

  IMPLICIT NONE

  private

  INTEGER , PARAMETER :: r8  = SELECTED_REAL_KIND(P=13,R=300)

  ! TODO ***********************************************************************************
  ! Usar spinup de execução de testes > do que de treinamento ?  (ex 7 para treino e 8 para execução)


  ! 259200  = 3 dias
  ! 1209600 = 14 dias  - usado para ensemble 2014-2015
  ! 4233600 = 49 dias - 3 dias a mais que o dia de inicio da avaliacao do script grads da chuva
  
  ! 20995200 = 243 dias - 01/09/2014 - usado para as rodadas do periodo seco
  ! 5097600  = 59 dias  - 01/03/2014 - usado para as rodadas do periodo umido
  INTEGER , PARAMETER :: dt_spinup = 1209600  ! Spin up time in seconds for microphysics HughMorrison NN
  INTEGER :: jdt_spinup ! Spin up count for microphysics HughMorrison


  !TODO ********************* pass as argument 
  INTEGER, PARAMETER :: nColsFixed = 5  ! pq = 5 ??? não lembro
  INTEGER, PARAMETER :: kMaxFixed = 64  ! TODO - Allocate cm kmax, pra não precisar compilar

  ! > 1 for LSTM, MLP_T, CNN_T, MLP_LSTM_NOLS
  INTEGER, PARAMETER :: window = 5
  
  INTEGER :: times_passed = 0

  integer :: estimated_out_NN_file_size
  real :: estimated_out_NN_file_size_real
  
  

  REAL(KIND=r8) :: sl_t(window, kMaxFixed)
  REAL(KIND=r8) :: Tc_t(window, 1:nColsFixed, 1:kMaxFixed)
  REAL(KIND=r8) :: qv_t(window, 1:nColsFixed, 1:kMaxFixed)
  REAL(KIND=r8) :: qc_t(window, 1:nColsFixed, 1:kMaxFixed)
  REAL(KIND=r8) :: qr_t(window, 1:nColsFixed, 1:kMaxFixed)
  REAL(KIND=r8) :: qi_t(window, 1:nColsFixed, 1:kMaxFixed)
  REAL(KIND=r8) :: qs_t(window, 1:nColsFixed, 1:kMaxFixed)
  REAL(KIND=r8) :: qg_t(window, 1:nColsFixed, 1:kMaxFixed)
  REAL(KIND=r8) :: ni_t(window, 1:nColsFixed, 1:kMaxFixed)
  REAL(KIND=r8) :: ns_t(window, 1:nColsFixed, 1:kMaxFixed)
  REAL(KIND=r8) :: nr_t(window, 1:nColsFixed, 1:kMaxFixed)
  REAL(KIND=r8) :: NG_t(window, 1:nColsFixed, 1:kMaxFixed)
  REAL(KIND=r8) :: NC_t(window, 1:nColsFixed, 1:kMaxFixed)
  REAL(KIND=r8) :: kzh_tke_t(window, 1:nColsFixed, 1:kMaxFixed)
  REAL(KIND=r8) :: omega_t(window, 1:nColsFixed, 1:kMaxFixed)
  REAL(KIND=r8) :: LSRAIN_t(window, 1:nColsFixed)
  REAL(KIND=r8) :: LSSNOW_t(window, 1:nColsFixed)


  public Init_Micro_HugMorr_NN, hugMorr_spinup_reached, RunMicro_HugMorr_NN, fill_vars_timesteps
  

CONTAINS


  SUBROUTINE Init_Micro_HugMorr_NN(dt)
    IMPLICIT NONE

    REAL(KIND=r8), INTENT(IN) :: dt  ! timestep  = delta time of model integration. Eg. 60, 360 seconds
    CHARACTER(LEN=*), PARAMETER     :: h='**(Init_Micro_HugMorr)**'

    jdt_spinup = dt_spinup / dt
    ! estimated_out_NN_file_size = kMax * window * 100  ! 
    estimated_out_NN_file_size = 0  
    estimated_out_NN_file_size_real = 0
    

  END SUBROUTINE Init_Micro_HugMorr_NN


   SUBROUTINE hugMorr_spinup_reached(jdt, reached, passed)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: jdt
    LOGICAL, INTENT(out) :: reached
    LOGICAL, INTENT(out) :: passed

    reached = jdt .gt. jdt_spinup

    if (reached) then
      times_passed = times_passed + 1
    endif
    passed = times_passed >= window

  END SUBROUTINE hugMorr_spinup_reached


  subroutine fill_vars_timesteps(nCols, kMax, sl, Tc, qv, qc, qr, qi, qs, qg, ni, ns, nr, NG, NC, kzh_tke, omega, LSRAIN, LSSNOW)
    implicit none

    INTEGER      , INTENT(IN) :: nCols
    INTEGER      , INTENT(IN) :: kMax
    REAL(KIND=r8), INTENT(IN) :: sl(kMax)
    REAL(KIND=r8), INTENT(IN) :: Tc(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(IN) :: qv(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(IN) :: qc(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(IN) :: qr(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(IN) :: qi(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(IN) :: qs(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(IN) :: qg(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(IN) :: ni(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(IN) :: ns(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(IN) :: nr(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(IN) :: NG(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(IN) :: NC(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(IN) :: kzh_tke(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(IN) :: omega(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(IN) :: LSRAIN(1:nCols)
    REAL(KIND=r8), INTENT(IN) :: LSSNOW(1:nCols)

    if ( times_passed > 0 .and. times_passed <= window ) then
      print*, 'times passed ===========', times_passed, sl(1:kMax)
      sl_t(times_passed, 1:kMax)             = sl(1:kMax)
      Tc_t(times_passed, 1:nCols, 1:kMax)    = Tc(1:nCols, 1:kMax)
      qv_t(times_passed, 1:nCols, 1:kMax)    = qv(1:nCols, 1:kMax)
      qc_t(times_passed, 1:nCols, 1:kMax)    = qc(1:nCols, 1:kMax)
      qr_t(times_passed, 1:nCols, 1:kMax)    = qr(1:nCols, 1:kMax)
      qi_t(times_passed, 1:nCols, 1:kMax)    = qi(1:nCols, 1:kMax)
      qs_t(times_passed, 1:nCols, 1:kMax)    = qs(1:nCols, 1:kMax)
      qg_t(times_passed, 1:nCols, 1:kMax)    = qg(1:nCols, 1:kMax)
      ni_t(times_passed, 1:nCols, 1:kMax)    = ni(1:nCols, 1:kMax)
      ns_t(times_passed, 1:nCols, 1:kMax)    = ns(1:nCols, 1:kMax)
      nr_t(times_passed, 1:nCols, 1:kMax)    = nr(1:nCols, 1:kMax)
      NG_t(times_passed, 1:nCols, 1:kMax)    = NG(1:nCols, 1:kMax)
      NC_t(times_passed, 1:nCols, 1:kMax)    = NC(1:nCols, 1:kMax)
      kzh_tke_t(times_passed, 1:nCols, 1:kMax)   = kzh_tke(1:nCols, 1:kMax)
      omega_t(times_passed, 1:nCols, 1:kMax) = omega(1:nCols, 1:kMax)
      LSRAIN_t(times_passed, 1:nCols)          = LSRAIN(1:nCols)
      LSSNOW_t(times_passed, 1:nCols)          = LSSNOW(1:nCols)
    else 
      ! only for windows > 1
      if (window > 1 .and. times_passed > window) then
        print*, 'times passed windows===========', times_passed, sl(1:kMax)
        sl_t(1:window - 1, 1:kMax)             = sl_t(2:window, 1:kMax)
        Tc_t(1:window - 1, 1:nCols, 1:kMax)    = Tc_t(2:window, 1:nCols, 1:kMax)
        qv_t(1:window - 1, 1:nCols, 1:kMax)    = qv_t(2:window, 1:nCols, 1:kMax)
        qc_t(1:window - 1, 1:nCols, 1:kMax)    = qc_t(2:window, 1:nCols, 1:kMax)
        qr_t(1:window - 1, 1:nCols, 1:kMax)    = qr_t(2:window, 1:nCols, 1:kMax)
        qi_t(1:window - 1, 1:nCols, 1:kMax)    = qi_t(2:window, 1:nCols, 1:kMax)
        qs_t(1:window - 1, 1:nCols, 1:kMax)    = qs_t(2:window, 1:nCols, 1:kMax)
        qg_t(1:window - 1, 1:nCols, 1:kMax)    = qg_t(2:window, 1:nCols, 1:kMax)
        ni_t(1:window - 1, 1:nCols, 1:kMax)    = ni_t(2:window, 1:nCols, 1:kMax)
        ns_t(1:window - 1, 1:nCols, 1:kMax)    = ns_t(2:window, 1:nCols, 1:kMax)
        nr_t(1:window - 1, 1:nCols, 1:kMax)    = nr_t(2:window, 1:nCols, 1:kMax)
        NG_t(1:window - 1, 1:nCols, 1:kMax)    = NG_t(2:window, 1:nCols, 1:kMax)
        NC_t(1:window - 1, 1:nCols, 1:kMax)    = NC_t(2:window, 1:nCols, 1:kMax)
        kzh_tke_t(1:window - 1, 1:nCols, 1:kMax)   = kzh_tke_t(2:window, 1:nCols, 1:kMax)
        omega_t(1:window - 1, 1:nCols, 1:kMax) = omega_t(2:window, 1:nCols, 1:kMax)
        LSRAIN_t(1:window - 1, 1:nCols)        = LSRAIN_t(2:window, 1:nCols)
        LSSNOW_t(1:window - 1, 1:nCols)        = LSSNOW_t(2:window, 1:nCols)

        print*, 'times passed windows===========', times_passed, sl(1:kMax)
        sl_t(window, 1:kMax)             = sl(1:kMax)
        Tc_t(window, 1:nCols, 1:kMax)    = Tc(1:nCols, 1:kMax)
        qv_t(window, 1:nCols, 1:kMax)    = qv(1:nCols, 1:kMax)
        qc_t(window, 1:nCols, 1:kMax)    = qc(1:nCols, 1:kMax)
        qr_t(window, 1:nCols, 1:kMax)    = qr(1:nCols, 1:kMax)
        qi_t(window, 1:nCols, 1:kMax)    = qi(1:nCols, 1:kMax)
        qs_t(window, 1:nCols, 1:kMax)    = qs(1:nCols, 1:kMax)
        qg_t(window, 1:nCols, 1:kMax)    = qg(1:nCols, 1:kMax)
        ni_t(window, 1:nCols, 1:kMax)    = ni(1:nCols, 1:kMax)
        ns_t(window, 1:nCols, 1:kMax)    = ns(1:nCols, 1:kMax)
        nr_t(window, 1:nCols, 1:kMax)    = nr(1:nCols, 1:kMax)
        NG_t(window, 1:nCols, 1:kMax)    = NG(1:nCols, 1:kMax)
        NC_t(window, 1:nCols, 1:kMax)    = NC(1:nCols, 1:kMax)
        kzh_tke_t(window, 1:nCols, 1:kMax)   = kzh_tke(1:nCols, 1:kMax)
        omega_t(window, 1:nCols, 1:kMax) = omega(1:nCols, 1:kMax)
        LSRAIN_t(window, 1:nCols)        = LSRAIN_t(window - 1, 1:nCols)  ! repeat last input - only out vars
        LSSNOW_t(window, 1:nCols)        = LSSNOW_t(window - 1, 1:nCols)  ! repeat last input - only out vars

      endif
    endif  

  end subroutine fill_vars_timesteps
    ! TODO - descomente abaixo para apenas extrair os dados CSV da microfísica original
    !  ret = .false.


  SUBROUTINE RunMicro_HugMorr_NN( &
       nCols       , &!INTEGER      , INTENT(IN   ) :: nCols
       kMax        , &!INTEGER      , INTENT(IN   ) :: kMax 
       sl          , &
       tc          , &!REAL(KIND=r8), INTENT(INOUT) :: Tc (1:nCols, 1:kMax)
       QV          , &!REAL(KIND=r8), INTENT(INOUT) :: qv (1:nCols, 1:kMax)
       QC          , &!REAL(KIND=r8), INTENT(INOUT) :: qc (1:nCols, 1:kMax)
       QR          , &!REAL(KIND=r8), INTENT(INOUT) :: qr (1:nCols, 1:kMax)
       QI          , &!REAL(KIND=r8), INTENT(INOUT) :: qi (1:nCols, 1:kMax)
       QS          , &!REAL(KIND=r8), INTENT(INOUT) :: qs (1:nCols, 1:kMax)
       QG          , &!REAL(KIND=r8), INTENT(INOUT) :: qg (1:nCols, 1:kMax)
       NI          , &!REAL(KIND=r8), INTENT(INOUT) :: ni (1:nCols, 1:kMax)
       NS          , &!REAL(KIND=r8), INTENT(INOUT) :: ns (1:nCols, 1:kMax)
       NR          , &!REAL(KIND=r8), INTENT(INOUT) :: nr (1:nCols, 1:kMax)
       NG          , &!REAL(KIND=r8), INTENT(INOUT) :: NG (1:nCols, 1:kMax)
       NC          , &!REAL(KIND=r8), INTENT(INOUT) :: NC (1:nCols, 1:kMax)
       kzh_tke         , &!REAL(KIND=r8), INTENT(IN   ) :: kzh_tke (1:nCols, 1:kMax)
       omega       , &!REAL(KIND=r8), INTENT(IN   ) :: omega  ! omega (Pa/s)
       LSRAIN      , &!REAL(KIND=r8), INTENT(OUT) :: LSRAIN(1:nCols)
       LSSNOW        )!REAL(KIND=r8), INTENT(OUT) :: LSSNOW(1:nCols)
      !  EFFCS       , &!REAL(KIND=r8), INTENT(OUT  ) :: EFFCS (1:nCols, 1:kMax)   ! EFFCS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
      !  EFFIS       , &!REAL(KIND=r8), INTENT(OUT  ) :: EFFIS (1:nCols, 1:kMax)   ! EFFIS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)       

    IMPLICIT NONE
    INTEGER      , INTENT(IN   ) :: nCols
    INTEGER      , INTENT(IN   ) :: kMax
    REAL(KIND=r8), INTENT(IN   ) :: sl(kMax)
    REAL(KIND=r8), INTENT(INOUT) :: Tc(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: qv(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: qc(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: qr(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: qi(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: qs(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: qg(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: ni(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: ns(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: nr(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: NG(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: NC(1:nCols, 1:kMax)
    ! REAL(KIND=r8), INTENT(IN   ) :: tke(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(IN   ) :: kzh_tke(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(IN   ) :: omega (1:nCols, 1:kMax) ! omega (Pa/s)
    ! REAL(KIND=r8), INTENT(OUT  ) :: EFFCS (1:nCols, 1:kMax)   ! EFFCS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
    ! note: effis not currently passed out of microphysics (no coupling of ice eff rad with radiation)
    ! REAL(KIND=r8), INTENT(OUT  ) :: EFFIS (1:nCols, 1:kMax)   ! EFFIS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
    REAL(KIND=r8), INTENT(OUT  ) :: LSRAIN(1:nCols)
    REAL(KIND=r8), INTENT(OUT  ) :: LSSNOW(1:nCols)


    REAL(KIND=r8) :: start, finish, start_python, finish_python
    INTEGER :: k, in_unit, out_unit, dummy_k, out_nn_size, t_nn
    logical :: out_nn_exists

    integer,dimension(8) :: t ! arguments for date_and_time
    integer :: ms1,ms2  ! start and end times [ms]
    real :: dt                ! desired sleep interval [ms]

    call cpu_time(start)
    in_unit = 112233
    out_unit = 445566


    print*, "INPUTED =============================================================================================="
    
        ! data for python NN predict
    open(unit=in_unit,file='./BAM1D_HUMO_in.csv',status='unknown')
    write(in_unit, '(A)') 'k,sl,Tc,qv,qc,qr,qi,qs,qg,ni,ns,nr,NG,NC,tke,omega,LSRAIN,LSSNOW'
    do t_nn=1, window
      do k=1, kMax
        print*, k, sl_t(t_nn, k), Tc_t(t_nn, 1,k), qv_t(t_nn, 1,k), qc_t(t_nn, 1,k), qr_t(t_nn, 1,k), qi_t(t_nn, 1,k), qs_t(t_nn, 1,k), qg_t(t_nn, 1,k), &
        ni_t(t_nn, 1,k), ns_t(t_nn, 1,k), nr_t(t_nn, 1,k), NG_t(t_nn, 1,k), NC_t(t_nn, 1,k), kzh_tke_t(t_nn, 1,k), omega_t(t_nn, 1,k), LSRAIN_t(t_nn, 1), LSSNOW_t(t_nn, 1)

        call csv_write(in_unit, k, .false.)
        call csv_write(in_unit, sl_t(t_nn, k), .false.)
        call csv_write(in_unit, Tc_t(t_nn, 1,k), .false.)
        call csv_write(in_unit, qv_t(t_nn, 1,k), .false.)
        call csv_write(in_unit, qc_t(t_nn, 1,k), .false.)
        call csv_write(in_unit, qr_t(t_nn, 1,k), .false.)
        call csv_write(in_unit, qi_t(t_nn, 1,k), .false.)
        call csv_write(in_unit, qs_t(t_nn, 1,k), .false.)
        call csv_write(in_unit, qg_t(t_nn, 1,k), .false.)
        call csv_write(in_unit, ni_t(t_nn, 1,k), .false.)
        call csv_write(in_unit, ns_t(t_nn, 1,k), .false.)
        call csv_write(in_unit, nr_t(t_nn, 1,k), .false.)
        call csv_write(in_unit, NG_t(t_nn, 1,k), .false.)
        call csv_write(in_unit, NC_t(t_nn, 1,k), .false.)
        ! call csv_write(in_unit, tke_t(t_nn, 1,k), .false.)
        call csv_write(in_unit, kzh_tke_t(t_nn, 1,k), .false.)
        call csv_write(in_unit, omega_t(t_nn, 1,k), .false.)
        call csv_write(in_unit, LSRAIN_t(t_nn, 1), .false.)
        call csv_write(in_unit, LSSNOW_t(t_nn, 1), .true.)
      enddo
      
    enddo
    close(in_unit)
    
    call cpu_time(start_python)
    ! PYTHON READS INPUT AND WRITE OUTPUT
    
    ! waits a little more time to calculate the first out file size
    if (estimated_out_NN_file_size_real == 0) then
      call sleep(1)
      inquire( file = "./BAM1D_HUMN_out.csv", exist = out_nn_exists, size = estimated_out_NN_file_size )
      estimated_out_NN_file_size_real = estimated_out_NN_file_size - (estimated_out_NN_file_size/5)
    endif

    print*, "waiting python NN outs ..."
    dt=10
    do 
      ! call sleep(1)
      inquire( file = "./BAM1D_HUMN_out.csv", exist = out_nn_exists, size = out_nn_size )
      print*, "ESTIMATED SIZE = ", estimated_out_NN_file_size_real, " FILE SIZE= ", out_nn_size
      ! TODO
      ! print*, "out_nn_exists .and. out_nn_size = ", out_nn_exists, out_nn_size
      if ( out_nn_exists .and. out_nn_size >= estimated_out_NN_file_size_real ) then
        exit
      else
        call date_and_time(values=t)
        ms1=(t(5)*3600+t(6)*60+t(7))*1000+t(8)
        do ! check time:
          call date_and_time(values=t)
          ms2=(t(5)*3600+t(6)*60+t(7))*1000+t(8)
          if(ms2-ms1>=dt)exit
        enddo
        ! exit
      endif
    enddo  
    
    call cpu_time(finish_python)
    print*, 'Python Time = ',finish_python-start_python

    open(unit=out_unit,file='./BAM1D_HUMN_out.csv',status='old', action='read')
    print*, "  =============================================================================================="
    do k=1, kMax
      read(out_unit, *) Tc(1,k), qv(1,k) , qc(1,k), qr(1,k), qi(1,k), qs(1,k), qg(1,k), ni(1,k), ns(1,k), nr(1,k), NG(1,k), NC(1,k), LSRAIN(1), LSSNOW(1)
      print*, Tc(1,k), qv(1,k) , qc(1,k), qr(1,k), qi(1,k), qs(1,k), qg(1,k), ni(1,k), ns(1,k), nr(1,k), NG(1,k), NC(1,k), LSRAIN(1), LSSNOW(1)
    enddo
    
    close(out_unit, status='delete')

    
    call cpu_time(finish)
    print*, 'Total Micro_NN Time = ',finish-start


  END SUBROUTINE RunMicro_HugMorr_NN


END MODULE Micro_HugMorr_NN
