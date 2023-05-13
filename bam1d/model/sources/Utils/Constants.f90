!
!  $Author: pkubota $
!  $Date: 2008/09/23 17:51:54 $
!  $Revision: 1.9 $
!
MODULE Constants

  IMPLICIT NONE

  ! Selecting Kinds
  INTEGER, PARAMETER :: r4 = SELECTED_REAL_KIND(6)  ! Kind for 32-bits Real Numbers
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)   ! Kind for 32-bits Integer Numbers
  INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers
  INTEGER, PARAMETER :: r16 = SELECTED_REAL_KIND(15)! Kind for 128-bits Real Numbers

  ! Dimension for Diagnostics Arrays
  INTEGER , PARAMETER :: ndavl = 261 ! Number of Available Diagnostics Fields
  INTEGER , PARAMETER :: ndrq  = 209 ! Number of Desired Diagnostics Fields
  INTEGER , PARAMETER :: ncdg  =  20 ! Number of Combined Field Components
  INTEGER , PARAMETER :: jxavl = 160 ! Maximum Number of Available plus Combined Fields
  INTEGER , PARAMETER :: jxcdg =  50 ! Maximum Number of Combined Field Components
  REAL (KIND=r8), PARAMETER :: undef = 1.0e53_r8
  ! Dimension for Grid History Arrays
  INTEGER , PARAMETER :: ngfld=86 ! Number of Available Grid History Fields

  ! Dimension for Units Arrays
  INTEGER , PARAMETER :: nunits = 230        ! Number of Defined Units
  INTEGER , PARAMETER :: numx   = nunits-1   ! Number of Defined U nits minus 1
  INTEGER , PARAMETER :: ngrp   = nunits/10  ! Number of Units Groups
  INTEGER , PARAMETER :: ngrmx  = ngrp-1     ! Number of Units Groups minus 1
  INTEGER , PARAMETER :: ncf    = 1000       ! Number of First Unit Conversion Coefficents
  INTEGER , PARAMETER :: ncf2   = 500        ! Number of Second Unit Conversion Coefficents
  ! Dimension for Vegetation Arrays
  INTEGER , PARAMETER :: ityp = 13 ! Number of Vegetation Types
  INTEGER , PARAMETER :: imon = 12 ! Number of Months with Defined Vegetation Types
  INTEGER , PARAMETER :: icg  = 2  ! Number of Vegetation Parameters
  INTEGER , PARAMETER :: iwv  = 3  ! Number of Radiation Wavelengths
  INTEGER , PARAMETER :: idp  = 3  ! Number of Soil Layer Parameters
  INTEGER , PARAMETER :: ibd  = 2  ! Number of Vegetation Stage

  ! Constants
  REAL (KIND=r8), PARAMETER   :: pai   =        3.14159265358979_r8! constant pi=3.1415926e0
  REAL (KIND=r8), PARAMETER   :: pihalf=                 pai/2.0_r8! constant pi / 2
  REAL (KIND=r8), PARAMETER   :: tbar  =                   300.0_r8! 300 K
  REAL (KIND=r8), PARAMETER   :: pscons=                  1000.0_r8! A constant value for surface pressure
  REAL (KIND=r8), PARAMETER   :: cp    =                  1004.6_r8! specific heat of air           (j/kg/k)
  REAL (KIND=r8), PARAMETER   :: gasr  =                  287.05_r8! gas constant of dry air        (j/kg/k)
  REAL (KIND=r8), PARAMETER   :: er    =               6370000.0_r8! earth's radius                 (m)
  REAL (KIND=r8), PARAMETER   :: eriv  =                  1.0_r8/er! earth's radius inverse         (1/m)
  REAL (KIND=r8), PARAMETER   :: ersqiv=               1.0_r8/er**2! earth's radius square inverse  (1/m**2)
  REAL (KIND=r8), PARAMETER   :: grav  =                   9.8e0_r8! gravity constant               (m/s**2)
  REAL (KIND=r8), PARAMETER   :: ga2   =               grav/(er*er)!
  REAL (KIND=r8), PARAMETER   :: rk    =                    gasr/cp!
  REAL (KIND=r8), PARAMETER   :: raa   =                 gasr/er**2!
  REAL (KIND=r8), PARAMETER   :: twomg =             1.458492e-4_r8!
  REAL (KIND=r8), PARAMETER   :: tdelt =                   100.0_r8! equator to pole surface radiative equilibrium 
  ! temperature difference (k)
  REAL (KIND=r8), PARAMETER   :: tsfc0 =                   300.0_r8! equator surface radiative equilibrium temperature (k)
  REAL (KIND=r8), PARAMETER   :: h0    =                   8.2e3_r8! h0 scale height of radiative equilibrium temperature
  ! assuming isothermal atmosphere (m)
  REAL (KIND=r8), PARAMETER   :: rlaps =                  6.5e-3_r8! rlaps radiative equilibrium temperature lapse rate (km)
  ! from surface to stratosphere
  REAL (KIND=r8), PARAMETER   :: tstrat=                   200.0_r8! tstrat    stratospheric radiative equilibrium 
  ! temperature (k)
  REAL (KIND=r8), PARAMETER   :: qmin  =                 1.0e-12_r8! minimum value relative humidity relhum
  ! constant qmin = 1.0e-12 
  REAL (KIND=r8)              :: root2                             ! SQRT(2.0)
  REAL (KIND=r8)              :: coriol                            !
  REAL (KIND=r8), PARAMETER   :: pie  =              3.1415926e0_r8! constant pi=3.1415926e0   
  REAL (KIND=r8), PARAMETER   :: pai12  =              pie/12.0_r8 ! pi/12
  REAL (KIND=r8), PARAMETER   :: hl  =                    2.52e6_r8! heat of evaporation of water     (j/kg) 
  REAL (KIND=r8), PARAMETER   :: stefan =                5.67e-8_r8! stefan Stefan Boltzman constant 
  REAL (KIND=r8), PARAMETER   :: solcon =              1365.00e0_r8! solar constant (wgne value)    (w/m**2) 
  REAL (KIND=r8), PARAMETER   :: rmwmd  =                0.622e0_r8! fracao molar entre a agua e o ar 
  REAL (KIND=r8), PARAMETER   :: rmwmdi =                 1.61e0_r8!
  REAL (KIND=r8), PARAMETER   :: e0c  =                   6.11e0_r8!
  REAL (KIND=r8), PARAMETER   :: t000   =                 299.e0_r8! Temp. global de ref. ao nivel do mar. (K)
  REAL (KIND=r8), PARAMETER   :: p000   =                1013.e0_r8! Pressao de referencia ao nivel do mar (mb)
  REAL (KIND=r8), PARAMETER   :: p00  =                  1000.e0_r8! Pressao de referencia ao nivel do mar (mb)
  REAL (KIND=r8), PARAMETER   :: snomel= 333624.2e0_r8 *1000.0e0_r8! Calor latente de fusao is expressed in (j m-1)
  REAL (KIND=r8), PARAMETER   :: tf  =                  273.16e0_r8! Temperatura de congelamento (K)=273.16e0
  REAL (KIND=r8), PARAMETER   :: epsfac =                0.622e0_r8! Constante 0.622 Razao entre as massas 
  !  moleculares do vapor e do ar seco
  REAL (KIND=r8), PARAMETER   :: athird =         1.0e0_r8/3.0e0_r8!
  REAL (KIND=r8), PARAMETER   :: clai   =   4.2e0_r8 *1000.0e0_r8 *0.2e0_r8! heat capacity of foliage    
  REAL (KIND=r8), PARAMETER   :: cw  =4.2e0_r8 *1000.0e0_r8 *1000.0e0_r8! liquid water heat capacity     (j/m**3)
  REAL (KIND=r8), PARAMETER   :: delq   =       0.608e0_r8! constant delq = 0.608e0
  REAL (KIND=r8), PARAMETER   :: tice   =      271.16e0_r8! constant tice
  REAL (KIND=r8), PARAMETER   :: icealv =         0.8e0_r8! constant icealv
  REAL (KIND=r8), PARAMETER   :: icealn =         0.4e0_r8! constant icealn
  REAL (KIND=r8), PARAMETER   :: oceald =      0.0419e0_r8! constant oceald
  REAL (KIND=r8), PARAMETER   :: z0ice  =       0.001e0_r8! 
  REAL (KIND=r8), PARAMETER   :: con_rd =    2.8705e+2_r8 ! gas constant air (J/kg/K)
  REAL (KIND=r8), PARAMETER   :: con_rv =    4.6150e+2_r8 ! gas constant H2O (J/kg/K)
  REAL (KIND=r8), PARAMETER   :: EPS    =    con_rd/con_rv        ! #
  REAL (KIND=r8), PARAMETER   :: EPSM1  =    con_rd/con_rv-1.0_r8 ! #

  ! Dobson Unit (DU): It is assumed, that the ozone of the column in pure
  ! form is present at normal pressure (1013.15 hPa) and at 0deg C (273.15 K).
  ! The thickness of this layer is than given in "milli-centimetre" (10-5 m).
  ! For example 300 DU correspond to a layer of 3 mm thickness under the
  ! conditions given above.
  ! Conversion:
  !      dp_ozone = -rho_ozone g dz_ozone
  !      dz_ozone =  -dp_ozone/g/rho_ozone
  !                          dp_air    R_ozone * To
  !               = -o3mix * ------ * -------------- = -o3mix * dp_air/g * gm2dob
  !                            g            Po
  !
  ! But R_ozone = R*/M_ozone = 8.316963 J/mol/K /(47.998 g/mol) = 173.28 W/kg/K
  ! Therefore, the conversion constant (kg/m2 to dobson):
  !                gm2dob = 46716.3694_r8
  !
  REAL(KIND=r8)   , PARAMETER :: gm2dob = 46716.3694_r8 ! hmjb
  ! 
  REAL (KIND=r8), ALLOCATABLE :: tdampr(:)                        !
  REAL (KIND=r8), ALLOCATABLE :: tdampf(:)                        !
  REAL (KIND=r8), ALLOCATABLE :: tov(:)                           ! 300 k at all levels (substituir tov por tbar)  

CONTAINS

  SUBROUTINE InitConstants (kMax) 

    IMPLICIT NONE

    INTEGER , INTENT(IN) :: kMax

    coriol = SQRT(2.0_r8/3.0_r8)*twomg
    root2  = SQRT(2.0_r8)
    ALLOCATE(tov(kMax))
    tov = tbar
    !tov = tbar + 6.0_r8
    ALLOCATE(tdampr(kMax))
    tdampr=10.0_r8
    ALLOCATE(tdampf(kMax))
    tdampf(1)=2.0_r8; tdampf(2)=5.0_r8; tdampf(3:kmax)=10.0_r8

  END SUBROUTINE InitConstants
END MODULE Constants
